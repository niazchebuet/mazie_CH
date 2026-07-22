import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction
import pandas as pd
import os

def is_exchange_reaction(rxn):
    """Check if a reaction is an exchange reaction"""
    # Check by ID patterns
    if any(prefix in rxn.id.upper() for prefix in ['EX_', 'DM_', 'SINK_', 'DEMAND_']):
        return True
    
    # Check if reaction has only one metabolite (typical for exchange reactions)
    if len(rxn.metabolites) == 1:
        return True
    
    # Check if all metabolites are in extracellular compartment
    if all('_e' in met.id or '_e0' in met.id or 'extracellular' in met.compartment.lower() 
          for met in rxn.metabolites):
        return True
    
    return False

def split_reversible_reactions(model):
    model_copy = model.copy()
    reactions_to_remove = []
    reactions_to_add = []
    
    # Different bounds for regular vs exchange reactions
    regular_ub = 1000.0
    exchange_ub = 100.0
    default_lb = 0.0
    
    for rxn in model_copy.reactions:
        # Determine appropriate upper bound
        ub = exchange_ub if is_exchange_reaction(rxn) else regular_ub
        
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            # Forward reaction
            fwd = Reaction(f"{rxn.id}_FWD")
            fwd.name = f"{rxn.name} (fwd)"
            fwd.lower_bound, fwd.upper_bound = default_lb, ub
            fwd.gene_reaction_rule = rxn.gene_reaction_rule
            fwd.subsystem = rxn.subsystem
            fwd.add_metabolites(rxn.metabolites)
            
            # Reverse reaction
            rev = Reaction(f"{rxn.id}_REV")
            rev.name = f"{rxn.name} (rev)"
            rev.lower_bound, rev.upper_bound = default_lb, ub
            rev.gene_reaction_rule = rxn.gene_reaction_rule
            rev.subsystem = rxn.subsystem
            rev.add_metabolites({met: -coeff for met, coeff in rxn.metabolites.items()})
            
            reactions_to_add.extend([fwd, rev])
            reactions_to_remove.append(rxn)
            
        elif rxn.lower_bound < 0 and rxn.upper_bound <= 0:
            # Convert reverse-only to forward
            new_rxn = Reaction(f"{rxn.id}_FWD")
            new_rxn.name = f"{rxn.name} (conv)"
            new_rxn.lower_bound, new_rxn.upper_bound = default_lb, ub
            new_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            new_rxn.subsystem = rxn.subsystem
            new_rxn.add_metabolites({met: -coeff for met, coeff in rxn.metabolites.items()})
            
            reactions_to_add.append(new_rxn)
            reactions_to_remove.append(rxn)
        
        elif rxn.lower_bound >= 0 and rxn.upper_bound > 0:
            # Already irreversible forward reaction - just update bounds
            rxn.lower_bound = default_lb
            rxn.upper_bound = ub
    
    model_copy.remove_reactions(reactions_to_remove)
    model_copy.add_reactions(reactions_to_add)
    
    # Final pass to ensure all reactions have correct bounds
    for rxn in model_copy.reactions:
        if rxn.lower_bound < 0:
            rxn.lower_bound = default_lb
        
        # Set appropriate upper bound based on reaction type
        target_ub = exchange_ub if is_exchange_reaction(rxn) else regular_ub
        if rxn.upper_bound != target_ub:
            rxn.upper_bound = target_ub
    
    return model_copy

def process_model(sbml_file):
    # Load and split
    model = read_sbml_model(sbml_file)
    split_model = split_reversible_reactions(model)
    
    # Set objective (find biomass or use current)
    biomass = next((r for r in model.reactions if r.objective_coefficient != 0), None)
    if not biomass:
        biomass = next((r for r in model.reactions if 'biomass' in r.id.lower()), None)
    if biomass:
        fwd_id = f"{biomass.id}_FWD"
        if fwd_id in [r.id for r in split_model.reactions]:
            split_model.objective = fwd_id
        else:
            split_model.objective = biomass.id
    
    # Run analyses
    fba = split_model.optimize()
    pfba = cobra.flux_analysis.pfba(split_model)
    
    # Create Excel data
    data = []
    for r in split_model.reactions:
        is_exchange = is_exchange_reaction(r)
        data.append({
            'Reaction_ID': r.id,
            'Name': r.name,
            'Subsystem': r.subsystem,
            'Genes': r.gene_reaction_rule,
            'Equation': r.reaction,
            'Lower_Bound': r.lower_bound,
            'Upper_Bound': r.upper_bound,
            'Is_Exchange': is_exchange,
            'FBA_Flux': fba.fluxes[r.id],
            'pFBA_Flux': pfba.fluxes[r.id],
            'Active': abs(pfba.fluxes[r.id]) > 1e-6
        })
    
    df = pd.DataFrame(data)
    
    # Save files
    base = os.path.splitext(sbml_file)[0]
    write_sbml_model(split_model, f"irreversible_{base}.xml")
    
    with pd.ExcelWriter(f"flux_{base}.xlsx") as writer:
        df.to_excel(writer, sheet_name='All_Fluxes', index=False)
        df[df['Active']].to_excel(writer, sheet_name='Active_Only', index=False)
        
        # Exchange reactions summary
        exchange_rxns = df[df['Is_Exchange'] == True]
        if not exchange_rxns.empty:
            exchange_rxns.to_excel(writer, sheet_name='Exchange_Reactions', index=False)
        
        # Regular reactions summary
        regular_rxns = df[df['Is_Exchange'] == False]
        if not regular_rxns.empty:
            regular_rxns.to_excel(writer, sheet_name='Regular_Reactions', index=False)
        
        pd.DataFrame({
            'Metric': ['Growth_Rate_FBA', 'Growth_Rate_pFBA', 'Total_Reactions', 'Active_Reactions', 
                      'Exchange_Reactions', 'Regular_Reactions'],
            'Value': [fba.objective_value, pfba.objective_value, len(split_model.reactions), 
                     df['Active'].sum(), len(exchange_rxns), len(regular_rxns)]
        }).to_excel(writer, sheet_name='Summary', index=False)
    
    print(f"Growth: {fba.objective_value:.3f} | Active: {df['Active'].sum()}/{len(split_model.reactions)}")
    print(f"Regular reactions upper bound: 1000")
    print(f"Exchange reactions upper bound: 100")
    print(f"Exchange reactions found: {len(df[df['Is_Exchange'] == True])}")

# Usage
process_model("cochliobolus_maize_leaf_community.xml")  # Change filename here