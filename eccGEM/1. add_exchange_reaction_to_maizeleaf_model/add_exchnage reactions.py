import cobra
from cobra import Reaction, Metabolite
import pandas as pd

def add_metabolite_if_missing(model, met_id, compartment):
    """Add metabolite to model if it doesn't exist"""
    try:
        return model.metabolites.get_by_id(met_id)
    except KeyError:
        # Create new metabolite
        metabolite = Metabolite(
            id=met_id,
            name=f"Metabolite {met_id.split('_')[0]}",
            compartment=compartment
        )
        model.add_metabolites([metabolite])
        print(f"Added missing metabolite: {met_id}")
        return metabolite

def add_exchange_and_transport_reactions(model, metabolite_ids):
    added_reactions = []
    reaction_counter = 3
    
    for met_id in metabolite_ids:
        extracellular_met_id = met_id  # "C00267_f"
        cytoplasmic_met_id = met_id.replace('_f', '_c')  # "C00267_c"
        
        # Add metabolites if they don't exist
        extracellular_met = add_metabolite_if_missing(model, extracellular_met_id, 'f')
        cytoplasmic_met = add_metabolite_if_missing(model, cytoplasmic_met_id, 'c')
        
        print(f"Processing: {extracellular_met.id} - {extracellular_met.name}")
            
        code = met_id.split('_')[0]
        exchange_id = f"R_ExBe{reaction_counter}"
        transport_id = f"R_TransBc{reaction_counter}"
        
        # Add exchange reaction with specified bounds
        if exchange_id not in model.reactions:
            exchange_reaction = Reaction(exchange_id)
            exchange_reaction.name = f"BSC exchange of {code}"
            exchange_reaction.lower_bound = -10.0  # Set to -10
            exchange_reaction.upper_bound = 100.0  # Set to 100
            
            exchange_reaction.add_metabolites({extracellular_met: -1.0})
            
            model.add_reactions([exchange_reaction])
            added_reactions.append(exchange_reaction.id)
            print(f"Added exchange: {exchange_id}")
        
        # Add transport reaction
        if transport_id not in model.reactions:
            transport_reaction = Reaction(transport_id)
            transport_reaction.name = f"BSC Transport of {code} in out of cell"
            transport_reaction.lower_bound = -1000.0
            transport_reaction.upper_bound = 1000.0
            
            transport_reaction.add_metabolites({
                extracellular_met: -1.0,
                cytoplasmic_met: 1.0
            })
            
            model.add_reactions([transport_reaction])
            added_reactions.append(transport_reaction.id)
            print(f"Added transport: {transport_id}")
        
        reaction_counter += 1
    
    return added_reactions

def get_all_exchange_reactions(model, solution=None):
    exchange_data = []
    
    for reaction in model.reactions:
        if (len(reaction.metabolites) == 1 and list(reaction.metabolites.values())[0] < 0) or "ExBe" in reaction.id:
            metabolite = list(reaction.metabolites.keys())[0]
            flux = solution.fluxes[reaction.id] if solution and solution.status == 'optimal' else 0
            
            exchange_data.append({
                'Reaction_ID': reaction.id,
                'Reaction_Name': reaction.name,
                'Type': 'Exchange',
                'Metabolite_ID': metabolite.id,
                'Metabolite_Name': metabolite.name,
                'Formula': metabolite.formula,
                'Compartment': metabolite.compartment,
                'Lower_Bound': reaction.lower_bound,
                'Upper_Bound': reaction.upper_bound,
                'Flux_Value': flux,
                'Reaction_String': reaction.reaction
            })
        
        elif "TransBc" in reaction.id and len(reaction.metabolites) == 2:
            reactants = [met for met, coef in reaction.metabolites.items() if coef < 0]
            products = [met for met, coef in reaction.metabolites.items() if coef > 0]
            
            if reactants and products:
                flux = solution.fluxes[reaction.id] if solution and solution.status == 'optimal' else 0
                
                exchange_data.append({
                    'Reaction_ID': reaction.id,
                    'Reaction_Name': reaction.name,
                    'Type': 'Transport',
                    'Metabolite_ID': reactants[0].id,
                    'Metabolite_Name': reactants[0].name,
                    'Formula': reactants[0].formula,
                    'Compartment': f"{reactants[0].compartment}→{products[0].compartment}",
                    'Lower_Bound': reaction.lower_bound,
                    'Upper_Bound': reaction.upper_bound,
                    'Flux_Value': flux,
                    'Reaction_String': reaction.reaction
                })
    
    return exchange_data

def main():
    model = cobra.io.read_sbml_model("izma6518_L.xml")
    print(f"Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")
    
    metabolite_ids = ["C00267_f", "C00022_f", "C00152_f", "C00049_f", "C00095_f", "C00124_f", "C00159_f"]
    
    initial_count = len(model.reactions)
    initial_met_count = len(model.metabolites)
    
    added_reactions = add_exchange_and_transport_reactions(model, metabolite_ids)
    
    final_count = len(model.reactions)
    final_met_count = len(model.metabolites)
    
    print(f"Added {len(added_reactions)} reactions. Total: {initial_count} -> {final_count}")
    print(f"Metabolites: {initial_met_count} -> {final_met_count}")
    
    solution = model.optimize()
    print(f"FBA Status: {solution.status}, Objective: {solution.objective_value:.6f}")
    
    exchange_data = get_all_exchange_reactions(model, solution)
    
    df = pd.DataFrame(exchange_data)
    df.to_excel("exchange_transport_reactions.xlsx", index=False)
    print(f"Excel file created with {len(exchange_data)} reactions")
    
    cobra.io.write_sbml_model(model, "izma6518_L_exchange_added.xml")
    print(f"Model saved: izma6518_L_exchange_added.xml")

if __name__ == "__main__":
    main()