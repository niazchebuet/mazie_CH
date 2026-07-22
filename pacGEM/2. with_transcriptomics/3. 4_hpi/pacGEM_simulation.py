import cobra
import pandas as pd
from optlang import Objective
# Load the community model
model = cobra.io.read_sbml_model("4_hpi.xml")
print("Model loaded successfully")

# Set solver to glpk or other desired commercial solver
model.solver = 'glpk'

# Load protein cost data
ch_data = pd.read_excel("community_model_per_unit_protein_data_CH.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])
zma_data = pd.read_excel("community_model_per_unit_protein_data_ZMA.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])
# Create protein cost dictionaries
model_reaction_ids = set(r.id for r in model.reactions)
ch_protein_costs = {}
for _, row in ch_data.iterrows():
    reaction_id = row["Reactions"]
    if reaction_id in model_reaction_ids:
        ch_protein_costs[reaction_id] = row["Protein_cost"]
zma_protein_costs = {}
for _, row in zma_data.iterrows():
    reaction_id = row["Reactions"]
    if reaction_id in model_reaction_ids:
        zma_protein_costs[reaction_id] = row["Protein_cost"]
print("Protein costs loaded")
# Define target reactions (excluding GROWTH)
other_target_reactions = [
    "MR00338_LSQBKT_d_RSQBKT_", "R02446_LSQBKT_p_RSQBKT_", "R02446_LSQBKT_q_RSQBKT_",
    "R01121_LSQBKT_c_RSQBKT_", "R01121_LSQBKT_d_RSQBKT_", "R05611_LSQBKT_p_RSQBKT__FWD",
    "R05884_LSQBKT_q_RSQBKT_", "R00289_LSQBKT_c_RSQBKT__REV", "R00289_LSQBKT_d_RSQBKT__REV",
    "R05884_LSQBKT_p_RSQBKT_", "R05611_LSQBKT_q_RSQBKT__REV", "R09248_LSQBKT_q_RSQBKT_",
    "R00086_LSQBKT_m_RSQBKT_", "R00086_LSQBKT_t_RSQBKT_", "R00086_LSQBKT_n_RSQBKT_",
    "R00086_LSQBKT_u_RSQBKT_", "R00086_LSQBKT_i_RSQBKT__REV", "R00086_LSQBKT_j_RSQBKT__REV"
]
# Target reaction objects
other_target_reaction_objects = []
for target in other_target_reactions:
    if target in model.reactions:
        other_target_reaction_objects.append(model.reactions.get_by_id(target))
# Validate GROWTH reaction
growth_reaction = None
if "GROWTH" in model.reactions:
    growth_reaction = model.reactions.get_by_id("GROWTH")
else:
    print("GROWTH reaction not found!")
    exit()
# Define exchange reactions (replace with full list as needed)
ch_exchange_ids = [
    "BOUNDARY_ACe", "BOUNDARY_CO2e", "BOUNDARY_DGLCe", "BOUNDARY_He",
    "BOUNDARY_H2O2e", "BOUNDARY_LACe", "BOUNDARY_GLYCOe", "BOUNDARY_CELLOBe",
    "BOUNDARY_O2e_FWD", "BOUNDARY_O2e_REV"
]
zma_exchange_ids = ["ExBe1_FWD", "ExBe1_REV", "ExBe10_FWD", "ExBe10_REV", "ExBe12_FWD", "ExBe12_REV", "ExBe13_FWD", "ExBe13_REV"]
exchange_reactions = []
ch_exchange_set, zma_exchange_set = set(ch_exchange_ids), set(zma_exchange_ids)
for reaction in model.reactions:
    if reaction.id in ch_exchange_set or reaction.id in zma_exchange_set:
        exchange_reactions.append(reaction)
if "ExBe15_REV" in model.reactions:
    model.reactions.get_by_id("ExBe15_REV").bounds = (10, model.reactions.get_by_id("ExBe15_REV").bounds[1])
print("Setting custom proteome constraints...")
# --- Main protein constraints for CH ---
ch_variables = {}
for reaction_id, cost in ch_protein_costs.items():
    reaction = model.reactions.get_by_id(reaction_id)
    ch_variables[reaction.forward_variable] = cost
# Linear expressions for ribo/other in terms of GROWTH variable (units [g/gDW])
ribo_expr = 0.084 + 0.226 * growth_reaction.forward_variable    # Example
# Instead of a single constraint, use two linear constraints for p_other (piecewise-line plateau)
met_expr = sum(coeff * var for var, coeff in ch_variables.items()) / 1000
from optlang.symbolics import Zero
p_other = model.problem.Variable("p_other", lb=0)
model.add_cons_vars(p_other)
# 1. Linear decay upper bound
linear_expr = 0.516 - 1.181 * growth_reaction.forward_variable
model.add_cons_vars(
    model.problem.Constraint(
        p_other - linear_expr, ub=0, name="CH_p_other_linear"
    )
)
# 2. Plateau minimum (lower bound)
plateau_min = 0.10
model.add_cons_vars(
    model.problem.Constraint(
        p_other, lb=plateau_min, name="CH_p_other_plateau"
    )
)
# 3. Total proteome constraint
proteome_constraint = model.problem.Constraint(
    met_expr + ribo_expr + p_other,
    ub=0.75,
    name="CH_protein_capacity"
)
model.add_cons_vars(proteome_constraint)
# ZMA protein constraint as before
chunk_size = 1000
reactions_list = list(zma_protein_costs.items())
num_chunks = (len(reactions_list) + chunk_size - 1) // chunk_size
chunk_aux_vars = []
for i in range(num_chunks):
    chunk_start, chunk_end = i * chunk_size, min((i + 1) * chunk_size, len(reactions_list))
    chunk = reactions_list[chunk_start:chunk_end]
    aux_var_name = f"zma_protein_chunk_{i+1}"
    chunk_aux_var = model.problem.Variable(aux_var_name, lb=0)
    model.add_cons_vars(chunk_aux_var)
    chunk_aux_vars.append(chunk_aux_var)
    chunk_expr = sum(
        model.reactions.get_by_id(reaction_id).forward_variable * cost
        for reaction_id, cost in chunk
    )
    chunk_constraint = model.problem.Constraint(
        chunk_aux_var - chunk_expr, lb=0, ub=0, name=f"zma_chunk_{i+1}_balance"
    )
    model.add_cons_vars(chunk_constraint)
total_zma_constraint = model.problem.Constraint(
    sum(chunk_aux_vars), lb=0, ub=600, name="ZMA_total_protein_capacity"
)
model.add_cons_vars(total_zma_constraint)
print("Protein capacity constraints applied")
# === Optimization steps ===
print("Step 1: Maximizing cochliobolus biomass growth rate")
model.objective = Objective(growth_reaction.forward_variable, direction='max')
solution1 = model.optimize()
if solution1.status == "optimal":
    max_growth = growth_reaction.flux
    # Fix growth at optimal
    growth_constraint = model.problem.Constraint(
        growth_reaction.forward_variable,
        lb=max_growth, ub=max_growth,
        name="fixed_growth"
    )
    model.add_cons_vars(growth_constraint)
else:
    print("Step 1 optimization failed")
    exit()
print("Step 2: Maximizing maize target reactions")
if other_target_reaction_objects:
    objective_expr = sum(reaction.forward_variable for reaction in other_target_reaction_objects)
    model.objective = Objective(objective_expr, direction='max')
    solution2 = model.optimize()
    if solution2.status == "optimal":
        max_other_targets_sum = solution2.objective_value
        print(f"Maize target objective value: {max_other_targets_sum:.5f}")
        other_targets_constraint = model.problem.Constraint(
            sum(reaction.forward_variable for reaction in other_target_reaction_objects),
            lb=max_other_targets_sum, ub=max_other_targets_sum,
            name="fixed_other_targets"
        )
        model.add_cons_vars(other_targets_constraint)
    else:
        print("Step 2 optimization failed")
        exit()
else:
    max_other_targets_sum = 0
    print(f"Maize target objective value: {max_other_targets_sum:.5f}")
print("Step 3: Minimizing total flux of exchange reactions for cochliobolus and maize")
if exchange_reactions:
    objective_expr = sum(reaction.forward_variable for reaction in exchange_reactions)
    model.objective = Objective(objective_expr, direction='min')
    solution3 = model.optimize()
    if solution3.status == "optimal":
        min_exchange_sum = solution3.objective_value
    else:
        print("Step 3 optimization failed")
        exit()
else:
    min_exchange_sum = 0
# === Protein pools reporting ===
ch_growth = growth_reaction.flux
p_met = sum(model.reactions.get_by_id(rxn).flux * cost for rxn, cost in ch_protein_costs.items()) / 1000   # [g/gDW]
p_ribo = 0.084 + 0.226 * ch_growth     # [g/gDW]
val_linear = 0.516 - 1.181 * ch_growth
p_other_val = plateau_min if val_linear < plateau_min else val_linear
zma_usage = sum(model.reactions.get_by_id(rxn).flux * cost for rxn, cost in zma_protein_costs.items()) / 1000   # [g/gDW]
print("\n-- PROTEIN REQUIREMENTS --")
print(f"Cochliobolus growth rate: {ch_growth:.3f} per hour")
print(f"Metabolic protein requirement: {p_met:.3f} g/gDW")
print(f"Ribosomal protein requirement: {p_ribo:.3f} g/gDW")
print(f"Other protein requirement: {p_other_val:.3f} g/gDW")
print(f"Total used: {(p_met+p_ribo+p_other_val):.3f} g/gDW (should be <= 0.75)")
print(f"ZMA protein usage: {zma_usage:.3f} g/gDW")
# === Export to Excel ===
with pd.ExcelWriter("CH_Results.xlsx") as writer:
    ch_flux_data = []
    for reaction_id in ch_protein_costs.keys():
        if reaction_id in model.reactions:
            flux = model.reactions.get_by_id(reaction_id).flux
            ch_flux_data.append({
                "Reaction_ID": reaction_id,
                "Flux (mmol/gDW/hr)": flux
            })
    ch_flux_df = pd.DataFrame(ch_flux_data)
    ch_flux_df.to_excel(writer, sheet_name='Flux_Distribution', index=False)
    ch_protein_data = []
    for reaction_id, cost in ch_protein_costs.items():
        flux = model.reactions.get_by_id(reaction_id).flux
        protein_required = (cost * flux) / 1000  # g/gDW
        ch_protein_data.append({
            "Reaction_ID": reaction_id,
            "Protein_Required (g/gDW)": protein_required
        })
    total_ch_protein = p_met
    ch_protein_data.append({
        "Reaction_ID": "TOTAL (metabolic only)",
        "Protein_Required (g/gDW)": total_ch_protein
    })
    ch_protein_df = pd.DataFrame(ch_protein_data)
    ch_protein_df.to_excel(writer, sheet_name='Protein_Requirements', index=False)
    summary_df = pd.DataFrame([
        {"Pool": "Metabolic", "g/gDW": p_met},
        {"Pool": "Ribosomal", "g/gDW": p_ribo},
        {"Pool": "Other", "g/gDW": p_other_val},
        {"Pool": "Total", "g/gDW": p_met + p_ribo + p_other_val}
    ])
    summary_df.to_excel(writer, sheet_name='Protein_Pools', index=False)
with pd.ExcelWriter("ZMA_Results.xlsx") as writer:
    zma_flux_data = []
    for reaction_id in zma_protein_costs.keys():
        if reaction_id in model.reactions:
            flux = model.reactions.get_by_id(reaction_id).flux
            zma_flux_data.append({
                "Reaction_ID": reaction_id,
                "Flux (mmol/gDW/hr)": flux
            })
    zma_flux_df = pd.DataFrame(zma_flux_data)
    zma_flux_df.to_excel(writer, sheet_name='Flux_Distribution', index=False)
    zma_protein_data = []
    for reaction_id, cost in zma_protein_costs.items():
        flux = model.reactions.get_by_id(reaction_id).flux
        protein_required = (cost * flux) / 1000  # g/gDW
        zma_protein_data.append({
            "Reaction_ID": reaction_id,
            "Protein_Required (g/gDW)": protein_required
        })
    total_zma_protein = zma_usage
    zma_protein_data.append({
        "Reaction_ID": "TOTAL",
        "Protein_Required (g/gDW)": total_zma_protein
    })
    zma_protein_df = pd.DataFrame(zma_protein_data)
    zma_protein_df.to_excel(writer, sheet_name='Protein_Requirements', index=False)
print("Results exported to CH_Results.xlsx and ZMA_Results.xlsx")