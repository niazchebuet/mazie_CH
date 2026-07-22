import cobra
import pandas as pd
from optlang import Objective
# Load the community model
model = cobra.io.read_sbml_model("irreversible_cochliobolus_maize_leaf_community.xml")
print("Model loaded successfully")
# Load protein cost data
ch_data = pd.read_excel("community_model_per_unit_protein_data_CH.xlsx",
                       usecols=[0, 1], names=["Reactions", "Protein_cost"])
zma_data = pd.read_excel("community_model_per_unit_protein_data_ZMA.xlsx",
                        usecols=[0, 1], names=["Reactions", "Protein_cost"])
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
# Validate target reactions
valid_other_targets = []
other_target_reaction_objects = []
for target in other_target_reactions:
    if target in model.reactions:
        valid_other_targets.append(target)
        other_target_reaction_objects.append(model.reactions.get_by_id(target))
# Validate GROWTH reaction
growth_reaction = None
if "GROWTH" in model.reactions:
    growth_reaction = model.reactions.get_by_id("GROWTH")
# Define exchange reactions
ch_exchange_ids = [
    "BOUNDARY_ACe", "BOUNDARY_CO2e", "BOUNDARY_DGLCe", "BOUNDARY_He",
    "BOUNDARY_H2O2e", "BOUNDARY_LACe", "BOUNDARY_GLYCOe", "BOUNDARY_CELLOBe",
    "BOUNDARY_O2e_FWD", "BOUNDARY_O2e_REV"
]
zma_exchange_ids = [
    "ExB_C00719", "ExBe23", "ExBe24", "ExBe25", "ExBe27", "ExBe31", "ExBe33",
    "ExBe35", "ExBe36", "ExBe37", "ExBe38", "ExBe42", "ExBe44", "ExBe46",
    "ExBe48a", "ExBe48b", "ExBe49", "ExBe5", "ExBe50", "ExBe54", "ExBe55",
    "ExBe58", "ExBe8", "ExLB_C00001", "ExLB_C00002", "ExLB_C00008", "ExLB_C00009",
    "ExLB_C00014", "ExLB_C00022", "ExLB_C00025", "ExLB_C00026", "ExLB_C00031",
    "ExLB_C00036", "ExLB_C00037", "ExLB_C00041", "ExLB_C00042", "ExLB_C00044",
    "ExLB_C00047", "ExLB_C00048", "ExLB_C00049", "ExLB_C00062", "ExLB_C00063",
    "ExLB_C00064", "ExLB_C00065", "ExLB_C00073", "ExLB_C00075", "ExLB_C00078",
    "ExLB_C00079", "ExLB_C00082", "ExLB_C00085", "ExLB_C00089", "ExLB_C00092",
    "ExLB_C00095", "ExLB_C00099", "ExLB_C00116", "ExLB_C00117", "ExLB_C00121",
    "ExLB_C00122", "ExLB_C00123", "ExLB_C00124", "ExLB_C00131", "ExLB_C00135",
    "ExLB_C00137", "ExLB_C00148", "ExLB_C00149", "ExLB_C00152", "ExLB_C00158",
    "ExLB_C00159", "ExLB_C00160", "ExLB_C00181", "ExLB_C00183", "ExLB_C00185",
    "ExLB_C00187", "ExLB_C00188", "ExLB_C00189", "ExLB_C00191", "ExLB_C00199",
    "ExLB_C00208", "ExLB_C00209", "ExLB_C00216", "ExLB_C00238", "ExLB_C00244",
    "ExLB_C00249", "ExLB_C00257", "ExLB_C00258", "ExLB_C00259", "ExLB_C00286",
    "ExLB_C00315", "ExLB_C00333", "ExLB_C00334", "ExLB_C00407", "ExLB_C00417",
    "ExLB_C00458", "ExLB_C00459", "ExLB_C00490", "ExLB_C00491", "ExLB_C00492",
    "ExLB_C00493", "ExLB_C00507", "ExLB_C00532", "ExLB_C00590", "ExLB_C00617",
    "ExLB_C00712", "ExLB_C00740", "ExLB_C00794", "ExLB_C00805", "ExLB_C00811",
    "ExLB_C00956", "ExLB_C01083", "ExLB_C01235", "ExLB_C01494", "ExLB_C01530",
    "ExLB_C01595", "ExLB_C01793", "ExLB_C02226", "ExLB_C02325", "ExLB_C02477",
    "ExLB_C02591", "ExLB_C02612", "ExLB_C02646", "ExLB_C03284", "ExLB_C06427",
    "ExLB_C08362", "ExLB_C15804", "ExLB_C15805", "ExLB_C15806", "ExLM_C00001",
    "ExLM_C00002", "ExLM_C00008", "ExLM_C00009", "ExLM_C00014", "ExLM_C00022",
    "ExLM_C00025", "ExLM_C00026", "ExLM_C00031", "ExLM_C00036", "ExLM_C00037",
    "ExLM_C00041", "ExLM_C00042", "ExLM_C00044", "ExLM_C00047", "ExLM_C00048",
    "ExLM_C00049", "ExLM_C00062", "ExLM_C00063", "ExLM_C00064", "ExLM_C00065",
    "ExLM_C00073", "ExLM_C00075", "ExLM_C00078", "ExLM_C00079", "ExLM_C00082",
    "ExLM_C00085", "ExLM_C00089", "ExLM_C00092", "ExLM_C00095", "ExLM_C00099",
    "ExLM_C00116", "ExLM_C00117", "ExLM_C00121", "ExLM_C00122", "ExLM_C00123",
    "ExLM_C00124", "ExLM_C00131", "ExLM_C00135", "ExLM_C00137", "ExLM_C00148",
    "ExLM_C00149", "ExLM_C00152", "ExLM_C00158", "ExLM_C00159", "ExLM_C00160",
    "ExLM_C00181", "ExLM_C00183", "ExLM_C00185", "ExLM_C00187", "ExLM_C00188",
    "ExLM_C00189", "ExLM_C00191", "ExLM_C00199", "ExLM_C00208", "ExLM_C00209",
    "ExLM_C00216", "ExLM_C00238", "ExLM_C00244", "ExLM_C00249", "ExLM_C00257",
    "ExLM_C00258", "ExLM_C00259", "ExLM_C00286", "ExLM_C00315", "ExLM_C00333",
    "ExLM_C00334", "ExLM_C00369", "ExLM_C00407", "ExLM_C00417", "ExLM_C00458",
    "ExLM_C00459", "ExLM_C00490", "ExLM_C00491", "ExLM_C00492", "ExLM_C00493",
    "ExLM_C00507", "ExLM_C00532", "ExLM_C00590", "ExLM_C00617", "ExLM_C00712",
    "ExLM_C00740", "ExLM_C00794", "ExLM_C00805", "ExLM_C00811", "ExLM_C00956",
    "ExLM_C01083", "ExLM_C01205", "ExLM_C01235", "ExLM_C01494", "ExLM_C01530",
    "ExLM_C01595", "ExLM_C01753", "ExLM_C01793", "ExLM_C02325", "ExLM_C02477",
    "ExLM_C02591", "ExLM_C02646", "ExLM_C03284", "ExLM_C05442", "ExLM_C06427",
    "ExLM_C08362", "ExLM_MC01042", "ExM_C00719", "ExMe_C12303", "ExMe23",
    "ExMe24", "ExMe25", "ExMe27", "ExMe31", "ExMe33", "ExMe35", "ExMe36",
    "ExMe37", "ExMe38", "ExMe42", "ExMe44", "ExMe46", "ExMe48a", "ExMe48b",
    "ExMe49", "ExMe5", "ExMe50", "ExMe54", "ExMe55", "ExMe58", "ExMe8",
    "ExBe1_FWD", "ExBe1_REV", "ExBe10_FWD", "ExBe10_REV", "ExBe12_FWD",
    "ExBe12_REV", "ExBe13_FWD", "ExBe13_REV", "ExBe15_FWD", "ExBe15_REV",
    "ExBe18_FWD", "ExBe18_REV", "ExBe19_FWD", "ExBe19_REV", "ExBe2_FWD",
    "ExBe2_REV", "ExBe21_FWD", "ExBe21_REV", "ExBe3_FWD", "ExBe3_REV",
    "ExBe30_FWD", "ExBe30_REV", "ExBe41_FWD", "ExBe41_REV", "ExBe43_FWD",
    "ExBe43_REV", "ExBe6_FWD", "ExBe6_REV", "ExBe9_FWD", "ExBe9_REV",
    "ExLB_C00698_FWD", "ExLB_C00698_REV", "ExLM_C00698_FWD", "ExLM_C00698_REV",
    "ExMe1_FWD", "ExMe1_REV", "ExMe12_FWD", "ExMe12_REV", "ExMe13_FWD",
    "ExMe13_REV", "ExMe15_FWD", "ExMe15_REV", "ExMe17_FWD", "ExMe17_REV",
    "ExMe18_FWD", "ExMe18_REV", "ExMe19_FWD", "ExMe19_REV", "ExMe2_FWD",
    "ExMe2_REV", "ExMe21_FWD", "ExMe21_REV", "ExMe3_FWD", "ExMe3_REV",
    "ExMe30_FWD", "ExMe30_REV", "ExMe4_FWD", "ExMe4_REV", "ExMe41_FWD",
    "ExMe41_REV", "ExMe43_FWD", "ExMe43_REV", "ExMe9_FWD", "ExMe9_REV"
]
# Get exchange reaction objects
ch_exchange_set = set(ch_exchange_ids)
zma_exchange_set = set(zma_exchange_ids)
exchange_reactions = []
for reaction in model.reactions:
    if reaction.id in ch_exchange_set or reaction.id in zma_exchange_set:
        exchange_reactions.append(reaction)
# Set initial constraints
model.reactions.get_by_id("ExBe15_REV").bounds = (10, model.reactions.get_by_id("ExBe15_REV").bounds[1])
# Create CH protein constraint
ch_variables = {}
for reaction_id, cost in ch_protein_costs.items():
    reaction = model.reactions.get_by_id(reaction_id)
    ch_variables[reaction.forward_variable] = cost
ch_protein_constraint = model.problem.Constraint(
    sum(coeff * var for var, coeff in ch_variables.items()),
    lb=0, ub=600, name="CH_protein_capacity"
)
model.add_cons_vars(ch_protein_constraint)
# Create ZMA protein constraint (chunked approach)
chunk_size = 1000
reactions_list = list(zma_protein_costs.items())
num_chunks = (len(reactions_list) + chunk_size - 1) // chunk_size
chunk_aux_vars = []
for i in range(num_chunks):
    chunk_start = i * chunk_size
    chunk_end = min((i + 1) * chunk_size, len(reactions_list))
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
# STEP 1: Maximize GROWTH (fungal biomass)
print("Step 1: Maximizing GROWTH")
if growth_reaction:
    model.objective = Objective(growth_reaction.forward_variable, direction='max')
    solution1 = model.optimize()
    if solution1.status == "optimal":
        max_growth = growth_reaction.flux
        # Fix GROWTH at maximum value
        growth_constraint = model.problem.Constraint(
            growth_reaction.forward_variable,
            lb=max_growth, ub=max_growth,
            name="fixed_growth"
        )
        model.add_cons_vars(growth_constraint)
    else:
        exit()
else:
    exit()
# STEP 2: Maximize sum of other target reactions
print("Step 2: Maximizing other targets")
if other_target_reaction_objects:
    objective_expr = sum(reaction.forward_variable for reaction in other_target_reaction_objects)
    model.objective = Objective(objective_expr, direction='max')
    solution2 = model.optimize()
    if solution2.status == "optimal":
        max_other_targets_sum = solution2.objective_value
        # Fix sum of other target reactions
        other_targets_constraint = model.problem.Constraint(
            sum(reaction.forward_variable for reaction in other_target_reaction_objects),
            lb=max_other_targets_sum, ub=max_other_targets_sum,
            name="fixed_other_targets"
        )
        model.add_cons_vars(other_targets_constraint)
    else:
        exit()
else:
    max_other_targets_sum = 0
# STEP 3: Minimize exchange reactions
print("Step 3: Minimizing exchange reactions")
if exchange_reactions:
    objective_expr = sum(reaction.forward_variable for reaction in exchange_reactions)
    model.objective = Objective(objective_expr, direction='min')
    solution3 = model.optimize()
    if solution3.status == "optimal":
        min_exchange_sum = solution3.objective_value
    else:
        exit()
else:
    min_exchange_sum = 0
# Calculate and display protein usage
ch_usage = sum(model.reactions.get_by_id(rxn).flux * cost for rxn, cost in ch_protein_costs.items())
zma_usage = sum(model.reactions.get_by_id(rxn).flux * cost for rxn, cost in zma_protein_costs.items())
print(f"CH protein usage: {ch_usage/1000:.3f} g/gDW")
print(f"ZMA protein usage: {zma_usage/1000:.3f} g/gDW")
# Export results to separate Excel files for CH and ZMA
with pd.ExcelWriter("CH_Results.xlsx") as writer:
    # Sheet 1: Flux Distribution
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
    # Sheet 2: Protein Requirements
    ch_protein_data = []
    for reaction_id, cost in ch_protein_costs.items():
        flux = model.reactions.get_by_id(reaction_id).flux
        protein_required = (cost * flux) / 1000  # Convert to g/gDW
        ch_protein_data.append({
            "Reaction_ID": reaction_id,
            "Protein_Required (g/gDW)": protein_required
        })
    # Add total row
    total_ch_protein = ch_usage / 1000
    ch_protein_data.append({
        "Reaction_ID": "TOTAL",
        "Protein_Required (g/gDW)": total_ch_protein
    })
    ch_protein_df = pd.DataFrame(ch_protein_data)
    ch_protein_df.to_excel(writer, sheet_name='Protein_Requirements', index=False)
with pd.ExcelWriter("ZMA_Results.xlsx") as writer:
    # Sheet 1: Flux Distribution
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
    # Sheet 2: Protein Requirements
    zma_protein_data = []
    for reaction_id, cost in zma_protein_costs.items():
        flux = model.reactions.get_by_id(reaction_id).flux
        protein_required = (cost * flux) / 1000  # Convert to g/gDW
        zma_protein_data.append({
            "Reaction_ID": reaction_id,
            "Protein_Required (g/gDW)": protein_required
        })
    # Add total row
    total_zma_protein = zma_usage / 1000
    zma_protein_data.append({
        "Reaction_ID": "TOTAL",
        "Protein_Required (g/gDW)": total_zma_protein
    })
    zma_protein_df = pd.DataFrame(zma_protein_data)
    zma_protein_df.to_excel(writer, sheet_name='Protein_Requirements', index=False)
print("Results exported to CH_Results.xlsx and ZMA_Results.xlsx")