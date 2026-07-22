from cobra import io
import pandas as pd

leaf_file = "izma6518_L_exchange_added.xml"
excel_output = "producing_metabolic_reactions.xlsx"

mets_of_interest = [
    "C04831_p", "C04831_q", "C00509_d", "C00509_c", "C00509_p", "C00509_q",
    "C00129_c", "C00129_d", "C00129_p", "C00129_q", "C00129_x",
    "C00029_c", "C00029_d", "C00029_q", "C00029_p", "C00029_m", "C00029_n", "C00029_r", "C00029_s", "C00029_x"
]

def is_metabolic(reaction):
    rxn_id = reaction.id.lower()
    non_metabolic_patterns = [
        'ex_', 'exchange', 'exbe', 'boundary', 'transporter', 'transport', 'tp_', 'sink', 'drain'
    ]
    return not any([pat in rxn_id for pat in non_metabolic_patterns])

print("Loading leaf model ...")
leaf_model = io.read_sbml_model(leaf_file)

# First sheet: metabolite-producing metabolic reactions
results = []
all_metabolic_rxn_ids = set()
for met_id in mets_of_interest:
    if met_id in leaf_model.metabolites:
        met = leaf_model.metabolites.get_by_id(met_id)
        producing_rxns = [
            rxn.id for rxn in met.reactions
            if rxn.get_coefficient(met) > 0 and is_metabolic(rxn)
        ]
        results.append({
            "Metabolite ID": met_id,
            "Number of Producing Metabolic Reactions": len(producing_rxns),
            "Metabolic Producing Reaction IDs": ", ".join(producing_rxns) if producing_rxns else "None"
        })
        all_metabolic_rxn_ids.update(producing_rxns)
    else:
        results.append({
            "Metabolite ID": met_id,
            "Number of Producing Metabolic Reactions": 0,
            "Metabolic Producing Reaction IDs": "Not found in model"
        })

df_prod = pd.DataFrame(results)

# Second sheet: all metabolic reactions in the model
all_metabolic_rxns = [
    rxn for rxn in leaf_model.reactions
    if is_metabolic(rxn)
]

df_all_metabolic = pd.DataFrame({
    "Metabolic Reaction ID": [rxn.id for rxn in all_metabolic_rxns],
    "Metabolic Reaction Name": [rxn.name for rxn in all_metabolic_rxns]
})

with pd.ExcelWriter(excel_output) as writer:
    df_prod.to_excel(writer, sheet_name="ProducingRxns", index=False)
    df_all_metabolic.to_excel(writer, sheet_name="AllMetabolicReactions", index=False)

print("Excel file created with two sheets:", excel_output)