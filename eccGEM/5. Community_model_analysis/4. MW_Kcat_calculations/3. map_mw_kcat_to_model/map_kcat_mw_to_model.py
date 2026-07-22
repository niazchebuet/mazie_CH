import cobra
import pandas as pd

# Load the COBRA model
model = cobra.io.read_sbml_model('irreversible_cochliobolus_maize_leaf_community.xml')

# Print model statistics
print("Model Statistics:")
print(f"Model ID: {model.id}")
print(f"Total reactions: {len(model.reactions)}")
print(f"Total metabolites: {len(model.metabolites)}")
print(f"Total genes: {len(model.genes)}")
print(f"Reversible reactions: {len([rxn for rxn in model.reactions if rxn.reversibility])}")
print(f"Irreversible reactions: {len([rxn for rxn in model.reactions if not rxn.reversibility])}")
print(f"Exchange reactions: {len(model.exchanges)}")
print(f"Demand reactions: {len(model.demands)}")
print(f"Sink reactions: {len(model.sinks)}")
print("-" * 50)

# Read the MW_kcat data from TurNuP results
kcat_data = pd.read_excel('kcat_mw_results_TurNuP.xlsx')

# Get all reactions from the model
all_reactions = [rxn.id for rxn in model.reactions]

# Process MW_kcat data - handle duplicates by taking minimum value for each reaction
kcat_dict = {}
for _, row in kcat_data.iterrows():
    reaction = row['Reactions']
    mw_kcat = row['MW_kcat']
    
    if reaction in kcat_dict:
        # If reaction already exists, take the minimum MW_kcat value
        kcat_dict[reaction] = min(kcat_dict[reaction], mw_kcat)
    else:
        kcat_dict[reaction] = mw_kcat

# Find global minimum MW_kcat value from the entire dataset
global_min_mw_kcat = kcat_data['MW_kcat'].min()

# Create result list
result_data = []
missing_reactions = 0

for reaction in all_reactions:
    if reaction in kcat_dict:
        # Use the minimum MW_kcat value for this reaction
        mw_kcat_value = kcat_dict[reaction]
    else:
        # Use global minimum for missing reactions
        mw_kcat_value = global_min_mw_kcat
        missing_reactions += 1
    
    result_data.append({
        'Reactions': reaction,
        'MW_kcat': mw_kcat_value
    })

# Create DataFrame
result_df = pd.DataFrame(result_data)

# Save to Excel file
result_df.to_excel('community_model_per_unit_protein_data.xlsx', index=False)

# Print summary information
print("Analysis complete!")
print(f"Total reactions in model: {len(all_reactions)}")
print(f"Reactions with MW_kcat data: {len(kcat_dict)}")
print(f"Reactions using global minimum: {missing_reactions}")
print(f"Coverage: {(len(kcat_dict)/len(all_reactions)*100):.1f}%")
print(f"Global minimum MW_kcat value: {global_min_mw_kcat}")
print(f"Output saved to: community_model_per_unit_protein_data.xlsx")