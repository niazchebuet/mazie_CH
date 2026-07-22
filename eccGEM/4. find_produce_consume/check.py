import cobra
from cobra.flux_analysis import pfba

# Load model
model = cobra.io.read_sbml_model('irreversible_cochliobolus_maize_leaf_community.xml')

# Set objective to GROWTH (maximize biomass)
model.objective = 'Bio_Nplus'

# Knock out reactions
reactions_to_knockout = [' ']
for reaction_id in reactions_to_knockout:
    if reaction_id in model.reactions:
        model.reactions.get_by_id(reaction_id).bounds = (0, 0)

# Original optimization
pfba_solution = model.optimize()
print(f"Original Growth: {pfba_solution.objective_value:.6f} h⁻¹")

# Find active AKGm reactions and print table
print(f"\n{'Reactions':<15} {'Equation':<60} {'Flux':<15}")
print("-" * 90)

for reaction in model.reactions:
    if any('C00149_d' in met.id for met in reaction.metabolites):
        flux = pfba_solution.fluxes[reaction.id]
        if abs(flux) > 1e-6:
            equation = reaction.reaction
            print(f"{reaction.id:<15} {equation:<60} {flux:<15.6f}")