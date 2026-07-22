from cobra import io, Model, Reaction, Metabolite

# Load the metabolic model from an SBML file
model = io.read_sbml_model('iCH384_id_change.xml')

# Run FBA
solution = model.optimize()

# Print the objective value
print("Objective value: ", solution.objective_value)