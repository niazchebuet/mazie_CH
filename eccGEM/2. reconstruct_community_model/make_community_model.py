from cobra import io, Model

# File paths
fungus_file = "iCH384_id_change.xml"
leaf_file = "izma6518_L_exchange_added.xml"
output_file = "cochliobolus_maize_leaf_community.xml"

# Fungus-specific boundary exchange reactions to remove
fungus_boundary_rxns_to_remove = [
    "BOUNDARY_GLCe", "BOUNDARY_PYRe", "BOUNDARY_ASNe", "BOUNDARY_ASPe",
    "BOUNDARY_FRUe", "BOUNDARY_GLACe", "BOUNDARY_MANe",
    "BOUNDARY_H2Oe", "BOUNDARY_HNO3e", "BOUNDARY_PIe",
    "BOUNDARY_SLFe", "BOUNDARY_NH3e"
]

# Leaf-specific exchange reactions to remove (update these names if needed to match your model!)
leaf_exchange_rxns_to_remove = [
    "R_ExBe3", "R_ExBe4", "R_ExBe5", "R_ExBe6", "R_ExBe7", "R_ExBe8", "R_ExBe9"
]

def safe_remove_rxns(model, rxn_ids, org):
    to_remove = []
    missing = []
    for rid in rxn_ids:
        if rid in model.reactions:
            to_remove.append(rid)
        else:
            missing.append(rid)
    model.remove_reactions(to_remove)
    if missing:
        print(f"[{org}] Skipped missing reactions: {missing}")

# --- Load models ---
print("Loading fungus model ...")
fungus = io.read_sbml_model(fungus_file)
print("Loading leaf model ...")
leaf = io.read_sbml_model(leaf_file)

# --- Remove specified exchanges ---
print("\nRemoving fungus-only exchanges ...")
safe_remove_rxns(fungus, fungus_boundary_rxns_to_remove, "Fungus")
print("Removing leaf-only exchanges ...")
safe_remove_rxns(leaf, leaf_exchange_rxns_to_remove, "Leaf")

# --- Merge ---
print("\nMerging models ...")
community = Model("cochliobolus_maize_leaf_community")
community.add_reactions(fungus.reactions)
community.add_reactions(leaf.reactions)

# --- Save merged model ---
print(f"\nSaving model as {output_file} ...")
io.write_sbml_model(community, output_file)
print("Saved community model:", output_file)