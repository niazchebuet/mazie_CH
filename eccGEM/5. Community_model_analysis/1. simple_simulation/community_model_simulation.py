from cobra import io
import pandas as pd
from cobra.flux_analysis import pfba

# Path to the metabolic community SBML file
community_model_file = "cochliobolus_maize_leaf_community.xml"
# Excel filename to store results
excel_output = "community_stage2_pfba_active_flux_distribution.xlsx"

# Reaction IDs for key steps (edit these if your model changes)
target_rxns = [
    "MR00338_LSQBKT_d_RSQBKT_",
    "R02446_LSQBKT_p_RSQBKT_",
    "R02446_LSQBKT_q_RSQBKT_",
    "R01121_LSQBKT_c_RSQBKT_",
    "R01121_LSQBKT_d_RSQBKT_",
    "R05611_LSQBKT_p_RSQBKT_",#FWD
    "R05884_LSQBKT_q_RSQBKT_",
    "R00289_LSQBKT_c_RSQBKT_",#REV
    "R00289_LSQBKT_d_RSQBKT_",#REV
    "R05884_LSQBKT_p_RSQBKT_",
    "R05611_LSQBKT_q_RSQBKT_",#REV
    "R09248_LSQBKT_q_RSQBKT_"
    "R00086_LSQBKT_m_RSQBKT_"
    "R00086_LSQBKT_t_RSQBKT_"
    "R00086_LSQBKT_n_RSQBKT_"
    "R00086_LSQBKT_u_RSQBKT_"
    "R00086_LSQBKT_i_RSQBKT_"#REV
    "R00086_LSQBKT_j_RSQBKT_"#REV
]
fungal_biomass_rxn = "GROWTH"     # Fungal biomass reaction
photon_uptake_rxn = "ExBe15"      # Photon uptake reaction
leaf_biomass_rxn = "Bio_Nplus"    # Leaf biomass reaction

# Load the SBML metabolic community model
community = io.read_sbml_model(community_model_file)

# --------------------------------------------------------
# 1. Stage 1: Maximize fungal biomass (GROWTH)
# Objective: max v_GROWTH
# Constraints: all model stoichiometry and bounds
# --------------------------------------------------------
community.objective = fungal_biomass_rxn
sol1 = community.optimize()
growth_value = sol1.fluxes[fungal_biomass_rxn]

# --------------------------------------------------------
# 2. Stage 2a: Fix fungal biomass and photon uptake;
#    maximize sum of leaf pathogenic defense metabolites reaction fluxes
# Constraints:
#   - v_GROWTH = optimal value from Stage 1
#   - v_ExBe15 = -100 (fixed photon import)
# Objective: max sum_i v_product_i
# --------------------------------------------------------
community.reactions.get_by_id(fungal_biomass_rxn).bounds = (growth_value, growth_value)
if photon_uptake_rxn in community.reactions:
    community.reactions.get_by_id(photon_uptake_rxn).bounds = (-100, -100)
target_found = [community.reactions.get_by_id(rid).flux_expression for rid in target_rxns if rid in community.reactions]
community.objective = sum(target_found)
sol2 = community.optimize()

# --------------------------------------------------------
# 3. Stage 2b: pFBA for realistic fluxes
# Constraints:
#   - v_GROWTH = optimal
#   - v_ExBe15 = -100
#   - v_product_i = optimal from Stage 2a for leaf pathogenic defense metabolites reaction fluxes
# Objective: minimize total flux sum over internal reactions (pFBA)
# --------------------------------------------------------
for rid in target_rxns:
    if rid in community.reactions:
        flux_val = sol2.fluxes[rid]
        community.reactions.get_by_id(rid).bounds = (flux_val, flux_val)
pfba_solution = pfba(community)

# --------------------------------------------------------
# Export only reactions carrying significant flux (|Flux| >= 0.00001)
# Excel contains Reaction ID, Flux, Reaction Name
# --------------------------------------------------------
filtered_flux_info = []
for rid, val in pfba_solution.fluxes.items():
    # Only include reactions with |flux| >= 1e-5
    if abs(val) >= 1e-5:
        rxn_name = community.reactions.get_by_id(rid).name if rid in community.reactions else ""
        filtered_flux_info.append((rid, val, rxn_name))
filtered_flux_df = pd.DataFrame(filtered_flux_info, columns=["Reaction ID", "Flux", "Reaction Name"])
filtered_flux_df.to_excel(excel_output, index=False)