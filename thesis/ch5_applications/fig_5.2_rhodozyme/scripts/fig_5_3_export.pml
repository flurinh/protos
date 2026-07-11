# Figure 5.3 — Export
# Run after setting views: @fig_5_3_export.pml
scene panel_A, recall
ray 2400, 1800
png ../figures/fig_5_3a_candidate_region.png, dpi=300

scene panel_B, recall
ray 2400, 1800
png ../figures/fig_5_3b_placed_theozyme.png, dpi=300

save ../figures/fig_5_3_theozyme_placement.pse
print("Export complete.")
