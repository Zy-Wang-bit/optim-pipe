# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure4C_sdab_partial_shifted_mutant.pml
reinitialize
set retain_order, 1
set ray_opaque_background, off
set orthoscopic, on
set antialias, 2
set depth_cue, off
set ambient, 0.42
set direct, 0.58
set spec_reflect, 0.12
set label_size, 14
set label_color, black
set stick_radius, 0.16
set sphere_scale, 0.42
bg_color white

load structures/sdab_interface/00_WT_parent.cif, wt
load structures/sdab_interface/02_plausible_partial_shifted_003413d2c0.cif, cand
align cand and chain B, wt and chain B
hide everything, all
show cartoon, wt and chain B
show cartoon, wt and chain A
show cartoon, cand and chain A
cartoon tube, wt and chain A+B
cartoon tube, cand and chain A
color gray80, wt and chain B
color lightblue, wt and chain A
color cyan, cand and chain A
set cartoon_transparency, 0.72, wt and chain A
set cartoon_transparency, 0.18, wt and chain B

select wt_ag_iface, byres (wt and chain B within 5 of (wt and chain A))
select cand_ag_iface, byres (cand and chain B within 5 of (cand and chain A))
select cand_ab_iface, byres (cand and chain A within 5 of (cand and chain B))
select cand_mut_sites, cand and chain A and resi 100+108+110
show sticks, wt_ag_iface or cand_ag_iface or cand_ab_iface or cand_mut_sites
show spheres, cand_mut_sites and name CA
color marine, wt_ag_iface
color cyan, cand_ag_iface
color cyan, cand_ab_iface
color yellow, cand_mut_sites

label cand and chain A and resi 100 and name CA, "Q100H"
label cand and chain A and resi 110 and name CA, "D110H"

orient wt_ag_iface or cand_ag_iface or cand_ab_iface
zoom wt_ag_iface or cand_ag_iface or cand_ab_iface or cand_mut_sites, 6
