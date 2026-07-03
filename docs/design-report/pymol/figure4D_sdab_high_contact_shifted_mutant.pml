# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure4D_sdab_high_contact_shifted_mutant.pml
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
load structures/sdab_interface/04_unsupported_high_score_interface_shift_04e3eaf858.cif, cand
align cand and chain B, wt and chain B
hide everything, all
show cartoon, wt and chain B
show cartoon, wt and chain A
show cartoon, cand and chain A
cartoon tube, wt and chain A+B
cartoon tube, cand and chain A
color gray80, wt and chain B
color lightblue, wt and chain A
color orange, cand and chain A
set cartoon_transparency, 0.72, wt and chain A
set cartoon_transparency, 0.18, wt and chain B

select wt_ag_iface, byres (wt and chain B within 5 of (wt and chain A))
select cand_ag_iface, byres (cand and chain B within 5 of (cand and chain A))
select cand_ab_iface, byres (cand and chain A within 5 of (cand and chain B))
select cand_mut_sites, cand and chain A and resi 72+100+110+111
show sticks, wt_ag_iface or cand_ag_iface or cand_ab_iface or cand_mut_sites
show spheres, cand_mut_sites and name CA
color marine, wt_ag_iface
color orange, cand_ag_iface
color orange, cand_ab_iface
color yellow, cand_mut_sites

label cand and chain A and resi 110 and name CA, "D110H"
label cand and chain A and resi 111 and name CA, "Y111H"

orient wt_ag_iface or cand_ag_iface or cand_ab_iface
zoom wt_ag_iface or cand_ag_iface or cand_ab_iface or cand_mut_sites, 6
