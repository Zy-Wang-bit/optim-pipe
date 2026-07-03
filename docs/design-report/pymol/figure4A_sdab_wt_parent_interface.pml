# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure4A_sdab_wt_parent_interface.pml
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
set sphere_scale, 0.38
bg_color white

load structures/sdab_interface/00_WT_parent.cif, wt
hide everything, all
show cartoon, wt and chain A+B
cartoon tube, wt and chain A+B
color lightblue, wt and chain A
color gray80, wt and chain B
set cartoon_transparency, 0.18, wt and chain B

select wt_ab_iface, byres (wt and chain A within 5 of (wt and chain B))
select wt_ag_iface, byres (wt and chain B within 5 of (wt and chain A))
show sticks, wt_ab_iface or wt_ag_iface
color marine, wt_ab_iface
color slate, wt_ag_iface

orient wt_ab_iface or wt_ag_iface
zoom wt_ab_iface or wt_ag_iface, 6
