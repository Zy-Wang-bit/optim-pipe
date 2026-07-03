# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure1_1e62_vl1_40_sites.pml
reinitialize
set retain_order, 1
set ray_opaque_background, off
set orthoscopic, on
set antialias, 2
set depth_cue, off
set ambient, 0.42
set direct, 0.58
set spec_reflect, 0.12
set label_size, 15
set label_color, black
set stick_radius, 0.16
set sphere_scale, 0.42
bg_color white

load structures/AF3-1E62-AeS-1.pdb, e62
hide everything, all
show cartoon, e62 and chain B+C
cartoon tube, e62 and chain B+C
color gray80, e62 and chain C
color lightblue, e62 and chain B
set cartoon_transparency, 0.62, e62 and chain C

select vl_window, e62 and chain B and resi 1-40
select direct_seed, e62 and chain B and resi 24+26+28+31+33+34
select control_seed, e62 and chain B and resi 35+38
select rescue_main, e62 and chain B and resi 27+32+36
select rescue_context, e62 and chain B and resi 1+3+12+14+17+18+22+25+29+30+40
select protected_site, e62 and chain B and resi 23
select site_focus, direct_seed or control_seed or rescue_main or rescue_context or protected_site

color wheat, vl_window
color orange, direct_seed
color purple, control_seed
color forest, rescue_main
color skyblue, rescue_context
color red, protected_site

show sticks, byres (direct_seed or control_seed or rescue_main or protected_site)
show spheres, (site_focus and name CA)
set sphere_scale, 0.46, (direct_seed or control_seed or rescue_main or protected_site) and name CA
set sphere_scale, 0.28, rescue_context and name CA

label e62 and chain B and resi 24 and name CA, "K24"
label e62 and chain B and resi 31 and name CA, "Y31"
label e62 and chain B and resi 35 and name CA, "Q35"
label e62 and chain B and resi 38 and name CA, "Y38"
label e62 and chain B and resi 23 and name CA, "L23C"

orient site_focus
zoom site_focus or (e62 and chain C within 7 of site_focus), 6
