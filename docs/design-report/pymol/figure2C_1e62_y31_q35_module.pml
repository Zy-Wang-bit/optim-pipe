# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure2C_1e62_y31_q35_module.pml
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
set stick_radius, 0.17
set sphere_scale, 0.48
bg_color white

load structures/AF3-1E62-AeS-1.pdb, e62
hide everything, all
show cartoon, e62 and chain B+C
cartoon tube, e62 and chain B+C
color gray82, e62 and chain C
color lightblue, e62 and chain B
set cartoon_transparency, 0.64, e62 and chain C

select seed_pair, e62 and chain B and resi 31+35
select rescue_set, e62 and chain B and resi 27+32+36
select local_context, e62 and chain B and resi 24+33+34+38
select site_focus, seed_pair or rescue_set or local_context

color orange, seed_pair
color forest, rescue_set
color gray55, local_context
show sticks, byres (seed_pair or rescue_set)
show spheres, (seed_pair or rescue_set or local_context) and name CA
set sphere_scale, 0.32, local_context and name CA

label e62 and chain B and resi 31 and name CA, "Y31H"
label e62 and chain B and resi 35 and name CA, "Q35H"
label e62 and chain B and resi 32 and name CA, "S32"
label e62 and chain B and resi 36 and name CA, "N36"

orient site_focus
zoom site_focus or (e62 and chain C within 7 of site_focus), 6
