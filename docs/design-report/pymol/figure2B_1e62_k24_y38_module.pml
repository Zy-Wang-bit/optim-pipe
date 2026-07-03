# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure2B_1e62_k24_y38_module.pml
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

select seed_pair, e62 and chain B and resi 24+38
select context_rescue, e62 and chain B and resi 1+3+12+14+17+18+22+25
select local_context, e62 and chain B and resi 27+31+32+35+36
select site_focus, seed_pair or context_rescue or local_context

color orange, seed_pair
color skyblue, context_rescue
color gray55, local_context
show sticks, byres (seed_pair)
show spheres, (seed_pair or context_rescue or local_context) and name CA
set sphere_scale, 0.30, context_rescue and name CA
set sphere_scale, 0.30, local_context and name CA

label e62 and chain B and resi 24 and name CA, "K24H"
label e62 and chain B and resi 38 and name CA, "Y38H"

orient site_focus
zoom site_focus or (e62 and chain C within 7 of site_focus), 6
