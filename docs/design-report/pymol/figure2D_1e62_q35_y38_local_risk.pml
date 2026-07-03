# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure2D_1e62_q35_y38_local_risk.pml
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

select seed_pair, e62 and chain B and resi 35+38
select local_risk, e62 and chain B and resi 32+33
select local_context, e62 and chain B and resi 24+27+31+34+36
select site_focus, seed_pair or local_risk or local_context

color purple, seed_pair
color red, local_risk
color gray55, local_context
show sticks, byres (seed_pair or local_risk)
show spheres, (seed_pair or local_risk or local_context) and name CA
set sphere_scale, 0.32, local_context and name CA

label e62 and chain B and resi 35 and name CA, "Q35H"
label e62 and chain B and resi 38 and name CA, "Y38H"
label e62 and chain B and resi 32 and name CA, "S32 risk"
label e62 and chain B and resi 33 and name CA, "S33 risk"

orient site_focus
zoom site_focus or (e62 and chain C within 7 of site_focus), 6
