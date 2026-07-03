# Run PyMOL from this directory:
#   cd docs/design-report/pymol
#   pymol figure3_sdab_vhh72_111_sites.pml
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

load structures/sdAb_mutation_frequency_bfactor.pdb, sdab
hide everything, all
show cartoon, sdab and chain A+B
cartoon tube, sdab and chain A+B
color gray80, sdab and chain B
color coral, sdab and chain A
set cartoon_transparency, 0.62, sdab and chain B

select vhh_window, sdab and chain A and resi 72-111
select primary_seed, sdab and chain A and resi 100+110+111
select risk_seed, sdab and chain A and resi 102+105+108
select rescue_main, sdab and chain A and resi 76+87+89+101+107
select rescue_context, sdab and chain A and resi 72+81+82+85+88+91+93
select protected_site, sdab and chain A and resi 96
select site_focus, primary_seed or risk_seed or rescue_main or rescue_context or protected_site

color wheat, vhh_window
color orange, primary_seed
color purple, risk_seed
color forest, rescue_main
color skyblue, rescue_context
color red, protected_site

show sticks, byres (primary_seed or risk_seed or rescue_main or protected_site)
show spheres, site_focus and name CA
set sphere_scale, 0.46, (primary_seed or risk_seed or rescue_main or protected_site) and name CA
set sphere_scale, 0.28, rescue_context and name CA

label sdab and chain A and resi 100 and name CA, "Q100"
label sdab and chain A and resi 110 and name CA, "D110"
label sdab and chain A and resi 111 and name CA, "Y111"
label sdab and chain A and resi 96 and name CA, "A96C"

orient site_focus
zoom site_focus or (sdab and chain B within 7 of site_focus), 6
