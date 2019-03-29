#make EDA analytical beams using Randall's IDL code

.comp /data/code/git/ben-astronomy/AAVS-1/Randall_IDL/antennas.pro

dipole_height_m = 0.3
d_m = 2*dipole_height_m


make_azza_arrays_fov,180,az,za,mask,180.0
field_x = az*0.0
p = where(za lt !pi/2)
ground_plane_effect = za*0.0
proj_x = sqrt(1-(sin(az[p])*sin(za[p]))^2)

freq_MHz = 160.0

lambda = 300.0/freq_MHz
d_in_lambda = d_m/lambda
gp_norm = end_fire_2element(d_in_lambda,0.0)
ground_plane_effect[p] = end_fire_2element(d_in_lambda,za[p])/gp_norm
field_x[p] = proj_x*ground_plane_effect[p]
; include a 1/cos(za) term to account for pixel area changing in sine projection
write_beam_fits,float(field_x^2)/cos(za),'model_xx.fits',0.0,0.0,160.0,0.0