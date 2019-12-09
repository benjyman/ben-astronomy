pro make_EDA_analytical_beams
   ; .comp /data/code/git/ben-astronomy/AAVS-1/Randall_IDL/antennas.pro
   
   dipole_height_m = 0.3
   d_m = 2*dipole_height_m
   
   
   make_azza_arrays_fov,180,az,za,mask,180.0
   field_x = az*0.0
   field_y = az*0.0
   p = where(za lt !pi/2)
   ground_plane_effect = za*0.0
   proj_x = sqrt(1-(sin(az[p])*sin(za[p]))^2)
   proj_y = sqrt(1-(cos(az[p])*sin(za[p]))^2)
   
   for freq_MHz = 50, 150 do begin
   
      beam_name_string_x = 'model_' + STRTRIM(freq_MHz, 2) + '_MHz_xx.fits'
      beam_name_string_no_cos_za_x = 'model_' + STRTRIM(freq_MHz, 2) + '_MHz_xx_no_cos_za.fits'
      beam_name_string_y = 'model_' + STRTRIM(freq_MHz, 2) + '_MHz_yy.fits'
      beam_name_string_no_cos_za_y = 'model_' + STRTRIM(freq_MHz, 2) + '_MHz_yy_no_cos_za.fits'
   
      lambda = 300.0/freq_MHz
      d_in_lambda = d_m/lambda
      gp_norm = end_fire_2element(d_in_lambda,0.0)
      ground_plane_effect[p] = end_fire_2element(d_in_lambda,za[p])/gp_norm
      field_x[p] = proj_x*ground_plane_effect[p]
      field_y[p] = proj_y*ground_plane_effect[p]
      ;xx
      ; include a 1/cos(za) term to account for pixel area changing in sine projection
      write_beam_fits,float(field_x^2)/cos(za),beam_name_string_x,0.0,0.0,freq_MHz,0.0
      ; write one without the 1/cos(za) term to use for the beam weights
      write_beam_fits,float(field_x^2),beam_name_string_no_cos_za_x,0.0,0.0,freq_MHz,0.0
      ;yy
      ; include a 1/cos(za) term to account for pixel area changing in sine projection
      write_beam_fits,float(field_y^2)/cos(za),beam_name_string_y,0.0,0.0,freq_MHz,0.0
      ; write one without the 1/cos(za) term to use for the beam weights
      write_beam_fits,float(field_y^2),beam_name_string_no_cos_za_y,0.0,0.0,freq_MHz,0.0
   endfor
end