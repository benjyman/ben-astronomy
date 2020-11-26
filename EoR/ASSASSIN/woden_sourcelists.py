#!/usr/bin/env python
'''Make hpx sky maps and woden skymodels from them for ASSASSIN sims
'''

def write_woden_sourcelists(hpx_fits_filename,freq_MHz,nside,time_string='',dipole_height_m=0.3,pol='X'):

   #try to get altaz working
   #source_name = '3C353'
   #source = SkyCoord.from_name(source_name)
   #altaz = source.transform_to(AltAz(obstime=time,location=eda2_loc))
   #print("source's Altitude = {0.alt:.2}".format(altaz))
   #sys.exit()
   
   wavelength = 300./freq_MHz
   name_base = hpx_fits_filename.split(".fits")[0] 
   
   source_list_name = '%s_sourcelist.txt' % name_base

   data = hp.read_map(hpx_fits_filename,nest=False)
   pix_inds = np.arange(hp.nside2npix(nside))
   l, b = hp.pix2ang(nside,pix_inds,lonlat=True)
   gal_coords = SkyCoord(l*u.deg, b*u.deg, frame='galactic')
   ra = gal_coords.icrs.ra.value
   dec = gal_coords.icrs.dec.value
   #need to change this for brightness temp not MJy/steradian
   #fluxes = data*hp.nside2pixarea(nside,degrees=False)*1e+6
   pix_area_sr = hp.nside2pixarea(nside,degrees=False)
   scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
   print("wavelength is %0.3f m, scale map by %s to get to Jy/pix" % (wavelength,scale))
   data_jy_per_pix = data * scale
   #fig = plt.figure(figsize=(10,10))
   #hp.mollview(log10(data), sub=(2,1,1), fig=fig,title='Galactic')
   #hp.mollview(log10(fluxes), sub=(2,1,2), fig=fig,title='Equatorial')
   #fig.savefig('%s_woden_map.png' % name_base, bbox_inches='tight')
   #plt.close()
   
   if time_string!='':
      year,month,day,hour,min,sec = time_string.split('_')
      time = Time('%d-%02d-%02d %02d:%02d:%02d' % (float(year),float(month),float(day),float(hour),float(min),float(sec)))
   
      angular_source_list_name = '%s_%s_pol_%s_angular_sourcelist.txt' % (name_base,time_string,pol)
      beam_source_list_name = '%s_%s_pol_%s_beam_sourcelist.txt' % (name_base,time_string,pol)
      global_foreground_source_list_name = '%s_%s_pol_%s_global_foreground_sourcelist.txt' % (name_base,time_string,pol)
   
      #beam stuff:
      altaz = gal_coords.transform_to(AltAz(obstime=time,location=eda2_loc))
      theta_array = 90. - altaz.alt.value
      phi_array = 90. - altaz.az.value
      beam_data_array = calc_beam_values(theta_array,phi_array,pol,dipole_height_m,wavelength)

      #work out global foreground signal and angular structure sourcelist 
      sky_with_beam = beam_data_array * data
      sum_of_beam_weights = np.nansum(beam_data_array)
      global_foreground_signal_value_K = np.nansum(sky_with_beam) /  sum_of_beam_weights  
      print("beam_weighted_av_sky at %0.3f MHz is %0.4E K" % (freq_MHz,global_foreground_signal_value_K))
      
      global_foreground_data_K = (data * 0.) + global_foreground_signal_value_K
      global_foreground_data_jy_per_pix = global_foreground_data_K * scale
      
      data_angular_only_K = data - global_foreground_signal_value_K
      data_angular_only_jy_per_pix = data_angular_only_K * scale
   
            
      #make some maps as a test:
      if np.isclose(freq_MHz,50.0):
          plt.clf()
          map_title="GSM from MWA at %s %0.3f MHz" % (time_string,freq_MHz)
          ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
          ##hp.orthview(map=sky_with_beam,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
          hp.mollview(map=data,coord='G',title=map_title,min=0, max=7000)
          fig_name="check_gsm_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (time_string,freq_MHz,pol)
          figmap = plt.gcf()
          figmap.savefig(fig_name,dpi=500)
          print("saved %s" % fig_name)

          plt.clf()
          map_title="beam from MWA at %s %0.3f MHz" % (time_string,freq_MHz)
          ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
          ##hp.orthview(map=sky_with_beam,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
          hp.mollview(map=beam_data_array,coord='G',title=map_title,min=0, max=1)
          fig_name="check_beam_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (time_string,freq_MHz,pol)
          figmap = plt.gcf()
          figmap.savefig(fig_name,dpi=500)
          print("saved %s" % fig_name)  
        
          
          plt.clf()
          map_title="sky with beam from MWA at %s %0.3f MHz" % (time_string,freq_MHz)
          ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
          ##hp.orthview(map=sky_with_beam,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
          hp.mollview(map=sky_with_beam,coord='G',title=map_title,min=0, max=7000)
          fig_name="check_sky_with_beam_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (time_string,freq_MHz,pol)
          figmap = plt.gcf()
          figmap.savefig(fig_name,dpi=500)
          print("saved %s" % fig_name)  
          
          
      source_ind = 0
      with open(global_foreground_source_list_name,'w') as outfile:
          for ind,val in enumerate(global_foreground_data_jy_per_pix):
              if source_ind == 0:
                  outfile.write('SOURCE pygsm P %d G 0 S 0 0\n' %len(global_foreground_data_jy_per_pix))
              outfile.write('COMPONENT POINT %.7f %.6f\n' %(ra[ind]/15.0,dec[ind]))
              outfile.write('LINEAR 150e+6 %.10f 0 0 0 0.0\n' % val)
              outfile.write('ENDCOMPONENT\n')
              source_ind += 1
          outfile.write('ENDSOURCE')
      print("wrote %s" % global_foreground_source_list_name)
       
      source_ind = 0
      with open(angular_source_list_name,'w') as outfile:
          for ind,val in enumerate(data_angular_only_jy_per_pix):
              if source_ind == 0:
                  outfile.write('SOURCE pygsm P %d G 0 S 0 0\n' %len(data_angular_only_jy_per_pix))
              outfile.write('COMPONENT POINT %.7f %.6f\n' %(ra[ind]/15.0,dec[ind]))
              outfile.write('LINEAR 150e+6 %.10f 0 0 0 0.0\n' % val)
              outfile.write('ENDCOMPONENT\n')
              source_ind += 1
          outfile.write('ENDSOURCE')
      print("wrote %s" % angular_source_list_name)
          
      source_ind = 0
      with open(beam_source_list_name,'w') as outfile:
          for ind,beam_val in enumerate(beam_data_array):
              if source_ind == 0:
                  outfile.write('SOURCE pygsm P %d G 0 S 0 0\n' %len(beam_data_array))
              outfile.write('COMPONENT POINT %.7f %.6f\n' %(ra[ind]/15.0,dec[ind]))
              outfile.write('LINEAR 150e+6 %.10f 0 0 0 0.0\n' % beam_val)
              outfile.write('ENDCOMPONENT\n')
              source_ind += 1
          outfile.write('ENDSOURCE')
      print("wrote %s" % beam_source_list_name)

   
   else:
      global_foreground_signal_value_K = np.nan
                
   source_ind = 0
   with open(source_list_name,'w') as outfile:
       for ind,data_val in enumerate(data_jy_per_pix):
           if source_ind == 0:
               outfile.write('SOURCE pygsm P %d G 0 S 0 0\n' %len(data_jy_per_pix))
           outfile.write('COMPONENT POINT %.7f %.6f\n' %(ra[ind]/15.0,dec[ind]))
           outfile.write('LINEAR 150e+6 %.10f 0 0 0 0.0\n' % data_val)
           outfile.write('ENDCOMPONENT\n')
           source_ind += 1
       outfile.write('ENDSOURCE')
   print("wrote %s" % source_list_name)   
       
   return(global_foreground_signal_value_K)

def get_beam_value(theta,phi,dipole_height_m,wavelength,pol):
   #theta = ZA, phi angle anticlockwise from East looking down on array
      #set to zero below horizon
   if theta > np.pi/2.:
      short_dipole_parallel_beam_value = 0
   else:   
      ##This is for YY dipole: !
      if (pol=='Y'):
         theta_parallel=np.arccos(np.sin(theta)*np.cos(phi))
      else:
      #This is for XX dipole!
         theta_parallel=np.arccos(np.sin(theta)*np.sin(phi)) 
      
      d_in_lambda = (2. * dipole_height_m)/wavelength
      gp_effect = 2.*np.sin(np.pi*d_in_lambda*np.cos(theta))
      voltage_parallel=np.sin(theta_parallel) * gp_effect
      short_dipole_parallel_beam_value = voltage_parallel**2       #* factor
   
   return(short_dipole_parallel_beam_value)
   
   #normalise to one at max
   beam_max = np.max(short_dipole_parallel_beam_map)
   short_dipole_parallel_beam_map = short_dipole_parallel_beam_map / beam_max
   
   return short_dipole_parallel_beam_map

def write_woden_skymodels(freq_MHz_list,nside,time_string,dipole_height_m,pol_list,fine_chan_khz=10):
   start_freq_MHz = float(freq_MHz_list[0])
   gsm = GlobalSkyModel()
   for pol in pol_list:
        global_foreground_value_list = []
        for band_num,freq_MHz in enumerate(freq_MHz_list):
           freq_MHz = float(freq_MHz)
           print("Freq %0.3f MHz" % freq_MHz)
           #name_base = "woden_map_centre_chan_%03d_band_%02d_freq_%0.3f_MHz_hpx" % (centre_chan,band_num,freq_MHz)
           #dont put freq as stuffs up naming on pawsey for array job
           name_base = "woden_map_start_freq_%0.3f_band_%s_hpx" % (start_freq_MHz,band_num)
           gsm_filename = "%s_gsm.fits" % name_base
           gsm_map = gsm.generate(freq_MHz)
           hp.write_map(gsm_filename,gsm_map,coord='G',nest=False,overwrite=True)
           print("saved %s" % gsm_filename)
           global_foreground_value = write_woden_sourcelists(gsm_filename,freq_MHz,nside,time_string,dipole_height_m,pol) 
           global_foreground_value_list.append(global_foreground_value)

        name_base = "woden_map_start_freq_%0.3f_hpx" % (start_freq_MHz)
        global_foreground_value_array = np.asarray(global_foreground_value_list)        
        global_foreground_value_array_filename = "%s_%s_pol_%s_global_foreground.npy" % (name_base,time_string,pol)
        np.save(global_foreground_value_array_filename,global_foreground_value_array)
        
        print("saved %s" % global_foreground_value_array_filename)          
   
   for band_num,freq_MHz in enumerate(freq_MHz_list):
        freq_MHz = float(freq_MHz)
        print("Freq %0.3f MHz" % freq_MHz)
        #name_base = "woden_map_centre_chan_%03d_band_%02d_freq_%0.3f_MHz_hpx" % (centre_chan,band_num,freq_MHz)
        #dont put freq as stuffs up naming on pawsey for array job
        name_base = "woden_map_start_freq_%0.3f_band_%s_hpx" % (start_freq_MHz,band_num)        
        EDGES_uniform_filename = "%s_EDGES_uniform.fits" % name_base
        unity_uniform_filename = "%s_unity_uniform.fits" % name_base
        freq_MHz_array = np.asarray([freq_MHz])
        s_21_array_EDGES = plot_S21_EDGES(nu_array=freq_MHz_array)
        s_21_array_EDGES_value = s_21_array_EDGES[0]
        global_EDGES_uniform_map = (gsm_map * 0.0) + s_21_array_EDGES_value
        hp.write_map(EDGES_uniform_filename,global_EDGES_uniform_map,coord='G',nest=False,overwrite=True)
        print("saved %s" % EDGES_uniform_filename)
        write_woden_sourcelists(EDGES_uniform_filename,freq_MHz,nside) 
        #print(global_EDGES_uniform_map)
        unity_uniform_sky_temp = 1.
        unity_map_uniform = (gsm_map * 0.0) + unity_uniform_sky_temp
        hp.write_map(unity_uniform_filename,unity_map_uniform,coord='G',nest=False,overwrite=True)
        print("saved %s" % unity_uniform_filename)
        write_woden_sourcelists(unity_uniform_filename,freq_MHz,nside) 

if __name__ == "__main__":
    import argparse
    
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)
            
            from argparse import RawTextHelpFormatter
            
    parser = argparse.ArgumentParser(description="Run this to make WODEN sky models on \
            supercomputer with array job",formatter_class=SmartFormatter)

    parser.add_argument('--band', default='0',
        help='Just does one band at a time. Frequency depends on other arguments such as --daniel e.g. --band_nums=1')
            
    parser.add_argument('--daniel', default=False, action='store_true',
        help='Frequencies correspond to Daniels sims. Band 0 is at 50 MHz and frequency goes up in increments of 1.28 MHz')
            
            
    args = parser.parse_args()
    
    if args.band: band = int(args.band)
    
    
    if (args.daniel):
       print('band %s daniel' % band)       
    else:
       print('band %s' % band) 
       
       
       
            