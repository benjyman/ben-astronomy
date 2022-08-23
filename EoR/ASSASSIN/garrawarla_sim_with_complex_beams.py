#!/usr/bin/env python3

import os,sys
import healpy as hp
import numpy as np
import math
from pygdsm import GlobalSkyModel
import h5py
from scipy.interpolate import interp2d,griddata,interp1d


mwa_latitude_deg = -26.70331940
mwa_longitude_deg = 116.670575
mwa_latitude_rad = float(mwa_latitude_deg/180.*np.pi)
mwa_longitude_rad = float(mwa_longitude_deg/180.*np.pi)

fine_chans_per_EDA2_chan = 27
fine_chan_width_Hz =  (400./512.) / 27. * 1.0e6 #(see Marcin email)
fine_chan_width_MHz = fine_chan_width_Hz / 1.0e6

gsm_dir = "/astro/mwaeor/bmckinley/EoR/EDA2/gsm_maps/"

def calc_uvw(x_diff,y_diff,z_diff,freq_MHz,hourangle=0.,declination=mwa_latitude_rad):
   #https://web.njit.edu/~gary/728/Lecture6.html
   print("calculating uvws")
   wavelength = 300./freq_MHz
   x_diff_lambda = x_diff/wavelength
   y_diff_lambda = y_diff/wavelength
   z_diff_lambda = z_diff/wavelength
   
   uu_array = (np.sin(hourangle)*x_diff_lambda + np.cos(hourangle)*y_diff_lambda + 0.*z_diff_lambda) #* (1./wavelength)
   vv_array = (-np.sin(declination)*np.cos(hourangle)*x_diff_lambda + np.sin(declination)*np.sin(hourangle)*y_diff_lambda + np.cos(declination)*z_diff_lambda) #* (1./wavelength)
   ww_array = (np.cos(declination)*np.cos(hourangle)*x_diff_lambda - np.cos(declination)*np.sin(hourangle)*y_diff_lambda + np.sin(declination)*z_diff_lambda) #* (1./wavelength)
   return(uu_array,vv_array,ww_array)
   
def calc_x_y_z_diff_from_E_N_U(E_pos_ant1,E_pos_ant2,N_pos_ant1,N_pos_ant2,U_pos_ant1,U_pos_ant2):
   #Look at instruction to go from E,N,U to u,v,w (https://web.njit.edu/~gary/728/Lecture6.html
   
   E_diff = E_pos_ant2 - E_pos_ant1
   N_diff = N_pos_ant2 - N_pos_ant1
   U_diff = U_pos_ant2 - U_pos_ant1
   
   x_diff = -1*N_diff*np.sin(mwa_latitude_rad) + U_diff*np.cos(mwa_latitude_rad)
   y_diff = E_diff 
   z_diff = N_diff*np.cos(mwa_latitude_rad) + U_diff*np.sin(mwa_latitude_rad)
                 
   #x_pos_ant2 = -1*N_pos_ant2*np.sin(mwa_latitude_rad) + U_pos_ant2*np.cos(mwa_latitude_rad)
   #y_pos_ant2 = E_pos_ant2
   #z_pos_ant2 = N_pos_ant2*np.cos(mwa_latitude_rad) + U_pos_ant2*np.sin(mwa_latitude_rad)

   return(x_diff,y_diff,z_diff)
   
def calculate_freq_MHz_fine_chan_subarray(EDA2_chan):
   freq_MHz_centre = (400./512.)*EDA2_chan
   fine_chan_index_array = np.arange(fine_chans_per_EDA2_chan) - 14
   freq_MHz_fine_chan_sub_array = fine_chan_index_array*fine_chan_width_MHz + freq_MHz_centre + (fine_chan_width_MHz/2.)
   return freq_MHz_fine_chan_sub_array
   
def calc_beam_normalisation(freq_MHz_fine_chan,LNA_impedance):
   mu_0 = 4.*np.pi*10**(-7)
   w = 2.*np.pi*freq_MHz_fine_chan*1e6; # angular frequency
   normalisation_factor = (-4.*np.pi*1j / (mu_0 * w)) * LNA_impedance
   return(normalisation_factor)
   
def calc_LNA_impedance_eda2(freq_MHz_fine_chan):
   #See Daniel email from 18/09/2018 on normalisation factor, updated email 25 Nov for EDA2 normalisation calculation
   #For EDA2 the LNA impedance is actually based on RLC lumped circuit model:
   #Series component
   Ls = 2.e-9 #nH
   #Parallel component
   Rp = 914. #ohms
   Lp = 450.e-9 #nH
   Cp = 3.2e-12 #pF
   #So the full impedance is calculated as follows:
   w = 2.*np.pi*freq_MHz_fine_chan*1e6; # angular frequency

   Z = 1j*w*Ls + (1./Rp + (1j*w*Cp) + 1./(1j*w*Lp))**(-1)
   return(Z)
   
def simulate_eda2_with_complex_beams(EDA2_chan,EDA2_chan_index,lst_hrs,EDA2_obs_time,freq_MHz_fine_chan_index,beam_dir='/astro/mwaeor/bmckinley/EoR/EDA2/',git_repo_dir='/astro/mwaeor/bmckinley/code/',nside=32,plot_from_saved=False,sim_unity=True,sim_pt_source=False,check_figs=False,fine_chan=True,base_freq_index_offset=0):

   antenna_layout_filename='/%sant_pos_eda2_combined_on_ground_sim.txt' % git_repo_dir
   
   os.chdir('sims')
   
   freq_MHz_fine_chan_subarray = calculate_freq_MHz_fine_chan_subarray(EDA2_chan)
   freq_MHz_fine_chan = float(freq_MHz_fine_chan_subarray[freq_MHz_fine_chan_index])
   
   test_n_ants = 256
   n_baselines_test = int(test_n_ants*(test_n_ants-1) / 2.)
   print("n_baselines from test %s ants: %s" % (test_n_ants, n_baselines_test))
   n_baselines_test_with_autos = int(test_n_ants*(test_n_ants-1) / 2. + test_n_ants)
   print("n_baselines from test %s ants with autos: %s" % (test_n_ants, n_baselines_test_with_autos))
   
   npix = hp.nside2npix(nside)

   with open(antenna_layout_filename) as f:
      lines = f.readlines()
      n_ants = len(lines) 
   n_baselines = n_ants*(n_ants-1) / 2.
   print("n_baselines from %s ants: %s" % (n_ants, n_baselines))

   E_pos_ant_index_0 = float(lines[0].split()[4])
   N_pos_ant_index_0 = float(lines[0].split()[5])
   U_pos_ant_index_0 = 0.
      
   E_pos_array_ant = np.asarray([ant_string.split()[4] for ant_string in lines[0:test_n_ants]],dtype=float)
   N_pos_array_ant = np.asarray([ant_string.split()[5] for ant_string in lines[0:test_n_ants]],dtype=float)
   U_pos_array_ant = 0. * N_pos_array_ant
   
  
   delta_x_array = E_pos_array_ant - E_pos_ant_index_0 #- np.mean(E_pos_array_ant) #
   delta_y_array = N_pos_array_ant - N_pos_ant_index_0 #- np.mean(N_pos_array_ant) #
   delta_z_array = U_pos_array_ant - U_pos_ant_index_0 # - np.mean(U_pos_array_ant) #

   azimuth_deg_array = np.flip(np.arange(361) + 90.)
   azimuth_deg_array[azimuth_deg_array >= 361] = azimuth_deg_array[azimuth_deg_array >= 361] - 361
   zenith_angle_deg_array = np.arange(91)
   
   #(az = phi)
   repetitions = 91
   az_ang_repeats_array = np.tile(azimuth_deg_array, (repetitions, 1))
   az_ang_repeats_array_rad = az_ang_repeats_array / 180. * np.pi
   az_ang_repeats_array_rad_transpose = np.transpose(az_ang_repeats_array_rad)
   #(zenith_ang = theta)
   repetitions = 361
   zenith_ang_repeats_array = np.tile(zenith_angle_deg_array, (repetitions, 1))
   zenith_ang_repeats_array_rad = zenith_ang_repeats_array / 180. * np.pi
   az_ang_repeats_array_flat_rad = az_ang_repeats_array_rad.flatten()    
   zenith_angle_repeats_array_flat_rad = zenith_ang_repeats_array_rad.flatten('F')
               
   #az_ang_repeats_array_rad_transpose_cos = np.cos(az_ang_repeats_array_rad_transpose)
   az_ang_repeats_array_rad_transpose_cos = np.cos(az_ang_repeats_array_rad_transpose)
   az_ang_repeats_array_rad_transpose_sin = np.sin(az_ang_repeats_array_rad_transpose)
   zenith_ang_repeats_array_rad_sin = np.sin(zenith_ang_repeats_array_rad)

   #hpx angle stuff
   hpx_pix_num_array = np.arange(npix)
   hpx_angles_rad = hp.pix2ang(nside,hpx_pix_num_array) 
   hpx_angles_rad_zenith_angle = hpx_angles_rad[0]
   hpx_angles_rad_azimuth = hpx_angles_rad[1]
               
   #element-by-element multiplication:
   k_x =  np.einsum('ij,ij->ij',az_ang_repeats_array_rad_transpose_cos, zenith_ang_repeats_array_rad_sin)
   k_y =  np.einsum('ij,ij->ij',az_ang_repeats_array_rad_transpose_sin, zenith_ang_repeats_array_rad_sin)
   
   uu_array = np.empty(n_baselines_test)
   vv_array = np.empty(n_baselines_test)
   ww_array = np.empty(n_baselines_test)
   baseline_number_array = np.empty(n_baselines_test)
     
   ######################array job - no for loop
   #simulated EEPs are half a finechan width off
   #this is just a temporary fudge ,  I am going to need to get Daniel to re-interpolate the beams
   number = (freq_MHz_fine_chan+(fine_chan_width_MHz/2.))*1000000
   freq_Hz_string = "%d" % ((np.round(number/5.) * 5.) - (5*EDA2_chan_index))
   wavelength = 300./freq_MHz_fine_chan
   k_0=2.*np.pi / wavelength 
   #adding in geometric delay
   phase_delta = k_0 * (np.einsum('i,jk->jki', delta_x_array, k_x) + np.einsum('i,jk->jki', delta_y_array, k_y))
   
   lst_deg = lst_hrs * 15.
   
   unity_sky_value = 1.
   
   #Do all the gsm sky stuff outside the pol loop
   #dont generate gsm maps here, do on desktop and copy across
   #gsm = GlobalSkyModel(freq_unit='MHz')
   #gsm_map_512 = gsm.generate(freq_MHz_fine_chan)
   
            
   gsm_fits_name = "%sgsm_map_%0.3f_MHz.fits" % (gsm_dir,freq_MHz_fine_chan)
   gsm_map_512 = hp.read_map(gsm_fits_name)
   print("Read %s" % gsm_fits_name)
   
   full_nside = hp.npix2nside(gsm_map_512.shape[0])
   
   gsm_auto_array = np.empty(test_n_ants)
   baseline_length_lambda_array = np.empty(n_baselines_test)
   gsm_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
   gsm_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0])
   
   #do all the rotation stuff on the full res maps ()
   zenith_dec = mwa_latitude_deg
   zenith_ra = lst_deg
   dec_rotate_gsm = zenith_dec-90. 
   ra_rotate_gsm = zenith_ra
   r_gsm_C = hp.Rotator(coord=['G','C'])
   #r_gsm_dec = hp.Rotator(rot=[0,dec_rotate_gsm], deg=True) #, coord=['C', 'C']
   #r_gsm_ra = hp.Rotator(rot=[ra_rotate_gsm,0], deg=True)
   r_gsm_C_ra_dec = hp.Rotator(coord=['G','C'],rot=[ra_rotate_gsm,dec_rotate_gsm], deg=True)
   r_beam = hp.Rotator(rot=[-90,0], deg=True)
      
   if sim_unity:
      unity_auto_array = np.empty(test_n_ants)
      unity_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
      unity_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0])
   if sim_pt_source:
      zenith_point_source_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
      zenith_point_source_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0]) 
      zenith_point_source_auto_array = np.empty(test_n_ants)
   
      #follow Jacks advice - 1 K (Jy?) point source at zenith
      #going to rotate sky not beam. Beam is in spherical coords, not RA/dec, but they are the same if you do 90 deg - dec for theta
      #SO beam centre 'zenith' is actually ra=0,dec=90-mwa_lat
      #see here: http://faraday.uwyo.edu/~admyers/ASTR5160/handouts/51609.pdf 
      
      #zenith_pixel = hp.ang2pix(512,np.radians(90.-(zenith_dec)),np.radians(zenith_ra))
      #point_source_at_zenith_sky_512[zenith_pixel] = 1.
      
      #put the 1 Jy point source at zenith in downgraded
      point_source_at_zenith_sky_nside = hp.ud_grade(gsm_map_512,nside) * 0.
      zenith_pixel = hp.ang2pix(nside,np.radians(90.-(zenith_dec)),np.radians(zenith_ra))
      point_source_at_zenith_sky_nside[zenith_pixel] = 1.0
      
      if check_figs:
         plt.clf()
         map_title=""
         #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         hp.mollview(map=point_source_at_zenith_sky_nside,rot=(0,90,0),title=map_title)
         fig_name="check1_pt_src_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
   
   
      r_pt_src_ra_dec = hp.Rotator(coord=['C'],rot=[ra_rotate_gsm,dec_rotate_gsm], deg=True)
      
      #plt.clf()
      #map_title=""
      ########hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      #hp.mollview(map=hp.ud_grade(rotated_point_source_at_zenith_sky_ra_512,nside),rot=(0,0,0),title=map_title)
      #fig_name="check2ra_pt_src_LST_%0.1f_%s_%0.3f_MHz.png" % (lst_deg,pol,freq_MHz_fine_chan)
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)      
      
      #rotated_point_source_at_zenith_sky_ra_dec_512 = r_gsm_dec.rotate_map(rotated_point_source_at_zenith_sky_ra_512)
      # b_1 = r_beam.rotate_map_alms(beam_1, use_pixel_weights=False)
      rotated_point_source_at_zenith_sky_ra_dec = r_pt_src_ra_dec.rotate_map_alms(point_source_at_zenith_sky_nside,use_pixel_weights=False)
      
      if check_figs:
         plt.clf()
         map_title=""
         ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         hp.mollview(map=rotated_point_source_at_zenith_sky_ra_dec,rot=(0,90,0),title=map_title)
         fig_name="check3ra_dec_pt_src_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)         
      
      #point_source_at_zenith_sky = hp.ud_grade(rotated_point_source_at_zenith_sky_ra_dec,nside)
      point_source_at_zenith_sky = rotated_point_source_at_zenith_sky_ra_dec
   
      if check_figs:
         plt.clf()
         map_title=""
         ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         hp.mollview(map=point_source_at_zenith_sky,rot=(0,90,0),title=map_title)
         fig_name="check4ra_dec_pt_src_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)    
     
      zenith_point_source_repeats_array = np.tile(point_source_at_zenith_sky, (test_n_ants,1))
      zenith_point_source_repeats_array = np.transpose(zenith_point_source_repeats_array)
   
   if check_figs:
      plt.clf()
      map_title=""
      ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      hp.mollview(map=gsm_map_512,coord="G",rot=(0,90,0),title=map_title)
      fig_name="check1_gsm.png" 
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)  
   

   
   #plt.clf()
   #map_title=""
   #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
   #hp.mollview(map=rotated_gsm_C_512,rot=(0,0,0),title=map_title)
   #fig_name="check2_gsm_celestial.png" 
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name)
   
   
   #plt.clf()
   #map_title=""
   #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
   #hp.mollview(map=rotated_gsm_C_ra_512,rot=(0,0,0),title=map_title)
   #fig_name="check3_gsm_rot_ra.png" 
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name)         

   rotated_gsm_C_ra_dec_512 = r_gsm_C_ra_dec.rotate_map_alms(gsm_map_512,use_pixel_weights=False)
   
   if check_figs:
      plt.clf()
      map_title=""
      hp.orthview(map=rotated_gsm_C_ra_dec_512,half_sky=False,rot=(0,90,0),title=map_title)
      #hp.mollview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512,nside),rot=(0,90,0),title=map_title)
      fig_name="check4_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
   
   rotated_gsm_C_ra_dec_512_extra_90 = r_beam.rotate_map_alms(rotated_gsm_C_ra_dec_512,use_pixel_weights=False)
   
   if check_figs:
      plt.clf()
      map_title=""
      hp.orthview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512_extra_90,nside),half_sky=False,rot=(0,90,0),title=map_title)
      ##hp.mollview(map=rotated_gsm_C_ra_dec_512_extra_90,rot=(0,90,0),title=map_title)
      fig_name="check4a_extra_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
   #downgrade res and scale (dont scale - leave as K)
   gsm_map = hp.ud_grade(rotated_gsm_C_ra_dec_512_extra_90,nside) #* scale
   
   if check_figs:
      plt.clf()
      map_title=""
      hp.orthview(map=gsm_map,half_sky=False,rot=(0,90,0),title=map_title)
      #hp.mollview(map=gsm_map,rot=(0,90,0),title=map_title)
      fig_name="check_gsm_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz_fine_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
   
   #make a cube with copies of tsky, dimensions (npix,test_n_ants)
   gsm_repeats_array = np.tile(gsm_map, (test_n_ants,1))
   gsm_repeats_array = np.transpose(gsm_repeats_array)
   
   unity_sky_repeats_array = gsm_repeats_array * 0. + unity_sky_value
   
   for pol_index1,pol1 in enumerate(['X','Y']):  #,'Y'
      if pol1=='X':
         daniel_pol1 = 'Y'
      else:
         daniel_pol1 = 'X'
      for pol_index2,pol2 in enumerate(['X','Y']):
         if pol2=='X':
            daniel_pol2 = 'Y'
         else:
            daniel_pol2 = 'X'         
         uu_array_filename = "uu_array_%0.3f.npy" % (freq_MHz_fine_chan) 
         vv_array_filename = "vv_array_%0.3f.npy" % (freq_MHz_fine_chan)
         ww_array_filename = "ww_array_%0.3f.npy" % (freq_MHz_fine_chan)    
         baseline_number_array_filename = "baseline_number_array_%0.3f.npy" % (freq_MHz_fine_chan)
      
         unity_cross_visibility_real_array_filename = "unity_cross_visibility_real_array_%s%s_%0.3f.npy" % (pol1,pol2,freq_MHz_fine_chan)
         unity_cross_visibility_imag_array_filename = "unity_cross_visibility_imag_array_%s%s_%0.3f.npy" % (pol1,pol2,freq_MHz_fine_chan)
         gsm_cross_visibility_real_array_filename = "gsm_cross_visibility_real_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
         gsm_cross_visibility_imag_array_filename = "gsm_cross_visibility_imag_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
         zenith_point_source_cross_visibility_real_array_filename = "zenith_point_source_cross_visibility_real_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
         zenith_point_source_cross_visibility_imag_array_filename = "zenith_point_source_cross_visibility_imag_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
           
         unity_auto_array_filename = "unity_auto_array_%s_%0.3f.npy" % (pol1,freq_MHz_fine_chan)
         gsm_auto_array_filename = "gsm_auto_array_%s_%s_%0.3f.npy" % (EDA2_obs_time,pol1,freq_MHz_fine_chan)
         zenith_point_source_auto_array_filename = "zenith_point_source_auto_array_%s_%s_%0.3f.npy" % (EDA2_obs_time,pol1,freq_MHz_fine_chan)
         
         baseline_length_lambda_array_filename = "baseline_length_lambda_array_%s_%0.3f.npy" % (pol1,freq_MHz_fine_chan)
   
         if not plot_from_saved:
            if fine_chan:
               EEP_name1 = '%sfreq_interp/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (beam_dir,freq_Hz_string,daniel_pol1)
               EEP_name2 = '%sfreq_interp/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (beam_dir,freq_Hz_string,daniel_pol2)
               #SITARA MATLAB beam file here: /md0/EoR/EDA2/EEPs/SITARA/chall_beam_Y.mat (1 MHz res) (70 - 200 MHz?)
               #beam_data1 = mat73.loadmat(EEP_name1)
               beam_data1 = h5py.File(EEP_name1, 'r')
                
               #keys = beam_data1_file.keys()
               #print(keys)
               print("loaded %s " % EEP_name1)
               
            
               beam_data2 = h5py.File(EEP_name2, 'r')
               print("loaded %s " % EEP_name2)

               E_phi_cube1 = np.transpose(np.array(beam_data1['Ephi']))[:,:,0:test_n_ants]
               E_theta_cube1 = np.transpose(np.array(beam_data1['Etheta']))[:,:,0:test_n_ants]
   
               E_phi_cube2 = np.transpose(np.array(beam_data2['Ephi']))[:,:,0:test_n_ants]
               E_theta_cube2 = np.transpose(np.array(beam_data2['Etheta']))[:,:,0:test_n_ants]              

               E_phi_cube1 = E_phi_cube1['real']+E_phi_cube1['imag']*1j
               E_theta_cube1 = E_theta_cube1['real']+E_theta_cube1['imag']*1j
               E_phi_cube2 = E_phi_cube2['real']+E_phi_cube2['imag']*1j
               E_theta_cube2 = E_theta_cube2['real']+E_theta_cube2['imag']*1j
               
            else:
               EEP_name1 = '%snew_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (beam_dir,freq_Hz_string,daniel_pol1)
               EEP_name2 = '%snew_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (beam_dir,freq_Hz_string,daniel_pol2)
               #SITARA MATLAB beam file here: /md0/EoR/EDA2/EEPs/SITARA/chall_beam_Y.mat (1 MHz res) (70 - 200 MHz?)
               beam_data1 = loadmat(EEP_name1)
               print("loaded %s " % EEP_name1)
               beam_data2 = loadmat(EEP_name2)
               print("loaded %s " % EEP_name2)
               
               E_phi_cube1 = beam_data1['Ephi'][:,:,0:test_n_ants]
               E_theta_cube1 = beam_data1['Etheta'][:,:,0:test_n_ants]
   
               E_phi_cube2 = beam_data2['Ephi'][:,:,0:test_n_ants]
               E_theta_cube2 = beam_data2['Etheta'][:,:,0:test_n_ants]
            
            #Daniels version of normalising the beam i.e. getting in correct unitless form  by getting rid on m squared
            #Doesn't have to be 1 at zenith, just need to make sure use same normalised beam for signal extraction
            
            lna_impedance = calc_LNA_impedance_eda2(freq_MHz_fine_chan)
            beam_norm = calc_beam_normalisation(freq_MHz_fine_chan,lna_impedance)
            print("normalising beam by %0.5f %0.5fj " % (beam_norm.real,beam_norm.imag))
            
            E_phi_cube1 = E_phi_cube1 * beam_norm
            E_theta_cube1 = E_theta_cube1 * beam_norm
            E_phi_cube2 = E_phi_cube2 * beam_norm
            E_theta_cube2 = E_theta_cube2 * beam_norm               
                
            E_phi_cube_plus_delta1 = E_phi_cube1 * np.exp(1j*phase_delta)
            E_theta_cube_plus_delta1 = E_theta_cube1 * np.exp(1j*phase_delta)
            E_phi_cube_plus_delta2 = E_phi_cube2 * np.exp(1j*phase_delta)
            E_theta_cube_plus_delta2 = E_theta_cube2 * np.exp(1j*phase_delta)     

            E_theta_cube1 = E_theta_cube_plus_delta1
            E_phi_cube1 = E_phi_cube_plus_delta1
            E_theta_cube2 = E_theta_cube_plus_delta2
            E_phi_cube2 = E_phi_cube_plus_delta2
                      
            E_phi_cube_slice1 = E_phi_cube1[:,:,0:0+test_n_ants]
            #E_phi_cube_slice = E_phi_cube[:,:,ant_index_1]
            #E_phi_cube_slice_flat = E_phi_cube_slice.flatten('F')
            #see https://stackoverflow.com/questions/18757742/how-to-flatten-only-some-dimensions-of-a-numpy-array
            flattened_size = int(E_phi_cube_slice1.shape[0]*E_phi_cube_slice1.shape[1])
            E_phi_cube_slice1 = E_phi_cube_slice1.transpose([1,0,2])
            E_phi_cube_slice_flat1 = E_phi_cube_slice1.reshape(flattened_size,E_phi_cube_slice1.shape[2])
            #E_phi_cube_slice_flat_1 = E_phi_cube_slice.reshape(-1, E_phi_cube_slice.shape[-1])
            regridded_to_hpx_E_phi_complex1 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_phi_cube_slice_flat1, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
            #do all the inter stuff first, then rotate at end (otherwise end up with hole in the middle)
      
            if (pol1==pol2):
               regridded_to_hpx_E_phi_complex2 = regridded_to_hpx_E_phi_complex1
            else:
               E_phi_cube_slice2 = E_phi_cube2[:,:,0:0+test_n_ants]
               #E_phi_cube_slice = E_phi_cube[:,:,ant_index_1]
               #E_phi_cube_slice_flat = E_phi_cube_slice.flatten('F')
               #see https://stackoverflow.com/questions/18757742/how-to-flatten-only-some-dimensions-of-a-numpy-array
               flattened_size = int(E_phi_cube_slice2.shape[0]*E_phi_cube_slice2.shape[1])
               E_phi_cube_slice2 = E_phi_cube_slice2.transpose([1,0,2])
               E_phi_cube_slice_flat2 = E_phi_cube_slice2.reshape(flattened_size,E_phi_cube_slice2.shape[2])
               #E_phi_cube_slice_flat_1 = E_phi_cube_slice.reshape(-1, E_phi_cube_slice.shape[-1])
               regridded_to_hpx_E_phi_complex2 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_phi_cube_slice_flat2, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
               #do all the inter stuff first, then rotate at end (otherwise end up with hole in the middle)
                        
            
            #repeat for E_theta:
            #E_theta_cube_slice = E_theta_cube[:,:,ant_index_1]
            #E_theta_cube_slice_flat = E_theta_cube_slice.flatten('F')
            E_theta_cube_slice1 = E_theta_cube1[:,:,0:0+test_n_ants]
            E_theta_cube_slice1 = E_theta_cube_slice1.transpose([1,0,2])
            #E_theta_cube_slice_flat = E_theta_cube_slice.reshape(-1, E_theta_cube_slice.shape[-1])
            E_theta_cube_slice_flat1 = E_theta_cube_slice1.reshape(flattened_size,E_theta_cube_slice1.shape[2])
            regridded_to_hpx_E_theta_complex1 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_theta_cube_slice_flat1, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
   
            if (pol1==pol2):
               regridded_to_hpx_E_theta_complex2 = regridded_to_hpx_E_theta_complex1
            else:
               E_theta_cube_slice2 = E_theta_cube2[:,:,0:0+test_n_ants]
               E_theta_cube_slice2 = E_theta_cube_slice2.transpose([1,0,2])
               E_theta_cube_slice_flat2 = E_theta_cube_slice2.reshape(flattened_size,E_theta_cube_slice1.shape[2])
               regridded_to_hpx_E_theta_complex2 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_theta_cube_slice_flat2, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
      
            
            #new array way:
            complex_beam_cube1 = np.empty((npix,2,test_n_ants), dtype=complex)
            complex_beam_cube1[:,0,:] = regridded_to_hpx_E_theta_complex1
            complex_beam_cube1[:,1,:] = regridded_to_hpx_E_phi_complex1
   
            complex_beam_cube2 = np.empty((npix,2,test_n_ants), dtype=complex)
            complex_beam_cube2[:,0,:] = regridded_to_hpx_E_theta_complex2
            complex_beam_cube2[:,1,:] = regridded_to_hpx_E_phi_complex2
                                           

            start_index = 0
            end_index = start_index + (test_n_ants-1)
            for ant_index_1 in range(0,test_n_ants):
               print(ant_index_1)
               #sitara-like baseline ant1,ant2 = (0,69) or 0,10
               #for ant_index_1 in range(68,69):
               #from jishnu via slack:
               #Tsky_12  = np.sum(b_12*(sky_val))/np.sum(np.abs(b_12)) 
               #Tsky_11 = np.abs(np.sum(b_1*(sky_val))/np.sum(np.abs(b_1)))
               start_cube_index = ant_index_1
               end_cube_index = test_n_ants
               power_pattern_cube = np.einsum('ij,ij...->i...', complex_beam_cube1[:,:,ant_index_1], np.conj(complex_beam_cube2[:,:,ant_index_1:test_n_ants]))
               power_pattern_cube = np.nan_to_num(power_pattern_cube)
               
               power_pattern_cube_mag = np.abs(power_pattern_cube)

               if sim_unity:
                  unity_sky_beam_cube = np.einsum('ij,ij->ij',unity_sky_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
                  unity_sky_beam_sum_array = np.einsum('ij->j',unity_sky_beam_cube)
                  unity_visibility_array = unity_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
                  unity_visibility_real_array = np.real(unity_visibility_array)
                  unity_visibility_imag_array = np.imag(unity_visibility_array)
                  unity_cross_visibility_real_array[start_index:end_index] = unity_visibility_real_array[1:]
                  unity_cross_visibility_imag_array[start_index:end_index] = unity_visibility_imag_array[1:]
                  unity_auto_array[ant_index_1] = unity_visibility_real_array[0]
               if sim_pt_source:
                  zenith_point_source_sky_beam_cube = np.einsum('ij,ij->ij',zenith_point_source_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
                  zenith_point_source_sky_beam_sum_array = np.einsum('ij->j',zenith_point_source_sky_beam_cube)
                  zenith_point_source_visibility_array = zenith_point_source_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
                  
                  if check_figs:
                     plt.clf()
                     map_title=""
                     hp.orthview(map=zenith_point_source_sky_beam_cube[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                     #hp.mollview(map=gsm_sky_beam_cube[:,0],rot=(0,90,0),title=map_title)
                     fig_name="check6_pt_source_sky_beam_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz_fine_chan)
                     figmap = plt.gcf()
                     figmap.savefig(fig_name)
                     print("saved %s" % fig_name) 
               
                  zenith_point_source_visibility_real_array = np.real(zenith_point_source_visibility_array)
                  zenith_point_source_visibility_imag_array = np.imag(zenith_point_source_visibility_array)
                  zenith_point_source_cross_visibility_real_array[start_index:end_index] = zenith_point_source_visibility_real_array[1:]
                  zenith_point_source_cross_visibility_imag_array[start_index:end_index] = zenith_point_source_visibility_imag_array[1:]               
                  zenith_point_source_auto_array[ant_index_1] = zenith_point_source_visibility_real_array[0]
                                    
               gsm_sky_beam_cube = np.einsum('ij,ij->ij',gsm_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
               gsm_sky_beam_sum_array = np.einsum('ij->j',gsm_sky_beam_cube)
               
               ######sanity check:
               #plt.clf()
               #map_title=""
               #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
               #hp.mollview(map=power_pattern_cube[:,0],title=map_title)
               #fig_name="check4_pattern_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz_fine_chan)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)    
               
               gsm_visibility_array = gsm_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
               gsm_visibility_real_array = np.real(gsm_visibility_array)
               gsm_visibility_imag_array = np.imag(gsm_visibility_array)
               gsm_cross_visibility_real_array[start_index:end_index] = gsm_visibility_real_array[1:]
               gsm_cross_visibility_imag_array[start_index:end_index] = gsm_visibility_imag_array[1:]
               gsm_auto_array[ant_index_1] = gsm_visibility_real_array[0]
   
   
               #sum_unity_sky_beam = np.nansum(unity_sky_beam)
               #sum_mag_beam = np.nansum(np.abs(rotated_power_pattern))
               #visibility = sum_unity_sky_beam / sum_mag_beam
              
               if check_figs:
                  #####sanity check:
                  plt.clf()
                  map_title=""
                  hp.orthview(map=power_pattern_cube_mag[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                  #hp.mollview(map=power_pattern_cube_mag[:,0],rot=(0,90,0),title=map_title)
                  fig_name="check4_complex_power_pattern_mag_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz_fine_chan)
                  figmap = plt.gcf()
                  figmap.savefig(fig_name)
                  print("saved %s" % fig_name)
               
   
               #rotated_beam = r_beam.rotate_map_alms(power_pattern_cube_mag[:,0],use_pixel_weights=False)
   
               #plt.clf()
               #map_title=""
               #hp.orthview(map=rotated_beam,half_sky=False,rot=(0,90,0),title=map_title)
               ##hp.mollview(map=rotated_beam,rot=(0,90,0),title=map_title)
               #fig_name="check4a_rot_complex_power_pattern_mag_%s_%s_%s_%0.3f_MHz.png" % (ant_index_1,'0',pol,freq_MHz_fine_chan)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)
                              
               if check_figs:
                  #####sanity check:
                  plt.clf()
                  map_title=""
                  hp.orthview(map=gsm_sky_beam_cube[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                  #hp.mollview(map=gsm_sky_beam_cube[:,0],rot=(0,90,0),title=map_title)
                  fig_name="check5_gsm_sky_beam_%0.1f_%s_%s_%s%s_%0.3f_MHz.png" % (lst_deg,ant_index_1,'0',pol1,pol2,freq_MHz_fine_chan)
                  figmap = plt.gcf()
                  figmap.savefig(fig_name)
                  print("saved %s" % fig_name)               
   
               ##phase sanity check?
               #power_pattern_cube_phase = np.angle(power_pattern_cube)
               #power_pattern_cube_phase[power_pattern_cube_phase < -np.pi] = power_pattern_cube_phase[power_pattern_cube_phase < -np.pi] + (2 * np.pi)  
               #power_pattern_cube_phase[power_pattern_cube_phase > np.pi] = power_pattern_cube_phase[power_pattern_cube_phase > np.pi] - (2 * np.pi) 
               #plt.clf()
               #map_title=""
               #hp.orthview(map=power_pattern_cube_phase[:,1],half_sky=True,rot=(0,90,0),title=map_title)
               ##hp.mollview(map=power_pattern_cube_phase[:,1],rot=(0,90,0),title=map_title)
               #fig_name="check4_complex_power_pattern_phase_%s_%s_%s_%0.3f_MHz.png" % (ant_index_1,'1',pol,freq_MHz_fine_chan)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)
                  
                        
               ### Add some text to the png
               #img = Image.open("%s" % fig_name)
               #draw = ImageDraw.Draw(img)
               #font = ImageFont.truetype('FreeSans.ttf',30)
               #draw.text((10, 10),"%0.3f MHz\n  %s %s " % (freq_MHz_fine_chan,ant_name_1,ant_name_2),(0,0,0),font=font)
               #img.save("%s" % fig_name)
   
               E_pos_ant1 = float(lines[ant_index_1].split()[4])
               N_pos_ant1 = float(lines[ant_index_1].split()[5])
               U_pos_ant1 = 0.
               
               E_pos_array_ant2 = np.asarray([ant_string.split()[4] for ant_string in lines[ant_index_1:test_n_ants]],dtype=float)
               N_pos_array_ant2 = np.asarray([ant_string.split()[5] for ant_string in lines[ant_index_1:test_n_ants]],dtype=float)
               U_pos_array_ant2 = N_pos_array_ant2 * 0.
               
               #due to Daniel's beam sims following different conventions for azimuth we need to rotate the E,N coords by 90 degrees
               #alternatively we could rotate the beams themselves (and the sky), but since the beams are cubes this is hard and slow
               
               x_diff,y_diff,z_diff = calc_x_y_z_diff_from_E_N_U(E_pos_ant1,E_pos_array_ant2,N_pos_ant1,N_pos_array_ant2,U_pos_ant1,U_pos_array_ant2)
               
               if (pol_index1==0 and pol_index2==0):
                  #now need delta_x and delta_y (difference in x and x pos for each antenna wrt antenna index 0)
                  uu_sub_array,vv_sub_array,ww_sub_array = calc_uvw(x_diff,y_diff,z_diff,freq_MHz_fine_chan,hourangle=0.,declination=mwa_latitude_rad)
                  #uvw calcs from: https://web.njit.edu/~gary/728/Lecture6.html
                  #see also https://www.atnf.csiro.au/research/radio-school/2014/talks/ATNF2014-Advanced.pdf
                  #also see CASA coord convention doc https://casa.nrao.edu/Memos/CoordConvention.pdf
                  ## Assign x->East and x-North. This is the local geographic csys
                  baseline_number_sub_array = (ant_index_1 * 256) + np.arange(ant_index_1,test_n_ants)
                  baseline_length_lambda = np.sqrt((uu_sub_array)**2 + (vv_sub_array)**2 + (ww_sub_array)**2) #u,v,w already in wavelength units / wavelength
                       

               if (pol_index1==0 and pol_index2==0 and EDA2_chan_index==0):
                  uu_array[start_index:end_index] = uu_sub_array[1:]
                  vv_array[start_index:end_index] = vv_sub_array[1:]
                  ww_array[start_index:end_index] = ww_sub_array[1:]
                  baseline_number_array[start_index:end_index] = baseline_number_sub_array[1:]
                  baseline_length_lambda_array[start_index:end_index] = baseline_length_lambda[1:]
   
               start_index += (test_n_ants-1-ant_index_1)
               end_index += (test_n_ants-2-ant_index_1)
            
            #save the arrays
            if sim_unity:
               np.save(unity_cross_visibility_real_array_filename,unity_cross_visibility_real_array)
               np.save(unity_cross_visibility_imag_array_filename,unity_cross_visibility_imag_array)
               np.save(unity_auto_array_filename,unity_auto_array)
            if sim_pt_source:
               np.save(zenith_point_source_cross_visibility_real_array_filename,zenith_point_source_cross_visibility_real_array)
               np.save(zenith_point_source_cross_visibility_imag_array_filename,zenith_point_source_cross_visibility_imag_array)                      
               np.save(zenith_point_source_auto_array_filename,zenith_point_source_auto_array)
            np.save(gsm_cross_visibility_real_array_filename,gsm_cross_visibility_real_array)
            np.save(gsm_cross_visibility_imag_array_filename,gsm_cross_visibility_imag_array)   
            np.save(gsm_auto_array_filename,gsm_auto_array)
             
            np.save(baseline_length_lambda_array_filename,baseline_length_lambda_array)
            
            if (pol_index1==0 and pol_index2==0):
               #print(uu_array)
               #print(vv_array)
               #print(ww_array)
               np.save(uu_array_filename,uu_array)
               np.save(vv_array_filename,vv_array)
               np.save(ww_array_filename,ww_array)
               np.save(baseline_number_array_filename,baseline_number_array)
               
         else:
            if sim_unity:
               unity_cross_visibility_real_array = np.load(unity_cross_visibility_real_array_filename)
               unity_cross_visibility_imag_array = np.load(unity_cross_visibility_imag_array_filename)
               unity_auto_array = np.load(unity_auto_array_filename)
            if sim_pt_source:
               zenith_point_source_cross_visibility_real_array = np.load(zenith_point_source_cross_visibility_real_array_filename)
               zenith_point_source_cross_visibility_imag_array = np.load(zenith_point_source_cross_visibility_imag_array_filename)                      
               zenith_point_source_auto_array = np.load(zenith_point_source_auto_array_filename)
            gsm_cross_visibility_real_array = np.load(gsm_cross_visibility_real_array_filename)
            gsm_cross_visibility_imag_array = np.load(gsm_cross_visibility_imag_array_filename) 
            gsm_auto_array = np.load(gsm_auto_array_filename)
              
            baseline_length_lambda_array = np.load(baseline_length_lambda_array_filename)
            baseline_number_array = np.load(baseline_number_array_filename)
            uu_array = np.load(uu_array_filename)
            vv_array = np.load(vv_array_filename)
            ww_array = np.load(ww_array_filename)
         
         if check_figs:
            if sim_unity:
               ##plot the real part of the visibilty - unity
               plt.clf()
               map_title="uniform response"
               plt.scatter(baseline_length_lambda_array,unity_cross_visibility_real_array,s=1)
               plt.xlim(0.0, 1)
               plt.ylabel("Real part of vis")
               plt.xlabel("Baseline length (wavelengths)")
               fig_name="unity_response_from_complex_beams_%s%s_%0.3f.png" % (pol1,pol2,freq_MHz_fine_chan)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
      
               ##plot the imaginary part of the visibilty - unity
               plt.clf()
               map_title="uniform response"
               plt.scatter(baseline_length_lambda_array,unity_cross_visibility_imag_array,s=1)
               plt.xlim(0.2, 2)
               plt.ylabel("Imag. part of vis")
               plt.xlabel("Baseline length (wavelengths)")
               fig_name="unity_response_from_complex_beams_imag_%s%s_%0.3f.png" % (pol1,pol2,freq_MHz_fine_chan)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)      
               
               #Autos!
               ant_index_array = range(0,test_n_ants)
               plt.clf()
               map_title="uniform auto response"
               plt.scatter(ant_index_array,unity_auto_array,s=1)
               #plt.xlim(0.2, 2)
               plt.ylabel("Auto vis")
               plt.xlabel("Antenna index")
               fig_name="unity_auto_response_from_complex_beams_%s_%0.3f.png" % (pol1,freq_MHz_fine_chan)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
                 
   
            if sim_pt_source:
               ##plot the real part of the visibilty - point_source
               plt.clf()
               map_title="zenith point source response"
               plt.scatter(baseline_length_lambda_array,zenith_point_source_cross_visibility_real_array,s=1)
               plt.xlim(0.2, 2)
               plt.ylabel("Real part of vis")
               plt.xlabel("Baseline length (wavelengths)")
               fig_name="zenith_point_source_response_from_complex_beams_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
      
               plt.clf()
               map_title="zenith point source response"
               plt.scatter(baseline_length_lambda_array,zenith_point_source_cross_visibility_imag_array,s=1)
               plt.xlim(0.2, 2)
               plt.ylabel("Imag part of vis")
               plt.xlabel("Baseline length (wavelengths)")
               fig_name="zenith_point_source_response_from_complex_beams_imag_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
         
            ##plot the real part of the visibilty 
            plt.clf()
            map_title="gsm response"
            plt.scatter(baseline_length_lambda_array,gsm_cross_visibility_real_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Real part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="gsm_response_from_complex_beams_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
    
            plt.clf()
            map_title="gsm response"
            plt.scatter(baseline_length_lambda_array,gsm_cross_visibility_imag_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Imag part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="gsm_response_from_complex_beams_imag_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz_fine_chan)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
    
            ##
            plt.clf()
            map_title="gsm auto response"
            plt.scatter(ant_index_array,gsm_auto_array,s=1)
            #plt.xlim(0.2, 2)
            plt.ylabel("Auto vis (K)")
            plt.xlabel("Intenna index")
            fig_name="gsm_auto_response_from_complex_beams_%s_%s_%0.3f.png" % (EDA2_obs_time,pol1,freq_MHz_fine_chan)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)  
                            

   os.chdir('..')

if __name__ == "__main__":
    import argparse
    
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)
            
            from argparse import RawTextHelpFormatter
            
    parser = argparse.ArgumentParser(description="Run this to make assassin2 sims and \
            run global signal analysis on Garrawarla with array jobs",formatter_class=SmartFormatter)

    parser.add_argument('--EDA2_chan', default='64',
        help='EDA2 coarse chan. e.g. --eda2_chan="64"')
        
    parser.add_argument('--EDA2_chan_index', default='0',
        help='index of the EDA2 coarse chan in the overall EDA2 chan list (i.e. EDA_chan 64 corresponds to EDA2_chan_index 0). This is needed to find the right beam file. e.g. --eda2_chan_index="0"')
        
    parser.add_argument('--EDA2_obs_time', default='20200303T000000',
        help='EDA2 obs time that corresponds to EDA2 coarse chan. \
             e.g. --eda2_obs_time_list="20200303T000000"')

    parser.add_argument('--lst_hrs', default='6.0',
        help='LST (in hours) that corresponds to the EDA2 obs time. e.g. --lst_hrs="6.0"')

    parser.add_argument('--freq_MHz_fine_chan_index', default=None,
        help='Index of fine channel subarray for this EDA2 chan (int, 0 - 26). e.g. --freq_MHz_fine_chan_index="0"')

    parser.add_argument('--beam_dir', default='/astro/mwaeor/bmckinley/EoR/EDA2/EEPs/',
        help='Directory that beams are stored in e.g. --lst_hrs="/astro/mwaeor/bmckinley/EoR/EDA2/EEPs/"')

    parser.add_argument('--git_repo_dir', default='/astro/mwaeor/bmckinley/code/',
        help='directory where my code is kept --git_repo_dir="/astro/mwaeor/bmckinley/code/"')


    args = parser.parse_args()
    
    if args.EDA2_chan:
       EDA2_chan = float(args.EDA2_chan)
       
    if args.EDA2_chan_index:
       EDA2_chan_index = int(args.EDA2_chan_index)
       
    if args.EDA2_obs_time:
       EDA2_obs_time = args.EDA2_obs_time
       
    if args.lst_hrs:
       lst_hrs = float(args.lst_hrs)

    if args.freq_MHz_fine_chan_index:
       freq_MHz_fine_chan_index = int(args.freq_MHz_fine_chan_index) 

    if args.beam_dir:
       beam_dir = args.beam_dir

    if args.git_repo_dir:
       git_repo_dir = args.git_repo_dir

    ##this is for LST 60.0 deg, will do rest of lsts later
    #year,month,day,hour,min,sec = 2015,11,29,15,40,29 #LST=60 deg
    #time_string = '%d_%02d_%02d_%02d_%02d_%02d' % (year,month,day,hour,min,sec)
   
    simulate_eda2_with_complex_beams(EDA2_chan=EDA2_chan,EDA2_chan_index=EDA2_chan_index,lst_hrs=lst_hrs,EDA2_obs_time=EDA2_obs_time,freq_MHz_fine_chan_index=freq_MHz_fine_chan_index,beam_dir=beam_dir,git_repo_dir=git_repo_dir)  


