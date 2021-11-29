#!/usr/bin/env python3
#ASSASSIN: (All Sky SignAl Short Spacing INterferometer)
#v2 uses complex beam patterns to simulate model vis in ms format

import matplotlib
matplotlib.use('Agg')
import os,sys
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
#from PIL import Image
#from PIL import ImageFont
#from PIL import ImageDraw 
import healpy as hp
from pygdsm import GSMObserver
from pygdsm import GlobalSkyModel
from pygdsm import GlobalSkyModel2016
from datetime import datetime, date
import time



EDA2_data_dir = '/md0/EoR/EDA2/20200303_data/'   #2020 paper
EDA2_chan_list = range(64,127)  #20200303:

def get_eda2_lst(eda_time_string):
   year, month, day, hour, minute, second = eda_time_string[0:4], eda_time_string[4:6],eda_time_string[6:8], eda_time_string[9:11],eda_time_string[11:13],eda_time_string[13:15]

   eda2_observer = ephem.Observer()
   eda2_observer.lon, eda2_observer.lat = mwa_longitude_ephem, mwa_latitude_ephem
   eda2_observer.date = '%s/%s/%s %s:%s:%s' % (year,month,day,hour,minute,second)
   eda2_obs_lst = (eda2_observer.sidereal_time()) 
   print("LST is")
   print(eda2_obs_lst)
   eda2_obs_lst_hrs = eda2_obs_lst / 2 / np.pi * 24.
   
   return(eda2_obs_lst_hrs)

def simulate_eda2_with_complex_beams(freq_MHz_list,lst_hrs,nside=512,antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt',plot_from_saved=False,EDA2_chan='None',EDA2_obs_time='None',n_obs_concat='None'):
   test_n_ants = 256
   n_baselines_test = int(test_n_ants*(test_n_ants-1) / 2.)
   print("n_baselines from test %s ants: %s" % (test_n_ants, n_baselines_test))
   n_baselines_test_with_autos = int(test_n_ants*(test_n_ants-1) / 2. + test_n_ants)
   print("n_baselines from test %s ants with autos: %s" % (test_n_ants, n_baselines_test_with_autos))
   
   npix = hp.nside2npix(nside)
   lst_hrs = float(lst_hrs)
   lst_deg = lst_hrs * 15.
   with open(antenna_layout_filename) as f:
      lines = f.readlines()
      n_ants = len(lines) 
   n_baselines = n_ants*(n_ants-1) / 2.
   print("n_baselines from %s ants: %s" % (n_ants, n_baselines))

   #Jy per pix conversion
   healpix_pixel_area_sr = 4*np.pi/npix
   #antenna coord stuff (for calculating UVWs)
   
   #need to put x,y,z transformations in here too
   
   #check if daniels positions are 1,2 or 4,5 index in this file
   #I think need to do this in the ant loop
   E_pos_ant_index_0 = float(lines[0].split()[4])
   N_pos_ant_index_0 = float(lines[0].split()[5])
   U_pos_ant_index_0 = 0.
      
   E_pos_array_ant = np.asarray([ant_string.split()[4] for ant_string in lines[0:test_n_ants]],dtype=float)
   N_pos_array_ant = np.asarray([ant_string.split()[5] for ant_string in lines[0:test_n_ants]],dtype=float)
   U_pos_array_ant = 0. * N_pos_array_ant
   
   #no this is wrong I think because it is all in sperical coors of az el
   ##delta_x_array,delta_y_array,delta_z_array = calc_x_y_z_diff_from_E_N_U(E_pos_ant_index_0,E_pos_array_ant,N_pos_ant_index_0,N_pos_array_ant,U_pos_ant_index_0,U_pos_array_ant)
   
   ##need to change this to mean position (don't think it matters)?
   
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
   
   ##need to split up into phi and theta components (this must be wrong)
   #k_x_phi = k_0 * az_ang_repeats_array_rad_transpose_cos
   #k_x_theta = k_0 * zenith_ang_repeats_array_rad_sin
   #k_y_phi = k_0 * az_ang_repeats_array_rad_transpose_sin
   #k_y_theta = k_0 * zenith_ang_repeats_array_rad_sin

   #hpx rotate stuff, rotate the complex beams to zenith at the required LST
   dec_rotate = 90. - float(mwa_latitude_ephem)
   ra_rotate = lst_deg
   r_beam_dec = hp.Rotator(rot=[0,dec_rotate], coord=['C', 'C'], deg=True) 
   r_beam_ra = hp.Rotator(rot=[ra_rotate,0], coord=['C', 'C'], deg=True)    
   #see Jishnu's SITARA paper 1, appendix B
   #first get it working for a single baseline, ant index 0 and 1
   #Get the right beams (EEPs from Daniel) and convert (interpolate) to healpix

   uu_array = np.empty(n_baselines_test)
   vv_array = np.empty(n_baselines_test)
   ww_array = np.empty(n_baselines_test)
   baseline_number_array = np.empty(n_baselines_test)
         
   for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
      freq_MHz = float(freq_MHz)
      freq_Hz_string = "%d" % (freq_MHz*1000000)
      wavelength = 300./freq_MHz
      k_0=2.*np.pi / wavelength 
      #adding in geometric delay
      phase_delta = k_0 * (np.einsum('i,jk->jki', delta_x_array, k_x) + np.einsum('i,jk->jki', delta_y_array, k_y))
      
      unity_sky_value = 1. 
      
      #need to either generate gsm map in MJy.sr and convert to Jy/pix, or convert manually:
      scale = (2. * k * 1.0e26 * healpix_pixel_area_sr) / (wavelength**2)
      print("scale map by %s to get to Jy/pix" % scale)
      
      unity_auto_array = np.empty(test_n_ants)
      gsm_auto_array = np.empty(test_n_ants)
      zenith_point_source_auto_array = np.empty(test_n_ants)
      baseline_length_lambda_array = np.empty(n_baselines_test)

      unity_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
      unity_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0])
      gsm_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
      gsm_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0])
      zenith_point_source_cross_visibility_real_array =  np.empty(baseline_length_lambda_array.shape[0])
      zenith_point_source_cross_visibility_imag_array =  np.empty(baseline_length_lambda_array.shape[0]) 
      
      #Do all the gsm sky stuff outside the pol loop
      #Now we have beams, need a sky!
      gsm = GlobalSkyModel(freq_unit='MHz')
      gsm_map_512 = gsm.generate(freq_MHz)
      full_nside = hp.npix2nside(gsm_map_512.shape[0])
      #print(full_nside)
      #point_source_at_zenith_sky_512 = gsm_map_512 * 0.
                 
      #follow Jacks advice - 1 K (Jy?) point source at zenith
      #going to rotate sky not beam. Beam is in spherical coords, not RA/dec, but they are the same if you do 90 deg - dec for theta
      #SO beam centre 'zenith' is actually ra=0,dec=90-mwa_lat
      #see here: http://faraday.uwyo.edu/~admyers/ASTR5160/handouts/51609.pdf
      zenith_dec = mwa_latitude_deg
      zenith_ra = lst_deg 
      
      #zenith_pixel = hp.ang2pix(512,np.radians(90.-(zenith_dec)),np.radians(zenith_ra))
      #point_source_at_zenith_sky_512[zenith_pixel] = 1.
      
      #put the 1 Jy point source at zenith in downgraded
      point_source_at_zenith_sky_nside = hp.ud_grade(gsm_map_512,nside) * 0.
      zenith_pixel = hp.ang2pix(nside,np.radians(90.-(zenith_dec)),np.radians(zenith_ra))
      point_source_at_zenith_sky_nside[zenith_pixel] = 1.0
      
      plt.clf()
      map_title=""
      #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      hp.mollview(map=point_source_at_zenith_sky_nside,rot=(0,90,0),title=map_title)
      fig_name="check1_pt_src_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
      
      #do all the rotation stuff on the full res maps ()
      dec_rotate_gsm = zenith_dec-90. 
      ra_rotate_gsm = zenith_ra
      
      #i've got my coord system messed up here a bit - dec (theta) should go first in the rotation ...
      #https://zonca.dev/2021/03/rotate-maps-healpy.html
      r_gsm_C = hp.Rotator(coord=['G','C'])
      #r_gsm_dec = hp.Rotator(rot=[0,dec_rotate_gsm], deg=True) #, coord=['C', 'C']
      #r_gsm_ra = hp.Rotator(rot=[ra_rotate_gsm,0], deg=True)
      r_gsm_C_ra_dec = hp.Rotator(coord=['G','C'],rot=[ra_rotate_gsm,dec_rotate_gsm], deg=True)
      r_pt_src_ra_dec = hp.Rotator(coord=['C'],rot=[ra_rotate_gsm,dec_rotate_gsm], deg=True)
      r_beam = hp.Rotator(rot=[-90,0], deg=True)
      #jishnu uses this rotator for the beam:
      #r_beam    = hp.Rotator(rot=[90, -90.0], coord=['C', 'C'], deg=True)
      
      #rotated_point_source_at_zenith_sky_ra_512 = r_gsm_ra.rotate_map(point_source_at_zenith_sky_512)      
   
      #plt.clf()
      #map_title=""
      ########hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      #hp.mollview(map=hp.ud_grade(rotated_point_source_at_zenith_sky_ra_512,nside),rot=(0,0,0),title=map_title)
      #fig_name="check2ra_pt_src_LST_%0.1f_%s_%0.3f_MHz.png" % (lst_deg,pol,freq_MHz)
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)      
      
      #rotated_point_source_at_zenith_sky_ra_dec_512 = r_gsm_dec.rotate_map(rotated_point_source_at_zenith_sky_ra_512)
      # b_1 = r_beam.rotate_map_alms(beam_1, use_pixel_weights=False)
      rotated_point_source_at_zenith_sky_ra_dec = r_pt_src_ra_dec.rotate_map_alms(point_source_at_zenith_sky_nside,use_pixel_weights=False)
      
      plt.clf()
      map_title=""
      ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      hp.mollview(map=rotated_point_source_at_zenith_sky_ra_dec,rot=(0,90,0),title=map_title)
      fig_name="check3ra_dec_pt_src_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)         
      
      #point_source_at_zenith_sky = hp.ud_grade(rotated_point_source_at_zenith_sky_ra_dec,nside)
      point_source_at_zenith_sky = rotated_point_source_at_zenith_sky_ra_dec
   
      plt.clf()
      map_title=""
      ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      hp.mollview(map=point_source_at_zenith_sky,rot=(0,90,0),title=map_title)
      fig_name="check4ra_dec_pt_src_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)    
        
   
      
      plt.clf()
      map_title=""
      ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      hp.mollview(map=gsm_map_512,coord="G",rot=(0,90,0),title=map_title)
      fig_name="check1_gsm.png" 
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)  
      
      
      ##convert to celestial coords
      #rotated_gsm_C_512 = r_gsm_C.rotate_map_alms(gsm_map_512,use_pixel_weights=False)
   
      #plt.clf()
      #map_title=""
      #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      #hp.mollview(map=rotated_gsm_C_512,rot=(0,0,0),title=map_title)
      #fig_name="check2_gsm_celestial.png" 
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)
      
      #rotate the sky insted of the beams (cause beams are a cube per ant)
      #rotated_gsm_C_ra_512 = r_gsm_ra.rotate_map_alms(rotated_gsm_C_512,use_pixel_weights=False)
   
      #plt.clf()
      #map_title=""
      #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
      #hp.mollview(map=rotated_gsm_C_ra_512,rot=(0,0,0),title=map_title)
      #fig_name="check3_gsm_rot_ra.png" 
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)         
      
      
      #rotated_gsm_C_ra_dec_512 = r_gsm_dec.rotate_map_alms(rotated_gsm_C_ra_512,use_pixel_weights=False)
      rotated_gsm_C_ra_dec_512 = r_gsm_C_ra_dec.rotate_map_alms(gsm_map_512,use_pixel_weights=False)
   
      plt.clf()
      map_title=""
      hp.orthview(map=rotated_gsm_C_ra_dec_512,half_sky=False,rot=(0,90,0),title=map_title)
      #hp.mollview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512,nside),rot=(0,90,0),title=map_title)
      fig_name="check4_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
      rotated_gsm_C_ra_dec_512_extra_90 = r_beam.rotate_map_alms(rotated_gsm_C_ra_dec_512,use_pixel_weights=False)
      
      plt.clf()
      map_title=""
      hp.orthview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512_extra_90,nside),half_sky=False,rot=(0,90,0),title=map_title)
      ##hp.mollview(map=rotated_gsm_C_ra_dec_512_extra_90,rot=(0,90,0),title=map_title)
      fig_name="check4a_extra_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
      #downgrade res and scale
      gsm_map = hp.ud_grade(rotated_gsm_C_ra_dec_512_extra_90,nside) * scale
      
      plt.clf()
      map_title=""
      hp.orthview(map=gsm_map,half_sky=False,rot=(0,90,0),title=map_title)
      #hp.mollview(map=gsm_map,rot=(0,90,0),title=map_title)
      fig_name="check_gsm_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      
      
      #make a cube with copies of tsky, dimensions (npix,test_n_ants)
      gsm_repeats_array = np.tile(gsm_map, (test_n_ants,1))
      gsm_repeats_array = np.transpose(gsm_repeats_array)
      
      zenith_point_source_repeats_array = np.tile(point_source_at_zenith_sky, (test_n_ants,1))
      zenith_point_source_repeats_array = np.transpose(zenith_point_source_repeats_array)
      
      unity_sky_repeats_array = gsm_repeats_array * 0 + unity_sky_value

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
            uu_array_filename = "uu_array_%0.3f.npy" % (freq_MHz) 
            vv_array_filename = "vv_array_%0.3f.npy" % (freq_MHz)
            ww_array_filename = "ww_array_%0.3f.npy" % (freq_MHz)    
            baseline_number_array_filename = "baseline_number_array_%0.3f.npy" % (freq_MHz)
         
            unity_cross_visibility_real_array_filename = "unity_cross_visibility_real_array_%s%s_%0.3f.npy" % (pol1,pol2,freq_MHz)
            unity_cross_visibility_imag_array_filename = "unity_cross_visibility_imag_array_%s%s_%0.3f.npy" % (pol1,pol2,freq_MHz)
            gsm_cross_visibility_real_array_filename = "gsm_cross_visibility_real_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            gsm_cross_visibility_imag_array_filename = "gsm_cross_visibility_imag_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            zenith_point_source_cross_visibility_real_array_filename = "zenith_point_source_cross_visibility_real_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            zenith_point_source_cross_visibility_imag_array_filename = "zenith_point_source_cross_visibility_imag_array_%s_%s%s_%0.3f.npy" % (EDA2_obs_time,pol1,pol2,freq_MHz)
              
            unity_auto_array_filename = "unity_auto_array_%s_%0.3f.npy" % (pol1,freq_MHz)
            gsm_auto_array_filename = "gsm_auto_array_%s_%s_%0.3f.npy" % (EDA2_obs_time,pol1,freq_MHz)
            zenith_point_source_auto_array_filename = "zenith_point_source_auto_array_%s_%s_%0.3f.npy" % (EDA2_obs_time,pol1,freq_MHz)
            
            baseline_length_lambda_array_filename = "baseline_length_lambda_array_%s_%0.3f.npy" % (pol1,freq_MHz)
   
            if not plot_from_saved:
               EEP_name1 = '/md0/EoR/EDA2/EEPs/new_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (freq_Hz_string,daniel_pol1)
               EEP_name2 = '/md0/EoR/EDA2/EEPs/new_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (freq_Hz_string,daniel_pol2)
               #SITARA MATLAB beam file here: /md0/EoR/EDA2/EEPs/SITARA/chall_beam_Y.mat (1 MHz res) (70 - 200 MHz?)
               beam_data1 = loadmat(EEP_name1)
               print("loaded %s " % EEP_name1)
               beam_data2 = loadmat(EEP_name2)
               print("loaded %s " % EEP_name2)
               
               E_phi_cube1 = beam_data1['Ephi'][:,:,0:test_n_ants]
               E_theta_cube1 = beam_data1['Etheta'][:,:,0:test_n_ants]

               E_phi_cube2 = beam_data2['Ephi'][:,:,0:test_n_ants]
               E_theta_cube2 = beam_data2['Etheta'][:,:,0:test_n_ants]
               
               #maybe this is not right, if the sims were done to normalise a zenith-pointed staion beam of 256 antenna (as for eda-1) then maybe we just need to normalise by sqrt(256)?
               #lna_impedance = calc_LNA_impedance_eda2(freq_MHz)
               #beam_norm = calc_beam_normalisation(freq_MHz,lna_impedance)
               beam_norm = np.sqrt(256)
               print("normalising beam by %0.5f %0.5fj " % (beam_norm.real,beam_norm.imag))
               E_phi_cube1 = E_phi_cube1 * beam_norm
               E_theta_cube1 = E_theta_cube1 * beam_norm
               E_phi_cube2 = E_phi_cube2 * beam_norm
               E_theta_cube2 = E_theta_cube2 * beam_norm               
                   
               #This is code from plotting the AAVS1 antennas phases, they need to be referenced to some antenna or location
               #We need to add in the geometric phase wrt to a single location for all antennas
               #done outside of pol loop
        
               #phase_delta[phase_delta < -np.pi] = phase_delta[phase_delta < -np.pi] + (2 * np.pi)  
               #phase_delta[phase_delta > np.pi] = phase_delta[phase_delta > np.pi] - (2 * np.pi) 
      
               E_phi_cube_plus_delta1 = E_phi_cube1 * np.exp(1j*phase_delta)
               E_theta_cube_plus_delta1 = E_theta_cube1 * np.exp(1j*phase_delta)
               E_phi_cube_plus_delta2 = E_phi_cube2 * np.exp(1j*phase_delta)
               E_theta_cube_plus_delta2 = E_theta_cube2 * np.exp(1j*phase_delta)     
               #phase wrap causing trouble?
               #E_phi_cube_plus_delta_mag = np.abs(E_phi_cube_plus_delta)
               #E_phi_cube_plus_delta_phase = np.angle(E_phi_cube_plus_delta)
               #E_phi_cube_plus_delta_phase[E_phi_cube_plus_delta_phase < -np.pi] = E_phi_cube_plus_delta_phase[E_phi_cube_plus_delta_phase < -np.pi] + (2 * np.pi)  
               #E_phi_cube_plus_delta_phase[E_phi_cube_plus_delta_phase > np.pi] = E_phi_cube_plus_delta_phase[E_phi_cube_plus_delta_phase > np.pi] - (2 * np.pi) 
               #E_phi_cube = E_phi_cube_plus_delta_mag * np.exp(1j*E_phi_cube_plus_delta_phase)
               #         
               #E_theta_cube_plus_delta_mag = np.abs(E_theta_cube_plus_delta)
               #E_theta_cube_plus_delta_phase = np.angle(E_theta_cube_plus_delta)
               #E_theta_cube_plus_delta_phase[E_theta_cube_plus_delta_phase < -np.pi] = E_theta_cube_plus_delta_phase[E_theta_cube_plus_delta_phase < -np.pi] + (2 * np.pi)  
               #E_theta_cube_plus_delta_phase[E_theta_cube_plus_delta_phase > np.pi] = E_theta_cube_plus_delta_phase[E_theta_cube_plus_delta_phase > np.pi] - (2 * np.pi) 
               #E_theta_cube = E_theta_cube_plus_delta_mag * np.exp(1j*E_theta_cube_plus_delta_phase)
               
               E_theta_cube1 = E_theta_cube_plus_delta1
               E_phi_cube1 = E_phi_cube_plus_delta1
               E_theta_cube2 = E_theta_cube_plus_delta2
               E_phi_cube2 = E_phi_cube_plus_delta2
                         
               #now add this phase to both E_theta and E_phi
               
               #phase_delta_1 = k_x * delta_x_array[1] + k_y * delta_y_array[1]
               #phase_delta_1[phase_delta_1 < -np.pi] = phase_delta_1[phase_delta_1 < -np.pi] + (2 * np.pi)  
               #phase_delta_1[phase_delta_1 > np.pi] = phase_delta_1[phase_delta_1 > np.pi] - (2 * np.pi) 
               
               
               ##again I don't think you can split them like this
               #phase_delta_array_phi_1 = delta_x_array[1] * k_x_phi + delta_y_array[1] * k_y_phi
               #phase_delta_array_theta_1 = delta_x_array[1] * k_x_theta + delta_y_array[1] * k_y_theta
               ##keep the phases from wrapping
               #phase_delta_array_phi_1[phase_delta_array_phi_1 < -np.pi] = phase_delta_array_phi_1[phase_delta_array_phi_1 < -np.pi] + (2 * np.pi)  
               #phase_delta_array_phi_1[phase_delta_array_phi_1 > np.pi] = phase_delta_array_phi_1[phase_delta_array_phi_1 > np.pi] - (2 * np.pi) 
               #phase_delta_array_theta_1[phase_delta_array_theta_1 < -np.pi] = phase_delta_array_theta_1[phase_delta_array_theta_1 < -np.pi] + (2 * np.pi)  
               #phase_delta_array_theta_1[phase_delta_array_theta_1 > np.pi] = phase_delta_array_theta_1[phase_delta_array_theta_1 > np.pi] - (2 * np.pi)
               ##now see if we can add this phase delta to E_phi and E_theta
               #E_phi_cube_slice_1 = E_phi_cube[:,:,1]
               #E_phi_cube_slice_1_mag = np.abs(E_phi_cube_slice_1)
               #E_phi_cube_slice_1_phase = np.angle(E_phi_cube_slice_1)
               #E_phi_cube_slice_1_phase_plus_delta = E_phi_cube_slice_1_phase + phase_delta_array_phi_1
               ##keep the phases from wrapping
               #E_phi_cube_slice_1_phase_plus_delta[E_phi_cube_slice_1_phase_plus_delta < -np.pi] = E_phi_cube_slice_1_phase_plus_delta[E_phi_cube_slice_1_phase_plus_delta < -np.pi] + (2 * np.pi)  
               #E_phi_cube_slice_1_phase_plus_delta[E_phi_cube_slice_1_phase_plus_delta > np.pi] = E_phi_cube_slice_1_phase_plus_delta[E_phi_cube_slice_1_phase_plus_delta > np.pi] - (2 * np.pi) 
               #E_theta_cube_slice_1 = E_theta_cube[:,:,1]
               #E_theta_cube_slice_1_mag = np.abs(E_theta_cube_slice_1)
               #E_theta_cube_slice_1_phase = np.angle(E_theta_cube_slice_1)
               #E_theta_cube_slice_1_phase_plus_delta = E_theta_cube_slice_1_phase + phase_delta_array_theta_1
               #keep the phases from wrapping
               #E_theta_cube_slice_1_phase_plus_delta[E_theta_cube_slice_1_phase_plus_delta < -np.pi] = E_theta_cube_slice_1_phase_plus_delta[E_theta_cube_slice_1_phase_plus_delta < -np.pi] + (2 * np.pi)  
               #E_theta_cube_slice_1_phase_plus_delta[E_theta_cube_slice_1_phase_plus_delta > np.pi] = E_theta_cube_slice_1_phase_plus_delta[E_theta_cube_slice_1_phase_plus_delta > np.pi] - (2 * np.pi) 
               #E_phi_cube_slice_1_with_phase_delta = E_phi_cube_slice_1_mag * np.exp(1j*E_phi_cube_slice_1_phase_plus_delta)
               #E_phi_cube[:,:,1] = E_phi_cube_slice_1_with_phase_delta
               #E_theta_cube_slice_1_with_phase_delta = E_theta_cube_slice_1_mag * np.exp(1j*E_theta_cube_slice_1_phase_plus_delta)
               #E_theta_cube[:,:,1] = E_theta_cube_slice_1_with_phase_delta
               
               ###########################
               
               #delta_x = ant_pos_x - ant002_pos_x_daniel
               #delta_y = ant_pos_y - ant002_pos_y_daniel
               
               ##For two antennas, the gains can be separated such that:
               ## g1g2* (e_ant1)T (e_ant2)*
               ##where g1 and g2 are the complex gain solutions and e_ant_1 and e_ant2 are the complex beam vectors with components in e_theta and e_phi
               ##So it is this quantity: (e_ant1)T (e_ant2)* that affects the gain solutions, hence the gain solution phase will correct for variations in this
               ##So we need to calculate (e_ant1)T (e_ant2)* and find the phase of it, call it the beam_gain
               #
               #beam_gain_complex = Etheta_complex_Ant002*Etheta_complex_ant.conjugate() + Ephi_complex_Ant002*Ephi_complex_ant.conjugate()
               #beam_gain_complex_phase_rad = np.angle(beam_gain_complex)
               #phase_ramp_phi_theta = np.zeros([361,91])
               #for phi in range(0,361):
               #   for theta in range(0,91):
               #      phi_rad = phi/180.*np.pi
               #      theta_rad = theta/180.*np.pi
               #      k_x = k_0 * np.cos(phi_rad)* np.sin(theta_rad)
               #      k_y = k_0 * np.sin(phi_rad)* np.sin(theta_rad)
               #      phase_delta = k_x * delta_x + k_y * delta_y
               #      while (phase_delta < -np.pi ):
               #         phase_delta += 2*np.pi
               #      while (phase_delta > np.pi):
               #         phase_delta -= 2*np.pi
               #      phase_ramp_phi_theta[phi,theta] = phase_delta
               #      beam_gain_complex_phase_rad[phi,theta] += phase_delta
               #      if (beam_gain_complex_phase_rad[phi,theta] >= np.pi):
               #         beam_gain_complex_phase_rad[phi,theta] -= 2*np.pi
               #      if (beam_gain_complex_phase_rad[phi,theta] <= -np.pi):
               #         beam_gain_complex_phase_rad[phi,theta] += 2*np.pi 
               #
               #phase_pattern_deg = beam_gain_complex_phase_rad/np.pi*180.
               #phase_ramp_phi_theta_deg = phase_ramp_phi_theta/np.pi*180.
               #np.save(phase_pattern_ant_filename,beam_gain_complex_phase_rad)
               #print "saved %s" % phase_pattern_ant_filename
               
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
               
                       
               ##repeat for phase ramp:
               #phase_delta_slice = phase_delta.transpose([1,0,2])
               #phase_delta_slice_flat = phase_delta_slice.reshape(flattened_size,phase_delta_slice.shape[2])
               #regridded_to_hpx_phase_delta = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), phase_delta_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
      
               #plt.clf()
               #map_title="rotated beam sim"
               #hp.orthview(map=regridded_to_hpx_phase_delta[:,2],half_sky=True,rot=(0,90,0),title=map_title)
               #fig_name="check3_phase_delta_%s_%0.3f_MHz.png" % (pol,freq_MHz)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)      
       
               #rotate appropriately:
               #rotated_E_phi_complex_1 = r_beam.rotate_map(regridded_to_hpx_E_phi_complex_1)      
               #rotated_E_theta_complex_1 = r_beam.rotate_map(regridded_to_hpx_E_theta_complex_1) 
                    
               ###sanity check power pattern:
               #power_pattern_1 = np.abs(regridded_to_hpx_E_phi_complex_1[:,1])**2 + np.abs(regridded_to_hpx_E_theta_complex_1[:,1])**2
               #rotated_power_pattern_1  = r_beam_dec.rotate_map(power_pattern_1)
               
               #plt.clf()
               #map_title="rotated beam sim"
               #hp.orthview(map=rotated_power_pattern_1,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
               #fig_name="check1_complex_power_pattern_%s_%s_%0.3f_MHz.png" % (ant_index_1,pol,freq_MHz)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)
               
               #Add some text to the png
               #img = Image.open("%s" % fig_name)
               #draw = ImageDraw.Draw(img)
               #font = ImageFont.truetype('FreeSans.ttf',30)
               #draw.text((10, 10),"%0.3f MHz\n  %s " % (freq_MHz,ant_name_1),(0,0,0),font=font)
               #img.save("%s" % fig_name)
               
               #old matrix way:
               #complex_beam_1 = np.matrix(np.empty((npix,2), dtype=complex))
               #complex_beam_1[:,0] = np.matrix(regridded_to_hpx_E_theta_complex_1).transpose()
               #complex_beam_1[:,1] = np.matrix(regridded_to_hpx_E_phi_complex_1).transpose()
               
               #new array way:
               complex_beam_cube1 = np.empty((npix,2,test_n_ants), dtype=complex)
               complex_beam_cube1[:,0,:] = regridded_to_hpx_E_theta_complex1
               complex_beam_cube1[:,1,:] = regridded_to_hpx_E_phi_complex1

               complex_beam_cube2 = np.empty((npix,2,test_n_ants), dtype=complex)
               complex_beam_cube2[:,0,:] = regridded_to_hpx_E_theta_complex2
               complex_beam_cube2[:,1,:] = regridded_to_hpx_E_phi_complex2
                                              
               ##DOn't need to repeat separately for ant 2 as it is all done in the cube abaove!
               ##repeat for ant 2
               #E_phi_cube_slice_2 = E_phi_cube[:,:,ant_index_2]
               #E_phi_cube_slice_flat_2 = E_phi_cube_slice_2.flatten('F')
               #regridded_to_hpx_E_phi_complex_2 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_phi_cube_slice_flat_2, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
               #      
               ##repeat for E_theta:
               #E_theta_cube_slice_2 = E_theta_cube[:,:,ant_index_2]
               #E_theta_cube_slice_flat_2 = E_theta_cube_slice_2.flatten('F')
               #regridded_to_hpx_E_theta_complex_2 = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_theta_cube_slice_flat_2, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
               
               #complex_beam_2 = np.matrix(np.empty((npix,2), dtype=complex))
               #complex_beam_2[:,0] = np.matrix(regridded_to_hpx_E_theta_complex_2).transpose()
               #complex_beam_2[:,1] = np.matrix(regridded_to_hpx_E_phi_complex_2).transpose()
               
               #complex_beam_cube_H = np.conj(complex_beam_cube.transpose([1,0,2]))
               
               #construct the power pattern using the diagonal as in SITARA 1 appendix B
               
               ###sanity check power pattern:
               #power_pattern = np.abs(np.array(complex_beam_1[:,0]))**2 + np.abs(np.array(complex_beam_2[:,1]))**2
               #power_pattern = power_pattern[:,0]
               #rotated_power_pattern = r_beam.rotate_map(power_pattern)
               #plt.clf()
               #map_title=""
               #hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
               #fig_name="check3_complex_power_pattern_%s_%s_%0.3f_MHz.png" % (ant_index_2,pol,freq_MHz)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)
               
               ### Add some text to the png
               #img = Image.open("%s" % fig_name)
               #draw = ImageDraw.Draw(img)
               #font = ImageFont.truetype('FreeSans.ttf',30)
               #draw.text((10, 10),"%0.3f MHz\n  %s " % (freq_MHz,ant_name_2),(0,0,0),font=font)
               #img.save("%s" % fig_name)
               
               #forget this costly matrix multiply stuff, do with arrays and only compute the terms that end up in the diagonal
               #https://stackoverflow.com/questions/14758283/is-there-a-numpy-scipy-dot-product-calculating-only-the-diagonal-entries-of-the
               #beam_matmul = np.matmul(complex_beam_1,complex_beam_2_H)
               #power_pattern = np.asarray(beam_matmul.diagonal())[0,:]
               
               #try for one slice to start
               #beam_matmul = np.dot(complex_beam_cube[:,:,0],complex_beam_cube_H[:,:,0])
               #beam_matmul = np.einsum('ij,kj',complex_beam_cube[:,:,0],np.conj(complex_beam_cube[:,:,0]))
               #gonna have to learn python einsum: https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
               #https://ajcr.net/Basic-guide-to-einsum/
               #https://stackoverflow.com/questions/64533775/calculation-diagonal-elements-of-matrices-multiplication-in-numpy-package
               #pretty sure the above step and the one below can be achieved in one go with einsum
               #power_pattern = np.einsum('ii->i', beam_matmul)
         

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
                  
                  #We cant normalise the beams like this as the cross pol (xy yx) beams will have different magnitude to the xx yy
                  #power_pattern_cube = power_pattern_cube / np.max(np.abs(power_pattern_cube))
                  
                  ################
                  ##add in the phase delta to ant index 1
                  #power_pattern_cube_mag = np.abs(power_pattern_cube)
                  #power_pattern_cube_phase = np.angle(power_pattern_cube)
                   
                  #power_pattern_cube_phase_plus_delta = power_pattern_cube_phase + regridded_to_hpx_phase_delta
                  
                  ##keep the phases from wrapping
                  #power_pattern_cube_phase_plus_delta[power_pattern_cube_phase_plus_delta < -np.pi] = power_pattern_cube_phase_plus_delta[power_pattern_cube_phase_plus_delta < -np.pi] + (2 * np.pi)  
                  #power_pattern_cube_phase_plus_delta[power_pattern_cube_phase_plus_delta > np.pi] = power_pattern_cube_phase_plus_delta[power_pattern_cube_phase_plus_delta > np.pi] - (2 * np.pi) 
                  #power_pattern_cube_slice_with_phase_delta = power_pattern_cube_mag * np.exp(1j*power_pattern_cube_phase_plus_delta)
                  #power_pattern_cube = power_pattern_cube_slice_with_phase_delta
                  #power_pattern_cube = np.nan_to_num(power_pattern_cube)
                  
                  #May need to rotate the sky map instead of the beams, I dont think healpy can rotate a cube, and 
                  #it makes sense anyway from an efficiency point of view since there is just one sky (ignoring frequency) but many beams...
                  #rotated_power_pattern = r_beam_dec.rotate_map_alms(power_pattern_cube[:,0],use_pixel_weights=False)
                  #rotated_power_pattern = r_beam_ra.rotate_map_alms(rotated_power_pattern,use_pixel_weights=False)
                  unity_sky_beam_cube = np.einsum('ij,ij->ij',unity_sky_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
                  gsm_sky_beam_cube = np.einsum('ij,ij->ij',gsm_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
                  zenith_point_source_sky_beam_cube = np.einsum('ij,ij->ij',zenith_point_source_repeats_array[:,ant_index_1:test_n_ants], power_pattern_cube)
               
                  ######sanity check:
                  #plt.clf()
                  #map_title=""
                  #######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
                  #hp.mollview(map=power_pattern_cube[:,0],title=map_title)
                  #fig_name="check4_pattern_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz)
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)    
                  
                  
                  unity_sky_beam_sum_array = np.einsum('ij->j',unity_sky_beam_cube)
                  gsm_sky_beam_sum_array = np.einsum('ij->j',gsm_sky_beam_cube)
                  zenith_point_source_sky_beam_sum_array = np.einsum('ij->j',zenith_point_source_sky_beam_cube)
                  
                  power_pattern_cube_mag = np.abs(power_pattern_cube)
                  
                  #what if we don't normalise for the beam weights? - don't if input hpx map is in Jy/pix
                  power_pattern_cube_mag_sum_array = np.einsum('ij->j',power_pattern_cube_mag)
                  unity_visibility_array = unity_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
                  gsm_visibility_array = gsm_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
                  zenith_point_source_visibility_array = zenith_point_source_sky_beam_sum_array #/ power_pattern_cube_mag_sum_array
   
                  #print(zenith_point_source_visibility_array)
   
                  #sum_unity_sky_beam = np.nansum(unity_sky_beam)
                  #sum_mag_beam = np.nansum(np.abs(rotated_power_pattern))
                  #visibility = sum_unity_sky_beam / sum_mag_beam
   
   
                  #####sanity check:
                  plt.clf()
                  map_title=""
                  hp.orthview(map=power_pattern_cube_mag[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                  #hp.mollview(map=power_pattern_cube_mag[:,0],rot=(0,90,0),title=map_title)
                  fig_name="check4_complex_power_pattern_mag_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz)
                  figmap = plt.gcf()
                  figmap.savefig(fig_name)
                  print("saved %s" % fig_name)
                  

                  #rotated_beam = r_beam.rotate_map_alms(power_pattern_cube_mag[:,0],use_pixel_weights=False)
   
                  #plt.clf()
                  #map_title=""
                  #hp.orthview(map=rotated_beam,half_sky=False,rot=(0,90,0),title=map_title)
                  ##hp.mollview(map=rotated_beam,rot=(0,90,0),title=map_title)
                  #fig_name="check4a_rot_complex_power_pattern_mag_%s_%s_%s_%0.3f_MHz.png" % (ant_index_1,'0',pol,freq_MHz)
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)
                                 
   
                  #####sanity check:
                  plt.clf()
                  map_title=""
                  hp.orthview(map=gsm_sky_beam_cube[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                  #hp.mollview(map=gsm_sky_beam_cube[:,0],rot=(0,90,0),title=map_title)
                  fig_name="check5_gsm_sky_beam_%0.1f_%s_%s_%s%s_%0.3f_MHz.png" % (lst_deg,ant_index_1,'0',pol1,pol2,freq_MHz)
                  figmap = plt.gcf()
                  figmap.savefig(fig_name)
                  print("saved %s" % fig_name)               
   
                  plt.clf()
                  map_title=""
                  hp.orthview(map=zenith_point_source_sky_beam_cube[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                  #hp.mollview(map=gsm_sky_beam_cube[:,0],rot=(0,90,0),title=map_title)
                  fig_name="check6_pt_source_sky_beam_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz)
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
                  #fig_name="check4_complex_power_pattern_phase_%s_%s_%s_%0.3f_MHz.png" % (ant_index_1,'1',pol,freq_MHz)
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)
                     
                           
                  ### Add some text to the png
                  #img = Image.open("%s" % fig_name)
                  #draw = ImageDraw.Draw(img)
                  #font = ImageFont.truetype('FreeSans.ttf',30)
                  #draw.text((10, 10),"%0.3f MHz\n  %s %s " % (freq_MHz,ant_name_1,ant_name_2),(0,0,0),font=font)
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
                     uu_sub_array,vv_sub_array,ww_sub_array = calc_uvw(x_diff,y_diff,z_diff,freq_MHz,hourangle=0.,declination=mwa_latitude_rad)
                     #uvw calcs from: https://web.njit.edu/~gary/728/Lecture6.html
                     #see also https://www.atnf.csiro.au/research/radio-school/2014/talks/ATNF2014-Advanced.pdf
                     #also see CASA coord convention doc https://casa.nrao.edu/Memos/CoordConvention.pdf
                     ## Assign x->East and x-North. This is the local geographic csys
                  
                     baseline_number_sub_array = (ant_index_1 * 256) + np.arange(ant_index_1,test_n_ants)
                  
                     baseline_length_lambda = np.sqrt((uu_sub_array)**2 + (vv_sub_array)**2 + (ww_sub_array)**2) #u,v,w already in wavelength units / wavelength
                          
                  #print(baseline_length_array[10])
                  #print(baseline_length_array[68])
                  
                  #need to add the tracking phase to phase the data to zenith
                  #gsm_visibility_array_phase_zenith = gsm_visibility_array * np.exp(2*np.pi*1j*ww_sub_array)
                  #unity_visibility_array_phase_zenith = unity_visibility_array * np.exp(2*np.pi*1j*ww_sub_array)
                  
                  ##print(visibility_list)
                  unity_visibility_real_array = np.real(unity_visibility_array)
                  unity_visibility_imag_array = np.imag(unity_visibility_array)
                  gsm_visibility_real_array = np.real(gsm_visibility_array)
                  gsm_visibility_imag_array = np.imag(gsm_visibility_array)
                  zenith_point_source_visibility_real_array = np.real(zenith_point_source_visibility_array)
                  zenith_point_source_visibility_imag_array = np.imag(zenith_point_source_visibility_array)
                                 
                  #already in wavelength units
                  ###baseline_length_lambda = baseline_length_array / wavelength
                  ##plot the real part of the visibilty
                  #plt.clf()
                  #map_title="uniform response"
                  #plt.scatter(baseline_length_lambda_array,visibility_real_array)
                  #plt.ylabel("Real part of vis")
                  #plt.xlabel("Baseline length (wavelengths)")
                  #fig_name="unity_response_from_complex_beams.png"
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)
                  
                  #print(visibility_real_array[10])
                  #baseline_index_10_unity_vis_list.append(visibility_real_array[10])
                  #visibility_real_array = np.roll(visibility_real_array,ant_index_1)
                  #add in visibility list stuff here:
                  unity_cross_visibility_real_array[start_index:end_index] = unity_visibility_real_array[1:]
                  unity_cross_visibility_imag_array[start_index:end_index] = unity_visibility_imag_array[1:]
                  gsm_cross_visibility_real_array[start_index:end_index] = gsm_visibility_real_array[1:]
                  gsm_cross_visibility_imag_array[start_index:end_index] = gsm_visibility_imag_array[1:]
                  zenith_point_source_cross_visibility_real_array[start_index:end_index] = zenith_point_source_visibility_real_array[1:]
                  zenith_point_source_cross_visibility_imag_array[start_index:end_index] = zenith_point_source_visibility_imag_array[1:]               
                  
                  if (pol_index1==0 and pol_index2==0 and freq_MHz_index==0):
                     uu_array[start_index:end_index] = uu_sub_array[1:]
                     vv_array[start_index:end_index] = vv_sub_array[1:]
                     ww_array[start_index:end_index] = ww_sub_array[1:]
                     baseline_number_array[start_index:end_index] = baseline_number_sub_array[1:]
                     baseline_length_lambda_array[start_index:end_index] = baseline_length_lambda[1:]
   
                  unity_auto_array[ant_index_1] = unity_visibility_real_array[0]
                  gsm_auto_array[ant_index_1] = gsm_visibility_real_array[0]
                  zenith_point_source_auto_array[ant_index_1] = zenith_point_source_visibility_real_array[0]
                  
               
                  start_index += (test_n_ants-1-ant_index_1)
                  end_index += (test_n_ants-2-ant_index_1)
               
               #save the arrays
               np.save(unity_cross_visibility_real_array_filename,unity_cross_visibility_real_array)
               np.save(unity_cross_visibility_imag_array_filename,unity_cross_visibility_imag_array)
               np.save(gsm_cross_visibility_real_array_filename,gsm_cross_visibility_real_array)
               np.save(gsm_cross_visibility_imag_array_filename,gsm_cross_visibility_imag_array)    
               np.save(zenith_point_source_cross_visibility_real_array_filename,zenith_point_source_cross_visibility_real_array)
               np.save(zenith_point_source_cross_visibility_imag_array_filename,zenith_point_source_cross_visibility_imag_array)                      
               np.save(baseline_length_lambda_array_filename,baseline_length_lambda_array)
               np.save(unity_auto_array_filename,unity_auto_array)
               np.save(gsm_auto_array_filename,gsm_auto_array)
               np.save(zenith_point_source_auto_array_filename,zenith_point_source_auto_array)
               if (pol_index1==0 and pol_index2==0 and freq_MHz_index==0):
                  #print(uu_array)
                  #print(vv_array)
                  #print(ww_array)
                  np.save(uu_array_filename,uu_array)
                  np.save(vv_array_filename,vv_array)
                  np.save(ww_array_filename,ww_array)
                  np.save(baseline_number_array_filename,baseline_number_array)
                  
            else:
               unity_cross_visibility_real_array = np.load(unity_cross_visibility_real_array_filename)
               unity_cross_visibility_imag_array = np.load(unity_cross_visibility_imag_array_filename)
               gsm_cross_visibility_real_array = np.load(gsm_cross_visibility_real_array_filename)
               gsm_cross_visibility_imag_array = np.load(gsm_cross_visibility_imag_array_filename)   
               zenith_point_source_cross_visibility_real_array = np.load(zenith_point_source_cross_visibility_real_array_filename)
               zenith_point_source_cross_visibility_imag_array = np.load(zenith_point_source_cross_visibility_imag_array_filename)                      
               baseline_length_lambda_array = np.load(baseline_length_lambda_array_filename)
               unity_auto_array = np.load(unity_auto_array_filename)
               gsm_auto_array = np.load(gsm_auto_array_filename)
               zenith_point_source_auto_array = np.load(zenith_point_source_auto_array_filename)
               baseline_number_array = np.load(baseline_number_array_filename)
               uu_array = np.load(uu_array_filename)
               vv_array = np.load(vv_array_filename)
               ww_array = np.load(ww_array_filename)
            
            ##plot the real part of the visibilty - unity
            plt.clf()
            map_title="uniform response"
            plt.scatter(baseline_length_lambda_array,unity_cross_visibility_real_array,s=1)
            plt.xlim(0.0, 1)
            plt.ylabel("Real part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="unity_response_from_complex_beams_%s%s_%0.3f.png" % (pol1,pol2,freq_MHz)
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
            fig_name="unity_response_from_complex_beams_imag_%s%s_%0.3f.png" % (pol1,pol2,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)        
   
            ##plot the real part of the visibilty - unity
            plt.clf()
            map_title="gsm response"
            plt.scatter(baseline_length_lambda_array,gsm_cross_visibility_real_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Real part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="gsm_response_from_complex_beams_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            plt.clf()
            map_title="gsm response"
            plt.scatter(baseline_length_lambda_array,gsm_cross_visibility_imag_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Imag part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="gsm_response_from_complex_beams_imag_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            ##plot the real part of the visibilty - point_source
            plt.clf()
            map_title="zenith point source response"
            plt.scatter(baseline_length_lambda_array,zenith_point_source_cross_visibility_real_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Real part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="zenith_point_source_response_from_complex_beams_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            plt.clf()
            map_title="zenith point source response"
            plt.scatter(baseline_length_lambda_array,zenith_point_source_cross_visibility_imag_array,s=1)
            plt.xlim(0.2, 2)
            plt.ylabel("Imag part of vis")
            plt.xlabel("Baseline length (wavelengths)")
            fig_name="zenith_point_source_response_from_complex_beams_imag_%s_%s%s_%0.3f.png" % (EDA2_obs_time,pol1,pol2,freq_MHz)
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
            fig_name="unity_auto_response_from_complex_beams_%s_%0.3f.png" % (pol1,freq_MHz)
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
            fig_name="gsm_auto_response_from_complex_beams_%s_%s_%0.3f.png" % (EDA2_obs_time,pol1,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)  
                            
            #now to extract global sky temp from gsm sim vis by comparing to unity vis - but
            #since we now use the complex beam pattern, we don't expect unity sky response to 
            #be purely real - not sure how to deal with that!
            #Then to making these vis into uvfits/miriad format, or at least phased properly (to zenith?)
            
            ##Lets look at some data at the same LST
            #uvfits_filename = "/md0/EoR/EDA2/20200303_data/%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            #print("%s" % uvfits_filename)
            #hdulist = fits.open(uvfits_filename)
            #uvtable = hdulist[0].data
            #uvtable_header = hdulist[0].header
            #hdulist.close()
            #visibilities = uvtable['DATA']
            #visibilities_shape = visibilities.shape
            #print("visibilities_shape")
            #print(visibilities_shape)

def write_to_miriad_vis(freq_MHz_list,lst_hrs,EDA2_chan='None',EDA2_obs_time='None',n_obs_concat='None',antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt',pt_source=False):
   lst_hrs = float(lst_hrs)
   lst_deg = lst_hrs * 15.
   with open(antenna_layout_filename) as f:
      lines = f.readlines()
      n_ants = len(lines) 
   n_baselines = n_ants*(n_ants-1) / 2.
   print("read %s" % antenna_layout_filename)
   print("n_baselines from %s ants: %s" % (n_ants, n_baselines))
   
   #time stuff
   year, month, day, hour, minute, second = EDA2_obs_time[0:4], EDA2_obs_time[4:6],EDA2_obs_time[6:8], EDA2_obs_time[9:11],EDA2_obs_time[11:13],EDA2_obs_time[13:15]
   eda2_astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (float(year), float(month), float(day), float(hour), float(minute), float(second))
   print(eda2_astropy_time_string)
   eda2_astropy_time = Time(eda2_astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
  
   
   for freq_MHz in freq_MHz_list:
      wavelength = 300./freq_MHz
      
      #initiate the miriad uv file outside the pol loop
      print("initiate the miriad uv file outside the pol loop")
      print("writing to miriad file")
      NFFT = 1
      SITE_LON = 116.670575
      FREQ_GHz = np.array([freq_MHz/1000.])#np.linspace(0.000, 0.250, NFFT);
      
      #timestr = time.strftime("%Y%m%d-%H%M%S")
      
      mir_file = "%s_%0.3f.vis" % (EDA2_obs_time,freq_MHz)
      mir_file_uvfits_name = "%s_%0.3f.uvfits" % (EDA2_obs_time,freq_MHz)
      mir_file_ms_name = "%s_%0.3f.ms" % (EDA2_obs_time,freq_MHz)
      # mir_file = "test.vis"
      cmd = "rm -rf %s" % mir_file
      print(cmd)
      os.system(cmd)
      
      print ("Writing out "+ mir_file)
      uv = a.miriad.UV(mir_file, 'new')
      uv['history'] = 'test file\n'
      
      uv.add_var('latitud','d')
      uv.add_var('npol','i')
      uv.add_var('nspect', 'i')
      uv.add_var('obsdec', 'd')
      uv.add_var('vsource', 'r')
      uv.add_var('ischan', 'i')
      uv.add_var('operator', 'a')
      uv.add_var('nants', 'i')
      uv.add_var('baseline', 'r')
      uv.add_var('sfreq', 'd')
      uv.add_var('inttime', 'r')
      uv.add_var('source', 'a')
      uv.add_var('epoch', 'r')
      uv.add_var('version', 'a')
      uv.add_var('ra', 'd')
      uv.add_var('restfreq', 'd')
      uv.add_var('nschan', 'i')
      uv.add_var('sdf', 'd')
      uv.add_var('corr', 'r')
      uv.add_var('freq', 'd')
      uv.add_var('longitu', 'd')
      uv.add_var('nchan', 'i')
      uv.add_var('tscale', 'r')
      uv.add_var('antpos', 'd')
      uv.add_var('telescop', 'a')
      uv.add_var('pol', 'i')
      uv.add_var('coord', 'd')
      uv.add_var('veldop', 'r')
      uv.add_var('lst', 'd')
      uv.add_var('time', 'd')
      uv.add_var('dec', 'd')
      uv.add_var('obsra', 'd')
      uv.add_var('jyperk', 'r')
      uv.add_var('systemp', 'r')
      
      uv['latitud'] = mwa_latitude_deg*np.pi/180.0
      uv['npol'] = 4
      uv['nspect'] = 1
      uv['obsdec'] = mwa_latitude_deg*np.pi/180.0 
      uv['vsource'] = 0.0
      uv['ischan'] = 0
      uv['operator'] = 'J'
      uv['nants'] = 256
      uv['baseline'] = 0.0
      uv['sfreq'] = FREQ_GHz[0]
      uv['inttime'] = 1.0
      uv['source'] = 'zenith'
      uv['epoch'] = 2000.0
      uv['version'] = 'A'
      uv['ra'] = lst_deg*np.pi/180.0
      uv['restfreq'] = 0.0
      uv['nschan'] = NFFT
      uv['sdf'] = 28/1000000. #FREQ_GHz[1]-FREQ_GHz[0]
      uv['corr'] = 0.0 
      uv['freq'] = FREQ_GHz[0]
      uv['longitu'] = mwa_longitude_deg*np.pi/180.0
      uv['nchan'] = NFFT
      uv['tscale'] = 0.0
      uv['antpos'] = 0.0
      uv['telescop'] = 'EDA2_sim'
      uv['pol'] = np.array([-5,-6,-7,-8])   #-5 is xx, -6 yy, -7 xy 8 yx
      uv['coord'] = 0.0
      uv['veldop'] = 0.0
      uv['lst'] = 0.0
      uv['time'] = 0.0
      uv['dec'] = mwa_latitude_deg*np.pi/180.0
      uv['obsra'] = lst_deg*np.pi/180.0
      uv['jyperk'] = 1.0
      uv['systemp'] = 1.0
      
      
      #need to set the x,y,z positions of the antennas properly and then use the generate uvw to get the uvws to write.
      beam = a.phs.Beam(FREQ_GHz)
      ants = []
      for ant_string in lines:
         E_pos = float(ant_string.split()[4])
         N_pos = float(ant_string.split()[5])
         U_pos = 0.
         x_pos,y_pos,z_pos = calc_x_y_z_pos_from_E_N_U(E_pos,N_pos,U_pos)
         ants.append(a.phs.Antenna(x_pos,y_pos,z_pos,beam,delay=0,offset=0))
      
      aa = a.phs.AntennaArray(ants=ants,location=("116:40:14.07","-26:42:11.23"))
      aa.set_jultime(eda2_astropy_time.jd)
      uv['lst'] = float(aa.sidereal_time())
      uv['antpos'] = np.arange(256*3,dtype='d')
      
      data_mask = np.zeros(NFFT)
      
      gsm_cross_visibility_real_array_filename_XX = "gsm_cross_visibility_real_array_%s_XX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_cross_visibility_imag_array_filename_XX = "gsm_cross_visibility_imag_array_%s_XX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_real_array_filename_XX = "zenith_point_source_cross_visibility_real_array_%s_XX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_imag_array_filename_XX = "zenith_point_source_cross_visibility_imag_array_%s_XX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_auto_array_filename_X = "gsm_auto_array_%s_X_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_auto_array_filename_XX = "zenith_point_source_auto_array_%s_XX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)

      gsm_cross_visibility_real_array_filename_YY = "gsm_cross_visibility_real_array_%s_YY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_cross_visibility_imag_array_filename_YY = "gsm_cross_visibility_imag_array_%s_YY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_real_array_filename_YY = "zenith_point_source_cross_visibility_real_array_%s_YY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_imag_array_filename_YY = "zenith_point_source_cross_visibility_imag_array_%s_YY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_auto_array_filename_Y = "gsm_auto_array_%s_Y_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_auto_array_filename_YY = "zenith_point_source_auto_array_%s_YY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)

      gsm_cross_visibility_real_array_filename_XY= "gsm_cross_visibility_real_array_%s_XY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_cross_visibility_imag_array_filename_XY = "gsm_cross_visibility_imag_array_%s_XY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_real_array_filename_XY = "zenith_point_source_cross_visibility_real_array_%s_XY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_imag_array_filename_XY = "zenith_point_source_cross_visibility_imag_array_%s_XY_%0.3f.npy" % (EDA2_obs_time,freq_MHz)

      gsm_cross_visibility_real_array_filename_YX = "gsm_cross_visibility_real_array_%s_YX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      gsm_cross_visibility_imag_array_filename_YX = "gsm_cross_visibility_imag_array_%s_YX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_real_array_filename_YX = "zenith_point_source_cross_visibility_real_array_%s_YX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
      zenith_point_source_cross_visibility_imag_array_filename_YX = "zenith_point_source_cross_visibility_imag_array_%s_YX_%0.3f.npy" % (EDA2_obs_time,freq_MHz)
     
      uu_array_filename = "uu_array_%0.3f.npy" % (freq_MHz) 
      vv_array_filename = "vv_array_%0.3f.npy" % (freq_MHz) 
      ww_array_filename = "ww_array_%0.3f.npy" % (freq_MHz) 
      baseline_number_array_filename = "baseline_number_array_%0.3f.npy" % (freq_MHz) 
      
      if not pt_source:
         gsm_cross_visibility_real_array_XX = np.load(gsm_cross_visibility_real_array_filename_XX)
         gsm_cross_visibility_imag_array_XX = np.load(gsm_cross_visibility_imag_array_filename_XX)   
         gsm_cross_visibility_real_array_YY = np.load(gsm_cross_visibility_real_array_filename_YY)
         gsm_cross_visibility_imag_array_YY = np.load(gsm_cross_visibility_imag_array_filename_YY)     
         gsm_cross_visibility_real_array_XY = np.load(gsm_cross_visibility_real_array_filename_XY)
         gsm_cross_visibility_imag_array_XY = np.load(gsm_cross_visibility_imag_array_filename_XY) 
         gsm_cross_visibility_real_array_YX = np.load(gsm_cross_visibility_real_array_filename_YX)
         gsm_cross_visibility_imag_array_YX = np.load(gsm_cross_visibility_imag_array_filename_YX) 
      else:
         zenith_point_source_cross_visibility_real_array_XX = np.load(zenith_point_source_cross_visibility_real_array_filename_XX)
         zenith_point_source_cross_visibility_imag_array_XX = np.load(zenith_point_source_cross_visibility_imag_array_filename_XX)              
         zenith_point_source_auto_array_XX = np.load(zenith_point_source_auto_array_filename_XX)
         zenith_point_source_cross_visibility_real_array_YY = np.load(zenith_point_source_cross_visibility_real_array_filename_YY)
         zenith_point_source_cross_visibility_imag_array_YY = np.load(zenith_point_source_cross_visibility_imag_array_filename_YY)              
         zenith_point_source_auto_array_YY = np.load(zenith_point_source_auto_array_filename_YY)
         zenith_point_source_cross_visibility_real_array_XY = np.load(zenith_point_source_cross_visibility_real_array_filename_XY)
         zenith_point_source_cross_visibility_imag_array_XY = np.load(zenith_point_source_cross_visibility_imag_array_filename_XY)              
         zenith_point_source_cross_visibility_real_array_YX = np.load(zenith_point_source_cross_visibility_real_array_filename_YX)
         zenith_point_source_cross_visibility_imag_array_YX = np.load(zenith_point_source_cross_visibility_imag_array_filename_YX)              
      
      baseline_number_array = np.load(baseline_number_array_filename)

      #miriad expects uvw in nanosecs so need to multiply by wavelength, divide by speed of light and times by 1e9
      #swapping uu and vv and making the new vv negative () fixes the rotation, but then the uvw ratios between data ad sime don't match so well (this way they match exactly)
      uu_array = np.load(uu_array_filename) * wavelength * 1.0e9 / c
      vv_array = np.load(vv_array_filename) * wavelength * 1.0e9 / c
      ww_array = np.load(ww_array_filename) * wavelength * 1.0e9 / c
      
      #print(uu_array_m)
      #print(vv_array_m)
      #print(ww_array_m)
      
      #phase centre zenith
      #print(np.max(ww_array_m))
      #print(gsm_cross_visibility_real_array[0:10])
      #print(gsm_cross_visibility_imag_array[0:10])
      
      ##GSM:
      if not pt_source:
         auto_array_X = np.load(gsm_auto_array_filename_X)
         cross_visibility_complex_array_XX = gsm_cross_visibility_real_array_XX + 1j*gsm_cross_visibility_imag_array_XX
         auto_array_Y = np.load(gsm_auto_array_filename_Y)
         cross_visibility_complex_array_YY = gsm_cross_visibility_real_array_YY + 1j*gsm_cross_visibility_imag_array_YY
         cross_visibility_complex_array_XY = gsm_cross_visibility_real_array_XY + 1j*gsm_cross_visibility_imag_array_XY
         cross_visibility_complex_array_YX = gsm_cross_visibility_real_array_YX + 1j*gsm_cross_visibility_imag_array_YX
      else:
         #Try point source instead
         auto_array_X = np.load(zenith_point_source_auto_array_filename_X)
         cross_visibility_complex_array_XX = zenith_point_source_cross_visibility_real_array_XX + 1j*zenith_point_source_cross_visibility_imag_array_XX
         auto_array_Y = np.load(zenith_point_source_auto_array_filename_Y)
         cross_visibility_complex_array_YY = zenith_point_source_cross_visibility_real_array_YY + 1j*zenith_point_source_cross_visibility_imag_array_YY
         cross_visibility_complex_array_XY = zenith_point_source_cross_visibility_real_array_XY + 1j*zenith_point_source_cross_visibility_imag_array_XY
         cross_visibility_complex_array_YX = zenith_point_source_cross_visibility_real_array_YX + 1j*zenith_point_source_cross_visibility_imag_array_YX
         
         
      cross_visibility_real_array_XX = np.real(cross_visibility_complex_array_XX)
      cross_visibility_imag_array_XX = np.imag(cross_visibility_complex_array_XX)
      cross_visibility_real_array_YY = np.real(cross_visibility_complex_array_YY)
      cross_visibility_imag_array_YY = np.imag(cross_visibility_complex_array_YY) 
      cross_visibility_real_array_XY = np.real(cross_visibility_complex_array_XY)
      cross_visibility_imag_array_XY = np.imag(cross_visibility_complex_array_XY)
      cross_visibility_real_array_YX = np.real(cross_visibility_complex_array_YX)
      cross_visibility_imag_array_YX = np.imag(cross_visibility_complex_array_YX)     
      

      #not needed for zenith?
      #cross_visibility_complex_array_add_phase = cross_visibility_complex_array * np.exp(1j*ww_array_m)
      #cross_visibility_real_array = np.real(cross_visibility_complex_array_add_phase)
      #cross_visibility_imag_array = np.imag(cross_visibility_complex_array_add_phase)
      
      #miriad wants u.v.w in units of seconds, but it divides by c itself! (does it also divide by wavelength?) ... nanosecs?
      #i am lost, but maybe there is hope here: https://github.com/HERA-Team/aipy/blob/main/doc/source/tutorial.rst

      #print(c)
      #print(uu_array)
      #print(vv_array)
      #print(ww_array)

      #from AIPY doco:write(preamble, data, flags=None)
      #Write the next data record. data must be a complex, masked array. preamble must be (uvw, t, (i,j)), where
      #uvw is an array of u,v,w, t is the Julian date, and (i,j) is an antenna pair.
      #add data
      
      #put in cross first then do autos
      #cant do x and y in a loop, can only go through the uvdata file once in order (miriad, am I right?)
      for baseline_number_index, baseline_number in enumerate(baseline_number_array):
         complex_cross_XX = np.asarray([cross_visibility_real_array_XX[baseline_number_index] + 1j*cross_visibility_imag_array_XX[baseline_number_index]])
         cross_vis_XX = np.ma.array(complex_cross_XX, mask=data_mask, dtype=np.complex64)
         complex_cross_YY = np.asarray([cross_visibility_real_array_YY[baseline_number_index] + 1j*cross_visibility_imag_array_YY[baseline_number_index]])
         cross_vis_YY = np.ma.array(complex_cross_YY, mask=data_mask, dtype=np.complex64)
         complex_cross_XY = np.asarray([cross_visibility_real_array_XY[baseline_number_index] + 1j*cross_visibility_imag_array_XY[baseline_number_index]])
         cross_vis_XY = np.ma.array(complex_cross_XY, mask=data_mask, dtype=np.complex64)
         complex_cross_YX = np.asarray([cross_visibility_real_array_YX[baseline_number_index] + 1j*cross_visibility_imag_array_YX[baseline_number_index]])
         cross_vis_YX = np.ma.array(complex_cross_YX, mask=data_mask, dtype=np.complex64)
         #print(cross_vis)
         uvw_array = [uu_array[baseline_number_index],vv_array[baseline_number_index],ww_array[baseline_number_index]]
         #print(uvw_array)
         ant1,ant2 = decode_baseline(baseline_number)
         #uvw = aa.gen_uvw(ant1,ant2)
         #print(uvw.shape)
         #print(uvw)
         #uvw_array = [uvw[0,0,0],uvw[1,0,0],uvw[2,0,0]]
         #print(baseline_number)
         #print(ant1,ant2)
         #LOOK AT AIPY DOCS, MIRIAD USES A DIFFERENT BASELINE INDEXING SCHEME! (no zero index antenna) this could be the problem function: ij2bl
         #I think I can just add one to the antenna indices ... nope
         uvw_12 = np.array(uvw_array, dtype=np.double)
         preamble = (uvw_12, eda2_astropy_time.jd, (ant1,ant2)) 
         #real data from 2018 paper is missing antenna 116 (among others) but can't have ant 255 cause it breaks casa
         if ant2<255 and ant1<255:
            #print("changing pol to -5 xx")
            uv['pol'] = -5   #-5 is xx, -6 yy, -7 xy 8 yx
            #put in the auto:
            if ant2==(ant1+1):
               auto = auto_array_X[ant1]
               auto_power = np.asarray([auto])
               auto_vis = np.ma.array(auto_power, mask=data_mask, dtype=np.complex64)
               uvw_array = [0,0,0]
               uvw_11 = np.array(uvw_array, dtype=np.double)
               auto_preamble = (uvw_11, eda2_astropy_time.jd, (ant1,ant1)) 
               uv.write(auto_preamble,auto_vis)
            uv.write(preamble,cross_vis_XX)
            #print("changing pol to -6 yy")
            uv['pol'] = -6   #-5 is xx, -6 yy, -7 xy 8 yx
            #put in the auto:
            if ant2==(ant1+1):
               auto = auto_array_Y[ant1]
               auto_power = np.asarray([auto])
               auto_vis = np.ma.array(auto_power, mask=data_mask, dtype=np.complex64)
               uvw_array = [0,0,0]
               uvw_11 = np.array(uvw_array, dtype=np.double)
               auto_preamble = (uvw_11, eda2_astropy_time.jd, (ant1,ant1)) 
               uv.write(auto_preamble,auto_vis)
            uv.write(preamble,cross_vis_YY)
            #print("changing pol to -7 xy")
            uv['pol'] = -7
            #put in the auto:
            if ant2==(ant1+1):
               auto = auto_array_X[ant1] * 0.
               auto_power = np.asarray([auto])
               auto_vis = np.ma.array(auto_power, mask=data_mask, dtype=np.complex64)
               uvw_array = [0,0,0]
               uvw_11 = np.array(uvw_array, dtype=np.double)
               auto_preamble = (uvw_11, eda2_astropy_time.jd, (ant1,ant1)) 
               uv.write(auto_preamble,auto_vis)
            uv.write(preamble,cross_vis_XY)
            #print("changing pol to -8 yx")
            uv['pol'] = -8
            #put in the auto:
            if ant2==(ant1+1):
               auto = auto_array_X[ant1] * 0.
               auto_power = np.asarray([auto])
               auto_vis = np.ma.array(auto_power, mask=data_mask, dtype=np.complex64)
               uvw_array = [0,0,0]
               uvw_11 = np.array(uvw_array, dtype=np.double)
               auto_preamble = (uvw_11, eda2_astropy_time.jd, (ant1,ant1)) 
               uv.write(auto_preamble,auto_vis)
            uv.write(preamble,cross_vis_YX)

      #for auto_index,auto in enumerate(auto_array):
      #   auto_power = np.asarray([auto])
      #   auto_vis = np.ma.array(auto_power, mask=data_mask, dtype=np.complex64)
      #   uvw_array = [0,0,0]
      #   uvw_11 = np.array(uvw_array, dtype=np.double)
      #   preamble = (uvw_11, eda2_astropy_time.jd, (auto_index,auto_index)) 
      #   uv.write(preamble,auto_vis)
      
      del(uv)
      
      for image_pol in ['xx','yy']:
      
         ##check by image
         map_name = 'test_eda_%0.3f_%s.image' % (freq_MHz,image_pol)
         beam_name = 'test_eda_%0.3f_%s.beam' % (freq_MHz,image_pol)
         map_name_clean = 'test_eda_%0.3f_%s_clean.image' % (freq_MHz,image_pol)
         map_name_restor = 'test_eda_%0.3f_%s_restor.image' % (freq_MHz,image_pol)
         map_name_fits = 'test_eda_%0.3f_%s.fits' % (freq_MHz,image_pol)
         cmd = "rm -rf %s %s %s %s %s" % (map_name,beam_name,map_name_fits,map_name_clean,map_name_restor)
         print(cmd)
         os.system(cmd)
         cmd = "invert vis=%s map=%s beam=%s imsize=512 cell=900 stokes=%s robust=0" % (mir_file,map_name,beam_name,image_pol)
         print(cmd)
         os.system(cmd) 
         cmd = "clean map=%s beam=%s niters=50 imsize=512 cell=900 stokes=%s out=%s" % (map_name,beam_name,image_pol,map_name_clean)
         print(cmd)
         os.system(cmd) 
         cmd = "restor map=%s beam=%s model=%s out=%s" % (map_name,beam_name,map_name_clean,map_name_restor)
         print(cmd)
         os.system(cmd)      
         cmd = "fits op=xyout in=%s out=%s" % (map_name_restor,map_name_fits)
         print(cmd)
         os.system(cmd)


      check_uv = a.miriad.UV(mir_file)
      #print(check_uv.items())
      #print(check_uv.vars())
      #print(check_uv['nchan'])
      #print(check_uv['antpos'])
      print(check_uv['pol'], a.miriad.pol2str[check_uv['pol']])

      
      ##export uvfits
      cmd = "rm -rf %s" % (mir_file_uvfits_name)
      print(cmd)
      os.system(cmd)
      cmd = "fits op=uvout in=%s out=%s" % (mir_file,mir_file_uvfits_name)
      print(cmd)
      os.system(cmd)
      
      #look at the uvfits file
      uvfits_filename = mir_file_uvfits_name
      print("%s" % uvfits_filename)
      hdulist = fits.open(uvfits_filename)
      #hdulist.info()
      info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print(info_string)
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      hdulist.close()
      sim_UU_s_array = uvtable['UU']
      sim_VV_s_array = uvtable['VV']
      sim_WW_s_array = uvtable['WW']
      sim_baselines = uvtable['baseline']
      sim_visibilities = uvtable['DATA']
      sim_visibilities_shape = sim_visibilities.shape
      #print('sim_UU_s_array')
      #print(sim_UU_s_array)
      #sim_UU_m_array = sim_UU_s_array * c
      #print('sim_UU_m_array')
      #print(sim_UU_m_array)
      #u_in_miriad = sim_UU_s_array[0]
      #ratio_u_in = u_in/u_in_miriad
      #print('ratio_u_in')
      #print(ratio_u_in)
      print("sim visibilities_shape")
      print(sim_visibilities_shape)
      print('sim_baselines')
      print(len(sim_baselines))
      print(np.max(sim_baselines))
      
      
      #what is going on with these baselines
      for baseline_num_index,baseline_num in enumerate(sim_baselines):
         #print(sim_baselines[index])
         if baseline_num > (2**16):
            print(baseline_num)
            ant1,ant2 = decode_baseline(baseline_num)
            print(ant1,ant2)
      
      #got it to work by not adding the last antenna i.e above:
      #if ant2<255 and ant1<255:
            #uv.write(preamble,cross_vis)
      #need a better solution! but will be interesting to see if ms calibrates/images or if proper antenna coords are needed (I think yes for calibrate .. but wsclean might be okay ....)
      
      
      
      #print(len(sim_baseline))
      #print(len(sim_UU_s_array))
      #convert from miriad baseline numbers
      #converted_sim_baseline_ant_nums = [aa.bl2ij(bl) for bl in sim_baselines]
      converted_sim_baseline_ant_nums = [decode_baseline(bl) for bl in sim_baselines]
      converted_sim_baselines = [(cal_standard_baseline_number(ant1,ant2)) for ant1,ant2 in converted_sim_baseline_ant_nums]
      converted_sim_baselines = np.asarray(converted_sim_baselines,dtype='f')
      #print(converted_sim_baselines)

      #try to read the data in to casa as a ms
      #verify with wsclean - too many antennas - 255 max?
      #read in the uvfits file
      casa_cmd_filename = 'import_uvfits.sh'
      cmd = "rm -rf %s" % (mir_file_ms_name)
      print(cmd)
      os.system(cmd)
            
      cmd = "importuvfits(fitsfile='%s',vis='%s')" % (mir_file_uvfits_name,mir_file_ms_name)
      print(cmd)
      os.system(cmd)
      
      with open(casa_cmd_filename,'w') as f:
         f.write(cmd)
           
      cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
      print(cmd)
      os.system(cmd)
      
      test_image_name = mir_file_ms_name.split('.ms')[0]
      #try imaging the ms:
      #if pol=='X':
      #   wsclean_imsize = '512'
      #   wsclean_scale = '900asec'
      #   cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx -data-column DATA  %s " % (test_image_name,wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
      #   print(cmd)
      #   os.system(cmd) 
      
      #try adding a model column
      #use kariukis ms_utils
      ms_table = table(mir_file_ms_name,readonly=False)
      model_data = get_data(ms_table, col="DATA")
      add_col(ms_table, "MODEL_DATA")
      put_col(ms_table, "MODEL_DATA", model_data)
      ms_table.close()

      #try imaging the model column of the ms:
      wsclean_imsize = '512'
      wsclean_scale = '900asec'
      cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column MODEL_DATA  %s " % (test_image_name+"_model",wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
      print(cmd)
      os.system(cmd)
      
      #try calibrate?
      gain_solutions_name = 'test_cal_%s_%s_calibrate_sols.bin' % (EDA2_chan,EDA2_obs_time)
      calibrate_options = ''
      cmd = "rm -rf %s" % (gain_solutions_name)
      print(cmd)
      os.system(cmd)  
      #calibrate
      cmd = "calibrate %s %s %s " % (calibrate_options,mir_file_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      #plot cal sols
      cmd = "aocal_plot.py %s  " % (gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "applysolutions %s %s  " % (mir_file_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      #test image the CORRECTED data'
      cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (test_image_name+"_corrected",wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
      print(cmd)
      os.system(cmd)
      
      #check cal
      #calibrate_with_complex_beam_model(model_ms_name,eda2_ms_name)
      
      
      #if pol=='X':
      #   wsclean_imsize = '512'
      #   wsclean_scale = '900asec'
      #   cmd = "wsclean -name test_eda_%0.3f_%s_wsclean -size %s %s -multiscale -niter 1000 -scale %s -pol xx  %s " % (freq_MHz,pol,wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
      #   print(cmd)
      #   os.system(cmd) 
      #
      #sys.exit()   

      #Lets look at some data at the same LST
      uvfits_filename = "/md0/EoR/EDA2/20200303_data/%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      ms_name = "/md0/EoR/EDA2/20200303_data/%s/av_chan_%s_%s_plus_%s_obs.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      print("%s" % uvfits_filename)
      hdulist = fits.open(uvfits_filename)
      #hdulist.info()
      info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print(info_string)
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      hdulist.close()
      data_UU_s_array = uvtable['UU']
      data_VV_s_array = uvtable['VV']
      data_WW_s_array = uvtable['WW']
      data_baselines = uvtable['baseline']
      data_visibilities = uvtable['DATA']
      data_visibilities_shape = data_visibilities.shape
      print("data visibilities_shape")
      print(data_visibilities_shape)
      #print(len(data_baseline))
      print('data_UU_s_array')
      print(data_UU_s_array)
      data_UU_m_array = data_UU_s_array * c
      print('data_UU_m_array')
      print(data_UU_m_array)
      #print(np.max(data_baselines))
      #print(np.min(data_baselines))

      #make wsclean image to compare (no need to do twice, just image xx,yy):
      #wsclean does not resore the cleaned image - clean is same as dirty!
      #if pol=='X':
      #   wsclean_imsize = '512'
      #   wsclean_scale = '900asec'
      #   cmd = "wsclean -name cal_chan_%s_%s_ms -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
      #   print(cmd)
      #   os.system(cmd) 
      
      #uvfits_filename = "/md0/EoR/EDA2/20200303_data/64/cal_chan_64_20200303T133741.uvfits" 
      #uvfits_vis_name = 'test_eda_actual.vis'
      #cmd = "rm -rf %s" % (uvfits_vis_name)
      #print(cmd)
      #os.system(cmd)
      #cmd = "fits op=uvin in=%s out=%s" % (uvfits_filename,uvfits_vis_name)
      #print(cmd)
      #os.system(cmd)
      ###check by image
      #map_name = 'test_eda_actual.image'
      #beam_name = 'test_eda_actual.beam'
      #cmd = "rm -rf %s %s" % (map_name,beam_name)
      #print(cmd)
      #os.system(cmd)
      #cmd = "invert vis=%s map=%s beam=%s imsize=512 cell=600 stokes=xx" % (uvfits_vis_name,map_name,beam_name)
      #print(cmd)
      #os.system(cmd)


      #make a sim uv file that only has baselines that match the data file. Usually data will have some flagged
      #antennas, so this alleviates the problem of the 255 limit for casa/ms
      #Cant work out this baseline number stuff, just look at u,v values
      #find common indices 

#took this from: https://stackoverflow.com/questions/16216078/test-for-membership-in-a-2d-numpy-array
def asvoid(arr):
    """
    Based on http://stackoverflow.com/a/16973510/190597 (Jaime, 2013-06)
    View the array as dtype np.void (bytes). The items along the last axis are
    viewed as one value. This allows comparisons to be performed on the entire row.
    """
    arr = np.ascontiguousarray(arr)
    if np.issubdtype(arr.dtype, np.floating):
        """ Care needs to be taken here since
        np.array([-0.]).view(np.void) != np.array([0.]).view(np.void)
        Adding 0. converts -0. to 0.
        """
        arr += 0.
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))


def inNd(a, b, assume_unique=False):
    a = asvoid(a)
    b = asvoid(b)
    return np.in1d(a, b, assume_unique)
         

def calibrate_with_complex_beam_model(model_ms_name,eda2_ms_name):
      #casacore stuff: https://casacore.github.io/python-casacore/casacore_tables.html#casacore.tables.tablesummary
      #Okay looks good, now lets try to calibrate the corresponding ms from the EDA2
      #eda2_ms_name = "/md0/EoR/EDA2/20200303_data/64/20200303_133733_eda2_ch32_ant256_midday_avg8140.ms" 
      
      ##This whole plan is bad:
      ##need to get rid of ant 255 from eda2 data (to match model)
      ##export as uvfits with antenna range
      ##write out the uvfits file
      #casa_cmd_filename = 'export_cropped_uvfits.sh'
      #cropped_uvits_filename = eda2_ms_name.split('/')[-1].split('.ms')[0] + "_cropped255.uvfits"
      #cropped_ms_filename = eda2_ms_name.split('/')[-1].split('.ms')[0] + "_cropped255.ms"
      #cmd = "rm -rf %s %s %s" % (cropped_uvits_filename,cropped_ms_filename,casa_cmd_filename)
      #print(cmd)
      #os.system(cmd)
      #
      #ints = np.arange(0,256)
      #string_ints = [str(int) for int in ints]
      #antenna_string = ",".join(string_ints)
      #
      #cmd = "exportuvfits(vis='%s',fitsfile='%s',datacolumn='data',overwrite=True,writestation=False)" % (eda2_ms_name,cropped_uvits_filename)
      #print(cmd)
      #os.system(cmd)
      #
      #with open(casa_cmd_filename,'w') as f:
      #   f.write(cmd)
      #     
      #cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
      #print(cmd)
      #os.system(cmd)
      #
      #now import the cropped one!

      eda2_ms_table = table(eda2_ms_name,readonly=False)
      eda2_table_summary = tablesummary(eda2_ms_name)
      #print(eda2_table_summary)
      eda2_data = get_data(eda2_ms_table)
      print(eda2_data.shape)
      eda2_uvw = get_uvw(eda2_ms_table)
      eda2_ant1, eda2_ant2 = get_ant12(eda2_ms_name)
      eda2_ants = np.vstack((eda2_ant1,eda2_ant2)).T
      #print(eda2_ants.shape)
      #print(eda2_ants[0:300,0])
      #print(eda2_ants[0:300,1])
      #eda2 data is missing some ants, but has autos
      #eda2_flags = get_flags(eda2_ms_table)
      #print(eda2_flags.shape)
      #print(eda2_flags)
      
      model_ms_table = table(model_ms_name,readonly=True)
      model_data = get_data(model_ms_table)
      model_uvw = get_uvw(model_ms_table)    
      model_ant1, model_ant2 = get_ant12(model_ms_name)
      model_ants = np.vstack((model_ant1,model_ant2)).T
      #print(model_ants.shape)      
      #print(model_ants[0:300,0])
      #print(model_ants[0:300,1])
      #model has no autos and has all ants except 255

      #do it by antenna number instead
      
      #find common indices: https://stackoverflow.com/questions/2333593/return-common-element-indices-between-two-numpy-arrays
      #This won't work because uvw is multi dimensional
      #model_ms_indices = np.nonzero(np.in1d(model_uvw, eda2_uvw))[0]
      model_ms_indices = inNd(model_ants, eda2_ants, assume_unique=False)
      n_common = np.count_nonzero(model_ms_indices) 
      print(n_common)
      eda2_ms_indices = inNd(eda2_ants, model_ants, assume_unique=False)
      n_common = np.count_nonzero(eda2_ms_indices)
      print(n_common)
      
      #ind = np.lexsort((b,a)) # Sort by a, then by b
      #common_eda2_uvw_sorted = np.sort(common_eda2_uvw,axis=0)
      
      eda2_common_ant1 = eda2_ant1[eda2_ms_indices]
      eda2_common_ant2 = eda2_ant2[eda2_ms_indices]
      eda2_common_sort_inds = np.lexsort((eda2_common_ant2,eda2_common_ant1)) # Sort by a, then by b

      model_common_ant1 = model_ant1[model_ms_indices]
      model_common_ant2 = model_ant2[model_ms_indices]
      model_common_sort_inds = np.lexsort((model_common_ant2,model_common_ant1)) # Sort by a, then by b
      
      #don't sort the ant arrays, these are what we use to sort the other arrays!
      #model_common_ant1_sorted = model_common_ant1[model_common_sort_inds]
      #model_common_ant2_sorted = model_common_ant2[model_common_sort_inds]
            
      common_eda2_uvw = eda2_uvw[eda2_ms_indices]
      common_eda2_uvw_sorted = common_eda2_uvw[eda2_common_sort_inds]
      
      common_eda2_data = eda2_data[eda2_ms_indices]
      common_eda2_data_sorted = common_eda2_data[eda2_common_sort_inds]
      #print(common_eda2_data_sorted.shape)
      #print(common_eda2_data_sorted[0:10,0,0])
      
      common_model_uvw = model_uvw[model_ms_indices]
      common_model_uvw_sorted = common_model_uvw[model_common_sort_inds]

      common_model_data = model_data[model_ms_indices]
      common_model_data_sorted = common_model_data[model_common_sort_inds]
      #print(common_model_data_sorted.shape)
      #print(common_model_data_sorted[0:10,0,0])      
      
      #for model go through ant2 array and data array and uvw array, insert zeros after wherever ant==254
      ant2_list=[]
      for ant2_index,ant2 in enumerate(model_common_ant2):
         if ant2==254:
            ant2_list.append(ant2_index)
      #print(len(ant2_list))
      #print(ant2_list)
      
      old_model_common_ant2 = model_common_ant2
      old_model_common_ant1 = model_common_ant1
      old_common_model_data_sorted = common_model_data_sorted
      counter=0
      
      #a = np.array([[1, 1], [2, 2], [3, 3]])
      #print(a)
      #b = np.insert(a, 1, 5, axis=0)
      #print(b)
      
      
      for ant2_index in ant2_list:
         ant2_index+=counter
         ant1_value = old_model_common_ant1[ant2_index-1]
         #new_data_value = np.zeros((1,4))
         #print(new_data_value.shape)
         #sys.exit()
         #print(old_array.shape)
         new_model_common_ant2 = np.insert(old_model_common_ant2,ant2_index+1,255)
         new_model_common_ant1 = np.insert(old_model_common_ant1,ant2_index+1,ant1_value)
         new_common_model_data_sorted = np.insert(old_common_model_data_sorted,ant2_index+1,0,axis=0)
         #print(old_common_model_data_sorted[ant2_index-2:ant2_index+5])
         #print(new_common_model_data_sorted[ant2_index-2:ant2_index+6])
         #print(new_array[ant2_index+1])
         old_model_common_ant2 = new_model_common_ant2
         old_model_common_ant1 = new_model_common_ant1
         old_common_model_data_sorted = new_common_model_data_sorted
         counter+=1
      new_model_common_ant2 = np.append(new_model_common_ant2,np.array([254,255,255]))
      new_model_common_ant1 = np.append(new_model_common_ant1,np.array([254,254,255]))
      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      
      #print(new_common_model_data_sorted.shape)
      repetitions = 32
      new_common_model_data_sorted_tile = np.tile(new_common_model_data_sorted, (repetitions, 1))
      #print(new_common_model_data_sorted_tile.shape)
      
      #print(len(new_model_common_ant2_sorted))
      #print((new_model_common_ant2_sorted[-200:-1]))
      #print(eda2_ant2.shape)
      #print(eda2_ant2[-200:-1])
      
      #print(common_model_data_sorted.shape)
      #print(common_model_data_sorted[-5:])
      #print(new_common_model_data_sorted.shape)
      #print((new_common_model_data_sorted[-5:]))

      #model has no ant 255 and no autos, data is missing some other antennas, but has 255 and autos
      #going to need to add autos into model, and to add in missing ant 255 correlations as dummy data, and then flag that ant in the data before calibration
      #arrrgghhhh
      try:
         add_col(eda2_ms_table, "MODEL_DATA")
      except:
         pass
      put_col(eda2_ms_table, "MODEL_DATA", new_common_model_data_sorted_tile)
      eda2_ms_table.close()     


      #try imaging the model column of the ms:
      wsclean_imsize = '512'
      wsclean_scale = '900asec'
      test_image_name = "complex_beam_test"
      cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column MODEL_DATA  %s " % (test_image_name+"_model",wsclean_imsize,wsclean_imsize,wsclean_scale,eda2_ms_name)
      print(cmd)
      os.system(cmd)
      
      #need to flag ant 255 .... also flux scale is about half (not treating polarisation correctly?)
      #use uv cutoff at half wavelength?
      #try cal with -ch 32
      
      #try calibrate?
      gain_solutions_name = 'complex_beam_test_calibrate_sols.bin' 
      calibrate_options = ''
      cmd = "rm -rf %s" % (gain_solutions_name)
      print(cmd)
      os.system(cmd)  
      #calibrate
      cmd = "calibrate -ch 32 %s %s %s " % (calibrate_options,eda2_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      #plot cal sols
      cmd = "aocal_plot.py %s  " % (gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "applysolutions %s %s  " % (eda2_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      #test image the CORRECTED data'
      cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (test_image_name+"_corrected",wsclean_imsize,wsclean_imsize,wsclean_scale,eda2_ms_name)
      print(cmd)
      os.system(cmd)



#times
EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan(EDA2_data_dir,EDA2_chan_list)
EDA2_obs_time_list_each_chan = EDA2_obs_time_list_each_chan[0:]
n_obs_concat_list = [len(obs_list) for obs_list in EDA2_obs_time_list_each_chan] 
EDA2_obs_time_list = [item[0] for item in EDA2_obs_time_list_each_chan] 

#chans
EDA2_chan_list_array = np.asarray(EDA2_chan_list)
freq_MHz_array = 400./512.*EDA2_chan_list_array
freq_MHz_list = freq_MHz_array[0:]
EDA2_chan_list = EDA2_chan_list[0:]

#LSTs
lst_hrs_list = []
for EDA2_obs_time_index,EDA2_obs_time in enumerate(EDA2_obs_time_list):
   #there might have been no obs:
   if EDA2_obs_time!=0:
      print(EDA2_obs_time)
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time)
      lst_hrs_list.append(lst_eda2_hrs)
   else:
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time_list[0])
      lst_hrs_list.append(lst_eda2_hrs)
      
      
      
      
      
      
      
      