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
import ephem
from scipy.io import loadmat
import aipy as a  #for making miriad files
from scipy.interpolate import interp2d,griddata,interp1d
from astropy.io import fits 
from casacore.tables import table, tablesummary
sys.path.append("/md0/code/git/ben-astronomy/ms")
from ms_utils import *
import statsmodels.api as sm
import numpy.polynomial.polynomial as poly

k = 1.38065e-23
c = 299792458.
fine_chan_width_Hz =  28935 #Hz should be:(400./512.)*1.0e6 

mwa_latitude_ephem = '-26.7'
mwa_longitude_ephem = '116.67'
mwa_latitude_deg = -26.70331940
mwa_longitude_deg = 116.670575
mwa_latitude_rad = float(mwa_latitude_deg/180.*np.pi)
mwa_longitude_rad = float(mwa_longitude_deg/180.*np.pi)
mwa_latitude_astropy = '-26.7d'
mwa_longitude_astropy = '116.67d'

EDA2_data_dir = '/md0/EoR/EDA2/20200303_data/'   #2020 paper
EDA2_chan_list = range(64,127)  #20200303:

def make_EDA2_obs_time_list_each_chan(base_dir,eda2_chan_list):
   obs_time_list_each_chan = []
   temp_txt_filename = 'uvfits_time_list.txt'
   for eda2_chan in eda2_chan_list:
      chan_obs_time_list = []
      chan_dir = "%s%s/" % (base_dir,eda2_chan)
      cmd = "ls -la %schan_%s_*.uvfits  > %s" % (chan_dir,eda2_chan,temp_txt_filename)
      os.system(cmd)
      with open(temp_txt_filename) as f:
         lines=f.readlines()
      for line in lines[1:-1]:
         obs_time = line.split('.uvfits')[0].split()[-1].split('_')[-1]
         chan_obs_time_list.append(obs_time)
      #there might be no obs, if so just put in a zero
      if len(chan_obs_time_list)!=0:
         obs_time_list_each_chan.append(chan_obs_time_list)
      else:
         obs_time_list_each_chan.append([0])
   return obs_time_list_each_chan
   
def get_eda2_lst(eda_time_string):
   year, month, day, hour, minute, second = eda_time_string[0:4], eda_time_string[4:6],eda_time_string[6:8], eda_time_string[9:11],eda_time_string[11:13],eda_time_string[13:15]
   eda2_observer = ephem.Observer()
   eda2_observer.lon, eda2_observer.lat = mwa_longitude_ephem, mwa_latitude_ephem
   eda2_observer.date = '%s/%s/%s %s:%s:%s' % (year,month,day,hour,minute,second)
   eda2_obs_lst = (eda2_observer.sidereal_time()) 
   #print("LST is")
   #print(eda2_obs_lst)
   eda2_obs_lst_hrs = eda2_obs_lst / 2 / np.pi * 24.
   
   return(eda2_obs_lst_hrs)

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

def calc_x_y_z_pos_from_E_N_U(E_pos_ant,N_pos_ant,U_pos_ant):
   #Look at instruction to go from E,N,U to u,v,w (https://web.njit.edu/~gary/728/Lecture6.html
   x_pos_ant = -1*N_pos_ant*np.sin(mwa_latitude_rad) + U_pos_ant*np.cos(mwa_latitude_rad)
   y_pos_ant = E_pos_ant
   z_pos_ant = N_pos_ant*np.cos(mwa_latitude_rad) + U_pos_ant*np.sin(mwa_latitude_rad)
                 
   return(x_pos_ant,y_pos_ant,z_pos_ant)

def decode_baseline(baseline_code):
   blcode = int(baseline_code)
   if baseline_code > 65536:
       baseline_code -= 65536
       a2 = int(baseline_code % 2048)
       a1 = int((baseline_code - a2) / 2048)
   else:
       a2 = int(baseline_code % 256)
       a1 = int((baseline_code - a2) / 256)
   return a1,a2
   
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
   
def cal_standard_baseline_number(ant1,ant2):
   baseline_number = (ant1 * 256) + ant2
   return baseline_number

def calc_LNA_impedance_eda2(freq_MHz):
   #See Daniel email from 18/09/2018 on normalisation factor, updated email 25 Nov for EDA2 normalisation calculation
   #For EDA2 the LNA impedance is actually based on RLC lumped circuit model:
   #Series component
   Ls = 2.e-9 #nH
   #Parallel component
   Rp = 914. #ohms
   Lp = 450.e-9 #nH
   Cp = 3.2e-12 #pF
   #So the full impedance is calculated as follows:
   w = 2.*np.pi*freq_MHz*1e6; # angular frequency

   Z = 1j*w*Ls + (1./Rp + (1j*w*Cp) + 1./(1j*w*Lp))**(-1)
   #You can double check your impedance calculation with Marcin’s paper (pg 4)
   return(Z)

def calc_beam_normalisation(freq_MHz,LNA_impedance):
   mu_0 = 4.*np.pi*10**(-7)
   w = 2.*np.pi*freq_MHz*1e6; # angular frequency
   normalisation_factor = (-4.*np.pi*1j / (mu_0 * w)) * LNA_impedance
   return(normalisation_factor)
   
def simulate_eda2_with_complex_beams(freq_MHz_list,lst_hrs,nside=512,antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt',plot_from_saved=False,EDA2_obs_time_list=[],sim_unity=True,sim_pt_source=False,check_figs=False):
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

   #don't convert to Jy, leave as K
   #Jy per pix conversion
   #healpix_pixel_area_sr = 4*np.pi/npix
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
      
      lst_hrs = float(lst_hrs_list[freq_MHz_index])
      lst_deg = lst_hrs * 15.

      #hpx rotate stuff, rotate the complex beams to zenith at the required LST
      dec_rotate = 90. - float(mwa_latitude_ephem)
      ra_rotate = lst_deg
      r_beam_dec = hp.Rotator(rot=[0,dec_rotate], coord=['C', 'C'], deg=True) 
      r_beam_ra = hp.Rotator(rot=[ra_rotate,0], coord=['C', 'C'], deg=True)    
      #see Jishnu's SITARA paper 1, appendix B
      #first get it working for a single baseline, ant index 0 and 1
      #Get the right beams (EEPs from Daniel) and convert (interpolate) to healpix
   
      EDA2_obs_time = EDA2_obs_time_list[freq_MHz_index]
   
      #don't convert to Jy, leave as K
      #need to either generate gsm map in MJy.sr and convert to Jy/pix, or convert manually:
      #scale = (2. * k * 1.0e26 * healpix_pixel_area_sr) / (wavelength**2)
      #print("scale map by %s to get to Jy/pix" % scale)
      
      unity_sky_value = 1.
      
      #Do all the gsm sky stuff outside the pol loop
      #Now we have beams, need a sky!
      gsm = GlobalSkyModel(freq_unit='MHz')
      gsm_map_512 = gsm.generate(freq_MHz)
      full_nside = hp.npix2nside(gsm_map_512.shape[0])
      #print(full_nside)
      #point_source_at_zenith_sky_512 = gsm_map_512 * 0.
      
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
            fig_name="check1_pt_src_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
      
         #i've got my coord system messed up here a bit - dec (theta) should go first in the rotation ...
         #https://zonca.dev/2021/03/rotate-maps-healpy.html

         r_pt_src_ra_dec = hp.Rotator(coord=['C'],rot=[ra_rotate_gsm,dec_rotate_gsm], deg=True)
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
         
         if check_figs:
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
      
         if check_figs:
            plt.clf()
            map_title=""
            ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
            hp.mollview(map=point_source_at_zenith_sky,rot=(0,90,0),title=map_title)
            fig_name="check4ra_dec_pt_src_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
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
   
      if check_figs:
         plt.clf()
         map_title=""
         hp.orthview(map=rotated_gsm_C_ra_dec_512,half_sky=False,rot=(0,90,0),title=map_title)
         #hp.mollview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512,nside),rot=(0,90,0),title=map_title)
         fig_name="check4_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
      
      rotated_gsm_C_ra_dec_512_extra_90 = r_beam.rotate_map_alms(rotated_gsm_C_ra_dec_512,use_pixel_weights=False)
      
      if check_figs:
         plt.clf()
         map_title=""
         hp.orthview(map=hp.ud_grade(rotated_gsm_C_ra_dec_512_extra_90,nside),half_sky=False,rot=(0,90,0),title=map_title)
         ##hp.mollview(map=rotated_gsm_C_ra_dec_512_extra_90,rot=(0,90,0),title=map_title)
         fig_name="check4a_extra_gsm_rot_ra_decLST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
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
         fig_name="check_gsm_dgrade_LST_%0.1f_%0.3f_MHz.png" % (lst_deg,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
      
      
      #make a cube with copies of tsky, dimensions (npix,test_n_ants)
      gsm_repeats_array = np.tile(gsm_map, (test_n_ants,1))
      gsm_repeats_array = np.transpose(gsm_repeats_array)
      
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
               
               
               #Daniels version of normalising the beam i.e. getting in correct unitless form  by getting rid on m squared
               #Doesn't have to be 1 at zenith, just need to make sure use same normalised beam for signal extraction
               
               lna_impedance = calc_LNA_impedance_eda2(freq_MHz)
               beam_norm = calc_beam_normalisation(freq_MHz,lna_impedance)
               print("normalising beam by %0.5f %0.5fj " % (beam_norm.real,beam_norm.imag))
               
               #find where the total power (Ephi*Ephi + E_theta*Etheta) is maximum and divide both E_phi and E_theta by their
               #values at that position (squared). This gives you a power pattern with max 1 for each beam map and retains the
               #relative vlaues between E_phi and E_theta and the ss,yy,xy,yx beams...
               
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
                  
                  power_pattern_cube_mag = np.abs(power_pattern_cube)
                  #what if we don't normalise for the beam weights? - don't if input hpx map is in Jy/pix
                  #power_pattern_cube_mag_sum_array = np.einsum('ij->j',power_pattern_cube_mag)
                  
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
                        fig_name="check6_pt_source_sky_beam_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz)
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
                  #fig_name="check4_pattern_%s_%s_%s%s_%0.3f_MHz.png" % (ant_index_1,'0',pol1,pol2,freq_MHz)
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
                                 
                  if check_figs:
                     #####sanity check:
                     plt.clf()
                     map_title=""
                     hp.orthview(map=gsm_sky_beam_cube[:,0],half_sky=False,rot=(0,90,0),title=map_title)
                     #hp.mollview(map=gsm_sky_beam_cube[:,0],rot=(0,90,0),title=map_title)
                     fig_name="check5_gsm_sky_beam_%0.1f_%s_%s_%s%s_%0.3f_MHz.png" % (lst_deg,ant_index_1,'0',pol1,pol2,freq_MHz)
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
                

                  if (pol_index1==0 and pol_index2==0 and freq_MHz_index==0):
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
                    
      
               if sim_pt_source:
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
            
               ##plot the real part of the visibilty 
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

def write_to_miriad_vis(freq_MHz_list,lst_hrs_list,EDA2_obs_time_list,antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt',input_sky="gsm",check_figs=False):
   with open(antenna_layout_filename) as f:
      lines = f.readlines()
      n_ants = len(lines) 
   n_baselines = n_ants*(n_ants-1) / 2.
   print("read %s" % antenna_layout_filename)
   print("n_baselines from %s ants: %s" % (n_ants, n_baselines))
   
   for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
      #time stuff
      EDA2_obs_time = EDA2_obs_time_list[freq_MHz_index]
      year, month, day, hour, minute, second = EDA2_obs_time[0:4], EDA2_obs_time[4:6],EDA2_obs_time[6:8], EDA2_obs_time[9:11],EDA2_obs_time[11:13],EDA2_obs_time[13:15]
      eda2_astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (float(year), float(month), float(day), float(hour), float(minute), float(second))
      print(eda2_astropy_time_string)
      eda2_astropy_time = Time(eda2_astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
  
      lst_hrs = float(lst_hrs_list[freq_MHz_index])
      lst_deg = lst_hrs * 15.
      
      wavelength = 300./freq_MHz
      #initiate the miriad uv file outside the pol loop
      print("initiate the miriad uv file outside the pol loop")
      print("writing to miriad file")
      NFFT = 1
      SITE_LON = 116.670575
      FREQ_GHz = np.array([freq_MHz/1000.])#np.linspace(0.000, 0.250, NFFT);
      
      #timestr = time.strftime("%Y%m%d-%H%M%S")
      if input_sky=="gsm":
         mir_file = "%s_%0.3f.vis" % (EDA2_obs_time,freq_MHz)
         mir_file_uvfits_name = "%s_%0.3f.uvfits" % (EDA2_obs_time,freq_MHz)
         mir_file_ms_name = "%s_%0.3f.ms" % (EDA2_obs_time,freq_MHz)
      else:
         mir_file = "%s_%0.3f.vis" % (input_sky,freq_MHz)
         mir_file_uvfits_name = "%s_%0.3f.uvfits" % (input_sky,freq_MHz)
         mir_file_ms_name = "%s_%0.3f.ms" % (input_sky,freq_MHz)     
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
      
      if input_sky=="unity":
         cross_visibility_real_array_filename_XX = "%s_cross_visibility_real_array_XX_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_imag_array_filename_XX = "%s_cross_visibility_imag_array_XX_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_real_array_filename_YY = "%s_cross_visibility_real_array_YY_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_imag_array_filename_YY = "%s_cross_visibility_imag_array_YY_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_real_array_filename_XY=  "%s_cross_visibility_real_array_XY_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_imag_array_filename_XY = "%s_cross_visibility_imag_array_XY_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_real_array_filename_YX = "%s_cross_visibility_real_array_YX_%0.3f.npy" % (input_sky,freq_MHz)
         cross_visibility_imag_array_filename_YX = "%s_cross_visibility_imag_array_YX_%0.3f.npy" % (input_sky,freq_MHz)      
         auto_array_filename_X = "%s_auto_array_X_%0.3f.npy" % (input_sky,freq_MHz)
         auto_array_filename_Y = "%s_auto_array_Y_%0.3f.npy" % (input_sky,freq_MHz)
      else:
         cross_visibility_real_array_filename_XX = "%s_cross_visibility_real_array_%s_XX_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_imag_array_filename_XX = "%s_cross_visibility_imag_array_%s_XX_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_real_array_filename_YY = "%s_cross_visibility_real_array_%s_YY_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_imag_array_filename_YY = "%s_cross_visibility_imag_array_%s_YY_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_real_array_filename_XY=  "%s_cross_visibility_real_array_%s_XY_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_imag_array_filename_XY = "%s_cross_visibility_imag_array_%s_XY_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_real_array_filename_YX = "%s_cross_visibility_real_array_%s_YX_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         cross_visibility_imag_array_filename_YX = "%s_cross_visibility_imag_array_%s_YX_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         auto_array_filename_X = "%s_auto_array_%s_X_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
         auto_array_filename_Y = "%s_auto_array_%s_Y_%0.3f.npy" % (input_sky,EDA2_obs_time,freq_MHz)
      
      uu_array_filename = "uu_array_%0.3f.npy" % (freq_MHz) 
      vv_array_filename = "vv_array_%0.3f.npy" % (freq_MHz) 
      ww_array_filename = "ww_array_%0.3f.npy" % (freq_MHz) 
      baseline_number_array_filename = "baseline_number_array_%0.3f.npy" % (freq_MHz) 
      

      cross_visibility_real_array_XX = np.load(cross_visibility_real_array_filename_XX)
      cross_visibility_imag_array_XX = np.load(cross_visibility_imag_array_filename_XX)   
      cross_visibility_real_array_YY = np.load(cross_visibility_real_array_filename_YY)
      cross_visibility_imag_array_YY = np.load(cross_visibility_imag_array_filename_YY)     
      cross_visibility_real_array_XY = np.load(cross_visibility_real_array_filename_XY)
      cross_visibility_imag_array_XY = np.load(cross_visibility_imag_array_filename_XY) 
      cross_visibility_real_array_YX = np.load(cross_visibility_real_array_filename_YX)
      cross_visibility_imag_array_YX = np.load(cross_visibility_imag_array_filename_YX) 
      
      auto_array_X = np.load(auto_array_filename_X)
      auto_array_Y = np.load(auto_array_filename_Y)
      
      #why did I do this?
      #cross_visibility_complex_array_XX = cross_visibility_real_array_XX + 1j*cross_visibility_imag_array_XX
      #cross_visibility_complex_array_YY = cross_visibility_real_array_YY + 1j*cross_visibility_imag_array_YY
      #cross_visibility_complex_array_XY = cross_visibility_real_array_XY + 1j*cross_visibility_imag_array_XY
      #cross_visibility_complex_array_YX = cross_visibility_real_array_YX + 1j*cross_visibility_imag_array_YX

      #cross_visibility_real_array_XX = np.real(cross_visibility_complex_array_XX)
      #cross_visibility_imag_array_XX = np.imag(cross_visibility_complex_array_XX)
      #cross_visibility_real_array_YY = np.real(cross_visibility_complex_array_YY)
      #cross_visibility_imag_array_YY = np.imag(cross_visibility_complex_array_YY) 
      #cross_visibility_real_array_XY = np.real(cross_visibility_complex_array_XY)
      #cross_visibility_imag_array_XY = np.imag(cross_visibility_complex_array_XY)
      #cross_visibility_real_array_YX = np.real(cross_visibility_complex_array_YX)
      #cross_visibility_imag_array_YX = np.imag(cross_visibility_complex_array_YX)     
             
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
      
      if check_figs:
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
      
      ##look at the uvfits file
      #uvfits_filename = mir_file_uvfits_name
      #print("%s" % uvfits_filename)
      #hdulist = fits.open(uvfits_filename)
      ##hdulist.info()
      #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      ##print(info_string)
      #uvtable = hdulist[0].data
      #uvtable_header = hdulist[0].header
      #hdulist.close()
      #sim_UU_s_array = uvtable['UU']
      #sim_VV_s_array = uvtable['VV']
      #sim_WW_s_array = uvtable['WW']
      #sim_baselines = uvtable['baseline']
      #sim_visibilities = uvtable['DATA']
      #sim_visibilities_shape = sim_visibilities.shape
      #print('sim_UU_s_array')
      #print(sim_UU_s_array)
      #sim_UU_m_array = sim_UU_s_array * c
      #print('sim_UU_m_array')
      #print(sim_UU_m_array)
      #u_in_miriad = sim_UU_s_array[0]
      #ratio_u_in = u_in/u_in_miriad
      #print('ratio_u_in')
      #print(ratio_u_in)
      #print("sim visibilities_shape")
      #print(sim_visibilities_shape)
      #print('sim_baselines')
      #print(len(sim_baselines))
      #print(np.max(sim_baselines))
      
      ##what is going on with these baselines
      #for baseline_num_index,baseline_num in enumerate(sim_baselines):
      #   #print(sim_baselines[index])
      #   if baseline_num > (2**16):
      #      print(baseline_num)
      #      ant1,ant2 = decode_baseline(baseline_num)
      #      print(ant1,ant2)
      
      #got it to work by not adding the last antenna i.e above:
      #if ant2<255 and ant1<255:
            #uv.write(preamble,cross_vis)
      #need a better solution! but will be interesting to see if ms calibrates/images or if proper antenna coords are needed (I think yes for calibrate .. but wsclean might be okay ....)
   
      #print(len(sim_baseline))
      #print(len(sim_UU_s_array))
      ##convert from miriad baseline numbers
      ##converted_sim_baseline_ant_nums = [aa.bl2ij(bl) for bl in sim_baselines]
      ##converted_sim_baseline_ant_nums = [decode_baseline(bl) for bl in sim_baselines]
      #converted_sim_baselines = [(cal_standard_baseline_number(ant1,ant2)) for ant1,ant2 in converted_sim_baseline_ant_nums]
      #converted_sim_baselines = np.asarray(converted_sim_baselines,dtype='f')
      ##print(converted_sim_baselines)
   
      #try to read the data in to casa as a ms
      #verify with wsclean - too many antennas - 255 max?
      #read in the uvfits file
      casa_cmd_filename = 'import_uvfits.sh'
      cmd = "rm -rf %s" % (mir_file_ms_name)
      print(cmd)
      os.system(cmd)
            
      cmd = "importuvfits(fitsfile='%s',vis='%s')" % (mir_file_uvfits_name,mir_file_ms_name)
 
      with open(casa_cmd_filename,'w') as f:
         f.write(cmd)
           
      cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
      print(cmd)
      os.system(cmd)
      
      if check_figs:
         #####verification stuff 
         test_image_name = mir_file_ms_name.split('.ms')[0]
         wsclean_imsize = '512'
         wsclean_scale = '900asec'
         cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column DATA  %s " % (test_image_name,wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
         print(cmd)
         os.system(cmd) 
         
         ##try adding a model column
         #use kariukis ms_utils
         #ms_table = table(mir_file_ms_name,readonly=False)
         #model_data = get_data(ms_table, col="DATA")
         #add_col(ms_table, "MODEL_DATA")
         #put_col(ms_table, "MODEL_DATA", model_data)
         #ms_table.close()
   
         ##try imaging the model column of the ms:
         #wsclean_imsize = '512'
         #wsclean_scale = '900asec'
         #cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column MODEL_DATA  %s " % (test_image_name+"_model",wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
         #print(cmd)
         #os.system(cmd)
         
         ##try calibrate?
         #gain_solutions_name = 'test_cal_%s_%s_calibrate_sols.bin' % (EDA2_chan,EDA2_obs_time)
         #calibrate_options = ''
         #cmd = "rm -rf %s" % (gain_solutions_name)
         #print(cmd)
         #os.system(cmd)  
         ##calibrate
         #cmd = "calibrate %s %s %s " % (calibrate_options,mir_file_ms_name,gain_solutions_name)
         #print(cmd)
         #os.system(cmd)
         ##plot cal sols
         #cmd = "aocal_plot.py %s  " % (gain_solutions_name)
         #print(cmd)
         #os.system(cmd)
         
         #cmd = "applysolutions %s %s  " % (mir_file_ms_name,gain_solutions_name)
         #print(cmd)
         #os.system(cmd)
         
         #test image the CORRECTED data'
         #cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (test_image_name+"_corrected",wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
         #print(cmd)
         #os.system(cmd)
         
         #check cal
         #calibrate_with_complex_beam_model(model_ms_name,eda2_ms_name)
         
         
         #if pol=='X':
         #   wsclean_imsize = '512'
         #   wsclean_scale = '900asec'
         #   cmd = "wsclean -name test_eda_%0.3f_%s_wsclean -size %s %s -multiscale -niter 1000 -scale %s -pol xx  %s " % (freq_MHz,pol,wsclean_imsize,wsclean_imsize,wsclean_scale,mir_file_ms_name)
         #   print(cmd)
         #   os.system(cmd) 
         #   
   
         ##Lets look at some data at the same LST
         #uvfits_filename = "/md0/EoR/EDA2/20200303_data/%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
         #ms_name = "/md0/EoR/EDA2/20200303_data/%s/av_chan_%s_%s_plus_%s_obs.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
         #print("%s" % uvfits_filename)
         #hdulist = fits.open(uvfits_filename)
         ##hdulist.info()
         #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
         ##print(info_string)
         #uvtable = hdulist[0].data
         #uvtable_header = hdulist[0].header
         #hdulist.close()
         #data_UU_s_array = uvtable['UU']
         #data_VV_s_array = uvtable['VV']
         #data_WW_s_array = uvtable['WW']
         #data_baselines = uvtable['baseline']
         #data_visibilities = uvtable['DATA']
         #data_visibilities_shape = data_visibilities.shape
         #print("data visibilities_shape")
         #print(data_visibilities_shape)
         ##print(len(data_baseline))
         #print('data_UU_s_array')
         #print(data_UU_s_array)
         #data_UU_m_array = data_UU_s_array * c
         #print('data_UU_m_array')
         #print(data_UU_m_array)
         ##print(np.max(data_baselines))
         ##print(np.min(data_baselines))
   
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
      #n_common = np.count_nonzero(model_ms_indices) 
      #print(n_common)
      eda2_ms_indices = inNd(eda2_ants, model_ants, assume_unique=False)
      #n_common = np.count_nonzero(eda2_ms_indices)
      #print(n_common)
      
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

def calibrate_with_complex_beam_model_time_av(EDA2_chan_list,lst_list=[],plot_cal=False,uv_cutoff=0,per_chan_cal=False,EDA2_data=True,sim_only_EDA2=[],run_aoflagger=True,rfi_strategy_name='/md0/code/git/ben-astronomy/EoR/ASSASSIN/rfi_strategy_new.rfis'):
   if len(sim_only_EDA2)!=0:
      sim_only=True
   #think about how to use coherence:
   #coherence = cross12_avg/(np.sqrt(auto11_avg*auto22_avg))
   
   print("concatenate and flag EDA2 obs before calibration")
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):  
      #freq_MHz = np.round(400./512.*float(EDA2_chan))
      freq_MHz = 400./512.*float(EDA2_chan)
      wavelength = 300./freq_MHz
      if len(EDA2_chan_list)==1:
         obs_time_list = EDA2_obs_time_list_each_chan[chan_num]
      else:
         obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
      first_obstime = obs_time_list[0]
      n_obs_concat = len(obs_time_list)
      #concat_ms_name = "%s/%s_plus_%s_obs_%0.3f_MHz_concat.ms" % (EDA2_chan,first_obstime,n_obs_concat,freq_MHz)
      concat_ms_name = "%s/concat_chan_%s.ms" % (EDA2_chan,EDA2_chan)
      obs_concat_list = []
      for EDA2_obs_time_index,EDA2_obs_time in enumerate(obs_time_list):
         eda2_ms_name = "%s/%s_%s_eda2_ch32_ant256_midday_avg8140.ms" % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
         print(eda2_ms_name)
         obs_concat_list.append(eda2_ms_name)
         
         #add in the model column to each ms
         model_ms_name = "%s_%0.3f.ms" % (first_obstime,freq_MHz) 
         print(model_ms_name)
   
         eda2_ms_table = table(eda2_ms_name,readonly=False)
         eda2_data = get_data(eda2_ms_table)
         eda2_uvw = get_uvw(eda2_ms_table)
         eda2_ant1, eda2_ant2 = get_ant12(eda2_ms_name)
         eda2_ants = np.vstack((eda2_ant1,eda2_ant2)).T
      
         model_ms_table = table(model_ms_name,readonly=True)
         model_data = get_data(model_ms_table)
         model_uvw = get_uvw(model_ms_table)    
         model_ant1, model_ant2 = get_ant12(model_ms_name)
         model_ants = np.vstack((model_ant1,model_ant2)).T
   
         model_ms_indices = inNd(model_ants, eda2_ants, assume_unique=False)
         #n_common = np.count_nonzero(model_ms_indices) 
         #print(n_common)
         eda2_ms_indices = inNd(eda2_ants, model_ants, assume_unique=False)
         #n_common = np.count_nonzero(eda2_ms_indices)
         #print(n_common)
              
         eda2_common_ant1 = eda2_ant1[eda2_ms_indices]
         eda2_common_ant2 = eda2_ant2[eda2_ms_indices]
         eda2_common_sort_inds = np.lexsort((eda2_common_ant2,eda2_common_ant1)) # Sort by a, then by b
   
         model_common_ant1 = model_ant1[model_ms_indices]
         model_common_ant2 = model_ant2[model_ms_indices]
         model_common_sort_inds = np.lexsort((model_common_ant2,model_common_ant1)) # Sort by a, then by b
               
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
         ###print(new_common_model_data_sorted.shape)
         repetitions = 32
         new_common_model_data_sorted_tile = np.tile(new_common_model_data_sorted, (repetitions, 1))
         #print(new_common_model_data_sorted_tile.shape)
      
        
         try:
            add_col(eda2_ms_table, "MODEL_DATA")
         except:
            pass
         put_col(eda2_ms_table, "MODEL_DATA", new_common_model_data_sorted_tile)
         eda2_ms_table.close() 
         
      obs_concat_list_string = "','".join(obs_concat_list)
      
      casa_cmd_filename = '%s/concat_ms.sh' % EDA2_chan
      
      cmd = "rm -rf %s %s" % (concat_ms_name,casa_cmd_filename)
      print(cmd)
      os.system(cmd)
            
      cmd = "concat(vis=['%s'],concatvis='%s',dirtol='5arcmin')" % (obs_concat_list_string,concat_ms_name)
      
      with open(casa_cmd_filename,'w') as f:
         f.write(cmd)
           
      cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
      print(cmd)
      os.system(cmd)
      
      if run_aoflagger:
         #flag data with aoflagger
         cmd = "aoflagger -strategy %s -column DATA %s" % (rfi_strategy_name,concat_ms_name)
         print(cmd)
         os.system(cmd)
   
      #calibrate
      if uv_cutoff==0:
         gain_solutions_name = '%s/%s_%s_complex_beam_cal_sols_av.bin' % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
         calibrate_options = ''
      else:
         gain_solutions_name = '%s/%s_%s_complex_beam_cal_sols_%0.3f_m_av.bin' % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15],uv_cutoff_m)
         calibrate_options = '-minuv %0.3f ' % uv_cutoff_m
     
      cmd = "rm -rf %s" % (gain_solutions_name)
      print(cmd)
      os.system(cmd)  
      #calibrate
      cmd = "calibrate -ch 1 %s %s %s " % (calibrate_options,concat_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      #plot cal sols
      if plot_cal:
         #Plot the cal solutions
         cmd = "aocal_plot.py  %s " % (gain_solutions_name)
         print(cmd)
         os.system(cmd)
          
      cmd = "applysolutions %s %s  " % (concat_ms_name,gain_solutions_name)
      print(cmd)
      os.system(cmd)
      
      ##test image the CORRECTED data'
      #wsclean_imsize = '512'
      #wsclean_scale = '900asec' 
      #cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (test_image_name+"_corrected_av",wsclean_imsize,wsclean_imsize,wsclean_scale,concat_ms_name)
      #print(cmd)
      #os.system(cmd)
   
      ###flag corrected data? - Yes! This worked well
      if run_aoflagger:
         cmd = "aoflagger -strategy %s -column CORRECTED_DATA %s" % (rfi_strategy_name,concat_ms_name)
         print(cmd)
         os.system(cmd)
      
      #calibrate again with revised flags? Image and then self-cal? These are all good idea - but look at individual time/freq chans first
      
      
      
      #Forget this averaging stuff - can't work out how to average over different obsids with CASA
      #just use concat ms and average in extract_tsky...
      ##average with CASA mstransform (or split?)
      #av_ms_name = "%s/av_chan_%s_%s_plus_concat.ms" % (EDA2_chan,EDA2_chan,first_obstime)
      #
      #casa_cmd_filename = '%s/av_ms.sh' % EDA2_chan
      #
      #cmd = "rm -rf %s %s" % (av_ms_name,casa_cmd_filename)
      #print(cmd)
      #os.system(cmd)
      #      
      ##cmd = "mstransform(vis='%s',outputvis='%s',timeaverage=True,timebin='10s',datacolumn='all')" % (concat_ms_name,av_ms_name)
      #
      # 
      #with open(casa_cmd_filename,'w') as f:
      #   f.write(cmd)
      #     
      #cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
      #print(cmd)
      #os.system(cmd)

      
def old_calibrate():
      ##run calibrate on flagged but not averaged data.   
      print("averaging EDA2 obs in time before calibration")
      ##specify uv_cutoff in wavelengths, convert to m for 'calibrate'
      ##pol = pol_list[0]
      #gsm  = GlobalSkyModel()
      #wsclean_imsize = '512'
      #wsclean_scale = '900asec'  
      # 
      #number_of_good_obs_list = []
      #number_of_good_obs_list_filename = "number_of_good_obs_list.txt"
      ##if not sim_only:
       
      #for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):  
      #   #freq_MHz = np.round(400./512.*float(EDA2_chan))
      #   freq_MHz = 400./512.*float(EDA2_chan)
      #   wavelength = 300./freq_MHz
      #   centre_freq = float(freq_MHz)
      #   fine_chan_width_MHz = fine_chan_width_Hz/1000000.  
      #   if uv_cutoff!=0:
      #      uv_cutoff_m = uv_cutoff * wavelength
      #   lst = lst_list[EDA2_chan_index]
      #   lst_deg = (float(lst)/24)*360.
      #   if len(EDA2_chan_list)==1:
      #      obs_time_list = EDA2_obs_time_list_each_chan[chan_num]
      #   else:
      #      obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
      #   first_obstime = obs_time_list[0]
      #   
      #   #guard against cases where there are no data for that channel
      #   if first_obstime==0:
      #      continue
      #   else:
      #      pass
   
      #   #Don't do the time averaging yet, try calibrating each obs individually so you can weed out the bad ones
      #   #No need for a multi-freq model cube for each obs, just use centre freq 
      #   wsclean_cal_ms_name_list = []
      #   for EDA2_obs_time_index,EDA2_obs_time in enumerate(obs_time_list):
      #      eda2_ms_name = "%s/%s_%s_eda2_ch32_ant256_midday_avg8140.ms" % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
      #      print(eda2_ms_name)
      #      #model_ms_name = "/md0/EoR/EDA2/EEPs/new_20210616/%s_%0.3f.ms" % (first_obstime,freq_MHz) 
      #      #now just in same directory as unity ms
      #      model_ms_name = "%s_%0.3f.ms" % (first_obstime,freq_MHz) 
      #      print(model_ms_name)
      #      
      #      eda2_ms_table = table(eda2_ms_name,readonly=False)
      #      #eda2_table_summary = tablesummary(eda2_ms_name)
      #      #print(eda2_table_summary)
      #      eda2_data = get_data(eda2_ms_table)
      #      #print(eda2_data.shape)
      #      eda2_uvw = get_uvw(eda2_ms_table)
      #      eda2_ant1, eda2_ant2 = get_ant12(eda2_ms_name)
      #      eda2_ants = np.vstack((eda2_ant1,eda2_ant2)).T
   
      #    
      #     model_ms_table = table(model_ms_name,readonly=True)
      #      model_data = get_data(model_ms_table)
      #      model_uvw = get_uvw(model_ms_table)    
      #      model_ant1, model_ant2 = get_ant12(model_ms_name)
      #      model_ants = np.vstack((model_ant1,model_ant2)).T
      #
      #      #model_ms_indices = np.nonzero(np.in1d(model_uvw, eda2_uvw))[0]
      #      model_ms_indices = inNd(model_ants, eda2_ants, assume_unique=False)
      #      #n_common = np.count_nonzero(model_ms_indices) 
      #      #print(n_common)
      #      eda2_ms_indices = inNd(eda2_ants, model_ants, assume_unique=False)
      #      #n_common = np.count_nonzero(eda2_ms_indices)
      #      #print(n_common)
      #      
      #      
      #      #ind = np.lexsort((b,a)) # Sort by a, then by b
      #      #common_eda2_uvw_sorted = np.sort(common_eda2_uvw,axis=0)
      #      
      #      eda2_common_ant1 = eda2_ant1[eda2_ms_indices]
      #      eda2_common_ant2 = eda2_ant2[eda2_ms_indices]
      #      eda2_common_sort_inds = np.lexsort((eda2_common_ant2,eda2_common_ant1)) # Sort by a, then by b
   
            
      #      model_common_ant1 = model_ant1[model_ms_indices]
      #      model_common_ant2 = model_ant2[model_ms_indices]
      #      model_common_sort_inds = np.lexsort((model_common_ant2,model_common_ant1)) # Sort by a, then by b
            
      #      #don't sort the ant arrays, these are what we use to sort the other arrays!
      #      #model_common_ant1_sorted = model_common_ant1[model_common_sort_inds]
      #      #model_common_ant2_sorted = model_common_ant2[model_common_sort_inds]
                  
      #      common_eda2_uvw = eda2_uvw[eda2_ms_indices]
      #      common_eda2_uvw_sorted = common_eda2_uvw[eda2_common_sort_inds]
            
      #      common_eda2_data = eda2_data[eda2_ms_indices]
      #      common_eda2_data_sorted = common_eda2_data[eda2_common_sort_inds]
      #      #print(common_eda2_data_sorted.shape)
      #      #print(common_eda2_data_sorted[0:10,0,0])
            
      #      common_model_uvw = model_uvw[model_ms_indices]
      #      common_model_uvw_sorted = common_model_uvw[model_common_sort_inds]
      
      #      common_model_data = model_data[model_ms_indices]
      #      common_model_data_sorted = common_model_data[model_common_sort_inds]
      #      #print(common_model_data_sorted.shape)
      #      #print(common_model_data_sorted[0:10,0,0])      
            
      #      #for model go through ant2 array and data array and uvw array, insert zeros after wherever ant==254
      #      ant2_list=[]
      #      for ant2_index,ant2 in enumerate(model_common_ant2):
      #         if ant2==254:
      #            ant2_list.append(ant2_index)
      #      #print(len(ant2_list))
      #      #print(ant2_list)
            
      #      old_model_common_ant2 = model_common_ant2
      #      old_model_common_ant1 = model_common_ant1
      #      old_common_model_data_sorted = common_model_data_sorted
      #      counter=0
            
      #      #a = np.array([[1, 1], [2, 2], [3, 3]])
      #      #print(a)
      #      #b = np.insert(a, 1, 5, axis=0)
      #      #print(b)
            
            
      #      for ant2_index in ant2_list:
      #         ant2_index+=counter
      #         ant1_value = old_model_common_ant1[ant2_index-1]
      #         #new_data_value = np.zeros((1,4))
      #         #print(new_data_value.shape)
      #         #print(old_array.shape)
      #         new_model_common_ant2 = np.insert(old_model_common_ant2,ant2_index+1,255)
      #         new_model_common_ant1 = np.insert(old_model_common_ant1,ant2_index+1,ant1_value)
      #         new_common_model_data_sorted = np.insert(old_common_model_data_sorted,ant2_index+1,0,axis=0)
      #         #print(old_common_model_data_sorted[ant2_index-2:ant2_index+5])
      #         #print(new_common_model_data_sorted[ant2_index-2:ant2_index+6])
      #         #print(new_array[ant2_index+1])
      #         old_model_common_ant2 = new_model_common_ant2
      #         old_model_common_ant1 = new_model_common_ant1
      #         old_common_model_data_sorted = new_common_model_data_sorted
      #          counter+=1
      #      new_model_common_ant2 = np.append(new_model_common_ant2,np.array([254,255,255]))
      #      new_model_common_ant1 = np.append(new_model_common_ant1,np.array([254,254,255]))
      #      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      #      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      #      new_common_model_data_sorted = np.append(new_common_model_data_sorted,np.array([[[0,0,0,0]]]),axis=0)
      #      ###print(new_common_model_data_sorted.shape)
      #      repetitions = 32
      #      new_common_model_data_sorted_tile = np.tile(new_common_model_data_sorted, (repetitions, 1))
      #      #print(new_common_model_data_sorted_tile.shape)
            
      #      #print(len(new_model_common_ant2_sorted))
      #      #print((new_model_common_ant2_sorted[-200:-1]))
      #      #print(eda2_ant2.shape)
      #      #print(eda2_ant2[-200:-1])
            
      #      #print(common_model_data_sorted.shape)
      #      #print(common_model_data_sorted[-5:])
      #      #print(new_common_model_data_sorted.shape)
      #      #print((new_common_model_data_sorted[-5:]))
      
      #      #model has no ant 255 and no autos, data is missing some other antennas, but has 255 and autos
      #      #going to need to add autos into model, and to add in missing ant 255 correlations as dummy data, and then flag that ant in the data before calibration
      #      #arrrgghhhh
      #      try:
      #         add_col(eda2_ms_table, "MODEL_DATA")
      #      except:
      #         pass
      #      put_col(eda2_ms_table, "MODEL_DATA", new_common_model_data_sorted_tile)
      #      eda2_ms_table.close()     
      
      #
      #      ##try imaging the model column of the ms:
      #      #wsclean_imsize = '512'
      #      #wsclean_scale = '900asec'
      #      test_image_name = "complex_beam_test"
      #      #cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column MODEL_DATA  %s " % (test_image_name+"_model",wsclean_imsize,wsclean_imsize,wsclean_scale,eda2_ms_name)
      #      #print(cmd)
      #      #os.system(cmd)
            
      #      #need to flag ant 255 .... also flux scale is about half (not treating polarisation correctly?)
      #      #use uv cutoff at half wavelength?
      #      #try cal with -ch 32
            
      #      #try calibrate
      #      if uv_cutoff==0:
      #         gain_solutions_name = '%s/%s_%s_complex_beam_cal_sols.bin' % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
      #         calibrate_options = ''
      #      else:
      #         gain_solutions_name = '%s/%s_%s_complex_beam_cal_sols_%0.3f_m.bin' % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15],uv_cutoff_m)
      #         calibrate_options = '-minuv %0.3f ' % uv_cutoff_m
      #  
      #      cmd = "rm -rf %s" % (gain_solutions_name)
      #      print(cmd)
      #      os.system(cmd)  
      #      #calibrate
      #      cmd = "calibrate -ch 32 %s %s %s " % (calibrate_options,eda2_ms_name,gain_solutions_name)
      #      print(cmd)
      #      os.system(cmd)
      #      #plot cal sols
            
      #      #plot the sols and 
      #      if (os.path.isfile(gain_solutions_name)):
      #         wsclean_cal_ms_name_list.append(eda2_ms_name)
      #       
      #         if plot_cal:
      #            #Plot the cal solutions
      #            cmd = "aocal_plot.py  %s " % (gain_solutions_name)
      #            print(cmd)
      #            #os.system(cmd)
      #             
      #         cmd = "applysolutions %s %s  " % (eda2_ms_name,gain_solutions_name)
      #         print(cmd)
      #         os.system(cmd)
      #      
      #         #test image the CORRECTED data'
      #         #cmd = "wsclean -name %s -size %s %s -multiscale -weight briggs 0 -niter 500 -scale %s -pol xx,yy -data-column CORRECTED_DATA  %s " % (test_image_name+"_corrected",wsclean_imsize,wsclean_imsize,wsclean_scale,eda2_ms_name)
      #         #print(cmd)
      #         #os.system(cmd)
      #   
      #   number_of_good_obs = len(wsclean_cal_ms_name_list)
      #   print("number of good obs used in chan %s is %s" % (EDA2_chan,number_of_good_obs)) 
      #   number_of_good_obs_list.append(str(number_of_good_obs))
      #    
      #   #print(wsclean_cal_ms_name_list)
      #   #now average
      #   av_ms_name = "%s/av_chan_%s_%s_plus_%s_obs_complex.ms" % (EDA2_chan,EDA2_chan,first_obstime,number_of_good_obs)
      #
      #   for ms_name_index,ms_name in enumerate(wsclean_cal_ms_name_list):    
      #      #need to use kariukes ms_utils instead (better anyway!)
      #      
      #      #if this is the first observation create a copy of the ms to be the sum ms
      #      if(ms_name_index==0):
      #         cmd = 'rm -rf %s' % (av_ms_name)
      #         print(cmd)
      #         os.system(cmd)
      #         
      #         cmd = 'cp -r %s %s' % (ms_name,av_ms_name)
      #         print(cmd)
      #         os.system(cmd)
      #         
      #         sum_ms_table = table(ms_name,readonly=True)
      #         sum_ms_data = get_data(sum_ms_table, col="DATA")
      #         #print(sum_ms_data.shape)
      #         sum_ms_table.close()
   

def extract_global_signal_from_ms_complex(EDA2_chan_list=[],lst_list=[],uvdist_thresh_lambda=0.5,n_fine_chans_included=27):
   #max umber of timesteps is set to eight 
   max_n_timesteps = 8
   n_fine_chans = int(n_fine_chans_included*len(EDA2_chan_list))
   t_sky_K_array = np.zeros(len(EDA2_chan_list))
   t_sky_error_K_array = np.zeros(len(EDA2_chan_list))
   t_sky_K_array_flagged = np.zeros(len(EDA2_chan_list))
   t_sky_error_K_array_flagged = np.zeros(len(EDA2_chan_list))
   
   t_sky_K_array_filename = "t_sky_K_array_eda2.npy"
   t_sky_error_K_array_filename = "t_sky_error_K_array_eda2.npy"
   t_sky_K_array_filename_flagged = "t_sky_K_array_eda2_flagged.npy"
   t_sky_error_K_array_filename_flagged = "t_sky_error_K_array_eda2_flagged.npy"   
   
   t_sky_K_array_fine_chan = np.empty([len(EDA2_chan_list),n_fine_chans_included,max_n_timesteps])
   t_sky_error_K_array_fine_chan = np.empty([len(EDA2_chan_list),n_fine_chans_included,max_n_timesteps])
   t_sky_K_array_flagged_fine_chan = np.empty([len(EDA2_chan_list),n_fine_chans_included,max_n_timesteps])
   t_sky_error_K_array_flagged_fine_chan = np.empty([len(EDA2_chan_list),n_fine_chans_included,max_n_timesteps])
   t_sky_K_array_fine_chan[:] = np.nan
   t_sky_error_K_array_fine_chan[:] = np.nan
   t_sky_K_array_flagged_fine_chan[:] = np.nan
   t_sky_error_K_array_flagged_fine_chan[:] = np.nan
 
   t_sky_K_array_filename_fine_chan = "t_sky_K_array_eda2_fine_chan.npy"
   t_sky_error_K_array_filename_fine_chan = "t_sky_error_K_array_eda2_fine_chan.npy"
   t_sky_K_array_filename_flagged_fine_chan = "t_sky_K_array_eda2_flagged_fine_chan.npy"
   t_sky_error_K_array_filename_flagged_fine_chan = "t_sky_error_K_array_eda2_flagged_fine_chan.npy"   
   
   eda2_chan_fine_chan_timestep_array = np.empty([len(EDA2_chan_list),n_fine_chans_included,max_n_timesteps])
   eda2_chan_fine_chan_timestep_array[:,:,:] = np.nan
   eda2_chan_fine_chan_timestep_array_fileame = "eda2_chan_fine_chan_timestep_array.npy"

   #number_of_good_obs_list_filename = "number_of_good_obs_list.txt"
   #with open(number_of_good_obs_list_filename) as f:
   #   number_of_good_obs_list_string = f.read()
   #number_of_good_obs_list = number_of_good_obs_list_string.split(',')
   
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):  
      freq_MHz = 400./512.*float(EDA2_chan)
      wavelength = c / (freq_MHz*1e6)
      print(wavelength)
      
      #This Jy to K has to be wrong, there is no solid angle ... I think just need to leave all the sims etc in units of K?
      #jy_to_K = (wavelength**2) / (2. * k * 1.0e26) 
      #print(jy_to_K)
      
      #number_of_good_obs = number_of_good_obs_list[EDA2_chan_index]
      lst = lst_list[EDA2_chan_index]
      lst_deg = (float(lst)/24)*360.
      
      if len(EDA2_chan_list)==1:
         obs_time_list = EDA2_obs_time_list_each_chan[chan_num]
      else:
         obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
      first_obstime = obs_time_list[0]
         
      #av_ms_name = "%s/av_chan_%s_%s_plus_concat.ms" % (EDA2_chan,EDA2_chan,first_obstime)
      av_ms_name = "%s/concat_chan_%s.ms" % (EDA2_chan,EDA2_chan)
      #concat_ms_name = "%s/concat_chan_%s.ms" % (EDA2_chan,EDA2_chan)
      #av_ms_name = "%s/av_chan_%s_%s_plus_5_obs_complex.ms" % (EDA2_chan,EDA2_chan,first_obstime)
      
      print("getting data from %s" % av_ms_name)
      unity_ms_name = "unity_%0.3f.ms" % (freq_MHz)
      print("and unity response from %s" % unity_ms_name)
      
      eda2_ms_table = table(av_ms_name,readonly=True)
      eda2_data_complex = get_data(eda2_ms_table, col="CORRECTED_DATA")
      n_vis_eda2 = eda2_data_complex.shape[0]
      
      eda2_flags = get_flags(eda2_ms_table)
      #make data where flags are true nans
      flagged_indices = np.where(eda2_flags)
      eda2_data_complex_flagged = np.copy(eda2_data_complex)
      eda2_data_complex_flagged[flagged_indices] = np.nan
      eda2_data_complex = eda2_data_complex_flagged
      
      #do this later on common data
      #for now, average 32 (central 28) chans to 1, look at only x pol and real
      #eda2_data_complex_freq_av = np.mean(eda2_data_complex[:,2:30,:],axis=1)
      #eda2_data_complex_freq_av_x = eda2_data_complex_freq_av[:,0]
      #eda2_data_complex_freq_av_x_real=eda2_data_complex_freq_av_x.real
      #print(eda2_data_complex_freq_av_x.shape)
      #print(eda2_data_complex_freq_av_x_real[1])
      #sys.exit()

      #now get expected unity response ms
      unity_ms_table = table(unity_ms_name,readonly=True)
      unity_data_complex = get_data(unity_ms_table, col="DATA")
      n_vis_1_timestep = unity_data_complex.shape[0]
      #do this later
      #unity_data_complex_freq_av_x = unity_data_complex[:,0,1]
      #unity_data_complex_freq_av_x_real=unity_data_complex_freq_av_x.real      
      
      #start_index = n_vis_1_timestep
      #end_index = n_vis_1_timestep*2
      #eda2_data_complex = eda2_data_complex[start_index:end_index,:,:]
      
      #sort out common baselines just like in calibrate_with_complex_beams
      eda2_uvw = get_uvw(eda2_ms_table)
      #eda2_uvw = eda2_uvw[start_index:end_index]
      eda2_ant1, eda2_ant2 = get_ant12(av_ms_name)
      #eda2_ant1 = eda2_ant1[start_index:end_index]
      #eda2_ant2 = eda2_ant2[start_index:end_index]
      
      #count number of ants in eda2
      n_ant_eda2 = 0
      for ant_num_index, ant_num in enumerate(eda2_ant1):
         n_ant_eda2 += 1
         if (eda2_ant1[ant_num_index] != eda2_ant1[ant_num_index+1]):
            break
      n_baselines_eda2 = int((n_ant_eda2*(n_ant_eda2-1) / 2) + n_ant_eda2)

      n_timesteps_eda2 = int(n_vis_eda2/n_baselines_eda2)
      print("n_timesteps_eda2 %s" % (n_timesteps_eda2))


      eda2_ant1 = eda2_ant1[0:n_baselines_eda2]
      eda2_ant2 = eda2_ant2[0:n_baselines_eda2]
      eda2_ants = np.vstack((eda2_ant1,eda2_ant2)).T
      eda2_uvw = eda2_uvw[0:n_baselines_eda2]
      
      unity_uvw = get_uvw(unity_ms_table)    
      unity_ant1, unity_ant2 = get_ant12(unity_ms_name)
      unity_ants = np.vstack((unity_ant1,unity_ant2)).T
      
      ##lets just plot unity unity real vs unity uv dist without any sorting or cutting
      #unity_uvw_dist = np.sqrt(unity_uvw[:,0]**2 + unity_uvw[:,1]**2 + unity_uvw[:,2]**2)
      #print(unity_uvw_dist.shape)
      #unity_data_complex_x = unity_data_complex[:,0,0]
      #unity_data_complex_x_real = unity_data_complex_x.real
      #print(unity_data_complex_x_real.shape)

      ##plot the expected unity response vs baseline length
      #plt.clf()
      #plt.scatter(unity_uvw_dist,unity_data_complex_x_real)            
      #map_title="Unity_response_vs_uvdist" 
      #plt.xlabel("UV distance (wavelengths)")
      #plt.ylabel("Real component of visibility X pol (Jy)")
      #plt.legend(loc=1)
      ##plt.text(x_pos, y_pos, fit_string)
      ##plt.ylim([0, 3.5])
      #fig_name= "test_unity_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #plt.close()
      #print("saved %s" % fig_name)
      #sys.exit()

      unity_ms_indices = inNd(unity_ants, eda2_ants, assume_unique=False)
      #n_common = np.count_nonzero(unity_ms_indices) 
      #print(n_common)
      eda2_ms_indices = inNd(eda2_ants, unity_ants, assume_unique=False)
      #n_common = np.count_nonzero(eda2_ms_indices)
      #print(n_common)

      #print_index_end = 256
      eda2_common_ant1 = eda2_ant1[eda2_ms_indices]
      eda2_common_ant2 = eda2_ant2[eda2_ms_indices]
      #print(eda2_common_ant1[0:print_index_end])
      #print(eda2_common_ant2[0:print_index_end])
      eda2_common_sort_inds = np.lexsort((eda2_common_ant2,eda2_common_ant1)) # Sort by a, then by b

      unity_common_ant1 = unity_ant1[unity_ms_indices]
      unity_common_ant2 = unity_ant2[unity_ms_indices]
      #print(unity_common_ant1[0:print_index_end])
      #print(unity_common_ant2[0:print_index_end])
      unity_common_sort_inds = np.lexsort((unity_common_ant2,unity_common_ant1)) # Sort by a, then by b
      
      common_eda2_uvw = eda2_uvw[eda2_ms_indices]
      common_eda2_uvw_sorted = common_eda2_uvw[eda2_common_sort_inds]
      #print(common_eda2_uvw_sorted[0:10])
      #print(common_eda2_uvw_sorted.shape)
      #print(common_eda2_uvw_sorted[0:10])

      
      common_unity_uvw = unity_uvw[unity_ms_indices]
      common_unity_uvw_sorted = common_unity_uvw[unity_common_sort_inds]
      #print(common_unity_uvw_sorted[0:10])
      #print(common_unity_uvw_sorted.shape)
      #print(common_unity_uvw_sorted[0:10])
      common_unity_data = unity_data_complex[unity_ms_indices]
      common_unity_data_sorted = common_unity_data[unity_common_sort_inds]
      #print(common_unity_data_sorted[0:10])
      #print(common_unity_data_sorted.shape)
      #print(common_unity_data_sorted[0:10,0,0])   

      ###lets just plot unity unity real vs unity uv dist of the matched data 
      #unity_uvw_dist = np.sqrt(common_unity_uvw_sorted[:,0]**2 + common_unity_uvw_sorted[:,1]**2 + common_unity_uvw_sorted[:,2]**2)
      #print(unity_uvw_dist.shape)
      #unity_data_complex_x = common_unity_data_sorted[:,0,0]
      #unity_data_complex_x_real = unity_data_complex_x.real
      #print(unity_data_complex_x_real.shape)

      ##plot the expected unity response vs baseline length
      #plt.clf()
      #plt.scatter(unity_uvw_dist,unity_data_complex_x_real)            
      #map_title="Unity_response_vs_uvdist" 
      #plt.xlabel("UV distance (wavelengths)")
      #plt.ylabel("Real component of visibility X pol (Jy)")
      #plt.legend(loc=1)
      ##plt.text(x_pos, y_pos, fit_string)
      ##plt.ylim([0, 3.5])
      #fig_name= "test_unity_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #plt.close()
      #print("saved %s" % fig_name)
 
 
      uvdist_array_unity_m = np.sqrt(common_unity_uvw_sorted[:,0]**2 + common_unity_uvw_sorted[:,1]**2 + common_unity_uvw_sorted[:,2]**2)
      uvdist_array_unity_lambda = uvdist_array_unity_m / wavelength
      #we don't want the auto included here
      uvdist_array_unity_lambda_no_zero = uvdist_array_unity_lambda[uvdist_array_unity_lambda>0.]
      uvdist_array_unity_inds = uvdist_array_unity_lambda_no_zero.argsort()
      uvdist_array_unity_lambda_sorted = uvdist_array_unity_lambda_no_zero[uvdist_array_unity_inds]

      #uvdist_array_unity_lambda_sorted_no_zero = uvdist_array_unity_lambda_sorted[uvdist_array_unity_lambda_sorted>0]
      #print(uvdist_array_unity_lambda_sorted_no_zero[615:630])
      uvdist_array_unity_lambda_sorted_cut = uvdist_array_unity_lambda_sorted[uvdist_array_unity_lambda_sorted < uvdist_thresh_lambda]
        
      n_baselines_included_unity = len(uvdist_array_unity_lambda_sorted_cut)
      print("number of baselines included unity %s" % n_baselines_included_unity)

      #just looking at real,x for now
      common_unity_data_complex_freq_av_x = common_unity_data_sorted[:,0,0]
      common_unity_data_complex_freq_av_x_real=common_unity_data_complex_freq_av_x.real
      
      #remove the autos
      common_unity_data_complex_freq_av_x_real_sorted_no_auto = common_unity_data_complex_freq_av_x_real[uvdist_array_unity_lambda>0]
      common_unity_data_complex_freq_av_x_real_sorted = common_unity_data_complex_freq_av_x_real_sorted_no_auto[uvdist_array_unity_inds]
      
      #Instead of averaging in time and frequency, do all times and freqs separately to look for bad times/chans in x_y_plots
      start_index = 0
      end_index = n_baselines_eda2
      for timestep in range(n_timesteps_eda2):
         timestep_data = eda2_data_complex[start_index:end_index,:,:]
         start_index += n_baselines_eda2
         end_index += n_baselines_eda2
         for fine_chan_included in range(n_fine_chans_included):
            fine_chan_index = fine_chan_included+2
            timestep_data_fine_chan = timestep_data[:,fine_chan_index,:]
            base_name = "%s_timestep_%s_chan_%03d_fine_%02d" % (first_obstime,timestep,EDA2_chan,fine_chan_index)
            common_eda2_data = timestep_data_fine_chan[eda2_ms_indices]
            common_eda2_data_sorted = common_eda2_data[eda2_common_sort_inds]
            common_eda2_data_complex_x = common_eda2_data_sorted[:,0]
            common_eda2_data_complex_x_real=common_eda2_data_complex_x.real
            #there will be nans in here, but not in the unity - what happens? I think fine, they get droppen later an in the fitting
            #common_eda2_data_complex_freq_av_x_real_no_nans = common_eda2_data_complex_freq_av_x_real[~np.isnan(common_eda2_data_complex_freq_av_x_real)]
      
      
            #Need to sort by baseline length (then only use short baselines)
            uvdist_array_eda2_m = np.sqrt(common_eda2_uvw_sorted[:,0]**2 + common_eda2_uvw_sorted[:,1]**2 + common_eda2_uvw_sorted[:,2]**2)
            uvdist_array_eda2_lambda = uvdist_array_eda2_m / wavelength
            #we don't want the auto included here
            uvdist_array_eda2_lambda_no_zero = uvdist_array_eda2_lambda[uvdist_array_eda2_lambda>0.]
            uvdist_array_eda2_inds = uvdist_array_eda2_lambda_no_zero.argsort()
            uvdist_array_eda2_lambda_sorted = uvdist_array_eda2_lambda_no_zero[uvdist_array_eda2_inds]
            #print(uvdist_array_eda2_lambda_sorted[0:300])
      
            #uvdist_array_eda2_lambda_sorted_no_zero = uvdist_array_eda2_lambda_sorted[uvdist_array_eda2_lambda_sorted>0]
            #print(uvdist_array_eda2_lambda_sorted_no_zero[615:630])
            uvdist_array_eda2_lambda_sorted_cut = uvdist_array_eda2_lambda_sorted[uvdist_array_eda2_lambda_sorted < uvdist_thresh_lambda]
            #print(uvdist_array_eda2_lambda_sorted_cut[0:300])
            n_baselines_included_eda2 = len(uvdist_array_eda2_lambda_sorted_cut)
            print("number of baselines included eda2 %s" % n_baselines_included_eda2)
            
            #remove the autos:
            common_eda2_data_complex_x_real_sorted_no_auto = common_eda2_data_complex_x_real[uvdist_array_eda2_lambda>0]
            common_eda2_data_complex_x_real_sorted = common_eda2_data_complex_x_real_sorted_no_auto[uvdist_array_eda2_inds]
            common_eda2_data_complex_x_real_sorted_cut = common_eda2_data_complex_x_real_sorted[0:n_baselines_included_eda2]
            
      
            #unity and eda2 differ by two - I am probably not calculating the uvw accurately enough in sims
            #just use the eda2 data uvdist cutoff
            uvdist_array_unity_lambda_sorted_cut = uvdist_array_unity_lambda_sorted[0:n_baselines_included_eda2]
       
            common_unity_data_complex_x_real_sorted_cut = common_unity_data_complex_freq_av_x_real_sorted[0:n_baselines_included_eda2]
            
            print(uvdist_array_eda2_lambda_sorted_cut.shape)
            print(common_eda2_data_complex_x_real_sorted_cut.shape)
            #plot the expected unity response vs baseline length
            plt.clf()
            plt.scatter(uvdist_array_eda2_lambda_sorted_cut,common_eda2_data_complex_x_real_sorted_cut)            
            map_title="Unity_response_vs_uvdist" 
            plt.xlabel("UV distance (wavelengths)")
            plt.ylabel("Real component of visibility X pol (Jy)")
            plt.legend(loc=1)
            #plt.text(x_pos, y_pos, fit_string)
            #plt.ylim([0, 3.5])
            fig_name= "eda2_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            plt.close()
            print("saved %s" % fig_name)
      
            #plot the expected unity response vs baseline length
            plt.clf()
            plt.scatter(uvdist_array_unity_lambda_sorted_cut,common_unity_data_complex_x_real_sorted_cut)            
            map_title="Unity_response_vs_uvdist" 
            plt.xlabel("UV distance (wavelengths)")
            plt.ylabel("Real component of visibility X pol (Jy)")
            plt.legend(loc=1)
            #plt.text(x_pos, y_pos, fit_string)
            #plt.ylim([0, 3.5])
            fig_name= "unity_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            plt.close()
            print("saved %s" % fig_name)
      
            print(len(common_eda2_data_complex_x_real_sorted_cut))
            print(len(common_unity_data_complex_x_real_sorted_cut))
      
            if (len(common_eda2_data_complex_x_real_sorted_cut)>2 and len(common_unity_data_complex_x_real_sorted_cut)>2):
               model = sm.OLS(common_eda2_data_complex_x_real_sorted_cut, common_unity_data_complex_x_real_sorted_cut,missing='drop')
               results = model.fit()
               parameters = results.params
               #print parameters
               t_sky_K = parameters[0]
               t_sky_error_K = results.bse[0]    
         
               print("t_sky_K is %0.4E +/- %0.04f K" % (t_sky_K,t_sky_error_K))
               fit_string = "y=%0.1fx" % t_sky_K         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
               print(fit_string)
               
               #t_sky_K_list.append(t_sky_K)
               #t_sky_error_K_list.append(t_sky_error_K)
         
               plt.clf()
               plt.plot(common_unity_data_complex_x_real_sorted_cut, common_eda2_data_complex_x_real_sorted_cut,linestyle='None',marker='.')
               #deal with the missing values in the fit
               common_unity_data_complex_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_x_real_sorted_cut))[:,0]
               #print(common_unity_data_complex_x_real_sorted_cut_missing_indices.shape)
               common_unity_data_complex_x_real_sorted_cut_missing = common_unity_data_complex_x_real_sorted_cut[common_unity_data_complex_x_real_sorted_cut_missing_indices]
               plt.plot(common_unity_data_complex_x_real_sorted_cut_missing, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
               map_title="Data and fit" 
               plt.xlabel("Expected global-signal response")
               plt.ylabel("Real component of visibility X pol (Jy)")
               plt.legend(loc=1)
               #plt.text(x_pos, y_pos, fit_string)
               #plt.ylim([0, 3.5])
               fig_name= "x_y_OLS_%s_%s_pol.png" % (base_name,"x")
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               plt.close()
               print("saved %s" % fig_name)
            
            else:
               t_sky_K = np.nan
               t_sky_error_K = np.nan
            
            #save to arrays
            t_sky_K_array_fine_chan[EDA2_chan_index,fine_chan_included,timestep] = t_sky_K
            t_sky_error_K_array_fine_chan[EDA2_chan_index,fine_chan_included,timestep] = t_sky_error_K
            
            ###FLAGGING bit 
            max_deviations = 5
            common_eda2_data_complex_x_real_sorted_cut_std_dev = np.nanstd(common_eda2_data_complex_x_real_sorted_cut)
            max_distance_from_model = max_deviations * common_eda2_data_complex_x_real_sorted_cut_std_dev
            
            ##need to deal with missing values from fit
            common_eda2_data_complex_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_x_real_sorted_cut))[:,0]
            common_eda2_data_complex_x_real_sorted_cut_missing = common_eda2_data_complex_x_real_sorted_cut[common_eda2_data_complex_x_real_sorted_cut_missing_indices]
            distance_from_model = np.abs(common_eda2_data_complex_x_real_sorted_cut_missing - results.fittedvalues)
            common_eda2_data_complex_x_real_sorted_cut_flagged = np.copy(common_eda2_data_complex_x_real_sorted_cut_missing)
            common_eda2_data_complex_x_real_sorted_cut_flagged[distance_from_model > max_distance_from_model] = np.nan
            #print(common_eda2_data_complex_x_real_sorted_cut_flagged.shape)
            #print(common_unity_data_complex_x_real_sorted_cut_missing.shape)
      
            print(len(common_eda2_data_complex_x_real_sorted_cut_flagged))
            print(len(common_unity_data_complex_x_real_sorted_cut_missing))
            
            if (len(common_eda2_data_complex_x_real_sorted_cut_flagged)>2 and len(common_unity_data_complex_x_real_sorted_cut_missing)>2):
               model = sm.OLS(common_eda2_data_complex_x_real_sorted_cut_flagged, common_unity_data_complex_x_real_sorted_cut_missing,missing='drop')
               results = model.fit()
               parameters = results.params
               #print parameters
               t_sky_K_flagged = parameters[0]
               t_sky_error_K_flagged = results.bse[0]    
               
               print("t_sky_K flagged is %0.4E +/- %0.04f K" % (t_sky_K_flagged,t_sky_error_K_flagged))
               fit_string = "y=%0.1fx" % t_sky_K_flagged        #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
               print(fit_string)
               
               plt.clf()
               plt.plot(common_unity_data_complex_x_real_sorted_cut_missing, common_eda2_data_complex_x_real_sorted_cut_flagged,linestyle='None',marker='.')
               #deal with the missing values in the fit
               common_unity_data_complex_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_x_real_sorted_cut_flagged))[:,0]
               #print(common_unity_data_complex_x_real_sorted_cut_missing_indices.shape)
               common_unity_data_complex_x_real_sorted_cut_missing = common_unity_data_complex_x_real_sorted_cut_missing[common_unity_data_complex_x_real_sorted_cut_missing_indices]
               
               
               plt.plot(common_unity_data_complex_x_real_sorted_cut_missing, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
               map_title="Data and fit flagged" 
               plt.xlabel("Expected global-signal response")
               plt.ylabel("Real component of visibility X pol (Jy)")
               plt.legend(loc=1)
               #plt.text(x_pos, y_pos, fit_string)
               #plt.ylim([0, 3.5])
               fig_name= "x_y_OLS_%s_%s_pol_flagged.png" % (base_name,"x")
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               plt.close()
               print("saved %s" % fig_name)
            else:
               t_sky_K_flagged = np.nan
               t_sky_error_K_flagged = np.nan  
               
            t_sky_K_array_flagged_fine_chan[EDA2_chan_index,fine_chan_included,timestep] = t_sky_K_flagged
            t_sky_error_K_array_flagged_fine_chan[EDA2_chan_index,fine_chan_included,timestep] = t_sky_error_K_flagged
      
   
      ###

      #average the eda2 timesteps - 
      sum_data_eda2 = np.zeros([n_baselines_eda2,32,4],dtype='complex')
      n_unflagged_obs_sum = np.zeros([n_baselines_eda2,32,4],dtype='complex')
      start_index = 0
      end_index = n_baselines_eda2
      for timestep in range(n_timesteps_eda2):
         timestep_data = eda2_data_complex[start_index:end_index,:,:]
         unflag_data_boolean = ~eda2_flags[start_index:end_index,:,:]
         unflag_data_int = unflag_data_boolean.astype(int)
         n_unflagged_obs_sum += unflag_data_int.astype(complex)
         sum_data_eda2 += timestep_data
         start_index += n_baselines_eda2
         end_index += n_baselines_eda2

      #av_eda2_data = sum_data_eda2 / n_timesteps_eda2
      n_unflagged_obs_sum_zero_indices = np.where(n_unflagged_obs_sum==0+0j)
      av_eda2_data = sum_data_eda2 / n_unflagged_obs_sum
      eda2_data_complex = av_eda2_data


      common_eda2_data = eda2_data_complex[eda2_ms_indices]
      common_eda2_data_sorted = common_eda2_data[eda2_common_sort_inds]
      #print(common_eda2_data_sorted[0:10])
      #print(common_eda2_data_sorted.shape)
      #print(common_eda2_data_sorted[0:10,0,0])
 

      #in calibrate_with_complex_beams I had to go through a convoluted process to add in ant 255 to the model (unity in this case) and to add in the autos...
      #but I think this is fine now? maybe because I am adding autos in the sims, not sure whats going on with 255
      
      #freq av and only use real for now, be sure to use nanmean so you dont make everything a nan.
      start_fine_chan = 2
      end_fine_chan = n_fine_chans_included+2
      common_eda2_data_complex_freq_av = np.nanmean(common_eda2_data_sorted[:,start_fine_chan:end_fine_chan,:],axis=1)
      common_eda2_data_complex_freq_av_x = common_eda2_data_complex_freq_av[:,0]
      common_eda2_data_complex_freq_av_x_real=common_eda2_data_complex_freq_av_x.real
      #there will be nans in here, but not in the unity - what happens? I think fine, they get droppen later an in the fitting
      #common_eda2_data_complex_freq_av_x_real_no_nans = common_eda2_data_complex_freq_av_x_real[~np.isnan(common_eda2_data_complex_freq_av_x_real)]


      #Need to sort by baseline length (then only use short baselines)
      uvdist_array_eda2_m = np.sqrt(common_eda2_uvw_sorted[:,0]**2 + common_eda2_uvw_sorted[:,1]**2 + common_eda2_uvw_sorted[:,2]**2)
      uvdist_array_eda2_lambda = uvdist_array_eda2_m / wavelength
      #we don't want the auto included here
      uvdist_array_eda2_lambda_no_zero = uvdist_array_eda2_lambda[uvdist_array_eda2_lambda>0.]
      uvdist_array_eda2_inds = uvdist_array_eda2_lambda_no_zero.argsort()
      uvdist_array_eda2_lambda_sorted = uvdist_array_eda2_lambda_no_zero[uvdist_array_eda2_inds]
      #print(uvdist_array_eda2_lambda_sorted[0:300])

      #uvdist_array_eda2_lambda_sorted_no_zero = uvdist_array_eda2_lambda_sorted[uvdist_array_eda2_lambda_sorted>0]
      #print(uvdist_array_eda2_lambda_sorted_no_zero[615:630])
      uvdist_array_eda2_lambda_sorted_cut = uvdist_array_eda2_lambda_sorted[uvdist_array_eda2_lambda_sorted < uvdist_thresh_lambda]
      #print(uvdist_array_eda2_lambda_sorted_cut[0:300])
      n_baselines_included_eda2 = len(uvdist_array_eda2_lambda_sorted_cut)
      print("number of baselines included eda2 %s" % n_baselines_included_eda2)
      
      #remove the autos:
      common_eda2_data_complex_freq_av_x_real_sorted_no_auto = common_eda2_data_complex_freq_av_x_real[uvdist_array_eda2_lambda>0]
      common_eda2_data_complex_freq_av_x_real_sorted = common_eda2_data_complex_freq_av_x_real_sorted_no_auto[uvdist_array_eda2_inds]
      common_eda2_data_complex_freq_av_x_real_sorted_cut = common_eda2_data_complex_freq_av_x_real_sorted[0:n_baselines_included_eda2]
      

      #unity and eda2 differ by two - I am probably not calculating the uvw accurately enough in sims
      #just use the eda2 data uvdist cutoff
      uvdist_array_unity_lambda_sorted_cut = uvdist_array_unity_lambda_sorted[0:n_baselines_included_eda2]
 
      common_unity_data_complex_freq_av_x_real_sorted_cut = common_unity_data_complex_freq_av_x_real_sorted[0:n_baselines_included_eda2]
      

      #plot the expected unity response vs baseline length
      plt.clf()
      plt.scatter(uvdist_array_eda2_lambda_sorted_cut,common_eda2_data_complex_freq_av_x_real_sorted_cut)            
      map_title="Unity_response_vs_uvdist" 
      plt.xlabel("UV distance (wavelengths)")
      plt.ylabel("Real component of visibility X pol (Jy)")
      plt.legend(loc=1)
      #plt.text(x_pos, y_pos, fit_string)
      #plt.ylim([0, 3.5])
      fig_name= "eda2_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name)

      #plot the expected unity response vs baseline length
      plt.clf()
      plt.scatter(uvdist_array_unity_lambda_sorted_cut,common_unity_data_complex_freq_av_x_real_sorted_cut)            
      map_title="Unity_response_vs_uvdist" 
      plt.xlabel("UV distance (wavelengths)")
      plt.ylabel("Real component of visibility X pol (Jy)")
      plt.legend(loc=1)
      #plt.text(x_pos, y_pos, fit_string)
      #plt.ylim([0, 3.5])
      fig_name= "unity_response_vs_uvdist_%0.3f_MHz_%s_pol.png" % (freq_MHz,"x")
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name)

      if (len(common_eda2_data_complex_freq_av_x_real_sorted_cut)>2 and len(common_unity_data_complex_freq_av_x_real_sorted_cut)>2):
         model = sm.OLS(common_eda2_data_complex_freq_av_x_real_sorted_cut, common_unity_data_complex_freq_av_x_real_sorted_cut,missing='drop')
         results = model.fit()
         parameters = results.params
         #print parameters
         #t_sky_jy = parameters[0]
         #t_sky_error_jy = results.bse[0]
         #t_sky_K = jy_to_K * t_sky_jy
         #t_sky_error_K = jy_to_K * t_sky_error_jy
         t_sky_K = parameters[0]
         t_sky_error_K = results.bse[0]    
   
         print("t_sky_K is %0.4E +/- %0.04f K" % (t_sky_K,t_sky_error_K))
         fit_string = "y=%0.1fx" % t_sky_K         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
         print(fit_string)
         
         #t_sky_K_list.append(t_sky_K)
         #t_sky_error_K_list.append(t_sky_error_K)
   
         plt.clf()
         plt.plot(common_unity_data_complex_freq_av_x_real_sorted_cut, common_eda2_data_complex_freq_av_x_real_sorted_cut,linestyle='None',marker='.')
         #deal with the missing values in the fit
         common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_freq_av_x_real_sorted_cut))[:,0]
         #print(common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices.shape)
         common_unity_data_complex_freq_av_x_real_sorted_cut_missing = common_unity_data_complex_freq_av_x_real_sorted_cut[common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices]
         plt.plot(common_unity_data_complex_freq_av_x_real_sorted_cut_missing, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
         map_title="Data and fit" 
         plt.xlabel("Expected global-signal response")
         plt.ylabel("Real component of visibility X pol (Jy)")
         plt.legend(loc=1)
         #plt.text(x_pos, y_pos, fit_string)
         #plt.ylim([0, 3.5])
         fig_name= "x_y_OLS_plot_%0.3f_MHz_all_%s_pol.png" % (freq_MHz,"x")
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         plt.close()
         print("saved %s" % fig_name)
      else:
         t_sky_K = np.nan
         t_sky_error_K = np.nan
               
      t_sky_K_array[EDA2_chan_index] = t_sky_K
      t_sky_error_K_array[EDA2_chan_index] = t_sky_error_K
      
      ##FLAGGING bit - shouldn't need this anymore now that rfi flagging is being done properly...
      ##now use the fit to identify outliers probably due to rfi
      max_deviations = 5
      common_eda2_data_complex_freq_av_x_real_sorted_cut_std_dev = np.nanstd(common_eda2_data_complex_freq_av_x_real_sorted_cut)
      max_distance_from_model = max_deviations * common_eda2_data_complex_freq_av_x_real_sorted_cut_std_dev
      
      #need to deal with missing values from fit
      common_eda2_data_complex_freq_av_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_freq_av_x_real_sorted_cut))[:,0]
      common_eda2_data_complex_freq_av_x_real_sorted_cut_missing = common_eda2_data_complex_freq_av_x_real_sorted_cut[common_eda2_data_complex_freq_av_x_real_sorted_cut_missing_indices]
      distance_from_model = np.abs(common_eda2_data_complex_freq_av_x_real_sorted_cut_missing - results.fittedvalues)
      common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged = np.copy(common_eda2_data_complex_freq_av_x_real_sorted_cut_missing)
      common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged[distance_from_model > max_distance_from_model] = np.nan
      print(common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged.shape)
      print(common_unity_data_complex_freq_av_x_real_sorted_cut_missing.shape)

      #common_unity_data_complex_freq_av_x_real_sorted_cut_flagged = common_unity_data_complex_freq_av_x_real_sorted_cut_missing[[distance_from_model < max_distance_from_model]]
      if (len(common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged)>0 and len(common_unity_data_complex_freq_av_x_real_sorted_cut_missing)>0):
         model = sm.OLS(common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged, common_unity_data_complex_freq_av_x_real_sorted_cut_missing,missing='drop')
         results = model.fit()
         parameters = results.params
         #print parameters
         #t_sky_jy = parameters[0]
         #t_sky_error_jy = results.bse[0]
         #t_sky_K = jy_to_K * t_sky_jy
         #t_sky_error_K = jy_to_K * t_sky_error_jy
         t_sky_K_flagged = parameters[0]
         t_sky_error_K_flagged = results.bse[0]    
         
         print("t_sky_K flagged is %0.4E +/- %0.04f K" % (t_sky_K_flagged,t_sky_error_K_flagged))
         fit_string = "y=%0.1fx" % t_sky_K_flagged        #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
         print(fit_string)
         
         plt.clf()
         plt.plot(common_unity_data_complex_freq_av_x_real_sorted_cut_missing, common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged,linestyle='None',marker='.')
         #deal with the missing values in the fit
         common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices = np.argwhere(~np.isnan(common_eda2_data_complex_freq_av_x_real_sorted_cut_flagged))[:,0]
         #print(common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices.shape)
         common_unity_data_complex_freq_av_x_real_sorted_cut_missing = common_unity_data_complex_freq_av_x_real_sorted_cut_missing[common_unity_data_complex_freq_av_x_real_sorted_cut_missing_indices]
         
         
         plt.plot(common_unity_data_complex_freq_av_x_real_sorted_cut_missing, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
         map_title="Data and fit flagged" 
         plt.xlabel("Expected global-signal response")
         plt.ylabel("Real component of visibility X pol (Jy)")
         plt.legend(loc=1)
         #plt.text(x_pos, y_pos, fit_string)
         #plt.ylim([0, 3.5])
         fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol_%s_flagged.png" % (freq_MHz,"x",EDA2_obs_time)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         plt.close()
         print("saved %s" % fig_name)
      else:
         t_sky_K_flagged = np.nan
         t_sky_error_K_flagged = np.nan 
         
      t_sky_K_array_flagged[EDA2_chan_index] = t_sky_K_flagged
      t_sky_error_K_array_flagged[EDA2_chan_index] = t_sky_error_K_flagged

   np.save(t_sky_K_array_filename,t_sky_K_array)     
   np.save(t_sky_error_K_array_filename,t_sky_error_K_array)  
   np.save(t_sky_K_array_filename_flagged,t_sky_K_array_flagged)     
   np.save(t_sky_error_K_array_filename_flagged,t_sky_error_K_array_flagged) 
   print("saved %s" % t_sky_K_array_filename)
   print("saved %s" % t_sky_error_K_array_filename)
   print("saved %s" % t_sky_K_array_filename_flagged)
   print("saved %s" % t_sky_error_K_array_filename_flagged)
   
   np.save(t_sky_K_array_filename_fine_chan,t_sky_K_array_fine_chan)     
   np.save(t_sky_error_K_array_filename_fine_chan,t_sky_error_K_array_fine_chan)  
   np.save(t_sky_K_array_filename_flagged_fine_chan,t_sky_K_array_flagged_fine_chan)     
   np.save(t_sky_error_K_array_filename_flagged_fine_chan,t_sky_error_K_array_flagged_fine_chan) 
   print("saved %s" % t_sky_K_array_filename_fine_chan)
   print("saved %s" % t_sky_error_K_array_filename_fine_chan)
   print("saved %s" % t_sky_K_array_filename_flagged_fine_chan)
   print("saved %s" % t_sky_error_K_array_filename_flagged_fine_chan)
      
   #return(t_sky_K_array,t_sky_error_K_array)

def plot_t_sky_and_fit_foregrounds(freq_MHz_list,t_sky_K_array_filename,t_sky_error_K_array_filename,poly_order=5):
   base_name = t_sky_K_array_filename.split(".npy")[0]
   t_sky_K_array = np.load(t_sky_K_array_filename)
   t_sky_error_K_array = np.load(t_sky_error_K_array_filename)
   
   #print(freq_MHz_list)
   #print(t_sky_K_array)
   #print(t_sky_error_K_array)
   plt.clf()
   plt.errorbar(freq_MHz_list,t_sky_K_array,yerr=t_sky_error_K_array)
   map_title="Global Tsky" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("EDA2 global Tsky (K)")
   #plt.legend(loc=1)
   #plt.text(x_pos, y_pos, fit_string)
   #plt.ylim([0, 3.5])
   fig_name= "%s.png"  % base_name
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name)
      
   #log fit and subtract
   log_sky_array = np.log10(t_sky_K_array)
   log_freq_MHz_array = np.log10(np.asarray(freq_MHz_list))
   coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
   ffit = poly.polyval(log_freq_MHz_array, coefs)
   ffit_linear = 10**ffit
   
   #log_residual = log_signal_array_short_baselines - log_ffit
   residual_of_log_fit = ffit_linear - t_sky_K_array
   
   rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
   print("rms_of_residuals is %0.3f K" % rms_of_residuals)
   
   max_abs_residuals = np.max(np.abs(residual_of_log_fit))
   y_max = 1.5 * max_abs_residuals
   y_min = 1.5 * -max_abs_residuals
   
   #print(freq_array_cut)
   #print(residual_of_log_fit)
   
   plt.plot(freq_MHz_list,residual_of_log_fit,label="residual")
   plt.text(50, max_abs_residuals, "rms=%1.2f K" % rms_of_residuals)
   
    
   map_title="Residual for log polynomial order %s fit " % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("Frequency (MHz)")
   #if len(model_type_list)>1:
   plt.legend(loc=1)
   #plt.legend(loc=1)
   plt.ylim([y_min, y_max])
   fig_name= "%s_log_fit_residual_poly_order_%s.png" % (base_name,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   plt.close()    

def plot_t_sky_waterfalls_timestep_finechan(EDA2_chan_list,t_sky_K_array_filename_fine_chan,t_sky_error_K_array_filename_fine_chan):
   freq_MHz_array = np.asarray(EDA2_chan_list)*400./512.
   base_name = t_sky_K_array_filename_fine_chan.split(".npy")[0]
   t_sky_K_array_timestep_finechan = np.load(t_sky_K_array_filename_fine_chan)
   t_sky_error_K_array_timestep_finechan = np.load(t_sky_error_K_array_filename_fine_chan)
   print(t_sky_K_array_timestep_finechan.shape)
   print(t_sky_error_K_array_timestep_finechan.shape)
   print(freq_MHz_array.shape)
   #make waterfall plots of T_skay and Tsky_error for each coarse chan
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):
      t_sky_K_array_EDA2_chan = t_sky_K_array_timestep_finechan[EDA2_chan_index,:,:]
      t_sky_error_K_array_EDA2_chan = t_sky_error_K_array_timestep_finechan[EDA2_chan_index,:,:]
      plt.clf()
      map_title="Tsky waterfall EDA2 chan %s" % EDA2_chan
      plt.imshow(t_sky_K_array_EDA2_chan)
      plt.title(map_title)
      bar=plt.colorbar()
      plt.xlabel('timestep')
      plt.ylabel('fine chan')
      bar.set_label('Tsky K')
      fig_name= "%s_waterfall_EDA2_chan_%03d.png" % (base_name,EDA2_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
   
      plt.clf()
      map_title="Tsky error waterfall EDA2 chan %s" % EDA2_chan
      plt.imshow(t_sky_error_K_array_EDA2_chan)
      plt.title(map_title)
      bar=plt.colorbar()
      plt.xlabel('timestep')
      plt.ylabel('finechan')
      bar.set_label('Tsky error K')
      fig_name= "%s_error_waterfall_EDA2_chan_%03d.png" % (base_name,EDA2_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)

def manual_flag_and_average_tsky(EDA2_chan_list,t_sky_K_array_filename_fine_chan,t_sky_error_K_array_filename_fine_chan):
   #just visually inspect each waterfall and manually flag bad fine chans and time steps here
   freq_MHz_array = np.asarray(EDA2_chan_list)*400./512.
   base_name = t_sky_K_array_filename_fine_chan.split(".npy")[0]
   t_sky_K_array_manual_flagged_averaged_filename = "t_sky_K_array_manual_flagged_averaged.npy"
   t_sky_K_array_manual_flagged_averaged_error_filename = "t_sky_K_array_manual_flagged_averaged_error.npy"
   t_sky_K_array_manual_flagged_averaged = np.empty(len(EDA2_chan_list))
   t_sky_K_array_manual_flagged_averaged[:] = np.nan
   t_sky_K_array_manual_flagged_averaged_error = np.empty(len(EDA2_chan_list))
   t_sky_K_array_manual_flagged_averaged_error[:] = np.nan
      
   t_sky_K_array_timestep_finechan = np.load(t_sky_K_array_filename_fine_chan)
   t_sky_error_K_array_timestep_finechan = np.load(t_sky_error_K_array_filename_fine_chan)
   print(t_sky_K_array_timestep_finechan.shape)
   print(t_sky_error_K_array_timestep_finechan.shape)
   print(freq_MHz_array.shape)
   #flag then make new waterfall plots of T_skay and Tsky_error for each coarse chan
   #indexing goes: EDA2_chan,fine_chan,timestep
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):
      #manual flagging here
      if EDA2_chan==71:
         bad_timestep_list = [0]
         bad_fine_chan_list = [13,14,18,19,20,21,22,23,24,25]
         for bad_timestep in bad_timestep_list:
            t_sky_K_array_timestep_finechan[EDA2_chan_index,:,bad_timestep] = np.nan
            t_sky_error_K_array_timestep_finechan[EDA2_chan_index,:,bad_timestep] = np.nan
         for bad_fine_chan in bad_fine_chan_list:
            t_sky_K_array_timestep_finechan[EDA2_chan_index,bad_fine_chan,:] = np.nan
            t_sky_error_K_array_timestep_finechan[EDA2_chan_index,bad_fine_chan,:] = np.nan         
      
      plt.clf()
      map_title="Tsky waterfall EDA2 chan %s" % EDA2_chan
      plt.imshow(t_sky_K_array_timestep_finechan[EDA2_chan_index,:,:])
      plt.title(map_title)
      bar=plt.colorbar()
      plt.xlabel('timestep')
      plt.ylabel('fine chan')
      bar.set_label('Tsky K')
      fig_name= "%s_waterfall_EDA2_chan_%03d_manual_flag.png" % (base_name,EDA2_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
   
      plt.clf()
      map_title="Tsky error waterfall EDA2 chan %s" % EDA2_chan
      plt.imshow(t_sky_error_K_array_timestep_finechan[EDA2_chan_index,:,:])
      plt.title(map_title)
      bar=plt.colorbar()
      plt.xlabel('timestep')
      plt.ylabel('finechan')
      bar.set_label('Tsky error K')
      fig_name= "%s_error_waterfall_EDA2_chan_%03d_manual_flag.png" % (base_name,EDA2_chan)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      
      t_sky_K_array_manual_flagged_averaged[EDA2_chan_index] = np.nanmean(t_sky_K_array_timestep_finechan[EDA2_chan_index,:,:])
      t_sky_K_array_manual_flagged_averaged_error[EDA2_chan_index] = np.nanstd(t_sky_K_array_timestep_finechan[EDA2_chan_index,:,:])
   
   np.save(t_sky_K_array_manual_flagged_averaged_filename,t_sky_K_array_manual_flagged_averaged)
   np.save(t_sky_K_array_manual_flagged_averaged_error_filename,t_sky_K_array_manual_flagged_averaged_error)
   print("saved %s" % t_sky_K_array_manual_flagged_averaged_filename)
   print("saved %s" % t_sky_K_array_manual_flagged_averaged_error_filename)
      
      
#times
EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan(EDA2_data_dir,EDA2_chan_list)
EDA2_obs_time_list_each_chan = EDA2_obs_time_list_each_chan[0:]
#n_obs_concat_list = [len(obs_list) for obs_list in EDA2_obs_time_list_each_chan] 
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
      #print(EDA2_obs_time)
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time)
      lst_hrs_list.append(lst_eda2_hrs)
   else:
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time_list[0])
      lst_hrs_list.append(lst_eda2_hrs)
      
#print("lst_hrs_list")
#print(lst_hrs_list)
       
      
##unity only sim takes 2 min with nside 32, 6 mins with nside 64, similar 
#chan_num = 0
freq_MHz_list = freq_MHz_list[0:10]
lst_hrs_list = lst_hrs_list[0:10]
EDA2_obs_time_list = EDA2_obs_time_list[0:10]
EDA2_chan_list = EDA2_chan_list[0:10]
plot_from_saved = False
sim_unity=True
sim_pt_source=False
check_figs=False
#simulate_eda2_with_complex_beams(freq_MHz_list,lst_hrs_list,nside=32,plot_from_saved=plot_from_saved,EDA2_obs_time_list=EDA2_obs_time_list,sim_unity=sim_unity,sim_pt_source=sim_pt_source,check_figs=check_figs)
#input_sky = "gsm"   #zenith_point_source  #unity
#write_to_miriad_vis(freq_MHz_list,lst_hrs_list,EDA2_obs_time_list=EDA2_obs_time_list,input_sky=input_sky,check_figs=check_figs)
#input_sky = "unity"
#write_to_miriad_vis(freq_MHz_list,lst_hrs_list,EDA2_obs_time_list=EDA2_obs_time_list,input_sky=input_sky,check_figs=check_figs)
#######
#just for initial testing of calibration
#model_ms_name= "20200303T133733_50.000.ms"
#eda2_ms_name = "/md0/EoR/EDA2/20200303_data/64/20200303_133733_eda2_ch32_ant256_midday_avg8140.ms" 
#calibrate_with_complex_beam_model(model_ms_name=model_ms_name,eda2_ms_name=eda2_ms_name)
#######
#now with time averaging
plot_cal = True
per_chan_cal = False
run_aoflagger=True
#calibrate_with_complex_beam_model_time_av(EDA2_chan_list=EDA2_chan_list,lst_list=lst_hrs_list,plot_cal=plot_cal,uv_cutoff=0,per_chan_cal=per_chan_cal,run_aoflagger=run_aoflagger)
#sys.exit()
#now need to extract the global signal using the complex beams
#extract_global_signal_from_ms_complex(EDA2_chan_list=EDA2_chan_list,lst_list=lst_hrs_list)
#t_sky_K_array_filename = "t_sky_K_array_eda2.npy"
#t_sky_error_K_array_filename = "t_sky_error_K_array_eda2.npy"
#plot_t_sky_and_fit_foregrounds(freq_MHz_list,t_sky_K_array_filename,t_sky_error_K_array_filename)
#t_sky_K_array_filename = "t_sky_K_array_eda2_flagged.npy"
#t_sky_error_K_array_filename = "t_sky_error_K_array_eda2_flagged.npy"
#plot_t_sky_and_fit_foregrounds(freq_MHz_list,t_sky_K_array_filename,t_sky_error_K_array_filename)
#t_sky_K_array_filename_fine_chan = "t_sky_K_array_eda2_fine_chan.npy"
#t_sky_error_K_array_filename_fine_chan = "t_sky_error_K_array_eda2_fine_chan.npy"
#plot_t_sky_waterfalls_timestep_finechan(EDA2_chan_list,t_sky_K_array_filename_fine_chan,t_sky_error_K_array_filename_fine_chan)
t_sky_K_array_filename_fine_chan = "t_sky_K_array_eda2_flagged_fine_chan.npy"
t_sky_error_K_array_filename_fine_chan = "t_sky_error_K_array_eda2_flagged_fine_chan.npy"
#plot_t_sky_waterfalls_timestep_finechan(EDA2_chan_list,t_sky_K_array_filename_fine_chan,t_sky_error_K_array_filename_fine_chan)
manual_flag_and_average_tsky(EDA2_chan_list,t_sky_K_array_filename_fine_chan,t_sky_error_K_array_filename_fine_chan)
t_sky_K_array_filename = "t_sky_K_array_manual_flagged_averaged.npy"
t_sky_error_K_array_filename = "t_sky_K_array_manual_flagged_averaged_error.npy"
plot_t_sky_and_fit_foregrounds(freq_MHz_list,t_sky_K_array_filename,t_sky_error_K_array_filename)      
      
      