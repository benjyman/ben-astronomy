#!/usr/bin/env python
#ASSASSIN: (All Sky SignAl Short Spacing INterferometer)
#Script to replicate Cath's global step simulations
#
#Daniel:
# The format of the data is in standard spherical coordinates (theta, phi).
# Where theta is the zenith angle and phi is anti-clockwise going from east to north looking down at the array. 
# The 'origin' will be on the top left corner. Going from top to bottom is increasing in phi and going left to right is increasing in theta. 
# The data is in 1 deg steps.
#
import matplotlib
matplotlib.use('Agg')
import os,sys
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
import healpy as hp
from pygsm import GSMObserver
from pygsm import GlobalSkyModel
from pygsm import GlobalSkyModel2016
from datetime import datetime, date

from reproject import reproject_from_healpix
import pyfits
from astropy.wcs import WCS
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from scipy import signal
import numpy.polynomial.polynomial as poly
from pyuvdata import UVData

import random
from astropy import units as u
from astropy.coordinates import SkyCoord
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
from sklearn.linear_model import LinearRegression
import seaborn as sns
from scipy.optimize import curve_fit

from sklearn.preprocessing import normalize
#avoid astropy time error
from astropy.utils.iers import conf
conf.auto_max_age = None
from astropy.utils import iers
iers.conf.auto_download = False  
#from astroplan import download_IERS_A
#download_IERS_A()


  
sun_flux_density = 50000.0   #borkowski et al 1982?
#mu_0 = 4.*np.pi*10**(-7)
c = 3.0e8
k = 1.38065e-23
sq_deg_in_1_sr = (180./math.pi)**2
#Need to specify the location of the observatory (MWA ... put in exact EDA later)
# Setup observatory location - in this case, MWA, Australia
#latitude_degrees=-26.70331940
#longitude_degrees=116.67081524
#elevation_m=377.83

mwa_latitude_pyephem = "-26:42.199"
mwa_longitude_pyephem = "116:40.2488"
mwa_elevation = 0

mwa_latitude_deg = -26.70331940

mwa_latitude_astropy = '-26.7d'
mwa_longitude_astropy = '116.67d'

#for uvgen systemp values, sky temp at 180 MHz, beta is spectral index (Frank's 180 at 180), efficiency eta
T_180 = 180.0
beta = -2.5
eta = 1


#pointing dec for uvgen (needs to be zenith at each lst)
pointing_dec = "-26.7"

dipole_height_m = 0.3

#global signal params
#From Cath feasibility study (look at edges papers?)
#S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
C = 0.08/1000.
A = -150./1000.
delta_nu = 10.
nu_c = 70.

#These should normally be true
do_cleanup_images_and_vis = True
use_analytic_beam = True
generate_new_hpx = True
generate_new_vis = True
do_image_and_subtr_from_simulated_data = False

generate_new_assassin_sky = True
plot_gsm_map_hpx = False
plot_global_signal_map_hpx = False

#This is for Daniels FEKO model beams
generate_new_average_beam = False
#apply_normalisation = True
plot_all_beams = False
plot_average_beam = False

#array_ant_locations = 'other'
array_ant_locations = 'eda2'
#array_ant_locations = 'eda2_sub48'
#array_ant_locations = 'eda2_sim_compact_seed_123'
#array_ant_locations = 'eda2_sim_compact_seed_124'

#EDA2
n_fine_chans = 32
n_chans_omitted_each_edge = 1
      
#eda2 parameters
#instantaneous bw and channel width in MHz
inst_bw = 30.
#this has to be 1 MHz:
chan_width = 1.0
#fine chan stuff
centre_chan_index = 16
fine_chan_width_Hz = 28935 #Hz
                     
if array_ant_locations == 'eda2':
   array_ant_locations_filename = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant'
   array_label = 'eda_model'
elif array_ant_locations == 'eda2_sub48':
   array_ant_locations_filename = '/md0/code/git/ben-astronomy/EoR/ASSASSIN/EDA2_sub_array_48.txt'
   array_label = 'eda_sub48_model'
elif array_ant_locations == 'eda2_sim_compact_seed_123':
   array_ant_locations_filename = 'sim_array_seed_123_diameter_17_m_n_ants_256_min_spacing_75cm.txt'
   array_label = 'eda2_sim_compact_seed_123'
elif array_ant_locations == 'eda2_sim_compact_seed_124':
   array_ant_locations_filename = 'sim_array_seed_124_diameter_17_m_n_ants_256_min_spacing_75cm.txt'
   array_label = 'eda2_sim_compact_seed_124'
else:
   array_ant_locations_filename = ''
   array_label = 'array'

#beam_image_dir = "/md0/EoR/ASSASSIN/beam_fits/"
beam_image_dir = "/md0/EoR/ASSASSIN/beam_fits/fine_chan/"

#eda2_data_dir = '/md0/EoR/ASSASSIN/data_eda2/eda2_sub48/'
#EDA2_data_dir = '/md0/EoR/EDA2/20191213_data/'
#eda2_data_uvfits_name_list = ['chan_204_20190611T024741.uvfits']
#extract_from_eda2_data_outbase_name = 'eda2_sub48_data_'

start_lst_hrs = 2.0
#int_time_hrs = 1.0
#8 mins
int_time_hrs = 0.1333333

end_lst_hrs = start_lst_hrs + int_time_hrs
lst_list_array = np.arange(start_lst_hrs,end_lst_hrs,0.066667)
#1_hr_lst_0:
#lst_list_array = np.arange(0,2,0.066667)
lst_list = lst_list_array.astype(str)

#lst_list = ['0.120']
#lst_list = ['12']
#pol_list = ['X','Y']
sky_model = 'gsm'
#sky_model = 'gsm2016'
#sky_model = 'gmoss'

if sky_model=='gsm':
   NSIDE=512 #(to match pygsm)
      
      
pol_list = ['X']
#pol_list = ['Y']
#can be any of these, except if can only have 'diffuse' if not diffuse_global or diffuse_angular
#signal_type_list=['global','global_EDGES','diffuse','noise','gain_errors','diffuse_global','diffuse_angular']
signal_type_list=['diffuse','noise']
#signal_type_list=['diffuse']
#signal_type_list=['global_unity']
#signal_type_list=['diffuse_global','noise']
#signal_type_list=['global_EDGES']
#gsm_smooth_poly_order = 5
#can be 5,6,or 7 for joint fitting
poly_order = 7
#freq_MHz_list = np.arange(50,200,1)
start_chan=50
n_chan=150
chan_step = chan_width
freq_MHz_list = np.arange(start_chan,start_chan+n_chan,chan_step)
#freq_MHz_list = [160.]
#n_chan = len(freq_MHz_list)
n_pol = len(pol_list)
#harange definition for each 'snapshot' where we assume same sky good_freq_MHz_array (pb * sky). (start,stop,step in hrs)
#beam models and gsm maps are ~1 deg resolution, so sky takes 4 mins to move 1 deg, so generate new sky model * beam every 4 mins
#4 min (.066667 deg) obs with 2 sec  (.0005556 deg) steps (can average to 30 res sec later on)
#harange_string = "0,0.0666667,0.0005556"
#actually, just generate 1 timestep every 4 mins since this is how often you change the sky model. also, you don't want to phase track for multiple times, you always want the data phased to zenith
#pyuvdata won't accept 1 timestep!
harange_string = "0,0.0667,0.03335"


template_imsize = 512
template_cell_size_asec = 180./template_imsize*60.*60.
#amount to zero pad beam image for FFT
#want this to be big so we have good resolution in Fourier space (needed for close to zeros spacings)
padding_factor = 25
uvdist_wavelength_cutoff = 1.0
min_imaging_uvdist_wavelength = 0.5
max_imaging_uvdist_wavelength = 35.0
zero_spacing_leakage_threshold = 0.5

outbase_name = 'lst_%0.2f_hr_int_%0.2f_hr' % (start_lst_hrs,int_time_hrs)


#from Wayth et al Table 2, third column divided by 256 (not last col which is dipole in isolation)
Aeff_list = np.array([970,950,914,874,832,771,707,638,568,498,435,377,329,288,252,222,196])/256.
Aeff_freqs = np.arange(60,230,10)

#fit polynomial 
z = np.polyfit(Aeff_freqs, Aeff_list, 3)
p = np.poly1d(z)
Aeff_for_freq_MHz_list = p(freq_MHz_list)

def get_eda2_lst(eda_time_string="20151202T171727"):
   #Hack to get LST! use old time (2015) astropy.utils.iers.iers.IERSRangeError: (some) times are outside of range covered by IERS table.
   #eda_time_string = "20191202T171727"
   year, month, day, hour, minute, second = eda_time_string[0:4], eda_time_string[4:6],eda_time_string[6:8], eda_time_string[9:11],eda_time_string[11:13],eda_time_string[13:15]
   eda2_astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (float(year), float(month), float(day), float(hour), float(minute), float(second))
   print(eda2_astropy_time_string)
   eda2_astropy_time = Time(eda2_astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
   #calculate the local sidereal time for the MWA at eda2_astropy_time
   eda2_obs_lst = eda2_astropy_time.sidereal_time('apparent')
   eda2_obs_lst_hrs = eda2_obs_lst.value
   return eda2_obs_lst_hrs

def rotation_matrix(axis, theta):
   """
   Return the rotation matrix associated with counterclockwise rotation about
   the given axis by theta radians.
   https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
   """
   axis = np.asarray(axis)
   axis = axis / math.sqrt(np.dot(axis, axis))
   #axis_array = axis_array / np.linalg.norm(axis_array,axis=1)
   #axis_array = normalize(axis_array, axis=1, norm='l1')
   #print axis_array
   a = np.cos(theta / 2.0)
   #b_c_d = - axis_array * np.sin(theta / 2.0)
   #b = b_c_d[:,0]
   #c = b_c_d[:,1]
   #d = b_c_d[:,2]
   b, c, d = - axis * np.sin(theta / 2.0)

   aa, bb, cc, dd = a * a, b * b, c * c, d * d
   bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

   
   return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                    [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                    [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def rotate_vector(axis, theta, vector):

   return np.transpose(np.dot(rotation_matrix(axis, theta), np.transpose(vector)))
   #print(np.dot(rotation_matrix(axis, theta), v)) 
   
#v = [[1, 0, 0],[0, 1, 0],[0, 0, 1]]
#print v
#axis = [0, 1, 0]
#theta = np.pi/2.
#rotated_vector = rotate_vector(axis,theta,v)
#print np.round(rotated_vector)
##print rotated_vector.shape
#sys.exit()

def rotate_map(hmap, rot_theta, rot_phi):
    """
    https://stackoverflow.com/questions/24636372/apply-rotation-to-healpix-map-in-healpy
    Take hmap (a healpix map array) and return another healpix map array 
    which is ordered such that it has been rotated in (theta, phi) by the 
    amounts given.
    """
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Define a rotator
    r = hp.Rotator(deg=False, rot=[rot_phi,rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot)

    return rot_map

def make_hpx_beam_test(NSIDE,pol,wavelength,dipole_height_m):
   n_pix = hp.nside2npix(NSIDE)
   hpx_index_array = np.arange(0,n_pix,1)

   phi_array,theta_array = np.radians(hp.pix2ang(NSIDE,hpx_index_array,lonlat=True))
   

   short_dipole_parallel_beam_map = theta_array * 0.
   short_dipole_parallel_beam_map[theta_array < np.pi/2.] = 1.
   
   return short_dipole_parallel_beam_map    
    

def make_hpx_beam(NSIDE,pol,wavelength,dipole_height_m):
   n_pix = hp.nside2npix(NSIDE)
   hpx_index_array = np.arange(0,n_pix,1)

   theta_array,phi_array = hp.pix2ang(NSIDE,hpx_index_array)
   
   ##This is for YY dipole: !
   if (pol=='Y'):
      theta_parallel_array=np.arccos(np.sin(theta_array)*np.cos(phi_array))
   else:
   #This is for XX dipole!
      theta_parallel_array=np.arccos(np.sin(theta_array)*np.sin(phi_array)) 
   
   #theta_parallel_array=np.arccos(np.sin(theta_array)*np.cos(phi_array))
   
   d_in_lambda = (2. * dipole_height_m)/wavelength
   gp_effect_array = 2.*np.sin(np.pi*d_in_lambda*np.cos(theta_array))
   voltage_parallel_array=np.sin(theta_parallel_array) * gp_effect_array
   short_dipole_parallel_beam_map = voltage_parallel_array**2       #* factor
   
   #set to zero below horizon
   short_dipole_parallel_beam_map[theta_array > np.pi/2.]=0.
   
   #normalise to one at max
   beam_max = np.max(short_dipole_parallel_beam_map)
   short_dipole_parallel_beam_map = short_dipole_parallel_beam_map / beam_max
   
   return short_dipole_parallel_beam_map
      
      
def cleanup_images_and_vis(array_label,lst,freq_MHz,pol):
   print("###################doing cleanup########################")
   #print("lst %s hrs" % lst)
   lst_deg = (float(lst)/24)*360.
   #print("lst %s deg" % lst_deg)
   year=2000
   month=1
   day=1
   hour=np.floor(float(lst))
   minute=np.floor((float(lst)-hour) * 60.)
   second=((float(lst)-hour) * 60. - minute) * 60.
      
   date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
  
   gsm_hpx_fits_name = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
   reprojected_gsm_fitsname = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
   reprojected_gsm_im_name = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)
   
   global_signal_hpx_fits_name = "global_signal_map_LST_%03d_%0.3f_MHz_hpx.fits" % (lst_deg,freq_MHz)
   reprojected_global_signal_fitsname = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (lst_deg,freq_MHz)
   reprojected_global_signal_im_name = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.im" % (lst_deg,freq_MHz)
   
   #lna_impedance_aavs1_filename = "/data/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
   
   base_vis_name = "%s_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_image_name = "%s_no_src_LST_%03d_%0.3f_MHz.map" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_beam_name = "%s_no_src_LST_%03d_%0.3f_MHz.beam" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_image_fits_name = "%s_no_src_LST_%03d_%0.3f_MHz.fits" % (array_label,lst_deg,freq_MHz)
   
   reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
   reprojected_global_signal_im_Jy_per_pix_name =  "global_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (date_time_string,freq_MHz)
   
   cmd = "rm -rf %s %s %s %s %s %s %s %s %s %s %s" % (reprojected_gsm_fitsname,reprojected_gsm_im_name,global_signal_hpx_fits_name,reprojected_global_signal_fitsname,reprojected_global_signal_im_name,base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name,reprojected_gsm_im_Jy_per_pix_name,reprojected_global_signal_im_Jy_per_pix_name)
   print(cmd)
   os.system(cmd)
   
   #don't get rid of this - need it for calibrating eda2 data
   #apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,int(freq_MHz))
   apparent_sky_im_Tb_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz_Tb.im" % (lst_deg,pol,freq_MHz)
   #apparent_sky_im_Tb_fits_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz_Tb.fits" % (lst_deg,pol,int(freq_MHz))
   apparent_global_signal_im_name = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
   cmd = "rm -rf %s %s " % (apparent_global_signal_im_name,apparent_sky_im_Tb_name)
   print(cmd)
   os.system(cmd)
   
   model_sky_vis = "%s_plus_sky_LST_%03d_%s_pol_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   model_global_signal_vis = "%s_plus_G_LST_%03d_%s_pol_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   eda_model_noise_vis_name = "%s_N_LST_%03d_%s_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   cmd = "rm -rf %s %s %s" % (model_sky_vis,model_global_signal_vis,eda_model_noise_vis_name)
   print(cmd)
   os.system(cmd)
         
   beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
   beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz)
   cmd = "rm -rf %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name)
   print(cmd)
   os.system(cmd)        

   one_jy_source_vis_name =  "%s_one_jy_source_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   one_jy_source_vis_name_assassin =  "%s_one_jy_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   cmd = "rm -rf %s %s " % (one_jy_source_vis_name,one_jy_source_vis_name_assassin)
   print(cmd)
   os.system(cmd)

def cleanup_images_and_vis_assassin(array_label,lst,freq_MHz,pol):
   print("###################doing cleanup########################")
   #print("lst %s hrs" % lst)
   lst_deg = (float(lst)/24)*360.
   #print("lst %s deg" % lst_deg)
   year=2000
   month=1
   day=1
   hour=np.floor(float(lst))
   minute=np.floor((float(lst)-hour) * 60.)
   second=((float(lst)-hour) * 60. - minute) * 60.
      
   date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
  
   
   base_vis_name = "%s_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_image_name = "%s_no_src_LST_%03d_%0.3f_MHz.map" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_beam_name = "%s_no_src_LST_%03d_%0.3f_MHz.beam" % (array_label,lst_deg,freq_MHz)
   eda_model_no_source_image_fits_name = "%s_no_src_LST_%03d_%0.3f_MHz.fits" % (array_label,lst_deg,freq_MHz)
      
   cmd = "rm -rf %s %s %s %s " % (base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
   print(cmd)
   os.system(cmd)
   

   model_sky_vis = "%s_plus_sky_LST_%03d_%s_pol_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   model_global_signal_vis = "%s_plus_G_LST_%03d_%s_pol_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   eda_model_noise_vis_name = "%s_N_LST_%03d_%s_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
   cmd = "rm -rf %s %s %s" % (model_sky_vis,model_global_signal_vis,eda_model_noise_vis_name)
   print(cmd)
   os.system(cmd)
              

   one_jy_source_vis_name =  "%s_one_jy_source_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   one_jy_source_vis_name_assassin =  "%s_one_jy_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
   cmd = "rm -rf %s %s " % (one_jy_source_vis_name,one_jy_source_vis_name_assassin)
   print(cmd)
   os.system(cmd)

def global_sig_func(nu_array,C=0.,A=1.,delta_nu=20.,nu_c=78.):
   S_21 = C*nu_array + A*np.exp(-(nu_array - nu_c)**2/(2*(delta_nu)**2))
   return S_21
   
#def global_sig_EDGES_func(nu_array,A_EDGES=0.52,tau_EDGES=6.5,nu_nought_EDGES=78.3,w_EDGES=20.7):
def global_sig_EDGES_func(nu_array,A_EDGES):
   tau_EDGES = 6.5
   nu_nought_EDGES = 78.3
   w_EDGES = 20.7
   
   B_EDGES = (4 * (nu_array - nu_nought_EDGES)**2 ) / (w_EDGES**2) * np.log(-(1/tau_EDGES)*np.log((1 + np.exp(-tau_EDGES))/2.))
   S_21_EDGES = - A_EDGES * ((1-np.exp(-tau_EDGES*np.exp(B_EDGES)))/(1-np.exp(-tau_EDGES)))
   return S_21_EDGES
   
   #def diffuse_global_foreground_func(nu_array,coeff_list):
def diffuse_global_foreground_func_order_5(nu_array,a0,a1,a2,a3,a4):
   #polynomial = np.poly1d(coeff_array)
   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5)
   return polynomial

def diffuse_global_foreground_func_order_6(nu_array,a0,a1,a2,a3,a4,a5):
   #polynomial = np.poly1d(coeff_array)
   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5) + a5 * nu_array**(5-2.5)
   return polynomial

def diffuse_global_foreground_func_order_7(nu_array,a0,a1,a2,a3,a4,a5,a6):
   #polynomial = np.poly1d(coeff_array)
   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5) + a5 * nu_array**(5-2.5) + a6 * nu_array**(6-2.5)
   return polynomial   

def global_sig_EDGES_and_diffuse_fg_func_order_5(nu_array,A_EDGES,a0,a1,a2,a3,a4):
   tau_EDGES = 6.5
   nu_nought_EDGES = 78.3
   w_EDGES = 20.7
   
   B_EDGES = (4 * (nu_array - nu_nought_EDGES)**2 ) / (w_EDGES**2) * np.log(-(1/tau_EDGES)*np.log((1 + np.exp(-tau_EDGES))/2.))
   S_21_EDGES = - A_EDGES * ((1-np.exp(-tau_EDGES*np.exp(B_EDGES)))/(1-np.exp(-tau_EDGES)))

   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5)

   total_signal = S_21_EDGES + polynomial
   return total_signal

def global_sig_EDGES_and_diffuse_fg_func_order_6(nu_array,A_EDGES,a0,a1,a2,a3,a4,a5):
   tau_EDGES = 6.5
   nu_nought_EDGES = 78.3
   w_EDGES = 20.7
   
   B_EDGES = (4 * (nu_array - nu_nought_EDGES)**2 ) / (w_EDGES**2) * np.log(-(1/tau_EDGES)*np.log((1 + np.exp(-tau_EDGES))/2.))
   S_21_EDGES = - A_EDGES * ((1-np.exp(-tau_EDGES*np.exp(B_EDGES)))/(1-np.exp(-tau_EDGES)))

   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5) + a5 * nu_array**(5-2.5)

   total_signal = S_21_EDGES + polynomial
   return total_signal


def global_sig_EDGES_and_diffuse_fg_func_order_7(nu_array,A_EDGES,a0,a1,a2,a3,a4,a5,a6):
   tau_EDGES = 6.5
   nu_nought_EDGES = 78.3
   w_EDGES = 20.7
   
   B_EDGES = (4 * (nu_array - nu_nought_EDGES)**2 ) / (w_EDGES**2) * np.log(-(1/tau_EDGES)*np.log((1 + np.exp(-tau_EDGES))/2.))
   S_21_EDGES = - A_EDGES * ((1-np.exp(-tau_EDGES*np.exp(B_EDGES)))/(1-np.exp(-tau_EDGES)))

   polynomial = a0 * nu_array**(0-2.5) + a1 * nu_array**(1-2.5) + a2 * nu_array**(2-2.5) + a3 * nu_array**(3-2.5) + a4 * nu_array**(4-2.5) + a5 * nu_array**(5-2.5) + a6 * nu_array**(6-2.5)

   total_signal = S_21_EDGES + polynomial
   return total_signal

def plot_iso_ant_int_response(plot_only=True):
   X_iso_parallel_array_filename = "X_iso_parallel_array.npy" 
   baseline_length_lambda_array_filename = "baseline_length_lambda_fig1_array.npy"
   if not plot_only:
      wavelength = 1.
      #baseline_length_m_array = np.arange(0,2,0.1)
      #baseline_length_lambda_array = baseline_length_m_array / wavelength
      baseline_length_lambda_array = np.arange(0,2,0.01)
      X_iso_parallel_array = np.full(len(baseline_length_lambda_array),np.nan)
      for baseline_length_lambda_index,baseline_length_lambda in enumerate(baseline_length_lambda_array):
         #for fig1 of paper
         n_pix = hp.nside2npix(NSIDE)
         pixel_solid_angle = (4.*np.pi) / n_pix
         hpx_index_array = np.arange(0,n_pix,1)
         iso_beam_map = np.full(hp.nside2npix(NSIDE),1.)             
         
         baseline_theta_rad = np.pi/2.
         baseline_phi_rad = 0.
         baseline_vector_for_dot_array = baseline_length_lambda * hp.ang2vec(baseline_theta_rad,baseline_phi_rad)
         
         #baseline_vector_for_dot_array_mag = np.linalg.norm(baseline_vector_for_dot_array)
         #print baseline_vector_for_dot_array_mag
         
         #Need to rotate all these vectors by -pi/2, just like the hpx beam map, since hp doesnt use alt/az, so theta=pi/2 is actually pointing at the zenith in orthographic proj
         #https://vpython.org/contents/docs/VisualIntro.html
         #rotate about x axis
         #forget about rotating the vectors its the phase angle hpx array you need to rotate!
         rot_axis = [0,1,0]
         rot_theta = 0. #np.pi / 4.
         
         #baseline_vector_array_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_unit)
         #print baseline_vector_array_unit_rotated[0,:]
      
         sky_vector_array = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array)))
         #sky_vector_array_unrotated = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array)))
         #sky_vector_array_unrotated_test = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array[1])))
         #print sky_vector_array_unrotated_test
         #sys.exit()
         #do the same rotation for the sky vector
         #sky_vector_array = rotate_vector(rot_axis,rot_theta,sky_vector_array_unrotated)
                     
                     
                     
         #b_dot_r_single = np.dot(baseline_vector_test_for_dot_array[4],sky_vector_array[4])
         #print b_dot_r_single
         #print baseline_vector_for_dot_array[0:3]
         b_dot_r_array = (baseline_vector_for_dot_array * sky_vector_array).sum(axis=1)
      
         phase_angle_array = 2.*np.pi*b_dot_r_array/wavelength
         
         element_iso_array = iso_beam_map*np.exp(-1j*phase_angle_array)
       
         
         #This one gives an answer that is almost exactly 6 PI too big .... (making this smaller makes Tsky smaller)
         #X_short_parallel = (1./(4.*np.pi)) * np.sum(element_short_parallel_array) * (2.*np.pi/float(n_pix))
         
         #this one gives approximately the right answer ....  no exactly!
         X_iso_parallel =  np.sum(element_iso_array) * pixel_solid_angle # (4.*np.pi/float(n_pix))
         
         print(baseline_length_lambda)
         print(X_iso_parallel.real)
         
         X_iso_parallel_array[baseline_length_lambda_index] = X_iso_parallel.real
      
      X_iso_parallel_array_max = np.nanmax(X_iso_parallel_array)
      X_iso_parallel_array_norm = X_iso_parallel_array / X_iso_parallel_array_max
         
      
      np.save(X_iso_parallel_array_filename,X_iso_parallel_array_norm) 
      np.save(baseline_length_lambda_array_filename,baseline_length_lambda_array)
   
   else:
      baseline_length_lambda_array = np.load(baseline_length_lambda_array_filename)
      X_iso_parallel_array_norm = np.load(X_iso_parallel_array_filename)
      
   plt.clf()
   map_title="X_iso_parallel"
   plt.plot(baseline_length_lambda_array,X_iso_parallel_array_norm)
   plt.ylabel("Interferometer response")
   plt.xlabel("Baseline length (wavelengths)")
   fig_name="X_iso_parallel.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   
   
def plot_S21(nu_array=None,C=0.,A=1.,delta_nu=20.,nu_c=78.):
   #Global signal
   #Use S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
   S_21 = global_sig_func(nu_array,C,A,delta_nu,nu_c)
   plt.clf()
   map_title="S_21 vs freq"
   plt.plot(nu_array,S_21)
   plt.ylabel("Tb (K)")
   plt.xlabel("freq (MHz)")
   fig_name="s_21_vs_freq.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   return S_21
   
   
def plot_S21_EDGES(nu_array,A_EDGES = 0.52,tau_EDGES = 6.5,nu_nought_EDGES = 78.3,w_EDGES = 20.7):
   
   #B_EDGES = (4 * (nu_array - nu_nought_EDGES)**2 ) / (w_EDGES**2) * np.log(-(1/tau_EDGES)*np.log((1 + np.exp(-tau_EDGES))/2.))
   #S_21_EDGES = - A_EDGES * ((1-np.exp(-tau_EDGES*np.exp(B_EDGES)))/(1-np.exp(-tau_EDGES)))
   
   S_21_EDGES = global_sig_EDGES_func(nu_array,A_EDGES)
   
   plt.clf()
   map_title="S_21 vs freq"
   plt.plot(nu_array,S_21_EDGES)
   plt.ylabel("Tb (K)")
   plt.xlabel("freq (MHz)")
   fig_name="s_21_EDGES_vs_freq.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   return S_21_EDGES


def write_beam_fits_sin(cart_image,fitsname,lst=None):
   if lst==None:
      lst=60
   else:
      lst=lst
   hdu = fits.PrimaryHDU(cart_image)
   hdul = fits.HDUList([hdu])
   hdul.writeto(fitsname,clobber=True)
   header = hdul[0].header
   header['CTYPE1']  = 'RA---SIN'
   header['CRPIX1']  =  91
   header['CRVAL1']  = 60
   header['CDELT1']  = -0.636620
   header['CUNIT1']  = 'deg     '
   header['CTYPE2']  = 'DEC--SIN'
   header['CRPIX2']  = 91
   header['CRVAL2']  = -27
   header['CDELT2']  = 0.636620
   header['CUNIT2']  = 'deg  '
   #target_header = fits.Header.fromstring("NAXIS   =                    2\nNAXIS1  =                  180\nNAXIS2  =                  180\nCTYPE1  = 'RA---SIN'\nCRPIX1  =                   90\nCRVAL1  =                   -4\nCDELT1  =                    1\nCUNIT1  =' deg     '\nCTYPE2  = 'DEC--SIN'\nCRPIX2  =                   90\nCRVAL2  =                  -27\nCDELT2  =                    1\nCUNIT2  ='deg     '\nCOORDSYS= 'icrs    '\n", sep='\n')
   #target_header = fits.Header
   #pyfits.writeto(Power_pattern_average_interp_fitsname,cart_image,clobber=True)
   fits.update(fitsname,cart_image,header=header)
   print("wrote image %s" %  fitsname)

def recompute_ffts(pol,freq_MHz):

   if pol=='X':
      beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
   else:
      beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
   
   beam_plot_basename = beam_image_name.split('.')[0]
   #fft_savename = beam_plot_basename + "_beam_fft2.npy"
   #fft_savename_shift = beam_plot_basename + "_beam_fft2_shift.npy"
   
   #fft fits name
   fits_name_real = beam_plot_basename + "_beam_fft2_real_shift.fits"
   fits_name_imag = beam_plot_basename + "_beam_fft2_imag_shift.fits"
   fits_name_abs = beam_plot_basename + "_beam_fft2_abs_shift.fits"
   
   with fits.open(beam_image_name) as beam_hdulist:
      #beam_hdulist = fits.open(beam_image_name)
      beam_image_data = beam_hdulist[0].data
      beam_image_header = beam_hdulist[0].header
   
   #the first row and column are -0 for some reason leading to a non-symetric beam - remove the first row and first column (then you also end up with odd numbers...)
   
   #beam_image_data = beam_image_data[1:,1:]
   
   #Fourier response:
   image_data_shape = beam_image_data.shape
   #print image_data_shape
   beam_image_length = image_data_shape[0]
   image_centre = int(beam_image_length/2)

   ##test_data_1d = np.array([3.,2.,1.,0.,1.,2.])
   ##test_data_1d_fft = np.fft.fft(test_data_1d)
   ##print test_data_1d
   ##print test_data_1d_fft
   
   #pix_from_centre = 2
   pix_from_centre = int(beam_image_length/2)
   
   ##image_data_test = beam_image_data[image_centre-pix_from_centre:image_centre+pix_from_centre,image_centre-pix_from_centre:image_centre+pix_from_centre]
   ##print image_data_test
   ##image_data_test_fft2 = np.fft.fft2(image_data_test)
   ##print image_data_test_fft2
   ##image_data_test_fft2_imag = np.imag(image_data_test_fft2)
   ##print image_data_test_fft2_imag
   ##print np.max(abs(image_data_test_fft2_imag))
 
   #beam_image_data = beam_image_data[image_centre-pix_from_centre:image_centre+pix_from_centre,image_centre-pix_from_centre:image_centre+pix_from_centre]
 
   #zero pad the beam image for later multiplication
   #pad_width = int(padding_factor/2. * image_data_shape[0] - (image_data_shape[0]/2.)) 
   #beam_image_data_padded = np.pad(beam_image_data,pad_width,'constant', constant_values=0)
   #print "beam_image_data_padded size %s" % beam_image_data_padded.shape[0]
   #image_centre_padded = int(beam_image_data_padded.shape[0]/2.)
   
   padding = padding_factor*np.asarray(image_data_shape)
   #padding -= 1
   #print 'padding'
   #print padding
   
   #####Ronniys recipe:
   #####beam_image = gaussian
   #####padded_beam = numpy.pad(beam_image, padding_length)
   #####shifted_beam = numpy.fftshift(padded_beam)
   #####shifted_FT_beam = numpy.fft2(shifted_beam)
   #####FT_beam = numpy.ifftshift(shifted_FT_beam)

   padded_beam = np.pad(beam_image_data, (padding,padding),'constant', constant_values=(0, 0))
   shifted_beam = np.fft.fftshift(padded_beam)
   shifted_beam_fft2 = np.fft.fft2(shifted_beam)
   beam_fft2 = np.fft.ifftshift(shifted_beam_fft2)
   ##beam_fft2 = np.fft.fft2(beam_image_data,s=[180,180])
   ##beam_fft2 = np.fft.fft2(beam_image_data,s=padding)
   ##beam_fft2 = np.fft.fft2(beam_image_data)
   ##beam_fft2_shift = np.fft.fftshift(beam_fft2)
   beam_fft2_real = np.real(beam_fft2)
   beam_fft2_imag = np.imag(beam_fft2)
   #print beam_fft2_imag
   #print np.max(abs(beam_fft2_imag))
   beam_fft2_abs = abs(beam_fft2)

   #beam_fft2_real_shift_norm = beam_fft2_real_shift/np.max(beam_fft2_real_shift)
   pyfits.writeto(fits_name_real,beam_fft2_real,clobber=True)
   print("saved %s" % fits_name_real)
   
   pyfits.writeto(fits_name_imag,beam_fft2_imag,clobber=True)
   print("saved %s" % fits_name_imag)

   pyfits.writeto(fits_name_abs,beam_fft2_abs,clobber=True)
   print("saved %s" % fits_name_abs)
   #save the fft? - these are like a gig each - dont save!
   #np.save(fft_savename,beam_fft2)
   #np.save(fft_savename_shift,beam_fft2_shift)
   #print "saved %s" % fft_savename_shift
   #print "saved %s" % fft_savename
      
def plot_beam_and_weights(pol,freq_MHz):
   wavelength = 300./freq_MHz
   profile_plot_length = 500
   
   if pol=='X':
      beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
   else:
      beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
   
   beam_plot_basename = beam_image_name.split('.')[0]
   #fft_savename = beam_plot_basename + "_beam_fft2.npy"
   #fft_savename_shift = beam_plot_basename + "_beam_fft2_shift.npy"
   
   with fits.open(beam_image_name) as beam_hdulist:
      #beam_hdulist = fits.open(beam_image_name)
      beam_image_data = beam_hdulist[0].data
      beam_image_header = beam_hdulist[0].header
   
   #print beam_image_data[90,:]
   
   
   plt.clf()
   map_title="beam"
   plt.imshow(beam_image_data)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)


   #Fourier response real part:
   fits_name_real = beam_plot_basename + "_beam_fft2_real_shift.fits"  
   beam_fft2_real_shift_hdu_list = fits.open(fits_name_real)
   #print "opened %s" % fits_name
   beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data
   beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
    
   fits_name_abs = beam_plot_basename + "_beam_fft2_abs_shift.fits"  
   beam_fft2_abs_shift_hdu_list = fits.open(fits_name_abs)
   #print "opened %s" % fits_name
   beam_fft2_abs_shift = beam_fft2_abs_shift_hdu_list[0].data
   beam_fft2_abs_shift_norm = beam_fft2_abs_shift/beam_fft2_abs_shift.max()  

   beam_image_length = beam_image_data.shape[0]
   beam_fft2_length = beam_fft2_real_shift.shape[0]
   #print beam_image_length
   #print beam_fft2_real_length
   fft_centre_padded = int(beam_fft2_length/2.)
   ##at phase centre angle is small sin(theta) = theta
   #sine projection so half beam image corresponds to sin(90 deg) = 1
   theta_step_rad = 1.0/(beam_image_length/2.)
   spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
   spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
   #rows same as cols as image is square
   ##spatial_frequencies_rows = np.fft.fftfreq(beam_fft2.shape[0],d=D_step_wavelengths)
 

   #real
   
   plt.clf()
   map_title="beam fft2 real"
   plt.imshow(beam_fft2_real_shift)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_real.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   #profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 real profile"
   #u or v?
   #beam_fft2_profile = beam_fft2_abs_shift_log[:,image_centre_padded]
   beam_fft2_profile = beam_fft2_real_shift_norm[fft_centre_padded,:]
   beam_fft2_profile_other = beam_fft2_real_shift_norm[:,fft_centre_padded]
   
   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile real inner" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile[fft_centre_padded:fft_centre_padded+profile_plot_length])
   ##plt.plot(UU_array_wavelengths[0:profile_plot_length],beam_fft2_profile[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_real_inner.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile log abs real inner" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   beam_fft2_profile_abs_log = np.log10(np.abs(beam_fft2_profile))
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_abs_log[fft_centre_padded:fft_centre_padded+profile_plot_length])
   
   #plt.plot(UU_array_wavelengths[0:profile_plot_length],beam_fft2_profile_abs_log[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylim([0, 4.5])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_real_log_abs_inner.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
      
   #other axis inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile real inner other" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_other[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_real_inner_other.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile log abs real inner other" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   beam_fft2_profile_abs_log_other = np.log10(np.abs(beam_fft2_profile_other))
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_abs_log_other[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylim([0, 4.5])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_real_log_abs_inner_other.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #abs
   
   plt.clf()
   map_title="beam fft2 abs"
   plt.imshow(beam_fft2_abs_shift)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_abs.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   #profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 abs profile"
   #u or v?
   #beam_fft2_profile = beam_fft2_abs_shift_log[:,image_centre_padded]
   beam_fft2_profile = beam_fft2_abs_shift_norm[fft_centre_padded,:]
   beam_fft2_profile_other = beam_fft2_abs_shift_norm[:,fft_centre_padded]
   
   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile abs inner" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   #plt.plot(spatial_frequencies_cols_fftshift_wavelength[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile[fft_centre_padded:fft_centre_padded+profile_plot_length])
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_abs_inner.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile log abs inner" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   beam_fft2_profile_abs_log = np.log10(np.abs(beam_fft2_profile))
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_abs_log[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylim([0, 4.5])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_abs_log_inner.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
      
   #other axis inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile abs inner other" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_other[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_abs_inner_other.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile log abs inner other" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   beam_fft2_profile_abs_log_other = np.log10(np.abs(beam_fft2_profile_other))
   plt.plot(spatial_frequencies_cols_fftshift[fft_centre_padded:fft_centre_padded+profile_plot_length],beam_fft2_profile_abs_log_other[fft_centre_padded:fft_centre_padded+profile_plot_length])
   #plt.ylim([0, 4.5])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_profile_abs_log_inner_other.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

#this should be quickest!
def bisection(array,value):
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl
         
   
def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return (index,array[index])

#Jacks 'find closest'
def find_closet_uv(u=None,v=None,u_range=None,v_range=None,resolution=None):
    '''Given the float coords u,v, find the closest point in the gridding range arrays
       u_range,v_range, which have resolution = resolution
       Returns the index of the visibility in u_range and v_range, as well
       as the offset between the actual u,v location and the grid point values'''
    u_offs = np.abs(u_range - u)
    ##Find out where in the gridded u coords the current u lives;
    ##This is a boolean array of length len(u_offs)
    u_true = u_offs < resolution/2.0
    ##Find the index so we can access the correct entry in the container
    u_ind = np.where(u_true == True)[0]

    ##Use the numpy abs because it's faster (np.abs)
    v_offs = np.abs(v_range - v)
    v_true = v_offs < resolution/2
    v_ind = np.where(v_true == True)[0]

    ##If the u or v coord sits directly between two grid points,
    ##just choose the first one ##TODO choose smaller offset?
    if len(u_ind) == 0:
        u_true = u_offs <= resolution/2
        u_ind = np.where(u_true == True)[0]
    if len(v_ind) == 0:
        v_true = v_offs <= resolution/2
        v_ind = np.where(v_true == True)[0]
    u_ind,v_ind = u_ind[0],v_ind[0]

    ##TODO is this -ve legit??? Seems so...
    u_offs = u_range - u
    v_offs = v_range - v

    u_off = -(u_offs[u_ind] / resolution)
    v_off = -(v_offs[v_ind] / resolution)

    return u_ind,v_ind,u_off,v_off

   
#def compute_weights(lst_deg,pol,freq_MHz,recompute_fits=False):
#Need to do this on the fly now cause the concat files are too big!
def compute_weights(lst_list,freq_MHz_list,pol_list,signal_type_list=['diffuse','global','noise','gain_errors'],sky_model='gsm'):
   n_chan = len(freq_MHz_list)
   n_lst = len(lst_list)
   n_pol = len(pol_list)
   freq_MHz_end = freq_MHz_list[-1]
   lst_end = lst_list[-1]
   lst_deg_end = (float(lst_end)/24.)*360.
   for pol_index,pol in enumerate(pol_list):
      model_vis_name_base = "%s_LST_%03d_%s_%s_MHz" % (array_label,lst_deg_end,pol,int(freq_MHz_end))
      if 'noise' in signal_type_list:
          model_vis_name_base += '_N'
      if 'diffuse' in signal_type_list:
          model_vis_name_base += '_D_%s' % sky_model
      if 'global' in signal_type_list:
          model_vis_name_base += '_G' 
      if 'gain_errors' in signal_type_list:
          model_vis_name_base += '_GE'
                                   
      uvfits_name = "%s_concat_lsts.uvfits" % model_vis_name_base
      #print uvfits_name           
      hdulist = fits.open(uvfits_name)
      #hdulist.info()
      #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print info_string
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      visibilities = uvtable['DATA']

         
      #sky_uvfits_name = "eda_model_plus_sky_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))
      #model_global_signal_uvfits = "eda_model_plus_global_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz)) 
      #eda_model_noise_uvfits_name = "eda_model_noise_LST_%03d_%s_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))  
      #print sky_uvfits_name
      #hdulist = fits.open(sky_uvfits_name)
      ##hdulist.info()
      ##info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      ##print info_string
      #uvtable = hdulist[0].data
      #uvtable_header = hdulist[0].header
      #visibilities = uvtable['DATA']
      #print visibilities.shape
      n_vis = visibilities.shape[0]
      #print "n_vis %s" % n_vis
      #abs_vis = abs(visibilities)
      #log_abs_vis = np.log10(abs_vis)
      UU_s_array = uvtable['UU']
      UU_m_array = UU_s_array * c   
      VV_s_array = uvtable['VV']
      VV_m_array = VV_s_array * c
   
      #print "max UU_m_array %s" % np.max(UU_m_array)
      #print "max VV_m_array %s" % np.max(VV_m_array)
   
      
      #initialise an array of nans for the weights of shape n_pol,n_chan_n_vis
      weights_array = np.full((n_vis,n_lst,n_chan,n_pol), np.nan)
      #print weights_array.shape
            
      #weights_array_filename =  "weights_LST_%03d_%s_%s_MHz.npy" % (lst_deg,pol,int(freq_MHz))
      weights_array_filename = "weights_%s.npy" % model_vis_name_base
   
      #for pol_index,pol in enumerate(pol_list):
         #for chan_index,chan in enumerate(range(0,n_chan,1)):
      
      for chan_index,freq_MHz in enumerate(freq_MHz_list):
         wavelength = 300./float(freq_MHz)
         for lst_index,lst in enumerate(lst_list):  
            lst_deg = (float(lst)/24.)*360.
      
      
            if pol=='X':
               beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
            else:
               beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
            
            beam_plot_basename = beam_image_name.split('.')[0]
            #fft_savename = beam_plot_basename + "_beam_fft2.npy"
            #fft_savename_shift = beam_plot_basename + "_beam_fft2_shift.npy"
            
            with fits.open(beam_image_name) as beam_hdulist:
               #beam_hdulist = fits.open(beam_image_name)
               beam_image_data = beam_hdulist[0].data
               beam_image_header = beam_hdulist[0].header
            
            #Fourier response:
            fits_name= beam_plot_basename + "_beam_fft2_real_shift.fits"  
         
            beam_fft2_real_shift_hdu_list = fits.open(fits_name)
            #print "opened %s" % fits_name
            beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data
               
            beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
                 
            #fits_name= beam_plot_basename + "_beam_fft2_abs_shift_log.fits"
            #pyfits.writeto(fits_name,beam_fft2_abs_shift_log,clobber=True)
            #print "saved %s" % fits_name        
                     
            beam_image_length = beam_image_data.shape[0]
            beam_fft2_real_length = beam_fft2_real_shift.shape[0]
            fft_centre_padded = int(beam_fft2_real_length/2.)
            ##at phase centre angle is small sin(theta) = theta
            #sine projection so half beam image corresponds to sin(90 deg) = 1
            theta_step_rad = 1.0/(beam_image_length/2.)
            spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
            spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
            
            UU_wavelength_array = UU_m_array / wavelength
            VV_wavelength_array = VV_m_array / wavelength   
            
            for UU_wavelength_index, UU_wavelength in enumerate(UU_wavelength_array):
               
               VV_wavelength = VV_wavelength_array[UU_wavelength_index]
               #if UU_wavelength_index==0:
               #   print "print U %0.3f V_%0.3f pol %s freq %s MHz %s wavelength m " % (UU_wavelength,VV_wavelength,pol,freq_MHz,wavelength)
               #print "baseline U %s wavelengths" % UU_wavelength
               #print "baseline (U,V) %s,%s wavelengths" % (UU_wavelength,VV_wavelength)
               
               #only proceed with all this compute-expensive stuff if the uvdist is small i.e. less than 2 wavelengths
               uvdist_wavelengths = np.sqrt(UU_wavelength**2 + VV_wavelength**2)
               if uvdist_wavelengths < uvdist_wavelength_cutoff:
               
                  #print "uvdist %s wavelengths, smaller than threshold %s, proceeding" % (uvdist_wavelengths,uvdist_wavelength_cutoff)
                  
                  plot_basename = "U_%0.3f_V_%0.3f_pol_%s_%s_MHz" % (UU_wavelength,VV_wavelength,pol,freq_MHz)
                  
                  #find the nearest U,V in our Fourier response image
                  #This is going to be way too coarse, better to do this weighting analytically ... but then can't use Daniel's beam models...
                  #find the nearest U value in the Fourier beam response
                  #print spatial_frequencies_cols_fftshift
                  nearest_UU_index, nearest_UU_value = find_nearest(spatial_frequencies_cols_fftshift,UU_wavelength)
                  #print "nearest_UU_index,nearest_UU_value %s,%s" % (nearest_UU_index, nearest_UU_value)
                  nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
                  #print "nearest_VV_index,nearest_VV_value %s,%s" % (nearest_VV_index, nearest_VV_value)
                  
                  
                  ######Don't need to do this - use cheat way (shifting theorem below) - did for testing and is consistent
                  
                  ##make delta function image at this U,U coord
                  #delta_function_UV = beam_fft2*0.
                  #delta_function_UV[nearest_VV_index,nearest_UU_index] = 1.
                  #
                  #delta_function_UV_fft = np.fft.fft2(delta_function_UV)
                  #delta_function_UV_fft_shift = np.fft.fftshift(delta_function_UV_fft)
                  #beam_footprint_mult_delta_image_space = beam_image_data_padded * delta_function_UV_fft_shift
                  #
                  #beam_footprint_mult_delta_image_space_inv_fft = np.fft.ifft2(beam_footprint_mult_delta_image_space)
                  #beam_footprint_mult_delta_image_space_inv_fft_abs = np.abs(beam_footprint_mult_delta_image_space_inv_fft)
                  ##normalise to max of one
                  #beam_footprint_mult_delta_image_space_inv_fft_abs_norm = beam_footprint_mult_delta_image_space_inv_fft_abs/np.max(beam_footprint_mult_delta_image_space_inv_fft_abs)
                  #beam_footprint_mult_delta_image_space_inv_fft_abs_log = np.log10(beam_footprint_mult_delta_image_space_inv_fft_abs)
                  
                  #plt.clf()
                  #map_title="beam footprint Fourier"
                  #plt.imshow(c)
                  #plt.ylabel("log abs(vis)")
                  #plt.xlabel("UU (lambda)")
                  #fig_name= plot_basename + "_beam_footprint_fourier.png"
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)
              
                  
                  ###convolve this with the Fourier beam response
                  ###This is way too slow - instead fft the delta function -  multiply in image space, and inverse fft the result
                  ###beam_footprint_at_UV = signal.convolve2d(delta_function_UV, beam_fft2_abs, boundary='symm', mode='same')
            
                  #Now look at the value at the origin (U=V=0), for now, abs, but should we be looking at the real?
                  #This is the amount of 'leakage' into the zero spacing mode, which we want to use as a weight
                  #zero_spacing_leakage = beam_footprint_mult_delta_image_space_inv_fft_abs_norm[image_centre_padded,image_centre_padded]
                  #print "zero_spacing_leakage using fft is %s " % zero_spacing_leakage
                  
                  
                  #don't have to actually do any convolution or fft, just use the shifting theorem and read off the correct pixel in the beam fft image!
                  UU_cheat_index = int((beam_fft2_real_length/2.)-(nearest_UU_index-(beam_fft2_real_length/2.)))
                  VV_cheat_index = int((beam_fft2_real_length/2.)-(nearest_VV_index-(beam_fft2_real_length/2.)))
                  
                  #print "UU_cheat_index is %s, VV_cheat_index is %s " % (UU_cheat_index,VV_cheat_index)
                  
                  #are these around the right way?
                  zero_spacing_leakage_cheat = beam_fft2_real_shift_norm[VV_cheat_index,UU_cheat_index]
                  ###zero_spacing_leakage_cheat = beam_fft2_real_shift_norm[UU_cheat_index,VV_cheat_index]
                  #print "zero_spacing_leakage using cheat is %s " % zero_spacing_leakage_cheat
                  
                  #put this check in the extract signal function (so we still record the smaller and negative weights and can plot them to check against the beam )
                  #if zero_spacing_leakage_cheat > zero_spacing_leakage_threshold:
                     #print "zero_spacing_leakage is %s, using this baseline (%s of %s, chan %s)" % (zero_spacing_leakage_cheat,UU_wavelength_index,n_vis,int(freq_MHz))
                     #don't write a FITS file unless doing a subset i.e. making images for a talk
                     #fits_name= plot_basename + "_beam_footprint_fourier.fits"
                     #pyfits.writeto(fits_name,beam_footprint_mult_delta_image_space_inv_fft_abs_norm,clobber=True)
                     #print "saved %s" % fits_name            
                     
                     #now work out what to do with these weights! save in an array for now
                  weights_array[UU_wavelength_index,lst_index,chan_index,pol_index] = zero_spacing_leakage_cheat
                  #else:
                     #print "zero_spacing_leakage is %s, skipping this baseline (%s of %s, chan %s)" % (zero_spacing_leakage_cheat,UU_wavelength_index,n_vis,int(freq_MHz))
               #else:
                  #print "uv_dist is %s wavelengths, skipping this baseline (%s of %s, chan %s)" % (uvdist_wavelengths,UU_wavelength_index,n_vis,int(freq_MHz))
 
   np.save(weights_array_filename,weights_array)
   

def model_signal(lst_list,freq_MHz_list,pol_list,signal_type_list,outbase_name,poly_order,sky_model,array_ant_locations_filename,array_label):
   concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
   concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
         if 'noise' in signal_type_list:
            concat_output_name_base_X += '_N'
            concat_output_name_base_Y += '_N'
         if 'diffuse' in signal_type_list:
            concat_output_name_base_X += '_D_%s' % sky_model
            concat_output_name_base_Y += '_D_%s' % sky_model
         if 'diffuse_global' in signal_type_list:
             concat_output_name_base_X += '_DG'
             concat_output_name_base_Y += '_DG'
         if 'diffuse_angular' in signal_type_list:
             concat_output_name_base_X += '_DA'
             concat_output_name_base_Y += '_DA'
         if 'global' in signal_type_list:
            concat_output_name_base_X += '_G' 
            concat_output_name_base_Y += '_G'
         if 'gain_errors' in signal_type_list:
            concat_output_name_base_X += '_GE'
            concat_output_name_base_Y += '_GE'
        
         if do_image_and_subtr_from_simulated_data:
            uv_type_list = ['original','subtracted']
         else:
            uv_type_list = ['original']
         
         for uv_type in uv_type_list: 
            if pol=='X':  
               if uv_type=='original':                          
                  model_vis_name_base = concat_output_name_base_X
               else:
                  model_vis_name_base = concat_output_name_base_X + '_sub'
            else:
               if uv_type=='original': 
                  model_vis_name_base = concat_output_name_base_Y
               else:
                  model_vis_name_base = concat_output_name_base_Y + '_sub'
            model_vis_name_base += "_thresh_%0.2f" % (zero_spacing_leakage_threshold)        
   
            #signal_array_short_baselines_filename = outbase_name + "%s_signal.npy" % (model_vis_name_base)
            #signal_array_short_baselines = np.load(signal_array_short_baselines_filename)
            #signal_array_all_baselines_filename = "%s_signal_all_baselines.npy" % (model_vis_name_base)
            #signal_array_all_baselines = np.load(signal_array_all_baselines_filename)
            #print signal_array_all_baselines
            #signal_array_all_baselines_filename_abs = outbase_name + "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
            #signal_array_all_baselines_abs = np.load(signal_array_all_baselines_filename_abs)
            #signal_array_short_baselines_weighted_filename = outbase_name + "%s_signal_weighted.npy" % (model_vis_name_base)
            #signal_array_short_baselines_weighted = np.load(signal_array_short_baselines_weighted_filename)
   
             
            signal_array_short_baselines_Tb_filename = "%s_signal_Tb.npy" % (model_vis_name_base)
            signal_array_short_baselines_Tb = np.load(signal_array_short_baselines_Tb_filename)
            signal_array_all_baselines_Tb_filename = "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
            signal_array_all_baselines_Tb = np.load(signal_array_all_baselines_Tb_filename)
            #signal_array_all_baselines_filename_abs_Tb = "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
            #signal_array_all_baselines_abs_Tb = np.load(signal_array_all_baselines_filename_abs_Tb)
            signal_array_short_baselines_weighted_Tb_filename = "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
            signal_array_short_baselines_weighted_Tb = np.load(signal_array_short_baselines_weighted_Tb_filename)
        
            if ('diffuse' or 'diffuse_global' in signal_type_list):
               sky_averaged_diffuse_array_no_beam_lsts_filename =  "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base_X
               sky_averaged_diffuse_array_beam_X_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base_X
               sky_averaged_diffuse_array_beam_Y_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base_Y
               if pol == 'X':
                  sky_averaged_diffuse_array_beam_lsts_filename = sky_averaged_diffuse_array_beam_X_lsts_filename
               else:
                  sky_averaged_diffuse_array_beam_lsts_filename = sky_averaged_diffuse_array_beam_Y_lsts_filename
               sky_averaged_diffuse_array_no_beam_lsts_filename = sky_averaged_diffuse_array_no_beam_lsts_filename
               sky_averaged_diffuse_array_beam_lsts = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
               #sky_averaged_diffuse_array_no_beam_lsts = np.load(sky_averaged_diffuse_array_no_beam_lsts_filename)
            
               ##########
               #Sky - average from sky model stuff (not simulated through assassin sim, just sky time beam and average)
               #with beam
               coefs = poly.polyfit(freq_MHz_array, sky_averaged_diffuse_array_beam_lsts, poly_order)
               ffit = poly.polyval(freq_MHz_array, coefs)
      
               residual = ffit - sky_averaged_diffuse_array_beam_lsts
      
               
               plt.clf()
               plt.plot(freq_MHz_array,residual,label='residual')
               map_title="Residual for polynomial order %s fit to diffuse" % poly_order
               plt.ylabel("Residual Tb (K)")
               plt.xlabel("freq (MHz)")
               plt.legend(loc=1)
               fig_name= "%s_residual_sky_average_beam_poly_%s.png" % (model_vis_name_base,poly_order)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
   
               ##no beam
               #coefs = poly.polyfit(freq_MHz_array, sky_averaged_diffuse_array_no_beam_lsts, poly_order)
               #ffit = poly.polyval(freq_MHz_array, coefs)
      
               #residual = ffit - sky_averaged_diffuse_array_no_beam_lsts
      
               
               #plt.clf()
               #plt.plot(freq_MHz_array,residual,label='residual')
               #map_title="Residual for polynomial order %s fit to diffuse" % poly_order
               #plt.ylabel("Residual Tb (K)")
               #plt.xlabel("freq (MHz)")
               #plt.legend(loc=1)
               #fig_name= "%s_residual_sky_average_beam_poly_%s.png" % (model_vis_name_base,poly_order)
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)
            
               #in log log space:
               log_sky_array = np.log10(sky_averaged_diffuse_array_beam_lsts)
               log_freq_MHz_array = np.log10(freq_MHz_array)
               coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
               ffit = poly.polyval(log_freq_MHz_array, coefs)
               ffit_linear = 10**ffit
               
               #log_residual = log_signal_array_short_baselines - log_ffit
               residual_of_log_fit = ffit_linear - sky_averaged_diffuse_array_beam_lsts
               
               plt.clf()
               plt.plot(freq_MHz_array,residual_of_log_fit,label='residual of log fit')
               map_title="Residual for log polynomial order %s fit to diffuse" % poly_order
               plt.ylabel("Residual Tb (K)")
               plt.xlabel("freq (MHz)")
               plt.legend(loc=1)
               fig_name= "%s_log_fit_residual_sky_average_poly_%s.png" % (model_vis_name_base,poly_order)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name) 

               plt.clf()
               plt.plot(freq_MHz_array,ffit_linear,label='fit')
               plt.plot(freq_MHz_array,sky_averaged_diffuse_array_beam_lsts,label=sky_model)
               map_title="Polynomial order %s fit %s" % (poly_order,sky_model)
               plt.ylabel("Tb (K)")
               plt.xlabel("freq (MHz)")
               plt.legend(loc=1)
               fig_name= "%s_log_fit_sky_average_poly_%s.png" % (model_vis_name_base,poly_order)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
         
         
            ###################
            #all baselines abs
            #can't use the nan values:
            
         
            
            #good_idx = np.isfinite(signal_array_all_baselines_abs_Tb[0,:]) & np.isfinite(signal_array_all_baselines_abs_Tb[0,:])
            ##just use the first 70 frequencies to avoid that nasty sign change of the weights at 130 MHz
            #good_signal_array_all_baselines_abs_Tb = signal_array_all_baselines_abs_Tb[0,:][good_idx][0:70]
            
   
            #good_freq_MHz_array = freq_MHz_array[good_idx][0:70]
            #
            #coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_all_baselines_abs_Tb, poly_order)
            #ffit = poly.polyval(good_freq_MHz_array, coefs)
            #
            ##in log log space:
            #log_good_signal_array_all_baselines_abs_Tb = np.log10(good_signal_array_all_baselines_abs_Tb)
            #log_freq_MHz_array = np.log10(good_freq_MHz_array)
            #coefs = poly.polyfit(log_freq_MHz_array, log_good_signal_array_all_baselines_abs_Tb, poly_order)
            #ffit = poly.polyval(log_freq_MHz_array, coefs)
            #ffit_linear = 10**ffit
            #
            ##log_residual = log_signal_array_short_baselines - log_ffit
            #residual_of_log_fit = ffit_linear - good_signal_array_all_baselines_abs_Tb
            #
            #plt.clf()
            #plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual of log fit')
            #map_title="All base abs residual for log polynomial order %s fit to diffuse" % poly_order
            #plt.ylabel("All base abs residual Tb (K)")
            #plt.xlabel("freq (MHz)")
            #plt.legend(loc=1)
            #fig_name= "%s_log_fit_residual_poly_%s_all_base_abs.png" % (model_vis_name_base,poly_order)
            #figmap = plt.gcf()
            #figmap.savefig(fig_name)
            #print("saved %s" % fig_name) 

            ####################
            ##all baselines 
            ##can't use the nan values:
            #
            #good_idx = np.isfinite(signal_array_all_baselines[0,:]) & np.isfinite(signal_array_all_baselines[0,:])
            ###just use the first 70 frequencies to avoid that nasty sign change of the weights at 130 MHz
            #good_signal_array_all_baselines = signal_array_all_baselines[0,:][good_idx]
            #
            #
            #good_freq_MHz_array = freq_MHz_array[good_idx]
            #
            #coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_all_baselines, poly_order)
            #ffit = poly.polyval(good_freq_MHz_array, coefs)
            #
            ##in log log space:
            #log_good_signal_array_all_baselines = np.log10(good_signal_array_all_baselines)
            #log_freq_MHz_array = np.log10(good_freq_MHz_array)
            #coefs = poly.polyfit(log_freq_MHz_array, log_good_signal_array_all_baselines, poly_order)
            #ffit = poly.polyval(log_freq_MHz_array, coefs)
            #ffit_linear = 10**ffit
            
            ##log_residual = log_signal_array_short_baselines - log_ffit
            #residual_of_log_fit = ffit_linear - good_signal_array_all_baselines
            
            #plt.clf()
            #plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual of log fit')
            #map_title="All base residual for log polynomial order %s fit to diffuse" % poly_order
            #plt.ylabel("All base residual Tb (K)")
            #plt.xlabel("freq (MHz)")
            #plt.legend(loc=1)
            #fig_name= "%s_log_fit_residual_poly_%s_all_base.png" % (model_vis_name_base,poly_order)
            #figmap = plt.gcf()
            #figmap.savefig(fig_name)
            #print("saved %s" % fig_name) 
                     
        
            #######Unweighted (real)
            
            #can't use the nan values:
            good_idx = np.isfinite(signal_array_short_baselines_Tb[0,:]) & np.isfinite(signal_array_short_baselines_Tb[0,:])
            good_signal_array_short_baselines_Tb = signal_array_short_baselines_Tb[0,:][good_idx]
            
            good_freq_MHz_array = freq_MHz_array[good_idx]
            
            coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
            ffit = poly.polyval(good_freq_MHz_array, coefs)
   
            #plt.clf()
            #plt.plot(good_freq_MHz_array,ffit,label='model fit')
            #plt.plot(freq_MHz_array,signal_array_short_baselines_Tb[0,:],label='data')
            #map_title="Polynomial order %s fit %s" % (poly_order,sky_model)
            #plt.ylabel("Extracted Tb (K)")
            #plt.xlabel("freq (MHz)")
            #plt.legend(loc=1)
            #fig_name= "%s_model_short_poly_%s.png" % (model_vis_name_base,poly_order)
            #figmap = plt.gcf()
            #figmap.savefig(fig_name)
            #print("saved %s" % fig_name)
   
            residual = ffit - good_signal_array_short_baselines_Tb
   
            
            plt.clf()
            plt.plot(good_freq_MHz_array,residual,label='residual')
            map_title="Residual for polynomial order %s fit to diffuse" % poly_order
            plt.ylabel("Residual Tb (K)")
            plt.xlabel("freq (MHz)")
            plt.legend(loc=1)
            fig_name= "%s_residual_short_poly_%s.png" % (model_vis_name_base,poly_order)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            #in log log space:
            log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
            log_freq_MHz_array = np.log10(good_freq_MHz_array)
            coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
            ffit = poly.polyval(log_freq_MHz_array, coefs)
            ffit_linear = 10**ffit
         
            #log_residual = log_signal_array_short_baselines - log_ffit
            residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
            
            
            #plt.clf()
            #plt.plot(log_freq_MHz_array,log_residual,label='residual of log fit')
            #map_title="Log Residual for polynomial order %s fit to diffuse" % poly_order
            #plt.ylabel("Log residual Tb (K)")
            #plt.xlabel("log freq (MHz)")
            #plt.legend(loc=1)
            #fig_name= "%s_log_fit_residual_poly_%s.png" % (model_vis_name_base,poly_order)
            #figmap = plt.gcf()
            #figmap.savefig(fig_name)
            #print("saved %s" % fig_name)         
            
            #convert back to linear
            #inverse_log_residual = 10.**log_residual
            
            plt.clf()
            plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual of log fit')
            map_title="Residual for log polynomial order %s fit to diffuse" % poly_order
            plt.ylabel("Residual Tb (K)")
            plt.xlabel("freq (MHz)")
            plt.legend(loc=1)
            fig_name= "%s_log_fit_residual_poly_%s.png" % (model_vis_name_base,poly_order)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
   
            #####Weighted
            #can't use the nan values:
            good_idx = np.isfinite(signal_array_short_baselines_weighted_Tb[0,:]) & np.isfinite(signal_array_short_baselines_weighted_Tb[0,:])
            good_signal_array_short_baselines_Tb = signal_array_short_baselines_weighted_Tb[0,:][good_idx]
            
            good_freq_MHz_array = freq_MHz_array[good_idx]
            
            coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
            ffit = poly.polyval(good_freq_MHz_array, coefs)
   
            plt.clf()
            plt.plot(good_freq_MHz_array,ffit,label='model fit')
            plt.plot(freq_MHz_array,signal_array_short_baselines_weighted_Tb[0,:],label='data')
            map_title="Polynomial order %s fit %s" % (poly_order,sky_model)
            plt.ylabel("Weighted extracted Tb (K)")
            plt.xlabel("freq (MHz)")
            plt.legend(loc=1)
            fig_name= "%s_model_short_weighted_poly_%s.png" % (model_vis_name_base,poly_order)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            residual = ffit - good_signal_array_short_baselines_Tb

         
            plt.clf()
            plt.plot(good_freq_MHz_array,residual,label='residual')
            map_title="Weighted residual for polynomial order %s fit to diffuse" % poly_order
            plt.ylabel("Residual Tb (K)")
            plt.xlabel("freq (MHz)")
            plt.legend(loc=1)
            fig_name= "%s_residual_short_weighted_poly_%s.png" % (model_vis_name_base,poly_order)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
   
            #in log log space:
            log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
            log_freq_MHz_array = np.log10(good_freq_MHz_array)
            coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
            ffit = poly.polyval(log_freq_MHz_array, coefs)
            ffit_linear = 10**ffit
            
            #residual = log_signal_array_short_baselines - log_ffit
            residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
            
            plt.clf()
            plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
            map_title="Weighted residual for log polynomial order %s fit to diffuse" % poly_order
            plt.ylabel("Residual Tb (K)")
            plt.xlabel("freq (MHz)")
            plt.legend(loc=1)
            fig_name= "%s_log_fit_residual_weighted_poly_%s.png" % (model_vis_name_base,poly_order)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
         
def model_signal_from_assassin(lst_list,freq_MHz_list,pol_list,signal_type_list,sky_model,outbase_name,n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length_m,zero_spacing_leakage_threshold,poly_order):
   n_baselines = 1 #(by definition of the assassin design, each uvfits is for just one baseline)
   final_lst = lst_list[-1]
   final_lst_deg = (float(final_lst)/24.)*360. 
   final_freq = freq_MHz_list[-1] 
   freq_MHz_array = np.asarray(freq_MHz_list)
   radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
   zero_spacing_leakage_threshold_pc = zero_spacing_leakage_threshold * 100.
   max_arm_length_cm = max_arm_length_m * 100.
   min_arm_length_cm = min_arm_length_m * 100.
   for pol in pol_list:
      extracted_signal_name_base = "assassin_%s_%s_%03d_%03d_%s_LST_%03d_thresh_%02d_pc" % (n_ants_per_m_of_circumference,n_circles,max_arm_length_cm,min_arm_length_cm,pol,final_lst_deg,zero_spacing_leakage_threshold_pc)
      final_output_name_base = "%s_poly_%s" % (extracted_signal_name_base,poly_order)
      
      signal_array_short_baselines_weighted_Tb_filename = "%s_signal_Tb.npy" % (extracted_signal_name_base)
      signal_array_short_baselines_weighted_Tb = np.load(signal_array_short_baselines_weighted_Tb_filename)     
           
      #can't use the nan values:
      good_idx = np.isfinite(signal_array_short_baselines_weighted_Tb) & np.isfinite(signal_array_short_baselines_weighted_Tb)
      good_signal_array_short_baselines_Tb = signal_array_short_baselines_weighted_Tb[good_idx]
      
      good_freq_MHz_array = freq_MHz_array[good_idx]
           
      #in log log space:
      log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
      log_freq_MHz_array = np.log10(good_freq_MHz_array)
      coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      
      #residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
      
      plt.clf()
      plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
      map_title="Weighted residual for log polynomial order %s fit to diffuse" % poly_order
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s_log_fit_residual_weighted_poly_%s.png" % (final_output_name_base,poly_order)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 

def model_tsky_from_saved_data_eda2(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,EDA2_chan,n_obs,fine_chan_index,model_type,include_angular_info=True):
   centre_freq = float(freq_MHz_list[freq_MHz_index])
   centre_wavelength = 300./centre_freq
   #jy_to_K = (centre_wavelength**2) / (2. * k * 1.0e26) 

   fine_chan_width_MHz = fine_chan_width_Hz/1000000.    
   EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)
   obs_time_list = EDA2_obs_time_list_each_chan[freq_MHz_index]
   
   
   sky_averaged_diffuse_array_beam_lsts_filename = "%seda_model_%s_lst_2.00_hr_int_0.13_hr_N_D_gsm_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,pol)
   diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   #just doing one freq at a time right now for EDA2, not sure how this works with fine chans
   diffuse_global_value = diffuse_global_value_array[0] 
   
   fine_chan_width_MHz = fine_chan_width_Hz/1000000.
   freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
   wavelength_fine_chan = 300./freq_MHz_fine_chan
               
   #do for each obs_time:
   t_sky_K_list = []
   t_sky_error_K_list = []
   t_sky_K_list_flagged = []
   t_sky_error_K_list_flagged = []
   
   for EDA2_obs_time in obs_time_list:
      baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol_%s.npy" % (freq_MHz_fine_chan,pol,EDA2_obs_time)             
      X_short_parallel_array_filename = "X_uniform_resp_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)               
      X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)               
      X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)            
      real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
      Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
      
      if not os.path.exists(X_short_parallel_array_filename):
         print("%s does not exist, setting T_sky to nan" % X_short_parallel_array_filename)
         t_sky_K_list.append(np.nan)
         t_sky_error_K_list.append(np.nan)
         t_sky_K_list_flagged.append(np.nan)
         t_sky_error_K_list_flagged.append(np.nan)
         continue
      else:
         X_short_parallel_array = np.load(X_short_parallel_array_filename)
         print("loaded %s" % X_short_parallel_array_filename)
      
         real_vis_data_sorted_array = np.load(real_vis_data_sorted_array_filename).real
         print("loaded %s" % real_vis_data_sorted_array_filename)
   
         #The EDA2 data calibrated with miriad has a bunch of zeros, replace these with nans
         real_vis_data_sorted_array[np.isclose(real_vis_data_sorted_array,0)] = np.nan
         
         baseline_length_array_lambda_sorted_cut = np.load(baseline_length_array_lambda_sorted_cut_filename)
         print("loaded %s" % baseline_length_array_lambda_sorted_cut_filename)
   
         real_or_simulated_string = "EDA2"
   
         X_short_parallel_array_pure_parallel = np.load(X_short_parallel_array_filename_pure_parallel).real
         print("loaded %s" % X_short_parallel_array_filename_pure_parallel) 
         X_short_parallel_array_pure_inline = np.load(X_short_parallel_array_filename_pure_inline).real
         print("loaded %s" % X_short_parallel_array_filename_pure_inline)    
   
         
         if include_angular_info:
            Y_short_parallel_angular_array = np.load(Y_short_parallel_angular_array_filename).real
            print("loaded %s" % Y_short_parallel_angular_array_filename)
      
            #plot a histogram of Y values
            plt.clf()
            n, bins, patches = plt.hist(Y_short_parallel_angular_array)
            map_title="Histogram of Y values (angular response)" 
            fig_name= "hist_Y_angular_%0.3f_MHz_%s_pol_%s.png" % (centre_freq,pol,EDA2_obs_time)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            plt.close()
            print("saved %s" % fig_name)  
      
   
         ##plot X_vs_uvdist for vis and X
         #normalise both X and real vis to max 1
         if len(X_short_parallel_array_pure_inline) > 0:
            X_short_parallel_array_max_pure_inline = np.nanmax(X_short_parallel_array_pure_inline)
            X_short_parallel_array_norm_pure_inline = X_short_parallel_array_pure_inline / X_short_parallel_array_max_pure_inline
            
            X_short_parallel_array_max_pure_parallel = np.nanmax(X_short_parallel_array_pure_parallel)
            X_short_parallel_array_norm_pure_parallel = X_short_parallel_array_pure_parallel / X_short_parallel_array_max_pure_inline  
         
            X_short_parallel_array_max = np.nanmax(X_short_parallel_array)
            X_short_parallel_array_norm = X_short_parallel_array / X_short_parallel_array_max_pure_inline
            
            real_vis_data_sorted_max = np.nanmax(real_vis_data_sorted_array)
            real_vis_data_sorted_array_norm = real_vis_data_sorted_array / real_vis_data_sorted_max
            
            real_vis_data_sorted_array_norm_scaled = real_vis_data_sorted_array_norm * (-0.4)    #2. / (2.*np.pi)
            
            
            #plot X and pure inline and parallel for fig 1 of paper
            
            ## plot X and real vis vs baseline length
            plt.clf()
            plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='EDA-2')
            plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_parallel,label='parallel',color='red')
            plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_inline,label='inline',color='green')
            #plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_norm_offset,s=1,label='real vis norm')
            #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
            map_title="Response to uniform sky vs baseline length" 
            plt.xlabel("Baseline length (wavelengths)")
            plt.ylabel("Normalised visibility amplitude")
            plt.legend(loc=1)
            #plt.ylim([0, 20])
            fig_name= "X_vs_uv_dist_%0.3f_MHz_%s_pol_%s.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
            ##
         
            ## plot X and real vis vs baseline length for fig2
            plt.clf()
            plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
            plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_norm_scaled,s=1,label='Scaled %s visibility amplitude' % real_or_simulated_string)
            #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
            map_title="Response to uniform sky vs baseline length data" 
            plt.xlabel("Baseline length (wavelengths)")
            plt.ylabel("Visibility amplitude")
            plt.legend(loc=1)
            #plt.ylim([0, 20])
            fig_name= "X_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol_%s.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
            
            
            #try just using centre wavelength (defined above)
            jy_to_K = (wavelength_fine_chan**2) / (2. * k * 1.0e26) 
            
            if include_angular_info:
               Y_short_parallel_array_norm = Y_short_parallel_angular_array / X_short_parallel_array_max_pure_inline
         
               #full response
               #need to convert between Jy and K
               Y_short_parallel_angular_array_Jy = Y_short_parallel_angular_array / jy_to_K
            
               #need to update full response to include fine chans
               full_response_Jy = ((diffuse_global_value * X_short_parallel_array) / jy_to_K) + Y_short_parallel_angular_array_Jy
         
         
      
      
            #also include Y and the sum of X plus Y
            plt.clf()
            #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
            plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array,s=1,label='%s visibility amplitude' % real_or_simulated_string)
            #plt.scatter(baseline_length_array_lambda_sorted_cut,Y_short_parallel_array_norm,s=1,label='Expected angular response')
            
            #need to update update full response to include fine chans
            plt.scatter(baseline_length_array_lambda_sorted_cut,full_response_Jy,s=1,label='Expected full response Jy')
            ##plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
            map_title="Response to uniform sky vs baseline length data" 
            plt.xlabel("Baseline length (wavelengths)")
            plt.ylabel("Visibility amplitude")
            plt.legend(loc=1)
            #plt.ylim([0, 20])
            fig_name= "X_Y_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol_%s.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name) 
            
            
            
         if np.nansum(np.abs(X_short_parallel_array) > 0):
            if model_type=='OLS_fixed_intercept':
               model = sm.OLS(real_vis_data_sorted_array, X_short_parallel_array,missing='drop')
               results = model.fit()
               parameters = results.params
               #print parameters
               t_sky_jy = parameters[0]
               t_sky_error_jy = results.bse[0]
            elif model_type=='OLS_fixed_int_subtr_Y':
               #subtract Y from the data before fitting (should get rid of the angular variations)
               real_vis_data_sorted_array_subtr_Y = real_vis_data_sorted_array - Y_short_parallel_angular_array_Jy
               model = sm.OLS(real_vis_data_sorted_array_subtr_Y, X_short_parallel_array,missing='drop')
               results = model.fit()
               ##print results.summary()
               parameters = results.params
               #print parameters
               t_sky_jy = parameters[0]
               t_sky_error_jy = results.bse[0]
            elif model_type=='OLS_with_intercept':
               X_short_parallel_array = sm.add_constant(X_short_parallel_array)
               model = sm.OLS(real_vis_data_sorted_array, X_short_parallel_array,missing='drop')
               results = model.fit()
               ##print results.summary()
               parameters = results.params
               ##print parameters
               t_sky_jy = parameters[1]
               t_sky_error_jy = results.bse[1]
         else:
            print("X_short_parallel_array all NaNs, returning Tsky NaN")
            return(np.nan,np.nan,np.nan,np.nan,freq_MHz_fine_chan)
         
         t_sky_K = jy_to_K * t_sky_jy
         t_sky_error_K = jy_to_K * t_sky_error_jy
         print("t_sky_K is %0.4E +/- %0.04f K" % (t_sky_K,t_sky_error_K))
         fit_string = "y=%0.1fx" % t_sky_jy         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
         
         print("diffuse_global_value is %0.4E" % diffuse_global_value) 
         
         ratio_in_out = diffuse_global_value / t_sky_K
         print("ratio between input and output T_sky is %0.4f" % ratio_in_out )
         
         y_pos = np.max(results.fittedvalues)
         x_pos = 1.2 * np.min(X_short_parallel_array)
         
          
         #get rid of nans
         real_vis_data_sorted_array_nonans = real_vis_data_sorted_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
         X_short_parallel_array_nonans = X_short_parallel_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
         
         if (real_vis_data_sorted_array_nonans.shape==results.fittedvalues.shape):
            #plot in Jy
            plt.clf()
            if model_type=='OLS_with_intercept':
               plt.plot(X_short_parallel_array_nonans[:,1], real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans[:,1], results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
            elif model_type=="OLS_fixed_int_subtr_Y":
               real_vis_data_sorted_array_subtr_Y_nonans = real_vis_data_sorted_array_subtr_Y[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
               plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_subtr_Y_nonans,label='%s data - Y' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
            elif model_type=='OLS_fixed_intercept':
               plt.scatter(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans[0:len(results.fittedvalues)], results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
            else:
               plt.scatter(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')         
            
            map_title="Data and fit" 
            plt.xlabel("Expected global-signal response")
            plt.ylabel("Real component of visibility (Jy)")
            plt.legend(loc=1)
            plt.text(x_pos, y_pos, fit_string)
            #plt.ylim([0, 3.5])
            fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol_%s_%s.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time,model_type)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            plt.close()
            print("saved %s" % fig_name)  
            
            t_sky_K_list.append(t_sky_K)
            t_sky_error_K_list.append(t_sky_error_K)
      
            #FLAGGING bit
            #now use the fit to identify outliers probably due to rfi
            #subtract the model from the data
         
            real_vis_data_sorted_array_subtr_model = real_vis_data_sorted_array_nonans - results.fittedvalues
            #take the mean 
            real_vis_data_sorted_array_subtr_model_mean = np.nanmean(real_vis_data_sorted_array_subtr_model)
            real_vis_data_sorted_array_subtr_model_std = np.nanstd(real_vis_data_sorted_array_subtr_model)
         
            #mask values greater than 5 sigma away from mean
            thresh = 5.* real_vis_data_sorted_array_subtr_model_std
            real_vis_data_sorted_array_flagged = np.copy(real_vis_data_sorted_array_nonans)
            real_vis_data_sorted_array_flagged[(np.abs(real_vis_data_sorted_array_subtr_model) > thresh)] = np.nan
            
            
         
            #get rid of nans
            #real_vis_data_sorted_array_flagged = real_vis_data_sorted_array_flagged[np.argwhere(np.logical_not(np.isnan(real_vis_data_sorted_array_flagged)))]
            #X_short_parallel_array_flagged = X_short_parallel_array_nonans[np.argwhere(np.logical_not(np.isnan(real_vis_data_sorted_array_flagged)))]
            
            if (X_short_parallel_array_nonans.shape[0]>0):
               model = sm.OLS(real_vis_data_sorted_array_flagged, X_short_parallel_array_nonans,missing='drop')
               results = model.fit()
               ##print results.summary()
               parameters = results.params
               print parameters
            
               t_sky_jy = parameters[0]
               t_sky_error_jy = results.bse[0]
               
               X_short_parallel_array_nonans_nonans = X_short_parallel_array_nonans[np.logical_not(np.isnan(real_vis_data_sorted_array_flagged))]
               
               #convert to K
               t_sky_K_flagged = jy_to_K * t_sky_jy
               t_sky_error_K_flagged = jy_to_K * t_sky_error_jy
               print("t_sky_K_flagged is %0.4E +/- %0.04f K" % (t_sky_K_flagged,t_sky_error_K_flagged))
               
               fit_string = "y=%0.1fx" % t_sky_jy         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
               fit_string_K = "y=%0.1fx" % t_sky_K_flagged
              
               print("diffuse_global_value is %0.4E" % diffuse_global_value) 
            
               ratio_in_out = diffuse_global_value / t_sky_K_flagged
               print("ratio between input and output T_sky_flagged is %0.4f" % ratio_in_out )
                
               fitted_values_K = results.fittedvalues * jy_to_K
                 
               y_pos = np.max(results.fittedvalues)
               x_pos = 1.2 * np.min(X_short_parallel_array)
            
               y_pos_K = 1300.0 # np.max(fitted_values_K)
               x_pos_K = 0.60 # 1.2 * np.min(X_short_parallel_array)
            
               #plot in Jy
               plt.clf()
               plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_flagged,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
            
               map_title="Flagged data and fit" 
               plt.xlabel("Expected global-signal response")
               plt.ylabel("Real component of visibility (Jy) flagged")
               plt.legend(loc=1)
               plt.text(x_pos, y_pos, fit_string)
               #plt.ylim([0, 3.5])
               fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_flagged.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time,model_type)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               plt.close()
               print("saved %s" % fig_name) 
               
               #plot in K
               
               real_vis_data_sorted_array_flagged_K = real_vis_data_sorted_array_flagged * jy_to_K
               
               
               plt.clf()
               plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_flagged_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
               plt.plot(X_short_parallel_array_nonans_nonans, fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
            
               map_title="Flagged data and fit" 
               plt.xlabel("Expected global-signal response")
               plt.ylabel("Real component of visibility (K)")
               plt.legend(loc=1)
               plt.text(x_pos_K, y_pos_K, fit_string_K)
               #plt.ylim([0, 3.5])
               fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_flagged_K.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time,model_type)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               plt.close()
               print("saved %s" % fig_name) 
               
               sys.exit()
               
            else:
               t_sky_jy = np.nan
               t_sky_error_jy = np.nan
               t_sky_K_flagged = np.nan
               t_sky_error_K_flagged = np.nan
               
             
            
   
            t_sky_K_list_flagged.append(t_sky_K_flagged)
            t_sky_error_K_list_flagged.append(t_sky_error_K_flagged)
         else:
            print("not including obs %s (chan %s), bad data identified in flagging " % (EDA2_obs_time,EDA2_chan))
            t_sky_K_list_flagged.append(np.nan)
            t_sky_error_K_list_flagged.append(np.nan)
            t_sky_K_list.append(np.nan)
            t_sky_error_K_list.append(np.nan)
            
   t_sky_K_array = np.asarray(t_sky_K_list)
   t_sky_error_K_array = np.asarray(t_sky_error_K_list)
   
   mean_t_sky_K = np.nanmean(t_sky_K_array)
   std_dev_t_sky_K = np.nanstd(t_sky_K_array)

   t_sky_K_array_flagged = np.asarray(t_sky_K_list_flagged)
   t_sky_error_K_array_flagged = np.asarray(t_sky_error_K_list_flagged)
   
   mean_t_sky_K_flagged = np.nanmean(t_sky_K_array_flagged)
   std_dev_t_sky_K_flagged = np.nanstd(t_sky_K_array_flagged)
   
   

   return(mean_t_sky_K.real,std_dev_t_sky_K,mean_t_sky_K_flagged.real,std_dev_t_sky_K_flagged,freq_MHz_fine_chan)

         
def model_tsky_from_saved_data(freq_MHz_list,freq_MHz_index,lst_hrs,pol,signal_type_list,sky_model,array_label,model_type,EDA2_data=False,EDA2_chan='None',n_obs_concat=1,fine_chan_index=0,edge_chan=False,wsclean=False,fast=False):
   freq_MHz = freq_MHz_list[freq_MHz_index]
   centre_freq = float(freq_MHz)
   fine_chan_width_MHz = fine_chan_width_Hz/1000000.
   
   if EDA2_data:
      n_edge_chans_omitted = 5 #two at start and 3 at end
      n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
   else:
      n_fine_chans_used = 1
      n_edge_chans_omitted = 0
   
   bandwidth = (n_fine_chans_used + 1) * fine_chan_width_MHz
   
   if EDA2_data:
      freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
      #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
      #freq_MHz_fine_chan = centre_freq + (bandwidth/2.) - (fine_chan_index)*fine_chan_width_MHz 
   else:
      freq_MHz_fine_chan = centre_freq     
   wavelength = 300./float(freq_MHz_fine_chan)  
   #wavelength = 300./float(centre_freq)  
   
   if edge_chan:
      return(np.nan,np.nan,np.nan,np.nan,freq_MHz_fine_chan) 
   
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''

   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
   if 'global_unity' in signal_type_list:
       signal_type_postfix += '_GU'
       concat_output_name_base += '_GU'    
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse_global at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DG_%s' % sky_model
          concat_output_name_base += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse_angular at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DA_%s' % sky_model
          concat_output_name_base += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print("cant have global and global_EDGES in signal_type_list")
          sys.exit()
       else:
          signal_type_postfix += '_G' 
          concat_output_name_base += '_G' 
   if 'global_EDGES' in signal_type_list:
       signal_type_postfix += '_ED' 
       concat_output_name_base += '_ED' 
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'
   
   if EDA2_data==True:
      EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)
   else:
      EDA2_chan_dir = ''
   #wavelength = 300./float(freq_MHz)
   
   #get the diffuse global diffuse value used in the simulation (from gsm)
   if EDA2_data==True:
      sky_averaged_diffuse_array_beam_lsts_filename = "%s%s_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,concat_output_name_base)
   else:
      sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % (concat_output_name_base)
   #sky_averaged_diffuse_array_no_beam_lsts_filename = "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base
   freq_MHz_index = int(freq_MHz - 50)
   diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   print("loaded %s" % sky_averaged_diffuse_array_beam_lsts_filename)
   if EDA2_data==True:
      #just doing one eda freq at a time for now
      ##
      diffuse_global_value = diffuse_global_value_array[0]
   else:
      diffuse_global_value = diffuse_global_value_array[freq_MHz_index]
   
   #in here put bit to read X from miriad_sim_uvfits
   if not fast:
      if EDA2_data:
         X_short_parallel_array_filename = "X_short_parallel_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
      else:
         X_short_parallel_array_filename = "X_short_parallel_array_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
   else:
      X_short_parallel_array_filename = "unity_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
      real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
      baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)  
         
   X_short_parallel_array = np.load(X_short_parallel_array_filename)
   print("loaded %s" % X_short_parallel_array_filename)
   
   real_vis_data_sorted_array = np.load(real_vis_data_sorted_array_filename).real
   print("loaded %s" % real_vis_data_sorted_array_filename)

   baseline_length_array_lambda_sorted_cut = np.load(baseline_length_array_lambda_sorted_cut_filename)
   print("loaded %s" % baseline_length_array_lambda_sorted_cut_filename)
   
   if fast:
      X_short_parallel_array = np.concatenate(X_short_parallel_array)
      real_vis_data_sorted_array = np.concatenate(real_vis_data_sorted_array)

   if EDA2_data==True:
      real_or_simulated_string = "EDA2"
   else:
      real_or_simulated_string = "simulated"
         
   if not fast:
      X_short_parallel_array_pure_parallel = np.load(X_short_parallel_array_filename_pure_parallel).real
      print("loaded %s" % X_short_parallel_array_filename_pure_parallel) 
      X_short_parallel_array_pure_inline = np.load(X_short_parallel_array_filename_pure_inline).real
      print("loaded %s" % X_short_parallel_array_filename_pure_inline)    
   
      if include_angular_info:
         Y_short_parallel_angular_array = np.load(Y_short_parallel_angular_array_filename).real
         print("loaded %s" % Y_short_parallel_angular_array_filename)
      
         #plot a histogram of Y values
         plt.clf()
         n, bins, patches = plt.hist(Y_short_parallel_angular_array)
         map_title="Histogram of Y values (angular response)" 
         fig_name= "hist_Y_angular_%0.3f_MHz_%s_pol%s.png" % (freq_MHz,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         plt.close()
         print("saved %s" % fig_name)  
   

      ##plot X_vs_uvdist for vis and X
      #normalise both X and real vis to max 1
      if len(X_short_parallel_array_pure_inline) > 0:
         X_short_parallel_array_max_pure_inline = np.max(X_short_parallel_array_pure_inline)
         #print("X_short_parallel_array_max_pure_inline %E" % X_short_parallel_array_max_pure_inline)
         X_short_parallel_array_norm_pure_inline = X_short_parallel_array_pure_inline / X_short_parallel_array_max_pure_inline
         
         X_short_parallel_array_max_pure_parallel = np.max(X_short_parallel_array_pure_parallel)
         #print("X_short_parallel_array_max_pure_parallel %E" % X_short_parallel_array_max_pure_parallel)
         X_short_parallel_array_norm_pure_parallel = X_short_parallel_array_pure_parallel / X_short_parallel_array_max_pure_inline  
      
         X_short_parallel_array_max = np.max(X_short_parallel_array)
         #print("X_short_parallel_array_max %E" % X_short_parallel_array_max)
         X_short_parallel_array_norm = X_short_parallel_array / X_short_parallel_array_max_pure_inline
         
         #[0:n_baselines_included]
         real_vis_data_sorted_max = np.max(real_vis_data_sorted_array)
         #print("real_vis_data_sorted_max %E" % real_vis_data_sorted_max)
         real_vis_data_sorted_array_norm = real_vis_data_sorted_array / real_vis_data_sorted_max
         
         #offset X and real_vis by some arbitrary amount 0.5?
         #real_vis_data_sorted_array_norm_offset = real_vis_data_sorted_array_norm + 0.5
         #need to multiply (scale) not add!
         #real_vis_data_sorted_array_norm_scaled = real_vis_data_sorted_array_norm * 2. / (2.*np.pi)
         #for fig4:
         real_vis_data_sorted_array_norm_scaled = X_short_parallel_array_norm * (-0.4)
         
         #plot X and pure inline and parallel for fig 1 of paper
         
         ## plot X and real vis vs baseline length
         plt.clf()
         plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='EDA-2')
         plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_parallel,label='parallel',color='red')
         plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_inline,label='inline',color='green')
         #plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_norm_offset,s=1,label='real vis norm')
         #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="Response to uniform sky vs baseline length" 
         plt.xlabel("Baseline length (wavelengths)")
         plt.ylabel("Normalised visibility amplitude")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "X_vs_uv_dist_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
         ##
      
         ## plot X and real vis vs baseline length for fig2
         
         plt.clf()
         plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Global response (unity sky)')
         plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_norm_scaled,s=1,label='Visibility amplitude (simulations)') #% real_or_simulated_string)
         #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="Response to uniform sky vs baseline length data" 
         plt.xlabel("Baseline length (wavelengths)")
         plt.ylabel("Visibility amplitude")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "X_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
         
         
         
         jy_to_K = (wavelength**2) / (2. * k * 1.0e26) 
         
         if include_angular_info:
            Y_short_parallel_array_norm = Y_short_parallel_angular_array / X_short_parallel_array_max_pure_inline
      
            #full response
            #need to convert between Jy and K
         
            Y_short_parallel_angular_array_Jy = Y_short_parallel_angular_array / jy_to_K
         
            X_short_parallel_array_diffuse_Jy = (diffuse_global_value * X_short_parallel_array) / jy_to_K
            
            #need to update full response to include fine chans
            full_response_Jy =  X_short_parallel_array_diffuse_Jy + Y_short_parallel_angular_array_Jy
      
      
   
   
         #also include Y and the sum of X plus Y
         plt.clf()
         #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
         plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array,s=1,label='%s visibility amplitude' % real_or_simulated_string)
         #plt.scatter(baseline_length_array_lambda_sorted_cut,Y_short_parallel_array_norm,s=1,label='Expected angular response')
         
         #need to update update full response to include fine chans
         plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_diffuse_Jy,s=1,label='Expected uniform diffuse response Jy')
         plt.scatter(baseline_length_array_lambda_sorted_cut,full_response_Jy,s=1,label='Expected full response Jy')
         ##plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="Response to uniform sky vs baseline length data" 
         plt.xlabel("Baseline length (wavelengths)")
         plt.ylabel("Visibility amplitude")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "X_Y_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
         
         #Repeat in K
         #also include Y and the sum of X plus Y
         real_vis_data_sorted_array_K = real_vis_data_sorted_array * jy_to_K
         X_short_parallel_array_diffuse_Jy_K =  X_short_parallel_array_diffuse_Jy * jy_to_K
         full_response_Jy_K = full_response_Jy * jy_to_K
         
         plt.clf()
         #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
         plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_K,s=1,label='%s visibility amplitude' % real_or_simulated_string)
         #plt.scatter(baseline_length_array_lambda_sorted_cut,Y_short_parallel_array_norm,s=1,label='Expected angular response')
         
         #need to update update full response to include fine chans
         plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_diffuse_Jy_K,s=1,label='Expected uniform diffuse response')
         plt.scatter(baseline_length_array_lambda_sorted_cut,full_response_Jy_K,s=1,label='Expected full response')
         ##plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="Response to uniform sky vs baseline length data" 
         plt.xlabel("Baseline length (wavelengths)")
         plt.ylabel("Visibility amplitude (K)")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "X_Y_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s_K.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
         
         
   ##Assign vis and X values a group according to the Y_angular response
   ##Dont need this anymore, is for mixedlm which didnt work.
   #Y_angular_group_array = X_short_parallel_array * 0
   #X_group_array = X_short_parallel_array * 0
   ###use the hist output to split up X according to Y values
   #for Y_bin_index,Y_bin in enumerate(n):
   #   start_value = bins[Y_bin_index]
   #   end_value = bins[Y_bin_index+1]
   #   Y_angular_group_array[(Y_short_parallel_angular_array >= start_value) & (Y_short_parallel_angular_array < end_value)] = Y_bin_index
    
   #for OLS_fixed_int_min_vis:
   #plot a histogram of X values
   plt.clf()
   n_X, X_bins, X_patches = plt.hist(X_short_parallel_array)
   X_bin_width = X_bins[1] - X_bins[0]
   #print X_bin_width
   map_title="Histogram of X values (angular response)" 
   fig_name= "hist_X_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name) 
   
   ##split X up into bins  (this didn't work....)
   ###use the hist output to split up real vis according to X values
   ## and just use the minimum value of the vis in each X_bin
   #real_vis_data_min_in_X_bin_array = np.empty(len(X_bins)-1)
   #X_bin_centres_array = np.empty(len(X_bins)-1) 
   #for X_bin_index,n_bin in enumerate(n_X):
   #   start_value = X_bins[X_bin_index]
   #   end_value = X_bins[X_bin_index+1]
   #   X_group_array[(X_short_parallel_array >= start_value) & (X_short_parallel_array < end_value)] = X_bin_index
   #   real_vis_data_in_X_bin = real_vis_data_sorted_array[(X_short_parallel_array >= start_value) & (X_short_parallel_array < end_value)]
   #   #check if all nans
   #   if len(real_vis_data_in_X_bin) > 0:
   #      real_vis_data_min_in_X_bin = np.nanmin(real_vis_data_in_X_bin)
   #   else:
   #      real_vis_data_min_in_X_bin = np.nan
   #   #print real_vis_data_in_X_bin
   #   if (X_bin_index < len(n_X)):
   #      real_vis_data_min_in_X_bin_array[X_bin_index] = real_vis_data_min_in_X_bin
   #      X_bin_centres_array[X_bin_index] = start_value + (X_bin_width / 2.)
 
   
   #make a data frame
   #data = {'X_global':X_short_parallel_array,'Y_angular':Y_short_parallel_angular_array,'real_vis':real_vis_data_sorted_array, 'Y_angular_group':Y_angular_group_array,'X_group':X_group_array}
   data = {'X_global':X_short_parallel_array,'real_vis':real_vis_data_sorted_array}
   
 
   
   df = pd.DataFrame(data) 
   
   X_short_parallel_array = np.asarray(df['X_global']).real
   real_vis_data_sorted_array = np.asarray(df['real_vis'])
   
   ##remove points where visibility is zero
   ##print(df.sample(5))
   #df_abs = df.abs()
   #smallest_vis = df_abs.where(df_abs.real_vis < 10000)
   #missing_values_count = smallest_vis.sum()
   #print(missing_values_count)
   #
   #sys.exit()
   

   if np.nansum(np.abs(X_short_parallel_array)) > 0 and np.nansum(np.abs(real_vis_data_sorted_array) > 0):
      #random intercept model may work if you split data up into different bins for the value of X
      #Actually, should calculate the angular response for each baseline, then divide X up according to similar values of the angular response
      #https://www.statsmodels.org/dev/examples/notebooks/generated/mixed_lm_example.html
      #tilde in mixedlm mean "predicted by"
      
      ##OLS with fixed intercept at zero
      ##X_short_parallel_array = sm.add_constant(X_short_parallel_array)
      #Lets use scikit-learn instead of statsmodels (cause there are more tutorials and Kaggle prefers)
      if model_type=='OLS_fixed_intercept':
         model = sm.OLS(real_vis_data_sorted_array, X_short_parallel_array,missing='drop')
         results = model.fit()
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[0]
         t_sky_error_jy = results.bse[0]
      if model_type=='OLS_fixed_int_min_vis':
         model = sm.OLS(real_vis_data_min_in_X_bin_array, X_bin_centres_array,missing='drop')
         results = model.fit()
         ##print results.summary()
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[0]
         t_sky_error_jy = results.bse[0]
      if model_type=='OLS_with_int_min_vis':
         X_bin_centres_array = sm.add_constant(X_bin_centres_array)
         model = sm.OLS(real_vis_data_min_in_X_bin_array, X_bin_centres_array,missing='drop')
         results = model.fit()
         ##print results.summary()
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[1]
         t_sky_error_jy = results.bse[1]
      elif model_type=='OLS_fixed_int_subtr_Y':
         #subtract Y from the data before fitting (should get rid of the angular variations)
         real_vis_data_sorted_array_subtr_Y = real_vis_data_sorted_array - Y_short_parallel_angular_array_Jy
         model = sm.OLS(real_vis_data_sorted_array_subtr_Y, X_short_parallel_array,missing='drop')
         results = model.fit()
         ##print results.summary()
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[0]
         t_sky_error_jy = results.bse[0]
      elif model_type=='OLS_with_intercept':
         X_short_parallel_array = sm.add_constant(X_short_parallel_array)
         model = sm.OLS(real_vis_data_sorted_array, X_short_parallel_array,missing='drop')
         results = model.fit()
         ##print results.summary()
         parameters = results.params
         ##print parameters
         t_sky_jy = parameters[1]
         t_sky_error_jy = results.bse[1]
      elif model_type=='mixedlm':
         #test_data = sm.datasets.get_rdataset('dietox', 'geepack').data
         #md = smf.mixedlm("Weight ~ Time", data, groups=data["Pig"])
         #mdf = md.fit()
         #print(mdf.summary())
         #print(test_data.head())
         md = smf.mixedlm("real_vis ~ X_global", data, groups=data["Y_angular_group"])
         results = md.fit()
         print(results.summary()) 
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[1]
         t_sky_error_jy = results.bse[1]
      elif model_type=='WLS':
         print('using WLS')
         sns_plot = sns.regplot(X_short_parallel_array,real_vis_data_sorted_array)
         fig = sns_plot.get_figure()
         figname = 'vis_vs_X_sns_regplot'
         fig.savefig(figname) 
         print('saved %s ' % figname)
         #sns_plot_Y = sns.regplot(Y_short_parallel_angular_array,real_vis_data_sorted_array)
         #fig2 = sns_plot_Y.get_figure()
         #figname = 'vis_vs_Y_sns_regplot'
         #fig2.savefig(figname) 
         #print('saved %s ' % figname)
         #errors are inversely proportional to X, so can just use X_short_parallel_array is the weights!
         #also they are inversely proportion to Y 
         weights_x = [x**2 for x in X_short_parallel_array]
         #weights = weights_x / Y_short_parallel_angular_array
         # reshape for compatibility
         X_short_parallel_array_reshape = X_short_parallel_array.reshape(-1, 1)
         WLS = LinearRegression(fit_intercept=False)
         results = WLS.fit(X_short_parallel_array_reshape, real_vis_data_sorted_array, sample_weight=weights)
         print(WLS.intercept_, WLS.coef_)         
         t_sky_jy = WLS.coef_
         t_sky_error_jy = 0.01
      elif model_type=='OLS_global_angular':
         model = LinearRegression(fit_intercept=False)
         model.fit(df[['X_global','Y_angular']],df['real_vis'])
         #print model.coef_
         t_sky_jy = model.coef_[0]
         t_sky_error_jy = 0.01
   else:
      print("X_short_parallel_array all NaNs, returning Tsky NaN")
      return(np.nan,np.nan,np.nan,np.nan,freq_MHz_fine_chan)
      
      
   
   jy_to_K = (wavelength**2) / (2. * k * 1.0e26)   # * 6*PI? (or 4 pi and then fiddle X again?) There are also the visibility weights that I have ignored .... a factor of 120.8
   print("jy_to_K %.4E Jy" % jy_to_K)
   t_sky_K = jy_to_K * t_sky_jy
   t_sky_error_K = jy_to_K * t_sky_error_jy
   print("t_sky_K is %0.4E +/- %0.04f K" % (t_sky_K,t_sky_error_K))
   fit_string = "y=%0.1fx" % t_sky_jy         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
  
   print("diffuse_global_value is %0.4E" % diffuse_global_value) 
   
   ratio_in_out = diffuse_global_value / t_sky_K
   print("ratio between input and output T_sky is %0.4f" % ratio_in_out )
   
   y_pos = np.max(results.fittedvalues)
   x_pos = 1.2 * np.min(X_short_parallel_array)
   
    
   #get rid of nans
   real_vis_data_sorted_array_nonans = real_vis_data_sorted_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
   X_short_parallel_array_nonans = X_short_parallel_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
   
   plt.clf()
   if model_type=='OLS_with_intercept':
      plt.plot(X_short_parallel_array_nonans[:,1], real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans[:,1], results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=='OLS_with_int_min_vis':
      plt.plot(X_bin_centres_array[:,1], real_vis_data_min_in_X_bin_array,label='%s data min in X bin' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_bin_centres_array[:,1], results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=='OLS_fixed_int_min_vis':
      plt.plot(X_bin_centres_array, real_vis_data_min_in_X_bin_array,label='%s data min in X bin' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_bin_centres_array, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=="OLS_fixed_int_subtr_Y":
      real_vis_data_sorted_array_subtr_Y_nonans = real_vis_data_sorted_array_subtr_Y[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_subtr_Y_nonans,label='%s data - Y' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   else:
      plt.scatter(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   

   map_title="Data and fit" 
   plt.xlabel("Expected global-signal response")
   plt.ylabel("Real component of visibility (Jy)")
   plt.legend(loc=1)
   plt.text(x_pos, y_pos, fit_string)
   #plt.ylim([0, 3.5])
   fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name)  
   

   
   #now use the fit to identify outliers probably due to rfi
   #subtract the model from the data
   real_vis_data_sorted_array_subtr_model = real_vis_data_sorted_array_nonans - results.fittedvalues
   #take the mean 
   real_vis_data_sorted_array_subtr_model_mean = np.nanmean(real_vis_data_sorted_array_subtr_model)
   real_vis_data_sorted_array_subtr_model_std = np.nanstd(real_vis_data_sorted_array_subtr_model)

   #mask values greater than 5 sigma away from mean
   thresh = 5.* real_vis_data_sorted_array_subtr_model_std
   real_vis_data_sorted_array_flagged = np.copy(real_vis_data_sorted_array_nonans)
   real_vis_data_sorted_array_flagged[(np.abs(real_vis_data_sorted_array_subtr_model) > thresh)] = np.nan
   
   

   #get rid of nans
   #real_vis_data_sorted_array_flagged = real_vis_data_sorted_array_flagged[np.argwhere(np.logical_not(np.isnan(real_vis_data_sorted_array_flagged)))]
   #X_short_parallel_array_flagged = X_short_parallel_array_nonans[np.argwhere(np.logical_not(np.isnan(real_vis_data_sorted_array_flagged)))]
   
   if (X_short_parallel_array_nonans.shape[0]>0):
      model = sm.OLS(real_vis_data_sorted_array_flagged, X_short_parallel_array_nonans,missing='drop')
      results = model.fit()
      ##print results.summary()
      parameters = results.params
      print parameters
   
      t_sky_jy = parameters[0]
      t_sky_error_jy = results.bse[0]
      
      X_short_parallel_array_nonans_nonans = X_short_parallel_array_nonans[np.logical_not(np.isnan(real_vis_data_sorted_array_flagged))]
      
      plt.clf()
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_flagged,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans_nonans, results.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   
      map_title="Flagged data and fit" 
      plt.xlabel("Expected global-signal response")
      plt.ylabel("Real component of visibility (Jy) flagged")
      plt.legend(loc=1)
      plt.text(x_pos, y_pos, fit_string)
      #plt.ylim([0, 3.5])
      fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_flagged.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name) 
      
      #convert to K
      t_sky_K_flagged = jy_to_K * t_sky_jy
      t_sky_error_K_flagged = jy_to_K * t_sky_error_jy
      print("t_sky_K_flagged is %0.4E +/- %0.04f K" % (t_sky_K_flagged,t_sky_error_K_flagged))
      fit_string_K = "y=%0.1fx" % t_sky_K_flagged      #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
   
      #plot in K
      real_vis_data_sorted_array_flagged_K = real_vis_data_sorted_array_flagged * jy_to_K
      fitted_values_K = results.fittedvalues * jy_to_K
      y_pos_K = 0 # * np.max(fitted_values_K)
      x_pos_K = 0.05 #1.4 * np.min(X_short_parallel_array)
   
      plt.clf()
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_flagged_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans_nonans, fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')

      map_title="Flagged data and fit" 
      plt.xlabel("Global response (unity sky)")
      plt.ylabel("Visibility amplitude (simulations)")
      plt.legend(loc=1)
      plt.text(x_pos_K, y_pos_K, fit_string_K)
      #plt.ylim([0, 3.5])
      fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_flagged_K.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name) 
      
   else:
      t_sky_jy = np.nan
      t_sky_error_jy = np.nan
      t_sky_K_flagged_K = np.nan
      t_sky_error_K_flagged = np.nan
    

   
   return t_sky_K,t_sky_error_K,t_sky_K_flagged,t_sky_error_K_flagged,freq_MHz_fine_chan

def extract_data_from_eda2_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,EDA2_chan,n_obs,calculate_uniform_response=False,include_angular_info=True):
   centre_freq = float(freq_MHz_list[freq_MHz_index])
   centre_wavelength = 300./centre_freq
   fine_chan_width_MHz = fine_chan_width_Hz/1000000.    

   lst_hrs = lst_hrs_list[0]
   lst_deg = (float(lst_hrs)/24.)*360.
   
   #this needs to change for each chan or you are wasting time
   #max_baselines_included = 200 #I think 1680 is the most baselines I've seen used for the current data at lowest freq (but that was for concat... think you only need 320 max at lowest freq)
         
   #get the diffuse global diffuse value used in the simulation (from gsm)
   EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)          
   sky_averaged_diffuse_array_beam_lsts_filename = "%seda_model_%s_lst_2.00_hr_int_0.13_hr_N_D_gsm_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,pol)
   diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   #just doing one freq at a time right now for EDA2, not sure how this works with fine chans
   diffuse_global_value = diffuse_global_value_array[0]
   
   obs_time_list = EDA2_obs_time_list_each_chan[freq_MHz_index]
   
     
   #open one to get the number of fine chans
   uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_list[0])
   
   
   hdulist = fits.open(uvfits_filename)
   hdulist.info()
   info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
   #print info_string
   uvtable = hdulist[0].data
   uvtable_header = hdulist[0].header
   #print(uvtable_header)
   hdulist.close()

   visibilities_single = uvtable['DATA']
   visibilities_shape = visibilities_single.shape
   print("visibilities_shape")
   print(visibilities_shape)
   
   #need to understand this timestep stuff, for some reason EDA2 vis have more rows that expected ....
   #n_timesteps = n_vis/n_baselines
   #print "n_vis is %s, n_baselines is %s, so n_timesteps %s " % (n_vis,n_baselines,n_timesteps)
   
   #for some reason there is an extra column in the wsclean uvfits files, probably because of the way CASA exports them..
   if wsclean:
      n_fine_chans = visibilities_single.shape[4]
   else:
      n_fine_chans = visibilities_single.shape[3]

   for EDA2_obs_time in obs_time_list:
         uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
         #check here
         if not os.path.exists(uvfits_filename):
            continue
         else:
            hdulist = fits.open(uvfits_filename)
            hdulist.info()
            info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
            #print info_string
            uvtable = hdulist[0].data
            uvtable_header = hdulist[0].header
            #print(uvtable_header)
            hdulist.close()
            
            visibilities = uvtable['DATA']
   
            UU_s_array = uvtable['UU']
            UU_m_array = UU_s_array * c   
            VV_s_array = uvtable['VV']
            VV_m_array = VV_s_array * c
   
       
            #Need to sort by baseline length (then only use short baselines)
            baseline_length_array_m = np.sqrt(UU_m_array**2 + VV_m_array**2)
            
            baseline_length_array_m_inds = baseline_length_array_m.argsort()
            baseline_length_array_m_sorted_orig = baseline_length_array_m[baseline_length_array_m_inds]
            
            UU_m_array_sorted_orig = UU_m_array[baseline_length_array_m_inds]
            VV_m_array_sorted_orig = VV_m_array[baseline_length_array_m_inds]
            #real_vis_data_sorted_orig = real_vis_data[baseline_length_array_m_inds]
   
             
            #eda2 data may have bad baselines where uu=vv=0 (or are these the autos?), dont use these
            baseline_length_array_m_sorted = baseline_length_array_m_sorted_orig[UU_m_array_sorted_orig>0]
            VV_m_array_sorted = VV_m_array_sorted_orig[UU_m_array_sorted_orig>0]
            #real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
   
            
            #leave this here!!
            UU_m_array_sorted = UU_m_array_sorted_orig[UU_m_array_sorted_orig>0]
            
            baseline_length_array_lambda_sorted = baseline_length_array_m_sorted / centre_wavelength
                  
            baseline_length_array_lambda_sorted_cut = baseline_length_array_lambda_sorted[baseline_length_array_lambda_sorted < baseline_length_thresh_lambda]
                  
            n_baselines_included = len(baseline_length_array_lambda_sorted_cut)
            print("n_baselines_included %s for obs %s, EDA2 chan %s" % (n_baselines_included,EDA2_obs_time,EDA2_chan))      
   
            if calculate_uniform_response:
            
               n_pix = hp.nside2npix(NSIDE)
               print('n_pix')
               print(n_pix)
               pixel_solid_angle = (4.*np.pi) / n_pix
               #print("pixel_solid_angle is %0.4E" % pixel_solid_angle)
               hpx_index_array = np.arange(0,n_pix,1)
   
               #need to rotate this beam, since pix2ang in make_hpx_beam gives the colatitude, which is taken from the 'north pole' as healpy does not use az/alt
               #rot_theta_beam = - np.pi / 2.
               rot_theta_sky = np.pi / 2.
               rot_phi_beam = 0.
   
               T_sky_rough = 180.*(centre_freq/180.)**(-2.5)
               print("180@180 t_sky chan %s is %s" % (EDA2_chan,T_sky_rough))
                  
               #need to update this for doing each fine chan (simulate needs to be changed as well for diffuse global value, which is currently just on centre freqs)
               if include_angular_info:            
                  print("calculating time for lst_deg %0.3f " % (lst_deg))
                  #check vale
                  year=2000
                  month=1
                  day=1
                  hour=0
                  minute=0
                  second=0
                  #hour=int(np.floor(float(lst_hrs)))
                  #minute=int(np.floor((float(lst_hrs)-hour) * 60.))
                  #second=int(((float(lst_hrs)-hour) * 60. - minute) * 60.)
                  
                  date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
                  #print(date_time_string)
                  latitude, longitude, elevation = mwa_latitude_pyephem, mwa_longitude_pyephem, mwa_elevation
                  #Parkes (pygsm example)
                  #latitude, longitude, elevation = '-32.998370', '148.263659', 100
                  ov = GSMObserver()
                  ov.lon = longitude
                  ov.lat = latitude
                  ov.elev = elevation
                  
                  #print time_string
                  
                  #set an initial time using astropy Time
                  
                  astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (year, month, day, hour, minute, second)
                  time_initial = Time(astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
                  
                  #calculate the local sidereal time for the MWA at date_obs_initial
                  lst_initial = time_initial.sidereal_time('apparent')
                  
                  lst_initial_days = (lst_initial.value / 24.) * (23.9344696 /24.)
                  
                  #want lst == lst_hrs, so add time so that this is true
                  delta_time = TimeDelta(lst_initial_days, format='jd') 
                  
                  desired_lst_days = float(lst_hrs) / 24.
                  
                  time_final = time_initial - delta_time + TimeDelta(desired_lst_days, format='jd') 
                  
                  print('final LST is: ')
                  print(time_final.sidereal_time('apparent'))
                  
                  #close enough......
                  
                  time_string = time_final.utc.iso
                  #time_string = "%02d_%02d_%02d" % (hour,minute,second)
                  
                  hour = int(time_string.split(' ')[1].split(':')[0])
                  minute = int(time_string.split(' ')[1].split(':')[1])
                  minute = int(np.floor(float(time_string.split(' ')[1].split(':')[2])))
                  
                  date_obs = datetime(year, month, day, hour, minute, second)
                  ov.date = date_obs
                  
                  gsm_map = ov.generate(centre_freq)
   
               #Need to update this to do each fine chan
               gsm_map_angular = gsm_map - diffuse_global_value
               gsm_map_angular = rotate_map(gsm_map_angular, rot_theta_sky, rot_phi_beam)
               
               # no point doing this for each fine chan, can move down further (beam and phase angle barely changes, makes no diff to derived t_sky)
               #fine_chan_index_array = range(n_fine_chans)
               #for fine_chan_index in fine_chan_index_array:
               #fine_chan_index = int(fine_chan_index)
               
               ##data coming out of the TPMs is reversed by coarse chan so for 20200303_data (and 20200304), need to change the freq calculation
               #freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
               #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               
               #wavelength = 300./float(freq_MHz_fine_chan)
               
               #print("fine_chan index,MHz,wavelength")
               #print(fine_chan_index)
               #print(freq_MHz_fine_chan)
               #print(wavelength)
                  
               #beam stuff in function
               short_dipole_parallel_beam_map = make_hpx_beam(NSIDE,pol,centre_wavelength,dipole_height_m)
               #short_dipole_parallel_beam_map = make_hpx_beam(NSIDE,pol,wavelength,dipole_height_m)
      
               baseline_phi_rad_array = np.arctan(UU_m_array_sorted/VV_m_array_sorted)
               
               if pol == 'Y':
                  baseline_phi_rad_array_pure_parallel = baseline_phi_rad_array * 0. + np.pi/2.
                  baseline_phi_rad_array_pure_inline = baseline_phi_rad_array * 0. 
               else:
                  baseline_phi_rad_array_pure_parallel = baseline_phi_rad_array * 0. 
                  baseline_phi_rad_array_pure_inline = baseline_phi_rad_array * 0. + np.pi/2.      
               
               #print baseline_phi_rad_array.shape
               baseline_theta_rad_array = baseline_phi_rad_array * 0. + np.pi/2. 
               
               baseline_vector_array_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array)
               #print baseline_vector_array_unit[0,:]
               baseline_vector_array_pure_parallel_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array_pure_parallel)
               baseline_vector_array_pure_inline_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array_pure_inline)
               
               
               #Need to rotate all these vectors by -pi/2, just like the hpx beam map, since hp doesnt use alt/az, so theta=pi/2 is actually pointing at the zenith in orthographic proj
               #https://vpython.org/contents/docs/VisualIntro.html
               #rotate about x axis
               #forget about rotating the vectors its the phase angle hpx array you need to rotate!
               rot_axis = [0,1,0]
               rot_theta = 0. #np.pi / 4.
               
               baseline_vector_array_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_unit)
               
               #print baseline_vector_array_unit_rotated[0,:]
               
               
               baseline_vector_array_pure_parallel_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_pure_parallel_unit)
               baseline_vector_array_pure_inline_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_pure_inline_unit)
               
               baseline_vector_array = baseline_vector_array_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
               baseline_vector_array_pure_parallel = baseline_vector_array_pure_parallel_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
               baseline_vector_array_pure_inline = baseline_vector_array_pure_inline_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
               
               sky_vector_array_unrotated = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array)))
               #sky_vector_array_unrotated_test = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array[1])))
               #print sky_vector_array_unrotated_test
               #sys.exit()
               #do the same rotation for the sky vector
               sky_vector_array = rotate_vector(rot_axis,rot_theta,sky_vector_array_unrotated)
               
               baseline_vector_array = baseline_vector_array
               
               X_short_parallel_array = np.full(len(baseline_vector_array),np.nan,dtype=complex)
               Y_short_parallel_angular_array = np.full(len(baseline_vector_array),0,dtype=complex)
               
               X_short_parallel_array_pure_parallel = np.full(len(baseline_vector_array),np.nan,dtype=complex)
               X_short_parallel_array_pure_inline = np.full(len(baseline_vector_array),np.nan,dtype=complex)
            
               #need do this bit just once for each chan
               for baseline_vector_index in range(0,n_baselines_included):
                  #Just try doing the integral (sum) all in one go with one baseline
                  #baseline_vector_test = baseline_vector_array[0]
                  
                  baseline_vector_for_dot_array = baseline_vector_array[baseline_vector_index,:]
                  
                  baseline_vector_for_dot_array_pure_parallel = baseline_vector_array_pure_parallel[baseline_vector_index,:]
                  baseline_vector_for_dot_array_pure_inline = baseline_vector_array_pure_inline[baseline_vector_index,:]
                  
                  b_dot_r_array = (baseline_vector_for_dot_array * sky_vector_array).sum(axis=1)
               
                  b_dot_r_array_pure_parallel = (baseline_vector_for_dot_array_pure_parallel * sky_vector_array).sum(axis=1)
                  b_dot_r_array_pure_inline = (baseline_vector_for_dot_array_pure_inline * sky_vector_array).sum(axis=1)
                  
                  
                  phase_angle_array = 2.*np.pi*b_dot_r_array/centre_wavelength
                  
                  phase_angle_array_pure_parallel = 2.*np.pi*b_dot_r_array_pure_parallel/centre_wavelength
                  phase_angle_array_pure_inline = 2.*np.pi*b_dot_r_array_pure_inline/centre_wavelength
                  
                  element_short_parallel_array = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array)
                  
                  element_short_parallel_array_pure_parallel = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array_pure_parallel)
                  element_short_parallel_array_pure_inline = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array_pure_inline)
                  
                  #for angular info
                  if include_angular_info:
                     element_short_parallel_angular_array = short_dipole_parallel_beam_map * gsm_map_angular * np.exp(-1j*phase_angle_array)
      
                  X_short_parallel =  np.sum(element_short_parallel_array) * pixel_solid_angle # (4.*np.pi/float(n_pix))
               
                  X_short_parallel_pure_parallel =  np.sum(element_short_parallel_array_pure_parallel) * pixel_solid_angle # (4.*np.pi/float(n_pix))
                  X_short_parallel_pure_inline =  np.sum(element_short_parallel_array_pure_inline) * pixel_solid_angle # (4.*np.pi/float(n_pix))
               
                  if include_angular_info:
                     Y_short_parallel_angular =  np.sum(element_short_parallel_angular_array) * pixel_solid_angle
      
                  X_short_parallel_array[baseline_vector_index] = X_short_parallel
                  
                  X_short_parallel_array_pure_parallel[baseline_vector_index] = X_short_parallel_pure_parallel
                  X_short_parallel_array_pure_inline[baseline_vector_index] = X_short_parallel_pure_inline
                  
                  
                  if include_angular_info:
                     Y_short_parallel_angular_array[baseline_vector_index] = Y_short_parallel_angular
                  
                  #now for each fine chan work out the frequency and wavelength and do the baseline cutoff and save stuff
                  #save for each obs separately
                  
               ##This stuff got moved up to do new beams and new phase angle for each fine chan, but made no diff to result
               fine_chan_index_array = range(n_fine_chans)
               for fine_chan_index in fine_chan_index_array:
                  fine_chan_index = int(fine_chan_index)
                  
                  #data coming out of the TPMs is reversed by coarse chan so for 20200303_data (and 20200304), need to change the freq calculation
                  freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
                  #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
                  
                  #wavelength = 300./float(freq_MHz_fine_chan)
                  
                  print("fine_chan index,MHz,wavelength")
                  print(fine_chan_index)
                  print(freq_MHz_fine_chan)
                  #print(wavelength)
                  
                  if wsclean:
                     real_vis_data = visibilities[:,0,0,0,fine_chan_index,0,0]
                  else:
                     real_vis_data = visibilities[:,0,0,fine_chan_index,0,0]
                  
                  ##UU_m_array_sorted_orig = UU_m_array[baseline_length_array_m_inds]
                  ##VV_m_array_sorted_orig = VV_m_array[baseline_length_array_m_inds]
                  real_vis_data_sorted_orig = real_vis_data[baseline_length_array_m_inds]
   
                  #eda2 data may have bad baselines where uu=vv=0 (or are these the autos?), dont use these
                  ##baseline_length_array_m_sorted = baseline_length_array_m_sorted_orig[UU_m_array_sorted_orig>0]
                  ##VV_m_array_sorted = VV_m_array_sorted_orig[UU_m_array_sorted_orig>0]
                  real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
                  
                  #leave this here!!
                  ##UU_m_array_sorted = UU_m_array_sorted_orig[UU_m_array_sorted_orig>0]
                  
                  #EDA2 data may also have visibilities where the cal solutions are zero, jump to here
                  #write out ALL calibrated vis as uvfits, then check n_vis and n_timesteps for each calibrated uvfits
                   

                  baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol_%s.npy" % (freq_MHz_fine_chan,pol,EDA2_obs_time)
                  np.save(baseline_length_array_lambda_sorted_cut_filename,baseline_length_array_lambda_sorted_cut)
                  print("saved %s" % baseline_length_array_lambda_sorted_cut_filename)
                  
                  X_short_parallel_array_filename = "X_uniform_resp_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
                  np.save(X_short_parallel_array_filename,X_short_parallel_array[0:n_baselines_included])
               
                  print("saved %s" % X_short_parallel_array_filename)
               
                  X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
                  np.save(X_short_parallel_array_filename_pure_inline,X_short_parallel_array_pure_inline[0:n_baselines_included])
                  print("saved %s" % X_short_parallel_array_filename_pure_inline)
                  
                  X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
                  np.save(X_short_parallel_array_filename_pure_parallel,X_short_parallel_array_pure_parallel[0:n_baselines_included])
                  print("saved %s" % X_short_parallel_array_filename_pure_parallel)
                  
                  real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
                  np.save(real_vis_data_sorted_array_filename,real_vis_data_sorted[0:n_baselines_included])
                  print("saved %s" % real_vis_data_sorted_array_filename)
                  
                  #update for fine chans
                  if include_angular_info:
                     Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_chan_%s_%0.3f_MHz_%s_pol_%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,EDA2_obs_time)
                     np.save(Y_short_parallel_angular_array_filename,Y_short_parallel_angular_array[0:n_baselines_included])
                     print("saved %s" % Y_short_parallel_angular_array_filename)
      
              
   return(diffuse_global_value)

def solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,include_angular_info=False,EDA2_data=False, EDA2_obs_time='None',EDA2_chan='None',n_obs_concat=1,wsclean=False,fast=False,calculate_uniform_response=True):
   freq_MHz = freq_MHz_list[freq_MHz_index]
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
   if 'global_unity' in signal_type_list:
       signal_type_postfix += '_GU'
       concat_output_name_base += '_GU'
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse_global at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DG_%s' % sky_model
          concat_output_name_base += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse_angular at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DA_%s' % sky_model
          concat_output_name_base += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print("cant have global and global_EDGES in signal_type_list")
          sys.exit()
       else:
          signal_type_postfix += '_G' 
          concat_output_name_base += '_G' 
   if 'global_EDGES' in signal_type_list:
       signal_type_postfix += '_ED' 
       concat_output_name_base += '_ED'
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'
   
   #if EDA2_data==True:
   #   signal_type_postfix = "_EDA2_data"
   #n_baselines_included
   #n_baselines_included = 20
   
   
   n_lsts = len(lst_hrs_list)
   
   #for EDA2
   n_ants = 256

   n_pix = hp.nside2npix(NSIDE)
   print('n_pix')
   print(n_pix)
   pixel_solid_angle = (4.*np.pi) / n_pix
   #print("pixel_solid_angle is %0.4E" % pixel_solid_angle)
   hpx_index_array = np.arange(0,n_pix,1)
   
   #need to understand this timestep stuff, for some reason EDA2 vis have more rows that expected ....
   #eda2 data includes autos?
   if EDA2_data:
      n_baselines = n_ants*(n_ants-1) / 2. + 256
   else:
      n_baselines = n_ants*(n_ants-1) / 2.
   
   
   #open one uvfits file to get n_timesteps
   lst_hrs = lst_hrs_list[0]
   lst_deg = (float(lst_hrs)/24.)*360.
   if EDA2_data:
      if fast:
         print('doing fast')   
         obs_time_list = EDA2_obs_time_list_each_chan[freq_MHz_index]
         #print(obs_time_list)
         
         #open one to get the number of fine chans
         uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_list[0])
         
         #read the cal uvfits, extract real vis uu and vv
         print("%s" % uvfits_filename)
         hdulist = fits.open(uvfits_filename)
         hdulist.info()
         info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
         #print info_string
         uvtable = hdulist[0].data
         uvtable_header = hdulist[0].header
         #print(uvtable_header)
         hdulist.close()

         visibilities_single = uvtable['DATA']
         visibilities_shape = visibilities_single.shape
         print("visibilities_shape")
         print(visibilities_shape)
         
         #need to understand this timestep stuff, for some reason EDA2 vis have more rows that expected ....
         #n_timesteps = n_vis/n_baselines
         #print "n_vis is %s, n_baselines is %s, so n_timesteps %s " % (n_vis,n_baselines,n_timesteps)
         
         #for some reason there is an extra column in the wsclean uvfits files, probably because of the way CASA exports them..
         if wsclean:
            n_fine_chans = visibilities_single.shape[4]
         else:
            n_fine_chans = visibilities_single.shape[3]
                  
         if EDA2_data:
            #print("EDA2 data. Omitting 1 edge chan each side, %s chans present, %s chans used" % (n_fine_chans,n_fine_chans-2))
            #fine_chan_index_array = range(n_fine_chans)[1:n_fine_chans-1]
            fine_chan_index_array = range(n_fine_chans)
            #print(fine_chan_index_array)
         else:
            fine_chan_index_array = np.asarray([0])
         centre_freq = float(freq_MHz)
         fine_chan_width_MHz = fine_chan_width_Hz/1000000.         
         
         
         #now do for each fine chan:
         for fine_chan_index in fine_chan_index_array:
            fine_chan_index = int(fine_chan_index)
            baseline_length_array_lambda_sorted_cut_list = []

            if EDA2_data:
               #data coming out of the TPMs is reversed by coarse chan so for 20200303_data (and 20200304), need to change the freq calculation
               freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
               #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               #freq_MHz_fine_chan = centre_freq - (fine_chan_index)*fine_chan_width_MHz
            else:
               freq_MHz_fine_chan = freq_MHz
            wavelength = 300./float(freq_MHz_fine_chan)
            
            print("fine_chan index,MHz,wavelength")
            print(fine_chan_index)
            print(freq_MHz_fine_chan)
            print(wavelength)

            unity_vis_data_sorted_list = []
            baseline_length_array_lambda_sorted_cut_list = []
            real_vis_data_sorted_list = []
            for obs_time_fast in obs_time_list:
               uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_fast)
               unity_uvfits_filename = "%s/unity_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_fast)
               
               #read the cal uvfits, extract real vis uu and vv
               print("%s" % uvfits_filename)
               hdulist = fits.open(uvfits_filename)
               hdulist.info()
               info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
               #print info_string
               uvtable = hdulist[0].data
               uvtable_header = hdulist[0].header
               #print(uvtable_header)
               hdulist.close()
   
               visibilities_single = uvtable['DATA']
               visibilities_shape = visibilities_single.shape
               print("visibilities_shape")
               print(visibilities_shape)
    
               if wsclean:
                  real_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,0]
                  imag_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,1]
                  weights_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,2]
               else:
                  real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
                  imag_vis_data = visibilities_single[:,0,0,fine_chan_index,0,1]
                  weights_vis_data = visibilities_single[:,0,0,fine_chan_index,0,2]
      
               UU_s_array = uvtable['UU']
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV']
               VV_m_array = VV_s_array * c
      
               #print(VV_s_array[0:100])
               print(UU_s_array[0:100])

               #Need to sort by baseline length (then only use short baselines)
               baseline_length_array_m = np.sqrt(UU_m_array**2 + VV_m_array**2)
               
               baseline_length_array_m_inds = baseline_length_array_m.argsort()
               baseline_length_array_m_sorted_orig = baseline_length_array_m[baseline_length_array_m_inds]
               
               #use only baselines shorter than threshold
               
            
               UU_m_array_sorted_orig = UU_m_array[baseline_length_array_m_inds]
               VV_m_array_sorted_orig = VV_m_array[baseline_length_array_m_inds]
               real_vis_data_sorted_orig = real_vis_data[baseline_length_array_m_inds]
               imag_vis_data_sorted_orig = imag_vis_data[baseline_length_array_m_inds]
               weights_vis_data_sorted_orig = weights_vis_data[baseline_length_array_m_inds]
               

               #eda2 data may have bad baselines where uu=vv=0 (or are these the autos?), dont use these
               baseline_length_array_m_sorted = baseline_length_array_m_sorted_orig[UU_m_array_sorted_orig>0]
               VV_m_array_sorted = VV_m_array_sorted_orig[UU_m_array_sorted_orig>0]
               real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
               imag_vis_data_sorted = imag_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
               weights_vis_data_sorted = weights_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
               
               #leave this here!!
               UU_m_array_sorted = UU_m_array_sorted_orig[UU_m_array_sorted_orig>0]
               
               #EDA2 data may also have visibilities where the cal solutions are zero, jump to here
               #write out ALL calibrated vis as uvfits, then check n_vis and n_timesteps for each calibrated uvfits
                
               baseline_length_array_lambda_sorted = baseline_length_array_m_sorted / wavelength
               
               
               baseline_length_array_lambda_sorted_cut = baseline_length_array_lambda_sorted[baseline_length_array_lambda_sorted < baseline_length_thresh_lambda]
               

               n_baselines_included = len(baseline_length_array_lambda_sorted_cut)
               print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included,obs_time_fast,fine_chan_index))
               
               baseline_length_array_lambda_sorted_cut_list.append(baseline_length_array_lambda_sorted_cut)
               real_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included])
            
            
               

                  
               #######################################################################################################
               #now repeat for unity sky to get X_short_parallel!
               print("%s" % unity_uvfits_filename)
               #TEST!
               #unity_uvfits_filename = '/md0/EoR/ASSASSIN/solve_for_tsky_weighted/global_unity/eda_model_LST_030_X_50_MHz_GU.uvfits'
               hdulist = fits.open(unity_uvfits_filename)
               hdulist.info()
               info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
               #print info_string
               uvtable = hdulist[0].data
               uvtable_header = hdulist[0].header
               #print(uvtable_header)
               hdulist.close()
      
               visibilities_single = uvtable['DATA']
               visibilities_shape = visibilities_single.shape
               print("visibilities_shape")
               print(visibilities_shape)
               
               UU_s_array = uvtable['UU']
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV']
               VV_m_array = VV_s_array * c
      

               #the uvfits files used to make the unity sky data had reversed fine channel ordering (true for 20200303 and 20200304) - marcin will fix this in later data
               #yes but this should not affect the unity uvfits
               ####fine_chan_index_array = fine_chan_index_array[::-1]
            

               #real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
               #TEST!
               real_vis_data = visibilities_single[:,0,0,0,0,0]
               
               print(real_vis_data)
               sys.exit()  
                            
               #Need to sort by baseline length (then only use short baselines)
               baseline_length_array_m = np.sqrt(UU_m_array**2 + VV_m_array**2)
               
               baseline_length_array_m_inds = baseline_length_array_m.argsort()
               baseline_length_array_m_sorted_orig = baseline_length_array_m[baseline_length_array_m_inds]
               
               #use only baselines shorter than threshold
               
            
               UU_m_array_sorted_orig = UU_m_array[baseline_length_array_m_inds]
               VV_m_array_sorted_orig = VV_m_array[baseline_length_array_m_inds]
               real_vis_data_sorted_orig = real_vis_data[baseline_length_array_m_inds]

               

               #eda2 data may have bad baselines where uu=vv=0 (or are these the autos?), dont use these
               baseline_length_array_m_sorted = baseline_length_array_m_sorted_orig[UU_m_array_sorted_orig>0]
               VV_m_array_sorted = VV_m_array_sorted_orig[UU_m_array_sorted_orig>0]
               real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]

               
               #leave this here!!
               UU_m_array_sorted = UU_m_array_sorted_orig[UU_m_array_sorted_orig>0]
               
               #EDA2 data may also have visibilities where the cal solutions are zero, jump to here
               #write out ALL calibrated vis as uvfits, then check n_vis and n_timesteps for each calibrated uvfits
                
               baseline_length_array_lambda_sorted = baseline_length_array_m_sorted / wavelength
               
               
               baseline_length_array_lambda_sorted_cut = baseline_length_array_lambda_sorted[baseline_length_array_lambda_sorted < baseline_length_thresh_lambda]
               

               n_baselines_included = len(baseline_length_array_lambda_sorted_cut)
               print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included,obs_time_fast,fine_chan_index))
               
        

               real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
 

               unity_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included])

            
            #uvfits_real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix,obs_time_fast)
            #uvfits_real_vis_data_sorted = np.load(uvfits_real_vis_data_sorted_array_filename)
            #print("loaded %s" % uvfits_real_vis_data_sorted_array_filename)
            
            
            #baseline_length_array_lambda_sorted_cut_filename_unity = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_unity_%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix,obs_time_fast)
            #np.save(baseline_length_array_lambda_sorted_cut_filename_unity,baseline_length_array_lambda_sorted_cut)
            #print("saved %s" % baseline_length_array_lambda_sorted_cut_filename_unity)
            
            unity_vis_data_sorted_array = np.asarray(unity_vis_data_sorted_list)
            baseline_length_array_lambda_sorted_cut_array = np.asarray(baseline_length_array_lambda_sorted_cut_list)
            real_vis_data_sorted_array = np.asarray(real_vis_data_sorted_list)
            
            
            
            unity_vis_data_sorted_array_filename = "unity_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
            np.save(unity_vis_data_sorted_array_filename,unity_vis_data_sorted_array)
            print("saved %s" % unity_vis_data_sorted_array_filename)         

            baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
            np.save(baseline_length_array_lambda_sorted_cut_filename,baseline_length_array_lambda_sorted_cut_array)
            print("saved %s" % baseline_length_array_lambda_sorted_cut_filename)
   
            real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
            np.save(real_vis_data_sorted_array_filename,real_vis_data_sorted_array)
            print("saved %s" % real_vis_data_sorted_array_filename)
           
      else:
         pass
      if n_obs_concat==1:
         if wsclean==True:
            uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
         else:
            uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
      else:
         if wsclean==True:
            uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
         else:
            #concat_chan_64_20191202T171525_n_obs_13.uvfits
            #uvfits_filename = "%s/av_chan_%s_%s_n_obs_%s_t_av_cal_freq_av.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
   else:
      uvfits_filename = "%s_LST_%03d_%s_%2d_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
   
   hdulist = fits.open(uvfits_filename)
   uvtable = hdulist[0].data
   visibilities = uvtable['DATA']
   n_vis = visibilities.shape[0]
   
   #need to understand this timestep stuff, for some reason EDA2 vis have more rows that expected ....
   n_timesteps = n_vis/n_baselines
   #print "n_vis is %s, n_baselines is %s, so n_timesteps %s " % (n_vis,n_baselines,n_timesteps)
   
   #for some reason there is an extra column in the wsclean uvfits files, probably because of the way CASA exports them..
   if wsclean:
      n_fine_chans = visibilities.shape[4]
   else:
      n_fine_chans = visibilities.shape[3]
   
   
   
   uvfits_filename_list = []
   if EDA2_data:
      if n_obs_concat==1:
         if wsclean==True:
            uvfits_filename = "%s/chan_%s_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            ms_filename = "%s/chan_%s_%s_ws_cal.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
         else:
            uvfits_filename = "%s/chan_%s_%s_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
      else:
         if wsclean==True:
            uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat) 
            ms_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat) 
         else:
            #uvfits_filename = "%s/av_chan_%s_%s_n_obs_%s_t_av_cal_freq_av.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)      
      uvfits_filename_list = [uvfits_filename]

   else:
      for lst_hrs in lst_hrs_list:
         lst_deg = (float(lst_hrs)/24.)*360.
         uvfits_filename = "%s_LST_%03d_%s_%2d_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
         uvfits_filename_list.append(uvfits_filename)
   
   #1. Get the u,v and visibilities for each fine chan 
   print("Solving for Tsky from:") 
   for uvfits_filename_index,uvfits_filename in enumerate(uvfits_filename_list):
      print("%s" % uvfits_filename)
      hdulist = fits.open(uvfits_filename)
      hdulist.info()
      info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print info_string
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      #print(uvtable_header)
      hdulist.close()
      #freq_hz = float(uvtable_header['crval4'])
      #freq_MHz = freq_hz/1000000.
      #print "freq is %0.1f " % freq_MHz
      
   
      visibilities_single = uvtable['DATA']
      visibilities_shape = visibilities_single.shape
      print("visibilities_shape")
      print(visibilities_shape)
      
      if EDA2_data:
         #print("EDA2 data. Omitting 1 edge chan each side, %s chans present, %s chans used" % (n_fine_chans,n_fine_chans-2))
         #fine_chan_index_array = range(n_fine_chans)[1:n_fine_chans-1]
         fine_chan_index_array = range(n_fine_chans)
         #print(fine_chan_index_array)
      else:
         fine_chan_index_array = np.asarray([0])
      centre_freq = float(freq_MHz)
      fine_chan_width_MHz = fine_chan_width_Hz/1000000.
      for fine_chan_index in fine_chan_index_array:
         fine_chan_index = int(fine_chan_index)
         if EDA2_data:
            #for 20200303 and 20200304 data fine chan order is reversed
            freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
            #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz
            #freq_MHz_fine_chan = centre_freq - (fine_chan_index)*fine_chan_width_MHz
         else:
            freq_MHz_fine_chan = freq_MHz
         wavelength = 300./float(freq_MHz_fine_chan)
         
         print("fine_chan index,MHz,wavelength")
         print(fine_chan_index)
         print(freq_MHz_fine_chan)
         print(wavelength)
         
         #dont predefine the size of these arrays, just take whatever comes in the concat uvfits files
         #real_vis_data_unweighted_single = visibilities_single[:,0,0,fine_chan_index,0,0]
         #imag_vis_data_single = visibilities_single[:,0,0,fine_chan_index,0,1]
         #weights_vis_data_single = visibilities_single[:,0,0,fine_chan_index,0,2]
         
         
         if wsclean:
            real_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,0]
            imag_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,1]
            weights_vis_data = visibilities_single[:,0,0,0,fine_chan_index,0,2]
         else:
            real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
            imag_vis_data = visibilities_single[:,0,0,fine_chan_index,0,1]
            weights_vis_data = visibilities_single[:,0,0,fine_chan_index,0,2]
         
         
         
         #real_vis_data_single = real_vis_data_unweighted_single #* weights_vis_data
         
         #print "real visibilities shape"
         #print real_vis_data.shape
   
         
         #get the UU and VV 
         #UU_s_array_single = uvtable['UU']
         #UU_m_array_single = UU_s_array_single * c   
         #VV_s_array_single = uvtable['VV']
         #VV_m_array_single = VV_s_array_single * c

         UU_s_array = uvtable['UU']
         UU_m_array = UU_s_array * c   
         VV_s_array = uvtable['VV']
         VV_m_array = VV_s_array * c

         
         
         #forget all this timestep stuff, just use all the data
         #start_index = int(uvfits_filename_index * n_baselines * n_timesteps)
         #end_index = int(start_index + (n_baselines * n_timesteps))
         #real_vis_data[start_index:end_index] = real_vis_data_single
         #imag_vis_data[start_index:end_index] = imag_vis_data_single
         #weights_vis_data[start_index:end_index] = weights_vis_data_single
         #UU_m_array[start_index:end_index] = UU_m_array_single
         #VV_m_array[start_index:end_index] = VV_m_array_single
      
      

         #Need to sort by baseline length (then only use short baselines)
         baseline_length_array_m = np.sqrt(UU_m_array**2 + VV_m_array**2)
         
         baseline_length_array_m_inds = baseline_length_array_m.argsort()
         baseline_length_array_m_sorted_orig = baseline_length_array_m[baseline_length_array_m_inds]
         
         #use only baselines shorter than threshold
         
      
         UU_m_array_sorted_orig = UU_m_array[baseline_length_array_m_inds]
         VV_m_array_sorted_orig = VV_m_array[baseline_length_array_m_inds]
         real_vis_data_sorted_orig = real_vis_data[baseline_length_array_m_inds]
         imag_vis_data_sorted_orig = imag_vis_data[baseline_length_array_m_inds]
         weights_vis_data_sorted_orig = weights_vis_data[baseline_length_array_m_inds]
         

            
         #eda2 data may have bad baselines where uu=vv=0 (or are these the autos?), dont use these
         baseline_length_array_m_sorted = baseline_length_array_m_sorted_orig[UU_m_array_sorted_orig>0]
         VV_m_array_sorted = VV_m_array_sorted_orig[UU_m_array_sorted_orig>0]
         real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
         imag_vis_data_sorted = imag_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
         weights_vis_data_sorted = weights_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
         
         #leave this here!!
         UU_m_array_sorted = UU_m_array_sorted_orig[UU_m_array_sorted_orig>0]
         
         #EDA2 data may also have visibilities where the cal solutions are zero, jump to here
         #write out ALL calibrated vis as uvfits, then check n_vis and n_timesteps for each calibrated uvfits
          
         baseline_length_array_lambda_sorted = baseline_length_array_m_sorted / wavelength
         
         baseline_length_array_lambda_sorted_cut = baseline_length_array_lambda_sorted[baseline_length_array_lambda_sorted < baseline_length_thresh_lambda]
        
         n_baselines_included = len(baseline_length_array_lambda_sorted_cut)
         print("n_baselines_included %s fine chan %s" % (n_baselines_included,fine_chan_index))
         
         if not fast:
            #save the baseline length array sorted
            baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
            np.save(baseline_length_array_lambda_sorted_cut_filename,baseline_length_array_lambda_sorted_cut)
            print("saved %s" % baseline_length_array_lambda_sorted_cut_filename)
         
         #some sort of beam solid angle ....
         #real_vis_data_sorted_kelvin = (wavelength**2 / (2.*k) ) * real_vis_data_sorted
          
         #sys.exit()
         #need dummy hpx map to get vectors? - nope
         #sky_map =  np.full(hp.nside2npix(NSIDE),0.)
         
         #iso_beam_map = np.full(hp.nside2npix(NSIDE),1.)
         #for hpx_index,beam_value in enumerate(short_dipole_parallel_beam_map):
         
         
         #beam stuff in function
         short_dipole_parallel_beam_map = make_hpx_beam(NSIDE,pol,wavelength,dipole_height_m)
         
         #need to rotate this beam, since pix2ang in make_hpx_beam gives the colatitude, which is taken from the 'north pole' as healpy does not use az/alt
         #rot_theta_beam = - np.pi / 2.
         rot_theta_sky = np.pi / 2.
         rot_phi_beam = 0.
         
         #Dont rotate the beam map and the phase angle map, rotate the sky only (quicker!)
         #short_dipole_parallel_beam_map = rotate_map(short_dipole_parallel_beam_map, rot_theta_beam, rot_phi_beam)
         
         #I think this is actually the wrong way round, should be UU/VV maybe ... to get the X/Y pols right.
         #TRY it!
         #baseline_phi_rad_array = np.arctan(VV_m_array_sorted/UU_m_array_sorted)
         baseline_phi_rad_array = np.arctan(UU_m_array_sorted/VV_m_array_sorted)
         
         if pol == 'Y':
            baseline_phi_rad_array_pure_parallel = baseline_phi_rad_array * 0. + np.pi/2.
            baseline_phi_rad_array_pure_inline = baseline_phi_rad_array * 0. 
         else:
            baseline_phi_rad_array_pure_parallel = baseline_phi_rad_array * 0. 
            baseline_phi_rad_array_pure_inline = baseline_phi_rad_array * 0. + np.pi/2.      
      
         #print baseline_phi_rad_array[0:10]
         #print baseline_phi_rad_array.shape
         baseline_theta_rad_array = baseline_phi_rad_array * 0. + np.pi/2. 
         
         baseline_vector_array_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array)
         #print baseline_vector_array_unit[0,:]
         baseline_vector_array_pure_parallel_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array_pure_parallel)
         baseline_vector_array_pure_inline_unit=hp.ang2vec(baseline_theta_rad_array,baseline_phi_rad_array_pure_inline)
         
         
         #Need to rotate all these vectors by -pi/2, just like the hpx beam map, since hp doesnt use alt/az, so theta=pi/2 is actually pointing at the zenith in orthographic proj
         #https://vpython.org/contents/docs/VisualIntro.html
         #rotate about x axis
         #forget about rotating the vectors its the phase angle hpx array you need to rotate!
         rot_axis = [0,1,0]
         rot_theta = 0. #np.pi / 4.
         
         baseline_vector_array_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_unit)
         
         #print baseline_vector_array_unit_rotated[0,:]
         
         
         baseline_vector_array_pure_parallel_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_pure_parallel_unit)
         baseline_vector_array_pure_inline_unit_rotated = rotate_vector(rot_axis,rot_theta,baseline_vector_array_pure_inline_unit)
         
         baseline_vector_array = baseline_vector_array_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
         baseline_vector_array_pure_parallel = baseline_vector_array_pure_parallel_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
         baseline_vector_array_pure_inline = baseline_vector_array_pure_inline_unit_rotated * np.transpose([baseline_length_array_m_sorted,baseline_length_array_m_sorted,baseline_length_array_m_sorted])
         
         sky_vector_array_unrotated = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array)))
         #sky_vector_array_unrotated_test = np.transpose(np.asarray(hp.pix2vec(NSIDE,hpx_index_array[1])))
         #print sky_vector_array_unrotated_test
         #sys.exit()
         #do the same rotation for the sky vector
         sky_vector_array = rotate_vector(rot_axis,rot_theta,sky_vector_array_unrotated)
         
         baseline_vector_array = baseline_vector_array[0:n_baselines_included]
         X_short_parallel_array = np.empty(len(baseline_vector_array),dtype=complex)
         Y_short_parallel_angular_array = np.empty(len(baseline_vector_array),dtype=complex)
         
         X_short_parallel_array_pure_parallel = np.empty(len(baseline_vector_array),dtype=complex)
         X_short_parallel_array_pure_inline = np.empty(len(baseline_vector_array),dtype=complex)
      
         #print baseline_length_array_lambda_sorted.shape
         #print real_vis_data_sorted.shape
         #plot amp vs uvdistance for simulated data
         plt.clf()
         plt.scatter(baseline_length_array_lambda_sorted,real_vis_data_sorted,s=1,label='real vis data')
         #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="real vis vs uvdistance" 
         plt.xlabel("uv-distance (lambda)")
         plt.ylabel("real vis")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name) 
         
         #First work out what things should look like
         #need to implement this for each lst - for now just use the last one (i.e. dont modify lst_deg)
         
         #Check value is reasonable
         #180 at 180 rule
         T_sky_rough = 180.*(freq_MHz_fine_chan/180.)**(-2.5)
         print("180@180 t_sky is %s" % T_sky_rough)
         
                  
         print("checking global value for lst_deg %s freq_MHz %s" % (lst_deg,freq_MHz))
         #check vale
         year=2000
         month=1
         day=1
         hour=0
         minute=0
         second=0
         #hour=int(np.floor(float(lst_hrs)))
         #minute=int(np.floor((float(lst_hrs)-hour) * 60.))
         #second=int(((float(lst_hrs)-hour) * 60. - minute) * 60.)
         
         date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
         #print(date_time_string)
         latitude, longitude, elevation = mwa_latitude_pyephem, mwa_longitude_pyephem, mwa_elevation
         #Parkes (pygsm example)
         #latitude, longitude, elevation = '-32.998370', '148.263659', 100
         ov = GSMObserver()
         ov.lon = longitude
         ov.lat = latitude
         ov.elev = elevation
         
         #print time_string
         
         #set an initial time using astropy Time
         
         astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (year, month, day, hour, minute, second)
         time_initial = Time(astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
         
         #calculate the local sidereal time for the MWA at date_obs_initial
         lst_initial = time_initial.sidereal_time('apparent')
         
         lst_initial_days = (lst_initial.value / 24.) * (23.9344696 /24.)
         
         #want lst == lst_hrs, so add time so that this is true
         delta_time = TimeDelta(lst_initial_days, format='jd') 
         
         desired_lst_days = float(lst_hrs) / 24.
         
         time_final = time_initial - delta_time + TimeDelta(desired_lst_days, format='jd') 
         
         print('final LST is: ')
         print(time_final.sidereal_time('apparent'))
         
         #close enough......
         
         time_string = time_final.utc.iso
         #time_string = "%02d_%02d_%02d" % (hour,minute,second)
         
         hour = int(time_string.split(' ')[1].split(':')[0])
         minute = int(time_string.split(' ')[1].split(':')[1])
         minute = int(np.floor(float(time_string.split(' ')[1].split(':')[2])))
         
         date_obs = datetime(year, month, day, hour, minute, second)
         ov.date = date_obs
         
         
         gsm_map = ov.generate(freq_MHz) * 0.
         
         #get the diffuse global diffuse value used in the simulation (from gsm)
         if EDA2_data==True:
            EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)          
            sky_averaged_diffuse_array_beam_lsts_filename = "%s%s_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,concat_output_name_base)       
         else:
            sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % (concat_output_name_base)
         #sky_averaged_diffuse_array_no_beam_lsts_filename = "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base
         freq_MHz_index = int(freq_MHz - 50)
         diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
         if EDA2_data==True:
            #just doing one freq at a time right now for EDA2, not sure how this works with fine chans
            diffuse_global_value = diffuse_global_value_array[0]
         else:
            diffuse_global_value = diffuse_global_value_array[freq_MHz_index]
         
         #ov1 = GlobalSkyModel()
         
         if 'global' in signal_type_list:
            if 'global_EDGES' in signal_type_list:
               print("cant have global and global edges in signal_type_list")
               sys.exit()
            else:
               global_signal_value = s_21_array[freq_MHz_index]
               gsm_map += global_signal_value
         if 'global_EDGES' in signal_type_list:
            if 'global' in signal_type_list:
               print("cant have global and global edges in signal_type_list")
               sys.exit()
            else:
               global_signal_value = s_21_array_EDGES[freq_MHz_index]
               gsm_map += global_signal_value
         if 'diffuse' in signal_type_list:
            gsm_map += ov.generate(freq_MHz)
         if 'diffuse_global' in signal_type_list:
            if 'diffuse' in signal_type_list:
               print("can't have diffuse and diffuse global at the same time.")
               sys.exit()
            else:
               gsm_map += diffuse_global_value
         if 'diffuse_angular' in signal_type_list:
            if 'diffuse' in signal_type_list:
               print("can't have diffuse and diffuse angular at the same time.")
               sys.exit()
            else:
               gsm_map += ov.generate(freq_MHz)
               gsm_map -= diffuse_global_value
               #print("diffuse_angular not implemented yet")
               #sys.exit()
         
         gsm_map = rotate_map(gsm_map, rot_theta_sky, rot_phi_beam)
         
      
      
         if include_angular_info:
            #Need to update this to do each fine chan
            gsm_map_angular = ov.generate(freq_MHz) - diffuse_global_value
            gsm_map_angular = rotate_map(gsm_map_angular, rot_theta_sky, rot_phi_beam)
            #gsm_map_angular = ov.generate(freq_MHz) - np.mean(ov.generate(freq_MHz))   #this is probly wrong
            ###gsm_map_angular_abs_max = np.max(np.abs(gsm_map_angular))
            ###gsm_map_angular_norm = gsm_map_angular / gsm_map_angular_abs_max
         
         #show the sky maps multiplied by the beam
         #gsm_map[theta_array > np.pi/2.]=np.nan
         
         #av_sky_no_beam = np.nanmean(gsm_map)
         #print("av_sky_no_beam %0.4E" % av_sky_no_beam)
         
         #dont need to rotate gsm - just had the time of observation wrong!
         #fixed it to work out time to have correct(ish) LST, still out by 2.5 mins for some reason...
            
         #sky_with_beam = gsm_map * short_dipole_parallel_beam_map
         sky_with_beam = gsm_map * short_dipole_parallel_beam_map
         #print(sky_with_beam)
         
         #Close! But you need the beam in celestial coords!!!!
         
         
         
         #make an image of this
         if np.isclose(freq_MHz, 70):
            #sky with beam
            plt.clf()
            map_title="GSM from MWA at %s hrs LST %0.3f MHz" % (lst_hrs,freq_MHz)
            ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=sky_with_beam,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
            #hp.mollview(map=sky_with_beam,coord='C',title=map_title,rot=(0,0,0),min=0, max=7000)
            fig_name="check_%s_with_beam_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (sky_model,lst_hrs,freq_MHz,pol)
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print("saved %s" % fig_name)
            
            #beam only
            plt.clf()
            map_title="GSM from MWA at %s hrs LST %0.3f MHz" % (lst_hrs,freq_MHz)
            ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=short_dipole_parallel_beam_map,half_sky=False,title=map_title,min=0, max=1)
            #hp.mollview(map=short_dipole_parallel_beam_map,coord='C',title=map_title,rot=(0,0,0),min=0, max=1)
            fig_name="check_beam_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (lst_hrs,freq_MHz,pol)
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print("saved %s" % fig_name)      
         
            #rotated beam only
            plt.clf()
            map_title="GSM from MWA at %s hrs LST %0.3f MHz" % (lst_hrs,freq_MHz)
            ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=short_dipole_parallel_beam_map,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=1)
            #hp.mollview(map=short_dipole_parallel_beam_map,coord='C',title=map_title,rot=(0,0,0),min=0, max=1)
            fig_name="check_rotated_beam_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (lst_hrs,freq_MHz,pol)
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print("saved %s" % fig_name)     
            
            
            #sky only
            plt.clf()
            map_title="GSM from MWA at %s hrs LST %0.3f MHz" % (lst_hrs,freq_MHz)
            ##hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=gsm_map,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
            #hp.mollview(map=gsm_map,coord='C',title=map_title,rot=(0,0,0),min=0, max=7000)
            fig_name="check_%s_%s_hrs_LST_%0.3f_MHz_%0.3f_pol.png" % (sky_model,lst_hrs,freq_MHz,pol)
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print("saved %s" % fig_name)    
                 
            ##sky only rotated
            #plt.clf()
            #map_title="GSM from MWA at %s hrs LST %s MHz" % (lst_hrs,int(freq_MHz))
            ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            #hp.orthview(map=rotated_gsm_map,half_sky=False,title=map_title,rot=(0,0,0),min=0, max=7000)
            ##hp.mollview(map=gsm_map,coord='C',title=map_title,rot=(0,0,0),min=0, max=7000)
            #fig_name="check_rotated_%s_%s_hrs_LST_%s_MHz_%s_pol.png" % (sky_model,lst_hrs,int(freq_MHz),pol)
            #figmap = plt.gcf()
            #figmap.savefig(fig_name,dpi=500)
            #print("saved %s" % fig_name)    
                 
         
         sum_of_beam_weights = np.nansum(short_dipole_parallel_beam_map)
         #print(sum_of_beam_weights)
         beam_weighted_av_sky = np.nansum(sky_with_beam) /  sum_of_beam_weights  #(2.*np.pi/float(n_pix)) #
         print("beam_weighted_av_sky at %0.3f MHz is %0.4E" % (freq_MHz,beam_weighted_av_sky))
         
         
         ##########
         
         #print sky_vector_array.shape
         #print baseline_vector_array.shape
         
         #n_element_for_dot_array = n_baselines * n_pix
         #print n_element_for_dot_array
         #array_size = n_element_for_dot_array * 8
         #array_size_Gb = array_size/1000000000.
         #print "%E Gb" % array_size_Gb
         #haha 800 Gb!
         #big_baseline_vector_for_dot_array = np.full(n_element_for_dot_array,0)
         
         #calculate the expected unity sky signal analytically
         if not fast:
            if calculate_uniform_response:
               print("calculating uniform response")
               for baseline_vector_index in range(0,n_baselines_included):
                  #Just try doing the integral (sum) all in one go with one baseline
                  #baseline_vector_test = baseline_vector_array[0]
                  
                  baseline_vector_for_dot_array = baseline_vector_array[baseline_vector_index,:]
                  
                  baseline_vector_for_dot_array_pure_parallel = baseline_vector_array_pure_parallel[baseline_vector_index,:]
                  baseline_vector_for_dot_array_pure_inline = baseline_vector_array_pure_inline[baseline_vector_index,:]
                  
                  #baseline_vector_for_dot_array_mag = np.linalg.norm(baseline_vector_for_dot_array)
                  #print baseline_vector_for_dot_array_mag
               
                  #b_dot_r_single = np.dot(baseline_vector_test_for_dot_array[4],sky_vector_array[4])
                  #print b_dot_r_single
                  #print baseline_vector_for_dot_array[0:3]
                  b_dot_r_array = (baseline_vector_for_dot_array * sky_vector_array).sum(axis=1)
             
                  b_dot_r_array_pure_parallel = (baseline_vector_for_dot_array_pure_parallel * sky_vector_array).sum(axis=1)
                  b_dot_r_array_pure_inline = (baseline_vector_for_dot_array_pure_inline * sky_vector_array).sum(axis=1)
                  
                  
                  phase_angle_array = 2.*np.pi*b_dot_r_array/wavelength
                  
                  phase_angle_array_pure_parallel = 2.*np.pi*b_dot_r_array_pure_parallel/wavelength
                  phase_angle_array_pure_inline = 2.*np.pi*b_dot_r_array_pure_inline/wavelength
                  
                  #Need to rotate the phase angle maps to match rotated beam and the sky
                  #might change this so that we just rotate the sky map instead - cause this rotation is in the baseline loop and 
                  #so it slows everything down a lot!
                  #phase_angle_array = rotate_map(phase_angle_array, rot_theta_beam, rot_phi_beam)
                  #phase_angle_array_pure_parallel = rotate_map(phase_angle_array_pure_parallel, rot_theta_beam, rot_phi_beam)
                  #phase_angle_array_pure_inline = rotate_map(phase_angle_array_pure_inline, rot_theta_beam, rot_phi_beam)
                  
                  #element_iso_array = iso_beam_map*np.exp(-1j*phase_angle_array)
                  element_short_parallel_array = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array)
                  
                  element_short_parallel_array_pure_parallel = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array_pure_parallel)
                  element_short_parallel_array_pure_inline = short_dipole_parallel_beam_map * np.exp(-1j*phase_angle_array_pure_inline)
                  
                  #for angular info
                  if include_angular_info:
                     element_short_parallel_angular_array = short_dipole_parallel_beam_map * gsm_map_angular * np.exp(-1j*phase_angle_array)
                     #print short_dipole_parallel_beam_map[0:100]
                     #print gsm_map_angular[:3]
                     #print phase_angle_array[0:3]
                     #print element_short_parallel_angular_array[0:3]
                  
                  #print element_iso_array[0:5]
                  #visibility_iso = np.sum(element_iso_array)
                  # Maybe this should be (1/2.*np.pi) as we are only looking at half the sky .... or 1/4pi as in singh et al?
                  
                  #This one gives an answer that is almost exactly 6 PI too big .... (making this smaller makes Tsky smaller)
                  #X_short_parallel = (1./(4.*np.pi)) * np.sum(element_short_parallel_array) * (2.*np.pi/float(n_pix))
                  
                  #this one gives approximately the right answer ....  no exactly!
                  X_short_parallel =  np.sum(element_short_parallel_array) * pixel_solid_angle # (4.*np.pi/float(n_pix))
            
                  X_short_parallel_pure_parallel =  np.sum(element_short_parallel_array_pure_parallel) * pixel_solid_angle # (4.*np.pi/float(n_pix))
                  X_short_parallel_pure_inline =  np.sum(element_short_parallel_array_pure_inline) * pixel_solid_angle # (4.*np.pi/float(n_pix))
            
            
                  if include_angular_info:
                     Y_short_parallel_angular =  np.sum(element_short_parallel_angular_array) * pixel_solid_angle
                  
                  #X_short_parallel = (1./(4.*np.pi)) * np.sum(element_short_parallel_array) * (float(n_pix))
                  #get rid of 1/4pi?, but divide by 2 cause only 1 polarisation (this is a fudge, I dunno about this unit conversion at all......)
                  #X_short_parallel =  (1./2.) * np.sum(element_short_parallel_array)
                  
                  #print visibility_iso
                  #only interested in the real component (that has the global signal)
                  X_short_parallel_array[baseline_vector_index] = X_short_parallel
                  
                  X_short_parallel_array_pure_parallel[baseline_vector_index] = X_short_parallel_pure_parallel
                  X_short_parallel_array_pure_inline[baseline_vector_index] = X_short_parallel_pure_inline
                  
                  
                  if include_angular_info:
                     Y_short_parallel_angular_array[baseline_vector_index] = Y_short_parallel_angular
                  
               #print X_short_parallel_array
               #save X_short_parallel_array
               X_short_parallel_array_filename = "X_short_parallel_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
               np.save(X_short_parallel_array_filename,X_short_parallel_array)
            
               print("saved %s" % X_short_parallel_array_filename)
            
               X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
               np.save(X_short_parallel_array_filename_pure_inline,X_short_parallel_array_pure_inline)
               print("saved %s" % X_short_parallel_array_filename_pure_inline)
               
               X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
               np.save(X_short_parallel_array_filename_pure_parallel,X_short_parallel_array_pure_parallel)
               print("saved %s" % X_short_parallel_array_filename_pure_parallel)
               
               #update for fine chans
               if include_angular_info:
                  Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
                  np.save(Y_short_parallel_angular_array_filename,Y_short_parallel_angular_array)
                  print("saved %s" % Y_short_parallel_angular_array_filename)
               
            real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
            np.save(real_vis_data_sorted_array_filename,real_vis_data_sorted[0:n_baselines_included])
            print("saved %s" % real_vis_data_sorted_array_filename)
            
         else:
            print("already got expected unity sky values, X, from miriad sims")              
            pass
            #rm miriad_X_vis_filename
            #cmd = "rm -rf %s tmp.vis" % miriad_X_vis_filename
            #print(cmd)
            #os.system(cmd)
            #
            ##Instead of this try using pyuvdata to read in ms and output a uvfits (to be read in using miriad!)
            ##
            ###read in uvfits filename 
            #cmd = "fits in=%s op=uvin out=tmp.vis" % (uvfits_filename)
            #print(cmd)
            #os.system(cmd)
        
            #cmd = "prthd in=tmp.vis "
            #print(cmd)
            #os.system(cmd)
            #
            #cmd = "uvlist vis=tmp.vis options=array,full"
            #print(cmd)
            #os.system(cmd)

            
            
            ##something is wrong when u do uvgen - try exporting a s auvfits and re-importing
            #cmd = "fits in=tmp.vis op=uvout out=tmp.uvfits" 
            #print(cmd)
            #os.system(cmd)
            

            
            #cmd = "fits in=tmp.uvfits op=uvin out=%s" % miriad_X_vis_filename
            #print(cmd)
            #os.system(cmd)
            #         
            #cmd = "rm -rf tmp.vis" 
            #print(cmd)
            #os.system(cmd)
            
            ###pyuvdata:
            #tmp_uv = UVData()
            #tmp_uv.read_ms(ms_filename)
            #
            #sys.exit()
            #nup
            ########################
    
   return(beam_weighted_av_sky)
   
         
         ##normalise both X and real vis to max 1
         ##Don't need to do this plot here anymore, it is done later (so that it can be done in plot_only mode)
         #X_short_parallel_array_max = np.max(X_short_parallel_array)
         #print("X_short_parallel_array_max %E" % X_short_parallel_array_max)
         #X_short_parallel_array_norm = X_short_parallel_array / X_short_parallel_array_max
         #
         #real_vis_data_sorted_max = np.max(real_vis_data_sorted[0:n_baselines_included])
         #print("real_vis_data_sorted_max %E" % real_vis_data_sorted_max)
         #real_vis_data_sorted_norm = real_vis_data_sorted[0:n_baselines_included] / real_vis_data_sorted_max
         #
         #plt.clf()
         #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='X norm')
         #plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_norm,s=1,label='real vis norm')
         ##plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         #map_title="X vs uvdistance" 
         #plt.xlabel("uv-distance (lambda)")
         #plt.ylabel("X")
         #plt.legend(loc=1)
         ##plt.ylim([0, 20])
         #fig_name= "X_vs_uv_dist_%d_MHz_%s_pol%s.png" % (freq_MHz,pol,signal_type_postfix)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name) 
         
         
         #plot X_short_parallel_array vs uvdist
          
         
         #Solve for Tsky
         #From Vedantham et al 2015 Moon and my model_moon.py
         #theta_hat2=np.matmul(np.matmul(np.linalg.inv(np.matmul(H2_matrix.T,H2_matrix)),H2_matrix.T),vec_D)
         #for us:
         # vis = X Tsky + noise
         
         ##Don't need all this stuff as you do it in model_t_sky function from saved data
         #X_transpose = np.transpose(X_short_parallel_array)
         #sum_of_squares = X_transpose.dot(X_short_parallel_array)
         #sum_of_squares_inverse = sum_of_squares**(-1) 
         ##print sum_of_squares_inverse
         #sum_of_squares_inverse_Xt = sum_of_squares_inverse * (X_transpose)
         #t_sky_jy = sum_of_squares_inverse_Xt.dot(real_vis_data_sorted[0:n_baselines_included])
         #print "t_sky is %0.2f Jy" % t_sky_jy
         ##convert jy to K
         ##jy_to_K = (2*k*(2*np.pi)) / (wavelength**2 * 10**(-26))
         ##jy_to_K = (2*k) / (10**(-26))
         ##from simulate where we convert to Jy/pix: scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
         #jy_to_K = (wavelength**2) / (2. * k * 1.0e26)   # * 6*PI? (or 4 pi and then fiddle X again?) There are also the visibility weights that I have ignored .... a factor of 120.8
         #print("jy_to_K %.4E" % jy_to_K)
         #t_sky_K = jy_to_K * t_sky_jy
         #print("t_sky_K is %0.4E" % t_sky_K)
        # 
         ##Residuals of fit
         ##print real_vis_data_sorted[0:n_baselines_included]
         ##print X_short_parallel_array
         #residuals_jy = real_vis_data_sorted[0:n_baselines_included] - (X_short_parallel_array * t_sky_jy)
         ##print residuals_jy
         #residuals_K = jy_to_K * residuals_jy
         #rms_residual_K = np.sqrt(np.mean(np.square(residuals_K)))
         ##print rms_residual_K
         #
         #parameter_covariance = rms_residual_K**2 * sum_of_squares_inverse
         #print "parameter_covariance"
         #print parameter_covariance
         
         #t_sky_error_K = np.sqrt(parameter_covariance)
          
         #print("t_sky_K = %0.4E +- %0.3f K for %s shortest baselines included" % (t_sky_K,t_sky_error_K,n_baselines_included))
      
      
      
         
         #ratio = beam_weighted_av_sky / t_sky_K
         #print("ratio between beam_weighted_av_sky %0.4E K and t_sky_K %0.4E +- %0.3f K is %0.3f" % (beam_weighted_av_sky,t_sky_K,t_sky_error_K,ratio)) 
         
         #ratio = beam_weighted_av_sky / t_sky_K
         #print("ratio between beam_weighted_av_sky %0.4E K and t_sky_K %0.4E +- %0.3f K is %0.3f" % (beam_weighted_av_sky,t_sky_K,t_sky_error_K,ratio)) 
      
         #Don't do this here, just save Y as above and then use random intercept LS regression later (Y is used to group X ....)
         #if include_angular_info:
         #   if n_baselines_included > 0:
         #      H_matrix = np.column_stack((X_short_parallel_array.real,Y_short_parallel_angular_array.real))
         #      theta_hat = np.matmul(np.matmul(np.linalg.inv(np.matmul(H_matrix.T,H_matrix)),H_matrix.T),real_vis_data_sorted[0:n_baselines_included])
         #      print "theta_hat"
         #      print theta_hat
         #      t_sky_global_K = theta_hat[0] * jy_to_K
         #      t_sky_ang_K = theta_hat[1] * jy_to_K
         #      print "t_sky_global_K is %E K" % t_sky_global_K
         #      print "t_sky_ang_K is %E K" % t_sky_ang_K
         #   else:
         #      t_sky_global_K = np.nan
         #      t_sky_ang_K = np.nan
         #else:
         #   t_sky_global_K = np.nan
         #   t_sky_ang_K = np.nan
         
         #don't return t_sky, error and n_baselines as these are derived from saved files later on oin model_t_sky
         #return(t_sky_K,t_sky_error_K,beam_weighted_av_sky,n_baselines_included)
      
         
         
         #make diagnostic plots
         
         #fitsname="sky_with_beam_h_m_s_%s.fits" % (time_string)
         #hp.write_map(fitsname,sky_with_beam,dtype=np.float32, overwrite=True)
         #plt.clf()
         #figname="sky_with_beam_h_m_s_%s.png" % (time_string)
         #title = "sky map"
         ##short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
         ##hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
         ##iso_beam_map = hp.read_map(beam_map_fitsname)
         #map_projected=hp.orthview(map=sky_with_beam,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=title,rot=(0,90,0))
         #figmap = plt.gcf()
         #figmap.savefig(figname,dpi=500,bbox_inches='tight') 
         #print "saved figure: %s" % figname
         #plt.close() 
      
         #fitsname="beam.fits" 
         #hp.write_map(fitsname,short_dipole_parallel_beam_map,dtype=np.float32, overwrite=True)
         #plt.clf()
         #figname="beam.png" 
         #title = "beam map"
         ##short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
         ##hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
         ##iso_beam_map = hp.read_map(beam_map_fitsname)
         #map_projected=hp.orthview(map=short_dipole_parallel_beam_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=title,rot=(0,90,0))
         #figmap = plt.gcf()
         #figmap.savefig(figname,dpi=500,bbox_inches='tight') 
         #print "saved figure: %s" % figname
         #plt.close()    

def joint_model_fit_t_sky_measured(lst_hrs_list,freq_MHz_list,pol_list,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,poly_order_list,global_signal_model,plot_only=False):

   poly_order_list_string = '_'.join(str(e) for e in poly_order_list)
   pol = pol_list[0]
   freq_MHz_array = np.asarray(freq_MHz_list)
   
   lst_string = '_'.join(lst_hrs_list)
   lst_hrs = lst_hrs_list[0]
   
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse global at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DG_%s' % sky_model
          concat_output_name_base += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse angular at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DA_%s' % sky_model
          concat_output_name_base += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print('cant have global_EDGES and global in signal_type_list')
          sys.exit()
       else:
          signal_type_postfix += '_G' 
          concat_output_name_base += '_G' 
   if 'global_EDGES' in signal_type_list:
       signal_type_postfix += '_ED' 
       concat_output_name_base += '_ED' 
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'


   t_sky_measured_array_filename = "t_sky_measured_array_lsts_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
   #t_sky_measured_error_array_filename = "t_sky_measured_error_array_lsts_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
   
   t_sky_measured_array = np.load(t_sky_measured_array_filename)
   #get rid of nans:
   #because nan != nan:
   okay_indices = t_sky_measured_array == t_sky_measured_array
   t_sky_measured_array_okay = t_sky_measured_array[okay_indices]
   freq_MHz_array_okay = freq_MHz_array[okay_indices]
   global_signal_model_okay = global_signal_model[okay_indices]
   
   joint_fit_list = []
   joint_fit_global_EDGES_list = []
   joint_fit_diffuse_global_foreground_list = []
   for poly_order in poly_order_list:
      if (('diffuse' in signal_type_list or 'diffuse_global' in signal_type_list) and ('global' not in signal_type_list and 'global_EDGES' not in signal_type_list)):
         print('fitting just diffuse model (no global signal)')
         #jointly fit model to data
         if poly_order==5:
            fitpars, covmat = curve_fit(diffuse_global_foreground_func_order_5,freq_MHz_array_okay,t_sky_measured_array_okay)
            joint_fit = diffuse_global_foreground_func_order_5(freq_MHz_array_okay,*fitpars)
         elif poly_order==7:
            fitpars, covmat = curve_fit(diffuse_global_foreground_func_order_7,freq_MHz_array_okay,t_sky_measured_array_okay)
            joint_fit = diffuse_global_foreground_func_order_7(freq_MHz_array_okay,*fitpars)
         elif poly_order==6:
            fitpars, covmat = curve_fit(diffuse_global_foreground_func_order_6,freq_MHz_array_okay,t_sky_measured_array_okay)
            joint_fit = diffuse_global_foreground_func_order_6(freq_MHz_array_okay,*fitpars)
   
         
      if (('global_EDGES' in signal_type_list) and ('diffuse' not in signal_type_list and 'diffuse_global' not in signal_type_list)):
         print('fitting just global_EDGES model (no diffuse foreground)')
         #jointly fit model to data
         fitpars, covmat = curve_fit(global_sig_EDGES_func,freq_MHz_array_okay,t_sky_measured_array_okay)
   
         joint_fit = global_sig_EDGES_func(freq_MHz_array_okay,*fitpars)   

      if (('global_EDGES' in signal_type_list) and ('diffuse' in signal_type_list or 'diffuse_global' in signal_type_list)):
         print('fitting global_EDGES model and diffuse foreground)')
         #jointly fit model to data
         if poly_order==5:
            fitpars, covmat = curve_fit(global_sig_EDGES_and_diffuse_fg_func_order_5,freq_MHz_array_okay,t_sky_measured_array_okay)
            joint_fit = global_sig_EDGES_and_diffuse_fg_func_order_5(freq_MHz_array_okay,*fitpars)
            joint_fit_diffuse_global_foreground = diffuse_global_foreground_func_order_5(freq_MHz_array_okay,fitpars[1],fitpars[2],fitpars[3],fitpars[4],fitpars[5])
   
         elif poly_order==6:
            fitpars, covmat = curve_fit(global_sig_EDGES_and_diffuse_fg_func_order_6,freq_MHz_array_okay,t_sky_measured_array_okay)
            joint_fit = global_sig_EDGES_and_diffuse_fg_func_order_6(freq_MHz_array_okay,*fitpars)
            joint_fit_diffuse_global_foreground = diffuse_global_foreground_func_order_6(freq_MHz_array_okay,fitpars[1],fitpars[2],fitpars[3],fitpars[4],fitpars[5],fitpars[6])
         
         elif poly_order==7:
            fitpars, covmat = curve_fit(global_sig_EDGES_and_diffuse_fg_func_order_7,freq_MHz_array_okay,t_sky_measured_array_okay)            
            joint_fit = global_sig_EDGES_and_diffuse_fg_func_order_7(freq_MHz_array_okay,*fitpars)     
            joint_fit_diffuse_global_foreground = diffuse_global_foreground_func_order_7(freq_MHz_array_okay,fitpars[1],fitpars[2],fitpars[3],fitpars[4],fitpars[5],fitpars[6],fitpars[7])
   
         
       
         joint_fit_global_EDGES = global_sig_EDGES_func(freq_MHz_array_okay,fitpars[0])
         
         joint_fit_list.append(joint_fit)
         joint_fit_global_EDGES_list.append(joint_fit_global_EDGES)
         joint_fit_diffuse_global_foreground_list.append(joint_fit_diffuse_global_foreground)
         
      #variances = covmat.diagonal()
      #std_devs = np.sqrt(variances)
      #print fitpars,std_devs
    
   #plot the recovered polynomial and recovered global signal separately 
   plt.clf()
   for joint_fit_global_EDGES_index,joint_fit_global_EDGES in enumerate(joint_fit_global_EDGES_list):
      plt.plot(freq_MHz_array_okay,joint_fit_global_EDGES,label='recovered order %s' % poly_order_list[joint_fit_global_EDGES_index])
   plt.plot(freq_MHz_array_okay,global_signal_model_okay,label='EDGES input',linestyle='dashed')
   map_title="t_sky gobal EDGES joint fit" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("T_sky (K)")
   plt.legend(loc="lower right")
   #plt.ylim([0, 20])
   fig_name= "t_sky_joint_fit_global_EDGES_LST_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order_list_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   plt.clf()
   for joint_fit_global_EDGES_index,joint_fit_global_EDGES in enumerate(joint_fit_global_EDGES_list):
      plt.plot(freq_MHz_array_okay,joint_fit_diffuse_global_foreground,label='Diffuse poly order %s' % poly_order_list[joint_fit_global_EDGES_index])
   map_title="t_sky diffuse joint fit" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("T_sky (K)")
   plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "t_sky_joint_fit_diffuse_global_foreground_LST_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order_list_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
   
   #residuals from fit:\
   plt.clf()
   for joint_fit_index,joint_fit in enumerate(joint_fit_list):
      residuals = t_sky_measured_array_okay - joint_fit
      rms_of_residuals = np.sqrt(np.mean(residuals**2))
      max_abs_residuals = np.max(np.abs(residuals)) * 0.9

      #plot residuals
      plt.plot(freq_MHz_array_okay,residuals,label='Residuals poly order %s' % poly_order_list[joint_fit_index] )
      
   map_title="residuals from joint fitting" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Residual T_sky (K)")
   plt.legend(loc=1)
   plt.text(50, max_abs_residuals, "rms=%0.3f" % rms_of_residuals)
   #plt.ylim([0, 20])
   fig_name= "t_sky_residuals_joint_fit_LST_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order_list_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
   
   
   #these two lines are equivalent
   #pl.plot(x, sinfunc(x, fitpars[0], fitpars[1]), 'r-')
   #pl.plot(x, sinfunc(x, *fitpars), 'r-')

   #plot the measured t sky and joint fit
   plt.clf()
   plt.plot(freq_MHz_array_okay,t_sky_measured_array_okay,label='Tsky measured')
   plt.plot(freq_MHz_array_okay,joint_fit,label='Tsky joint fit')
   map_title="t_sky and joint fit" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("T_sky (K)")
   plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "t_sky_measured_and_joint_fit_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)


def plot_tsky_for_multiple_freqs(lst_hrs_list,freq_MHz_list,pol_list,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,poly_order,plot_only=False,include_angular_info=False,model_type_list=['OLS_fixed_intercept'],EDA2_data=False,EDA2_chan_list='None',n_obs_concat_list=[],wsclean=False,fast=False,no_modelling=False,calculate_uniform_response=True):
   #for plot_expected_rms_noise_eda2 below
   int_time = 0.28 * 5.
   bw_Hz = 27./32. * fine_chan_width_Hz

   pol = pol_list[0]
   freq_MHz_array = np.asarray(freq_MHz_list)

   if EDA2_data:
      n_edge_chans_omitted = 5 #two at start and 3 at end
      n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
   else:
      n_fine_chans_used = 1
      n_edge_chans_omitted = 0
   
   #lst_string = '_'.join(lst_hrs_list)
   #just use the first LST for LST string
   lst_string = lst_hrs_list[0]
   lst_hrs = lst_hrs_list[0]
   
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   
   
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
   if 'global_unity' in signal_type_list:
       signal_type_postfix += '_GU'
       concat_output_name_base += '_GU'
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse global at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DG_%s' % sky_model
          concat_output_name_base += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse angular at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DA_%s' % sky_model
          concat_output_name_base += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print('cant have global_EDGES and global in signal_type_list')
          sys.exit()
       else:
          signal_type_postfix += '_G' 
          concat_output_name_base += '_G'
   if 'global_EDGES' in signal_type_list:
       signal_type_postfix += '_ED' 
       concat_output_name_base += '_ED' 
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'
   
   if EDA2_data==True:
      signal_type_postfix = "_EDA2_data"
   
   
   t_sky_theoretical_array = np.full(len(freq_MHz_list),np.nan)
   #don't use this as you have one for each fine chan now
   n_baselines_used_array = np.full(len(freq_MHz_list),np.nan)
   
   #for including angular info:
   #t_sky_measured_global_array = np.full(len(freq_MHz_list),np.nan)
   #t_sky_measured_angular_array = np.full(len(freq_MHz_list),np.nan)

   
   t_sky_theoretical_array_filename = "t_sky_theoretical_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   n_baselines_used_array_filename = "n_baselines_used_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   
   ##for including angular info:
   t_sky_measured_global_array_filename = "t_sky_measured_global_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   t_sky_measured_angular_array_filename = "t_sky_measured_angular_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   
   #freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   
   
   if not plot_only:
      for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
         #using global variable chan num so I can run lots in parallel (yes I know this is bad programming)
         if EDA2_data==True:
            if len(freq_MHz_list)==1:
               EDA2_chan = EDA2_chan_list[0]
               EDA2_obs_time = EDA2_obs_time_list[chan_num]
            else:
               EDA2_chan = EDA2_chan_list[freq_MHz_index]
               EDA2_obs_time = EDA2_obs_time_list[freq_MHz_index]
            if len(n_obs_concat_list) > 0:
               if len(freq_MHz_list)==1:
                  n_obs_concat = n_obs_concat_list[chan_num]
               else:
                  n_obs_concat = n_obs_concat_list[freq_MHz_index]
            else:
               n_obs_concat = 1
         else:
            EDA2_chan = None
            EDA2_obs_time = None
            n_obs_concat = 1
         
         #only the theoretical beam weighted av is taken from this function now,, rest is derived from saved files
         #t_sky_measured,t_sky_measured_error,t_sky_theoretical,n_baselines_used = solve_for_tsky_from uvfits(freq_MHz,lst_hrs_list,pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,include_angular_info=include_angular_info,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat)
         
         ###THis is the one that works for the SIMS! at 20200422
         if EDA2_data:
            if wsclean:
               t_sky_theoretical = solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,include_angular_info=include_angular_info,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,wsclean=wsclean,fast=fast,calculate_uniform_response=calculate_uniform_response)
            else:
               t_sky_theoretical = extract_data_from_eda2_uvfits(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs_list=lst_hrs_list,pol=pol,EDA2_chan=EDA2_chan,n_obs=n_obs_concat,calculate_uniform_response=calculate_uniform_response)
         else:
            t_sky_theoretical = solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,include_angular_info=include_angular_info,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,wsclean=wsclean,fast=fast,calculate_uniform_response=calculate_uniform_response)
         #t_sky_measured_array[freq_MHz_index] = t_sky_measured
         #t_sky_measured_error_array[freq_MHz_index] = t_sky_measured_error
         t_sky_theoretical_array[freq_MHz_index] = t_sky_theoretical
         #n_baselines_used_array[freq_MHz_index] = n_baselines_used  
      #np.save(t_sky_measured_array_filename,t_sky_measured_array)
      #np.save(t_sky_measured_error_array_filename,t_sky_measured_error_array)
      np.save(t_sky_theoretical_array_filename,t_sky_theoretical_array)
      #np.save(n_baselines_used_array_filename,n_baselines_used_array)
      #if include_angular_info:
      #   np.save(t_sky_measured_global_array_filename,t_sky_measured_global_array)
      #   np.save(t_sky_measured_angular_array_filename,t_sky_measured_angular_array)
   
   #this array needs to be bigger now as we are using the fine channels of EDA2 data
   
   
   if not no_modelling:
   
      t_sky_array_length = int(len(freq_MHz_list) * n_fine_chans_used)
      t_sky_measured_array = np.full(t_sky_array_length,np.nan)
      t_sky_measured_error_array = np.full(t_sky_array_length,np.nan)
      t_sky_measured_array_flagged = np.full(t_sky_array_length,np.nan)
      t_sky_measured_error_array_flagged = np.full(t_sky_array_length,np.nan)
      freq_MHz_fine_array = np.full(t_sky_array_length,np.nan)
      
      for model_type in model_type_list:
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
         t_sky_measured_array_filename_flagged = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename_flagged = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
   
         #this replaces all the matrix stuff you do in model_tsky_from_saved_data
         
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            if EDA2_data==True:
               EDA2_chan = EDA2_chan_list[freq_MHz_index]
               if len(n_obs_concat_list) > 0:
                  #n_obs_concat = n_obs_concat_list[freq_MHz_index]
                  n_obs_concat = n_obs_concat_list[0]
               else:
                  n_obs_concat = 1
               #fine_chan_index_array = range(n_fine_chans)[n_chans_omitted_each_edge:n_fine_chans-n_chans_omitted_each_edge]
               #omit 2 fine chans at the low end and 3 edge chans at the high end and you are sweet, no gaps
               #fine_chan_index_array = n_fine_chans_used - np.arange(n_fine_chans-5)+2
               fine_chan_index_array = np.arange(n_fine_chans-5)+2
               for fine_chan_index_index,fine_chan_index in enumerate(fine_chan_index_array):
                  freq_MHz_index_fine = freq_MHz_index*n_fine_chans_used + fine_chan_index_index
                  ######freq_MHz_index_fine = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
                  #oversampled PFB - dont need all this edge chan stuff
                  #channel_remainder = freq_MHz_index_fine % n_fine_chans
                  #print channel_remainder
                  #if (channel_remainder < n_chans_omitted_each_edge or channel_remainder >= (n_fine_chans-n_chans_omitted_each_edge)):
                  #   edge_chan=True
                  #else:
                  #   edge_chan=False
                  edge_chan=False 
                  
                  if EDA2_data:
                     if wsclean:
                        t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs=lst_hrs,pol=pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,model_type=model_type,EDA2_data=EDA2_data,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,fine_chan_index=fine_chan_index,edge_chan=edge_chan,wsclean=wsclean,fast=fast)
                     else:
                        t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data_eda2(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs_list=lst_hrs_list,pol=pol,EDA2_chan=EDA2_chan,n_obs=n_obs_concat,fine_chan_index=fine_chan_index,model_type=model_type)
                  else:
                     t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs=lst_hrs,pol=pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,model_type=model_type,EDA2_data=EDA2_data,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,fine_chan_index=fine_chan_index,edge_chan=edge_chan,wsclean=wsclean,fast=fast)
                  
                  print(freq_MHz_fine)
                  print(freq_MHz_index_fine)
                  
                  
                  t_sky_measured_array[freq_MHz_index_fine] = t_sky_measured
                  t_sky_measured_error_array[freq_MHz_index_fine] = t_sky_measured_error
                  t_sky_measured_array_flagged[freq_MHz_index_fine] = t_sky_measured_flagged
                  t_sky_measured_error_array_flagged[freq_MHz_index_fine] = t_sky_measured_error_flagged
                  freq_MHz_fine_array[freq_MHz_index_fine] = freq_MHz_fine
         
               
            else:
               EDA2_chan = None
               edge_chan=False 
               n_obs_concat = 1
               fine_chan_index=0
               t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs=lst_hrs,pol=pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,model_type=model_type,EDA2_data=EDA2_data,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,fine_chan_index=fine_chan_index,edge_chan=edge_chan,wsclean=wsclean,fast=fast)
               t_sky_measured_array[freq_MHz_index] = t_sky_measured
               t_sky_measured_error_array[freq_MHz_index] = t_sky_measured_error
               t_sky_measured_array_flagged[freq_MHz_index] = t_sky_measured_flagged
               t_sky_measured_error_array_flagged[freq_MHz_index] = t_sky_measured_error_flagged
               freq_MHz_fine_array = freq_MHz_array    
               
         np.save(t_sky_measured_array_filename,t_sky_measured_array)
         np.save(t_sky_measured_error_array_filename,t_sky_measured_error_array) 
         np.save(t_sky_measured_array_filename_flagged,t_sky_measured_array_flagged)
         np.save(t_sky_measured_error_array_filename_flagged,t_sky_measured_error_array_flagged)       
         np.save(freq_MHz_fine_array_filename,freq_MHz_fine_array)
  
      
      
      
   #Dont need to load these...  
   #yes you do now.....below!     
   ###t_sky_measured_array = np.load(t_sky_measured_array_filename)
   ####t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
   t_sky_theoretical_array = np.load(t_sky_theoretical_array_filename)
   #n_baselines_used_array = np.load(n_baselines_used_array_filename)
   #if include_angular_info:
   #   t_sky_measured_global_array = np.load(t_sky_measured_global_array_filename)
   #   t_sky_measured_angular_array = np.load(t_sky_measured_angular_array_filename)

   
   #make a plot of diffuse global input  
   if 'diffuse_global' in signal_type_list:
      sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base
      sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base
      diffuse_global_value_array_X = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
      #diffuse_global_value_array_Y = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
      
      plt.clf()
      plt.plot(freq_MHz_array,diffuse_global_value_array_X,label='sim input X')
      #plt.plot(freq_MHz_array,diffuse_global_value_array_Y,label='sim input Y')
      map_title="t_sky beam averaged input" 
      plt.xlabel("freq (MHz)")
      plt.ylabel("t_sky (K)")
      plt.legend(loc=1)
      #plt.ylim([0, 20])
      fig_name= "t_sky_beam_av_input_%s%s.png" % (lst_string,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
    
      #subtract a polynomial fit just do for X for now
      #in log log space:
      sky_array = diffuse_global_value_array_X[t_sky_measured_array>0.]
      log_sky_array = np.log10(sky_array)
      freq_array_cut = freq_MHz_array[t_sky_measured_array>0.]
      log_freq_MHz_array = np.log10(freq_array_cut)
      coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      
      #log_residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - sky_array
      
      rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
      print("rms_of_residuals for X is %0.3f K" % rms_of_residuals)
      
      max_abs_residuals = np.max(np.abs(residual_of_log_fit))
      y_max = 1.5 * max_abs_residuals
      y_min = 1.5 * -max_abs_residuals
      
      plt.clf()
      plt.plot(freq_array_cut,residual_of_log_fit,label='residual of log fit X')
      map_title="Residual for log polynomial order %s fit " % poly_order
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      plt.text(50, max_abs_residuals, "rms=%0.3f" % rms_of_residuals)
      plt.ylim([y_min, y_max])
      fig_name= "eda2_log_fit_residual_tsy_input_poly_%s_lst_%s%s.png" % (poly_order,lst_string,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
   
   
   #Get the theoretical diffuse value from the chan folders instead (in case you have run each chan separately for solve_from_uvfits) 
   if EDA2_data==True:
      label2='predicted GSM'
      t_sky_theoretical_list = []
      EDA2_chan_list_input = EDA2_chan_list[0:len(freq_MHz_list)]
      for EDA2_chan in EDA2_chan_list_input:
            EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)          
            sky_averaged_diffuse_array_beam_lsts_filename = "%s%s_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,concat_output_name_base)       
            diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
            diffuse_global_value = diffuse_global_value_array[0]   
            t_sky_theoretical_list.append(diffuse_global_value)
      t_sky_theoretical_array = np.asarray(t_sky_theoretical_list)      
   else:
      label2 = 'input'
   
   #for plotting a subset of chans:
   if EDA2_data:
      fine_chans_per_EDA2_chan = n_fine_chans-5
      length_freq_MHz_fine_chan_to_plot = int(len(EDA2_chan_list_input) * fine_chans_per_EDA2_chan)
   else:
      length_freq_MHz_fine_chan_to_plot = len(freq_MHz_list)
   
   #unflagged
   plt.clf()
   for model_type in model_type_list:
      #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
      if model_type=='OLS_fixed_intercept':
         if EDA2_data:
            label1='measured sky temp.'
         else:
            label1='ignore angular response'
      elif  model_type=='OLS_fixed_int_subtr_Y':
         label1='subtract angular response'
      else:
         label1='recovered'
         
      t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
      t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
      freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
    
      
      t_sky_measured_array = np.load(t_sky_measured_array_filename)
      t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
      t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
      t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
      freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
      freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
      #print(t_sky_measured_array)
      #print(t_sky_measured_array.shape)
      #print(t_sky_measured_error_array)
      #print(t_sky_measured_error_array.shape)
      #print(freq_MHz_fine_array)
      #print(freq_MHz_fine_array.shape)
      # 
      #print(t_sky_theoretical_array)
      #print(freq_MHz_list)
      
      plt.errorbar(freq_MHz_fine_array,t_sky_measured_array,yerr=t_sky_measured_error_array,label=label1)
   if len(freq_MHz_list)==1:
      plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
   else:
      plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2)
 
   #if 'diffuse_global' in signal_type_list:
   #   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
   #if include_angular_info:
   #   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      print(signal_type_list)
      plt.legend(loc='upper right')
   else:
      plt.legend(loc='lower right')
   if EDA2_data:
      plt.ylim([500, 5000])
   else:
      plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
  
   ###Also plot the average measurement for each EDA2 coarse chan
   t_sky_measure_av_per_EDA2_chan = np.full(len(freq_MHz_list),np.nan)
   t_sky_measure_av_per_EDA2_chan_err = np.full(len(freq_MHz_list),np.nan)
   t_sky_measure_av_per_EDA2_chan_weighted = np.full(len(freq_MHz_list),np.nan)
   t_sky_measure_av_per_EDA2_chan_err_weighted = np.full(len(freq_MHz_list),np.nan)
   
   
   plt.clf()
   
   for model_type in model_type_list:
      #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
      if model_type=='OLS_fixed_intercept':
         if EDA2_data:
            label1='measured sky temp.'
         else:
            label1='ignore angular response'
         label3='ignore angular resp. weighted'
      elif  model_type=='OLS_fixed_int_subtr_Y':
         label1='subtract angular response'
         label3='subtract angular resp. weighted'
      else:
         label1='recovered'
         label3='recovered weighted'
      t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
      t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
      freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
    
      t_sky_measured_array = np.load(t_sky_measured_array_filename)
      t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
      t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
      t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
      freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
      freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
      
       
      if EDA2_data==True:
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            freq_range_min = freq_MHz - (centre_chan_index * fine_chan_width_Hz/1000000.)
            freq_range_max = freq_MHz + (centre_chan_index * fine_chan_width_Hz/1000000.)
            #print(freq_range_min)
            #print(freq_range_max)
            
            indices = np.where(np.logical_and(freq_MHz_fine_array>=freq_range_min,freq_MHz_fine_array<=freq_range_max))
            t_sky_measured_EDA2_chan = t_sky_measured_array[indices]
            t_sky_measured_error_chan = t_sky_measured_error_array[indices]
            #print(t_sky_measured_EDA2_chan)
            #print(freq_MHz_fine_array[indices])
            t_sky_measure_av_per_EDA2_chan[freq_MHz_index] = np.nanmean(t_sky_measured_EDA2_chan)
            t_sky_measure_av_per_EDA2_chan_err[freq_MHz_index] = np.nanstd(t_sky_measured_EDA2_chan)
            
            #weighted average:
            weight_array = 1./(t_sky_measured_error_chan**2)
            weighted_array = t_sky_measured_EDA2_chan * weight_array
            weighted_sum = np.sum(weighted_array)
            sum_of_weights = np.sum(weight_array)
            weighted_mean = weighted_sum / sum_of_weights
            st_dev_w_mean = np.sqrt(1./(sum_of_weights))
            
            t_sky_measure_av_per_EDA2_chan_weighted[freq_MHz_index] = weighted_mean
            t_sky_measure_av_per_EDA2_chan_err_weighted[freq_MHz_index] = st_dev_w_mean
            
      else:
         t_sky_measure_av_per_EDA2_chan = t_sky_measured_array
         t_sky_measure_av_per_EDA2_chan_err = t_sky_measured_error_array
         t_sky_measure_av_per_EDA2_chan_weighted = t_sky_measured_array
         t_sky_measure_av_per_EDA2_chan_err_weighted = t_sky_measured_error_array
         
      plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan,yerr=t_sky_measure_av_per_EDA2_chan_err,label=label1)
      #if EDA2_data:
      #   plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan_weighted,yerr=t_sky_measure_av_per_EDA2_chan_err_weighted,label=label3)
   if len(freq_MHz_list)==1:
      plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
   else:
      plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2)
   
   #if 'diffuse_global' in signal_type_list:
   #   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
   #if include_angular_info:
   #   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      #print(signal_type_list)
      if len(freq_MHz_list)==1:
         plt.legend(loc='lower right')
      else:
         plt.legend(loc='upper right')
   else:
      plt.legend(loc='lower right')
   if EDA2_data:
      plt.ylim([500, 5000])
   else:
      plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s_per_chan_av.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   ####Flagged
   plt.clf()
   for model_type in model_type_list:
      #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
      if model_type=='OLS_fixed_intercept':
         if EDA2_data:
            label1='measured sky temp.'
         else:
            label1='ignore angular response'
      elif  model_type=='OLS_fixed_int_subtr_Y':
         label1='subtract angular response'
      else:
         label1='recovered'
         
      t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
    
    
      t_sky_measured_array = np.load(t_sky_measured_array_filename)
      t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
      t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
      t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
      freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
      freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
       
      plt.errorbar(freq_MHz_fine_array,t_sky_measured_array,yerr=t_sky_measured_error_array,label=label1)
   if len(freq_MHz_list)==1:
      plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
   else:
      plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2)
 
   #if 'diffuse_global' in signal_type_list:
   #   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
   #if include_angular_info:
   #   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured flagged" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      print(signal_type_list)
      plt.legend(loc='upper right')
   else:
      plt.legend(loc='lower right')
   if EDA2_data:
      plt.ylim([500, 5000])
   else:
      plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s_flagged.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
  
   ###Also plot the average measurement for each EDA2 coarse chan
   t_sky_measure_av_per_EDA2_chan = np.full(len(freq_MHz_list),np.nan)
   t_sky_measure_av_per_EDA2_chan_err = np.full(len(freq_MHz_list),np.nan)
   
   plt.clf()
   
   for model_type in model_type_list:
      #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
      if model_type=='OLS_fixed_intercept':
         if EDA2_data:
            label1='measured sky temp.'
         else:
            label1='ignore angular response'
         label3='ignore angular resp. weighted'
      elif  model_type=='OLS_fixed_int_subtr_Y':
         label1='subtract angular response'
         label3='subtract angular resp. weighted'
      else:
         label1='recovered'
         label3='recovered weighted'
      t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
    
      t_sky_measured_array = np.load(t_sky_measured_array_filename)
      t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
      t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
      t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
      freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
      freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
       
      
          
      if EDA2_data==True:
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            freq_range_min = freq_MHz - (centre_chan_index * fine_chan_width_Hz/1000000.)
            freq_range_max = freq_MHz + (centre_chan_index * fine_chan_width_Hz/1000000.)
            #print(freq_range_min)
            #print(freq_range_max)
            
            indices = np.where(np.logical_and(freq_MHz_fine_array>=freq_range_min,freq_MHz_fine_array<=freq_range_max))
            t_sky_measured_EDA2_chan = t_sky_measured_array[indices]
            t_sky_measured_error_chan = t_sky_measured_error_array[indices]
            #print(t_sky_measured_EDA2_chan)
            #print(freq_MHz_fine_array[indices])
            t_sky_measure_av_per_EDA2_chan[freq_MHz_index] = np.nanmean(t_sky_measured_EDA2_chan)
            t_sky_measure_av_per_EDA2_chan_err[freq_MHz_index] = np.nanstd(t_sky_measured_EDA2_chan)
            
            #weighted average:
            weight_array = 1./(t_sky_measured_error_chan**2)
            weighted_array = t_sky_measured_EDA2_chan * weight_array
            weighted_sum = np.sum(weighted_array)
            sum_of_weights = np.sum(weight_array)
            weighted_mean = weighted_sum / sum_of_weights
            st_dev_w_mean = np.sqrt(1./(sum_of_weights))
            
            t_sky_measure_av_per_EDA2_chan_weighted[freq_MHz_index] = weighted_mean
            t_sky_measure_av_per_EDA2_chan_err_weighted[freq_MHz_index] = st_dev_w_mean
            
      else:
         t_sky_measure_av_per_EDA2_chan = t_sky_measured_array
         t_sky_measure_av_per_EDA2_chan_err = t_sky_measured_error_array
         t_sky_measure_av_per_EDA2_chan_weighted = t_sky_measured_array
         t_sky_measure_av_per_EDA2_chan_err_weighted = t_sky_measured_error_array
         
      plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan,yerr=t_sky_measure_av_per_EDA2_chan_err,label=label1)
      #if EDA2_data:
      #   plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan_weighted,yerr=t_sky_measure_av_per_EDA2_chan_err_weighted,label=label3)
   if len(freq_MHz_list)==1:
      plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
   else:
      plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2)
   
   #if 'diffuse_global' in signal_type_list:
   #   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
   #if include_angular_info:
   #   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      #print(signal_type_list)
      if len(freq_MHz_list)==1:
         plt.legend(loc='lower right')
      else:
         plt.legend(loc='upper right')
   else:
      plt.legend(loc='lower right')
   if EDA2_data:
      plt.ylim([500, 5000])
   else:
      plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s_per_chan_av_flagged.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
 
        
   
   
      

   
   #repeat for using angular info
   #if include_angular_info:
   #   #subtract a polynomial fit
   #   #in log log space:
   #   sky_array = t_sky_measured_array[t_sky_measured_array>0.]
   #   log_sky_array = np.log10(sky_array)
   #   freq_array_cut = freq_MHz_array[t_sky_measured_array>0.]
   #   log_freq_MHz_array = np.log10(freq_array_cut)
   #   coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
   #   ffit = poly.polyval(log_freq_MHz_array, coefs)
   #   ffit_linear = 10**ffit
   #   
   #   #log_residual = log_signal_array_short_baselines - log_ffit
   #   residual_of_log_fit_glob_ang = ffit_linear - sky_array
   #   
   #   rms_of_residuals_glob_ang = np.sqrt(np.mean(residual_of_log_fit_glob_ang**2))
   #   print("rms_of_residuals_glob_ang is %0.3f K" % rms_of_residuals_glob_ang)
   

   
   #no error bar plot what we want
   plt.clf()
   for model_type in model_type_list:
      #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
      if model_type=='OLS_fixed_intercept':
         if EDA2_data:
            label1='measured sky temp.'
         else:
            label1='ignore angular response'
         y_offset=1
         colour='tab:blue'
      elif  model_type=='OLS_fixed_int_subtr_Y':
         label1='subtract angular response'
         y_offset=0
         colour='tab:orange'
      else:
         label1='recovered'
         
      t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
      freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
    
    
      t_sky_measured_array = np.load(t_sky_measured_array_filename)
      t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
      t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
      t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
      freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
      freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
       
       
      #subtract a polynomial fit
      #in log log space:
      sky_array = t_sky_measured_array[t_sky_measured_array>0.]
      log_sky_array = np.log10(sky_array)
      if n_fine_chans_used==1:
         freq_array_cut = freq_MHz_array[t_sky_measured_array>0.]
      else:
         freq_array_cut = freq_MHz_fine_array[t_sky_measured_array>0.]
      log_freq_MHz_array = np.log10(freq_array_cut)
      coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      
      #log_residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - sky_array
      
      rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
      print("rms_of_residuals is %0.3f K" % rms_of_residuals)
      
      max_abs_residuals = np.max(np.abs(residual_of_log_fit))
      y_max = 1.5 * max_abs_residuals
      y_min = 1.5 * -max_abs_residuals
   
      #temporary just for paper:
      y_max = 100
      y_min = -100
   
      plt.plot(freq_array_cut,residual_of_log_fit,label=label1)
      #plt.text(50, max_abs_residuals + y_offset, "rms=%0.3f" % rms_of_residuals,{'color': colour})
      plt.text(50, 75, "rms=%0.3f" % rms_of_residuals,{'color': colour})
   
      
   map_title="Residual for log polynomial order %s fit " % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("freq (MHz)")
   if len(model_type_list)>1:
      plt.legend(loc=1)
   
   plt.ylim([y_min, y_max])
   fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_no_e_bars.png" % (poly_order,lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   plt.close()         
            
            
   plt.clf()
   plt.errorbar(freq_array_cut,residual_of_log_fit,yerr=t_sky_measured_error_array[t_sky_measured_array>0.],label='residual of log fit')
   #if include_angular_info:
   #   plt.plot(freq_array_cut,residual_of_log_fit_glob_ang,label='residual of glob ang')
   map_title="Residual for log polynomial order %s fit " % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("freq (MHz)")
   if len(model_type_list)>1:
      plt.legend(loc=1)
   plt.text(50, max_abs_residuals, "rms=%0.3f" % rms_of_residuals)
   plt.ylim([y_min, y_max])
   fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_%s.png" % (poly_order,lst_string,signal_type_postfix,model_type)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   plt.close() 
   
   
   
   
   #repeat for per chan av
   if EDA2_data:
      plt.clf()
      for model_type in model_type_list:
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='ignore angular response'
            else:
               label1='ignore angular response'
            y_offset=1
            colour='tab:blue'
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular response'
            y_offset=0
            colour='tab:orange'
         else:
            label1='recovered'
            
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s_%s.npy" % (lst_string,signal_type_postfix,model_type)
       
       
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
         
         print(length_freq_MHz_fine_chan_to_plot)
         print(t_sky_measured_array)
         print(freq_MHz_fine_array)
        
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            freq_range_min = freq_MHz - (centre_chan_index * fine_chan_width_Hz/1000000.)
            freq_range_max = freq_MHz + (centre_chan_index * fine_chan_width_Hz/1000000.)
            print(freq_range_min)
            print(freq_range_max)
            
            indices = np.where(np.logical_and(freq_MHz_fine_array>=freq_range_min,freq_MHz_fine_array<=freq_range_max))
            t_sky_measured_EDA2_chan = t_sky_measured_array[indices]
            t_sky_measured_error_chan = t_sky_measured_error_array[indices]
            print(t_sky_measured_EDA2_chan)
            print(freq_MHz_fine_array[indices])
            t_sky_measure_av_per_EDA2_chan[freq_MHz_index] = np.nanmean(t_sky_measured_EDA2_chan)
            t_sky_measure_av_per_EDA2_chan_err[freq_MHz_index] = np.nanstd(t_sky_measured_EDA2_chan)
          
         #subtract a polynomial fit
         #in log log space:
         sky_array = t_sky_measure_av_per_EDA2_chan[t_sky_measure_av_per_EDA2_chan>0.]
         log_sky_array = np.log10(sky_array)
         freq_array_cut = freq_MHz_array[t_sky_measure_av_per_EDA2_chan>0.]
         t_sky_theoretical_array_cut = t_sky_theoretical_array[t_sky_measure_av_per_EDA2_chan>0.]

         log_freq_MHz_array = np.log10(freq_array_cut)
         coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
         ffit = poly.polyval(log_freq_MHz_array, coefs)
         ffit_linear = 10**ffit
         
         #log_residual = log_signal_array_short_baselines - log_ffit
         residual_of_log_fit = ffit_linear - sky_array
         
         rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
         print("rms_of_residuals is %0.3f K" % rms_of_residuals)
         
         max_abs_residuals = np.max(np.abs(residual_of_log_fit))
         y_max = 1.5 * max_abs_residuals
         y_min = 1.5 * -max_abs_residuals
         
         print(freq_array_cut)
         print(residual_of_log_fit)
         
         plt.plot(freq_array_cut,residual_of_log_fit,label=label1)
         plt.text(50, max_abs_residuals + y_offset, "rms=%0.3f" % rms_of_residuals,{'color': colour})
         
         
         #include expected noise estimate:
         #expected_noise = plot_expected_rms_noise_eda2(freq_MHz_list=freq_array_cut,t_sky_theoretical_array=t_sky_theoretical_array_cut,int_time=int_time,bandwidth_Hz=bw_Hz)
         #plt.plot(freq_array_cut,expected_noise,label="expected rms noise")
      
      
      map_title="Residual for log polynomial order %s fit " % poly_order
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("freq (MHz)")
      #if len(model_type_list)>1:
      #   plt.legend(loc=1)
      
      plt.legend(loc=1)
      plt.ylim([y_min, y_max])
      fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_no_e_bars_per_chan_av.png" % (poly_order,lst_string,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      plt.close()    
   
               
   #n baselines plot
   plt.clf()
   plt.plot(freq_MHz_list,n_baselines_used_array)
   map_title="n_baselines_included" 
   plt.xlabel("freq (MHz)")
   plt.ylabel("number of baselines")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "n_baselines_included_lst_%s%s_%s.png" % (lst_string,signal_type_postfix,model_type)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   plt.close()  
           
def plot_expected_rms_noise_eda2(freq_MHz_list,t_sky_theoretical_array,int_time,bandwidth_Hz):
   print("plotting expected noise")
   #from EDA1 paper wayth et al 2017, table 2:
   A_eff_array_per_dipole = np.asarray([970.,950.,914.,874.,832.,771.,707.,638.,568.,498.,435.,377.,329.,288.,252.,222.,196.]) / 256.
   freq_MHz_for_A_eff_array_per_dipole = np.asarray([60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.])    
   #fit a line:
   coefs = poly.polyfit(freq_MHz_for_A_eff_array_per_dipole, A_eff_array_per_dipole, 1)
   ffit = poly.polyval(freq_MHz_for_A_eff_array_per_dipole, coefs)
   
   #plot it to check 
   plt.clf()
   plt.plot(freq_MHz_for_A_eff_array_per_dipole,A_eff_array_per_dipole)
   plt.plot(freq_MHz_for_A_eff_array_per_dipole,ffit)
   map_title="A_eff EDA2 dipoles" 
   plt.xlabel("freq (MHz)")
   plt.ylabel("A_eff (m)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "A_eff_EDA2.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   plt.close()  
   
   A_eff_for_calc_array = poly.polyval(freq_MHz_list, coefs)
   
   T_rms1 = t_sky_theoretical_array / np.sqrt(int_time * bandwidth_Hz)
   
   #T_rms2 = ((c**2 / (freq_MHz_list*1000000.)**2) / (2*A_eff_for_calc_array)) * T_rms1
   
   #T_rms3 = T_rms1 / (2*A_eff_for_calc_array)
   
   
   plt.clf()
   plt.plot(freq_MHz_list,T_rms1,label='T_rms1')
   #plt.plot(freq_MHz_list,T_rms2,label='T_rms2')
   #plt.plot(freq_MHz_list,T_rms3,label='T_rms3')
   map_title="T_rms EDA2" 
   plt.xlabel("freq (MHz)")
   plt.ylabel("T_rms (K)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "T_rms_EDA2.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   plt.close() 
   
   return T_rms1  
   

def extract_signal_from_sims(lst_list,freq_MHz_list,pol_list,signal_type_list,outbase_name,sky_model,array_ant_locations_filename,array_label):
   with open(array_ant_locations_filename,'r') as f:
      lines = f.readlines()
   n_ants = int(len(lines))
   n_baselines = n_ants*(n_ants-1) / 2. 
   concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
   concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
         if 'noise' in signal_type_list:
            concat_output_name_base_X += '_N'
            concat_output_name_base_Y += '_N'
         if 'diffuse' in signal_type_list:
            concat_output_name_base_X += '_D_%s' % sky_model
            concat_output_name_base_Y += '_D_%s' % sky_model
         if 'diffuse_global' in signal_type_list:
             concat_output_name_base_X += '_DG'
             concat_output_name_base_Y += '_DG'
         if 'diffuse_angular' in signal_type_list:
             concat_output_name_base_X += '_DA'
             concat_output_name_base_Y += '_DA'
         if 'global' in signal_type_list:
            concat_output_name_base_X += '_G' 
            concat_output_name_base_Y += '_G'
         if 'gain_errors' in signal_type_list:
            concat_output_name_base_X += '_GE'
            concat_output_name_base_Y += '_GE'
         
         if do_image_and_subtr_from_simulated_data:
            uv_type_list = ['original','subtracted']
         else:
            uv_type_list = ['original']
         
         for uv_type in uv_type_list:
            if pol=='X':  
               if uv_type=='original':                      
                  uvfits_name = "%s_concat_lsts.uvfits" % concat_output_name_base_X
                  model_vis_name_base = concat_output_name_base_X
               else:
                  uvfits_name = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_X
                  model_vis_name_base = concat_output_name_base_X + '_sub'
            else:
               if uv_type=='original':
                  uvfits_name = "%s_concat_lsts.uvfits" % concat_output_name_base_Y
                  model_vis_name_base = concat_output_name_base_Y
               else:
                  uvfits_name = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_Y
                  model_vis_name_base = concat_output_name_base_Y + '_sub'
            model_vis_name_base += "_thresh_%0.2f" % (zero_spacing_leakage_threshold)
            print(uvfits_name)
            hdulist = fits.open(uvfits_name)
            #hdulist.info()
            info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
            #print info_string
            uvtable = hdulist[0].data
            uvtable_header = hdulist[0].header
            visibilities = uvtable['DATA']
            #print visibilities.shape
            visibilities_shape = visibilities.shape
            print("visibilities shape")
            print(visibilities_shape)
            
            #get the UU and VV so we can check whether we are using short baselines
            UU_s_array = uvtable['UU']
            UU_m_array = UU_s_array * c   
            VV_s_array = uvtable['VV']
            VV_m_array = VV_s_array * c
            
            n_vis = visibilities.shape[0]
            n_timesteps = n_vis/n_baselines
            #print "n_timesteps %s " % n_timesteps
            timestep_array = np.arange(0,n_timesteps,1)
            
            #No longer can do this, need to determine weights on the fly
            ##weights_array_filename = "weights_LST_%03d_%s_%s_MHz.npy" % (lst_deg,pol,int(freq_MHz))
            #weights_array_filename = "weights_%s.npy" % model_vis_name_base
            #weights_array = np.load(weights_array_filename)
         
               
            signal_array_short_baselines = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_filename = "%s_signal.npy" % (model_vis_name_base)
            signal_array_short_baselines_Tb = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_Tb_filename = "%s_signal_Tb.npy" % (model_vis_name_base)
            signal_array_all_baselines = np.full((n_pol,n_chan),0.0)
            signal_array_all_baselines_filename = "%s_signal_all_baselines.npy" % (model_vis_name_base)
            signal_array_all_baselines_Tb = np.full((n_pol,n_chan),0.0)
            signal_array_all_baselines_Tb_filename = "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
            signal_array_all_baselines_abs = np.full((n_pol,n_chan),0.0)
            signal_array_all_baselines_filename_abs = "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
            signal_array_all_baselines_abs_Tb = np.full((n_pol,n_chan),0.0)
            signal_array_all_baselines_filename_abs_Tb = "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
            signal_array_short_baselines_weighted = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_weighted_filename = "%s_signal_weighted.npy" % (model_vis_name_base)
            signal_array_short_baselines_weighted_Tb = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_weighted_Tb_filename = "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
            number_baselines_used_array = np.full((n_pol,n_chan),0.0)
            number_baselines_used_array_filename = "%s_number_baselines_used.npy" % (model_vis_name_base)
            sum_of_weights_all_baselines_array = np.full((n_pol,n_chan),np.nan)
            sum_of_weights_all_baselines_array_filename = "%s_sum_of_weights_all_baselines.npy" % (model_vis_name_base)
            sum_of_weights_short_baselines_array = np.full((n_pol,n_chan),np.nan)
            sum_of_weights_short_baselines_array_filename = "%s_sum_of_weights_short_baselines.npy" % (model_vis_name_base)
         
            #chan_sum_of_weights = 0.0
            #chan_vis_real_sum_weighted = 0.0
            for chan_index,freq_MHz in enumerate(freq_MHz_list):
               wavelength = 300./freq_MHz
               
               UU_wavelength_array = UU_m_array / wavelength
               VV_wavelength_array = VV_m_array / wavelength
               
               if pol=='X':
                  beam_image_name = "%smodel_%0.3f_MHz_xx.fits" % (beam_image_dir,freq_MHz)
               else:
                  beam_image_name = "%smodel_%0.3f_MHz_yy.fits" % (beam_image_dir,freq_MHz)
               
               #cmd = "cp %s%s ." % (beam_image_dir,beam_image_name)
               #print(cmd)
               #os.system(cmd)
               
               beam_plot_basename = beam_image_name.split('.')[0].split('/')[-1]
               
               with fits.open(beam_image_name) as beam_hdulist:
                  beam_image_data = beam_hdulist[0].data
                  beam_image_header = beam_hdulist[0].header
               
               #Fourier response:
               fits_name = beam_image_dir + beam_plot_basename + "_beam_fft2_real_shift.fits"  
               
               #cmd = "cp %s%s ." % (beam_image_dir,fits_name)
               #print(cmd)
               #os.system(cmd)
               
               beam_fft2_real_shift_hdu_list = fits.open(fits_name)
               beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data   
               beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
                    
               beam_image_length = beam_image_data.shape[0]
               beam_fft2_real_length = beam_fft2_real_shift.shape[0]
               fft_centre_padded = int(beam_fft2_real_length/2.)
               ##at phase centre angle is small sin(theta) = theta
               #sine projection so half beam image corresponds to sin(90 deg) = 1
               theta_step_rad = 1.0/(beam_image_length/2.)
               spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
               spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
               
               u_max = np.abs(spatial_frequencies_cols_fftshift[0])
               v_max = np.abs(spatial_frequencies_cols_fftshift[0])
               resolution = np.abs(spatial_frequencies_cols_fftshift[0] - spatial_frequencies_cols_fftshift[1])
               u_range = np.arange(-u_max,u_max,resolution)
               v_range = np.arange(-v_max,v_max,resolution)
                           
               number_baselines_used_sum = 0.0
               chan_sum_of_weights_all_baselines = 0.0
               chan_sum_of_weights_short_baselines = 0.0
               chan_vis_real_sum_weighted = 0.0
               chan_vis_real_sum = 0.0
               chan_vis_real_sum_all_baselines = 0.0
               chan_vis_abs_sum_all_baselines = 0.0
               #lst_vis_real_sum = 0.0
               #lst_vis_real_sum_weighted = 0.0
               chan_vis_used = 0.0
               for lst_index,lst in enumerate(lst_list):  
                  lst_deg = (float(lst)/24.)*360.         
                  print("lst_deg %s freq_MHz %s" % (lst_deg,freq_MHz))
                  #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
                  timestep_vis_real_sum = 0.
                  timestep_vis_real_sum_weighted = 0.
                  for timestep in timestep_array:
                     timestep_baselines_used = 0.
                     vis_start_index = int(timestep*n_baselines)
                     timestep_visibilities = visibilities[vis_start_index:int(vis_start_index+n_baselines),0,0,0,chan_index,0,:]
                     #print timestep_visibilities.shape
                     for timestep_visibility_index,timestep_visibility_real in enumerate(timestep_visibilities[:,0]):
                        visibility_index = vis_start_index + timestep_visibility_index
                        #only proceed with all this compute-expensive stuff if the uvdist is small i.e. less than 2 wavelengths
                        UU_wavelength = UU_wavelength_array[visibility_index]
                        VV_wavelength = VV_wavelength_array[visibility_index]
                        uvdist_wavelengths = np.sqrt(UU_wavelength**2 + VV_wavelength**2)
                        if uvdist_wavelengths < uvdist_wavelength_cutoff:
                           #print timestep_visibility_index
                           timestep_visibility_imag = timestep_visibilities[timestep_visibility_index,1]
                           ###visibility_weight = visibilities[visibility_index,0,0,0,2]
                           ##print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
                           complex_vis = np.complex(timestep_visibility_real,timestep_visibility_imag)
                           abs_vis = abs(complex_vis)
                           ##print "complex_vis %s" % complex_vis                   
                           #weighted: (this is equivalent to gridding), finding the value that the visibility would be it u,v=0 after convolving with a Fourier beam kernel:
                           
                           #all baselines:
                           chan_vis_real_sum_all_baselines += timestep_visibility_real
                           chan_vis_abs_sum_all_baselines += abs_vis
                           
                           #only short baselines
                           #Need to work this out now...(take stuff from compute weights function)
                           #uv_zero_weighting = weights_array[visibility_index,lst_index,chan_index,pol_index]
   
                           #print "uvdist %s wavelengths, smaller than threshold %s, proceeding" % (uvdist_wavelengths,uvdist_wavelength_cutoff)
                           #plot_basename = "U_%0.3f_V_%0.3f_pol_%s_%s_MHz" % (UU_wavelength,VV_wavelength,pol,freq_MHz)
                        
                           ###find the nearest U value in the Fourier beam response (need to interpolate to make this better.....) 
                           #use Jacks, quicker?
                           #bisection should be quickest
                           nearest_UU_index, nearest_UU_value = find_nearest(spatial_frequencies_cols_fftshift,UU_wavelength)
                           #nearest_UU_index = bisection(spatial_frequencies_cols_fftshift,UU_wavelength)
                           #nearest_UU_value = spatial_frequencies_cols_fftshift[nearest_UU_index]
                           #print "nearest_UU_index,nearest_UU_value %s,%s" % (nearest_UU_index, nearest_UU_value)
                           #print "nearest_UU_index %s" % (nearest_UU_index)
                           nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
                           #nearest_VV_index = bisection(spatial_frequencies_cols_fftshift,VV_wavelength)
                           #nearest_VV_value = spatial_frequencies_cols_fftshift[nearest_VV_index]
                           #print "nearest_VV_index,nearest_VV_value %s,%s" % (nearest_VV_index, nearest_VV_value)
                           #print "nearest_VV_index %s" % (nearest_VV_index)
                           U_offset = nearest_UU_value - UU_wavelength
                           V_offset = nearest_VV_value - VV_wavelength
   
                           #print "U offset is %s, V offset is %s" % (U_offset,V_offset)
                           
                           #u_ind,v_ind,u_off,v_off = find_closet_uv(u=UU_wavelength,v=VV_wavelength,u_range=u_range,v_range=v_range,resolution=resolution)
                           #print "jacks u_ind,v_in,u_off,v_off %s %s %s %s (offset in pix?)" % (u_ind,v_ind,u_off,v_off)
                           
                           #using mine:
                           UU_cheat_index = int((beam_fft2_real_length/2.)-(nearest_UU_index-(beam_fft2_real_length/2.)))
                           VV_cheat_index = int((beam_fft2_real_length/2.)-(nearest_VV_index-(beam_fft2_real_length/2.)))
                           
                           #using jacks:
                           #UU_cheat_index = int((beam_fft2_real_length/2.)-(u_ind-(beam_fft2_real_length/2.)))
                           #VV_cheat_index = int((beam_fft2_real_length/2.)-(v_ind-(beam_fft2_real_length/2.)))
   
                           #print "UU_cheat_index is %s, VV_cheat_index is %s " % (UU_cheat_index,VV_cheat_index)
                           
                           
                           #are these around the right way?
                           uv_zero_weighting = beam_fft2_real_shift_norm[VV_cheat_index,UU_cheat_index]
                           #uv_zero_weighting_symmetric = beam_fft2_real_shift_norm[nearest_VV_index,nearest_UU_index]
   
                           #print "uv_zero_weighting %s" % uv_zero_weighting
                           #print "uv_zero_weighting symmetric %s" % uv_zero_weighting_symmetric
                           
                           #all baselines signal add to weights array                        
                           chan_sum_of_weights_all_baselines += uv_zero_weighting
                           if uv_zero_weighting > zero_spacing_leakage_threshold:
                              #print "uv_zero_weighting %s for baseline %s"  % (uv_zero_weighting,visibility_index)
                              #complex_vis_weighted = uv_zero_weighting*complex_vis
                              #abs_vis_weighted = uv_zero_weighting*abs_vis
                              #abs_vis_weighted = np.abs(complex_vis_weighted)
                              #print "complex_vis_weighted %s" % complex_vis_weighted
                              #signal_array_weighted[pol_index,chan_index,0] += complex_vis_weighted
                              
                              #timestep_vis_real_sum += timestep_visibility_real
                              #timestep_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                              chan_vis_real_sum += timestep_visibility_real
                              chan_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                              #timestep_baselines_used += 1.
                              chan_vis_used += 1.
                              number_baselines_used_sum += 1.
                              chan_sum_of_weights_short_baselines += uv_zero_weighting
                              #if timestep_visibility_real<0:
                              #print "timestep_visibility_real %s at LST %s %s MHz" % (timestep_visibility_real,lst,freq_MHz)
                              #print "for baseline U %s V %s at visibility index %s and uv_zero_weighting %s" % (UU_m_array[visibility_index],VV_m_array[visibility_index],visibility_index,uv_zero_weighting)
                     #
                         
                     #timestep_vis_real_sum_weighted_norm = timestep_vis_real_sum_weighted*chan_sum_of_weights**2
                  
                          
                  #timestep_vis_real_average = timestep_vis_real_sum/n_timesteps
                  
                  #timestep_vis_real_sum_weighted_norm_average = timestep_vis_real_sum_weighted_norm/n_timesteps
                  #timestep_vis_real_average_baseline_norm = timestep_vis_real_average/timestep_baselines_used
                  #lst_vis_real_sum += timestep_vis_real_average
                  #lst_vis_real_sum_weighted += timestep_vis_real_sum_weighted_norm_average
                      
               #lst_vis_real_average = lst_vis_real_sum/len(lst_list)
               #lst_vis_real_average_weighted = lst_vis_real_sum_weighted/len(lst_list)
               
               #all_baselines:
               signal_array_all_baselines[pol_index,chan_index] = chan_vis_real_sum_all_baselines
               signal_array_all_baselines_abs[pol_index,chan_index] = chan_vis_abs_sum_all_baselines
               
               #weights:
               #if not np.isclose(chan_sum_of_weights,0.0):
               sum_of_weights_short_baselines_array[pol_index,chan_index] = chan_sum_of_weights_short_baselines
               sum_of_weights_all_baselines_array[pol_index,chan_index] = chan_sum_of_weights_all_baselines
               if chan_vis_used>0.:
                  #chan_vis_real_average = chan_vis_real_sum/chan_vis_used
                  print("chan_vis_real_sum")
                  print(chan_vis_real_sum)
                  #chan_vis_real_weighted_average = chan_vis_real_sum_weighted/chan_sum_of_weights
               
                  number_baselines_used_average =   number_baselines_used_sum / (n_timesteps*len(lst_list)) 
                  print("av number baselines used for chan %s is %s" % (chan_index,number_baselines_used_average))
               
               
                  #signal_array_short_baselines[pol_index,chan_index] = lst_vis_real_average
                  #signal_array_short_baselines_weighted[pol_index,chan_index] = lst_vis_real_average_weighted
                  
                  #using average?
                  #signal_array_short_baselines[pol_index,chan_index] = chan_vis_real_average
                  #signal_array_short_baselines_weighted[pol_index,chan_index] = chan_vis_real_weighted_average
               
                  signal_array_short_baselines[pol_index,chan_index] = chan_vis_real_sum
                  signal_array_short_baselines_weighted[pol_index,chan_index] = chan_vis_real_sum_weighted
                  
                  #number_baselines_used_array[pol_index,chan_index] = chan_vis_used            
                  number_baselines_used_array[pol_index,chan_index] = number_baselines_used_average 
               
                  #sum_of_weights += uv_zero_weighting
                  #no_baselines_weighted += 1
               
                  #print "sum_of_weights %s" % sum_of_weights
                  
                  #print abs(signal_array[pol_index,chan_index,0])
               
                  #normalise by dividing by the sum of the weights_array
                  #print "sum_of_weights %s at freq %s MHz with %s baselines" % (sum_of_weights,freq_MHz,no_baselines_weighted)
                  #signal_array_weighted_norm = signal_array_weighted/sum_of_weights
                    
                  #plt.clf()
                  #map_title="unweighted real vis vs freq x pol"
                  #plt.plot(freq_MHz_list,np.real(signal_array_unweighted[0,:,0]))
                  #plt.ylabel("unweighted sum vis real Jy")
                  #plt.xlabel("freq (MHz)")
                  #fig_name= plot_basename + "_real_vis_vs_freq_unweighted.png"
                  #figmap = plt.gcf()
                  #figmap.savefig(fig_name)
                  #print("saved %s" % fig_name)       
         
            #convert to brightness temp
            wavelength_array = 300./freq_MHz_array
            #print sum_of_weights_array
            jy_to_tb_all_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_all_baselines_array) 
            jy_to_tb_all_baselines_abs = (wavelength_array**2) / (2. * k * 1.0e26 * abs(sum_of_weights_all_baselines_array)) 
            jy_to_tb_short_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_short_baselines_array) 
            #jy_to_tb = (wavelength_array**2) / (2. * k * 1.0e26)
            signal_array_all_baselines_Tb = signal_array_all_baselines * jy_to_tb_all_baselines
            signal_array_short_baselines_Tb = signal_array_short_baselines * jy_to_tb_short_baselines
            signal_array_all_baselines_abs_Tb = signal_array_all_baselines_abs * jy_to_tb_all_baselines_abs
            signal_array_short_baselines_weighted_Tb = signal_array_short_baselines_weighted * jy_to_tb_short_baselines
               
            np.save(signal_array_short_baselines_filename,signal_array_short_baselines)
            np.save(signal_array_short_baselines_Tb_filename,signal_array_short_baselines_Tb)
            np.save(signal_array_all_baselines_filename,signal_array_all_baselines)
            np.save(signal_array_all_baselines_Tb_filename,signal_array_all_baselines_Tb)
            np.save(signal_array_all_baselines_filename_abs,signal_array_all_baselines_abs)
            np.save(signal_array_all_baselines_filename_abs_Tb,signal_array_all_baselines_abs_Tb)
            np.save(signal_array_short_baselines_weighted_filename,signal_array_short_baselines_weighted)
            np.save(signal_array_short_baselines_weighted_Tb_filename,signal_array_short_baselines_weighted_Tb)
            np.save(number_baselines_used_array_filename,number_baselines_used_array)
            np.save(sum_of_weights_all_baselines_array_filename,sum_of_weights_all_baselines_array)
            np.save(sum_of_weights_short_baselines_array_filename,sum_of_weights_short_baselines_array)
            
def extract_signal_from_sims_multi_coherent(lst_list,freq_MHz_list,pol_list,signal_type_list,outbase_name,sky_model,array_ant_locations_filename,array_label,n_arrays,min_spacing_m):
   min_spacing_cm = min_spacing_m * 100.
   concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
   concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
         if 'noise' in signal_type_list:
            concat_output_name_base_X += '_N'
            concat_output_name_base_Y += '_N'
         if 'diffuse' in signal_type_list:
            concat_output_name_base_X += '_D_%s' % sky_model
            concat_output_name_base_Y += '_D_%s' % sky_model
         if 'global' in signal_type_list:
            concat_output_name_base_X += '_G' 
            concat_output_name_base_Y += '_G'
         if 'gain_errors' in signal_type_list:
            concat_output_name_base_X += '_GE'
            concat_output_name_base_Y += '_GE'
         
         if do_image_and_subtr_from_simulated_data:
            uv_type_list = ['original','subtracted']
         else:
            uv_type_list = ['original']
         
         for uv_type in uv_type_list:
            if pol=='X':  
               if uv_type=='original':                      
                  model_vis_name_base = concat_output_name_base_X
               else:
                  model_vis_name_base = concat_output_name_base_X + '_sub'
            else:
               if uv_type=='original':
                  model_vis_name_base = concat_output_name_base_Y
               else:
                  model_vis_name_base = concat_output_name_base_Y + '_sub'
            model_vis_name_base += "_thresh_%0.2f" % (zero_spacing_leakage_threshold)

            signal_array_short_baselines_weighted = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_weighted_Tb = np.full((n_pol,n_chan),0.0)
            signal_array_short_baselines_weighted_Tb_filename = "%s_signal_weighted_Tb_multi.npy" % (model_vis_name_base)
            number_baselines_used_array = np.full((n_pol,n_chan),0.0)
            number_baselines_used_array_filename = "%s_number_baselines_used_multi.npy" % (model_vis_name_base)
            sum_of_weights_short_baselines_array = np.full((n_pol,n_chan),np.nan)
            sum_of_weights_short_baselines_array_filename = "%s_sum_of_weights_short_baselines_multi.npy" % (model_vis_name_base)
            
            for chan_index,freq_MHz in enumerate(freq_MHz_list):
               wavelength = 300./freq_MHz
               
               if pol=='X':
                  beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
               else:
                  beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
               
               beam_plot_basename = beam_image_name.split('.')[0]
               
               with fits.open(beam_image_name) as beam_hdulist:
                  beam_image_data = beam_hdulist[0].data
                  beam_image_header = beam_hdulist[0].header
               
               #Fourier response:
               fits_name= beam_plot_basename + "_beam_fft2_real_shift.fits"  
               
               beam_fft2_real_shift_hdu_list = fits.open(fits_name)
               beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data   
               beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
                    
               beam_image_length = beam_image_data.shape[0]
               beam_fft2_real_length = beam_fft2_real_shift.shape[0]
               fft_centre_padded = int(beam_fft2_real_length/2.)
               ##at phase centre angle is small sin(theta) = theta
               #sine projection so half beam image corresponds to sin(90 deg) = 1
               theta_step_rad = 1.0/(beam_image_length/2.)
               spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
               spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
               
               u_max = np.abs(spatial_frequencies_cols_fftshift[0])
               v_max = np.abs(spatial_frequencies_cols_fftshift[0])
               resolution = np.abs(spatial_frequencies_cols_fftshift[0] - spatial_frequencies_cols_fftshift[1])
               u_range = np.arange(-u_max,u_max,resolution)
               v_range = np.arange(-v_max,v_max,resolution)
                           
               number_baselines_used_sum = 0.0
               chan_sum_of_weights_all_baselines = 0.0
               chan_sum_of_weights_short_baselines = 0.0
               chan_vis_real_sum_weighted = 0.0
               chan_vis_real_sum = 0.0
               chan_vis_real_sum_all_baselines = 0.0
               chan_vis_abs_sum_all_baselines = 0.0
               chan_vis_used = 0.0
               
               
               for array_number in range(n_arrays):
                  array_label = 'array%03d_sep_%03d' % (array_number,min_spacing_cm)
                  concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
                  concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
                  if 'noise' in signal_type_list:
                     concat_output_name_base_X += '_N'
                     concat_output_name_base_Y += '_N'
                  if 'diffuse' in signal_type_list:
                     concat_output_name_base_X += '_D_%s' % sky_model
                     concat_output_name_base_Y += '_D_%s' % sky_model
                  if 'global' in signal_type_list:
                     concat_output_name_base_X += '_G' 
                     concat_output_name_base_Y += '_G'
                  if 'gain_errors' in signal_type_list:
                     concat_output_name_base_X += '_GE'
                     concat_output_name_base_Y += '_GE'
            
                  if pol=='X':  
                     if uv_type=='original':                      
                        uvfits_name = "%s_concat_lsts.uvfits" % concat_output_name_base_X
                     else:
                        uvfits_name = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_X
                  else:
                     if uv_type=='original':
                        uvfits_name = "%s_concat_lsts.uvfits" % concat_output_name_base_Y
                     else:
                        uvfits_name = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_Y
      
                  print(uvfits_name)           
                  hdulist = fits.open(uvfits_name)
                  #hdulist.info()
                  info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
                  #print info_string
                  uvtable = hdulist[0].data
                  uvtable_header = hdulist[0].header
                  visibilities = uvtable['DATA']
                  #print visibilities.shape
                  visibilities_shape = visibilities.shape
                  print("visibilities shape")
                  print(visibilities_shape)
                  
                  #get the UU and VV so we can check whether we are using short baselines
                  UU_s_array = uvtable['UU']
                  UU_m_array = UU_s_array * c   
                  VV_s_array = uvtable['VV']
                  VV_m_array = VV_s_array * c
                  
                  n_vis = visibilities.shape[0]
                  n_timesteps = n_vis/n_baselines
                  #print "n_timesteps %s " % n_timesteps
                  timestep_array = np.arange(0,n_timesteps,1)
            
                  UU_wavelength_array = UU_m_array / wavelength
                  VV_wavelength_array = VV_m_array / wavelength
                  
                  for lst_index,lst in enumerate(lst_list):  
                     lst_deg = (float(lst)/24.)*360.         
                     print("lst_deg %s freq_MHz %s" % (lst_deg,freq_MHz))
                     #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
                     timestep_vis_real_sum = 0.
                     timestep_vis_real_sum_weighted = 0.
                     for timestep in timestep_array:
                        timestep_baselines_used = 0.
                        vis_start_index = int(timestep*n_baselines)
                        timestep_visibilities = visibilities[vis_start_index:int(vis_start_index+n_baselines),0,0,0,chan_index,0,:]
                        #print timestep_visibilities.shape
                        for timestep_visibility_index,timestep_visibility_real in enumerate(timestep_visibilities[:,0]):
                           visibility_index = vis_start_index + timestep_visibility_index
                           #only proceed with all this compute-expensive stuff if the uvdist is small i.e. less than 2 wavelengths
                           UU_wavelength = UU_wavelength_array[visibility_index]
                           VV_wavelength = VV_wavelength_array[visibility_index]
                           uvdist_wavelengths = np.sqrt(UU_wavelength**2 + VV_wavelength**2)
                           if uvdist_wavelengths < uvdist_wavelength_cutoff:
                              #print timestep_visibility_index
                              timestep_visibility_imag = timestep_visibilities[timestep_visibility_index,1]
                              ###visibility_weight = visibilities[visibility_index,0,0,0,2]
                              ##print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
                              complex_vis = np.complex(timestep_visibility_real,timestep_visibility_imag)
                              abs_vis = abs(complex_vis)                 
   
                              nearest_UU_index, nearest_UU_value = find_nearest(spatial_frequencies_cols_fftshift,UU_wavelength)
                              nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
   
                              U_offset = nearest_UU_value - UU_wavelength
                              V_offset = nearest_VV_value - VV_wavelength
                                 
                              #using mine:
                              UU_cheat_index = int((beam_fft2_real_length/2.)-(nearest_UU_index-(beam_fft2_real_length/2.)))
                              VV_cheat_index = int((beam_fft2_real_length/2.)-(nearest_VV_index-(beam_fft2_real_length/2.)))                           
                              
                              #are these around the right way?
                              uv_zero_weighting = beam_fft2_real_shift_norm[VV_cheat_index,UU_cheat_index]
                              #uv_zero_weighting_symmetric = beam_fft2_real_shift_norm[nearest_VV_index,nearest_UU_index]
      
                              #print "uv_zero_weighting %s" % uv_zero_weighting
                              #print "uv_zero_weighting symmetric %s" % uv_zero_weighting_symmetric
                              
                              #all baselines signal add to weights array                        
                              chan_sum_of_weights_all_baselines += uv_zero_weighting
                              if uv_zero_weighting > zero_spacing_leakage_threshold:                              
                                 #timestep_vis_real_sum += timestep_visibility_real
                                 #timestep_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                                 chan_vis_real_sum += timestep_visibility_real
                                 chan_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                                 #timestep_baselines_used += 1.
                                 chan_vis_used += 1.
                                 number_baselines_used_sum += 1.
                                 chan_sum_of_weights_short_baselines += uv_zero_weighting
               
               #weights:
               sum_of_weights_short_baselines_array[pol_index,chan_index] = chan_sum_of_weights_short_baselines
               if chan_vis_used>0.:
                  #chan_vis_real_average = chan_vis_real_sum/chan_vis_used
                  print("chan_vis_real_sum")
                  print(chan_vis_real_sum)
                  #chan_vis_real_weighted_average = chan_vis_real_sum_weighted/chan_sum_of_weights
               
                  number_baselines_used_average =   number_baselines_used_sum / (n_timesteps*len(lst_list)) 
                  print("av number baselines used for chan %s is %s" % (chan_index,number_baselines_used_average))
               
                  signal_array_short_baselines_weighted[pol_index,chan_index] = chan_vis_real_sum_weighted
                  
                  #number_baselines_used_array[pol_index,chan_index] = chan_vis_used            
                  number_baselines_used_array[pol_index,chan_index] = number_baselines_used_average 
                  
         
            #convert to brightness temp
            wavelength_array = 300./freq_MHz_array
            jy_to_tb_short_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_short_baselines_array) 
            #jy_to_tb = (wavelength_array**2) / (2. * k * 1.0e26)

            signal_array_short_baselines_weighted_Tb = signal_array_short_baselines_weighted * jy_to_tb_short_baselines
               
            np.save(signal_array_short_baselines_weighted_Tb_filename,signal_array_short_baselines_weighted_Tb)
            np.save(number_baselines_used_array_filename,number_baselines_used_array)
            np.save(sum_of_weights_short_baselines_array_filename,sum_of_weights_short_baselines_array)     

def extract_signal_from_assassin(lst_list,freq_MHz_list,pol_list,signal_type_list,sky_model,outbase_name,n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length_m,zero_spacing_leakage_threshold):
   n_baselines = 1 #(by definition of the assassin design, each uvfits is for just one baseline)
   final_lst = lst_list[-1]
   final_lst_deg = (float(final_lst)/24.)*360. 
   final_freq = freq_MHz_list[-1] 
   freq_MHz_array = np.asarray(freq_MHz_list)
   radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
   zero_spacing_leakage_threshold_pc = zero_spacing_leakage_threshold * 100.
   max_arm_length_cm = max_arm_length_m * 100.
   min_arm_length_cm = min_arm_length_m * 100.
   for pol in pol_list:
      final_output_name_base = "assassin_%s_%s_%03d_%03d_%s_LST_%03d_thresh_%02d_pc" % (n_ants_per_m_of_circumference,n_circles,max_arm_length_cm,min_arm_length_cm,pol,final_lst_deg,zero_spacing_leakage_threshold_pc)
      signal_array_short_baselines_weighted = np.full(n_chan,0.0)
      signal_array_short_baselines_weighted_Tb = np.full(n_chan,0.0)
      number_baselines_used_array = np.full(n_chan,0.0)
      sum_of_weights_short_baselines_array = np.full(n_chan,np.nan)
      signal_array_short_baselines_weighted_Tb_filename = "%s_signal_Tb.npy" % (final_output_name_base)
      number_baselines_used_array_filename = "%s_n_baselines_used.npy" % (final_output_name_base)
      sum_of_weights_short_baselines_array_filename = "%s_sum_of_weights.npy" % (final_output_name_base)
      for chan_index,freq_MHz in enumerate(freq_MHz_list):
         wavelength = 300./freq_MHz
         
         chan_vis_real_sum_weighted = 0.0
         chan_vis_used = 0.0
         
         if pol=='X':
            beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
         else:
            beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
         
         beam_plot_basename = beam_image_name.split('.')[0]
         
         with fits.open(beam_image_name) as beam_hdulist:
            beam_image_data = beam_hdulist[0].data
            beam_image_header = beam_hdulist[0].header
         
         #Fourier response:
         fits_name= beam_plot_basename + "_beam_fft2_real_shift.fits"  
         
         beam_fft2_real_shift_hdu_list = fits.open(fits_name)
         beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data   
         beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
              
         beam_image_length = beam_image_data.shape[0]
         beam_fft2_real_length = beam_fft2_real_shift.shape[0]
         fft_centre_padded = int(beam_fft2_real_length/2.)
         ##at phase centre angle is small sin(theta) = theta
         #sine projection so half beam image corresponds to sin(90 deg) = 1
         theta_step_rad = 1.0/(beam_image_length/2.)
         spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
         spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
         
         u_max = np.abs(spatial_frequencies_cols_fftshift[0])
         v_max = np.abs(spatial_frequencies_cols_fftshift[0])
         resolution = np.abs(spatial_frequencies_cols_fftshift[0] - spatial_frequencies_cols_fftshift[1])
         u_range = np.arange(-u_max,u_max,resolution)
         v_range = np.arange(-v_max,v_max,resolution)
            
         number_baselines_used_sum = 0.0
         chan_sum_of_weights_short_baselines = 0.0
    
         #then for each sub-array 
         for circle_number in range(0,n_circles):
            #circle 1 is the smallest
            radius = (min_arm_length_m + circle_number * radial_spacing) 
            diameter = radius * 2.
            diameter_cm = int(diameter*100.)
            #work out circumference 
            circumference = math.pi * diameter
            #print diameter
            n_angles = int(round(circumference * n_ants_per_m_of_circumference))
            angle_increment = (2.*math.pi)/n_angles
            #(remember only need half of them!)
            angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
            #angle_array_deg = angle_array_rad / math.pi * 180.
            #print angle_array_deg
            for angle_rad in angle_array_rad:
               angle_deg = int(angle_rad/np.pi*180.)
               sub_array_string = "sub_%03d_%03d" % (diameter_cm,angle_deg)
               
               model_vis_name_base = "%s_LST_%03d_%s_%s_MHz" % (sub_array_string,final_lst_deg,pol,int(final_freq))
               
               if 'noise' in signal_type_list:
                  model_vis_name_base += '_N'
               if 'diffuse' in signal_type_list:
                  model_vis_name_base += '_D_%s' % sky_model
               if 'global' in signal_type_list:
                  model_vis_name_base += '_G' 
               if 'gain_errors' in signal_type_list:
                  model_vis_name_base += '_GE'
                 
               uvfits_name = "%s_concat_lsts.uvfits" % model_vis_name_base
               print(uvfits_name)
                          
               hdulist = fits.open(uvfits_name)
               #hdulist.info()
               info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
               #print info_string
               uvtable = hdulist[0].data
               uvtable_header = hdulist[0].header
               visibilities = uvtable['DATA']
               #print visibilities.shape
               visibilities_shape = visibilities.shape
               print("visibilities shape")
               print(visibilities_shape)
               
               #get the UU and VV so we can check whether we are using short baselines
               UU_s_array = uvtable['UU']
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV']
               VV_m_array = VV_s_array * c
               
               n_vis = visibilities.shape[0]
               n_timesteps = n_vis/n_baselines
               #print "n_timesteps %s " % n_timesteps
               timestep_array = np.arange(0,n_timesteps,1)
         
               UU_wavelength_array = UU_m_array / wavelength
               VV_wavelength_array = VV_m_array / wavelength
               
               for lst_index,lst in enumerate(lst_list):  
                  #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
                  timestep_vis_real_sum = 0.
                  timestep_vis_real_sum_weighted = 0.
                  for timestep in timestep_array:
                     timestep_baselines_used = 0.
                     vis_start_index = int(timestep*n_baselines)
                     timestep_visibilities = visibilities[vis_start_index:int(vis_start_index+n_baselines),0,0,0,chan_index,0,:]
                     #print timestep_visibilities.shape
                     for timestep_visibility_index,timestep_visibility_real in enumerate(timestep_visibilities[:,0]):
                        visibility_index = vis_start_index + timestep_visibility_index
                        #only proceed with all this compute-expensive stuff if the uvdist is small i.e. less than 2 wavelengths
                        UU_wavelength = UU_wavelength_array[visibility_index]
                        VV_wavelength = VV_wavelength_array[visibility_index]
                        uvdist_wavelengths = np.sqrt(UU_wavelength**2 + VV_wavelength**2)
                        if uvdist_wavelengths < uvdist_wavelength_cutoff:
                           #print timestep_visibility_index
                           timestep_visibility_imag = timestep_visibilities[timestep_visibility_index,1]
                           ###visibility_weight = visibilities[visibility_index,0,0,0,2]
                           ##print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
                           complex_vis = np.complex(timestep_visibility_real,timestep_visibility_imag)
                           abs_vis = abs(complex_vis)                 
   
                           nearest_UU_index, nearest_UU_value = find_nearest(spatial_frequencies_cols_fftshift,UU_wavelength)
                           nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
   
                           U_offset = nearest_UU_value - UU_wavelength
                           V_offset = nearest_VV_value - VV_wavelength
                              
                           #using mine:
                           UU_cheat_index = int((beam_fft2_real_length/2.)-(nearest_UU_index-(beam_fft2_real_length/2.)))
                           VV_cheat_index = int((beam_fft2_real_length/2.)-(nearest_VV_index-(beam_fft2_real_length/2.)))                           
                           
                           #are these around the right way?
                           uv_zero_weighting = beam_fft2_real_shift_norm[VV_cheat_index,UU_cheat_index]
                           #uv_zero_weighting_symmetric = beam_fft2_real_shift_norm[nearest_VV_index,nearest_UU_index]
      
                           #print "uv_zero_weighting %s" % uv_zero_weighting
                           #print "uv_zero_weighting symmetric %s" % uv_zero_weighting_symmetric
                                                 
                           if uv_zero_weighting > zero_spacing_leakage_threshold:                              
                              chan_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                              chan_vis_used += 1.
                              number_baselines_used_sum += 1.
                              chan_sum_of_weights_short_baselines += uv_zero_weighting
               
         #weights:
         sum_of_weights_short_baselines_array[chan_index] = chan_sum_of_weights_short_baselines
         if chan_vis_used>0.:
            number_baselines_used_average =   number_baselines_used_sum / (n_timesteps*len(lst_list)) 
            print("av number baselines used for chan %s is %s" % (chan_index,number_baselines_used_average))
         
            signal_array_short_baselines_weighted[chan_index] = chan_vis_real_sum_weighted
                       
            number_baselines_used_array[chan_index] = number_baselines_used_average 
                  
         
      #convert to brightness temp
      wavelength_array = 300./freq_MHz_array
      jy_to_tb_short_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_short_baselines_array) 
      #jy_to_tb = (wavelength_array**2) / (2. * k * 1.0e26)

      signal_array_short_baselines_weighted_Tb = signal_array_short_baselines_weighted * jy_to_tb_short_baselines
         
      np.save(signal_array_short_baselines_weighted_Tb_filename,signal_array_short_baselines_weighted_Tb)
      np.save(number_baselines_used_array_filename,number_baselines_used_array)
      np.save(sum_of_weights_short_baselines_array_filename,sum_of_weights_short_baselines_array)     
            
            
            
            
def extract_signal_from_multiple(uvfits_list_filename):
   #read the txt file containing list of uvfits files
   uvfits_filename_list = []
   with open(uvfits_list_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      uvfits_filename_list.append(line.strip())
   for uvfits_filename in uvfits_filename_list:
      #print uvfits_filename           
      hdulist = fits.open(uvfits_filename)
      #hdulist.info()
      info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print info_string
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      #print uvtable_header
      #get info from header
      n_chan = int(uvtable_header['NAXIS4'])
      ref_pix = int(uvtable_header['CRPIX4'])
      centre_pix = int(n_chan/2)
      cdelt4 = int(uvtable_header['CDELT4'])
      crval4 = float(uvtable_header['CRVAL4'])
      #print "number of chans %s, ref pix %s, delt %s, crval %s" % (n_chan,ref_pix,cdelt4,crval4)
      if ref_pix == 1:
         #case for sims
         start_freq = crval4
      elif (ref_pix==centre_pix):
         #NEED TO TEST THIS WITH DATA!
         start_freq = crval4 - (centre_pix*cdelt4)
      else:
         print("invalid ref pixel for freq")
         sys.exit()
      
      start_freq_MHz = start_freq/1000000.
      chan_width_MHz = float(cdelt4/1000000.)
      end_freq_MHz = start_freq_MHz + (chan_width_MHz * n_chan)
      
      visibilities = uvtable['DATA']
      #print visibilities.shape
      visibilities_shape = visibilities.shape
      print("visibilities shape")
      print(isibilities_shape)
      
      n_pol = visibilities_shape[-2]
      
      print("number of chans %s, start freq %s MHz, chan width %s MHz, npol %s" % (n_chan,start_freq_MHz,chan_width_MHz,n_pol))
     
      freq_MHz_list = np.arange(start_freq_MHz,end_freq_MHz,chan_width_MHz)
      #print freq_MHz_list
   
      #get the UU and VV so we can check whether we are using short baselines
      UU_s_array = uvtable['UU']
      UU_m_array = UU_s_array * c   
      VV_s_array = uvtable['VV']
      VV_m_array = VV_s_array * c
            
      n_vis = visibilities.shape[0]
      n_timesteps = n_vis/n_baselines
      #print "n_timesteps %s " % n_timesteps
      timestep_array = np.arange(0,n_timesteps,1)
            
      for pol in range(0,n_pol):
         model_vis_name_base = uvfits_list_filename.split('.')[0]
         model_vis_name_base += "_thresh_%0.2f_pol_%s" % (zero_spacing_leakage_threshold,int(pol))
         #print model_vis_name_base
         signal_array_short_baselines = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_filename = "%s_signal.npy" % (model_vis_name_base)
         signal_array_short_baselines_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_Tb_filename = "%s_signal_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename = "%s_signal_all_baselines.npy" % (model_vis_name_base)
         signal_array_all_baselines_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_Tb_filename = "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines_abs = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename_abs = "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
         signal_array_all_baselines_abs_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename_abs_Tb = "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
         signal_array_short_baselines_weighted = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_weighted_filename = "%s_signal_weighted.npy" % (model_vis_name_base)
         signal_array_short_baselines_weighted_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_weighted_Tb_filename = "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
         number_baselines_used_array = np.full((n_pol,n_chan),0.0)
         number_baselines_used_array_filename = "%s_number_baselines_used.npy" % (model_vis_name_base)
         sum_of_weights_all_baselines_array = np.full((n_pol,n_chan),np.nan)
         sum_of_weights_all_baselines_array_filename = "%s_sum_of_weights_all_baselines.npy" % (model_vis_name_base)
         sum_of_weights_short_baselines_array = np.full((n_pol,n_chan),np.nan)
         sum_of_weights_short_baselines_array_filename = "%s_sum_of_weights_short_baselines.npy" % (model_vis_name_base)
         
         #chan_sum_of_weights = 0.0
         #chan_vis_real_sum_weighted = 0.0
         for chan_index,freq_MHz in enumerate(freq_MHz_list):
            wavelength = 300./freq_MHz
            
            UU_wavelength_array = UU_m_array / wavelength
            VV_wavelength_array = VV_m_array / wavelength
            
            if pol==0:
               beam_image_name = "model_%0.3f_MHz_xx.fits" % freq_MHz
            else:
               beam_image_name = "model_%0.3f_MHz_yy.fits" % freq_MHz
            
            beam_plot_basename = beam_image_name.split('.')[0]
            
            with fits.open(beam_image_name) as beam_hdulist:
               beam_image_data = beam_hdulist[0].data
               beam_image_header = beam_hdulist[0].header
            
            #Fourier response:
            fits_name= beam_plot_basename + "_beam_fft2_real_shift.fits"  
            
            beam_fft2_real_shift_hdu_list = fits.open(fits_name)
            beam_fft2_real_shift = beam_fft2_real_shift_hdu_list[0].data   
            beam_fft2_real_shift_norm = beam_fft2_real_shift/beam_fft2_real_shift.max()
                 
            beam_image_length = beam_image_data.shape[0]
            beam_fft2_real_length = beam_fft2_real_shift.shape[0]
            fft_centre_padded = int(beam_fft2_real_length/2.)
            ##at phase centre angle is small sin(theta) = theta
            #sine projection so half beam image corresponds to sin(90 deg) = 1
            theta_step_rad = 1.0/(beam_image_length/2.)
            spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_real_shift.shape[1],d=theta_step_rad)
            spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
            
            u_max = np.abs(spatial_frequencies_cols_fftshift[0])
            v_max = np.abs(spatial_frequencies_cols_fftshift[0])
            resolution = np.abs(spatial_frequencies_cols_fftshift[0] - spatial_frequencies_cols_fftshift[1])
            u_range = np.arange(-u_max,u_max,resolution)
            v_range = np.arange(-v_max,v_max,resolution)
                        
            number_baselines_used_sum = 0.0
            chan_sum_of_weights_all_baselines = 0.0
            chan_sum_of_weights_short_baselines = 0.0
            chan_vis_real_sum_weighted = 0.0
            chan_vis_real_sum = 0.0
            chan_vis_real_sum_all_baselines = 0.0
            chan_vis_abs_sum_all_baselines = 0.0
            #lst_vis_real_sum = 0.0
            #lst_vis_real_sum_weighted = 0.0
            chan_vis_used = 0.0
            
            #sys.exit()
            
            for lst_index,lst in enumerate(lst_list):  
               lst_deg = (float(lst)/24.)*360.         
               #print "lst_deg %s freq_MHz %s" % (lst_deg,freq_MHz)
               #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
               timestep_vis_real_sum = 0.
               timestep_vis_real_sum_weighted = 0.
               for timestep in timestep_array:
                  timestep_baselines_used = 0.
                  vis_start_index = int(timestep*n_baselines)
                  timestep_visibilities = visibilities[vis_start_index:int(vis_start_index+n_baselines),0,0,0,chan_index,0,:]
                  #print timestep_visibilities.shape
                  for timestep_visibility_index,timestep_visibility_real in enumerate(timestep_visibilities[:,0]):
                     visibility_index = vis_start_index + timestep_visibility_index
                     #only proceed with all this compute-expensive stuff if the uvdist is small i.e. less than 2 wavelengths
                     UU_wavelength = UU_wavelength_array[visibility_index]
                     VV_wavelength = VV_wavelength_array[visibility_index]
                     uvdist_wavelengths = np.sqrt(UU_wavelength**2 + VV_wavelength**2)
                     if uvdist_wavelengths < uvdist_wavelength_cutoff:
                        #print timestep_visibility_index
                        timestep_visibility_imag = timestep_visibilities[timestep_visibility_index,1]
                        ###visibility_weight = visibilities[visibility_index,0,0,0,2]
                        ##print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
                        complex_vis = np.complex(timestep_visibility_real,timestep_visibility_imag)
                        abs_vis = abs(complex_vis)
                        ##print "complex_vis %s" % complex_vis                   
                        #weighted: (this is equivalent to gridding), finding the value that the visibility would be it u,v=0 after convolving with a Fourier beam kernel:
                        
                        #all baselines:
                        chan_vis_real_sum_all_baselines += timestep_visibility_real
                        chan_vis_abs_sum_all_baselines += abs_vis
                        
                        #only short baselines
                        #Need to work this out now...(take stuff from compute weights function)
                        #uv_zero_weighting = weights_array[visibility_index,lst_index,chan_index,pol_index]
   
                        #print "uvdist %s wavelengths, smaller than threshold %s, proceeding" % (uvdist_wavelengths,uvdist_wavelength_cutoff)
                        #plot_basename = "U_%0.3f_V_%0.3f_pol_%s_%s_MHz" % (UU_wavelength,VV_wavelength,pol,freq_MHz)
                     
                        ###find the nearest U value in the Fourier beam response (need to interpolate to make this better.....) 
                        #use Jacks, quicker?
                        #bisection should be quickest
                        nearest_UU_index, nearest_UU_value = find_nearest(spatial_frequencies_cols_fftshift,UU_wavelength)
                        #nearest_UU_index = bisection(spatial_frequencies_cols_fftshift,UU_wavelength)
                        #nearest_UU_value = spatial_frequencies_cols_fftshift[nearest_UU_index]
                        #print "nearest_UU_index,nearest_UU_value %s,%s" % (nearest_UU_index, nearest_UU_value)
                        #print "nearest_UU_index %s" % (nearest_UU_index)
                        nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
                        #nearest_VV_index = bisection(spatial_frequencies_cols_fftshift,VV_wavelength)
                        #nearest_VV_value = spatial_frequencies_cols_fftshift[nearest_VV_index]
                        #print "nearest_VV_index,nearest_VV_value %s,%s" % (nearest_VV_index, nearest_VV_value)
                        #print "nearest_VV_index %s" % (nearest_VV_index)
                        U_offset = nearest_UU_value - UU_wavelength
                        V_offset = nearest_VV_value - VV_wavelength
   
                        #print "U offset is %s, V offset is %s" % (U_offset,V_offset)
                        
                        #u_ind,v_ind,u_off,v_off = find_closet_uv(u=UU_wavelength,v=VV_wavelength,u_range=u_range,v_range=v_range,resolution=resolution)
                        #print "jacks u_ind,v_in,u_off,v_off %s %s %s %s (offset in pix?)" % (u_ind,v_ind,u_off,v_off)
                        
                        #using mine:
                        UU_cheat_index = int((beam_fft2_real_length/2.)-(nearest_UU_index-(beam_fft2_real_length/2.)))
                        VV_cheat_index = int((beam_fft2_real_length/2.)-(nearest_VV_index-(beam_fft2_real_length/2.)))
                        
                        #using jacks:
                        #UU_cheat_index = int((beam_fft2_real_length/2.)-(u_ind-(beam_fft2_real_length/2.)))
                        #VV_cheat_index = int((beam_fft2_real_length/2.)-(v_ind-(beam_fft2_real_length/2.)))
   
                        #print "UU_cheat_index is %s, VV_cheat_index is %s " % (UU_cheat_index,VV_cheat_index)
                        
                        
                        #are these around the right way?
                        uv_zero_weighting = beam_fft2_real_shift_norm[VV_cheat_index,UU_cheat_index]
                        #uv_zero_weighting_symmetric = beam_fft2_real_shift_norm[nearest_VV_index,nearest_UU_index]
   
                        #print "uv_zero_weighting %s" % uv_zero_weighting
                        #print "uv_zero_weighting symmetric %s" % uv_zero_weighting_symmetric
                        
                        #all baselines signal add to weights array                        
                        chan_sum_of_weights_all_baselines += uv_zero_weighting
                        if uv_zero_weighting > zero_spacing_leakage_threshold:
                           #print "uv_zero_weighting %s for baseline %s"  % (uv_zero_weighting,visibility_index)
                           #complex_vis_weighted = uv_zero_weighting*complex_vis
                           #abs_vis_weighted = uv_zero_weighting*abs_vis
                           #abs_vis_weighted = np.abs(complex_vis_weighted)
                           #print "complex_vis_weighted %s" % complex_vis_weighted
                           #signal_array_weighted[pol_index,chan_index,0] += complex_vis_weighted
                           
                           #timestep_vis_real_sum += timestep_visibility_real
                           #timestep_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                           chan_vis_real_sum += timestep_visibility_real
                           chan_vis_real_sum_weighted += timestep_visibility_real * uv_zero_weighting
                           #timestep_baselines_used += 1.
                           chan_vis_used += 1.
                           number_baselines_used_sum += 1.
                           chan_sum_of_weights_short_baselines += uv_zero_weighting
                           #if timestep_visibility_real<0:
                           #print "timestep_visibility_real %s at LST %s %s MHz" % (timestep_visibility_real,lst,freq_MHz)
                           #print "for baseline U %s V %s at visibility index %s and uv_zero_weighting %s" % (UU_m_array[visibility_index],VV_m_array[visibility_index],visibility_index,uv_zero_weighting)
                  #
                      
                  #timestep_vis_real_sum_weighted_norm = timestep_vis_real_sum_weighted*chan_sum_of_weights**2
               
                       
               #timestep_vis_real_average = timestep_vis_real_sum/n_timesteps
               
               #timestep_vis_real_sum_weighted_norm_average = timestep_vis_real_sum_weighted_norm/n_timesteps
               #timestep_vis_real_average_baseline_norm = timestep_vis_real_average/timestep_baselines_used
               #lst_vis_real_sum += timestep_vis_real_average
               #lst_vis_real_sum_weighted += timestep_vis_real_sum_weighted_norm_average
                   
            #lst_vis_real_average = lst_vis_real_sum/len(lst_list)
            #lst_vis_real_average_weighted = lst_vis_real_sum_weighted/len(lst_list)
            
            #all_baselines:
            signal_array_all_baselines[pol_index,chan_index] = chan_vis_real_sum_all_baselines
            signal_array_all_baselines_abs[pol_index,chan_index] = chan_vis_abs_sum_all_baselines
            
            #weights:
            #if not np.isclose(chan_sum_of_weights,0.0):
            sum_of_weights_short_baselines_array[pol_index,chan_index] = chan_sum_of_weights_short_baselines
            sum_of_weights_all_baselines_array[pol_index,chan_index] = chan_sum_of_weights_all_baselines
            if chan_vis_used>0.:
               #chan_vis_real_average = chan_vis_real_sum/chan_vis_used
               #print "chan_vis_real_sum"
               #print chan_vis_real_sum
               #chan_vis_real_weighted_average = chan_vis_real_sum_weighted/chan_sum_of_weights
            
               number_baselines_used_average =   number_baselines_used_sum / (n_timesteps*len(lst_list)) 
               #print "av number baselines used for chan %s is %s" % (chan_index,number_baselines_used_average)
            
            
               #signal_array_short_baselines[pol_index,chan_index] = lst_vis_real_average
               #signal_array_short_baselines_weighted[pol_index,chan_index] = lst_vis_real_average_weighted
               
               #using average?
               #signal_array_short_baselines[pol_index,chan_index] = chan_vis_real_average
               #signal_array_short_baselines_weighted[pol_index,chan_index] = chan_vis_real_weighted_average
            
               signal_array_short_baselines[pol_index,chan_index] = chan_vis_real_sum
               signal_array_short_baselines_weighted[pol_index,chan_index] = chan_vis_real_sum_weighted
               
               #number_baselines_used_array[pol_index,chan_index] = chan_vis_used            
               number_baselines_used_array[pol_index,chan_index] = number_baselines_used_average 
            
               #sum_of_weights += uv_zero_weighting
               #no_baselines_weighted += 1
            
               #print "sum_of_weights %s" % sum_of_weights
               
               #print abs(signal_array[pol_index,chan_index,0])
            
               #normalise by dividing by the sum of the weights_array
               #print "sum_of_weights %s at freq %s MHz with %s baselines" % (sum_of_weights,freq_MHz,no_baselines_weighted)
               #signal_array_weighted_norm = signal_array_weighted/sum_of_weights
                 
               #plt.clf()
               #map_title="unweighted real vis vs freq x pol"
               #plt.plot(freq_MHz_list,np.real(signal_array_unweighted[0,:,0]))
               #plt.ylabel("unweighted sum vis real Jy")
               #plt.xlabel("freq (MHz)")
               #fig_name= plot_basename + "_real_vis_vs_freq_unweighted.png"
               #figmap = plt.gcf()
               #figmap.savefig(fig_name)
               #print("saved %s" % fig_name)       
         
         #convert to brightness temp
         wavelength_array = 300./freq_MHz_array
         #print sum_of_weights_array
         jy_to_tb_all_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_all_baselines_array) 
         jy_to_tb_all_baselines_abs = (wavelength_array**2) / (2. * k * 1.0e26 * abs(sum_of_weights_all_baselines_array)) 
         jy_to_tb_short_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_short_baselines_array) 
         #jy_to_tb = (wavelength_array**2) / (2. * k * 1.0e26)
         signal_array_all_baselines_Tb = signal_array_all_baselines * jy_to_tb_all_baselines
         signal_array_short_baselines_Tb = signal_array_short_baselines * jy_to_tb_short_baselines
         signal_array_all_baselines_abs_Tb = signal_array_all_baselines_abs * jy_to_tb_all_baselines_abs
         signal_array_short_baselines_weighted_Tb = signal_array_short_baselines_weighted * jy_to_tb_short_baselines
            
         np.save(signal_array_short_baselines_filename,signal_array_short_baselines)
         np.save(signal_array_short_baselines_Tb_filename,signal_array_short_baselines_Tb)
         np.save(signal_array_all_baselines_filename,signal_array_all_baselines)
         np.save(signal_array_all_baselines_Tb_filename,signal_array_all_baselines_Tb)
         np.save(signal_array_all_baselines_filename_abs,signal_array_all_baselines_abs)
         np.save(signal_array_all_baselines_filename_abs_Tb,signal_array_all_baselines_abs_Tb)
         np.save(signal_array_short_baselines_weighted_filename,signal_array_short_baselines_weighted)
         np.save(signal_array_short_baselines_weighted_Tb_filename,signal_array_short_baselines_weighted_Tb)
         np.save(number_baselines_used_array_filename,number_baselines_used_array)
         np.save(sum_of_weights_all_baselines_array_filename,sum_of_weights_all_baselines_array)
         np.save(sum_of_weights_short_baselines_array_filename,sum_of_weights_short_baselines_array)
            
            
def plot_from_uvfits(uvfits_name, freq_MHz):
               
   plot_basename = uvfits_name.split('.')[0]
   
   freq_MHz = float(freq_MHz)
   wavelength = 300.0 / freq_MHz
   
   hdulist = fits.open(uvfits_name)
   hdulist.info()
   info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
   #print info_string
   uvtable = hdulist[0].data
   uvtable_header = hdulist[0].header
   #print uvtable_header
   
   #super_threshold_indices = a > thresh
   #a[super_threshold_indices] = 0
   
   #print uvtable['DATA'].shape
   real_visibilities_at_freq_MHz1 = uvtable['DATA'][:,0,0,0,0,0]
   abs_vis = abs(real_visibilities_at_freq_MHz1)
   log_abs_vis = np.log10(abs_vis)
   UU_s = uvtable['UU']
   UU_m = UU_s * c
   VV_s = uvtable['VV']
   VV_m = VV_s * c
   
   uvdist_m = np.sqrt(UU_m**2 + VV_m**2)
   
   #print "min abs uvdist_m"
   #print np.abs(uvdist_m).min()
   #print "max abs uvdist_m"
   #print np.abs(uvdist_m).max()
   
   
   UU_lambda = UU_m/wavelength
   uvdist_lambda = uvdist_m/wavelength
   plt.clf()
   map_title="abs vis vs UU"
   plt.plot(UU_lambda,log_abs_vis,linestyle='none',marker='.')
   plt.ylabel("log abs(vis)")
   plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "log_abs_vis_vs_uu.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   plt.clf()
   map_title="real vis vs UU"
   plt.plot(uvdist_lambda,real_visibilities_at_freq_MHz1,linestyle='none',marker='.')
   plt.ylabel("real vis Jy")
   plt.xlabel("uvdist (lambda)")
   fig_name= plot_basename + "real_vis_vs_uvdist.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   #print UU_lambda.shape
   uvdist_lambda_short_indices = np.where(abs(uvdist_lambda) < 2.)
   uvdist_lambda_short = uvdist_lambda[uvdist_lambda_short_indices]
   real_visibilities_at_freq_MHz1_short = real_visibilities_at_freq_MHz1[uvdist_lambda_short_indices]
   #print real_visibilities_at_freq_MHz1_short.shape
   plt.clf()
   map_title="real vis vs uvdist zoom"
   plt.plot(uvdist_lambda_short,real_visibilities_at_freq_MHz1_short,linestyle='none',marker='.')
   plt.ylabel("real vis Jy")
   plt.xlabel("uvdist (lambda)")
   fig_name= plot_basename + "real_vis_vs_uvdist_zoom.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   
   #for signal extraction only interested in short baselines (work out what this threshold should be based on U=0 leakage value later)
   threshold = 4.0
   threshold_indices = UU_m > threshold
   #short baselines only - as a function of frequency

   #for now, just take the max value of the visibilities at each freq
   #max_vis_vs_freq_list = []
   #for chan_index in np.arange(0,150):
   #   max_vis = np.max(uvtable['DATA'][:,0,0,0,chan_index,0,0])
   #   max_vis_vs_freq_list.append(max_vis)
   #max_vis_vs_freq_array = np.asarray(max_vis_vs_freq_list)
   #abs_max_vis_vs_freq_array = abs(max_vis_vs_freq_array)
   #log_abs_max_vis_vs_freq_array = np.log10(max_vis_vs_freq_array)
   #plt.clf()
   #map_title="abs vis vs freq for short baselines"
   #plt.plot(freq_MHz_list,log_abs_max_vis_vs_freq_array,linestyle='none',marker='.')
   #plt.ylabel("log abs(max vis)")
   #plt.xlabel("freq (MHz)")
   #fig_name= plot_basename + "log_abs_max_vis_vs_freq.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name)   
 
   #plt.clf()
   #map_title="log abs vis vs freq"
   #plt.plot(freq_MHz_array,abs_vis,linestyle='none',marker='.')
   #fig_name="abs_vis_vs_uu.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name)


def concat_uvfits(uvfits_list,outname,delete=False):
   tmp_uv = UVData()
   first_filename = uvfits_list[0]
   
   ##checking to see why there are three timesteps... had to get harange right
   #hdulist = fits.open(first_filename)
   ##hdulist.info()
   ##info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
   ##print info_string
   #uvtable = hdulist[0].data
   #print first_filename
   #print uvtable.shape
   #n_timesteps = uvtable.shape[0]/n_baselines
   #print "n_timesteps %s " % n_timesteps

   
   fits.setval(first_filename, 'TELESCOP', value='MWA' )
   fits.setval(first_filename, 'OBJECT', value='zenith' )
   fits.setval(first_filename, 'CDELT4', value=1.000E+06)
   fits.setval(first_filename, 'DATE', value=1.0000)
   fits.setval(first_filename, '_DATE', value=0.5000)
   fits.setval(first_filename, 'INTTIM', value=240)
   #print first_filename
   tmp_uv.read_uvfits(first_filename)
   #run_check=False,check_extra=False,run_check_acceptability=False
   for uvfits_name in uvfits_list[1:]: 
      #open the uvfits file and put in MWA as telescope
      fits.setval(uvfits_name, 'TELESCOP', value='MWA')
      fits.setval(first_filename, 'OBJECT', value='zenith' )
      fits.setval(uvfits_name, 'CDELT4', value=1.000E+06)
      fits.setval(uvfits_name, 'DATE', value=1.0000)
      fits.setval(uvfits_name, '_DATE', value=0.5000)
      fits.setval(uvfits_name, 'INTTIM', value=240)
      new_uv = UVData()
      new_uv.read_uvfits(uvfits_name)
      #print(new_uv.Nfreqs)
      #print tmp_uv.antenna_positions[0:5]
      #print new_uv.antenna_positions[0:5]
      tmp_uv = tmp_uv + new_uv
      #print(tmp_uv.Nfreqs)
   tmp_uv.write_uvfits(outname)
   print("wrote %s" % outname)
   #delete all the temporary uvfits files
   if delete:
      for uvfits_name in uvfits_list:
         cmd = "rm -rf %s" % (uvfits_name)
         print(cmd)
         os.system(cmd)
         #also delete the miriad vis file
         vis_name = uvfits_name.split('.')[0] + '.vis'
         cmd = "rm -rf %s" % (vis_name)
         print(cmd)
         os.system(cmd)
         
         
def concat_uvfits_casa(uvfits_list,outname,delete=False):
   #write a bash file for casa to read
   for uvfits_name in uvfits_list:
      vis_name = uvfits_name.split('.')[0] + '.ms'
      casa_string = "importuvfits(fitsfile='%s',vis='%s')" % (uvfits_name,vis_name)
      casa_filename = 'casa_concat.sh'
      with open(casa_filename,'w') as f:
         f.write(casa_string)
      cmd = "casa -c --nocrashreport %s" % casa_filename
      print(cmd)
      os.system(cmd)

   if delete:
      for uvfits_name in uvfits_list:
         cmd = "rm -rf %s" % (uvfits_name)
         print(cmd)
         os.system(cmd)
         #also delete the miriad vis file
         vis_name = uvfits_name.split('.')[0] + '.vis'
         cmd = "rm -rf %s" % (vis_name)
         print(cmd)
         os.system(cmd)
         #also delete the miriad ms file
         vis_name = uvfits_name.split('.')[0] + '.ms'
         cmd = "rm -rf %s" % (vis_name)
         print(cmd)
         os.system(cmd)

         
  
def concat_vis(vis_list,outname):
   tmp_list = []
   for vis_index,vis in enumerate(vis_list): 
      tmp_name = "t%s.v" % vis_index
      tmp_list.append(tmp_name)
      cmd = "mv  %s %s" % (vis,tmp_name)
      print(cmd)
      os.system(cmd)

   vis_list_string = ','.join(tmp_list)
   cmd = "uvaver vis=%s out=%s line=channel,%s,1,1,1" % (vis_list_string,outname,n_chan)
   print(cmd)
   os.system(cmd)
   #cmd = "rm -rf t*.v" 
   #print(cmd)
   #os.system(cmd)
   
def polar2cartesian(r, t, grid, x, y, order=3):

    X, Y = np.meshgrid(x, y)

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(X, Y)

    ir = interp1d(r, np.arange(len(r)), bounds_error=False)
    it = interp1d(t, np.arange(len(t)))

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    return map_coordinates(grid, np.array([new_ir, new_it]),
                            order=order).reshape(new_r.shape)
                          
#Now for the beam (EDA)!
def generate_new_average_beam():
   if freq_MHz==160.:
      EDA_chan=125
   else:
      print("no EDA beam data for this frequency (%s MHz)" % freq_MHz)
   for freq_MHz in freq_MHz_list:
      for pol in pol_list:   
         Etheta_mag_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Etheta_mag_ch%s.fits" % (EDA_chan,pol,EDA_chan)
         Etheta_phase_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Etheta_phase_ch%s.fits" % (EDA_chan,pol,EDA_chan)
         Ephi_mag_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Ephi_mag_ch%s.fits" % (EDA_chan,pol,EDA_chan)
         Ephi_phase_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Ephi_phase_ch%s.fits" % (EDA_chan,pol,EDA_chan)
   
         Ephi_mag_HDUlist = pyfits.open(Ephi_mag_filename)
         Ephi_mag_header = Ephi_mag_HDUlist[0].header
      
         Ephi_mag_HDUlist.verify('fix')
         Ephi_mag_data = Ephi_mag_HDUlist[0].data
         
         Ephi_phase_HDUlist = pyfits.open(Ephi_phase_filename)
         Ephi_phase_header = Ephi_phase_HDUlist[0].header
         Ephi_phase_HDUlist.verify('fix')
         Ephi_phase_data = Ephi_phase_HDUlist[0].data   
      
         Etheta_mag_HDUlist = pyfits.open(Etheta_mag_filename)
         Etheta_mag_header = Etheta_mag_HDUlist[0].header
         Etheta_mag_HDUlist.verify('fix')
         Etheta_mag_data = Etheta_mag_HDUlist[0].data
       
         Etheta_phase_HDUlist = pyfits.open(Etheta_phase_filename)
         Etheta_phase_header = Etheta_phase_HDUlist[0].header
         Etheta_phase_HDUlist.verify('fix')
         Etheta_phase_data = Etheta_phase_HDUlist[0].data  
         
         #Want an average beam pattern
         average_beam_array = np.zeros([361,91])
         
         for ant_index in range(0,n_ants-1):
            #ant_name_randall = line.split()[0][3:6]
            ant_name_randall = "%03d" % (ant_index+1)
            #print("Ant %s" % ant_name_randall)
           
            Power_pattern_ant_interp_title = 'Antenna %s Power Pattern %s Pol %s MHz' % (ant_name_randall,pol,int(freq_MHz))
            Power_pattern_ant_interp_figname = 'Power_pattern_Ant%s_%s_%s_MHz_interp.png' % (ant_name_randall,pol,int(freq_MHz))
            power_pattern_ant_filename = '/mnt/md0/AAVS-1/beam_pattern_tests/new_data/%s_MHz/AAVS1_power_pattern_linear_Ant%s_%s_theta_phi.npy' % (int(freq_MHz),ant_name_randall,pol)
            
            Ephi_phase_data_ant = Ephi_phase_data[ant_index]
            Ephi_mag_data_ant = Ephi_mag_data[ant_index]
            Etheta_phase_data_ant = Etheta_phase_data[ant_index]
            Etheta_mag_data_ant = Etheta_mag_data[ant_index]
            
            power_pattern_ant = Etheta_mag_data_ant**2 + Ephi_mag_data_ant**2
            
            #print(power_pattern_ant.shape)
            #print(power_pattern_ant.max())
            #need to normalise (can do at end?)
            #include a 1/cos(za) term to account for pixel area changing in sine projection
            #need to check I have the order of ZA right...
            for za_deg in range(0,91):
               za_rad = (za_deg/180.)*np.pi
               #print("zenith angle is %s deg, divide by factor cos(za)=%s" % (za_deg,np.cos(za_rad)))
               power_pattern_ant[:,za_deg] = power_pattern_ant[:,za_deg]/np.cos(za_rad)
            
            average_beam_array += power_pattern_ant
            
            #Ephi_real_ant = Ephi_mag_data_ant * np.cos(Ephi_phase_data_ant)
            #Ephi_imag_ant = Ephi_mag_data_ant * np.sin(Ephi_phase_data_ant)
            #Ephi_complex_ant = Ephi_real_ant+Ephi_imag_ant*1j
            #Etheta_real_ant = Etheta_mag_data_ant * np.cos(Etheta_phase_data_ant)
            #Etheta_imag_ant = Etheta_mag_data_ant * np.sin(Etheta_phase_data_ant)
            #Etheta_complex_ant = Etheta_real_ant+Etheta_imag_ant*1j
         
            #if apply_normalisation:
            #            
            #   Ephi_complex_ant = Ephi_complex_ant * normalisation_factor
            #   Etheta_complex_ant = Etheta_complex_ant * normalisation_factor
            
            
            if plot_all_beams:
               #Cartesian grid:
               power_patter_ant_log = 10*np.log10(power_pattern_ant)
               
               # Define original polar grid
               
               nr = 91
               nt = 361
               r = np.linspace(1, 91, nr)
               t = np.linspace(-np.pi, np.pi, nt)
               z = np.swapaxes(power_patter_ant_log,0,1)
               
               nx = 180
               ny = 180
               x = np.linspace(-90, 90., nx)
               y = np.linspace(-90., 90., ny)
               
               # Interpolate polar grid to cartesian grid (nearest neighbor)
               
               #blank out below horizon:
               centre_x = nx/2
               centre_y = ny/2
               #centre_x = 3
               #centre_y = 3
               y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
               mask = x**2+y**2 <= centre_x**2
               mask = 1*mask.astype(float)
               
               
               fig = plt.figure()
               ax = fig.add_subplot(111)
               cart_image = polar2cartesian(r, t, z, x, y, order=3)
               cart_image = cart_image*mask
               cart_image[cart_image==-0.] = np.nan
               ###get east and west right
               ##cart_image = np.flip(cart_image,axis=0)
               ###sort out N-S (same as rotating 180
               ##cart_image = np.flip(cart_image,axis=1)
               cart_image = np.rot90(cart_image, k=1, axes=(0,1))
               
               #max_ant001_db = 10*np.log10(max_ant001)
               #min_ant001_db = 10*np.log10(min_ant001)
               #
   
               img = ax.imshow(cart_image, interpolation='nearest', vmax=0,vmin=-25)
               #fig.savefig('test4.png')
            
               cbar = plt.colorbar(img, ax=ax)
               cbar.set_label('Power', rotation=270, labelpad=10)
               #plt.colorbar(img, ax=ax)
               plt.title(Power_pattern_ant_interp_title)
               plt.savefig(Power_pattern_ant_interp_figname)
               plt.close()
               print("saved %s" % Power_pattern_ant_interp_figname)
      
         #calculate average
         average_beam_array = average_beam_array/n_ants
         #normalise the average beam
         # max of average beam array:
         #print average_beam_array
         average_array_max = average_beam_array.max()
         #print("%s" % average_array_max)
         average_beam_array_norm = average_beam_array/average_beam_array.max()
      
         #sin projection sampling 
         #Need to do a better job than this simple nearest neighbour interpolation - but fine for now.
         #print("av beam array shape is")
         #print average_beam_array.shape
         sin_projected_beam_array=average_beam_array*0
         step=1./91.
         sin_za_array=np.arange(0,1,step)
         #print sin_za_array
         sin_proj_sampling_rad = np.arcsin(sin_za_array)
         sin_proj_sampling_deg = sin_proj_sampling_rad/np.pi*180.
         sin_proj_sampling_deg_round = np.around(sin_proj_sampling_deg)
         sin_proj_sampling_deg_indexes = sin_proj_sampling_deg_round.astype(int)
         #print sin_proj_sampling_deg_indexes.shape
         
         for az in range(0,361):
            power_vs_za_array_proj = sin_proj_sampling_deg_indexes*0.
            for orig_za_index in range(0,91):
               new_za_index = sin_proj_sampling_deg_indexes[orig_za_index]
               power_vs_za_array_proj[orig_za_index] = average_beam_array[az,new_za_index]
            sin_projected_beam_array[az,:] = power_vs_za_array_proj
            
         sin_projected_beam_array_norm = sin_projected_beam_array/sin_projected_beam_array.max()
            
         #   sin_za_array = np.sin(za_array_rad)
         #   print za_array_rad
         #   print sin_za_array
      
         if plot_average_beam:
            power_pattern_average_interp_title = 'Average Antenna Power Pattern %s Pol %s MHz' % (pol,int(freq_MHz))
            power_pattern_average_interp_figname = 'power_pattern_average_%s_%s_MHz_interp.png' % (pol,int(freq_MHz))
            power_pattern_average_interp_fitsname =  'power_pattern_average_%s_%s_MHz_interp.fits' % (pol,int(freq_MHz))
            
            
            #Cartesian grid:
            #power_pattern_ant_log = 10*np.log10(average_beam_array_norm)
            
            # Define original polar grid
            
            nr = 91
            nt = 361
            r = np.linspace(1, 91, nr)
            t = np.linspace(-np.pi, np.pi, nt)
            z = np.swapaxes(average_beam_array_norm,0,1)
            
            nx = 180
            ny = 180
            x = np.linspace(-90, 90., nx)
            y = np.linspace(-90., 90., ny)
            
            # Interpolate polar grid to cartesian grid (nearest neighbor)
            
            #blank out below horizon:
            centre_x = nx/2
            centre_y = ny/2
            #centre_x = 3
            #centre_y = 3
            y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
            mask = x**2+y**2 <= centre_x**2
            mask = 1*mask.astype(float)
            
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cart_image = polar2cartesian(r, t, z, x, y, order=3)
            cart_image = cart_image*mask
            cart_image[cart_image==-0.] = np.nan
            ###get east and west right
            ##cart_image = np.flip(cart_image,axis=0)
            ###sort out N-S (same as rotating 180
            ##cart_image = np.flip(cart_image,axis=1)
            cart_image = np.rot90(cart_image, k=1, axes=(0,1))
            
            #max_ant001_db = 10*np.log10(max_ant001)
            #min_ant001_db = 10*np.log10(min_ant001)
            #
        
            img = ax.imshow(cart_image, interpolation='nearest')
            #fig.savefig('test4.png')
            
            cbar = plt.colorbar(img, ax=ax)
            cbar.set_label('Power', rotation=270, labelpad=10)
            #plt.colorbar(img, ax=ax)
            plt.title(power_pattern_average_interp_title)
            plt.savefig(power_pattern_average_interp_figname)
            plt.close()
            print("saved %s" % power_pattern_average_interp_figname)
            
            pyfits.writeto(power_pattern_average_interp_fitsname,cart_image,clobber=True)
            print("wrote fits image %s" %  power_pattern_average_interp_fitsname)
   
   
            ##SIN projected!
            power_pattern_average_interp_sin_title = 'Average Antenna Power Pattern SIN Proj %s Pol %s MHz' % (pol,int(freq_MHz))
            power_pattern_average_interp_sin_figname = 'power_pattern_average_%s_%s_MHz_interp_sin.png' % (pol,int(freq_MHz))
            power_pattern_average_interp_sin_fitsname =  'power_pattern_average_%s_%s_MHz_interp_sin.fits' % (pol,int(freq_MHz))
            
            
            #Cartesian grid:
            #power_pattern_ant_log = 10*np.log10(average_beam_array_norm)
            
            # Define original polar grid
            
            nr = 91
            nt = 361
            r = np.linspace(1, 91, nr)
            t = np.linspace(-np.pi, np.pi, nt)
            z = np.swapaxes(sin_projected_beam_array_norm,0,1)
            
            nx = 180
            ny = 180
            x = np.linspace(-90, 90., nx)
            y = np.linspace(-90., 90., ny)
            
            # Interpolate polar grid to cartesian grid (nearest neighbor)
            
            #blank out below horizon:
            centre_x = nx/2
            centre_y = ny/2
            #centre_x = 3
            #centre_y = 3
            y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
            mask = x**2+y**2 <= centre_x**2
            mask = 1*mask.astype(float)
            
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cart_image_sin = polar2cartesian(r, t, z, x, y, order=3)
            cart_image_sin = cart_image_sin*mask
            cart_image_sin[cart_image_sin==-0.] = np.nan
            ###get east and west right
            ##cart_image = np.flip(cart_image,axis=0)
            ###sort out N-S (same as rotating 180
            ##cart_image = np.flip(cart_image,axis=1)
            cart_image_sin = np.rot90(cart_image_sin, k=1, axes=(0,1))
            
            #max_ant001_db = 10*np.log10(max_ant001)
            #min_ant001_db = 10*np.log10(min_ant001)
            #
        
            img = ax.imshow(cart_image_sin, interpolation='nearest')
            #fig.savefig('test4.png')
            
            cbar = plt.colorbar(img, ax=ax)
            cbar.set_label('Power', rotation=270, labelpad=10)
            #plt.colorbar(img, ax=ax)
            plt.title(power_pattern_average_interp_sin_title)
            plt.savefig(power_pattern_average_interp_sin_figname)
            plt.close()
            print("saved %s" % power_pattern_average_interp_sin_figname)
            
        
         ##need to write as a fits file with the correct header
   
         write_beam_fits_sin(cart_image_sin,power_pattern_average_interp_sin_fitsname)
         
         #Okay, I think this looks alright. Have a SIN projected beam and a SIN projected gsm

def generate_smooth_hpx_cube(sky_model,freq_MHz_list,gsm_smooth_poly_order):
   freq_MHz_array = np.asarray(freq_MHz_list)
   log_freq_MHz_array = np.log10(freq_MHz_array)
   if sky_model == 'gsm2016':
      gsm = GlobalSkyModel2016() 
   elif  sky_model == 'gsm':
      gsm = GlobalSkyModel()
   else:
      print('unknown sky_model')
      sys.exit()
        
   gsm_hpx_cube_smooth_filename = "%s_cube_hpx.npy" % (sky_model)
    
   cmd = "rm -rf %s" % (gsm_hpx_cube_smooth_filename)
   print(cmd)
   os.system(cmd)

   template_gsm = gsm.generate(100)
   npix = template_gsm.shape[0]
   n_freq = freq_MHz_array.shape[0]
   gsm_cube = np.full((npix,n_freq), np.nan)
   gsm_cube_smooth = np.full((npix,n_freq), np.nan)
   
   for freq_MHz_index,freq_MHz in enumerate(freq_MHz_array):
      gsm_map = gsm.generate(freq_MHz)
      gsm_cube[:,freq_MHz_index] = gsm_map
    
   #for each pixel in the cube, fit a polynomial in log space across frequency
   for pix_index,pix in enumerate(gsm_cube[:,0][0:1]):
      spectrum = gsm_cube[pix_index,:]
      log_spectrum = np.log10(spectrum)
      coefs = poly.polyfit(log_freq_MHz_array, log_spectrum, gsm_smooth_poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      gsm_cube_smooth[pix_index,:] = ffit_linear
      
      plt.clf()
      plt.plot(freq_MHz_array,ffit_linear,label='model fit')
      plt.plot(freq_MHz_array,spectrum,label='spectrum')
      map_title = "%s fit pix %s" % (sky_model,pix)
      plt.ylabel("Tb (K)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name = "%s_cube_hpx_fit_pix_%s.png" % (sky_model,pix_index)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      
      residual = spectrum - ffit_linear
      
      plt.clf()
      plt.plot(freq_MHz_array,residual,label='residual')
      map_title = "%s residual %s" % (sky_model,pix)
      plt.ylabel("residual Tb (K)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name = "%s_cube_hpx_fit_residual_pix_%s.png" % (sky_model,pix_index)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      
      
   np.save(gsm_hpx_cube_smooth_filename,gsm_cube_smooth)
   print("saved %s" % (gsm_hpx_cube_smooth_filename))

#generate_smooth_hpx_cube(sky_model,freq_MHz_list,gsm_smooth_poly_order)

def generate_apparent_sky_model(pol,lst_hrs,freq_MHz):
   print("lst %s hrs" % lst)
   lst_deg = (float(lst)/24.)*360.
   apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)

   print(freq_MHz)
   #and datetime of observation (eventually do this for many dates)
   year=2000
   month=1
   day=1
   hour=np.floor(float(lst))
   minute=np.floor((float(lst)-hour) * 60.)
   second=((float(lst)-hour) * 60. - minute) * 60.
   
   date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
   #print date_time_string
   #for miriad time just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
   day_frac_plus_one = float(lst)/24. + 1
   miriad_uvgen_time_string = '00JAN%1.3f' % day_frac_plus_one
   #print miriad_uvgen_time_string
   wavelength = 300./freq_MHz
   freq_GHz = freq_MHz/1000.
   
   gsm_hpx_fits_name = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
   reprojected_gsm_fitsname = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
   reprojected_gsm_im_name = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)
   
   base_vis_name = "%s_LST_%03d_%s_MHz.vis" % (array_label,lst_deg,int(freq_MHz))
   #eda_model_uvfits_name = "eda_model_LST_%03d_%s_MHz.uvfits" % (lst_deg,int(freq_MHz))
   eda_model_no_source_image_name = "%s_no_src_LST_%03d_%s_MHz.map" % (array_label,lst_deg,int(freq_MHz))
   eda_model_no_source_beam_name = "%s_no_src_LST_%03d_%s_MHz.beam" % (array_label,lst_deg,int(freq_MHz))
   eda_model_no_source_image_fits_name = "%s_no_src_LST_%03d_%s_MHz.fits" % (array_label,lst_deg,int(freq_MHz))
 


   #Need a pb-attenuated gsm image in slant-orthographic projection
   # - obtain gsm image at desired frequency and datetime (healpix)
   if sky_model == 'gsm2016':
      #forget about all this observer stuff, just generate the all-sky map and let reproject take care of the rest 
      #ov = GSMObserver2016()
      gsm = GlobalSkyModel2016() 
   elif  sky_model == 'gsm':
      #ov = GSMObserver()
      gsm = GlobalSkyModel()
   elif sky_model == 'gmoss':
      gmoss_data_text_file_name = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/GMOSS_sky_spectra.txt"
      #from Mayuri: The data is in HEALPIX Nested Scheme with NSIDE 16 and 5 degrees resolution. 
      #It contains GMOSS spectra generated at the 3072 pixels. The first row gives the frequency [GHz] 
      #going from 0.04 GHz to 0.2 GHz in steps of 0.001 GHz. Sky spectra are in Kelvin units
      sky_model_lines_list = []
      with open(gmoss_data_text_file_name, 'r') as f:
         sky_model_lines = f.readlines()
      for line in sky_model_lines:
         sky_model_lines_list.append(line.split())
      sky_model_data_array = np.full((len(sky_model_lines_list),len(sky_model_lines_list[0])),np.nan)
      for line_values_index,line_vales in enumerate(sky_model_lines_list):
         sky_model_data_array[line_values_index,:] = [float(i.strip()) for i in sky_model_lines_list[line_values_index]]
      gmoss_freqs_GHz = sky_model_data_array[0]
      #gmoss_freqs_GHz = np.asarray(gmoss_freqs_GHz)
      gmoss_freqs_MHz = gmoss_freqs_GHz*1000.
   else:
      print('unknown sky_model')
      sys.exit()
   #ov.lon = longitude_degrees
   #ov.lat = latitude_degrees
   #ov.elev = elevation_m
   #
   #ov.date = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
   
   
   #1. get 'blank' visibilities in miriad and make an image (1 Jy pt source) to use as a template for the output projection
   if generate_new_vis:
   
      cmd = "rm -rf %s" % base_vis_name
      print(cmd)
      os.system(cmd)
     
      #need to change this so that each uvgen call has phase centre at zenith (RA = LST)
      cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (array_ant_locations_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,base_vis_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "prthd in=tmp.vis "
      print(cmd)
      os.system(cmd)
            
      cmd = "uvlist vis=tmp.vis options=array,full"
      print(cmd)
      os.system(cmd)

      
      
      #output the no source uvfits file for checking
      cmd = "fits in=%s out=%s op=uvout" % (eda_model_vis_name,eda_model_uvfits_name)
      print(cmd)
      os.system(cmd)
      
      sys.exit()
      
      cmd = "rm -rf %s %s %s " % (eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "invert vis=%s map=%s beam=%s imsize=%s cell=%s stokes=xx options=mfs" % (base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,template_imsize,template_cell_size_asec)
      print(cmd)
      os.system(cmd)
   
      cmd = "fits in=%s out=%s op=xyout" % (eda_model_no_source_image_name,eda_model_no_source_image_fits_name)
      print(cmd)
      os.system(cmd)
      
      
   if generate_new_hpx:
      cmd = "rm -rf %s %s" % (gsm_hpx_fits_name,reprojected_gsm_im_name)
      print(cmd)
      os.system(cmd)
      if sky_model == 'gmoss':
         gmoss_freq_index, gmoss_freq_value = find_nearest(gmoss_freqs_MHz,freq_MHz)
         #print "gmoss_freq_index, gmoss_freq_value %s %s" % (gmoss_freq_index, gmoss_freq_value)
         sky_model_data_nested = sky_model_data_array[:,gmoss_freq_index][1:]
         hp.write_map(gsm_hpx_fits_name,sky_model_data_nested,coord='G',nest=True)
         gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
      else:
         gsm_map = gsm.generate(freq_MHz)
         hp.write_map(gsm_hpx_fits_name,gsm_map,coord='G')
   else:
      if sky_model == 'gmoss':
         sky_model_data_nested = hp.read_map(gsm_hpx_fits_name,nest=True)
         gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
      else:
         gsm_map = hp.read_map(gsm_hpx_fits_name)
   
   #plot?
   if plot_gsm_map_hpx:
      #plot
      plt.clf()
      map_title="GSM from MWA at %.0f:%.0f:%.0f and %0.3f MHz" % (hour,minute,second,freq_MHz)
      #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
      hp.orthview(map=gsm_map,half_sky=False,title=map_title)
      #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
      fig_name="%s_map_%s_%0.3f_MHz.png" % (sky_model,date_time_string,freq_MHz)
      figmap = plt.gcf()
      figmap.savefig(fig_name,dpi=500)
      print("saved %s" % fig_name)
   
   #Miriad doesn't seem to be able to import the hpx file
   #Try using reproject_from_healpix
   #output_projection is the header of the pt source made above
   if os.path.isfile(eda_model_no_source_image_fits_name) and os.access(eda_model_no_source_image_fits_name, os.R_OK):
      hdulist = pyfits.open(eda_model_no_source_image_fits_name)
   else:
      print("Either file %s is missing or is not readable" % eda_model_no_source_image_fits_name)
      #continue        
   
   data=hdulist[0].data[0,0,:,:]
   no_source_header=hdulist[0].header
   pix_size_deg = float(no_source_header['CDELT1'])
   pix_area_deg_sq = pix_size_deg*pix_size_deg
   pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
   
   #needs cooordsys keyword
   #pt_source_header['COORDSYS'] = 'icrs'
   
   #print(pt_source_header)
   
   #del pt_source_header[8]
   #del pt_source_header[8]
   del no_source_header['history']
                   
   #print(pt_source_header)
   
   target_wcs = WCS(no_source_header)
   
   target_wcs=target_wcs.dropaxis(2)
   target_wcs=target_wcs.dropaxis(2)
                   
   #hdu_hpx = pyfits.open(gsm_hpx_fits_name)[1]
   ##hdu_hpx.info()
   #hpx_header = hdu_hpx.header
   #print(hpx_header)
   
   hdu_gsm = fits.open(gsm_hpx_fits_name)[1]
   #hdu_gsm.info()
   #data_gsm = hdu_gsm.data
   
   reprojected_gsm_map,footprint = reproject_from_healpix(hdu_gsm, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
   #reprojected_gsm_map,footprint = reproject_from_healpix((data_gsm,'galactic'), target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',nested=False)
   
   #calc sky-average temp, no beam
   sky_average_temp_no_beam = np.nanmean(reprojected_gsm_map)
   #print sky_average_temp_no_beam
   sky_averaged_diffuse_array_no_beam_lsts[freq_MHz_index] += sky_average_temp_no_beam
   
   #write the map to fits
   pyfits.writeto(reprojected_gsm_fitsname,reprojected_gsm_map,clobber=True)
   pyfits.update(reprojected_gsm_fitsname,reprojected_gsm_map,header=no_source_header)
   print("wrote image %s" %  reprojected_gsm_fitsname)
   
   #Do this GSM map stuff here as it doesn't depend on pol
   cmd = "rm -rf %s" % (reprojected_gsm_im_name)
   print(cmd)
   os.system(cmd)     
   
   cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_fitsname,reprojected_gsm_im_name)
   print(cmd)
   os.system(cmd)
   
   #uvmodel requires the model to be in Jy/pix
   #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
   
   scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
   print("scale map by %s to get to Jy/pix" % scale)
   
   reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
   
   cmd = "rm -rf %s" % reprojected_gsm_im_Jy_per_pix_name
   print(cmd)
   os.system(cmd)
         
   cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name)
   print(cmd)
   os.system(cmd)

   model_vis_name_base = "%s_LST_%03d_%s_%s_MHz" % (array_label,lst_deg,pol,int(freq_MHz))
   
   if use_analytic_beam:
      if pol=='X':
         beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
      else:
         beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
   else:
         beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
   
   #power_pattern_average_interp_sin_im_name = 'power_pattern_average_%s_%s_MHz_interp_sin.im' % (pol,int(freq_MHz))
   #power_pattern_average_interp_sin_regrid_gsm_im_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.im' % (pol,int(freq_MHz))
   #power_pattern_average_interp_sin_regrid_gsm_fits_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.fits' % (pol,int(freq_MHz))
        
   beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
   beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
   beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz)
   
   cmd = "rm -rf %s %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name, beam_image_sin_projected_puthd_fits_name)
   print(cmd)
   os.system(cmd)
   
   cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
   print(cmd)
   os.system(cmd)
   
   
   #put in the correct ra in the header (ra = lst for zenith) 
   #puthd in="$beam/crval1" value=$lst_degs
   cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
   print(cmd)
   os.system(cmd)         
      
   #write out as a fits file to check header
   cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
   print(cmd)
   os.system(cmd)
   
   #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
   cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name)
   print(cmd)
   os.system(cmd)   
   
   #Sweet! finally have a gsm and a beam.
   #now multiply them to get the apparent sky and put that into uvmodel, 
   
   apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)

   cmd = "rm -rf %s" % (apparent_sky_im_name)
   print(cmd)
   os.system(cmd)
   
   cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
   print(cmd)
   os.system(cmd)

           
#change simulate around a bit for ASSASSIN to speed it up (assume same sky for the multiple arrays)
def simulate_assassin(lst_list,freq_MHz_list,pol_list,signal_type_list,sky_model,outbase_name,n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length_m):
   if generate_new_assassin_sky:
      ##2. Sky model
      #Need a pb-attenuated gsm image in slant-orthographic projection
      # - obtain gsm image at desired frequency and datetime (healpix)
      if sky_model == 'gsm2016':
         #forget about all this observer stuff, just generate the all-sky map and let reproject take care of the rest 
         #ov = GSMObserver2016()
         gsm = GlobalSkyModel2016() 
      elif  sky_model == 'gsm':
         #ov = GSMObserver()
         gsm = GlobalSkyModel()
      elif sky_model == 'gmoss':
         gmoss_data_text_file_name = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/GMOSS_sky_spectra.txt"
         #from Mayuri: The data is in HEALPIX Nested Scheme with NSIDE 16 and 5 degrees resolution. 
         #It contains GMOSS spectra generated at the 3072 pixels. The first row gives the frequency [GHz] 
         #going from 0.04 GHz to 0.2 GHz in steps of 0.001 GHz. Sky spectra are in Kelvin units
         sky_model_lines_list = []
         with open(gmoss_data_text_file_name, 'r') as f:
            sky_model_lines = f.readlines()
         for line in sky_model_lines:
            sky_model_lines_list.append(line.split())
         sky_model_data_array = np.full((len(sky_model_lines_list),len(sky_model_lines_list[0])),np.nan)
         for line_values_index,line_vales in enumerate(sky_model_lines_list):
            sky_model_data_array[line_values_index,:] = [float(i.strip()) for i in sky_model_lines_list[line_values_index]]
         gmoss_freqs_GHz = sky_model_data_array[0]
         #gmoss_freqs_GHz = np.asarray(gmoss_freqs_GHz)
         gmoss_freqs_MHz = gmoss_freqs_GHz*1000.
      else:
         print('unknown sky_model')
         sys.exit()
               
      for lst in lst_list:
         print("lst %s hrs" % lst)
         lst_deg = (float(lst)/24.)*360.
         print("lst %s deg" % lst_deg)
         
         year=2000
         month=1
         day=1
         hour=np.floor(float(lst))
         minute=np.floor((float(lst)-hour) * 60.)
         second=((float(lst)-hour) * 60. - minute) * 60.
         
         date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
         print(date_time_string)
         #for miriad time just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
         day_frac_plus_one = float(lst)/24. + 1
         miriad_uvgen_time_string = '00JAN%1.3f' % day_frac_plus_one
         print(miriad_uvgen_time_string)
         
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            print(freq_MHz)
            wavelength = 300./freq_MHz
            freq_GHz = freq_MHz/1000.
            
            gsm_hpx_fits_name = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_fitsname = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_im_name = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)
   
            global_signal_hpx_fits_name = "global_signal_map_LST_%03d_%0.3f_MHz_hpx.fits" % (lst_deg,freq_MHz)
            reprojected_global_signal_fitsname = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (lst_deg,freq_MHz)
            reprojected_global_signal_im_name = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.im" % (lst_deg,freq_MHz)
      
            if generate_new_hpx:
               cmd = "rm -rf %s %s" % (gsm_hpx_fits_name,reprojected_gsm_im_name)
               print(cmd)
               os.system(cmd)
               if sky_model == 'gmoss':
                  gmoss_freq_index, gmoss_freq_value = find_nearest(gmoss_freqs_MHz,freq_MHz)
                  #print "gmoss_freq_index, gmoss_freq_value %s %s" % (gmoss_freq_index, gmoss_freq_value)
                  sky_model_data_nested = sky_model_data_array[:,gmoss_freq_index][1:]
                  hp.write_map(gsm_hpx_fits_name,sky_model_data_nested,coord='G',nest=True)
                  gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
               else:
                  gsm_map = gsm.generate(freq_MHz)
                  hp.write_map(gsm_hpx_fits_name,gsm_map,coord='G')
            else:
               if sky_model == 'gmoss':
                  sky_model_data_nested = hp.read_map(gsm_hpx_fits_name,nest=True)
                  gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
               else:
                  gsm_map = hp.read_map(gsm_hpx_fits_name)
      
            #plot?
            if plot_gsm_map_hpx:
               #plot
               plt.clf()
               map_title="GSM from MWA at %.0f:%.0f:%.0f and %0.3f MHz" % (hour,minute,second,freq_MHz)
               #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
               hp.orthview(map=gsm_map,half_sky=False,title=map_title)
               #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
               fig_name="%s_map_%s_%0.3fMHz.png" % (sky_model,date_time_string,freq_MHz)
               figmap = plt.gcf()
               figmap.savefig(fig_name,dpi=500)
               print("saved %s" % fig_name)
      
            ##########
            #assassin-specific bit:
            #angles are anticlockwise from East
            
            #need an image to use to make the reprojected sky maps, can use the same image (as we are assuming the same LST for all arrays)
            #so get your array of angles and diameters first, then take the first one to use for the initial image
            radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
            circle_array = range(0,n_circles)
            
            circle_number = circle_array[0]
   
            #circle 1 is the smallest
            radius = (min_arm_length_m + circle_number * radial_spacing) 
            diameter = radius * 2.
            diameter_cm = int(diameter*100.)
            #work out circumference 
            circumference = math.pi * diameter
            #print diameter
            n_angles = int(round(circumference * n_ants_per_m_of_circumference))
            angle_increment = (2.*math.pi)/n_angles
            #(remember only need half of them!)
            angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment    
            angle_rad = angle_array_rad[0]
            angle_deg = int(angle_rad/np.pi*180.)
   
            antenna_position_filename = "assassin_baseline_%03d_deg_%03d_cm.txt" % (angle_deg,diameter_cm)
            ant_1_x_offset = radius * np.cos(angle_rad)
            ant_1_y_offset = radius * np.sin(angle_rad)
            ant_1_position_string = "%0.3f   %0.3f   0\n" % (ant_1_x_offset,ant_1_y_offset)
            ant_2_angle = angle_rad + np.pi
            ant_2_x_offset = radius * np.cos(ant_2_angle)
            ant_2_y_offset = radius * np.sin(ant_2_angle)
            ant_2_position_string = "%0.3f   %0.3f   0" % (ant_2_x_offset,ant_2_y_offset)         
            with open(antenna_position_filename,'w') as f:
               f.write(ant_1_position_string)
               f.write(ant_2_position_string)
            print("wrote %s" % antenna_position_filename)
            
            sub_array_string = "sub_%03d_%03d" % (diameter_cm,angle_deg)
               
            base_vis_name = "%s_LST_%03d_%s_MHz.vis" % (sub_array_string,lst_deg,int(freq_MHz))
            eda_model_no_source_image_name = "no_src_LST_%03d_%s_MHz.map" % (lst_deg,int(freq_MHz))
            eda_model_no_source_beam_name = "no_src_LST_%03d_%s_MHz.beam" % (lst_deg,int(freq_MHz))
            eda_model_no_source_image_fits_name = "no_src_LST_%03d_%s_MHz.fits" % (lst_deg,int(freq_MHz))
            
            #put all the apparent sky stuff here and just do it once (for each freq)
            cmd = "rm -rf %s" % base_vis_name
            print(cmd)
            os.system(cmd)
            
            #need to change this so that each uvgen call has phase centre at zenith (RA = LST)
            cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (antenna_position_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,base_vis_name)
            print(cmd)
            os.system(cmd)
            
            #output the no source uvfits file for checking
            #cmd = "fits in=%s out=%s op=uvout" % (eda_model_vis_name,eda_model_uvfits_name)
            #print(cmd)
            #os.system(cmd)
            
            cmd = "rm -rf %s %s %s " % (eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "invert vis=%s map=%s beam=%s imsize=%s cell=%s stokes=xx options=mfs" % (base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,template_imsize,template_cell_size_asec)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyout" % (eda_model_no_source_image_name,eda_model_no_source_image_fits_name)
            print(cmd)
            os.system(cmd)
   
            #Miriad doesn't seem to be able to import the hpx file
            #Try using reproject_from_healpix
            #output_projection is the header of the pt source made above
            if os.path.isfile(eda_model_no_source_image_fits_name) and os.access(eda_model_no_source_image_fits_name, os.R_OK):
               hdulist = pyfits.open(eda_model_no_source_image_fits_name)
            else:
               print("Either file %s is missing or is not readable" % eda_model_no_source_image_fits_name)
               #continue        
            
            data=hdulist[0].data[0,0,:,:]
            no_source_header=hdulist[0].header
            pix_size_deg = float(no_source_header['CDELT1'])
            pix_area_deg_sq = pix_size_deg*pix_size_deg
            pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
                   
            #print(pt_source_header)
            
            #del pt_source_header[8]
            #del pt_source_header[8]
            del no_source_header['history']
                            
            #print(pt_source_header)
            
            target_wcs = WCS(no_source_header)
            
            target_wcs=target_wcs.dropaxis(2)
            target_wcs=target_wcs.dropaxis(2)
            
            hdu_gsm = fits.open(gsm_hpx_fits_name)[1]
            #hdu_gsm.info()
            #data_gsm = hdu_gsm.data
            
            reprojected_gsm_map,footprint = reproject_from_healpix(hdu_gsm, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
           
            
            #write the map to fits
            pyfits.writeto(reprojected_gsm_fitsname,reprojected_gsm_map,clobber=True)
            pyfits.update(reprojected_gsm_fitsname,reprojected_gsm_map,header=no_source_header)
            print("wrote image %s" %  reprojected_gsm_fitsname)
            
            #Do this GSM map stuff here as it doesn't depend on pol
            cmd = "rm -rf %s" % (reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)     
            
            cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_fitsname,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)
            
            #uvmodel requires the model to be in Jy/pix
            #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
            
            scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
            print("scale map by %s to get to Jy/pix" % scale)
            
            reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
            
            cmd = "rm -rf %s" % reprojected_gsm_im_Jy_per_pix_name
            print(cmd)
            os.system(cmd)
                  
            cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name)
            print(cmd)
            os.system(cmd)
            
            ########
            #Repeat the above for the global signal
            #make an all sky global signal map here too:
            s_21_value = s_21_array[freq_MHz_index]
            s_21_hpx_map = gsm_map * 0.0 + s_21_value
            
            cmd = "rm -rf %s %s" % (global_signal_hpx_fits_name,reprojected_global_signal_im_name)
            print(cmd)
            os.system(cmd)
               
            hp.write_map(global_signal_hpx_fits_name,s_21_hpx_map,coord='C')
            
            #plot?
            if plot_global_signal_map_hpx:
               #plot
               plt.clf()
               map_title="Global signal from MWA at %.0f:%.0f:%.0f and %0.3f MHz" % (hour,minute,second,freq_MHz)
               #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
               hp.orthview(map=s_21_hpx_map,half_sky=False,title=map_title)
               #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
               fig_name="global_signal_map_%s_%0.3f_MHz.png" % (date_time_string,freq_MHz)
               figmap = plt.gcf()
               figmap.savefig(fig_name)
               print("saved %s" % fig_name)
            
            reprojected_global_signal_map,footprint = reproject_from_healpix(global_signal_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
            
            #write the map to fits
            pyfits.writeto(reprojected_global_signal_fitsname,reprojected_global_signal_map,clobber=True)
            pyfits.update(reprojected_global_signal_fitsname,reprojected_global_signal_map,header=no_source_header)
            print("wrote image %s" %  reprojected_global_signal_fitsname)
            
            cmd = "rm -rf %s" % (reprojected_global_signal_im_name)
            print(cmd)
            os.system(cmd)     
            
            cmd = "fits in=%s out=%s op=xyin" % (reprojected_global_signal_fitsname,reprojected_global_signal_im_name)
            print(cmd)
            os.system(cmd)
            
            #uvmodel requires the model to be in Jy/pix
            #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
            
            reprojected_global_signal_im_Jy_per_pix_name =  "global_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (date_time_string,freq_MHz)
            
            cmd = "rm -rf %s" % reprojected_global_signal_im_Jy_per_pix_name
            print(cmd)
            os.system(cmd)
                  
            cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_global_signal_im_name,reprojected_global_signal_im_Jy_per_pix_name)
            print(cmd)
            os.system(cmd)
             
            #do for all pols:
            for pol in pol_list:      
               if use_analytic_beam:
                  if pol=='X':
                     beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                  else:
                     beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
               else:
                     beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
                    
               beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
               beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
               beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz)
               
               cmd = "rm -rf %s %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name, beam_image_sin_projected_puthd_fits_name)
               print(cmd)
               os.system(cmd)
               
               cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
               print(cmd)
               os.system(cmd)
               
               
               #put in the correct ra in the header (ra = lst for zenith) 
               #puthd in="$beam/crval1" value=$lst_degs
               cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
               print(cmd)
               os.system(cmd)         
               
               #write out as a fits file to check header
               cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
               print(cmd)
               os.system(cmd)
               
               #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
               cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name)
               print(cmd)
               os.system(cmd)   
            
            
               #Sweet! finally have a gsm and a beam.
               #now multiply them to get the apparent sky and put that into uvmodel, 
            
               apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
               
               cmd = "rm -rf %s " % (apparent_sky_im_name)
               print(cmd)
               os.system(cmd)
            
               cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
               print(cmd)
               os.system(cmd)
               
               #Repeat for global signal
               apparent_global_signal_im_name = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
               apparent_global_signal_fits_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.fits" % (lst_deg,pol,freq_MHz)
               apparent_global_signal_im_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz_cropped.im" % (lst_deg,pol,freq_MHz)
            
               cmd = "rm -rf %s %s %s" % (apparent_global_signal_im_name,apparent_global_signal_fits_name_cropped,apparent_global_signal_im_name_cropped)
               print(cmd)
               os.system(cmd)
            
               cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_global_signal_im_Jy_per_pix_name,apparent_global_signal_im_name)
               print(cmd)
               os.system(cmd)
                     
   ######           
   #then for each sub-array (and each lst,freq,pol)    
   radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
   for circle_number in range(0,n_circles):
      #circle 1 is the smallest
      radius = (min_arm_length_m + circle_number * radial_spacing) 
      diameter = radius * 2.
      diameter_cm = int(diameter*100.)
      #work out circumference 
      circumference = math.pi * diameter
      #print diameter
      n_angles = int(round(circumference * n_ants_per_m_of_circumference))
      angle_increment = (2.*math.pi)/n_angles
      #(remember only need half of them!)
      angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
      #angle_array_deg = angle_array_rad / math.pi * 180.
      #print angle_array_deg
      for angle_rad in angle_array_rad:
         angle_deg = int(angle_rad/np.pi*180.)
         sub_array_string = "sub_%03d_%03d" % (diameter_cm,angle_deg)
         antenna_position_filename = "assassin_baseline_%03d_deg_%03d_cm.txt" % (angle_deg,diameter_cm)
         ant_1_x_offset = radius * np.cos(angle_rad)
         ant_1_y_offset = radius * np.sin(angle_rad)
         ant_1_position_string = "%0.3f   %0.3f   0\n" % (ant_1_x_offset,ant_1_y_offset)
         ant_2_angle = angle_rad + np.pi
         ant_2_x_offset = radius * np.cos(ant_2_angle)
         ant_2_y_offset = radius * np.sin(ant_2_angle)
         ant_2_position_string = "%0.3f   %0.3f   0" % (ant_2_x_offset,ant_2_y_offset)         
         with open(antenna_position_filename,'w') as f:
            f.write(ant_1_position_string)
            f.write(ant_2_position_string)
         print("wrote %s" % antenna_position_filename)
         
         for pol in pol_list:
            model_vis_uvfits_list_lsts = []
            for lst in lst_list:
               model_vis_uvfits_list = []
               
               print("lst %s hrs" % lst)
               lst_deg = (float(lst)/24.)*360.
               print("lst %s deg" % lst_deg)
               
               year=2000
               month=1
               day=1
               hour=np.floor(float(lst))
               minute=np.floor((float(lst)-hour) * 60.)
               second=((float(lst)-hour) * 60. - minute) * 60.
               
               date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
               #print date_time_string
               #for miriad time just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
               day_frac_plus_one = float(lst)/24. + 1
               miriad_uvgen_time_string = '00JAN%1.3f' % day_frac_plus_one
               #print miriad_uvgen_time_string 
            
               for freq_MHz in freq_MHz_list:
                  print(freq_MHz)
                  wavelength = 300./freq_MHz
                  freq_GHz = freq_MHz/1000.
                  
                  #for gain errors:
                  one_jy_source_vis_name =  "%s_one_jy_LST_%03d_%s_MHz.vis" % (sub_array_string,lst_deg,int(freq_MHz))
                  gain_solutions_name_amp = "%s_gains_%03d_%s_MHz_amp.txt" % (sub_array_string,lst_deg,int(freq_MHz))
                  gain_solutions_name_phase = "%s_gains_%03d_%s_MHz_phase.txt" % (sub_array_string,lst_deg,int(freq_MHz))
                  ###########
               
                                
                  n_chan = len(freq_MHz_list)
                  n_lsts = len(lst_list)
               
                  sky_averaged_diffuse_array_beam_lsts = np.full((n_chan),0.0)
                  sky_averaged_diffuse_array_no_beam_lsts = np.full((n_chan),0.0)
   
                  model_vis_name_base = "%s_LST_%03d_%s_%s_MHz" % (sub_array_string,lst_deg,pol,int(freq_MHz))

                  apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)

                  apparent_global_signal_fits_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.fits" % (lst_deg,pol,freq_MHz)
                  apparent_global_signal_im_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz_cropped.im" % (lst_deg,pol,freq_MHz)
         
                  base_vis_name = "%s_LST_%03d_%s_%s_MHz.vis" % (sub_array_string,lst_deg,pol,int(freq_MHz))
                  
                  cmd = "rm -rf %s" % base_vis_name
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (antenna_position_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,base_vis_name)
                  print(cmd)
                  os.system(cmd)
                  
                  
                  #then put into the vis 
                        
                  if 'noise' in signal_type_list:
                     model_vis_name_base += '_N'
                     out_vis_name = model_vis_name_base + '.vis'
                     
                     ##Put in the uvgen for the noise here (so we get different noise for X and Y pol (i.e. remove noise uvgen from above!!))
                     #generate noise only visibilites
                     #for systemp, just use sky temp 180 at 180, for JyperK, use Aeff from Table 2 of Wayth et al 2017 (EDA paper)
                     systemp = T_180*(freq_MHz/180.0)**beta
                     print('systemp %s K' % systemp)
                     A_eff = Aeff_for_freq_MHz_list[freq_MHz_index]
                     print('A_eff %s K' % A_eff)
                     JperK = (2.0*k*10**26)/(eta*A_eff)
                     print('JperK %s ' % JperK)
                     SEFD = systemp * JperK
                     print('SEFD %s ' % SEFD)
                  
                     #generate the noise-only uvfits file
                     cmd = "rm -rf %s" % out_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 systemp=%s jyperk=%s out=%s stokes=xx  " % (antenna_position_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec, harange_string, systemp, JperK, out_vis_name)
                     print(cmd)
                     os.system(cmd)
                     
                     base_vis_name = out_vis_name
                     
                  if 'diffuse' in signal_type_list:
                  
                     model_vis_name_base += '_D_%s' % sky_model
                     out_vis_name = model_vis_name_base + '.vis'
                     
                     cmd = "rm -rf %s" % out_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     
                      
                     cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_sky_im_name,out_vis_name)
                     print(cmd)
                     os.system(cmd)   
                     
                     
                     #don't need the previous basevis any more
                     cmd = "rm -rf %s" % base_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     base_vis_name = out_vis_name
                     
                  if 'global' in signal_type_list:
                  
                     model_vis_name_base += '_G'
                     out_vis_name = model_vis_name_base + '.vis'
                  
                     cmd = "rm -rf %s" % out_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_global_signal_im_name,out_vis_name)
                     print(cmd)
                     os.system(cmd)
                     
                     #remove the previous base_vis
                     cmd = "rm -rf %s" % base_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     base_vis_name = out_vis_name
           
                  if 'gain_errors' in signal_type_list:
                     cmd = "rm -rf %s" % one_jy_source_vis_name
                     print(cmd)
                     os.system(cmd)
            
                     # randall's recipe for gain errors
                     #rm -rf "$gainvis" "$noisevis"
                     ## make 1Jy source with antenna-based amp/phase errors, no noise
                     #$MIRCAT/point.source
                     #uvgen ant=antfile.txt baseunit=-3.33564 corr=$corr freq=$freq radec=$lst,$dec harange=0,$ha_max,0.0005556 lat=-26.7 source=/tmp/source1Jy.txt out="$gainvis" time=00JAN1$fakedayfrac gnoise=10 pnoise=90 | grep noise
                     cmd = "uvgen source=$MIRCAT/point.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx gnoise=10 pnoise=90  " % (antenna_position_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,one_jy_source_vis_name)
                     print(cmd)
                     os.system(cmd)
                     ## selfcal this to solve for the gain errors
                     #selfcal vis="$gainvis" flux=1 options=amplitude,noscale
                     cmd = "selfcal vis=%s flux=1 options=amplitude,noscale" % (one_jy_source_vis_name)
                     print(cmd)
                     os.system(cmd)
                     ## write out the solutions
                     #gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_amp.txt
                     #gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_pha.txt yaxis=phase
                     
                     cmd = "rm -rf %s %s" % (gain_solutions_name_amp,gain_solutions_name_phase)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "gpplt vis=%s log=%s" % (one_jy_source_vis_name,gain_solutions_name_amp)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "gpplt vis=%s log=%s yaxis=phase" % (one_jy_source_vis_name,gain_solutions_name_phase)
                     print(cmd)
                     os.system(cmd)
                     #collect these up and plot them later on
                     
                     
                     model_vis_name_base += '_GE'
                     out_vis_name = model_vis_name_base + '.vis'
                     
                     cmd = "rm -rf %s" % out_vis_name
                     print(cmd)
                     os.system(cmd)
                              
                     cmd = "gpcopy vis=%s out=%s options=relax" % (one_jy_source_vis_name,base_vis_name)
                     print(cmd)
                     os.system(cmd)
                     #use uvcat so that the solutions are applied and become the 'data' (not corrected data)
                     #uvcat vis=/tmp/tempvis.uv out="$outname"
                     cmd = "uvcat vis=%s out=%s" % (base_vis_name,out_vis_name)
                     print(cmd)
                     os.system(cmd)
                     
                     #remove the previous base_vis
                     cmd = "rm -rf %s" % base_vis_name
                     print(cmd)
                     os.system(cmd)
                     
                     base_vis_name = out_vis_name
                        
                  out_vis_uvfits_name = model_vis_name_base + '.uvfits' 
   
                  
                  cmd = "rm -rf %s" % (out_vis_uvfits_name)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "fits in=%s out=%s op=uvout options=nocal,nopol,nopass" % (out_vis_name,out_vis_uvfits_name)
                  print(cmd)
                  os.system(cmd)
                  
                  sky_averaged_diffuse_array_no_beam_lsts_filename =  "%s_sky_averaged_diffuse_no_beam.npy" % model_vis_name_base
                  sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % model_vis_name_base
                  
                  output_concat_vis_pyuvdata_name_lsts = "%s_concat_lsts.vis" % model_vis_name_base
                  output_concat_uvfits_pyuvdata_name_lsts = "%s_concat_lsts.uvfits" % model_vis_name_base 
                  
                  
                  #add visname to concat list
                  model_vis_uvfits_list.append(out_vis_uvfits_name)
                  
                  #delete all the intermediate images and vis that are no longer required
                  #sys.exit()
                  if do_cleanup_images_and_vis:
                     cleanup_images_and_vis_assassin(sub_array_string,lst,freq_MHz,pol)
      
      
               #print model_vis_uvfits_list
               
               #cleanup_images_and_vis(sub_array_string,lst,freq_MHz,pol)      
         
               output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
               output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
               
               cmd = "rm -rf %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name)
               print(cmd)
               os.system(cmd)
         
               concat_uvfits(model_vis_uvfits_list,output_concat_uvfits_pyuvdata_name)

               model_vis_uvfits_list_lsts.append(output_concat_uvfits_pyuvdata_name)
            
            cmd = 'rm -rf %s' % output_concat_uvfits_pyuvdata_name_lsts
            print(cmd)
            os.system(cmd)
            
            concat_uvfits(model_vis_uvfits_list_lsts,output_concat_uvfits_pyuvdata_name_lsts)
   
   #Delete all those intermediate files!
   for circle_number in range(0,n_circles):
      #circle 1 is the smallest
      radius = (min_arm_length_m + circle_number * radial_spacing) 
      diameter = radius * 2.
      diameter_cm = int(diameter*100.)
      #work out circumference 
      circumference = math.pi * diameter
      #print diameter
      n_angles = int(round(circumference * n_ants_per_m_of_circumference))
      angle_increment = (2.*math.pi)/n_angles
      #(remember only need half of them!)
      angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
      #angle_array_deg = angle_array_rad / math.pi * 180.
      #print angle_array_deg
      for angle_rad in angle_array_rad:
        angle_deg = int(angle_rad/np.pi*180.)
        sub_array_string = "sub_%03d_%03d" % (diameter_cm,angle_deg)
        for pol in pol_list:
          for lst in lst_list:
            lst_deg = (float(lst)/24.)*360.
            for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
               model_vis_name_base = "%s_LST_%03d_%s_%s_MHz" % (sub_array_string,lst_deg,pol,int(freq_MHz))
               if 'noise' in signal_type_list:
                  model_vis_name_base += '_N'
               if 'diffuse' in signal_type_list:
                  model_vis_name_base += '_D_%s' % sky_model
               if 'global' in signal_type_list:
                  model_vis_name_base += '_G' 
               if 'gain_errors' in signal_type_list:
                  model_vis_name_base += '_GE'
               
               out_vis_name = model_vis_name_base + '.vis'
               out_vis_uvfits_name = model_vis_name_base + '.uvfits'
               
               cmd = "rm -rf %s %s " % (out_vis_name,out_vis_uvfits_name)
               print(cmd)
               os.system(cmd)
               
            output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
            output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
               
            cmd = "rm -rf %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name)
            print(cmd)
            os.system(cmd)         

#main program
def simulate(lst_list,freq_MHz_list,pol_list,signal_type_list,sky_model,outbase_name,array_ant_locations_filename,array_label,EDA2_data=False):
   if EDA2_data:
      lst_list = [lst_list[0]]
   ##split up into chunks of 30 MHz
   #remainder = len(freq_MHz_list) % inst_bw
   #assert remainder == 0, "n_chan (%s) must be a multiple of the instantaneous BW (%s MHz)" % (len(freq_MHz_list),inst_bw)
   #n_freq_chunks = int(len(freq_MHz_list)/inst_bw)
   #
   #for chunk in range(n_freq_chunks):
   ##centre freq defined as the freq of the 16th chan (zero-based index 15)
   #   centre_freq_index = int(chunk*inst_bw + 15)
   #   centre_freq = freq_MHz_list[centre_freq_index]
   #   print centre_freq


   concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
   concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
   if 'noise' in signal_type_list:
       concat_output_name_base_X += '_N'
       concat_output_name_base_Y += '_N'
   if 'diffuse' in signal_type_list:
       concat_output_name_base_X += '_D_%s' % sky_model
       concat_output_name_base_Y += '_D_%s' % sky_model
   if 'global_unity' in signal_type_list:
       concat_output_name_base_X += '_GU'
       concat_output_name_base_Y += '_GU'
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("can't have diffuse and diffuse global at the same time.")
          sys.exit()
       else:
          concat_output_name_base_X += '_DG_%s' % sky_model
          concat_output_name_base_Y += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("can't have diffuse and diffuse angular at the same time.")
          sys.exit()
       else:
         concat_output_name_base_X += '_DA_%s' % sky_model
         concat_output_name_base_Y += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print("cant have global_EGES and global in signal_type_list")
          sys.exit()
       else:
          concat_output_name_base_X += '_G' 
          concat_output_name_base_Y += '_G'
   if 'global_EDGES' in signal_type_list:
       concat_output_name_base_X += '_ED' 
       concat_output_name_base_Y += '_ED'
   if 'gain_errors' in signal_type_list:
       concat_output_name_base_X += '_GE'
       concat_output_name_base_Y += '_GE'

 
   sky_averaged_diffuse_array_no_beam_lsts_filename =  "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base_X
   sky_averaged_diffuse_array_beam_X_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base_X
   sky_averaged_diffuse_array_beam_Y_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base_Y
   
   output_concat_vis_pyuvdata_name_lsts_X = "%s_concat_lsts.vis" % concat_output_name_base_X
   output_concat_uvfits_pyuvdata_name_lsts_X = "%s_concat_lsts.uvfits" % concat_output_name_base_X   
   output_concat_vis_pyuvdata_name_lsts_Y = "%s_concat_lsts.vis" % concat_output_name_base_Y
   output_concat_uvfits_pyuvdata_name_lsts_Y = "%s_concat_lsts.uvfits" % concat_output_name_base_Y 
   
   output_concat_vis_pyuvdata_name_lsts_X_sub = "%s_concat_lsts_sub.vis" % concat_output_name_base_X
   output_concat_uvfits_pyuvdata_name_lsts_X_sub = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_X   
   output_concat_vis_pyuvdata_name_lsts_Y_sub = "%s_concat_lsts_sub.vis" % concat_output_name_base_Y
   output_concat_uvfits_pyuvdata_name_lsts_Y_sub = "%s_concat_lsts_sub.uvfits" % concat_output_name_base_Y 
   
   n_chan = len(freq_MHz_list)
   n_lsts = len(lst_list)
   #Do this stuff for each lst and each freq:     
   #model_sky_uvfits_list_X_lsts = []
   #model_global_signal_uvfits_list_X_lsts = []
   #model_noise_uvfits_list_X_lsts = []
   #model_sky_uvfits_list_Y_lsts = []
   #model_global_signal_uvfits_list_Y_lsts = []
   #model_noise_uvfits_list_Y_lsts = []
   model_vis_uvfits_list_X_lsts = []
   model_vis_uvfits_list_Y_lsts = []
   model_vis_uvfits_list_X_lsts_sub = []
   model_vis_uvfits_list_Y_lsts_sub = []
   
   sky_averaged_diffuse_array_beam_X_lsts = np.full((n_chan),0.0)
   sky_averaged_diffuse_array_beam_Y_lsts = np.full((n_chan),0.0)
   sky_averaged_diffuse_array_no_beam_lsts = np.full((n_chan),0.0)
   
   for lst in lst_list:
      #model_sky_vis_list_X = []
      #model_global_signal_vis_list_X = []
      #model_sky_vis_list_Y = []
      #model_global_signal_vis_list_Y = []
      #model_sky_uvfits_list_X = []
      #model_global_signal_uvfits_list_X = []
      #model_noise_uvfits_list_X = []
      #model_sky_uvfits_list_Y = []
      #model_global_signal_uvfits_list_Y = []
      #model_noise_uvfits_list_Y = []
      model_vis_uvfits_list_X = []
      model_vis_uvfits_list_Y = []
      model_vis_uvfits_list_X_sub = []
      model_vis_uvfits_list_Y_sub = []
      
      print("lst %s hrs" % lst)
      lst_deg = (float(lst)/24.)*360.
      
      print("lst %s deg" % lst_deg)
      for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
         print(freq_MHz)
         #and datetime of observation (eventually do this for many dates)
         year=2000
         month=1
         day=1
         hour=np.floor(float(lst))
         minute=np.floor((float(lst)-hour) * 60.)
         second=((float(lst)-hour) * 60. - minute) * 60.
         
         date_time_string = '%s_%02d_%02d_%02d_%02d_%02d' % (year,float(month),float(day),hour,minute,second)
         print(date_time_string)
         #for miriad time just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
         day_frac_plus_one = float(lst)/24. + 1
         miriad_uvgen_time_string = '00JAN%1.3f' % day_frac_plus_one
         print(miriad_uvgen_time_string)
         wavelength = 300./freq_MHz
         freq_GHz = freq_MHz/1000.
         
         gsm_hpx_fits_name = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
         reprojected_gsm_fitsname = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
         reprojected_gsm_im_name = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)

               
         global_signal_hpx_fits_name = "global_signal_map_LST_%03d_%0.3f_MHz_hpx.fits" % (lst_deg,freq_MHz)
         unity_sky_hpx_fits_name = "unity_sky_map_LST_%03d_%0.3f_MHz_hpx.fits" % (lst_deg,freq_MHz)
         reprojected_global_signal_fitsname = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (lst_deg,freq_MHz)
         reprojected_unity_sky_fitsname = "unity_sky_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (lst_deg,freq_MHz)
         reprojected_global_signal_im_name = "global_signal_map_LST_%03d_%0.3f_MHz_reprojected.im" % (lst_deg,freq_MHz)
         reprojected_unity_sky_im_name = "unity_sky_map_LST_%03d_%0.3f_MHz_reprojected.im" % (lst_deg,freq_MHz)
         
         #lna_impedance_aavs1_filename = "/md0/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
         
         base_vis_name = "%s_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
         unity_sky_base_vis_name = "%s_unity_sky_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
         unity_sky_out_vis_name = "%s_unity_sky_out_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
         unity_sky_out_uvfits_name = "%s_unity_sky_out_LST_%03d_%0.3f_MHz.uvfits" % (array_label,lst_deg,freq_MHz)
         #eda_model_uvfits_name = "eda_model_LST_%03d_%s_MHz.uvfits" % (lst_deg,int(freq_MHz))
         eda_model_no_source_image_name = "%s_no_src_LST_%03d_%0.3f_MHz.map" % (array_label,lst_deg,freq_MHz)
         eda_model_no_source_beam_name = "%s_no_src_LST_%03d_%0.3f_MHz.beam" % (array_label,lst_deg,freq_MHz)
         eda_model_no_source_image_fits_name = "%s_no_src_LST_%03d_%0.3f_MHz.fits" % (array_label,lst_deg,freq_MHz)
         
         #for gain errors:
         one_jy_source_vis_name =  "%s_one_jy_source_LST_%03d_%0.3f_MHz.vis" % (array_label,lst_deg,freq_MHz)
         gain_solutions_name_amp = "%s_gains_%03d_%0.3f_MHz_amp.txt" % (array_label,lst_deg,freq_MHz)
         gain_solutions_name_phase = "%s_gains_%03d_%0.3f_MHz_phase.txt" % (array_label,lst_deg,freq_MHz)
         
         #if apply_normalisation:
         #   #get the lna impedance for this frequency
         #   with open(lna_impedance_aavs1_filename) as f:
         #      lines = f.readlines()
         #   for line in lines:
         #      freq = line.split()[0]
         #      if float(freq) == float(freq_MHz):
         #         lna_z_real = float(line.split()[1])
         #         lna_z_imag = float(line.split()[2])
         #
         #   #print "%s %s" % (lna_z_real,lna_z_imag)
         #   lna_z_complex=complex(lna_z_real,lna_z_imag)
         #   print "LNA impedance is %s" % lna_z_complex
         #
         #   #w = 2pi*freq
         #   w = 2.*np.pi*float(freq_MHz)*10**(6)
         #
         #   #N_factor = (-j4pi/(w u0))*Z_lna
         #   normalisation_factor = (4.*np.pi*1j / (mu_0 * w)) * lna_z_complex
         
         #Generate model visibilities
         #1. Get a sky model (healpix map of GSM2017) - easy
         #2. Get a beam model in healpix (not entirely straightforward, but I think I have done this before) just use short dipole model.
         #3. Generate a point source model (I don't understand how this is done)
         #4. Generate a global EoR signal (I think I understand how to do this from Caths eqn 3 etc
          
         #Get EDA antenna positions (you have these)
         
         #Use zenith pointing centre and phase centre of RA = 4h, Dec= -27 degrees
         
         #Generate simulated visibilities ... with miriad to start?
         #uvgen to make a blank visibility set (no sources, but with correct observing parameters)
         #uvmodel with the 'add' option to take the pb attenuated gsm image and add it to the visibility file from uvgen
         
         #2. Sky model
         #Need a pb-attenuated gsm image in slant-orthographic projection
         # - obtain gsm image at desired frequency and datetime (healpix)
         if sky_model == 'gsm2016':
            #forget about all this observer stuff, just generate the all-sky map and let reproject take care of the rest 
            #ov = GSMObserver2016()
            gsm = GlobalSkyModel2016() 
         elif  sky_model == 'gsm':
            #ov = GSMObserver()
            gsm = GlobalSkyModel()
         elif sky_model == 'gmoss':
            gmoss_data_text_file_name = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/GMOSS_sky_spectra.txt"
            #from Mayuri: The data is in HEALPIX Nested Scheme with NSIDE 16 and 5 degrees resolution. 
            #It contains GMOSS spectra generated at the 3072 pixels. The first row gives the frequency [GHz] 
            #going from 0.04 GHz to 0.2 GHz in steps of 0.001 GHz. Sky spectra are in Kelvin units
            sky_model_lines_list = []
            with open(gmoss_data_text_file_name, 'r') as f:
               sky_model_lines = f.readlines()
            for line in sky_model_lines:
               sky_model_lines_list.append(line.split())
            sky_model_data_array = np.full((len(sky_model_lines_list),len(sky_model_lines_list[0])),np.nan)
            for line_values_index,line_vales in enumerate(sky_model_lines_list):
               sky_model_data_array[line_values_index,:] = [float(i.strip()) for i in sky_model_lines_list[line_values_index]]
            gmoss_freqs_GHz = sky_model_data_array[0]
            #gmoss_freqs_GHz = np.asarray(gmoss_freqs_GHz)
            gmoss_freqs_MHz = gmoss_freqs_GHz*1000.
         else:
            print('unknown sky_model')
            sys.exit()
         #ov.lon = longitude_degrees
         #ov.lat = latitude_degrees
         #ov.elev = elevation_m
         #
         #ov.date = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
         
         
         #1. get 'blank' visibilities in miriad and make an image (1 Jy pt source) to use as a template for the output projection
         if generate_new_vis:
         
            cmd = "rm -rf %s %s" % (unity_sky_base_vis_name,base_vis_name)
            print(cmd)
            os.system(cmd)
           
            #need to change this so that each uvgen call has phase centre at zenith (RA = LST)
            cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (array_ant_locations_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,base_vis_name)
            print(cmd)
            os.system(cmd)
         
            #and for unity vis (for use later on in solving for T_sky)
            cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (array_ant_locations_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,unity_sky_base_vis_name)
            print(cmd)
            os.system(cmd)
            
                        
            #output the no source uvfits file for checking
            #cmd = "fits in=%s out=%s op=uvout" % (eda_model_vis_name,eda_model_uvfits_name)
            #print(cmd)
            #os.system(cmd)
            
            cmd = "rm -rf %s %s %s " % (eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "invert vis=%s map=%s beam=%s imsize=%s cell=%s stokes=xx options=mfs" % (base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,template_imsize,template_cell_size_asec)
            print(cmd)
            os.system(cmd)
         
            cmd = "fits in=%s out=%s op=xyout" % (eda_model_no_source_image_name,eda_model_no_source_image_fits_name)
            print(cmd)
            os.system(cmd)
            
            if 'gain_errors' in signal_type_list:
            
               cmd = "rm -rf %s" % one_jy_source_vis_name
               print(cmd)
               os.system(cmd)
   
               # randall's recipe for gain errors
               #rm -rf "$gainvis" "$noisevis"
               ## make 1Jy source with antenna-based amp/phase errors, no noise
               #$MIRCAT/point.source
               #uvgen ant=antfile.txt baseunit=-3.33564 corr=$corr freq=$freq radec=$lst,$dec harange=0,$ha_max,0.0005556 lat=-26.7 source=/tmp/source1Jy.txt out="$gainvis" time=00JAN1$fakedayfrac gnoise=10 pnoise=90 | grep noise
               cmd = "uvgen source=$MIRCAT/point.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx gnoise=10 pnoise=90  " % (array_ant_locations_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,one_jy_source_vis_name)
               print(cmd)
               os.system(cmd)
               ## selfcal this to solve for the gain errors
               #selfcal vis="$gainvis" flux=1 options=amplitude,noscale
               cmd = "selfcal vis=%s flux=1 options=amplitude,noscale" % (one_jy_source_vis_name)
               print(cmd)
               os.system(cmd)
               ## write out the solutions
               #gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_amp.txt
               #gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_pha.txt yaxis=phase
               
               
               cmd = "rm -rf %s %s" % (gain_solutions_name_amp,gain_solutions_name_phase)
               print(cmd)
               os.system(cmd)
               
               cmd = "gpplt vis=%s log=%s" % (one_jy_source_vis_name,gain_solutions_name_amp)
               print(cmd)
               os.system(cmd)
               
               cmd = "gpplt vis=%s log=%s yaxis=phase" % (one_jy_source_vis_name,gain_solutions_name_phase)
               print(cmd)
               os.system(cmd)
            
               #collect these up and plot them later on
            
         if generate_new_hpx:
            cmd = "rm -rf %s %s" % (gsm_hpx_fits_name,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)
            if sky_model == 'gmoss':
               gmoss_freq_index, gmoss_freq_value = find_nearest(gmoss_freqs_MHz,freq_MHz)
               #print "gmoss_freq_index, gmoss_freq_value %s %s" % (gmoss_freq_index, gmoss_freq_value)
               sky_model_data_nested = sky_model_data_array[:,gmoss_freq_index][1:]
               hp.write_map(gsm_hpx_fits_name,sky_model_data_nested,coord='G',nest=True)
               gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
            else:
               gsm_map = gsm.generate(freq_MHz)
               #print np.max(gsm_map)
               hp.write_map(gsm_hpx_fits_name,gsm_map,coord='G')
               
            #Now repeat this for each of the EDA2s 32 fine chans
            centre_freq = float(freq_MHz)
            fine_chan_width_MHz = fine_chan_width_Hz/1000000.
            for fine_chan_index in range(0,32):
               if EDA2_data:
                  freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
                  #dont reverse the fine chan index here - that is done later in calibrate_eda2...
                  #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               else:
                  freq_MHz_fine_chan = centre_freq     
               wavelength_fine_chan = 300./float(freq_MHz_fine_chan)  
               
               gsm_hpx_fits_name_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz_fine_chan)
               reprojected_gsm_fitsname_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz_fine_chan)
               reprojected_gsm_im_name_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz_fine_chan)
               
               
               cmd = "rm -rf %s %s" % (gsm_hpx_fits_name_fine_chan,reprojected_gsm_im_name_fine_chan)
               print(cmd)
               os.system(cmd)
               
               #only gsm supported
               gsm_map = gsm.generate(freq_MHz_fine_chan)
               #print np.max(gsm_map)
               hp.write_map(gsm_hpx_fits_name_fine_chan,gsm_map,coord='G')
               print("wrote %s" % gsm_hpx_fits_name_fine_chan)
               
            
         else:
            if sky_model == 'gmoss':
               sky_model_data_nested = hp.read_map(gsm_hpx_fits_name,nest=True)
               gsm_map = hp.reorder(sky_model_data_nested, n2r=True)
            else:
               gsm_map = hp.read_map(gsm_hpx_fits_name)
         
         
         
         #plot?
         if plot_gsm_map_hpx:
            #plot
            plt.clf()
            map_title="GSM from MWA at %.0f:%.0f:%.0f and %0.3f MHz" % (hour,minute,second,freq_MHz)
            #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=gsm_map,half_sky=False,title=map_title)
            #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
            fig_name="%s_map_%s_%0.3f_MHz.png" % (sky_model,date_time_string,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print("saved %s" % fig_name)
         
         
         #Miriad doesn't seem to be able to import the hpx file
         #Try using reproject_from_healpix
         #output_projection is the header of the pt source made above
         if os.path.isfile(eda_model_no_source_image_fits_name) and os.access(eda_model_no_source_image_fits_name, os.R_OK):
            hdulist = pyfits.open(eda_model_no_source_image_fits_name)
         else:
            print("Either file %s is missing or is not readable" % eda_model_no_source_image_fits_name)
            #continue        
         
         data=hdulist[0].data[0,0,:,:]
         no_source_header=hdulist[0].header
         pix_size_deg = float(no_source_header['CDELT1'])
         pix_area_deg_sq = pix_size_deg*pix_size_deg
         pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
         
         #needs cooordsys keyword
         #pt_source_header['COORDSYS'] = 'icrs'
         
         #print(pt_source_header)
         
         #del pt_source_header[8]
         #del pt_source_header[8]
         del no_source_header['history']
                         
         #print(pt_source_header)
         
         target_wcs = WCS(no_source_header)
         
         target_wcs=target_wcs.dropaxis(2)
         target_wcs=target_wcs.dropaxis(2)
                         
         #hdu_hpx = pyfits.open(gsm_hpx_fits_name)[1]
         ##hdu_hpx.info()
         #hpx_header = hdu_hpx.header
         #print(hpx_header)
         
         #for centre chan (original behaviour)
         
         hdu_gsm = fits.open(gsm_hpx_fits_name)[1]
         #print hdu_gsm.header
         #see bottom for accessing table data https://python4astronomers.github.io/astropy/fits.html
         input_tuple = (gsm_map,'G')
         
         reprojected_gsm_map,footprint = reproject_from_healpix(hdu_gsm, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
         #reprojected_gsm_map,footprint = reproject_from_healpix((data_gsm,'galactic'), target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',nested=False)
         #reprojected_gsm_map,footprint = reproject_from_healpix(input_tuple, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',nested=False)
         

         
         ##calc sky-average temp, no beam - 
         #Dont use this as input for visibilities, use the beam weighted on calculated below in the beamy bit
         sky_average_temp_no_beam = np.nanmean(reprojected_gsm_map)
         #print sky_average_temp_no_beam
         sky_averaged_diffuse_array_no_beam_lsts[freq_MHz_index] += sky_average_temp_no_beam
         
         
         
         
         
         #write the reprojected gsm maps to fits
         pyfits.writeto(reprojected_gsm_fitsname,reprojected_gsm_map,clobber=True)
         pyfits.update(reprojected_gsm_fitsname,reprojected_gsm_map,header=no_source_header)
         print("wrote image %s" %  reprojected_gsm_fitsname)

         
                  
         #Do this GSM map stuff here as it doesn't depend on pol
         cmd = "rm -rf %s" % (reprojected_gsm_im_name)
         print(cmd)
         os.system(cmd)     
         
         cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_fitsname,reprojected_gsm_im_name)
         print(cmd)
         os.system(cmd)


                 
         #uvmodel requires the model to be in Jy/pix
         #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
         
         scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
         print("scale map by %s to get to Jy/pix" % scale)
         
         reprojected_gsm_im_Jy_per_pix_name =  "%s_%s_%0.3f_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)

         
         cmd = "rm -rf %s" % (reprojected_gsm_im_Jy_per_pix_name)
         print(cmd)
         os.system(cmd)
               
         cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name)
         print(cmd)
         os.system(cmd)

        
         ########
         #Repeat the above for the global signal
         #make an all sky global signal map here too:
         if 'global_EDGES' in signal_type_list:
            s_21_value = s_21_array_EDGES[freq_MHz_index]
         else:
            s_21_value = s_21_array[freq_MHz_index]
         s_21_hpx_map = gsm_map * 0.0 + s_21_value
         
         jy_to_K = (wavelength**2) / (2. * k * 1.0e26) 
         
         unity_sky_value = 1. * jy_to_K
         #What value do you actually need to put in here to get the desired result .... I think it is 1 / Jy_to_k?
         #play until you get a gradient of one in x_y_plot!
         
         unity_sky_hpx_map = gsm_map * 0.0 + unity_sky_value
         
         cmd = "rm -rf %s %s %s %s" % (global_signal_hpx_fits_name,reprojected_global_signal_im_name,unity_sky_hpx_fits_name,reprojected_unity_sky_im_name)
         print(cmd)
         os.system(cmd)
            
         hp.write_map(global_signal_hpx_fits_name,s_21_hpx_map,coord='C')
         hp.write_map(unity_sky_hpx_fits_name,unity_sky_hpx_map,coord='C')
         
         #plot?
         if plot_global_signal_map_hpx:
            #plot
            plt.clf()
            map_title="Global signal from MWA at %.0f:%.0f:%.0f and %0.3f MHz" % (hour,minute,second,freq_MHz)
            #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=s_21_hpx_map,half_sky=False,title=map_title)
            #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
            fig_name="global_signal_map_%s_%0.3f_MHz.png" % (date_time_string,freq_MHz)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print("saved %s" % fig_name)
              
         reprojected_global_signal_map,footprint = reproject_from_healpix(global_signal_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
         reprojected_unity_sky_map,footprint = reproject_from_healpix(unity_sky_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
         
         
         #write the map to fits
         pyfits.writeto(reprojected_global_signal_fitsname,reprojected_global_signal_map,clobber=True)
         pyfits.update(reprojected_global_signal_fitsname,reprojected_global_signal_map,header=no_source_header)
         print("wrote image %s" %  reprojected_global_signal_fitsname)

         pyfits.writeto(reprojected_unity_sky_fitsname,reprojected_unity_sky_map,clobber=True)
         pyfits.update(reprojected_unity_sky_fitsname,reprojected_unity_sky_map,header=no_source_header)
         print("wrote image %s" %  reprojected_unity_sky_fitsname)
                  
         cmd = "rm -rf %s %s" % (reprojected_global_signal_im_name,reprojected_unity_sky_im_name)
         print(cmd)
         os.system(cmd)     
         
         cmd = "fits in=%s out=%s op=xyin" % (reprojected_global_signal_fitsname,reprojected_global_signal_im_name)
         print(cmd)
         os.system(cmd)
         
         cmd = "fits in=%s out=%s op=xyin" % (reprojected_unity_sky_fitsname,reprojected_unity_sky_im_name)
         print(cmd)
         os.system(cmd)
         
         #uvmodel requires the model to be in Jy/pix
         #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
         
         reprojected_global_signal_im_Jy_per_pix_name =  "global_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (date_time_string,freq_MHz)
         reprojected_unity_sky_im_Jy_per_pix_name =  "unity_sky_%s_%0.3f_MHz_reproj_Jy_pix.im" % (date_time_string,freq_MHz)
         
         
         cmd = "rm -rf %s %s" % (reprojected_global_signal_im_Jy_per_pix_name,reprojected_unity_sky_im_Jy_per_pix_name)
         print(cmd)
         os.system(cmd)
               
         cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_global_signal_im_name,reprojected_global_signal_im_Jy_per_pix_name)
         print(cmd)
         os.system(cmd)

         cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_unity_sky_im_name,reprojected_unity_sky_im_Jy_per_pix_name)
         print(cmd)
         os.system(cmd)
                   
         #do for all pols:
         for pol in pol_list:
            #model_sky_vis = "eda_model_plus_sky_LST_%03d_%s_pol_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
            #model_global_signal_vis = "eda_model_plus_global_LST_%03d_%s_pol_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
            #model_sky_uvfits = "eda_model_plus_sky_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))
            #model_global_signal_uvfits = "eda_model_plus_global_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz)) 
            #eda_model_noise_vis_name = "eda_model_noise_LST_%03d_%s_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
            #eda_model_noise_uvfits_name = "eda_model_noise_LST_%03d_%s_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))        
            model_vis_name_base = "%s_LST_%03d_%s_%0.3f_MHz" % (array_label,lst_deg,pol,freq_MHz)

               
            if use_analytic_beam:
               if pol=='X':
                  beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                  beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
               else:
                  beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                  beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
            else:
                  beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
                  
    
            #power_pattern_average_interp_sin_im_name = 'power_pattern_average_%s_%s_MHz_interp_sin.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_im_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_fits_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.fits' % (pol,int(freq_MHz))

            cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
            print(cmd)
            os.system(cmd)

            cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname_no_cos_za)
            print(cmd)
            os.system(cmd)
                                         
            beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
            beam_image_sin_projected_im_name_no_cos_za = 'beam_image_sin_projected_%s_%0.3f_MHz_no_cos_za.im' % (pol,freq_MHz)
            beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_fits_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.fits' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_im_name_no_cos_za =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid_no_cos_za.im' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid_no_cos_za.fits' % (pol,freq_MHz)
           
           
            cmd = "rm -rf %s %s %s %s %s %s" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name, beam_image_sin_projected_puthd_fits_name,beam_image_sin_projected_regrid_gsm_im_name_no_cos_za, beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za, beam_image_sin_projected_im_name_no_cos_za)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname_no_cos_za,beam_image_sin_projected_im_name_no_cos_za)
            print(cmd)
            os.system(cmd)
            
            #put in the correct ra in the header (ra = lst for zenith) 
            #puthd in="$beam/crval1" value=$lst_degs
            cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
            print(cmd)
            os.system(cmd)         

            cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name_no_cos_za,lst_deg)
            print(cmd)
            os.system(cmd) 
            
            #write out as a fits file to check header
            cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
            print(cmd)
            os.system(cmd)
            
            #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
            cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)   

            cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name_no_cos_za,beam_image_sin_projected_regrid_gsm_im_name_no_cos_za,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)
                        
         
            #Sweet! finally have a gsm and a beam.
            
            #diffuse
            #now multiply them to get the apparent sky and put that into uvmodel, 
         
            apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
            apparent_sky_fits_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.fits" % (lst_deg,pol,freq_MHz)
            apparent_sky_im_Tb_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz_Tb.im" % (lst_deg,pol,freq_MHz)
            apparent_sky_im_Tb_fits_name = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz_Tb.fits" % (lst_deg,pol,freq_MHz)
            
            cmd = "rm -rf %s %s %s %s " % (apparent_sky_im_name,apparent_sky_im_Tb_name,apparent_sky_im_Tb_fits_name,apparent_sky_fits_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_name,apparent_sky_fits_name)
            print(cmd)
            os.system(cmd)
            
            #repeat the apparent sky simulation for each of the 32 EDA2 fine chans. Beam doesnt change much so just use same beam
            #########################################
            #########################################
            for fine_chan_index in range(0,32):
               if EDA2_data:
                  freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
                  #dont reverse the fine chan index here - that is done later in calibrate_eda2...
                  #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               else:
                  freq_MHz_fine_chan = centre_freq     
               wavelength_fine_chan = 300./float(freq_MHz_fine_chan) 
               
               gsm_hpx_fits_name_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz_fine_chan)
               reprojected_gsm_fitsname_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz_fine_chan)
               reprojected_gsm_im_name_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz_fine_chan)
               
               hdu_gsm_fine_chan = fits.open(gsm_hpx_fits_name_fine_chan)[1]
               #print hdu_gsm.header
               #see bottom for accessing table data https://python4astronomers.github.io/astropy/fits.html
               gsm_map_fine_chan = hp.read_map(gsm_hpx_fits_name_fine_chan)
               input_tuple_fine_chan = (gsm_map_fine_chan,'G')
               
               reprojected_gsm_map_fine_chan,footprint_fine_chan = reproject_from_healpix(hdu_gsm_fine_chan, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
               #reprojected_gsm_map,footprint = reproject_from_healpix((data_gsm,'galactic'), target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',nested=False)
               #reprojected_gsm_map,footprint = reproject_from_healpix(input_tuple, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',nested=False)
               
               #write the reprojected gsm maps to fits
               pyfits.writeto(reprojected_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,clobber=True)
               pyfits.update(reprojected_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,header=no_source_header)
               print("wrote image %s" %  reprojected_gsm_fitsname_fine_chan)
             
               #Do this GSM map stuff here as it doesn't depend on pol
               cmd = "rm -rf %s" % (reprojected_gsm_im_name_fine_chan)
               print(cmd)
               os.system(cmd)     
               
               cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_fitsname_fine_chan,reprojected_gsm_im_name_fine_chan)
               print(cmd)
               os.system(cmd)
         
               #uvmodel requires the model to be in Jy/pix
               #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
               
               scale_fine_chan = (2. * k * 1.0e26 * pix_area_sr) / (wavelength_fine_chan**2)
               print("scale map by %s to get to Jy/pix" % scale_fine_chan)
               
               reprojected_gsm_im_Jy_per_pix_name_fine_chan =  "%s_%s_%0.3f_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz_fine_chan)

               
               cmd = "rm -rf %s" % (reprojected_gsm_im_Jy_per_pix_name_fine_chan)
               print(cmd)
               os.system(cmd)
                     
               cmd = "maths exp=%s*%s out=%s " % (scale_fine_chan,reprojected_gsm_im_name_fine_chan,reprojected_gsm_im_Jy_per_pix_name_fine_chan)
               print(cmd)
               os.system(cmd)
      
               #do for all pols:
               for pol in pol_list:
                     
                  if use_analytic_beam:
                     if pol=='X':
                        beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                        beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                     else:
                        beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                        beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
                  else:
                        beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
                                                   
                  beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
                  
                  beam_image_sin_projected_regrid_gsm_im_name_fine_chan =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz_fine_chan)
                  beam_image_sin_projected_regrid_gsm_fits_name_fine_chan =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.fits' % (pol,freq_MHz_fine_chan)

                 
                  cmd = "rm -rf %s " % (beam_image_sin_projected_regrid_gsm_im_name_fine_chan)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name_fine_chan,reprojected_gsm_im_name_fine_chan)
                  print(cmd)
                  os.system(cmd)   
      
                  #now multiply them to get the apparent sky and put that into uvmodel, 
               
                  apparent_sky_im_name_fine_chan = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz_fine_chan)
                  apparent_sky_fits_name_fine_chan = "apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.fits" % (lst_deg,pol,freq_MHz_fine_chan)

                  
                  cmd = "rm -rf %s %s" % (apparent_sky_im_name_fine_chan,apparent_sky_fits_name_fine_chan)
                  print(cmd)
                  os.system(cmd)
               
                  cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name_fine_chan,reprojected_gsm_im_Jy_per_pix_name_fine_chan,apparent_sky_im_name_fine_chan)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_name_fine_chan,apparent_sky_fits_name_fine_chan)
                  print(cmd)
                  os.system(cmd)
            
            
            #########################################
            #########################################
            #THis is done above per fine chan
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name,apparent_sky_im_Tb_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_Tb_name,apparent_sky_im_Tb_fits_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name)
            print(cmd)
            os.system(cmd)

            cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_regrid_gsm_im_name_no_cos_za,beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za)
            print(cmd)
            os.system(cmd)
              
            #just for sky-averaged temp plotting:
            with fits.open(apparent_sky_im_Tb_fits_name) as hdu_list:
               data = hdu_list[0].data
            
            #this should be beam-weighted average
            #sky_average_temp_beam = np.nanmean(data)
            with fits.open(beam_image_sin_projected_regrid_gsm_fits_name) as hdu_list:
            #with fits.open(beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za) as hdu_list:
               beam_image_sin_projected_regrid_gsm = hdu_list[0].data
            #these probably have a different number of pixels blanked with nan and that is why you get the 
            #wrong answer for the beam weighted average. yes matter ....
            
            #This gives the wrong answer..... I think it has to do with the cos(az) term in the beam or something....
            #I think the way to go is to make the beams fresh in hpx (as I do in solve_for_tsky_from uvfits), then multiply by 
            #the healpix gsm, then reproject with reproject from healpy ....
            
            sum_of_beam_weights = np.nansum(beam_image_sin_projected_regrid_gsm)
            
            sky_average_temp_beam_weighted = np.nansum(data) / sum_of_beam_weights

            if pol=='X':
               sky_averaged_diffuse_array_beam_X_lsts[freq_MHz_index] += sky_average_temp_beam_weighted
            else:
               sky_averaged_diffuse_array_beam_Y_lsts[freq_MHz_index] += sky_average_temp_beam_weighted

            
            #Need this above for diffuse global, so move that bit here:
            #Not sure if I am introducing double beam effects now .....
            ##
            gsm_global_hpx_fits_name = "%s_global_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_global_fitsname = "%s_global_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_global_im_name = "%s_global_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)
            
            gsm_angular_hpx_fits_name = "%s_angular_map_LST_%03d_%0.3f_MHz_hpx.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_angular_fitsname = "%s_angular_map_LST_%03d_%0.3f_MHz_reprojected.fits" % (sky_model,lst_deg,freq_MHz)
            reprojected_gsm_angular_im_name = "%s_angular_map_LST_%03d_%0.3f_MHz_reprojected.im" % (sky_model,lst_deg,freq_MHz)
            
            reprojected_gsm_global_map = reprojected_gsm_map * 0. + sky_average_temp_beam_weighted
            reprojected_gsm_angular_map = reprojected_gsm_map - sky_average_temp_beam_weighted
            
            pyfits.writeto(reprojected_gsm_global_fitsname,reprojected_gsm_global_map,clobber=True)
            pyfits.update(reprojected_gsm_global_fitsname,reprojected_gsm_global_map,header=no_source_header)
            print("wrote image %s" %  reprojected_gsm_global_fitsname)
            
            pyfits.writeto(reprojected_gsm_angular_fitsname,reprojected_gsm_angular_map,clobber=True)
            pyfits.update(reprojected_gsm_angular_fitsname,reprojected_gsm_angular_map,header=no_source_header)
            print("wrote image %s" %  reprojected_gsm_angular_fitsname)
            
            cmd = "rm -rf %s %s" % (reprojected_gsm_global_im_name,reprojected_gsm_angular_im_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_global_fitsname,reprojected_gsm_global_im_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_angular_fitsname,reprojected_gsm_angular_im_name)
            print(cmd)
            os.system(cmd)
            
            reprojected_gsm_global_im_Jy_per_pix_name =  "%s_DG_%s_%0.3f_MHz_%s_pol_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz,pol)
            reprojected_gsm_angular_im_Jy_per_pix_name =  "%s_DA_%s_%0.3f_MHz_%s_pol_reprojected_Jy_pix.im" % (sky_model,date_time_string,freq_MHz,pol)
            
            cmd = "rm -rf %s %s" % (reprojected_gsm_global_im_Jy_per_pix_name,reprojected_gsm_angular_im_Jy_per_pix_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_global_im_name,reprojected_gsm_global_im_Jy_per_pix_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_angular_im_name,reprojected_gsm_angular_im_Jy_per_pix_name)
            print(cmd)
            os.system(cmd)         
            ##
                 
            
            #Repeat for global signal
            apparent_global_signal_im_name = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
            apparent_global_signal_fits_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz.fits" % (lst_deg,pol,freq_MHz)
            apparent_global_signal_im_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%0.3f_MHz_cropped.im" % (lst_deg,pol,freq_MHz)
         
            cmd = "rm -rf %s %s %s" % (apparent_global_signal_im_name,apparent_global_signal_fits_name_cropped,apparent_global_signal_im_name_cropped)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_global_signal_im_Jy_per_pix_name,apparent_global_signal_im_name)
            print(cmd)
            os.system(cmd)

            #Repeat for unity sky
            apparent_unity_sky_im_name = "apparent_unity_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)

            cmd = "rm -rf %s" % (apparent_unity_sky_im_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_unity_sky_im_Jy_per_pix_name,apparent_unity_sky_im_name)
            print(cmd)
            os.system(cmd)
                        
            ##crop the image
            ##output a cropped fits file
            #cmd = "fits in=%s out=%s region=quarter op=xyout" % (apparent_global_signal_im_name,apparent_global_signal_fits_name_cropped)
            #print(cmd)
            #os.system(cmd)            
            #
            ##read it back in
            #cmd = "fits in=%s out=%s region=quarter op=xyin" % (apparent_global_signal_fits_name_cropped,apparent_global_signal_im_name_cropped)
            #print(cmd)
            #os.system(cmd)
            
            #repeat for diffuse global
            apparent_sky_diffuse_global_im_name = "apparent_DG_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
 
            cmd = "rm -rf %s " % (apparent_sky_diffuse_global_im_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_global_im_Jy_per_pix_name,apparent_sky_diffuse_global_im_name)
            print(cmd)
            os.system(cmd)

            #repeat for diffuse angular
            apparent_sky_diffuse_angular_im_name = "apparent_angular_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
 
            cmd = "rm -rf %s " % (apparent_sky_diffuse_angular_im_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_angular_im_Jy_per_pix_name,apparent_sky_diffuse_angular_im_name)
            print(cmd)
            os.system(cmd)                      
               
            #then put into the vis 
            
            if 'noise' in signal_type_list:
            
               out_vis_name =  "%s_N_LST_%03d_%s_%0.3f_MHz.vis" % (array_label,lst_deg,pol,freq_MHz)
               model_vis_name_base += '_N'
               
               
               ##Put in the uvgen for the noise here (so we get different noise for X and Y pol (i.e. remove noise uvgen from above!!))
               #generate noise only visibilites
               #for systemp, just use sky temp 180 at 180, for JyperK, use Aeff from Table 2 of Wayth et al 2017 (EDA paper)
               systemp = T_180*(freq_MHz/180.0)**beta
               print('systemp %s K' % systemp)
               A_eff = Aeff_for_freq_MHz_list[freq_MHz_index]
               print('A_eff %s K' % A_eff)
               JperK = (2.0*k*10**26)/(eta*A_eff)
               print('JperK %s ' % JperK)
               SEFD = systemp * JperK
               print('SEFD %s ' % SEFD)
           
               #generate the noise-only uvfits file
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
            
               cmd = "uvgen source=$MIRCAT/no.source ant='%s' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 systemp=%s jyperk=%s out=%s stokes=xx  " % (array_ant_locations_filename,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec, harange_string, systemp, JperK, out_vis_name)
               print(cmd)
               os.system(cmd)
            
               #cmd = "rm -rf %s" % eda_model_noise_uvfits_name
               #print(cmd)
               #os.system(cmd)
               
               #cmd = "fits in=%s out=%s op=uvout" % (eda_model_noise_vis_name,eda_model_noise_uvfits_name)
               #print(cmd)
               #os.system(cmd)  
               
               base_vis_name = out_vis_name
            
            if 'diffuse' in signal_type_list:
            
               model_vis_name_base += '_D_%s' % sky_model
               out_vis_name = model_vis_name_base + '.vis'
               
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
            
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_sky_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
   
               #cmd = "rm -rf %s" % model_sky_uvfits
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "fits in=%s out=%s op=uvout" % (model_sky_vis,model_sky_uvfits)
               #print(cmd)
               #os.system(cmd)            
               
               base_vis_name = out_vis_name

            if 'diffuse_global' in signal_type_list:
            
               model_vis_name_base += '_DG_%s' % sky_model
               out_vis_name = model_vis_name_base + '.vis'
               
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
            
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_sky_diffuse_global_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
   
               #cmd = "rm -rf %s" % model_sky_uvfits
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "fits in=%s out=%s op=uvout" % (model_sky_vis,model_sky_uvfits)
               #print(cmd)
               #os.system(cmd)            
               
               base_vis_name = out_vis_name

            if 'diffuse_angular' in signal_type_list:
            
               model_vis_name_base += '_DA_%s' % sky_model
               out_vis_name = model_vis_name_base + '.vis'
               
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
            
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_sky_diffuse_angular_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
   
               #cmd = "rm -rf %s" % model_sky_uvfits
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "fits in=%s out=%s op=uvout" % (model_sky_vis,model_sky_uvfits)
               #print(cmd)
               #os.system(cmd)            
               
               base_vis_name = out_vis_name
                                             
            if 'global' in signal_type_list:
            
               model_vis_name_base += '_G'
               out_vis_name = model_vis_name_base + '.vis'
            
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
               
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_global_signal_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
               
               #cmd = "rm -rf %s" % model_global_signal_uvfits
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "fits in=%s out=%s op=uvout" % (model_global_signal_vis,model_global_signal_uvfits)
               #print(cmd)
               #os.system(cmd)
   
               
               #remove the previous base_vis
               cmd = "rm -rf %s" % base_vis_name
               print(cmd)
               os.system(cmd)
               
               base_vis_name = out_vis_name

            if 'global_EDGES' in signal_type_list:
            
               model_vis_name_base += '_ED'
               out_vis_name = model_vis_name_base + '.vis'
            
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
               
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_global_signal_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
               
               #cmd = "rm -rf %s" % model_global_signal_uvfits
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "fits in=%s out=%s op=uvout" % (model_global_signal_vis,model_global_signal_uvfits)
               #print(cmd)
               #os.system(cmd)
   
               
               #remove the previous base_vis
               cmd = "rm -rf %s" % base_vis_name
               print(cmd)
               os.system(cmd)
               
               base_vis_name = out_vis_name

            if 'global_unity' in signal_type_list:
            
               model_vis_name_base += '_GU'
               out_vis_name = model_vis_name_base + '.vis'
            
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
               
               cmd = "uvmodel vis=%s model=%s options=add,mfs out=%s" % (base_vis_name,apparent_unity_sky_im_name,out_vis_name)
               print(cmd)
               os.system(cmd)
               
               #remove the previous base_vis
               cmd = "rm -rf %s" % base_vis_name
               print(cmd)
               os.system(cmd)
               
               base_vis_name = out_vis_name              
               
               #cmd = "prthd in=%s " % out_vis_name
               #print(cmd)
               #os.system(cmd)
               #      
               #cmd = "uvlist vis=%s options=array,full" % out_vis_name
               #print(cmd)
               #os.system(cmd)
               
            if 'gain_errors' in signal_type_list:
               model_vis_name_base += '_GE'
               out_vis_name = model_vis_name_base + '.vis'
               
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
                     
               cmd = "gpcopy vis=%s out=%s options=relax" % (one_jy_source_vis_name,base_vis_name)
               print(cmd)
               os.system(cmd)
               #use uvcat so that the solutions are applied and become the 'data' (not corrected data)
               #uvcat vis=/tmp/tempvis.uv out="$outname"
               cmd = "uvcat vis=%s out=%s" % (base_vis_name,out_vis_name)
               print(cmd)
               os.system(cmd)
               
               #remove the previous base_vis
               cmd = "rm -rf %s" % base_vis_name
               print(cmd)
               os.system(cmd)
               
               base_vis_name = out_vis_name
               
            out_vis_uvfits_name = model_vis_name_base + '.uvfits' 
            out_image_subtr_vis_name = model_vis_name_base + '_sub.vis' 
            out_image_subtr_vis_uvfits_name = model_vis_name_base + '_sub.uvfits' 
            
            cmd = "rm -rf %s %s %s" % (out_vis_uvfits_name,out_image_subtr_vis_name,out_image_subtr_vis_uvfits_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=uvout options=nocal,nopol,nopass" % (out_vis_name,out_vis_uvfits_name)
            print(cmd)
            os.system(cmd)
            
            #add visname to concat list
            if pol=='X':
               #model_sky_vis_list_X.append(model_sky_vis)
               #model_global_signal_vis_list_X.append(model_global_signal_vis)
               #model_sky_uvfits_list_X.append(model_sky_uvfits)
               #model_global_signal_uvfits_list_X.append(model_global_signal_uvfits)
               #model_noise_uvfits_list_X.append(eda_model_noise_uvfits_name)
               model_vis_uvfits_list_X.append(out_vis_uvfits_name)
               model_vis_uvfits_list_X_sub.append(out_image_subtr_vis_uvfits_name)
            else:
               #model_sky_vis_list_Y.append(model_sky_vis)
               #model_global_signal_vis_list_Y.append(model_global_signal_vis)    
               #model_sky_uvfits_list_Y.append(model_sky_uvfits)
               #model_global_signal_uvfits_list_Y.append(model_global_signal_uvfits)      
               #model_noise_uvfits_list_Y.append(eda_model_noise_uvfits_name)  
               model_vis_uvfits_list_Y.append(out_vis_uvfits_name)
               model_vis_uvfits_list_Y_sub.append(out_image_subtr_vis_uvfits_name)

            #need to image the datasets here before concatenation, as miriad can't re-import the uv files after pyuvdata has messed with them in concatenation
            if do_image_and_subtr_from_simulated_data:
            
               image_basename = model_vis_name_base + '_sim_imaged'
            
               cmd = "rm -rf %s_xx.map %s.beam %s_xx.model %s_xx.restor %s_xx_sub.map %s_sub.beam %s_xx_sub.model %s_xx_sub.restor" % (image_basename,image_basename,image_basename,image_basename,image_basename,image_basename,image_basename,image_basename)
               print(cmd)
               os.system(cmd)
               
               #only image with long baselines (don't want to be subtracting global signal!)
               #uvrange(uvmin,uvmax) Select visibilities with uv radius between uvmin and uvmax (in kilowavelenghts). If only one value is given, uvmin is taken as zero
               uvmin = min_imaging_uvdist_wavelength/1000.
               uvmax = max_imaging_uvdist_wavelength/1000.
               
               cmd = "invert vis=%s map=%s_xx.map beam=%s.beam options=double imsize=512 stokes=xx robust=-0.5 cell=1800 select=uvrange\(%0.6f,%0.6f\)" % (out_vis_name,image_basename,image_basename,uvmin,uvmax)
               print(cmd)
               os.system(cmd)
               
               cmd = "clean map=%s_xx.map beam=%s.beam out=%s_xx.model niters=2000" % (image_basename,image_basename,image_basename)
               print(cmd)
               os.system(cmd)
               
               cmd = "restor model=%s_xx.model  beam=%s.beam map=%s_xx.map out=%s_xx.restor " % (image_basename,image_basename,image_basename,image_basename)
               print(cmd)
               os.system(cmd)
               
               #only need the .model one
               cmd = "rm -rf %s_xx.map %s.beam %s_xx.restor " % (image_basename,image_basename,image_basename)
               print(cmd)
               os.system(cmd)
               
               cmd = "uvmodel vis=%s model=%s_xx.model options=subtract,mfs out=%s" % (out_vis_name,image_basename,out_image_subtr_vis_name)
               print(cmd)
               os.system(cmd)
               
               cmd = "fits in=%s out=%s op=uvout options=nocal,nopol,nopass" % (out_image_subtr_vis_name,out_image_subtr_vis_uvfits_name)
               print(cmd)
               os.system(cmd)
            
               #image the subtr vis to see if it worked
               #cmd = "invert vis=%s map=%s_xx_sub.map beam=%s_sub.beam options=double imsize=512 stokes=xx robust=-0.5 cell=1800 select=uvrange\(%0.6f,%0.6f\)" % (out_image_subtr_vis_name,image_basename,image_basename,uvmin,uvmax)
               #print(cmd)
               #os.system(cmd)
               
               #cmd = "clean map=%s_xx_sub.map beam=%s_sub.beam out=%s_xx_sub.model niters=2000" % (image_basename,image_basename,image_basename)
               #print(cmd)
               #os.system(cmd)
            
               #cmd = "restor model=%s_xx_sub.model  beam=%s_sub.beam map=%s_xx_sub.map out=%s_xx_sub.restor " % (image_basename,image_basename,image_basename,image_basename)
               #print(cmd)
               #os.system(cmd)
               
               
            #delete all the intermediate images and vis that are no longer required
            #sys.exit()
            if do_cleanup_images_and_vis:
               cleanup_images_and_vis(array_label,lst,freq_MHz,pol)
      
      ###concat vis both pols for each freq
      #Don't need to do this concat anymore....
      for pol in pol_list:
         output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
         output_concat_vis_pyuvdata_name_sub = "%s_concat_freqs_sub.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name_sub = "%s_concat_freqs_sub.uvfits" % model_vis_name_base
         
         cmd = "rm -rf %s %s %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name,output_concat_vis_pyuvdata_name_sub,output_concat_uvfits_pyuvdata_name_sub)
         print(cmd)
         os.system(cmd)
         if pol=='X':
             concat_uvfits(model_vis_uvfits_list_X,output_concat_uvfits_pyuvdata_name)
             if do_image_and_subtr_from_simulated_data:
                concat_uvfits(model_vis_uvfits_list_X_sub,output_concat_uvfits_pyuvdata_name_sub)
         else:
             concat_uvfits(model_vis_uvfits_list_Y,output_concat_uvfits_pyuvdata_name)
             if do_image_and_subtr_from_simulated_data:
                concat_uvfits(model_vis_uvfits_list_Y_sub,output_concat_uvfits_pyuvdata_name_sub)
      
      #   model_sky_vis_list_string = ','.join(model_sky_vis_list)
      #   cmd = "uvaver vis=%s out=%s" % (model_sky_vis_list_string,output_concat_model_vis_name)
      #   print(cmd)
      #   os.system(cmd)      
      
      #   #gsm
      #   output_concat_model_pyuvdata_name = "eda_model_LST_%03d_%s_pol_concat.vis" % (lst_deg,pol)
      #   output_concat_model_uvfits_name = "eda_model_LST_%03d_%s_pol_concat.uvfits" % (lst_deg,pol)
      #   
      #   cmd = "rm -rf %s" % output_concat_model_uvfits_name
      #   print(cmd)
      #   os.system(cmd)
      #   cmd = "rm -rf %s" % output_concat_model_pyuvdata_name
      #   print(cmd)
      #   os.system(cmd)
      #   if pol=='X':
      #       concat_uvfits(model_sky_uvfits_list_X,output_concat_model_uvfits_name)
      #   else:
      #       concat_uvfits(model_sky_uvfits_list_Y,output_concat_model_uvfits_name)
      #   model_sky_vis_list_string = ','.join(model_sky_vis_list)
      #   cmd = "uvaver vis=%s out=%s" % (model_sky_vis_list_string,output_concat_model_vis_name)
      #   print(cmd)
      #   os.system(cmd)
      #
      #  #global signal
      #  global_output_concat_model_vis_name = "global_eda_model_LST_%03d_%s_pol_concat.vis" % (lst_deg,pol)
      #  global_output_concat_model_uvfits_name = "global_eda_model_LST_%03d_%s_pol_concat.uvfits" % (lst_deg,pol)
      #  cmd = "rm -rf %s" % global_output_concat_model_uvfits_name
      #  print(cmd)
      #  os.system(cmd)
      #  cmd = "rm -rf %s" % global_output_concat_model_vis_name
      #  print(cmd)
      #  os.system(cmd)
      #  if pol=='X':
      #      concat_uvfits(model_global_signal_uvfits_list_X,global_output_concat_model_uvfits_name)
      #  else:
      #      concat_uvfits(model_global_signal_uvfits_list_Y,global_output_concat_model_uvfits_name)
      #      
      #  Noise:
      #  noise_output_concat_model_vis_name = "noise_eda_model_LST_%03d_%s_pol_concat.vis" % (lst_deg,pol)
      #  noise_output_concat_model_uvfits_name = "noise_eda_model_LST_%03d_%s_pol_concat.uvfits" % (lst_deg,pol)
      #  cmd = "rm -rf %s" % noise_output_concat_model_uvfits_name
      # print(cmd)
      # os.system(cmd)
      # cmd = "rm -rf %s" % noise_output_concat_model_vis_name
      # print(cmd)
      # os.system(cmd)
      # if pol=='X':
      #      concat_uvfits(model_noise_uvfits_list_X,noise_output_concat_model_uvfits_name)
      #  else:
      #      concat_uvfits(model_noise_uvfits_list_Y,noise_output_concat_model_uvfits_name)
      #   
         if pol=='X':
            #model_sky_uvfits_list_X_lsts.append(output_concat_model_uvfits_name)
            #model_global_signal_uvfits_list_X_lsts.append(global_output_concat_model_uvfits_name)
            #model_noise_uvfits_list_X_lsts.append(noise_output_concat_model_uvfits_name)
            model_vis_uvfits_list_X_lsts.append(output_concat_uvfits_pyuvdata_name)
            if do_image_and_subtr_from_simulated_data:
               model_vis_uvfits_list_X_lsts_sub.append(output_concat_uvfits_pyuvdata_name_sub)
         else:   
            #model_sky_uvfits_list_Y_lsts.append(output_concat_model_uvfits_name)
            #model_global_signal_uvfits_list_Y_lsts.append(global_output_concat_model_uvfits_name)      
            #model_noise_uvfits_list_Y_lsts.append(noise_output_concat_model_uvfits_name)
            model_vis_uvfits_list_Y_lsts.append(output_concat_uvfits_pyuvdata_name)
            if do_image_and_subtr_from_simulated_data:
               model_vis_uvfits_list_Y_lsts_sub.append(output_concat_uvfits_pyuvdata_name_sub)
      #     
         ####################################################################################     
      #   model_global_vis_list_string = ','.join(model_global_signal_vis_list)
      #   cmd = "uvaver vis=%s out=%s" % (model_global_vis_list_string,global_output_concat_model_vis_name)
      #   print(cmd)
      #   os.system(cmd)
         
         #Combine gsm and global
         
         ##write out uvfits file
         #cmd = "fits in=%s out=%s op=uvout" % (output_concat_model_vis_name,output_concat_model_uvfits_name)
         #print(cmd)
         #os.system(cmd)
   
         #cmd = "fits in=%s out=%s op=uvout" % (global_output_concat_model_vis_name,global_output_concat_model_uvfits_name)
         #print(cmd)
         #os.system(cmd)
   
   
         #then can finally try to repeat Caths stuff 
   
         #plot_freq_GHz = '0.050'
         
         #gsm
         #plot Amp vs U(lambda)
   
   
         #uv_dist_plot_name="eda_model_LST_%s_%s_pol_amp_vs_uvdist.png" % (lst,pol)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uvdist,amp options=nobase select=-auto' % (uv_dist_plot_name,output_concat_model_vis_name)
         #print(cmd)
         #os.system(cmd)
         #uv_dist_plot_name="eda_model_LST_%s_%s_pol_amp_vs_uc_%s_GHz.png" % (lst,pol,plot_freq_GHz)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uc,amp options=nobase select=frequency\(%s\)' % (uv_dist_plot_name,output_concat_model_vis_name,plot_freq_GHz)
         #print(cmd)
         #os.system(cmd)
         #amp_freq_plot_name="eda_model_LST_%s_%s_pol_amp_vs_freq.png" % (lst,pol)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=freq,amp options=nobase ' % (amp_freq_plot_name,output_concat_model_vis_name)
         #print(cmd)
         #os.system(cmd)
         #real_v_uc_plot_name="eda_model_LST_%s_%s_pol_real_vs_uc_%s_GHz.png" % (lst,pol,plot_freq_GHz)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uc,real options=nobase select=frequency\(%s\)' % (real_v_uc_plot_name,output_concat_model_vis_name,plot_freq_GHz)
         #print(cmd)
         #os.system(cmd)
   
         ##global signal
         #uv_dist_plot_name="eda_model_global_LST_%s_%s_pol_amp_vs_uvdist.png" % (lst,pol)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uvdist,amp options=nobase select=-auto' % (uv_dist_plot_name,global_output_concat_model_vis_name)
         #print(cmd)
         #os.system(cmd)
         #uv_dist_plot_name="eda_model_global_LST_%s_%s_pol_amp_vs_uc_%s_GHz.png" % (lst,pol,plot_freq_GHz)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uc,amp options=nobase select=frequency\(%s\)' % (uv_dist_plot_name,global_output_concat_model_vis_name,plot_freq_GHz)
         #print(cmd)
         #os.system(cmd)
         #amp_freq_plot_name="eda_model_global_LST_%s_%s_pol_amp_vs_freq.png" % (lst,pol)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=freq,amp options=nobase ' % (amp_freq_plot_name,global_output_concat_model_vis_name)
         #print(cmd)
         #os.system(cmd)
         #real_v_uc_plot_name="eda_model_global_LST_%s_%s_pol_real_vs_uc_%s_GHz.png" % (lst,pol,plot_freq_GHz)
         #cmd = 'uvplt device="%s/png" vis=%s  axis=uc,real options=nobase select=frequency\(%s\)' % (real_v_uc_plot_name,global_output_concat_model_vis_name,plot_freq_GHz)
         #print(cmd)
         #os.system(cmd)    
   
         #save data to a file with uvlist instead of plotting 
   
         #Going to need to dive into the uvfits files and plot stuff myself without miriad
         #see https://mail.python.org/pipermail/astropy/2013-December/002681.html
   
   #Can do this concat now!
   #DON'T Concat into a 2-hour chunk (this fails anyway because of pyuvdata and different antenna positions for different lsts for some reason)
   #Data will likely come in a number of 'snapshot' uvfits files so might as well just digest these individually in extract signal etc
   
   #get the average diffuse over all lsts ... NO, just do for one lst at a time.
   #sky_averaged_diffuse_array_no_beam_lsts = sky_averaged_diffuse_array_no_beam_lsts/n_lsts
   np.save(sky_averaged_diffuse_array_no_beam_lsts_filename,sky_averaged_diffuse_array_no_beam_lsts)
   for pol in pol_list:
      if pol == 'X':
         #sky_averaged_diffuse_array_beam_X_lsts = sky_averaged_diffuse_array_beam_X_lsts/n_lsts
         np.save(sky_averaged_diffuse_array_beam_X_lsts_filename,sky_averaged_diffuse_array_beam_X_lsts)
      else:
         #sky_averaged_diffuse_array_beam_Y_lsts = sky_averaged_diffuse_array_beam_Y_lsts/n_lsts
         np.save(sky_averaged_diffuse_array_beam_Y_lsts_filename,sky_averaged_diffuse_array_beam_Y_lsts)
   
   #concat the lsts together  
   cmd = "rm -rf %s %s %s %s" % (output_concat_vis_pyuvdata_name_lsts_X,output_concat_vis_pyuvdata_name_lsts_Y,output_concat_vis_pyuvdata_name_lsts_X_sub,output_concat_vis_pyuvdata_name_lsts_Y_sub)
   print(cmd)
   os.system(cmd)
      
   for pol in pol_list:
      if pol=='X':
          concat_uvfits(model_vis_uvfits_list_X_lsts,output_concat_uvfits_pyuvdata_name_lsts_X)
          if do_image_and_subtr_from_simulated_data:
             concat_uvfits(model_vis_uvfits_list_X_lsts_sub,output_concat_uvfits_pyuvdata_name_lsts_X_sub)
      else:
          concat_uvfits(model_vis_uvfits_list_Y_lsts,output_concat_uvfits_pyuvdata_name_lsts_Y)
          if do_image_and_subtr_from_simulated_data:
             concat_uvfits(model_vis_uvfits_list_Y_lsts_sub,output_concat_uvfits_pyuvdata_name_lsts_Y_sub)
       
   #remove the intermediate uvfits and concat freq uvfits
   for lst in lst_list:
      lst_deg = (float(lst)/24.)*360.
      for pol in pol_list:
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            model_vis_name_base = "%s_LST_%03d_%s_%0.3f_MHz" % (array_label,lst_deg,pol,freq_MHz)
            if 'noise' in signal_type_list:
               model_vis_name_base += '_N'
            if 'diffuse' in signal_type_list:
               model_vis_name_base += '_D_%s' % sky_model
            if 'global' in signal_type_list:
               model_vis_name_base += '_G' 
            if 'global_EDGES' in signal_type_list:
               model_vis_name_base += '_ED' 
            if 'gain_errors' in signal_type_list:
               model_vis_name_base += '_GE'
            
            out_vis_name = model_vis_name_base + '.vis'
            out_vis_uvfits_name = model_vis_name_base + '.uvfits'
            out_vis_name_sub = model_vis_name_base + '_sub.vis'
            out_vis_uvfits_name_sub = model_vis_name_base + '_sub.uvfits'
            
            #Don't remove the individual uvfits files (may not need concat files in future)
            #cmd = "rm -rf %s %s %s %s" % (out_vis_name,out_vis_uvfits_name,out_vis_name_sub,out_vis_uvfits_name_sub)
            cmd = "rm -rf %s %s" % (out_vis_name,out_vis_name_sub)
            print(cmd)
            os.system(cmd)
            
         output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
         output_concat_vis_pyuvdata_name_sub = "%s_concat_freqs_sub.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name_sub = "%s_concat_freqs_sub.uvfits" % model_vis_name_base
            
         cmd = "rm -rf %s %s %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name,output_concat_vis_pyuvdata_name_sub,output_concat_uvfits_pyuvdata_name_sub)
         print(cmd)
         os.system(cmd)


            
            
def plot_signal(lst_list,freq_MHz_list,pol_list,signal_type_list,outbase_name,sky_model,array_ant_locations_filename,array_label):
   concat_output_name_base_X = "%s_X_%s" % (array_label,outbase_name)
   concat_output_name_base_Y = "%s_Y_%s" % (array_label,outbase_name)
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
      if 'noise' in signal_type_list:
         concat_output_name_base_X += '_N'
         concat_output_name_base_Y += '_N'
      if 'diffuse' in signal_type_list:
         concat_output_name_base_X += '_D_%s' % sky_model
         concat_output_name_base_Y += '_D_%s' % sky_model
      if 'diffuse_global' in signal_type_list:
          concat_output_name_base_X += '_DG'
          concat_output_name_base_Y += '_DG'
      if 'diffuse_angular' in signal_type_list:
          concat_output_name_base_X += '_DA'
          concat_output_name_base_Y += '_DA'
      if 'global' in signal_type_list:
         concat_output_name_base_X += '_G' 
         concat_output_name_base_Y += '_G'
      if 'gain_errors' in signal_type_list:
         concat_output_name_base_X += '_GE'
         concat_output_name_base_Y += '_GE'
               
      if do_image_and_subtr_from_simulated_data:
         uv_type_list = ['original','subtracted']
      else:
         uv_type_list = ['original']
      
      for uv_type in uv_type_list:      
         if pol=='X':  
            if uv_type=='original':                          
               model_vis_name_base = concat_output_name_base_X
            else:
               model_vis_name_base = concat_output_name_base_X + '_sub'
         else:
            if uv_type=='original':
               model_vis_name_base = concat_output_name_base_Y
            else:
               model_vis_name_base = concat_output_name_base_Y + '_sub'
         model_vis_name_base += "_thresh_%0.2f" % (zero_spacing_leakage_threshold)
   
         signal_array_short_baselines_filename = "%s_signal.npy" % (model_vis_name_base)
         signal_short_baselines = np.load(signal_array_short_baselines_filename)
         signal_array_short_baselines_Tb_filename = "%s_signal_Tb.npy" % (model_vis_name_base)
         signal_short_baselines_Tb = np.load(signal_array_short_baselines_Tb_filename)
         signal_array_all_baselines_filename = "%s_signal_all_baselines.npy" % (model_vis_name_base)
         signal_all_baselines = np.load(signal_array_all_baselines_filename)
         signal_array_all_baselines_Tb_filename = "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
         signal_all_baselines_Tb = np.load(signal_array_all_baselines_Tb_filename)
         signal_array_all_baselines_filename_abs = "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
         signal_all_baselines_abs = np.load(signal_array_all_baselines_filename_abs)
         signal_array_all_baselines_filename_abs_Tb = "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
         signal_all_baselines_abs_Tb = np.load(signal_array_all_baselines_filename_abs_Tb)
         signal_array_short_baselines_weighted_filename = "%s_signal_weighted.npy" % (model_vis_name_base)
         signal_short_baselines_weighted = np.load(signal_array_short_baselines_weighted_filename)
         signal_array_short_baselines_weighted_Tb_filename = "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
         signal_short_baselines_weighted_Tb = np.load(signal_array_short_baselines_weighted_Tb_filename)
         number_baselines_used_array_filename = "%s_number_baselines_used.npy" % (model_vis_name_base)
         number_baselines_used_array = np.load(number_baselines_used_array_filename)
         sum_of_weights_all_baselines_array_filename = "%s_sum_of_weights_all_baselines.npy" % (model_vis_name_base)
         sum_of_weights_all_baselines_array = np.load(sum_of_weights_all_baselines_array_filename)
         sum_of_weights_short_baselines_array_filename = "%s_sum_of_weights_short_baselines.npy" % (model_vis_name_base)
         sum_of_weights_short_baselines_array = np.load(sum_of_weights_short_baselines_array_filename)
         signal_short_baselines_log = np.log10(abs(signal_short_baselines[0,:]))
         signal_all_baselines_log = np.log10(abs(signal_all_baselines[0,:]))
         signal_all_baselines_abs_log = np.log10(signal_all_baselines_abs[0,:])
         signal_short_baselines_weighted_log = np.log10(abs(signal_short_baselines_weighted[0,:]))
      
      
         ##all baselines Tb real
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines_Tb[0,:])
         #map_title="all baselines real Tb vs freq %s pol" % (pol)
         #plt.ylabel("sum vis real Tb")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s%s_real_vis_vs_freq_all_baselines_Tb.png" % (outbase_name,model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)
         
         #short baselines real Tb
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines_Tb[0,:])
         map_title="short baselines real Tb vs freq %s pol" % (pol)
         plt.ylabel("sum vis real Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_real_vis_vs_freq_short_baselines_Tb.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #short baselines real weighted Tb
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines_weighted_Tb[0,:])
         map_title="short baselines weighted real Tb vs freq %s pol" % (pol)
         plt.ylabel("sum vis real Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_real_vis_vs_freq_short_baselines_weighted_Tb.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
      
         ##all baselines abs Tb
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines_abs_Tb[0,:])
         #map_title="all baselines abs Tb vs freq %s pol" % (pol)
         #plt.ylabel("sum vis real Tb")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_abs_vis_vs_freq_all_baselines_Tb.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)
   
         #short baselines real
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines[0,:])
         map_title="short baselines real vis vs freq %s pol" % (pol)
         plt.ylabel("sum vis real (Jy)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_real_vis_vs_freq_short_baselines.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #short baselines real log
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines_log)
         map_title="short baselines log abs real vis vs freq %s pol" % (pol)
         plt.ylabel("log sum abs real vis (Jy)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_log_abs_real_vis_vs_freq_short_baseline.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
      
         #short baselines real weighted
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines_weighted[0,:])
         map_title="short baselines weighted real vis vs freq %s pol" % (pol)
         plt.ylabel("sum vis real (Jy)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_real_vis_vs_freq_short_baselines_weighted.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
       
         #short baselines real weighted log
         plt.clf()
         plt.plot(freq_MHz_list,signal_short_baselines_weighted_log)
         map_title="short baselines weighted real vis vs freq %s pol" % (pol)
         plt.ylabel("log sum abs real vis (Jy)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_log_real_vis_vs_freq_short_baselines_weighted.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)    
   
         ##all baselines real
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines[0,:])
         #map_title="all baselines real vis vs freq %s pol" % (pol)
         #plt.ylabel("sum vis real Jy")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_real_vis_vs_freq_all_baselines.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)
         
         #all baselines real log
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines_log)
         #map_title="all baselines log abs real vis vs freq %s pol" % (pol)
         #plt.ylabel("log sum abs real vis Jy")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_log_abs_real_vis_vs_freq_all_baseline.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)

         #all baselines abs
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines_abs[0,:])
         #map_title="all baselines abs vis vs freq %s pol" % (pol)
         #plt.ylabel("sum vis real Jy")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_abs_vis_vs_freq_all_baselines.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)
         
         #all baselines abs log
         #plt.clf()
         #plt.plot(freq_MHz_list,signal_all_baselines_abs_log)
         #map_title="all baselines log abs vis vs freq %s pol" % (pol)
         #plt.ylabel("log sum abs vis Jy")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_log_abs_vis_vs_freq_all_baseline.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)
         
         #number of baselines
         plt.clf()
         plt.plot(freq_MHz_list,number_baselines_used_array[0,:])
         map_title="Number of short baselines used"
         plt.ylabel("Number of baselines")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_number_of_baselines_used.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #sum of weights all baselines
         #plt.clf()
         #plt.plot(freq_MHz_list,sum_of_weights_all_baselines_array[0,:])
         #map_title="Sum of u,v=0 weights all baselines"
         #plt.ylabel("Sum of weights")
         #plt.xlabel("freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s_sum_of_weights_all_baselines.png" % (model_vis_name_base)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)

         #sum of weights short baselines
         plt.clf()
         plt.plot(freq_MHz_list,sum_of_weights_short_baselines_array[0,:])
         map_title="Sum of u,v=0 weights short baselines"
         plt.ylabel("Sum of weights")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s_sum_of_weights_short_baselines.png" % (model_vis_name_base)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
    
def extract_signal_from_eda2_data(eda2_data_uvfits_name_list,outbase_name,array_label):   
   outbase_name += "_thresh_%0.2f" % (zero_spacing_leakage_threshold)
   #marcin put some data in: bmckinley@bighorns.pawsey.org.au:/raid/data/eda/eda2/2019_06_11_Sun_test
   if array_label == 'eda2_sub48':
      n_ants = 48
   elif (array_label == 'eda2'):
      n_ants = 256
   n_baselines = n_ants*(n_ants-1) / 2. 
   for uvfits_name_index,uvfits_name in enumerate(eda2_data_uvfits_name_list):                         
         uvfits_filename = "%s%s" % (eda2_data_dir,uvfits_name)
         #print uvfits_filename           
         hdulist = fits.open(uvfits_filename)
         #hdulist.info()
         info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
         #print info_string
         uvtable = hdulist[0].data
         uvtable_header = hdulist[0].header
         visibilities = uvtable['DATA']
         visibilities_shape = visibilities.shape
         #print "visibilities shape"
         #print visibilities_shape
         
         n_pol = visibilities_shape[4]
         n_freq = visibilities_shape[3]
         
         #get the UU and VV so we can check whether we are using short baselines
         UU_s_array = uvtable['UU']
         UU_m_array = UU_s_array * c   
         VV_s_array = uvtable['VV']
         VV_m_array = VV_s_array * c
         
         n_vis = visibilities.shape[0]
         n_timesteps = n_vis/n_baselines
         #print "n_timesteps %s " % n_timesteps
         timestep_array = np.arange(0,n_timesteps,1)
            
         print("file %s has %s visibilities, %s timesteps, %s pols and %s freq chans" % (uvfits_filename,n_vis,n_timesteps,n_pol,n_freq))

def plot_EDA2_cal_sols(EDA2_chan,EDA2_obs_time,fine_chan_index,phase_sol_filename,amp_sol_filename,n_ants=256):
   #read the cal log into python and extract data for plotting. use pandas?
   #amp 
   #df = pd.read_csv(amp_sol_filename,header=2,index_col=False,engine='python',error_bad_lines=True)
   #print(df.head)
   #or look at randalls magic code....
   #append uncommented lines to list and note times
   
   
    
   with open(amp_sol_filename) as f:
      #looks like readlines ignores commented # lines
      lines = f.readlines()
   
   timestep_count = 1
   freq_bin_count = 0
   
   current_cal_time = 0
   sol_values_list = []
   
   #go through once to work out how many freq_bins and timesteps 
   for line_index,line in enumerate(lines):
      line_string_list = line.split()
      if line_string_list[0] != '#':
         if len(line_string_list)==8:
            sol_values = line_string_list[2:]
            for value in sol_values:  
               sol_values_list.append(float(value))
            timestep_count += 1
            new_cal_time = int("%s%s%s" % (line_string_list[1].strip()[0:2],line_string_list[1].strip()[3:5],line_string_list[1].strip()[6:8]))
            if line_index==0:
               current_cal_time = new_cal_time
            if new_cal_time > current_cal_time:
               #print new_cal_time
               #print(line_string_list)   
               current_cal_time = new_cal_time
            else:
               freq_bin_count += 1
               current_cal_time = new_cal_time
               timestep_count = 1
         else:
            sol_values = line_string_list
            for value in sol_values:  
               sol_values_list.append(float(value))
      

   data_frame = np.full([freq_bin_count,timestep_count,n_ants],np.nan)
   #print data_frame.shape
   
   #now go through again to put the data in the right place in the data frame
   for freq_bin in np.arange(freq_bin_count):
      for timestep in np.arange(timestep_count): 
         start_index = int(freq_bin*timestep)
         data_frame[freq_bin,timestep,:] = np.asarray(sol_values_list[start_index:start_index+n_ants])
   
   #print data_frame
   #you beauty, now plot it!
   
   
   for ant in np.arange(n_ants):
      plt.clf()
      for timestep in np.arange(timestep_count):
         cal_sols = data_frame[:,timestep,ant]
         #print cal_sols
         plt.plot(np.arange(freq_bin_count),cal_sols,label='timestep %s' % timestep)
      map_title="Cal sols ants"
      plt.ylabel("Amp")
      plt.xlabel("freq_bin")
      plt.legend(loc=1)
      fig_name= "cal_plot_amp_%s_%s_fine_chan_%s_ant_%s.png" % (EDA2_chan,EDA2_obs_time,fine_chan_index,ant)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)  
   
   
   
   
def calibrate_eda2_data(EDA2_chan_list,obs_type='night',lst_list=[],pol_list=[],sky_model_im_name='',n_obs_concat_list=[],concat=False,wsclean=False,plot_cal=False,uv_cutoff=0,per_chan_cal=False):
   #specify uv_cutoff in wavelengths, convert to m for 'calibrate'
   pol = pol_list[0]
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):  
       if len(n_obs_concat_list) > 0:
          if len(EDA2_chan_list)==1:
             n_obs_concat = n_obs_concat_list[chan_num]
          else:
             n_obs_concat = n_obs_concat_list[EDA2_chan_index]
       else:
          n_obs_concat = 1
       #freq_MHz = np.round(400./512.*float(EDA2_chan))
       freq_MHz = 400./512.*float(EDA2_chan)
       wavelength = 300./freq_MHz
       if uv_cutoff!=0:
          uv_cutoff_m = uv_cutoff * wavelength
       lst = lst_list[EDA2_chan_index]
       lst_deg = (float(lst)/24)*360.
       if len(EDA2_chan_list)==1:
          obs_time_list = EDA2_obs_time_list_each_chan[chan_num]
       else:
          obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
       first_obstime = obs_time_list[0]
       
       #gaurd against cases where there are no data for that channel
       if first_obstime==0:
          continue
       else:
          pass
      
       concat_vis_name = "%s/concat_chan_%s_%s_n_obs_%s.vis" % (EDA2_chan,EDA2_chan,first_obstime,n_obs_concat)
       concat_uvfits_name = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,first_obstime,n_obs_concat)

       concat_ms_name_wsclean_cal = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.ms" % (EDA2_chan,EDA2_chan,first_obstime,n_obs_concat)
       concat_uvfits_name_wsclean_cal = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,first_obstime,n_obs_concat)
              
       miriad_cal_vis_name_list = []
       miriad_cal_uvfits_name_list = []
       wsclean_cal_ms_name_list = []
       
       
       
       for EDA2_obs_time in obs_time_list:    
            #EDA2_obs_time = EDA2_obs_time_list[EDA2_chan_index]
            simulation_uvfits_name = "%s/eda_model_LST_%03d_%s_%0.3f_MHz_D_gsm.uvfits" % (EDA2_chan,lst_deg,pol,freq_MHz)
            simulation_ms_name = "%s/eda_model_LST_%03d_%s_%0.3f_MHz_D_gsm.ms" % (EDA2_chan,lst_deg,pol,freq_MHz)
            
            uvfits_filename = "%s/chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            uvfits_vis_filename = "%s/chan_%s_%s.vis" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            unity_sky_uvfits_filename = "%s/unity_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            unity_sky_vis_filename = "%s/unity_chan_%s_%s.vis" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            zero_sky_vis_filename = "%s/unity_chan_%s_%s.vis" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            #miriad selfcal
            calibrated_uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            #wsclean predict calibrate
            calibrated_uvfits_filename_wsclean = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            ms_name = "%s/%s_%s_eda2_ch32_ant256_midday_avg8140.ms" % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
            miriad_vis_name = uvfits_filename.split('.')[0] + '.vis'
            

            if not wsclean:
               cmd = "rm -rf %s %s" % (miriad_vis_name,calibrated_uvfits_filename)
               print(cmd)
               os.system(cmd)
               
               cmd = "fits in=%s out=%s op=uvin" % (uvfits_filename,miriad_vis_name)
               print(cmd)
               os.system(cmd)
            
               ##flag three edge chans for each coarse band
               #cmd = "uvflag vis=%s edge=3,3 flagval=flag options=brief,noquery" % (miriad_vis_name)
               #print(cmd)
               #os.system(cmd)
               
               ##check data, see how many chans there are in uvfits
               #with fits.open(uvfits_filename) as hdulist:
               #    uvtable = hdulist[0].data
               #    data = uvtable['DATA']
               #    shape=data.shape
               #    print shape
               #    visibilities_real_1chan_x = data[:,0,0,0,0,0]
               #    visibilities_imag_1chan_x = data[:,0,0,0,0,1]
               #    weights_1chan_x = data[:,0,0,0,0,2]
               #    
               #    
               #    end_print_index=10
               #    print visibilities_real_1chan_x[0:end_print_index]
               #    print visibilities_imag_1chan_x[0:end_print_index]
               #    print weights_1chan_x[0:end_print_index]

                
                #plt.clf()
                #n, bins, patches = plt.hist(visibilities_real_1chan_x)
                #map_title="Histogram of real vis" 
                #fig_name= "hist_real_vis.png"
                #figmap = plt.gcf()
                #figmap.savefig(fig_name)
                #plt.close()
                #print("saved %s" % fig_name)  
   
            
            if obs_type=='sun':
               
               cmd = "mfcal vis=%s flux=%s,0.150,0 refant=2 interval=1 stokes=xx,yy,xy,yx" % (miriad_vis_name,sun_flux_density)
               print(cmd)
               os.system(cmd)
               
               gain_solutions_name_amp = 'sun_%s_amp.txt' % (uvfits_filename.split('.')[0])
               gain_solutions_name_phase = 'sun_%s_ph.txt' % (uvfits_filename.split('.')[0])
               
               cmd = "rm -rf %s %s" % (gain_solutions_name_amp,gain_solutions_name_phase)
               print(cmd)
               os.system(cmd)
               
               cmd = "gpplt vis=%s device=output.png/png yaxis=phase options=bandpass nxy=7,7 log=%s" % (miriad_vis_name,gain_solutions_name_amp)
               print(cmd)
               os.system(cmd)
   
               gain_solutions_name_amp = 'sun_%s_amp.txt' % (uvfits_filename.split('.')[0])
               gain_solutions_name_phase = 'sun_%s_ph.txt' % (uvfits_filename.split('.')[0])
               
               cmd = "selfcal vis=%s flux=%s options=amplitude,noscale" % (miriad_vis_name,sun_flux_density)
               print(cmd)
               os.system(cmd)
               #### write out the solutions
               ##gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_amp.txt
               ##gpplt vis=$gainvis log=/tmp/GAINS/gains_${freq}MHz_lst${lst}_pha.txt yaxis=phase
               
               
               cmd = "rm -rf %s %s" % (gain_solutions_name_amp,gain_solutions_name_phase)
               print(cmd)
               os.system(cmd)
               
               cmd = "gpplt vis=%s log=%s" % (miriad_vis_name,gain_solutions_name_amp)
               print(cmd)
               os.system(cmd)
               
               cmd = "gpplt vis=%s log=%s yaxis=phase" % (miriad_vis_name,gain_solutions_name_phase)
               print(cmd)
               os.system(cmd)
               
               
            elif(obs_type=='night'):
               print("calibrating using sky model")
               #generate_apparent_sky_model(pol=pol,lst_hrs=lst,freq_MHz=freq)
               apparent_sky_im_name = "%s/apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (EDA2_chan,lst_deg,pol,freq_MHz)
               apparent_sky_fits_name = "%s/apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.fits" % (EDA2_chan,lst_deg,pol,freq_MHz)
               apparent_sky_model_fits_prefix = "%s/apparent_sky_LST_%03d_%s_pol_%0.3f_MHz" % (EDA2_chan,lst_deg,pol,freq_MHz)
               apparent_sky_casa_imname = "%s/apparent_sky_LST_%03d_%s_pol_%0.3f_MHz.model" % (EDA2_chan,lst_deg,pol,freq_MHz)
               
               uncal_ms_image_prefix = "uncal_chan_%s_%s_ms" % (EDA2_chan,EDA2_obs_time)
               uncal_ms_image_name = "%s-image.fits" % uncal_ms_image_prefix
               
               
               gsm_hpx_fits_name = "%s/%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (EDA2_chan,sky_model,lst_deg,freq_MHz)
               reprojected_to_wsclean_gsm_prefix = "%s/%s_map_LST_%03d_%0.3f_MHz_hpx_reprojected_wsclean" % (EDA2_chan,sky_model,lst_deg,freq_MHz)
               reprojected_to_wsclean_gsm_fitsname = "%s.fits" % (reprojected_to_wsclean_gsm_prefix)
               reprojected_to_wsclean_gsm_fitsname_Jy_per_pix = "%s_Jy_per_pix.fits" % (reprojected_to_wsclean_gsm_prefix)
               reprojected_to_wsclean_gsm_im_name_Jy_per_pix = "%s_map_LST_%03d_%0.3f_MHz_hpx_reprojected_wsclean_Jy_per_pix.im" % (sky_model,lst_deg,freq_MHz)

               #print apparent_sky_im_name
               
               if wsclean==True:
                  #reset pol back to the value input originally (gets changed in loop below!)
                  pol = pol_list[0]
                  #Even though using ms for wsclean, need to still use the uvfits file to generate unity sky visibilities
                  ############
                  #try uvgen on the uncal uvfits 
                  apparent_unity_sky_im_name = "%s/apparent_unity_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (EDA2_chan,lst_deg,pol,freq_MHz)
                  apparent_unity_sky_im_name_copy = "apparent_unity_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
                  #uv_dist_plot_name = "test_uvdist.png"
                  
                  cmd = "cp -r %s %s" % (apparent_unity_sky_im_name,apparent_unity_sky_im_name_copy)
                  print(cmd)
                  os.system(cmd)
                  
                  
                  apparent_zero_sky_im_name = "apparent_zero_sky_LST_%03d_%s_pol_%0.3f_MHz.im" % (lst_deg,pol,freq_MHz)
                  
                  cmd = "rm -rf %s %s %s %s %s" % (unity_sky_uvfits_filename,unity_sky_vis_filename,zero_sky_vis_filename,uvfits_vis_filename,apparent_zero_sky_im_name)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "maths exp=%s*0.0 out=%s " % (apparent_unity_sky_im_name_copy,apparent_zero_sky_im_name)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "fits in=%s op=uvin out=%s" % (uvfits_filename,uvfits_vis_filename)
                  print(cmd)
                  os.system(cmd)

                  cmd = "uvmodel vis=%s model=%s options=replace,mfs out=%s" % (uvfits_vis_filename,apparent_zero_sky_im_name,zero_sky_vis_filename)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "uvmodel vis=%s model=%s options=replace,mfs out=%s" % (zero_sky_vis_filename,apparent_zero_sky_im_name,unity_sky_vis_filename)
                  print(cmd)
                  os.system(cmd)
                  

                  
                  #cmd = 'uvplt device="%s/png" vis=%s  axis=uvdist,amp options=nobase select=-auto' % (uv_dist_plot_name,out_vis)
                  #print(cmd)
                  #os.system(cmd)  
         
                  cmd = "fits in=%s op=uvout options=nocal,nopol,nopass out=%s" % (unity_sky_vis_filename,unity_sky_uvfits_filename)
                  print(cmd)
                  os.system(cmd)
                  
                  wsclean_imsize = '512'
                  wsclean_scale = '900asec'
                  
                  print("cal using wsclean predict and calibrate / CASA bandpass")
                  
                  #Cant import the EDA sim uvfits file into ms - too many antennas - what if you make an mwa 32T simulated uvfits file instead!? then import to ms and image at low res (do this in simulate())
                  
                  #make a wsclean image of the uncalibrated ms just to get an image header to reproject to:
                  cmd = "wsclean -name %s -size %s %s -scale %s -pol xx -data-column DATA %s " % (uncal_ms_image_prefix,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                  print(cmd)
                  os.system(cmd)
                  
                  
                  
                  ###
                  if os.path.isfile(uncal_ms_image_name) and os.access(uncal_ms_image_name, os.R_OK):
                     hdulist = pyfits.open(uncal_ms_image_name)
                  else:
                     print("Either file %s is missing or is not readable" % uncal_ms_image_name)
                     #continue        
                  
                  data=hdulist[0].data[0,0,:,:]
                  new_header=hdulist[0].header
                  
                  pix_size_deg = float(new_header['CDELT1'])
                  pix_area_deg_sq = pix_size_deg*pix_size_deg
                  pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
                  
                  #needs cooordsys keyword
                  #pt_source_header['COORDSYS'] = 'icrs'
                  
                  #print(pt_source_header)
                  
                  del new_header[8]
                  del new_header[8]
                  del new_header['history']
                  new_header['bunit'] ='Jy/Pixel'                                  
                  #print(pt_source_header)
                  
                  target_wcs = WCS(new_header)
            
                  target_wcs=target_wcs.dropaxis(2)
                  target_wcs=target_wcs.dropaxis(2)
                                  
                  #hdu_hpx = pyfits.open(gsm_hpx_fits_name)[1]
                  ##hdu_hpx.info()
                  #hpx_header = hdu_hpx.header
                  #print(hpx_header)
                  
                  ##########################
                  ##########################
                  #repeat the below for each fine chan to get each correctly-name apparent sky
                  ##########################
                  ##########################
                  
                  
                  hdu_gsm = fits.open(gsm_hpx_fits_name)[1]
                  print(hdu_gsm.header)
                  #see bottom for accessing table data https://python4astronomers.github.io/astropy/fits.html
                  
                  reprojected_gsm_map,footprint = reproject_from_healpix(hdu_gsm, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
                  
                  #write the reprojected gsm maps to fits
                  pyfits.writeto(reprojected_to_wsclean_gsm_fitsname,reprojected_gsm_map,clobber=True)
                  #print new_header
                  pyfits.update(reprojected_to_wsclean_gsm_fitsname,reprojected_gsm_map,header=new_header)
                  print("wrote image %s" %  reprojected_to_wsclean_gsm_fitsname)
 
                  #model needs to be in Jy/pix
                  #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
         
                  scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
                  print("scale map by %s to get to Jy/pix" % scale)
                  
                  ###
                  #Dont need this now as can use reprojected to wsclean image
                  ##rename apparent sky fits image
                  #cmd = "cp %s %s-model.fits" % (apparent_sky_fits_name,apparent_sky_model_fits_prefix)
                  #print(cmd)
                  #os.system(cmd)                  
                  #
                  
                  #check the model image for non-finite values and scale to Jy per pix:
                  with fits.open("%s" % (reprojected_to_wsclean_gsm_fitsname)) as hdu_list:
                     data = hdu_list[0].data
                  #replace nans with zeros
                  data_new = np.nan_to_num(data)
                  data_new_jy_per_pix = data_new * scale
                  
                  #write out a new fits file
                  fits.writeto("%s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix),data_new_jy_per_pix,clobber=True)
                  pyfits.update(reprojected_to_wsclean_gsm_fitsname_Jy_per_pix,data_new_jy_per_pix,header=new_header)
                  print("saved %s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix))
                  
                  for pol in ['X','Y']:
                     if use_analytic_beam:
                        if pol=='X':
                           beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                           beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                        else:
                           beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                           beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
                     else:
                        beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
     
     
                     cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
                     print(cmd)
                     os.system(cmd)
                     
                     #Need to regrid the beam to the reproject gsm
                     beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
                     beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
                     beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.im' % (pol,freq_MHz)
                     beam_image_sin_projected_regrid_gsm_fits_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.fits' % (pol,freq_MHz)
                     
                     cmd = "rm -rf %s %s %s %s %s" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name,beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
                     print(cmd)
                     os.system(cmd)
                     
                     
                     #put in the correct ra in the header (ra = lst for zenith) 
                     #puthd in="$beam/crval1" value=$lst_degs
                     cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
                     print(cmd)
                     os.system(cmd) 
                     
                     #write out as a fits file to check header
                     cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "fits in=%s out=%s op=xyin" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix,reprojected_to_wsclean_gsm_im_name_Jy_per_pix)
                     print(cmd)
                     os.system(cmd)
                     
                     #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
                     cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix)
                     print(cmd)
                     os.system(cmd)  
                     
                     #Now have a gsm and a beam. multiply 'em'
                     apparent_sky_fits_name_prefix = "apparent_sky_LST_%03d_%0.3f_MHz_wsclean" % (lst_deg,freq_MHz)
                     apparent_sky_im_name = "%s.im" % (apparent_sky_fits_name_prefix)
                     apparent_sky_fits_name = "apparent_sky_LST_%03d_%0.3f_MHz_wsclean-%s%s-model.fits" % (lst_deg,freq_MHz,pol,pol)
                     
                     cmd = "rm -rf %s %s " % (apparent_sky_im_name,apparent_sky_fits_name)
                     print(cmd)
                     os.system(cmd)

                     cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix,apparent_sky_im_name)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_name,apparent_sky_fits_name)
                     print(cmd)
                     os.system(cmd) 
            
            
                     #check the model image for non-finite values 
                     with fits.open("%s" % (apparent_sky_fits_name)) as hdu_list:
                        data = hdu_list[0].data
                     #replace nans with zeros
                     data_new = np.nan_to_num(data)
                     
                     #write out a new fits file
                     fits.writeto("%s" % (apparent_sky_fits_name),data_new,clobber=True)
                     pyfits.update(apparent_sky_fits_name,data_new,header=new_header)
                     print("saved %s" % (apparent_sky_fits_name))
                  
                  
                     #check the model image for non-finite values:
                     #with fits.open("%s-model.fits" % (apparent_sky_model_fits_prefix)) as hdu_list:
                     #   data = hdu_list[0].data
                     #print(np.max(data))  
                     #print(np.min(data))
                     
                  ##########################
                  ##########################
                  for fine_chan_index in range(0,32):
                     centre_freq = float(freq_MHz)
                     fine_chan_width_MHz = fine_chan_width_Hz/1000000.   
                     
                     #dont reverse chan order
                     freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
                     #freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz
                     
                     wavelength_fine_chan = 300./float(freq_MHz_fine_chan)
                     
                     gsm_hpx_fits_name_fine_chan = "%s/%s_map_LST_%03d_%0.3f_MHz_hpx.fits" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan)
               
                     reprojected_to_wsclean_gsm_prefix_fine_chan = "%s/%s_map_LST_%03d_%0.3f_MHz_hpx_reprojected_wsclean" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan)
                     reprojected_to_wsclean_gsm_fitsname_fine_chan = "%s.fits" % (reprojected_to_wsclean_gsm_prefix_fine_chan)
                     reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan = "%s_Jy_per_pix.fits" % (reprojected_to_wsclean_gsm_prefix_fine_chan)
                     reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan = "%s_map_LST_%03d_%0.3f_MHz_hpx_reprojected_wsclean_Jy_per_pix.im" % (sky_model,lst_deg,freq_MHz_fine_chan)


                     hdu_gsm_fine_chan = fits.open(gsm_hpx_fits_name_fine_chan)[1]
                     #print(hdu_gsm_fine_chan.header)
                     #see bottom for accessing table data https://python4astronomers.github.io/astropy/fits.html
                     
                     reprojected_gsm_map_fine_chan,footprint_fine_chan = reproject_from_healpix(hdu_gsm_fine_chan, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
                     
                     #write the reprojected gsm maps to fits
                     pyfits.writeto(reprojected_to_wsclean_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,clobber=True)
                     #print new_header
                     pyfits.update(reprojected_to_wsclean_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,header=new_header)
                     print("wrote image %s" %  reprojected_to_wsclean_gsm_fitsname_fine_chan)
    
                     #model needs to be in Jy/pix
                     #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
            
                     scale_fine_chan = (2. * k * 1.0e26 * pix_area_sr) / (wavelength_fine_chan**2)
                     print("scale map by %s to get to Jy/pix" % scale_fine_chan)
                     
                     ###
                     #Dont need this now as can use reprojected to wsclean image
                     ##rename apparent sky fits image
                     #cmd = "cp %s %s-model.fits" % (apparent_sky_fits_name,apparent_sky_model_fits_prefix)
                     #print(cmd)
                     #os.system(cmd)                  
                     #
                     
                     #check the model image for non-finite values and scale to Jy per pix:
                     with fits.open("%s" % (reprojected_to_wsclean_gsm_fitsname_fine_chan)) as hdu_list_fine_chan:
                        data_fine_chan = hdu_list_fine_chan[0].data
                     #replace nans with zeros
                     data_new_fine_chan = np.nan_to_num(data_fine_chan)
                     data_new_jy_per_pix_fine_chan = data_new_fine_chan * scale_fine_chan
                     
                     #write out a new fits file
                     fits.writeto("%s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan),data_new_jy_per_pix_fine_chan,clobber=True)
                     pyfits.update(reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan,data_new_jy_per_pix_fine_chan,header=new_header)
                     print("saved %s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan))
                     
                     for pol in ['X','Y']:
                        if use_analytic_beam:
                           if pol=='X':
                              beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                              beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                           else:
                              beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                              beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
                        else:
                           beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
        
        
                        cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
                        print(cmd)
                        os.system(cmd)
                        
                        #Need to regrid the beam to the reproject gsm
                        beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
                        beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
                        beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.im' % (pol,freq_MHz)
                        beam_image_sin_projected_regrid_gsm_fits_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.fits' % (pol,freq_MHz)
                        
                        cmd = "rm -rf %s %s %s %s %s" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name,beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan)
                        print(cmd)
                        os.system(cmd)
                        
                        cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
                        print(cmd)
                        os.system(cmd)
                        
                        
                        #put in the correct ra in the header (ra = lst for zenith) 
                        #puthd in="$beam/crval1" value=$lst_degs
                        cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
                        print(cmd)
                        os.system(cmd) 
                        
                        #write out as a fits file to check header
                        cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
                        print(cmd)
                        os.system(cmd)
                        
                        cmd = "fits in=%s out=%s op=xyin" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan)
                        print(cmd)
                        os.system(cmd)
                        
                        
                        #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
                        cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan)
                        print(cmd)
                        os.system(cmd)  
                        
                        
                        
                        #Now have a gsm and a beam. multiply 'em'
                        apparent_sky_fits_name_prefix_fine_chan = "apparent_sky_LST_%03d_%0.3f_MHz_wsclean" % (lst_deg,freq_MHz)
                        apparent_sky_im_name_fine_chan = "apparent_sky_LST_%03d_%0.3f_MHz_wsclean-%04d.im" % (lst_deg,freq_MHz,fine_chan_index)
                        apparent_sky_fits_name_fine_chan = "%s-%04d-%s%s-model.fits" % (apparent_sky_fits_name_prefix_fine_chan,fine_chan_index,pol,pol)
                        
                        cmd = "rm -rf %s %s" % (apparent_sky_im_name_fine_chan,apparent_sky_fits_name_fine_chan)
                        print(cmd)
                        os.system(cmd)
   
                        cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan,apparent_sky_im_name_fine_chan)
                        print(cmd)
                        os.system(cmd)
                        
                        cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_name_fine_chan,apparent_sky_fits_name_fine_chan)
                        print(cmd)
                        os.system(cmd) 
               
                        print("wrote %s" % apparent_sky_fits_name_fine_chan)
                        
                        
                        #check the model image for non-finite values 
                        with fits.open("%s" % (apparent_sky_fits_name_fine_chan)) as hdu_list:
                           data = hdu_list[0].data
                        #replace nans with zeros
                        data_new = np.nan_to_num(data)
                        
                        #write out a new fits file
                        fits.writeto("%s" % (apparent_sky_fits_name_fine_chan),data_new,clobber=True)
                        pyfits.update(apparent_sky_fits_name_fine_chan,data_new,header=new_header)
                        print("saved %s" % (apparent_sky_fits_name_fine_chan))
                        
                        
                        
                  ###########################
                  ###########################
                  
                  
                  
                  ## predict a model onto the ms for calibration
                  ##cmd = "wsclean -predict -name %s -size 512 512 -scale 1800asec -pol xx %s " % (apparent_sky_fits_name, EDA2_chan,EDA2_obs_time,ms_name)
                  #cmd = "wsclean -predict -name %s -size %s %s -scale %s -pol xx,yy %s " % (apparent_sky_fits_name_prefix,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                  #print(cmd)
                  #os.system(cmd)                  
                  
                  ###make an image to check 
                  #cmd = "wsclean -name model_col_chan_%s_%s_ms -size %s %s -scale %s -pol xx -data-column MODEL_DATA %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                  #print(cmd)
                  #os.system(cmd)
                  
                  # predict a multi-channel model
                  cmd = "wsclean -predict -name %s -size %s %s -scale %s -pol xx,yy -channels-out 32 %s " % (apparent_sky_fits_name_prefix_fine_chan,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                  print(cmd)
                  os.system(cmd)
                  
                  ##make  images to check 
                  #cmd = "wsclean -name model_col_chan_%s_%s_ms -size %s %s -scale %s -pol xx -data-column MODEL_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                  #print(cmd)
                  #os.system(cmd)
                  
                  #crikey I think it actually works
                  
                  #hmmm seemed to actually work! We'll see ...
                  if uv_cutoff==0:
                     gain_solutions_name = 'cal_%s_%s_calibrate_sols.bin' % (EDA2_chan,EDA2_obs_time)
                     calibrate_options = ''
                  else:
                     gain_solutions_name = 'cal_%s_%s_calibrate_sols_uvcutoff_%0.3f_m.bin' % (EDA2_chan,EDA2_obs_time,uv_cutoff_m)
                     calibrate_options = '-minuv %0.3f ' % uv_cutoff_m
                  
                  cmd = "rm -rf %s" % (gain_solutions_name)
                  print(cmd)
                  os.system(cmd)
                  
                  #calibrate
                  
                  #calibrate on all 32 chans to increase SNR (poor results if you don't do this)
                  cmd = "calibrate  -ch 32 %s %s %s " % (calibrate_options,ms_name,gain_solutions_name)
                  print(cmd)
                  os.system(cmd)
                  
                  #cmd = "calibrate  %s %s %s " % (calibrate_options,ms_name,gain_solutions_name)
                  #print(cmd)
                  #os.system(cmd)
                  
                  
                  
                  #plot the sols and 
                  if (os.path.isfile(gain_solutions_name)):
                     if plot_cal:
                        #Plot the cal solutions
                        cmd = "aocal_plot.py  %s " % (gain_solutions_name)
                        print(cmd)
                        os.system(cmd)
                        
                     #write the calibrated uvfits file out ?
                     #cmd = "fits in=%s out=%s op=uvout" % (miriad_vis_name,calibrated_uvfits_filename)
                     #print(cmd)
                     #os.system(cmd)
                     
                     wsclean_cal_ms_name_list.append(ms_name) 

                     cmd = "applysolutions %s %s " % (ms_name,gain_solutions_name)
                     print(cmd)
                     os.system(cmd)
                     
                     ###make an image to check (both pols) 32 chans, skip for now, takes ages
                     #cmd = "wsclean -name cal_chan_%s_%s_ms -size %s %s -auto-threshold 5 -scale %s -pol xx,yy -data-column CORRECTED_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
                     #print(cmd)
                     #os.system(cmd) 
                     
                     #write out the uvfits file
                     casa_cmd_filename = 'export_individual_uvfits.sh'
                     cmd = "rm -rf %s %s" % (calibrated_uvfits_filename_wsclean,casa_cmd_filename)
                     print(cmd)
                     os.system(cmd)
                           
                     cmd = "exportuvfits(vis='%s',fitsfile='%s',datacolumn='corrected',overwrite=True,writestation=False)" % (ms_name,calibrated_uvfits_filename_wsclean)
                     print(cmd)
                     os.system(cmd)
        
                     with open(casa_cmd_filename,'w') as f:
                        f.write(cmd)
                          
                     cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
                     print(cmd)
                     os.system(cmd)
                     
                     
                  else:
                     print("no cal solutions for %s" % (ms_name))
                     continue
                           
                  
                               #############
                  #experiment with not concatenating and whether miriad will accpt an exported uvfits file and write it back out
                  #out_uvfits = 'test.uvfits'
                  #test_vis = 'test.vis'
                  #casa_cmd_filename = "casa_command.py"
                  #
                  #cmd = "rm -rf %s " % (casa_cmd_filename)
                  #print(cmd)
                  #os.system(cmd)
                  #
                  #cmd = "exportuvfits(vis='%s',fitsfile='%s',datacolumn='corrected',overwrite=True,writestation=False,multisource=False)" % (ms_name,out_uvfits)
                  #print(cmd)
                  #os.system(cmd)
                  #
                  #with open(casa_cmd_filename,'w') as f:
                  #   f.write(cmd)
                  #
                  #cmd = "casa --nohead --nogui -c %s" % casa_cmd_filename
                  #print(cmd)
                  #os.system(cmd)
                  #
                  #cmd = "fits in=%s op=uvin out=%s" % (out_uvfits,test_vis)
                  #print(cmd)
                  #os.system(cmd)
                  #
                  #
                  #cmd = "prthd in=%s " % test_vis
                  #print(cmd)
                  #os.system(cmd)
                  #
                  #cmd = "uvlist vis=%s options=array,full" % test_vis
                  #print(cmd)
                  #os.system(cmd)
                  
                  ##something is wrong when u do uvgen - try exporting a s auvfits and re-importing
                  #cmd = "fits in=%s op=uvout out=tmp.uvfits" % (test_vis)
                  #print(cmd)
                  #os.system(cmd)
                  
                  
             
                  #image with wsclean to get a model image
                  
                  #wsclean -predict the model into the real ms
                  
                  #use calibrate to ... calibrate
                  
                  #cmd = "wsclean -predict %s " % (ms_name)
                  #print(cmd)
                  #os.system(cmd)
                  
                  #if solutions were found:
                     #plot cal solutions
                  
                     #apply solutions applysolutions()
                  
                     #write out to uvfits with casa : exportuvfits(data_column='CORRECTED_DATA)
                     
                     #casa_cal_vis_name_list.append(ms_name)
                  
                  ##import apparent sky image in Jy/pix
                  #casa_cmd_filename = "casa_cmd.py"
                  #casa_cmd = "importfits(fitsimage='%s',imagename='%s',overwrite=True)" % (apparent_sky_fits_name,apparent_sky_casa_imname)
                  #with open(casa_cmd_filename,'w') as f:
                  #   f.write(casa_cmd)
                  #print("wrote %s " % (casa_cmd_filename))
                  #cmd = "casa --nohead --nogui -c %s " % casa_cmd_filename
                  #print(cmd)
                  #os.system(cmd)
                  
                  #use setjy with a sky model setjy(vis=myinput,field=myfluxcal,spw=myspw,modimage=myfluxmod,standard='Perley-Butler 2010',scalebychan=True)
                  #this is extremely slow or doesnt work - going to have to use wsclean -predict
                  #casa_cmd = "setjy(vis='%s',modimage='%s',standard='manual',usescratch=True)" % (ms_name,apparent_sky_casa_imname)
                  #with open(casa_cmd_filename,'w') as f:
                  #   f.write(casa_cmd)
                  #print("wrote %s " % (casa_cmd_filename))
                  #cmd = "casa -c %s " % casa_cmd_filename
                  #print(cmd)
                  #os.system(cmd)
                  
                  
                  
                  
                  #calibrate using bandpass: bandpass(vis=myinput,caltable=myfluxbpasstable,field=myfluxcal,refant=myrefant,solint='inf',combine='scan',solnorm=False, # for calibrationbandtype='B',gaintable=[myfluxbpphasetable,myfluxbpamptable],gaincurve=True)
                  
                  #if solutions were found:
                     #plot cal solutions
                  
                     #apply solutions applycal()
                  
                     #write out to uvfits: exportuvfits(data_column='CORRECTED_DATA)
                     
                     #casa_cal_vis_name_list.append(ms_name)
                  
               else:
                  print("cal using miriad selfcal")
                  
                  ##GO back to just using selfcal on the whole 32 fine chans ... all this splitting etc is not practical, this should be
                  #okay when the eda2 is fixed and the huge phase ramps are gone.
                  
                  #how many chans
                  #cmd = "uvlist vis=%s " % (miriad_vis_name)
                  #print(cmd) 
                  #os.system(cmd)
                  
                  with fits.open(uvfits_filename) as hdulist:
                     uvtable = hdulist[0].data
                     visibilities = uvtable['DATA']
                  
                  n_fine_chans = visibilities.shape[3]
                  
                  centre_freq = float(freq_MHz)
                  
                  #need to split off each chan and cal individually since selfcal does not apply bandpass corrections
                  #Noooo this is a dumb idea, just go back to using the average cal solutions across the band
                  #miriad_cal_vis_name_list_fine_chans = []
                  #miriad_cal_uvfits_name_list_fine_chans = []
                  #concat_vis_freq_name = "chan_%s_%s_all_f_ch.vis" % (EDA2_chan,EDA2_obs_time)
                  #concat_uvfits_freq_name = "chan_%s_%s_all_f_ch.uvfits" % (EDA2_chan,EDA2_obs_time)
                  #
                  #for fine_chan_index in range(0,n_fine_chans):
                  #   centre_freq = float(freq_MHz)
                  #
                  #   fine_chan_width_MHz = fine_chan_width_Hz/1000000.   
                  #   fine_chan_index = int(fine_chan_index)
                  #   freq_MHz_fine_chan = freq_MHz + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
                  #  
                  #   start_freq_GHz = (freq_MHz_fine_chan - fine_chan_width_MHz/2.) / 1000.
                  #   end_freq_GHz = (freq_MHz_fine_chan + fine_chan_width_MHz/2.) / 1000.
                  #   
                  #   #uvsplit has its own naming convention for outputs (read help)
                  #   orig_split_vis_name = "eda2cal.%s" % int(freq_MHz)
                  #   
                  #   cmd = "rm -rf %s" % (orig_split_vis_name)
                  #   print(cmd) 
                  #   os.system(cmd)
                  #   
                  #   #cmd = "uvsplit vis=%s select=frequency\(%0.6f,%0.6f\)" % (miriad_vis_name,start_freq_GHz,end_freq_GHz)
                  #   cmd = "uvaver vis=%s line=channel,1,%s,1,1 out=%s " % (miriad_vis_name,int(fine_chan_index+1),orig_split_vis_name)
                  #   print(cmd) 
                  #   os.system(cmd)
                  #   
                  #   ##check how many chans in new vis
                  #   #cmd = "uvlist vis=%s" % (orig_split_vis_name)
                  #   cmd = "prthd in=%s" % (orig_split_vis_name)
                  #   print(cmd) 
                  #   os.system(cmd)
                     
                     
                     
                  #   #rename the vis produced
                  #   new_split_vis_name = "fine_ch_%s.vis" % (fine_chan_index)
                  #   new_split_uvfits_name = "fine_ch_%s.uvfits" % (fine_chan_index)
                  #   
                  #   cmd = "rm -rf %s %s" % (new_split_vis_name,new_split_uvfits_name)
                  #   print(cmd) 
                  #   os.system(cmd)

                  #   cmd = "mv  %s %s" % (orig_split_vis_name,new_split_vis_name)
                  #   print(cmd) 
                  #   os.system(cmd)
                     
                     
                      
                     #select=uvrange(uvmin,uvmax) in kilolambda
                     #also test with gpscal ... not sure what selfcal is doing with polarisation/stokes....
                  if uv_cutoff!=0:
                     cmd = "selfcal vis=%s interval=1 model=%s line=channel,1,1,1,1 options=amplitude,noscale,mfs select=uvrange\(%0.5f,50\)" % (miriad_vis_name,apparent_sky_im_name,uv_cutoff)
                     gain_solutions_name_amp = 'cal_%s_%s_%0.3f_MHz_amp_uvcut_%0.5f.txt' % (EDA2_chan,EDA2_obs_time,centre_freq,uv_cutoff)
                     gain_solutions_name_phase = 'cal_%s_%s_%0.3f_MHz_ph_uvcut_%0.5f.txt' % (EDA2_chan,EDA2_obs_time,centre_freq,uv_cutoff)
                  else:
                     cmd = "selfcal vis=%s interval=1 model=%s options=amplitude,noscale,mfs" % (miriad_vis_name,apparent_sky_im_name)
                     gain_solutions_name_amp = 'cal_%s_%s_%0.3f_MHz_amp.txt' % (EDA2_chan,EDA2_obs_time,centre_freq)
                     gain_solutions_name_phase = 'cal_%s_%s_%0.3f_MHz_ph.txt' % (EDA2_chan,EDA2_obs_time,centre_freq)
                  print(cmd) 
                  os.system(cmd)
                     
                  
                  cmd = "rm -rf %s %s" % (gain_solutions_name_amp,gain_solutions_name_phase)
                  print(cmd)
                  os.system(cmd)
                  
                  #this doesnt work - need to just run gpplt and then plot manually
                  #cmd = "gpplt vis=%s device=output.png/png yaxis=amplitude options=bandpass nxy=7,7 log=%s" % (miriad_vis_name,gain_solutions_name_amp)
                  cmd = "gpplt vis=%s yaxis=amplitude log=%s" % (miriad_vis_name,gain_solutions_name_amp)
                  print(cmd)
                  os.system(cmd)
      
                  cmd = "gpplt vis=%s yaxis=phase log=%s" % (miriad_vis_name,gain_solutions_name_phase)
                  print(cmd)
                  os.system(cmd)               
                  
                  
                  #plot the sols and make in image to check
                  if (os.path.isfile(gain_solutions_name_phase) and os.path.isfile(gain_solutions_name_amp)):
                     if plot_cal:
                        plot_EDA2_cal_sols(EDA2_chan,EDA2_obs_time,fine_chan_index,gain_solutions_name_phase,gain_solutions_name_amp)
                     #write the calibrated uvfits file out
                     #cmd = "fits in=%s out=%s op=uvout" % (miriad_vis_name,calibrated_uvfits_filename)
                     #print(cmd)
                     #os.system(cmd)
                     
                     miriad_cal_vis_name_list.append(miriad_vis_name) 
                     miriad_cal_uvfits_name_list.append(calibrated_uvfits_filename) 
                     
                     #write out the uvfits file (this also applies the calibration)
                     cmd = "fits in=%s out=%s op=uvout" % (miriad_vis_name,calibrated_uvfits_filename)
                     print(cmd)
                     os.system(cmd)
                     
                     image_basename="miriad_cal_chan_%s_%s" % (EDA2_chan,EDA2_obs_time)
      
                     cmd = "rm -rf %s_xx.map %s.beam %s_xx.model %s_xx.restor" % (image_basename,image_basename,image_basename,image_basename)
                     print(cmd)
                     os.system(cmd)
                     
                     cmd = "invert vis=%s map=%s_xx.map beam=%s.beam options=double,mfs imsize=512 stokes=xx robust=-0.5 cell=1800 " % (miriad_vis_name,image_basename,image_basename)
                     print(cmd)
                     os.system(cmd)
                     cmd = "clean map=%s_xx.map beam=%s.beam out=%s_xx.model niters=2000" % (image_basename,image_basename,image_basename)
                     print(cmd)
                     os.system(cmd)
                     cmd = "restor model=%s_xx.model  beam=%s.beam map=%s_xx.map out=%s_xx.restor " % (image_basename,image_basename,image_basename,image_basename)
                     print(cmd)
                     os.system(cmd)
                        
                     cmd = "fits in=%s_xx.map out=%s_xx.fits op=xyout" % (image_basename,image_basename)
                     print(cmd)
                     os.system(cmd)
         
                  else:
                     print("no cal solutions for %s" % (miriad_vis_name))
                     continue
                   
                  ##concat them in freq here
                  #print("concatenating calibrated data in freq with pyuvdata")
             
                  #cmd = "rm -rf %s %s" % (concat_vis_freq_name,concat_uvfits_freq_name)
                  #print(cmd)
                  #os.system(cmd)
                  
                  #concat_uvfits(miriad_cal_uvfits_name_list_fine_chans,concat_uvfits_freq_name)
                  ##vis_string = ','.join(miriad_cal_vis_name_list_fine_chans)
                 
                  #sys.exit()
                  
                  #cmd = "uvaver vis=%s out=%s options=nopol" % (vis_string,concat_vis_freq_name)
                  #print(cmd)
                  #os.system(cmd)
                  ##then
                  
                  ##check how many chans in new vis
                  ##cmd = "uvlist vis=%s" % (concat_vis_freq_name)
                  #cmd = "prthd in=%s" % (concat_vis_freq_name)
                  #print(cmd) 
                  #os.system(cmd)
                  
                  ##pot uv data to see whats going on
                  #cmd = "uvplt vis=%s axis=freq,amplitude device=uvplt_amp_v_freq.png/png" % (concat_vis_freq_name)
                  #print(cmd)
                  #os.system(cmd)

                  #miriad_cal_vis_name_list.append(concat_vis_freq_name)
                  #miriad_cal_uvfits_name_list.append(concat_uvfits_freq_name)
       
                  
       #print miriad_cal_vis_name_list
       
       if concat==True:
          if wsclean==True:
             print("concatenating using CASA")
             
             casa_cmd_filename = "casa_command_concat.py"
             
             cmd = "rm -rf %s %s %s" % (concat_ms_name_wsclean_cal,concat_uvfits_name_wsclean_cal,casa_cmd_filename)
             print(cmd)
             os.system(cmd)
             
             vis_string = "','".join(wsclean_cal_ms_name_list)
             
             #need to write these commmands to files and then run 'casa -c --nohead --nogui file.py'
             cmd = "concat(vis=['%s'],concatvis='%s',dirtol='3600arcsec')" % (vis_string,concat_ms_name_wsclean_cal)
             print(cmd)
             os.system(cmd)

             with open(casa_cmd_filename,'w') as f:
                f.write(cmd)
                  
             cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
             print(cmd)
             os.system(cmd)


             #calibrate per chan on concat ms
             if per_chan_cal==True:
                if uv_cutoff==0:
                   gain_solutions_name = 'cal_%s_%s_calibrate_sols_per_chan.bin' % (EDA2_chan,EDA2_obs_time)
                   calibrate_options = ''
                else:
                   gain_solutions_name = 'cal_%s_%s_calibrate_sols_uvcutoff_%0.3f_m_per_chan.bin' % (EDA2_chan,EDA2_obs_time,uv_cutoff_m)
                   calibrate_options = '-minuv %0.3f ' % uv_cutoff_m

                #calibrate on each chan and use the corrected data
                cmd = "calibrate  -ch 1 -datacolumn CORRECTED_DATA %s %s %s " % (calibrate_options,concat_ms_name_wsclean_cal,gain_solutions_name)
                print(cmd)
                os.system(cmd)
                
                #apply sols
                cmd = "applysolutions -datacolumn CORRECTED_DATA %s %s " % (concat_ms_name_wsclean_cal,gain_solutions_name)
                print(cmd)
                os.system(cmd)
                
             casa_cmd_filename = "casa_command_exportuvfits.py"   

             cmd = "rm -rf %s " % (casa_cmd_filename)
             print(cmd)
             os.system(cmd)
                   
             cmd = "exportuvfits(vis='%s',fitsfile='%s',datacolumn='corrected',overwrite=True,writestation=False)" % (concat_ms_name_wsclean_cal,concat_uvfits_name_wsclean_cal)
             print(cmd)
             os.system(cmd)

             with open(casa_cmd_filename,'w') as f:
                f.write(cmd)
                  
             cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
             print(cmd)
             os.system(cmd)
             
             #casa still causing /tmp to fill up with crash reports
             cmd = "rm -rf /tmp/*" 
             print(cmd)
             os.system(cmd)


             #image the concat file to check spectral behaviour
             ##make an image to check (both pols) 32 chans
             cmd = "wsclean -name concat_cal_chan_%s_ms -size %s %s -auto-threshold 5 -niter 10000 -scale %s -pol xx,yy -data-column CORRECTED_DATA -channels-out 32 %s " % (EDA2_chan,wsclean_imsize,wsclean_imsize,wsclean_scale,concat_ms_name_wsclean_cal)
             print(cmd)
             os.system(cmd)
                          
          else:
             print("concatenating calibrated data into one vis and uvfits with miriad")
             
             cmd = "rm -rf %s %s" % (concat_vis_name,concat_uvfits_name)
             print(cmd)
             os.system(cmd)
             
             vis_string = ','.join(miriad_cal_vis_name_list)
         
             cmd = "uvaver vis=%s out=%s options=nopol" % (vis_string,concat_vis_name)
             print(cmd)
             os.system(cmd)
                  
             cmd = "fits in=%s out=%s op=uvout" % (concat_vis_name,concat_uvfits_name)
             print(cmd)
             os.system(cmd)
               
             #image the concat vis to seee if the cal worked   

               
          
def image_eda2_data(eda2_data_uvfits_name_list):
   for uvfits_name_index,uvfits_name in enumerate(eda2_data_uvfits_name_list):                         
         uvfits_filename = "%s%s" % (eda2_data_dir,uvfits_name)
         miriad_vis_name = uvfits_name.split('.')[0] + '.vis'

         image_basename=uvfits_name.split('.')[0]
      
         cmd = "rm -rf %s_xx.map %s.beam %s_xx.model %s_xx.restor" % (image_basename,image_basename,image_basename,image_basename)
         print(cmd)
         os.system(cmd)
         
         cmd = "invert vis=%s map=%s_xx.map beam=%s.beam options=double imsize=512 stokes=xx robust=-0.5 cell=1800 " % (miriad_vis_name,image_basename,image_basename)
         print(cmd)
         os.system(cmd)
         cmd = "clean map=%s_xx.map beam=%s.beam out=%s_xx.model niters=2000" % (image_basename,image_basename,image_basename)
         print(cmd)
         os.system(cmd)
         cmd = "restor model=%s_xx.model  beam=%s.beam map=%s_xx.map out=%s_xx.restor " % (image_basename,image_basename,image_basename,image_basename)
         print(cmd)
         os.system(cmd)
            
         cmd = "fits in=%s_xx.map out=%s_xx.fits op=xyout" % (image_basename,image_basename)
         print(cmd)
         os.system(cmd)


def plot_baseline_length_counts(array_layout_filename,freq_MHz,lambda_threshold=0.5,m_threshold=3.0):
   
   plot_basename = array_layout_filename.split('/')[-1].split('.')[0]
   wavelength = 300./float(freq_MHz)
   
   antenna_position_x_list=[]
   antenna_position_y_list=[]
   
   baseline_length_list = []
   
   with open(array_layout_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[0])
      antenna_position_y = float(line.strip().split()[1])
      antenna_position_x_list.append(antenna_position_x)
      antenna_position_y_list.append(antenna_position_y)   
   
   antenna_position_x_m_array = np.asarray(antenna_position_x_list)
   antenna_position_y_m_array = np.asarray(antenna_position_y_list)

   n_ants = len(antenna_position_x_m_array)
   
   n_baselines_predicted = (n_ants * (n_ants-1))/2
   #print(n_baselines_predicted)
   
   for ant1 in range(n_ants/2):
      for ant2 in range(n_ants):
         if ant1!=ant2:
            x_length = antenna_position_x_m_array[ant1] - antenna_position_x_m_array[ant2]
            y_length = antenna_position_y_m_array[ant1] - antenna_position_y_m_array[ant2]
            baseline_length = np.sqrt(x_length**2 + y_length**2)
            baseline_length_list.append(baseline_length)
            
   #print(len(baseline_length_list))

   
   baseline_length_array = np.asarray(baseline_length_list)
   #histogram of baseline lengths
   #plot a histogram of Y values
   plt.clf()
   n, bins, patches = plt.hist(baseline_length_array)
   map_title="Histogram of baseline lengths" 
   fig_name= "hist_baseline_lenghts_%s.png" % (plot_basename)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name)  
   
   baseline_length_array_m_inds = baseline_length_array.argsort()
   baseline_length_array_sorted = baseline_length_array[baseline_length_array_m_inds]
   baseline_length_array_sorted_m_cut = baseline_length_array_sorted[baseline_length_array_sorted<m_threshold]
   
   baseline_length_array_sorted_lambda = baseline_length_array_sorted / wavelength
   baseline_length_array_sorted_lambda_cut = baseline_length_array_sorted_lambda[baseline_length_array_sorted_lambda < lambda_threshold]
   
   n_m_cut = len(baseline_length_array_sorted_m_cut)
   n_wavelength_cut = len(baseline_length_array_sorted_lambda_cut)
   
   print("n baselines theoretical is %s" % n_baselines_predicted)
   print("n baselines from list length is %s" % len(baseline_length_list))
   print("n baseline less than %s m is %s" % (m_threshold,n_m_cut))
   print("n baseline less than %s wavelengths is %s" % (lambda_threshold,n_wavelength_cut))
   
   
def plot_antenna_array(array_layout_filename,ylim=18):
   # txt file needs to be miriad format (E W U)
   antenna_layout_basename = array_layout_filename.split('/')[-1].split('.')[0]
   antenna_name_list = range(1,257)

   antenna_position_x_list=[]
   antenna_position_y_list=[]
   with open(array_layout_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[0])
      antenna_position_y = float(line.strip().split()[1])
      antenna_position_x_list.append(antenna_position_x)
      antenna_position_y_list.append(antenna_position_y)   
   
   antenna_position_x_m = np.asarray(antenna_position_x_list)
   antenna_position_y_m = np.asarray(antenna_position_y_list)
   
   #Plot antenna positions
   antenna_position_plot_figname = "%s_layout.png" % (antenna_layout_basename)
   antenna_position_plot_title = 'Antenna Positions AAVS-1 Tests'
   
   fig, ax = plt.subplots()
   ax.scatter(antenna_position_x_m,antenna_position_y_m, marker='+')
   
   #don't have the names
   #for i, name in enumerate(antenna_name_list):
   #   ax.annotate(str(name), (antenna_position_x_m[i],antenna_position_y_m[i]),size=5, color='grey')  

   plt.xlabel('X offset from centre (m) ')
   plt.ylabel('Y offset from centre (m) ')
   plt.xlim((-ylim,ylim))
   plt.ylim((-ylim,ylim))
   x0,x1 = ax.get_xlim()
   y0,y1 = ax.get_ylim()
   ax.set_aspect(abs(x1-x0)/abs(y1-y0))
   plt.title(antenna_position_plot_title)
   plt.savefig(antenna_position_plot_figname,dpi = 900)
   print("saved %s" % antenna_position_plot_figname)
   plt.close()

     
def grow_new_array(seed,n_ants,min_spacing_m):
   #want to pack them in as tight as possible, so work out diameter (4/pi comes from ration of area of a square to a circle)
   diameter = (4./math.pi) * min_spacing_m * math.sqrt(n_ants)
   print("diameter of array with %s ants is %0.2f m" % (int(n_ants),diameter))
   min_spacing_cm = int(min_spacing_m*100.)
   new_ant_array_locations_filename = 'sim_array_seed_%03d_s_%03d_a_%03d.txt' % (seed,min_spacing_cm,n_ants)
   x_offset_list = []
   y_offset_list = []
   radius=diameter/2.
   random.seed(seed)
   #print 'random seed %s' % (seed)
   antennas_placed = 0
   while (antennas_placed<n_ants):
      x_offset_array = np.asarray(x_offset_list)
      y_offset_array = np.asarray(y_offset_list)
      #pick a random x,y offset from centre
      x_offset = (random.random()-0.5) * diameter
      y_offset = (random.random()-0.5) * diameter
      dist_from_centre = math.sqrt(x_offset**2+y_offset**2)
      x_offsets_from_placed_ants = x_offset_array - x_offset
      y_offsets_from_placed_ants = y_offset_array - y_offset
      distances_to_placed_ants = np.sqrt(x_offsets_from_placed_ants**2 + y_offsets_from_placed_ants**2)
      if antennas_placed > 0:
         min_distance_between_ants = np.min(distances_to_placed_ants)
      else:
         min_distance_between_ants = diameter
      while(dist_from_centre>radius or min_distance_between_ants<min_spacing_m):
         x_offset = (random.random()-0.5) * diameter
         y_offset = (random.random()-0.5) * diameter
         dist_from_centre = math.sqrt(x_offset**2+y_offset**2)
         x_offsets_from_placed_ants = x_offset_array - x_offset
         y_offsets_from_placed_ants = y_offset_array - y_offset
         distances_to_placed_ants = np.sqrt(x_offsets_from_placed_ants**2 + y_offsets_from_placed_ants**2)
         if antennas_placed > 0:
            min_distance_between_ants = np.min(distances_to_placed_ants)
         else:
            min_distance_between_ants = diameter   
      x_offset_list.append(x_offset)
      y_offset_list.append(y_offset)
      antennas_placed += 1
      print('ant %s dist_from_centre %s, min distance to another ant %s' % (antennas_placed,dist_from_centre,min_distance_between_ants))
   with open(new_ant_array_locations_filename,'w') as f:
      for x_offset_index,x_offset in enumerate(x_offset_list):
         y_offset = y_offset_list[x_offset_index]
         line = "%0.3f   %0.3f   0\n" % (x_offset,y_offset)
         f.write(line)
   print("wrote %s" % new_ant_array_locations_filename)
   return new_ant_array_locations_filename
   #print x_offset_list
   #print y_offset_list
   
def model_average_from_arrays(n_arrays,min_spacing_m,zero_spacing_leakage_threshold,addition_type,n_ants):
   #addition_type can be 'simple' or 'coherent'
   #simple just averages each array signal that has been extracted seperately
   #coherent uses the outputs from cohently adding the weighted signal previously using extract_signal_from_sims_multi
   min_spacing_cm = int(min_spacing_m*100.)
   freq_MHz_array = np.asarray(freq_MHz_list)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   model_vis_name_base = 'test_%03d_arrays_%03d_s_%03d_a_%03d_pc_thresh_%s' % (n_arrays,min_spacing_cm,n_ants,zero_spacing_leakage_threshold_percent,addition_type)
   if addition_type == 'simple':
      #signal_filename_01 = 'array000_sep_%03d_%03d_ants_X_lst_2.00_hr_int_0.13_hr_N_D_gmoss_G_thresh_%0.2f_signal_weighted_Tb.npy' % (min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      #signal_filename_01 = "array000_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_N_D_gmoss_G_thresh_%0.2f_signal_weighted_Tb.npy" % (min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      signal_filename_01 = "array000_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy" % (min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      signal_array_01 = np.load(signal_filename_01)
      signal_array_average = signal_array_01 * 0.0
      for subarray in range(0,n_arrays):
         #signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_N_D_gmoss_G_thresh_%0.2f_signal_weighted_Tb.npy' % (subarray,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
         signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (subarray,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
         signal_array = np.load(signal_filename)
         signal_array_average += signal_array
      signal_array_average = signal_array_average/n_arrays
   elif addition_type == 'coherent':
      #signal_filename = 'array%03d_sep_%03d_%03d_ants_X_lst_2.00_hr_int_0.13_hr_N_D_gmoss_G_thresh_%0.2f_signal_weighted_Tb_multi.npy' % (n_arrays-1,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      #signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_N_D_gmoss_G_thresh_%0.2f_signal_weighted_Tb_multi.npy' % (n_arrays-1,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb_multi.npy' % (n_arrays-1,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      signal_array_average = np.load(signal_filename)
   
   #can't use the nan values:
   good_idx = np.isfinite(signal_array_average[0,:]) & np.isfinite(signal_array_average[0,:])
   good_signal_array_short_baselines_Tb = signal_array_average[0,:][good_idx]
   
   good_freq_MHz_array = freq_MHz_array[good_idx]
   
   coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
   ffit = poly.polyval(good_freq_MHz_array, coefs)
   
   #in log log space:
   log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
   log_freq_MHz_array = np.log10(good_freq_MHz_array)
   coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
   ffit = poly.polyval(log_freq_MHz_array, coefs)
   ffit_linear = 10**ffit
   
   #residual = log_signal_array_short_baselines - log_ffit
   residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
   
   rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
   
   plt.clf()
   plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("freq (MHz)")
   plt.legend(loc=1)
   fig_name= "%s_log_fit_residual_weighted_poly_%s.png" % (model_vis_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   return rms_of_residuals

def model_and_plot_random_array_layout_residuals(n_arrays,min_spacing_m,zero_spacing_leakage_threshold,n_ants):
   min_spacing_cm = int(min_spacing_m*100.)
   freq_MHz_array = np.asarray(freq_MHz_list)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   final_plot_name_base = 'test_%03d_layouts_%03d_s_%03d_a_%03d_pc_thresh' % (n_arrays,min_spacing_cm,n_ants,zero_spacing_leakage_threshold_percent)
   rms_of_residuals_list = []
   random_array_number_list = []
   for subarray in range(0,n_arrays):
      plot_name_base = 'test_%03d_layout_%03d_s_%03d_a_%03d_pc_thresh' % (subarray,min_spacing_cm,n_ants,zero_spacing_leakage_threshold_percent)
      signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (subarray,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      signal_array = np.load(signal_filename)

      #can't use the nan values:
      good_idx = np.isfinite(signal_array[0,:]) & np.isfinite(signal_array[0,:])
      good_signal_array_short_baselines_Tb = signal_array[0,:][good_idx]
   
      good_freq_MHz_array = freq_MHz_array[good_idx]
   
      coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
      ffit = poly.polyval(good_freq_MHz_array, coefs)
   
      #in log log space:
      log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
      log_freq_MHz_array = np.log10(good_freq_MHz_array)
      coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
   
      #residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
   
      rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
      rms_of_residuals_list.append(rms_of_residuals)
      random_array_number_list.append(subarray)
      
   
      plt.clf()
      plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
      map_title="Weighted residual for log polynomial order %s fit" % poly_order
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s_log_fit_residual_weighted_poly_%s.png" % (plot_name_base,poly_order)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
   
   rms_of_residuals_array = np.asarray(rms_of_residuals_list)
   random_array_number_array = np.asarray(random_array_number_list)
   
   plt.clf()
   plt.plot(random_array_number_array,rms_of_residuals_array,label='residual from log fit')

   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.xlabel("random array layout number")
   plt.ylabel("rms residuals (K)")
   plt.legend(loc=1)
   #plt.ylim([0, 3.5])
   fig_name= "%s_n_layouts_vs_residuals_poly_%s.png" % (final_plot_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   #now go through and selectively add arrays
   #find the lowest rms array:
   #index_min = np.argmin(rms_of_residuals_array)
   #random_array_number_of_lowest_rms = random_array_number_array[index_min]
   
   #sort the rms array and array number array by order of ascending rms
   #zip the two arrays:
   rms_array_number_zip = zip(rms_of_residuals_array,random_array_number_array)
   #sort ascending according to rms
   rms_array_number_zip.sort()
   #print rms_array_number_zip
   #array_number_sorted = [x for y, x in rms_array_number_zip]
   rms_of_residuals_array_sorted = zip(*rms_array_number_zip)[0]
   random_array_number_array_sorted = zip(*rms_array_number_zip)[1]
   #print rms_of_residuals_array_sorted
   #print random_array_number_array_sorted
  
   random_array_number_of_lowest_rms = int(random_array_number_array_sorted[0])
   lowest_rms = rms_of_residuals_array_sorted[0]
   
   rms_of_residuals_selective_list = []
   random_array_numbers_used_selective_list = []
   
   first_signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (random_array_number_of_lowest_rms,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
   first_signal_array = np.load(first_signal_filename)
   selective_sum_signal_array = first_signal_array*0.0
   
   #can't use the nan values:
   good_idx = np.isfinite(first_signal_array[0,:]) & np.isfinite(first_signal_array[0,:])
   good_signal_array_short_baselines_Tb = first_signal_array[0,:][good_idx]
   
   good_freq_MHz_array = freq_MHz_array[good_idx]
   
   coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
   ffit = poly.polyval(good_freq_MHz_array, coefs)
   
   #in log log space:
   log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
   log_freq_MHz_array = np.log10(good_freq_MHz_array)
   coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
   ffit = poly.polyval(log_freq_MHz_array, coefs)
   ffit_linear = 10**ffit
   
   #residual = log_signal_array_short_baselines - log_ffit
   residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
   
   rms_of_residuals_selective = np.sqrt(np.mean(residual_of_log_fit**2))
   rms_of_residuals_selective_list.append(rms_of_residuals_selective)
   random_array_numbers_used_selective_list.append(random_array_number_of_lowest_rms)
   
   selective_sum_signal_array += first_signal_array
   current_rms = rms_of_residuals_selective
   current_n_arrays_used = 1.0
   
   #Then step through all the rest of the arrays:
   for array_number in random_array_number_array_sorted[1:]:
      new_signal_filename = 'array%03d_s_%03d_a_%03d_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (array_number,min_spacing_cm,n_ants,zero_spacing_leakage_threshold)
      new_signal_array = np.load(new_signal_filename)
      new_selective_sum_signal_array = selective_sum_signal_array + new_signal_array
      new_selective_average_signal_array = new_selective_sum_signal_array / (current_n_arrays_used + 1.0)
      
      #can't use the nan values:
      good_idx = np.isfinite(new_selective_average_signal_array[0,:]) & np.isfinite(new_selective_average_signal_array[0,:])
      good_signal_array_short_baselines_Tb = new_selective_average_signal_array[0,:][good_idx]
      
      good_freq_MHz_array = freq_MHz_array[good_idx]
      
      coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
      ffit = poly.polyval(good_freq_MHz_array, coefs)
      
      #in log log space:
      log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
      log_freq_MHz_array = np.log10(good_freq_MHz_array)
      coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      
      #residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
      new_rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
      
      #only append if rms has been lowered
      if (new_rms_of_residuals < current_rms):
         current_rms = new_rms_of_residuals
         current_n_arrays_used += 1.0
         selective_sum_signal_array = new_selective_sum_signal_array
         rms_of_residuals_selective_list.append(current_rms)
         random_array_numbers_used_selective_list.append(array_number)
   
   rms_of_residuals_selective_array = np.asarray(rms_of_residuals_selective_list)
   #print rms_of_residuals_selective_list
   #print random_array_numbers_used_selective_list
         
   #expect noise to go down as sqrt of n_ants? - probly not
   n_arrays_array = np.arange(1,int(len(random_array_numbers_used_selective_list))+1)
   expected_residuals = rms_of_residuals_selective_list[0]/np.sqrt(n_arrays_array)
   
   plt.clf()
   plt.plot(n_arrays_array,rms_of_residuals_selective_array,label='residual from log fit')
   plt.plot(n_arrays_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.xlabel("number of arrays used")
   plt.ylabel("rms residuals (K)")
   plt.legend(loc=1)
   plt.ylim([0, lowest_rms+0.2])
   fig_name= "%s_n_arrays_vs_residuals_poly_%s_selective.png" % (final_plot_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)     

def model_and_plot_assassin_residuals_selective(n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length_m,zero_spacing_leakage_threshold):
   freq_MHz_array = np.asarray(freq_MHz_list)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   max_arm_length_cm = (max_arm_length_m * 100.)
   min_arm_length_cm = (min_arm_length_m * 100.)
   final_plot_name_base = 'assassin_%02d_ants_per_m_%02d_circles_%03d_to_%03d_length_%03d_pc_thresh' % (n_ants_per_m_of_circumference,n_circles,min_arm_length_cm,max_arm_length_cm,zero_spacing_leakage_threshold_percent)
   rms_of_residuals_list = []
   array_number_list = []
   array_length_list = []
   array_angle_list = []
   #angles are anticlockwise from East
   radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
   array_number = 0
   for circle_number in range(0,n_circles):
      #circle 1 is the smallest
      radius = (min_arm_length_m + circle_number * radial_spacing) 
      diameter = radius * 2.
      diameter_cm = int(diameter*100.)
      #work out circumference 
      circumference = math.pi * diameter
      #print diameter
      n_angles = int(round(circumference * n_ants_per_m_of_circumference))
      angle_increment = (2.*math.pi)/n_angles
      #(remember only need half of them!)
      angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
      #angle_array_deg = angle_array_rad / math.pi * 180.
      #print angle_array_deg
      for angle_rad in angle_array_rad:
         array_number += 1
         angle_deg = int(angle_rad/np.pi*180.)
         signal_filename = 'assassin_%03d_deg_%03d_cm_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (angle_deg,diameter_cm,zero_spacing_leakage_threshold)
         signal_array = np.load(signal_filename)
   
         #can't use the nan values:
         good_idx = np.isfinite(signal_array[0,:]) & np.isfinite(signal_array[0,:])
         good_signal_array_short_baselines_Tb = signal_array[0,:][good_idx]
      
         good_freq_MHz_array = freq_MHz_array[good_idx]
      
         coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
         ffit = poly.polyval(good_freq_MHz_array, coefs)
      
         #in log log space:
         log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
         log_freq_MHz_array = np.log10(good_freq_MHz_array)
         coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
         ffit = poly.polyval(log_freq_MHz_array, coefs)
         ffit_linear = 10**ffit
      
         #residual = log_signal_array_short_baselines - log_ffit
         residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
      
         rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
         rms_of_residuals_list.append(rms_of_residuals)
         array_length_list.append(diameter_cm)
         array_angle_list.append(angle_deg)
         array_number_list.append(array_number)
         

   rms_of_residuals_array = np.asarray(rms_of_residuals_list)
   array_number_array = np.asarray(array_number_list)
   array_length_array = np.asarray(array_length_list)
   array_angle_array = np.asarray(array_angle_list)
   
   array_length_angle_zip = zip(array_length_array,array_angle_array)
   
   #now go through and selectively add arrays
   #find the lowest rms array:
   #index_min = np.argmin(rms_of_residuals_array)
   #random_array_number_of_lowest_rms = random_array_number_array[index_min]
   
   #sort the rms array and array number array by order of ascending rms
   #zip the two arrays:
   rms_array_number_zip = zip(rms_of_residuals_array,array_number_array)
   rms_array_length_zip = zip(rms_of_residuals_array,array_length_array)
   rms_array_angle_zip = zip(rms_of_residuals_array,array_angle_array)
   
   #sort ascending according to rms
   rms_array_number_zip.sort()
   rms_array_length_zip.sort()
   rms_array_angle_zip.sort()
   
   #print rms_array_number_zip
   #array_number_sorted = [x for y, x in rms_array_number_zip]
   rms_of_residuals_array_sorted = zip(*rms_array_number_zip)[0]
   array_number_array_sorted = zip(*rms_array_number_zip)[1]
   array_length_array_sorted = zip(*rms_array_length_zip)[1]
   array_angle_array_sorted = zip(*rms_array_angle_zip)[1]

   #print rms_of_residuals_array_sorted
   #print array_number_array_sorted
   #print array_length_array_sorted
   #print array_angle_array_sorted
   
   #print rms_of_residuals_array_sorted
   #print random_array_number_array_sorted
  
   array_number_of_lowest_rms = int(array_number_array_sorted[0])
   array_length_of_lowest_rms = int(array_length_array_sorted[0])
   array_angle_of_lowest_rms = int(array_angle_array_sorted[0])

   lowest_rms = rms_of_residuals_array_sorted[0]
   
   rms_of_residuals_selective_list = []
   array_numbers_used_selective_list = []
   array_lengths_used_selective_list = []
   array_angles_used_selective_list = []

   first_signal_filename = 'assassin_%03d_deg_%03d_cm_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (array_angle_of_lowest_rms,array_length_of_lowest_rms,zero_spacing_leakage_threshold)
   
   first_signal_array = np.load(first_signal_filename)
   selective_sum_signal_array = first_signal_array*0.0
   
   #can't use the nan values:
   good_idx = np.isfinite(first_signal_array[0,:]) & np.isfinite(first_signal_array[0,:])
   good_signal_array_short_baselines_Tb = first_signal_array[0,:][good_idx]
   
   good_freq_MHz_array = freq_MHz_array[good_idx]
   
   coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
   ffit = poly.polyval(good_freq_MHz_array, coefs)
   
   #in log log space:
   log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
   log_freq_MHz_array = np.log10(good_freq_MHz_array)
   coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
   ffit = poly.polyval(log_freq_MHz_array, coefs)
   ffit_linear = 10**ffit
   
   #residual = log_signal_array_short_baselines - log_ffit
   residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
   
   best_residual_of_log_fit = residual_of_log_fit
   
   rms_of_residuals_selective = np.sqrt(np.mean(residual_of_log_fit**2))
   rms_of_residuals_selective_list.append(rms_of_residuals_selective)
   array_numbers_used_selective_list.append(array_number_of_lowest_rms)
   array_lengths_used_selective_list.append(array_length_of_lowest_rms)
   array_angles_used_selective_list.append(array_angle_of_lowest_rms)
   
   selective_sum_signal_array += first_signal_array
   current_rms = rms_of_residuals_selective
   current_n_arrays_used = 1.0
   
   print("first signal:")
   print(first_signal_filename)
   
   #print current_rms
   
   #Then step through all the rest of the arrays:
   for array_number_index,array_number in enumerate(array_number_array_sorted[1:]):
      array_length = array_length_array_sorted[1:][array_number_index]
      array_angle = array_angle_array_sorted[1:][array_number_index]
      new_signal_filename = 'assassin_%03d_deg_%03d_cm_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (array_angle,array_length,zero_spacing_leakage_threshold)
      #print new_signal_filename
      new_signal_array = np.load(new_signal_filename)
      new_selective_sum_signal_array = selective_sum_signal_array + new_signal_array
      new_selective_average_signal_array = new_selective_sum_signal_array / (current_n_arrays_used + 1.0)
      
      #can't use the nan values:
      good_idx = np.isfinite(new_selective_average_signal_array[0,:]) & np.isfinite(new_selective_average_signal_array[0,:])
      good_signal_array_short_baselines_Tb = new_selective_average_signal_array[0,:][good_idx]
      
      good_freq_MHz_array = freq_MHz_array[good_idx]
      
      coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
      ffit = poly.polyval(good_freq_MHz_array, coefs)
      
      #in log log space:
      log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
      log_freq_MHz_array = np.log10(good_freq_MHz_array)
      coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
      ffit = poly.polyval(log_freq_MHz_array, coefs)
      ffit_linear = 10**ffit
      
      #residual = log_signal_array_short_baselines - log_ffit
      residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
      new_rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
      
      #print new_rms_of_residuals
      
      #only append if rms has been lowered
      if (new_rms_of_residuals < current_rms):
         best_residual_of_log_fit = residual_of_log_fit
         current_rms = new_rms_of_residuals
         current_n_arrays_used += 1.0
         selective_sum_signal_array = new_selective_sum_signal_array
         rms_of_residuals_selective_list.append(current_rms)
         array_numbers_used_selective_list.append(array_number)
   
   rms_of_residuals_selective_array = np.asarray(rms_of_residuals_selective_list)
   #print rms_of_residuals_selective_list
   #print random_array_numbers_used_selective_list
         
   #expect noise to go down as sqrt of n_ants? - probly not
   n_arrays_array = np.arange(1,int(len(array_numbers_used_selective_list))+1)
   expected_residuals = rms_of_residuals_selective_list[0]/np.sqrt(n_arrays_array)
   
   plt.clf()
   plt.plot(n_arrays_array,rms_of_residuals_selective_array,label='residual from log fit')
   plt.plot(n_arrays_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.xlabel("number of arrays used")
   plt.ylabel("rms residuals (K)")
   plt.legend(loc=1)
   plt.ylim([0, lowest_rms+0.2])
   fig_name= "%s_n_arrays_vs_residuals_poly_%s_selective.png" % (final_plot_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   good_idx = np.isfinite(best_residual_of_log_fit) & np.isfinite(best_residual_of_log_fit)
   good_best_residual_of_log_fit = best_residual_of_log_fit[good_idx]  
   good_freq_MHz_array = freq_MHz_array[good_idx]
   
   #plot the final best residual plot 
   plt.clf()
   plt.plot(good_freq_MHz_array,best_residual_of_log_fit,label='residual from log fit')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("freq (MHz)")
   plt.legend(loc=1)
   fig_name= "%s_best_log_fit_residual_poly_%s_selective.png" % (final_plot_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

def model_and_plot_assassin_residuals(n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length_m,zero_spacing_leakage_threshold):
   freq_MHz_array = np.asarray(freq_MHz_list)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   max_arm_length_cm = (max_arm_length_m * 100.)
   min_arm_length_cm = (min_arm_length_m * 100.)
   final_plot_name_base = 'assassin_%02d_ants_per_m_%02d_circles_%03d_to_%03d_length_%03d_pc_thresh' % (n_ants_per_m_of_circumference,n_circles,min_arm_length_cm,max_arm_length_cm,zero_spacing_leakage_threshold_percent)
   rms_of_residuals_list = []
   array_number_list = []
   array_length_list = []
   array_angle_list = []
   #angles are anticlockwise from East
   radial_spacing = (max_arm_length_m - min_arm_length_m) / (n_circles-1)
   array_number = 0
   
   signal_array_sum = np.zeros((1,n_chan))
   signal_array_sum_nan_to_num = np.zeros((1,n_chan))
   arrays_used_per_freq_sum = np.zeros((1,n_chan))
   
   
   for circle_number in range(0,n_circles):
      #circle 1 is the smallest
      radius = (min_arm_length_m + circle_number * radial_spacing) 
      diameter = radius * 2.
      diameter_cm = int(diameter*100.)
      #work out circumference 
      circumference = math.pi * diameter
      #print diameter
      n_angles = int(round(circumference * n_ants_per_m_of_circumference))
      angle_increment = (2.*math.pi)/n_angles
      #(remember only need half of them!)
      angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
      #angle_array_deg = angle_array_rad / math.pi * 180.
      #print angle_array_deg
      for angle_rad in angle_array_rad:
         array_number += 1
         angle_deg = int(angle_rad/np.pi*180.)
         signal_filename = 'assassin_%03d_deg_%03d_cm_X_lst_2.00_hr_int_0.13_hr_D_gmoss_thresh_%0.2f_signal_weighted_Tb.npy' % (angle_deg,diameter_cm,zero_spacing_leakage_threshold)
         signal_array = np.load(signal_filename)
         
         arrays_used_per_freq = signal_array*0
         arrays_used_per_freq[signal_array>=0] = 1
         arrays_used_per_freq[signal_array<0] = 1
         arrays_used_per_freq_nan_to_num = np.nan_to_num(arrays_used_per_freq)
         arrays_used_per_freq_sum += arrays_used_per_freq_nan_to_num
         #
         signal_array_nan_to_num = np.nan_to_num(signal_array)
         signal_array_sum_nan_to_num += signal_array_nan_to_num
         
         signal_array_sum += signal_array
         array_length_list.append(diameter_cm)
         array_angle_list.append(angle_deg)
         array_number_list.append(array_number)
         
   n_arrays = float(len(array_number_list))
   ##signal_array_average = signal_array_sum/n_arrays
   signal_array_average = signal_array_sum_nan_to_num/arrays_used_per_freq_sum
   
   #print signal_array_average
   
   array_number_array = np.asarray(array_number_list)
   array_length_array = np.asarray(array_length_list)
   array_angle_array = np.asarray(array_angle_list)
   
   
   #can't use the nan values:
   good_idx = np.isfinite(signal_array_average[0,:]) & np.isfinite(signal_array_average[0,:])
   good_signal_array_short_baselines_Tb = signal_array_average[0,:][good_idx]
   
   good_freq_MHz_array = freq_MHz_array[good_idx]
   
   coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
   ffit = poly.polyval(good_freq_MHz_array, coefs)
   
   #in log log space:
   log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
   log_freq_MHz_array = np.log10(good_freq_MHz_array)
   coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
   ffit = poly.polyval(log_freq_MHz_array, coefs)
   ffit_linear = 10**ffit
   
   #residual = log_signal_array_short_baselines - log_ffit
   residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
   
   rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))

   print('rms of residuals %0.3f ' % rms_of_residuals)
   
   #plot the final best residual plot 
   plt.clf()
   plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.ylabel("Residual Tb (K)")
   plt.xlabel("freq (MHz)")
   plt.legend(loc=1)
   fig_name= "%s_log_fit_residual_poly_%s_all.png" % (final_plot_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
def plot_rms_residuals_for_arrays(n_arrays,min_spacing_m,zero_spacing_leakage_threshold,addition_type,n_ants):
   #addition_type can be 'simple' or 'coherent'
   min_spacing_cm = int(min_spacing_m*100.)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   model_vis_name_base = 'test_%03d_arrays_combined_%03d_s_%03d_a_%03d_pc_thresh_%s' % (n_arrays,min_spacing_cm,n_ants,zero_spacing_leakage_threshold_percent,addition_type)
   rms_of_residuals_list = []
   n_arrays_list = []
   for n_subarrays in range(1,n_arrays+1):
      rms_of_residuals = model_average_from_arrays(n_subarrays,min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,addition_type=addition_type,n_ants=n_ants)
      rms_of_residuals_list.append(rms_of_residuals)
      n_arrays_list.append(n_subarrays)
      rms_of_residuals_array = np.asarray(rms_of_residuals_list)
      n_arrays_array = np.asarray(n_arrays_list)
    
    
   #expect noise to go down as sqrt of n_arrays
   expected_residuals = rms_of_residuals_array[0]/np.sqrt(n_arrays_array)
   
   plt.clf()
   plt.plot(n_arrays_array,rms_of_residuals_array,label='residual from log fit')
   plt.plot(n_arrays_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.xlabel("number of sub arrays")
   plt.ylabel("rms residuals (K)")
   plt.legend(loc=1)
   plt.ylim([0, 3.5])
   fig_name= "%s_n_arrays_vs_residuals_poly_%s.png" % (model_vis_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name) 


def plot_rms_residuals_for_n_ants(n_ant_min,n_ant_max,min_spacing_m,zero_spacing_leakage_threshold):
   min_spacing_cm = int(min_spacing_m*100.)
   zero_spacing_leakage_threshold_percent = zero_spacing_leakage_threshold*100.
   model_vis_name_base = 'test_%03d_min_ant_%03d_max_ant_%03d_s_%03d_pc_thresh' % (n_ant_min,n_ant_max,min_spacing_cm,zero_spacing_leakage_threshold_percent)
   rms_of_residuals_list = []
   n_ants_list = []
   for n_ants in range(n_ant_min,n_ant_max+1):
      rms_of_residuals = model_average_from_arrays(1,min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,addition_type='simple',n_ants=n_ants)
      rms_of_residuals_list.append(rms_of_residuals)
      n_ants_list.append(n_ants)
   rms_of_residuals_array = np.asarray(rms_of_residuals_list)
   n_ants_array = np.asarray(n_ants_list)
    
    
   #expect noise to go down as sqrt of n_ants? - probly not
   #expected_residuals = rms_of_residuals_array[0]/np.sqrt(n_arrays_array)
   
   plt.clf()
   plt.plot(n_ants_array,rms_of_residuals_array,label='residual from log fit')
   #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
   map_title="Weighted residual for log polynomial order %s fit" % poly_order
   plt.xlabel("number of ants in subarray")
   plt.ylabel("rms residuals (K)")
   plt.legend(loc=1)
   plt.ylim([0, 20])
   fig_name= "%s_n_ants_vs_residuals_poly_%s.png" % (model_vis_name_base,poly_order)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 

def simulate_and_extract_assassin_baselines(n_ants_per_m_of_circumference,n_circles,max_arm_length_m,min_arm_length):
   #angles are anticlockwise from East
   radial_spacing = (max_arm_length_m - min_arm_length) / (n_circles-1)
   for circle_number in range(0,n_circles):
      #circle 1 is the smallest
      radius = (min_arm_length + circle_number * radial_spacing) 
      diameter = radius * 2.
      diameter_cm = int(diameter*100.)
      #work out circumference 
      circumference = math.pi * diameter
      #print diameter
      n_angles = int(round(circumference * n_ants_per_m_of_circumference))
      angle_increment = (2.*math.pi)/n_angles
      #(remember only need half of them!)
      angle_array_rad = np.arange(1,n_angles/2+1) * angle_increment
      #angle_array_deg = angle_array_rad / math.pi * 180.
      #print angle_array_deg
      for angle_rad in angle_array_rad:
         angle_deg = int(angle_rad/np.pi*180.)
         antenna_position_filename = "assassin_baseline_%03d_deg_%03d_cm.txt" % (angle_deg,diameter_cm)
         ant_1_x_offset = radius * np.cos(angle_rad)
         ant_1_y_offset = radius * np.sin(angle_rad)
         ant_1_position_string = "%0.3f   %0.3f   0\n" % (ant_1_x_offset,ant_1_y_offset)
         ant_2_angle = angle_rad + np.pi
         ant_2_x_offset = radius * np.cos(ant_2_angle)
         ant_2_y_offset = radius * np.sin(ant_2_angle)
         ant_2_position_string = "%0.3f   %0.3f   0" % (ant_2_x_offset,ant_2_y_offset)         
         with open(antenna_position_filename,'w') as f:
            f.write(ant_1_position_string)
            f.write(ant_2_position_string)
         print("wrote %s" % antenna_position_filename)
      
         array_label = "assassin_%03d_deg_%03d_cm" % (angle_deg,diameter_cm)
         plot_antenna_array(array_layout_filename=antenna_position_filename,ylim=2)
         simulate(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=antenna_position_filename,array_label=array_label)
         extract_signal_from_sims(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=antenna_position_filename,array_label=array_label)
    
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
            
def make_image_movie_from_ds9(EDA2_chan_list,n_obs_list,out_movie_name):
   #Need to actually plug in linux box monitor and use native screen (cant do over ssh from Mac)
   image_counter = 0
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):
      obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
      for obs_time_index,obs_time in enumerate(obs_time_list):
         image_name = "cal_chan_%s_%s_ms-XX-image.fits" % (EDA2_chan,obs_time)
         if os.path.exists(image_name):
            image_counter += 1
            print(image_name)
            freq_MHz = float(EDA2_chan) * (50./64)
             
            out_png_image_name="cal_image_%03d_ds9.png" % (image_counter)
            
            #cmd = "ds9 -fits %s_xx_restor.fits -scale limits 1 100 -cmap invert yes -colorbar yes -grid yes -grid axes type exterior -export jpg %s_xx_restor.jpg -exit " % (image_basename,image_basename)
            #"-scale limits %s %s" % (chan_scale_min_list_xx[chan_index],chan_scale_max_list_xx[chan_index])
            cmd = "ds9 %s -invert -colorbar yes -view buttons yes -view panner yes -view magnifier yes -view info yes -view filename yes  -zoom 1 -width 512 -height 512 -grid  -saveimage png %s -quit" % (image_name,out_png_image_name)
            print(cmd)
            os.system(cmd)
            
            #Add some text to the png
            img = Image.open("%s" % out_png_image_name)
            draw = ImageDraw.Draw(img)
            ##font = ImageFont.truetype(<font-file>, <font-size>)
            ##font = ImageFont.truetype("sans-serif.ttf", 16)
            font = ImageFont.truetype('FreeSans.ttf',30)
            ##draw.text((x, y),"Sample Text",(r,g,b))
            draw.text((10, 10),"EDA2 chan %s (%.0f MHz)\nTime %s" % (EDA2_chan,freq_MHz,obs_time),(0,0,0),font=font)
            ##draw.text((256, 256),"Channel %s" % EDA2_chan,(0,0,0))
            img.save("%s" % out_png_image_name)
   
   cmd = "ffmpeg -framerate 2 -i cal_image_%03d_ds9.png -c:v libx264 -r 30 -pix_fmt yuv420p out_movie.mp4" 
   print(cmd)
   os.system(cmd)
   
   cmd = "mv out_movie.mp4 %s" % out_movie_name
   print(cmd)
   os.system(cmd)


#SIMS

#calculate the global 21cm signal:
s_21_array = plot_S21(nu_array=freq_MHz_list,C=C,A=A,delta_nu=delta_nu,nu_c=nu_c)
s_21_array_EDGES = plot_S21_EDGES(nu_array=freq_MHz_list)

#assassin:
#simulate_and_extract_assassin_baselines(n_ants_per_m_of_circumference=2,n_circles=5,max_arm_length_m=1.5,min_arm_length=0.325)
#sys.exit()
#model_and_plot_assassin_residuals(n_ants_per_m_of_circumference=2,n_circles=5,max_arm_length_m=1.5,min_arm_length_m=0.325,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold)


#simulate_assassin(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,n_ants_per_m_of_circumference=2,n_circles=5,max_arm_length_m=1.5,min_arm_length_m=0.325)
#extract_signal_from_assassin(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,n_ants_per_m_of_circumference=2,n_circles=5,max_arm_length_m=1.5,min_arm_length_m=0.325,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold)
#model_signal_from_assassin(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,n_ants_per_m_of_circumference=2,n_circles=5,max_arm_length_m=1.5,min_arm_length_m=0.325,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,poly_order=poly_order)
#sys.exit()

#EDA2 sims to compare to Singh et al 2015 


#EDA2 sims and new data extraction

#############################################################
#Welcome back for 2020, remember the latest sim results are in: /md0/EoR/ASSASSIN/solve_for_tsky_weighted/
#and the EDA data work is in: /md0/EoR/EDA2/20191213_data/90
#next thing average EDA2 data in time
#############################################################



#lst_hrs_list = ['2.2','2.4','2.6']
#lst_hrs_list = ['2']

#freq = 400/512 * chan number
#freq_MHz_list=[50.]

#lst_hrs_list = ['2.0','2.2','2.4','2.6']
lst_hrs_list = ['2']

EDA2_data = True

#EDA2_filenames = ["chan_64_20191202T171525_calibrated.uvfits","chan_77_20191202T171629_calibrated.uvfits","chan_90_20191202T171727_calibrated.uvfits","chan_103_20191202T171830_calibrated.uvfits","chan_116_20191202T171928_calibrated.uvfits","chan_129_20191202T172027_calibrated.uvfits"]

#20190929 data:
#EDA2_chan_list = [64,77,90,103,116,129]

#20200217 data (dont use 63 / 49 MHz dont have beam model! missing 127 and 128 so just do to incl 126)
#20200303:
EDA2_chan_list = range(64,127)
#20200304:

#test progress 
#EDA2_chan_list = range(64,130)
#EDA2_chan_list = range(64,74)
#middle time for each chan
#EDA2_obs_time_list = ["20191202T171525","20191202T171629","20191202T171727","20191202T171830","20191202T171928","20191202T172027"]

##all times for each chan
#EDA2_obs_time_list_each_chan = [
#   ['20191202T171525','20191202T171530','20191202T171535','20191202T171540','20191202T171545','20191202T171550','20191202T171554','20191202T171559','20191202T171604','20191202T171609','20191202T171614','20191202T171618','20191202T171624'],
#   ['20191202T171629','20191202T171633','20191202T171638','20191202T171643','20191202T171648','20191202T171653','20191202T171658','20191202T171702','20191202T171707','20191202T171712','20191202T171717','20191202T171722'],
#   ['20191202T171727','20191202T171732','20191202T171737','20191202T171741','20191202T171746','20191202T171751','20191202T171756','20191202T171801','20191202T171805','20191202T171810','20191202T171815','20191202T171820','20191202T171825'],
#   ['20191202T171830','20191202T171835','20191202T171840','20191202T171845','20191202T171849','20191202T171854','20191202T171859','20191202T171904','20191202T171909','20191202T171913','20191202T171918','20191202T171923'],
#   ['20191202T171928','20191202T171933','20191202T171938','20191202T171943','20191202T171948','20191202T171952','20191202T171957','20191202T172002','20191202T172007','20191202T172012','20191202T172017','20191202T172021'],
#   ['20191202T172027','20191202T172032','20191202T172036','20191202T172041','20191202T172046','20191202T172051','20191202T172056','20191202T172100','20191202T172105','20191202T172110','20191202T172115','20191202T172120','20191202T172125']
#   ]

#default:
#SOme of these do not calibrate, so exclude them (best results so far):
#EDA2_obs_time_list_each_chan = [
#   ['20191202T171525','20191202T171530','20191202T171535','20191202T171540','20191202T171545','20191202T171550','20191202T171554','20191202T171559','20191202T171604','20191202T171609','20191202T171614','20191202T171618'],
#   ['20191202T171629','20191202T171633','20191202T171638','20191202T171643','20191202T171648','20191202T171653','20191202T171658','20191202T171702','20191202T171707','20191202T171712','20191202T171717','20191202T171722'],
#   ['20191202T171727','20191202T171732','20191202T171737','20191202T171741','20191202T171746','20191202T171751','20191202T171756','20191202T171801','20191202T171805','20191202T171810','20191202T171815','20191202T171820','20191202T171825'],
#   ['20191202T171830','20191202T171835','20191202T171840','20191202T171845','20191202T171849','20191202T171854','20191202T171859','20191202T171904','20191202T171909','20191202T171913','20191202T171918','20191202T171923'],
#   ['20191202T171928','20191202T171933','20191202T171938','20191202T171943','20191202T171948','20191202T171952','20191202T171957','20191202T172002','20191202T172007','20191202T172012','20191202T172017','20191202T172021'],
#   ['20191202T172027','20191202T172032','20191202T172041','20191202T172115']
#   ]

#Use older data at 50 MHz from Marcin. UTC chan_64_20190929T213132.uvfits is the same LST (5.8666 hrs)
#they are space 5 mins apart, so just use 2 for now
#20190929T213132
#20190929T213632
#
#20190929_213132_eda2_ch32_ant256_midday_avg8140.ms
#20190929_213632_eda2_ch32_ant256_midday_avg8140.ms

#20190929 data (sub 20190929 data for 50MHz):
#EDA2_obs_time_list_each_chan = [
#   ['20190929T213132','20190929T213632'],
#   ['20191202T171629','20191202T171633','20191202T171638','20191202T171643','20191202T171648','20191202T171653','20191202T171658','20191202T171702','20191202T171707','20191202T171712','20191202T171717','20191202T171722'],
#   ['20191202T171727','20191202T171732','20191202T171737','20191202T171741','20191202T171746','20191202T171751','20191202T171756','20191202T171801','20191202T171805','20191202T171810','20191202T171815','20191202T171820','20191202T171825'],
#   ['20191202T171830','20191202T171835','20191202T171840','20191202T171845','20191202T171849','20191202T171854','20191202T171859','20191202T171904','20191202T171909','20191202T171913','20191202T171918','20191202T171923'],
#   ['20191202T171928','20191202T171933','20191202T171938','20191202T171943','20191202T171948','20191202T171952','20191202T171957','20191202T172002','20191202T172007','20191202T172012','20191202T172017','20191202T172021'],
#   ['20191202T172027','20191202T172032','20191202T172041','20191202T172115']
#   ]
   


##for use in testing when you just want 1 obs per freq
#EDA2_obs_time_list_each_chan = [
#   ['20191202T171525','20191202T171530'],
#   ['20191202T171629'],
#   ['20191202T171727'],
#   ['20191202T171830'],
#   ['20191202T171928'],
#   ['20191202T172027'],
#   ]




#EDA2_data_dir = '/md0/EoR/EDA2/20191213_data/'
EDA2_data_dir = '/md0/EoR/EDA2/20200303_data/'
#inv fine chans
#EDA2_data_dir = '/md0/EoR/EDA2/inv_uvfits/20200303_213605/'
#20200217 data (don't use 'edge' times):
#EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan("/md0/EoR/EDA2/20200303_data/",EDA2_chan_list)
#EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan("/md0/EoR/EDA2/20200304_data/",EDA2_chan_list)
EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan(EDA2_data_dir,EDA2_chan_list)


EDA2_obs_time_list_each_chan = EDA2_obs_time_list_each_chan[0:]

n_obs_concat_list = [len(obs_list) for obs_list in EDA2_obs_time_list_each_chan] 

EDA2_obs_time_list = [item[0] for item in EDA2_obs_time_list_each_chan] 


#just for SIMS for EDA2 data (REMOVE!!!), have to manually cd to chan dir and run just simulate:
#EDA2_chan_list = [129]
#EDA2_obs_time_list = ['20191202T172027']

EDA2_chan_list_array = np.asarray(EDA2_chan_list)

#EDA2_obs_time = '20191202T171727'
#freq_MHz_array = np.round(400./512.*EDA2_chan_list_array)
freq_MHz_array = 400./512.*EDA2_chan_list_array
freq_MHz_list = freq_MHz_array[0:]
EDA2_chan_list = EDA2_chan_list[0:]

print(freq_MHz_list)
print(EDA2_chan_list)



#eda2_data_filename = "chan_%s_%s.uvfits" % (int(EDA2_chan),EDA2_obs_time)
#eda2_data_filename = "chan_90_20191202T171727.uvfits"
#eda2_data_uvfits_name_list=[eda2_data_filename]




##need to use an old year ...
lst_hrs_list = []
for EDA2_obs_time_index,EDA2_obs_time in enumerate(EDA2_obs_time_list):
   #there might have been no obs:
   if EDA2_obs_time!=0:
      print(EDA2_obs_time)
      lst_eda2_hrs = "%0.1f" % get_eda2_lst("2015%s" % EDA2_obs_time[4::])
      #print(lst_eda2_hrs)
      lst_hrs_list.append(lst_eda2_hrs)
   else:
      #just use the first LST
      lst_eda2_hrs = "%0.1f" % get_eda2_lst("2015%s" % EDA2_obs_time_list[0][4::])
      #print(lst_eda2_hrs)
      lst_hrs_list.append(lst_eda2_hrs)


#EDA2
#plot_antenna_array(array_layout_filename=array_ant_locations_filename)
#sys.exit()
#plot_baseline_length_counts(array_layout_filename = array_ant_locations_filename,freq_MHz=50.)
#sys.exit()

#Step 1: simulate


#cd into each chan dir separately and run simulate to get apparent sky images (don't worry that it crashes on concat freq step)
#need to fix this so you can just run like the other functoins below for multiple eda2 chans
#for EDA2 chans [64,77,90,103,116,129], freq_MHz_list = [ 50.  60.  70.  80.  91. 101.]
#sims:
#freq_MHz_list=[50.]
#lst_hrs_list = ['2']
#simulate(lst_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label,EDA2_data=False)

#DATA: (repeat twice with 'diffuse' then 'global_unity')
#pol_list = ['Y']
#chan_num = 0
#freq_MHz_list = [freq_MHz_array[chan_num]]
###if FAST: for data need to simulate with 'global_unity' and then separately 'diffuse' (only if fast)
#for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
#   if len(freq_MHz_list)==1:
#      EDA2_chan = EDA2_chan_list[chan_num]
#   else:
#      EDA2_chan = EDA2_chan_list[freq_MHz_index]
#   new_dir = "./%s" % EDA2_chan
#   os.chdir(new_dir)
#   freq_MHz_input_list = [freq_MHz]
#   lst_hrs_list_input = [lst_hrs_list[freq_MHz_index]]
#   simulate(lst_list=lst_hrs_list_input,freq_MHz_list=freq_MHz_input_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label,EDA2_data=True)
#   #the concat step is causing /tmp to fill with casa crash reports
#   cmd = "rm -rf /tmp/*" 
#   print(cmd)
#   os.system(cmd)
#   os.chdir('./..')
##sys.exit()

#Step 2: calibrate

##calibrate each individually first and concat
##do this outside chan dir
##if doing individual chans:
#chan_num = 0
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]
#plot_cal = False
#wsclean = False
#wsclean = True
#concat=True
#per_chan_cal = True
#calibrate_eda2_data(EDA2_chan_list=EDA2_chan_list,obs_type='night',lst_list=lst_hrs_list,pol_list=pol_list,n_obs_concat_list=n_obs_concat_list,concat=concat,wsclean=wsclean,plot_cal=plot_cal,uv_cutoff=0,per_chan_cal=per_chan_cal)
#sys.exit()

#Need to plug in monitor to namorrodor, can't do this with nohup or remotely
#make_image_movie_from_ds9(EDA2_chan_list,n_obs_concat_list,'20200303_data.mp4')
#sys.exit()
 
#after calibration:
#no cal solutions for 64/chan_64_20191202T171624.vis
#no cal solutions for 129/chan_129_20191202T172036.vis
#no cal solutions for 129/chan_129_20191202T172046.vis
#no cal solutions for 129/chan_129_20191202T172051.vis
#no cal solutions for 129/chan_129_20191202T172056.vis
#no cal solutions for 129/chan_129_20191202T172100.vis
#no cal solutions for 129/chan_129_20191202T172105.vis
#no cal solutions for 129/chan_129_20191202T172110.vis
#no cal solutions for 129/chan_129_20191202T172120.vis
#no cal solutions for 129/chan_129_20191202T172125.vis

#then concat
##concat data 
#concat_EDA2_data(EDA2_chan_list=EDA2_chan_list,EDA2_obs_time_list_each_chan=EDA2_obs_time_list_each_chan,pol='X')
#sys.exit()
#


#plot_EDA2_cal_sols('cal_64_ph.txt','cal_64_amp.txt')
#sys.exit()

#model_type = 'OLS_with_intercept'
#model_type = 'mixedlm'

#model_type_list = ['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
model_type_list = ['OLS_fixed_intercept']
#model_type = 'OLS_fixed_int_subtr_Y'

#model_type = 'OLS_fixed_int_min_vis'
#model_type = 'OLS_with_int_min_vis'

#model_type = 'WLS'
#model_type = 'OLS_global_angular'


#look here: https://towardsdatascience.com/when-and-how-to-use-weighted-least-squares-wls-models-a68808b1a89d

#to save time:
#poor man's parallel programming:
#run this 5 seperate times, in 5 seperate screens, one for each freq with plot_only = False to get saved data products, then run with plot_only=True once with all freqs
#dont need to change lst_hrs_list, for eda2 data processing only the first lst ios used
#freq_MHz_list = [freq_MHz_list[5]]
#EDA2_chan_list = [EDA2_chan_list[5]]


#also need to change EDA2_obs_time_list_each_chan[x:] above, this is just so the correct first obsid is selected (only good for concat data, more than 1 obs)
#if doing individual chans:
#EDA2 data:
#freq_MHz_list = [freq_MHz_array[0]]
#EDA2_chan_list = [EDA2_chan_list[0]]

#for sims:
#freq_MHz_list = np.arange(start_chan,start_chan+n_chan,chan_step)
#freq_MHz_array = np.asarray(freq_MHz_list)
#lst_hrs_list=['2']
#poly_order_list=[5,6,7]
#poly_order=7

#plot_iso_ant_int_response()
#sys.exit()

plot_only = True
baseline_length_thresh_lambda = 0.5
include_angular_info = True


#up to here with plot_only = False
#chan_num = 90 - 64
#chan_num = 0
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]
#freq_MHz_list = freq_MHz_array[chan_num:chan_num+35]
#EDA2_chan_list = EDA2_chan_list[chan_num:chan_num+35]
#wsclean=False # for sims or miriad cal
#sim for paper plot 1 
wsclean=True # for data
fast=False
no_modelling=True
calculate_uniform_response=False
plot_tsky_for_multiple_freqs(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order=poly_order,plot_only=plot_only,include_angular_info=include_angular_info,model_type_list=model_type_list, EDA2_data=EDA2_data,EDA2_chan_list=EDA2_chan_list,n_obs_concat_list=n_obs_concat_list,wsclean=wsclean,fast=fast,no_modelling=no_modelling,calculate_uniform_response=calculate_uniform_response)



###obs seem to underestimate the global temp - Y sub improves it a bit.
#I reckon because orig sims are with a much quieter patch of sky - this is the GP overhead!
#need to sim the actual LST of the obs. Wait for new data from Marcin (current 12/2/2020)


#joint_model_fit_t_sky_measured(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order_list=poly_order_list,global_signal_model=s_21_array_EDGES,plot_only=plot_only)


#model_tsky_from_saved_data(freq_MHz=50,lst_hrs=2,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label)
#take a look at MLE https://towardsdatascience.com/a-gentle-introduction-to-maximum-likelihood-estimation-9fbff27ea12f
#especially just looking at your data
#try OLS and WLS from here:
# https://www.statsmodels.org/dev/examples/notebooks/generated/ols.html
#mixed LM? no since all baselines see all sky vectors the matrices get too big!
#https://www.statsmodels.org/dev/examples/notebooks/generated/mixed_lm_example.html
# https://www.statsmodels.org/dev/examples/notebooks/generated/wls.html


       
#sys.exit()



#EDA2 Sims
#simulate(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label)
#extract_signal_from_sims(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label)
#plot_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label)
#model_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,poly_order=poly_order,sky_model=sky_model,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label)


#calibrate simulated data with noise and errors
#calibrate_eda2_data(eda2_data_uvfits_name_list=['eda_model_X_lst_2.00_hr_int_0.13_hr_noise_diffuse_gmoss_global_gain_errors_concat_lsts.uvfits'],obs_type='sun',lst_list=[2],freq_list=[],pol_list=[])

#arrays:
#new_array = grow_new_array(seed=123,diameter=35,n_ants=256,min_spacing_m=1.5)
#plot_antenna_array(array_layout_filename=new_array)
#what if you halve the array size and the min spacing?
#new_array = grow_new_array(seed=124,diameter=17,n_ants=256,min_spacing_m=0.75)
#plot_antenna_array(array_layout_filename=new_array)


#EDA2
#plot_antenna_array(array_layout_filename=array_ant_locations_filename)

#DATA

#extract_signal_from_eda2_data(eda2_data_uvfits_name_list=eda2_data_uvfits_name_list,outbase_name=extract_from_eda2_data_outbase_name,array_label='eda2_sub48')

#calibrate_eda2_data(eda2_data_uvfits_name_list=eda2_data_uvfits_name_list,obs_type='sun')
#image_eda2_data(eda2_data_uvfits_name_list=eda2_data_uvfits_name_list)
#plot_antenna_array(array_layout_filename=new_array)



##
#with miriad 512 is too many ants! Doesn't matter - just make two arrays of 256 instead
#what if double the number of ants (have to increase diameter to 23m)?
#new_array = grow_new_array(seed=123,diameter=23,n_ants=512,min_spacing_m=0.75)

##############################
#Test to see how residuals behave as you increase the number of antenna stations, using 16 ant stations

##n_ants = 16
#n_ant_min = 16
#n_ant_max = 16
##n_ant_min = 16
##n_ant_max = 16
##n_arrays = 200
#n_arrays = 20
#min_spacing_m = 0.75
#min_spacing_cm = min_spacing_m * 100.
#
#
#for n_ants in range(n_ant_min,n_ant_max+1):
#   for subarray in range(0,n_arrays):
#      array_label = 'array%03d_s_%03d_a_%03d' % (subarray,min_spacing_cm,n_ants)
#      new_array = grow_new_array(seed=subarray,n_ants=n_ants,min_spacing_m=min_spacing_m)
#      plot_antenna_array(array_layout_filename=new_array)
#      simulate(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=new_array,array_label=array_label)
#      extract_signal_from_sims(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=new_array,array_label=array_label)
#      #extract_signal_from_sims_multi_coherent(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=new_array,array_label=array_label,n_arrays=subarray+1,min_spacing_m=min_spacing_m)
#      plot_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,sky_model=sky_model,array_ant_locations_filename=new_array,array_label=array_label)
#      model_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name=outbase_name,poly_order=poly_order,sky_model=sky_model,array_ant_locations_filename=new_array,array_label=array_label)
#   
#   plot_rms_residuals_for_arrays(n_arrays=n_arrays,min_spacing_m=min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,addition_type='simple',n_ants=n_ants)
#
#   model_and_plot_random_array_layout_residuals(n_arrays=n_arrays,min_spacing_m=min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,n_ants=n_ants)



#before trying to add coherently, just average the 8 ones that worked
#model_average_from_arrays(n_arrays=8)


#plot_rms_residuals_for_n_ants(n_ant_min=4,n_ant_max=n_ant_max-1,min_spacing_m=min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold)

#plot_rms_residuals_for_arrays(n_arrays=n_arrays,min_spacing_m=min_spacing_m,zero_spacing_leakage_threshold=zero_spacing_leakage_threshold,addition_type='coherent',n_ants=n_ants)


#extract_signal_from_multiple('test_uvfits_list.txt')


#######
#no longer used:
#compute_weights(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model)


#Grid seach, maybe I want to do this:
#from sklearn.grid_search import ParameterGrid
#param_grid = {'param1': [value1, value2, value3], 'paramN' : [value1, value2, valueM]}
#grid = ParameterGrid(param_grid)
#for params in grid:
#    your_function(params['param1'], params['param2'])

#Or see: https://stackoverflow.com/questions/13370570/elegant-grid-search-in-python-numpy

#Actually - this is pretty helpful for cure fitting:
#http://keflavich.github.io/astr2600_notebooks/Lecture23_DataFitting.html#/10



