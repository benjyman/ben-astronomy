#!/usr/bin/env python3
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
#from PIL import Image
#from PIL import ImageFont
#from PIL import ImageDraw 
import healpy as hp
from pygdsm import GSMObserver
from pygdsm import GlobalSkyModel
from pygdsm import GlobalSkyModel2016
from datetime import datetime, date
import time

from reproject import reproject_from_healpix
import pyfits
from astropy.wcs import WCS
from astropy.io import fits 
from scipy.interpolate import interp2d,griddata,interp1d
from scipy.ndimage import map_coordinates
from scipy import signal
import numpy.polynomial.polynomial as poly
#from pyuvdata import UVData

import random
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
#from sklearn.linear_model import LinearRegression
#import seaborn as sns
from scipy.optimize import curve_fit

#from sklearn.preprocessing import normalize
#avoid astropy time error 
from astropy.utils.iers import conf
conf.auto_max_age = None
from astropy.utils import iers
iers.conf.auto_download = False  
#from astroplan import download_IERS_A
#download_IERS_A()

#i dont think pyrap works under python 3 (old) - use Kariukis ms_utils instead
#import pyrap.tables as pt
#from pyrap.tables import *
#import pyrap
sys.path.append("/md0/code/git/ben-astronomy/ms")
from ms_utils import *

import ephem

from scipy.io import loadmat

#for making miriad files:
import aipy as a

from casacore.tables import table, tablesummary




#color defs for color blindness contrast
#from https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
color_black = '#000000'
color_orange='#E69F00'
color_light_blue = '#56B4E9'
color_green = '#009E73'
color_yellow = '#F0E442'
color_dark_blue = '#0072B2'
color_orange_red = '#D55E00'
color_pink = '#CC79A7'


  
sun_flux_density = 50000.0   #borkowski et al 1982?
#mu_0 = 4.*np.pi*10**(-7)
c = 3.0e8
k = 1.38065e-23
sq_deg_in_1_sr = (180./math.pi)**2
#Need to specify the location of the observatory (MWA ... put in exact EDA later)
# Setup observatory location - in this case, MWA, Australia
mwa_latitude_degrees=-26.70331940
mwa_longitude_degrees=116.67081524
mwa_latitude_rad = float(mwa_latitude_degrees/180.*np.pi)
mwa_longitude_rad = float(mwa_longitude_degrees/180.*np.pi)
#elevation_m=377.83

mwa_latitude_pyephem = "-26:42.199"
mwa_longitude_pyephem = "116:40.2488"
mwa_elevation = 0

mwa_latitude_deg = -26.70331940
mwa_longitude_deg = 116.670575

mwa_latitude_astropy = '-26.7d'
mwa_longitude_astropy = '116.67d'

mwa_latitude_ephem = '-26.7'
mwa_longitude_ephem = '116.67'


eda2_loc = EarthLocation(lat=-26.70*u.deg, lon=116.67*u.deg, height=0*u.m)
utcoffset = 8*u.hour  # Western Australian Standard Time
   
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
fine_chan_width_Hz =  28935 #Hz should be:(400./512.)*1.0e6 
                     
if array_ant_locations == 'eda2':
   #array_ant_locations_filename = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant'
   #miriad want north, east, up
   array_ant_locations_filename = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_NEU.ant'
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
pol_list = ['X','Y']
sky_model = 'gsm'
#sky_model = 'gsm2016'
#sky_model = 'gmoss'

if sky_model=='gsm':
   NSIDE=512 #(to match pygsm)
      
#can be any of these, except if can only have 'diffuse' if not diffuse_global or diffuse_angular
#signal_type_list=['global','global_EDGES','diffuse','noise','gain_errors','diffuse_global','diffuse_angular']
signal_type_list=['diffuse','noise'] #fig9, 10b?
#signal_type_list=['single_point'] #tests with Jack and WODEN
#signal_type_list=['diffuse_global','noise'] #fig7
#signal_type_list=['diffuse_global','diffuse_angular']
#signal_type_list=['diffuse']
#signal_type_list=['global_unity']
#signal_type_list=['diffuse_global','noise','global_EDGES'] #fig8b
#signal_type_list=['global_EDGES','noise'] #fig6b
#signal_type_list=['global_EDGES'] #fig5
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

def print_eda2_systemp(chan_num):
   freq_MHz = freq_MHz_list[chan_num]
   print(freq_MHz)
   systemp = T_180*(freq_MHz/180.0)**beta
   print('systemp %s K' % systemp)
   A_eff = Aeff_for_freq_MHz_list[chan_num]
   print('A_eff %s K' % A_eff)
   JperK = (2.0*k*10**26)/(eta*A_eff)
   print('JperK %s ' % JperK)
   SEFD = systemp * JperK
   print('SEFD %s ' % SEFD)
   
#print_eda2_systemp(0)
#sys.exit() 

#def get_eda2_lst(eda_time_string="20151202T171727"):
#   ###Change this to use ephem instead of astropy
#   #Hack to get LST! use old time (2015) astropy.utils.iers.iers.IERSRangeError: (some) times are outside of range covered by IERS table.
#   #eda_time_string = "20191202T171727"
#   year, month, day, hour, minute, second = eda_time_string[0:4], eda_time_string[4:6],eda_time_string[6:8], eda_time_string[9:11],eda_time_string[11:13],eda_time_string[13:15]
#   eda2_astropy_time_string = '%4d-%02d-%02d %02d:%02d:%02.1d' % (float(year), float(month), float(day), float(hour), float(minute), float(second))
#   print(eda2_astropy_time_string)
#   eda2_astropy_time = Time(eda2_astropy_time_string, scale='utc', location=(mwa_longitude_astropy, mwa_latitude_astropy))
#   #calculate the local sidereal time for the MWA at eda2_astropy_time
#   eda2_obs_lst = eda2_astropy_time.sidereal_time('apparent')
#   eda2_obs_lst_hrs = eda2_obs_lst.value
#   return eda2_obs_lst_hrs

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

def calc_beam_values(theta_array_deg,phi_array_deg,pol,dipole_height_m,wavelength):
   #need angles in radians!!
   theta_array = theta_array_deg/180. * np.pi
   phi_array = phi_array_deg/180. * np.pi
   ##This is for YY dipole: !
   if (pol=='Y'):
      theta_parallel_array=np.arccos(np.sin(theta_array)*np.cos(phi_array))
   else:
   #This is for XX dipole!
      theta_parallel_array=np.arccos(np.sin(theta_array)*np.sin(phi_array))
      
   d_in_lambda = (2. * dipole_height_m)/wavelength
   gp_effect_array = 2.*np.sin(np.pi*d_in_lambda*np.cos(theta_array))
   voltage_parallel_array=np.sin(theta_parallel_array) * gp_effect_array
   short_dipole_parallel_beam_map = voltage_parallel_array**2    
   
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
   
   reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
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
   plt.xlabel("Frequency (MHz)")
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
   plt.xlabel("Frequency (MHz)")
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
               plt.xlabel("Frequency (MHz)")
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
               #plt.xlabel("Frequency (MHz)")
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
               plt.xlabel("Frequency (MHz)")
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
               plt.xlabel("Frequency (MHz)")
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
            #plt.xlabel("Frequency (MHz)")
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
            #plt.xlabel("Frequency (MHz)")
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
            #plt.xlabel("Frequency (MHz)")
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
            plt.xlabel("Frequency (MHz)")
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
            #plt.xlabel("log Frequency (MHz)")
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
            plt.xlabel("Frequency (MHz)")
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
            plt.xlabel("Frequency (MHz)")
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
            plt.xlabel("Frequency (MHz)")
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
            plt.xlabel("Frequency (MHz)")
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
      plt.xlabel("Frequency (MHz)")
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
            plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='EDA-2',color='#377eb8',marker='.')
            plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_parallel,label='parallel',color='#ff7f00',linestyle='--')
            plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_inline,label='inline',color='#4daf4a',linestyle='-')
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
            fig_name= "x_y_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol_%s.png" % (freq_MHz_fine_chan,pol,EDA2_obs_time)
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
               print(parameters)
            
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
   lst_deg = (float(lst_hrs)/24.)*360.
   #print(freq_MHz_list)
   print(freq_MHz_index)
   start_freq = freq_MHz_list[0]
   freq_MHz = freq_MHz_list[freq_MHz_index]
   #print(freq_MHz)
   centre_freq = float(freq_MHz)
   fine_chan_width_MHz = fine_chan_width_Hz/1000000.
   
   if EDA2_data:
      n_edge_chans_omitted = 5 #two at start and 3 at end
      n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
   else:
      n_fine_chans_used = 1
      n_edge_chans_omitted = 0
   
   bandwidth = (n_fine_chans_used + 1) * fine_chan_width_MHz
   
   #maybe it is actually saved in the correct order now...., so don't reverse here?
   #this makes no difference .. why not ?
   if EDA2_data:  #hack 2
      if not reverse_fine_chans:
         freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
      else:
         #try reversed
         freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
         #freq_MHz_fine_chan = centre_freq + (bandwidth/2.) - (fine_chan_index)*fine_chan_width_MHz 
   else:
      freq_MHz_fine_chan = centre_freq     
   wavelength = 300./float(freq_MHz_fine_chan)  
   #wavelength = 300./float(centre_freq)  
   
   #how does (synthesised) beam solid angle come into this?
   jy_to_K = (wavelength**2) / (2. * k * 1.0e26)   # * 6*PI? (or 4 pi and then fiddle X again?) There are also the visibility weights that I have ignored .... a factor of 120.8
   print("jy_to_K %.4E Jy" % jy_to_K)
   
   if edge_chan:
      return(np.nan,np.nan,np.nan,np.nan,np.nan,freq_MHz_fine_chan) 
   
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
      sky_averaged_diffuse_array_beam_lsts_filename = "%s/sky_av_input_cal_%s_LST_%0.3f_%0.3f_MHz_pol_%s.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan,pol)
      #sky_averaged_diffuse_array_beam_lsts_filename = "%s%s_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,concat_output_name_base)
      #sky_averaged_diffuse_array_beam_lsts_filename =  "woden_map_start_freq_%0.3f_hpx_%s_%s_%s_%s_%s_%s_pol_%s_global_foreground.npy" % (start_freq,year,month,day,hour,min,sec,pol)
   else:
      #sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % (concat_output_name_base)
      sky_averaged_diffuse_array_beam_lsts_filename =  "woden_map_start_freq_%0.3f_hpx_%s_%s_%s_%s_%s_%s_pol_%s_global_foreground.npy" % (start_freq,year,month,day,hour,min,sec,pol)
   #sky_averaged_diffuse_array_no_beam_lsts_filename = "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base
   #freq_MHz_index = int(freq_MHz - 50)
   diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   print("loaded %s" % sky_averaged_diffuse_array_beam_lsts_filename)
   if EDA2_data==True:
      #just doing one eda freq at a time for now
      ##
      diffuse_global_value = diffuse_global_value_array[0]
   else:
      #test
      #diffuse_global_value = diffuse_global_value_array[0]
      diffuse_global_value = diffuse_global_value_array[freq_MHz_index]
   
   #in here put bit to read X from miriad_sim_uvfits
   if not fast:
      if EDA2_data:
         X_short_parallel_array_filename = "X_short_parallel_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         if include_angular_info:
            Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
         baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
      else:
         X_short_parallel_array_filename = "X_short_parallel_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         if include_angular_info:
            Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
         real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         #Bad temporary hack just to get figure done for paper:
         #X_short_parallel_array_filename = "X_short_parallel_array_%s_MHz_%s_pol.npy" % (int(freq_MHz_fine_chan),pol)
         #X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_%s_MHz_%s_pol.npy" % (int(freq_MHz_fine_chan),pol)
         #X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_%s_MHz_%s_pol.npy" % (int(freq_MHz_fine_chan),pol)
         #if include_angular_info:
         #   Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_%s_MHz_%s_pol%s.npy" % (int(freq_MHz_fine_chan),pol,signal_type_postfix)
         #real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%s_MHz_%s_pol%s.npy" % (int(freq_MHz_fine_chan),pol,signal_type_postfix)
         #baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%s_MHz_%s_pol%s.npy" % (int(freq_MHz_fine_chan),pol,signal_type_postfix)
   else:
      X_short_parallel_array_filename = "unity_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
      real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
      sim_vis_data_sorted_array_filename = "sim_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
      baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)  
      if include_angular_info:
         Y_short_parallel_angular_array_filename = "angular_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
                  
   X_short_parallel_array = np.load(X_short_parallel_array_filename)
   print("loaded %s" % X_short_parallel_array_filename)
   
   real_vis_data_sorted_array = np.load(real_vis_data_sorted_array_filename).real
   print("loaded %s" % real_vis_data_sorted_array_filename)

   sim_vis_data_sorted_array = np.load(sim_vis_data_sorted_array_filename).real
   print("loaded %s" % sim_vis_data_sorted_array_filename)

   baseline_length_array_lambda_sorted_cut = np.load(baseline_length_array_lambda_sorted_cut_filename)
   print("loaded %s" % baseline_length_array_lambda_sorted_cut_filename)
   
   if fast:
      X_short_parallel_array = np.concatenate(X_short_parallel_array)
      #X_short_parallel_array = X_short_parallel_array * jy_to_K
      real_vis_data_sorted_array = np.concatenate(real_vis_data_sorted_array)
      sim_vis_data_sorted_array = np.concatenate(sim_vis_data_sorted_array)
      
   if EDA2_data==True:
      real_or_simulated_string = "EDA2"
   else:
      real_or_simulated_string = "simulated"

   if include_angular_info:
      Y_short_parallel_angular_array = np.load(Y_short_parallel_angular_array_filename).real
      print("loaded %s" % Y_short_parallel_angular_array_filename)
      
      if fast:
         Y_short_parallel_angular_array = np.concatenate(Y_short_parallel_angular_array)
         #   Y_short_parallel_angular_array = Y_short_parallel_angular_array

      #plot a histogram of Y values
      plt.clf()
      n, bins, patches = plt.hist(Y_short_parallel_angular_array)
      map_title="Histogram of Y values (angular response)" 
      fig_name= "hist_Y_angular_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name)
      
               
   if not fast:
      X_short_parallel_array_pure_parallel = np.load(X_short_parallel_array_filename_pure_parallel).real
      print("loaded %s" % X_short_parallel_array_filename_pure_parallel) 
      X_short_parallel_array_pure_inline = np.load(X_short_parallel_array_filename_pure_inline).real
      print("loaded %s" % X_short_parallel_array_filename_pure_inline)    
   
  
   

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
         #for fig4a:
         real_vis_data_sorted_array_norm_scaled = X_short_parallel_array_norm * (-0.4)
         
         #plot X and pure inline and parallel for fig 1 of paper
         
         #print(len(X_short_parallel_array_norm))
         #print(len(baseline_length_array_lambda_sorted_cut))
         
         #This is paper 1, fig3, at 70 MHz, /md0/EoR/ASSASSIN/solve_for_tsky_weighted/global_EDGES/x_pol
         #thresh 2.0 lambda
         #make sure using new numbering sys for freq e.g. 70.000_MHz not 70_MHz
         ## plot X and real vis vs baseline length
         print(baseline_length_array_lambda_sorted_cut.shape)
         print(X_short_parallel_array_norm.shape)
         plt.clf()
         plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='EDA-2',color=color_light_blue,marker='.')
         plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_parallel,label='parallel',color=color_orange,linestyle='--')
         plt.plot(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm_pure_inline,label='inline',color=color_green,linestyle='-')
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
      
      
      
         ## plot X and real vis vs baseline length for fig4a, paper 1
         
         plt.clf()
         #just plot evry 5th point [::5]
         plt.scatter(baseline_length_array_lambda_sorted_cut[::10],X_short_parallel_array_norm[::10],label='Global response (unity sky)',color=color_dark_blue,marker='s',s=3)
         plt.scatter(baseline_length_array_lambda_sorted_cut[::10],real_vis_data_sorted_array_norm_scaled[::10],label='Visibility amplitude (simulations)',color=color_orange_red,marker='>',s=3) #% real_or_simulated_string)
         #plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
         map_title="Response to uniform sky vs baseline length data" 
         plt.xlabel("Baseline length (wavelengths)")
         plt.ylabel("Visibility amplitude")
         plt.legend(loc=1)
         #plt.ylim([0, 20])
         fig_name= "X_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name) #,dpi=1000 for paper fig4a
         print("saved %s" % fig_name) 
         
   #also do for fast   
   X_short_parallel_array_diffuse_Jy = (diffuse_global_value * X_short_parallel_array) / jy_to_K
   
   if include_angular_info:
      #Y_short_parallel_array_norm = Y_short_parallel_angular_array / X_short_parallel_array_max_pure_inline
   
      #full response
      #need to convert between Jy and K
   
      Y_short_parallel_angular_array_Jy = Y_short_parallel_angular_array / jy_to_K

      #need to update full response to include fine chans
      full_response_Jy =  X_short_parallel_array_diffuse_Jy + Y_short_parallel_angular_array_Jy
   
   
   
   
   #also include Y and the sum of X plus Y
   plt.clf()
   #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
   #fist attempt:
   #plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array,s=1,label='%s visibility amplitude' % real_or_simulated_string)
   plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array,label='%s visibility amplitude' % real_or_simulated_string,color=color_yellow,marker='o',s=3)
   #plt.scatter(baseline_length_array_lambda_sorted_cut,Y_short_parallel_array_norm,s=1,label='Expected angular response')
   
   #need to update update full response to include fine chans
   plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_diffuse_Jy,label='Expected uniform diffuse response Jy',color=color_dark_blue,marker='s',s=3)
   if include_angular_info:
      plt.scatter(baseline_length_array_lambda_sorted_cut,full_response_Jy,label='Expected full response Jy',color=color_orange_red,marker='>',s=3)
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
   
   
   
   #fig13 paper1
   #Repeat in K
   #also include Y and the sum of X plus Y
   real_vis_data_sorted_array_K = real_vis_data_sorted_array * jy_to_K
   X_short_parallel_array_diffuse_Jy_K =  X_short_parallel_array_diffuse_Jy * jy_to_K
   if include_angular_info:
      full_response_Jy_K = full_response_Jy * jy_to_K
   
   plt.clf()
   #plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_norm,s=1,label='Expected uniform sky response')
   plt.scatter(baseline_length_array_lambda_sorted_cut,real_vis_data_sorted_array_K,label='%s visibility amplitude' % real_or_simulated_string,color=color_yellow,marker='o',s=3)
   #plt.scatter(baseline_length_array_lambda_sorted_cut,Y_short_parallel_array_norm,s=1,label='Expected angular response')
   
   #need to update update full response to include fine chans
   plt.scatter(baseline_length_array_lambda_sorted_cut,X_short_parallel_array_diffuse_Jy_K,label='Expected uniform diffuse response',color=color_dark_blue,marker='s',s=3)
   if include_angular_info:
      plt.scatter(baseline_length_array_lambda_sorted_cut,full_response_Jy_K,label='Expected full response',color=color_orange_red,marker='>',s=3)
   ##plt.plot(n_ants_array,expected_residuals,label='sqrt(n_arrays)',linestyle=':')
   map_title="Response to uniform sky vs baseline length data" 
   plt.xlabel("Baseline length (wavelengths)")
   plt.ylabel("Visibility amplitude (K)")
   plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "X_Y_and_real_vis_vs_uv_dist_%0.3f_MHz_%s_pol%s_K.png" % (freq_MHz_fine_chan,pol,signal_type_postfix)
   figmap = plt.gcf()
   # fig13:
   #figmap.savefig(fig_name,dpi=1000)
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
   data = {'X_global':X_short_parallel_array,'real_vis':real_vis_data_sorted_array,'sim_vis':sim_vis_data_sorted_array}
   
   print(X_short_parallel_array.shape)
   print(real_vis_data_sorted_array.shape)
   print(sim_vis_data_sorted_array.shape)
   
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
   

   if (np.nansum(np.abs(X_short_parallel_array)) > 0 and np.nansum(np.abs(real_vis_data_sorted_array)) > 0 and np.nansum(np.abs(sim_vis_data_sorted_array)) > 0):
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
         
         #repeat for sim
         model_sim = sm.OLS(sim_vis_data_sorted_array, X_short_parallel_array,missing='drop')
         results_sim = model_sim.fit()
         parameters_sim = results_sim.params
         #print parameters
         t_sky_sim_jy = parameters_sim[0]
         t_sky_sim_error_jy = results_sim.bse[0]         
         
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
         if fast and woden:
            real_vis_data_sorted_array_subtr_Y = real_vis_data_sorted_array - Y_short_parallel_angular_array
         elif fast and EDA2_data:
            real_vis_data_sorted_array_subtr_Y = real_vis_data_sorted_array - Y_short_parallel_angular_array
         else:
            real_vis_data_sorted_array_subtr_Y = real_vis_data_sorted_array - Y_short_parallel_angular_array_Jy
         real_vis_data_sorted_array_subtr_Y_K = real_vis_data_sorted_array_subtr_Y * jy_to_K
         model = sm.OLS(real_vis_data_sorted_array_subtr_Y, X_short_parallel_array,missing='drop')
         results = model.fit()
         ##print results.summary()
         parameters = results.params
         #print parameters
         t_sky_jy = parameters[0]
         t_sky_error_jy = results.bse[0]


         #repeat for sim
         if fast and EDA2_data:
            sim_vis_data_sorted_array_subtr_Y = sim_vis_data_sorted_array - Y_short_parallel_angular_array
            ##plot the subtr_y and Y_short array
            plt.clf()
            plt.plot(X_short_parallel_array, sim_vis_data_sorted_array_subtr_Y,label='%s Y' % 'sim',linestyle='None',marker='.')
            ##plt.plot(X_short_parallel_array_nonans, sim_vis_data_sorted_array_nonans,label='%s data' % 'sim',linestyle='None',marker='.')
            ##plt.plot(X_short_parallel_array_nonans, results_sim.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')    
         
            map_title="Data and fit" 
            plt.xlabel("Expected global-signal response")
            plt.ylabel("Real comp. Y (Jy)")
            plt.legend(loc=1)
            ##plt.text(x_pos_sim, y_pos_sim, fit_string_sim)
            ##plt.ylim([0, 3.5])
            fig_name= "check_Y_x_y_OLS_sim_plot_%0.3f_MHz_%s_pol%s_%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            plt.close()
            print("saved %s" % fig_name)
         else:
            print("can't do sims")
         #sim_vis_data_sorted_array_subtr_Y_K = sim_vis_data_sorted_array_subtr_Y * jy_to_K
         model_sim = sm.OLS(sim_vis_data_sorted_array_subtr_Y, X_short_parallel_array,missing='drop')
         results_sim = model_sim.fit()
         ##print results.summary()
         parameters_sim = results_sim.params
         #print parameters
         t_sky_sim_jy = parameters_sim[0]
         t_sky_sim_error_jy = results_sim.bse[0]   
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
      return(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,freq_MHz_fine_chan)
   
   #print("jy_to_K is %0.4E" % jy_to_K)
   
   t_sky_K =   t_sky_jy #* jy_to_K
   t_sky_error_K =  t_sky_error_jy #* jy_to_K 
   print("t_sky_K is %0.4E +/- %0.04f K" % (t_sky_K,t_sky_error_K))
   fit_string = "y=%0.1fx" % t_sky_jy         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
   fit_string_K = "y=%0.1fx" % t_sky_K
   
   print("diffuse_global_value is %0.4E" % diffuse_global_value) 
   
   ratio_in_out = diffuse_global_value / t_sky_K
   print("ratio between input and output T_sky is %0.4f" % ratio_in_out )

   y_pos = np.max(results.fittedvalues)
   x_pos = 1.2 * np.min(X_short_parallel_array)
   
   fitted_values_K = results.fittedvalues * jy_to_K
   
   y_pos_K = np.max(fitted_values_K)
   x_pos = 1.2 * np.min(X_short_parallel_array)
   
   #repeat for sims
   t_sky_sim_K =   t_sky_sim_jy #* jy_to_K
   t_sky_sim_error_K =  t_sky_sim_error_jy #* jy_to_K 
   print("t_sky_sim_K is %0.4E +/- %0.04f K" % (t_sky_sim_K,t_sky_sim_error_K))
   fit_string_sim = "y=%0.1fx" % t_sky_sim_jy         #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
   fit_string_K_sim = "y=%0.1fx" % t_sky_sim_K
   
   ratio_in_out_sim = diffuse_global_value / t_sky_sim_K
   print("ratio between input and output T_sky_sim is %0.4f" % ratio_in_out_sim )
   
   #if model_type=='OLS_fixed_int_subtr_Y':
   #   sys.exit()
   
   y_pos_sim = np.max(results_sim.fittedvalues)
   x_pos_sim = 1.2 * np.min(X_short_parallel_array)
   
   fitted_values_K_sim = results_sim.fittedvalues * jy_to_K
   
   y_pos_K_sim = np.max(fitted_values_K_sim)
   x_pos_sim = 1.2 * np.min(X_short_parallel_array)
    
   #get rid of nans
   real_vis_data_sorted_array_nonans = real_vis_data_sorted_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
   #sim_vis_data_sorted_array_nonans = sim_vis_data_sorted_array[(np.logical_not(np.isnan(sim_vis_data_sorted_array)))]
   X_short_parallel_array_nonans = X_short_parallel_array[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
   real_vis_data_sorted_array_nonans_K = real_vis_data_sorted_array_nonans * jy_to_K
   #sim_vis_data_sorted_array_nonans_K = sim_vis_data_sorted_array_nonans * jy_to_K
   
   
   #in Jy
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
   
   #in K
   plt.clf()
   if model_type=='OLS_with_intercept':
      plt.plot(X_short_parallel_array_nonans[:,1], real_vis_data_sorted_array_nonans_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans[:,1], fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=='OLS_with_int_min_vis':
      plt.plot(X_bin_centres_array[:,1], real_vis_data_min_in_X_bin_array,label='%s data min in X bin' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_bin_centres_array[:,1], fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=='OLS_fixed_int_min_vis':
      plt.plot(X_bin_centres_array, real_vis_data_min_in_X_bin_array,label='%s data min in X bin' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_bin_centres_array, fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
   elif model_type=="OLS_fixed_int_subtr_Y":
      real_vis_data_sorted_array_subtr_Y_nonans_K = real_vis_data_sorted_array_subtr_Y_K[(np.logical_not(np.isnan(real_vis_data_sorted_array)))]
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_subtr_Y_nonans_K,label='%s data - Y' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
   else:
      plt.scatter(X_short_parallel_array_nonans, real_vis_data_sorted_array_nonans_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array_nonans, fitted_values_K, 'r--.', label="OLS fit",linestyle='--',marker='None')
   

   map_title="Data and fit" 
   plt.xlabel("Expected global-signal response")
   plt.ylabel("Real component of visibility (K)")
   plt.legend(loc=1)
   plt.text(x_pos, y_pos_K, fit_string_K)
   #plt.ylim([0, 3.5])
   fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_K.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name) 
   
   #repeat for sims
   #in Jy
   plt.clf()
   if model_type=="OLS_fixed_int_subtr_Y":
      #sim_vis_data_sorted_array_subtr_Y_nonans = sim_vis_data_sorted_array_subtr_Y[(np.logical_not(np.isnan(sim_vis_data_sorted_array)))]
      plt.plot(X_short_parallel_array, sim_vis_data_sorted_array_subtr_Y,label='%s data - Y' % 'sim',linestyle='None',marker='.')
      plt.plot(X_short_parallel_array, sim_vis_data_sorted_array,label='%s data' % 'sim',linestyle='None',marker='.')
      plt.plot(X_short_parallel_array, results_sim.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   else:
      plt.scatter(X_short_parallel_array, sim_vis_data_sorted_array,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      plt.plot(X_short_parallel_array, results_sim.fittedvalues, 'r--.', label="OLS fit",linestyle='--',marker='None')
   
   #maybe the problem has to do with the nonans stuff

   map_title="Data and fit" 
   plt.xlabel("Expected global-signal response")
   plt.ylabel("Real comp. sim (Jy)")
   plt.legend(loc=1)
   plt.text(x_pos_sim, y_pos_sim, fit_string_sim)
   #plt.ylim([0, 3.5])
   fig_name= "x_y_OLS_sim_plot_%0.3f_MHz_%s_pol%s_%s.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   plt.close()
   print("saved %s" % fig_name) 
   
   #flagging
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
   #only do for eda2 data ATM:
   
   if (X_short_parallel_array_nonans.shape[0]>0):
    if EDA2_data:
      model = sm.OLS(real_vis_data_sorted_array_flagged, X_short_parallel_array_nonans,missing='drop')
      results = model.fit()
      ##print results.summary()
      parameters = results.params
      #print parameters
   
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
      t_sky_K_flagged =  t_sky_jy #* jy_to_K
      t_sky_error_K_flagged = t_sky_error_jy  # * jy_to_K 
      print("t_sky_K_flagged is %0.4E +/- %0.04f K" % (t_sky_K_flagged,t_sky_error_K_flagged))
      #fig6b and 10b:
      #fit_string_K = "y=%0.1fx --> T(70 MHz)=-0.4 K" % t_sky_K_flagged      #t_sky_K=%0.6f K" % (t_sky_jy,t_sky_K)
      #6b
      #fit_string_K = "y=%0.1fx" % t_sky_K_flagged 
      #10b
      fit_string_K = "y=%sx" % int(t_sky_K_flagged)
   
      #plot in K
      real_vis_data_sorted_array_flagged_K = real_vis_data_sorted_array_flagged  #* jy_to_K
      fitted_values_K = results.fittedvalues # * jy_to_K
      #fig10b
      y_pos_K = 1200 # * np.max(fitted_values_K) 
      x_pos_K = 0.45 #1.4 * np.min(X_short_parallel_array)
      #fig6b
      #y_pos_K = -0.03 # * np.max(fitted_values_K) 
      #x_pos_K = 0.44 #1.4 * np.min(X_short_parallel_array)
   
   
      plt.clf()
      #loosley dashed (0, (5, 10)))
      plt.plot(X_short_parallel_array_nonans, real_vis_data_sorted_array_flagged_K,label='%s data' % real_or_simulated_string,linestyle='None',marker='.')
      #10b
      plt.plot(X_short_parallel_array_nonans_nonans, fitted_values_K,color='red', label="OLS fit",linestyle=(0, (3, 3)),marker='None')
      #6b:
      #plt.plot(X_short_parallel_array_nonans_nonans, fitted_values_K,color='red', label="OLS fit",linestyle=(0, (4, 6)),marker='None')
      
 
      #fig10 (and 5b and 6b for sims i think)
      map_title="Flagged data and fit" 
      plt.xlabel("Global response (unity sky)")
      plt.ylabel("Visibility amplitude (K)")
      plt.legend(loc=1)
      plt.text(x_pos_K, y_pos_K,fit_string_K,color='red')
      #plt.ylim([0, 3.5])
      fig_name= "x_y_OLS_plot_%0.3f_MHz_%s_pol%s_%s_flagged_K.png" % (freq_MHz_fine_chan,pol,signal_type_postfix,model_type)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      plt.close()
      print("saved %s" % fig_name) 
    
    else:
      t_sky_K_flagged = np.nan
      t_sky_error_K_flagged = np.nan  
      
   else:
      t_sky_jy = np.nan
      t_sky_error_jy = np.nan
      t_sky_K_flagged = np.nan
      t_sky_error_K_flagged = np.nan
    
   return t_sky_K,t_sky_error_K,t_sky_K_flagged,t_sky_error_K_flagged,t_sky_sim_K,t_sky_sim_error_K,diffuse_global_value,freq_MHz_fine_chan

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
                  if not reverse_fine_chans:
                     freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
                  else:
                     freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
                  
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

def solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,include_angular_info=False,EDA2_data=False, EDA2_obs_time='None',EDA2_chan='None',n_obs_concat=1,wsclean=False,fast=False,calculate_uniform_response=True,woden=True,noise_coupling=True):
   if pol=='X':
      pol_index = 0
   elif pol=='Y':
      pol_index = 1
   else:
      print('pol %s not recognised' % pol)
      sys.exit()
   freq_MHz = freq_MHz_list[freq_MHz_index]
   start_freq = freq_MHz_list[0]
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
       if woden:
          type = "gsm"
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
       if woden:
          type = "EDGES_uniform"
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'
   
   #if EDA2_data==True:
   #   signal_type_postfix = "_EDA2_data"
   #n_baselines_included
   #n_baselines_included = 20
   
   
   n_lsts = len(lst_hrs_list)
   

   #lst_hrs = lst_hrs_list[0]
   lst_hrs = lst_hrs_list[freq_MHz_index]
   lst_deg = (float(lst_hrs)/24.)*360.
   
   #for EDA2
   #n_ants = 256

   n_pix = hp.nside2npix(NSIDE)
   print('n_pix')
   print(n_pix)
   pixel_solid_angle = (4.*np.pi) / n_pix
   #print("pixel_solid_angle is %0.4E" % pixel_solid_angle)
   hpx_index_array = np.arange(0,n_pix,1)
   
   #need to understand this timestep stuff, for some reason EDA2 vis have more rows that expected ....
   #eda2 data includes autos?
   if EDA2_data:
      #n_baselines = n_ants*(n_ants-1) / 2. + 256
      #take the centre chan av temp as input
      sky_averaged_temp_cal_input_filename = "%s/sky_av_input_cal_%s_LST_%0.3f_%0.3f_MHz_pol_%s.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz,pol)
      sky_averaged_temp_cal_input_array = np.load(sky_averaged_temp_cal_input_filename)
      print("loaded %s " % sky_averaged_temp_cal_input_filename)
      beam_weighted_av_sky = sky_averaged_temp_cal_input_array
      print("beam_weighted_av_sky is %E" % beam_weighted_av_sky)
   #else:
      #pass
      #n_baselines = n_ants*(n_ants-1) / 2.
   
   if fast:
      print('doing fast')  
      if EDA2_data:
         #obs_time_list = EDA2_obs_time_list_each_chan[freq_MHz_index]
         #print(obs_time_list)
         #open one to get the number of fine chans
         #uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_list[0])
         uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
         
         
      elif woden:
         if noise_coupling:
            uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
         else:
            uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
         obs_time_list = ['1']
      else:
         uvfits_filename = "%s_LST_%03d_%s_%0.3f_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
         obs_time_list = ['1']
      #read the cal uvfits, extract real vis uu and vv
      print("%s" % uvfits_filename)
      hdulist = fits.open(uvfits_filename)
      #hdulist.info()
      #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
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
    
      #only output a file for the central 27 chans
      if EDA2_data:
         n_edge_chans_omitted = 5 #two at start and 3 at end
         n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
      else:
         n_fine_chans_used = 1
         n_edge_chans_omitted = 0
          
      if EDA2_data:
         #print("EDA2 data. Omitting 1 edge chan each side, %s chans present, %s chans used" % (n_fine_chans,n_fine_chans-2))
         #fine_chan_index_array = range(n_fine_chans)[1:n_fine_chans-1]
         #changing this to test (15 Feb 2021)
         fine_chan_index_array = np.arange(n_fine_chans_used) + 2
         #print(fine_chan_index_array)
      else:
         fine_chan_index_array = np.asarray([0])
      centre_freq = float(freq_MHz)
      fine_chan_width_MHz = fine_chan_width_Hz/1000000.         
      
      
      #now do for each fine chan:
      for fine_chan_index in fine_chan_index_array:
         fine_chan_index = int(fine_chan_index)
         baseline_length_array_lambda_sorted_cut_list = []

         if EDA2_data: #hack 1
            #data coming out of the TPMs is reversed by coarse chan so for 20200303_data (and 20200304), need to change the freq calculation
            #2020 paper:
            if not reverse_fine_chans:
               freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
            else:
               freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               #freq_MHz_fine_chan = centre_freq - (fine_chan_index)*fine_chan_width_MHz
         else:
            freq_MHz_fine_chan = freq_MHz
         wavelength = 300./float(freq_MHz_fine_chan)
         
         print("fine_chan index,MHz,wavelength")
         print(fine_chan_index)
         print(freq_MHz_fine_chan)
         print(wavelength)
         
         unity_vis_data_sorted_list = []
         angular_vis_data_sorted_list = []
         sim_vis_data_sorted_list = []
         baseline_length_array_lambda_sorted_cut_list = []
         real_vis_data_sorted_list = []
         for obs_time_fast in [EDA2_obs_time]:
            if EDA2_data:
               #uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_fast)
               uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
               #this is old pre woden miriad stuff
               unity_uvfits_filename = "%s/unity_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,obs_time_fast,pol)
               angular_uvfits_filename = "%s/angular_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,obs_time_fast,pol)
               gsm_uvfits_filename = "%s/gsm_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,obs_time_fast,pol)
               
               #this is woden - havent run these yet
               #unity_uvfits_filename = "%s/woden_LST_%0.3f_unity_uniform_start_freq_%0.3f_band%02d.uvfits" % (EDA2_chan,lst_deg,start_freq,freq_MHz_index) 
               #angular_uvfits_filename = "%s/woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular_band%02d.uvfits" % (EDA2_chan,lst_deg,start_freq,pol,freq_MHz_index) 
            elif woden:
               if noise_coupling:
                  uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
               else:
                  uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) 
               unity_uvfits_filename = "woden_LST_%0.3f_unity_uniform_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,start_freq,freq_MHz_index) 
               angular_uvfits_filename = "woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular_band%02d.uvfits" % (lst_deg,start_freq,pol,freq_MHz_index) 
            
            #open the unity fits file to get the matching baselines
            with fits.open(unity_uvfits_filename) as hdulist:
               uvtable = hdulist[0].data
               unity_baseline_array = uvtable['BASELINE']
               #print(unity_baseline_array)
               
               
            #read the cal uvfits, extract real vis uu and vv
            print("%s" % uvfits_filename)
            with fits.open(uvfits_filename) as hdulist:
               #hdulist.info()
               #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
               #print info_string
               uvtable = hdulist[0].data
               uvtable_header = hdulist[0].header
               #print(uvtable_header)
               #hdulist.close()

               data_baseline_array = uvtable['BASELINE']
               print(data_baseline_array.shape)
               print(data_baseline_array.min())
               print(data_baseline_array.max())
               #incl autos so b in quad eqn is -3
               n_baselines = data_baseline_array.shape[0]
               n_ants = int((-1 + math.sqrt(1 + 4*2*n_baselines) ) / 2)
               print("n_ants is %s" % n_ants)
               
               full_auto_baseline_number_array = np.arange(1,257)*256 + np.arange(1,257)
               
               ##find the indices where the baseline numbers in the data match the expected full autos (should match n_ants)
               #auto_baselines_in_data_mask = np.in1d(full_auto_baseline_number_array,data_baseline_array)
               #auto_baselines_in_data_count = np.count_nonzero(auto_baselines_in_data_mask)
               #print("auto_baselines_in_data %s" % auto_baselines_in_data_count)
               
               #ant1,ant2 = split_baseline(np.asarray([257,258,259]))
               #print(ant1,ant2)
               #sys.exit()
               
               #find indices of unity file where the baseline value is the same as the data baseline value
               #To find the indices of the elements in A that are present in B, I can do
               unity_common_inds = np.in1d(unity_baseline_array,data_baseline_array)
               count = np.count_nonzero(unity_common_inds)
               print('count of unity_common_inds')
               print(count)
               
               #also get the indices of the data file that are common with the unity file (cause data file has autos and sims dont)
               data_common_inds = np.in1d(data_baseline_array,unity_baseline_array)
               count2 = np.count_nonzero(data_common_inds)
               print('count of unity_common_inds')
               print(count2)

               #get the autocorrelation baseline indices. done
               #coherence = cross_a1_a2 / np.sqrt(auto_a1 * auto_a2)
    
               


               #only use dat where there is a common index with unity file
               visibilities_single = uvtable['DATA'][data_common_inds]
               visibilities_shape = visibilities_single.shape
               print("visibilities_shape")
               print(visibilities_shape)
               
               if wsclean:
                  real_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,0]
                  imag_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,1]
                  weights_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,2]
               else:
                  real_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,0]
                  imag_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,1]
                  weights_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,2]
         
               UU_s_array = uvtable['UU'][data_common_inds]
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV'][data_common_inds]
               VV_m_array = VV_s_array * c
      
               #print(VV_s_array[0:100])
               #print(UU_s_array[0:100])

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
            

            n_baselines_included_data = len(baseline_length_array_lambda_sorted_cut)
            print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included_data,obs_time_fast,fine_chan_index))
            
            baseline_length_array_lambda_sorted_cut_list.append(baseline_length_array_lambda_sorted_cut)
            real_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included_data])

               
            #######################################################################################################
            #now repeat for unity sky to get X_short_parallel!
            #note there are separate uvfits for each fine chan
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
            
            #only use the values where the baselines are common
            visibilities_single = uvtable['DATA'][unity_common_inds]
            visibilities_shape = visibilities_single.shape
            print("visibilities_shape")
            print(visibilities_shape)
            
            
            UU_s_array = uvtable['UU'][unity_common_inds]
            UU_m_array = UU_s_array * c   
            VV_s_array = uvtable['VV'][unity_common_inds]
            VV_m_array = VV_s_array * c
            
            #fine chan index is zero since we have a new uvfits for each fine chan
            unity_fine_chan_index = 0
            #if reverse_fine_chans:
            #   unity_fine_chan_index = 31 - fine_chan_index
            #else:
            #   unity_fine_chan_index = fine_chan_index
            #the uvfits files used to make the unity sky data had reversed fine channel ordering (true for 20200303 and 20200304) - marcin will fix this in later data
            #yes but this should not affect the unity uvfits
            ####fine_chan_index_array = fine_chan_index_array[::-1]
         
            print(visibilities_single.shape)
            #real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
            #TEST!
            real_vis_data = visibilities_single[:,0,0,unity_fine_chan_index,0,0]
                         
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
            print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included,obs_time_fast,unity_fine_chan_index))
            

            real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
 
            #hack just to get it working for a check Feb 5 2021
            #unity_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included])
            unity_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included_data])
  
            #sys.exit()
  
            #######
            #######################################################################################################
            #now repeat for angular sky to get Y_short_parallel!
            print("%s" % angular_uvfits_filename)
            #TEST!
            #unity_uvfits_filename = '/md0/EoR/ASSASSIN/solve_for_tsky_weighted/global_unity/eda_model_LST_030_X_50_MHz_GU.uvfits'
            hdulist = fits.open(angular_uvfits_filename)
            hdulist.info()
            info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
            #print info_string
            uvtable = hdulist[0].data
            uvtable_header = hdulist[0].header
            #print(uvtable_header)
            hdulist.close()
      
            #only use the values where the baselines are common
            visibilities_single = uvtable['DATA'][unity_common_inds]
            visibilities_shape = visibilities_single.shape
            print("visibilities_shape")
            print(visibilities_shape)
            
            UU_s_array = uvtable['UU'][unity_common_inds]
            UU_m_array = UU_s_array * c   
            VV_s_array = uvtable['VV'][unity_common_inds]
            VV_m_array = VV_s_array * c
      

            #the uvfits files used to make the unity sky data had reversed fine channel ordering (true for 20200303 and 20200304) - marcin will fix this in later data
            #yes but this should not affect the unity uvfits
            ####fine_chan_index_array = fine_chan_index_array[::-1]
         

            #real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
            #TEST!
            real_vis_data = visibilities_single[:,0,0,unity_fine_chan_index,pol_index,0]
            
            #print(real_vis_data)
            #sys.exit()  
                         
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
            print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included,obs_time_fast,unity_fine_chan_index))
            
      

            real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
 

            angular_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included_data])

            #######################################################################################################
            #now repeat for full gsm sky to get sim_vis_data_sorted!
            print("%s" % gsm_uvfits_filename)
            #TEST!
            #unity_uvfits_filename = '/md0/EoR/ASSASSIN/solve_for_tsky_weighted/global_unity/eda_model_LST_030_X_50_MHz_GU.uvfits'
            hdulist = fits.open(gsm_uvfits_filename)
            hdulist.info()
            info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
            #print info_string
            uvtable = hdulist[0].data
            uvtable_header = hdulist[0].header
            #print(uvtable_header)
            hdulist.close()
      
            visibilities_single = uvtable['DATA'][unity_common_inds]
            visibilities_shape = visibilities_single.shape
            print("visibilities_shape")
            print(visibilities_shape)
            
            UU_s_array = uvtable['UU'][unity_common_inds]
            UU_m_array = UU_s_array * c   
            VV_s_array = uvtable['VV'][unity_common_inds]
            VV_m_array = VV_s_array * c
      

            #the uvfits files used to make the unity sky data had reversed fine channel ordering (true for 20200303 and 20200304) - marcin will fix this in later data
            #yes but this should not affect the unity uvfits
            ####fine_chan_index_array = fine_chan_index_array[::-1]
         

            #real_vis_data = visibilities_single[:,0,0,fine_chan_index,0,0]
            #TEST!
            real_vis_data = visibilities_single[:,0,0,unity_fine_chan_index,pol_index,0]
            
            #print(real_vis_data)
            #sys.exit()  
                         
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
            print("n_baselines_included %s for obs %s, fine chan %s" % (n_baselines_included,obs_time_fast,unity_fine_chan_index))
            
      

            real_vis_data_sorted = real_vis_data_sorted_orig[UU_m_array_sorted_orig>0]
 

            sim_vis_data_sorted_list.append(real_vis_data_sorted[0:n_baselines_included_data])


         
         #######
            
         #uvfits_real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix,obs_time_fast)
         #uvfits_real_vis_data_sorted = np.load(uvfits_real_vis_data_sorted_array_filename)
         #print("loaded %s" % uvfits_real_vis_data_sorted_array_filename)
         
         
         #baseline_length_array_lambda_sorted_cut_filename_unity = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_unity_%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix,obs_time_fast)
         #np.save(baseline_length_array_lambda_sorted_cut_filename_unity,baseline_length_array_lambda_sorted_cut)
         #print("saved %s" % baseline_length_array_lambda_sorted_cut_filename_unity)
         
         unity_vis_data_sorted_array = np.asarray(unity_vis_data_sorted_list)
         angular_vis_data_sorted_array = np.asarray(angular_vis_data_sorted_list)
         sim_vis_data_sorted_array = np.asarray(sim_vis_data_sorted_list)
         baseline_length_array_lambda_sorted_cut_array = np.asarray(baseline_length_array_lambda_sorted_cut_list)
         real_vis_data_sorted_array = np.asarray(real_vis_data_sorted_list)
         
         
         
         unity_vis_data_sorted_array_filename = "unity_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
         np.save(unity_vis_data_sorted_array_filename,unity_vis_data_sorted_array)
         print(unity_vis_data_sorted_array.shape)
         print(np.max(unity_vis_data_sorted_array))
         print("saved %s" % unity_vis_data_sorted_array_filename)         

         

         angular_vis_data_sorted_array_filename = "angular_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
         np.save(angular_vis_data_sorted_array_filename,angular_vis_data_sorted_array)
         print("saved %s" % angular_vis_data_sorted_array_filename)  

         sim_vis_data_sorted_array_filename = "sim_vis_data_sorted_array_%0.3f_MHz_%s_pol_fast.npy" % (freq_MHz_fine_chan,pol)
         np.save(sim_vis_data_sorted_array_filename,sim_vis_data_sorted_array)
         print("saved %s" % sim_vis_data_sorted_array_filename)  
                  
         baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         np.save(baseline_length_array_lambda_sorted_cut_filename,baseline_length_array_lambda_sorted_cut_array)
         print("saved %s" % baseline_length_array_lambda_sorted_cut_filename)
   
         real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s_fast.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
         np.save(real_vis_data_sorted_array_filename,real_vis_data_sorted_array)
         print(real_vis_data_sorted_array.shape)
         print("saved %s" % real_vis_data_sorted_array_filename)

         
      #what the hell was this else here for?     
      #else:
      #   pass
      
      #   Dont think we need or want any of this stuff - repeated below?
      #if n_obs_concat==1:
      #   if wsclean==True:
      #      uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
      #   elif woden:
      #      if noise_coupling:
      #         uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
      #      else:
      #         uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index)
      #   else:
      #      uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
      #else:
      #   if wsclean==True:
      #      #uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      #      #use calibrated average data
      #      uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      #   else:
      #      #concat_chan_64_20191202T171525_n_obs_13.uvfits
      #      #uvfits_filename = "%s/av_chan_%s_%s_n_obs_%s_t_av_cal_freq_av.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      #      uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
   else:
      print("not doing fast")
      if EDA2_data:
         if n_obs_concat==1:
            if wsclean==True:
               uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            else:
               uvfits_filename = "%s/cal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
         else:
            if wsclean==True:
               #uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
               #use calibrated average data
               uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            else:
               uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
      else:
         if woden:
            if noise_coupling:
               uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
            else: 
               uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
         else:
            uvfits_filename = "%s_LST_%03d_%s_%0.3f_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
   
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

   #only output a file for the central 27 chans
   if EDA2_data:
      n_edge_chans_omitted = 5 #two at start and 3 at end
      n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
   else:
      n_fine_chans_used = 1
      n_edge_chans_omitted = 0
         
   uvfits_filename_list = []
   if EDA2_data:
      if n_obs_concat==1:
         if wsclean==True:
            #uvfits_filename = "%s/chan_%s_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
            #use calibrated average data
            uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            #ms_filename = "%s/chan_%s_%s_ws_cal.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
         else:
            uvfits_filename = "%s/chan_%s_%s_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
      else:
         if wsclean==True:
            #uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat) 
            #ms_filename = "%s/concat_chan_%s_%s_n_obs_%s_ws_cal.ms" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat) 
            #use calibrated average data
            uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
         else:
            #uvfits_filename = "%s/av_chan_%s_%s_n_obs_%s_t_av_cal_freq_av.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
            uvfits_filename = "%s/concat_chan_%s_%s_n_obs_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)      
      uvfits_filename_list = [uvfits_filename]

   else:
      for lst_hrs in lst_hrs_list:
         lst_deg = (float(lst_hrs)/24.)*360.
         if woden:
            if noise_coupling:
               uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
            else:
               uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) 
         else:
            uvfits_filename = "%s_LST_%03d_%s_%0.3f_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
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
         fine_chan_index_array = np.arange(n_fine_chans_used) + 2
         #print(fine_chan_index_array)
      else:
         fine_chan_index_array = np.asarray([0])
      centre_freq = float(freq_MHz)
      fine_chan_width_MHz = fine_chan_width_Hz/1000000.
      for fine_chan_index in fine_chan_index_array:
         fine_chan_index = int(fine_chan_index)
         if EDA2_data:
            #for 20200303 and 20200304 data fine chan order is reversed?
            if not reverse_fine_chans:
               freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
            else:
               freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz
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
            real_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,0]
            imag_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,1]
            weights_vis_data = visibilities_single[:,0,0,0,fine_chan_index,pol_index,2]
         else:
            real_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,0]
            imag_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,1]
            weights_vis_data = visibilities_single[:,0,0,fine_chan_index,pol_index,2]
         
         
         
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
         

            
         #eda2 data have autos (where uu=vv=0), dont use these
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
            if EDA2_data:
               baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
            else:
               baseline_length_array_lambda_sorted_cut_filename = "baseline_length_array_lambda_sorted_cut_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
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
         #lst_initial = time_initial.sidereal_time('apparent')
         
         #lst_initial_days = (lst_initial.value / 24.) * (23.9344696 /24.)
         
         ##want lst == lst_hrs, so add time so that this is true
         #delta_time = TimeDelta(lst_initial_days, format='jd') 
         
         #desired_lst_days = float(lst_hrs) / 24.
         
         #time_final = time_initial - delta_time + TimeDelta(desired_lst_days, format='jd') 
         
         #print('final LST is: ')
         #print(time_final.sidereal_time('apparent'))
         #
         #print('final time is: ')
         #print(time_final)
         
         #close enough......
         
         #time_string = time_final.utc.iso
         #time_string = "%02d_%02d_%02d" % (hour,minute,second)
         
         #hour = int(time_string.split(' ')[1].split(':')[0])
         #minute = int(time_string.split(' ')[1].split(':')[1])
         #minute = int(np.floor(float(time_string.split(' ')[1].split(':')[2])))
         
         #date_obs = datetime(year, month, day, hour, minute, second)
         print("WARNING: using hard-coded time of date_obs = datetime(2015, 11, 29, 15, 33, 43) to match woden sims for testing angular subtraction" )
         year = 2015
         month = 11
         day = 29
         hour = 15
         min = 33
         sec = 43
         date_obs = datetime(year, month, day, hour, min, sec)
         ov.date = date_obs
         
         
         gsm_map = ov.generate(freq_MHz) * 0.
         
         #get the diffuse global diffuse value used in the simulation (from gsm)
         if EDA2_data==True:
            EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)          
            sky_averaged_diffuse_array_beam_lsts_filename = "%s%s_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,concat_output_name_base)       
            #sky_averaged_diffuse_array_beam_lsts_filename =  "woden_map_start_freq_%0.3f_hpx_%s_%s_%s_%s_%s_%s_pol_%s_global_foreground.npy" % (start_freq,year,month,day,hour,min,sec,pol)
         else:
            #sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % (concat_output_name_base)
            sky_averaged_diffuse_array_beam_lsts_filename =  "woden_map_start_freq_%0.3f_hpx_%s_%s_%s_%s_%s_%s_pol_%s_global_foreground.npy" % (start_freq,year,month,day,hour,min,sec,pol)
         #sky_averaged_diffuse_array_no_beam_lsts_filename = "%s_sky_averaged_diffuse_no_beam.npy" % concat_output_name_base
         #freq_MHz_index = int(freq_MHz - 50)
         diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
         if EDA2_data==True:
            #just doing one freq at a time right now for EDA2, not sure how this works with fine chans
            diffuse_global_value = diffuse_global_value_array[0]
         else:
            #diffuse_global_value = diffuse_global_value_array[0]
            diffuse_global_value = diffuse_global_value_array[freq_MHz_index]
         
         #ov1 = GlobalSkyModel()
         
         if 'global' in signal_type_list:
            if 'global_EDGES' in signal_type_list:
               print("cant have global and global edges in signal_type_list")
               sys.exit()
            else:
               freq_MHz_array = np.asarray([freq_MHz_fine_chan])
               s_21_array = plot_S21(nu_array=freq_MHz_array)
               global_signal_value = s_21_array[0]
               #global_signal_value = s_21_array[freq_MHz_index]
               gsm_map += global_signal_value
         if 'global_EDGES' in signal_type_list:
            if 'global' in signal_type_list:
               print("cant have global and global edges in signal_type_list")
               sys.exit()
            else:
               freq_MHz_array = np.asarray([freq_MHz_fine_chan])
               s_21_array_EDGES = plot_S21_EDGES(nu_array=freq_MHz_array)
               global_signal_value = s_21_array_EDGES[0]
               ##global_signal_value = s_21_array_EDGES[freq_MHz_index]
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
            
            #d = ov.view(logged=True)
            #fig_name="testing_1999-12-31.png"
            #figmap = plt.gcf()
            #figmap.savefig(fig_name,dpi=500)
            #print("saved %s" % fig_name)
            
            #2015-11-29T15:33:43
            #date_obs = datetime(2015, 11, 29, 15, 33, 43)
            #ov.date = date_obs
            #ov.generate(freq_MHz)
            #d = ov.view(logged=True)
            #fig_name="testing_2015-11-29.png"
            #figmap = plt.gcf()
            #figmap.savefig(fig_name,dpi=500)
            #print("saved %s" % fig_name)
            #sys.exit()
                  
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
            fig_name="check_%s_%s_hrs_LST_%0.3f_MHz_%s_pol.png" % (sky_model,lst_hrs,freq_MHz,pol)
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
         if not EDA2_data:
            beam_weighted_av_sky = np.nansum(sky_with_beam) /  sum_of_beam_weights  #(2.*np.pi/float(n_pix)) #
            print("beam_weighted_av_sky at %0.3f MHz is %0.4E for %s pol" % (freq_MHz,beam_weighted_av_sky,pol))
         
         
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
               X_short_parallel_array_filename = "X_short_parallel_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
               np.save(X_short_parallel_array_filename,X_short_parallel_array)
            
               print("saved %s" % X_short_parallel_array_filename)
            
               X_short_parallel_array_filename_pure_inline = "X_short_parallel_array_pure_inline_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
               np.save(X_short_parallel_array_filename_pure_inline,X_short_parallel_array_pure_inline)
               print("saved %s" % X_short_parallel_array_filename_pure_inline)
                
               X_short_parallel_array_filename_pure_parallel = "X_short_parallel_array_pure_parallel_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
               np.save(X_short_parallel_array_filename_pure_parallel,X_short_parallel_array_pure_parallel)
               print("saved %s" % X_short_parallel_array_filename_pure_parallel)
               
               #update for fine chans
               if include_angular_info:
                  Y_short_parallel_angular_array_filename = "Y_short_parallel_angular_array_%0.3f_MHz_%s_pol.npy" % (freq_MHz_fine_chan,pol)
                  np.save(Y_short_parallel_angular_array_filename,Y_short_parallel_angular_array)
                  print("saved %s" % Y_short_parallel_angular_array_filename)
             
            if EDA2_data:   
               real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_chan_%s_%0.3f_MHz_%s_pol%s.npy" % (EDA2_chan,freq_MHz_fine_chan,pol,signal_type_postfix)
            else:
               real_vis_data_sorted_array_filename = "real_vis_data_sorted_array_%0.3f_MHz_%s_pol%s.npy" % (freq_MHz_fine_chan,pol,signal_type_postfix)
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
    
   return(n_baselines_included)
   
         
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

def joint_model_fit_t_sky_measured(lst_hrs_list,freq_MHz_list,pol_list,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,poly_order_list,global_signal_model,plot_only=False,model_type='OLS_fixed_intercept'):

   poly_order_list_string = '_'.join(str(e) for e in poly_order_list)
   pol = pol_list[0]
   freq_MHz_array = np.asarray(freq_MHz_list)
   
   linestyle_list = ['--','-.',':']
   
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


   t_sky_measured_array_filename = "t_sky_measured_array_lsts_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
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
   #fig8a paper1
   plt.clf()
   for joint_fit_global_EDGES_index,joint_fit_global_EDGES in enumerate(joint_fit_global_EDGES_list):
      plt.plot(freq_MHz_array_okay,joint_fit_global_EDGES,label='recovered order %s' % poly_order_list[joint_fit_global_EDGES_index],linestyle=linestyle_list[joint_fit_global_EDGES_index])
   plt.plot(freq_MHz_array_okay,global_signal_model_okay,label='EDGES input',linestyle='-',alpha=0.7)
   map_title="t_sky gobal EDGES joint fit" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   plt.legend(loc="lower right")
   #plt.ylim([0, 20])
   fig_name= "t_sky_joint_fit_global_EDGES_LST_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order_list_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
   plt.clf()
   for joint_fit_global_EDGES_index,joint_fit_global_EDGES in enumerate(joint_fit_global_EDGES_list):
      plt.plot(freq_MHz_array_okay,joint_fit_diffuse_global_foreground,label='Diffuse poly order %s' % poly_order_list[joint_fit_global_EDGES_index],linestyle=linestyle_list[joint_fit_global_EDGES_index])
   map_title="t_sky diffuse joint fit" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "t_sky_joint_fit_diffuse_global_foreground_LST_%s%s_order_%s.png" % (lst_string,signal_type_postfix,poly_order_list_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
   #fig8b paper 1
   #residuals from fit:\
   plt.clf()
   for joint_fit_index,joint_fit in enumerate(joint_fit_list):
      residuals = t_sky_measured_array_okay - joint_fit
      rms_of_residuals = np.sqrt(np.mean(residuals**2))
      max_abs_residuals = np.max(np.abs(residuals)) * 0.9

      #plot residuals
      plt.plot(freq_MHz_array_okay,residuals,label='Residuals poly order %s' % poly_order_list[joint_fit_index],linestyle=linestyle_list[joint_fit_index])
      
   map_title="residuals from joint fitting" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Residual Tb (K)")
   plt.legend(loc=1)
   plt.text(55, max_abs_residuals, "rms=%0.3f K" % rms_of_residuals,color='green')
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


def plot_tsky_for_multiple_freqs(lst_hrs_list,freq_MHz_list,pol_list,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,poly_order,plot_only=False,include_angular_info=False,model_type_list=['OLS_fixed_intercept'],EDA2_data=False,EDA2_chan_list='None',n_obs_concat_list=[],wsclean=False,fast=False,no_modelling=False,calculate_uniform_response=True,woden=True,noise_coupling=True):
   #for plot_expected_rms_noise_eda2 below
   int_time = 4.*60. #0.28 * 5.
   bw_Hz = 1000000 #27. * fine_chan_width_Hz

   #pol = pol_list[0]
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
   
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   
   
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
   if 'global_unity' in signal_type_list:
       signal_type_postfix += '_GU'
   if 'diffuse_global' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse global at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DG_%s' % sky_model
   if 'diffuse_angular' in signal_type_list:
       if 'diffuse' in signal_type_list:
          print("cant have diffuse and diffuse angular at same time")
          sys.exit()
       else:
          signal_type_postfix += '_DA_%s' % sky_model
   if 'global' in signal_type_list:
       if 'global_EDGES' in signal_type_list:
          print('cant have global_EDGES and global in signal_type_list')
          sys.exit()
       else:
          signal_type_postfix += '_G' 
   if 'global_EDGES' in signal_type_list:
       signal_type_postfix += '_ED'  
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'

   if EDA2_data==True:
      signal_type_postfix = "_EDA2_data"
   

   if not plot_only:
      for pol in pol_list:
         #t_sky_theoretical_array = np.full(len(freq_MHz_list),np.nan)
         n_baselines_used_array = np.full(len(freq_MHz_list),np.nan)
         n_baselines_used_array_filename = "n_baselines_included_lst_%s_pol_%s%s.npy" % (lst_string,pol,signal_type_postfix)
         #t_sky_theoretical_array_filename = "t_sky_theoretical_array_lst_%s_pol_%s%s.npy" % (lst_string,pol,signal_type_postfix)

         print(EDA2_obs_time_list)
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
                  n_baselines_used = solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,include_angular_info=include_angular_info,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,wsclean=wsclean,fast=fast,calculate_uniform_response=calculate_uniform_response,woden=woden,noise_coupling=noise_coupling)
               else:
                  t_sky_theoretical = extract_data_from_eda2_uvfits(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs_list=lst_hrs_list,pol=pol,EDA2_chan=EDA2_chan,n_obs=n_obs_concat,calculate_uniform_response=calculate_uniform_response)
            else:
               t_sky_theoretical,n_baselines_used = solve_for_tsky_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,include_angular_info=include_angular_info,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,wsclean=wsclean,fast=fast,calculate_uniform_response=calculate_uniform_response,woden=woden,noise_coupling=noise_coupling)
            #t_sky_measured_array[freq_MHz_index] = t_sky_measured
            #t_sky_measured_error_array[freq_MHz_index] = t_sky_measured_error
            #t_sky_theoretical_array[freq_MHz_index] = t_sky_theoretical
            n_baselines_used_array[freq_MHz_index] = n_baselines_used  

         #np.save(t_sky_theoretical_array_filename,t_sky_theoretical_array)
         np.save(n_baselines_used_array_filename,n_baselines_used_array)

   
   t_sky_array_length = int(len(freq_MHz_list) * n_fine_chans_used)
   
   freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   
   

   freq_MHz_fine_array = np.full(t_sky_array_length,np.nan)
   
   if not no_modelling:
      for pol_index,pol in enumerate(pol_list):
         t_sky_beam_wtd_av_array_filename = "t_sky_beam_wtd_av_lst_%s_pol_%s.npy" % (lst_string,pol)
         t_sky_beam_wtd_av_array = np.full(t_sky_array_length,np.nan)
         
         t_sky_measured_array = np.full(t_sky_array_length,np.nan)
         t_sky_measured_error_array = np.full(t_sky_array_length,np.nan)
         t_sky_sim_measured_array = np.full(t_sky_array_length,np.nan)
         t_sky_sim_measured_error_array = np.full(t_sky_array_length,np.nan)
         t_sky_measured_array_flagged = np.full(t_sky_array_length,np.nan)
         t_sky_measured_error_array_flagged = np.full(t_sky_array_length,np.nan)
         


         for model_type in model_type_list:
            t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
            t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
            t_sky_sim_measured_array_filename = "t_sky_sim_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
            t_sky_sim_measured_error_array_filename = "t_sky_sim_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
            t_sky_measured_array_filename_flagged = "t_sky_measured_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
            t_sky_measured_error_array_filename_flagged = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
            

            #this replaces all the matrix stuff you do in model_tsky_from_saved_data
            
            for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
               if EDA2_data==True:
                  EDA2_chan = EDA2_chan_list[freq_MHz_index]
                  lst_hrs = lst_hrs_list[freq_MHz_index]
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
                     #what if you change this to fine_chan_index?
                     #freq_MHz_index_fine = freq_MHz_index*n_fine_chans_used + fine_chan_index_index
                     freq_MHz_index_fine = freq_MHz_index*n_fine_chans_used + fine_chan_index
                     
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
                           t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,t_sky_sim_measured,t_sky_sim_measured_error,beam_wtd_av,freq_MHz_fine = model_tsky_from_saved_data(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs=lst_hrs,pol=pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,model_type=model_type,EDA2_data=EDA2_data,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,fine_chan_index=fine_chan_index,edge_chan=edge_chan,wsclean=wsclean,fast=fast)
                        else:
                           t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data_eda2(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs_list=lst_hrs_list,pol=pol,EDA2_chan=EDA2_chan,n_obs=n_obs_concat,fine_chan_index=fine_chan_index,model_type=model_type)
                     else:
                        t_sky_measured,t_sky_measured_error,t_sky_measured_flagged,t_sky_measured_error_flagged,freq_MHz_fine = model_tsky_from_saved_data(freq_MHz_list=freq_MHz_list,freq_MHz_index=freq_MHz_index,lst_hrs=lst_hrs,pol=pol,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,model_type=model_type,EDA2_data=EDA2_data,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,fine_chan_index=fine_chan_index,edge_chan=edge_chan,wsclean=wsclean,fast=fast)
                     
                     print(freq_MHz_fine)
                     print(freq_MHz_index_fine)
                     
                     
                     t_sky_measured_array[freq_MHz_index_fine-2] = t_sky_measured
                     t_sky_measured_error_array[freq_MHz_index_fine-2] = t_sky_measured_error
                     t_sky_sim_measured_array[freq_MHz_index_fine-2] = t_sky_sim_measured
                     t_sky_sim_measured_error_array[freq_MHz_index_fine-2] = t_sky_sim_measured_error
                     t_sky_measured_array_flagged[freq_MHz_index_fine-2] = t_sky_measured_flagged
                     t_sky_measured_error_array_flagged[freq_MHz_index_fine-2] = t_sky_measured_error_flagged
                     t_sky_beam_wtd_av_array[freq_MHz_index_fine-2] = beam_wtd_av
                     
                     freq_MHz_fine_array[freq_MHz_index_fine-2] = freq_MHz_fine
                        
            
                  
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
            np.save(t_sky_sim_measured_array_filename,t_sky_sim_measured_array)
            np.save(t_sky_sim_measured_error_array_filename,t_sky_sim_measured_error_array)             
            np.save(t_sky_measured_array_filename_flagged,t_sky_measured_array_flagged)
            np.save(t_sky_measured_error_array_filename_flagged,t_sky_measured_error_array_flagged)      
         
         np.save(t_sky_beam_wtd_av_array_filename,t_sky_beam_wtd_av_array)    
         np.save(freq_MHz_fine_array_filename,freq_MHz_fine_array)
         
         #print(freq_MHz_fine_array)
         #sys.exit()
         
   #make a plot of diffuse global input  
   #dont need this at the moment
   #if 'diffuse_global' in signal_type_list:
   #   sky_averaged_diffuse_array_beam_lsts_filename = "%s_sky_averaged_diffuse_beam.npy" % concat_output_name_base
   #   diffuse_global_value_array_X = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   #   #diffuse_global_value_array_Y = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
   #   
   #   plt.clf()
   #   plt.plot(freq_MHz_array,diffuse_global_value_array_X,label='sim input X')
   #   #plt.plot(freq_MHz_array,diffuse_global_value_array_Y,label='sim input Y')
   #   map_title="t_sky beam averaged input" 
   #   plt.xlabel("Frequency (MHz)")
   #   plt.ylabel("t_sky (K)")
   #   plt.legend(loc=1)
   #   #plt.ylim([0, 20])
   #   fig_name= "t_sky_beam_av_input_%s%s.png" % (lst_string,signal_type_postfix)
   #   figmap = plt.gcf()
   #   figmap.savefig(fig_name)
   #   print("saved %s" % fig_name)
   # 
   #   #subtract a polynomial fit just do for X for now
   #   #in log log space:
   #   sky_array = t_sky_theoretical_array[t_sky_measured_array>0.]
   #   log_sky_array = np.log10(sky_array)
   #   freq_array_cut = freq_MHz_array[t_sky_measured_array>0.]
   #   log_freq_MHz_array = np.log10(freq_array_cut)
   #   coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
   #   ffit = poly.polyval(log_freq_MHz_array, coefs)
   #   ffit_linear = 10**ffit
   #   
   #   #log_residual = log_signal_array_short_baselines - log_ffit
   #   residual_of_log_fit = ffit_linear - sky_array
   #   
   #   rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
   #   print("rms_of_residuals for X is %0.3f K" % rms_of_residuals)
   #   
   #   max_abs_residuals = np.max(np.abs(residual_of_log_fit))
   #   y_max = 1.5 * max_abs_residuals
   #   y_min = 1.5 * -max_abs_residuals
   #   
   #   plt.clf()
   #   plt.plot(freq_array_cut,residual_of_log_fit,label='residual of log fit X')
   #   map_title="Residual for log polynomial order %s fit " % poly_order
   #   plt.ylabel("Residual Tb (K)")
   #   plt.xlabel("Frequency (MHz)")
   #   plt.legend(loc=1)
   #   plt.text(50, max_abs_residuals, "rms=%0.3f" % rms_of_residuals)
   #   plt.ylim([y_min, y_max])
   #   fig_name= "eda2_log_fit_residual_tsy_input_poly_%s_lst_%s%s.png" % (poly_order,lst_string,signal_type_postfix)
   #   figmap = plt.gcf()
   #   figmap.savefig(fig_name)
   #   print("saved %s" % fig_name) 
   
   #Get the theoretical diffuse value from the chan folders instead (in case you have run each chan separately for solve_from_uvfits) 
   linestyle_list = ['dashdot','-','dotted','--','solid','dotted']
   color_list=[color_dark_blue,color_orange_red,color_green,color_orange,color_yellow,color_pink,color_yellow]
   plotting_index = 0
   #unflagged
   plt.clf()
   fig, ax1 = plt.subplots()
   freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
   
   t_sky_monopole_list = []
   t_sky_monopole_array_to_plot_filename = "t_sky_monopole_array_plot.npy"
            
   for pol_index,pol in enumerate(pol_list):
      n_baselines_used_array_filename = "n_baselines_included_lst_%s_pol_%s%s.npy" % (lst_string,pol,signal_type_postfix)  
      if EDA2_data==True:
         label2 = 'beam wtd av GSM %s' % pol
         label4 = 'monopole GSM'
         #t_sky_theoretical_list = []
         t_sky_theoretical_array_to_plot_filename = "t_sky_theoretical_array_plot_pol_%s.npy" % pol

         EDA2_chan_list_input = EDA2_chan_list[0:len(freq_MHz_list)]
         for index,EDA2_chan in enumerate(EDA2_chan_list_input):
               freq_MHz = 400./512.*float(EDA2_chan)
               EDA2_chan_dir = "%s%s/" % (EDA2_data_dir,EDA2_chan)  
               EDA2_chan_index = int(EDA2_chan - 64)
               lst = lst_hrs_list[EDA2_chan_index]
               lst_deg = (float(lst)/24)*360.        
               #sky_averaged_diffuse_array_beam_lsts_filename =  "woden_map_start_freq_%0.3f_hpx_%s_%s_%s_%s_%s_%s_pol_%s_global_foreground.npy" % (start_freq,year,month,day,hour,min,sec,pol)      
               #sky_averaged_diffuse_array_beam_lsts_filename = "%seda_model_%s_lst_2.00_hr_int_0.13_hr_N_D_gsm_sky_averaged_diffuse_beam.npy" % (EDA2_chan_dir,pol)
               #sky_averaged_diffuse_array_beam_lsts_filename = "t_sky_theoretical_array_lst_%s_pol_%s%s.npy" % (lst_string,pol,signal_type_postfix)
               #sky_averaged_diffuse_array_beam_lsts_filename = "%s/sky_av_input_cal_%s_LST_%0.3f_%0.3f_MHz_pol_%s.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz,pol)
               #diffuse_global_value_array = np.load(sky_averaged_diffuse_array_beam_lsts_filename)
               #diffuse_global_value = diffuse_global_value_array[0]   
               #t_sky_theoretical_list.append(diffuse_global_value)
               
               #monopole is same for both pols
               if pol_index==0:
                  sky_global_monopole_temp_filename = "%s/sky_mono_temp_%s_LST_%0.3f_%0.3f_MHz.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz)
                  sky_global_monopole_array = np.load(sky_global_monopole_temp_filename)
                  sky_global_monopole_value = sky_global_monopole_array[0]
                  t_sky_monopole_list.append(sky_global_monopole_value)
                  
         #t_sky_theoretical_array = np.asarray(t_sky_theoretical_list)        
         #np.save(t_sky_theoretical_array_to_plot_filename,t_sky_theoretical_array)
         
         #if pol_index==0:
         #   t_sky_monopole_array = np.asarray(t_sky_monopole_list)
         #   np.save(t_sky_monopole_array_to_plot_filename,t_sky_monopole_array)
      else:
         label2 = 'input %s' % pol
         t_sky_theoretical_array_filename = "t_sky_theoretical_array_lst_%s_pol_%s%s.npy" % (lst_string,pol,signal_type_postfix)
         t_sky_theoretical_array = np.load(t_sky_theoretical_array_filename)
      
      #for plotting a subset of chans:
      if EDA2_data:
         fine_chans_per_EDA2_chan = n_fine_chans-5
         length_freq_MHz_fine_chan_to_plot = int(len(EDA2_chan_list_input) * fine_chans_per_EDA2_chan)
      else:
         length_freq_MHz_fine_chan_to_plot = len(freq_MHz_list)
      
      for model_type_index,model_type in enumerate(model_type_list):
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='measured Tsky %s' % pol 
               label5='sim measured Tsky %s' % pol
            else:
               #label1='ignore angular response'
               #fig5: and fig6a
               #label1='recovered'
               #fig9a:
               label1='ignore angular %s' % pol
               label2='input %s' % pol
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular %s' % pol
         else:
            label1='recovered %s' % pol
            
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_beam_wtd_av_array_filename = "t_sky_beam_wtd_av_lst_%s_pol_%s.npy" % (lst_string,pol)
         
         t_sky_beam_wtd_av_array = np.load(t_sky_beam_wtd_av_array_filename)
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
         n_baselines_used_array = np.load(n_baselines_used_array_filename)
         
         #fig5 and fig6a, paper 1 
         ax1.set_xlabel("Frequency (MHz)")
         ax1.set_ylabel("Sky temperature (K)", color='black')
         ax1.tick_params(axis='y', labelcolor='black')
         ax1.errorbar(freq_MHz_fine_array,t_sky_measured_array,yerr=t_sky_measured_error_array,label=label1,color=color_list[plotting_index],linestyle=linestyle_list[plotting_index],alpha=0.7)
         plotting_index += 1
         #plt.errorbar(freq_MHz_fine_array,t_sky_measured_array,yerr=t_sky_measured_error_array,label=label1,color=color_list[model_type_index],linestyle=model_type_linestyle_list[model_type_index+1],alpha=0.7)
      if len(freq_MHz_list)==1:
         #ax1.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
         ax1.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2)
         plotting_index += 1
         #ax1.legend(loc=1)
      else:
         #ax1.plot(freq_MHz_list,t_sky_theoretical_array,label=label2,color=color_list[plotting_index],linestyle=linestyle_list[plotting_index])
         ax1.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2,color=color_list[plotting_index],linestyle=linestyle_list[plotting_index])
         plotting_index += 1
         #ax1.legend(loc=1)
         #fig5
         #ax2 = ax1.twinx()
         #ax2.set_ylabel('no. baselines used', color=color_green)  # we already handled the x-label with ax1
         #ax2.tick_params(axis='y', labelcolor=color_green)
         #ax2.plot(freq_MHz_list,n_baselines_used_array,label='no. baselines used',color=color_green,linestyle='--')
         #ax2.legend(loc=0)
         lines, labels = ax1.get_legend_handles_labels()
         #lines2, labels2 = ax2.get_legend_handles_labels()
         #ax1.legend(lines + lines2, labels + labels2, loc=7)
         #fig9:
         #plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2,color=color_green,linestyle=':')
      #if 'diffuse_global' in signal_type_list:
      #   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
      #if include_angular_info:
      #   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured" 
   #plt.xlabel("Frequency (MHz)")
   #plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      print(signal_type_list)
      plt.legend(loc='upper right')
   else:
      #comment out for fig5
      plt.legend(loc='lower right')
   #if EDA2_data:
   #   plt.ylim([500, 5000])
   #else:
      #plt.ylim([-1, 0.5])
      #commented out for fig5
   fig_name= "t_sky_measured_lst_%s%s.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   #percentage_diff = ((t_sky_theoretical_array - t_sky_measured_array) / t_sky_theoretical_array) * 100.
   #print(percentage_diff)
   #percentage_diff_max_index = np.nanargmax(percentage_diff)
   #print(percentage_diff_max_index)
   #percentage_diff_max_freq = freq_MHz_list[percentage_diff_max_index]
   #print(percentage_diff_max_freq)
   #print(percentage_diff[percentage_diff_max_index])

   plt.clf()
   
   for pol_index,pol in enumerate(pol_list):
      ###Also plot the average measurement for each EDA2 coarse chan
      t_sky_measure_av_per_EDA2_chan = np.full(len(freq_MHz_list),np.nan)
      t_sky_measure_av_per_EDA2_chan_err = np.full(len(freq_MHz_list),np.nan)
      t_sky_measure_av_per_EDA2_chan_weighted = np.full(len(freq_MHz_list),np.nan)
      t_sky_measure_av_per_EDA2_chan_err_weighted = np.full(len(freq_MHz_list),np.nan)
      for model_type in model_type_list:
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='measured Tsky %s' % pol
               label2 = 'beam wtd av GSM %s' % pol
            else:
               label1='ignore angular %s' % pol
            label3='ignore angular resp. weighted'
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular %s' % pol
            label3='subtract angular resp. weighted'
         else:
            label1='recovered %s' % pol
            label3='recovered weighted'
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
         
       
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[0:length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[0:length_freq_MHz_fine_chan_to_plot]
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[0:length_freq_MHz_fine_chan_to_plot]
         
         #if EDA2_data:
         #   t_sky_theoretical_array_to_plot_filename = "t_sky_theoretical_array_plot_pol_%s.npy" % pol
         #   t_sky_theoretical_array = np.load(t_sky_theoretical_array_to_plot_filename)
          
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
         #plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
         plt.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2)
         
         #plt.scatter(freq_MHz_list,t_sky_monopole_array,label=label3)
      else:
         #plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2)
         plt.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2)
         #plt.plot(freq_MHz_list,t_sky_monopole_array,label=label3)
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
   #else:
      #plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s_per_chan_av.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   ####Flagged fig11 paper1
   plt.clf()
   
   for pol_index,pol in enumerate(pol_list):
      for model_type in model_type_list:
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='measured Tsky %s' % pol
               label5='sim measured Tsky %s' % pol
               label4='monopole T_sky'
               label2='beam wtd av %s' % pol
            else:
               label1='ignore angular %s' % pol
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular %s' % pol
            label5='sim subtr. ang. %s' % pol
         else:
            label1='recovered %s' % pol
            
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)

         t_sky_sim_measured_array_filename = "t_sky_sim_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_sim_measured_error_array_filename = "t_sky_sim_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
                
         start_fine_chan = int(fine_chans_per_EDA2_chan * chan_num)
         
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_sim_measured_array = np.load(t_sky_sim_measured_array_filename)
         t_sky_sim_measured_array = t_sky_sim_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_sim_measured_error_array = np.load(t_sky_sim_measured_error_array_filename)
         t_sky_sim_measured_error_array = t_sky_sim_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
                  
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
          
         if EDA2_data:
            length_freq_MHz_fine_chan_to_plot_theoretical = length_freq_MHz_fine_chan_to_plot / fine_chans_per_EDA2_chan
            #t_sky_theoretical_array_to_plot_filename = "t_sky_theoretical_array_plot_pol_%s.npy" % pol
            #t_sky_theoretical_array = np.load(t_sky_theoretical_array_to_plot_filename)
            #t_sky_theoretical_array = t_sky_theoretical_array[0:length_freq_MHz_fine_chan_to_plot_theoretical]
            
            t_sky_monopole_array_to_plot_filename = "t_sky_monopole_array_plot.npy"
            t_sky_monopole_array = np.load(t_sky_monopole_array_to_plot_filename)
            t_sky_monopole_array = t_sky_monopole_array[0:length_freq_MHz_fine_chan_to_plot_theoretical]
            freq_MHz_list_coarse = freq_MHz_list[0:length_freq_MHz_fine_chan_to_plot_theoretical]
          
         plt.errorbar(freq_MHz_fine_array,t_sky_measured_array,yerr=t_sky_measured_error_array,label=label1,linestyle='-',alpha=0.7)
      
         #and sims
         plt.errorbar(freq_MHz_fine_array,t_sky_sim_measured_array,yerr=t_sky_sim_measured_error_array,label=label5,alpha=0.7)         
         
      #theoretical fine chan
      ax1.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2,linestyle=':')
         
      ##hack to remove theoretical line for smoothness test
      #if len(freq_MHz_list)==1:
      #   plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
      #   #if pol_index==0:
      #   #   plt.scatter(freq_MHz_list,t_sky_monopole_array,label=label4)
      #else: 
      #   plt.plot(freq_MHz_list_coarse,t_sky_theoretical_array,label=label2,linestyle=':')
      #   #if pol_index==0:
      #   #   plt.plot(freq_MHz_list_coarse,t_sky_monopole_array,label=label4)
      #   #plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2,linestyle=':')
      
      
      
      ###if 'diffuse_global' in signal_type_list:
      ###   plt.plot(freq_MHz_list,diffuse_global_value_array,label='input')
      ###if include_angular_info:
      ###   plt.plot(freq_MHz_list,t_sky_measured_global_array,label='with ang info')

   map_title="t_sky measured flagged" 
   plt.xlabel("Frequency (MHz)")
   plt.ylabel("Sky temperature (K)")
   if ('diffuse_global' in signal_type_list or 'diffuse' in signal_type_list or 'diffuse_angular' in signal_type_list):
      print(signal_type_list)
      plt.legend(loc='upper right')
   else:
      plt.legend(loc='lower right')
   #if EDA2_data:
   #   plt.ylim([500, 5000])
      #check first 5 chans:
      #plt.ylim([3500, 4500]) #hack
   #else:
      #plt.ylim([0, 4000])
   fig_name= "t_sky_measured_lst_%s%s_flagged.png" % (lst_string,signal_type_postfix)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name) 
   
   plt.clf()
   
   for pol_index,pol in enumerate(pol_list):
      ###Also plot the average measurement for each EDA2 coarse chan
      t_sky_measure_av_per_EDA2_chan = np.full(len(freq_MHz_list),np.nan)
      t_sky_measure_av_per_EDA2_chan_err = np.full(len(freq_MHz_list),np.nan)
      for model_type in model_type_list:
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='measured Tsky %s' % pol
            else:
               label1='ignore angular %s' % pol
            label3='ignore angular resp. weighted'
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular %s' % pol
            label3='subtract angular resp. weighted'
         else:
            label1='recovered %s' % pol
            label3='recovered weighted'
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
       
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
        
         #if EDA2_data:
         #   t_sky_theoretical_array_to_plot_filename = "t_sky_theoretical_array_plot_pol_%s.npy" % pol
         #   t_sky_theoretical_array = np.load(t_sky_theoretical_array_to_plot_filename)
             
         if EDA2_data==True:
            for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list_coarse):
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
            
         plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan,yerr=t_sky_measure_av_per_EDA2_chan_err,label=label1,linestyle='-',alpha=0.7)
         #if EDA2_data:
         #   plt.errorbar(freq_MHz_list,t_sky_measure_av_per_EDA2_chan_weighted,yerr=t_sky_measure_av_per_EDA2_chan_err_weighted,label=label3)
      
      #theoretical fine chan
      ax1.plot(freq_MHz_fine_array,t_sky_beam_wtd_av_array,label=label2,linestyle=':')
      #if len(freq_MHz_list)==1:
      #   plt.scatter(freq_MHz_list,t_sky_theoretical_array,label=label2)
      #else:
      #   plt.plot(freq_MHz_list,t_sky_theoretical_array,label=label2,linestyle=':')
   
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
   #if EDA2_data:
   #plt.ylim([500, 5000])
   #else:
      #plt.ylim([0, 4000])
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
   linestyle_list = ['dashdot','solid','dotted','dashed','dashdot','solid','dotted','dashed']
   color_list=[color_dark_blue,color_orange_red,color_green,color_orange,color_black,color_light_blue,color_pink,color_yellow]
   plotting_index = 0
   #unflagged
   plt.clf()
   fig, ax1 = plt.subplots()
   #also take the average of the X and Y pol residuals
   for pol_index,pol in enumerate(pol_list):
      plot_log = False
      for model_type_index,model_type in enumerate(model_type_list):
         #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
         if model_type=='OLS_fixed_intercept':
            if EDA2_data:
               label1='measured Tsky %s' % pol
               label2='sim measured Tsky %s' % pol
            else:
               #fig9
               label1='ignore angular %s' % pol
               #fig7:
               #label1='residual from log fit'
            y_offset=1
            #colour='tab:blue'
            colour=color_dark_blue
         elif  model_type=='OLS_fixed_int_subtr_Y':
            label1='subtract angular %s' % pol
            label2='sim subtr. ang. %s' % pol
            y_offset=0
            #colour='tab:orange'
            colour=color_orange_red
         else:
            label1='recovered %s' % pol
         
         #for fig9b (2020) change this to unflagged
         #old, no pol
         #t_sky_measured_array_filename = "t_sky_measured_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
         #t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s%s_%s_flagged.npy" % (lst_string,signal_type_postfix,model_type)
        
         t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)

         t_sky_sim_measured_array_filename = "t_sky_sim_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         t_sky_sim_measured_error_array_filename = "t_sky_sim_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)

         #for sims use unflagged (need to update model_from_saved_data to use subtr_Y data properly for flagged ATM does nothing!):
         #t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         #t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s.npy" % (lst_string,pol,signal_type_postfix,model_type)
         
         freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
         
         t_sky_measured_array = np.load(t_sky_measured_array_filename)
         t_sky_measured_array = t_sky_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
         t_sky_measured_error_array = t_sky_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         
         t_sky_beam_wtd_av_array = t_sky_beam_wtd_av_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         
         t_sky_sim_measured_array = np.load(t_sky_sim_measured_array_filename)
         t_sky_sim_measured_array = t_sky_sim_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         t_sky_sim_measured_error_array = np.load(t_sky_sim_measured_error_array_filename)
         t_sky_sim_measured_error_array = t_sky_sim_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
                 
         freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
         freq_MHz_fine_array = freq_MHz_fine_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
         n_baselines_used_array = np.load(n_baselines_used_array_filename)
         n_baselines_used_array = n_baselines_used_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
          

         #subtract a polynomial fit
         #in log log space:
         sky_array = t_sky_measured_array[t_sky_measured_array>0.]
         sky_sim_array = t_sky_sim_measured_array[t_sky_measured_array>0.]
         sky_beam_wtd_av_array = t_sky_beam_wtd_av_array[t_sky_measured_array>0.]
         
         if not EDA2_data:
            t_sky_theoretical_array_cut = t_sky_theoretical_array[t_sky_measured_array>0.]
            n_baselines_used_array_cut = n_baselines_used_array[t_sky_measured_array>0.]
            t_sky_measured_error_array_cut = t_sky_measured_error_array[t_sky_measured_array>0.]
            
         log_sky_array = np.log10(sky_array)
         log_sky_sim_array = np.log10(sky_sim_array)
         log_sky_beam_wtd_av_array = np.log10(sky_beam_wtd_av_array)
         
         if n_fine_chans_used==1:
            freq_array_cut = freq_MHz_array[t_sky_measured_array>0.]
         else:
            freq_array_cut = freq_MHz_fine_array[t_sky_measured_array>0.]
  
         log_freq_MHz_array = np.log10(freq_array_cut)
         
         if len(log_sky_array) != 0:
            plot_log = True
            coefs = poly.polyfit(log_freq_MHz_array, log_sky_array, poly_order)
            ffit = poly.polyval(log_freq_MHz_array, coefs)
            ffit_linear = 10**ffit
            
            #log_residual = log_signal_array_short_baselines - log_ffit
            residual_of_log_fit = ffit_linear - sky_array
            if pol_index==0 and model_type=='OLS_fixed_intercept': 
               both_pol_residual_sum_array = np.zeros(len(freq_array_cut))
            if pol_index==0 and model_type=='OLS_fixed_int_subtr_Y':
               both_pol_residual_sum_array_subtr_Y = np.zeros(len(freq_array_cut))
            if model_type=='OLS_fixed_intercept':
               if len(both_pol_residual_sum_array)==len(residual_of_log_fit):
                  both_pol_residual_sum_array += residual_of_log_fit
            if model_type=='OLS_fixed_int_subtr_Y':
               both_pol_residual_sum_array_subtr_Y += residual_of_log_fit
            rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
            print("rms_of_residuals data %s %s is %0.3f K" % (model_type,pol,rms_of_residuals))
            
            max_abs_residuals = np.max(np.abs(residual_of_log_fit))
            #
            #y_max = 1.5 * max_abs_residuals
            #y_min = 1.5 * -max_abs_residuals
            if model_type!='OLS_fixed_int_subtr_Y':
               if pol!='Y':
                  y_max = 1.5 * max_abs_residuals
                  y_min = 1.5 * -max_abs_residuals
               else:
                  y_max = 1.5 * max_abs_residuals
                  y_min = 1.5 * -max_abs_residuals
            
            ##temporary just for paper fig12a:
            #y_max = 100
            #y_min = -100
            
            #print(sky_array)
            #print(residual_of_log_fit)
            
            #fig7b:
            #PUT THIS BACK IN JUST TEMP TO CHECK SIMS plt.plot(freq_array_cut,residual_of_log_fit,label=label1,linestyle=linestyle_list[plotting_index],color=color_list[plotting_index])
            #plt.plot(freq_array_cut,residual_of_log_fit,label=label1,linestyle=linestyle_list[model_type_index])
            #plt.text(50, max_abs_residuals + y_offset, "%srms=%1.2f K" % (linestyle_list[model_type_index],rms_of_residuals),{'color': colour})
            #plt.text(50, 75, "rms=%2.1f K" % rms_of_residuals,{'color': colour})
            y_offset = 10. + plotting_index * 1.3
            plt.text(50, y_offset, "rms=%0.2f K" % rms_of_residuals,{'color': color_list[plotting_index]})
            plotting_index += 1
            
            #comment out for fig9b
            if not EDA2_data:
               if not woden:
                  if model_type_index==0:
                     expected_noise = plot_expected_rms_noise_eda2(freq_MHz_list=freq_array_cut,t_sky_theoretical_array=t_sky_theoretical_array_cut,n_baselines_used_array=n_baselines_used_array_cut,int_time=int_time,bandwidth_Hz=bw_Hz)
                     plt.plot(freq_array_cut,expected_noise,label="expected rms noise",color='red',linestyle='--')
                     #for referee comments on fig7b:
                     plt.plot(freq_array_cut,t_sky_measured_error_array_cut,label="OLS fit error",color='green',linestyle='-.')
                     #print(expected_noise)
         #EDA2_sims
         if EDA2_data:
            if len(log_sky_sim_array) != 0:
               plot_log = True
               coefs = poly.polyfit(log_freq_MHz_array, log_sky_sim_array, poly_order)
               ffit = poly.polyval(log_freq_MHz_array, coefs)
               ffit_linear = 10**ffit
               
               #log_residual = log_signal_array_short_baselines - log_ffit
               residual_of_log_fit = ffit_linear - sky_sim_array
               if pol_index==0 and model_type=='OLS_fixed_intercept': 
                  both_pol_residual_sum_array = np.zeros(len(freq_array_cut))
               if pol_index==0 and model_type=='OLS_fixed_int_subtr_Y':
                  both_pol_residual_sum_array_subtr_Y = np.zeros(len(freq_array_cut))
               if model_type=='OLS_fixed_intercept':
                  if len(both_pol_residual_sum_array)==len(residual_of_log_fit):
                     both_pol_residual_sum_array += residual_of_log_fit
               if model_type=='OLS_fixed_int_subtr_Y':
                  both_pol_residual_sum_array_subtr_Y += residual_of_log_fit
               rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
               print("rms_of_residuals sims %s %s is %0.3f K" % (model_type,pol,rms_of_residuals))
               
               max_abs_residuals = np.max(np.abs(residual_of_log_fit))
               #
               #y_max = 1.5 * max_abs_residuals
               #y_min = 1.5 * -max_abs_residuals
               if model_type!='OLS_fixed_int_subtr_Y':
                  if pol!='Y':
                     y_max = 1.5 * max_abs_residuals
                     y_min = 1.5 * -max_abs_residuals
                  else:
                     y_max = 1.5 * max_abs_residuals
                     y_min = 1.5 * -max_abs_residuals
               
               ##temporary just for paper fig12a:
               #y_max = 100
               #y_min = -100
               
               #print(sky_array)
               #print(residual_of_log_fit)
               
               #fig7b:
               #print("plotting index %s" % plotting_index)
               plt.plot(freq_array_cut,residual_of_log_fit,label=label2,linestyle=linestyle_list[plotting_index],color=color_list[plotting_index])
               #plt.plot(freq_array_cut,residual_of_log_fit,label=label1,linestyle=linestyle_list[model_type_index])
               #plt.text(50, max_abs_residuals + y_offset, "%srms=%1.2f K" % (linestyle_list[model_type_index],rms_of_residuals),{'color': colour})
               #plt.text(50, 75, "rms=%2.1f K" % rms_of_residuals,{'color': colour})
               y_offset = 10. + plotting_index * 1.3
               plt.text(50, y_offset, "rms=%0.2f K" % rms_of_residuals,{'color': color_list[plotting_index]})
               plotting_index += 1
                        
         #Check that the beam weighted average input is in fact smooth

         #beam wtd average fine chan
         if EDA2_data:
            if len(log_sky_beam_wtd_av_array) != 0:
               plot_log = True
               coefs = poly.polyfit(log_freq_MHz_array, log_sky_beam_wtd_av_array, poly_order)
               ffit = poly.polyval(log_freq_MHz_array, coefs)
               ffit_linear = 10**ffit
               
               #log_residual = log_signal_array_short_baselines - log_ffit
               residual_of_log_fit = ffit_linear - sky_beam_wtd_av_array
               if pol_index==0 and model_type=='OLS_fixed_intercept': 
                  both_pol_residual_sum_array = np.zeros(len(freq_array_cut))
               if pol_index==0 and model_type=='OLS_fixed_int_subtr_Y':
                  both_pol_residual_sum_array_subtr_Y = np.zeros(len(freq_array_cut))
               if model_type=='OLS_fixed_intercept':
                  if len(both_pol_residual_sum_array)==len(residual_of_log_fit):
                     both_pol_residual_sum_array += residual_of_log_fit
               if model_type=='OLS_fixed_int_subtr_Y':
                  both_pol_residual_sum_array_subtr_Y += residual_of_log_fit
               rms_of_residuals = np.sqrt(np.mean(residual_of_log_fit**2))
               print("rms_of_residuals beam wtd av %s %s is %0.3f K" % (model_type,pol,rms_of_residuals))
               
               max_abs_residuals = np.max(np.abs(residual_of_log_fit))
               #
               #y_max = 1.5 * max_abs_residuals
               #y_min = 1.5 * -max_abs_residuals
               if model_type!='OLS_fixed_int_subtr_Y':
                  if pol!='Y':
                     y_max = 1.5 * max_abs_residuals
                     y_min = 1.5 * -max_abs_residuals
                  else:
                     y_max = 1.5 * max_abs_residuals
                     y_min = 1.5 * -max_abs_residuals
               
               ##temporary just for paper fig12a:
               #y_max = 100
               #y_min = -100
               
               #print(sky_array)
               #print(residual_of_log_fit)
               
               #fig7b:
               #print("plotting index %s" % plotting_index)
               plt.plot(freq_array_cut,residual_of_log_fit,label="beam_wtd_av %s" % pol)
               #plt.plot(freq_array_cut,residual_of_log_fit,label=label1,linestyle=linestyle_list[model_type_index])
               #plt.text(50, max_abs_residuals + y_offset, "%srms=%1.2f K" % (linestyle_list[model_type_index],rms_of_residuals),{'color': colour})
               #plt.text(50, 75, "rms=%2.1f K" % rms_of_residuals,{'color': colour})
               y_offset = 10. + plotting_index * 1.3
               plt.text(50, y_offset, "rms=%0.2f K" % rms_of_residuals)
               #plotting_index += 1
                        
         #Check that the beam weighted average input is in fact smooth

                
   #good for most EDA2 plots:
   y_min, y_max = -20 , 20
   
   if model_type=='OLS_fixed_intercept': 
      both_pol_residual_av_array = both_pol_residual_sum_array / float(len(pol_list))
      rms_of_residuals_combined = np.sqrt(np.mean(both_pol_residual_av_array**2))
   if 'OLS_fixed_int_subtr_Y' in model_type_list:
      both_pol_residual_av_array_subtr_Y = both_pol_residual_sum_array_subtr_Y / float(len(pol_list))
      rms_of_residuals_combined_subtr_Y = np.sqrt(np.mean(both_pol_residual_av_array_subtr_Y**2))
   
   if plot_log == True:
      #fig9b paper1 (#and 7b)
      map_title="Residual for log polynomial order %s fit " % poly_order
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("Frequency (MHz)")
      #if len(model_type_list)>1:
      #   plt.legend(loc=1)
      plt.legend(loc=1)
      #plt.ylim([y_min, y_max])
      fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_no_e_bars.png" % (poly_order,lst_string,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      plt.close()         

      plt.clf()
      fig, ax1 = plt.subplots()
      if model_type=='OLS_fixed_intercept': 
         plt.plot(freq_array_cut,both_pol_residual_av_array,color=color_dark_blue,label='comb. pol. ignore ang.')
         plt.text(50, 2.5, "rms=%0.5f K" % rms_of_residuals_combined,color=color_dark_blue)
      if 'OLS_fixed_int_subtr_Y' in model_type_list:
         plt.plot(freq_array_cut,both_pol_residual_av_array_subtr_Y,color=color_orange,label='comb. pol. subtr. ang.')
      if 'OLS_fixed_int_subtr_Y' in model_type_list:
         plt.text(50, 2, "rms=%0.5f K" % rms_of_residuals_combined_subtr_Y,color=color_orange)
      map_title="Residual for combined pols "
      plt.ylabel("Residual Tb (K)")
      plt.xlabel("Frequency (MHz)")
      #if len(model_type_list)>1:
      #   plt.legend(loc=1)
      plt.legend(loc=1)
      plt.ylim([y_min, y_max])
      fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_no_e_bars_combined_pol.png" % (poly_order,lst_string,signal_type_postfix)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name)
      plt.close()  
            
      #plt.clf()
      #plt.errorbar(freq_array_cut,residual_of_log_fit,yerr=t_sky_measured_error_array[t_sky_measured_array>0.],label='residual of log fit')
      ##if include_angular_info:
      ##   plt.plot(freq_array_cut,residual_of_log_fit_glob_ang,label='residual of glob ang')
      #map_title="Residual for log polynomial order %s fit " % poly_order
      #plt.ylabel("Residual Tb (K)")
      #plt.xlabel("Frequency (MHz)")
      #if len(model_type_list)>1:
      #   plt.legend(loc=1)
      #plt.text(50, max_abs_residuals, "rms=%0.3f K" % rms_of_residuals)
      #plt.ylim([y_min, y_max])
      #fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_%s.png" % (poly_order,lst_string,signal_type_postfix,model_type)
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name) 
      #plt.close() 

      #repeat for per chan av
      if EDA2_data:
         plt.clf()
         for model_type in model_type_list:
            t_sky_measure_av_per_EDA2_chan_comb_pol = np.full(len(freq_MHz_list),np.nan)
            for pol_index,pol in enumerate(pol_list):
               t_sky_measure_av_per_EDA2_chan = np.full(len(freq_MHz_list),np.nan)
               #['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
               if model_type=='OLS_fixed_intercept':
                  if EDA2_data:
                     label1='ignore angular %s' % pol
                     label3='ignore ang comb pol'
                  else:
                     label1='ignore angular %s' % pol
                  y_offset=1
                  colour='tab:blue'
               elif  model_type=='OLS_fixed_int_subtr_Y':
                  label1='subtract angular %s' % pol
                  y_offset=0
                  colour='tab:orange'
                  label3='subtr ang comb pol'
               else:
                  label1='recovered %s' % pol
                  
               t_sky_measured_array_filename = "t_sky_measured_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
               t_sky_measured_error_array_filename = "t_sky_measured_error_array_lst_%s_pol_%s%s_%s_flagged.npy" % (lst_string,pol,signal_type_postfix,model_type)
               freq_MHz_fine_array_filename = "freq_MHz_fine_array_lst_%s%s.npy" % (lst_string,signal_type_postfix)
               
             
               t_sky_measured_array = np.load(t_sky_measured_array_filename)
               t_sky_measured_array = t_sky_measured_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
               
               if pol_index == 0:
                  t_sky_measured_comb_pol_sum_array = np.load(t_sky_measured_array_filename)
                  t_sky_measured_comb_pol_sum_array = t_sky_measured_comb_pol_sum_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
               
               t_sky_measured_error_array = np.load(t_sky_measured_error_array_filename)
               t_sky_measured_error_array = t_sky_measured_error_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
               
               if pol_index != 0:
                  t_sky_measured_comb_pol_sum_array += t_sky_measured_error_array
               freq_MHz_fine_array = np.load(freq_MHz_fine_array_filename)
               freq_MHz_fine_array = freq_MHz_fine_array[start_fine_chan:start_fine_chan+length_freq_MHz_fine_chan_to_plot]
               n_baselines_used_array = np.load(n_baselines_used_array_filename)
              
               for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list_coarse):
                  freq_range_min = freq_MHz - (centre_chan_index * fine_chan_width_Hz/1000000.)
                  freq_range_max = freq_MHz + (centre_chan_index * fine_chan_width_Hz/1000000.)
                  #print(freq_range_min)
                  #print(freq_range_max)
                  
                  indices = np.where(np.logical_and(freq_MHz_fine_array>=freq_range_min,freq_MHz_fine_array<=freq_range_max))
                  t_sky_measured_EDA2_chan = t_sky_measured_array[indices]
                  t_sky_measured_error_chan = t_sky_measured_error_array[indices]
                  t_sky_measure_av_per_EDA2_chan[freq_MHz_index] = np.nanmean(t_sky_measured_EDA2_chan)
                  t_sky_measure_av_per_EDA2_chan_err[freq_MHz_index] = np.nanstd(t_sky_measured_EDA2_chan)
                
               #subtract a polynomial fit
               #in log log space:
               sky_array = t_sky_measure_av_per_EDA2_chan[t_sky_measure_av_per_EDA2_chan>0.]
               log_sky_array = np.log10(sky_array)
               freq_array_cut = freq_MHz_array[t_sky_measure_av_per_EDA2_chan>0.]
               if not EDA2_data:
                  t_sky_theoretical_array_cut = t_sky_theoretical_array[t_sky_measure_av_per_EDA2_chan>0.]
                  n_baselines_used_array_cut = n_baselines_used_array[t_sky_measure_av_per_EDA2_chan>0.]
               
               
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
               
               #print(freq_array_cut)
               #print(residual_of_log_fit)
               
               plt.plot(freq_array_cut,residual_of_log_fit,label=label1)
               plt.text(50, max_abs_residuals + y_offset, "rms=%1.2f K" % rms_of_residuals,{'color': colour})
               
               
               #include expected noise estimate:
               if not EDA2_data:
                  expected_noise = plot_expected_rms_noise_eda2(freq_MHz_list=freq_array_cut,t_sky_theoretical_array=t_sky_theoretical_array_cut,n_baselines_used_array=n_baselines_used_array_cut,int_time=int_time,bandwidth_Hz=bw_Hz)
                  plt.plot(freq_array_cut,expected_noise,label="expected rms noise")
         
            if len(pol_list) > 1:
               t_sky_measured_comb_pol_av_array = t_sky_measured_comb_pol_sum_array / float(len(pol_list))
               for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
                  freq_range_min = freq_MHz - (centre_chan_index * fine_chan_width_Hz/1000000.)
                  freq_range_max = freq_MHz + (centre_chan_index * fine_chan_width_Hz/1000000.)
                  #print(freq_range_min)
                  #print(freq_range_max)
                  
                  indices = np.where(np.logical_and(freq_MHz_fine_array>=freq_range_min,freq_MHz_fine_array<=freq_range_max))
                  t_sky_measured_EDA2_chan = t_sky_measured_comb_pol_av_array[indices]
                  #t_sky_measured_error_chan = t_sky_measured_error_array[indices]
                  t_sky_measure_av_per_EDA2_chan_comb_pol[freq_MHz_index] = np.nanmean(t_sky_measured_EDA2_chan)
                  #t_sky_measure_av_per_EDA2_chan_err[freq_MHz_index] = np.nanstd(t_sky_measured_EDA2_chan)
               #subtract a polynomial fit
               #in log log space:
               sky_array = t_sky_measure_av_per_EDA2_chan_comb_pol[t_sky_measure_av_per_EDA2_chan_comb_pol>0.]
               log_sky_array = np.log10(sky_array)
               freq_array_cut = freq_MHz_array[t_sky_measure_av_per_EDA2_chan_comb_pol>0.]
               if not EDA2_data:
                  t_sky_theoretical_array_cut = t_sky_theoretical_array[t_sky_measure_av_per_EDA2_chan_comb_pol>0.]
                  n_baselines_used_array_cut = n_baselines_used_array[t_sky_measure_av_per_EDA2_chan_comb_pol>0.]
               
               
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
               
               #print(freq_array_cut)
               #print(residual_of_log_fit)
               
               plt.plot(freq_array_cut,residual_of_log_fit,label=label3)
               plt.text(50, max_abs_residuals + y_offset, "rms=%1.2f K" % rms_of_residuals,{'color': colour})
         
         
         map_title="Residual for log polynomial order %s fit " % poly_order
         plt.ylabel("Residual Tb (K)")
         plt.xlabel("Frequency (MHz)")
         #if len(model_type_list)>1:
         plt.legend(loc=1)
         
         #plt.legend(loc=1)
         plt.ylim([y_min, y_max])
         fig_name= "eda2_log_fit_residual_tsy_measured_poly_%s_lst_%s%s_no_e_bars_per_chan_av.png" % (poly_order,lst_string,signal_type_postfix)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         plt.close()    
      
               
   #n baselines plot
   if not EDA2_data:
      plt.clf()
      plt.plot(freq_MHz_list,n_baselines_used_array)
      map_title="n_baselines_included" 
      plt.xlabel("Frequency (MHz)")
      plt.ylabel("number of baselines")
      #plt.legend(loc=1)
      #plt.ylim([0, 20])
      fig_name= "n_baselines_included_lst_%s%s_%s.png" % (lst_string,signal_type_postfix,model_type)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      plt.close()  
           
def plot_expected_rms_noise_eda2(freq_MHz_list,t_sky_theoretical_array,n_baselines_used_array,int_time,bandwidth_Hz):
   print("plotting expected noise")
   #for good explanation from LWA: http://lwa.phys.unm.edu/obsstatus/obsstatus006.html#toc14
   #from EDA1 paper wayth et al 2017, table 2:
   A_eff_array_per_dipole = np.asarray([970.,950.,914.,874.,832.,771.,707.,638.,568.,498.,435.,377.,329.,288.,252.,222.,196.]) / 256.
   freq_MHz_for_A_eff_array_per_dipole = np.asarray([60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.])    
   #fit a line:
   coefs = poly.polyfit(freq_MHz_for_A_eff_array_per_dipole, A_eff_array_per_dipole, 1)
   ffit = poly.polyval(freq_MHz_for_A_eff_array_per_dipole, coefs)
   
   #plot it to check 
   #plt.clf()
   #plt.plot(freq_MHz_for_A_eff_array_per_dipole,A_eff_array_per_dipole)
   #plt.plot(freq_MHz_for_A_eff_array_per_dipole,ffit)
   #map_title="A_eff EDA2 dipoles" 
   #plt.xlabel("Frequency (MHz)")
   #plt.ylabel("A_eff (m)")
   ##plt.legend(loc=1)
   ##plt.ylim([0, 20])
   #fig_name= "A_eff_EDA2.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name) 
   #plt.close()  
   
   A_eff_for_calc_array = poly.polyval(freq_MHz_list, coefs)
   
   T_rms1 = t_sky_theoretical_array / np.sqrt(2 * n_baselines_used_array * int_time * bandwidth_Hz)
   #T_rms1 = t_sky_theoretical_array / np.sqrt(n_baselines_used_array * int_time * bandwidth_Hz)
   #T_rms2 = ((c**2 / (freq_MHz_list*1000000.)**2) / (2*A_eff_for_calc_array)) * T_rms1
   
   #This is the one in the original paper submitted, with the extra factor of 2 at the bottom....remove it for referee response and consistency with TMS p232
   #T_rms3 =  ((c**2 / (freq_MHz_list*1000000.)**2) / (2*A_eff_for_calc_array)) * T_rms1
   T_rms3 =  ((c**2 / (freq_MHz_list*1000000.)**2) / (A_eff_for_calc_array)) * T_rms1
   
   
   #plt.clf()
   ##plt.plot(freq_MHz_list,T_rms1,label='T_rms1')
   ###plt.plot(freq_MHz_list,T_rms2,label='T_rms2')
   #plt.plot(freq_MHz_list,T_rms3,label='T_rms3')
   #map_title="T_rms EDA2" 
   #plt.xlabel("Frequency (MHz)")
   #plt.ylabel("T_rms (K)")
   ###plt.legend(loc=1)
   ###plt.ylim([0, 20])
   #fig_name= "T_rms_EDA2.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name) 
   #plt.close() 
   
   #print(freq_MHz_list)
   #print(n_baselines_used_array)
   #print(t_sky_theoretical_array)
   
   #print(T_rms3)
   
   return T_rms3  
   

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
                  #plt.xlabel("Frequency (MHz)")
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
               #plt.xlabel("Frequency (MHz)")
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
   #plt.xlabel("Frequency (MHz)")
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
      plt.xlabel("Frequency (MHz)")
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
      plt.xlabel("Frequency (MHz)")
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
   
   reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
   
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
            
            reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%0.3f_MHz_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)
            
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
                   
def simulate_sitara(lst_hrs,nside=512):
   n_ants = 2
   freq_MHz_array = np.arange(70,70+151)
   npix = hp.nside2npix(nside)
   lst_hrs = float(lst_hrs)
   lst_deg = lst_hrs * 15.
   unity_sky_value = 1.
   
   #hpx rotate stuff, rotate the complex beams to zenith at the required LST
   dec_rotate = 90. - float(mwa_latitude_ephem)
   ra_rotate = lst_deg
   r_beam_dec = hp.Rotator(rot=[0,dec_rotate], coord=['C', 'C'], deg=True) 
   r_beam_ra = hp.Rotator(rot=[ra_rotate,0], coord=['C', 'C'], deg=True)    
   
   EEP_name = '/md0/EoR/EDA2/EEPs/SITARA/chall_beam_Y.mat'
   #SITARA MATLAB beam file here: /md0/EoR/EDA2/EEPs/SITARA/chall_beam_Y.mat (1 MHz res) (70 - 200 MHz?)
   beam_data = loadmat(EEP_name)
   print("loaded %s " % EEP_name)
   #print(beam_data.keys())
   
   E_phi_element_1_cube = beam_data['Ephi_elem1']
   E_phi_element_2_cube = beam_data['Ephi_elem2']
   E_theta_element_1_cube = beam_data['Etheta_elem1']
   E_theta_element_2_cube = beam_data['Etheta_elem2']

   #angle interpolation stuff
   azimuth_deg_array = np.flip(np.arange(361) + 90.)
   azimuth_deg_array[azimuth_deg_array >= 361] = azimuth_deg_array[azimuth_deg_array >= 361] - 361
   zenith_angle_deg_array = np.arange(91)
   
   repetitions = 91
   ang_repeats_array = np.tile(azimuth_deg_array, (repetitions, 1))
   az_ang_repeats_array_flat = ang_repeats_array.flatten()
   
   repetitions = 361
   ang_repeats_array = np.tile(zenith_angle_deg_array, (repetitions, 1))
   
   zenith_angle_repeats_array_flat = ang_repeats_array.flatten('F')
   zenith_angle_repeats_array_flat_rad = zenith_angle_repeats_array_flat / 180. * np.pi
   az_ang_repeats_array_flat_rad = (az_ang_repeats_array_flat) / 180. * np.pi
       
   hpx_pix_num_array = np.arange(npix)
   hpx_angles_rad = hp.pix2ang(nside,hpx_pix_num_array) 
   hpx_angles_rad_zenith_angle = hpx_angles_rad[0]
   hpx_angles_rad_azimuth = hpx_angles_rad[1]
   
   test_n_freqs = 151

   E_phi_element_1_cube_slice = E_phi_element_1_cube[:,:,0:0+test_n_freqs]
   flattened_size = int(E_phi_element_1_cube_slice.shape[0]*E_phi_element_1_cube_slice.shape[1])
   E_phi_element_1_cube_slice = E_phi_element_1_cube_slice.transpose([1,0,2])
   E_phi_element_1_cube_slice_flat = E_phi_element_1_cube_slice.reshape(flattened_size,E_phi_element_1_cube_slice.shape[2])
   regridded_to_hpx_E_phi_element_1_complex = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_phi_element_1_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')

   E_phi_element_2_cube_slice = E_phi_element_2_cube[:,:,0:0+test_n_freqs]
   E_phi_element_2_cube_slice = E_phi_element_2_cube_slice.transpose([1,0,2])
   E_phi_element_2_cube_slice_flat = E_phi_element_2_cube_slice.reshape(flattened_size,E_phi_element_2_cube_slice.shape[2])
   regridded_to_hpx_E_phi_element_2_complex = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_phi_element_2_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')

   E_theta_element_1_cube_slice = E_theta_element_1_cube[:,:,0:0+test_n_freqs]
   E_theta_element_1_cube_slice = E_theta_element_1_cube_slice.transpose([1,0,2])
   E_theta_element_1_cube_slice_flat = E_theta_element_1_cube_slice.reshape(flattened_size,E_theta_element_1_cube_slice.shape[2])
   regridded_to_hpx_E_theta_element_1_complex = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_theta_element_1_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
   
   E_theta_element_2_cube_slice = E_theta_element_2_cube[:,:,0:0+test_n_freqs]
   E_theta_element_2_cube_slice = E_theta_element_2_cube_slice.transpose([1,0,2])
   E_theta_element_2_cube_slice_flat = E_theta_element_2_cube_slice.reshape(flattened_size,E_theta_element_2_cube_slice.shape[2])
   regridded_to_hpx_E_theta_element_2_complex = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), E_theta_element_2_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method='cubic')
   

   
   #element 1
   complex_beam_cube_1 = np.empty((npix,2,test_n_freqs), dtype=complex)
   complex_beam_cube_1[:,0,:] = regridded_to_hpx_E_theta_element_1_complex
   complex_beam_cube_1[:,1,:] = regridded_to_hpx_E_phi_element_1_complex

   #element 2
   complex_beam_cube_2 = np.empty((npix,2,test_n_freqs), dtype=complex)
   complex_beam_cube_2[:,0,:] = regridded_to_hpx_E_theta_element_2_complex
   complex_beam_cube_2[:,1,:] = regridded_to_hpx_E_phi_element_2_complex

   ####sanity check power pattern:
   #power_pattern_1 = np.abs(complex_beam_cube_1[:,0,0])**2 + np.abs(complex_beam_cube_1[:,1,0])**2
   #rotated_power_pattern_1  = r_beam_dec.rotate_map(power_pattern_1)
   #plt.clf()
   #map_title="rotated beam sim"
   #hp.orthview(map=rotated_power_pattern_1,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
   #fig_name="check1_complex_power_pattern_sitara.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name)

   #Now we have beams, need a sky!
   gsm = GlobalSkyModel()
   gsm_map_512_cube = gsm.generate(freq_MHz_array[0:test_n_freqs])
   gsm_map_cube = hp.ud_grade(gsm_map_512_cube,nside)
   
   #(don't)rotate gsm - rotate the beam, easier for SITARA since there are only 2
   #dec_rotate_gsm = 90. - float(mwa_latitude_ephem)
   #r_gsm_dec = hp.Rotator(rot=[0,dec_rotate_gsm], coord=['C', 'C'], deg=True)
   #r_gsm_ra = hp.Rotator(rot=[-lst_deg,0], coord=['C', 'C'], deg=True)
   ##convert to celestial coords
   #r_gsm_C = hp.Rotator(coord=['G','C'])
   #rotated_gsm_C = r_gsm_C.rotate_map(gsm_map)
   ##rotate the sky insted of the beams (cause beams are a cube)
   #rotated_gsm_C_dec = r_gsm_dec.rotate_map(rotated_gsm_C)
   #rotated_gsm_C_dec_ra = r_gsm_ra.rotate_map(rotated_gsm_C_dec)
   #unity_sky_repeats_array = gsm_map_cube * 0 + unity_sky_value
   #unity_sky_repeats_array = np.transpose(unity_sky_repeats_array)
   
   unity_sky_repeats_array = np.zeros((npix,test_n_freqs)) + unity_sky_value
   
   power_pattern_cube = np.einsum('ij...,ij...->i...', complex_beam_cube_1, np.conj(complex_beam_cube_2))
   power_pattern_cube = np.nan_to_num(power_pattern_cube)
   unity_sky_beam_cube = np.einsum('ij,ij->ij',unity_sky_repeats_array, power_pattern_cube)
   
   ######sanity check:
   #plt.clf()
   #map_title=""
   ######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
   #hp.mollview(map=unity_sky_beam_cube[:,2],title=map_title)
   #fig_name="check3_complex_power_pattern_sitara.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print("saved %s" % fig_name) 

   unity_sky_beam_sum_array = np.einsum('ij->j',unity_sky_beam_cube)
   power_pattern_cube_mag = np.abs(power_pattern_cube)
   power_pattern_cube_mag_sum_array = np.einsum('ij->j',power_pattern_cube_mag)
   visibility_array = unity_sky_beam_sum_array / power_pattern_cube_mag_sum_array


   #plot uniform response vs freq
   plt.clf()
   map_title="uniform response sitara"
   plt.plot(freq_MHz_array[0:test_n_freqs],visibility_array)
   plt.ylabel("Real part of vis")
   plt.xlabel("Frequency (MHz)")
   fig_name="uniform_response_from_complex_beams_sitara.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
  
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
   if 'single_point' in signal_type_list:
       concat_output_name_base_X += '_SP'
       concat_output_name_base_Y += '_SP'
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
        
         #get rid of unwanted tmp files in /tmp
         cmd = "rm -rf /tmp/tmp*"
         print(cmd)
         os.system(cmd)
         
      
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
         
         reprojected_gsm_im_Jy_per_pix_name =  "%s_%s_%0.3f_MHz_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz)

         
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
            s_21_array = plot_S21(nu_array=freq_MHz_array)
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
                  #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
               else:
                  beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                  #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
            else:
                  beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
                  
    
            #power_pattern_average_interp_sin_im_name = 'power_pattern_average_%s_%s_MHz_interp_sin.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_im_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_fits_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.fits' % (pol,int(freq_MHz))

            cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
            print(cmd)
            os.system(cmd)

            #cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname_no_cos_za)
            #print(cmd)
            #os.system(cmd)
                                         
            beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz)
            #beam_image_sin_projected_im_name_no_cos_za = 'beam_image_sin_projected_%s_%0.3f_MHz_no_cos_za.im' % (pol,freq_MHz)
            beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.im' % (pol,freq_MHz)
            beam_image_sin_projected_regrid_gsm_fits_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid.fits' % (pol,freq_MHz)
            #beam_image_sin_projected_regrid_gsm_im_name_no_cos_za =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid_no_cos_za.im' % (pol,freq_MHz)
            #beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_regrid_no_cos_za.fits' % (pol,freq_MHz)
           
           
            cmd = "rm -rf %s %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name, beam_image_sin_projected_puthd_fits_name)
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
            print(cmd)
            os.system(cmd)
            
            #cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname_no_cos_za,beam_image_sin_projected_im_name_no_cos_za)
            #print(cmd)
            #os.system(cmd)
            
            #put in the correct ra in the header (ra = lst for zenith) 
            #puthd in="$beam/crval1" value=$lst_degs
            cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name,lst_deg)
            print(cmd)
            os.system(cmd)         

            #cmd = 'puthd in="%s/crval1" value=%0.4f,"degrees"' % (beam_image_sin_projected_im_name_no_cos_za,lst_deg)
            #print(cmd)
            #os.system(cmd) 
            
            #write out as a fits file to check header
            cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name)
            print(cmd)
            os.system(cmd)
            
            #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
            cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)   

            #cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name_no_cos_za,beam_image_sin_projected_regrid_gsm_im_name_no_cos_za,reprojected_gsm_im_name)
            #print(cmd)
            #os.system(cmd)
                        
         
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
               
               reprojected_gsm_im_Jy_per_pix_name_fine_chan =  "%s_%s_%0.3f_MHz_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz_fine_chan)

               
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
                        beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz_fine_chan,'xx')
                        #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                     else:
                        beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz_fine_chan,'yy')
                        #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
                  else:
                        beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
                                                    
                  cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
                  print(cmd)
                  os.system(cmd)
                  
                  time.sleep(1)
                  
                  #if not EDA2_data:
                  #   if use_analytic_beam:
                  #      if pol=='X':
                  #         alt_beam_image_sin_projected_fitsname = "model_%s_MHz_%s.fits" % (int(freq_MHz),'xx')
                  #      else:
                  #         alt_beam_image_sin_projected_fitsname = "model_%s_MHz_%s.fits" % (int(freq_MHz),'yy')
                  #   cmd = "cp %s %s " % (beam_image_sin_projected_fitsname,alt_beam_image_sin_projected_fitsname)
                  #   print(cmd)
                  #   os.system(cmd)
                  
                  beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz_fine_chan)
                  
                  cmd = "rm -rf %s " % (beam_image_sin_projected_im_name)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
                  print(cmd)
                  os.system(cmd)
                  
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

            #cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_regrid_gsm_im_name_no_cos_za,beam_image_sin_projected_regrid_gsm_fits_name_no_cos_za)
            #print(cmd)
            #os.system(cmd)
              
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
            
            reprojected_gsm_global_im_Jy_per_pix_name =  "%s_DG_%s_%0.3f_MHz_%s_pol_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz,pol)
            reprojected_gsm_angular_im_Jy_per_pix_name =  "%s_DA_%s_%0.3f_MHz_%s_pol_reproj_Jy_pix.im" % (sky_model,date_time_string,freq_MHz,pol)
            
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
            if 'single_point' in signal_type_list:
               model_vis_name_base += '_SP'
               out_vis_name = model_vis_name_base + '.vis'
               
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
               
               pointing_dec_SP = "-26.70331940"
               array_ant_locations_filename_255 = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_255_NEU.ant'
               #point_jack.source
               #flux,dra,ddec,bmaj,bmin,bpa,iflux,ipa,vflux
               #     1.0000    3600.0000    -4668.05016    0.0000    0.0000    0.0000    0.0000    0.0000
               #day_frac_plus_one_SP = float(lst)/24. + 1
               #miriad_uvgen_time_string_SP = '00JUN%1.3f' % day_frac_plus_one_SP
               cmd = "uvgen source=$MIRCAT/point_jack.source ant='%s' baseunit=-3.33564 corr='32,1,0,0.93' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (array_ant_locations_filename_255,miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec_SP, harange_string, out_vis_name)
               print(cmd)
               os.system(cmd)
               
               base_vis_name = out_vis_name
               
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
            #if do_cleanup_images_and_vis:
             #  cleanup_images_and_vis(array_label,lst,freq_MHz,pol)
      
      ###concat vis both pols for each freq
      #Don't need to do this concat anymore....
      #for pol in pol_list:
      #   output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
      #   output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
      #   output_concat_vis_pyuvdata_name_sub = "%s_concat_freqs_sub.vis" % model_vis_name_base
      #   output_concat_uvfits_pyuvdata_name_sub = "%s_concat_freqs_sub.uvfits" % model_vis_name_base
      #   
      #   cmd = "rm -rf %s %s %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name,output_concat_vis_pyuvdata_name_sub,output_concat_uvfits_pyuvdata_name_sub)
      #   print(cmd)
      #   os.system(cmd)
      #   
      #   print(model_vis_uvfits_list_X)
      #   
      #   if pol=='X':
      #       concat_uvfits(model_vis_uvfits_list_X,output_concat_uvfits_pyuvdata_name)
      #       if do_image_and_subtr_from_simulated_data:
      #          concat_uvfits(model_vis_uvfits_list_X_sub,output_concat_uvfits_pyuvdata_name_sub)
      #   else:
      #       concat_uvfits(model_vis_uvfits_list_Y,output_concat_uvfits_pyuvdata_name)
      #       if do_image_and_subtr_from_simulated_data:
      #          concat_uvfits(model_vis_uvfits_list_Y_sub,output_concat_uvfits_pyuvdata_name_sub)
      
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
      #   if pol=='X':
      #      #model_sky_uvfits_list_X_lsts.append(output_concat_model_uvfits_name)
      ##      #model_global_signal_uvfits_list_X_lsts.append(global_output_concat_model_uvfits_name)
      #      #model_noise_uvfits_list_X_lsts.append(noise_output_concat_model_uvfits_name)
      #      model_vis_uvfits_list_X_lsts.append(output_concat_uvfits_pyuvdata_name)
      #      if do_image_and_subtr_from_simulated_data:
      #         model_vis_uvfits_list_X_lsts_sub.append(output_concat_uvfits_pyuvdata_name_sub)
      #   else:   
      #      #model_sky_uvfits_list_Y_lsts.append(output_concat_model_uvfits_name)
      #      #model_global_signal_uvfits_list_Y_lsts.append(global_output_concat_model_uvfits_name)      
      #      #model_noise_uvfits_list_Y_lsts.append(noise_output_concat_model_uvfits_name)
      #      model_vis_uvfits_list_Y_lsts.append(output_concat_uvfits_pyuvdata_name)
      #      if do_image_and_subtr_from_simulated_data:
      #         model_vis_uvfits_list_Y_lsts_sub.append(output_concat_uvfits_pyuvdata_name_sub)
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
   #cmd = "rm -rf %s %s %s %s" % (output_concat_vis_pyuvdata_name_lsts_X,output_concat_vis_pyuvdata_name_lsts_Y,output_concat_vis_pyuvdata_name_lsts_X_sub,output_concat_vis_pyuvdata_name_lsts_Y_sub)
   #print(cmd)
   #os.system(cmd)
   #   
   #for pol in pol_list:
   #   if pol=='X':
   #       concat_uvfits(model_vis_uvfits_list_X_lsts,output_concat_uvfits_pyuvdata_name_lsts_X)
   #       if do_image_and_subtr_from_simulated_data:
   #          concat_uvfits(model_vis_uvfits_list_X_lsts_sub,output_concat_uvfits_pyuvdata_name_lsts_X_sub)
   #   else:
   #       concat_uvfits(model_vis_uvfits_list_Y_lsts,output_concat_uvfits_pyuvdata_name_lsts_Y)
   #       if do_image_and_subtr_from_simulated_data:
   #          concat_uvfits(model_vis_uvfits_list_Y_lsts_sub,output_concat_uvfits_pyuvdata_name_lsts_Y_sub)
       
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
         #plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
         #plt.xlabel("Frequency (MHz)")
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
         plt.xlabel("Frequency (MHz)")
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
   
def inspect_cross_auto_power(lst_hrs_list,freq_MHz_list,pol_list,signal_type_list,sky_model,array_label,baseline_length_thresh_lambda,poly_order,EDA2_data=False,EDA2_chan_list='None',n_obs_concat_list=[],wsclean=False,fast=False,woden=True,noise_coupling=True,baseline_number=0,unity=False,angular=False):
   n_fine_chans = 32
   if EDA2_data:
      n_edge_chans_omitted = 5 #two at start and 3 at end
      n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
   else:
      n_fine_chans_used = 1
   n_edge_chans_omitted = 0
   
   
         
   ant1, ant2 = decode_baseline(baseline_number)
   for pol in pol_list:
      print("inspecting cross and/or power(s) for baseline %s, ant %s and ant %s " % (int(baseline_number), ant1, ant2))
      cross_power_array = np.full(int(len(freq_MHz_list)*n_fine_chans_used),np.nan)
      auto_power_ant1_array = np.full(int(len(freq_MHz_list)*n_fine_chans_used),np.nan)
      auto_power_ant2_array = np.full(int(len(freq_MHz_list)*n_fine_chans_used),np.nan)
      fine_chan_array = np.full(int(len(freq_MHz_list)*n_fine_chans_used),np.nan)
      unity_sky_value_K_array = np.full(int(len(freq_MHz_list)*n_fine_chans_used),np.nan)
      
      cross_power_array_filename = 'cross_power_%s_%s.npy' % (ant1,ant2)
      auto_power_ant1_array_filename = 'auto_power_%s.npy' % (ant1)
      auto_power_ant2_array_filename = 'auto_power_%s.npy' % (ant2)
      fine_chan_array_filename = 'fine_chan_array_inspect_auto_cross_power.npy' 
      unity_sky_value_K_array_filename = 'unity_sky_value_K_array.npy' 
      
      
      #print(EDA2_obs_time_list)
      for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
         start_fine_index =  freq_MHz_index * n_fine_chans_used 
         end_fine_index = (freq_MHz_index * n_fine_chans_used) + n_fine_chans_used
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
         
            auto_power_ant1_chan_array,fine_chan_array_chan,unity_sky_value_K_array_chan = get_visibility_power_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list,sky_model=sky_model,array_label=array_label,EDA2_data=EDA2_data,EDA2_obs_time=EDA2_obs_time,EDA2_chan=EDA2_chan,n_obs_concat=n_obs_concat,wsclean=wsclean,fast=fast,woden=woden,noise_coupling=noise_coupling,baseline_number=baseline_number,unity=unity,angular=angular)
        
        
         auto_power_ant1_array[start_fine_index:end_fine_index] = auto_power_ant1_chan_array
         fine_chan_array[start_fine_index:end_fine_index] = fine_chan_array_chan
         unity_sky_value_K_array[start_fine_index:end_fine_index] = unity_sky_value_K_array_chan
         
         
      np.save(auto_power_ant1_array_filename,auto_power_ant1_array)
      np.save(fine_chan_array_filename,fine_chan_array)
     
      print(auto_power_ant1_array.shape)
      print(auto_power_ant1_array)
      print(fine_chan_array.shape)
      print(fine_chan_array)
 

      plt.clf()
      
      plt.plot(unity_sky_value_K_array)
      
      map_title="unity sky value K" 
      #plt.xlabel("Frequency (MHz)")
      plt.ylabel("unity sky value K")
      plt.legend(loc='lower right')
      #if EDA2_data:
         #plt.ylim([500, 5000])
      fig_name= "unity_sky_values_K.png" 
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
    
      plt.clf()
      
      plt.plot(fine_chan_array)
      
      map_title="Frequency" 
      #plt.xlabel("Frequency (MHz)")
      plt.ylabel("Freq MHz")
      plt.legend(loc='lower right')
      #if EDA2_data:
         #plt.ylim([500, 5000])
      fig_name= "fine_chan_array_plot.png" 
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
      #plt.clf()
      #
      #plt.plot(auto_power_ant1_array)
      #
      #map_title="Cross power" 
      ##plt.xlabel("Frequency (MHz)")
      #plt.ylabel("Cross power (Jy)")
      #plt.legend(loc='lower right')
      ##if EDA2_data:
      #   #plt.ylim([500, 5000])
      #fig_name= "cross_power_baseline_%s_ants_%s_%s_pol_%s_noxlabel.png" % (baseline_number,ant1,ant2,pol)
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name) 
         
      plt.clf()
      
      plt.plot(fine_chan_array,auto_power_ant1_array)
      
      map_title="Cross power" 
      plt.xlabel("Frequency (MHz)")
      plt.ylabel("Cross power (Jy)")
      plt.legend(loc='lower right')
      #if EDA2_data:
         #plt.ylim([500, 5000])
      fig_name= "cross_power_baseline_%s_ants_%s_%s_pol_%s.png" % (baseline_number,ant1,ant2,pol)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print("saved %s" % fig_name) 
      
   
def get_visibility_power_from_uvfits(freq_MHz_list,freq_MHz_index,lst_hrs_list,pol,signal_type_list,sky_model,array_label,EDA2_data=False,EDA2_obs_time='None',EDA2_chan='None',n_obs_concat=1,wsclean=False,fast=False,woden=True,noise_coupling=True,baseline_number=0,unity=False,angular=False):
   if pol=='X':
      pol_index = 0
   elif pol=='Y':
      pol_index = 1
   else:
      print('pol %s not recognised' % pol)
      sys.exit()
   freq_MHz = freq_MHz_list[freq_MHz_index]
   start_freq = freq_MHz_list[0]
   concat_output_name_base = "%s_%s_%s" % (array_label,pol,outbase_name)
   output_prefix = "%s" % (array_label)
   signal_type_postfix = ''
   n_fine_chans = 32
   if 'noise' in signal_type_list:
       signal_type_postfix += '_N'
       concat_output_name_base += '_N'
   if 'diffuse' in signal_type_list:
       signal_type_postfix += '_D_%s' % sky_model
       concat_output_name_base += '_D_%s' % sky_model
       if woden:
          type = "gsm"
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
       if woden:
          type = "EDGES_uniform"
   if 'gain_errors' in signal_type_list:
       signal_type_postfix += '_GE'
       concat_output_name_base += '_GE'
   
   n_lsts = len(lst_hrs_list)

   lst_hrs = lst_hrs_list[freq_MHz_index]
   lst_deg = (float(lst_hrs)/24.)*360.
   
   if fast:
      if EDA2_data:
         if unity:
            uvfits_filename = "%s/unity_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,centre_chan_index,EDA2_obs_time,pol)
         elif angular:
            uvfits_filename = "%s/angular_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,centre_chan_index,EDA2_obs_time,pol)
         else:
            uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat) 
      elif woden:
         if noise_coupling:
            uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
         else:
            uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
         obs_time_list = ['1']
      else:
         uvfits_filename = "%s_LST_%03d_%s_%0.3f_MHz%s.uvfits" % (output_prefix,lst_deg,pol,freq_MHz,signal_type_postfix)
         obs_time_list = ['1']
      #read the cal uvfits, extract real vis uu and vv
      print("%s" % uvfits_filename)
      hdulist = fits.open(uvfits_filename)
      #hdulist.info()
      #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print info_string
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      #print(uvtable_header)
      hdulist.close()

      visibilities_single = uvtable['DATA']
      visibilities_shape = visibilities_single.shape
      print("visibilities_shape")
      print(visibilities_shape)

      UU_s = uvtable['UU'][baseline_number]
      UU_m = UU_s * c   
      VV_s = uvtable['VV'][baseline_number]
      VV_m = VV_s * c

      baseline_length_m = np.sqrt(UU_m**2 + VV_m**2)
            
      print("baseline length is %0.3f m" % baseline_length_m)
      
      if EDA2_data:
         n_edge_chans_omitted = 5 #two at start and 3 at end
         n_fine_chans_used = n_fine_chans - n_edge_chans_omitted
      else:
         n_fine_chans_used = 1
         n_edge_chans_omitted = 0
       
      #only return the central n_fine_chans_used 
      if EDA2_data:
         fine_chan_index_array = np.arange(n_fine_chans_used) + 2
      else:
         fine_chan_index_array = np.asarray([0])
      
      centre_freq = float(freq_MHz)
      fine_chan_width_MHz = fine_chan_width_Hz/1000000.         
       
      vis_power_list = []
      freq_MHz_fine_chan_list = []
      unity_sky_value_list = []
      #now do for each fine chan:
      for fine_chan_index in fine_chan_index_array:
         fine_chan_index = int(fine_chan_index)
         baseline_length_array_lambda_sorted_cut_list = []

         if EDA2_data: #hack 1
            #data coming out of the TPMs is reversed by coarse chan so for 20200303_data (and 20200304), need to change the freq calculation
            #2020 paper:
            if not reverse_fine_chans:
               freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz 
            else:
               freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz 
               #freq_MHz_fine_chan = centre_freq - (fine_chan_index)*fine_chan_width_MHz
         else:
            freq_MHz_fine_chan = freq_MHz
         wavelength = 300./float(freq_MHz_fine_chan)
         
         jy_to_K = (wavelength**2) / (2. * k * 1.0e26) 
         unity_sky_value = 1. * jy_to_K
         unity_sky_value_list.append(unity_sky_value)
         #print("fine_chan index,MHz,wavelength")
         #print(fine_chan_index)
         #print(freq_MHz_fine_chan)
         #print(wavelength)
         
         freq_MHz_fine_chan_list.append(freq_MHz_fine_chan)

         for obs_time_fast in [EDA2_obs_time]:
            if EDA2_data:
               #uvfits_filename = "%s/wscal_chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,obs_time_fast)
               #this is old pre woden miriad stuff
               if unity:
                  uvfits_filename = "%s/unity_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,obs_time_fast,pol)
               elif angular:
                  uvfits_filename = "%s/angular_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,obs_time_fast,pol)
               else:
                  uvfits_filename = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time,n_obs_concat)
               #this is woden - havent run these yet
               #unity_uvfits_filename = "%s/woden_LST_%0.3f_unity_uniform_start_freq_%0.3f_band%02d.uvfits" % (EDA2_chan,lst_deg,start_freq,freq_MHz_index) 
               #angular_uvfits_filename = "%s/woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular_band%02d.uvfits" % (EDA2_chan,lst_deg,start_freq,pol,freq_MHz_index) 
            elif woden:
               if noise_coupling:
                  uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d_nc.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) #woden_LST_60.000_gsm_start_freq_50.000_band99.uvfits 
               elif unity:
                  uvfits_filename = "woden_LST_%0.3f_unity_uniform_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,start_freq,freq_MHz_index) 
               elif angular:
                  uvfits_filename = "woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular_band%02d.uvfits" % (lst_deg,start_freq,pol,freq_MHz_index) 
               else:
                  uvfits_filename = "woden_LST_%0.3f_%s_start_freq_%0.3f_band%02d.uvfits" % (lst_deg,type,start_freq,freq_MHz_index) 

            if unity or angular: 
               print("%s" % uvfits_filename)
               with fits.open(uvfits_filename) as hdulist:
                  uvtable = hdulist[0].data
                  uvtable_header = hdulist[0].header
                  visibilities_single = uvtable['DATA']
                  visibilities_shape = visibilities_single.shape
                  print(visibilities_shape)
                  
                  real_vis_data = visibilities_single[baseline_number,0,0,0,pol_index,0]
                  imag_vis_data = visibilities_single[baseline_number,0,0,0,pol_index,1]
                  weights_vis_data = visibilities_single[baseline_number,0,0,0,pol_index,2]
            
                  total_vis_power = np.sqrt(real_vis_data**2 + imag_vis_data**2)
                  print(total_vis_power)
                  vis_power_list.append(total_vis_power)
  
            else:
               #read the cal uvfits, extract real vis uu and vv
               print("%s" % uvfits_filename)
               with fits.open(uvfits_filename) as hdulist:
                  #hdulist.info()
                  #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
                  #print info_string
                  uvtable = hdulist[0].data
                  uvtable_header = hdulist[0].header
                  #print(uvtable_header)
                  #hdulist.close()
   
                  visibilities_single = uvtable['DATA']
                  visibilities_shape = visibilities_single.shape
                  #print("visibilities_shape")
                  #print(visibilities_single)
                  #print(visibilities_shape)
                  
                  if wsclean:
                     real_vis_data = visibilities_single[baseline_number,0,0,0,fine_chan_index,pol_index,0]
                     imag_vis_data = visibilities_single[baseline_number,0,0,0,fine_chan_index,pol_index,1]
                     weights_vis_data = visibilities_single[baseline_number,0,0,0,fine_chan_index,pol_index,2]
                  else:
                     real_vis_data = visibilities_single[baseline_number,0,0,fine_chan_index,pol_index,0]
                     imag_vis_data = visibilities_single[baseline_number,0,0,fine_chan_index,pol_index,1]
                     weights_vis_data = visibilities_single[baseline_number,0,0,fine_chan_index,pol_index,2]
            
                  total_vis_power = np.sqrt(real_vis_data**2 + imag_vis_data**2)
                  print(total_vis_power)
                  vis_power_list.append(total_vis_power)

               
      vis_power_array = np.asarray(vis_power_list)
      freq_MHz_fine_chan_array = np.asarray(freq_MHz_fine_chan_list)
      unity_sky_value_array = np.asarray(unity_sky_value_list) 
      #print(vis_power_array)
      #print(vis_power_array.shape)
      #print(freq_MHz_fine_chan_array)
      #print(freq_MHz_fine_chan_array.shape)
      
      #plt.clf()
      #plt.plot(vis_power_array)
      #map_title="Cross power" 
      ##plt.xlabel("Frequency (MHz)")
      #plt.ylabel("Cross power (Jy)")
      #plt.legend(loc='lower right')
      ##if EDA2_data:
      #   #plt.ylim([500, 5000])
      #fig_name= "cross_power_baseline_%s_pol_%s_noxlabel_edachan_%s.png" % (baseline_number,pol,EDA2_chan)
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name) 
      
      return vis_power_array,freq_MHz_fine_chan_array,unity_sky_value_array
      
def calibrate_eda2_data_time_av(EDA2_chan_list,obs_type='night',lst_list=[],pol_list=[],sky_model_im_name='',n_obs_concat_list=[],concat=False,wsclean=False,plot_cal=False,uv_cutoff=0,per_chan_cal=False,EDA2_data=True,sim_only_EDA2=[]):
   if len(sim_only_EDA2)!=0:
      sim_only=True
   #think about how to use coherence:
   #coherence = cross12_avg/(np.sqrt(auto11_avg*auto22_avg))
   print("averaging EDA2 obs in time before calibration")
   #specify uv_cutoff in wavelengths, convert to m for 'calibrate'
   #pol = pol_list[0]
   gsm  = GlobalSkyModel()
   wsclean_imsize = '512'
   wsclean_scale = '900asec'  
   
   #if not sim_only:
     
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
       centre_freq = float(freq_MHz)
       fine_chan_width_MHz = fine_chan_width_Hz/1000000.  
       if uv_cutoff!=0:
          uv_cutoff_m = uv_cutoff * wavelength
       lst = lst_list[EDA2_chan_index]
       lst_deg = (float(lst)/24)*360.
       if len(EDA2_chan_list)==1:
          obs_time_list = EDA2_obs_time_list_each_chan[chan_num]
       else:
          obs_time_list = EDA2_obs_time_list_each_chan[EDA2_chan_index]
       first_obstime = obs_time_list[0]
       
       #guard against cases where there are no data for that channel
       if first_obstime==0:
          continue
       else:
          pass
       
       
       #Don't do the time averaging yet, try calibrating each obs individually so you can weed out the bad ones
       #No need for a multi-freq model cube for each obs, just use centre freq 
       wsclean_cal_ms_name_list = []
       for EDA2_obs_time_index,EDA2_obs_time in enumerate(obs_time_list):
          print(EDA2_obs_time)
          gsm_hpx_fits_name_chan = "%s/%s_map_LST_%0.3f_%0.3f_MHz_hpx.fits" % (EDA2_chan,sky_model,lst_deg,freq_MHz)
          unity_hpx_fits_name_chan = "%s/unity_sky_map_LST_%0.3f_%0.3f_MHz_hpx.fits" % (EDA2_chan,lst_deg,freq_MHz)
       
          reprojected_to_wsclean_gsm_prefix_chan = "%s/%s_map_LST_%0.3f_%0.3f_MHz_hpx_reprojected_wsclean" % (EDA2_chan,sky_model,lst_deg,freq_MHz)
          reprojected_to_wsclean_gsm_fitsname_chan = "%s.fits" % (reprojected_to_wsclean_gsm_prefix_chan)
          reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_chan = "%s_Jy_per_pix.fits" % (reprojected_to_wsclean_gsm_prefix_chan)
          reprojected_to_wsclean_gsm_im_name_Jy_per_pix_chan = "%s_map_LST_%0.3f_%0.3f_MHz_hpx_reproj_Jy_per_pix.im" % (sky_model,lst_deg,freq_MHz)

          # remake the gsm files it is not hard and save as hpx fits (so dont need to run simulate anymore)
       
          gsm_map = gsm.generate(freq_MHz)
       
          hp.write_map(gsm_hpx_fits_name_chan,gsm_map,coord='G',overwrite=True)
          print("wrote %s" % gsm_hpx_fits_name_chan)
                        
          hdu_gsm_chan = fits.open(gsm_hpx_fits_name_chan)[1]
  
          uncal_ms_image_prefix = "uncal_chan_%s_%s_ms" % (EDA2_chan,EDA2_obs_time)
          uncal_ms_image_name = "%s-image.fits" % uncal_ms_image_prefix

          print("cal using wsclean predict and calibrate / CASA bandpass")
          ms_name = "%s/%s_%s_eda2_ch32_ant256_midday_avg8140.ms" % (EDA2_chan,EDA2_obs_time[0:8],EDA2_obs_time[9:15])
          
          #make a wsclean image of the uncalibrated ms just to get an image header to reproject to:
          cmd = "wsclean -name %s -size %s %s -scale %s -pol xx -data-column DATA %s " % (uncal_ms_image_prefix,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
          print(cmd)
          os.system(cmd)
       
          ###
          if os.path.isfile(uncal_ms_image_name) and os.access(uncal_ms_image_name, os.R_OK):
             hdulist = fits.open(uncal_ms_image_name)
          else:
             print("Either file %s is missing or is not readable" % uncal_ms_image_name)
             #continue        

       
          data=hdulist[0].data[0,0,:,:]
          new_header=hdulist[0].header
       
          pix_size_deg = float(new_header['CDELT1']) #* 1.5      #hack - try a bigger pix area and check the effect
          pix_area_deg_sq = pix_size_deg*pix_size_deg
          pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
       
          del new_header[8]
          del new_header[8]
          del new_header['history']
          new_header['bunit'] ='Jy/Pixel'                                  
          #print(pt_source_header)
          
          target_wcs = WCS(new_header)
          
          target_wcs=target_wcs.dropaxis(2)
          target_wcs=target_wcs.dropaxis(2)
                     
          reprojected_gsm_map_chan,footprint_chan = reproject_from_healpix(hdu_gsm_chan, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)

          #write the reprojected gsm maps to fits
          fits.writeto(reprojected_to_wsclean_gsm_fitsname_chan,reprojected_gsm_map_chan,clobber=True)
          #print new_header
          fits.update(reprojected_to_wsclean_gsm_fitsname_chan,reprojected_gsm_map_chan,header=new_header)
          print("wrote image %s" %  reprojected_to_wsclean_gsm_fitsname_chan)
               
          #model needs to be in Jy/pix
          scale_chan = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
          print("scale map by %s to get to Jy/pix" % scale_chan)
       

          #check the model image for non-finite values and scale to Jy per pix:
          with fits.open("%s" % (reprojected_to_wsclean_gsm_fitsname_chan)) as hdu_list_chan:
             data_chan = hdu_list_chan[0].data
          #replace nans with zeros
          data_new_chan = np.nan_to_num(data_chan)
          data_new_jy_per_pix_chan = data_new_chan * scale_chan
          
          #write out a new fits file
          fits.writeto("%s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_chan),data_new_jy_per_pix_chan,clobber=True)
          fits.update(reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_chan,data_new_jy_per_pix_chan,header=new_header)
          print("saved %s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_chan))
       
          for pol in ['X','Y']:
             if use_analytic_beam:
                if pol=='X':
                   beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
                   #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                else:
                   beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                   #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
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
             
             cmd = "rm -rf %s %s %s %s %s" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name,beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_chan)
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
          
             cmd = "fits in=%s out=%s op=xyin" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_chan,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_chan)
             print(cmd)
             os.system(cmd)
                          
             #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
             cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_chan)
             print(cmd)
             os.system(cmd)  
             
             #write out as a fits file for av sky temp calc
             cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name)
             print(cmd)
             os.system(cmd)
          
             #Now have a gsm and a beam. multiply 'em'
             apparent_sky_fits_name_prefix_chan = "apparent_sky_LST_%0.3f_%0.3f_MHz_wsclean" % (lst_deg,freq_MHz)
             apparent_sky_im_name_chan = "apparent_sky_LST_%0.3f_%0.3f_MHz_wsclean.im" % (lst_deg,freq_MHz)
             apparent_sky_fits_name_chan = "%s-%s%s-model.fits" % (apparent_sky_fits_name_prefix_chan,pol,pol)
             
             cmd = "rm -rf %s %s" % (apparent_sky_im_name_chan,apparent_sky_fits_name_chan)
             print(cmd)
             os.system(cmd)
         
             cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_chan,apparent_sky_im_name_chan)
             print(cmd)
             os.system(cmd)
             
             cmd = "fits in=%s out=%s op=xyout" % (apparent_sky_im_name_chan,apparent_sky_fits_name_chan)
             print(cmd)
             os.system(cmd) 
             
             print("wrote %s" % apparent_sky_fits_name_chan)
             
             #check the model image for non-finite values 
             with fits.open("%s" % (apparent_sky_fits_name_chan)) as hdu_list:
                data = hdu_list[0].data
             #replace nans with zeros
             data_new = np.nan_to_num(data)
             
             #write out a new fits file
             fits.writeto("%s" % (apparent_sky_fits_name_chan),data_new,clobber=True)
             fits.update(apparent_sky_fits_name_chan,data_new,header=new_header)
             print("saved %s" % (apparent_sky_fits_name_chan))
   
             cmd = "rm -rf %s" % (apparent_sky_im_name_chan)
             print(cmd)
             os.system(cmd)              
   
             #now predict and calibrate

          # predict a model
          cmd = "wsclean -predict -name %s -size %s %s -scale %s -pol xx,yy  %s " % (apparent_sky_fits_name_prefix_chan,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
          print(cmd)
          os.system(cmd)
          
          ##make  images to check 
          #cmd = "wsclean -name model_col_chan_%s_%s_ms -size %s %s -scale %s -pol xx -data-column MODEL_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
          #print(cmd)
          #os.system(cmd)
          
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
          
          #plot the sols and 
          if (os.path.isfile(gain_solutions_name)):
          
             wsclean_cal_ms_name_list.append(ms_name)
          
             if plot_cal:
                #Plot the cal solutions
                cmd = "aocal_plot.py  %s " % (gain_solutions_name)
                print(cmd)
                os.system(cmd)
                
             #write the calibrated uvfits file out ?
             #cmd = "fits in=%s out=%s op=uvout" % (miriad_vis_name,calibrated_uvfits_filename)
             #print(cmd)
             #os.system(cmd)
             
             #no need to apply sols
             #cmd = "applysolutions %s %s " % (ms_name,gain_solutions_name)
             #print(cmd)
             #os.system(cmd)
             
             ###make an image to check (both pols) 32 chans, skip for now, takes ages
             #cmd = "wsclean -name cal_chan_%s_%s_ms -size %s %s -auto-threshold 5 -scale %s -pol xx,yy -data-column CORRECTED_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
             #print(cmd)
             #os.system(cmd) 
    
          else:
             print("no cal solutions for %s" % (ms_name))
             continue

       number_of_good_obs = len(wsclean_cal_ms_name_list)
       print("number of good obs used in chan %s is %s" % (EDA2_chan,number_of_good_obs)) 
       
       #now average
       av_uvfits_name = "%s/av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,first_obstime,number_of_good_obs)
       av_ms_name = "%s/av_chan_%s_%s_plus_%s_obs.ms" % (EDA2_chan,EDA2_chan,first_obstime,number_of_good_obs)
       sum_ms_name = "%s/sum_chan_%s_%s_plus_%s_obs.ms" % (EDA2_chan,EDA2_chan,first_obstime,number_of_good_obs)
       calibrated_uvfits_filename_wsclean = "%s/cal_av_chan_%s_%s_plus_%s_obs.uvfits" % (EDA2_chan,EDA2_chan,first_obstime,number_of_good_obs)
       
       
       
       for ms_name_index,ms_name in enumerate(wsclean_cal_ms_name_list):    
          t = pt.table(ms_name, readonly=True, ack=False)
          print('Total rows in  %s  = %s' % (ms_name,str(t.nrows())))
          intTime = t.getcell("INTERVAL", 0)
          print('Integration time:\t%f sec' % (intTime))

          # open the antenna and spectral window subtables
          tant = pt.table(t.getkeyword('ANTENNA'), readonly=True, ack=False)
          tsp = pt.table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False)
          numChannels = len(tsp.getcell('CHAN_FREQ',0))
          print('Number of channels:\t%d' % (numChannels))
          print('Reference frequency:\t%5.2f MHz' % (tsp.getcell('REF_FREQUENCY',0)/1.e6))
          
          t_autos = t.query('ANTENNA1 = ANTENNA2') #autos
          t_cross = t.query('ANTENNA1 != ANTENNA2')
          t_all = t.query('ANTENNA1 != ANTENNA2 OR ANTENNA1 = ANTENNA2').DATA[0]
          #print(t_all)

          #This bit is just to learn about antenna tables so I can do the WODEN sims right on pawsey
          #ant_table = pt.table(ms_name+"::ANTENNA", readonly=True)
          #pos = ant_table.getcol("POSITION")
          #print(pos)
          #flag_row = ant_table.getcol("FLAG_ROW") 
          #print(flag_row)
          #sys.exit()
          

          #if this is the first observation create a copy of the the ms to be the sum ms
          if(ms_name_index==0):
             cmd = 'rm -rf %s %s' % (av_ms_name,sum_ms_name)
             print(cmd)
             os.system(cmd)
             
             pt.taql('select from %s where ANTENNA1 != ANTENNA2 OR ANTENNA1 = ANTENNA2 giving %s as plain' % (ms_name,sum_ms_name))

    
          #update the sum ms by adding the current ms data
          if(ms_name_index!=0):
             pt.taql('update %s t1, %s t2 set DATA = t2.DATA+t1.DATA ' % (sum_ms_name,ms_name))
          
          #otherwise gets stuck?
          t.unlock()
       #e.g Put CORRECTED_DATA from that MS into DATA column of this MS
       #It requires that both tables have the same number of rows
       #update this.ms, that.ms t2 set DATA = t2.CORRECTED_DATA!
       
       pt.taql('select from %s where ANTENNA1 != ANTENNA2 OR ANTENNA1 = ANTENNA2 giving %s as plain' % (sum_ms_name,av_ms_name))
       
       pt.taql('update %s t1 set DATA=t1.DATA/%s ' % (av_ms_name,float(number_of_good_obs)))
      
       cmd = 'rm -rf %s' % (sum_ms_name)
       print(cmd)
       os.system(cmd)
       
       #check:
       #t_av = pt.table(av_ms_name, readonly=True, ack=False)
       #t_av_data = t_av.query('ANTENNA1 != ANTENNA2 OR ANTENNA1 = ANTENNA2').DATA[0]
       #print(t_av_data)
             
       #we have an averaged ms

       ##########################
       ##########################
       apparent_unity_sky_im_name_fine_chan_list_X = []
       apparent_unity_sky_im_name_fine_chan_list_Y = []
       apparent_angular_sky_im_name_fine_chan_list_X = []
       apparent_angular_sky_im_name_fine_chan_list_Y = []
       max_gsm_list = []
       scale_fine_chan_list = []
       
       for fine_chan_index in range(0,32):
          #reversing or not reversing the order of the input images in the model appears to make zero difference to the result!
          #dont reverse chan order (2020 paper)
          if not reverse_fine_chans:
             freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
          else:
             #do reverse chan order (time av cal 2021 miriad)
             freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz
          
          wavelength_fine_chan = 300./float(freq_MHz_fine_chan)
          
          #(these get made in here now, no need to run simulate
          gsm_hpx_fits_name_fine_chan = "%s/%s_map_LST_%0.3f_%0.3f_MHz_hpx.fits" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan)
          unity_hpx_fits_name_fine_chan = "%s/unity_sky_map_LST_%0.3f_%0.3f_MHz_hpx.fits" % (EDA2_chan,lst_deg,centre_freq)
          
       
          reprojected_to_wsclean_gsm_prefix_fine_chan = "%s/%s_map_LST_%0.3f_%0.3f_MHz_hpx_reprojected_wsclean" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan)
          reprojected_to_wsclean_gsm_fitsname_fine_chan = "%s.fits" % (reprojected_to_wsclean_gsm_prefix_fine_chan)
          reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan = "%s_Jy_per_pix.fits" % (reprojected_to_wsclean_gsm_prefix_fine_chan)
          reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan = "%s_map_LST_%0.3f_%0.3f_MHz_hpx_reproj_Jy_per_pix.im" % (sky_model,lst_deg,freq_MHz_fine_chan)


          #going to make unity sky ones too (so dont need to run simulate)
          reprojected_to_wsclean_unity_prefix_fine_chan = "%s/unity_map_LST_%0.3f_%0.3f_MHz_hpx_reprojected_wsclean" % (EDA2_chan,lst_deg,freq_MHz_fine_chan)
          reprojected_to_wsclean_unity_fitsname_fine_chan = "%s.fits" % (reprojected_to_wsclean_unity_prefix_fine_chan)
          reprojected_to_wsclean_unity_fitsname_Jy_per_pix_fine_chan = "%s_Jy_per_pix.fits" % (reprojected_to_wsclean_unity_prefix_fine_chan)
          reprojected_to_wsclean_unity_im_name_Jy_per_pix_fine_chan = "unity_map_LST_%0.3f_%0.3f_MHz_hpx_reproj_Jy_per_pix.im" % (lst_deg,freq_MHz_fine_chan)

          jy_to_K = (wavelength_fine_chan**2) / (2. * k * 1.0e26) 
          unity_sky_value = 1. #* jy_to_K (i think this should be divide to get to Jy from K?)
          
          gsm_map = gsm.generate(freq_MHz_fine_chan) 
          unity_map = gsm_map * 0. + unity_sky_value
          
          #here get the true predicted global monolpole (not beam averaged)
          sky_global_monopole_temp_filename = "%s/sky_mono_temp_%s_LST_%0.3f_%0.3f_MHz.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan)
          monopole_temp = np.mean(gsm_map)   
          
          monopole_temp_array = np.asarray([monopole_temp])
          np.save(sky_global_monopole_temp_filename,monopole_temp_array)
          print("saved %s" % (sky_global_monopole_temp_filename))
          print("only saving monopole temps, see line 11262 and 11646" )
          
          #HACK TO ONLY SAVE monopole temp
          #continue
          
          max_gsm = np.max(gsm_map)
          max_gsm_list.append(max_gsm)
          
          hp.write_map(gsm_hpx_fits_name_fine_chan,gsm_map,coord='G',overwrite=True)
          print("wrote %s" % gsm_hpx_fits_name_fine_chan)
                    
          hp.write_map(unity_hpx_fits_name_fine_chan,unity_map,coord='G',overwrite=True)
          print("wrote %s" % unity_hpx_fits_name_fine_chan)         
          
          hdu_gsm_fine_chan = fits.open(gsm_hpx_fits_name_fine_chan)[1]
          hdu_unity_fine_chan = fits.open(unity_hpx_fits_name_fine_chan)[1]
  
          uncal_ms_image_prefix = "uncal_chan_%s_%s_ms" % (EDA2_chan,first_obstime)
          uncal_ms_image_name = "%s-image.fits" % uncal_ms_image_prefix
          wsclean_imsize = '512'
          wsclean_scale = '900asec'
          
          print("cal using wsclean predict and calibrate / CASA bandpass")
          
          #Cant import the EDA sim uvfits file into ms - too many antennas - what if you make an mwa 32T simulated uvfits file instead!? then import to ms and image at low res (do this in simulate())
          
          #make a wsclean image of the uncalibrated ms just to get an image header to reproject to:
          cmd = "wsclean -name %s -size %s %s -scale %s -pol xx -data-column DATA %s " % (uncal_ms_image_prefix,wsclean_imsize,wsclean_imsize,wsclean_scale,av_ms_name)
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
          
          pix_size_deg = float(new_header['CDELT1']) #* 1.5      #hack - try a bigger pix area and check the effect
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
                     
          reprojected_gsm_map_fine_chan,footprint_fine_chan = reproject_from_healpix(hdu_gsm_fine_chan, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
          reprojected_unity_map_fine_chan,footprint_fine_chan = reproject_from_healpix(hdu_unity_fine_chan, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear',field=0)
          
          
          #hdu_gsm_fine_chan.close()
          #unity_gsm_fine_chan.close()
          
          #write the reprojected gsm maps to fits
          pyfits.writeto(reprojected_to_wsclean_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,clobber=True)
          #print new_header
          pyfits.update(reprojected_to_wsclean_gsm_fitsname_fine_chan,reprojected_gsm_map_fine_chan,header=new_header)
          print("wrote image %s" %  reprojected_to_wsclean_gsm_fitsname_fine_chan)
          
          #untiy (miriad)
          pyfits.writeto(reprojected_to_wsclean_unity_fitsname_fine_chan,reprojected_unity_map_fine_chan,clobber=True)
          #print new_header
          pyfits.update(reprojected_to_wsclean_unity_fitsname_fine_chan,reprojected_unity_map_fine_chan,header=new_header)
          print("wrote image %s" %  reprojected_to_wsclean_unity_fitsname_fine_chan)
    
          
          #model needs to be in Jy/pix
          #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
       
          scale_fine_chan = (2. * k * 1.0e26 * pix_area_sr) / (wavelength_fine_chan**2)  
          print("scale map by %s to get to Jy/pix" % scale_fine_chan)
          
          scale_fine_chan_list.append(scale_fine_chan)

          #check the model image for non-finite values and scale to Jy per pix:
          with fits.open("%s" % (reprojected_to_wsclean_gsm_fitsname_fine_chan)) as hdu_list_fine_chan:
             data_fine_chan = hdu_list_fine_chan[0].data
          #replace nans with zeros
          data_new_fine_chan = np.nan_to_num(data_fine_chan)
          ##
          #hack for test:
        
          #scale_fine_chan_modified = scale_fine_chan * 1.2
          #scale_fine_chan_modified = scale_fine_chan * (1.2 - float(fine_chan_index)*0.0125)
          #print("scale_fine_chan_modified is ")
          #print(scale_fine_chan_modified)
          #print("scale_fine_chan is ")
          #print(scale_fine_chan)
          #This does make a difference - you can change the slope! So why does reversing the channels above not do anything?
          #data_new_jy_per_pix_fine_chan = data_new_fine_chan * scale_fine_chan_modified  # hack
          
          #hack to remove scaling to jy from brightness temp
          
          data_new_jy_per_pix_fine_chan = data_new_fine_chan * scale_fine_chan 
          
          #write out a new fits file
          fits.writeto("%s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan),data_new_jy_per_pix_fine_chan,clobber=True)
          pyfits.update(reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan,data_new_jy_per_pix_fine_chan,header=new_header)
          print("saved %s" % (reprojected_to_wsclean_gsm_fitsname_Jy_per_pix_fine_chan))
          
          #unity miriad
          #check the model image for non-finite values and scale to Jy per pix:
          with fits.open("%s" % (reprojected_to_wsclean_unity_fitsname_fine_chan)) as hdu_list_fine_chan:
             data_fine_chan = hdu_list_fine_chan[0].data
          #replace nans with zeros
          data_new_fine_chan = np.nan_to_num(data_fine_chan)
          data_new_jy_per_pix_fine_chan = data_new_fine_chan * scale_fine_chan
          
          #write out a new fits file
          fits.writeto("%s" % (reprojected_to_wsclean_unity_fitsname_Jy_per_pix_fine_chan),data_new_jy_per_pix_fine_chan,clobber=True)
          pyfits.update(reprojected_to_wsclean_unity_fitsname_Jy_per_pix_fine_chan,data_new_jy_per_pix_fine_chan,header=new_header)
          print("saved %s" % (reprojected_to_wsclean_unity_fitsname_Jy_per_pix_fine_chan))
          
          for pol in ['X','Y']:
             if use_analytic_beam:
                if pol=='X':
                   beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz_fine_chan,'xx')
                   #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                else:
                   beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz_fine_chan,'yy')
                   #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
             else:
                beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
       
             cmd = "cp %s%s . " % (beam_image_dir,beam_image_sin_projected_fitsname)
             print(cmd)
             os.system(cmd)
             
             #Need to regrid the beam to the reproject gsm
             beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%0.3f_MHz.im' % (pol,freq_MHz_fine_chan)
             beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%0.3f_MHz_puthd.fits' % (pol,freq_MHz_fine_chan)
             beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.im' % (pol,freq_MHz_fine_chan)
             beam_image_sin_projected_regrid_gsm_fits_name =  'beam_image_sin_projected_%s_%0.3f_MHz_gsm_wsclean_regrid.fits' % (pol,freq_MHz_fine_chan)
             
             cmd = "rm -rf %s %s %s %s %s %s" % (beam_image_sin_projected_im_name,beam_image_sin_projected_puthd_fits_name,beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan,reprojected_to_wsclean_unity_im_name_Jy_per_pix_fine_chan)
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

             #unity miriad
             cmd = "fits in=%s out=%s op=xyin" % (reprojected_to_wsclean_unity_fitsname_Jy_per_pix_fine_chan,reprojected_to_wsclean_unity_im_name_Jy_per_pix_fine_chan)
             print(cmd)
             os.system(cmd)
                          
             #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
             cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan)
             print(cmd)
             os.system(cmd)  
             
             #write out as a fits file for av sky temp calc
             cmd = "fits in=%s out=%s op=xyout" % (beam_image_sin_projected_regrid_gsm_im_name,beam_image_sin_projected_regrid_gsm_fits_name)
             print(cmd)
             os.system(cmd)
             
             
             #Now have a gsm and a beam. multiply 'em'
             apparent_sky_fits_name_prefix_fine_chan = "apparent_sky_LST_%0.3f_%0.3f_MHz_wsclean" % (lst_deg,freq_MHz)
             apparent_sky_im_name_fine_chan = "apparent_sky_LST_%0.3f_%0.3f_MHz_wsclean-%04d.im" % (lst_deg,freq_MHz,fine_chan_index)
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
                
             #get the correct input average sky temp value (before removing nans)
             sky_averaged_temp_cal_input_filename = "%s/sky_av_input_cal_%s_LST_%0.3f_%0.3f_MHz_pol_%s.npy" % (EDA2_chan,sky_model,lst_deg,freq_MHz_fine_chan,pol)
             with fits.open("%s" % (apparent_sky_fits_name_fine_chan)) as hdu_list:
                sky_with_beam = hdu_list[0].data
                sky_with_beam_K = sky_with_beam / scale_fine_chan
             with fits.open("%s" % (beam_image_sin_projected_regrid_gsm_fits_name)) as hdu_list:
                beam_data = hdu_list[0].data
             beam_weighted_av_sky = np.nansum(sky_with_beam_K) / np.nansum(beam_data)
             beam_weighted_av_sky_Jy = np.nansum(sky_with_beam) / np.nansum(beam_data)
             
             print("beam_weighted_av_sky_Jy is %0.4E" % beam_weighted_av_sky_Jy)
             
             #check the model image for non-finite values 
             with fits.open("%s" % (apparent_sky_fits_name_fine_chan)) as hdu_list:
                data = hdu_list[0].data
             #replace nans with zeros
             data_new = np.nan_to_num(data)
             
             #write out a new fits file
             fits.writeto("%s" % (apparent_sky_fits_name_fine_chan),data_new,clobber=True)
             pyfits.update(apparent_sky_fits_name_fine_chan,data_new,header=new_header)
             print("saved %s" % (apparent_sky_fits_name_fine_chan))
             

             
             print("beam_weighted_av_sky is %E K (%E Jy) from %s and %s" % (beam_weighted_av_sky,beam_weighted_av_sky_Jy,apparent_sky_fits_name_fine_chan,beam_image_sin_projected_regrid_gsm_fits_name))
             
             beam_weighted_av_sky_array = np.asarray([beam_weighted_av_sky])

             np.save(sky_averaged_temp_cal_input_filename,beam_weighted_av_sky_array)
             print("saved %s" % (sky_averaged_temp_cal_input_filename))
             
             ####unity miriad (and angular)
             #Now have a unity map and a beam. multiply 'em'
             #apparent_unity_sky_im_name_fine_chan = "u_%0.3f_%s.im" % (freq_MHz_fine_chan,pol)
             #apparent_angular_sky_im_name_fine_chan = "a_%0.3f_%s.im" % (freq_MHz_fine_chan,pol)
             apparent_unity_sky_im_name_fine_chan = "u_%0.3f_%0.3f_MHz-%04d-%s.im" % (lst_deg,freq_MHz,fine_chan_index,pol)
             apparent_angular_sky_im_name_fine_chan_no_beam = "a_nb_%0.3f_%0.3f_MHz_%04d_%s.im" % (lst_deg,freq_MHz,fine_chan_index,pol)
             apparent_angular_sky_im_name_fine_chan = "a_%0.3f_%0.3f_MHz_%04d_%s.im" % (lst_deg,freq_MHz,fine_chan_index,pol)
             
             if pol=='X':
                apparent_unity_sky_im_name_fine_chan_list_X.append(apparent_unity_sky_im_name_fine_chan)
                apparent_angular_sky_im_name_fine_chan_list_X.append(apparent_angular_sky_im_name_fine_chan)
             else:
                apparent_unity_sky_im_name_fine_chan_list_Y.append(apparent_unity_sky_im_name_fine_chan)
                apparent_angular_sky_im_name_fine_chan_list_Y.append(apparent_angular_sky_im_name_fine_chan)
             
             cmd = "rm -rf %s %s" % (apparent_unity_sky_im_name_fine_chan,apparent_angular_sky_im_name_fine_chan_no_beam)
             print(cmd)
             os.system(cmd)
   
             cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_to_wsclean_unity_im_name_Jy_per_pix_fine_chan,apparent_unity_sky_im_name_fine_chan)
             print(cmd)
             os.system(cmd)
            
             #originally I was subtracting the global value from apparent gsm sky to get angular only
             #What I actually need to do (as in WODEN sims) is to subtract the global value from the image 
             #unattenuated by the beam and then apply the beam!
             
             maths_sky_im_name_fine_chan_no_beam = "nb_LST_%0.3f_%0.3f_MHz_wsclean_%04d.im" % (lst_deg,freq_MHz,fine_chan_index)
             maths_beam_image_sin_projected_im_name = 'beam_%s_%0.3f_MHz.im' % (pol,freq_MHz_fine_chan)

             cmd = "rm -rf %s %s" % (maths_sky_im_name_fine_chan_no_beam,maths_beam_image_sin_projected_im_name)
             print(cmd)
             os.system(cmd)
             
             cmd = "cp -r %s %s" % (reprojected_to_wsclean_gsm_im_name_Jy_per_pix_fine_chan,maths_sky_im_name_fine_chan_no_beam)
             print(cmd)
             os.system(cmd)     

             cmd = "cp -r %s %s" % (beam_image_sin_projected_regrid_gsm_im_name,maths_beam_image_sin_projected_im_name)
             print(cmd)
             os.system(cmd)  
                          
             cmd = "maths exp=%s-%0.12f out=%s " % (maths_sky_im_name_fine_chan_no_beam,beam_weighted_av_sky_Jy,apparent_angular_sky_im_name_fine_chan_no_beam)
             print(cmd)
             os.system(cmd)

             cmd = "rm -rf %s" % (apparent_angular_sky_im_name_fine_chan)
             print(cmd)
             os.system(cmd)              
             
             cmd = "maths exp=%s*%s out=%s " % (maths_beam_image_sin_projected_im_name,apparent_angular_sky_im_name_fine_chan_no_beam,apparent_angular_sky_im_name_fine_chan)
             print(cmd)
             os.system(cmd)             

             cmd = "rm -rf %s %s" % (maths_sky_im_name_fine_chan_no_beam,maths_beam_image_sin_projected_im_name)
             print(cmd)
             os.system(cmd)      
             
             #will use imcat to make cubes for uvmodel
             
             
             #instead of using imcat, write out a new uvfits for each fine chan
       
             #
             #This next bit used to be outside of the fine chan loop,  moving it in and getting rid of multi-freq unity image (just use central chan) - do this for angular too TODO
             
             ####
             #stuff for FAST - uses uvfits files even though calibrating ms data
             #need to do for both pols
             
             EDA2_obs_time = obs_time_list[0]
             uvfits_filename = "%s/chan_%s_%s.uvfits" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
             uvfits_vis_filename = "%s/chan_%s_%s.vis" % (EDA2_chan,EDA2_chan,EDA2_obs_time)
             uvfits_filename_fine_chan = "%s/chan_%s_fine_%02d_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time)
             uvfits_vis_filename_fine_chan = "%s/chan_%s_fine_%02d_%s.vis" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time)
             
             unity_sky_uvfits_filename = "%s/unity_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol)
             unity_sky_vis_filename = "%s/unity_chan_%s_%02d_%s_pol_%s.vis" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol) 
             angular_sky_uvfits_filename = "%s/angular_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol)
             angular_sky_vis_filename = "%s/angular_chan_%s_fine_%02d_%s_pol_%s.vis" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol) 
             #also uvmodel the full gsm, then you can compare the data to the sims
             gsm_sky_uvfits_filename = "%s/gsm_chan_%s_fine_%02d_%s_pol_%s.uvfits" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol)
             gsm_sky_vis_filename = "%s/gsm_chan_%s_fine_%02d_%s_pol_%s.vis" % (EDA2_chan,EDA2_chan,fine_chan_index,EDA2_obs_time,pol) 
             
             
                       
             #apparent_unity_sky_im_cube_name = "apparent_unity_sky_LST_%0.3f_%0.3f_MHz_cube_pol_%s.im" % (lst_deg,freq_MHz,pol)    
             #apparent_angular_sky_im_cube_name = "apparent_angular_sky_LST_%0.3f_%0.3f_MHz_cube_pol_%s.im" % (lst_deg,freq_MHz,pol)
             #
             #if pol=='X':
             #   unity_model_im_list_string = ','.join(apparent_unity_sky_im_name_fine_chan_list_X)
             #   angular_model_im_list_string = ','.join(apparent_angular_sky_im_name_fine_chan_list_X)
             #else:
             #   unity_model_im_list_string = ','.join(apparent_unity_sky_im_name_fine_chan_list_Y)
             #   angular_model_im_list_string = ','.join(apparent_angular_sky_im_name_fine_chan_list_Y)
             #   
             #cmd = "rm -rf %s %s %s %s %s %s %s" % (unity_sky_uvfits_filename,unity_sky_vis_filename,uvfits_vis_filename,apparent_unity_sky_im_cube_name,angular_sky_uvfits_filename,angular_sky_vis_filename,apparent_angular_sky_im_cube_name)
             cmd = "rm -rf %s %s %s %s %s %s %s %s" % (unity_sky_uvfits_filename,unity_sky_vis_filename,uvfits_filename_fine_chan,uvfits_vis_filename_fine_chan,angular_sky_uvfits_filename,angular_sky_vis_filename,gsm_sky_uvfits_filename,gsm_sky_vis_filename)
             print(cmd)
             os.system(cmd)
             
             #cmd = "imcat in=%s out=%s options=relax" % (unity_model_im_list_string,apparent_unity_sky_im_cube_name)
             #print(cmd)
             #os.system(cmd)
             #
             ##imcat angular
             #cmd = "imcat in=%s out=%s options=relax" % (angular_model_im_list_string,apparent_angular_sky_im_cube_name)
             #print(cmd)
             #os.system(cmd)
             
             #cmd = "fits in=%s op=uvin out=%s" % (uvfits_filename,uvfits_vis_filename)
             #print(cmd)
             #os.system(cmd)
             
             #make a copy of the uvfits_filename so we can modify the header with the fine chan freq
             cmd = "cp -r %s %s" % (uvfits_filename,uvfits_filename_fine_chan)
             print(cmd)
             os.system(cmd)
             
             #modify the uvfits_filename_fine_chan header here
             with fits.open("%s" % (uvfits_filename_fine_chan),'update') as hdu_list:
                #data = hdu_list[0].data
                for hdu in hdu_list:
                   hdu.header['CRVAL4'] = freq_MHz_fine_chan * 1.e6
                   hdu.header['FREQ'] = freq_MHz_fine_chan * 1.e6
           
             ##check  
             #with fits.open("%s" % (uvfits_filename_fine_chan),readonly=True) as hdu_list:
             #   #data = hdu_list[0].data
             #   header = hdu_list[0].header
             #   print(header)
                  
             cmd = "fits in=%s op=uvin out=%s" % (uvfits_filename_fine_chan,uvfits_vis_filename_fine_chan)
             print(cmd)
             os.system(cmd)
            
             
             #cmd = "puthd in=%s/crval4 value=%0.3f" % (uvfits_vis_filename_fine_chan,freq_MHz_fine_chan)
             #print(cmd)
             #os.system(cmd)
             
             #apparent_unity_sky_im_name_centre_chan = "u_%0.3f_%s.im" % (centre_freq,pol)
             #apparent_angular_sky_im_name_centre_chan = "a_%0.3f_%s.im" % (centre_freq,pol)
             
             cmd = "uvmodel vis=%s model=%s options=replace out=%s" % (uvfits_vis_filename_fine_chan,apparent_unity_sky_im_name_fine_chan,unity_sky_vis_filename)
             print(cmd)
             os.system(cmd)
             
             #angular
             cmd = "uvmodel vis=%s model=%s options=replace out=%s" % (uvfits_vis_filename_fine_chan,apparent_angular_sky_im_name_fine_chan,angular_sky_vis_filename)
             print(cmd)
             os.system(cmd)          
             
             #full gsm
             cmd = "uvmodel vis=%s model=%s options=replace out=%s" % (uvfits_vis_filename_fine_chan,apparent_sky_im_name_fine_chan,gsm_sky_vis_filename)
             print(cmd)
             os.system(cmd)             
             
             
             #cmd = 'uvplt device="%s/png" vis=%s  axis=uvdist,amp options=nobase select=-auto' % (uv_dist_plot_name,out_vis)
             #print(cmd)
             #os.system(cmd)  
             
             cmd = "fits in=%s op=uvout options=nocal,nopol,nopass out=%s" % (unity_sky_vis_filename,unity_sky_uvfits_filename)
             print(cmd)
             os.system(cmd)
             
             cmd = "fits in=%s op=uvout options=nocal,nopol,nopass out=%s" % (angular_sky_vis_filename,angular_sky_uvfits_filename)
             print(cmd)
             os.system(cmd)

             cmd = "fits in=%s op=uvout options=nocal,nopol,nopass out=%s" % (gsm_sky_vis_filename,gsm_sky_uvfits_filename)
             print(cmd)
             os.system(cmd)
                    
             #get rid of the copied fine chan uvfits and the unity vis
             cmd = "rm -rf %s %s %s %s" % (uvfits_filename_fine_chan,unity_sky_vis_filename,angular_sky_vis_filename,gsm_sky_vis_filename)
             print(cmd)
             os.system(cmd)
             
             ####
       #second hack to only do monopole signal   
       #continue
          
       ###########################
       ###########################
       #print(max_gsm_list)
       #print(scale_fine_chan_list)
        
       #sys.exit()
         
       # predict a multi-channel model
       cmd = "wsclean -predict -name %s -size %s %s -scale %s -pol xx,yy -channels-out 32 %s " % (apparent_sky_fits_name_prefix_fine_chan,wsclean_imsize,wsclean_imsize,wsclean_scale,av_ms_name)
       print(cmd)
       os.system(cmd)
       
       
       ##make  images to check 
       #cmd = "wsclean -name model_col_chan_%s_%s_ms -size %s %s -scale %s -pol xx -data-column MODEL_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
       #print(cmd)
       #os.system(cmd)
       
       #crikey I think it actually works
       
       #hmmm seemed to actually work! We'll see ...
       if uv_cutoff==0:
          #gain_solutions_name = 'cal_%s_%s_calibrate_sols.bin' % (EDA2_chan,EDA2_obs_time)
          gain_solutions_name = "av_chan_%s_%s_plus_%s_obs_cal_sols.bin" % (EDA2_chan,first_obstime,len(obs_time_list))
          calibrate_options = ''
       else:
          #gain_solutions_name = 'cal_%s_%s_calibrate_sols_uvcutoff_%0.3f_m.bin' % (EDA2_chan,EDA2_obs_time,uv_cutoff_m)
          gain_solutions_name = "av_chan_%s_%s_plus_%s_obs_cal_sols_uvcutoff_%0.3f_m.bin" % (EDA2_chan,first_obstime,len(obs_time_list),uv_cutoff_m)
          calibrate_options = '-minuv %0.3f ' % uv_cutoff_m
       
       cmd = "rm -rf %s" % (gain_solutions_name)
       print(cmd)
       os.system(cmd)
       
       #calibrate
       
       #replace this with cal using autos / coherence?
       #########
       #try without on av data set - didn't seem to work
       #calibrate on all 32 chans to increase SNR (poor results if you don't do this)
       cmd = "calibrate  -ch 32 -datacolumn DATA %s %s %s " % (calibrate_options,av_ms_name,gain_solutions_name)
       print(cmd)
       os.system(cmd)
       

       #plot the sols and 
       if (os.path.isfile(gain_solutions_name)):
          if plot_cal:
             #Plot the cal solutions
             cmd = "aocal_plot.py  %s " % (gain_solutions_name)
             print(cmd)
             os.system(cmd)
             
          #indent this if doing precal steps above   

          cmd = "applysolutions %s %s " % (av_ms_name,gain_solutions_name)
          print(cmd)
          os.system(cmd)
       
          ###make an image to check (both pols) 32 chans, skip for now, takes ages 
          #cmd = "wsclean -name cal_chan_%s_%s_ms -size %s %s -scale %s -pol xx,yy -data-column CORRECTED_DATA -channels-out 32 %s " % (EDA2_chan,EDA2_obs_time,wsclean_imsize,wsclean_imsize,wsclean_scale,ms_name)
          #print(cmd)
          #os.system(cmd) 
          
    
          
          #############
   
          #calibrate per chan on time av ms
          if per_chan_cal==True:
             if uv_cutoff==0:
                gain_solutions_name = 'cal_%s_%s_calibrate_sols_per_chan.bin' % (EDA2_chan,first_obstime)
                calibrate_options = ''
             else:
                gain_solutions_name = 'cal_%s_%s_calibrate_sols_uvcutoff_%0.3f_m_per_chan.bin' % (EDA2_chan,first_obstime,uv_cutoff_m)
                calibrate_options = '-minuv %0.3f ' % uv_cutoff_m
   
             #calibrate on each chan and use the corrected data
             cmd = "calibrate  -ch 1 -datacolumn CORRECTED_DATA %s %s %s " % (calibrate_options,av_ms_name,gain_solutions_name)
             #for av ms name use data column
             #cmd = "calibrate  -ch 1 %s %s %s " % (calibrate_options,av_ms_name,gain_solutions_name)
             print(cmd)
             os.system(cmd)
             
             #apply sols
             cmd = "applysolutions -datacolumn CORRECTED_DATA %s %s " % (av_ms_name,gain_solutions_name)
             #cmd = "applysolutions  %s %s " % (av_ms_name,gain_solutions_name)
             print(cmd)
             os.system(cmd)
             
             if plot_cal:
                #Plot the cal solutions
                cmd = "aocal_plot.py  %s " % (gain_solutions_name)
                print(cmd)
                os.system(cmd)
             
                      
          #write out the uvfits file
          casa_cmd_filename = 'export_individual_uvfits.sh'
          cmd = "rm -rf %s %s" % (calibrated_uvfits_filename_wsclean,casa_cmd_filename)
          print(cmd)
          os.system(cmd)
                 
          cmd = "exportuvfits(vis='%s',fitsfile='%s',datacolumn='corrected',overwrite=True,writestation=False)" % (av_ms_name,calibrated_uvfits_filename_wsclean)
          print(cmd)
          os.system(cmd)
          
          with open(casa_cmd_filename,'w') as f:
             f.write(cmd)
               
          cmd = "casa --nohead --nogui --nocrashreport -c %s" % casa_cmd_filename
          print(cmd)
          os.system(cmd)
          
          #end indent
       
          #remove those nasty tmp files that CASA writes
          cmd = "rm -rf /tmp/tmp*"
          print(cmd)
          os.system(cmd)
                    
       else:
          print("no cal solutions for %s" % (av_ms_name))
          continue
   

       
       
       
       
       
   
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
       
       #guard against cases where there are no data for that channel
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
                   
                  cmd = "rm -rf %s %s %s %s %s" % (unity_sky_uvfits_filename,unity_sky_vis_filename,zero_sky_vis_filename,uvfits_vis_filename)
                  print(cmd)
                  os.system(cmd)
                  
                  cmd = "fits in=%s op=uvin out=%s" % (uvfits_filename,uvfits_vis_filename)
                  print(cmd)
                  os.system(cmd)
               
                  cmd = "uvmodel vis=%s model=%s options=replace,mfs out=%s" % (uvfits_vis_filename,apparent_unity_sky_im_name_copy,unity_sky_vis_filename)
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
                           #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                        else:
                           beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                           #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
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
                     
                     if not reverse_fine_chans:
                        #dont reverse chan order
                        freq_MHz_fine_chan = centre_freq + (fine_chan_index - centre_chan_index)*fine_chan_width_MHz
                     else:
                        freq_MHz_fine_chan = centre_freq - (fine_chan_index - centre_chan_index + 1)*fine_chan_width_MHz
                     
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
                              #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'xx')
                           else:
                              beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'yy')
                              #beam_image_sin_projected_fitsname_no_cos_za = "model_%0.3f_MHz_%s_no_cos_za.fits" % (freq_MHz,'yy')
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
   print("txt file needs to be miriad format (North, East, Up)")
   antenna_layout_basename = array_layout_filename.split('/')[-1].split('.')[0]
   antenna_name_list = range(1,257)

   antenna_position_x_list=[]
   antenna_position_y_list=[]
   with open(array_layout_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[1])
      antenna_position_y = float(line.strip().split()[0])
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
   plt.xlabel("Frequency (MHz)")
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
      plt.xlabel("Frequency (MHz)")
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
   plt.xlabel("Frequency (MHz)")
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
   plt.xlabel("Frequency (MHz)")
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

def compare_uvfits(uvfitsname1,uvfitsname2):
   print(uvfitsname1)
   print(uvfitsname2)
   uvfitsname1_base = uvfitsname1.split('.uvfits')[0]
   uvfitsname2_base = uvfitsname2.split('.uvfits')[0]
   hdulist1 = fits.open(uvfitsname1)
   #hdulist1.info()
   #info_string1 = [(x,x.data.shape,x.data.dtype.names) for x in hdulist1]
   #print info_string1
   uvtable1 = hdulist1[0].data
   uvtable_header1 = hdulist1[0].header
   visibilities1 = uvtable1['DATA'][:,0,0,0,0,0]
   n_vis1 = visibilities1.shape[0]

   UU_s_array1 = uvtable1['UU']
   UU_m_array1 = UU_s_array1 * c   
   VV_s_array1 = uvtable1['VV']
   VV_m_array1 = VV_s_array1 * c
   print(n_vis1)
   #print(len(VV_m_array1))

   UU_m_array1_one_timestep = UU_m_array1[0:32385]
   VV_m_array1_one_timestep = VV_m_array1[0:32385]
   visibilities1_one_timestep = visibilities1[0:32385]
   
   hdulist2 = fits.open(uvfitsname2)
   #hdulist2.info()
   #info_string2 = [(x,x.data.shape,x.data.dtype.names) for x in hdulist2]
   #print info_string2
   uvtable2 = hdulist2[0].data
   uvtable_header2 = hdulist2[0].header
   visibilities2 = uvtable2['DATA'][:,0,0,0,0,0]
   n_vis2 = visibilities2.shape[0]

   UU_s_array2 = uvtable2['UU']
   UU_m_array2 = UU_s_array2 * c   
   VV_s_array2 = uvtable2['VV']
   VV_m_array2 = VV_s_array2 * c
   #print(n_vis2)
   #print(len(VV_m_array2))
   
   baseline_length_array_m1 = np.sqrt(UU_m_array1_one_timestep**2 + VV_m_array1_one_timestep**2)
   baseline_length_array_m_inds1 = baseline_length_array_m1.argsort()
   baseline_length_array_m_sorted_orig1 = baseline_length_array_m1[baseline_length_array_m_inds1]
   UU_m_array_sorted_orig1 = UU_m_array1_one_timestep[baseline_length_array_m_inds1]
   VV_m_array_sorted_orig1 = VV_m_array1_one_timestep[baseline_length_array_m_inds1]
   visibilities1_one_timestep_sorted = visibilities1_one_timestep[baseline_length_array_m_inds1]
   
   baseline_length_array_m2 = np.sqrt(UU_m_array2**2 + VV_m_array2**2)
   baseline_length_array_m_inds2 = baseline_length_array_m2.argsort()
   baseline_length_array_m_sorted_orig2 = baseline_length_array_m2[baseline_length_array_m_inds2]
   UU_m_array_sorted_orig2 = UU_m_array2[baseline_length_array_m_inds2]
   VV_m_array_sorted_orig2 = VV_m_array2[baseline_length_array_m_inds2] 
   visibilities2_sorted = visibilities2[baseline_length_array_m_inds2]
   #print(baseline_length_array_m_sorted_orig1[0:10])
   #print(baseline_length_array_m_sorted_orig2[0:10])
   
   #print(VV_m_array_sorted_orig1[0:10])
   #print(UU_m_array_sorted_orig2[0:10])  
   
   #make a uv plot:
   plt.clf()
   plt.scatter(UU_m_array_sorted_orig1,VV_m_array_sorted_orig1,s=1,marker='.')
   map_title="UV footprint" 
   plt.xlabel("UU (m)")
   plt.ylabel("VV (m)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "UV_plot_%s.png" % (uvfitsname1_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)    

   plt.clf()
   plt.scatter(UU_m_array_sorted_orig2,VV_m_array_sorted_orig2,s=1,marker='.')
   map_title="UV footprint" 
   plt.xlabel("UU (m)")
   plt.ylabel("VV (m)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "UV_plot_%s.png" % (uvfitsname2_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)   
   
   #now plot the visibility values as a function of baseline length
   #real
   plt.clf()
   plt.plot(baseline_length_array_m_sorted_orig1,visibilities1_one_timestep_sorted.real)
   map_title="Real vis vs baseline length" 
   plt.xlabel("Baseline length (m)")
   plt.ylabel("Visibility amplitude real (Jy)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "real_vis_vs_baseline_length_%s.png" % (uvfitsname1_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)    
   
   #imag
   plt.clf()
   plt.plot(baseline_length_array_m_sorted_orig1,visibilities1_one_timestep_sorted.imag)
   map_title="Imag vis vs baseline length" 
   plt.xlabel("Baseline length (m)")
   plt.ylabel("Visibility amplitude real (Jy)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "imag_vis_vs_baseline_length_%s.png" % (uvfitsname1_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)

   #real
   plt.clf()
   plt.plot(baseline_length_array_m_sorted_orig2,visibilities2_sorted.real)
   map_title="Real vis vs baseline length" 
   plt.xlabel("Baseline length (m)")
   plt.ylabel("Visibility amplitude real (Jy)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "real_vis_vs_baseline_length_%s.png" % (uvfitsname2_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)    
   
   #imag
   plt.clf()
   plt.plot(baseline_length_array_m_sorted_orig2,visibilities2_sorted.imag)
   map_title="Imag vis vs baseline length" 
   plt.xlabel("Baseline length (m)")
   plt.ylabel("Visibility amplitude real (Jy)")
   #plt.legend(loc=1)
   #plt.ylim([0, 20])
   fig_name= "imag_vis_vs_baseline_length_%s.png" % (uvfitsname2_base)
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print("saved %s" % fig_name)
   
   
def plot_internal_noise_coupling(frequency_MHz_array,mnm_even_filename,antenna_positions_filename,plot=False):
   n_freqs = len(frequency_MHz_array) 
   mnm_even_array = np.load(mnm_even_filename)
   #mnm_even_array_real = mnm_even_array.real
   
   ##deleting Daniels antenna 237 (zero indexed) as it has no match in WODEN sims
   ##need to rerun woden sims with different 255 antennas (i.e. deleting ant x;7.000, y:3.00 (AAVS1_loc_uvgen_match_daniel_255_NEU.ant))
   ##print(mnm_even_array.shape)
   #new_mnm_even_filename = mnm_even_filename.split(".npy")[0].split("/")[-1] + "_255.npy"
   #new_mnm_even_array = np.copy(mnm_even_array)
   #new_mnm_even_array = np.delete(new_mnm_even_array, 237, 1)
   #new_mnm_even_array = np.delete(new_mnm_even_array, 237, 2)
   #print(new_mnm_even_array.shape)
   #np.save(new_mnm_even_filename,new_mnm_even_array)
   #print("saved %s" % new_mnm_even_filename)
   
   #for one freq, plot the real and abs values of the internal noise as a function of baseline length
   #For zenith, u,v is just the diff of the E and the diff of the N
   antenna_position_x_list = []
   antenna_position_y_list = []
   with open(antenna_positions_filename) as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[1])
      antenna_position_y = float(line.strip().split()[0])
      antenna_position_x_list.append(antenna_position_x)
      antenna_position_y_list.append(antenna_position_y)   
   
   antenna_position_x_m = np.asarray(antenna_position_x_list)
   antenna_position_y_m = np.asarray(antenna_position_y_list)
   
   #sort by descending E coord just to compare uv coverage with previous sims
   antenna_position_x_m_ascending_inds = antenna_position_x_m.argsort()
   antenna_position_x_m_sorted = antenna_position_x_m[antenna_position_x_m_ascending_inds][::-1]
   antenna_position_y_m_sorted = antenna_position_y_m[antenna_position_x_m_ascending_inds][::-1]
   
   n_ant = antenna_position_x_m_sorted.shape[0]
   n_baselines = int((n_ant * (n_ant-1.))/2.)
   print("n_ant %s" % n_ant)
   
   #calculate u and v for each baseline
   #uu_list_sorted = []
   #vv_list_sorted = []
   #
   #for ant_1_index in range(0,n_ant):
   #   for ant_2_index in range(0,n_ant):
   #      if (ant_2_index >= ant_1_index):
   #         uu = antenna_position_x_m_sorted[ant_1_index] - antenna_position_x_m_sorted[ant_2_index]
   #         vv = antenna_position_y_m_sorted[ant_1_index] - antenna_position_y_m_sorted[ant_2_index]
   #         uu_list_sorted.append(uu)
   #         vv_list_sorted.append(vv)

   #uu_array_sorted = np.asarray(uu_list_sorted)
   #vv_array_sorted = np.asarray(vv_list_sorted)
   
   #plt.clf()
   #plot_filename = "uv_plot_sorted_test_eda2_daniel.png"
   #plt.scatter(uu_array_sorted,vv_array_sorted,s=1,marker='.')
   #plt.gcf()
   #plt.savefig(plot_filename)
   #print("save %s" % plot_filename)
   
   #okay looks exactly as expected!
   uv_correlation_array = np.zeros((n_baselines,n_freqs+2),dtype=complex)
   for freq_index in range(0,n_freqs):
      print(freq_index)
      #now do the same thing without the sorting so that you can allocate the correct correlation to each baseline  
      #calculate u and v for each baseline
      uu_list = []
      vv_list = []
      correlation_list = []
      
      for ant_1_index in range(0,n_ant):
         for ant_2_index in range(0,n_ant):
            #only uniques baselines and dont include autos
            if (ant_2_index > ant_1_index):
               uu = antenna_position_x_m[ant_1_index] - antenna_position_x_m[ant_2_index]
               vv = antenna_position_y_m[ant_1_index] - antenna_position_y_m[ant_2_index]
               uu_list.append(uu)
               vv_list.append(vv)
               correlation = mnm_even_array[freq_index,ant_1_index,ant_2_index]
               correlation_list.append(correlation)
               
   
      uu_array = np.asarray(uu_list)
      vv_array = np.asarray(vv_list)
      correlation_array = np.asarray(correlation_list) 
      print(correlation_array)
      #also make a file with the u,v and the correlation
      
      #populate uv_correlation_array
      n_baselines = len(uu_array)
      uv_correlation_array[:,0] = uu_array
      uv_correlation_array[:,1] = vv_array
      uv_correlation_array[:,2+freq_index] = correlation_array
   
   uv_correlation_array_filename = "uv_correlation.npy"  
   np.save(uv_correlation_array_filename,uv_correlation_array)
   print("saved %s" % uv_correlation_array_filename)
      
   if plot:
      plt.clf()
      plot_filename = "uv_plot_unsorted_eda2_daniel.png"
      plt.scatter(uu_array,vv_array,s=1,marker='.')
      plt.gcf()
      plt.savefig(plot_filename)
      print("save %s" % plot_filename)
      
      n_baselines = uu_array.shape[0]
      print("n_baselines %s" % n_baselines)
      
      baseline_length_array_m = np.sqrt(uu_array**2 + vv_array**2)
      baseline_length_array_m_inds = baseline_length_array_m.argsort()
      baseline_length_array_m_sorted = baseline_length_array_m[baseline_length_array_m_inds]
      correlation_array_sorted = correlation_array[baseline_length_array_m_inds]
      
      length_list = [35,10,5,2]
      #real
      for length in length_list:
         plt.clf()
         plot_filename = "correlation_real_vs_baseline_length_eda2_daniel_cutoff_%s.png" % length
         plt.plot(baseline_length_array_m_sorted[baseline_length_array_m_sorted<length],correlation_array_sorted[baseline_length_array_m_sorted<length].real)
         plt.gcf()
         plt.savefig(plot_filename)
         print("save %s" % plot_filename)
      
      #abs
      for length in length_list:
         plt.clf()
         plot_filename = "correlation_abs_vs_baseline_length_eda2_daniel_cutoff_%s.png" % length
         plt.plot(baseline_length_array_m_sorted[baseline_length_array_m_sorted<length],abs(correlation_array_sorted[baseline_length_array_m_sorted<length]))
         plt.gcf()
         plt.savefig(plot_filename)
         print("save %s" % plot_filename)   
          
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
        

   
          
def write_woden_sims_sbatch_file(freq_MHz_list,time_string,pol_list):
   #array_layout = "/astro/mwaeor/bmckinley/code/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_match_daniel_255.ant"
   array_layout = "/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_255.txt"
   n_bands = len(freq_MHz_list)
   lowest_channel_freq = freq_MHz_list[0]
   year,month,day,hour,min,sec = time_string.split('_')
   time_formatted = '%d-%02d-%02dT%02d:%02d:%02d' % (float(year),float(month),float(day),float(hour),float(min),float(sec))
   #2015-11-29T15:33:43
   #Jack changed it so you don't need a metafits and you can do 1 MHz-space coarse chans
   type_list = ["gsm","gsm_uniform","EDGES_uniform","unity_uniform","angular"]
   #type_list = ["gsm"]
   #LST_deg = 58.13223745343605    #from template metafits RA=LST
   LST_deg = 60.0
   for type in type_list:
         name_base = "woden_eda2_sbatch_%s_lst_%0.3f" % (type,LST_deg)  
         sbatch_filename = "%s.sh" % name_base
         if (type=="gsm" or type=="EDGES_uniform" or type=="unity_uniform"):
            sourcelist_name = "woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_%s_sourcelist.txt" % (lowest_channel_freq,type)
            output_uvfits_prepend = "/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/data/woden_LST_%0.3f_%s_start_freq_%0.3f" % (LST_deg,type,lowest_channel_freq)
            with open('%s' % sbatch_filename,'w') as outfile:
               outfile.write("#!/bin/bash --login\n#SBATCH --nodes=1\n#SBATCH --partition=gpuq\n#SBATCH --gres=gpu:1\n")
               outfile.write("#SBATCH --time=00:20:00\n#SBATCH --account=mwaeor\n#SBATCH --nodes=1\n#SBATCH --mem=10gb\n")
               outfile.write("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --array=0-%0.0f\n\n" % n_bands)
      
               outfile.write("module swap gcc gcc/5.5.0\nmodule use /pawsey/mwa/software/python3/modulefiles\nmodule load erfa/1.7.0\n")
               outfile.write("module load json-c/0.14\nmodule load hdf5/1.10.5\nmodule load cfitsio/3.48\nmodule load cmake/3.15.0\n")
               outfile.write("module load cuda/10.2\nmodule load pal/0.9.8\nmodule load python/3.8.2\nmodule load astropy/4.0.1.post1\n\n")
      
               outfile.write("source /astro/mwaeor/jline/software/WODEN/build/init_WODEN.sh\n")
               outfile.write("export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH\n\n")
               
               outfile.write("cd /astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN\n")
               outfile.write("mkdir -p data\n")
               outfile.write("time python /astro/mwaeor/jline/software/WODEN/build/run_woden.py \\\n")
               outfile.write("   --ra0=%0.5f --dec0=-26.70 \\\n" % LST_deg)
               outfile.write("   --num_freq_channels=1 --num_time_steps=1 \\\n")
               outfile.write("   --freq_res=10e+3 --time_res=0.28 \\\n")
               outfile.write("   --lowest_channel_freq=%0.0fe+6 \\\n" % lowest_channel_freq)
               outfile.write("   --coarse_band_width=1e+6 \\\n")
               outfile.write("   --cat_filename=/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/%s \\\n" % sourcelist_name)
               #outfile.write("   --metafits_filename=/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/WODEN/centre_chan_%03d_metafits_ppds.fits \\\n" % centre_chan)
               outfile.write("   --output_uvfits_prepend=%s \\\n" % output_uvfits_prepend)
               outfile.write("   --sky_crop_components \\\n")
               #outfile.write("   --EDA2_sim \\\n")
               outfile.write("   --primary_beam=EDA2 \\\n")
               outfile.write("   --array_layout=%s \\\n" % array_layout)
               outfile.write("   --band_nums=${SLURM_ARRAY_TASK_ID} \\\n")
               outfile.write("   --date=%s \\\n" % time_formatted)
               outfile.write("   --chunking_size=2500\n")
               
            print("wrote %s" % sbatch_filename)
         
         for pol in pol_list:
            sbatch_filename = "%s_pol_%s.sh" % (name_base,pol)
            if type =="angular":
               sourcelist_name = "woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_gsm_%s_pol_%s_angular_sourcelist.txt" % (lowest_channel_freq,time_string,pol)
               output_uvfits_prepend = "/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/data/woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular" % (LST_deg,lowest_channel_freq,pol)
            elif type =="gsm_uniform":
               sourcelist_name = "woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_gsm_uniform_%s_pol_%s_sourcelist.txt" % (lowest_channel_freq,time_string,pol)
               output_uvfits_prepend = "/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/data/woden_LST_%0.3f_gsm_uniform_start_freq_%0.3f_pol_%s" % (LST_deg,lowest_channel_freq,pol)            
            else:
               continue
            with open('%s' % sbatch_filename,'w') as outfile:
               outfile.write("#!/bin/bash --login\n#SBATCH --nodes=1\n#SBATCH --partition=gpuq\n#SBATCH --gres=gpu:1\n")
               outfile.write("#SBATCH --time=00:20:00\n#SBATCH --account=mwaeor\n#SBATCH --nodes=1\n#SBATCH --mem=10gb\n")
               outfile.write("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --array=0-%0.0f\n\n" % n_bands)
         
               outfile.write("module swap gcc gcc/5.5.0\nmodule use /pawsey/mwa/software/python3/modulefiles\nmodule load erfa/1.7.0\n")
               outfile.write("module load json-c/0.14\nmodule load hdf5/1.10.5\nmodule load cfitsio/3.48\nmodule load cmake/3.15.0\n")
               outfile.write("module load cuda/10.2\nmodule load pal/0.9.8\nmodule load python/3.8.2\nmodule load astropy/4.0.1.post1\n\n")
         
               outfile.write("source /astro/mwaeor/jline/software/WODEN/build/init_WODEN.sh\n")
               outfile.write("export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH\n\n")
               
               outfile.write("cd /astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN\n")
               outfile.write("mkdir -p data\n")
               outfile.write("time python /astro/mwaeor/jline/software/WODEN/build/run_woden.py \\\n")
               outfile.write("   --ra0=%0.5f --dec0=-26.70 \\\n" % LST_deg)
               outfile.write("   --num_freq_channels=1 --num_time_steps=1 \\\n")
               outfile.write("   --freq_res=10e+3 --time_res=0.28 \\\n")
               outfile.write("   --lowest_channel_freq=%0.0fe+6 \\\n" % lowest_channel_freq)
               outfile.write("   --coarse_band_width=1e+6 \\\n")
               outfile.write("   --cat_filename=/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/%s \\\n" % sourcelist_name)
               #outfile.write("   --metafits_filename=/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/WODEN/centre_chan_%03d_metafits_ppds.fits \\\n" % centre_chan)
               outfile.write("   --output_uvfits_prepend=%s \\\n" % output_uvfits_prepend)
               outfile.write("   --sky_crop_components \\\n")
               #outfile.write("   --EDA2_sim \\\n")
               outfile.write("   --primary_beam=EDA2 \\\n")
               outfile.write("   --array_layout=%s \\\n" % array_layout)
               outfile.write("   --band_nums=${SLURM_ARRAY_TASK_ID} \\\n")
               outfile.write("   --date=%s \\\n" % time_formatted)
               outfile.write("   --chunking_size=2500\n")
                  
            print("wrote %s" % sbatch_filename)       
         
def find_missing_antennas(antenna_positions_filename_1,antenna_positions_filename_2,tolerance_m):
   mapped_antenna_positions_filename = "ant_list_1_mapped_to_ant_list_2.txt"
   antenna_position_x_list_1 = []
   antenna_position_y_list_1 = []
   antenna_position_x_list_2 = []
   antenna_position_y_list_2 = []
   map_ant_1_to_ant_2_x_list = []
   map_ant_1_to_ant_2_y_list = []
   map_ant_1_to_ant_2_map_index_list = []
   with open(antenna_positions_filename_1) as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[1])
      antenna_position_y = float(line.strip().split()[0])
      antenna_position_x_list_1.append(antenna_position_x)
      antenna_position_y_list_1.append(antenna_position_y)   
   with open(antenna_positions_filename_2) as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[1])
      antenna_position_y = float(line.strip().split()[0])
      antenna_position_x_list_2.append(antenna_position_x)
      antenna_position_y_list_2.append(antenna_position_y)  
      
   antenna_position_x_m_1 = np.asarray(antenna_position_x_list_1)
   antenna_position_y_m_1 = np.asarray(antenna_position_y_list_1)
   antenna_position_x_m_2 = np.asarray(antenna_position_x_list_2)
   antenna_position_y_m_2 = np.asarray(antenna_position_y_list_2)  
   
   #go through each line of first ant file and find the corresponding antenna in second ant file
   for x_1_index,x_1 in enumerate(antenna_position_x_m_1):
      match = False
      y_1 = antenna_position_y_m_1[x_1_index]
      for x_2_index,x_2 in enumerate(antenna_position_x_m_2):
         y_2 = antenna_position_y_m_2[x_2_index]
         if (np.isclose(x_1,x_2,atol=float(tolerance_m)) and np.isclose(y_1,y_2,atol=float(tolerance_m))):
            match = True
            map_ant_1_to_ant_2_x_list.append(x_1) 
            map_ant_1_to_ant_2_y_list.append(y_1)    
            map_ant_1_to_ant_2_map_index_list.append(x_2_index)
            #print("antenna_list_1 ant number %s, x:%0.3f, y:%0.3f matches antenna_list_2 ant number %s, x:%0.3f, y:%0.3f" % (x_1_index,x_1,y_1,x_2_index,x_2,y_2))
      if match is False:
         print("NO MATCH FOUND for antenna_list_1 ant number %s (zero indexed), x:%0.3f, y:%0.3f" % (x_1_index,x_1,y_1))

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

def cal_standard_baseline_number(ant1,ant2):
   baseline_number = (ant1 * 256) + ant2
   return baseline_number
           
def add_noise_coupling_to_sim_uvfits(uvfits_filename,uv_correlation_array_filename_x,uv_correlation_array_filename_y):
   #testing for 50 MHz:
   freq_index = 0
   freq_MHz = freq_index * 1.28 + 50.
   wavelength = 300. / freq_MHz 
   jy_to_K = (wavelength**2) / (2. * k * 1.0e26)
   
   # get the values from the noise coupling array created using plot_internal_noise_coupling
   uv_correlation_array_x = np.load(uv_correlation_array_filename_x)
   uv_correlation_array_y = np.load(uv_correlation_array_filename_y)

   uvfits_name_base = uvfits_filename.split(".uvfits")[0].split("/")[-1]
   new_uvfits_name = uvfits_name_base + "_nc.uvfits"
   
   #from Jack:
   #with fits.open('woden_LST_60.000_EDGES_uniform_chan_063_band01.uvfits') as hdulist:
   #    data = hdulist[0].data.data
   #    print(data.shape)
   #    noise = np.random.normal(0,1,data.shape[0])
   #
   #    print("real XX before",data[0,0,0,0,0,0])
   #    print("real YY before",data[0,0,0,0,0,1])
   # 
   #    data[:,0,0,0,0,0] += noise
   #    data[:,0,0,0,0,1] += noise
   # 
   #    hdulist.writeto("test.uvfits",overwrite=True)
   #
   # 
   # with fits.open('test.uvfits') as hdulist:
   #     data = hdulist[0].data.data
   # 
   #     print("real XX after",data[0,0,0,0,0,0])
   #     print("real YY after",data[0,0,0,0,0,1])
   
   print("adding noise coupling to %s" % uvfits_filename)
   with fits.open(uvfits_filename) as hdulist:
      #hdulist.info()
      #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      #print info_string
      uvtable = hdulist[0].data
      visibilities_single = uvtable['DATA']
      visibilities_shape = visibilities_single.shape
      print("visibilities_shape")
      print(visibilities_shape)
   
      UU_s_array = uvtable['UU']
      UU_m_array = UU_s_array * c   
      VV_s_array = uvtable['VV']
      VV_m_array = VV_s_array * c
   
      #baselines = uvtable['BASELINE']
      #for baseline_index in range(250,260):
      #   ant1,ant2 = decode_baseline(baselines[baseline_index])
      #   print(ant1)
      #   print(ant2)
      #sys.exit()
   
      #uv_ratio = UU_m_array / uv_correlation_array_x[:,0]
      #do a better check than this using - isclose()
      #plot uv diff and ratio to check - ordering is right!
      #uv_diff = UU_m_array - uv_correlation_array[:,0]
      
      #plt.clf()
      #plt.plot(uv_ratio)
      #plt.ylabel("ratio woden uv to Daniel uv")
      #plt.xlabel("index")
      #plt.ylim((-2,2))
      #fig_name="uv_ratio_check.png"
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)
   
      #plt.clf()
      #plt.plot(uv_diff)
      #plt.ylabel("diff woden uv to Daniel uv (m)")
      #plt.xlabel("index")
      #plt.ylim((-1,1))
      #fig_name="uv_diff_check.png"
      #figmap = plt.gcf()
      #figmap.savefig(fig_name)
      #print("saved %s" % fig_name)
   
      data = hdulist[0].data.data
      
      ####X pol
      pol_index = 0
 
      internal_noise_real = uv_correlation_array_x[:,2+freq_index].real
      internal_noise_imag = uv_correlation_array_x[:,2+freq_index].imag
      
      internal_noise_real_jy = internal_noise_real / jy_to_K
      internal_noise_imag_jy = internal_noise_imag / jy_to_K
      #always WODEN not wsclean
      data[:,0,0,0,pol_index,0] += internal_noise_real_jy
      data[:,0,0,0,pol_index,1] += internal_noise_imag_jy
       
      ####Y pol
      pol_index = 1
      internal_noise_real = uv_correlation_array_y[:,2+freq_index].real
      internal_noise_imag = uv_correlation_array_y[:,2+freq_index].imag
      
      internal_noise_real_jy = internal_noise_real / jy_to_K
      internal_noise_imag_jy = internal_noise_imag / jy_to_K
      
      data[:,0,0,0,pol_index,0] += internal_noise_real_jy
      data[:,0,0,0,pol_index,1] += internal_noise_imag_jy
 
      #now write out new uvfits file:
      hdulist.writeto(new_uvfits_name,overwrite=True)
      print("saved %s" % (new_uvfits_name))

def split_baseline(baseline, shift=None):
    """
    Given a baseline, split it into it consistent stand ID numbers.
    """
    
    if baseline.any() >= 65536:
        ant1 = np.floor_divide((baseline - 65536),2048)
        ant1 = np.mod((baseline - 65536),2048)
        #ant1 = (baseline - 65536) // 2048
        #ant2 = (baseline - 65536) % 2048
    else:
        ant1 = np.floor_divide(baseline,256)
        ant2 = np.mod(baseline,256)
        
    return ant1,ant2

def convert_matlab_EEPs_by_ant(ant_index,freq_MHz_list,nside,method='cubic',analytic_beam_dir='/md0/EoR/ASSASSIN/beam_fits/fine_chan/',antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt'):
   #method = 'nearest', 'linear', 'cubic'
   with open(antenna_layout_filename) as f:
      lines = f.readlines()
   ant_name = lines[ant_index].split()[0]
   npix = hp.nside2npix(nside)
   print("npix is %s" % npix)
   #['kx', 'geo_ph', '__header__', '__globals__', 'Ephi', 'kz', 'ky', '__version__', 'Etheta']
   #(<type 'numpy.ndarray'>, (361, 91, 256))
   #(<type 'numpy.ndarray'>, (361, 91, 256))
   print("converting MATLAB EEPs")
   for pol in ['X','Y']:
      for freq_MHz in freq_MHz_list:
         wavelength = 300./ freq_MHz
         freq_Hz_string = "%d" % (freq_MHz*1000000)
         EEP_name = '/md0/EoR/EDA2/EEPs/new_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (freq_Hz_string,pol)
         beam_data = loadmat(EEP_name)
         print("loaded %s " % EEP_name)
         #print(beam_data.keys())
         #print(type(beam_data['Ephi']),beam_data['Ephi'].shape)
         #print(type(beam_data['Etheta']),beam_data['Etheta'].shape)
         #These are complex arrays
         #print(E_phi.dtype)
         #print(E_phi.shape)
         #print(E_phi[360][1])
         #print(E_phi[0][1])
         #data for az 0 and az 360 is same as expected, discard angle 360
         #(361, 91, 256)
         E_phi = beam_data['Ephi'][0:361][:]
         E_theta = beam_data['Etheta'][0:361][:]
         # The format of the data is in standard spherical coordinates (theta, phi).
         # Where theta is the zenith angle and phi is anti-clockwise going from east to north looking down at the array. 
         # The 'origin' will be on the top left corner. Going from top to bottom is increasing in phi and going left to right is increasing in theta. 
         
         #daniels azimuth different from analytic beam - need to revcerse it and add 90 deg
         #azimuth_deg_array = np.arange(361)
         azimuth_deg_array = np.flip(np.arange(361) + 90.)
         azimuth_deg_array[azimuth_deg_array >= 361] = azimuth_deg_array[azimuth_deg_array >= 361] - 361
         
         zenith_angle_deg_array = np.arange(91)
         
         beam_power_cube = np.abs(E_theta)**2 + np.abs(E_phi)**2
         #print(beam_power_cube.shape)
         
         beam_power_cube_slice = beam_power_cube[:,:,ant_index]
         #print(beam_power_cube_slice.shape)
         #print(beam_power_cube_slice[:,0])
         #sys.exit()
         #beam_power_cube_slice_transpose = np.transpose(beam_power_cube_slice)
         beam_power_cube_slice_flat = beam_power_cube_slice.flatten('F')
         
         #see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html
         
         repetitions = 91
         
         ang_repeats_array = np.tile(azimuth_deg_array, (repetitions, 1))
         #print(az_repeats_array[0:2])
         #print(az_repeats_array.shape)
         az_ang_repeats_array_flat = ang_repeats_array.flatten()
         #print(az_ang_repeats_array_flat.shape)
         #print(az_ang_repeats_array_flat[0:750])
         
         repetitions = 361
         ang_repeats_array = np.tile(zenith_angle_deg_array, (repetitions, 1))
         
         zenith_angle_repeats_array_flat = ang_repeats_array.flatten('F')
         #print(zenith_angle_repeats_array_flat.shape)
         #print(zenith_angle_repeats_array_flat[0:750])
         
         zenith_angle_repeats_array_flat_rad = zenith_angle_repeats_array_flat / 180. * np.pi
         az_ang_repeats_array_flat_rad = (az_ang_repeats_array_flat) / 180. * np.pi
         
         #print(zenith_angle_repeats_array_flat_rad.shape)
         #print(az_ang_repeats_array_flat_rad.shape)
         #print(beam_power_cube_slice_flat.shape)
         
         #beam_interp = interp2d(zenith_angle_repeats_array_flat_rad,az_ang_repeats_array_flat_rad, beam_power_cube_slice_flat)
             
         
         hpx_pix_num_array = np.arange(npix)
         hpx_angles_rad = hp.pix2ang(nside,hpx_pix_num_array) 
         hpx_angles_rad_zenith_angle = hpx_angles_rad[0]
         hpx_angles_rad_azimuth = hpx_angles_rad[1]
         
         #print(az_ang_repeats_array_flat_rad[0:750])
         #print(zenith_angle_repeats_array_flat_rad[0:750])
         #print(beam_power_cube_slice_flat[0:750])

         regridded_to_hpx_beam_power = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), beam_power_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method=method)
         regridded_to_hpx_beam_power_norm = regridded_to_hpx_beam_power / np.nanmax(regridded_to_hpx_beam_power)
         regridded_to_hpx_beam_power_norm_log = np.log10(regridded_to_hpx_beam_power_norm)
         
         #rotate appropriately:
         dec_rotate = 90. - float(mwa_latitude_ephem)
         r_beam = hp.Rotator(rot=[0,dec_rotate], coord=['C', 'C'], deg=True)  
         rotated_beam = r_beam.rotate_map(regridded_to_hpx_beam_power_norm)
         
         #sim_beam_hpx_fits_name = 'sim_beam_rotated_hpx_%s_%0.3f_MHz.fits' % (pol,freq_MHz)
         #hp.write_map(sim_beam_hpx_fits_name,rotated_beam,coord='C',nest=False,overwrite=True)
         
         
         ##print(regridded_to_hpx_beam_power.shape)
         ##view regridded beam
         plt.clf()
         map_title="rotated beam sim"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=rotated_beam,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="rotated_sim_beam_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)

         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s " % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
         
         
         #This is somewhat tricky, the analytic beam images are in sin projection, the sims are in spherical coords, need to reproject
         #1. can generate hpx_analytic beams and then treat them same as the sims
         #2. need to check that when you take an analytic hpx beam and reproject it to sin projection it matches the beam images we are using
         
         #compare to : make_hpx_beam(NSIDE,pol,wavelength,dipole_height_m):
         
         #azimuth angle is defined differently between analytic beams and sime - add 90 deg somewhere
         
         analytic_beam_hpx = make_hpx_beam(nside,pol,wavelength,dipole_height_m)
         
         #also rotate analytic beam: 
         rotated_beam_analytic = r_beam.rotate_map(analytic_beam_hpx)
         
         analytic_beam_hpx_fits_name = 'analytic_beam_rotated_hpx_%s_%0.3f_MHz.fits' % (pol,freq_MHz)
         hp.write_map(analytic_beam_hpx_fits_name,rotated_beam_analytic,coord='C',nest=False,overwrite=True)
         
         
         plt.clf()
         map_title="rotated analytic beam hpx"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=rotated_beam_analytic,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="rotated_analytic_beam_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s" % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
                  
         #beam_image_sin_projected_fitsname = "model_%0.3f_MHz_%s.fits" % (freq_MHz,'xx')
         
         #I think we need to use adaptive resampling: https://reproject.readthedocs.io/en/stable/celestial.html#adaptive-resampling
         #(or not use reproject_from_healpix, just do everything in hpx)
         #So I think what I was doing originally is sort of okay because the interpolation of the gsm was done in surface brightess units (Temp)
         #Still some inaccuracy as the same pixel size in sr was assumed for all pix, which is not the case but maybe close enough for near 
         #zenith where beam is peaked. Acutally, is this roughly taken into account by the sin projection term in the beams - think so.
         
         
         
         #Forget the reproection
         #try using reproject_from_hpx to reproject to your RA 0, Dec=latitude beam image in sin projection
         #model_beam_image_sin_proj_name = '/md0/EoR/ASSASSIN/beam_fits/model_100_MHz_xx.fits'
         #with fits.open(model_beam_image_sin_proj_name) as hdulist:
         #   model_header = hdulist[0].header
         #
         #reproj_imsize = int(model_header['naxis1'])
         #print(reproj_imsize)
         #target_wcs = WCS(model_header)
         #
         #hdu_sim_beam = fits.open(sim_beam_hpx_fits_name)[1]
         #reprojected_sim_beam,footprint = reproject_from_healpix(hdu_sim_beam,target_wcs,shape_out=(reproj_imsize,reproj_imsize), order='bilinear',field=0)
         #write the map to fits
         #reprojected_sim_beam_fitsname = "reprojected_sim_beam_%s_%0.3f_MHz.fits" % (pol,freq_MHz)
         #fits.writeto(reprojected_sim_beam_fitsname,reprojected_sim_beam,overwrite=True)
         #fits.update(reprojected_sim_beam_fitsname,reprojected_sim_beam,header=model_header)
         #print("wrote image %s" %  reprojected_sim_beam_fitsname)
         #
         #repeat for analytic:
         #hdu_analytic_beam = fits.open(analytic_beam_hpx_fits_name)[1]
         #reprojected_analytic_beam,footprint = reproject_from_healpix(hdu_analytic_beam,target_wcs,shape_out=(reproj_imsize,reproj_imsize), order='bilinear',field=0)
         #write the map to fits
         #reprojected_analytic_beam_fitsname = "reprojected_analytic_beam_%s_%0.3f_MHz.fits" % (pol,freq_MHz)
         #fits.writeto(reprojected_analytic_beam_fitsname,reprojected_analytic_beam,overwrite=True)
         #fits.update(reprojected_analytic_beam_fitsname,reprojected_analytic_beam,header=model_header)
         #print("wrote image %s" %  reprojected_analytic_beam_fitsname)
   
   
         #Do our beam comparisons in healpix
         diff_analytic_sim_beams = rotated_beam_analytic - rotated_beam
         plt.clf()
         map_title="diff analytic sim beams"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=diff_analytic_sim_beams,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="diff_analytic_sim_beams_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s" % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
         
      #make movies:
      movie_name = "beam_diff_movie_ant_index_%s_%s.mp4" % (ant_index,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "diff_analytic_sim_beams_ant_index_%s_%s_*_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (ant_index,pol,movie_name)
      print(cmd)
      os.system(cmd)

      movie_name = "analytic_beam_movie_ant_index_%s_%s.mp4" % (ant_index,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "rotated_analytic_beam_ant_index_%s_%s_*_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (ant_index,pol,movie_name)
      print(cmd)
      os.system(cmd)

      movie_name = "sim_beam_movie_ant_index_%s_%s.mp4" % (ant_index,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "rotated_sim_beam_ant_index_%s_%s_*_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (ant_index,pol,movie_name)
      print(cmd)
      os.system(cmd)
         
         ##print(regridded_to_hpx_beam_power.shape)
         ##view regridded beam
         #plt.clf()
         #map_title="regridded beam sim"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         #hp.orthview(map=regridded_to_hpx_beam_power_norm_log,half_sky=True,rot=(0,90,0),title=map_title)
         ##hp.mollview(map=regridded_to_hpx_beam_power,coord='C',title=map_title,rot=(0,0,0)) #,min=0, max=7000)
         #fig_name="regridded_sim_beam_log_%s_%0.3f_MHz.png" % (pol,freq_MHz)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print("saved %s" % fig_name)        
         
         ##cube 256 stuff - probably wrong
         #beam_power_cube_flat = beam_power_cube.reshape(-1, beam_power_cube.shape[-1])
         #print(beam_power_cube_flat.shape)
         #
         #an_array = azimuth_deg_array
         #repetitions = 91
         #repeats_array = np.transpose([an_array] * repetitions)
         #print(repeats_array.shape)
 
         #cube_reps = 256
         #phi_theta_array_cube = np.repeat(repeats_array[:, :, np.newaxis], cube_reps, axis=2)
         #print(phi_theta_array_cube.shape)
 
         #phi_theta_array_cube_flat = phi_theta_array_cube.reshape(-1, phi_theta_array_cube.shape[-1])
         #print(phi_theta_array_cube_flat.shape)
         
         #beam_interp = interp2d(phi_theta_array_cube_flat, azimuth_deg_array, beam_power_cube)
         #print(beam_interp.shape)
         ####
         
         #look at code/AAVS-1/dipole_beam_tests.py for what you did before, antenna ordering etc
         #from daniel email to Jishnu: 
         
         ##Below is from Jishnu (Load_beam_feko_mollview_1.ipynb): 
         # "Etheta_elem1 = beam_data['Etheta_elem1']\n",
         #"Etheta_elem2 = beam_data['Etheta_elem2']\n",
         #"\n",
         #"Ephi_elem1   = beam_data['Ephi_elem1']\n",
         #"Ephi_elem2   = beam_data['Ephi_elem2']\n",
         #"\n",
         #"Az_pt = beam_data['Az_pt'][0,:]\n",
         #"ZA_pt = beam_data['ZA_pt'][0,:]\n",
         #"\n",
         #"fullbeam_1 = np.abs(Etheta_elem1)**2 + np.abs(Ephi_elem1)**2\n",
         #"fullbeam_2 = np.abs(Etheta_elem2)**2 + np.abs(Ephi_elem2)**2\n",
         #"\n",
         #"def calc_beam_values(n_side, n_pix, nu):\n",
         #"    \n",
         #"    beam_index = int(round(nu-70))\n",
         #"    \n",
         #"    if beam_index < 0:\n",
         #"        beam_index = 0\n",
         #"    if beam_index > 150:\n",
         #"        beam_index = -1\n",
         #"    \n",
         #"    beam_interp = interpolate.interp2d(ZA_pt, Az_pt, fullbeam_1[:,:,beam_index])\n",
         #"    dipole_beam_map = np.zeros(n_pix, dtype=float)\n",
         #"    \n",
         #"    for pix_index in range(n_pix):   \n",
         #"        ZA, Az = hp.pix2ang(n_side, pix_index, nest=False, lonlat=False)   \n",
         #"        ZA = np.pi-ZA\n",
         #"        if ZA<np.pi/2:\n",
         #"            dipole_beam_map[pix_index] = beam_interp(ZA*180/np.pi, Az*180/np.pi)\n",
         #"\n",
         #"    beam_max = np.max(dipole_beam_map)\n",
         #"    dipole_beam_map = dipole_beam_map / beam_max \n",
         #"    \n",
         #"    return dipole_beam_map"  
         #     
         #     beam_map = calc_beam_values(n_side=64, n_pix=49152, nu=100)\n",
         #     r_beam = hp.Rotator(rot=[-90,-90], coord=['C', 'C'], deg=True)  #rot=[ra_deg, dec_deg]
         #     antenna_beam = r_beam.rotate_map_pixel(calc_beam_values(n_side=64, n_pix=49152, nu=150))

def convert_matlab_EEPs_by_freq(freq_MHz_list,freq_MHz_list_index,nside,method='cubic',analytic_beam_dir='/md0/EoR/ASSASSIN/beam_fits/fine_chan/',antenna_layout_filename='/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt'):
   #method = 'nearest', 'linear', 'cubic'
   npix = hp.nside2npix(nside)
   
   with open(antenna_layout_filename) as f:
      lines = f.readlines()
      n_ants = len(lines) 

   freq_MHz = freq_MHz_list[freq_MHz_list_index]
   wavelength = 300./ freq_MHz
   freq_Hz_string = "%d" % (freq_MHz*1000000)
         
   print("converting MATLAB EEPs at %0.3f MHz freq for all antennas in %s" % (freq_MHz,antenna_layout_filename))
   for pol in ['X','Y']:
      EEP_name = '/md0/EoR/EDA2/EEPs/new_20210616/FEKO_EDA2_256_elem_%sHz_%spol.mat' % (freq_Hz_string,pol)
      beam_data = loadmat(EEP_name)
      print("loaded %s " % EEP_name)

      E_phi = beam_data['Ephi'][0:361][:]
      E_theta = beam_data['Etheta'][0:361][:]
      # The format of the data is in standard spherical coordinates (theta, phi).
      # Where theta is the zenith angle and phi is anti-clockwise going from east to north looking down at the array. 
      # The 'origin' will be on the top left corner. Going from top to bottom is increasing in phi and going left to right is increasing in theta. 
      
      #daniels azimuth different from analytic beam - need to revcerse it and add 90 deg
      #azimuth_deg_array = np.arange(361)
      azimuth_deg_array = np.flip(np.arange(361) + 90.)
      azimuth_deg_array[azimuth_deg_array >= 361] = azimuth_deg_array[azimuth_deg_array >= 361] - 361
      
      zenith_angle_deg_array = np.arange(91)
      
      beam_power_cube = np.abs(E_theta)**2 + np.abs(E_phi)**2
      #print(beam_power_cube.shape)
      for ant_index in range(n_ants):
         ant_name = lines[ant_index].split()[0]
         beam_power_cube_slice = beam_power_cube[:,:,ant_index]
         #print(beam_power_cube_slice.shape)
         #print(beam_power_cube_slice[:,0])
         #sys.exit()
         #beam_power_cube_slice_transpose = np.transpose(beam_power_cube_slice)
         beam_power_cube_slice_flat = beam_power_cube_slice.flatten('F')
         
         #see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html
         
         repetitions = 91
         
         ang_repeats_array = np.tile(azimuth_deg_array, (repetitions, 1))
         #print(az_repeats_array[0:2])
         #print(az_repeats_array.shape)
         az_ang_repeats_array_flat = ang_repeats_array.flatten()
         #print(az_ang_repeats_array_flat.shape)
         #print(az_ang_repeats_array_flat[0:750])
         
         repetitions = 361
         ang_repeats_array = np.tile(zenith_angle_deg_array, (repetitions, 1))
         
         zenith_angle_repeats_array_flat = ang_repeats_array.flatten('F')
         #print(zenith_angle_repeats_array_flat.shape)
         #print(zenith_angle_repeats_array_flat[0:750])
         
         zenith_angle_repeats_array_flat_rad = zenith_angle_repeats_array_flat / 180. * np.pi
         az_ang_repeats_array_flat_rad = (az_ang_repeats_array_flat) / 180. * np.pi
         
         #print(zenith_angle_repeats_array_flat_rad.shape)
         #print(az_ang_repeats_array_flat_rad.shape)
         #print(beam_power_cube_slice_flat.shape)
         
         #beam_interp = interp2d(zenith_angle_repeats_array_flat_rad,az_ang_repeats_array_flat_rad, beam_power_cube_slice_flat)
             
         hpx_pix_num_array = np.arange(npix)
         hpx_angles_rad = hp.pix2ang(nside,hpx_pix_num_array) 
         hpx_angles_rad_zenith_angle = hpx_angles_rad[0]
         hpx_angles_rad_azimuth = hpx_angles_rad[1]
         
         #print(az_ang_repeats_array_flat_rad[0:750])
         #print(zenith_angle_repeats_array_flat_rad[0:750])
         #print(beam_power_cube_slice_flat[0:750])

         regridded_to_hpx_beam_power = griddata((az_ang_repeats_array_flat_rad,zenith_angle_repeats_array_flat_rad), beam_power_cube_slice_flat, (hpx_angles_rad_azimuth,hpx_angles_rad_zenith_angle), method=method)
         regridded_to_hpx_beam_power_norm = regridded_to_hpx_beam_power / np.nanmax(regridded_to_hpx_beam_power)
         regridded_to_hpx_beam_power_norm_log = np.log10(regridded_to_hpx_beam_power_norm)
         
         #rotate appropriately:
         dec_rotate = 90. - float(mwa_latitude_ephem)
         r_beam = hp.Rotator(rot=[0,dec_rotate], coord=['C', 'C'], deg=True)  
         rotated_beam = r_beam.rotate_map(regridded_to_hpx_beam_power_norm)
              
         ##print(regridded_to_hpx_beam_power.shape)
         ##view regridded beam
         plt.clf()
         map_title="rotated beam sim"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=rotated_beam,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="rotated_sim_beam_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)

         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s " % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
         
         
         #This is somewhat tricky, the analytic beam images are in sin projection, the sims are in spherical coords, need to reproject
         #1. can generate hpx_analytic beams and then treat them same as the sims
         #2. need to check that when you take an analytic hpx beam and reproject it to sin projection it matches the beam images we are using
         
         #compare to : make_hpx_beam(NSIDE,pol,wavelength,dipole_height_m):
         
         #azimuth angle is defined differently between analytic beams and sime - add 90 deg somewhere
         
         analytic_beam_hpx = make_hpx_beam(nside,pol,wavelength,dipole_height_m)
         
         #also rotate analytic beam: 
         rotated_beam_analytic = r_beam.rotate_map(analytic_beam_hpx)
         
         analytic_beam_hpx_fits_name = 'analytic_beam_rotated_hpx_%s_%0.3f_MHz.fits' % (pol,freq_MHz)
         hp.write_map(analytic_beam_hpx_fits_name,rotated_beam_analytic,coord='C',nest=False,overwrite=True)
         
         
         plt.clf()
         map_title="rotated analytic beam hpx"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=rotated_beam_analytic,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="rotated_analytic_beam_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s" % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
         
         #Do our beam comparisons in healpix
         diff_analytic_sim_beams = rotated_beam_analytic - rotated_beam
         plt.clf()
         map_title="diff analytic sim beams"
         ###hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=diff_analytic_sim_beams,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
         fig_name="diff_analytic_sim_beams_ant_index_%s_%s_%0.3f_MHz.png" % (ant_index,pol,freq_MHz)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print("saved %s" % fig_name)
         
         #Add some text to the png
         img = Image.open("%s" % fig_name)
         draw = ImageDraw.Draw(img)
         font = ImageFont.truetype('FreeSans.ttf',30)
         draw.text((10, 10),"%0.3f MHz\n  %s" % (freq_MHz,ant_name),(0,0,0),font=font)
         img.save("%s" % fig_name)
            
         
      #make movies:
      movie_name = "beam_diff_movie_%0.3f_MHz_%s.mp4" % (freq_MHz,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "diff_analytic_sim_beams_ant_index_*_%s_%0.3f_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (pol,freq_MHz,movie_name)
      print(cmd)
      os.system(cmd)

      movie_name = "analytic_beam_movie_%0.3f_MHz_%s.mp4" % (freq_MHz,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "rotated_analytic_beam_ant_index_*_%s_%0.3f_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (pol,freq_MHz,movie_name)
      print(cmd)
      os.system(cmd)

      movie_name = "sim_beam_movie_%0.3f_MHz_%s.mp4" % (freq_MHz,pol)
      cmd = 'ffmpeg -framerate 2 -pattern_type glob -i "rotated_sim_beam_ant_index_*_%s_%0.3f_MHz.png" -c:v libx264 -r 30 -pix_fmt yuv420p -nostdin -y %s' % (pol,freq_MHz,movie_name)
      print(cmd)
      os.system(cmd)


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

#from Jack
def enh2xyz(east,north,height,latitiude):
    sl = sin(latitiude)
    cl = cos(latitiude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl
    return X,Y,Z
    
def get_uvw_freq(x_length=None,y_length=None,z_length=None,dec=None,ha=None,freq=None):
    '''Takes the baseline length in metres and uses the frequency'''
    
    scale = freq / VELC
    X = x_length * scale
    Y = y_length * scale
    Z = z_length * scale
    
    u = sin(ha)*X + cos(ha)*Y
    v = -sin(dec)*cos(ha)*X + sin(dec)*sin(ha)*Y + cos(dec)*Z
    w = cos(dec)*cos(ha)*X - cos(dec)*sin(ha)*Y + sin(dec)*Z
    
    return u,v,w
    
    
               
def get_antenna_table_from_uvfits(uvfits_name):
   print("getting antenna table from %s " % uvfits_name)
   with fits.open(uvfits_name) as hdulist:
      uvtable = hdulist[0].data
      uvtable_header = hdulist[0].header
      hdulist.info()
      info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
      print(info_string)
      #print(uvtable_header)
      #see https://pyuvdata.readthedocs.io/en/v1.5/_modules/pyuvdata/uvfits.html
      ant_hdu = hdulist['AIPS AN']
      print(ant_hdu.data)

def cross_match_eda2_ant_pos_with_sims(ant_pos_on_ground_filename,ant_pos_sim_filename,combined_ant_pos_name_filename):
   ##Open up one of Daniel's files and take a look:
   #Xpol Ephi Phase:
   antenna_name_list = range(1,257)

   #Daniels sims:
   antenna_position_x_list=[]
   antenna_position_y_list=[]
   with open(ant_pos_sim_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      antenna_position_x = float(line.strip().split()[1])
      antenna_position_y = float(line.strip().split()[2])
      antenna_position_x_list.append(antenna_position_x)
      antenna_position_y_list.append(antenna_position_y)   
   
   antenna_position_x_m = np.asarray(antenna_position_x_list)
   antenna_position_y_m = np.asarray(antenna_position_y_list)
   
   
   #Plot antenna positions
   antenna_position_plot_figname = 'ant_pos_eda2_sims.png'
   antenna_position_plot_title = 'Antenna Positions EDA2 Beam Sims '
   
   fig, ax = plt.subplots()
   ax.scatter(antenna_position_x_m,antenna_position_y_m, marker='.')
   
   for i, name in enumerate(antenna_name_list):
       ax.annotate(str(name), (antenna_position_x_m[i],antenna_position_y_m[i]),size=5)
       
   plt.xlabel('X offset from centre (m) ')
   plt.ylabel('Y offset from centre (m) ')
   plt.title(antenna_position_plot_title)
   plt.savefig(antenna_position_plot_figname,dpi = 300)
   print('saved %s ' % antenna_position_plot_figname)
   
   #Randalls:
   antenna_position_x_list_on_ground=[]
   antenna_position_y_list_on_ground=[]
   antenna_names_list_on_ground=[]
   new_antenna_position_list=[]
   with open(ant_pos_on_ground_filename,'r') as f:
      lines = f.readlines()
   for line in lines:
      antenna_name_on_ground = line.strip().split('\t')[0]
      antenna_position_x_on_ground = float(line.strip().split()[1])
      antenna_position_y_on_ground = float(line.strip().split()[2])
      antenna_names_list_on_ground.append(antenna_name_on_ground)
      antenna_position_x_list_on_ground.append(antenna_position_x_on_ground)
      antenna_position_y_list_on_ground.append(antenna_position_y_on_ground)
         
   
   antenna_position_x_m_on_ground = np.asarray(antenna_position_x_list_on_ground)
   antenna_position_y_m_on_ground = np.asarray(antenna_position_y_list_on_ground)

   position_tolerance = 0.1
   for pos_index,x_pos_on_ground in enumerate(antenna_position_x_m_on_ground):
      y_pos_on_ground = antenna_position_y_m_on_ground[pos_index]
      antenna_name_on_ground = antenna_names_list_on_ground[pos_index]
      for x_pos_sim_index,x_pos_sim in enumerate(antenna_position_x_m):
         y_pos_sim = antenna_position_y_m[x_pos_sim_index]
         antenna_name_sim = str(x_pos_sim_index+1)
         x_diff = abs(x_pos_on_ground-x_pos_sim)
         y_diff = abs(y_pos_on_ground-y_pos_sim)
         if ((x_diff < position_tolerance) and (y_diff < position_tolerance)):
            new_position_line = "%s %s %s %s %s %s" % (antenna_name_on_ground,x_pos_on_ground,y_pos_on_ground,antenna_name_sim,x_pos_sim,y_pos_sim)
            new_antenna_position_list.append(new_position_line)
   
   with open (combined_ant_pos_name_filename, 'w') as outfile:
       outfile.write("\n".join(new_antenna_position_list))
   print('wrote %s' % combined_ant_pos_name_filename)
   
   
   #Plot antenna positions
   antenna_position_plot_figname = 'antenna_positions_aavs1_beam_tests_randall.png'
   antenna_position_plot_title = 'Antenna Positions AAVS-1 Beam Tests Randall'
   
   fig, ax = plt.subplots()
   ax.scatter(antenna_position_x_m,antenna_position_y_m, marker='.')
   
   for i, name in enumerate(antenna_name_list):
       ax.annotate(str(name), (antenna_position_x_m[i],antenna_position_y_m[i]),size=5)
   
   plt.xlabel('X offset from centre (m) ')
   plt.ylabel('Y offset from centre (m) ')
   plt.title(antenna_position_plot_title)
   plt.savefig(antenna_position_plot_figname,dpi = 300)
   print('save %s ' % antenna_position_plot_figname)



         
#antenna_positions_filename_1 = "/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_255_NEU.ant"  
#antenna_positions_filename_2 = "/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_NEU.ant"        
#antenna_positions_filename_1 = "/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_match_daniel_255_NEU.ant"
#antenna_positions_filename_1 = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_NEU.txt"
#antenna_positions_filename_2 = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_NEU_255.txt"
#tol_m = 0.1
#find_missing_antennas(antenna_positions_filename_1,antenna_positions_filename_2,tol_m)
#sys.exit() 
#NO MATCH FOUND for antenna_list_1 ant number 237, x:-0.300, y:-0.254

#modified plot_internal_noise_coupling to write a new internal_noise_matrix_filename, without Daniels antenna 237 (zero indexed)
#run this just once to get new mnm_even_255 with daniels ant 237 removed
#internal_noise_matrix_filename = "/md0/EoR/ASSASSIN/noise_coupling/mnm_even_eda2_y.npy"
#antenna_positions_filename = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_NEU.txt"
##Then run this for correct uv_correlation file
#internal_noise_matrix_filename = "/md0/EoR/ASSASSIN/noise_coupling/mnm_even_eda2_255_x.npy"
#antenna_positions_filename = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_NEU_255.txt"
#frequency_MHz_array_mnm = (np.arange(0,218) * 1.28 ) + 50
#plot_internal_noise_coupling(frequency_MHz_array_mnm,internal_noise_matrix_filename,antenna_positions_filename)
#sys.exit()

#year,month,day,hour,min,sec = 2015,11,29,15,40,29
#time_string = '%d%02d%02dT%02d%02d%02d' % (year,month,day,hour,min,sec)
#lst_for_woden_hrs = get_eda2_lst(eda_time_string=time_string)
#lst_for_woden_deg = lst_for_woden_hrs*15.
#diff_with_60_hrs = (60. / 15.) - lst_for_woden_hrs
#print(lst_for_woden_deg)
#print(diff_with_60_hrs)
#diff_with_60_min = diff_with_60_hrs * 60.
#print(diff_with_60_min)
#diff_with_60_sec = diff_with_60_min * 60.
#print(diff_with_60_sec)
#sys.exit()

#uvfits_filename = "/md0/EoR/ASSASSIN/WODEN/data/global_edges/woden_LST_60.000_EDGES_uniform_chan_063_band01.uvfits"  
#uvfits_filename = "/md0/EoR/ASSASSIN/noise_coupling/woden_LST_60.000_EDGES_uniform_start_freq_50.000_band00.uvfits"  
#uv_correlation_array_filename_x = "uv_correlation.npy"
#add_noise_coupling_to_sim_uvfits(uvfits_filename,uv_correlation_array_filename_x,uv_correlation_array_filename_x)
#sys.exit()





#SIMS

#calculate the global 21cm signal:
#s_21_array = plot_S21(nu_array=freq_MHz_list,C=C,A=A,delta_nu=delta_nu,nu_c=nu_c)
#s_21_array_EDGES = plot_S21_EDGES(nu_array=freq_MHz_list)

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
#lst_hrs_list = ['2']

EDA2_data = True

#EDA2_filenames = ["chan_64_20191202T171525_calibrated.uvfits","chan_77_20191202T171629_calibrated.uvfits","chan_90_20191202T171727_calibrated.uvfits","chan_103_20191202T171830_calibrated.uvfits","chan_116_20191202T171928_calibrated.uvfits","chan_129_20191202T172027_calibrated.uvfits"]

#20190929 data:
#EDA2_chan_list = [64,77,90,103,116,129]

#20200217 data (dont use 63 / 49 MHz dont have beam model! missing 127 and 128 so just do to incl 126)
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
EDA2_data_dir = '/md0/EoR/EDA2/20200303_data/'   #2020 paper
#EDA2_data_dir = '/md0/EoR/EDA2/20210813/'

#20200303:
EDA2_chan_list = range(64,127) 
#EDA2_chan_list = [142]

#inv fine chans
#EDA2_data_dir = '/md0/EoR/EDA2/inv_uvfits/20200303_213605/'
#20200217 data (don't use 'edge' times):
#EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan("/md0/EoR/EDA2/20200303_data/",EDA2_chan_list)
#EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan("/md0/EoR/EDA2/20200304_data/",EDA2_chan_list)
EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan(EDA2_data_dir,EDA2_chan_list)


EDA2_obs_time_list_each_chan = EDA2_obs_time_list_each_chan[0:]

n_obs_concat_list = [len(obs_list) for obs_list in EDA2_obs_time_list_each_chan] 

EDA2_obs_time_list = [item[0] for item in EDA2_obs_time_list_each_chan] 

#test just to see when Galaxy rises now
#EDA2_obs_time_list = ['20210621T150000']

#just for SIMS for EDA2 data (REMOVE!!!), have to manually cd to chan dir and run just simulate:
#EDA2_chan_list = [129]
#EDA2_obs_time_list = ['20191202T172027']

EDA2_chan_list_array = np.asarray(EDA2_chan_list)

#EDA2_obs_time = '20191202T171727'
#freq_MHz_array = np.round(400./512.*EDA2_chan_list_array)
freq_MHz_array = 400./512.*EDA2_chan_list_array
freq_MHz_list = freq_MHz_array[0:]
EDA2_chan_list = EDA2_chan_list[0:]

#print(freq_MHz_list)
#print(EDA2_chan_list)

#eda2_data_filename = "chan_%s_%s.uvfits" % (int(EDA2_chan),EDA2_obs_time)
#eda2_data_filename = "chan_90_20191202T171727.uvfits"
#eda2_data_uvfits_name_list=[eda2_data_filename]



lst_hrs_list = []
for EDA2_obs_time_index,EDA2_obs_time in enumerate(EDA2_obs_time_list):
   #there might have been no obs:
   if EDA2_obs_time!=0:
      print(EDA2_obs_time)
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time)
      #lst_eda2_hrs = "%0.1f" % get_eda2_lst("2015%s" % EDA2_obs_time[4::])
      #print(lst_eda2_hrs)
      lst_hrs_list.append(lst_eda2_hrs)
   else:
      #just use the first LST
      #lst_eda2_hrs = "%0.1f" % get_eda2_lst("2015%s" % EDA2_obs_time_list[0][4::])
      lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time_list[0])
      #print(lst_eda2_hrs)
      lst_hrs_list.append(lst_eda2_hrs)

#print(lst_hrs_list)


#EDA2
#array_layout_filename = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_255_NEU.ant'
#array_layout_filename = "/md0/code/git/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_NEU.txt"
#plot_antenna_array(array_layout_filename=array_layout_filename)
#sys.exit()
#plot_baseline_length_counts(array_layout_filename = array_ant_locations_filename,freq_MHz=50.)
#sys.exit()

#Step 1: simulate

#need to re-simulate diffuse_global_noise?


#cd into each chan dir separately and run simulate to get apparent sky images (don't worry that it crashes on concat freq step)
#need to fix this so you can just run like the other functoins below for multiple eda2 chans
#for EDA2 chans [64,77,90,103,116,129], freq_MHz_list = [ 50.  60.  70.  80.  91. 101.]
#sims:
#freq_MHz_list=[70.]
#lst_hrs_list = ['0']
#freq_MHz_list = np.arange(start_chan,start_chan+n_chan,chan_step)
#freq_MHz_array = np.asarray(freq_MHz_list)
#lst_hrs_list=['2']
#do it here: /md0/EoR/ASSASSIN/solve_for_tsky_weighted/jack_tests/single_point
#simulate(lst_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label,EDA2_data=False)
#uvfitsname1 = "eda_model_LST_000_X_70.000_MHz_SP.uvfits"
#uvfitsname2 = "ben_test_band01.uvfits"
#uvfitsname3 = "chan_94_20200303T135201.uvfits"
#uvfitsname4 = "eda_model_LST_030_X_70_MHz_DG_gsm.uvfits"
#compare_uvfits(uvfitsname1,uvfitsname2)
#sys.exit()

#DATA: (repeat twice with 'diffuse' then 'global_unity')
#pol_list = ['Y']
#chan_num = 0
#freq_MHz_list = [freq_MHz_array[chan_num]]
####if FAST: for data need to simulate with 'global_unity' and then separately 'diffuse' (only if fast)
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
#sys.exit()

#Step 2: calibrate





##do this outside chan dir
##if doing individual chans:
#chan_num = 1
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]
#freq_MHz_list = freq_MHz_array[0:2]
#EDA2_chan_list = EDA2_chan_list[0:2]
plot_cal = False
#wsclean = False
wsclean = True
concat=True
per_chan_cal = True
reverse_fine_chans = False   #this should always be false!

#uvfits_name = "64/chan_64_20200303T133729.uvfits"
#get_antenna_table_from_uvfits(uvfits_name)
#sys.exit()
##Daniels sims and EDA2 positions from Marcin and antenna numbering match!!!!
#ant_pos_on_ground_filename = '/md0/code/git/ben-astronomy/EoR/ASSASSIN/antenna_locations_marcin.txt'
#ant_pos_sim_filename = '/md0/code/git/ben-astronomy/EoR/ASSASSIN/EDA2_loc_190605.txt'
combined_ant_pos_name_filename = '/md0/code/git/ben-astronomy/EoR/ASSASSIN/ant_pos_eda2_combined_on_ground_sim.txt'
#cross_match_eda2_ant_pos_with_sims(ant_pos_on_ground_filename,ant_pos_sim_filename,combined_ant_pos_name_filename)
#sys.exit()

#simulate_sitara(lst_hrs_list[0],nside=32)
#sys.exit()

##unity only sim takes 2 min with nside 32, 6 mins with nside 64, similar 
chan_num = 0
plot_from_saved = False
simulate_eda2_with_complex_beams([freq_MHz_list[chan_num]],lst_hrs_list[chan_num],nside=32,plot_from_saved=plot_from_saved,EDA2_chan=EDA2_chan_list[chan_num],EDA2_obs_time=EDA2_obs_time_list[chan_num],n_obs_concat=n_obs_concat_list[chan_num])
pt_source=False
#currently hacked so xy and yx are zero ... gives correct flux scale, if they are non zero, corrected image is roughly twice as bright
#need to read stara paper comments and refs from refereee on polarisation!
#try different baseline length cuts for calibration and see how that affects the flux scale ...
write_to_miriad_vis([freq_MHz_list[chan_num]],lst_hrs_list[chan_num],EDA2_chan=EDA2_chan_list[chan_num],EDA2_obs_time=EDA2_obs_time_list[chan_num],n_obs_concat=n_obs_concat_list[chan_num],pt_source=pt_source)
model_ms_name= "20200303T133733_50.000.ms"
eda2_ms_name = "/md0/EoR/EDA2/20200303_data/64/20200303_133733_eda2_ch32_ant256_midday_avg8140.ms" 
calibrate_with_complex_beam_model(model_ms_name=model_ms_name,eda2_ms_name=eda2_ms_name)
sys.exit()

#ant_index = 0
#convert_matlab_EEPs_by_ant(ant_index,freq_MHz_list,nside=512,method='cubic')
#freq_MHz_list_index = 0
#convert_matlab_EEPs_by_freq(freq_MHz_list,freq_MHz_list_index,nside=512)


#don't need to cal separately each pol anymore, using wsclean predict and calibrate!
#for pol in pol_list:
pol_list_input = []
chan_num = 0
#New cal Jan 2021 - try to average data in time first before cal
#2 Feb try withinitial full BW cal
#calibrate_eda2_data_time_av(EDA2_chan_list=EDA2_chan_list,obs_type='night',lst_list=lst_hrs_list,pol_list=pol_list_input,n_obs_concat_list=n_obs_concat_list,concat=concat,wsclean=wsclean,plot_cal=plot_cal,uv_cutoff=0,per_chan_cal=per_chan_cal)
#sys.exit()

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
#model_type_list = ['OLS_fixed_intercept']
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
#pol_list = ['X','Y']
#chan_num = 1
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]


#for sims:
#freq_MHz_list = np.arange(start_chan,start_chan+n_chan,chan_step)
#freq_MHz_array = np.asarray(freq_MHz_list)

#lst_hrs_list=['4']
#poly_order_list=[5,6,7]
#poly_order=5

#plot_iso_ant_int_response()
#sys.exit()
 
#plot_only = False
#baseline_length_thresh_lambda = 0.5
#include_angular_info = True

###WODEN sims signal extraction
#woden=False
#wsclean=True
#fast=True
#no_modelling=False
#calculate_uniform_response=False
#noise_coupling=False
##woden sims from 50 to 193 MHz
##dont need this any more
##woden: centre_chans_number_list = [63,87,111,135,159,183]
#
##make sourcelist and sbatch files on namorrodor:      
##new centre chans list for 1 MHz wide chans
##UTC time: woden sims for LST 60 zenith "2015-11-29T15:33:43"
##year,month,day,hour,min,sec = 2020,03,03,15,00,00
#freq_MHz_list = np.arange(0,218)*1.28 + 50.
#freq_MHz_list = freq_MHz_list[0:100]
#print(freq_MHz_list)
#year,month,day,hour,min,sec = 2015,11,29,15,40,29 #LST=60 deg
##year,month,day,hour,min,sec = 2015,11,29,15,33,43
#time_string = '%d_%02d_%02d_%02d_%02d_%02d' % (year,month,day,hour,min,sec)
##write_woden_skymodels(freq_MHz_list,NSIDE,time_string,dipole_height_m,pol_list)
##write_woden_sims_sbatch_file(freq_MHz_list,time_string,pol_list)
##sys.exit()
#
##freq_MHz_list = np.asarray(woden_chan_list) + (np.arange(0,24)-13) 
###just use the orig freq list from 50 to 199 MHz
#outbase_name = 'lst_%0.2f_hr' % (float(lst_hrs_list[0]))
#print(freq_MHz_list)
###simulate to get the theoretical beam weighted global signal 
#simulate(lst_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,outbase_name=outbase_name,array_ant_locations_filename=array_ant_locations_filename,array_label=array_label,EDA2_data=False)

#plot_tsky_for_multiple_freqs(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list_input,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order=poly_order,plot_only=plot_only,include_angular_info=include_angular_info,model_type_list=model_type_list, EDA2_data=EDA2_data,EDA2_chan_list=EDA2_chan_list,n_obs_concat_list=n_obs_concat_list,wsclean=wsclean,fast=fast,no_modelling=no_modelling,calculate_uniform_response=calculate_uniform_response,woden=woden,noise_coupling=noise_coupling)

#sys.exit()


#############################################################################################
model_type_list = ['OLS_fixed_intercept','OLS_fixed_int_subtr_Y']
#model_type_list = ['OLS_fixed_intercept']
pol_list_input = ['X','Y']
poly_order=5
plot_only = True
baseline_length_thresh_lambda = 0.5
include_angular_info = True
woden=False
wsclean=True
fast=True
no_modelling=True
calculate_uniform_response=False
noise_coupling=False

#have a look at the auto and cross power as a function of frequency for selected baselines and see if continuous
#do this also in calibtate loop above for sanity check - am I calibrating twice on corrected data?? (you are, but does it matter?)#baseline_number = 10
#pick angular or unity or neither (not both)
#baseline_number = 10
#unity = False
#angular=False
#inspect_cross_auto_power(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list_input,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order=poly_order,EDA2_data=EDA2_data,EDA2_chan_list=EDA2_chan_list,n_obs_concat_list=n_obs_concat_list,wsclean=wsclean,fast=fast,woden=woden,noise_coupling=noise_coupling,baseline_number=baseline_number,unity=unity,angular=angular)
#sys.exit()

#
#EDA 2 data:
#Plot subset of freqs?
#chan_num = 90 - 64 #90 = 70MHz
#this is the start coarse chan
chan_num = 0
n_coarse_chans_to_plot = 63     #127 - 64 = 63 (all coarse chans)
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]
freq_MHz_list = freq_MHz_array[chan_num:chan_num+n_coarse_chans_to_plot]
EDA2_chan_list = EDA2_chan_list[chan_num:chan_num+n_coarse_chans_to_plot]
#
plot_tsky_for_multiple_freqs(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list_input,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order=poly_order,plot_only=plot_only,include_angular_info=include_angular_info,model_type_list=model_type_list, EDA2_data=EDA2_data,EDA2_chan_list=EDA2_chan_list,n_obs_concat_list=n_obs_concat_list,wsclean=wsclean,fast=fast,no_modelling=no_modelling,calculate_uniform_response=calculate_uniform_response,woden=woden,noise_coupling=noise_coupling)
sys.exit()


#sims - new idea: just always do sim plots when analysing data  - I think we already have the miriad sims after the calibrate stage!

#plot_only = False
#include_angular_info = True
#no_modelling=False
#up to here with plot_only = False
#chan_num = 90 - 64 #90 = 70MHz
#chan_num = 20
#freq_MHz_list = [freq_MHz_array[chan_num]]
#EDA2_chan_list = [EDA2_chan_list[chan_num]]
#freq_MHz_list = freq_MHz_array[chan_num:chan_num+35]
#EDA2_chan_list = EDA2_chan_list[chan_num:chan_num+35]
#wsclean=False # for sims or miriad cal
#sim for paper plot 1 
#wsclean=True # for data
#fast=False
#calculate_uniform_response=False


#plot_tsky_for_multiple_freqs(lst_hrs_list=lst_hrs_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model,array_label=array_label,baseline_length_thresh_lambda=baseline_length_thresh_lambda,poly_order=poly_order,plot_only=plot_only,include_angular_info=include_angular_info,model_type_list=model_type_list, EDA2_data=EDA2_data,EDA2_chan_list=EDA2_chan_list,n_obs_concat_list=n_obs_concat_list,wsclean=wsclean,fast=fast,no_modelling=no_modelling,calculate_uniform_response=calculate_uniform_response,woden=woden,noise_coupling=noise_coupling)

#Need to change colors of plots throughout so they are suitable for color blindness and use dotted, dashed, or dot dashed lines instead of just colours (and different symbols in scatter plots)
#See: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40     for colors to use

#https://gist.github.com/thriveth/8560036  
#CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
                  
                  
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

#beam normalisation
freq_MHz_list = np.arange(50,300)
LNA_impedance_list =[]
beam_norm_list = []
for freq_MHz in freq_MHz_list:
   LNA_impedance = calc_LNA_impedance_eda2(freq_MHz)
   LNA_impedance_list.append(LNA_impedance)
   beam_norm = calc_beam_normalisation(freq_MHz,LNA_impedance)
   beam_norm_list.append(beam_norm)
   #print(LNA_impedance)

LNA_impedance_list = np.asarray(LNA_impedance_list,dtype=complex)
beam_norm_list = np.asarray(beam_norm_list,dtype=complex)

plt.clf()
map_title=""
#######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
plt.plot(freq_MHz_list,beam_norm_list.real)
fig_name="beam_norm_real.png" 
figmap = plt.gcf()
figmap.savefig(fig_name)
print("saved %s" % fig_name) 

plt.clf()
map_title=""
#######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
plt.plot(freq_MHz_list,beam_norm_list.imag)
fig_name="beam_norm_imag.png" 
figmap = plt.gcf()
figmap.savefig(fig_name)
print("saved %s" % fig_name) 

#compare to Sokolowski et al 2017 MWA beam paper
plt.clf()
map_title=""
#######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
plt.plot(freq_MHz_list,LNA_impedance_list.real)
fig_name="lna_impedance_real.png" 
figmap = plt.gcf()
figmap.savefig(fig_name)
print("saved %s" % fig_name) 

plt.clf()
map_title=""
#######hp.orthview(map=rotated_power_pattern,half_sky=True,rot=(0,float(mwa_latitude_ephem),0),title=map_title)
plt.plot(freq_MHz_list,LNA_impedance_list.imag)
fig_name="lna_impedance_imag.png" 
figmap = plt.gcf()
figmap.savefig(fig_name)
print("saved %s" % fig_name) 



