#!/usr/bin/env python
#ASSASSIN: (All Sky SigAl Short SPacing INterferometer)
#Script to replicate Cath's global step simulations

#Daniel:
# The format of the data is in standard spherical coordinates (theta, phi).
# Where theta is the zenith angle and phi is anti-clockwise going from east to north looking down at the array. 
# The 'origin' will be on the top left corner. Going from top to bottom is increasing in phi and going left to right is increasing in theta. 
# The data is in 1 deg steps.

import numpy as np
import healpy as hp
#from pygsm import GSMObserver2016
#from pygsm import GSMObserver
from pygsm import GlobalSkyModel
from pygsm import GlobalSkyModel2016
from datetime import datetime, date
import matplotlib.pyplot as plt
import os,sys
from reproject import reproject_from_healpix
import pyfits
from astropy.wcs import WCS
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from scipy import signal
import numpy.polynomial.polynomial as poly
from pyuvdata import UVData
import math


#mu_0 = 4.*np.pi*10**(-7)
c = 3.0e8
k = 1.38065e-23
sq_deg_in_1_sr = (180./math.pi)**2
#Need to specify the location of the observatory (MWA ... put in exact EDA later)
# Setup observatory location - in this case, MWA, Australia
latitude_degrees=-26.70331940
longitude_degrees=116.67081524
elevation_m=377.83
#for uvgen systemp values, sky temp at 180 MHz, beta is spectral index (Frank's 180 at 180), efficiency eta
T_180 = 180.0
beta = -2.5
eta = 1

#pointing dec for uvgen (needs to be zenith at each lst)
pointing_dec = "-26.7"


#global signal params
#From Cath feasibility study (look at edges papers?)
#S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
C = 0.08/1000.
A = -150./1000.
delta_nu = 10.
nu_c = 70.

use_analytic_beam = True
generate_new_hpx = True
plot_gsm_map_hpx = False
generate_new_vis = True
plot_global_signal_map_hpx = False

#This is for Daniels FEKO model beams
generate_new_average_beam = False
#apply_normalisation = True
plot_all_beams = False
plot_average_beam = False

#specify date
#date_time_string = '2018_09_25_00_00_00'
#lst_list = ['00', '03', '06', '09', '12', '15', '18', '21']
#lst_list = ['00']
lst_list_array = np.arange(0,2,0.066667)
#lst_list_array = np.arange(0,0.066667,0.066667)
#lst_list_array = np.arange(0,0.1333333,0.066667)
lst_list = lst_list_array.astype(str)
#lst_list = ['0.120']
#lst_list = ['12']
#pol_list = ['X','Y']
sky_model = 'gsm'
#sky_model = 'gsm2016'
pol_list = ['X']
#signal_type_list=['global','diffuse','noise','gain_errors']
signal_type_list=['diffuse']
gsm_smooth_poly_order = 5
poly_order = 8
#freq_MHz_list = np.arange(50,200,1)
start_chan=50
n_chan=150
chan_step = 1
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

n_ants = 256
n_baselines = n_ants*(n_ants-1) / 2.
template_imsize = 512
template_cell_size_asec = 180./template_imsize*60.*60.
#amount to zero pad beam image for FFT
#want this to be big so we have good resolution in Fourier space (needed for close to zeros spacings)
padding_factor = 25
uvdist_wavelength_cutoff = 1.0
zero_spacing_leakage_threshold = 0.1


#from Wayht et al Table 2, third column divided by 256 (not last col which is dipole in isolation)
Aeff_list = np.array([970,950,914,874,832,771,707,638,568,498,435,377,329,288,252,222,196])/256.
Aeff_freqs = np.arange(60,230,10)

#fit polynomial 
z = np.polyfit(Aeff_freqs, Aeff_list, 3)
p = np.poly1d(z)
Aeff_for_freq_MHz_list = p(freq_MHz_list)

def cleanup_images_and_vis(lst,freq_MHz,pol):
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
  
   gsm_hpx_fits_name = "%s_map_LST_%03d_%s_MHz_hpx.fits" % (sky_model,lst_deg,int(freq_MHz))
   reprojected_gsm_fitsname = "%s_map_LST_%03d_%s_MHz_reprojected.fits" % (sky_model,lst_deg,int(freq_MHz))
   reprojected_gsm_im_name = "%s_map_LST_%03d_%s_MHz_reprojected.im" % (sky_model,lst_deg,int(freq_MHz))
   
   global_signal_hpx_fits_name = "global_signal_map_LST_%03d_%s_MHz_hpx.fits" % (lst_deg,int(freq_MHz))
   reprojected_global_signal_fitsname = "global_signal_map_LST_%03d_%s_MHz_reprojected.fits" % (lst_deg,int(freq_MHz))
   reprojected_global_signal_im_name = "global_signal_map_LST_%03d_%s_MHz_reprojected.im" % (lst_deg,int(freq_MHz))
   
   #lna_impedance_aavs1_filename = "/data/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
   
   base_vis_name = "eda_model_LST_%03d_%s_MHz.vis" % (lst_deg,int(freq_MHz))
   eda_model_no_source_image_name = "eda_no_source_LST_%03d_%s_MHz.map" % (lst_deg,int(freq_MHz))
   eda_model_no_source_beam_name = "eda_no_source_LST_%03d_%s_MHz.beam" % (lst_deg,int(freq_MHz))
   eda_model_no_source_image_fits_name = "eda_no_source_LST_%03d_%s_MHz.fits" % (lst_deg,int(freq_MHz))
   
   reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%s_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,int(freq_MHz))
   reprojected_global_signal_im_Jy_per_pix_name =  "global_map_%s_%s_MHz_reproj_Jy_pix.im" % (date_time_string,int(freq_MHz))
   
   cmd = "rm -rf %s %s %s %s %s %s %s %s %s %s %s %s" % (gsm_hpx_fits_name,reprojected_gsm_fitsname,reprojected_gsm_im_name,global_signal_hpx_fits_name,reprojected_global_signal_fitsname,reprojected_global_signal_im_name,base_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name,reprojected_gsm_im_Jy_per_pix_name,reprojected_global_signal_im_Jy_per_pix_name)
   print cmd
   os.system(cmd)
   
   apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%s_MHz.im" % (lst_deg,pol,int(freq_MHz))
   apparent_global_signal_im_name = "apparent_global_signal_LST_%03d_%s_pol_%s_MHz.im" % (lst_deg,pol,int(freq_MHz))
   cmd = "rm -rf %s %s" % (apparent_sky_im_name,apparent_global_signal_im_name)
   print cmd
   os.system(cmd)
   
   model_sky_vis = "eda_model_plus_sky_LST_%03d_%s_pol_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
   model_global_signal_vis = "eda_model_plus_global_LST_%03d_%s_pol_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
   eda_model_noise_vis_name = "eda_model_noise_LST_%03d_%s_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
   cmd = "rm -rf %s %s %s" % (model_sky_vis,model_global_signal_vis,eda_model_noise_vis_name)
   print cmd
   os.system(cmd)
         
   beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%s_MHz.im' % (pol,int(freq_MHz))
   beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%s_MHz_gsm_regrid.im' % (pol,int(freq_MHz))
   cmd = "rm -rf %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name)
   print cmd
   os.system(cmd)        

   one_jy_source_vis_name =  "one_jy_source_LST_%03d_%s_MHz.vis" % (lst_deg,int(freq_MHz))
   cmd = "rm -rf %s " % (one_jy_source_vis_name)
   print cmd
   os.system(cmd)

#for lst in lst_list: 
#   cleanup_images_and_vis(lst)
#sys.exit()   

      
def plot_S21(nu_array=None,C=0.,A=1.,delta_nu=20.,nu_c=78.):
   #Global signal
   #Use S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
   S_21 = C*nu_array + A*np.exp(-(nu_array - nu_c)**2/(2*(delta_nu)**2))
   plt.clf()
   map_title="S_21 vs freq"
   plt.plot(nu_array,S_21)
   plt.ylabel("Tb (K)")
   plt.xlabel("freq (MHz)")
   fig_name="s_21_vs_freq.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   return S_21

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
   print "wrote image %s" %  fitsname

def recompute_ffts(pol,freq_MHz):

   if pol=='X':
      beam_image_name = "model_%s_MHz_xx.fits" % int(freq_MHz)
   else:
      beam_image_name = "model_%s_MHz_yy.fits" % int(freq_MHz)
   
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
   print beam_fft2_imag
   print np.max(abs(beam_fft2_imag))
   beam_fft2_abs = abs(beam_fft2)

   #beam_fft2_real_shift_norm = beam_fft2_real_shift/np.max(beam_fft2_real_shift)
   pyfits.writeto(fits_name_real,beam_fft2_real,clobber=True)
   print "saved %s" % fits_name_real
   
   pyfits.writeto(fits_name_imag,beam_fft2_imag,clobber=True)
   print "saved %s" % fits_name_imag

   pyfits.writeto(fits_name_abs,beam_fft2_abs,clobber=True)
   print "saved %s" % fits_name_abs 
   #save the fft? - these are like a gig each - dont save!
   #np.save(fft_savename,beam_fft2)
   #np.save(fft_savename_shift,beam_fft2_shift)
   #print "saved %s" % fft_savename_shift
   #print "saved %s" % fft_savename
      
def plot_beam_and_weights(pol,freq_MHz):
   wavelength = 300./freq_MHz
   profile_plot_length = 500
   
   if pol=='X':
      beam_image_name = "model_%s_MHz_xx.fits" % int(freq_MHz)
   else:
      beam_image_name = "model_%s_MHz_yy.fits" % int(freq_MHz)
   
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
   print "saved %s" % fig_name


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
   print "saved %s" % fig_name
   
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
   print "saved %s" % fig_name

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
   print "saved %s" % fig_name
      
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
   print "saved %s" % fig_name

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
   print "saved %s" % fig_name

   #abs
   
   plt.clf()
   map_title="beam fft2 abs"
   plt.imshow(beam_fft2_abs_shift)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= beam_plot_basename + "_beam_fft2_abs.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   
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
   print "saved %s" % fig_name

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
   print "saved %s" % fig_name
      
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
   print "saved %s" % fig_name

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
   print "saved %s" % fig_name
   
   
def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return (index,array[index])
   
   
#def compute_weights(lst_deg,pol,freq_MHz,recompute_fits=False):
def compute_weights(lst_list,freq_MHz_list,pol_list,signal_type_list=['diffuse','global','noise','gain_errors'],sky_model='gsm'):
   n_chan = len(freq_MHz_list)
   n_lst = len(lst_list)
   n_pol = len(pol_list)
   freq_MHz_end = freq_MHz_list[-1]
   lst_end = lst_list[-1]
   lst_deg_end = (float(lst_end)/24.)*360.
   for pol_index,pol in enumerate(pol_list):
      model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg_end,pol,int(freq_MHz_end))
      if 'noise' in signal_type_list:
          model_vis_name_base += '_noise'
      if 'diffuse' in signal_type_list:
          model_vis_name_base += '_diffuse_%s' % sky_model
      if 'global' in signal_type_list:
          model_vis_name_base += '_global' 
      if 'gain_errors' in signal_type_list:
          model_vis_name_base += '_gain_errors'
                                   
      uvfits_name = "%s_concat_lsts.uvfits" % model_vis_name_base
      print uvfits_name           
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
               beam_image_name = "model_%s_MHz_xx.fits" % int(freq_MHz)
            else:
               beam_image_name = "model_%s_MHz_yy.fits" % int(freq_MHz)
            
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
                  #print "saved %s" % fig_name
              
                  
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
                  
                  #put this check in the extract_signal function (so we still record the smaller and negative weights and can plot them to check against the beam )
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
   

#for lst in lst_list[0:2]:
#   lst_deg = (float(lst)/24.)*360.
#   for pol in pol_list:
#      for freq_MHz in freq_MHz_list[0:2]:
#         compute_weights(lst_deg,pol,freq_MHz)

def model_signal(lst_list,freq_MHz_list,pol_list,signal_type_list=['diffuse','global','noise','gain_errors'],outbase_name='extracted_1_',poly_order=4,sky_model='gsm'):
   #for signal_type in signal_type_list:
   freq_MHz_end = freq_MHz_list[-1]
   lst_end = lst_list[-1]
   lst_deg_end = (float(lst_end)/24.)*360.
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
         model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg_end,pol,int(freq_MHz_end))
         if 'noise' in signal_type_list:
             model_vis_name_base += '_noise'
         if 'diffuse' in signal_type_list:
             model_vis_name_base += '_diffuse_%s' % sky_model
         if 'global' in signal_type_list:
             model_vis_name_base += '_global'
         if 'gain_errors' in signal_type_list:
             model_vis_name_base += '_gain_errors'         

         #signal_array_short_baselines_filename = outbase_name + "%s_signal.npy" % (model_vis_name_base)
         #signal_array_short_baselines = np.load(signal_array_short_baselines_filename)
         #signal_array_all_baselines_filename = outbase_name + "%s_signal_all_baselines.npy" % (model_vis_name_base)
         #signal_array_all_baselines = np.load(signal_array_all_baselines_filename)
         #signal_array_all_baselines_filename_abs = outbase_name + "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
         #signal_array_all_baselines_abs = np.load(signal_array_all_baselines_filename_abs)
         #signal_array_short_baselines_weighted_filename = outbase_name + "%s_signal_weighted.npy" % (model_vis_name_base)
         #signal_array_short_baselines_weighted = np.load(signal_array_short_baselines_weighted_filename)

          
         signal_array_short_baselines_Tb_filename = outbase_name + "%s_signal_Tb.npy" % (model_vis_name_base)
         signal_array_short_baselines_Tb = np.load(signal_array_short_baselines_Tb_filename)
         signal_array_all_baselines_Tb_filename = outbase_name + "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines_Tb = np.load(signal_array_all_baselines_Tb_filename)
         signal_array_all_baselines_filename_abs_Tb = outbase_name + "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines_abs_Tb = np.load(signal_array_all_baselines_filename_abs_Tb)
         signal_array_short_baselines_weighted_Tb_filename = outbase_name + "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
         signal_array_short_baselines_weighted_Tb = np.load(signal_array_short_baselines_weighted_Tb_filename)
        
         #######Unweighted
         
         #can't use the nan values:
         good_idx = np.isfinite(signal_array_short_baselines_Tb[0,:]) & np.isfinite(signal_array_short_baselines_Tb[0,:])
         good_signal_array_short_baselines_Tb = signal_array_short_baselines_Tb[0,:][good_idx]
         
         good_freq_MHz_array = freq_MHz_array[good_idx]
         
         coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
         ffit = poly.polyval(good_freq_MHz_array, coefs)

         plt.clf()
         plt.plot(good_freq_MHz_array,ffit,label='model fit')
         plt.plot(freq_MHz_array,signal_array_short_baselines_Tb[0,:],label='data')
         map_title="Polynomial order %s fit %s" % (poly_order,sky_model)
         plt.ylabel("Extracted Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_model_short_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name

         residual = ffit - good_signal_array_short_baselines_Tb

         
         plt.clf()
         plt.plot(good_freq_MHz_array,residual,label='residual')
         map_title="Residual for polynomial order %s fit to diffuse" % poly_order
         plt.ylabel("Residual Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_residual_short_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name

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
         #fig_name= "%s%s_log_fit_residual_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print "saved %s" % fig_name         
         
         #convert back to linear
         #inverse_log_residual = 10.**log_residual
         
         plt.clf()
         plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual of log fit')
         map_title="Residual for log polynomial order %s fit to diffuse" % poly_order
         plt.ylabel("Residual Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_log_fit_residual_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name 

         #####Weighted
         #can't use the nan values:
         good_idx = np.isfinite(signal_array_short_baselines_weighted_Tb[0,:]) & np.isfinite(signal_array_short_baselines_Tb[0,:])
         good_signal_array_short_baselines_Tb = signal_array_short_baselines_weighted_Tb[0,:][good_idx]
         
         good_freq_MHz_array = freq_MHz_array[good_idx]
         
         coefs = poly.polyfit(good_freq_MHz_array, good_signal_array_short_baselines_Tb, poly_order)
         ffit = poly.polyval(good_freq_MHz_array, coefs)

         plt.clf()
         plt.plot(good_freq_MHz_array,ffit,label='model fit')
         plt.plot(freq_MHz_array,signal_array_short_baselines_weighted_Tb[0,:],label='data')
         map_title="Polynomial order %s fit %s" % (poly_order,sky_model)
         plt.ylabel("Eighted extracted Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_model_short_weighted_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name

         residual = ffit - good_signal_array_short_baselines_Tb

         
         plt.clf()
         plt.plot(good_freq_MHz_array,residual,label='residual')
         map_title="Weighted residual for polynomial order %s fit to diffuse" % poly_order
         plt.ylabel("Residual Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_residual_short_weighted_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name

         #in log log space:
         log_signal_array_short_baselines = np.log10(good_signal_array_short_baselines_Tb)
         log_freq_MHz_array = np.log10(good_freq_MHz_array)
         coefs = poly.polyfit(log_freq_MHz_array, log_signal_array_short_baselines, poly_order)
         ffit = poly.polyval(log_freq_MHz_array, coefs)
         ffit_linear = 10**ffit
         
         #residual = log_signal_array_short_baselines - log_ffit
         residual_of_log_fit = ffit_linear - good_signal_array_short_baselines_Tb
         
         #plt.clf()
         #plt.plot(log_freq_MHz_array,log_residual,label='log residual')
         #map_title="Weighted log residual for polynomial order %s fit to diffuse" % poly_order
         #plt.ylabel("Log residual Tb (K)")
         #plt.xlabel("log freq (MHz)")
         #plt.legend(loc=1)
         #fig_name= "%s%s_log_residual_weighted_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         #figmap = plt.gcf()
         #figmap.savefig(fig_name)
         #print "saved %s" % fig_name         
         
         #convert back to linear
         #inverse_log_residual = 10.**log_residual
         
         plt.clf()
         plt.plot(good_freq_MHz_array,residual_of_log_fit,label='residual from log fit')
         map_title="Weighted residual for log polynomial order %s fit to diffuse" % poly_order
         plt.ylabel("Residual Tb (K)")
         plt.xlabel("freq (MHz)")
         plt.legend(loc=1)
         fig_name= "%s%s_log_fit_residual_weighted_poly_%s.png" % (outbase_name,model_vis_name_base,poly_order)
         figmap = plt.gcf()
         figmap.savefig(fig_name)
         print "saved %s" % fig_name 
         
         
def extract_signal_from_sims(lst_list,freq_MHz_list,pol_list,signal_type_list=['diffuse','global','noise','gain_errors'],outbase_name='extracted_1_',sky_model='gsm'):
   #for signal_type in signal_type_list:
   freq_MHz_end = freq_MHz_list[-1]
   lst_end = lst_list[-1]
   lst_deg_end = (float(lst_end)/24.)*360.
   freq_MHz_array = np.asarray(freq_MHz_list)
   for pol_index,pol in enumerate(pol_list):
         model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg_end,pol,int(freq_MHz_end))
         if 'noise' in signal_type_list:
             model_vis_name_base += '_noise'
         if 'diffuse' in signal_type_list:
             model_vis_name_base += '_diffuse_%s' % sky_model
         if 'global' in signal_type_list:
             model_vis_name_base += '_global' 
         if 'gain_errors' in signal_type_list:
             model_vis_name_base += '_gain_errors'
                                      
         uvfits_name = "%s_concat_lsts.uvfits" % model_vis_name_base
         print uvfits_name           
         hdulist = fits.open(uvfits_name)
         #hdulist.info()
         info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
         #print info_string
         uvtable = hdulist[0].data
         uvtable_header = hdulist[0].header
         visibilities = uvtable['DATA']
         print visibilities.shape

               
         signal_array_short_baselines = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_filename = outbase_name + "%s_signal.npy" % (model_vis_name_base)
         signal_array_short_baselines_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_Tb_filename = outbase_name + "%s_signal_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename = outbase_name + "%s_signal_all_baselines.npy" % (model_vis_name_base)
         signal_array_all_baselines_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_Tb_filename = outbase_name + "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
         signal_array_all_baselines_abs = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename_abs = outbase_name + "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
         signal_array_all_baselines_abs_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_all_baselines_filename_abs_Tb = outbase_name + "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
         signal_array_short_baselines_weighted = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_weighted_filename = outbase_name + "%s_signal_weighted.npy" % (model_vis_name_base)
         signal_array_short_baselines_weighted_Tb = np.full((n_pol,n_chan),0.0)
         signal_array_short_baselines_weighted_Tb_filename = outbase_name + "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
         number_baselines_used_array = np.full((n_pol,n_chan),0.0)
         number_baselines_used_array_filename = outbase_name + "%s_number_baselines_used.npy" % (model_vis_name_base)
         sum_of_weights_all_baselines_array = np.full((n_pol,n_chan),np.nan)
         sum_of_weights_all_baselines_array_filename = outbase_name + "%s_sum_of_weights_all_baselines.npy" % (model_vis_name_base)
         sum_of_weights_short_baselines_array = np.full((n_pol,n_chan),np.nan)
         sum_of_weights_short_baselines_array_filename = outbase_name + "%s_sum_of_weights_short_baselines.npy" % (model_vis_name_base)
         
         #chan_sum_of_weights = 0.0
         #chan_vis_real_sum_weighted = 0.0
         for chan_index,freq_MHz in enumerate(freq_MHz_list):
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
               #if signal_type == 'sky':
               #   uvfits_name = "eda_model_plus_sky_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))
               #elif signal_type == 'global':
               #   uvfits_name = "eda_model_plus_global_LST_%03d_%s_pol_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz)) 
               #elif signal_type == 'noise':
               #   uvfits_name = "eda_model_noise_LST_%03d_%s_%s_MHz.uvfits" % (lst_deg,pol,int(freq_MHz))  
               #else:
               #   print "unrecognised signal type, can not extract signal"
               #   sys.exit()

               
               #get the UU and VV so we can check whether we are using short baselines
               UU_s_array = uvtable['UU']
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV']
               VV_m_array = VV_s_array * c
               
               #print "VV_m_array.shape" 
               #print VV_m_array.shape
               
               #print visibilities.shape
               n_vis = visibilities.shape[0]
               n_timesteps = n_vis/n_baselines
               #print "n_timesteps %s " % n_timesteps
               timestep_array = np.arange(0,n_timesteps,1)
               #weights_array_filename = "weights_LST_%03d_%s_%s_MHz.npy" % (lst_deg,pol,int(freq_MHz))
               weights_array_filename = "weights_%s.npy" % model_vis_name_base
               weights_array = np.load(weights_array_filename)
                      
               #for pol_index,pol in enumerate(pol_list):
               #   for chan_index,freq_MHz in enumerate(freq_MHz_list):
               
               wavelength = 300./freq_MHz
               #print freq_MHz
               #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
               timestep_vis_real_sum = 0.
               timestep_vis_real_sum_weighted = 0.
               for timestep in timestep_array:
                  timestep_baselines_used = 0.
                  vis_start_index = int(timestep*n_baselines)
                  timestep_visibilities = visibilities[vis_start_index:int(vis_start_index+n_baselines),0,0,0,chan_index,0,:]
                  #print timestep_visibilities.shape
                  for timestep_visibility_index,timestep_visibility_real in enumerate(timestep_visibilities[:,0]):
                     #this is wrong as the visibilities are ungridded
                     ##print visibility_real
                     #print timestep_visibilities.shape
                     timestep_visibility_imag = timestep_visibilities[timestep_visibility_index,1]
                     ###visibility_weight = visibilities[visibility_index,0,0,0,2]
                     ##print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
                     complex_vis = np.complex(timestep_visibility_real,timestep_visibility_imag)
                     abs_vis = abs(complex_vis)
                     ##print "complex_vis %s" % complex_vis
                     ##signal_array_unweighted[pol_index,chan_index,0] += complex_vis
                     #if signal_type == 'sky':
                     #   #sky_signal_array_all_baselines[pol_index,chan_index,0] += abs_vis
                     #elif signal_type == 'global':
                     #   global_signal_array_all_baselines[pol_index,chan_index,0] += abs_vis
                     #elif signal_type == 'noise':
                     #   noise_array_all_baselines[pol_index,chan_index,0] += abs_vis
                     #else:
                     #   print "unrecognised signal type, cannot extract signal"
                     #   sys.exit()
                     
                     
                     #weighted: (this is equivalent to gridding), finding the value that the visibility would be it u,v=0 after convolving with a Fourier beam kernel:
                     visibility_index = vis_start_index + timestep_visibility_index
                     #print "visibility_index" 
                     #print visibility_index
                     uv_zero_weighting = weights_array[visibility_index,lst_index,chan_index,pol_index]
                     #print "weights_array shape" 
                     #print weights_array.shape
                     #all baselines:
                     chan_vis_real_sum_all_baselines += timestep_visibility_real
                     chan_vis_abs_sum_all_baselines += abs_vis
                     #only short baselines
                     if not np.isnan(uv_zero_weighting):
                        #all baselines signal add to weights array                        
                        chan_sum_of_weights_all_baselines += uv_zero_weighting
                        if uv_zero_weighting > 0.1:
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
                           #  print "timestep_visibility_real %s at LST %s %s MHz" % (timestep_visibility_real,lst,freq_MHz)
                           #  print "for baseline U %s V %s at visibility index %s and uv_zero_weighting %s" % (UU_m_array[visibility_index],VV_m_array[visibility_index],visibility_index,uv_zero_weighting)
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
               print "chan_vis_real_sum"
               print chan_vis_real_sum
               #chan_vis_real_weighted_average = chan_vis_real_sum_weighted/chan_sum_of_weights
            
               number_baselines_used_average =   number_baselines_used_sum / (n_timesteps*len(lst_list)) 
               print "av number baselines used for chan %s is %s" % (chan_index,number_baselines_used_average)
            
            
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
               #print "saved %s" % fig_name       
         
         #convert to brightness temp
         wavelength_array = 300./freq_MHz_array
         #print sum_of_weights_array
         jy_to_tb_all_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_all_baselines_array) 
         jy_to_tb_short_baselines = (wavelength_array**2) / (2. * k * 1.0e26 * sum_of_weights_short_baselines_array) 
         #jy_to_tb = (wavelength_array**2) / (2. * k * 1.0e26)
         signal_array_all_baselines_Tb = signal_array_all_baselines * jy_to_tb_all_baselines
         signal_array_short_baselines_Tb = signal_array_short_baselines * jy_to_tb_short_baselines
         signal_array_all_baselines_abs_Tb = signal_array_all_baselines_abs * jy_to_tb_all_baselines
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
   print info_string
   uvtable = hdulist[0].data
   uvtable_header = hdulist[0].header
   #print uvtable_header
   
   #super_threshold_indices = a > thresh
   #a[super_threshold_indices] = 0
   
   print uvtable['DATA'].shape
   real_visibilities_at_freq_MHz1 = uvtable['DATA'][:,0,0,0,0,0]
   abs_vis = abs(real_visibilities_at_freq_MHz1)
   log_abs_vis = np.log10(abs_vis)
   UU_s = uvtable['UU']
   UU_m = UU_s * c
   VV_s = uvtable['VV']
   VV_m = VV_s * c
   
   uvdist_m = np.sqrt(UU_m**2 + VV_m**2)
   
   print "min abs uvdist_m"
   print np.abs(uvdist_m).min()
   print "max abs uvdist_m"
   print np.abs(uvdist_m).max()
   
   
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
   print "saved %s" % fig_name

   plt.clf()
   map_title="real vis vs UU"
   plt.plot(uvdist_lambda,real_visibilities_at_freq_MHz1,linestyle='none',marker='.')
   plt.ylabel("real vis Jy")
   plt.xlabel("uvdist (lambda)")
   fig_name= plot_basename + "real_vis_vs_uvdist.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   
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
   print "saved %s" % fig_name

   
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
   #print "saved %s" % fig_name   
 
   #plt.clf()
   #map_title="log abs vis vs freq"
   #plt.plot(freq_MHz_array,abs_vis,linestyle='none',marker='.')
   #fig_name="abs_vis_vs_uu.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print "saved %s" % fig_name


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
   #sys.exit()
   
   fits.setval(first_filename, 'TELESCOP', value='MWA' )
   fits.setval(first_filename, 'OBJECT', value='zenith' )
   fits.setval(first_filename, 'CDELT4', value=1.000E+06)
   fits.setval(first_filename, 'DATE', value=1.0000)
   fits.setval(first_filename, '_DATE', value=0.5000)
   fits.setval(first_filename, 'INTTIM', value=240)
   print first_filename
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
   print "wrote %s" % outname
   #delete all the temporary uvfits files
   if delete:
      for uvfits_name in uvfits_list:
         cmd = "rm -rf %s" % (uvfits_name)
         print cmd
         os.system(cmd)
         #also delete the miriad vis file
         vis_name = uvfits_name.split('.')[0] + '.vis'
         cmd = "rm -rf %s" % (vis_name)
         print cmd
         os.system(cmd)
         
         
def concat_uvfits_casa(uvfits_list,outname,delete=False):
   #write a bash file for casa to read
   for uvfits_name in uvfits_list:
      vis_name = uvfits_name.split('.')[0] + '.ms'
      casa_string = "importuvfits(fitsfile='%s',vis='%s')" % (uvfits_name,vis_name)
      casa_filename = 'casa_concat.sh'
      with open(casa_filename,'w') as f:
         f.write(casa_string)
      cmd = "casa -c %s" % casa_filename
      print cmd
      os.system(cmd)

   if delete:
      for uvfits_name in uvfits_list:
         cmd = "rm -rf %s" % (uvfits_name)
         print cmd
         os.system(cmd)
         #also delete the miriad vis file
         vis_name = uvfits_name.split('.')[0] + '.vis'
         cmd = "rm -rf %s" % (vis_name)
         print cmd
         os.system(cmd)
         #also delete the miriad ms file
         vis_name = uvfits_name.split('.')[0] + '.ms'
         cmd = "rm -rf %s" % (vis_name)
         print cmd
         os.system(cmd)

         
  
def concat_vis(vis_list,outname):
   tmp_list = []
   for vis_index,vis in enumerate(vis_list): 
      tmp_name = "t%s.v" % vis_index
      tmp_list.append(tmp_name)
      cmd = "mv  %s %s" % (vis,tmp_name)
      print cmd
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
            print "wrote fits image %s" %  power_pattern_average_interp_fitsname
   
   
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
      print 'unknown sky_model'
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
      print "saved %s" % fig_name
      
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
      print "saved %s" % fig_name
      
      
   np.save(gsm_hpx_cube_smooth_filename,gsm_cube_smooth)
   print "saved %s" % (gsm_hpx_cube_smooth_filename)

#generate_smooth_hpx_cube(sky_model,freq_MHz_list,gsm_smooth_poly_order)
#sys.exit()

#main program
def simulate(lst_list,freq_MHz_list,pol_list,signal_type_list,sky_model='gsm'):
   #Do this stuff for each lst and each freq:     
   #model_sky_uvfits_list_X_lsts = []
   #model_global_signal_uvfits_list_X_lsts = []
   #model_noise_uvfits_list_X_lsts = []
   #model_sky_uvfits_list_Y_lsts = []
   #model_global_signal_uvfits_list_Y_lsts = []
   #model_noise_uvfits_list_Y_lsts = []
   model_vis_uvfits_list_X_lsts = []
   model_vis_uvfits_list_Y_lsts = []
   
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
         print date_time_string
         #for miriad time just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
         day_frac_plus_one = float(lst)/24. + 1
         miriad_uvgen_time_string = '00JAN%1.3f' % day_frac_plus_one
         print miriad_uvgen_time_string
         wavelength = 300./freq_MHz
         freq_GHz = freq_MHz/1000.
         
         gsm_hpx_fits_name = "%s_map_LST_%03d_%s_MHz_hpx.fits" % (sky_model,lst_deg,int(freq_MHz))
         reprojected_gsm_fitsname = "%s_map_LST_%03d_%s_MHz_reprojected.fits" % (sky_model,lst_deg,int(freq_MHz))
         reprojected_gsm_im_name = "%s_map_LST_%03d_%s_MHz_reprojected.im" % (sky_model,lst_deg,int(freq_MHz))
         
         global_signal_hpx_fits_name = "global_signal_map_LST_%03d_%s_MHz_hpx.fits" % (lst_deg,int(freq_MHz))
         reprojected_global_signal_fitsname = "global_signal_map_LST_%03d_%s_MHz_reprojected.fits" % (lst_deg,int(freq_MHz))
         reprojected_global_signal_im_name = "global_signal_map_LST_%03d_%s_MHz_reprojected.im" % (lst_deg,int(freq_MHz))
         
         #lna_impedance_aavs1_filename = "/data/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
         
         base_vis_name = "eda_model_LST_%03d_%s_MHz.vis" % (lst_deg,int(freq_MHz))
         #eda_model_uvfits_name = "eda_model_LST_%03d_%s_MHz.uvfits" % (lst_deg,int(freq_MHz))
         eda_model_no_source_image_name = "eda_no_source_LST_%03d_%s_MHz.map" % (lst_deg,int(freq_MHz))
         eda_model_no_source_beam_name = "eda_no_source_LST_%03d_%s_MHz.beam" % (lst_deg,int(freq_MHz))
         eda_model_no_source_image_fits_name = "eda_no_source_LST_%03d_%s_MHz.fits" % (lst_deg,int(freq_MHz))
         
         #for gain errors:
         one_jy_source_vis_name =  "one_jy_source_LST_%03d_%s_MHz.vis" % (lst_deg,int(freq_MHz))
         gain_solutions_name_amp = "gains_%03d_%s_MHz_amp.txt" % (lst_deg,int(freq_MHz))
         gain_solutions_name_phase = "gains_%03d_%s_MHz_phase.txt" % (lst_deg,int(freq_MHz))
         
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
         else:
            print 'unknown sky_model'
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
            cmd = "uvgen source=$MIRCAT/no.source ant='/data/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx  " % (miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,base_vis_name)
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
               cmd = "uvgen source=$MIRCAT/point.source ant='/data/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 out=%s stokes=xx gnoise=10 pnoise=90  " % (miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec,harange_string,one_jy_source_vis_name)
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

            gsm_map = gsm.generate(freq_MHz)

            hp.write_map(gsm_hpx_fits_name,gsm_map,coord='C')
         else:
            gsm_map = hp.read_map(gsm_hpx_fits_name)
         
         #plot?
         if plot_gsm_map_hpx:
            #plot
            plt.clf()
            map_title="GSM from MWA at %.0f:%.0f:%.0f and %s MHz" % (hour,minute,second,int(freq_MHz))
            #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=gsm_map,half_sky=False,title=map_title)
            #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
            fig_name="%s_map_%s_%sMHz.png" % (sky_model,date_time_string,int(freq_MHz))
            figmap = plt.gcf()
            figmap.savefig(fig_name,dpi=500)
            print "saved %s" % fig_name
         
         #Miriad doesn't seem to be able to import the hpx file
         #Try using reproject_from_healpix
         #output_projection is the header of the pt source made above
         if os.path.isfile(eda_model_no_source_image_fits_name) and os.access(eda_model_no_source_image_fits_name, os.R_OK):
            hdulist = pyfits.open(eda_model_no_source_image_fits_name)
         else:
            print "Either file %s is missing or is not readable" % eda_model_no_source_image_fits_name
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
         
         reprojected_gsm_map,footprint = reproject_from_healpix(gsm_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
         
         #write the map to fits
         pyfits.writeto(reprojected_gsm_fitsname,reprojected_gsm_map,clobber=True)
         pyfits.update(reprojected_gsm_fitsname,reprojected_gsm_map,header=no_source_header)
         print "wrote image %s" %  reprojected_gsm_fitsname
   
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
         print "scale map by %s to get to Jy/pix" % scale
         
         reprojected_gsm_im_Jy_per_pix_name =  "%s_map_%s_%s_MHz_reprojected_Jy_pix.im" % (sky_model,date_time_string,int(freq_MHz))
         
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
            map_title="Global signal from MWA at %.0f:%.0f:%.0f and %s MHz" % (hour,minute,second,int(freq_MHz))
            #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
            hp.orthview(map=s_21_hpx_map,half_sky=False,title=map_title)
            #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
            fig_name="global_signal_map_%s_%sMHz.png" % (date_time_string,int(freq_MHz))
            figmap = plt.gcf()
            figmap.savefig(fig_name)
            print "saved %s" % fig_name
              
         reprojected_global_signal_map,footprint = reproject_from_healpix(global_signal_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
         
         #write the map to fits
         pyfits.writeto(reprojected_global_signal_fitsname,reprojected_global_signal_map,clobber=True)
         pyfits.update(reprojected_global_signal_fitsname,reprojected_global_signal_map,header=no_source_header)
         print "wrote image %s" %  reprojected_global_signal_fitsname
         
         cmd = "rm -rf %s" % (reprojected_global_signal_im_name)
         print(cmd)
         os.system(cmd)     
         
         cmd = "fits in=%s out=%s op=xyin" % (reprojected_global_signal_fitsname,reprojected_global_signal_im_name)
         print(cmd)
         os.system(cmd)
         
         #uvmodel requires the model to be in Jy/pix
         #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
         
         reprojected_global_signal_im_Jy_per_pix_name =  "global_map_%s_%s_MHz_reproj_Jy_pix.im" % (date_time_string,int(freq_MHz))
         
         cmd = "rm -rf %s" % reprojected_global_signal_im_Jy_per_pix_name
         print(cmd)
         os.system(cmd)
               
         cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_global_signal_im_name,reprojected_global_signal_im_Jy_per_pix_name)
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
            model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg,pol,int(freq_MHz))

               
            if use_analytic_beam:
               if pol=='X':
                  beam_image_sin_projected_fitsname = "model_%s_MHz_%s.fits" % (int(freq_MHz),'xx')
               else:
                  beam_image_sin_projected_fitsname = "model_%s_MHz_%s.fits" % (int(freq_MHz),'yy')
            else:
                  beam_image_sin_projected_fitsname = "power_pattern_average_%s_%s_MHz_sin_regrid.fits" % (pol,int(freq_MHz))
    
            #power_pattern_average_interp_sin_im_name = 'power_pattern_average_%s_%s_MHz_interp_sin.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_im_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.im' % (pol,int(freq_MHz))
            #power_pattern_average_interp_sin_regrid_gsm_fits_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.fits' % (pol,int(freq_MHz))
                 
            beam_image_sin_projected_im_name = 'beam_image_sin_projected_%s_%s_MHz.im' % (pol,int(freq_MHz))
            beam_image_sin_projected_puthd_fits_name = 'beam_image_sin_projected_%s_%s_MHz_puthd.fits' % (pol,int(freq_MHz))
            beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%s_MHz_gsm_regrid.im' % (pol,int(freq_MHz))
            
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
         
            apparent_sky_im_name = "apparent_sky_LST_%03d_%s_pol_%s_MHz.im" % (lst_deg,pol,int(freq_MHz))
         
            cmd = "rm -rf %s " % apparent_sky_im_name
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
            print(cmd)
            os.system(cmd)
   
            #Repeat for global signal
            apparent_global_signal_im_name = "apparent_global_signal_LST_%03d_%s_pol_%s_MHz.im" % (lst_deg,pol,int(freq_MHz))
            apparent_global_signal_fits_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%s_MHz.fits" % (lst_deg,pol,int(freq_MHz))
            apparent_global_signal_im_name_cropped = "apparent_global_signal_LST_%03d_%s_pol_%s_MHz_cropped.im" % (lst_deg,pol,int(freq_MHz))
         
            cmd = "rm -rf %s %s %s" % (apparent_global_signal_im_name,apparent_global_signal_fits_name_cropped,apparent_global_signal_im_name_cropped)
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_global_signal_im_Jy_per_pix_name,apparent_global_signal_im_name)
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
                              
         
            #then put into the vis 
            
            if 'noise' in signal_type_list:
            
               out_vis_name =  "eda_model_noise_LST_%03d_%s_%s_MHz.vis" % (lst_deg,pol,int(freq_MHz))
               model_vis_name_base += '_noise'
               
               
               ##Put in the uvgen for the noise here (so we get different noise for X and Y pol (i.e. remove noise uvgen from above!!))
               #generate noise only visibilites
               #for systemp, just use sky temp 180 at 180, for JyperK, use Aeff from Table 2 of Wayth et al 2017 (EDA paper)
               systemp = T_180*(freq_MHz/180.0)**beta
               print 'systemp %s K' % systemp
               A_eff = Aeff_for_freq_MHz_list[freq_MHz_index]
               print 'A_eff %s K' % A_eff
               JperK = (2.0*k*10**26)/(eta*A_eff)
               print 'JperK %s ' % JperK
               SEFD = systemp * JperK
               print 'SEFD %s ' % SEFD
           
               #generate the noise-only uvfits file
               cmd = "rm -rf %s" % out_vis_name
               print(cmd)
               os.system(cmd)
            
               cmd = "uvgen source=$MIRCAT/no.source ant='/data/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant' baseunit=-3.33564 corr='1,1,0,1' time=%s freq=%.4f,0.0 radec='%2.3f,%s' harange=%s lat=-26.70331940 systemp=%s jyperk=%s out=%s stokes=xx  " % (miriad_uvgen_time_string,freq_GHz,float(lst),pointing_dec, harange_string, systemp, JperK, out_vis_name)
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
            
               model_vis_name_base += '_diffuse_%s' % sky_model
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
               
            if 'global' in signal_type_list:
            
               model_vis_name_base += '_global'
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
        
            if 'gain_errors' in signal_type_list:
               model_vis_name_base += '_gain_errors'
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
            
            cmd = "rm -rf %s" % out_vis_uvfits_name
            print(cmd)
            os.system(cmd)
            
            cmd = "fits in=%s out=%s op=uvout" % (out_vis_name,out_vis_uvfits_name)
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
            else:
               #model_sky_vis_list_Y.append(model_sky_vis)
               #model_global_signal_vis_list_Y.append(model_global_signal_vis)    
               #model_sky_uvfits_list_Y.append(model_sky_uvfits)
               #model_global_signal_uvfits_list_Y.append(model_global_signal_uvfits)      
               #model_noise_uvfits_list_Y.append(eda_model_noise_uvfits_name)  
               model_vis_uvfits_list_Y.append(out_vis_uvfits_name)
            #Export each individually to uvfits and then combine with pyuvdata (miriad concat doesn't update the header and so you can't properly export a concat'd vis set to uvfits (i think))
            
            #(Don't) Skip the concat step - it takes too long, we can just cycle through each individual uvfits file (1 for each freq) in the analysis steps    
            
            #delete all the intermediate images and vis that are no longer required
            cleanup_images_and_vis(lst,freq_MHz,pol)
      
      ###concat vis both pols for each freq
      for pol in pol_list:
         output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
         
         cmd = "rm -rf %s" % output_concat_vis_pyuvdata_name
         print(cmd)
         os.system(cmd)
         cmd = "rm -rf %s" % output_concat_uvfits_pyuvdata_name
         print(cmd)
         os.system(cmd)
         if pol=='X':
             concat_uvfits(model_vis_uvfits_list_X,output_concat_uvfits_pyuvdata_name)
         else:
             concat_uvfits(model_vis_uvfits_list_Y,output_concat_uvfits_pyuvdata_name)
      
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
         else:   
            #model_sky_uvfits_list_Y_lsts.append(output_concat_model_uvfits_name)
            #model_global_signal_uvfits_list_Y_lsts.append(global_output_concat_model_uvfits_name)      
            #model_noise_uvfits_list_Y_lsts.append(noise_output_concat_model_uvfits_name)
            model_vis_uvfits_list_Y_lsts.append(output_concat_uvfits_pyuvdata_name)
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
   #Data will likely come in a number of 'snapshot' uvfits files so might as well just digest these individually in extract_signal etc
   
   #concat the lsts together in 2-hour chunk
   for pol in pol_list:
      output_concat_vis_pyuvdata_name_lsts = "%s_concat_lsts.vis" % model_vis_name_base
      output_concat_uvfits_pyuvdata_name_lsts = "%s_concat_lsts.uvfits" % model_vis_name_base   
      cmd = "rm -rf %s" % output_concat_vis_pyuvdata_name_lsts
      print(cmd)
      os.system(cmd)
      cmd = "rm -rf %s" % output_concat_uvfits_pyuvdata_name_lsts
      print(cmd)
      os.system(cmd)
      if pol=='X':
          concat_uvfits(model_vis_uvfits_list_X_lsts,output_concat_uvfits_pyuvdata_name_lsts)
      else:
          concat_uvfits(model_vis_uvfits_list_Y_lsts,output_concat_uvfits_pyuvdata_name_lsts)
       
   #remove the intermediate uvfits and concat freq uvfits
   for lst in lst_list:
      lst_deg = (float(lst)/24.)*360.
      for pol in pol_list:
         for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
            model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg,pol,int(freq_MHz))
            if 'noise' in signal_type_list:
               model_vis_name_base += '_noise'
            if 'diffuse' in signal_type_list:
               model_vis_name_base += '_diffuse_%s' % sky_model
            if 'global' in signal_type_list:
               model_vis_name_base += '_global' 
            if 'gain_errors' in signal_type_list:
               model_vis_name_base += '_gain_errors'
            
            out_vis_name = model_vis_name_base + '.vis'
            out_vis_uvfits_name = model_vis_name_base + '.uvfits'
            
            cmd = "rm -rf %s %s" % (out_vis_name,out_vis_uvfits_name)
            print(cmd)
            os.system(cmd)
            
         output_concat_vis_pyuvdata_name = "%s_concat_freqs.vis" % model_vis_name_base
         output_concat_uvfits_pyuvdata_name = "%s_concat_freqs.uvfits" % model_vis_name_base
            
         cmd = "rm -rf %s %s" % (output_concat_vis_pyuvdata_name,output_concat_uvfits_pyuvdata_name)
         print(cmd)
         os.system(cmd)
            
            
            
def plot_signal(lst_list,freq_MHz_list,pol_list,signal_type_list=['diffuse','global','noise','gain_errors'],outbase_name='extracted_1_',sky_model='gsm'):
   freq_MHz_end = freq_MHz_list[-1]
   lst_end = lst_list[-1]
   lst_deg_end = (float(lst_end)/24.)*360.
   for pol_index,pol in enumerate(pol_list):
      model_vis_name_base = "eda_model_LST_%03d_%s_%s_MHz" % (lst_deg_end,pol,int(freq_MHz_end))
      if 'noise' in signal_type_list:
          model_vis_name_base += '_noise'
      if 'diffuse' in signal_type_list:
          model_vis_name_base += '_diffuse_%s' % sky_model
      if 'global' in signal_type_list:
          model_vis_name_base += '_global' 
      if 'gain_errors' in signal_type_list:
          model_vis_name_base += '_gain_errors'

      signal_array_short_baselines_filename = outbase_name + "%s_signal.npy" % (model_vis_name_base)
      signal_short_baselines = np.load(signal_array_short_baselines_filename)
      signal_array_short_baselines_Tb_filename = outbase_name + "%s_signal_Tb.npy" % (model_vis_name_base)
      signal_short_baselines_Tb = np.load(signal_array_short_baselines_Tb_filename)
      signal_array_all_baselines_filename = outbase_name + "%s_signal_all_baselines.npy" % (model_vis_name_base)
      signal_all_baselines = np.load(signal_array_all_baselines_filename)
      signal_array_all_baselines_Tb_filename = outbase_name + "%s_signal_all_baselines_Tb.npy" % (model_vis_name_base)
      signal_all_baselines_Tb = np.load(signal_array_all_baselines_Tb_filename)
      signal_array_all_baselines_filename_abs = outbase_name + "%s_signal_all_baselines_abs.npy" % (model_vis_name_base)
      signal_all_baselines_abs = np.load(signal_array_all_baselines_filename_abs)
      signal_array_all_baselines_filename_abs_Tb = outbase_name + "%s_signal_all_baselines_abs_Tb.npy" % (model_vis_name_base)
      signal_all_baselines_abs_Tb = np.load(signal_array_all_baselines_filename_abs_Tb)
      signal_array_short_baselines_weighted_filename = outbase_name + "%s_signal_weighted.npy" % (model_vis_name_base)
      signal_short_baselines_weighted = np.load(signal_array_short_baselines_weighted_filename)
      signal_array_short_baselines_weighted_Tb_filename = outbase_name + "%s_signal_weighted_Tb.npy" % (model_vis_name_base)
      signal_short_baselines_weighted_Tb = np.load(signal_array_short_baselines_weighted_Tb_filename)
      number_baselines_used_array_filename = outbase_name + "%s_number_baselines_used.npy" % (model_vis_name_base)
      number_baselines_used_array = np.load(number_baselines_used_array_filename)
      sum_of_weights_all_baselines_array_filename = outbase_name + "%s_sum_of_weights_all_baselines.npy" % (model_vis_name_base)
      sum_of_weights_all_baselines_array = np.load(sum_of_weights_all_baselines_array_filename)
      sum_of_weights_short_baselines_array_filename = outbase_name + "%s_sum_of_weights_short_baselines.npy" % (model_vis_name_base)
      sum_of_weights_short_baselines_array = np.load(sum_of_weights_short_baselines_array_filename)
      signal_short_baselines_log = np.log10(abs(signal_short_baselines[0,:]))
      signal_all_baselines_log = np.log10(abs(signal_all_baselines[0,:]))
      signal_all_baselines_abs_log = np.log10(signal_all_baselines_abs[0,:])
      signal_short_baselines_weighted_log = np.log10(abs(signal_short_baselines_weighted[0,:]))
      
      #all baselines Tb real
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines_Tb[0,:])
      map_title="all baselines real Tb vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Tb")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_all_baselines_Tb.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #short baselines real Tb
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines_Tb[0,:])
      map_title="short baselines real Tb vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Tb")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_short_baselines_Tb.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #short baselines real weighted Tb
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines_weighted_Tb[0,:])
      map_title="short baselines weighted real Tb vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Tb")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_short_baselines_weighted_Tb.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #all baselines abs Tb
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines_abs_Tb[0,:])
      map_title="all baselines abs Tb vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Tb")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_abs_vis_vs_freq_all_baselines_Tb.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name

      #short baselines real
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines[0,:])
      map_title="short baselines real vis vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_short_baselines.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #short baselines real log
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines_log)
      map_title="short baselines log abs real vis vs freq %s pol" % (pol)
      plt.ylabel("log sum abs real vis Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_log_abs_real_vis_vs_freq_short_baseline.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #short baselines real weighted
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines_weighted[0,:])
      map_title="short baselines weighted real vis vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_short_baselines_weighted.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
    
      #short baselines real weighted log
      plt.clf()
      plt.plot(freq_MHz_list,signal_short_baselines_weighted_log)
      map_title="short baselines weighted real vis vs freq %s pol" % (pol)
      plt.ylabel("log sum abs real vis (Jy)")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_log_real_vis_vs_freq_short_baselines_weighted.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name    

      #all baselines real
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines[0,:])
      map_title="all baselines real vis vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_real_vis_vs_freq_all_baselines.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #all baselines real log
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines_log)
      map_title="all baselines log abs real vis vs freq %s pol" % (pol)
      plt.ylabel("log sum abs real vis Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_log_abs_real_vis_vs_freq_all_baseline.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name

      #all baselines abs
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines_abs[0,:])
      map_title="all baselines abs vis vs freq %s pol" % (pol)
      plt.ylabel("sum vis real Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_abs_vis_vs_freq_all_baselines.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #all baselines abs log
      plt.clf()
      plt.plot(freq_MHz_list,signal_all_baselines_abs_log)
      map_title="all baselines log abs vis vs freq %s pol" % (pol)
      plt.ylabel("log sum abs vis Jy")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_log_abs_vis_vs_freq_all_baseline.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #number of baselines
      plt.clf()
      plt.plot(freq_MHz_list,number_baselines_used_array[0,:])
      map_title="Number of short baselines used"
      plt.ylabel("Number of baselines")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_number_of_baselines_used.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      
      #sum of weights all baselines
      plt.clf()
      plt.plot(freq_MHz_list,sum_of_weights_all_baselines_array[0,:])
      map_title="Sum of u,v=0 weights all baselines"
      plt.ylabel("Sum of weights")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_sum_of_weights_all_baselines.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name

      #sum of weights short baselines
      plt.clf()
      plt.plot(freq_MHz_list,sum_of_weights_short_baselines_array[0,:])
      map_title="Sum of u,v=0 weights short baselines"
      plt.ylabel("Sum of weights")
      plt.xlabel("freq (MHz)")
      plt.legend(loc=1)
      fig_name= "%s%s_sum_of_weights_short_baselines.png" % (outbase_name,model_vis_name_base)
      figmap = plt.gcf()
      figmap.savefig(fig_name)
      print "saved %s" % fig_name
      

#calculate the global 21cm signal:
s_21_array = plot_S21(nu_array=freq_MHz_list,C=C,A=A,delta_nu=delta_nu,nu_c=nu_c)

simulate(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model)
#sys.exit() 

compute_weights(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,sky_model=sky_model)
#sys.exit()

extract_signal_from_sims(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name='test_concat_',sky_model=sky_model)
#sys.exit()

plot_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name='test_concat_',sky_model=sky_model)

model_signal(lst_list=lst_list,freq_MHz_list=freq_MHz_list,pol_list=pol_list,signal_type_list=signal_type_list,outbase_name='test_concat_',poly_order=poly_order,sky_model=sky_model)



#list = ['eda_model_plus_global_LST_000_X_pol_50_MHz.uvfits','eda_model_noise_LST_000_X_50_MHz.uvfits','eda_model_plus_sky_LST_000_X_pol_50_MHz.uvfits']  
#list = ['eda_model_plus_global_LST_000_X_pol_50_MHz.uvfits','eda_model_plus_global_LST_001_X_pol_50_MHz.uvfits']
#outname = 'test_concat_LST_000_001_X_pol_50_MHz.uvfits'
#concat_uvfits(list,outname)
#sys.exit()

#for pol in pol_list:
#   for freq_MHz in freq_MHz_list:
#      recompute_ffts(pol,freq_MHz)
#      #plot_beam_and_weights(pol,freq_MHz)
#
#for lst in lst_list:
#   lst_deg = (float(lst)/24.)*360.
#   for pol in pol_list:
#      for freq_MHz in freq_MHz_list:
#         #cleanup_images_and_vis(lst,freq_MHz,pol)
#         compute_weights(lst_deg,pol,freq_MHz)

