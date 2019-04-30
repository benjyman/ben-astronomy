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
from pygsm import GSMObserver2016
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

#global signal params
#From Cath feasibility study (look at edges papers?)
#S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
C = 0.08/1000.
A = -150./1000.
delta_nu = 10.
nu_c = 70.

recompute_ffts = False

add_diffuse_model = True
add_global_signal = True

use_analytic_beam = True
generate_new_hpx = True
plot_gsm_map_hpx = False
generate_new_vis = True
plot_global_signal_map_hpx = True

#This is for Daniels FEKO model beams
generate_new_average_beam = False
#apply_normalisation = True
plot_all_beams = False
plot_average_beam = False

#specify date
#date_time_string = '2018_09_25_00_00_00'
#lst_list = ['00', '03', '06', '09', '12', '15', '18', '21']
lst_list = ['00']
#lst_list = ['12']
#pol_list = ['X','Y']
pol_list = ['X']
#freq_MHz_list = np.arange(50,200,1)
start_chan=50
n_chan=150
chan_step = 1
freq_MHz_list = np.arange(start_chan,start_chan+n_chan,chan_step)
#freq_MHz_list = [160.]
#n_chan = len(freq_MHz_list)
n_pol = len(pol_list)


n_ants = 256
template_imsize = 512
template_cell_size_asec = 180./template_imsize*60.*60.
#amount to zero pad beam image for FFT
#want this to be big so we have good resolution in Fourier space (needed for close to zeros spacings)
padding_factor = 50
uvdist_wavelength_cutoff = 10.0
zero_spacing_leakage_threshold = 0.1

def plot_S21(nu_array=None,C=0.,A=1.,delta_nu=20.,nu_c=78.):
   #Global signal
   #Use S_21 = C*nu + A*exp((nu - nu_c)**2/(2*delta_nu)**2)
   S_21 = C*nu_array + A*np.exp(-(nu_array - nu_c)**2/(2*(delta_nu)**2))
   plt.clf()
   map_title="S_21 vs freq"
   plt.plot(nu_array,S_21)
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

def plot_from_beam_image(beam_image_name):
   profile_plot_length = 500
   plot_basename = beam_image_name.split(".")[0]
   fft_savename = plot_basename + "_beam_fft2.npy"
   fft_savename_shift = plot_basename + "_beam_fft2_shift.npy"
   wavelength = 300./70.
   hdulist = fits.open(beam_image_name)
   image_data = hdulist[0].data
   image_header = hdulist[0].header
   plt.clf()
   map_title="beam"
   plt.imshow(image_data)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "_beam.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   
   #Fourier response:
   #zero pd to get better resolution in Fourier space:
   image_data_shape = image_data.shape
   padding = padding_factor*np.asarray(image_data_shape)
   
   if recompute_ffts:
      beam_fft2 = np.fft.fft2(image_data,s=padding)
      beam_fft2_shift = np.fft.fftshift(beam_fft2)
      #save the fft
      np.save(fft_savename,beam_fft2)
      np.save(fft_savename_shift,beam_fft2_shift)
   else:
      beam_fft2 = np.load(fft_savename)
      beam_fft2_shift = np.load(fft_savename_shift)
      
   
   beam_fft2_abs = np.abs(beam_fft2)
   beam_fft2_abs_shift = np.abs(beam_fft2_shift)
   beam_fft2_abs_shift_log = np.log10(beam_fft2_abs_shift)
   
   #FreqCompRows = np.fft.fftfreq(FFTData.shape[0],d=2)
   #FreqCompCols = np.fft.fftfreq(FFTData.shape[1],d=2)
   
   #no pad
   beam_image_length = image_data.shape[0]
   #image_centre = int(beam_image_length/2)
   
   image_centre_padded = int(beam_fft2_shift.shape[0]/2.)
   
   ##at phase centre angle is small sin(theta) = theta
   ##l_step = 1.0/(beam_image_length/2.)
   theta_step_rad = 1.0/(beam_image_length/2.)
   D_step_m = wavelength * theta_step_rad
   D_step_wavelengths = D_step_m/wavelength
   spatial_frequencies_cols = np.fft.fftfreq(beam_fft2.shape[1],d=D_step_wavelengths)
   #spatial_frequencies_rows = np.fft.fftfreq(beam_fft2.shape[0],d=D_step_wavelengths)
   spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
   
   #pad - this is wrong I think:
   #padded_beam_image_length = padding[0]
   #padded_theta_step_rad = 1.0/(padded_beam_image_length/2.)
   #D_step_m = wavelength * padded_theta_step_rad
   #D_step_wavelengths = D_step_m/wavelength
   #spatial_frequencies_cols = np.fft.fftfreq(beam_fft2.shape[1],d=D_step_wavelengths)
   
   
   plt.clf()
   map_title="beam fft2"
   plt.imshow(beam_fft2_abs_shift_log)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "_beam_fft2.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   
   #fits
   fits_name= plot_basename + "_beam_fft2.fits"
   pyfits.writeto(fits_name,beam_fft2_abs_shift_log,clobber=True)
   print "saved %s" % fits_name

   #profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile"
   #u or v?
   #beam_fft2_profile = beam_fft2_abs_shift_log[:,image_centre_padded]
   beam_fft2_profile = beam_fft2_abs_shift_log[image_centre_padded,:]
   
   #print beam_fft2_profile
   plt.plot(spatial_frequencies_cols_fftshift,beam_fft2_profile)
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "_beam_fft2_profile.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name

   #inner profile of beam response at v=something
   plt.clf()
   map_title="beam fft2 profile inner" 
   #print 'image centre %s' % image_centre
   #print 'profile plot length %s' % 100
   plt.plot(spatial_frequencies_cols_fftshift[image_centre_padded:image_centre_padded+profile_plot_length],beam_fft2_profile[image_centre_padded:image_centre_padded+profile_plot_length])
   #plt.ylabel("log abs(vis)")
   #plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "_beam_fft2_profile_inner.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name

def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return (index,array[index])
   
   
def compute_weights(uvfits_name):
   print uvfits_name
   hdulist = fits.open(uvfits_name)
   hdulist.info()
   info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
   #print info_string
   uvtable = hdulist[0].data
   uvtable_header = hdulist[0].header
   visibilities = uvtable['DATA']
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
   weights_array = np.full((n_pol,n_chan, n_vis),np.nan)
   #print weights_array.shape
         
   plot_basename = uvfits_name.split('.')[0]
   weights_array_filename = plot_basename + "_weights.npy"
   
   for pol_index,pol in enumerate(pol_list):
      for chan_index,chan in enumerate(range(0,n_chan,1)):
   
         #chan = 0
         print chan
         freq_MHz = 50.+chan*1.0
         wavelength = 300./freq_MHz
         
         if pol=='X':
            beam_image_name = "model_%s_MHz_xx.fits" % int(freq_MHz)
         else:
            beam_image_name = "model_%s_MHz_yy.fits" % int(freq_MHz)
         
         beam_plot_basename = beam_image_name.split('.')[0]
         #fft_savename = beam_plot_basename + "_beam_fft2.npy"
         #fft_savename_shift = beam_plot_basename + "_beam_fft2_shift.npy"
            
         beam_hdulist = fits.open(beam_image_name)
         beam_image_data = beam_hdulist[0].data
         beam_image_header = beam_hdulist[0].header
         
         #Fourier response:
         image_data_shape = beam_image_data.shape
         beam_image_length = image_data_shape[0]
         #image_centre = int(beam_image_length/2)
            
         #zero pad the beam image for later multiplication
         pad_width = int(padding_factor/2. * image_data_shape[0] - (image_data_shape[0]/2.)) 
         beam_image_data_padded = np.pad(beam_image_data,pad_width,'constant', constant_values=0)
         #print "beam_image_data_padded size %s" % beam_image_data_padded.shape[0]
         image_centre_padded = int(beam_image_data_padded.shape[0]/2.)
   
         padding = padding_factor*np.asarray(image_data_shape)
         
         fits_name= beam_plot_basename + "_beam_fft2_abs_shift_norm.fits"
         
         if recompute_ffts:
            beam_fft2 = np.fft.fft2(beam_image_data,s=padding)
            beam_fft2_abs = np.abs(beam_fft2)
            beam_fft2_abs_shift = np.abs(beam_fft2_shift)
            beam_fft2_abs_shift_norm = beam_fft2_abs_shift/np.max(beam_fft2_abs_shift)
            pyfits.writeto(fits_name,beam_fft2_abs_shift_norm,clobber=True)
            print "saved %s" % fits_name 
            #save the fft? - these are like a gig each - dont save!
            #np.save(fft_savename,beam_fft2)
            #np.save(fft_savename_shift,beam_fft2_shift)
            #print "saved %s" % fft_savename_shift
            #print "saved %s" % fft_savename
         else:
            #beam_fft2 = np.load(fft_savename)
            beam_fft2_abs_shift_norm_hdu_list = fits.open(fits_name)
            print "opened %s" % fits_name
            beam_fft2_abs_shift_norm = beam_fft2_abs_shift_norm_hdu_list[0].data
            
         beam_fft2_abs_shift_log = np.log10(beam_fft2_abs_shift_norm)
              

         #fits_name= beam_plot_basename + "_beam_fft2_abs_shift_log.fits"
         #pyfits.writeto(fits_name,beam_fft2_abs_shift_log,clobber=True)
         #print "saved %s" % fits_name        
                  
         beam_image_length = beam_image_data.shape[0]
         beam_fft2_abs_length = beam_fft2_abs_shift_norm.shape[0]
         ##at phase centre angle is small sin(theta) = theta
         ##l_step = 1.0/(beam_image_length/2.)
         theta_step_rad = 1.0/(beam_image_length/2.)
         D_step_m = wavelength * theta_step_rad
         D_step_wavelengths = D_step_m/wavelength
         spatial_frequencies_cols = np.fft.fftfreq(beam_fft2_abs_shift_norm.shape[1],d=D_step_wavelengths)
         spatial_frequencies_cols_fftshift = np.fft.fftshift(spatial_frequencies_cols)
         #spatial_frequencies_rows = np.fft.fftfreq(beam_fft2.shape[0],d=D_step_wavelengths)
         
   
         UU_wavelength_array = UU_m_array / wavelength
         VV_wavelength_array = VV_m_array / wavelength   
         
         for UU_wavelength_index, UU_wavelength in enumerate(UU_wavelength_array):
            
            #print "baseline U %s wavelengths" % UU_wavelength
            VV_wavelength = VV_wavelength_array[UU_wavelength_index]
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
               print "nearest_UU_index,nearest_UU_value %s,%s" % (nearest_UU_index, nearest_UU_value)
               nearest_VV_index, nearest_VV_value = find_nearest(spatial_frequencies_cols_fftshift,VV_wavelength)
               print "nearest_VV_index,nearest_VV_value %s,%s" % (nearest_VV_index, nearest_VV_value)
               
               
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
               UU_cheat_index = int((beam_fft2_abs_length/2.)-(nearest_UU_index-(beam_fft2_abs_length/2.)))
               VV_cheat_index = int((beam_fft2_abs_length/2.)-(nearest_VV_index-(beam_fft2_abs_length/2.)))
               
               print "UU_cheat_index is %s, VV_cheat_index is %s " % (UU_cheat_index,VV_cheat_index)
               
               zero_spacing_leakage_cheat = beam_fft2_abs_shift_norm[VV_cheat_index,UU_cheat_index]
               print "zero_spacing_leakage using cheat is %s " % zero_spacing_leakage_cheat
               
               
               if zero_spacing_leakage_cheat > zero_spacing_leakage_threshold:
                  print "zero_spacing_leakage is %s, using this baseline (%s of %s, chan %s)" % (zero_spacing_leakage_cheat,UU_wavelength_index,n_vis,int(freq_MHz))
                  #don't write a FITS file unless doing a subset i.e. making images for a talk
                  #fits_name= plot_basename + "_beam_footprint_fourier.fits"
                  #pyfits.writeto(fits_name,beam_footprint_mult_delta_image_space_inv_fft_abs_norm,clobber=True)
                  #print "saved %s" % fits_name            
                  
                  #now work out what to do with these weights! save in an array for now
                  weights_array[pol_index,chan_index,UU_wavelength_index] = zero_spacing_leakage_cheat
               #else:
                  #print "zero_spacing_leakage is %s, skipping this baseline (%s of %s, chan %s)" % (zero_spacing_leakage,UU_wavelength_index,n_vis,int(freq_MHz))
            #else:
               #print "uv_dist is %s, skipping this baseline (%s of %s, chan %s)" % (uvdist_wavelengths,UU_wavelength_index,n_vis,int(freq_MHz))
       
   np.save(weights_array_filename,weights_array)

#plot_from_beam_image("model_70_MHz_xx.fits")
#sys.exit()
compute_weights('global_eda_model_LST_00_X_pol_concat.uvfits')
sys.exit()

def extract_signal(uvfits_name):
   print uvfits_name
   hdulist = fits.open(uvfits_name)
   hdulist.info()
   info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
   print info_string
   uvtable = hdulist[0].data
   uvtable_header = hdulist[0].header
   visibilities = uvtable['DATA']
   print visibilities.shape
   n_vis = visibilities.shape[0]
   
   plot_basename = uvfits_name.split('.')[0]
   weights_array_filename = plot_basename + "_weights.npy" 
   weights_array = np.load(weights_array_filename)
   
   signal_array_unweighted = np.full((n_pol,n_chan,1),0,dtype=complex)
   signal_array_weighted = np.full((n_pol,n_chan,1),0,dtype=complex)
   signal_array_unweighted_filename = plot_basename + "_signal_unweighted.npy" 
   signal_array_weighted_filename = plot_basename + "_signal_weighted.npy" 
   
   for pol_index,pol in enumerate(pol_list):
      for chan_index,freq_MHz in enumerate(freq_MHz_list):
   
         wavelength = 300./freq_MHz
         #print freq_MHz
         #keep in mind pol_index may be wrong (need to rerun sims with pol=xx,yy only, and not sure what the last index is even for ... real imag weight?)
         sum_of_weights = 0
         no_baselines_weighted = 0
         for visibility_index,visibility_real in enumerate(visibilities[:,0,0,0,chan_index,pol_index,0]):
            #print visibility_real
            visibility_imag = visibilities[visibility_index,0,0,0,chan_index,pol_index,1]
            ##visibility_weight = visibilities[visibility_index,0,0,0,chan_index,pol_index,2]
            #print " %s %s i, weight:%s " % (visibility_real,visibility_imag,visibility_weight)
            complex_vis = np.complex(visibility_real,visibility_imag)
            #print "complex_vis %s" % complex_vis
            signal_array_unweighted[pol_index,chan_index,0] += complex_vis
            
            #now for weighted:
            uv_zero_weighting = weights_array[pol_index,chan_index,visibility_index]
            if not np.isnan(uv_zero_weighting):
               #print "uv_zero_weighting %s for baseline %s"  % (uv_zero_weighting,visibility_index)
               complex_vis_weighted = uv_zero_weighting*complex_vis
               #print "complex_vis_weighted %s" % complex_vis_weighted
               signal_array_weighted[pol_index,chan_index,0] += complex_vis_weighted
               sum_of_weights += uv_zero_weighting
               no_baselines_weighted += 1
               #print "sum_of_weights %s" % sum_of_weights
            
         #print abs(signal_array[pol_index,chan_index,0])
         
         #normalise by dividing by the sum of the weights_array
         print "sum_of_weights %s at freq %s MHz with %s baselines" % (sum_of_weights,freq_MHz,no_baselines_weighted)
         signal_array_weighted_norm = signal_array_weighted/sum_of_weights
        
   plt.clf()
   map_title="unweighted real vis vs freq x pol"
   plt.plot(freq_MHz_list,np.real(signal_array_unweighted[0,:,0]))
   plt.ylabel("unweighted sum vis real Jy")
   plt.xlabel("freq (MHz)")
   fig_name= plot_basename + "_real_vis_vs_freq_unweighted.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name       
   
   plt.clf()
   map_title="unweighted abs vis vs freq x pol"
   plt.plot(freq_MHz_list,np.abs(signal_array_unweighted[0,:,0]))
   plt.ylabel("unweighted sum vis real Jy")
   plt.xlabel("freq (MHz)")
   fig_name= plot_basename + "_abs_vis_vs_freq_unweighted.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name 
   
   plt.clf()
   map_title="weighted real vis vs freq x pol"
   plt.plot(freq_MHz_list,np.real(signal_array_weighted_norm[0,:,0]))
   plt.ylabel("weighted sum vis real Jy")
   plt.xlabel("freq (MHz)")
   fig_name= plot_basename + "_real_vis_vs_freq_weighted.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name       
   
   plt.clf()
   map_title="weighted abs vis vs freq x pol"
   plt.plot(freq_MHz_list,np.abs(signal_array_weighted_norm[0,:,0]))
   plt.ylabel("weighted sum vis real Jy")
   plt.xlabel("freq (MHz)")
   fig_name= plot_basename + "_abs_vis_vs_freq_weighted.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name 
   
#extract_signal('global_eda_model_LST_00_X_pol_concat.uvfits')
#sys.exit()
   
def plot_from_uvfits(uvfits_name):
   
   plot_basename = uvfits_name.split('.')[0]
   
   chan = 0
   freq_MHz1 = 50.+chan*1.0
   wavelength1 = 300./freq_MHz1
   #freq_MHz_array = np.arange(50,200,1)
   freq_MHz_list = np.arange(50,200,1)
   
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
   visibilities_at_freq_MHz1 = uvtable['DATA'][:,0,0,0,0,0,0]
   abs_vis = abs(visibilities_at_freq_MHz1)
   log_abs_vis = np.log10(abs_vis)
   UU_s = uvtable['UU']
   UU_m = UU_s * c
   
   UU_lambda = UU_m/wavelength1
   plt.clf()
   map_title="abs vis vs UU"
   plt.plot(UU_lambda,log_abs_vis,linestyle='none',marker='.')
   plt.ylabel("log abs(vis)")
   plt.xlabel("UU (lambda)")
   fig_name= plot_basename + "log_abs_vis_vs_uu.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name
   
   #for signal extraction only interested in short baselines (work out what this threshold should be based on U=0 leakage value later)
   threshold = 4.0
   threshold_indices = UU_m > threshold
   #short baselines only - as a function of frequency

   #for now, just take the max value of the visibilities at each freq
   max_vis_vs_freq_list = []
   for chan_index in np.arange(0,150):
      max_vis = np.max(uvtable['DATA'][:,0,0,0,chan_index,0,0])
      max_vis_vs_freq_list.append(max_vis)
   max_vis_vs_freq_array = np.asarray(max_vis_vs_freq_list)
   abs_max_vis_vs_freq_array = abs(max_vis_vs_freq_array)
   log_abs_max_vis_vs_freq_array = np.log10(max_vis_vs_freq_array)
   plt.clf()
   map_title="abs vis vs freq for short baselines"
   plt.plot(freq_MHz_list,log_abs_max_vis_vs_freq_array,linestyle='none',marker='.')
   plt.ylabel("log abs(max vis)")
   plt.xlabel("freq (MHz)")
   fig_name= plot_basename + "log_abs_max_vis_vs_freq.png"
   figmap = plt.gcf()
   figmap.savefig(fig_name)
   print "saved %s" % fig_name   

   #plt.clf()
   #map_title="log abs vis vs freq"
   #plt.plot(freq_MHz_array,abs_vis,linestyle='none',marker='.')
   #fig_name="abs_vis_vs_uu.png"
   #figmap = plt.gcf()
   #figmap.savefig(fig_name)
   #print "saved %s" % fig_name
   
   

#plot_from_uvfits('global_eda_model_LST_00_X_pol_concat.uvfits')
#plot_from_uvfits('eda_model_LST_00_X_pol_concat.uvfits')
#sys.exit()





def concat_uvfits(uvfits_list,outname):
   tmp_uv = UVData()
   first_filename = uvfits_list[0]
   fits.setval(first_filename, 'TELESCOP', value='MWA')
   fits.setval(first_filename, 'CDELT4', value=1.000E+06)
   #fits.setval(first_filename, 'DATE', value=2)
   #fits.setval(first_filename, '_DATE', value=0.5)
   #fits.setval(first_filename, 'INTTIM', value=2)
   print first_filename
   tmp_uv.read_uvfits(first_filename)
   #run_check=False,check_extra=False,run_check_acceptability=False
   for uvfits_name in uvfits_list[1:]: 
      #open the uvfits file and put in MWA as telescope
      fits.setval(uvfits_name, 'TELESCOP', value='MWA')
      fits.setval(uvfits_name, 'CDELT4', value=1.000E+06)
      #fits.setval(uvfits_name, 'DATE', value=2)
      #fits.setval(uvfits_name, '_DATE', value=0.5)
      #fits.setval(uvfits_name, 'INTTIM', value=2)
      new_uv = UVData()
      new_uv.read_uvfits(uvfits_name)
      print(new_uv.Nfreqs)
      tmp_uv = tmp_uv + new_uv
      print(tmp_uv.Nfreqs)
   tmp_uv.write_uvfits(outname)
   print "wrote %s" % outname
  
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

#calculate the global 21cm signal:
s_21_array = plot_S21(nu_array=freq_MHz_list,C=C,A=A,delta_nu=delta_nu,nu_c=nu_c)
                          
#Do this stuff for each lst and each freq:                   
for lst in lst_list:
   #model_sky_vis_list_X = []
   #model_global_signal_vis_list_X = []
   #model_sky_vis_list_Y = []
   #model_global_signal_vis_list_Y = []
   model_sky_uvfits_list_X = []
   model_global_signal_uvfits_list_X = []
   model_sky_uvfits_list_Y = []
   model_global_signal_uvfits_list_Y = []
   print("lst %s hrs" % lst)
   lst_deg = (float(lst)/24)*360.
   print("lst %s deg" % lst_deg)
   for freq_MHz_index,freq_MHz in enumerate(freq_MHz_list):
      print(freq_MHz)
      #for date just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
      date_time_string = '1900_01_01_%02d_00_00' % (float(lst))
      print date_time_string
      wavelength = 300./freq_MHz
      freq_GHz = freq_MHz/1000.
      
      gsm_hpx_fits_name = "gsm_map_LST_%s_%s_MHz_hpx.fits" % (lst,int(freq_MHz))
      reprojected_gsm_fitsname = "gsm_map_LST_%s_%s_MHz_reprojected.fits" % (lst,int(freq_MHz))
      reprojected_gsm_im_name = "gsm_map_LST_%s_%s_MHz_reprojected.im" % (lst,int(freq_MHz))
      
      global_signal_hpx_fits_name = "global_signal_map_LST_%s_%s_MHz_hpx.fits" % (lst,int(freq_MHz))
      reprojected_global_signal_fitsname = "global_signal_map_LST_%s_%s_MHz_reprojected.fits" % (lst,int(freq_MHz))
      reprojected_global_signal_im_name = "global_signal_map_LST_%s_%s_MHz_reprojected.im" % (lst,int(freq_MHz))
      
      #lna_impedance_aavs1_filename = "/data/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
      
      eda_model_vis_name = "eda_model_LST_%s_%s_MHz.vis" % (lst,int(freq_MHz))
      eda_model_no_source_image_name = "eda_no_source_LST_%s_%s_MHz.map" % (lst,int(freq_MHz))
      eda_model_no_source_beam_name = "eda_no_source_LST_%s_%s_MHz.beam" % (lst,int(freq_MHz))
      eda_model_no_source_image_fits_name = "eda_no_source_LST_%s_%s_MHz.fits" % (lst,int(freq_MHz))
      
      
      
      #if apply_normalisation:
      #   #get the lna impedance for this frequenccy
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

      ov = GSMObserver2016()
      ov.lon = longitude_degrees
      ov.lat = latitude_degrees
      ov.elev = elevation_m
      
      #and datetime of observation (eventually do this for many dates)
      year=int(date_time_string.split("_")[0])
      month=int(date_time_string.split("_")[1])
      day=int(date_time_string.split("_")[2])
      hour=int(date_time_string.split("_")[3])
      minute=int(date_time_string.split("_")[4])
      second=int(date_time_string.split("_")[5])
      
      ov.date = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
      
      
      #1. get 'blank' visibilities in miriad and make an image (1 Jy pt source) to use as a template for the output projection
      if generate_new_vis:
      
         cmd = "rm -rf %s" % eda_model_vis_name
         print(cmd)
         os.system(cmd)
      
         cmd = "uvgen source=$MIRCAT/no.source ant='/data/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant' baseunit=-3.33564 corr='1,1,0,1' time=19SEP26:00:00:00.0 freq=%.4f,0.0 radec='4,-27' harange=0,0.001,0.0005 lat=-26.70331940 out=%s stokes=xx,yy  " % (freq_GHz,eda_model_vis_name)
         print(cmd)
         os.system(cmd)
         
         cmd = "rm -rf %s %s %s " % (eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
         print(cmd)
         os.system(cmd)
         
         cmd = "invert vis=%s map=%s beam=%s imsize=%s cell=%s stokes=i options=mfs" % (eda_model_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,template_imsize,template_cell_size_asec)
         print(cmd)
         os.system(cmd)
      
         cmd = "fits in=%s out=%s op=xyout" % (eda_model_no_source_image_name,eda_model_no_source_image_fits_name)
         print(cmd)
         os.system(cmd)
         
      
             
      if generate_new_hpx:
         cmd = "rm -rf %s %s" % (gsm_hpx_fits_name,reprojected_gsm_im_name)
         print(cmd)
         os.system(cmd)
         gsm_map = ov.generate(freq_MHz)
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
         fig_name="gsm_map_%s_%sMHz.png" % (date_time_string,int(freq_MHz))
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
      
      reprojected_gsm_im_Jy_per_pix_name =  "gsm_map_%s_%s_MHz_reprojected_Jy_pix.im" % (date_time_string,int(freq_MHz))
      
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
      
      
      ########
      
      
      
      #Now for the beam (EDA)!
      if generate_new_average_beam:
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
            
       
      #do for all pols:
      for pol in pol_list:
         model_sky_vis = "eda_model_plus_sky_LST_%s_%s_pol_%s_MHz.vis" % (lst,pol,int(freq_MHz))
         model_global_signal_vis = "eda_model_plus_global_LST_%s_%s_pol_%s_MHz.vis" % (lst,pol,int(freq_MHz))
         model_sky_uvfits = "eda_model_plus_sky_LST_%s_%s_pol_%s_MHz.uvfits" % (lst,pol,int(freq_MHz))
         model_global_signal_uvfits = "eda_model_plus_global_LST_%s_%s_pol_%s_MHz.uvfits" % (lst,pol,int(freq_MHz))         
         if pol=='X':
            #model_sky_vis_list_X.append(model_sky_vis)
            #model_global_signal_vis_list_X.append(model_global_signal_vis)
            model_sky_uvfits_list_X.append(model_sky_uvfits)
            model_global_signal_uvfits_list_X.append(model_global_signal_uvfits)
         else:
            #model_sky_vis_list_Y.append(model_sky_vis)
            #model_global_signal_vis_list_Y.append(model_global_signal_vis)    
            model_sky_uvfits_list_Y.append(model_sky_uvfits)
            model_global_signal_uvfits_list_Y.append(model_global_signal_uvfits)      
         
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
         beam_image_sin_projected_regrid_gsm_im_name =  'beam_image_sin_projected_%s_%s_MHz_gsm_regrid.im' % (pol,int(freq_MHz))
         
         cmd = "rm -rf %s %s " % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name)
         print(cmd)
         os.system(cmd)
      
         cmd = "fits in=%s out=%s op=xyin" % (beam_image_sin_projected_fitsname,beam_image_sin_projected_im_name)
         print(cmd)
         os.system(cmd)
      
         #put in the correct ra in the header (ra = lst for zenith) 
         #puthd in="$beam/crval1" value=$lst_degs
         cmd = 'puthd in="%s/crval1" value=%0.2f' % (beam_image_sin_projected_im_name,lst_deg)
         print(cmd)
         os.system(cmd)         
      
         #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
         cmd = "regrid in=%s out=%s tin=%s tol=0" % (beam_image_sin_projected_im_name,beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_name)
         print(cmd)
         os.system(cmd)   
         
      
         #Sweet! finally have a gsm and a beam.
         #now multiply them to get the apparent sky and put that into uvmodel, 
      
         apparent_sky_im_name = "apparent_sky_LST_%s_%s_pol_%s_MHz.im" % (lst,pol,int(freq_MHz))
      
         cmd = "rm -rf %s " % apparent_sky_im_name
         print(cmd)
         os.system(cmd)
      
         cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
         print(cmd)
         os.system(cmd)

         #Repeat for global signal
         apparent_global_signal_im_name = "apparent_global_signal_LST_%s_%s_pol_%s_MHz.im" % (lst,pol,int(freq_MHz))
      
         cmd = "rm -rf %s " % apparent_global_signal_im_name
         print(cmd)
         os.system(cmd)
      
         cmd = "maths exp=%s*%s out=%s " % (beam_image_sin_projected_regrid_gsm_im_name,reprojected_global_signal_im_Jy_per_pix_name,apparent_global_signal_im_name)
         print(cmd)
         os.system(cmd)
         
               
      
         #then put into the vis 

         
         if add_diffuse_model == True:
            cmd = "rm -rf %s" % model_sky_vis
            print(cmd)
            os.system(cmd)
         
            cmd = "uvmodel vis=%s model=%s options=add out=%s" % (eda_model_vis_name,apparent_sky_im_name,model_sky_vis)
            print(cmd)
            os.system(cmd)

            cmd = "rm -rf %s" % model_sky_uvfits
            print(cmd)
            os.system(cmd)
         
            cmd = "fits in=%s out=%s op=uvout" % (model_sky_vis,model_sky_uvfits)
            print(cmd)
            os.system(cmd)            
            
        
         if add_global_signal:
            cmd = "rm -rf %s" % model_global_signal_vis
            print(cmd)
            os.system(cmd)
            
            cmd = "uvmodel vis=%s model=%s options=add out=%s" % (eda_model_vis_name,apparent_global_signal_im_name,model_global_signal_vis)
            print(cmd)
            os.system(cmd)
            
            cmd = "rm -rf %s" % model_global_signal_uvfits
            print(cmd)
            os.system(cmd)
         
            cmd = "fits in=%s out=%s op=uvout" % (model_global_signal_vis,model_global_signal_uvfits)
            print(cmd)
            os.system(cmd)   

            #delete the vis!
            
         
            #Export each individually to uvfits and then combine with pyuvdata (miriad concat doesn't update the header and so you can't properly export a concat'd vis set to uvfits (i think))
         
         
         
         ##do this out of freq loop
         ##then concat the model viss for each lst into one with multiple freqs
         ##if freq_MHz_index == 0:  
         ##  new_tmp_concat_vis = "tmp_eda_model_LST_%s_concat_%s_pol_%s.vis" % (lst,pol,freq_MHz_index)
         ##   cmd = "rm -rf %s" % new_tmp_concat_vis
         ##   print(cmd)
         ##   os.system(cmd)
         ##   cmd = "uvaver vis=%s out=%s" % (model_sky_vis,new_tmp_concat_vis)
         ##   print(cmd)
         ##   os.system(cmd)
         ##else:
         ##   prev_tmp_concat_vis = "tmp_eda_model_LST_%s_concat_%s_pol_%s.vis" % (lst,pol,freq_MHz_index-1)
         ##   new_tmp_concat_vis = "tmp_eda_model_LST_%s_concat_%s_pol_%s.vis" % (lst,pol,freq_MHz_index)
         ##   cmd = "rm -rf %s" % new_tmp_concat_vis
         ##   print(cmd)
         ##   os.system(cmd)
         ##   cmd = "uvaver vis=%s,%s out=%s" % (prev_tmp_concat_vis,eda_model_vis_name,new_tmp_concat_vis)
         ##   print(cmd)
         ##   os.system(cmd)
         ##   #remove the previous tmp_concat_vis
         ##   cmd = "rm -rf %s" % prev_tmp_concat_vis
         ##   print(cmd)
         ##   os.system(cmd)
         
            
   ##rename the concat vis both pols
   for pol in pol_list:
      #gsm
      output_concat_model_pyuvdata_name = "eda_model_LST_%s_%s_pol_concat.vis" % (lst,pol)
      output_concat_model_uvfits_name = "eda_model_LST_%s_%s_pol_concat.uvfits" % (lst,pol)
      cmd = "rm -rf %s" % output_concat_model_uvfits_name
      print(cmd)
      os.system(cmd)
      cmd = "rm -rf %s" % output_concat_model_pyuvdata_name
      print(cmd)
      os.system(cmd)
      if pol=='X':
          concat_uvfits(model_sky_uvfits_list_X,output_concat_model_uvfits_name)
      else:
          concat_uvfits(model_sky_uvfits_list_Y,output_concat_model_uvfits_name)
   #   model_sky_vis_list_string = ','.join(model_sky_vis_list)
   #   cmd = "uvaver vis=%s out=%s" % (model_sky_vis_list_string,output_concat_model_vis_name)
   #   print(cmd)
   #   os.system(cmd)

      #global signal
      global_output_concat_model_vis_name = "global_eda_model_LST_%s_%s_pol_concat.vis" % (lst,pol)
      global_output_concat_model_uvfits_name = "global_eda_model_LST_%s_%s_pol_concat.uvfits" % (lst,pol)
      cmd = "rm -rf %s" % global_output_concat_model_uvfits_name
      print(cmd)
      os.system(cmd)
      cmd = "rm -rf %s" % global_output_concat_model_vis_name
      print(cmd)
      os.system(cmd)
      if pol=='X':
          concat_uvfits(model_global_signal_uvfits_list_X,global_output_concat_model_uvfits_name)
      else:
          concat_uvfits(model_global_signal_uvfits_list_Y,global_output_concat_model_uvfits_name)
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



       
       
       
    
 
 
 

