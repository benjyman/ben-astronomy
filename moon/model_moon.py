#python code to model moon spectrum, remove rfi component ets.
import matplotlib
matplotlib.use('Agg')
import pyfits
import pylab as py
import numpy as np
from scipy import signal
from scipy import fftpack
import statsmodels.api as sm
#from sklearn import linear_model
import matplotlib.pyplot as plt
import os
import os.path
import math
from scipy import optimize
from pygsm import GlobalSkyModel2016
import healpy as hp


#Choose one or the other: LFSM or GSM (GSM2016)
use_model_gsm=True
use_model_lfsm=False

new_model_galaxy=False
generate_new_beam_maps=False
#set if you just want to plot already saved data
plot_only=True
#RFI broadening or not:
do_RFI_broadening=True


if use_model_gsm:
   NSIDE=1024 #nside for healpix maps (pygsm default is 512 for pygsm and 1024 for pygsm2016 (checked))
   beam_map_filename_base="/data/moon/primary_beam/gaussian_beam_gsm"
if use_model_lfsm:
   NSIDE=128
   beam_map_filename_base="/data/moon/primary_beam/gaussian_beam_lfsm"

#downloaded pixel transformation fits file from (https://lambda.gsfc.nasa.gov/toolbox/tb_pixelcoords.cfm)
#for pygsm
#pixel_coordinate_file="/data/code/healpix/pixel_coords_map_ring_celestial_res9.fits"
#for pygsm2016:
pixel_coordinate_file_gsm2016="/data/code/healpix/pixel_coords_map_ring_celestial_res10.fits"
pixel_coordinate_file_lfsm="/data/code/healpix/pixel_coords_map_ring_celestial_res7.fits"

#beamformer delays for obsid chan 69:1127321744 :  24,18,12,6,22,16,10,4,20,14,8,2,18,12,6,0

#range of RA DEC (in degrees) to compute the mean galactic signal for (i.e. the patch of sky that the moon is in and that the cropped images cover
#CenA: ra_min=13*15., ra_max=14*15., dec_min=-53., dec_max=-33.
#paper 1 moon pos:  22:52:17.79, dec =  -4:51:28.7,
#need the pointing position! RA,DEC: 333.00, -7.70
ra_pointing=333.00
dec_pointing=-7.70
#only used for old (wrong) way:
ra_min=ra_pointing-20.0
ra_max=ra_pointing+20.0
dec_min=dec_pointing-20.0
dec_max=dec_pointing+20.0
#Need to change the way I do this: the mean (global) signal is the 'beam-weighted average' of the entire sky so need to generate an 
#average beam map (including all obsids) at each frequency

if use_model_lfsm:
   average_temp_filename="average_lfsm_temp_for_moon_region.npy"
   gsm_fractional_error=0.15
if use_model_gsm:
   average_temp_filename="average_gsm_temp_for_moon_region.npy"
   gsm_fractional_error=0.15

#set if using already cropped images
use_cropped_images=False

#Set if using the small images. 'small images are 2048 x 20148' 
#But still need to crop as especially at the higher freqs, the primary beam is a lot narrower and PB errors affect the difference images.
small_images=True
ionpeeled_images=False

#Portion of difference image used to measure rms (diff is from cropped image - middle quarter)
#[1100:1600,1100:1600]
rms_start_x=0
rms_end_x=450
rms_start_y=0
rms_end_y=450

plot_images=False
#speed of light
c=299792458.0
#Boltzmann's constant
k=1.380648*10**(-23)
moon_radius_deg=0.2774
#Solid angle subtended by the Moon in steradians
#Omega=6.67*10.0**(-5)
#For a 33.29 arcmin moon e.g.: 20150926: http://aa.usno.navy.mil/imagery/disk?body=moon&year=2015&month=9&day=26&hour=16&minute=00
Omega=7.365*10.0**(-5)
T_moon=230 #K
#Moon RA 22:49:02 DEC -5:34:25 (2015-09-26) (RA 342.45 DEC -5.573)

T_refl_gal_150=25.26 # temp of reflected galactic emisssion at 150 from reflection modeling (GSM2016) (takes into account albedo of 0.07)
refl_gal_alpha=-2.50
#RFI model mask radius in arcmin - limits the extent of the gausssian reflection model
#rfi_model_mask_size=2.0
rfi_model_mask_size=(3.75/2)

#filename for big numpy arrays
big_Smoon_average_stddev_spectrum_npy_filename="big_Smoon_average_stddev_spectrum.npy"
big_Srfi_average_stddev_spectrum_npy_filename="big_Srfi_average_stddev_spectrum.npy"

big_Smoon_predicted_npy_filename="big_Smoon_predicted.npy"
big_Ssky_predicted_npy_filename="big_Ssky_predicted.npy"


#functions
def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    

#these aren't real beam maps, just a temporary solution until marchin builds his code to get the FEKO sim maps in healpix
#for now they are just Gaussians in hpx
def generate_beam_maps(frequency_array):
   if use_model_lfsm:
         pixel_coordinate_file=pixel_coordinate_file_lfsm
   if use_model_gsm:
         pixel_coordinate_file=pixel_coordinate_file_gsm2016
   #read in the coordinate transform file to convert pixel to ra dec (celestial)
   pixel_coordinate_transform_map_ra=hp.read_map(pixel_coordinate_file,field=0)
   pixel_coordinate_transform_map_dec=hp.read_map(pixel_coordinate_file,field=1)
   for freq_index,freq_MHz in enumerate(frequency_array):
      beam_map_filename="%s_%.0fMHz.npy" % (beam_map_filename_base,freq_MHz)
      beam_map_figname="%s_%.0fMHz.png" % (beam_map_filename_base,freq_MHz)
      print "generating beam map %s" % beam_map_filename
      beam_map=np.empty(hp.nside2npix(NSIDE))
      #print "beam_map.shape %s" % beam_map.shape
      fwhm=25.0/(freq_MHz/150)
      #print lfsm_map.shape
      for hpx_pixel,value in enumerate(beam_map):
         pix_ra_deg=pixel_coordinate_transform_map_ra[hpx_pixel]
         pix_dec_deg=pixel_coordinate_transform_map_dec[hpx_pixel]  
         #print "pix_ra_deg-ra_pointing %s" % str(pix_ra_deg-ra_pointing)
         #print "pix_dec_deg-dec_pointing %s" % str(pix_dec_deg-dec_pointing)
         
         beam_value=np.exp(-4*np.log(2) * ((pix_ra_deg-ra_pointing)**2 + (pix_dec_deg-dec_pointing)**2) / fwhm**2)
         beam_map[hpx_pixel]=beam_value

      #save the  average temp array
      np.save(beam_map_filename,beam_map)

      ##save the beam map as png           
      plt.clf()
      map_title="Gaussian beam %.0fMHz" % freq_MHz
      #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
      hp.mollview(map=beam_map,coord='G',xsize=400,title=map_title, max=1.0, min=0.001)
      figmap = plt.gcf()
      figmap.savefig(beam_map_figname,dpi=100) 
      plt.close()

#model=either lfsm or gsm
def model_galaxy(frequency_array, model_data):
   print "Doing Galaxy modeling with %s" % model_data
   average_temp_array=np.empty(len(frequency_array))
   if (model_data=='gsm'):
      gsm = GlobalSkyModel2016(freq_unit='MHz')
   for freq_index,freq_MHz in enumerate(frequency_array):
      beam_map_filename="%s_%.0fMHz.npy" % (beam_map_filename_base,freq_MHz)
      beam_map=np.load(beam_map_filename)
      if (model_data=='lfsm'):
         #Just read in maps already generated using realizelfsm.py
         lfsm_data_filename="/data/moon/low_freq_sky_model_lwa/output-%.0f.dat" % freq_MHz
         galaxy_map=np.loadtxt(lfsm_data_filename)      
      if (model_data=='gsm'):
         galaxy_map=gsm.generate(freq_MHz)
      
      mean_galaxy=np.mean(galaxy_map)
      print "mean_galaxy %s" % mean_galaxy
      #calculate the beam-weighted mean of the lfsm
      beam_weighted_mean_galaxy=np.sum(galaxy_map*beam_map)/np.sum(beam_map)
      print "beam_weighted_mean_galaxy %s freq %.0f MHz" % (beam_weighted_mean_galaxy,freq_MHz)
      average_temp_array[freq_index]=beam_weighted_mean_galaxy

   #print average_temp_array
   #save the  average temp array
   np.save(average_temp_filename,average_temp_array)

#no longer used:
def model_gsm(frequency_array):
   average_temp_array=np.empty(len(frequency_array))
   pixel_coordinate_file=pixel_coordinate_file_gsm2016
   gsm = GlobalSkyModel2016(freq_unit='MHz')
   #ream in the coordinate transform file to convert pixel to ra dec (celestial)
   pixel_coordinate_transform_map_ra=hp.read_map(pixel_coordinate_file,field=0)
   pixel_coordinate_transform_map_dec=hp.read_map(pixel_coordinate_file,field=1)
   #Do for each freq
   for freq_index,freq_MHz in enumerate(frequency_array):
      gsm_map=gsm.generate(freq_MHz)
      for hpx_pixel,gsm_value in enumerate(gsm_map):
         pix_ra_deg=pixel_coordinate_transform_map_ra[hpx_pixel]
         pix_dec_deg=pixel_coordinate_transform_map_dec[hpx_pixel]  
         if not (ra_min < pix_ra_deg < ra_max and dec_min < pix_dec_deg < dec_max):
             #print pix_ra_deg
             #print pix_dec_deg
             gsm_map[hpx_pixel]=np.nan
      #Calculate the mean value
      average_temp=np.nanmean(gsm_map)
      average_temp_array[freq_index]=average_temp

   #print average_temp_array
   #save the  average temp array
   np.save(average_temp_filename,average_temp_array)

#To get a spectrum across the whole MWA observing band, we need to repeat this for all 5 bands observed

#make one big freaking array to cover each 1.28 MHz subband
#band 1 - centre 69: 57:80
#band 2 - centre 93: 81:104 
#band 3 - centre 121: 109:132 (ORBCOMM skipped 4 chans 105:108)
#band 4 - centre 145: 133:156
#band 5 - centre 169: 157:180
number_coarse_chans=int(5*24+4)
big_chan_array=np.arange(number_coarse_chans)+57
big_freq_array=big_chan_array*1.28

 
if generate_new_beam_maps:
    generate_beam_maps(big_freq_array)

#initialise the other big arrays
big_Smoon_average_stddev_spectrum=np.empty([number_coarse_chans,2])
big_Srfi_average_stddev_spectrum=np.zeros([number_coarse_chans,2])
big_Ssky_predicted_values=np.zeros(number_coarse_chans)
big_Smoon_predicted_values=np.zeros(number_coarse_chans)
big_Tb_value_error=np.zeros([number_coarse_chans,2])


band_centre_chans=[69,93,121,145,169]
#band_centre_chans=[169]

if (not plot_only):

 for centre_chan in band_centre_chans:
   
   #set rms_thresholds for difference images
   rms_threshold=1.0
   
   #read the info files to find out how many observations there are and what the on and off-moon obsids are:
   on_moon_filename="20150926_moon_%s.txt" % (str(centre_chan))
   off_moon_filename="20150929_off_moon1_%s.txt" % (str(centre_chan))
   
   print on_moon_filename
   
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]
   
   with open(on_moon_filename,'r') as on_moon_file:
      on_moon_lines = on_moon_file.readlines()
   for line in on_moon_lines:
      on_moon_obsid_list.append(line.strip())
   
   with open(off_moon_filename,'r') as off_moon_file:
      off_moon_lines = off_moon_file.readlines()
   for line in off_moon_lines:
      off_moon_obsid_list.append(line.strip())
   
   n_obs=len(on_moon_obsid_list)
   n_chans=24 
   
   #multiple linear regresion (see http://stackoverflow.com/questions/11479064/multiple-linear-regression-in-python/14971531#14971531)
   #using statsmodels
       
       
   #read in moon images
   #for each obsid
   #for each channel
   if (use_cropped_images): 
      xstart_moon,xend_moon,ystart_moon,yend_moon=0,520,0,520
      xstart_psf,xend_psf,ystart_psf,yend_psf=0,520,0,520
   elif (small_images):
      #small images are 2048 x 20148
      #xstart_moon,xend_moon,ystart_moon,yend_moon=0,2048,0,2048
      #xstart_psf,xend_psf,ystart_psf,yend_psf=0,2048,0,2048    
      #Crop to only use the inner quarter of the images
      xstart_moon,xend_moon,ystart_moon,yend_moon=512,1536,512,1536
      xstart_psf,xend_psf,ystart_psf,yend_psf=512,1536,512,1536 
   else:
      #xstart,xend,ystart,yend=2000,3120,2000,3120
      #crop images to make FFTs quicker
      #69, 93 and 121 are all 5120 x 5120 pix
      if (centre_chan <= 121): 
         xstart_moon,xend_moon,ystart_moon,yend_moon=2300,2820,2300,2820
         xstart_psf,xend_psf,ystart_psf,yend_psf=2299,2819,2301,2821
      #145,169 are 3840 x 3840
      else:
         xstart_moon,xend_moon,ystart_moon,yend_moon=1660,2180,1660,2180
         xstart_psf,xend_psf,ystart_psf,yend_psf=1659,2179,1659,2179
   
   #initialise arrays of length n_chans to store Smoon and Srfi values
   Smoon_spectrum_values=np.empty([n_chans,n_obs])
   Smoon_spectrum_values[:]=np.nan 
   Srfi_spectrum_values=np.empty([n_chans,n_obs])
   Srfi_spectrum_values[:]=np.nan
   Smoon_average_stddev_spectrum=np.zeros([n_chans,2])
   Srfi_average_stddev_spectrum=np.zeros([n_chans,2])
   
   
   #need to do all this stuff for each frequency channel (and MF) and each obsid
   channel_list=range(n_chans)
   channel_list.append('MFS')
   for chan in channel_list:
      for obsid_index in range(n_obs):
         #sort out filenames 
         on_moon_obsid=on_moon_obsid_list[obsid_index]
         off_moon_obsid=off_moon_obsid_list[obsid_index] 
         if (chan != 'MFS'):
            chan_string='%.04d' % chan
         else:
            chan_string='MFS'
         if (use_cropped_images): 
            moon_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            off_moon_fitsname="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I_cropped.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            on_moon_psf_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            moon_difference_fitsname="images/difference_%s_%s_cotter_20150929_moon_%s_peeled-%s-I_cropped.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
            
            #moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s_dirty-I_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            #off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s_dirty-I_cropped.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            #psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s-psf_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            #moon_difference_fitsname="/data/moon/2017/difference_%s_%s_on_off_moon_%s-%s-I_cropped.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
         else:
            #moon_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001_dirty_applied-I.fits"
            #moon_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            #off_moon_fitsname="1127572080_cotter_20150929_moon_93_track_off_moon_paired_1127313592_peeled-0001_dirty_applied-I.fits"
            #off_moon_fitsname="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            #psf_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001-psf.fits"
            #psf_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            #difference image filename
            #moon_difference_fitsname="images/difference_%s_%s_cotter_20150929_moon_%s_peeled-%s-I.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
            if (ionpeeled_images):
               moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               on_moon_psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               off_moon_psf_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s-psf.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               moon_difference_fitsname="/data/moon/2017/difference_%s_%s_on_off_moon_%s_peeled-%s-I.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
               psf_difference_fitsname="/data/moon/2017/difference_%s_%s_psf_%s_peeled-%s-psf.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string) 
            else:
               moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s_dirty-I.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s_dirty-I.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               on_moon_psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s-psf.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               off_moon_psf_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s-psf.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               moon_difference_fitsname="/data/moon/2017/difference_%s_%s_on_off_moon_%s-%s-I.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
               psf_difference_fitsname="/data/moon/2017/difference_%s_%s_psf_%s-%s-psf.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)

         print "############################"
         print moon_fitsname  
         print off_moon_fitsname
         if os.path.isfile(moon_fitsname) and os.access(moon_fitsname, os.R_OK):
            moon_hdulist = pyfits.open(moon_fitsname)
         else:
            print "Either file %s is missing or is not readable" % moon_fitsname
            continue        
         moon_data=moon_hdulist[0].data[0,0,:,:]
         moon_data=np.nan_to_num(moon_data)
         moon_header=moon_hdulist[0].header
         moon_zoom=moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]
         pix_size_deg = np.abs(float(moon_header['cdelt1']))
         moon_radius_pix = np.round(moon_radius_deg/pix_size_deg)
        
         max_moon_value=np.max(moon_zoom)
         min_moon_value=np.min(moon_zoom)  
         #print max value of moon image in Jy/beam to check against rfi + moon total flux
         #print "Max moon is %s Jy/beam" % (np.max(moon_data))
         print "Max moon is %s Jy/beam" % max_moon_value
         print "Min moon is %s Jy/beam" % min_moon_value
   
         if (max_moon_value==min_moon_value):
            print "Moon max and min are the same - something is wrong with the image, discarding."
            continue
         #Need to have the images in Jy per pixel for the equations to make sense
         #know the pix area in degrees^2 = pix_size_deg x pix_size_deg
         pix_area_deg_sq = pix_size_deg * pix_size_deg
         #beam area (gaussian restoring beam)
         bmaj_deg=np.abs(float(moon_header['bmaj'])) 
         bmin_deg=np.abs(float(moon_header['bmin']))   
         #beamarea for 2d gussian 2 pi major minor / (8 ln 2) = 1.133 maj min ?
         beam_area_deg_sq=1.133 * bmaj_deg * bmin_deg
         n_pixels_per_beam=beam_area_deg_sq/pix_area_deg_sq
         #to convert from Jy/beam to Jy per pis, divide by n_pixels_per_beam
         
         #read in off moon images
         if os.path.isfile(off_moon_fitsname) and os.access(off_moon_fitsname, os.R_OK):
             off_moon_hdulist = pyfits.open(off_moon_fitsname)
         else:
             print "Either file %s is missing or is not readable" % off_moon_fitsname
             continue
         off_moon_hdulist = pyfits.open(off_moon_fitsname)
         off_moon_data=off_moon_hdulist[0].data[0,0,:,:]
         off_moon_data=np.nan_to_num(off_moon_data)
         off_moon_header=off_moon_hdulist[0].header
         off_moon_zoom=off_moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]
         
         max_off_moon_value=np.max(off_moon_zoom)
         min_off_moon_value=np.min(off_moon_zoom) 
         
         print "Max off_moon is %s Jy/beam" % max_off_moon_value
         print "Min off_moon is %s Jy/beam" % min_off_moon_value
   
         if (max_off_moon_value==min_off_moon_value):
            print "Off_moon max and min are the same - something is wrong with the image, discarding."
            continue
         
         #difference the images
         moon_minus_sky=moon_zoom-off_moon_zoom
         #convert to Jy per pix
         moon_minus_sky_jyppix=moon_minus_sky/n_pixels_per_beam

         #print out the rms of the difference image ( corner only )
         difference_rms=np.sqrt(np.mean(np.square(moon_minus_sky[rms_start_x:rms_end_x,rms_start_y:rms_end_y])))
         
         #difference rms in Jy/pixel:
         difference_rms_jy_per_pixel=np.sqrt(np.mean(np.square(moon_minus_sky_jyppix[rms_start_x:rms_end_x,rms_start_y:rms_end_y])))
         print "rms of difference map %s is %s in Jy/pixel." % (moon_difference_fitsname,difference_rms_jy_per_pixel)
         
         #If the rms is too high then just put a nan in the arrays and move on to next obsid
         if (difference_rms>rms_threshold or math.isnan(difference_rms)==True):
            print "rms of difference map %s is %s in Jy/beam. Discarding difference image." % (moon_difference_fitsname,difference_rms)
            #place values in the arrays
            if (chan != 'MFS'):
               Smoon_spectrum_values[chan,obsid_index]=np.nan
               Srfi_spectrum_values[chan,obsid_index]=np.nan
            continue  
         else:
            print "rms of difference map %s is %s in Jy/beam. Writing out difference image." % (moon_difference_fitsname,difference_rms)
            #write out the difference image
            pyfits.writeto(moon_difference_fitsname,moon_minus_sky,clobber=True)
            pyfits.update(moon_difference_fitsname,moon_minus_sky,header=moon_header)
            print "wrote difference image %s" %  moon_difference_fitsname

         #read in psf images
         if os.path.isfile(on_moon_psf_fitsname) and os.access(on_moon_psf_fitsname, os.R_OK):         
            psf_hdulist = pyfits.open(on_moon_psf_fitsname)
         else:
            print "Either file %s is missing or is not readable" % on_moon_psf_fitsname
            continue
         psf_data=psf_hdulist[0].data[0,0,:,:]
         psf_data=np.nan_to_num(psf_data)
         psf_header=psf_hdulist[0].header
         psf_zoom=psf_data[xstart_psf:xend_psf,ystart_psf:yend_psf]
         psf_zoom=np.require(psf_zoom, dtype=np.float32)
         #convert to Jy per pix
         psf_zoom_jyppix=psf_zoom/n_pixels_per_beam
 
         #if os.path.isfile(off_moon_psf_fitsname) and os.access(off_moon_psf_fitsname, os.R_OK):         
         #   off_moon_psf_hdulist = pyfits.open(off_moon_psf_fitsname)
         #else:
         #   print "Either file %s is missing or is not readable" % off_moon_psf_fitsname
         #   continue
         #off_moon_psf_data=off_moon_psf_hdulist[0].data[0,0,:,:]
         #off_moon_psf_header=off_moon_psf_hdulist[0].header
         #off_moon_psf_zoom=off_moon_psf_data[xstart_psf:xend_psf,ystart_psf:yend_psf]
         #off_moon_psf_zoom=np.require(off_moon_psf_data, dtype=np.float32)
         #convert to Jy per pix
         #off_moon_psf_zoom_jyppix=off_moon_psf_zoom/n_pixels_per_beam
         #
         #make the psf difference image and save
         ##difference the images
         #psf_difference=psf_zoom-off_moon_psf_zoom
   
         #write out the difference image
         #pyfits.writeto(psf_difference_fitsname,psf_difference,clobber=True)
         #pyfits.update(psf_difference_fitsname,psf_difference,header=moon_header)
         #print "wrote psf difference image %s" %  psf_difference_fitsname
         
         
         image_length=moon_zoom.shape[0]
         image_height=moon_zoom.shape[1]
         moon_mask = np.zeros((image_length,image_height))
   
         #define the centre of the disk (moon position)
         a,b = (image_length/2)-1, (image_height/2)-1
   
         y,x = np.ogrid[-a:image_length-a, -b:image_height-b]
         mask = x*x + y*y <= moon_radius_pix*moon_radius_pix
         moon_mask[mask]=1

         #make rfi model mask
         rfi_radius_deg=rfi_model_mask_size/60.0
         rfi_radius_pix = np.round(rfi_radius_deg/pix_size_deg)
         print "rfi radius in pix is %s " % rfi_radius_pix
         rfi_model_mask = np.zeros((image_length,image_height))
         rfi_mask = x*x + y*y <= rfi_radius_pix*rfi_radius_pix
         rfi_model_mask[rfi_mask]=1

         #do some maths
         #following Harish's notation:
         #D = moon_zoom
         #M = moon_mask
         #P= psf_zoom
         #G= M convol P

         ##Rather than a convolution, take the FFT of both images and multiply them, then inverse FFT
         #moon_mask_fft=fftpack.fft2(moon_mask)
         #psf_fft=fftpack.fft2(psf_zoom)
         #G_fourier=moon_mask_fft*psf_fft
         #G_image=np.real(fftpack.ifft2(G_fourier))
         #G_image_shift=fftpack.fftshift(G_image)


         #G=signal.convolve2d(moon_mask,psf_zoom)
         G_image_shift=signal.fftconvolve(moon_mask,psf_zoom_jyppix,mode='same')

         #Vectorise
         vec_D=moon_minus_sky_jyppix.flatten('F')
         vec_G=G_image_shift.flatten('F')
         vec_PSF=psf_zoom_jyppix.flatten('F')

         H=[vec_G,vec_PSF]

         #using statsmodels
         X = np.array(H).T
         X_const=sm.add_constant(X)
         #results = sm.OLS(endog=vec_D, exog=X_const).fit()
         #beta_hat=results.params
         #print results.summary()


         ##Try again using sklearn
         #reg = linear_model.LinearRegression(fit_intercept=True)
         #reg.fit(X,vec_D)
         ##print reg.coef_
         ##print reg.intercept_

         # Now do just with numpy:
         #beta_hat = np.linalg.lstsq(X_const,vec_D)[0]
         #print beta_hat

         #Do it using the MLE from Vedantham et al 2015  (eqn 17)
         #H_matrix=np.column_stack((vec_G,vec_PSF))
         #theta_hat=np.matmul(np.matmul(np.linalg.inv(np.matmul(H_matrix.T,H_matrix)),H_matrix.T),vec_D)
         #print "theta_hat is %s" % str(theta_hat)


         #Use whichever one for now
         #S_moon=beta_hat[1]
         #S_RFI=beta_hat[2]
         #S_moon=theta_hat[0]
         #S_RFI=theta_hat[1]


         #print "First test: S_moon is %s Jy and S_RFI is %s Jy" % (S_moon, S_RFI)

         #Woohoo! all give basically the same answer and it is about what we expect, the Moon disk is close to zero (constant offset is negative so that is good) and the RF is about 29 Jy

         #Create reconstructed moon dirty image with RFI removed.
         #reconstructed_moon=S_moon*G_image_shift

         #Show earthshine (RFI) alone:
         #RFI_alone=S_RFI*psf_zoom_jyppix

         #Make a residual map (vedantham et al 2015 eqn 18)
         #theta_hat=[[S_moon],[S_RFI]]
         #H_theta_hat=S_moon*vec_G + S_RFI*vec_PSF
         #vec_R=vec_D - H_theta_hat
         #R=np.reshape(vec_R, psf_zoom.shape, order='F')
         #equivalent to just:
         #R=moon_minus_sky_jyppix-reconstructed_moon-RFI_alone
         
         ##############
         #What I've done isn't right - need to re-do the maths with this new convolution?
         #residuals show we are oversubtracting in the middle and under-subtracting 
         #at the edges so the RFI must be broader than the PSF (as expected)
         #So instead of using the psf * S_RFI as the model, use a Gaussian of width 3.86 arcmin (Vedantham et al 2015)
         #convolved with the PSF
         #gaussian_fwhm_deg=3.86/60.0
         gaussian_fwhm_deg=3.75/60.0
         gaussian_fwhm_pix = np.round(gaussian_fwhm_deg/pix_size_deg)
         
         #For the broadened RFI, try a disk instead of a Gaussian 
         #old:
         #RFI_broadening=makeGaussian(image_length, fwhm = gaussian_fwhm_pix, center=(a+1,b+1))
         #new:
         RFI_broadening=rfi_model_mask
         #limit extent of broadened rfi model
         #RFI_broadening=RFI_broadening*rfi_model_mask
         ##Convolve this with the PSF
         #RFI_broadening_fft=fftpack.fft2(RFI_broadening)
         #RFI_convolved_fourier=psf_fft*RFI_broadening_fft
         #RFI_convolved_image=np.real(fftpack.ifft2(RFI_convolved_fourier))
         #shift and normalise - this is the new PSF for RFI
         #RFI_convolved_shift=fftpack.fftshift(RFI_convolved_image)/np.max(RFI_convolved_image)
         
         RFI_convolved_image=signal.fftconvolve(psf_zoom,RFI_broadening,mode='same')
         ##Don't need this normalising step
         #RFI_convolved_shift=(RFI_convolved_image)/np.max(RFI_convolved_image)
         RFI_convolved_shift=RFI_convolved_image
         #convert Jy per pix
         RFI_convolved_shift_jyppix=RFI_convolved_shift/n_pixels_per_beam   
         
         vec_RFI=RFI_convolved_shift_jyppix.flatten('F')
         
         ##Use this for broadened RFI
         #H2=[vec_G,vec_RFI]
         if not do_RFI_broadening:
            #Go back to pt source RFI for test
            vec_RFI=vec_PSF
         
         H2=[vec_G,vec_RFI]
         
         
         #using statsmodels
         X2 = np.array(H2).T
         X2_const=sm.add_constant(X2)
         
         # Now do just with numpy:
         #beta_hat2 = np.linalg.lstsq(X2_const,vec_D)[0]

         #using vedantham MLE as above
         H2_matrix=np.column_stack((vec_G,vec_RFI))
         theta_hat2=np.matmul(np.matmul(np.linalg.inv(np.matmul(H2_matrix.T,H2_matrix)),H2_matrix.T),vec_D)
         print "theta_hat2 is %s Jy/pix" % str(theta_hat2)

        
         #S_moon2=beta_hat2[1]
         #S_RFI2=beta_hat2[2]
         S_moon2=theta_hat2[0]
         S_RFI2=theta_hat2[1]         

         #Convert to total flux in Jy for Moon (assume for now RFI is a point source therefore flux in Jy/pix is same as total flux is Jy)
         S_moon2_tot_Jy=np.sum(S_moon2*moon_mask)
         print "Initial total moon flux density is %s Jy" % S_moon2_tot_Jy
         if not do_RFI_broadening:
            #for non-point source version:
            RFI_alone2=S_RFI2
            print "Initial total RFI flux density assuming pt source is %s Jy" % RFI_alone2
         else:
            RFI_alone2=np.sum(S_RFI2*RFI_broadening)
         #enforce Srfi is always positive
         if (S_RFI2 <0.0):
            print "Srfi is negative, but this is unphysical, enforcing Srfi to be positive"
            vec_RFI=-1.0*vec_RFI
         #enforce S_moon to be negative if we are in coarse chan 69 or 93 
         if (S_moon2 >= 0.0 and centre_chan <= 93):
            #repeat the analysis above but reverse the sign of G (i.e. equivalent to making the mask M negative?)
            print "Smoon is positive, but frequency is below 150 MHz, enforcing Smoon to be negative"
            vec_G=-1.0*vec_G
         #enforce S_moon to be positive if we are in coarse chan 145 or 169 
         if (S_moon2 <= 0.0 and centre_chan >= 145):
            print "Smoon is negative, but frequency is above 150 MHz, enforcing Smoon to be positive"
            vec_G=-1.0*vec_G

         H2_matrix=np.column_stack((vec_G,vec_RFI))
         theta_hat2=np.matmul(np.matmul(np.linalg.inv(np.matmul(H2_matrix.T,H2_matrix)),H2_matrix.T),vec_D)
         print "New theta_hat2 is %s Jy/pix for negative Smoon" % str(theta_hat2) 
         S_moon2=theta_hat2[0]
         #Need to conserve flux density: check RFI and moon flux density against max Moon
         S_RFI2=theta_hat2[1]

         #
         
         if not do_RFI_broadening:
            RFI_alone2=S_RFI2
         else:
            RFI_alone2=np.sum(S_RFI2*RFI_broadening) 
         
         #Calculate the errors:
         covariance_matrix_theta=(difference_rms_jy_per_pixel**2)*(np.linalg.inv(np.matmul(H2_matrix.T,H2_matrix)))
         #print "covariance_matrix_theta"
         #print covariance_matrix_theta
         S_moon2_error=np.sqrt(covariance_matrix_theta[0,0])
         
         #Convert to total flux in Jy for Moon (assume for now RFI is a point source therefore flux in Jy/pix is same as total flux is Jy)
         S_moon2_tot_Jy=np.sum(S_moon2*moon_mask)
         #Error from rms in Jy/pixel
         S_moon2_tot_Jy_error=np.sum(S_moon2_error*moon_mask)
         #print "S_moon2_tot_Jy_error"
         #print S_moon2_tot_Jy_error
         
         #print "Final Smoon, with RFI broadening: is %s Jy/pix and S_RFI is %s Jy/pix for moon obsid %s chan %s" % (S_moon2, S_RFI2,on_moon_obsid,chan_string)
         if not do_RFI_broadening:
            print "Final total Moon flux density, without RFI broadening: is %s Jy +- %s Jy and S_RFI is %s Jy for moon obsid %s chan %s" % (S_moon2_tot_Jy,S_moon2_tot_Jy_error, RFI_alone2,on_moon_obsid,chan_string)
         else:
            print "Final total Moon flux density, with RFI broadening: is %s Jy +- %s Jy and S_RFI is %s Jy for moon obsid %s chan %s" % (S_moon2_tot_Jy,S_moon2_tot_Jy_error, RFI_alone2,on_moon_obsid,chan_string)
         #print "Final total Moon flux density, with RFI broadening: is %s Jy and S_RFI is %s Jy for moon obsid %s chan %s" % (S_moon2_tot_Jy, RFI_alone2,on_moon_obsid,chan_string)
         
         reconstructed_moon2=S_moon2*G_image_shift
         R2=moon_minus_sky_jyppix-reconstructed_moon2-RFI_alone2
         
         #place values in the arrays
         if (chan != 'MFS'):
             Smoon_spectrum_values[chan,obsid_index]=S_moon2_tot_Jy
             Srfi_spectrum_values[chan,obsid_index]=RFI_alone2    

   #work out average values and standard deviations
   #Smoon_spectrum_values=np.empty([n_chans,n_obs])
   #Srfi_spectrum_values=np.empty([n_chans,n_obs])
   #Smoon_average_stddev_spectrum=np.empty([n_chans,2])
   #Srfi_average_stddev_spectrum=np.empty([n_chans,2])
   #print Smoon_spectrum_values
   Smoon_average_stddev_spectrum[:,0]=np.nanmean(Smoon_spectrum_values, axis=1)
   std_dev_Smoon=np.nanstd(Smoon_spectrum_values, axis=1)
   #count the number of non-nan data points
   #~ inverts the boolean matrix returned from np.isnan
   check_non_zero_array=~np.isnan(Smoon_spectrum_values)
   Smoon_data_points=np.sum(check_non_zero_array,axis=1)
   print "Number of Smoon data points is %s" % Smoon_data_points
   std_dev_of_mean_Smoon=std_dev_Smoon/np.sqrt(Smoon_data_points)
   
   #Go back to using std dev rather than std dev of mean
   #Smoon_average_stddev_spectrum[:,1]=std_dev_of_mean_Smoon
   Smoon_average_stddev_spectrum[:,1]=std_dev_Smoon
   #print "S_moon average and std dev of mean for each chan:" 
   #print Smoon_average_stddev_spectrum
   Srfi_average_stddev_spectrum[:,0]=np.nanmean(Srfi_spectrum_values, axis=1)
   std_dev_Srfi=np.nanstd(Srfi_spectrum_values, axis=1)
   check_non_zero_array_rfi=~np.isnan(Srfi_spectrum_values)
   Srfi_data_points=np.sum(check_non_zero_array_rfi,axis=1)
   std_dev_of_mean_Srfi=std_dev_Srfi/np.sqrt(Srfi_data_points)

   #Go back to using std dev rather than std dev of mean 
   #Srfi_average_stddev_spectrum[:,1]=std_dev_of_mean_Srfi
   Srfi_average_stddev_spectrum[:,1]=std_dev_Srfi
   #print "S_rfi average and std dev of mean for each chan:"
   #print Srfi_average_stddev_spectrum
   
   #what freqs?
   if centre_chan==69:
      freq_array=(np.arange(24)+57)*1.28
   if centre_chan==93:
      freq_array=(np.arange(24)+57+24)*1.28
   if centre_chan==121:
      freq_array=(np.arange(24)+57+24+24+4)*1.28
   if centre_chan==145:
      freq_array=(np.arange(24)+57+24+24+24+4)*1.28
   if centre_chan==169:
      freq_array=(np.arange(24)+57+24+24+24+24+4)*1.28

   #Now calculate the inferred background temperature in kilokelvin
   #Tb = T_moon + 160.0*(freq_array/60.0)**(-2.24) - ((10.0**(-26))*(c**2)*Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(freq_array*10**6)**2)
   Tb = T_moon + T_refl_gal_150*(freq_array/150.0)**(refl_gal_alpha) - (10.0**(-26)*c**2*Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(freq_array*10**6)**2)
   Tb_error=abs( (10.0**(-26)*c**2*Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
   #print "Tb in K is"
   #print Tb
   #print Tb_error
   
   #put in to the big arrays:
   if (centre_chan==69):
      big_Smoon_average_stddev_spectrum[57-57:57-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[57-57:57-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[57-57:57-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[57-57:57-57+24,1]=Srfi_average_stddev_spectrum[:,1]
   if (centre_chan==93):
      big_Smoon_average_stddev_spectrum[81-57:81-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[81-57:81-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[81-57:81-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[81-57:81-57+24,1]=Srfi_average_stddev_spectrum[:,1]
   if (centre_chan==121):
      big_Smoon_average_stddev_spectrum[109-57:109-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[109-57:109-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[109-57:109-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[109-57:109-57+24,1]=Srfi_average_stddev_spectrum[:,1]
   if (centre_chan==145):
      big_Smoon_average_stddev_spectrum[133-57:133-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[133-57:133-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[133-57:133-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[133-57:133-57+24,1]=Srfi_average_stddev_spectrum[:,1]
   if (centre_chan==169):
      big_Smoon_average_stddev_spectrum[157-57:157-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[157-57:157-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[157-57:157-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[157-57:157-57+24,1]=Srfi_average_stddev_spectrum[:,1]


   ##now plot the  Smoon with std dev as error bars for each subband
   plt.clf()
   Smoon_plot=plt.figure(1)
   plt.errorbar(freq_array,Smoon_average_stddev_spectrum[:,0],yerr=Smoon_average_stddev_spectrum[:,1])
   ##plt.plot(freq_array,Smoon_average_stddev_spectrum[:,0])
   plt.title('Moon flux density vs frequency for MWA')
   plt.ylabel('Mean Moon Flux Density (Smoon in Jy)')
   plt.xlabel('Frequency (MHz)')
   Smoon_plot.savefig('Smoon_plot.png')

   
   #Plot the inferred and predicted background temp
   plt.clf()
   Tb_plot=plt.figure(2)
   plt.errorbar(freq_array,Tb,yerr=Tb_error)
   plt.plot(freq_array,T_sky_predicted)
   plt.title('Inferred backroud temperature vs frequency for MWA')
   plt.ylabel('Inferred background temperature (Tb in K)')
   plt.xlabel('Frequency (MHz)')
   Tb_plot.savefig('Tb_plot.png')

   print freq_array
   print Smoon_average_stddev_spectrum[:,0]

 #save the big array for later use
 np.save(big_Smoon_average_stddev_spectrum_npy_filename, big_Smoon_average_stddev_spectrum)
 np.save(big_Srfi_average_stddev_spectrum_npy_filename,big_Srfi_average_stddev_spectrum)
 np.save(big_Smoon_predicted_npy_filename,big_Smoon_predicted_values)
 np.save(big_Ssky_predicted_npy_filename,big_Ssky_predicted_values)

#Now plot the big array
if (plot_only):
   big_Srfi_average_stddev_spectrum=np.load(big_Srfi_average_stddev_spectrum_npy_filename)
   big_Smoon_average_stddev_spectrum=np.load(big_Smoon_average_stddev_spectrum_npy_filename)
   big_Ssky_predicted_values=np.load(big_Ssky_predicted_npy_filename)
   big_Smoon_predicted_values=np.load(big_Smoon_predicted_npy_filename)
else:
   pass

######################################################################
#Do all the fitting stuff first before plotting anything!
#Now calculate the inferred background temperature in kilokelvin
#Tb = T_moon + 160.0*(big_freq_array/60.0)**(-2.24) - ((10.0**(-26))*(c**2)*big_Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(big_freq_array*10**6)**2)
Tb = T_moon + T_refl_gal_150*(big_freq_array/150.0)**(refl_gal_alpha) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(big_freq_array*10**6)**2)
#Tb_error=abs(0.07*Tsky_150*(big_freq_array/150.0)**(-2.5) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
#Tb_error=abs(0.07*160*(big_freq_array/60.0)**(-2.24) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
Tb_error=abs((10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(big_freq_array*10**6)**2))
#print "Tb in K is"
#print Tb
#print Tb_error

#Now we have worked out the Tb and error we want to model it as power law
#We want the amplitude to be the amp at 150 MHz so use
#logx = np.log10(freq_array)
logx=np.log10(big_freq_array/150.0)
logy = np.log10(Tb)
Tb_error[Tb_error == 0] = 200.
print Tb
print Tb_error
logyerr = Tb_error / Tb

powerlaw2 = lambda x, amp, index: amp * (x**index)

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit,
                       args=(logx, logy, logyerr), full_output=1)

pfinal = out[0]
covar = out[1]
#print pfinal
#print covar

alpha_gal_measured = pfinal[1]
T_150_measured = 10.0**pfinal[0]

alpha_gal_measured_err = np.sqrt( covar[1][1] )
T_150_measured_err = np.sqrt( covar[0][0] ) * T_150_measured

T_sky_measured_fit = powerlaw2(big_freq_array/150.0, T_150_measured, alpha_gal_measured)


################

#Save and plot the inferred and predicted background temp
tb_filename="inferred_sky_background_temp_and_error.npy"
big_Tb_value_error[:,0]=Tb
big_Tb_value_error[:,1]=Tb_error
np.save(tb_filename,big_Tb_value_error)

plt.clf()
Tb_plot=plt.figure(6)
plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
#plt.plot(big_freq_array,T_sky_predicted,label="Predicted")
plt.plot(big_freq_array,T_sky_measured_fit,label="Powerlaw fit")
plt.title('Inferred background temperature vs frequency for MWA')
plt.ylabel('Inferred background temperature (Tb in K)')
plt.xlabel('Frequency (MHz)')
plt.text(150, 500, 'Temp_150MHz = %5.2f +/- %5.3f' % (T_150_measured, T_150_measured_err))
plt.text(150, 400, 'Index = %5.2f +/- %5.3f' % (alpha_gal_measured, alpha_gal_measured_err))
plt.legend(loc=1)
Tb_plot.savefig('big_inferred_Tb_plot.png')

#predicted sky temp - now using model GSM (pyGSM)
#T_sky_predicted=Tsky_150*((big_freq_array/150.0)**(-2.5))
### Could try out healpix_util (https://github.com/esheldon/healpix_util)...
#for now just use what you did for gsm_reflection...
#FOR TESTING 
########test_freq_array=np.array([150,160])
########big_freq_array=test_freq_array


#print pixel_coordinate_transform_map_ra[0]
#print pixel_coordinate_transform_map_dec[0]
if new_model_galaxy:
   if use_model_gsm:   
      model_galaxy(big_freq_array,'gsm') 
   if use_model_lfsm:
      model_galaxy(big_freq_array,'lfsm')

   
#load the average temp file
average_gsm_temp_data=np.load(average_temp_filename)
temp_error=gsm_fractional_error*average_gsm_temp_data
#We want the amplitude to be the amp at 150 MHz so use
#logx = np.log10(freq_array)
logx=np.log10(big_freq_array/150.0)
logy = np.log10(average_gsm_temp_data)
logyerr = temp_error / average_gsm_temp_data

powerlaw = lambda x, amp, index: amp * (x**index)

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit,
                       args=(logx, logy, logyerr), full_output=1)

pfinal = out[0]
covar = out[1]
#print pfinal
#print covar

alpha_gal = pfinal[1]
T_150 = 10.0**pfinal[0]

alpha_gal_err = np.sqrt( covar[1][1] )
T_150_err = np.sqrt( covar[0][0] ) * T_150

T_sky_fit = powerlaw(big_freq_array/150.0, T_150, alpha_gal)


#plot this stuff:
plt.clf()
if use_model_gsm:
   plot_filename="best_fit_powerlaw_gsm.png"
   plot_title='Best Fit Power Law GSM2016'
   plot_ylabel='Average GSM Temp (K)'
if use_model_lfsm:  
   plot_filename="best_fit_powerlaw_lfsm.png"
   plot_title='Best Fit Power Law LFSM'
   plot_ylabel='Average LFSM Temp (K)'
fit_plot=plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(big_freq_array, powerlaw(big_freq_array/150.0, T_150, alpha_gal))     # Fit
plt.errorbar(big_freq_array, average_gsm_temp_data, yerr=temp_error, fmt='k.')  # Data
plt.text(150, 1000, 'Temp_150MHz = %5.2f +/- %5.2f' % (T_150, T_150_err))
plt.text(150, 700, 'Index = %5.2f +/- %5.2f' % (alpha_gal, alpha_gal_err))
plt.title(plot_title)
plt.xlabel('Frequency (MHz)')
plt.ylabel(plot_ylabel)
plt.xlim(70, 240)

plt.subplot(2, 1, 2)
plt.loglog(big_freq_array, powerlaw(big_freq_array/150.0, T_150, alpha_gal))
plt.errorbar(big_freq_array, average_gsm_temp_data, yerr=temp_error, fmt='k.')  # Data
plt.xlabel('Frequency (log scale)')
plt.ylabel('Temp (log scale)')
plt.xlim(70, 240)

fit_plot.savefig(plot_filename)

##save the new gsm map            
#plt.clf()
#map_title="GSM map masked"
#hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
##hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
#fig_name="gsm_map_%sMHz.png" % (str(freq_MHz))
#figmap = plt.gcf()
#figmap.savefig(fig_name,dpi=100) 


###############################
#Now work out predicted values for stuff:
###Predicted values:
####predicted sky temp - now using model GSM
big_T_sky_predicted=T_150*((big_freq_array/150.0)**(alpha_gal))
#print "Predicted Tsky is:"
#print T_sky_predicted

#predicted sky flux density for moon-disk area on sky
big_Ssky_predicted_values=(2.0*k*big_T_sky_predicted*Omega)/((300.0/big_freq_array)**(2)*10.0**(-26))
#print "Predicted Ssky is:"
#print S_sky_predicted

big_T_moon_predicted=np.zeros(len(big_freq_array)) + T_moon +  (T_refl_gal_150)*(big_freq_array/150.0)**(refl_gal_alpha)

#predicted sky flux density for moon-disk area on sky
big_S_moon_predicted=(2.0*k*big_T_moon_predicted*Omega)/((300.0/big_freq_array)**(2)*10.0**(-26))
#print "Predicted Smoon is:"
#print S_moon_predicted

#Assume that we know the Galactic temp (e.g.from modelling above) and estimate the brightness temp of the Moon!
#Convert to Moon brightness temperature units 
Tmoon_measured_average= (big_Smoon_average_stddev_spectrum[:,0]*(300.0/big_freq_array)**2*1e-26) / (2.0*k*Omega) + big_T_sky_predicted
Tmoon_measured_stddev= (big_Smoon_average_stddev_spectrum[:,1]*(300.0/big_freq_array)**2*1e-26) / (2.0*k*Omega)


#predicted moon-sky difference
big_predicted_moon_sky_difference=big_Smoon_predicted_values-big_Ssky_predicted_values


##############
plt.clf()
big_Smoon_plot=plt.figure(2)
plt.errorbar(big_freq_array,big_Smoon_average_stddev_spectrum[:,0],yerr=big_Smoon_average_stddev_spectrum[:,1],label="Measured")
plt.errorbar(big_freq_array,big_predicted_moon_sky_difference,label="Predicted")
plt.title('Moon Flux Density vs Frequency for MWA')
plt.ylabel('Mean Moon Flux Density (Jy)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
big_Smoon_plot.savefig('big_Smoon_and_predicted_moon_sky.png')


plt.clf()
big_Tmoon_plot=plt.figure(3)
plt.errorbar(big_freq_array,Tmoon_measured_average,yerr=Tmoon_measured_stddev,label="Measured")
plt.errorbar(big_freq_array,big_T_moon_predicted,label="Predicted")
plt.title('Moon Brightness Temperature vs Frequency for MWA')
plt.ylabel('Estimated Mean Moon Brightness Temperature (K)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
big_Tmoon_plot.savefig('big_Tmoon_and_predicted_moon_sky_T.png')

plt.clf()
big_Tmoon_plot=plt.figure(4)
plt.errorbar(big_freq_array,Tmoon_measured_average,yerr=Tmoon_measured_stddev,label="Measured")
plt.errorbar(big_freq_array,big_T_moon_predicted,label="Predicted")
axes = plt.gca()
#axes.set_xlim([xmin,xmax])
axes.set_ylim([-100,500])
plt.title('Moon Brightness Temperature vs Frequency for MWA')
plt.ylabel('Mean Moon Brightness Temperature (K)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
big_Tmoon_plot.savefig('big_Tmoon_and_predicted_moon_sky_T_zoom.png')

#Plot the RFI
plt.clf()
big_Smoon_plot=plt.figure(5)
plt.errorbar(big_freq_array,big_Srfi_average_stddev_spectrum[:,0],yerr=big_Smoon_average_stddev_spectrum[:,1],label="Measured")
#plt.errorbar(big_freq_array,big_predicted_moon_sky_difference,label="Predicted")
plt.title('Reflected Moon RFI Flux Density vs Frequency for MWA')
plt.ylabel('Mean ReflectedcRFI Flux Density (Jy)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=1)
big_Smoon_plot.savefig('big_Srfi.png')




#bask in glory

