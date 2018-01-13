#python code to model moon spectrum, remove rfi component ets.
import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
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
from pygsm import GlobalSkyModel
import healpy as hp

def model_moon(options):
   stokes=options.stokes
   #Choose one only: LFSM or GSM ( or GSM2016)
   use_model_gsm=True
   use_model_gsm2016=False
   use_model_lfsm=False
   
   new_model_galaxy=False
   use_mwa_beams=True
   use_gaussian_beams=False
   generate_new_beam_maps=False
   #set if you just want to plot already saved data
   plot_only=False
   #RFI broadening or not:
   do_RFI_broadening=True
   #Set true to force Srfi to always be positive (or should we just set it to zero if found negative?)
   enforce_Srfi_positive=True
   #Set true to enforce Smoon to be negative below 150 MHz - don't do this as Smoon (disk component of modelling) can be positive 
   #below 150 MHz due to diffuse component of RFI reflection
   enforce_Smoon_negative=False
   
   
   if use_model_gsm2016:
      NSIDE=1024 #nside for healpix maps (pygsm default is 512 for pygsm and 1024 for pygsm2016 (checked))
      if use_gaussian_beams:
         beam_map_filename_base="/data/moon/primary_beam/gaussian_beam_gsm2016"
      if use_mwa_beams:
         beam_map_filename_base="/data/moon/primary_beam/1443286527_antpat" 
   if use_model_gsm:
      NSIDE=512 #nside for healpix maps (pygsm default is 512 for pygsm and 1024 for pygsm2016 (checked))
      if use_gaussian_beams:
         beam_map_filename_base="/data/moon/primary_beam/gaussian_beam_gsm"
      if use_mwa_beams:
         beam_map_filename_base="/data/moon/primary_beam/1443286527_antpat"  
   if use_model_lfsm:
      NSIDE=128
      if use_gaussian_beams:
         beam_map_filename_base="/data/moon/primary_beam/gaussian_beam_lfsm"
      if use_mwa_beams:
         beam_map_filename_base="/data/moon/primary_beam/1443286527_antpat"
   #downloaded pixel transformation fits file from (https://lambda.gsfc.nasa.gov/toolbox/tb_pixelcoords.cfm)
   #for pygsm
   #pixel_coordinate_file="/data/code/healpix/pixel_coords_map_ring_celestial_res9.fits"
   #for pygsm2016:
   pixel_coordinate_file_gsm2016="/data/code/healpix/pixel_coords_map_ring_celestial_res10.fits"
   pixel_coordinate_file_gsm="/data/code/healpix/pixel_coords_map_ring_celestial_res9.fits"
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
      global_average_temp_filename="global_average_lfsm_temp.npy"
      fixed_beam_average_temp_filename="fixed_beam_average_lfsm_temp.npy"
      chromaticity_correction_filename="chromaticity_correction_lfsm.npy"
      gsm_fractional_error=0.15
   if use_model_gsm2016:
      average_temp_filename="average_gsm2016_temp_for_moon_region.npy"
      global_average_temp_filename="global_average_gsm2016_temp.npy"
      fixed_beam_average_temp_filename="fixed_beam_average_gsm2016_temp.npy"
      chromaticity_correction_filename="chromaticity_correction_gsm2016.npy"
      gsm_fractional_error=0.15
   if use_model_gsm:
      average_temp_filename="average_gsm_temp_for_moon_region.npy"
      global_average_temp_filename="global_average_gsm_temp.npy"
      fixed_beam_average_temp_filename="fixed_beam_average_gsm_temp.npy"
      chromaticity_correction_filename="chromaticity_correction_gsm.npy"
      gsm_fractional_error=0.15
      
   #set if using already cropped images
   use_cropped_images=False
   
   #Set if using the small images. 'small images are 2048 x 20148' 
   #But still need to crop as especially at the higher freqs, the primary beam is a lot narrower and PB errors affect the difference images.
   small_images=True #i.e. the standard 2048x2048 images I am using now
   ionpeeled_images=False

   #xstart_moon_small,xend_moon_small,ystart_moon_small,yend_moon_small=512,1536,512,1536
   #The moon (noticeably the pt source of specular RFI) is smeared (by about 1 pixel) by the motion of the Moon over the 2-min snapshot (in 1 min the Moon moves ~.0089 deg relative to the background, our pixel size is .0085)
   #So on average the centroid of the specular RFI component and the Moon will be 1 pixel away from the centre. SO crop the image so that the 
   #average position of the Moon is in the centre)
   xstart_moon_small,xend_moon_small,ystart_moon_small,yend_moon_small=511,1535,513,1537
   
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
   T_moon=230. #K
   #Moon RA 22:49:02 DEC -5:34:25 (2015-09-26) (RA 342.45 DEC -5.573)
   
   T_refl_gal_150=25.26 # temp of reflected galactic emisssion at 150 from reflection modeling (GSM2016) (takes into account albedo of 0.07)
   refl_gal_alpha=-2.50
   #RFI model mask radius in arcmin - limits the extent of the gausssian reflection model
   #rfi_model_mask_size=2.0
   rfi_model_mask_size=(3.75/2)
   
   #filename for big numpy arrays
   big_Smoon_average_stddev_spectrum_npy_filename="big_Smoon_average_stddev_spectrum.npy"
   big_Srfi_average_stddev_spectrum_npy_filename="big_Srfi_average_stddev_spectrum.npy"
   big_rms_residual_images_npy_filename="big_rms_residual_images.npy"
   
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
       
   def set_shared_ylabel(a, ylabel, labelpad = 0.01):
    """Set a y label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    ylabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""

    f = a[0].get_figure()
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    top = a[0].get_position().y1
    bottom = a[-1].get_position().y0

    # get the coordinates of the left side of the tick labels 
    x0 = 1
    for at in a:
        at.set_ylabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.yaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.inverse_transformed(f.transFigure)
        xt = bboxes.x0
        if xt < x0:
            x0 = xt
    tick_label_left = x0

    # set position of label
    a[-1].set_ylabel(ylabel)
    a[-1].yaxis.set_label_coords(tick_label_left - labelpad,(bottom + top)/2, transform=f.transFigure)


   #these aren't real beam maps, just a temporary solution until marcin builds his code to get the FEKO sim maps in healpix
   #for now they are just Gaussians in hpx
   def generate_beam_maps(frequency_array):
      if use_model_lfsm:
            pixel_coordinate_file=pixel_coordinate_file_lfsm
      if use_model_gsm2016:
            pixel_coordinate_file=pixel_coordinate_file_gsm2016
      if use_model_gsm:
            pixel_coordinate_file=pixel_coordinate_file_gsm
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
   
         #save the beam map array
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
      #for beam-weighted average
      average_temp_array=np.empty(len(frequency_array))
      #for global (isotropic antenna) average
      global_average_temp_array=np.empty(len(frequency_array))
      #for fixed beam (150 MHz) average
      fixed_beam_average_temp_array=np.empty(len(frequency_array))
      #chromaticity correction
      chromaticity_correction_array=np.empty(len(frequency_array))
      if (model_data=='gsm2016'):
         gsm = GlobalSkyModel2016(freq_unit='MHz')
      if (model_data=='gsm'):
         gsm = GlobalSkyModel(freq_unit='MHz')
      for freq_index,freq_MHz in enumerate(frequency_array):
         if use_gaussian_beams:
            beam_map_filename="%s_%.0fMHz.npy" % (beam_map_filename_base,freq_MHz)
            beam_map=np.load(beam_map_filename)
         if use_mwa_beams:
            beam_map_filename="%s_%07.2fMHz.fits" % (beam_map_filename_base,freq_MHz)
            beam_map=hp.read_map(beam_map_filename)
            #for fixed beam map 
            fixed_freq_MHz=170*1.28
            fixed_beam_map_filename="%s_%07.2fMHz.fits" % (beam_map_filename_base,fixed_freq_MHz)
            #print 'fixed_beam_map_filename %s' % fixed_beam_map_filename
            fixed_beam_map=hp.read_map(fixed_beam_map_filename)            
            if (NSIDE!=512):
               beam_map=hp.ud_grade(beam_map,nside_out=NSIDE)
               fixed_beam_map=hp.ud_grade(fixed_beam_map,nside_out=NSIDE)
         if (model_data=='lfsm'):
            #Just read in maps already generated using realizelfsm.py
            lfsm_data_filename="/data/moon/low_freq_sky_model_lwa/output-%.0f.dat" % freq_MHz
            galaxy_map=np.loadtxt(lfsm_data_filename)      
         if (model_data=='gsm2016'):
            galaxy_map=gsm.generate(freq_MHz)
            fixed_galaxy_map=gsm.generate(fixed_freq_MHz)
         if (model_data=='gsm'):
            #galaxy_map=gsm.generate(freq_MHz)
            #try reading the angelica data gsm maps from Marcin
            galaxy_map_filename="/data/moon/bighorns/sky_maps_angelica/1443286527_angelica_%07.2fMHz.fits" % freq_MHz
            print "using angelica map %s" % galaxy_map_filename
            galaxy_map=hp.read_map(galaxy_map_filename)
            fixed_galaxy_map=gsm.generate(fixed_freq_MHz)
         #global average
         mean_sky=np.mean(galaxy_map)
         print "global average %s" % mean_sky
         global_average_temp_array[freq_index]=mean_sky
         #calculate the beam-weighted mean of the sky
         sky_seen_by_MWA=galaxy_map*beam_map
         beam_weighted_mean_sky=np.sum(sky_seen_by_MWA)/np.sum(beam_map)
         print "beam_weighted_mean_sky %s freq %.0f MHz" % (beam_weighted_mean_sky,freq_MHz)
         average_temp_array[freq_index]=beam_weighted_mean_sky
         
         #calculate the FIXED beam-weighted mean of the sky
         fixed_sky_seen_by_MWA=fixed_galaxy_map*fixed_beam_map
         fixed_beam_weighted_mean_sky=np.sum(fixed_sky_seen_by_MWA)/np.sum(fixed_beam_map)
         sky_seen_by_MWA_with_fixed_sky_only=fixed_galaxy_map*beam_map
         beam_weighted_mean_with_fixed_sky_only=np.sum(sky_seen_by_MWA_with_fixed_sky_only)/np.sum(beam_map)
         print "Fixed (150 MHz) beam_weighted_mean_sky %s freq %.0f MHz" % (fixed_beam_weighted_mean_sky,fixed_freq_MHz)
         fixed_beam_average_temp_array[freq_index]=fixed_beam_weighted_mean_sky
         chromaticity_correction_array[freq_index]=beam_weighted_mean_with_fixed_sky_only/fixed_beam_weighted_mean_sky
         #chromaticity_correction_array[freq_index]=np.sum(sky_seen_by_MWA_with_fixed_sky_only)/np.sum(fixed_sky_seen_by_MWA)
         
         #for printing log
         #stay positive!
         galaxy_map=np.where(galaxy_map > 0.0, galaxy_map, np.nan)
         sky_above_horizon=np.where(beam_map > 0.0, galaxy_map, np.nan)
         sky_seen_by_MWA=np.where(sky_seen_by_MWA > 0.0, sky_seen_by_MWA, np.nan)   
         fixed_sky_seen_by_MWA=np.where(fixed_sky_seen_by_MWA > 0.0, fixed_sky_seen_by_MWA, np.nan)   
         beam_map=np.where(beam_map > 0.0, beam_map, np.nan)       
         fixed_beam_map=np.where(fixed_beam_map > 0.0, fixed_beam_map, np.nan) 
         
         #print min(galaxy_map)
         #print min(beam_map)
         #print min(sky_seen_by_MWA)
    
         ##save the beam multiplied by the sky map as png  
         min_temp=1
         max_temp=np.nanmax(sky_seen_by_MWA) 
         plt.clf()
         map_title="Sky multiplied by MWA beam %.0fMHz" % freq_MHz
         sky_seen_by_MWA_figname="sky_multiplied_by_MWA_beam_%.0fMHz.png" % freq_MHz
         #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         hp.mollview(map=sky_seen_by_MWA,coord='G',xsize=400,title=map_title,norm='linear',min=min_temp, max=max_temp) #norm="log"
         figmap = plt.gcf()
         figmap.savefig(sky_seen_by_MWA_figname,dpi=100) 
         plt.close()
   
         #save the galaxy
         plt.clf()
         map_title="sky_%.0fMHz" % freq_MHz
         sky_figname="sky_%.0fMHz.png" % freq_MHz
         #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         hp.mollview(map=galaxy_map,coord='G',xsize=400,title=map_title,norm='linear',min=400,max=50000)
         figmap = plt.gcf()
         figmap.savefig(sky_figname,dpi=100)
         plt.close()
         
         ##save the sky above horizon  as png           
         plt.clf()
         map_title="Sky above horizon %.0fMHz" % freq_MHz
         sky_above_horizon_figname="sky_above_horizon_%.0fMHz.png" % freq_MHz
         #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         hp.mollview(map=sky_above_horizon,coord='G',xsize=400,title=map_title,norm='log',min=100,max=50000) #norm="log"
         figmap = plt.gcf()
         figmap.savefig(sky_above_horizon_figname,dpi=100)
         plt.close()
   
         #save the beam
         plt.clf()
         map_title="MWA beam %.0fMHz" % freq_MHz
         beam_MWA_figname="MWA_beam_%.0fMHz.png" % freq_MHz
         #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         hp.mollview(map=beam_map,coord='G',xsize=400,title=map_title,norm='log')
         figmap = plt.gcf()
         figmap.savefig(beam_MWA_figname,dpi=100)
         plt.close()
   
   
      #print average_temp_array
      #save the average temp array
      np.save(average_temp_filename,average_temp_array)
      np.save(global_average_temp_filename,global_average_temp_array)
      np.save(fixed_beam_average_temp_filename,fixed_beam_average_temp_array)
      np.save(chromaticity_correction_filename,chromaticity_correction_array)
      
   #no longer used:
   def model_gsm(frequency_array):
      #for beam-weighted average
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
   
   def fit_powerlaw(temp_data, temp_data_error, freq_array):
      #returns 'T_sky_fit,'alpha','alpha_err','T_150','T_150_err'
      #We want the amplitude to be the amp at 150 MHz so use
      #logx = np.log10(freq_array)
      logx=np.log10(freq_array/150.0)
      logy = np.log10(temp_data)
      logyerr = temp_data_error / temp_data
   
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
   
      alpha = pfinal[1]
      T_150 = 10.0**pfinal[0]
   
      alpha_err = np.sqrt( covar[1][1] )
      T_150_err = np.sqrt( covar[0][0] ) * T_150
   
      T_sky_fit = powerlaw(freq_array/150.0, T_150, alpha)
   
      return {'T_sky_fit':T_sky_fit,'alpha':alpha,'alpha_err':alpha_err,'T_150':T_150,'T_150_err':T_150_err}
   
   def fit_line(data, data_error, freq_array):
      #returns a fit to y=ax +b
      def func(x, a, b):
         return a*x + b
      
      popt, pcov = optimize.curve_fit(func, freq_array, data, p0=None, sigma=data_error, absolute_sigma=False, check_finite=True) 
      fit_line=func(freq_array,popt[0],popt[1])
      return {'line':fit_line,'a':popt[0],'b':popt[1]}
  
  
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
   big_rms_residual_images=np.zeros(number_coarse_chans)
   big_Smoon_average_stddev_spectrum=np.zeros([number_coarse_chans,2])
   big_Srfi_average_stddev_spectrum=np.zeros([number_coarse_chans,2])
   big_Ssky_predicted_values=np.zeros(number_coarse_chans)
   big_Smoon_predicted_values=np.zeros(number_coarse_chans)
   big_Tb_value_error=np.zeros([number_coarse_chans,2])
   
   
   band_centre_chans=[69,93,121,145,169]
   #band_centre_chans=[69,93,121]
   
   if (not plot_only):
   
    for centre_chan in band_centre_chans:
      
      #set rms_thresholds for difference images
      rms_threshold=1.0
      
      #read the info files to find out how many observations there are and what the on and off-moon obsids are:
      on_moon_filename="/data/moon/2017/20150926_moon_%s.txt" % (str(centre_chan))
      off_moon_filename="/data/moon/2017/20150929_off_moon1_%s.txt" % (str(centre_chan))
      
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
         xstart_moon,xend_moon,ystart_moon,yend_moon=xstart_moon_small,xend_moon_small,ystart_moon_small,yend_moon_small
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
      rms_residual_images=np.zeros([n_chans,n_obs])
      rms_residual_images[:]=np.nan 
      Smoon_spectrum_values=np.zeros([n_chans,n_obs])
      Smoon_spectrum_values[:]=np.nan 
      Srfi_spectrum_values=np.zeros([n_chans,n_obs])
      Srfi_spectrum_values[:]=np.nan
      Smoon_average_stddev_spectrum=np.zeros([n_chans,2])
      Srfi_average_stddev_spectrum=np.zeros([n_chans,2])
      rms_residual_images_average=np.zeros(n_chans)
      
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
            moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s_dirty-%s.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string,stokes)
            off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s_dirty-%s.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string,stokes)
            on_moon_psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s-psf.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            off_moon_psf_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s-psf.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            #output fitsnames:
            moon_difference_fitsname="difference_%s_%s_on_off_moon_%s-%s-%s.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            psf_difference_fitsname="difference_%s_%s_psf_%s-%s-psf.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
            moon_zoom_fitsname="moon_zoom_%s_%s-%s-%s.fits" % (on_moon_obsid,str(centre_chan),chan_string,stokes)
            off_moon_zoom_fitsname="off_moon_zoom_%s_paired_with_%s_%s-%s-%s.fits" % (off_moon_obsid,on_moon_obsid,str(centre_chan),chan_string,stokes)
            rfi_modelled_fitsname="rfi_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            moon_modelled_fitsname="moon_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            residual_modelled_fitsname="residual_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
   
            #if (use_cropped_images): 
               #moon_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)
               #off_moon_fitsname="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I_cropped.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               #on_moon_psf_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)
               #moon_difference_fitsname="images/difference_%s_%s_cotter_20150929_moon_%s_peeled-%s-I_cropped.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
               
               #moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s_dirty-I_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               #off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s-%s_dirty-I_cropped.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               #psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon-%s-psf_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
               #moon_difference_fitsname="/data/moon/2017/difference_%s_%s_on_off_moon_%s-%s-I_cropped.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
            #else:
               #moon_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001_dirty_applied-I.fits"
               #moon_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (on_moon_obsid,str(centre_chan),chan_string)
               #off_moon_fitsname="1127572080_cotter_20150929_moon_93_track_off_moon_paired_1127313592_peeled-0001_dirty_applied-I.fits"
               #off_moon_fitsname="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
               #psf_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001-psf.fits"
               #psf_fitsname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (on_moon_obsid,str(centre_chan),chan_string)
               #difference image filename
               #moon_difference_fitsname="images/difference_%s_%s_cotter_20150929_moon_%s_peeled-%s-I.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
               #if (ionpeeled_images):
                  #moon_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
                  #off_moon_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
                  #on_moon_psf_fitsname="/data/moon/2017/20150926/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
                  #off_moon_psf_fitsname="/data/moon/2017/20150929/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s-psf.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
                  #moon_difference_fitsname="/data/moon/2017/difference_%s_%s_on_off_moon_%s_peeled-%s-I.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string)
                  #psf_difference_fitsname="/data/moon/2017/difference_%s_%s_psf_%s_peeled-%s-psf.fits" % (on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string) 
                  
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
            moon_zoom=moon_data[ystart_moon:yend_moon,xstart_moon:xend_moon]
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
            #beam area for 2d gussian 2 pi major minor / (8 ln 2) = 1.133 maj min ?
            beam_area_deg_sq=1.133 * bmaj_deg * bmin_deg
            n_pixels_per_beam=beam_area_deg_sq/pix_area_deg_sq
            #to convert from Jy/beam to Jy per pix, divide by n_pixels_per_beam
            
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
            off_moon_zoom=off_moon_data[ystart_moon:yend_moon,xstart_moon:xend_moon]
            
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
               print "rms of difference map %s is %s in Jy/beam. Discarding images." % (moon_difference_fitsname,difference_rms)
               #place values in the arrays
               if (chan != 'MFS'):
                  Smoon_spectrum_values[chan,obsid_index]=np.nan
                  Srfi_spectrum_values[chan,obsid_index]=np.nan
               continue  
            else:
               print "rms of difference map %s is %s in Jy/beam. Writing out images." % (moon_difference_fitsname,difference_rms)
               #write out the difference image
               pyfits.writeto(moon_difference_fitsname,moon_minus_sky,clobber=True)
               pyfits.update(moon_difference_fitsname,moon_minus_sky,header=moon_header)
               print "wrote  image %s" %  moon_difference_fitsname
               #write out the difference image in Jy per pix
               moon_difference_jyperpix_fitsname=moon_difference_fitsname.split('.')[0]+'_jyperpix.fits'
               pyfits.writeto(moon_difference_jyperpix_fitsname,moon_minus_sky_jyppix,clobber=True)
               #should really change header so it says it is in Jy pper pix....later
               pyfits.update(moon_difference_jyperpix_fitsname,moon_minus_sky_jyppix,header=moon_header)
               print "wrote  image %s" %  moon_difference_jyperpix_fitsname
               
               
               
               #write out the moon image
               pyfits.writeto(moon_zoom_fitsname,moon_zoom,clobber=True)
               pyfits.update(moon_zoom_fitsname,moon_zoom,header=moon_header)
               print "wrote  image %s" %  moon_zoom_fitsname
               #write out the off moon image
               pyfits.writeto(off_moon_zoom_fitsname,off_moon_zoom,clobber=True)
               pyfits.update(off_moon_zoom_fitsname,off_moon_zoom,header=moon_header)
               print "wrote  image %s" %  off_moon_zoom_fitsname
   
            #read in psf images
            if os.path.isfile(on_moon_psf_fitsname) and os.access(on_moon_psf_fitsname, os.R_OK):         
               psf_hdulist = pyfits.open(on_moon_psf_fitsname)
            else:
               print "Either file %s is missing or is not readable" % on_moon_psf_fitsname
               continue
            psf_data=psf_hdulist[0].data[0,0,:,:]
            psf_data=np.nan_to_num(psf_data)
            psf_header=psf_hdulist[0].header
            psf_zoom=psf_data[ystart_psf:yend_psf,xstart_psf:xend_psf]
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
            gaussian_fwhm_deg=3.86/60.0
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
   
            #Convert to total flux in Jy for Moon (assume for now RFI is a point source therefore flux in Jy/pix is same as total flux in Jy)
            S_moon2_tot_Jy=np.sum(S_moon2*moon_mask)
            print "Initial total moon flux density is %s Jy" % S_moon2_tot_Jy
            if not do_RFI_broadening:
               #for non-point source version:
               RFI_alone2=S_RFI2
               print "Initial total RFI flux density assuming pt source is %s Jy" % RFI_alone2
            else:
               RFI_alone2=np.sum(S_RFI2*RFI_broadening)
            #enforce Srfi is always positive
            #take this out and rerun just for kicks
            if enforce_Srfi_positive:
               if (S_RFI2 <0.0):
                  print "Srfi is negative, but this is unphysical, enforcing Srfi to be positive"
                  vec_RFI=-1.0*vec_RFI
            if enforce_Smoon_negative:
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
            
            ###Don't think I have done this correctly
            #I think for RFI we need to do (S_RFI2x1 pix in centre) convolved with psf
            # and for moon: (S_moon2xmoon_mask) convolved with psf
            if do_RFI_broadening:
               reconstructed_RFI=S_RFI2*RFI_convolved_shift_jyppix
            else:
               reconstructed_RFI=S_RFI2*psf_zoom_jyppix
               #instead try this:
               ############RFI_pt=np.zeros((image_length,image_height))
               ############RFI_pt[image_length/2-1,image_height/2-1]=1
               ############reconstructed_RFI=signal.fftconvolve(psf_zoom_jyppix,(S_RFI2*RFI_pt),mode='same')
               
            reconstructed_moon2=S_moon2*G_image_shift
            #instead try this:
            ##########reconstructed_moon2=signal.fftconvolve(psf_zoom_jyppix,(S_moon2*moon_mask),mode='same') 
            
            #residual=moon_minus_sky_jyppix-reconstructed_moon2-RFI_alone2
            residual=moon_minus_sky_jyppix-reconstructed_moon2-reconstructed_RFI
            residual_rms=np.sqrt(np.mean(np.square(residual)))
            
            
            #place values in the arrays
            if (chan != 'MFS'):
                Smoon_spectrum_values[chan,obsid_index]=S_moon2_tot_Jy
                Srfi_spectrum_values[chan,obsid_index]=RFI_alone2 
                rms_residual_images[chan,obsid_index]=residual_rms
            
            #write out the modelled rfi image
            pyfits.writeto(rfi_modelled_fitsname,reconstructed_RFI,clobber=True)
            pyfits.update(rfi_modelled_fitsname,reconstructed_RFI,header=moon_header)
            print "wrote image %s" %  rfi_modelled_fitsname
            
            pyfits.writeto(moon_modelled_fitsname,reconstructed_moon2,clobber=True)
            pyfits.update(moon_modelled_fitsname,reconstructed_moon2,header=moon_header)
            print "wrote image %s" %  moon_modelled_fitsname
   
            pyfits.writeto(residual_modelled_fitsname,residual,clobber=True)
            pyfits.update(residual_modelled_fitsname,residual,header=moon_header)
            print "wrote image %s" %  residual_modelled_fitsname         
            
            
            
            
      #work out average values and standard deviations
      #Smoon_spectrum_values=np.empty([n_chans,n_obs])
      #Srfi_spectrum_values=np.empty([n_chans,n_obs])
      #Smoon_average_stddev_spectrum=np.empty([n_chans,2])
      #Srfi_average_stddev_spectrum=np.empty([n_chans,2])
      #print Smoon_spectrum_values
      
      #rms
      rms_residual_images_average=np.nanmean(rms_residual_images,axis=1)
      
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
      
      #put in to the big arrays:
      if (centre_chan==69):
         big_Smoon_average_stddev_spectrum[57-57:57-57+24,0]=Smoon_average_stddev_spectrum[:,0]
         big_Smoon_average_stddev_spectrum[57-57:57-57+24,1]=Smoon_average_stddev_spectrum[:,1]
         big_Srfi_average_stddev_spectrum[57-57:57-57+24,0]=Srfi_average_stddev_spectrum[:,0]
         big_Srfi_average_stddev_spectrum[57-57:57-57+24,1]=Srfi_average_stddev_spectrum[:,1]
         big_rms_residual_images[57-57:57-57+24]=rms_residual_images_average
      if (centre_chan==93):
         big_Smoon_average_stddev_spectrum[81-57:81-57+24,0]=Smoon_average_stddev_spectrum[:,0]
         big_Smoon_average_stddev_spectrum[81-57:81-57+24,1]=Smoon_average_stddev_spectrum[:,1]
         big_Srfi_average_stddev_spectrum[81-57:81-57+24,0]=Srfi_average_stddev_spectrum[:,0]
         big_Srfi_average_stddev_spectrum[81-57:81-57+24,1]=Srfi_average_stddev_spectrum[:,1]
         big_rms_residual_images[81-57:81-57+24]=rms_residual_images_average
      if (centre_chan==121):
         big_Smoon_average_stddev_spectrum[109-57:109-57+24,0]=Smoon_average_stddev_spectrum[:,0]
         big_Smoon_average_stddev_spectrum[109-57:109-57+24,1]=Smoon_average_stddev_spectrum[:,1]
         big_Srfi_average_stddev_spectrum[109-57:109-57+24,0]=Srfi_average_stddev_spectrum[:,0]
         big_Srfi_average_stddev_spectrum[109-57:109-57+24,1]=Srfi_average_stddev_spectrum[:,1]
         big_rms_residual_images[109-57:109-57+24]=rms_residual_images_average
      if (centre_chan==145):
         big_Smoon_average_stddev_spectrum[133-57:133-57+24,0]=Smoon_average_stddev_spectrum[:,0]
         big_Smoon_average_stddev_spectrum[133-57:133-57+24,1]=Smoon_average_stddev_spectrum[:,1]
         big_Srfi_average_stddev_spectrum[133-57:133-57+24,0]=Srfi_average_stddev_spectrum[:,0]
         big_Srfi_average_stddev_spectrum[133-57:133-57+24,1]=Srfi_average_stddev_spectrum[:,1]
         big_rms_residual_images[133-57:133-57+24]=rms_residual_images_average
      if (centre_chan==169):
         big_Smoon_average_stddev_spectrum[157-57:157-57+24,0]=Smoon_average_stddev_spectrum[:,0]
         big_Smoon_average_stddev_spectrum[157-57:157-57+24,1]=Smoon_average_stddev_spectrum[:,1]
         big_Srfi_average_stddev_spectrum[157-57:157-57+24,0]=Srfi_average_stddev_spectrum[:,0]
         big_Srfi_average_stddev_spectrum[157-57:157-57+24,1]=Srfi_average_stddev_spectrum[:,1]
         big_rms_residual_images[157-57:157-57+24]=rms_residual_images_average
   
   
      ##now plot the  Smoon with std dev as error bars for each subband
      #plt.clf()
      #Smoon_plot=plt.figure(1)
      #plt.errorbar(freq_array,Smoon_average_stddev_spectrum[:,0],yerr=Smoon_average_stddev_spectrum[:,1])
      ###plt.plot(freq_array,Smoon_average_stddev_spectrum[:,0])
      #plt.title('Moon flux density vs frequency for MWA')
      #plt.ylabel('Mean Moon Flux Density (Smoon in Jy)')
      #plt.xlabel('Frequency (MHz)')
      #Smoon_plot.savefig('Smoon_plot.png')
   
      
      ##Plot the inferred and predicted background temp
      #plt.clf()
      #Tb_plot=plt.figure(2)
      #plt.errorbar(freq_array,Tb,yerr=Tb_error)
      #plt.plot(freq_array,T_sky_predicted)
      #plt.title('Inferred backroud temperature vs frequency for MWA')
      #plt.ylabel('Inferred background temperature (Tb in K)')
      #plt.xlabel('Frequency (MHz)')
      #Tb_plot.savefig('Tb_plot.png')
   
      #print freq_array
      #print Smoon_average_stddev_spectrum[:,0]
   
    #save the big array for later use
    np.save(big_Smoon_average_stddev_spectrum_npy_filename, big_Smoon_average_stddev_spectrum)
    np.save(big_Srfi_average_stddev_spectrum_npy_filename,big_Srfi_average_stddev_spectrum)
    np.save(big_Srfi_average_stddev_spectrum_npy_filename,big_Srfi_average_stddev_spectrum)
    np.save(big_rms_residual_images_npy_filename,big_rms_residual_images)
   
   #Now plot the big array
   if (plot_only):
      big_Srfi_average_stddev_spectrum=np.load(big_Srfi_average_stddev_spectrum_npy_filename)
      big_Smoon_average_stddev_spectrum=np.load(big_Smoon_average_stddev_spectrum_npy_filename)
      big_Srfi_average_stddev_spectrum[big_Srfi_average_stddev_spectrum < 1e-10] = 0.
      big_Srfi_average_stddev_spectrum[big_Srfi_average_stddev_spectrum > 1e+10] = 0.
      big_Smoon_average_stddev_spectrum[big_Smoon_average_stddev_spectrum < -1e+20] = 0.
      big_Smoon_average_stddev_spectrum[big_Smoon_average_stddev_spectrum > 1e+20] = 0.
      big_rms_residual_images=np.load(big_rms_residual_images_npy_filename)
   else:
      pass
   
   
   #work out how to incorporate the addtional rms errors from the residual images - what units are they in?

   #the averages are v small - must be in Jy/ per pix - what uncertainty in moon flux density does this translate too? Multiply by how many pix there are in a Moon!
   #no - if you do that the numbers are too big - must be in Jy/beam ... but then they are too small!
   ##remake the moon_mask
   #image_length=xend_moon_small-xstart_moon_small
   #image_height=yend_moon_small-ystart_moon_small
   #moon_mask = np.zeros((image_length,image_height))

   ##define the centre of the disk (moon position)
   #a,b = (image_length/2)-1, (image_height/2)-1

   #y,x = np.ogrid[-a:image_length-a, -b:image_height-b]
   #pix_size_deg = 0.0085
   #print "using pix_size_deg %s to calculate moon mask for Smoon error" % pix_size_deg
   #moon_radius_pix = np.round(moon_radius_deg/pix_size_deg)
   #mask = x*x + y*y <= moon_radius_pix*moon_radius_pix
   #moon_mask[mask]=1
   #print "np.sum(moon_mask) %s" %  str(np.sum(moon_mask))
   #big_rms_residual_images_convert_to_Smoon_error=np.sum(moon_mask)*big_rms_residual_images

   #print big_rms_residual_images
   #print big_Smoon_average_stddev_spectrum[:,1]
   # Smoon and plot as a 4 panel with rfi, and subtr linear fit, and S_moon_rfi_subtracted
   
   big_Smoon_average_stddev_spectrum[:,1][big_Smoon_average_stddev_spectrum[:,1] == 0] = 20.

   line_fit_result=fit_line(big_Smoon_average_stddev_spectrum[:,0],big_Smoon_average_stddev_spectrum[:,1],big_freq_array)
   
   S_moon_minus_linear_fit=big_Smoon_average_stddev_spectrum[:,0]-line_fit_result['line']
   
   #Work out what the percentage of the peak RFI from specular reflection is contaminating the rest of the disk (use the FM band only)
   #try doing on the whole band and see if there is a freq dependence
   #RFI_peak_frac_in_disk=big_Smoon_average_stddev_spectrum[:,0]/big_Srfi_average_stddev_spectrum[:,0]
   #just use FM band and extrapolate a linear fit!
   RFI_peak_frac_in_disk_fm_band=S_moon_minus_linear_fit[12:27]/big_Srfi_average_stddev_spectrum[12:27,0]
   line_fit_result_RFI=fit_line(RFI_peak_frac_in_disk_fm_band,RFI_peak_frac_in_disk_fm_band*.1,big_freq_array[12:27])
   big_RFI_line_fit=big_freq_array*line_fit_result_RFI['a']+line_fit_result_RFI['b']
   av_RFI_peak_frac_in_disk_fm_band=np.mean(RFI_peak_frac_in_disk_fm_band)
   
   #wavelength dependence fit normalised to 1 at 100 MHz (99.84)
   big_freq_array_100MHz_index=21
   print 'big freq_array at index %s is %s ' % (big_freq_array_100MHz_index,big_freq_array[big_freq_array_100MHz_index])
   #normalise the evans ratio instead, to the value at 100 Mhz
   #big_RFI_line_fit_normalised_100MHz=big_RFI_line_fit/big_RFI_line_fit[big_freq_array_100MHz_index]
   
   #wavelength dependence (evans 1969, p223, eqns 31,32, power law) normalised to RFI_peak_frac_in_disk_fm_band at 100 MHz
   ratio_evans=(big_freq_array/big_freq_array[big_freq_array_100MHz_index])**0.58
   ratio_evans_normalised=ratio_evans/ratio_evans[big_freq_array_100MHz_index]*big_RFI_line_fit[big_freq_array_100MHz_index]
   
   #This theoretical wavelength dependence is clearly wrong (diffuse to specular should increase more quickly with freq)
   #Evans measurements are at the wrong freq anyway, but we know that it is a powerlaw dependence anyway, not linear
   #So instead try fitting a power law!
   #only do this fitting if RFI_peak_frac_in_disk_fm_band is all positive:
   if any(t < 0 for t in RFI_peak_frac_in_disk_fm_band):
      print "RFI_peak_frac_in_disk_fm_band has negative elements, skipping power law fit"
   else:
      #normalisation is at 100 MHz
      logx=np.log10(big_freq_array[12:27]/100.0)
      logy = np.log10(RFI_peak_frac_in_disk_fm_band)

      RFI_peak_frac_in_disk_fm_band_error=RFI_peak_frac_in_disk_fm_band*0.1
      logyerr = RFI_peak_frac_in_disk_fm_band_error / RFI_peak_frac_in_disk_fm_band
   
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
   
      alpha = pfinal[1]
      ratio_100 = 10.0**pfinal[0]
   
      alpha_error = np.sqrt( covar[1][1] )
      ratio_100_error = np.sqrt( covar[0][0] ) * ratio_100
   
      ratio_fit = powerlaw2(big_freq_array[12:27]/100.0, ratio_100, alpha)
      ratio_fit_big = powerlaw2(big_freq_array/100.0, ratio_100, alpha)
      
   #try using evans ratio  (normalised to powerlaw fit - makes it worse
   ####ratio_evans_normalised=ratio_evans/ratio_evans[big_freq_array_100MHz_index]*big_RFI_line_fit[big_freq_array_100MHz_index]
   
   #Okay so now we can remove av_RFI_peak_frac_in_disk_fm_band*RFI from each chan in S_moon
   #S_moon_RFI_subtracted=big_Smoon_average_stddev_spectrum[:,0]-(big_RFI_line_fit*big_Srfi_average_stddev_spectrum[:,0])
   #Use the powerlaw fit instead (even though it makes no difference whatsoever to the final result)
   S_moon_RFI_subtracted=big_Smoon_average_stddev_spectrum[:,0]-(ratio_fit_big*big_Srfi_average_stddev_spectrum[:,0])


   #add in the residual image rms (this is small, makes no difference) and the RFI rms (actually I think this is wrong, double dipping - plus makes the spectrum flatter!
   #Aha - need to also multiply the rfi error by the linear fit used for the subtracftion I think
   #big_Srfi_error_total_for_Smoon=big_Srfi_average_stddev_spectrum[:,1]*big_RFI_line_fit
   #S_moon_error_total=np.sqrt(big_rms_residual_images**2+big_Smoon_average_stddev_spectrum[:,1]**2+big_Srfi_error_total_for_Smoon**2)
   #S_moon_error_total=np.sqrt(big_rms_residual_images**2+big_Smoon_average_stddev_spectrum[:,1]**2+big_Srfi_average_stddev_spectrum[:,1]**2)*big_RFI_line_fit
   #Use powerlaw fit to scale error instead
   print big_rms_residual_images
   #The rms of the residual images is in Jy per pixel so in order to convert this to an uncertainty on Smoon, 
   #you need to convert to (got this from the mwa magellenic cloud paper):
   #1. convert residual image units back to Jy/beam. 
   #2. measure the rms (in Jy/beam)
   #3. sigma_integratedfluxdensity = rms x sqrt(N), where N is the integrated area in units of synthesised beam area
   # so for us N=moon_area / synthesised beam area = pi*(moon_radius_deg)^2 / 1.33*bmaj*bmin ? 
   # bmaj and bmin change with freq! so need to get these as a function of freq. so root N goes up with freq, which is good as we need
   # the errors to be larger at higher freq!
   
   
   S_moon_error_total=np.sqrt(big_rms_residual_images**2+big_Smoon_average_stddev_spectrum[:,1]**2+(big_Srfi_average_stddev_spectrum[:,1]*ratio_fit_big)**2)
   
   #S_moon_error_total=big_Smoon_average_stddev_spectrum[:,1] 
   #print S_moon_error_total

   
   plot_filename='big_S_moon_and_rfi.png'
   plot_title='S_moon and RFI'
   plt.clf()

   f,a = plt.subplots(3, sharex=True, gridspec_kw={'hspace':0})


   a[0].plot(big_freq_array,big_Smoon_average_stddev_spectrum[:,0],label="S_moon initial")
   #a[0].plot(big_freq_array,line_fit_result['line'],label="linear fit",linestyle='dashed')

   #a[0].text(150, 100, 'y=%5.2fx +%5.2f ' % (line_fit_result['a'],line_fit_result['a'])) 
   a[0].legend(loc=1)

   a[1].plot(big_freq_array,big_Srfi_average_stddev_spectrum[:,0],label="reflected RFI")
   a[1].legend(loc=1)
   
   #a[2].plot(big_freq_array,S_moon_minus_linear_fit,label="S_moon-linear fit")
   a[2].plot(big_freq_array,S_moon_RFI_subtracted,label="S_moon RFI excised")
   a[2].legend(loc=4)
   
   #a[3].plot(big_freq_array,S_moon_RFI_subtracted,label="S_moon_minus_linear_fit_RFI_subtracted")
   ##a[3].set_ylim([0,5])
   #a[3].legend(loc=2)
   plt.xlabel('Frequency (MHz)')
  
   
   set_shared_ylabel(a, 'Flux density (Jy)')
   
   f.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   plt.close()


   #####
   ####Plot the theoretical lambda dependence (evans 1969, p223, eqns 31,32, power law) on top of the linear dependence from FM band
   plot_filename='big_specular_vs_diffuse_refl_ratio.png'
   plot_title='Ratio of diffuse to specular reflection'
   plt.clf()

   ratio_plot=plt.figure()


   plt.plot(big_freq_array,big_RFI_line_fit,label="Linear fit FM")
   plt.plot(big_freq_array,ratio_evans_normalised,label="ratio Evans 1969")
   #Also add in the data used to get the linear fit:
   plt.plot(big_freq_array[12:27],RFI_peak_frac_in_disk_fm_band,label="FM ratio data")
   #add in the powerlaw fit result
   plt.plot(big_freq_array,ratio_fit_big,label="Powerlaw fit FM")

   #a[0].text(150, 100, 'y=%5.2fx +%5.2f ' % (line_fit_result['a'],line_fit_result['a'])) 
   plt.legend(loc=1)

   #a[0].plot(big_freq_array,big_Srfi_average_stddev_spectrum[:,0],label="reflected RFI")

   
   plt.xlabel('Frequency (MHz)')
   plt.title(plot_title)
   
   plt.ylabel('Flux density (Jy)')
   
   ratio_plot.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   plt.close()
   
   
   ######################################################################
   #Do all the fitting stuff first before plotting anything!
   #Use the RFI subtracted S_moon!!!
   #Now calculate the inferred background temperature 
   #Tb = T_moon + 160.0*(big_freq_array/60.0)**(-2.24) - ((10.0**(-26))*(c**2)*big_Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(big_freq_array*10**6)**2)
   Tb = T_moon + T_refl_gal_150*(big_freq_array/150.0)**(refl_gal_alpha) - (10.0**(-26)*c**2*S_moon_RFI_subtracted)/(2*k*Omega*(big_freq_array*10**6)**2) - 2.725 #remove Tcmb to get Galactic spectrum
   #Tb_error=abs(0.07*Tsky_150*(big_freq_array/150.0)**(-2.5) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
   #Tb_error=abs(0.07*160*(big_freq_array/60.0)**(-2.24) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
   Tb_error=abs((10.0**(-26)*c**2*S_moon_error_total)/(2*k*Omega*(big_freq_array*10**6)**2))
   #print "Tb in K is"
   #print Tb
   #print Tb_error
   
   #Now we have worked out the Tb and error we want to model it as power law
   #We want the amplitude to be the amp at 150 MHz so use
   #logx = np.log10(freq_array)
   
   #only do this fitting if Tb is all positive:
   if any(t < 0 for t in Tb):
      print "Tb has negative elements, skipping power law fit"
   else:
      logx=np.log10(big_freq_array/150.0)
      logy = np.log10(Tb)
      Tb_error[Tb_error == 0] = 200.
      Tb_error[Tb_error < 1e-10] = 200.
      Tb_error[Tb_error > 1e+10] = 200.
      #print Tb
      #print Tb_error
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
   
      #work out the residuals - what is the rms?
      residuals_big = Tb - T_sky_measured_fit
      #print len(residuals)
      residuals=residuals_big[48:]
      residual_rms=np.sqrt(np.mean(np.square(residuals)))
      #print Tb_error
      weights= np.reciprocal(Tb_error[48:].astype(np.float32))
      norm_weights=weights/np.sum(weights)
      #print weights
      weighted_residual_rms=np.sqrt(np.sum(norm_weights*np.square(residuals)))
      #print residual_rms
      #print weighted_residual_rms
      ################
   
      #Save and plot the inferred and predicted background temp
      tb_filename="inferred_sky_background_temp_and_error.npy"
      big_Tb_value_error[:,0]=Tb
      big_Tb_value_error[:,1]=Tb_error
      np.save(tb_filename,big_Tb_value_error)
   
      plt.clf()
      Tb_plot=plt.figure(6)
      plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
      plt.plot(big_freq_array,T_sky_measured_fit,label="Powerlaw fit")
      plt.title('Inferred background temperature vs frequency for MWA')
      plt.ylabel('Inferred background temperature (Tb in K)')
      plt.xlabel('Frequency (MHz)')
      plt.text(150, 1500, 'Temp_150MHz = %5.1f +/- %3.1f' % (T_150_measured, T_150_measured_err))
      plt.text(150, 1300, 'Index = %5.2f +/- %4.2f' % (alpha_gal_measured, alpha_gal_measured_err))
      plt.legend(loc=1)
      Tb_plot.savefig('big_inferred_Tb_plot.png')
      print "saved figure big_inferred_Tb_plot.png"
 
      #Plot subplot with log log space
      loglog_plot_filename='big_inferred_Tb_plot_loglog.png'
      plot_title='Inferred background temperature vs frequency for MWA'
      plot_ylabel='Inferred background \n temperature (Tb in K)'
      fit_plot=plt.figure(7)
      plt.subplot(2, 1, 1)
      plt.plot(big_freq_array, T_sky_measured_fit)     # Fit
      plt.errorbar(big_freq_array, Tb, yerr=Tb_error, fmt='k.')  # Data
      plt.text(150, 1000, 'Temp_150MHz = %5.2f +/- %5.2f' % (T_150_measured, T_150_measured_err))
      plt.text(150, 700, 'Index = %5.2f +/- %5.2f' % (alpha_gal_measured, alpha_gal_measured_err))
      plt.title(plot_title)
      plt.xlabel('Frequency (MHz)')
      plt.ylabel(plot_ylabel)
      plt.xlim(70, 240)
      plt.tight_layout()
   
      plt.subplot(2, 1, 2)
      plt.loglog(big_freq_array, T_sky_measured_fit)
      plt.errorbar(big_freq_array, Tb, yerr=Tb_error, fmt='k.')  # Data
      plt.xlabel('Frequency (log scale)')
      plt.ylabel('Temp (log scale)')
      plt.xlim(70, 240)
      plt.tight_layout()
   
      fit_plot.savefig(loglog_plot_filename)
      print "saved figure %s" % loglog_plot_filename
   
      plt.clf()
      plot_name="residual_plot.png"
      Tb_plot=plt.figure(8)
      #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
      plt.plot(big_freq_array[48:],residuals,label="Residuals")
      plt.title('Fit residuals  vs frequency for MWA')
      plt.ylabel('Residual temperature (K)')
      plt.xlabel('Frequency (MHz)')
      plt.text(180, -40, 'Residual RMS = %5.2f ' % (residual_rms))
      plt.text(180, -50, 'Weighted Residual RMS = %5.2f' % (weighted_residual_rms))
      plt.legend(loc=1)
      Tb_plot.savefig(plot_name)
      print "saved figure %s" %  plot_name
   
   
   
   
   #FOR TESTING 
   ########test_freq_array=np.array([150,160])
   ########big_freq_array=test_freq_array
   
   
   #print pixel_coordinate_transform_map_ra[0]
   #print pixel_coordinate_transform_map_dec[0]
   if new_model_galaxy:
      if use_model_gsm2016:   
         model_galaxy(big_freq_array,'gsm2016') 
      if use_model_gsm:   
         model_galaxy(big_freq_array,'gsm') 
      if use_model_lfsm:
         model_galaxy(big_freq_array,'lfsm')
   
      
   #load the average temp files
   average_gsm_temp_data=np.load(average_temp_filename)
   global_average_gsm_temp_data=np.load(global_average_temp_filename)
   fixed_beam_average_gsm_temp_data=np.load(fixed_beam_average_temp_filename)
   chromaticity_correction=np.load(chromaticity_correction_filename)
   #print 'chromaticity_correction:'
   #print chromaticity_correction
   #remove 2.727K for the CMB temp
   average_gsm_temp_data-=2.727
   
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
   
   
   ###############
   #Look at the difference between the global, and fixed beam, mean sky and the frequency-dependent beam-weighted mean sky
   ###############
   global_mean_difference= global_average_gsm_temp_data - average_gsm_temp_data
   fixed_beam_mean_difference= fixed_beam_average_gsm_temp_data - average_gsm_temp_data
   global_diff_residual_rms=np.sqrt(np.mean(np.square(global_mean_difference)))
   fixed_beam_diff_residual_rms=np.sqrt(np.mean(np.square(fixed_beam_mean_difference)))
   
   #plot stuff:
   plt.clf()
   if use_model_gsm:
      plot_filename="best_fit_powerlaw_gsm.png"
      plot_title='Best Fit Power Law GSM'
      plot_ylabel='Average GSM Temp (K)'
   if use_model_gsm2016:
      plot_filename="best_fit_powerlaw_gsm2016.png"
      plot_title='Best Fit Power Law GSM2016'
      plot_ylabel='Average GSM2016 Temp (K)'
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
   plt.tight_layout()
   
   plt.subplot(2, 1, 2)
   plt.loglog(big_freq_array, powerlaw(big_freq_array/150.0, T_150, alpha_gal))
   plt.errorbar(big_freq_array, average_gsm_temp_data, yerr=temp_error, fmt='k.')  # Data
   plt.xlabel('Frequency (log scale)')
   plt.ylabel('Temp (log scale)')
   plt.xlim(70, 240)
   plt.tight_layout()
   
   fit_plot.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   ##save the new gsm map            
   #plt.clf()
   #map_title="GSM map masked"
   #hp.orthview(map=gsm_map,coord='G',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
   ##hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
   #fig_name="gsm_map_%sMHz.png" % (str(freq_MHz))
   #figmap = plt.gcf()
   #figmap.savefig(fig_name,dpi=100) 
   
   plt.clf()
   if use_model_gsm:
      global_av_plot_filename="global_av_gsm.png"
      global_av_plot_title='Global average GSM'
      global_av_plot_ylabel='Temperature (K)'
      log_global_av_plot_filename="global_av_gsm_log.png"
      log_global_av_plot_title='Log Global average GSM'
      log_global_av_plot_ylabel='Temperature (K)'
      global_diff_plot_filename="global_diff_gsm.png"
      global_diff_plot_title='Residuals global - beam-weighted average GSM'
      global_diff_plot_ylabel='Temperature difference (K)'
      log_global_diff_plot_filename="global_diff_gsm_log.png"
      log_global_diff_plot_title='Log Residuals global - beam-weighted average GSM'
      log_global_diff_plot_ylabel='Temperature difference Log (K)'
      fixed_beam_diff_plot_filename="fixed_beam_diff_gsm.png"
      fixed_beam_diff_plot_title='Residuals fixed beam - beam-weighted average GSM'
      fixed_beam_diff_plot_ylabel='Temperature difference (K)'
   if use_model_gsm2016:
      global_av_plot_filename="global_av_gsm2016.png"
      global_av_plot_title='Global average GSM2016'
      global_av_plot_ylabel='Temperature (K)'
      log_global_av_plot_filename="global_av_gsm2016_log.png"
      log_global_av_plot_title='Log Global average GSM2016'
      log_global_av_plot_ylabel='Temperature (K)'
      global_diff_plot_filename="global_diff_gsm2016.png"
      global_diff_plot_title='Residuals global - beam-weighted average GSM2016'
      global_diff_plot_ylabel='Temperature difference (K)'
      log_global_diff_plot_filename="global_diff_gsm2016_log.png"
      log_global_diff_plot_title='Log Residuals global - beam-weighted average GSM2016'
      log_global_diff_plot_ylabel='Temperature difference Log (K)'
      fixed_beam_diff_plot_filename="fixed_beam_diff_gsm2016.png"
      fixed_beam_diff_plot_title='Residuals fixed beam - beam-weighted average GSM2016'
      fixed_beam_diff_plot_ylabel='Temperature difference (K)'
   if use_model_lfsm:  
      global_av_plot_filename="global_av_lfsm.png"
      global_av_plot_title='Global average LFSM'
      global_av_plot_ylabel='Temperature (K)'
      log_global_av_plot_filename="global_av_lfsm_log.png"
      log_global_av_plot_title='Log Global average LFSM'
      log_global_av_plot_ylabel='Temperature (K)'
      global_diff_plot_filename="global_diff_lfsm.png"
      global_diff_plot_title='Residuals global - beam-weighted average LFSM'
      global_diff_plot_ylabel='Temperature difference (K)'
      log_global_diff_plot_filename="global_diff_lfsm_log.png"
      log_global_diff_plot_title='Log Residuals global - beam-weighted average LFSM'
      log_global_diff_plot_ylabel='Temperature difference Log (K)'
      fixed_beam_diff_plot_filename="fixed_beam_diff_lfsm.png"
      fixed_beam_diff_plot_title='Residuals fixed beam - beam-weighted average LFSM'
      fixed_beam_diff_plot_ylabel='Temperature difference (K)'
   
   plt.clf()
   global_diff_plot=plt.figure()
   #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
   plt.plot(big_freq_array,global_mean_difference,label="Residuals")
   plt.title(global_diff_plot_title)
   plt.ylabel(global_diff_plot_ylabel)
   plt.xlabel('Frequency (MHz)')
   plt.text(180, -700, 'Residual RMS = %5.2f ' % (global_diff_residual_rms))
   plt.legend(loc=1)
   global_diff_plot.savefig(global_diff_plot_filename)
   print "saved figure %s" %  global_diff_plot_filename

   plt.clf()
   global_diff_plot=plt.figure()
   #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
   plt.loglog(big_freq_array,global_mean_difference,label="Residuals")
   plt.title(log_global_diff_plot_title)
   plt.ylabel(log_global_diff_plot_ylabel)
   plt.xlabel('Frequency (MHz)')
   plt.text(180, 10, 'Residual RMS = %5.2f ' % (global_diff_residual_rms))
   plt.legend(loc=1)
   global_diff_plot.savefig(log_global_diff_plot_filename)
   print "saved figure %s" %  log_global_diff_plot_filename
   
   #Plot the global average with no pb stuff
   #and a powerlaw fit
   global_average_gsm_temp_data_error=10
   powerlaw_fit_result = fit_powerlaw(global_average_gsm_temp_data,global_average_gsm_temp_data_error,big_freq_array)
   
   plt.clf()
   global_av_plot=plt.figure()
   #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
   plt.plot(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed') 
   plt.plot(big_freq_array,global_average_gsm_temp_data,label="Global average")
   plt.text(120, 1500, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   plt.text(120, 1300, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   plt.title(global_av_plot_title)
   plt.ylabel(global_av_plot_ylabel)
   plt.xlabel('Frequency (MHz)')
   #plt.text(180, -700, 'Residual RMS = %5.2f ' % (global_av_residual_rms))
   plt.legend(loc=1)
   global_av_plot.savefig(global_av_plot_filename)
   print "saved figure %s" %  global_av_plot_filename

   plt.clf()
   global_av_plot=plt.figure()
   #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
   plt.loglog(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed') 
   plt.loglog(big_freq_array,global_average_gsm_temp_data,label="Global average")
   plt.title(log_global_av_plot_title)
   plt.ylabel(log_global_av_plot_ylabel)
   plt.xlabel('Frequency (MHz)')
   #plt.text(180, 10, 'Residual RMS = %5.2f ' % (global_av_residual_rms))
   plt.legend(loc=1)
   global_av_plot.savefig(log_global_av_plot_filename)
   print "saved figure %s" %  log_global_av_plot_filename
   
   ################################################
   ################################################
   #Fit power law to difference plot, work out residuals
   
   logx=np.log10(big_freq_array/150.0)
   logy = np.log10(global_mean_difference)

   global_diff_error=10
   logyerr = global_diff_error / global_mean_difference
   
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
   
   alpha_global_diff = pfinal[1]
   T_150_global_diff = 10.0**pfinal[0]
   
   alpha_global_diff_err = np.sqrt( covar[1][1] )
   T_150_global_diff_err = np.sqrt( covar[0][0] ) * T_150_global_diff
   
   T_global_diff_fit = powerlaw2(big_freq_array/150.0, T_150_global_diff, alpha_global_diff)
   
   #work out the residuals - what is the rms?
   residuals = global_mean_difference - T_global_diff_fit
   #print len(residuals)
   residual_rms=np.sqrt(np.mean(np.square(residuals)))
   #print residual_rms
   ################
   
   plt.clf()
   Tb_plot=plt.figure()
   plt.errorbar(big_freq_array,global_mean_difference,yerr=global_diff_error,label="Global diff ")
   plt.plot(big_freq_array,T_global_diff_fit,label="Powerlaw fit")
   plt.title('Global Diff vs frequency for MWA')
   plt.ylabel('Global diff temperature (Tb in K)')
   plt.xlabel('Frequency (MHz)')
   plt.text(150, 500, 'Temp_150MHz = %5.2f +/- %5.3f' % (T_150_global_diff, T_150_global_diff_err))
   plt.text(150, 400, 'Index = %5.2f +/- %5.3f' % (alpha_global_diff, alpha_global_diff_err))
   plt.legend(loc=1)
   Tb_plot.savefig('global_diff_fit_plot.png')
   
   plt.clf()
   #Plot subplot with log log space
   loglog_plot_filename='global_diff_Tb_plot_loglog.png'
   plot_title='Global diff temperature vs frequency for MWA'
   plot_ylabel='Global diff \n temperature (Tb in K)'
   fit_plot=plt.figure()
   plt.subplot(2, 1, 1)
   plt.plot(big_freq_array, T_global_diff_fit)     # Fit
   plt.errorbar(big_freq_array, global_mean_difference, yerr=global_diff_error, fmt='k.')  # Data
   plt.text(150, 1000, 'Temp_150MHz = %5.2f +/- %5.2f' % (T_150_global_diff, T_150_global_diff_err))
   plt.text(150, 700, 'Index = %5.2f +/- %5.2f' % (alpha_global_diff, alpha_global_diff_err))
   plt.title(plot_title)
   plt.xlabel('Frequency (MHz)')
   plt.ylabel(plot_ylabel)
   plt.xlim(70, 240)
   plt.tight_layout()
   
   plt.subplot(2, 1, 2)
   plt.loglog(big_freq_array, T_global_diff_fit)
   plt.errorbar(big_freq_array, global_mean_difference, yerr=global_diff_error, fmt='k.')  # Data
   plt.xlabel('Frequency (log scale)')
   plt.ylabel('Temp (log scale)')
   plt.xlim(70, 240)
   plt.tight_layout()
   
   fit_plot.savefig(loglog_plot_filename)
   print "saved figure %s" % loglog_plot_filename
   
   plt.clf()
   plot_name="global_diff_residual_plot.png"
   diff_plot=plt.figure()
   #plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
   plt.plot(big_freq_array,residuals,label="Residuals")
   plt.title('Global Diff Fit residuals  vs frequency for MWA')
   plt.ylabel('Residual temperature (K)')
   plt.xlabel('Frequency (MHz)')
   plt.text(70, 10, 'Residual RMS = %5.2f K' % (residual_rms))
   plt.legend(loc=2)
   diff_plot.savefig(plot_name)
   print "saved figure %s" %  plot_name
   
   
   ###Plot Global as two panel
   plot_filename='global_difference_comparison_2panel.png'
   plot_title='Beam-weighted average - global average temperature'
   plt.clf()
   #fit_plot=plt.figure()

   f,a = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0})

   #plt.subplot(2,1,1)
   a[0].errorbar(big_freq_array,global_mean_difference,yerr=global_diff_error*0,label="Global diff ")
   a[0].plot(big_freq_array,T_global_diff_fit,label="Power law fit",linestyle='dashed')
   #plt.title('Global Diff vs frequency for MWA')
   #plt.ylabel('Global diff temperature (Tb in K)')
   #plt.xlabel('Frequency (MHz)')
   a[0].text(150, 500, 'T_150 = %5.2f +/- %5.3f' % (T_150_global_diff, T_150_global_diff_err))
   a[0].text(150, 400, 'Alpha = %5.2f +/- %5.3f' % (alpha_global_diff, alpha_global_diff_err))
   a[0].legend(loc=1)
   #plt.tight_layout()
   
   #plt.subplot(2,1,2)
   a[1].plot(big_freq_array,residuals,label="Residuals")
   #plt.title('Global Diff Fit residuals vs frequency for MWA')
   #plt.ylabel('Residual temperature (K)')
   plt.xlabel('Frequency (MHz)')
   a[1].text(175, -5, 'Residual RMS = %5.2f K' % (residual_rms))
   a[1].legend(loc=1)
   #plt.tight_layout()
   
   set_shared_ylabel(a, 'Temperature (K)')
   
   f.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   ###
   
   ######
   ###Plot fixed beam as two panel
   fixed_beam_mean_difference_error=10
   #powerlaw_fit_result = fit_powerlaw(fixed_beam_mean_difference,fixed_beam_mean_difference_error,big_freq_array)
   
   #work out the residuals - what is the rms?
   residuals = fixed_beam_mean_difference - powerlaw_fit_result['T_sky_fit']
   #print len(residuals)
   residual_rms=np.sqrt(np.mean(np.square(residuals)))
   
   plot_filename='fixed_beam_difference_comparison_2panel.png'
   plot_title='Beam-weighted average - fixed beam average temperature'
   plt.clf()
   #fit_plot=plt.figure()

   f,a = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0})

   #plt.subplot(2,1,1)
   a[0].errorbar(big_freq_array,fixed_beam_mean_difference,yerr=fixed_beam_mean_difference_error*0,label="Fixed beam diff ")
   #a[0].plot(big_freq_array,powerlaw_fit_result['T_sky_fit'],label="Power law fit",linestyle='dashed')
   #plt.title('Global Diff vs frequency for MWA')
   #plt.ylabel('Global diff temperature (Tb in K)')
   #plt.xlabel('Frequency (MHz)')
   #a[0].text(150, 500, 'T_150 = %5.2f +/- %5.3f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   #a[0].text(150, 400, 'Alpha = %5.2f +/- %5.3f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha']))
   a[0].legend(loc=1)
   #plt.tight_layout()
   
   #plt.subplot(2,1,2)
   a[1].plot(big_freq_array,residuals,label="Residuals")
   #plt.title('Global Diff Fit residuals vs frequency for MWA')
   #plt.ylabel('Residual temperature (K)')
   plt.xlabel('Frequency (MHz)')
   a[1].text(175, -5, 'Residual RMS = %5.2f K' % (residual_rms))
   a[1].legend(loc=1)
   #plt.tight_layout()
   
   set_shared_ylabel(a, 'Temperature (K)')
   
   f.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   ###
   
   
   #################################################
   ##################################################
   
   ###############################
   #Now work out predicted values for stuff:
   ###Predicted values:
   ####predicted sky temp - now using model GSM
   ##Predicted values of mine are wrong because I don't include the horizon and 
   #set an ambient temperature for below the horizon like Marcin does
   #So for this section use the old GSM model (pretty much the same) and Marcin's values eg:
   #alpha_gal=-2.5
   #T_150=278.
   #Actually, Marcin's is still wrong cause the GSM is wrong! (LFSM is worse!)
   big_T_sky_predicted=T_150*((big_freq_array/150.0)**(alpha_gal))
   #print "Predicted Tsky is:"
   #print T_sky_predicted
   
   #predicted sky flux density for moon-disk area on sky
   big_Ssky_predicted_values=(2.0*k*big_T_sky_predicted*Omega)/((300.0/big_freq_array)**(2)*10.0**(-26))
   #print "Predicted Ssky is:"
   #print S_sky_predicted
   
   big_T_moon_predicted=np.zeros(len(big_freq_array)) + T_moon +  (T_refl_gal_150)*(big_freq_array/150.0)**(refl_gal_alpha)
   
   #predicted sky flux density for moon-disk area on sky
   big_Smoon_predicted_values=(2.0*k*big_T_moon_predicted*Omega)/((300.0/big_freq_array)**(2)*10.0**(-26))
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
   plt.ylabel('Mean Reflected RFI Flux Density (Jy)')
   plt.xlabel('Frequency (MHz)')
   plt.legend(loc=1)
   big_Smoon_plot.savefig('big_Srfi.png')
   
   #Plot the RFI condensed scale
   plt.clf()
   big_Smoon_plot=plt.figure(5)
   plt.errorbar(big_freq_array,big_Srfi_average_stddev_spectrum[:,0],yerr=big_Smoon_average_stddev_spectrum[:,1],label="Measured")
   #plt.errorbar(big_freq_array,big_predicted_moon_sky_difference,label="Predicted")
   axes = plt.gca()
   axes.set_ylim([-2,10])
   plt.title('Reflected Moon RFI Flux Density vs Frequency for MWA')
   plt.ylabel('Mean Reflected RFI Flux Density (Jy)')
   plt.xlabel('Frequency (MHz)')
   plt.legend(loc=1)
   big_Smoon_plot.savefig('big_Srfi_zoom.png')


   

   ######################################################
   #Make plots for paper with all three models:
   #linear scale with fit:
   plot_filename='sky_model_comparison.png'
   plot_title='Beam-weighted average temperature vs frequency for MWA'
   plt.clf()
   
   plt.figure()
   #fit_plot=plt.figure()
   #https://stackoverflow.com/questions/32633322/changing-aspect-ratio-of-subplots-in-matplotlib
   #https://stackoverflow.com/questions/23528477/share-axes-in-matplotlib-for-only-part-of-the-subplots
   nrows, ncols = 3, 2
   dx, dy = 1.2, 1
   label_posx=110
   label1_posy=1350
   label2_posy=1150
   figsize = plt.figaspect(float(dy * nrows) / float(dx * ncols)) * 1.8
   fig,axes = plt.subplots(nrows,ncols, gridspec_kw={'hspace':0,'wspace':0.2}, figsize=figsize)

   average_temp_filename="average_gsm_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #plt.subplot(3, 1, 1)
   #axes.flat[0].get_shared_x_axes().join(axes.flat[0], axes.flat[4])
   axes.flat[0].plot(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[0].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0, label='GSM')  # Data
   axes.flat[0].text(label_posx, label1_posy, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   axes.flat[0].text(label_posx, label2_posy+50, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.title(plot_title)
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[0].legend(loc=1)
   #plt.tight_layout()
     
   average_temp_filename="average_gsm2016_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #plt.subplot(3, 1, 2)
   axes.flat[2].get_shared_x_axes().join(axes.flat[2], axes.flat[0])
   axes.flat[2].plot(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[2].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0,label='GSM2017')  # Data
   axes.flat[2].text(label_posx, label1_posy, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   axes.flat[2].text(label_posx, label2_posy+50, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[2].legend(loc=1)
   #plt.tight_layout()


   average_temp_filename="average_lfsm_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #plt.subplot(3, 1, 3)
   axes.flat[4].get_shared_x_axes().join(axes.flat[4], axes.flat[0])
   axes.flat[4].plot(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[4].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0,label='LFSM')  # Data
   axes.flat[4].text(label_posx, label1_posy, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   axes.flat[4].text(label_posx, label2_posy, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.xlabel('Frequency (MHz)')
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[4].legend(loc=1)
   #plt.tight_layout()
   set_shared_ylabel([axes.flat[0],axes.flat[2],axes.flat[4]], plot_ylabel)
   #f.savefig(plot_filename)
   #print "saved figure %s" % plot_filename
   
   #LOG scale with fit and no errors shown:
   #plot_filename='sky_model_comparison_log.png'
   #plot_title='Beam-weighted average temperature vs frequency for MWA'
   #plt.clf()
   #fit_plot=plt.figure()
   #f,a = plt.subplots(3, sharex=True, gridspec_kw={'hspace':0})
   
   
   average_temp_filename="average_gsm_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #ax1=plt.subplot(311)
   #plt.setp(ax1.get_xticklabels(), visible=False)
   #axes.flat[1].get_shared_x_axes().join(axes.flat[1], axes.flat[5])
   axes.flat[1].loglog(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[1].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0,label='GSM')  # Data
   #plt.text(150, 1000, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   #plt.text(150, 600, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.title(plot_title)
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[1].legend(loc=1)
   #plt.tight_layout()
     
   average_temp_filename="average_gsm2016_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #ax2=plt.subplot(312, sharex = ax1)
   #plt.setp(ax2.get_xticklabels(), visible=False)
   axes.flat[3].get_shared_x_axes().join(axes.flat[3], axes.flat[1])
   axes.flat[3].loglog(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[3].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0,label='GSM2017')  # Data
   #plt.text(150, 1000, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   #plt.text(150, 600, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[3].legend(loc=1)
   #plt.tight_layout()


   average_temp_filename="average_lfsm_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   powerlaw_fit_result = fit_powerlaw(average_temp_data,average_temp_data_error,big_freq_array)
   plot_ylabel='Sky Temp (K)'
   #ax3=plt.subplot(313, sharex = ax1)
   axes.flat[5].get_shared_x_axes().join(axes.flat[5], axes.flat[1])
   axes.flat[5].loglog(big_freq_array, powerlaw_fit_result['T_sky_fit'],label='fit',linestyle='dashed')     # Fit
   axes.flat[5].errorbar(big_freq_array, average_temp_data, yerr=average_temp_data_error*0,label='LFSM')  # Data
   #plt.text(150, 1000, 'T_150 = %5.2f +/- %5.2f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   #plt.text(150, 600, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   #plt.xlabel('Frequency (MHz)')
   #plt.ylabel(plot_ylabel)
   #plt.xlim(70, 240)
   axes.flat[5].legend(loc=1)
   #plt.tight_layout()
   
   #fit_plot.savefig(plot_filename)
   #print "saved figure %s" % plot_filename
   #set_shared_ylabel(a, plot_ylabel)
   fig.text(0.5, 0.04, 'Frequency (MHz)', ha='center')
   fig.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   #LOG
   #plt.subplot(2, 1, 2)
   #plt.loglog(big_freq_array, T_global_diff_fit)
   #plt.errorbar(big_freq_array, global_mean_difference, yerr=global_diff_error, fmt='k.')  # Data
   #plt.xlabel('Frequency (log scale)')
   #plt.ylabel('Temp (log scale)')
   #plt.xlim(70, 240)
   #plt.tight_layout()
   
   ####################################################################
   ######
   ###Plot corrected T_sky and correction factor as 2 panel
   ####################################
   ## T_sky corrected for chromatic beam effects
   ## T_meas = C(nu) * Tsky_corrected(nu)
   ## where C=beam_weighted_average_temp/beam_weighted_average_temp_at_150MHz
   
   #only do if Tb all positive
   if any(t < 0 for t in Tb):
      print "Tb has negative elements - skipping chromaticity correction"
   else:
      T_sky_corrected=Tb/chromaticity_correction
   
      T_sky_corrected_error=Tb_error/chromaticity_correction                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
      powerlaw_fit_result = fit_powerlaw(T_sky_corrected,T_sky_corrected_error,big_freq_array)
   
      plot_filename='chromaticity_correction_Tb_2panel.png'
      plot_title='T_sky corrected for beam chromaticity'
      plt.clf()
      #fit_plot=plt.figure()

      f,a = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0})

      #plt.subplot(2,1,1)
      T_sky_corrected_error=T_sky_corrected*0.15
      a[0].errorbar(big_freq_array,T_sky_corrected,yerr=T_sky_corrected_error*0,label="T_sky corrected ")
      a[0].plot(big_freq_array,powerlaw_fit_result['T_sky_fit'],label="Power law fit",linestyle='dashed')
      #plt.title('Global Diff vs frequency for MWA')
      #plt.ylabel('Global diff temperature (Tb in K)')
      #plt.xlabel('Frequency (MHz)')
      a[0].text(150, 1000, 'T_150 = %5.1f +/- %3.1f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
      a[0].text(150, 800, 'Alpha = %4.2f +/- %4.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
      a[0].legend(loc=1)
      a[0].set_ylabel("Temperature (K)")
      #plt.tight_layout()
   
      #plt.subplot(2,1,2)
      a[1].plot(big_freq_array,chromaticity_correction,label="chromaticity factor")
      #plt.title('Global Diff Fit residuals vs frequency for MWA')
      #plt.ylabel('Residual temperature (K)')
      plt.xlabel('Frequency (MHz)')
      #a[1].text(175, -5, 'Residual RMS = %5.2f K' % (residual_rms))
      a[1].legend(loc=1)
      a[1].set_ylabel("Chromaticity factor")
      #plt.tight_layout()
   
      #set_shared_ylabel(a, 'Temperature (K)')
   
      f.savefig(plot_filename)
      print "saved figure %s" % plot_filename
   
   
   ####################################################################
   ######
   ###Plot corrected SIMULATED (for IAU333 paper)  and correction factor as 2 panel
   ####################################
   ## T_sky corrected for chromatic beam effects
   ## T_meas = C(nu) * Tsky_corrected(nu)
   ## where C=beam_weighted_average_temp/beam_weighted_average_temp_at_150MHz
   average_temp_filename="average_gsm_temp_for_moon_region.npy"
   fractional_error=0.15
   average_temp_data=np.load(average_temp_filename)
   average_temp_data_error=average_temp_data*fractional_error
   
   sim_T_sky_corrected=average_temp_data/chromaticity_correction
   
   sim_T_sky_corrected_error=average_temp_data_error/chromaticity_correction                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
   powerlaw_fit_result = fit_powerlaw(sim_T_sky_corrected,sim_T_sky_corrected_error,big_freq_array)
   
   plot_filename='chromaticity_correction_Tb_2panel_SIM.png'
   plot_title='GSM sky temperature corrected for beam chromaticity'
   plt.clf()
   #fit_plot=plt.figure()

   f,a = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0})

   #plt.subplot(2,1,1)
   #T_sky_corrected_error=T_sky_corrected*0.15
   a[0].errorbar(big_freq_array,sim_T_sky_corrected,yerr=sim_T_sky_corrected_error*0,label="T_sky corrected ")
   a[0].plot(big_freq_array,powerlaw_fit_result['T_sky_fit'],label="Power law fit",linestyle='dashed')
   #plt.title('Global Diff vs frequency for MWA')
   a[0].set_ylabel('Sky temperature (K)')
   #plt.xlabel('Frequency (MHz)')
   a[0].text(150, 1000, 'T_150 = %5.1f +/- %5.1f' % (powerlaw_fit_result['T_150'], powerlaw_fit_result['T_150_err']))
   a[0].text(150, 820, 'Alpha = %5.2f +/- %5.2f' % (powerlaw_fit_result['alpha'], powerlaw_fit_result['alpha_err']))
   a[0].legend(loc=1)
   #plt.tight_layout()
   
   #plt.subplot(2,1,2)
   a[1].plot(big_freq_array,chromaticity_correction,label="chromaticity factor")
   #plt.title('Global Diff Fit residuals vs frequency for MWA')
   a[1].set_ylabel('Correction factor')
   plt.xlabel('Frequency (MHz)')
   #a[1].text(175, -5, 'Residual RMS = %5.2f K' % (residual_rms))
   a[1].legend(loc=1)
   #plt.tight_layout()
   
   #set_shared_ylabel(a, 'Temperature (K)')
   
   f.savefig(plot_filename)
   print "saved figure %s" % plot_filename
   
   
   
   
   ############################################################3
   
#bask in glory
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: model_moon.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--stokes',type='string', dest='stokes',default='I',help='Stokes parameter of Moon images. I,Q,U,V or linear (sqrt(Q^2+U^2)) e.g. --stokes="Q" [default=%default]')

(options, args) = parser.parse_args()

model_moon(options)
