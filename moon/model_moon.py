#python code to model moon spectrum, remove rfi component ets.

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

#set if you just want to plot already saved data
plot_only=False

#set if using already cropped images
use_cropped_images=False

#Set if using the small 500 px images from namorrodor
small_images=True
ionpeeled_images=True

plot_images=False
#speed of light
c=299792458.0
#Boltzmann's constant
k=1.380648*10**(-23)
#Solid angle subtended by the Moon in steradians
Omega=6.67*10.0**(-5)
T_moon=230 #K
#Moon RA 22:49:02 DEC -5:34:25 (2015-09-26) (RA 342.45 DEC -5.573)
Tsky_150=245 #K (Landecker and wielebinski 1970)
T_refl_gal_150=14.6 # temp of reflected galactic emisssion at 150 from reflection modeling (takes into account albedo of 0.07)
refl_gal_alpha=-2.50
#RFI model mask radius in arcmin - limits the extent of the gausssian reflection model
#rfi_model_mask_size=2.0
rfi_model_mask_size=(3.75/2)

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
    
#filename for big numpy arrays
big_Smoon_average_stddev_spectrum_npy_filename="big_Smoon_average_stddev_spectrum.npy"
big_Srfi_average_stddev_spectrum_npy_filename="big_Srfi_average_stddev_spectrum.npy"

big_Smoon_predicted_npy_filename="big_Smoon_predicted.npy"
big_Ssky_predicted_npy_filename="big_Ssky_predicted.npy"

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


#initialise the other big arrays
big_Smoon_average_stddev_spectrum=np.zeros([number_coarse_chans,2])
big_Srfi_average_stddev_spectrum=np.zeros([number_coarse_chans,2])
big_Ssky_predicted_values=np.zeros(number_coarse_chans)
big_Smoon_predicted_values=np.zeros(number_coarse_chans)
big_Tb_value_error=np.zeros([number_coarse_chans,2])


band_centre_chans=[69,93,121,145,169]
#band_centre_chans=[121]

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
      xstart_moon,xend_moon,ystart_moon,yend_moon=0,2048,0,2048
      xstart_psf,xend_psf,ystart_psf,yend_psf=0,2048,0,2048      
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
         moon_radius_pix = np.round(0.25/pix_size_deg)

         #print max value of moon image in Jy/beam to check against rfi + moon total flux
         print "Max moon is %s Jy/beam" % (np.max(moon_data))
         print "Max moon is %s Jy/beam" % (np.max(moon_zoom))
   
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
         
         #difference the images
         moon_minus_sky=moon_zoom-off_moon_zoom
         #convert to Jy per pix
         moon_minus_sky_jyppix=moon_minus_sky/n_pixels_per_beam

         #print out the rms of the difference image ( corner only )
         difference_rms=np.sqrt(np.mean(np.square(moon_minus_sky[0:150,0:150])))
         

         #If the rms is too high then just put a nan in the arrays and move on to next obsid
         if (difference_rms>rms_threshold or math.isnan(difference_rms)==True):
            print "rms of difference map %s is %s. Discarding difference image." % (moon_difference_fitsname,difference_rms)
            #place values in the arrays
            if (chan != 'MFS'):
               Smoon_spectrum_values[chan,obsid_index]=np.nan
               Srfi_spectrum_values[chan,obsid_index]=np.nan
            continue  
         else:
            print "rms of difference map %s is %s. Writing out difference image." % (moon_difference_fitsname,difference_rms)
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
         
         #For the broadened RFI, try a disk insted of a Gaussian 
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
         #Go back to pt source RFI for test
         vec_RFI=vec_PSF
         H2=[vec_G,vec_RFI]
         
         
         #using statsmodels
         X2 = np.array(H2).T
         X2_const=sm.add_constant(X2)
         
         # Now do just with numpy:
         beta_hat2 = np.linalg.lstsq(X2_const,vec_D)[0]

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
         #for non-point source version:
         #RFI_alone2=np.sum(S_RFI2*RFI_broadening)
         RFI_alone2=S_RFI2
         print "Initial total RFI flux density is %s Jy" % RFI_alone2

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
         #RFI_alone2=np.sum(S_RFI2*RFI_broadening)
         RFI_alone2=S_RFI2
         
         #Convert to total flux in Jy for Moon (assume for now RFI is a point source therefore flux in Jy/pix is same as total flux is Jy)
         S_moon2_tot_Jy=np.sum(S_moon2*moon_mask)
         
         print "Final Smoon, with RFI broadening: is %s Jy/pix and S_RFI is %s Jy/pix for moon obsid %s chan %s" % (S_moon2, S_RFI2,on_moon_obsid,chan_string)
         print "Final total Moon flux density, with RFI broadening: is %s Jy and S_RFI is %s Jy for moon obsid %s chan %s" % (S_moon2_tot_Jy, RFI_alone2,on_moon_obsid,chan_string)
         
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
   #print "Number of Smoon data points is %s" % Smoon_data_points
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

   #predicted sky temp
   T_sky_predicted=248.0*((freq_array/150.0)**(-2.5))
   #print "Predicted Tsky is:"
   #print T_sky_predicted

   #predicted sky flux density for moon-disk area on sky
   S_sky_predicted=(2.0*k*T_sky_predicted*Omega)/((300.0/freq_array)**(2)*10.0**(-26))
   #print "Predicted Ssky is:"
   #print S_sky_predicted
   
   #predicted sky flux density for moon-disk area on sky
   S_moon_predicted=(2.0*k*T_moon*Omega)/((300.0/freq_array)**(2)*10.0**(-26))
   #print "Predicted Smoon is:"
   #print S_moon_predicted
   
   
   
   #put in to the big arrays:
   if (centre_chan==69):
      big_Smoon_average_stddev_spectrum[57-57:57-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[57-57:57-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[57-57:57-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[57-57:57-57+24,1]=Srfi_average_stddev_spectrum[:,1]
      big_Ssky_predicted_values[57-57:57-57+24]=S_sky_predicted
      big_Smoon_predicted_values[57-57:57-57+24]=S_moon_predicted
   if (centre_chan==93):
      big_Smoon_average_stddev_spectrum[81-57:81-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[81-57:81-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[81-57:81-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[81-57:81-57+24,1]=Srfi_average_stddev_spectrum[:,1]
      big_Ssky_predicted_values[81-57:81-57+24]=S_sky_predicted
      big_Smoon_predicted_values[81-57:81-57+24]=S_moon_predicted
   if (centre_chan==121):
      big_Smoon_average_stddev_spectrum[109-57:109-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[109-57:109-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[109-57:109-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[109-57:109-57+24,1]=Srfi_average_stddev_spectrum[:,1]
      big_Ssky_predicted_values[109-57:109-57+24]=S_sky_predicted
      big_Smoon_predicted_values[109-57:109-57+24]=S_moon_predicted
   if (centre_chan==145):
      big_Smoon_average_stddev_spectrum[133-57:133-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[133-57:133-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[133-57:133-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[133-57:133-57+24,1]=Srfi_average_stddev_spectrum[:,1]
      big_Ssky_predicted_values[133-57:133-57+24]=S_sky_predicted
      big_Smoon_predicted_values[133-57:133-57+24]=S_moon_predicted
   if (centre_chan==169):
      big_Smoon_average_stddev_spectrum[157-57:157-57+24,0]=Smoon_average_stddev_spectrum[:,0]
      big_Smoon_average_stddev_spectrum[157-57:157-57+24,1]=Smoon_average_stddev_spectrum[:,1]
      big_Srfi_average_stddev_spectrum[157-57:157-57+24,0]=Srfi_average_stddev_spectrum[:,0]
      big_Srfi_average_stddev_spectrum[157-57:157-57+24,1]=Srfi_average_stddev_spectrum[:,1]
      big_Ssky_predicted_values[157-57:157-57+24]=S_sky_predicted
      big_Smoon_predicted_values[157-57:157-57+24]=S_moon_predicted

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
   
#predicted moon-sky difference
big_predicted_moon_sky_difference=big_Smoon_predicted_values-big_Ssky_predicted_values

plt.clf()
big_Smoon_plot=plt.figure(2)
plt.errorbar(big_freq_array,big_Smoon_average_stddev_spectrum[:,0],yerr=big_Smoon_average_stddev_spectrum[:,1],label="Measured")
plt.errorbar(big_freq_array,big_predicted_moon_sky_difference,label="Predicted")
plt.title('Moon Flux Density vs Frequency for MWA')
plt.ylabel('Mean Moon Flux Density (Jy)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
big_Smoon_plot.savefig('big_Smoon_and_predicted_moon_sky.png')

#Convert to brightness temperature units (Got to add in the Galactic average temp + plus the reflected galactic emission in the predicted signal (7% albedo?)
Tmoon_measured_average= (big_Smoon_average_stddev_spectrum[:,0]*(300.0/big_freq_array)**2*1e-26) / (2.0*k*Omega) + Tsky_150*(big_freq_array/150.0)**(-2.5)
Tmoon_measured_stddev= (big_Smoon_average_stddev_spectrum[:,1]*(300.0/big_freq_array)**2*1e-26) / (2.0*k*Omega)
predicted_Tmoon_big=np.zeros(len(big_freq_array)) + T_moon +  (T_refl_gal_150)*(big_freq_array/150.0)**(-2.5)

plt.clf()
big_Tmoon_plot=plt.figure(3)
plt.errorbar(big_freq_array,Tmoon_measured_average,yerr=Tmoon_measured_stddev,label="Measured")
plt.errorbar(big_freq_array,predicted_Tmoon_big,label="Predicted")
plt.title('Moon Brightness Temperature vs Frequency for MWA')
plt.ylabel('Mean Moon Brightness Temperature (K)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
big_Tmoon_plot.savefig('big_Tmoon_and_predicted_moon_sky_T.png')

plt.clf()
big_Tmoon_plot=plt.figure(4)
plt.errorbar(big_freq_array,Tmoon_measured_average,yerr=Tmoon_measured_stddev,label="Measured")
plt.errorbar(big_freq_array,predicted_Tmoon_big,label="Predicted")
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
plt.legend(loc=4)
big_Smoon_plot.savefig('big_Srfi.png')

#Now calculate the inferred background temperature in kilokelvin
#Tb = T_moon + 160.0*(big_freq_array/60.0)**(-2.24) - ((10.0**(-26))*(c**2)*big_Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(big_freq_array*10**6)**2)
Tb = T_moon + T_refl_gal_150*(big_freq_array/150.0)**(refl_gal_alpha) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,0])/(2*k*Omega*(big_freq_array*10**6)**2)
#Tb_error=abs(0.07*Tsky_150*(big_freq_array/150.0)**(-2.5) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
#Tb_error=abs(0.07*160*(big_freq_array/60.0)**(-2.24) - (10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(freq_array*10**6)**2))
Tb_error=abs((10.0**(-26)*c**2*big_Smoon_average_stddev_spectrum[:,1])/(2*k*Omega*(big_freq_array*10**6)**2))
#print "Tb in K is"
#print Tb
#print Tb_error

#predicted sky temp
T_sky_predicted=Tsky_150*((big_freq_array/150.0)**(-2.5))

#Save and plot the inferred and predicted background temp
tb_filename="inferred_sky_background_temp_and_error.npy"
big_Tb_value_error[:,0]=Tb
big_Tb_value_error[:,1]=Tb_error
np.save(tb_filename,big_Tb_value_error)

plt.clf()
Tb_plot=plt.figure(6)
plt.errorbar(big_freq_array,Tb,yerr=Tb_error,label="Measured")
plt.plot(big_freq_array,T_sky_predicted,label="Predicted")
plt.title('Inferred backgroud temperature vs frequency for MWA')
plt.ylabel('Inferred background temperature (Tb in K)')
plt.xlabel('Frequency (MHz)')
plt.legend(loc=4)
Tb_plot.savefig('big_inferred_Tb_plot.png')

if (plot_images):
   #Plot images:
   py.figure(1)
   py.clf()
   py.title("Dirty Moon Image")
   py.imshow( ( moon_zoom ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("moon_image.png")
   
   py.figure(2)
   py.clf()
   py.title("Moon-Sky Difference Image")
   py.imshow( ( moon_minus_sky ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("moon_minus_sky.png")

   py.figure(3)
   py.clf()
   py.title("Moon Mask")
   py.imshow( ( moon_mask ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_mask.png")

   py.figure(4)
   py.clf()
   py.title("Moon Mask Convolved with PSF")
   py.imshow( ( G_image_shift ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_mask_psf_convol.png")


   py.figure(5)
   py.clf()
   py.title("Reconstructed Moon Dirty Image")
   py.imshow( ( reconstructed_moon ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_reconstructed.png")


   py.figure(6)
   py.clf()
   py.title("Reconstructed RFI-only Dirty Image")
   py.imshow( ( RFI_alone ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_only.png")

   py.figure(7)
   py.clf()
   py.title("Dirty Off-Moon Image")
   py.imshow( ( off_moon_zoom ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("off_moon_image.png")

   py.figure(8)
   py.clf()
   py.title("PSF Image")
   py.imshow( ( psf_zoom ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("PSF_image.png")


   py.figure(9)
   py.clf()
   py.title("Residual Image")
   py.imshow( ( R ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("residual_image.png")

   #broadened RFI
   py.figure(10)
   py.clf()
   py.title("RFI Broadening Image")
   py.imshow( ( RFI_broadening ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_broadening_image.png")

   py.figure(11)
   py.clf()
   py.title("RFI Broadened Convolved Image")
   py.imshow( ( RFI_convolved_shift ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_broadened_convol_image.png")

   py.figure(12)
   py.clf()
   py.title("Broadened Reconstructed RFI-only Dirty Image")
   py.imshow( ( RFI_alone2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_only_broadened.png")

   py.figure(13)
   py.clf()
   py.title("Residual Image with Broadened RFI")
   py.imshow( ( R2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("residual_image_broadened_rfi.png")

   py.figure(14)
   py.clf()
   py.title("Reconstructed Moon Dirty RFI Broadened Image")
   py.imshow( ( reconstructed_moon2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_reconstructed_rfi_broadened.png")

#bask in glory

