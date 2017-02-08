#python code to model moon spectrum, remove rfi component ets.

import pyfits
import pylab as py
import numpy as np
from scipy import signal
from scipy import fftpack
import statsmodels.api as sm
from sklearn import linear_model

plot_images=False

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
    
#read the info files to find out how many observations there are and what the on and off-moon obsids are:
on_moon_filename="20150926_moon_93_test.txt"
off_moon_filename="20150929_off_moon1_93_test.txt"

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
#xstart,xend,ystart,yend=2000,3120,2000,3120
xstart_moon,xend_moon,ystart_moon,yend_moon=2300,2820,2300,2820
xstart_psf,xend_psf,ystart_psf,yend_psf=2299,2819,2301,2821

#initialise arrays of length n_chans to store Smoon and Srfi values
Smoon_spectrum_values=np.empty([n_chans,n_obs])
Srfi_spectrum_values=np.empty([n_chans,n_obs])
Smoon_average_stddev_spectrum=np.empty([n_chans,2])
Srfi_average_stddev_spectrum=np.empty([n_chans,2])


#need to do all this stuff for each frequency channel and each obsid
for chan in range(n_chans):
   for obsid_index in range(n_obs):
      #sort out filenames 
      on_moon_obsid=on_moon_obsid_list[obsid_index]
      off_moon_obsid=off_moon_obsid_list[obsid_index] 
      chan_string='%.04d' % chan
      #moon_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001_dirty_applied-I.fits"
      moon_fitsname="images/%s_cotter_20150926_moon_93_trackmoon_peeled-%s_dirty_applied-I.fits" % (on_moon_obsid,chan_string)
      #off_moon_fitsname="1127572080_cotter_20150929_moon_93_track_off_moon_paired_1127313592_peeled-0001_dirty_applied-I.fits"
      off_moon_fitsname="images/%s_cotter_20150929_moon_93_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (off_moon_obsid,on_moon_obsid,chan_string)
      #psf_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001-psf.fits"
      psf_fitsname="images/%s_cotter_20150926_moon_93_trackmoon_peeled-%s-psf.fits" % (on_moon_obsid,chan_string)

      moon_hdulist = pyfits.open(moon_fitsname)
      moon_data=moon_hdulist[0].data[0,0,:,:]
      moon_header=moon_hdulist[0].header
      moon_zoom=moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]
      pix_size_deg = np.abs(float(moon_header['CDELT1']))
      moon_radius_pix = np.round(0.25/pix_size_deg)


      #read in off moon images
      off_moon_hdulist = pyfits.open(off_moon_fitsname)
      off_moon_data=off_moon_hdulist[0].data[0,0,:,:]
      off_moon_header=off_moon_hdulist[0].header
      off_moon_zoom=off_moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]

      #difference the images
      moon_minus_sky=moon_zoom-off_moon_zoom

      #read in psf images
      psf_hdulist = pyfits.open(psf_fitsname)
      psf_data=psf_hdulist[0].data[0,0,:,:]
      psf_header=psf_hdulist[0].header
      psf_zoom=psf_data[xstart_psf:xend_psf,ystart_psf:yend_psf]
      psf_zoom=np.require(psf_zoom, dtype=np.float32)

      image_length=moon_zoom.shape[0]
      image_height=moon_zoom.shape[1]
      moon_mask = np.zeros((image_length,image_height))

      #define the centre of the disk (moon position)
      a,b = (image_length/2)-1, (image_height/2)-1

      y,x = np.ogrid[-a:image_length-a, -b:image_height-b]
      mask = x*x + y*y <= moon_radius_pix*moon_radius_pix
      moon_mask[mask]=1


      #do some maths
      #following Harish's notation:
      #D = moon_zoom
      #M = moon_mask
      #P= psf_zoom
      #G= M convol P

      #Rather than a convolution, take the FFT of both images and multiply them, then inverse FFT
      moon_mask_fft=fftpack.fft2(moon_mask)
      psf_fft=fftpack.fft2(psf_zoom)
      G_fourier=moon_mask_fft*psf_fft
      G_image=np.real(fftpack.ifft2(G_fourier))
      G_image_shift=fftpack.fftshift(G_image)


      #G=signal.convolve2d(moon_mask,psf_zoom)

      #Vectorise
      vec_D=moon_minus_sky.flatten('F')
      vec_G=G_image_shift.flatten('F')
      vec_PSF=psf_zoom.flatten('F')

      H=[vec_G,vec_PSF]

      #using statsmodels
      X = np.array(H).T
      X_const=sm.add_constant(X)
      results = sm.OLS(endog=vec_D, exog=X_const).fit()
      #print results.summary()

      #Try again using sklearn
      reg = linear_model.LinearRegression(fit_intercept=True)
      reg.fit(X,vec_D)
      #print reg.coef_
      #print reg.intercept_

      # Now do just with numpy:
      beta_hat = np.linalg.lstsq(X_const,vec_D)[0]
      #print beta_hat

      #Use whichever one for now
      S_moon=beta_hat[1]
      S_RFI=beta_hat[2]

      #print "First test: S_moon is %s Jy and S_RFI is %s Jy" % (S_moon, S_RFI)

      #Woohoo! all give basically the same answer and it is about what we expect, the Moon disk is close to zero (constant offset is negative so that is good) and the RF is about 29 Jy

      #Create reconstructed moon dirty image with RFI removed.
      reconstructed_moon=S_moon*G_image_shift

      #Show earthshine (RFI) alone:
      RFI_alone=S_RFI*psf_zoom

      #Make a residual map (vedantham et al 2015 eqn 18)
      #theta_hat=[[S_moon],[S_RFI]]
      #H_theta_hat=S_moon*vec_G + S_RFI*vec_PSF
      #vec_R=vec_D - H_theta_hat
      #R=np.reshape(vec_R, psf_zoom.shape, order='F')
      #equivalent to just:
      R=moon_minus_sky-reconstructed_moon-RFI_alone

      ###############
      #What I've done isn't right - need to re-do the maths with this new convolution?
      #residuals show we are oversubtracting in the middle and under-subtracting 
      #at the edges so the RFI must be broader than the PSF (as expected)
      #So instead of using the psf * S_RFI as the model, use a Gaussian of width 3.86 arcmin (Vedantham et al 2015)
      #convolved with the PSF
      #gaussian_fwhm_deg=3.86/60.0
      gaussian_fwhm_deg=3.75/60.0
      gaussian_fwhm_pix = np.round(gaussian_fwhm_deg/pix_size_deg)
      RFI_broadening=makeGaussian(image_length, fwhm = gaussian_fwhm_pix, center=(a+1,b+1))
      #Convolve this with the PSF
      RFI_broadening_fft=fftpack.fft2(RFI_broadening)
      RFI_convolved_fourier=psf_fft*RFI_broadening_fft
      RFI_convolved_image=np.real(fftpack.ifft2(RFI_convolved_fourier))
      #shift and normalise - this is the new PSF for RFI
      RFI_convolved_shift=fftpack.fftshift(RFI_convolved_image)/np.max(RFI_convolved_image)

      vec_RFI=RFI_convolved_shift.flatten('F')

      H2=[vec_G,vec_RFI]

      #using statsmodels
      X2 = np.array(H2).T
      X2_const=sm.add_constant(X2)

      # Now do just with numpy:
      beta_hat2 = np.linalg.lstsq(X2_const,vec_D)[0]

      S_moon2=beta_hat2[1]
      S_RFI2=beta_hat2[2]
      RFI_alone2=S_RFI2*RFI_convolved_shift

      print "With RFI broadening: S_moon is %s Jy and S_RFI is %s Jy for moon obsid %s chan %s" % (S_moon2, S_RFI2,on_moon_obsid,chan_string)

      reconstructed_moon2=S_moon2*G_image_shift
      R2=moon_minus_sky-reconstructed_moon2-RFI_alone2

      #place values in the arrays
      Smoon_spectrum_values[chan,obsid_index]=S_moon2
      Srfi_spectrum_values[chan,obsid_index]=S_RFI2      

#work out average values and standard deviations
#Smoon_spectrum_values=np.empty([n_chans,n_obs])
#Srfi_spectrum_values=np.empty([n_chans,n_obs])
#Smoon_average_stddev_spectrum=np.empty([n_chans,2])
#Srfi_average_stddev_spectrum=np.empty([n_chans,2])
#print Smoon_spectrum_values
Smoon_average_stddev_spectrum[:,0]=np.mean(Smoon_spectrum_values, axis=1)
Smoon_average_stddev_spectrum[:,1]=np.std(Smoon_spectrum_values, axis=1)
print "S_moon average and std dev for each chan:" 
print Smoon_average_stddev_spectrum
Srfi_average_stddev_spectrum[:,0]=np.mean(Srfi_spectrum_values, axis=1)
Srfi_average_stddev_spectrum[:,1]=np.std(Srfi_spectrum_values, axis=1)
print "S_rfi average and std dev for each chan:"
print Srfi_average_stddev_spectrum

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

