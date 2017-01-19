#python code to model moon spectrum, remove rfi component ets.

import pyfits
import pylab as py
import numpy as np
from scipy import signal
from scipy import fftpack
import statsmodels.api as sm
from sklearn import linear_model


#multiple linear regresion (see http://stackoverflow.com/questions/11479064/multiple-linear-regression-in-python/14971531#14971531)
#using statsmodels
    
    
#read in moon images
#for each obsid
#for each channel
#xstart,xend,ystart,yend=2000,3120,2000,3120
xstart_moon,xend_moon,ystart_moon,yend_moon=2300,2820,2300,2820
xstart_psf,xend_psf,ystart_psf,yend_psf=2299,2819,2301,2821
moon_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001_dirty_applied-I.fits"
moon_hdulist = pyfits.open(moon_fitsname)
moon_data=moon_hdulist[0].data[0,0,:,:]
moon_header=moon_hdulist[0].header
moon_zoom=moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]
pix_size_deg = np.abs(float(moon_header['CDELT1']))
moon_radius_pix = np.round(0.25/pix_size_deg)


#read in off moon images
off_moon_fitsname="1127572080_cotter_20150929_moon_93_track_off_moon_paired_1127313592_peeled-0001_dirty_applied-I.fits"
off_moon_hdulist = pyfits.open(off_moon_fitsname)
off_moon_data=off_moon_hdulist[0].data[0,0,:,:]
off_moon_header=off_moon_hdulist[0].header
off_moon_zoom=off_moon_data[xstart_moon:xend_moon,ystart_moon:yend_moon]

#difference the images
moon_minus_sky=moon_zoom-off_moon_zoom


#read in psf images
psf_fitsname="1127313592_cotter_20150926_moon_93_trackmoon_peeled-0001-psf.fits"
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
#print X
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

print "S_moon is %s Jy and S_RFI is %s Jy" % (S_moon, S_RFI)

#Woohoo! all give basically the same answer and it is about what we expect, the Moon disk is close to zero (constant offset is negative so that is good) and the RF is about 29 Jy

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

#Create reconstructed moon dirty image with RFI removed.
reconstructed_moon=S_moon*G_image_shift

py.figure(5)
py.clf()
py.title("Reconstructed Moon Dirty Image")
py.imshow( ( reconstructed_moon ), cmap=py.cm.Greys,origin='lower')
py.colorbar()
py.savefig("moon_reconstructed.png")

#Show earthshine (RFI) alone:
RFI_alone=S_RFI*psf_zoom

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

#Make a residual map (vedantham et al 2015 eqn 18)
#theta_hat=[[S_moon],[S_RFI]]
#H_theta_hat=S_moon*vec_G + S_RFI*vec_PSF
#vec_R=vec_D - H_theta_hat
#R=np.reshape(vec_R, psf_zoom.shape, order='F')
R=moon_minus_sky-reconstructed_moon-RFI_alone

py.figure(9)
py.clf()
py.title("Residual Image")
py.imshow( ( R ), cmap=py.cm.Greys,origin='lower')
py.colorbar()
py.savefig("residual_image.png")

print psf_data.shape
print psf_zoom.shape
print moon_mask.shape


#bask in glory

