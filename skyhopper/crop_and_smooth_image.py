#script to crop and smooth large image to skyhopper FoV and resolution

from astropy.io import fits
from astropy import wcs
import numpy as np
import scipy.ndimage as ndimage
from astropy.convolution import convolve,convolve_fft, Gaussian2DKernel, Tophat2DKernel
from astropy.modeling.models import Gaussian2D
from scipy.misc import imresize


final_image_size=512
#vista pix res = .0056 arcmin/pix
#4'' PSF = 11.9 = 12 pix

#Yfilter:
#large_image_name="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Y/New/vista_ngp_Y_new.fits"
#cropped_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Y/New/vista_ngp_Y_new_cropped.fits"
#cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Y/New/vista_ngp_Y_new_cropped_smoothed.fits"
#resized_cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Y/New/vista_ngp_Y_new_cropped_smoothed_resized.fits"

#Z filter (ngp field)
large_image_name="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Z/New/vista_ngp_Z_new.fits"
cropped_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Z/New/vista_ngp_Z_new_cropped.fits"
cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Z/New/vista_ngp_Z_new_cropped_smoothed.fits"
resized_cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/Z/New/vista_ngp_Z_new_cropped_smoothed_resized.fits"


#J:
#large_image_name="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/J/New/vista_ngp_J_new.fits"
#cropped_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/J/New/vista_ngp_J_new_cropped.fits"
#cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/J/New/vista_ngp_J_new_cropped_smoothed.fits"
#resized_cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/J/New/vista_ngp_J_new_cropped_smoothed_resized.fits"

#H:
#large_image_name="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/H/New/vista_ngp_H_new.fits"
#cropped_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/H/New/vista_ngp_H_new_cropped.fits"
#cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/H/New/vista_ngp_H_new_cropped_smoothed.fits"
#resized_cropped_smoothed_image_outname="/Users/benjaminmckinley/Documents/Astronomy/skyhopper/VISTA_images/H/New/vista_ngp_H_new_cropped_smoothed_resized.fits"


f = fits.open(large_image_name)
w = wcs.WCS(f[0].header)
newf = fits.PrimaryHDU()

#imsize is 2680 pix (.0056 arcmin per pix)

#Y:
#image=f[0].data[8950:11630,3550:6230]
#Z:
image=f[0].data[8950:11630,3550:6230]
#J:
#image=f[0].data[8950:11630,3550:6230]
#H:
#image=f[0].data[8950:11630,3550:6230]


newf.data = image
newf.header = f[0].header
newf.header.update(w[8950:11630,3550:6230].to_header())

fits.writeto(cropped_image_outname,newf.data,clobber=True)
fits.update(cropped_image_outname, newf.data, newf.header)

#smooth:
smoothed_image_gauss = ndimage.gaussian_filter(image, sigma=11.9, order=0)
##astropy - tapers edge of image ... otherwise same as above - use ndimage
#gauss_kernel = Gaussian2DKernel(12)
#smoothed_image_gauss = convolve_fft(image, gauss_kernel)

fits.writeto(cropped_smoothed_image_outname,smoothed_image_gauss,clobber=True)
fits.update(cropped_smoothed_image_outname, smoothed_image_gauss, newf.header)

#degrade the resolution of the images to 512x512 pix
resized_image=imresize(smoothed_image_gauss, [512, 512])

fits.writeto(resized_cropped_smoothed_image_outname,resized_image,clobber=True)



