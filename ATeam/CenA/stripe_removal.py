#!/usr/bin/env python
#find and remove horizontal stripes in image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
import numpy as np
from astropy.io import fits
from scipy import fftpack

image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar.fits"
out_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar_filtered.fits"

crop_x1, crop_x2, crop_y1, crop_y2 = 1750, 2000, 3100, 3350


hdulist = fits.open(image_filename)
img = hdulist[0].data
#img = img[0,0,:,:]

#Add 0.5 to image so there are no negatives (be sure to undo this later)
img += 0.5
 

#needs to be in Uint 8 for Canny
#Will need to do some thresholding/clipping and or select a 'quiet' part of the image 
#to be able to get sufficient dynamic range to detect the lines
#img_max = np.max(img)
#print img_max
#img_min = np.min(img)
#print img_min
#img_norm = img/img_max

img_clipped = np.clip(img,a_min=0,a_max=1)

img_clipped_max = np.max(img_clipped)
img_clipped_norm = img_clipped/img_clipped_max

img_uint8 = (img_clipped_norm*255).astype(np.uint8)

#check conversion
test_filename = 'test.fits'
fits.writeto(test_filename,img_uint8,clobber=True)
print "saved %s" % test_filename

#select a small region to find lines (250 pix is about 1 degree on the image)
img_uint8_cropped = img_uint8[crop_y1:crop_y2,crop_x1:crop_x2]
test_filename = 'cropped.fits'
fits.writeto(test_filename,img_uint8_cropped,clobber=True)
print "saved %s" % test_filename

subtr_mean_image = np.full((img_uint8_cropped.shape[0],img_uint8_cropped.shape[1]),np.nan)

#fit a (low-order) polynomial to each column and subtract it.
column_pixel_numbers = np.arange(0,img_uint8_cropped.shape[0])
for column_index in column_pixel_numbers:
   column_data = img_uint8_cropped[:,int(column_index)]
   #find the mean and std dev and clip anything more than 3 sigma
   sigma = np.std(column_data)
   mean = np.mean(column_data)
   column_data[column_data > (mean+3*sigma)] = mean
   
   #fit polynomial 
   z = np.polyfit(column_pixel_numbers, column_data, 8)
   p = np.poly1d(z)
   fit_line = p(column_pixel_numbers)
   
   #figname="poly_fit_column_%04d.png" % column_index
   #plt.clf()
   #plt.plot(column_pixel_numbers,fit_line,label='fit')
   #plt.plot(column_pixel_numbers,column_data,label='data')
   #plt.xlabel("pixel number")
   #plt.ylabel("data value")
   #plt.savefig(figname)
   #print("saved %s" % figname)
   #plt.close()
   
   #subtract fit line
   new_column_data = column_data - fit_line

   subtr_mean_image[:,int(column_index)] = new_column_data

filename = 'subtr_mean_image.fits'
fits.writeto(filename,subtr_mean_image,clobber=True)
print "saved %s" % test_filename


# Take the fourier transform of the image.
F1 = fftpack.fft2(subtr_mean_image)

print F1

# Now shift the quadrants around so that low spatial frequencies are in
# the center of the 2D fourier transformed image.
F2 = fftpack.fftshift(F1)
# Calculate a 2D power spectrum
psd2D = np.abs(F2)**2
F2_real = F2.real
F2_imag = F2.imag

#print psd2D

fitsname = "fft_psd_polar_image_mean_subtr.fits"
fits.writeto(fitsname,psd2D,clobber=True)
print "saved %s" % fitsname

fitsname = "fft_real_polar_image_mean_subtr.fits"
fits.writeto(fitsname,F2_real,clobber=True)
print "saved %s" % fitsname

fitsname = "fft_image_polar_image_mean_subtr.fits"
fits.writeto(fitsname,F2_imag,clobber=True)
print "saved %s" % fitsname

#find the indices of the four highest values
max_val_1_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
psd2D[max_val_1_index] = 0
max_val_2_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
psd2D[max_val_2_index] = 0
max_val_3_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
psd2D[max_val_3_index] = 0
max_val_4_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
psd2D[max_val_4_index] = 0

fft_filter = np.full(psd2D.shape,0.)
fft_filter[max_val_1_index] = 1
fft_filter[max_val_2_index] = 1
fft_filter[max_val_3_index] = 1
fft_filter[max_val_4_index] = 1

fitsname = "fft_filter.fits"
fits.writeto(fitsname,fft_filter,clobber=True)
print "saved %s" % fitsname

#fft_shift the filter and fft it to see if the stripes match
fft_filter_shift = fftpack.ifftshift(fft_filter)

fitsname = "fft_filter_shift.fits"
fits.writeto(fitsname,fft_filter_shift,clobber=True)
print "saved %s" % fitsname

fft_filter_shift_inverse = fftpack.ifft2(fft_filter_shift)
fft_filter_shift_inverse_real = fft_filter_shift_inverse.real
fft_filter_shift_inverse_imag = fft_filter_shift_inverse.imag
fft_filter_shift_inverse_abs = np.abs(fft_filter_shift_inverse)

fitsname = "fft_filter_shift_inverse_real.fits"
fits.writeto(fitsname,fft_filter_shift_inverse_real,clobber=True)
print "saved %s" % fitsname

fitsname = "fft_filter_shift_inverse_imag.fits"
fits.writeto(fitsname,fft_filter_shift_inverse_imag,clobber=True)
print "saved %s" % fitsname

fitsname = "fft_filter_shift_inverse_abs.fits"
fits.writeto(fitsname,fft_filter_shift_inverse_abs,clobber=True)
print "saved %s" % fitsname



#Maybe we don't need all this hough transform stufff
#img = img.astype(np.float32)

#img = cv2.imread(image_filename)
#gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
#edges = cv2.Canny(gray, 75, 150)

#edges = cv2.Canny(img_uint8_cropped, 0, 20,apertureSize = 3,L2gradient=True)

#test_filename = 'edges.fits'
#fits.writeto(test_filename,edges,clobber=True)
#print "saved %s" % test_filename


#lines = cv2.HoughLinesP(edges, 1, np.pi/180, 10, maxLineGap=3)

#print lines 
#for line in lines:
#   x1, y1, x2, y2 = line[0]
#   cv2.line(img_uint8_cropped, (x1, y1), (x2, y2), (0, 0, 128), 1)
#cv2.imshow("linesEdges", edges)
#cv2.imshow("linesDetected", img_uint8_cropped)
#cv2.waitKey(10000)
#cv2.destroyAllWindows()

