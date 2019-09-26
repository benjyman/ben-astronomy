#!/usr/bin/env python
#find and remove horizontal stripes in image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
import numpy as np
from astropy.io import fits
from scipy import fftpack,ndimage
from PIL import Image
import sys

do_subtract_poly = False

# define normalized 2D gaussian
def gaus2d(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    super_power=4
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.))**super_power)
    
image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar.fits"
out_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar_filtered.fits"
polar_rotated_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb_polar_rotated.fits"

crop_x1, crop_x2, crop_y1, crop_y2 = 1750, 2000, 3100, 3350


hdulist = fits.open(image_filename)
img = hdulist[0].data
#img = img[0,0,:,:]

##Add 0.5 to image so there are no negatives (be sure to undo this later)
#img += 0.5
# 
#
##needs to be in Uint 8 for Canny
##Will need to do some thresholding/clipping and or select a 'quiet' part of the image 
##to be able to get sufficient dynamic range to detect the lines
##img_max = np.max(img)
##print img_max
##img_min = np.min(img)
##print img_min
##img_norm = img/img_max
#
#img_clipped = np.clip(img,a_min=0,a_max=1)
#
#img_clipped_max = np.max(img_clipped)
#img_clipped_norm = img_clipped/img_clipped_max
#
#img_uint8 = (img_clipped_norm*255).astype(np.uint8)
#
##check conversion
#test_filename = 'test.fits'
#fits.writeto(test_filename,img_uint8,clobber=True)
#print "saved %s" % test_filename
#
##select a small region to find lines (250 pix is about 1 degree on the image)
#img_uint8_cropped = img_uint8[crop_y1:crop_y2,crop_x1:crop_x2]
#
#
#subtr_mean_image = np.full((img_uint8_cropped.shape[0],img_uint8_cropped.shape[1]),np.nan)
#
##fit a (low-order) polynomial to each column and subtract it.
#column_pixel_numbers = np.arange(0,img_uint8_cropped.shape[0])
#for column_index in column_pixel_numbers:
#   column_data = img_uint8_cropped[:,int(column_index)]
#   #find the mean and std dev and clip anything more than 3 sigma
#   sigma = np.std(column_data)
#   mean = np.mean(column_data)
#   column_data[column_data > (mean+3*sigma)] = mean
#   
#   #fit polynomial 
#   z = np.polyfit(column_pixel_numbers, column_data, 8)
#   p = np.poly1d(z)
#   fit_line = p(column_pixel_numbers)
#   
#   #figname="poly_fit_column_%04d.png" % column_index
#   #plt.clf()
#   #plt.plot(column_pixel_numbers,fit_line,label='fit')
#   #plt.plot(column_pixel_numbers,column_data,label='data')
#   #plt.xlabel("pixel number")
#   #plt.ylabel("data value")
#   #plt.savefig(figname)
#   #print("saved %s" % figname)
#   #plt.close()
#   
#   #subtract fit line
#   new_column_data = column_data - fit_line
#
#   subtr_mean_image[:,int(column_index)] = new_column_data
#
#filename = 'subtr_mean_image.fits'
#fits.writeto(filename,subtr_mean_image,clobber=True)
#print "saved %s" % test_filename
#
#
## Take the fourier transform of the image.
#F1 = fftpack.fft2(subtr_mean_image)
#
#print F1
#
## Now shift the quadrants around so that low spatial frequencies are in
## the center of the 2D fourier transformed image.
#F2 = fftpack.fftshift(F1)
## Calculate a 2D power spectrum
#psd2D = np.abs(F2)**2
#F2_real = F2.real
#F2_imag = F2.imag
#
##print psd2D
#
#fitsname = "fft_psd_polar_image_mean_subtr.fits"
#fits.writeto(fitsname,psd2D,clobber=True)
#print "saved %s" % fitsname
#
#fitsname = "fft_real_polar_image_mean_subtr.fits"
#fits.writeto(fitsname,F2_real,clobber=True)
#print "saved %s" % fitsname
#
#fitsname = "fft_image_polar_image_mean_subtr.fits"
#fits.writeto(fitsname,F2_imag,clobber=True)
#print "saved %s" % fitsname
#
##find the indices of the four highest values
#max_val_1_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
#psd2D[max_val_1_index] = 0
#max_val_2_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
#psd2D[max_val_2_index] = 0
#max_val_3_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
#psd2D[max_val_3_index] = 0
#max_val_4_index = np.unravel_index(np.argmax(psd2D, axis=None), psd2D.shape)
#psd2D[max_val_4_index] = 0
#
#fft_filter = np.full(psd2D.shape,0.)
#fft_filter[max_val_1_index] = 1
#fft_filter[max_val_2_index] = 1
#fft_filter[max_val_3_index] = 1
#fft_filter[max_val_4_index] = 1
#
#fitsname = "fft_filter.fits"
#fits.writeto(fitsname,fft_filter,clobber=True)
#print "saved %s" % fitsname
#
##fft_shift the filter and fft it to see if the stripes match
#fft_filter_shift = fftpack.ifftshift(fft_filter)
#
#fitsname = "fft_filter_shift.fits"
#fits.writeto(fitsname,fft_filter_shift,clobber=True)
#print "saved %s" % fitsname
#
#fft_filter_shift_inverse = fftpack.ifft2(fft_filter_shift)
#fft_filter_shift_inverse_real = fft_filter_shift_inverse.real
#fft_filter_shift_inverse_imag = fft_filter_shift_inverse.imag
#fft_filter_shift_inverse_abs = np.abs(fft_filter_shift_inverse)
#
#fitsname = "fft_filter_shift_inverse_real.fits"
#fits.writeto(fitsname,fft_filter_shift_inverse_real,clobber=True)
#print "saved %s" % fitsname
#
#fitsname = "fft_filter_shift_inverse_imag.fits"
#fits.writeto(fitsname,fft_filter_shift_inverse_imag,clobber=True)
#print "saved %s" % fitsname
#
#fitsname = "fft_filter_shift_inverse_abs.fits"
#fits.writeto(fitsname,fft_filter_shift_inverse_abs,clobber=True)
#print "saved %s" % fitsname
#
#######################
#

fitsname = "img.fits"
fits.writeto(fitsname,img,clobber=True)
print "saved %s" % fitsname

img = img.astype(float)
img_FFT = fftpack.fft2(img)
img_FFT_shift = fftpack.fftshift(img_FFT)
img_FFT_shift_abs = np.abs(img_FFT_shift)

fitsname = "img_FFT_shift_abs.fits"
fits.writeto(fitsname,img_FFT_shift_abs,clobber=True)
print "saved %s" % fitsname

if do_subtract_poly:
   #fit a (low-order) polynomial to each column and subtract it.
   subtr_poly_image = np.full(img.shape,np.nan)
   column_pixel_numbers = np.arange(0,img.shape[0])
   for column_index in column_pixel_numbers:
      column_data = img[:,int(column_index)]
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
   
      subtr_poly_image[:,int(column_index)] = new_column_data
   
   filename = 'subtr_poly_full_image.fits'
   fits.writeto(filename,subtr_poly_image,clobber=True)
   print "saved %s" % filename
   
   
   # Take the fourier transform of the image.
   F1 = fftpack.fft2(subtr_poly_image)
   
   #print F1

   ## Now shift the quadrants around so that low spatial frequencies are in
   ## the center of the 2D fourier transformed image.
   F2 = fftpack.fftshift(F1)
   ## Calculate a 2D power spectrum
   psd2D = np.abs(F2)**2
   F2_real = F2.real
   F2_imag = F2.imag
   
   ##print psd2D
   
   filename = "fft_psd_polar_image_poly_subtr_full.fits"
   fits.writeto(filename,psd2D,clobber=True)
   print "saved %s" % filename



###messing around with filters
centre_pix_num = int(img.shape[0]/2.)
#By eye from FT of poly_subtr image, can clearly see where artifacts are
#make filter
filter_width = 200 #5
filter_middle_part =  1000 #252
filter_length = 4096 #700



#stoopid basic filter
#horiz_line_filter[centre_pix_num+filter_middle_part:centre_pix_num+filter_middle_part+filter_length,centre_pix_num-filter_width:centre_pix_num+filter_width] = 0.
#horiz_line_filter[centre_pix_num-filter_middle_part-filter_length:centre_pix_num-filter_middle_part,centre_pix_num-filter_width:centre_pix_num+filter_width] = 0.

#use 2d gaussian
#https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html
centre_pix = img.shape[0]/2 

x = np.arange(-centre_pix, centre_pix,1)
y = np.arange(-centre_pix, centre_pix,1)
#x = np.linspace(-5, 5)
#y = np.linspace(-5, 5)
xx, yy = np.meshgrid(x, y, sparse=True) 
centre_x = 0
centre_y = 0
horiz_line_filter_1 = gaus2d(x=xx, y=yy, mx=centre_x, my=centre_y, sx=filter_width, sy=filter_length)
horiz_line_filter_1 = horiz_line_filter_1/np.max(horiz_line_filter_1)

horiz_line_filter_2 = gaus2d(x=xx, y=yy, mx=centre_x, my=centre_y, sx=filter_width, sy=filter_middle_part)
horiz_line_filter_2 = horiz_line_filter_2/np.max(horiz_line_filter_2)


horiz_line_filter = horiz_line_filter_1 - horiz_line_filter_2

horiz_line_filter = 1. - horiz_line_filter

###
#horiz_line_filter = np.full(img.shape,1.)
#horiz_line_filter[:,centre_pix-filter_width/2:centre_pix+filter_width/2] = 0
###

#
horiz_line_filter = gaus2d(x=xx, y=yy, mx=0, my=0, sx=150, sy=200)
horiz_line_filter = horiz_line_filter / np.max(horiz_line_filter)

#


fitsname = "horiz_line_filter.fits"
fits.writeto(fitsname,horiz_line_filter,clobber=True)
print "saved %s" % fitsname


#invert the filter to have a look in image space
horiz_line_filter_invert = 1. - horiz_line_filter
horiz_line_filter_shift = fftpack.ifftshift(horiz_line_filter_invert)
horiz_line_filter_shift_inverse = fftpack.ifft2(horiz_line_filter_shift)
horiz_line_filter_shift_inverse_real = horiz_line_filter_shift_inverse.real

fitsname = "horiz_line_filter_shift_inverse_real.fits"
fits.writeto(fitsname,horiz_line_filter_shift_inverse_real,clobber=True)
print "saved %s" % fitsname


#Apply filter to image FFT
#filtered_image_FFT_shift = F2 * horiz_line_filter
filtered_image_FFT_shift = horiz_line_filter * img_FFT_shift
filtered_image_FFT_shift_abs = np.abs(filtered_image_FFT_shift)

fitsname = "filtered_image_FFT_shift_abs.fits"
fits.writeto(fitsname,filtered_image_FFT_shift_abs,clobber=True)
print "saved %s" % fitsname

#shift back
filtered_image_FFT_shift_back = fftpack.ifftshift(filtered_image_FFT_shift)

#FFT back
filtered_image_FFT_shift_back_inverse = fftpack.ifft2(filtered_image_FFT_shift_back).real

fitsname = "filtered_image_FFT_shift_back_inverse.fits"
fits.writeto(fitsname,filtered_image_FFT_shift_back_inverse,clobber=True)
print "saved %s" % fitsname


#regrid back!
value = np.sqrt(((img.shape[0]/2.0)**2.0)+((img.shape[1]/2.0)**2.0))
polar_image_undone = cv2.linearPolar(filtered_image_FFT_shift_back_inverse,(img.shape[0]/2, img.shape[1]/2), value, cv2.WARP_INVERSE_MAP)
polar_undone_filename = "filtered_image_polar_undone.fits"

fits.writeto(polar_undone_filename,polar_image_undone,clobber=True)
print "saved %s" % polar_undone_filename


#lowpass = ndimage.gaussian_filter(img, 5)
#gauss_highpass = img - lowpass

#fitsname = "lowpass.fits"
#fits.writeto(fitsname,lowpass,clobber=True)
#print "saved %s" % fitsname

#fitsname = "gauss_highpass.fits"
#fits.writeto(fitsname,gauss_highpass,clobber=True)
#print "saved %s" % fitsname

#print gauss_highpass.shape
#gauss_highpass_FFT =  fftpack.fft2(gauss_highpass)

#gauss_highpass_FFT_shift =  fftpack.fftshift(gauss_highpass_FFT)
#gauss_highpass_FFT_shift_abs = np.abs(gauss_highpass_FFT_shift)

#fitsname = "gauss_highpass_FFT_shift_abs.fits"
#fits.writeto(fitsname,gauss_highpass_FFT_shift_abs,clobber=True)
#print "saved %s" % fitsname


##FFT the original cropped image, find those four peaks and set them to be the average of the surrounding pix
#img_uint8_cropped_FFT = fftpack.fft2(img_uint8_cropped)
##img_uint8_cropped_FFT_inverse = fftpack.ifft2(img_uint8_cropped_FFT).real
##img_uint8_cropped_FFT_inverse_real = img_uint8_cropped_FFT_inverse.real
#
#
#img_uint8_cropped_FFT_shift =  fftpack.fftshift(img_uint8_cropped_FFT)
#
#
##test_filename = 'cropped.fits'
##fits.writeto(test_filename,img_uint8_cropped,clobber=True)
##print "saved %s" % test_filename
#
#
##img_uint8_cropped_FFT_shift_psd2D = np.abs(img_uint8_cropped_FFT_shift)**2
#
#
##fitsname = "img_uint8_cropped_FFT_shift_psd2D.fits"
##fits.writeto(fitsname,img_uint8_cropped_FFT_shift_psd2D,clobber=True)
##print "saved %s" % fitsname
##max_value_index_list = [max_val_1_index,max_val_2_index,max_val_3_index,max_val_4_index]
##for max_val_index in max_value_index_list:
##
##   #print img_uint8_cropped_FFT[max_val_index[0]-1:img_uint8_cropped_FFT[max_val_index[0]+2
##   sum_surrounding_values = np.sum(img_uint8_cropped_FFT[max_val_index[0]-1:max_val_index[0]+2,max_val_index[1]-1:max_val_index[1]+2] - img_uint8_cropped_FFT[max_val_index])
##   mean_surrounding_values = sum_surrounding_values/8.
##   img_uint8_cropped_FFT_shift[max_val_index] = 0
#
#img_uint8_cropped_FFT_shift[:,125] = 0
#img_uint8_cropped_FFT_shift[:,126] = 0
#img_uint8_cropped_FFT_shift[:,124] = 0
#
#img_uint8_cropped_FFT_shiftback =  fftpack.ifftshift(img_uint8_cropped_FFT_shift)
#
#img_uint8_cropped_FFT_shiftback_inverse = fftpack.ifft2(img_uint8_cropped_FFT_shiftback)
#img_uint8_cropped_FFT_shiftback_inverse_real = img_uint8_cropped_FFT_shiftback_inverse.real
#
#fitsname = "img_uint8_cropped_FFT_shiftback_inverse_real.fits"
#fits.writeto(fitsname,img_uint8_cropped_FFT_shiftback_inverse_real,clobber=True)
#print "saved %s" % fitsname
#
##Try with full image
#
##FFT the original cropped image, find those four peaks and set them to be the average of the surrounding pix
#img = img.astype(float)
#img_FFT = fftpack.fft2(img)
#
#img_FFT_shift =  fftpack.fftshift(img_FFT)
#
#print img_FFT_shift.shape
#
#img_FFT_shift[:,2048] = 0
##img_FFT_shift[:,2047] = 0
##img_FFT_shift[:,2049] = 0
#
#
#img_FFT_shiftback =  fftpack.ifftshift(img_FFT_shift)
#
#img_FFT_shiftback_inverse = fftpack.ifft2(img_FFT_shiftback)
#img_FFT_shiftback_inverse_real = img_FFT_shiftback_inverse.real
#
#fitsname = "img_FFT_shiftback_inverse_real.fits"
#fits.writeto(fitsname,img_FFT_shiftback_inverse_real,clobber=True)
#print "saved %s" % fitsname


##rotate the image so that stripes are horizontal
#
#pil_image = Image.fromarray(polar_image)
#rotated_image = np.asarray(pil_image.rotate(45))
#fits.writeto(polar_rotated_filename,rotated_image,clobber=True)
#print "saved %s" % polar_rotated_filename


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

