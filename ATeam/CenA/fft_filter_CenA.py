from scipy import fftpack
import pyfits
import numpy as np
import pylab as py

#use the original image not the one that has been regridded to J2000 - NaNs!
#CenA_file_name="cenA_SPASS_i_orig.fits"
CenA_file_name="FFT/cropped/cenA_SPASS_J2000_regridMWA_cropped_new-I.fits"

hdulist=pyfits.open(CenA_file_name)
image = hdulist[0].data
orig_header=hdulist[0].header
image=np.require(image, dtype=np.float32)
print image.shape
# Take the fourier transform of the image.
F1 = fftpack.fft2(image)

# Now shift the quadrants around so that low spatial frequencies are in
# the center of the 2D fourier transformed image.
F2 = fftpack.fftshift( F1 )

# Calculate a 2D power spectrum
psd2D = np.abs( F2 )**2


#make a high pass filter
#smallest baseline approx 7.5m, smallest u = D/lambda = 7.5/2 = 3.75
#extent of SPASS image=30 deg (in x direction) so Lmax=sin(30/2) * 2 = .258 *2 =0.517
#so SPASS u res=1/0.517 = 1.93 
#So we want to filter out 2 (or 3.75/1.93) pixels either side of zero for u
#similarly for v: L-max is sin(20/2)*2 , v_resolution is 2.88 so we want to blank 3.75/2.88 = 1.3 pixels ie just one

#So make a filter 601*401 pixels which is 1 for all values except a 5x3 grid in the centre (and don't forget to FFT shift before inverse FFT)
high_pass_filter1=np.ones([865,865])
high_pass_filter1[432,432]=0

high_pass_filter2=np.zeros([401,601]) + 1
high_pass_filter2[199,300]=0
high_pass_filter2[200,300]=0
high_pass_filter2[201,300]=0
high_pass_filter2[199,299]=0
high_pass_filter2[200,299]=0
high_pass_filter2[201,299]=0
high_pass_filter2[199,298]=0
high_pass_filter2[200,298]=0
high_pass_filter2[201,298]=0
high_pass_filter2[199,301]=0
high_pass_filter2[200,301]=0
high_pass_filter2[201,301]=0
high_pass_filter2[199,302]=0
high_pass_filter2[200,302]=0
high_pass_filter2[201,302]=0

high_pass_filter3=np.zeros([401,601]) + 1
x_hanning=np.hanning(10)
y_hanning=np.hanning(6)
hanning_2D=np.sqrt(np.outer(y_hanning,x_hanning))
#high_pass_filter3[197:203,295:305]=high_pass_filter3[197:203,295:305]*(1-hanning_2D)
#high_pass_filter3[198:204,295:305]=high_pass_filter3[198:204,295:305]*(1-hanning_2D[0:6,0:10])
#high_pass_filter3[199:203,296:304]=high_pass_filter3[199:203,296:304]*(1-hanning_2D[1:5,1:9])
high_pass_filter3[199:202,298:303]=high_pass_filter3[199:202,298:303]*(1-hanning_2D[1:4,3:8])
#print high_pass_filter3[199:202,298:303].shape
#print hanning_2D

#filtered_FFT3=F2*high_pass_filter3
#filtered_FFT_psd3= np.abs( filtered_FFT3 )**2
##shift back the filtered FFT, and inverse FFT
#filtered_FFT_shiftback3=fftpack.ifftshift( filtered_FFT3 )
#filtered_image3 = np.real(fftpack.ifft2(filtered_FFT_shiftback3))

#filtered_FFT2=F2*high_pass_filter2
#filtered_FFT_psd2= np.abs( filtered_FFT2 )**2
##shift back the filtered FFT, and inverse FFT
#filtered_FFT_shiftback2=fftpack.ifftshift( filtered_FFT2 )
#filtered_image2 = np.real(fftpack.ifft2(filtered_FFT_shiftback2))

filtered_FFT1=F2*high_pass_filter1
filtered_FFT_psd1= np.abs( filtered_FFT1 )**2
#shift back the filtered FFT, and inverse FFT
filtered_FFT_shiftback1=fftpack.ifftshift( filtered_FFT1 )
filtered_image1 = np.real(fftpack.ifft2(filtered_FFT_shiftback1))

# Now plot up  image, FFT, filter, filtered FFT, HPF image
py.figure(1)
py.clf()
py.title("Original Image")
py.imshow( ( image ), cmap=py.cm.Greys,vmin=0.0, vmax=2,origin='lower')
py.savefig("FFT/Original Image")
 
 
py.figure(2)
py.clf()
py.title("PSD of Original Image")
py.imshow( (  np.log10( psd2D )),origin='lower')
py.savefig("FFT/PSD of Original Image")

#py.figure(3)
#py.clf()
#py.title("Hanning High Pass Filter")
#py.imshow(high_pass_filter3,cmap=py.cm.Greys,origin='lower')
#py.savefig("FFT/Hanning High Pass Filter")

#py.figure(4)
#py.clf()
#py.title("Hanning High Pass Filtered Image")
#py.imshow( ( filtered_image3 ), cmap=py.cm.Greys,vmin=0.0, vmax=2,origin='lower')
#py.savefig("FFT/Hanning High Pass Filtered Image")
#hanning_filtered_fits_name="FFT/hanning_hpf_filtered_CenA_SPASS.fits"
#pyfits.writeto(hanning_filtered_fits_name,filtered_image3)
#pyfits.update(hanning_filtered_fits_name, filtered_image3, orig_header)


py.figure(5)
py.clf()
py.title("5x3 High Pass Filter")
py.imshow(high_pass_filter2,cmap=py.cm.Greys,origin='lower')
py.savefig("FFT/5x3 High Pass Filter")

#py.figure(6)
#py.clf()
#py.title("5x3 High Pass Filtered Image")
#py.imshow( ( filtered_image2 ), cmap=py.cm.Greys,vmin=0.0, vmax=2,origin='lower')
#py.savefig("FFT/5x3 High Pass Filtered Image")
#five_by_three_filtered_fits_name="FFT/five_by_three_hpf_filtered_CenA_SPASS.fits"
#pyfits.writeto(five_by_three_filtered_fits_name,filtered_image2)
#pyfits.update(five_by_three_filtered_fits_name, filtered_image2, orig_header)


#py.figure(7)
#py.clf()
#py.title("DC High Pass Filter")
#py.imshow(high_pass_filter1,cmap=py.cm.Greys,origin='lower')
#py.savefig("FFT/DC High Pass Filter")

py.figure(8)
py.clf()
py.title("DC High Pass Filtered Image")
py.imshow( ( filtered_image1 ), cmap=py.cm.Greys,vmin=0.0, vmax=2,origin='lower')
py.savefig("FFT/DC High Pass Filtered Image")
DC_filtered_fits_name="FFT/DC_hpf_filtered_CenA_SPASS_cropped_new.fits"
pyfits.writeto(DC_filtered_fits_name,filtered_image1)
pyfits.update(DC_filtered_fits_name, filtered_image1, orig_header)


hdulist.close()
#py.figure(6)
#py.clf()
#py.imshow( hanning_2D )

#py.show()