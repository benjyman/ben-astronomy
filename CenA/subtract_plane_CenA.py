from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pyfits

#surface_fits_name='cenA_SPASS_J2000-I_bkg.fits'
#Read in the surface fits file produce by BANE.py
#surface_fits_name='cenA_SPASS_i_orig_bkg.fits'
#orig_image_fits_name='cenA_SPASS_i_orig.fits'
#subtr_output_image_name='cenA_SPASS_i_plane_subtr.fits'

surface_fits_name='CenA_Feain_1400_smooth_to_SPASS_first_then_regrid_SPASS_bkg.fits'
orig_image_fits_name='CenA_Feain_1400_smooth_to_SPASS_first_then_regrid_SPASS.fits'
subtr_output_image_name='CenA_Feain_1400_smooth_to_SPASS_first_then_regrid_SPASS_plane_subtr.fits'



hdulist = pyfits.open(surface_fits_name)
image_data=hdulist[0].data
image_header=hdulist[0].header
print image_data.shape


#m = 256 #size of the matrix
m_naxis_1=image_data.shape[0]
m_naxis_2=image_data.shape[1]


X1, X2 = np.mgrid[:m_naxis_1, :m_naxis_2]

fig = plt.figure()
ax = fig.add_subplot(3,1,1, projection='3d')
jet = plt.get_cmap('jet')



##generation of the surface
#F = 3        
#i = np.minimum(X1, m-X1-1)
#j = np.minimum(X2, m-X2-1)
#H = np.exp(-.5*(np.power(i, 2)  +  np.power(j, 2)   )/(F*F))
#Y = np.real(  np.fft.ifft2   (H  *  np.fft.fft2(  np.random.randn(m, m))))
#a = 0.0005; b = 0.0002; #parameters of the tilted plane
#Y = Y + (a*X1 + b*X2); #adding the plane
#Y = (Y - np.min(Y)) / (np.max(Y) - np.min(Y)) #data scaling

Y=image_data

#plot the initial topological surface
ax.plot_surface(X1,X2,Y, rstride = 1, cstride = 1, cmap = jet, linewidth = 0)

#plt.show()

#Regression
X = np.hstack(   ( np.reshape(X1, (m_naxis_1*m_naxis_2, 1)) , np.reshape(X2, (m_naxis_1*m_naxis_2, 1)) ) )
X = np.hstack(   ( np.ones((m_naxis_1*m_naxis_2, 1)) , X ))
YY = np.reshape(Y, (m_naxis_1*m_naxis_2, 1))

theta = np.dot(np.dot( np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()), YY)

plane = np.reshape(np.dot(X, theta), (m_naxis_1, m_naxis_2));

ax = fig.add_subplot(3,1,2, projection='3d')
ax.plot_surface(X1,X2,plane)
ax.plot_surface(X1,X2,Y, rstride = 1, cstride = 1, cmap = jet, linewidth = 0)


#Subtraction
Y_sub = Y - plane
ax = fig.add_subplot(3,1,3, projection='3d')
ax.plot_surface(X1,X2,Y_sub, rstride = 1, cstride = 1, cmap = jet, linewidth = 0)

#read in the original (pre-BANE.py) image:
hdulist_orig = pyfits.open(orig_image_fits_name)
image_data_orig=hdulist_orig[0].data
hdulist_orig_header=hdulist_orig[0].header



orig_image_sub=image_data_orig - plane

#create a new fits file to store the output 
pyfits.writeto(subtr_output_image_name,orig_image_sub)

pyfits.update(subtr_output_image_name, orig_image_sub, hdulist_orig_header)

hdulist_orig.close()
hdulist.close()
plt.show()




