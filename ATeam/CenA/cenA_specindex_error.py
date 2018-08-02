
import matplotlib.pyplot as plt
import numpy as np
import pyfits
np.seterr(divide='ignore')

specindex_fits_name='CenA_spec_index_154_2307_DC_HPF_SPASS_cropped_new.fits'
S1_fits_name='1086439360_121_cotter_CenA_posfix_cropped_new_smooth_to_SPASS.fits'
S2_fits_name='DC_hpf_filtered_CenA_SPASS_cropped_new.fits'
error_fits_name='CenA_spec_index_154_2307_DC_HPF_SPASS_cropped_new_error.fits'

delta_S1=0.25
delta_S2=0.05

specindex_hdulist = pyfits.open(specindex_fits_name)
specindex_data=specindex_hdulist[0].data[0,0,:,:]
specindex_header=specindex_hdulist[0].header
print specindex_data.shape

S1_hdulist = pyfits.open(S1_fits_name)
S1_data=S1_hdulist[0].data[0,0,:,:]
S1_header=S1_hdulist[0].header
print S1_data.shape

S2_hdulist = pyfits.open(S2_fits_name)
S2_data=S2_hdulist[0].data
S2_header=S2_hdulist[0].header
print S2_data.shape


#Calculate spec index error:

delta_alpha=np.abs((1.0/np.log(154.0/2307.0)) * (S2_data/S1_data) * (specindex_data)) * np.sqrt((delta_S1/S1_data)**2 + (delta_S1/S1_data)**2)


#create a new fits file to store the output 
pyfits.writeto(error_fits_name,delta_alpha)
pyfits.update(error_fits_name, delta_alpha, S2_header)

specindex_hdulist.close()
S1_hdulist.close()
S2_hdulist.close()





