#!/usr/bin/env python
'''Add internal noise to eda2 sims
'''
import numpy as np
from astropy.io import fits

c = 299792458.
k = 1.38065e-23

def add_noise_coupling_to_sim_uvfits(band,daniel=True,uv_correlation_array_filename_x='uv_correlation_eda2_255_daniel_x.npy',uv_correlation_array_filename_y='uv_correlation_eda2_255_daniel_y.npy'):
   # get the values from the noise coupling array created using plot_internal_noise_coupling
   uv_correlation_array_x = np.load(uv_correlation_array_filename_x)
   uv_correlation_array_y = np.load(uv_correlation_array_filename_y)
   freq_index = int(band)
   start_freq_MHz = 50.0
   freq_MHz = freq_index * 1.28 + start_freq_MHz
   wavelength = 300. / freq_MHz 
   jy_to_K = (wavelength**2) / (2. * k * 1.0e26)
   print("jy_to_K %E" % jy_to_K)
   if daniel:
      type_list = ["gsm","gsm_uniform","EDGES_uniform","unity_uniform","angular"]
      pol_list = ['X','Y']
      LST_deg = 60.0
      for pol in pol_list:
         for type in type_list:   
            if type =="angular":
               uvfits_prepend = "data/woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular" % (LST_deg,start_freq_MHz,pol)
               input_uvfits_name = "%s_band%02d.uvfits" % (uvfits_prepend,band)
               output_uvfits_name = "%s_band%02d_nc.uvfits" % (uvfits_prepend,band)
            if type =="gsm_uniform":
               uvfits_prepend = "data/woden_LST_%0.3f_gsm_uniform_start_freq_%0.3f_pol_%s_global_foreground" % (LST_deg,start_freq_MHz,pol)            
               input_uvfits_name = "%s_band%02d.uvfits" % (uvfits_prepend,band)
               output_uvfits_name = "%s_band%02d_nc.uvfits" % (uvfits_prepend,band)  
            if pol=='X':
               if (type=="gsm" or type=="EDGES_uniform" or type=="unity_uniform"):
                  uvfits_prepend = "data/woden_LST_%0.3f_%s_start_freq_%0.3f" % (LST_deg,type,start_freq_MHz)
                  input_uvfits_name = "%s_band%02d.uvfits" % (uvfits_prepend,band)
                  output_uvfits_name = "%s_band%02d_nc.uvfits" % (uvfits_prepend,band)
   
            print("adding noise coupling to %s" % input_uvfits_name)
            
            with fits.open(input_uvfits_name) as hdulist:
               #hdulist.info()
               #info_string = [(x,x.data.shape,x.data.dtype.names) for x in hdulist]
               #print info_string
               uvtable = hdulist[0].data
               visibilities_single = uvtable['DATA']
               visibilities_shape = visibilities_single.shape
               print("visibilities_shape")
               print(visibilities_shape)
            
               UU_s_array = uvtable['UU']
               UU_m_array = UU_s_array * c   
               VV_s_array = uvtable['VV']
               VV_m_array = VV_s_array * c
            
               data = hdulist[0].data.data
               
               ####X pol 
               pol_index = 0
          
               internal_noise_real = uv_correlation_array_x[:,2+freq_index].real
               print(internal_noise_real[0:20])
               internal_noise_imag = uv_correlation_array_x[:,2+freq_index].imag
               
               internal_noise_real_jy = internal_noise_real / jy_to_K
               internal_noise_imag_jy = internal_noise_imag / jy_to_K
               print(internal_noise_real_jy[0:20])
               #always WODEN not wsclean
               data[:,0,0,0,pol_index,0] += internal_noise_real_jy
               data[:,0,0,0,pol_index,1] += internal_noise_imag_jy
                
               ####Y pol
               pol_index = 1
               internal_noise_real = uv_correlation_array_y[:,2+freq_index].real
               internal_noise_imag = uv_correlation_array_y[:,2+freq_index].imag
               
               internal_noise_real_jy = internal_noise_real / jy_to_K
               internal_noise_imag_jy = internal_noise_imag / jy_to_K
               
               data[:,0,0,0,pol_index,0] += internal_noise_real_jy
               data[:,0,0,0,pol_index,1] += internal_noise_imag_jy
          
               #now write out new uvfits file:
               hdulist.writeto(output_uvfits_name,overwrite=True)
               print("saved %s" % (output_uvfits_name))
   else:
      print("can only do with daniel freqs (1.28 MHz spacing)")   
      
if __name__ == "__main__":
    import argparse
    
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)
            
            from argparse import RawTextHelpFormatter
            
    parser = argparse.ArgumentParser(description="Run this to make WODEN sky models on \
            supercomputer with array job",formatter_class=SmartFormatter)

    parser.add_argument('--band', default='0',
        help='Just does one band at a time. Frequency depends on other arguments such as --daniel e.g. --band_nums=1')
            
    parser.add_argument('--daniel', default=False, action='store_true',
        help='Frequencies correspond to Daniels sims. Band 0 is at 50 MHz and frequency goes up in increments of 1.28 MHz')
               
    args = parser.parse_args()
    
    if args.band:
       band = int(args.band)

    
    
    add_noise_coupling_to_sim_uvfits(band=band,daniel=args.daniel) 
    
    
    