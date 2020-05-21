#!/usr/bin/env python
#code to analyse CenA data including flux-scale setting and measurements for paper

import os,sys

def get_scaling_factor_from_core(image_name,freq_MHz,alpha):
   #get this from image using masking etc eventually, for now just use kvis
   measured_new_flux_density_core = 1053.0
   model_freq_MHz = 843
   #from NED
   model_flux_density_core = 392  #(from Israel 1998)
   predicted_new_flux_density_core = model_flux_density_core * (freq_MHz / model_freq_MHz)**alpha
   
   scaling_factor = predicted_new_flux_density_core / measured_new_flux_density_core
   
   print("predicted core flux density is %0.3f Jy" % predicted_new_flux_density_core)
   print("measured core flux density is %0.3f Jy" % measured_new_flux_density_core)
   print("Scaling factor is %0.3f" % scaling_factor)
   
   
   
def source_find_image(image_name):
   output_base = image_name.split('/')[-1].split('.fits')[0]
   output_table = output_base+'-aegean_sources.fits'
   
   cmd = "BANE %s" % image_name
   print(cmd)
   os.system(cmd)

   cmd = "aegean --telescope MWA --autoload --table %s %s" % (output_table,image_name)
   print(cmd)
   os.system(cmd)   

def regrid_optical(template_imagename,input_imagename):
   imported_template_imagename = template_imagename.split('.fits')[0]+'.image'
   imported_input_imagename = input_imagename.split('.fits')[0]+'.image'
   regridded_imagename = input_imagename.split('.fits')[0]+'_regridded.image'
   regridded_fitsname = input_imagename.split('.fits')[0]+'_regridded.fits'
   
   casa_string = "importfits(fitsimage='%s',imagename='%s')" % (template_imagename,imported_template_imagename)
   casa_filename = 'casa_import_fits.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa -c --nocrashreport %s" % casa_filename
   print(cmd)
   os.system(cmd)

   casa_string = "importfits(fitsimage='%s',imagename='%s')" % (input_imagename,imported_input_imagename)
   casa_filename = 'casa_import_fits2.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa -c --nocrashreport %s" % casa_filename
   print(cmd)
   os.system(cmd)
    
   casa_string = "imregrid(imagename='%s',output='%s',template='%s')" % (imported_input_imagename,regridded_imagename,imported_template_imagename)
   casa_filename = 'casa_imregrid.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa -c --nocrashreport %s" % casa_filename
   print(cmd)
   os.system(cmd)

   casa_string = "exportfits(imagename='%s',fitsimage='%s')" % (regridded_imagename,regridded_fitsname)
   casa_filename = 'casa_export_fits.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa -c --nocrashreport %s" % casa_filename
   print(cmd)
   os.system(cmd)   


#image_name = "CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits"  
#get_scaling_factor_from_core(image_name,185.,-0.7)

#source_find_image(image_name)

template_imagename = 'rband_1sec_tr_fl_geo_ha.fits'
regrid_optical(template_imagename,input_imagename)
