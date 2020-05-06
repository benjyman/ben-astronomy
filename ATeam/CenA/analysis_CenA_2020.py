#!/usr/bin/env python
#code to analyse CenA data including flux-scale setting and measurements for paper

import os,sys

def get_scaling_factor_from_core(image_name,freq_MHz):
   print('')
   
def source_find_image(image_name):
   output_base = image_name.split('/')[-1].split('.fits')[0]
   output_table = output_base+'-aegean_sources.fits'
   
   cmd = "BANE %s" % image_name
   print(cmd)
   os.system(cmd)

   cmd = "aegean --telescope MWA --autoload --table %s %s" % (output_table,image_name)
   print(cmd)
   os.system(cmd)   


image_name = "CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits"  
source_find_image(image_name)
   