#!/usr/bin/env python
#code to analyse CenA data including flux-scale setting and measurements for paper

import os,sys
from astropy.io import fits
import numpy as np

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
   #optical stuff from Connor etc is difficult as the WSCS from PixelInsight doesnt seem to match what 
   #kvis or ds9 expects. To get around it I did this:
   #Detele  the CD_1, CD1_2, CD2_1, CD2_2 rotation parameters
   #leave the CROTA and B
   #Change the CDELT parameters to both be positive
   #Then can view with ds9 with coordinate grid in publication mode and can rotate the image to match what you want 
   #can also do overlays. explore ds9 command line options again
   #The below works for just the NML
   #...or make a very high res radio image just for the WCS and regrid with that!
   #wsclean -name CenA_optical_template -size 4500 4500 -scale 0.0003135 -pol xx -interval 1 3 -channel-range 1 3 1121321488.ms 
   imported_template_imagename = template_imagename.split('.fits')[0]+'.image'
   imported_input_imagename = input_imagename.split('.fits')[0]+'.image'
   regridded_imagename = input_imagename.split('.fits')[0]+'_regridded.image'
   regridded_fitsname = input_imagename.split('.fits')[0]+'_regridded.fits'
   
   casa_string = "importfits(fitsimage='%s',imagename='%s',overwrite=True)" % (template_imagename,imported_template_imagename)
   casa_filename = 'casa_import_fits.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)

   casa_string = "importfits(fitsimage='%s',imagename='%s',overwrite=True)" % (input_imagename,imported_input_imagename)
   casa_filename = 'casa_import_fits2.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)
    
   casa_string = "imregrid(imagename='%s',output='%s',template='%s',overwrite=True)" % (imported_input_imagename,regridded_imagename,imported_template_imagename)
   casa_filename = 'casa_imregrid.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)

   casa_string = "exportfits(imagename='%s',fitsimage='%s',overwrite=True)" % (regridded_imagename,regridded_fitsname)
   casa_filename = 'casa_export_fits.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)   

def edit_optical_header(optical_image,edhead_image_output_name):
   #if the regridding doesnt seem to get the RA or DEC going the right way it is 
   #probably because the cdelt is the wrong sign. multiply by -1. below
   n_decimals = 10
   with fits.open(optical_image) as hdulist:
      data = hdulist[0].data
      header1 = hdulist[0].header
   while('HISTORY' in header1):
      del header1['HISTORY']
   del header1['CD1_1']
   del header1['CD1_2']
   del header1['CD2_1']
   del header1['CD2_2']   
   cdelt_old = header1['CDELT1']  
   cdelt_new = np.around(np.abs(float(cdelt_old)),decimals=n_decimals)
   #only for Connors, not Rolfs: RA 
   header1['CDELT1'] = cdelt_new * -1.
   cdelt_new = np.around(np.abs(float(header1['CDELT2'])),decimals=n_decimals)
   header1['CDELT2'] = cdelt_new *-1. #only for Connors, not Rolfs:
   #actually these all need to be in the right numeric format
   old_val = np.around(float(header1['CRPIX1']),decimals=5)
   header1['CRPIX1'] = old_val
   old_val = np.around(float(header1['CRPIX2']),decimals=5)
   header1['CRPIX2'] = old_val
   old_val = np.around(float(header1['CRVAL1']),decimals=n_decimals)
   header1['CRVAL1'] = old_val   
   old_val = np.around(float(header1['CRVAL2']),decimals=n_decimals)
   header1['CRVAL2'] = old_val            
   old_val = np.around(float(header1['CROTA1']),decimals=n_decimals)
   header1['CROTA1'] = old_val    
   old_val = np.around(float(header1['CROTA2']),decimals=n_decimals)
   header1['CROTA2'] = old_val     
   
   #print header1
   #sys.exit()
   
   #write new fits file
   fits.writeto(edhead_image_output_name,data,clobber=True)
   fits.update(edhead_image_output_name,data,header=header1)
   print("wrote image %s" %  edhead_image_output_name)     

def regrid_to_J2000(input_imagename):
   imported_input_imagename = input_imagename.split('.fits')[0]+'.image'
   regridded_imagename = input_imagename.split('.fits')[0]+'_regridJ2000.image'
   regridded_fitsname = input_imagename.split('.fits')[0]+'_regridJ2000.fits'

   casa_string = "importfits(fitsimage='%s',imagename='%s',overwrite=True)" % (input_imagename,imported_input_imagename)
   casa_filename = 'casa_import_fits2.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)
    
   casa_string = "imregrid(imagename='%s',output='%s',template='J2000',overwrite=True)" % (imported_input_imagename,regridded_imagename)
   casa_filename = 'casa_imregrid.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)

   casa_string = "exportfits(imagename='%s',fitsimage='%s',overwrite=True)" % (regridded_imagename,regridded_fitsname)
   casa_filename = 'casa_export_fits.sh'
   with open(casa_filename,'w') as f:
      f.write(casa_string)
   cmd = "casa --nocrashreport -c %s" % casa_filename
   print(cmd)
   os.system(cmd)   

def spectral_index_map(image_1_name,image_2_name,freq_MHz_low,freq_MHz_high,output_name_base,mask_level):
   spec_index_im_name = "%s_spec_index.im" % output_name_base
   spec_index_fits_name = "%s_spec_index.fits" % output_name_base
   #use image_1 as template, image_1 to be lowest freq
   
   #low freq
   hdulist = fits.open("%s" % (image_1_name))
   image_header = hdulist[0].header
   image_data = hdulist[0].data 
   target_bmaj = float(image_header['BMAJ']) *60.*60. 
   target_bmin = float(image_header['BMIN']) *60.*60.
   target_bpa = float(image_header['BPA'])    
   #freq_Hz_low = float(image_header['CRVAL3'])
   #freq_MHz_low = freq_Hz_low/1000000.
   hdulist.close()

   #hdulist = fits.open("%s" % (image_2_name))
   #image_header = hdulist[0].header
   #image_data = hdulist[0].data 
   #freq_Hz_high = float(image_header['CRVAL3'])
   #freq_MHz_high = freq_Hz_high/1000000.
   #hdulist.close()
   
   #read image_1 into miriad
   image_name_base_low = 'low'
   image_name_base_high = 'high'
   
   #read in images to miriad
   im_name_low = "%s.im" % image_name_base_low
   im_name_high = "%s.im" % image_name_base_high
   im_name_high_regrid = "%s_regrid.im" % image_name_base_high
   im_name_high_regrid_convol = "%s_regrid_convol.im" % image_name_base_high
   
   cmd = "rm -rf %s %s %s %s" % (im_name_low,im_name_high,im_name_high_regrid,im_name_high_regrid_convol)
   print(cmd)
   os.system(cmd)
   
   cmd = "fits in=%s out=%s op=xyin" % (image_1_name,im_name_low)
   print(cmd)
   os.system(cmd)

   cmd = "fits in=%s out=%s op=xyin" % (image_2_name,im_name_high)
   print(cmd)
   os.system(cmd)

   #regrid im2 to im1
   cmd = "regrid in=%s out=%s tin=%s" % (im_name_high,im_name_high_regrid,im_name_low)
   print(cmd)
   os.system(cmd) 
   
   #smooth im2 down
   cmd = "convol map=%s fwhm=%4f,%4f pa=%4f options=final out=%s " % (im_name_high_regrid,target_bmaj,target_bmin,target_bpa,im_name_high_regrid_convol)
   print(cmd)
   os.system(cmd)
   
   ##############
   #OPTIONAL FITS output
   #output the regrid convol
   regrid_convol_fitsname = "rosat_diffuse_north_regrid_convol_mwa_2deg.fits"
   cmd = "rm -rf %s" % (regrid_convol_fitsname)
   print(cmd)
   os.system(cmd)  

   cmd = "fits in=%s out=%s op=xyout" % (im_name_high_regrid_convol,regrid_convol_fitsname)
   print(cmd)
   os.system(cmd) 
 
   #make the spec index map
   cmd = "rm -rf %s %s" % (spec_index_im_name,spec_index_fits_name)
   print(cmd)
   os.system(cmd)
  
   cmd = "maths exp=log\(%s/%s\)/log\(%0.2f/%0.2f\) mask=%s>%0.5f  out=%s" % (im_name_high_regrid_convol,im_name_low,freq_MHz_high,freq_MHz_low,im_name_high_regrid_convol,mask_level,spec_index_im_name)
   print(cmd)
   os.system(cmd) 
   
   cmd = "fits in=%s out=%s op=xyout" % (spec_index_im_name,spec_index_fits_name)
   print(cmd)
   os.system(cmd)

def regrid_concvol(image_1_name,image_2_name_list,target_bmaj_deg,target_bmin_deg,target_bpa_deg,output_name_base):
   output_im_name = "%s_mosaic.im" % output_name_base
   output_fits_name = "%s_mosaic.fits" % output_name_base
   #use image_1 as template
   
   target_bmaj = float(target_bmaj_deg) *60.*60. 
   target_bmin = float(target_bmin_deg) *60.*60.
   target_bpa = float(target_bpa_deg)    
   
   #read image_1 into miriad
   image_name_base_1 = 'template'
   im_name_1 = "%s.im" % image_name_base_1

   cmd = "rm -rf %s %s %s " % (im_name_1,output_im_name,output_fits_name)
   print(cmd)
   os.system(cmd)
   
   cmd = "fits in=%s out=%s op=xyin" % (image_1_name,im_name_1)
   print(cmd)
   os.system(cmd)
   
   linmos_image_list = []
   for image_2_name in image_2_name_list:
      #fix header bits of images to be regridded
      hdulist = fits.open("%s" % (image_2_name))
      image_header = hdulist[0].header
      image_data = hdulist[0].data     
      if image_header['CTYPE1'] == 'EQU--CAR':
         image_header['CTYPE1'] = 'RA--CAR'
      if image_header['CTYPE2'] == 'EQU--CAR':
         image_header['CTYPE2'] = 'DEC--CAR'
      try:
         freq = image_header['CTYPE3']
      except:
         image_header['CTYPE3'] = 'FREQ'
         image_header['CRVAL3'] = 184955000.
         image_header['CRPIX3'] = 1.
         image_header['CDELT3'] = 30720000.
         image_header['CUNIT3'] = 'Hz'
      #resolution of PSPC on XRT on Rosat about 20''
      try:
         bmaj = image_header['BMAJ']
      except:
         image_header['BMAJ'] = 20.
         image_header['BMIN'] = 20.    
         image_header['BPA'] = 0.  
         image_header['TELESCOP'] = 'ATCA'
      #print(image_header)
      #sys.exit()
      #freq_Hz_low = float(image_header['CRVAL3'])
      #freq_MHz_low = freq_Hz_low/1000000.
      hdulist.close()
   
      #write new fits file
      fits.writeto(image_2_name,image_data,clobber=True)
      fits.update(image_2_name,image_data,header=image_header)
      print("wrote image %s" %  image_2_name)  
   
   
      image_name_base_2 = image_2_name.split('.fits')[0]
      output_im_2_name = "%s_regrid_convol.im" % image_name_base_2
      output_im_2_fits_name = "%s_regrid_convol.fits" % image_name_base_2
      
      #read in images to miriad
      im_name_2 = "%s.im" % image_name_base_2
      im_name_2_regrid = "%s_regrid.im" % image_name_base_2
      
      cmd = "rm -rf %s %s %s %s" % (im_name_2,output_im_2_name,output_im_2_fits_name,im_name_2_regrid)
      print(cmd)
      os.system(cmd)
   
      cmd = "fits in=%s out=%s op=xyin" % (image_2_name,im_name_2)
      print(cmd)
      os.system(cmd)
   
      #regrid im2 to im1
      cmd = "regrid in=%s out=%s tin=%s" % (im_name_2,im_name_2_regrid,im_name_1)
      print(cmd)
      os.system(cmd) 
      
      #smooth im2 down
      #cmd = "convol map=%s fwhm=%4f,%4f pa=%4f options=final out=%s " % (im_name_2_regrid,target_bmaj,target_bmin,target_bpa,output_im_2_name)
      #rosat (no final):
      cmd = "convol map=%s fwhm=%4f,%4f pa=%4f out=%s " % (im_name_2_regrid,target_bmaj,target_bmin,target_bpa,output_im_2_name)
      print(cmd)
      os.system(cmd)
      
      linmos_image_list.append(output_im_2_name)
      
      ##############
      #OPTIONAL FITS output
      #output the regrid convol
      cmd = "fits in=%s out=%s op=xyout" % (output_im_2_name,output_im_2_fits_name)
      print(cmd)
      os.system(cmd)
   
   #use miriad linmos to combine images on to the template image grid.
   linmos_image_list_string = ','.join(linmos_image_list)

   cmd = "linmos in=%s out=%s" % (linmos_image_list_string,output_im_name)
   print(cmd)
   os.system(cmd)

#spectral index ASKAP MWA:
#spectral_index_map('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits','CenA_i.fits',185,1400,'CenA_185_1400_MHz',0.3)


#regridding x ray rosat from galacto:
#regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits','g000p00r1b120pm.fits',0.5,0.5,0,'CenA_rosat_north_to_mwa')
cen_A_rosat_p10_list = ['932527p-p10.fits','932429p-p10.fits']
regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits',cen_A_rosat_p10_list,0.24,0.24,0,'CenA_middle_rosat_postage_to_mwa_lowband')
sys.exit()

#image_name = "CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits"  
#get_scaling_factor_from_core(image_name,185.,-0.7)

#source_find_image(image_name)

#regrid VLA image to J2000:
#image_name = "NGC_5128_I_21cm_chs1996.fits"
#regrid_to_J2000(image_name)
#sys.exit()

#template_imagename = 'rband_1sec_tr_fl_geo_ha.fits'
#template_imagename = 'CenA_optical_template-image.fits'
template_imagename = 'Ha_cont_subtracted_via_rband_scaling.fits'

#input_imagename = 'CenA_WCS_edhead.fits'
#input_imagename = 'CenA_WCS_Ha_edhead.fits'
#regrid_optical(template_imagename,input_imagename)

#this works well.
#regrid all of the 'new connor' images so that we can clearly show the new filament, that it is not a HII region,
# and how it lines up with the radio (and x-ray!) and is a similar distance from the core to the inner filament in the north
#in benjamin@namorrodor:/md0/ATeam/CenA/paper_2020/optical/new_connor/new_connor
#first edit the headers    #Detele  the CD1_1, CD1_2, CD2_1, CD2_2 rotation parameters, leave the CROTA and B, Change the CDELT parameters to both be positive
#cp /md0/ATeam/CenA/paper_2020/optical/CenA_optical_template-image.fits .
#new connor:
#input_name_list = ['1_Stacked_Image.fits','2_Gradient_Removal.fits','3_Separate_HII_regions_from_Ha.fits','4_Noise_Reduction.fits','5_Combined_Ha_with_RGB.fits','6_Histogram_Stretch.fits','7_Artifact_fixing_final_image.fits']
#mike sidonio:
#input_name_list = ['CenA_Ha_1050min_median.fits','CenA_lum_1470min_median_grad.fits']
#input_name_list = ['CenA_Ha_1050min_median_solved.fits','CenA_Lum_1470min_median_grad_solved.fits']
input_name_list = ['3_Separate_HII_regions_from_Ha.fits']
#i know this one works:
#input_name_list = ['CenA_WCS.fits']
#input_name_list = ['CenA_WCS_Ha.fits']
for input_name in input_name_list:
   edhead_name = input_name.split('.fits')[0]+'_edhead.fits'
   edit_optical_header(input_name,edhead_name)
   regrid_optical(template_imagename,edhead_name)
   #try no edhead for mikes:
   #regrid_optical(template_imagename,input_name)
   
#benjamin@namorrodor:/md0/ATeam/CenA/paper_2020/optical/new_connor
#kvis ../../../CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08.fits *_edhead_regridded.fits ../../x_ray/XMM_2001.fits 
#or on Mac: ~/Documents/Astronomy/ATeam/CenA/2020_paper/draft/overlays/optical/new_connor/regridded/
#kvis *_edhead_regridded.fits  ~/Documents/Astronomy/ATeam/CenA/CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits  ../../../x_ray/XMM_2001.fits  
#radio contours 0.1 0.2 0.4 0.8 1 2 4 8 16 32 64
#x xay contours 5 6 7 8 9
#look at new filament at 13:24:34.9, -43:09:11.09
#show doesnt get subtracted by HII regions, matches exactly with radio 'bulge'
#dealing with HII regions: http://www.arciereceleste.it/tutorial-pixinsight/cat-tutorial-eng/85-enhance-galaxy-ha-eng











