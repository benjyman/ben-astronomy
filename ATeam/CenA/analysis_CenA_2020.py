#!/usr/bin/env python
#code to analyse CenA data including flux-scale setting and measurements for paper

import os,sys
from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
import math

G = 6.674e-11
M_solar = 1.98847e30 
c = 299792458.
M_bh_cena = 1.0e8 * M_solar
pc_m = 3.0857e16 
yr_sec = 3.1557600e7
mp = 1.6726219e-27    #kg
kb = 1.38064852e-23
eV_J = 1.602176634e-19
erg_J = 1.0e-7
dyn_cm2_Pa = 0.1


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

def regrid_optical(template_imagename,input_imagename,smooth=0):
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
   regridded_fitsname_smooth = input_imagename.split('.fits')[0]+'_regridded_smooth.fits'
   
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
   
   if smooth!=0:
      hdulist = fits.open("%s" % (regridded_fitsname))
      image_header = hdulist[0].header
      image_data = np.nan_to_num(hdulist[0].data)
      hdulist.close()

      data_smooth = ndimage.gaussian_filter(image_data, sigma=(smooth, smooth), order=0)
   
      #write to fits:
      fits.writeto(regridded_fitsname_smooth,data_smooth,clobber=True)
      fits.update(regridded_fitsname_smooth,data_smooth,header=image_header)
      print("wrote image %s" %  regridded_fitsname_smooth) 
      
    

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
   #regrid_convol_fitsname = "rosat_diffuse_north_regrid_convol_mwa_2deg.fits"
   #cmd = "rm -rf %s" % (regrid_convol_fitsname)
   #print(cmd)
   #os.system(cmd)  
   #
   #cmd = "fits in=%s out=%s op=xyout" % (im_name_high_regrid_convol,regrid_convol_fitsname)
   #print(cmd)
   #os.system(cmd) 
 
   #make the spec index map
   cmd = "rm -rf %s %s" % (spec_index_im_name,spec_index_fits_name)
   print(cmd)
   os.system(cmd)
  
   cmd = "maths exp=log\(%s/%s\)/log\(%0.2f/%0.2f\) mask=%s.gt.%0.5f  out=%s" % (im_name_high_regrid_convol,im_name_low,freq_MHz_high,freq_MHz_low,im_name_high_regrid_convol,mask_level,spec_index_im_name)
   print(cmd)
   os.system(cmd) 
   
   cmd = "fits in=%s out=%s op=xyout" % (spec_index_im_name,spec_index_fits_name)
   print(cmd)
   os.system(cmd)

def regrid_concvol(image_1_name,image_2_name_list,target_bmaj_deg,target_bmin_deg,target_bpa_deg,output_name_base):
   mosaic_im_name = "%s_mosaic.im" % output_name_base
   mosaic_fits_name = "%s_mosaic.fits" % output_name_base
   mosaic_im_name_smooth = "%s_mosaic_smooth.im" % output_name_base
   mosaic_fits_name_smooth = "%s_mosaic_smooth.fits" % output_name_base
   #use image_1 as template
   
   #read template image to get data grid array:
   hdulist = fits.open("%s" % (image_1_name))
   image_header_1 = hdulist[0].header
   image_data_1 = hdulist[0].data  
   pix_size_deg = abs(float(image_header_1['CDELT1']))
   hdulist.close()
   
   target_bmaj = float(target_bmaj_deg) *60.*60. 
   target_bmin = float(target_bmin_deg) *60.*60.
   target_bpa = float(target_bpa_deg)    
   
   #read image_1 into miriad
   image_name_base_1 = 'template'
   im_name_1 = "%s.im" % image_name_base_1

   cmd = "rm -rf %s %s %s %s %s " % (im_name_1,mosaic_im_name,mosaic_fits_name,mosaic_im_name_smooth,mosaic_fits_name_smooth)
   print(cmd)
   os.system(cmd)
   
   cmd = "fits in=%s out=%s op=xyin" % (image_1_name,im_name_1)
   print(cmd)
   os.system(cmd)
   
   sum_image_data = image_data_1 * 0
   count_image_data_sum = image_data_1 * 0
   for image_2_name in image_2_name_list:
      #fix header bits of images to be regridded
      hdulist = fits.open("%s" % (image_2_name))
      image_header = hdulist[0].header
      image_data = hdulist[0].data     
      if image_header['CTYPE1'] == 'EQU--CAR':
         image_header['CTYPE1'] = 'RA--CAR'
      if image_header['CTYPE2'] == 'EQU--CAR':
         image_header['CTYPE2'] = 'DEC--CAR'
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
      fits_name_2_regrid = "%s_regrid.fits" % image_name_base_2
      
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
      
      cmd = "fits in=%s out=%s op=xyout" % (im_name_2_regrid,fits_name_2_regrid)
      print(cmd)
      os.system(cmd)
      
      ##smooth im2 down
      cmd = "convol map=%s fwhm=%4f,%4f pa=%4f options=final out=%s " % (im_name_2_regrid,target_bmaj,target_bmin,target_bpa,output_im_2_name)
      print(cmd)
      os.system(cmd)
      
      ##############
      #FITS output
      #output the regrid convol
      cmd = "fits in=%s out=%s op=xyout" % (output_im_2_name,output_im_2_fits_name)
      print(cmd)
      os.system(cmd)
   
      #sum up the un convol images
      #read in fits file and get data array
      hdulist = fits.open("%s" % (fits_name_2_regrid))
      image_header_convol = hdulist[0].header
      image_data_convol = np.nan_to_num(hdulist[0].data)
      hdulist.close()
      
      count_image_data = image_data_convol * 0
      count_image_data[image_data_convol>0] = 1
      
      sum_image_data += image_data_convol
      count_image_data_sum += count_image_data
      
      
   
   av_image_data = sum_image_data / count_image_data_sum
   av_image_data[np.isnan(av_image_data)] = 0.0
   av_image_data[np.isinf(av_image_data)] = 0.0
   
   #smooth
   av_image_data_smooth = ndimage.gaussian_filter(av_image_data, sigma=(3, 3), order=0)
   
   #write to fits:
   fits.writeto(mosaic_fits_name,av_image_data,clobber=True)
   fits.update(mosaic_fits_name,av_image_data,header=image_header_1)
   print("wrote image %s" %  mosaic_fits_name) 
   #cmd = "linmos in=%s out=%s" % (linmos_image_list_string,output_im_name)
   #print(cmd)
   #os.system(cmd)

   fits.writeto(mosaic_fits_name_smooth,av_image_data_smooth,clobber=True)
   fits.update(mosaic_fits_name_smooth,av_image_data_smooth,header=image_header_1)
   print("wrote image %s" %  mosaic_fits_name_smooth) 
   #cmd = "linmos in=%s out=%s" % (linmos_image_list_string,output_im_name)
   #print(cmd)
   #os.system(cmd)

def calculate_outflow_properties():
#G17 scale defs
   r_s_m = 2. * G * M_bh_cena / c**2
   r_s_pc = r_s_m / pc_m
   print("r_s_pc is %0.9f = %0.3f kpc" % (r_s_pc,r_s_pc/1000.0))
   print("5 * r_s_pc is %0.9f pc = %0.3f kpc" % (r_s_pc*5.,r_s_pc*5./1000.) )
   print("Micro: 1000 * r_s_pc is %0.9f pc = %0.3f kpc" % (r_s_pc*1000.,r_s_pc*1000./1000.) )
   print("Meso: 1e6 * r_s_pc is %0.9f pc = %0.3f kpc" % (r_s_pc*1.e6,r_s_pc*1e6/1000.) )
   print("Macro: 1e9 * r_s_pc is %0.9f pc = %0.3f kpc" % (r_s_pc*1.e9,r_s_pc*1e9/1000.) )
   print("r_vir 1e10 * r_s_pc = r_vir is %0.9f pc = %0.3f kpc" % (r_s_pc*1.e10,r_s_pc*1.e10/1000.) )
   
   #israel et al 2017 inner 500 pc outflow
   vel_km_s = 60.
   vel_m_s = vel_km_s * 1000.
   mass_out_kg = 1.2e7 * M_solar
   projected_length_pc = 174.
   line_of_sight_theta_deg = 24.
   line_of_sight_theta_rad = line_of_sight_theta_deg/180. * math.pi
   actual_length_pc = projected_length_pc / math.sin(line_of_sight_theta_rad)
   actual_velocity_m_s = vel_m_s / math.cos(line_of_sight_theta_rad)
   outflow_time_s = (actual_length_pc * pc_m) / actual_velocity_m_s
   outflow_time_Myr = outflow_time_s / (1e6 * yr_sec) #(1e6 * 365. * 24 * 60. * 60.)
   mass_outflow_rate_kg_s = mass_out_kg / outflow_time_s
   mass_outflow_rate_solarmass_yr = mass_outflow_rate_kg_s / M_solar * yr_sec
   
   print("Israel et al 2017 inner 500pc ouflow time is %0.5f Myr" % outflow_time_Myr)
   print("mass outflow rate %0.5f solar mass/yr " % mass_outflow_rate_solarmass_yr)
   print("actual velocity %0.5f km/s " % (actual_velocity_m_s/1000.))
   
   #krol et al 2020 kpc cold outflow
   line_of_sight_theta_deg = 45.
   vel_km_s = 800.
   #WEST:
   
   vel_m_s = vel_km_s * 1000.
   volume_kpc3 = 0.43
   density_n_cm3 = 0.5
   mass_out_kg = 5.5e6 * M_solar
   pressure_dyn_cm2 = 1.6e-10
   internal_energy_erg = 3e54
   projected_length_pc = 2.7 * 1e3
   
   
   line_of_sight_theta_rad = line_of_sight_theta_deg/180. * math.pi
   actual_length_pc = projected_length_pc / math.sin(line_of_sight_theta_rad)
   print(actual_length_pc)
   actual_velocity_m_s = vel_m_s / math.cos(line_of_sight_theta_rad)
   print(actual_velocity_m_s)
   outflow_time_s = (actual_length_pc * pc_m) / actual_velocity_m_s
   print(outflow_time_s)
   outflow_time_Myr = outflow_time_s / (1e6 * yr_sec) #(1e6 * 365. * 24 * 60. * 60.)
   mass_outflow_rate_kg_s = mass_out_kg / outflow_time_s
   mass_outflow_rate_solarmass_yr = mass_outflow_rate_kg_s / M_solar * yr_sec

   print("Krol et al 2017 2.7 kpc West ouflow time is %0.5f Myr" % outflow_time_Myr)
   print("Krol et al 2017 2.7 kpc West mass outflow rate %0.5f solar mass/yr " % mass_outflow_rate_solarmass_yr)

   #EAST:
   vel_m_s = vel_km_s * 1000.
   volume_kpc3 = 0.23
   density_n_cm3 = 0.3
   mass_out_kg = 1.6e6 * M_solar
   pressure_dyn_cm2 = 1.2e-10
   internal_energy_erg = 1e54
   projected_length_pc = 2.7 * 1e3
   line_of_sight_theta_rad = line_of_sight_theta_deg/180. * math.pi
   actual_length_pc = projected_length_pc / math.sin(line_of_sight_theta_rad)
   actual_velocity_m_s = vel_m_s / math.cos(line_of_sight_theta_rad)
   outflow_time_s = (actual_length_pc * pc_m) / actual_velocity_m_s
   outflow_time_Myr = outflow_time_s / (1e6 * yr_sec) #(1e6 * 365. * 24 * 60. * 60.)
   mass_outflow_rate_kg_s = mass_out_kg / outflow_time_s
   mass_outflow_rate_solarmass_yr = mass_outflow_rate_kg_s / M_solar * yr_sec

   print("Krol et al 2017 2.7 kpc East ouflow time is %0.5f Myr" % outflow_time_Myr)
   print("Krol et al 2017 2.7 kpc East mass outflow rate %0.5f solar mass/yr " % mass_outflow_rate_solarmass_yr)

   #from G17 unified (3.1):
   #adiabatic index
   gamma = 5./3.
   mu = 0.62 #is the average atomic weight for a fully ionized plasma with 25% He in mass
   Tx_ev = 0.5 * 1000.
   Tx_J = (eV_J * Tx_ev) / kb
   #adiabatic sound speed
   cs_ms = (gamma * kb * Tx_J / (mu * mp))**(0.5)
   
   #constant = ((gamma * kb) / (mu * mp))**(0.5)
   #print(constant)
   #cs_ms = 1.5e4 * (Tx_ev)**(0.5)
   #the factor 100 is there I think because G17 is wrong in its 10**4 value in sec 3.1 eqn
   cs_kms = cs_ms / 1000. 
   print("cs is %E km per sec" % cs_kms)
   
   #G17 stuff:
   #from MWA image: 
   r_thermalisation_kpc = 75.
   #r_thermalisation_kpc = 55. * (T_x_7p4)**(3./2.)  (EQN 26)
   #T_x_7p4 = (r_thermalisation_kpc/55.0)**(2./3.)
   #print("T_x_7p4 is %E keV" % T_x_7p4)
   
   #EQN 18:
   #v_out_km_s = 1.4e4 * T_x_7p4**(0.5)
   #print("inner vel v_out_km_s is %E " % v_out_km_s)
   
   #(EQN 23 )
   #entrainment_factor_hot = 40. * 1.**(2./3.) / T_x_7p4
   #print("entrainment_factor_hot is %E " % entrainment_factor_hot)
   
   #outer vel (EQN 21)
   #V_OUT_km_s = entrainment_factor_hot**(-0.5) * v_out_km_s
   #print("macro vel V_OUT_km_s is %E " % V_OUT_km_s)
   
   
   #EQN23 
   
   #sigma_steradians = 4.*np.pi
   #for plasma density:  Look at (from Max) OSullivan et al 2003 or Bogdan et al. 2010 for Chandra, plus Borkar 2020 for inner 180 pc
   #from borkar 2020 / kraft 2003 - use n(r) = n0 (r/r0)**-1.5, with n0 = 3.7e-2 cm-3, and r0=10 pc
   #hot_plasma_density = 0.4 * mp * 1.0e3               # 1.e-25 #
   #print("hot_plasma_density is %E" % hot_plasma_density)
   #ref_radius_kpc = 2.7
   #constant_eqn_23 = (sigma_steradians/(4.*np.pi)) * (hot_plasma_density/1.0e-25)* 40.
   #print("constant_23 = %0.5f" % constant_eqn_23)
   
   #neff NML
   #radius = 30.
   #width = 24.
   #sigma_steradians = 2. * 2. * np.pi * (1. - np.cos((math.atan(width/radius) / 2)))  # see wikepedia page on solid angle! ...extra factor of 2 because bipolar
   #print(sigma_steradians)
   #frac_sphere = sigma_steradians / (4.*np.pi)
   #print(frac_sphere)
   #constant_eqn_23 = (sigma_steradians/(4.*np.pi)) * (hot_plasma_density/1.0e-25)* 40.
   #print("constant_23 neff = %0.5f" % constant_eqn_23)
   
   
   #Gaspari eqn 26 (assumes omega is 4pi)
   #try using T_x_kraft_2009
   #sigma = 4.*np.pi
   radius = 63. #adj
   width = 90.  #opp = width/2
   print("########from MWA image estimate that the outflow has spread to a width of %s kpc at a distance of %s kpc from the core" % (radius,width))
   #tan(theta) = opp/adj
   theta = math.atan((width/2.)/radius)
   print("theta is %0.3f rad" % theta) 
   theta_deg = theta / np.pi * 180.
   print("#####theta_deg is %0.3f deg" % theta_deg) 
   opening_angle_deg = 2.*theta_deg
   sigma = 2. * np.pi * (1. - np.cos(theta))  # see wikepedia page on solid angle!
   print("#####opening solid angle sigma is %0.3f sr" % sigma)
   #T_x_kraft_2009 = 0.35 #keV
   #T_x_7p4 = T_x_kraft_2009 / 2.2
   const_eqn_26 = 55. / (sigma*2 / (4.*np.pi)) #kpc
   print("######const_eqn_26 is %E" % const_eqn_26)
   #r_th = const_eqn_26 * T_x_7p4**(3./2.)
   #print("r_th is %E kpc" % r_th)
   
   
   #Reversing this so we have r_th measured from radio:
   r_th = 75.
   const_eqn_26 = 55. / (sigma / (4.*np.pi)) #kpc
   T_x_7p4 = (r_th / const_eqn_26)**(2./3.)
   print("######T_x_7p4 is %E keV" % T_x_7p4)
   T_x_unscaled = T_x_7p4 * 2.2
   print("######T_x_unscaled is %E keV" % T_x_unscaled)

   #MWA south (arcmin):
   #radius = 63.
   #width = 90.
   #sigma_steradians = 2. * 2. * np.pi * (1. - np.cos((math.atan(width/radius) / 2)))  # see wikepedia page on solid angle! ...extra factor of 2 because bipolar
   ##sigma_steradians = 4.*np.pi
   #print(sigma_steradians)
   #sigma_degrees = sigma_steradians / np.pi * 180.
   #print("opening angle is south is %0.2f degrees" % sigma_degrees)
   #frac_sphere = sigma_steradians / (4.*np.pi)
   #print(frac_sphere)
   
   constant_eqn_23 = 40. * (sigma/(4.*np.pi))**(2./3.) #* (hot_plasma_density/1.0e-25) * (ref_radius_kpc/1.)
   print("########constant_23 MWA south = %0.5f" % constant_eqn_23)
   
   ##eqn 26
   #const_eqn_26 = 55. * (sigma_steradians / (4.*np.pi))**(2./3.) * (hot_plasma_density/1.0e-25)**(2./3.) * (ref_radius_kpc/1.)**(2./3.)
   #print("const_eqn_26 = %0.5f" % const_eqn_26)
   
   ##r_thermalisation_kpc = const_eqn_26 * (T_x_7p4)**(3./2.)  (EQN 26)
   #T_x_7p4 = (r_thermalisation_kpc/const_eqn_26)**(2./3.)
   #if we take the kraft et al value of Tx:
   ##T_x_7p4 = 0.35
   ##print("T_x_7p4 is %E keV" % T_x_7p4)
   
   #Tx_unscaled = T_x_7p4 * 2.2
   #print("Tx_unscaled is %E keV" % Tx_unscaled)
   
   
   #EQN 23
   entrainment_factor_hot = constant_eqn_23 * 1.**(2./3.) / T_x_7p4
   print("########entrainment_factor_hot is %E " % entrainment_factor_hot)
   
   #EQN 18: (inner outflow)
   v_out_km_s = 1.4e4 * T_x_7p4**(0.5)
   print("#########inner vel v_out_km_s is %E " % v_out_km_s)
   
   #outer outflow vel (EQN 21)
   V_OUT_km_s = entrainment_factor_hot**(-0.5) * v_out_km_s
   print("#########macro vel V_OUT_km_s is %E " % V_OUT_km_s) 
   
   #try equation 20 G17
   #ro_nought = 1. #kpc
   #M_OUT = sigma * ro_nought * r_nought * r_th * V_OUT_km_s
   
   
   ##EQN 20 outer mass outflow 
   #M_OUT = (sigma_steradians / (4.*np.pi)) *  (hot_plasma_density/1.0e-25) * (ref_radius_kpc/1.) * V_OUT_km_s
   #print("Don't know units! macro mass outflow rate M_OUT is %E " % M_OUT)
   ##EQN 19 outer mass outflow 
   
   #ENQN 7 Macro mass inflow rate M_cool (also look at sun2012)
   #T_x_7p4 = 0.8
   M_cool = 1.1 * T_x_7p4**2
   print("#########Macro mass inflow rate M_cool is %E solar mass / yr" % M_cool) 
   
   #micro outflow rate is approximately equal to the inflow rate (only a few percent of the mass is actually able to fall into the BH)
   #EQN 16
   m_out = M_cool
   print("#########micro mass outflow rate m_out is %E solar mass / year (roughly same as M_cool)" % m_out) 
   
   M_OUT = m_out * entrainment_factor_hot 
   print("#########Macro mass outflow rate (hot) m_out is %E solar mass / year" % M_OUT) 
   
   #The final link! G17 Eqn 12 - the micro mass inflow rate
   M_dot_in_micro = 0.03 * M_cool * T_x_7p4 
   print("#########micro mass inflow rate M_dot_in_micro is %E solar mass / year" % M_dot_in_micro) 
   
   #check if power of wind can create the x-ray clouds consistent with Kaft et al 2009:
   #G17 eqn 15:
   #M_OUT = 0.5
   #V_OUT_km_s = 2000.
   M_OUT_kg_s = M_OUT * M_solar / yr_sec
   print("M_OUT_kg_s is %E" % M_OUT_kg_s)
   power_out_macro_W = 0.5 * M_OUT_kg_s * (V_OUT_km_s * 1000.) ** 2.0
   print("power_out_macro_W is %E" % power_out_macro_W)
   power_out_macro_erg_s = power_out_macro_W / erg_J
   print("power_out_macro_erg_s is %E" % power_out_macro_erg_s)
   
   #NML power from kraft et al 2009:
   kraft_nml_radius_m = 35. * 1000. * pc_m
   kraft_nml_volume_m3 = 4./3. * np.pi * kraft_nml_radius_m**3
   kraft_nml_pressure_dyn_cm2 = 8.5e-13
   kraft_nml_pressure_N_m2 = kraft_nml_pressure_dyn_cm2 * 0.1
   E_mech_J = 4. * kraft_nml_pressure_N_m2 * kraft_nml_volume_m3
   E_mech_erg = E_mech_J / erg_J
   print("Kraft nml E_mech_erg is %E erg" % E_mech_erg)
   time_for_power_out_macro_s = E_mech_erg / power_out_macro_erg_s
   time_for_power_out_macro_Myr = time_for_power_out_macro_s / yr_sec / 1.e6
   print("time_for_power_out_macro_Myr is %E" % time_for_power_out_macro_Myr)
   
   
   
   ##Israel hot ionised  outflow
   #this number density is just for e-, dont know hot gas density at 1kpc - maybe could use krol at 2.7?
   #number_density_cm3 = 10.
   #density_g_cm3 = number_density_cm3 * mp #* mu
   #print("density_g_cm3 is %E" % density_g_cm3)
   #constant_eqn_23 = (sigma_steradians/(4.*np.pi)) * (density_g_cm3/1.0e-25)* 40.
   #print("constant_23 MWA south and israel CND density = %0.5f" % constant_eqn_23)
   
   
   
   #stuff from Tombesi 2014 (table 1)
   CenA_v_out_c = 0.005
   CenA_v_out_kms = CenA_v_out_c * c / 1000.
   s3C445_v_out_kms = 0.034 * c / 1000.
   #print("CenA UFO velocity is %0.5f km/sec for %0.5f c" % (CenA_v_out_kms,CenA_v_out_c))
   #print("compare to 3C445 UFO velocity is %0.5f km/sec" % s3C445_v_out_kms)
   
   #radio jet
   v_in_c = .1
   v_kms = v_in_c * c / 1000.
   #print("CenA UFO velocity is %0.5f km/sec for %0.5f c" % (v_kms,v_in_c))
   v_in_c = .45
   v_kms = v_in_c * c / 1000.
   #print("CenA UFO velocity is %0.5f km/sec for %0.5f c" % (v_kms,v_in_c))
   
   #temperature conversions
   temp_kT_k_eV = 0.16
   temp_kT_k_J = 1000. * eV_J * temp_kT_k_eV
   temp_T_K = temp_kT_k_J / kb
   #print("temp_kT %0.5f keV = temp_K %4E K " % (temp_kT_k_eV,temp_T_K))
   
   #HI stuff
   total_mass = 1.5e8 * M_solar
   
   #Eilek density: 
   B_field_micro = 1
   n_p = 3e-7 * B_field_micro**2
   #print("np frowm Eilek with 1 micro Gauss field is %E" % n_p)
   
   
   #Neff et al 2015b stuff:
   thermal_pressure_ism_neff_dyn_cm2 = 8.5e-13
   thermal_pressure_ism_neff_Pa = thermal_pressure_ism_neff_dyn_cm2 * 0.1
   power_1_neff_erg_s = 1.0e42
   power_1_neff_J_s = power_1_neff_erg_s * 1e-7
   #velocity_wind_1_neff_ms = 2000.0 * 1000.
   
   power_2_neff_erg_s = 1.0e43
   power_2_neff_J_s = power_2_neff_erg_s * erg_J
   #velocity_wind_2_neff_ms = 6300.0 * 1000.
   
   radius_shock_neff_1_m = np.sqrt(power_1_neff_J_s / (4. * np.pi * thermal_pressure_ism_neff_Pa))
   radius_shock_neff_1_kpc = (radius_shock_neff_1_m / pc_m) / 1000.
   #print("radius_shock_neff_1_kpc is %E kpc" % radius_shock_neff_1_kpc)

   radius_shock_neff_2_m = np.sqrt(power_2_neff_J_s / (4. * np.pi * thermal_pressure_ism_neff_Pa))
   radius_shock_neff_2_kpc = (radius_shock_neff_2_m / pc_m) / 1000.
   #print("radius_shock_neff_2_kpc is %E kpc" % radius_shock_neff_2_kpc)
   

def image_comparison_feain():
   feain_peak_core_mJy = 1441.68
   feain_rm_1_deg = 17.
   dynamic_range_1_deg =    feain_peak_core_mJy / feain_rm_1_deg
   
   print(dynamic_range_1_deg)
   
   mwa_dynamic_range = 50000.
   improvement_factor = mwa_dynamic_range / dynamic_range_1_deg
   print(improvement_factor)
   
#image_comparison_feain()
#sys.exit()

calculate_outflow_properties()
sys.exit()


#mask_level=0.1
#spectral index ASKAP MWA:
#spectral_index_map('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits','CenA_i.fits',185,1400,'CenA_185_1400_MHz',0.3)
#spectral_index_map('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits','image_cenauc13_image_tt0.fits',185,1400,'CenA_uc13_185_1400_MHz',0.3)
#spectral_index_map("CenA_2015_2018_joint_145_robust0_image_pb_8_ims_0000_08_weighted.fits","CenA_2015_2018_joint_145_robust0_image_pb_8_ims_0009_08_weighted.fits",171.135,198.755,'CenA_171_198_MHz',mask_level)


sys.exit()

#regridding x ray rosat from galacto:
#regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits','g000p00r1b120pm.fits',0.5,0.5,0,'CenA_rosat_north_to_mwa')
#cen_A_rosat_p10_list = ['932527p-p10.fits','932429p-p10.fits']
cen_A_rosat_p10_list = ['932428p-p10.fits','932429p-p10.fits','932430p-p10.fits','932526p-p10.fits','932527p-p10.fits','932528p-p10.fits','932623p-p10.fits','932624p-p10.fits','932625p-p10.fits']
cen_A_rosat_p20_list = ['932428p-p20.fits','932429p-p20.fits','932430p-p20.fits','932526p-p20.fits','932527p-p20.fits','932528p-p20.fits','932623p-p20.fits','932624p-p20.fits','932625p-p20.fits']
cen_A_rosat_p30_list = ['932428p-p30.fits','932429p-p30.fits','932430p-p30.fits','932526p-p30.fits','932527p-p30.fits','932528p-p30.fits','932623p-p30.fits','932624p-p30.fits','932625p-p30.fits']


#regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits',cen_A_rosat_p10_list,0.24,0.24,0,'rosat_low')
#regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits',cen_A_rosat_p20_list,0.24,0.24,0,'rosat_mid')
#regrid_concvol('CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits',cen_A_rosat_p30_list,0.08,0.08,0,'rosat_high_5arcmin')
#sys.exit()




#image_name = "CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits"  
#get_scaling_factor_from_core(image_name,185.,-0.7)

#source_find_image(image_name)

#regrid VLA image to J2000:
#image_name = "NGC_5128_I_21cm_chs1996.fits"
#regrid_to_J2000(image_name)
#sys.exit()

#template_imagename = 'rband_1sec_tr_fl_geo_ha.fits'
template_imagename = 'CenA_optical_template-image.fits'
#template_imagename = 'Ha_cont_subtracted_via_rband_scaling.fits'

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
input_name_list = ['CenA_2015_2018_joint_145_robust0_image_pb_8_ims_08_weighted.fits']
#i know this one works:
#input_name_list = ['CenA_WCS.fits']
#input_name_list = ['CenA_WCS_Ha.fits']
for input_name in input_name_list:
   #edhead_name = input_name.split('.fits')[0]+'_edhead.fits'
   #edit_optical_header(input_name,edhead_name)
   regrid_optical(template_imagename,input_name)    #edhead_name,smooth=2)
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











