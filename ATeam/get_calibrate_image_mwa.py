#!/usr/bin/env python
#does what it says

import os,sys
from astropy.io import fits
import numpy as np

#new ones (close to old from tfranzne): 2015: 1122198832 1122112400 1122025976 1121939544 1121853112 1121766680
#new ones  (close to old from msok) 2018: 1201123376 1201296272 1201382720 1201469160  1201555608 1201814952


obsid_list_2015 = ['1122198832','1122112400','1122025976','1121939544','1121853112','1121766680']
obsid_list_2018 = ['1201123376','1201296272','1201382720','1201469160','1201555608','1201814952']

def download_obs(obsid_list,dest_dir,timeres=4,freqres=40,ms=True):
   #write the csv file
   csv_filename = "manta_ray_download.csv"
   csv_line_list = []
   for obsid in obsid_list:
      if ms==True:
         csv_line = "obs_id=%s, job_type=c, timeres=%s, freqres=%s, edgewidth=80, conversion=ms, calibrate=false, allowmissing=true, flagdcchannels=true " % (obsid,timeres,freqres)
         csv_line_list.append(csv_line)
   
   csv_line_list_string = '\n'.join(csv_line_list)     
   with open(csv_filename,'w') as f:
      f.write(csv_line_list_string)
   
   cmd = "mkdir -p %s" % dest_dir
   print(cmd)
   os.system(cmd)
    
   #run manta ray command
   cmd = "mwa_client -c %s -d %s " % (csv_filename,dest_dir)
   print(cmd)
   os.system(cmd)

def unzip_obs(obsid_list,ms=True):
   for obsid in obsid_list:
      if ms==True:
         cmd = 'unzip %s_ms.zip' % obsid
         print(cmd)
         os.system(cmd)

def create_circular_mask(h, w, centre=None, radius=None):
    #give centre and array of two (x and y centre)
    if centre is None: # use the middle of the image
        centre = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(centre[0], centre[1], w - centre[0], h - centre[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_centre = np.sqrt((X - centre[0])**2 + (Y - centre[1])**2)
    mask = dist_from_centre <= radius
    return mask*1.

def mask_model_image(input_image_name,subtraction_radius_deg=0):
   output_image_name_base = input_image_name.split('.fits')[0].split('/')[-1]
   beam_name = "%s_beam" % output_image_name_base
   
   masked_fits_image_filename = "%s_masked_%0.3f_deg-I.fits" % (output_image_name_base,subtraction_radius_deg)
   
   hdulist = fits.open(input_image_name)
   image_header = hdulist[0].header
   image_data = hdulist[0].data[0,0,:,:]
   
   pixel_scale_deg = float(image_header['CDELT1'])
   imsize = int(image_header['NAXIS1'])
   image_centre = imsize/2.
   
   radius_pix = abs(subtraction_radius_deg/pixel_scale_deg)
   mask = create_circular_mask(imsize,imsize,centre=[image_centre,image_centre],radius=radius_pix)
   masked_image_data = mask * image_data

   fits.writeto(masked_fits_image_filename,masked_image_data,clobber=True)
   fits.update(masked_fits_image_filename,masked_image_data,header=image_header)
   print "wrote  image %s" %  masked_fits_image_filename
   
   return(masked_fits_image_filename)
   



def calibrate_obs(obsid_list,model_image='',generate_new_beams=True,metafits_dir='',ms=True):
   for obsid in obsid_list:
      metafits_filename = "%s/%s.metafits" % (metafits_dir,obsid)
      beam_name = "%s_beam" % obsid
      #get latest with wget?
      if generate_new_beams:
         if os.path.isfile(metafits_filename) and os.access(metafits_filename, os.R_OK):
            cmd = "beam -2016 -proto %s -name %s -m %s " % (model_image,beam_name,metafits_filename)
            print cmd
            os.system(cmd)
         else:
            print("Error: either file %s is missing or is not readable" % metafits_filename)
            sys.exit()
            
      if ms==True:
         #calibrate them
         if model_image!='':
            ms_name = "%s/%s.ms" % (metafits_dir,obsid)
            check_model_image_name = "check_model_%s" % obsid
            check_cal_image_name = "check_cal_%s" % obsid
            
            with fits.open(model_image) as hdulist:
               image_header = hdulist[0].header
               image_data = hdulist[0].data
            
            predict_pixel_scale_deg = float(image_header['CDELT1'])
            predict_imsize = int(image_header['NAXIS1'])
            check_scale_deg = 0.025
            check_imsize = 1024
            
            cal_sol_filename = "%s.bin"
            
            pbuncorrect_input_name = model_image.split('-I.fits')[0]

            uncorrected_image_name = "%s_%s_uncorr" % (pbuncorrect_input_name,obsid)
   
            #uncorrect image with beam
            cmd = "pbcorrect -uncorrect %s model.fits %s %s" % (uncorrected_image_name,beam_name,pbuncorrect_input_name)
            print cmd
            os.system(cmd)

            #predict in the model (wsclean -predict)
            cmd = "wsclean -predict -name %s -size %s %s -scale %.4f -pol xx,xy,yx,yy %s " % (uncorrected_image_name,int(predict_imsize),int(predict_imsize),abs(float(predict_pixel_scale_deg)),ms_name)
            print cmd
            os.system(cmd) 
         
            #image the  model just to see if it worked
            cmd = "wsclean -name %s -size %s %s -niter 0 -data-column MODEL_DATA -scale %s -small-inversion -pol xx  %s" % (check_model_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
            print cmd
            os.system(cmd)
            
            #calibrate
            cmd = "calibrate -minuv 60 %s %s " % (ms_name,cal_sol_filename)
            print cmd
            os.system(cmd)
            
            cmd = "applysolutions %s %s " % (ms_name,cal_sol_filename)
            print cmd
            os.system(cmd)
                     
            #image to see if it worked
            cmd = "wsclean -name %s -size %s %s -niter 0 -data-column CORRECTED_DATA -scale %s -small-inversion -pol xx  %s" % (check_cal_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
            print cmd
            os.system(cmd)
            





#2015 data:
#download_obs(obsid_list_2015,dest_dir='2015',timeres=8,freqres=80,ms=True)
#2018 data:
#download_obs(obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)  


#2015 data:
#unzip_obs(obsid_list_2015,ms=True)
#2018 data:
#unzip_obs(obsid_list_2018,ms=True)  
 
model_image_name = "/md0/ATeam/CenA/CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"
masked_fits_image_filename = mask_model_image(model_image_name,subtraction_radius_deg=4)

calibrate_obs(obsid_list_2015,model_image=masked_fits_image_filename,generate_new_beams=False,metafits_dir="/md0/ATeam/CenA/image_2/2015")

#use uniform weighting and auto thresholds for initial imaging and selfcal then one final robust 0 clean, will probably need to run first, see where it goes non-linear and adjust the niter












