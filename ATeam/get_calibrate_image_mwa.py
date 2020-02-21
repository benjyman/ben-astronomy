#!/usr/bin/env python
#does what it says

import os,sys
from astropy.io import fits
import numpy as np

#new ones (close to old from tfranzne): 2015: 1122198832 1122112400 1122025976 1121939544 1121853112 1121766680
#new ones  (close to old from msok) 2018: 1201123376 1201296272 1201382720 1201469160  1201555608 1201814952

#image_2:
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
   



def calibrate_obs(obsid_list,model_image='',model_wsclean_txt='',generate_new_beams=True,ms_dir='',ms=True):
   for obsid in obsid_list:
      metafits_filename = "%s/%s.metafits" % (ms_dir,obsid)       
      if ms==True:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         check_model_image_name = "check_model_%s" % obsid
         check_cal_image_name = "check_cal_%s" % obsid
         cal_sol_filename = "%s_sol1.bin" % (obsid)
         check_scale_deg = 0.025
         check_imsize = 1024
         #calibrate them
         if model_image!='':
            beam_name = "%s_beam" % obsid
            #get latest metafits with wget?
            if generate_new_beams:
               if os.path.isfile(metafits_filename) and os.access(metafits_filename, os.R_OK):
                  cmd = "beam -2016 -proto %s -name %s -m %s " % (model_image,beam_name,metafits_filename)
                  print(cmd)
                  os.system(cmd)
               else:
                  print("Error: either file %s is missing or is not readable" % metafits_filename)
                  sys.exit()
                  
            with fits.open(model_image) as hdulist:
               image_header = hdulist[0].header
               image_data = hdulist[0].data
            
            predict_pixel_scale_deg = float(image_header['CDELT1'])
            predict_imsize = int(image_header['NAXIS1'])

            
            pbuncorrect_input_name = model_image.split('-I.fits')[0]

            uncorrected_image_name = "%s_%s_uncorr" % (pbuncorrect_input_name,obsid)
            
            #need to change to Jy/pix at some point....
            
            #uncorrect image with beam
            cmd = "pbcorrect -uncorrect %s model.fits %s %s" % (uncorrected_image_name,beam_name,pbuncorrect_input_name)
            print(cmd)
            os.system(cmd)

            #predict in the model (wsclean -predict)
            cmd = "wsclean -predict -name %s -size %s %s -scale %.4f -pol xx,xy,yx,yy %s " % (uncorrected_image_name,int(predict_imsize),int(predict_imsize),abs(float(predict_pixel_scale_deg)),ms_name)
            print(cmd)
            os.system(cmd) 
         
            #image the  model just to see if it worked
            cmd = "wsclean -name %s -size %s %s -niter 0 -data-column MODEL_DATA -scale %s -small-inversion -pol xx  %s" % (check_model_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
            print(cmd)
            os.system(cmd)
            
            #calibrate
            cmd = "calibrate -minuv 60 %s %s " % (ms_name,cal_sol_filename)
            print(cmd)
            os.system(cmd)
            
         #calibrate from txt model output of wsclean
         elif model_wsclean_txt!='':    
            
            #calibrate
            cmd = "calibrate -minuv 60  -applybeam -m %s %s %s " % (model_wsclean_txt,ms_name,cal_sol_filename)
            print(cmd)
            os.system(cmd)            
            
            cmd = "applysolutions %s %s " % (ms_name,cal_sol_filename)
            print(cmd)
            os.system(cmd)
                     
            #image to see if it worked
            cmd = "wsclean -name %s -size %s %s -niter 0 -data-column CORRECTED_DATA -scale %s -small-inversion -pol xx  %s" % (check_cal_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
            print(cmd)
            os.system(cmd)
            

def jointly_deconvolve_idg(obsid_list,ms_dir_list,outname,wsclean_options):
   #wsclean -name  -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight %s -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 %s  " % (ms_string)
   ms_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)

   ms_string = ' '.join(ms_list)
   cmd = "wsclean -name %s %s %s  " % (outname, wsclean_options, ms_string)
   print(cmd)
   os.system(cmd)

def self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number):
   ms_list = []
   solution_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         solution_name = "%s/%s_selfcal_%02d.bin" % (ms_dir,obsid,self_cal_number)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)
            solution_list.append(solution_name)
            
   for ms_index,ms in enumerate(ms_list):
      solution_name = solution_list[ms_index]
      cmd = "calibrate %s %s %s" % (calibrate_options,ms,solution_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "applysolutions %s %s" % (ms,solution_name)
      print(cmd)
      os.system(cmd)
      
      #Plot the cal solutions
      cmd = "aocal_plot.py  %s " % (solution_name)
      print(cmd)
      os.system(cmd)
   
def flag_bad_ants(ms_name,ant_string):
   cmd = "flagantennae %s %s " % (ms_name,ant_string)
   print(cmd)
   os.system(cmd)

def ms_qa(obsid_list,ms_dir_list):
   ms_list = []
   aoqplot_savename_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         aoqplot_savename = "aoqplot_stddev_%s"
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)
            aoqplot_savename_list.append(aoqplot_savename)
            
   for ms_name_index,ms_name in enumerate(ms_list):
      aoqplot_savename = aoqplot_savename_list[ms_name_index]
      cmd = "aoquality collect %s; aoqplot -save %s StandardDeviation %s" % (ms_name,aoqplot_savename,ms_name)
      print(cmd)
      os.system(cmd)    

def average_images(base_dir,image_list,output_name_base):
   n_images = len(image_list)
   output_name = "%s_%02d.fits" % (output_name_base,n_images)
   for image_name_index,image_name in enumerate(image_list):
      hdulist = fits.open("%s%s" % (base_dir,image_name))
      image_header = hdulist[0].header
      image_data = hdulist[0].data[0,0,:,:]
      
      if image_name_index==0:
         image_data_sum = image_data*0.
         image_data_sum += image_data
      else:
         image_data_sum += image_data
   image_data_average = image_data_sum / float(n_images)
   fits.writeto(output_name,image_data_average,clobber=True)
   fits.update(output_name,image_data_average,header=image_header)
   print "wrote  image %s" %  output_name

      


#2015 data:
#download_obs(obsid_list_2015,dest_dir='2015',timeres=8,freqres=80,ms=True)
#2018 data:
#download_obs(obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)  


#2015 data:
#unzip_obs(obsid_list_2015,ms=True)
#2018 data:
#unzip_obs(obsid_list_2018,ms=True)  

#Do this first!!!!!
#ms_qa(obsid_list,ms_dir_list)
 
model_image_name = "/md0/ATeam/CenA/CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"
#masked_fits_image_filename = mask_model_image(model_image_name,subtraction_radius_deg=4)

#masked_fits_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb_masked_4.000_deg-I.fits"
#calibrate_obs(obsid_list_2015,model_image=masked_fits_image_filename,generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2015")
#calibrate_obs(obsid_list_2018,model_image=masked_fits_image_filename,generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2018")

#The above strategy doesn't work cause you are using a deconvolved image in Jy/beam!
#Use the same strategy as before with the first good 145 image i.e. the .txt model 
#calibrate_obs(obsid_list_2015,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2015")
#calibrate_obs(obsid_list_2018,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2018")


#use uniform weighting and auto thresholds for initial imaging and selfcal then one final robust 0 clean, will probably need to run first, see where it goes non-linear and adjust the niter
obsid_list = obsid_list_2015 + obsid_list_2018
ms_dir_list=["/md0/ATeam/CenA/image_2/2015","/md0/ATeam/CenA/image_2/2018"]

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_01_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=1)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=2)

#looks like there is something fishing going on, bad artefacts in middle of image. Looks like there is a bad reciever at least for 1122198832
#look closely at cal plots and also with aoquality collect / aoqplot
#flag_bad_ants("2015/1122198832.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1121766680.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1121853112.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1121939544.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1122025976.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1122112400.ms","49 56 60 104 105 106 107 108 109 110 111")
#flag_bad_ants("2015/1122198832.ms","49 56 60 104 105 106 107 108 109 110 111")
#potentially ant 49 and 10x ants on all 2015 obs

#now 2018
#flag_bad_ants("2018/1201123376.ms","84 76")
#flag_bad_ants("2018/1201296272.ms","84 76")
#flag_bad_ants("2018/1201382720.ms","84 76")
#flag_bad_ants("2018/1201469160.ms","84 76")
#flag_bad_ants("2018/1201555608.ms","84 76")
#flag_bad_ants("2018/1201814952.ms","84 76")

#should probably go back and start again from the start now the flagging is done properly

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_model_cal_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=1)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_01_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=2)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=3)


#looking good, but still some artefacts compared to first image - one more selfcal
#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_image2"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=4)

#wsclean -name CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2
#check output of last imaging and adjust threshold an niter accordingly

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image2"
#wsclean_options = "-size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
#based on first cleaning attempt above (goes non-linear)
#wsclean_options = "-size 4096 4096 -niter 350000  -threshold 0.015  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)

base_dir = "/md0/ATeam/CenA/"
image_list = ["image_2/CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image2-MFS-image-pb.fits","CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"]
output_name_base = "CenA_2015_2018_joint_145_robust0_image_pb"
average_images(base_dir=base_dir,image_list=image_list,output_name_base=output_name_base)

#looks good! now need more images to average!!! maybe re-image old data with uniform weighting and self-cal, then redo robust 0 cleaning (twice) this is because wsclean has been updated a lot since that image was made (new image looks better for sure)






