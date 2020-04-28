#!/usr/bin/env python
#does what it says

import os,sys
from astropy.io import fits
import numpy as np

#new ones (close to old from tfranzne): 2015: 1122198832 1122112400 1122025976 1121939544 1121853112 1121766680
#new ones  (close to old from msok) 2018: 1201123376 1201296272 1201382720 1201469160  1201555608 1201814952

#image_2:
#obsid_list_2015 = ['1122198832','1122112400','1122025976','1121939544','1121853112','1121766680']
#obsid_list_2018 = ['1201123376','1201296272','1201382720','1201469160','1201555608','1201814952']

#image_1 (original)
#obsid_list_2015 = ['1117031728','1121334536','1121420968','1121507392','1121593824','1121680256']
#obsid_list_2018 = ['1199663088','1200604688','1200691136','1200777584','1200864032','1200950480']

#image_3
#obsid_list_2015 = ['1122285264','1122371696','1122458120','1122544552','1122630984','1122717416']
#obsid_list_2018 = ['1201895248','1201900584','1201901392','1201981408','1201986752','1201987840']

#image_4
#obsid_list_2015 = ['1112806040','1112892200','1114782984','1114869144','1114955312','1115041472']
#obsid_list_2018 = ['1202239904','1202326064','1202410528','1202411608','1202418952','1202672864']

#image 6
#obsid_list_2015 = ['1115307768','1116343632','1116429792','1116603776','1116604056','1116773232']
#obsid_list_2018 = ['1202679384','1202756888','1202764984','1202843048','1202851152','1203015384']

#image 7
#the 2018/2019 data are actuall pointed off CenA (credit Natasha), use chgcentre....
#obsid_list_2015 = ['1120470256','1120815968','1120902392','1120988824','1121075248','1121161680']
#obsid_list_2018 = ['1244807072','1244376256','1243859272','1243341384','1243169056','1242739136'] #actually 2019


#image 8:
#the 2018/2019 data are actuall pointed off CenA (credit Natasha), use chgcentre....
obsid_list_2015 = ['1121248104','1121334536','1121420968','1121507392','1121593824','1121680256']
obsid_list_2018 = ['1242479744','1242221248','1242049824','1241705168','1236791704','1234292944'] #actually 2019
      
#image 11:
#obsid_list_2015 = ['1086347312','1086355112','1086433480','1086441280','1086519640','1086527440']
#obsid_list_2018 = ['1200604688','1200691136','1200777584','1200864032','1200950480','1201895248']

#image 12
#obsid_list_2015 = ['1086605808','1086613608','1086691968','1086699768','1086778136','1086785936']
#obsid_list_2018 = ['1201123376','1201296272','1201382720','1201469160','1201555608','1201814952']

def download_obs(obsid_list,dest_dir,timeres=4,freqres=40,ms=True):
   #write the csv file
   csv_filename = "%s/manta_ray_download.csv" % dest_dir
   csv_line_list = []
   
   cmd = "mkdir -p %s" % dest_dir
   print(cmd)
   os.system(cmd)
   
   for obsid in obsid_list:
      if ms==True:
         csv_line = "obs_id=%s, job_type=c, timeres=%s, freqres=%s, edgewidth=80, conversion=ms, calibrate=false, allowmissing=true, flagdcchannels=true " % (obsid,timeres,freqres)
         csv_line_list.append(csv_line)
   
   csv_line_list_string = '\n'.join(csv_line_list)     
   with open(csv_filename,'w') as f:
      f.write(csv_line_list_string)
   
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
            
            #do aocal_plot.py
            #Plot the cal solutions
            cmd = "aocal_plot.py  %s " % (cal_sol_filename)
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

def chgcentre_ms(obsid_list,ms_dir_list,target="CenA"):
   ms_list = []
   if target=="CenA":
      new_ra = "13h25m27.6152s"
      new_dec = "-43d01m08.805s"
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)
   for ms_index,ms in enumerate(ms_list):
      cmd = "chgcentre %s %s %s" % (ms,new_ra,new_dec)
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

def ms_qa(obsid_list,ms_dir_list,data_column="DATA"):
   ms_list = []
   aoqplot_savename_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         aoqplot_savename = "aoqplot_stddev_%s_%s" % (obsid,data_column)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)
            aoqplot_savename_list.append(aoqplot_savename)
            
   for ms_name_index,ms_name in enumerate(ms_list):
      aoqplot_savename = aoqplot_savename_list[ms_name_index]
      cmd = "aoquality collect -d %s %s; aoqplot -save %s StandardDeviation %s" % (data_column,ms_name,aoqplot_savename,ms_name)
      print(cmd)
      os.system(cmd)    

def average_images(base_dir,image_list,output_name_base,weighted_average=True):
   n_images = len(image_list)
   output_name = "%s_%02d.fits" % (output_name_base,n_images)
   if weighted_average==True:
      output_name_weighted = "%s_%02d_weighted.fits" % (output_name_base,n_images)
   
   #convolve the images to the same resolution first (miriad convol function)
   bmaj_list = []
   bmin_list = []
   bpa_list = []
   for image_name_index,image_name in enumerate(image_list):
      hdulist = fits.open("%s%s" % (base_dir,image_name))
      image_header = hdulist[0].header
      image_data = hdulist[0].data[0,0,:,:] 
      bmaj = float(image_header['BMAJ'])
      bmin = float(image_header['BMIN'])
      bpa = float(image_header['BPA'])
      bmaj_list.append(bmaj)
      bmin_list.append(bmin)
      bpa_list.append(bpa)

   #find the index of the biggest bpa
   max_ind = np.argmax(bmaj_list)
   bmaj_convol = bmaj_list[max_ind]
   bmin_convol = bmin_list[max_ind]
   bpa_convol = bpa_list[max_ind]
   
   
   max_flux_list = []
   weights_list = []
   for image_name_index,image_name in enumerate(image_list):
      image_name_base = image_name.split('.fits')[0]
      
      #read in image to miriad
      im_name = "%s.im" % image_name_base
      convol_imname = "%s_convol.im" % image_name_base
      convol_fitsname = "%s_convol.fits" % image_name_base
      
      cmd = "rm -rf %s %s %s" % (im_name,convol_imname,convol_fitsname)
      print(cmd)
      os.system(cmd)
      
      cmd = "fits in=%s out=%s op=xyin" % (image_name,im_name)
      print(cmd)
      os.system(cmd)
      
      cmd = "convol map=%s fwhm=%4f,%4f pa=%4f out=%s " % (im_name,bmaj_convol,bmin_convol,bpa_convol,convol_imname)
      print(cmd)
      os.system(cmd)
      
      
      cmd = "fits in=%s out=%s op=xyout" % (convol_imname,convol_fitsname)
      print(cmd)
      os.system(cmd)
      
      hdulist = fits.open("%s" % (convol_fitsname))
      image_header = hdulist[0].header
      image_data = hdulist[0].data[0,0,:,:]
      
      #need all images to be on same flux scale (have same max flux)
      data_max = np.max(image_data)
      max_flux_list.append(data_max)
      
      #choose a region to measure the rms in for weighting. 
      region_x_start,region_x_end,region_y_start,region_y_end = 2662,2778,1387,1491
      data_region = image_data[region_y_start:region_y_end,region_x_start:region_x_end]
      
      #write out the image region as fits just to check
      fits.writeto('region_check.fits',data_region,clobber=True)
      
      if weighted_average==True:
         #get rms of region
         rms = np.sqrt(np.mean(np.square(data_region)))
         variance_weight = 1./rms**2
         weights_list.append(variance_weight)
         #print "rms %3f" % rms
      
   
      
   scale_max_flux = np.max(max_flux_list)
   
   for image_name_index,image_name in enumerate(image_list):
      image_name_base = image_name.split('.fits')[0]
      image_name_scaled = "%s_scaled.fits" % image_name_base
      
      hdulist = fits.open("%s" % (image_name))
      image_header = hdulist[0].header
      image_data = hdulist[0].data[0,0,:,:]
   
      #work out scaling ratio:
      image_data_max = np.max(image_data)
      scaling = scale_max_flux / image_data_max
      scaled_data = scaling * image_data
      
      #write out scaled images to check
      fits.writeto(image_name_scaled,scaled_data,clobber=True)
      fits.update(image_name_scaled,scaled_data,header=image_header)
      print("wrote  image %s" %  image_name_scaled)
      
      #also do the weighted average
      if weighted_average==True:
         variance_weight = weights_list[image_name_index]
         print "weighting image by %s" % (variance_weight)
         
      if image_name_index==0:
         image_data_sum = scaled_data*0.
         image_data_sum += scaled_data
         if weighted_average==True:
            image_data_sum_weighted = scaled_data*0.
            image_data_sum_weighted += scaled_data*variance_weight
      else:
         image_data_sum += scaled_data
         if weighted_average==True:
            image_data_sum_weighted += scaled_data*variance_weight
         
         
   image_data_average = image_data_sum / float(n_images)
   fits.writeto(output_name,image_data_average,clobber=True)
   fits.update(output_name,image_data_average,header=image_header)
   print("wrote  image %s" %  output_name)
   

   if weighted_average==True:
      sum_of_weights = np.sum(weights_list)
      image_data_average_weighted = image_data_sum_weighted / sum_of_weights
      fits.writeto(output_name_weighted,image_data_average_weighted,clobber=True)
      fits.update(output_name_weighted,image_data_average_weighted,header=image_header)
      print("wrote  image %s" %  output_name_weighted)
      



#run in separate screens
#2015 data:
#download_obs(obsid_list_2015,dest_dir='2015',timeres=8,freqres=80,ms=True)
#2018 data:
#download_obs(obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)  
#sys.exit()

#2015 data:
#unzip_obs(obsid_list_2015,ms=True)
#2018 data:
#unzip_obs(obsid_list_2018,ms=True)  
#sys.exit()

#can do QA on namorodror by mounting /fred/oz048/bmckinle to namorrodor
#Do this first!!!!!
#the repeat after calibration
#ms_qa(obsid_list_2018,['/md0/ATeam/CenA/image_4/2018'],data_column="DATA")
#ms_qa(obsid_list_2018,['/md0/ozstar/ATeam/CenA/image4/2018'],data_column="DATA")
#qa image_4 
#ms_qa(obsid_list_2015,['/md0/ATeam/CenA/image_4/2015'],data_column="DATA")
#ms_qa(obsid_list_2015,['/md0/ozstar/ATeam/CenA/image4/2015'],data_column="DATA")
#sys.exit()


#flagged all image_3 49 (2015)
#flagged 76 80 for image_3 2018
#flag_bad_ants("2015/1122544552.ms","49 56")
#flag_bad_ants("2015/1122458120.ms","49 56")
#sys.exit()

#model_image_name = "/md0/ATeam/CenA/CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"
#masked_fits_image_filename = mask_model_image(model_image_name,subtraction_radius_deg=4)

####Dont use:!!!!!
#####masked_fits_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb_masked_4.000_deg-I.fits"
#####calibrate_obs(obsid_list_2015,model_image=masked_fits_image_filename,generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2015")
#####calibrate_obs(obsid_list_2018,model_image=masked_fits_image_filename,generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_2/2018")

#yes:
#The above strategy doesn't work cause you are using a deconvolved image in Jy/beam!
#Use the same strategy as before with the first good 145 image i.e. the .txt model 
#calibrate_obs(obsid_list_2015,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_7/2015")
#calibrate_obs(obsid_list_2018,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_7/2018")

#ms_qa(obsid_list_2015,['/md0/ATeam/CenA/image_7/2015'],data_column="CORRECTED_DATA")
#ms_qa(obsid_list_2018,['/md0/ATeam/CenA/image_7/2018'],data_column="CORRECTED_DATA")

#sys.exit()

#chgcentre for obs pointed off cena (image 7 and 8)
#chgcentre_ms(obsid_list_2015,['/md0/ATeam/CenA/image_8/2015'],target='CenA')
#chgcentre_ms(obsid_list_2018,['/md0/ATeam/CenA/image_8/2018'],target='CenA')

#sys.exit()


#use uniform weighting and auto thresholds for initial imaging and selfcal then one final robust 0 clean, will probably need to run first, see where it goes non-linear and adjust the niter
obsid_list = obsid_list_2015 + obsid_list_2018
#ms_dir_list=["/md0/ATeam/CenA/image_2/2015","/md0/ATeam/CenA/image_2/2018"]
#image_3:
ms_dir_list=["/md0/ATeam/CenA/image_3/2015","/md0/ATeam/CenA/image_3/2018"]

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

#base_dir = "/md0/ATeam/CenA/"
#image_list = ["image_2/CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image2-MFS-image-pb.fits","CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"]
#output_name_base = "CenA_2015_2018_joint_145_robust0_image_pb"
#average_images(base_dir=base_dir,image_list=image_list,output_name_base=output_name_base)

#looks good! now need more images to average!!! maybe re-image old data with uniform weighting and self-cal, then redo robust 0 cleaning (twice) this is because wsclean has been updated a lot since that image was made (new image looks better for sure)
#image_1_redo

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_image1"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=4)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image1"
#wsclean_options = "-size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
#based on first cleaning attempt above (goes non-linear)
#wsclean_options = "-size 4096 4096 -niter 350000  -threshold 0.015  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)



#In the meantime get new data for image 3

#calibrate image_3
#calibrate_obs(obsid_list_2015,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_3/2015")
#calibrate_obs(obsid_list_2018,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_3/2018")

#image 4
#calibrate_obs(obsid_list_2015,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_4/2015")
#calibrate_obs(obsid_list_2018,model_wsclean_txt='/md0/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',generate_new_beams=True,ms_dir="/md0/ATeam/CenA/image_4/2018")
#sys.exit()


#the repeat after calibration

#ms_qa(obsid_list_2018,['/md0/ATeam/CenA/image_3/2018'],data_column="CORRECTED_DATA")
#ms_qa(obsid_list_2015,['/md0/ATeam/CenA/image_3/2015'],data_column="CORRECTED_DATA")
#sys.exit()

#the repeat after calibration
#ms_qa(obsid_list_2015,['/md0/ATeam/CenA/image_4/2015'],data_column="CORRECTED_DATA")
#sys.exit()
#WHOA some whacky stuff going on with  1112806040 - do some serious flagging and re-cal



######################################################################
#Need to move all this stuff to ozstar /fred/oz048/bmckinle/ATeam
##image 3:
#go ahead with imaging and self cal
#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_model_cal_image3"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=1)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_01_image3"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=2)

#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_image3"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=3)

#looking good, but still some artefacts compared to first image - one more selfcal
#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_image3"
#wsclean_options = " -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
#calibrate_options = "-minuv 60"
#self_cal(obsid_list,ms_dir_list,calibrate_options,self_cal_number=4)


#image_outname = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image3"
###wsclean_options = "-size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
###based on first cleaning attempt above (goes non-linear)
#wsclean_options = "-size 4096 4096 -niter 350000  -threshold 0.015  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2"
#jointly_deconvolve_idg(obsid_list=obsid_list,ms_dir_list=ms_dir_list,outname=image_outname,wsclean_options=wsclean_options)
######################################################################



#############
#averaging images bit (keep at bottom)
base_dir = "/md0/ATeam/CenA/"
image_list = ["CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image1-MFS-image-pb.fits","image_2/CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image2-MFS-image-pb.fits","image_3/CenA_2015_2018_joint_idg_12_obs_145_selfcal_04_robust0_image3-MFS-image-pb.fits","image_4/image_04_selfcal_04_robust0-MFS-image-pb.fits","image_5/image_05_selfcal_05_robust0-MFS-image-pb.fits","image_6/image_06_selfcal_04_robust0-MFS-image-pb.fits","image_7/image_07_selfcal_05_robust0-MFS-image-pb.fits","image_8/image_08_selfcal_05_robust0-MFS-image-pb.fits"]
output_name_base = "CenA_2015_2018_joint_145_robust0_image_pb_8_ims"
average_images(base_dir=base_dir,image_list=image_list,output_name_base=output_name_base)






