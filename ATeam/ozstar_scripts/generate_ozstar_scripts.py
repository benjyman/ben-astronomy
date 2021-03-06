#!/usr/bin/env python
#generate scripts to lauch on ozstar for data processing

import os,sys
from astropy.io import fits
import numpy as np

def initiate_script(filename,time_hrs,partition_string='skylake',array_string=''):
   if partition_string=='skylake-gpu':
      gpu_string = '#SBATCH --gres=gpu:1' 
   else: 
      gpu_string = ''
   header_text = """#!/bin/bash -l
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=100000 
#SBATCH --account=oz048 
#SBATCH --time=%02d:00:00 
#SBATCH --partition=%s 
%s
%s

module load gcc/6.4.0 openmpi/3.0.0 
module load fftw/3.3.7 
module load gsl/2.4 
module load cfitsio/3.420 
module load boost/1.67.0-python-2.7.14 
module load hdf5/1.10.1 
module load openblas/0.2.20 
module load cuda/9.0.176 


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/MWA/CODE/lib/

"""        % (time_hrs,partition_string,array_string,gpu_string)
   with open(filename,'w') as f:
      f.write(header_text)

        
def generate_download(obsid_list,dest_dir,timeres=4,freqres=40,ms=True):
   print('generating download script') 
   
   cmd = "mkdir -p %s;\ncd %s" % (dest_dir,dest_dir)
   print(cmd)
   os.system(cmd)
   
   out_filename = '%s/download_obs.sh' % dest_dir
   initiate_script(out_filename,time_hrs=12)
   
   #write the csv file
   csv_filename = "%s/manta_ray_download.csv" % dest_dir
   csv_line_list = []
   for obsid in obsid_list:
      if ms==True:
         csv_line = "obs_id=%s, job_type=c, timeres=%s, freqres=%s, edgewidth=80, conversion=ms, calibrate=false, allowmissing=true, flagdcchannels=true " % (obsid,timeres,freqres)
         csv_line_list.append(csv_line)
   
   csv_line_list_string = '\n'.join(csv_line_list)     
   with open(csv_filename,'w') as f:
      f.write(csv_line_list_string)
   
   #run manta ray command
   cmd = "time python /fred/oz048/MWA/CODE/local_python/mantaray/scripts/mwa_client.py -c %s -d %s " % (csv_filename,dest_dir)
   with open(out_filename,'a') as f:
      f.write(cmd)
   
   print("wrote %s" % (out_filename))
   return(out_filename)
   
   
def generate_unzip(obsid_list,dest_dir):
   out_filename = '%s/unzip_obs.sh' % dest_dir
   initiate_script(out_filename,time_hrs=2) 
   
   cmd = "cd %s;\nunzip -o '*.zip'\nrm *.zip " % dest_dir
   with open(out_filename,'a') as f:
      f.write(cmd)
   
   print("wrote %s" % (out_filename))
   return(out_filename)
   
def generate_model_cal(obsid_list,model_wsclean_txt='',dest_dir=''):
   #individual job script for each obsid
   print('generating model calibration script')
   out_filename_list = []
   for obsid in obsid_list:
      out_filename = '%s/model_cal_%s.sh' % (dest_dir,obsid)
      initiate_script(out_filename,time_hrs=2) 
   
      cmd = "cd %s\n" % dest_dir
      with open(out_filename,'a') as f:
         f.write(cmd)
      
      metafits_filename = "%s/%s.metafits" % (dest_dir,obsid)       
      ms_name = "%s/%s.ms" % (dest_dir,obsid)
      check_model_image_name = "check_model_%s" % obsid
      check_cal_image_name = "check_cal_%s" % obsid
      cal_sol_filename = "%s_sol1.bin" % (obsid)
      check_scale_deg = 0.025
      check_imsize = 1024
      #calibrate them
      cmd_list = []
      
      
      cmd = "module load boost/1.66.0-python-2.7.14 \n"
      cmd_list.append(cmd)
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/calibrate -minuv 60  -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -applybeam -m %s %s %s \n" % (model_wsclean_txt,ms_name,cal_sol_filename)
      cmd_list.append(cmd)
      
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/applysolutions %s %s \n" % (ms_name,cal_sol_filename)
      cmd_list.append(cmd)
               
      #image to see if it worked
      cmd = "module load boost/1.67.0-python-2.7.14 \n"
      cmd_list.append(cmd)
      cmd = "/fred/oz048/MWA/CODE/bin/wsclean -name %s -size %s %s -niter 0 -data-column CORRECTED_DATA -scale %s -small-inversion -pol xx %s \n" % (check_cal_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
      cmd_list.append(cmd)
      
      #do aocal_plot.py
      #Plot the cal solutions
      #cmd = "aocal_plot.py  %s \n" % (cal_sol_filename)
      #cmd_list.append(cmd)

      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      out_filename_list.append(out_filename)
      print("wrote %s" % (out_filename))
      
   return(out_filename_list)
   
def generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base,wsclean_options,dest_dir,self_cal_number=0):
   print('generating wsclean script')
   out_filename = '%s/wsclean_selfcal_%02d.sh' % (dest_dir,self_cal_number)
   initiate_script(out_filename,time_hrs=20,partition_string='skylake-gpu')
   ms_list = []
   cmd_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)

   ms_string = ' '.join(ms_list)
   cmd = "time /fred/oz048/MWA/CODE/bin/wsclean -name %s %s %s  \n" % (out_image_name_base, wsclean_options, ms_string)
   cmd_list.append(cmd)
   
   with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
   print("wrote %s" % (out_filename))
   return(out_filename)
   
def generate_selfcal(obsid_list,ms_dir_list,calibrate_options,self_cal_number,dest_dir):
   print('generating selfcal scripts')

   ms_list = []
   solution_list = []
   for ms_dir in ms_dir_list:
      for obsid in obsid_list:
         ms_name = "%s/%s.ms" % (ms_dir,obsid)
         solution_name = "%s/%s_selfcal_%02d.bin" % (ms_dir,obsid,self_cal_number)
         if os.path.isdir(ms_name):
            ms_list.append(ms_name)
            solution_list.append(solution_name)
   out_filename_list = []        
   for ms_index,ms in enumerate(ms_list):
      out_filename = '%s/selfcal_%02d_ms_%s.sh' % (dest_dir,self_cal_number,ms_index)
      initiate_script(out_filename,time_hrs=2,partition_string='skylake')
      cmd_list = []
      cmd = "module load boost/1.66.0-python-2.7.14 \n"
      cmd_list.append(cmd)
   
      solution_name = solution_list[ms_index]
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/calibrate %s %s %s \n" % (calibrate_options,ms,solution_name)
      cmd_list.append(cmd)
      
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/applysolutions %s %s \n" % (ms,solution_name)
      cmd_list.append(cmd)
      
      #Plot the cal solutions ##need to get this working
      #cmd = "aocal_plot.py  %s \n" % (solution_name)
      #cmd_list.append(cmd)
      
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      out_filename_list.append(out_filename)
   return(out_filename_list)
   
def generate_sbatch_script_CenA(image_number,n_selfcals,download=False,model_cal=False):
   print('generating sbatch script for CenA')
   obsid_list_2015,obsid_list_2018,ms_dir_list = get_obsid_list(image_number)
   obsid_list = obsid_list_2015 + obsid_list_2018
   if download:
      out_filename = 'sbatch_launch_download_image_%02d.sh' % (image_number)
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
   
      out_filename_2015 = generate_download(obsid_list=obsid_list_2015,dest_dir='2015',timeres=4,freqres=40,ms=True)
      out_filename_2018 = generate_download(obsid_list=obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)   
      
      # first commands no dependencies
      cmd = 'jid1=$(sbatch %s | cut -f 4 -d " ") \n' % out_filename_2015
      cmd_list.append(cmd)
      cmd = 'jid2=$(sbatch %s | cut -f 4 -d " ") \n' % out_filename_2018
      cmd_list.append(cmd)
      
      #unzip
      out_filename_2015 = generate_unzip(obsid_list=obsid_list_2015,dest_dir='2015')
      out_filename_2018 = generate_unzip(obsid_list=obsid_list_2018,dest_dir='2018')  
      
      #unzip commands depend on downloads
      cmd = 'jid3=$(sbatch --dependency=afterok:$jid1 %s | cut -f 4 -d " ") \n' % out_filename_2015
      cmd_list.append(cmd)
      cmd = 'jid4=$(sbatch --dependency=afterok:$jid2 %s | cut -f 4 -d " ") \n' % out_filename_2018
      cmd_list.append(cmd)
   
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      
      cmd = 'chmod +x %s' % out_filename
      print(cmd)
      os.system(cmd)
      
      return(out_filename)
      
   #do qa manually after download and unzip then proceed
   elif model_cal:
      out_filename = 'sbatch_launch_model_cal_image_%02d.sh' % (image_number)
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
   
      model_cal_out_filename_2015_list = generate_model_cal(obsid_list_2015,model_wsclean_txt='/fred/oz048/bmckinle/code/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[0])
      model_cal_out_filename_2018_list  = generate_model_cal(obsid_list_2018,model_wsclean_txt='/fred/oz048/bmckinle/code/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[1])
      
      model_cal_out_filename_list = model_cal_out_filename_2015_list + model_cal_out_filename_2018_list
       
      for model_cal_out_filename_index,model_cal_out_filename in enumerate(model_cal_out_filename_list):
         cmd = 'jid%s=$(sbatch %s | cut -f 4 -d " ") \n' % (model_cal_out_filename_index,model_cal_out_filename)
         cmd_list.append(cmd)

      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      
      cmd = 'chmod +x %s' % out_filename
      print(cmd)
      os.system(cmd)
      
      return(out_filename)

   #do qa again after model cal then proceed
   else:
      out_filename = 'sbatch_launch_wsclean_and_selfcal_image_%02d.sh' % (image_number)
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
      job_id_list = []
      
      for selfcal in range(0,n_selfcals,1):
         wsclean_job_index = int(selfcal + len(job_id_list)*selfcal)
         wsclean_options_1 = " -size 4096 4096 -j 8 -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam -idg-mode hybrid -pb-undersampling 4 -channels-out 8 -join-channels -fit-spectral-pol 2"
         calibrate_options_1 = "-minuv 60"
         
         out_image_name_base = "image_%02d_selfcal_%02d_uniform" % (image_number,selfcal)
         
         wsclean_out_filename = generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base=out_image_name_base,wsclean_options=wsclean_options_1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image%s' % image_number,self_cal_number=selfcal)
         
         if selfcal==0:
            cmd = 'jid%s=$(sbatch %s | cut -f 4 -d " ") \n' % (wsclean_job_index,wsclean_out_filename)
         else:
            cmd = 'jid%s=$(sbatch --dependency=afterok:%s %s | cut -f 4 -d " ") \n' % (wsclean_job_index,job_id_list_string,wsclean_out_filename)
         cmd_list.append(cmd)
         
         selfcal_out_filename_list = generate_selfcal(obsid_list,ms_dir_list,calibrate_options=calibrate_options_1,self_cal_number=selfcal+1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image%s' % image_number)
         
         job_id_list = []
         for selfcal_out_filename_index,selfcal_out_filename in enumerate(selfcal_out_filename_list):
            selfcal_job_index = selfcal*len(selfcal_out_filename_list) + selfcal_out_filename_index + 1 + selfcal
            cmd = 'jid%s=$(sbatch --dependency=afterok:$jid%s %s | cut -f 4 -d " ") \n' % (selfcal_job_index,wsclean_job_index,selfcal_out_filename)
            cmd_list.append(cmd)
            job_id_list.append('$jid%s' % selfcal_job_index)
         
         job_id_list_string = ':'.join(job_id_list)
         
      #final robust 0 clean
      out_image_name_base = "image_%02d_selfcal_%02d_robust0" % (image_number,selfcal+1)
      
      wsclean_options_2 = "-size 4096 4096 -j 8 -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam -idg-mode hybrid -pb-undersampling 4 -channels-out 10 -join-channels -fit-spectral-pol 2"
      wsclean_out_filename = generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base=out_image_name_base,wsclean_options=wsclean_options_2,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image%s' % image_number,self_cal_number=n_selfcals)
      cmd = 'id%s=$(sbatch --dependency=afterok:$jid%s %s | cut -f 4 -d " ") \n' % (n_selfcals,job_id_list_string,wsclean_out_filename)
      cmd_list.append(cmd)
    
      #FIRST=$(sbatch testrun.sh | cut -f 4 -d' ')
      #echo $FIRST
      #SECOND=$(sbatch -d afterany:$FIRST testrun.sh | cut -f 4 -d' ')
      #echo $SECOND
      #THIRD=$(sbatch -d afterany:$SECOND testrun.sh | cut -f 4 -d' ')
        
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      
      cmd = 'chmod +x %s' % out_filename
      print(cmd)
      os.system(cmd)

      return(out_filename)
   
def get_obsid_list(image_number):
   #best results are from images 4 and 5 where there are only 2 lsts in 205 data - so maybe it performs better if there are less (single?) lsts, 
   if image_number==1:
      #these ones are bad (from corrected_data qa): '1117031728' , '1199663088'
      lst_list_2015   = ['13.381','13.455','13.528','13.602','13.677']
      lst_list_2018   = ['13.154','13.234','13.313','13.392','13.471']
      obsid_list_2015 = ['1121334536','1121420968','1121507392','1121593824','1121680256']
      obsid_list_2018 = ['1200604688','1200691136','1200777584','1200864032','1200950480']
   elif image_number==2:
      #try removing 1201123376 1201296272  - some weirdness in phase sols
      lst_list_2015   = ['14.12','14.046','13.973','13.899','13.824','13.749']
      lst_list_2018   = ['13.629','13.787','13.866','13.943','14.022','14.259']
      obsid_list_2015 = ['1122198832','1122112400','1122025976','1121939544','1121853112','1121766680']
      obsid_list_2018 = ['1201123376','1201296272','1201382720','1201469160','1201555608','1201814952']
   elif image_number==3:
      #try removing 1201895248 - some weirdness in phase sols
      lst_list_2015   = ['14.195','14.27','14.342','14.417','14.491','14.566']
      lst_list_2018   = ['12.625','14.111','14.336','12.624','14.112','14.415']
      obsid_list_2015 = ['1122285264','1122371696','1122458120','1122544552','1122630984','1122717416']
      obsid_list_2018 = ['1201895248','1201900584','1201901392','1201981408','1201986752','1201987840']
   elif image_number==4:
      #leave out 1112806040, just going to make things worse
      #leave out 1112892200 (lst '13.867') as well, to see if just having one lst of 2015 data improves the image
      #including or leaving out 1112892200 made little difference, try using 2018 data all from the same LST
      lst_list_2015   = ['12.523','12.521','12.522','12.521','13.867']
      lst_list_2018   = ['12.625','12.624','12.15','12.451','14.496','13.221']
      obsid_list_2015 = ['1114782984','1114869144','1114955312','1115041472','1112892200']
      obsid_list_2018 = ['1202239904','1202326064','1202410528','1202411608','1202418952','1202672864']
   elif image_number==5:
      lst_list_2015   = ['14.694','12.522','14.695','12.521','14.694','12.522']
      lst_list_2018   = ['13.314','13.381','13.448','13.515','13.582','14.783']
      obsid_list_2015 = ['1115049272','1115127640','1115135440','1115213800','1115221600','1115299968']
      obsid_list_2018 = ['1202673200','1202673440','1202673680','1202673920','1202674160','1202678472']
   elif image_number==6:
      lst_list_2015   = ['14.695','15.223','15.222','15.683','15.761','14.883']
      lst_list_2018   = ['15.037','12.625','14.88','12.623','14.881','12.626']
      obsid_list_2015 = ['1115307768','1116343632','1116429792','1116603776','1116604056','1116773232']
      obsid_list_2018 = ['1202679384','1202756888','1202764984','1202843048','1202851152','1203015384']
   elif image_number==7:
      lst_list_2015   = ['12.646','12.94','13.012','13.087','13.159','13.234']
      lst_list_2018   = ['13.212','13.213','13.213','12.961','12.961','13.212']
      obsid_list_2015 = ['1120470256','1120815968','1120902392','1120988824','1121075248','1121161680']
      obsid_list_2018 = ['1244807072','1244376256','1243859272','1243341384','1243169056','1242739136'] #actually 2019
   elif image_number==8:
      #mabe leave ouut 1234292944.ms  ... higher rms in corrected data, looks a bit dodgy,'1234292944', ,'13.677'
      lst_list_2015   = ['13.306','13.381','13.455','13.528','13.602']
      lst_list_2018   = ['12.962','12.961','13.212','13.213','12.625','12.624']
      obsid_list_2015 = ['1121248104','1121334536','1121420968','1121507392','1121593824','1121680256']
      obsid_list_2018 = ['1242479744','1242221248','1242049824','1241705168','1236791704'] #actually 2019
   elif image_number==9:
      #bad: 1121321488?
      obsid_list_2015 = ['1121323288','1121325984','1123754200','1123668136','1123581872']
      obsid_list_2018 = ['1244815176','1244384352','1243867368','1243349480','1243177152','1242747232']
   elif image_number==10:
      #leaving out '1234301040' high std dev
      obsid_list_2015 = ['1123495280','1123409544','1123322416','1123237216','1123235984','1122976704']
      obsid_list_2018 = ['1242487840','1242229352','1242057920','1241713264','1236799800']
   elif image_number==11: 
      obsid_list_2015 = ['1086347312','1086355112','1086433480','1086441280','1086519640','1086527440']
      obsid_list_2018 = []
   elif image_number==12:
      #not used
      obsid_list_2015 = ['1086605808','1086613608','1086691968','1086699768','1086778136','1086785936']
      obsid_list_2018 = []

      
   

   ms_dir_list=["/fred/oz048/bmckinle/ATeam/CenA/image%s/2015" % int(image_number),"/fred/oz048/bmckinle/ATeam/CenA/image%s/2018" % int(image_number)]
 
   #return(obsid_list_2015,obsid_list_2018,lst_list_2015,lst_list_2018,ms_dir_list)
   return(obsid_list_2015,obsid_list_2018,ms_dir_list)

#image11 2015:
#1086433480 - some dodginess around ant 120
#same with 1086441280. flag 'em     ants 112-119 incl     
#image12: 1086605808 112-119 again, same 1086613608, same 1086699768,  1086778136
 

image_number = 9
n_selfcals = 5
generate_sbatch_script_CenA(image_number=image_number,n_selfcals=n_selfcals,download=True,model_cal=False)
generate_sbatch_script_CenA(image_number=image_number,n_selfcals=n_selfcals,download=False,model_cal=True)
generate_sbatch_script_CenA(image_number=image_number,n_selfcals=n_selfcals)

#generate_download(obsid_list=obsid_list_2015,dest_dir='2015',timeres=8,freqres=80,ms=True)
#generate_download(obsid_list=obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)   
#   
#generate_unzip(obsid_list=obsid_list_2015,dest_dir='2015')
#generate_unzip(obsid_list=obsid_list_2018,dest_dir='2018')  
#   
#generate_model_cal(obsid_list_2015,model_wsclean_txt='/fred/oz048/bmckinle/code/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[0])
#generate_model_cal(obsid_list_2018,model_wsclean_txt='/fred/oz048/bmckinle/code/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[1])
#
#wsclean_options_1 = " -size 4096 4096 -j 8 -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam -idg-mode hybrid -pb-undersampling 4 -channels-out 8 -join-channels -fit-spectral-pol 2"
#generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base='test1',wsclean_options=wsclean_options_1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4',self_cal_number=0)
#
#calibrate_options_1 = "-minuv 60"
#generate_selfcal(obsid_list,ms_dir_list,calibrate_options=calibrate_options_1,self_cal_number=1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4')

