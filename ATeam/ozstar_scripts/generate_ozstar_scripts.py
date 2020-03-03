#!/usr/bin/env python
#generate scripts to lauch on ozstar for data processing

import os,sys
from astropy.io import fits
import numpy as np

def initiate_script(filename,time_hrs,partition_string='skylake'):
   header_text = """#!/bin/bash -l
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=100000 
#SBATCH --account=oz048 
#SBATCH --time=%02d:00:00 
#SBATCH --gres=gpu:1 
#SBATCH --partition=%s 

module load gcc/6.4.0 openmpi/3.0.0 
module load fftw/3.3.7 
module load gsl/2.4 
module load cfitsio/3.420 
module load boost/1.67.0-python-2.7.14 
module load hdf5/1.10.1 
module load openblas/0.2.20 
module load cuda/9.0.176 


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/MWA/CODE/lib/

"""        % (time_hrs,partition_string)
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
   print('generating model calibration script')
   out_filename = '%s/model_cal.sh' % dest_dir
   initiate_script(out_filename,time_hrs=4) 
   
   cmd = "cd %s\n" % dest_dir
   with open(out_filename,'a') as f:
      f.write(cmd)
   
   for obsid in obsid_list:
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
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/calibrate -minuv 60  -applybeam -m %s %s %s \n" % (model_wsclean_txt,ms_name,cal_sol_filename)
      cmd_list.append(cmd)
      
      cmd = "/fred/oz048/MWA/CODE/mwa-reduce/build/applysolutions %s %s \n" % (ms_name,cal_sol_filename)
      cmd_list.append(cmd)
               
      #image to see if it worked
      cmd = "module load boost/1.67.0-python-2.7.14 \n"
      cmd_list.append(cmd)
      cmd = "wsclean -name %s -size %s %s -niter 0 -data-column CORRECTED_DATA -scale %s -small-inversion -pol xx %s \n" % (check_cal_image_name,int(check_imsize),int(check_imsize),check_scale_deg,ms_name)
      cmd_list.append(cmd)
      
      #do aocal_plot.py
      #Plot the cal solutions
      #cmd = "aocal_plot.py  %s \n" % (cal_sol_filename)
      #cmd_list.append(cmd)

      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]

   print("wrote %s" % (out_filename))
   return(out_filename)
   
def generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base,wsclean_options,dest_dir,self_cal_number=0):
   print('generating wsclean script')
   out_filename = '%s/wsclean_selfcal_%02d.sh' % (dest_dir,self_cal_number)
   initiate_script(out_filename,time_hrs=15,partition_string='skylake-gpu')
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
   print('generating selfcal script')
   out_filename = '%s/selfcal_%02d.sh' % (dest_dir,self_cal_number)
   initiate_script(out_filename,time_hrs=2,partition_string='skylake')
   cmd_list = []
   cmd = "module load boost/1.66.0-python-2.7.14 \n"
   cmd_list.append(cmd)
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
      cmd = "calibrate %s %s %s \n" % (calibrate_options,ms,solution_name)
      cmd_list.append(cmd)
      
      cmd = "applysolutions %s %s \n" % (ms,solution_name)
      cmd_list.append(cmd)
      
      #Plot the cal solutions ##need to get this working
      #cmd = "aocal_plot.py  %s \n" % (solution_name)
      #cmd_list.append(cmd)
      
   with open(out_filename,'a') as f:
      [f.write(cmd) for cmd in cmd_list]
   print("wrote %s" % (out_filename))
   return(out_filename)
   
def generate_sbatch_script_CenA(obsid_list,ms_dir_list,n_selfcals,download=False,model_cal=False):
   print('generating sbatch script for CenA')

   if download:
      out_filename = 'sbatch_launch_download.sh'
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
   
      out_filename_2015 = generate_download(obsid_list=obsid_list_2015,dest_dir='2015',timeres=4,freqres=40,ms=True)
      out_filename_2018 = generate_download(obsid_list=obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)   
      
      # first commands no dependencies
      cmd = 'jid1=$(sbatch %s) \n' % out_filename_2015
      cmd_list.append(cmd)
      cmd = 'jid2=$(sbatch %s) \n' % out_filename_2018
      cmd_list.append(cmd)
      
      #unzip
      out_filename_2015 = generate_unzip(obsid_list=obsid_list_2015,dest_dir='2015')
      out_filename_2018 = generate_unzip(obsid_list=obsid_list_2018,dest_dir='2018')  
      
      #unzip commands depend on downloads
      cmd = 'jid3=$(sbatch --dependency=afterok:$jid1 %s) \n' % out_filename_2015
      cmd_list.append(cmd)
      cmd = 'jid4=$(sbatch --dependency=afterok:$jid2 %s) \n' % out_filename_2018
      cmd_list.append(cmd)
   
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      return(out_filename)
      
   #do qa manually after download and unzip then proceed
   elif model_cal:
      out_filename = 'sbatch_launch_model_cal.sh'
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
   
      out_filename_2015 = generate_model_cal(obsid_list_2015,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[0])
      out_filename_2018 = generate_model_cal(obsid_list_2018,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[1])
      
      cmd = 'jid1=$(sbatch %s) \n' % out_filename_2015
      cmd_list.append(cmd)
      cmd = 'jid2=$(sbatch %s) \n' % out_filename_2018
      cmd_list.append(cmd)
      
      
      
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      return(out_filename)

   #do qa again after model cal then proceed
   else:
      out_filename = 'sbatch_launch_wsclean_and_selfcal.sh'
      cmd = '#! /bin/bash \n'
      with open(out_filename,'w') as f:
         f.write(cmd)
      cmd_list = []
      
      for selfcal in range(n_selfcals*2,2):
         print selfcal
         wsclean_options_1 = " -size 4096 4096 -j 8 -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam -idg-mode hybrid -pb-undersampling 4 -channels-out 8 -join-channels -fit-spectral-pol 2"
         calibrate_options_1 = "-minuv 60"
         
         out_filename = generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base='test1',wsclean_options=wsclean_options_1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4',self_cal_number=selfcal)
   
         cmd = 'jid%s=$(sbatch %s) \n' % (selfcal,out_filename)
         cmd_list.append(cmd)
         
         out_filename = generate_selfcal(obsid_list,ms_dir_list,calibrate_options=calibrate_options_1,self_cal_number=1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4')

         cmd = 'jid%s=$(sbatch %s) \n' % (selfcal+1,out_filename)
         cmd_list.append(cmd)
      
      
   
      with open(out_filename,'a') as f:
         [f.write(cmd) for cmd in cmd_list]
      print("wrote %s" % (out_filename))
      return(out_filename)
   


obsid_list_2015 = ['1112806040','1112892200','1114782984','1114869144','1114955312','1115041472']
obsid_list_2018 = ['1202239904','1202326064','1202410528','1202411608','1202418952','1202672864']
obsid_list = obsid_list_2015 + obsid_list_2018

ms_dir_list=["/fred/oz048/bmckinle/ATeam/CenA/image4/2015","/fred/oz048/bmckinle/ATeam/CenA/image4/2018"]

generate_sbatch_script_CenA(obsid_list,ms_dir_list,n_selfcals=4,download=True,model_cal=False)
generate_sbatch_script_CenA(obsid_list,ms_dir_list,n_selfcals=4,download=False,model_cal=True)
generate_sbatch_script_CenA(obsid_list,ms_dir_list,n_selfcals=4)

#generate_download(obsid_list=obsid_list_2015,dest_dir='2015',timeres=8,freqres=80,ms=True)
#generate_download(obsid_list=obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)   
#   
#generate_unzip(obsid_list=obsid_list_2015,dest_dir='2015')
#generate_unzip(obsid_list=obsid_list_2018,dest_dir='2018')  
#   
#generate_model_cal(obsid_list_2015,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[0])
#generate_model_cal(obsid_list_2018,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[1])
#
#wsclean_options_1 = " -size 4096 4096 -j 8 -mwa-path /fred/oz048/MWA/CODE/MWA_Tools/mwapy/data -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam -idg-mode hybrid -pb-undersampling 4 -channels-out 8 -join-channels -fit-spectral-pol 2"
#generate_wsclean_image(obsid_list,ms_dir_list,out_image_name_base='test1',wsclean_options=wsclean_options_1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4',self_cal_number=0)
#
#calibrate_options_1 = "-minuv 60"
#generate_selfcal(obsid_list,ms_dir_list,calibrate_options=calibrate_options_1,self_cal_number=1,dest_dir='/fred/oz048/bmckinle/ATeam/CenA/image4')

