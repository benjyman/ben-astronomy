#!/usr/bin/env python
#generate scripts to lauch on ozstar for data processing

import os,sys
from astropy.io import fits
import numpy as np

def initiate_script(filename,time_hrs):
   header_text = """#!/bin/bash -l
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=100000 
#SBATCH --account=oz048 
#SBATCH --time=%02d:00:00 
#SBATCH --gres=gpu:1 
#SBATCH --partition=skylake-gpu 

module load gcc/6.4.0 openmpi/3.0.0 
module load fftw/3.3.7 
module load gsl/2.4 
module load cfitsio/3.420 
module load boost/1.67.0-python-2.7.14 
module load hdf5/1.10.1 
module load openblas/0.2.20 
module load cuda/9.0.176 


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/MWA/CODE/lib/

"""        % (time_hrs)
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
   
def generate_unzip(obsid_list,dest_dir):
   out_filename = '%s/unzip_obs.sh' % dest_dir
   initiate_script(out_filename,time_hrs=2) 
   
   cmd = "cd %s;\nunzip -o '*.zip'\nrm *.zip " % dest_dir
   with open(out_filename,'a') as f:
      f.write(cmd)
   
   print("wrote %s" % (out_filename))
   
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
   
   
def generate_wsclean_image():
   print('generating wsclean script')
   
def generate_selfcal():
   print('generating selfcal script')
   
def generate_final_sbatch_script():
   print('generating final sbatch script')
   






obsid_list_2015 = ['1112806040','1112892200','1114782984','1114869144','1114955312','1115041472']
obsid_list_2018 = ['1202239904','1202326064','1202410528','1202411608','1202418952','1202672864']
obsid_list = obsid_list_2015 + obsid_list_2018

ms_dir_list=["/fred/oz048/bmckinle/ATeam/CenA/image4/2015","/fred/oz048/bmckinle/ATeam/CenA/image4/2018"]
   
generate_download(obsid_list=obsid_list_2015,dest_dir='2015',timeres=4,freqres=40,ms=True)
generate_download(obsid_list=obsid_list_2018,dest_dir='2018',timeres=4,freqres=40,ms=True)   
   
generate_unzip(obsid_list=obsid_list_2015,dest_dir='2015')
generate_unzip(obsid_list=obsid_list_2018,dest_dir='2018')  
   
generate_model_cal(obsid_list_2015,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[0])
generate_model_cal(obsid_list_2018,model_wsclean_txt='/fred/oz048/bmckinle/code/git/ben-astronomy/ATeam/CenA/models/CenA_core_wsclean_model.txt',dest_dir=ms_dir_list[1])





