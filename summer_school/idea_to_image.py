#!/usr/bin/env python
#Python script - idea to image for SKA Summer School Shanghai
import os
#1. Idea! (EoR!)
#2. Find some observations.
#- website (http://mwa-metadata01.pawsey.org.au MWA-guest guest)
#- find_observations.py

os.environ["ASVO_USER"] = "bmckinley"
os.environ["ASVO_PASS"] = "Bridgey2014"

def idea_to_image(options):

   #data_dir = "/md0/summer_school/"
   data_dir="/data/"

   #first lets download a calibrator observation - bright point source, cause this is easy and one way to calibrate the data
   obs_filename="calibrator.txt"
   cmd = 'find_observations.py --proj=G0009 --start=1075635272 --stop=1075647960 --obsname=high_PictorA* > %s' % obs_filename
   print cmd
   os.system(cmd)
   cal_obs_list=[]
   with open(obs_filename,'r') as f:
      lines = f.readlines()
   for line in lines[1:2]:
      obs=line.strip().split()[0]
      cal_obs_list.append(obs)  


   #now for the target data
   n_obs_to_image = 3 
   obs_filename = "eor2_high_season1_obs.txt"
   cmd = 'find_observations.py --proj=G0009 --racenter=10:20:00 --deccenter=-10:00:00 --start=1075652472 --stop=1075653448 --obsname=high_season1* > %s' % obs_filename
   print cmd
   #os.system(cmd)
   #read the obs file
   with open(obs_filename,'r') as f:
      lines = f.readlines()
   obs_list=[]
   for line in lines[1:n_obs_to_image+1]:
      obs=line.strip().split()[0]
      obs_list.append(obs)

   #make csv file
   csv_filename="manta_ray_csv_download.csv"
   manta_ray_options_string = "job_type=c, timeres=10, freqres=80, edgewidth=80, conversion=ms, allowmissing=true, flagdcchannels=true"
   csv_file_list=[]
   obs_list.append(cal_obs_list[0])
   for obs in obs_list:
      print obs
      csv_line_string = "obs_id=%s, %s" % (obs,manta_ray_options_string)
      csv_file_list.append(csv_line_string)

   #print csv_file_list
   with open(csv_filename,'w') as f:
      f.write("\n".join(csv_file_list))

   cmd = "mwa_client -c %s -d %s" % (csv_filename,data_dir)
   print cmd
   os.system(cmd)

   #this gives you zip files
   for obs in obs_list:
      cmd = "unzip %s_ms.zip" % obs
      print cmd
      os.system(cmd)
      #cmd = "rm %s_ms.zip" % obs
      #print cmd
      #os.system(cmd)
      #also need the metafits files, get the most up-to-date:
      cmd = "wget -O %s_metafits_ppds.fits http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s" % (obs,obs)
      print cmd
      os.system(cmd)

   #now you have your measurment sets! what has been done toi the data? see slides
  
   #These data are uncalibrated! 
   #To calibrate you do it two ways: 1. use calibrator with a simple sky model (single source), or 2. use a full model of the sky for your pointing
   
   #Method 2: full sky model:
   for obs in obs_list[1:]:
      print obs
      cmd="/srclists/srclist_by_beam.py -x --aocalibrate -m %s_metafits_ppds.fits -n 1000 -s /srclists/srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt" % (obs)
      print cmd
      #os.system(cmd)   
      #once you have a sourclist, use it to calibrate
      cmd="calibrate -minuv 60  -m srclist_pumav3_EoR0aegean_EoR1pietro+ForA_%s_aocal1000.txt -applybeam  %s.ms  %s_selfcal.bin" % (obs,obs,obs)
      print cmd
      os.system(cmd)


import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: model_moon.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--sky_model',type='string', dest='sky_model',default='gsm',help='Sky model to use for predicted disk-averaged sky (gsm,gsm2016 or lfsm) e.g. --sky_model="gsm2016" [default=%default]')
parser.add_option('--pre_cropped_images',action='store_true',dest='pre_cropped_images',default=False,help='set if using images pre-cropped in imaging step (i.e. with --crop_images e.g. --pre_cropped_images')


(options, args) = parser.parse_args()

idea_to_image(options)
