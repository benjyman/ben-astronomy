#!/usr/bin/env python
#Python script - idea to image for SKA Summer School Shanghai
import os
#1. Idea! (EoR!)
#2. Find some observations.
#- website (http://mwa-metadata01.pawsey.org.au MWA-guest guest)
#- find_observations.py

def idea_to_image(options):

   data_dir = "/md0/summer_school/"
   n_obs_to_image = 3 
   obs_filename = "eor2_high_season1_obs.txt"
   cmd = 'find_observations.py --proj=G0009 --racenter=10:20:00 --deccenter=-10:00:00 --limit=5 --obsname=high_season1* > %s' % obs_filename
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
   for obs in obs_list:
      print obs
      csv_line_string = "obs_id=%s, %s" % (obs,manta_ray_options_string)
      csv_file_list.append(csv_line_string)

   with open(csv_filename,'w') as f:
      f.write("\n".join(csv_file_list))

   cmd = "mwa_client -c %s -d %s" % (csv_filename,data_dir)
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
