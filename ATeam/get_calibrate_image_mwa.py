#!/usr/bin/env python
#does what it says

import os,sys

#new ones (close to old from tfranzne): 2015: 1122198832 1122112400 1122025976 1121939544 1121853112 1121766680
#new ones  (close to old from msok) 2018: 1201123376 1201296272 1201382720 1201469160  1201555608 1201814952


obsid_list_2015 = ['1122198832','1122112400','1122025976','1121939544','1121853112','1121766680']
obsid_list_2018 = ['1201123376','1201296272','1201382720','1201469160','1201555608','1201814952']

def download_obs(obsid_list,timeres=4,freqres=40,ms=True):
   #write the csv file
   csv_filename = "manta_ray_download.csv"
   for obsid in obsid_list:
      if ms==True:
         csv_line = "obs_id=%s, job_type=c, timeres=%s, freqres=%s, edgewidth=80, conversion=ms, calibrate=False, allowmissing=False, flagdcchannels=true" % (obsid,timeres,freqres)
         with open(csv_filename,'w') as f:
            f.write(csv_line)
   
   #run manta ray command
   cmd = "mwa_client -c %s -d . " % (csv_filename)
   print(cmd)
   #os.system(cmd)