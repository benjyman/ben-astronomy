#!/usr/bin/env python
import os
import subprocess
from datetime import datetime, timedelta, date, time
import string
import fileinput

#takes L-files and makes uvfits files for AAVS-1
#l-files are made such that the integration time is 1.0 s

#make a copy of the header file and replace
#N_SCANS   2     # number of scans (time instants) in correlation products
#FREQCENT  159.375 # observing center freq in MHz
#HA_HRS    -1.22787  # the HA at the *start* of the scan. (hours)
#RA_HRS    1.4888621   # the RA of the desired phase centre (hours)
#DEC_DEGS  9.3533632 # the DEC of the desired phase centre (degs)
#DATE      20180414  # YYYYMMDD
#TIME      030000    # HHMMSS

#base_dir = "/data/data_1/mark_jeff/26_04_2018/2018-04-26 00:00:02/"
#base_dir_fix = "/data/data_1/mark_jeff/26_04_2018/2018-04-26\ 00:00:02/"
base_dir = "/data/data_1/mark_jeff/26_04_2018/2018-04-26 20:00:00/"
base_dir_fix = "/data/data_1/mark_jeff/26_04_2018/2018-04-26\ 20:00:00/"


#for each chan
chan_list = range(15,415,20)
chan_list = [15]
lfile_basename = "test0_16-32chan"
uvfits_basename = "sun"
antenna_locations_filename = "~/rwayth/antenna_locations.txt"
instr_config_filename = "~/rwayth/instr_config.txt"

for chan in chan_list:
   #Get the time of the observation
   channel_dir = "%schannel_%s/" % (base_dir,chan)
   channel_dir_fix = "%schannel_%s/" % (base_dir_fix,chan)
   LACSPC_filename = "%s%s.LACSPC" % (channel_dir_fix,lfile_basename)
   LCCSPC_filename = "%s%s.LCCSPC" % (channel_dir_fix,lfile_basename)


   hd5_filename = os.listdir(channel_dir)[0]
   print hd5_filename

   #use h5ls to get the timestamp
   command = "h5ls -v %s%s" % (channel_dir_fix,hd5_filename)
   output = subprocess.check_output(command, shell=True) 
   #print output
   #output_list = output.split('\n')
   #timestamp = float(output_list[46].split('Data:')[1].strip())

   start_time = float(output.split('ts_start scalar')[1].split('Data:')[1].split('\n')[0].strip())
   #print "start time from hd5 file is %s" % start_time
   #end_time = float(output.split('ts_end scalar')[1].split('Data:')[1].split('\n')[0].strip())
   #print end_time

   #print timestamp_float
   #convert timestamp todate and a time in  YYYYMMDD and  HHMMSS UTC
   start_date_unix_utc=datetime.utcfromtimestamp(start_time)
   print "start time from hd5 file is UTC %s " % start_date_unix_utc
   #end_date_unix_utc=datetime.utcfromtimestamp(start_time)
   #print end_date_unix_utc

   #find the middle of the observation
   #Get the observation length
   #obslength=float((header)['EXPOSURE'])
   #print "EXPOSURE is %s s" % obslength
   obslength = 2.0*60*60 
   half_obslength=obslength/2.
   half_obslength_timedelta=timedelta(seconds=half_obslength)
   middle_of_observation=start_date_unix_utc+half_obslength_timedelta
   #print "middle_of_observation is %s" % middle_of_observation

   middle_of_observation_formatted=middle_of_observation.strftime('%Y-%m-%dT%H:%M:%S')
   print "middle_of_observation_formatted %s " % middle_of_observation_formatted

   #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
   new_date=string.replace(string.replace(middle_of_observation_formatted,'-','/'),'T',' ')
   #print new_date
   #find position of sun
   #src_file_name='%ssrc_file_moon_%s.txt' % (on_moon_directory,on_moon_obsid)
   cmd="~/rwayth/print_src.py --date='%s'" % (new_date)
   output = subprocess.check_output(cmd, shell=True) 
   sun_ra_dec=output.split('Sun')[1].split('\n')[0]
   sun_ra = sun_ra_dec.split('ra = ')[1].split(',')[0].strip()
   sun_ra_hours = float(sun_ra.split(':')[0])
   sun_ra_mins = float(sun_ra.split(':')[1])
   sun_ra_secs = float(sun_ra.split(':')[2])
   sun_ra_timedelta = timedelta(hours=sun_ra_hours, minutes=sun_ra_mins, seconds=sun_ra_secs)
   sun_ra_hours = sun_ra_timedelta.total_seconds() / 60 / 60
   sun_ra_hours_string = "%.7f" % sun_ra_hours
 
   sun_dec = sun_ra_dec.split('dec = ')[1].split(',')[0].strip()
   #print sun_dec
   sun_dec_degrees_only = float(sun_dec.split(':')[0])
   sun_dec_mins_only = float(sun_dec.split(':')[1])
   sun_dec_seconds_only = float(sun_dec.split(':')[2])
   sun_dec_decimal_mins = (sun_dec_seconds_only/60) + sun_dec_mins_only
   sun_dec_decimal_degrees = (sun_dec_decimal_mins/60) + sun_dec_degrees_only
   sun_dec_decimal_degrees_string = "%.7f" % sun_dec_decimal_degrees
   #print sun_dec_decimal_degrees_string

   #work out the sun HA
   #for HA in header need the START time of the scan, not middle
   #get LST
   start_of_observation_formatted=start_date_unix_utc.strftime('%Y-%m-%dT%H:%M:%S')
   print "start_of_observation_formatted %s " % start_of_observation_formatted

   #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
   new_date_start=string.replace(string.replace(start_of_observation_formatted,'-','/'),'T',' ')
   cmd="~/rwayth/print_src.py --date='%s'" % (new_date_start)
   new_output = subprocess.check_output(cmd, shell=True)

   lst = float(new_output.split('LST')[1].split('=')[1].split(')')[0].strip())
   sun_ha_hrs = lst - sun_ra_hours
   #print sun_ha_hrs

   #copy the header template and replace RA and DEC
   new_header_file = channel_dir_fix + "header.txt"
   cmd = "cp ~/rwayth/header_common.txt %s" % (new_header_file)
   #print cmd 
   os.system(cmd)
   new_header_file = channel_dir + "header.txt"
   new_ra_line = "RA_HRS    %s   # the RA of the desired phase centre (hours)" % sun_ra_hours_string
   new_dec_line = "DEC_DEGS  %s  # the DEC of the desired phase centre (degs)" % sun_dec_decimal_degrees_string 
   #get time in right format for header
   new_date_line = "DATE     %s%s%s  # YYYYMMDD" % (new_date.split('/')[0],new_date.split('/')[1],new_date.split('/')[2].split(' ')[0])
   #print new_date_line
   new_time_line = "TIME     %s%s%s  # HHMMSS" % (new_date.split(' ')[1].split(':')[0],new_date.split(':')[1],new_date.split(':')[2].split('\n')[0])
   #print new_time_line
   new_ha_line = "HA_HRS    %.7f  # the HA at the *start* of the scan. (hours)" % sun_ha_hrs
   for line in fileinput.input(new_header_file, inplace=True):
      if "RA_HRS" in line:
         print new_ra_line
      elif "DEC_DEGS" in line:
         print new_dec_line
      elif "DATE" in line:
         print new_date_line
      elif "TIME" in line:
         print new_time_line
      elif "HA_HRS" in line:
         print new_ha_line
      else:
         print line.strip('\n')

   
   #Get the frequency of the observation
   #N_SCANS?

   uvfits_filename = "%s%s_%s%s%s_%s.uvfits" % (channel_dir_fix,uvfits_basename,new_date.split('/')[0],new_date.split('/')[1],new_date.split('/')[2].split(' ')[0],chan)   
   print uvfits_filename
   #convert to uvfits
   cmd = "corr2uvfits -a %s -c %s -o %s -S %s -I %s -H %s -f 0" % (LACSPC_filename,LCCSPC_filename,uvfits_filename,antenna_locations_filename,instr_config_filename,new_header_file)
   print cmd
   os.system(cmd)
