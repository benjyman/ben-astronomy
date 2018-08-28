#!/usr/bin/env python 
#python script to cotter moon obs

from astropy.io import fits
import string
import os.path
import cmd
import datetime
import numpy as np
import subprocess


def manta_ray(options):
   machine=options.machine
   if machine=='namorrodor':
      mwa_dir = '/md0/moon/data/MWA/'
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'
       
   obsid=options.obsid
   if options.sister_obsid:
      sister_obsid=options.sister_obsid
   
   time_res=(options.time_freq_res).split(',')[0]
   freq_res=(options.time_freq_res).split(',')[1]
   
   
   if (options.track_off_moon):
      track_off_moon_string= options.track_off_moon_list
   else:
      track_off_moon_string=' '

   
   if options.sister_obsid:   
      print "downloading ms obsid %s paired with %s with time resolution %s and freq resolution %s" % (obsid,sister_obsid,time_res,freq_res)
   else:
      print "downloading ms obsid %s with time resolution %s and freq resolution %s" % (obsid,time_res,freq_res)

   if options.track_off_moon:
      track_off_moon_list=track_off_moon_string.split(',')

      track_off_moon_paired_obsid=track_off_moon_list[0].strip()

      track_off_moon_new_RA=track_off_moon_list[1].strip()

      track_off_moon_new_DEC=track_off_moon_list[2].strip()

      print "obs paired with %s centering at RA:%s DEC:%s" % (track_off_moon_paired_obsid,track_off_moon_new_RA,track_off_moon_new_DEC)
 

   data_dir='%s%s/' % (mwa_dir,obsid)
   
   #tagname=options.tagname
   epoch_ID=options.epoch_ID
   #base_name=obsid+'_cotter_'+tagname
   base_name= "%s_%s" % (obsid,epoch_ID)
   
   track_moon_string=' '
   
   #if (options.track_moon):
   #   base_name+='_trackmoon'
   #
   #if (options.track_off_moon):
   #   base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
   #
   #ms_name=data_dir+base_name+'.ms'
   
   if (options.track_moon):
      on_moon_basename=base_name + '_trackmoon'
      on_moon_ms_name=data_dir+on_moon_basename+'.ms'
      off_moon_base_name="%s_%s_track_off_moon_paired_%s" % (sister_obsid,epoch_ID,obsid)
      off_moon_ms_name=mwa_dir+sister_obsid+'/'+off_moon_base_name+'.ms'
      ms_name=on_moon_ms_name
   elif (options.track_off_moon):
      on_moon_basename="%s_%s_trackmoon" % (sister_obsid,epoch_ID) 
      on_moon_ms_name=mwa_dir+sister_obsid+'/'+on_moon_basename+'.ms'
      off_moon_base_name=base_name+'_track_off_moon_paired_' + sister_obsid
      off_moon_ms_name=data_dir+off_moon_base_name+'.ms'
      ms_name=off_moon_ms_name
   else:
      if (machine == "magnus" or machine == "galaxy"):
         ms_name=base_name+'.ms'
      else:
         ms_name=data_dir+base_name+'.ms'

   #Check if the ms already exists and is the right size then just exit - no need to make it again
   if options.use_existing_ms:
      if epoch_ID[0:5]=="2015B":
         #min size in Kb
         min_ms_size=800000.0
      else:
         min_ms_size=100000.0
   
      ms_size = float(subprocess.check_output(['du','-s', ms_name]).split()[0].decode('utf-8').split("'")[0])
      if (os.path.exists(ms_name) and ms_size >= min_ms_size):
         print "%s already exists and has size %s (greater than min size %s), exiting cotter_moon.py" % (ms_name, ms_size, min_ms_size)
         sys.exit(0)
         #exit  

   #metafits_filename=data_dir+obsid+'.metafits'
   metafits_filename="%s%s_metafits_ppds.fits" % (data_dir,obsid)
   
   
   #flag_file_name=data_dir+obsid+'_flags.zip'

   #check if metafits file exists, if not, download it
   #No always download new fits file
   #if not os.path.isfile(metafits_filename):
   #   #cmd="make_metafits.py -o %s -g %s" % (metafits_filename,obsid)
   #cmd="wget -O %s http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s" % (metafits_filename,obsid)
   #print cmd
   #os.system(cmd)

   ##unpack the flagfile
   #cmd="unzip %s -d %s " % (flag_file_name,data_dir)
   #print cmd
   #os.system(cmd)

   #if tracking moon, find moon position
   if (options.track_moon):
      #get the date and time of the observation from the metafits file
      try:
         HDU_list = fits.open(metafits_filename)
      except IOError, err:
         'Cannot open metadata file %s\n' % str(options.input_file)
      header=HDU_list[0].header 
      date_unix_utc_string = (header)['GOODTIME']
      #datetime.datetime.utcfromtimestamp(posix_time).strftime('%Y-%m-%dT%H:%M:%SZ')
      date_unix_utc=datetime.datetime.utcfromtimestamp(date_unix_utc_string)
      print "start time (GOODTIME) of observation is %s" % (date_unix_utc)
      
      #Get the observation length
      obslength=float((header)['EXPOSURE'])
      print "EXPOSURE is %s s" % obslength
      half_obslength=obslength/2.
      half_obslength_timedelta=datetime.timedelta(seconds=half_obslength)
      middle_of_observation=date_unix_utc+half_obslength_timedelta
      print "middle_of_observation is %s" % middle_of_observation
      
      middle_of_observation_formatted=middle_of_observation.strftime('%Y-%m-%dT%H:%M:%S')
      print "middle_of_observation_formatted %s " % middle_of_observation_formatted
      
      #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
      new_date=string.replace(string.replace(middle_of_observation_formatted,'-','/'),'T',' ')
      print new_date
      #find position of moon
      src_file_name='src_file_moon_%s.txt' % obsid 
      cmd="print_src.py --date='"+new_date+"' > "+ src_file_name
      print cmd
      os.system(cmd)
      #read src_file.txt to get Moon ra and dec 
      with open(src_file_name, "r") as infile:
         lines=infile.readlines()
         moon_ra=lines[6].split()[3].split(',')[0]
         moon_dec=lines[6].split()[6].split(',')[0]
         #print moon_ra
         #print moon_dec
      #get ra and dec in right format for cotter
      new_moon_dec=string.replace(moon_dec,":",".")
      track_moon_string=' %s %s ' % (moon_ra,new_moon_dec)
      print track_moon_string


   #if tracking off moon, shift to required moon  position from a previous time
   if (options.track_off_moon):
      off_moon_ra=track_off_moon_new_RA
      off_moon_dec=track_off_moon_new_DEC
      track_moon_string=' %s %s ' % (off_moon_ra,off_moon_dec)
      print track_moon_string
   
   ##flagfiles for this observation
   #flagfiles_string= " %s%s_%s.mwaf " % (data_dir,obsid,"%%")
   #print flagfiles_string

   ##flagfiles for paired sister observation
   #flagfiles_string_sister= " %s%s/%s_%s.mwaf " % (mwa_dir,sister_obsid,sister_obsid,"%%")
   #print flagfiles_string_sister
   
   #Don't do this. Instead, make the cotter file then merge the flag columns of the ms s later....
   ##combine the flagfiles - for each flag_file (coarse chan):
   #for coarse_chan in range(1,2):
   #   flag_file_1_name="%s%s/%s_%02d.mwaf" % (mwa_dir,obsid,obsid,coarse_chan)
   #   print flag_file_1_name
   #   flag_file_1_HDU_list=fits.open(flag_file_1_name)
   #   flag_file_1_HDU_list.info()
   #   header1=flag_file_1_HDU_list[0].header
   #   print header1
   #   data1=flag_file_1_HDU_list[0].data
   #   print data1
   #   header2=flag_file_1_HDU_list[1].header
   #   print header2
   #   data2=flag_file_1_HDU_list[1].data
      #print data2
   
   #need to put user options on the time and freq resolution - what is best for the new long baseline obs?
   #can probably get away with just halving each (double baselines)
   
   #make csv file
   csv_filename="manta_ray_csv_download.csv"
   manta_ray_options_string = "obs_id=%s, job_type=c, timeres=%s, freqres=%s, edgewidth=80, conversion=ms, allowmissing=false, flagdcchannels=true, noantennapruning=true" % (obsid,time_res,freq_res)

   with open(csv_filename,'w') as f:
      f.write(manta_ray_options_string)

   cmd = "mwa_client -c %s -d %s" % (csv_filename,data_dir)
   print cmd
   os.system(cmd)

   #this gives you zip files - unzip!
   cmd = "unzip %s%s_ms.zip" % (data_dir,obsid)
   print cmd
   os.system(cmd)
   
   #rename ms:
   cmd = "mv %s.ms %s" % (data_dir,obsid,ms_name)
   print cmd
   os.system(cmd)
   
   #change the phase centre of the ms
   cmd = "chgcentre %s %s" % (ms_name,track_moon_string)
   print cmd
   os.system(cmd)
   
   #cmd='cotter -flagfiles %s -norfi -noantennapruning -o %s %s -m %s -timeres %s -freqres %s %s*gpubox*.fits' % (flagfiles_string,ms_name,track_moon_string,metafits_filename,time_res,freq_res,data_dir)
   #print cmd
   #os.system(cmd)

   if (options.flag_ants):
      flag_ants_cmd_string="flagantennae %s %s " % (ms_name,options.flag_ants)
      print flag_ants_cmd_string
      os.system(flag_ants_cmd_string)
   
   #if (options.cleanup and os.path.exists(ms_name) and ms_size >= min_ms_size):
   #   cmd="rm -rf  %s*gpubox*.fits" % (data_dir)
   #   print cmd
   #   os.system(cmd)

   #do all this stuff in the selfcal stage instead - so you are sure that both ms have already meen created
   #here merge the flag columns for both the on_moon and off_moon ms
   #only do this step for track_off_moon, (always need to make on_moon ms first)
   #if (options.track_off_moon):
   #   on_moon_table=table(on_moon_ms_name,readonly=False)
   #   #print on_moon_table
   #   #on_moon_table_UVW=tablecolumn(on_moon_table,'UVW')
   #   #on_moon_table_time=tablecolumn(on_moon_table,'TIME')
   #   on_moon_table_flag=tablecolumn(on_moon_table,'FLAG')
   #
   #   off_moon_table=table(off_moon_ms_name,readonly=False)
   #   #print off_moon_table
   #   #off_moon_table_UVW=tablecolumn(off_moon_table,'UVW')
   #   #off_moon_table_time=tablecolumn(off_moon_table,'TIME')
   #   off_moon_table_flag=tablecolumn(off_moon_table,'FLAG')
   #
   #   new_off_moon_table_flag=np.logical_not(off_moon_table_flag)
   #
   # 
   #   #look in obs_list_generator.py for how to compare LSTs
   #
   #   #Not sure if 2018A obs are well enough LST matched
   #   #Continue anyway (come back and look at 2015A) 
   #
   #   #OR the two flag columns
   #   new_flag_column=np.logical_or(on_moon_table_flag,off_moon_table_flag)
   #   #print new_flag_column[20100:20101]
   # 
   #   on_moon_table.putcol('FLAG',new_flag_column)
   #   off_moon_table.putcol('FLAG',new_flag_column)
   #
   #   on_moon_table.close()
   #   off_moon_table.close()

   
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: manta_ray.py [obsid] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
#parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms file e.g. --tagname="" [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='List of antennas (space separated) to flag after cottering (andre indexing as for rfigui etc)  e.g. --flag_ants="56,60" [default=%default]')
parser.add_option('--use_existing_ms',action='store_true',dest='use_existing_ms',default=False,help='If there is already a measurment set dont remake it [default=%default]') 
parser.add_option('--obsid',type='string', dest='obsid',default='',help='obsid to be cottered e.g. --obsid="1199394880" [default=%default]')
parser.add_option('--sister_obsid',type='string', dest='sister_obsid',default='',help='sister_obsid e.g. --sister_obsid="1199396880" [default=%default]')
parser.add_option('--track_off_moon_list',type='string', dest='track_off_moon_list',default='',help='When track_off_moon is True. Details of on-moon pairing on_moon_obsid,RA,DEC of on_moon paired obs e.g. --track_off_moon_list="1199394880,13.0,0.56" [default=%default]')
parser.add_option('--time_freq_res',type='string', dest='time_freq_res',default='8,80',help='Time and then frequency resolution, comma separated e.g. --time_freq_res="8,80" [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')


(options, args) = parser.parse_args()

#obsid = args[0]
#if (options.track_off_moon):
#   track_off_moon_string= args[1]
#else:
#   track_off_moon_string=' '

manta_ray(options)


