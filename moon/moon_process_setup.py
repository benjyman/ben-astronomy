#script to set up Moon processing in bulk (on magnus)!
#just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date
import sys,os
import cmd
from datetime import datetime, timedelta, date, time
from astropy.io import fits
import string
import glob
from casacore.tables import table,tablecolumn,tablerow
import numpy as np
import subprocess

sidereal_day_sec = 86164
sourcelist='srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt'
havebeam_string = ''
namorrodor_image_dir = "/md0/moon/images/"

#get filepaths using gator_read_table_only_peeled.py (need to fix this script for peeled images)
#or just rsync -a -v bmckinley@magnus.pawsey.org.au:/astro/mwaeor/MWA/data/11273*/2018-08-*/*cropped*I.fits .

def check_metafits_tile_flags(on_moon_metafits_file_name,off_moon_metafits_file_name):
      #on moon
   try:
      print on_moon_metafits_file_name
      HDU_list = fits.open(on_moon_metafits_file_name)
   except IOError, err:
      'Cannot open metadata file %s\n' % on_moon_metafits_file_name
   data = HDU_list[1].data
   tile_flags = data['FLAG']
   tile_flagging_indices = data['ANTENNA']
   tile_name = data['TILENAME']
   #print tile_flags
   
   flag_antenna_indices_list = []
   flag_antenna_tilenames_list = []
   for index,tile_flag in enumerate(tile_flags):
      if (tile_flag==1):
         print tile_name[index]

def copy_images_from_magnus(filename):
   for filepath in open(filename):
      cmd = "rsync -a bmckinley@magnus.pawsey.org.au:%s*I.fits %s" % (filepath.rstrip(),namorrodor_image_dir)
      print cmd
      os.system(cmd)
      cmd = "rsync -a bmckinley@magnus.pawsey.org.au:%s*psf*.fits %s" % (filepath.rstrip(),namorrodor_image_dir)
      print cmd
      os.system(cmd)
   
#Functions
def write_obs_list(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir):
   on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
   #print on_off_moon_string
   if on_off_moon_string=='on_moon':
      start_time_string="%s 00:00:00" % on_moon_date
      start_time_datetime=datetime.strptime(start_time_string, '%Y-%m-%d %H:%M:%S')
      stop_time_datetime=start_time_datetime + timedelta(hours=24)
      stop_time_string=stop_time_datetime.strftime('%Y-%m-%d %H:%M:%S')
      find_observations_filename="%s%s_on_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   elif on_off_moon_string=='off_moon':
      start_time_string="%s 00:00:00" % off_moon_date
      start_time_datetime=datetime.strptime(start_time_string, '%Y-%m-%d %H:%M:%S')
      stop_time_datetime=start_time_datetime + timedelta(hours=24)
      stop_time_string=stop_time_datetime.strftime('%Y-%m-%d %H:%M:%S')
      find_observations_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   else:
      print "Bad values for on/off moon"
               
   cmd="find_observations.py --proj='G0017' --chan=%s --start=%s --stop=%s --obsname=EoRMoon_%s -q > %s" % (chan,start_time_string,stop_time_string,chan,find_observations_filename)
   print cmd
   os.system(cmd)
   return find_observations_filename
   
def write_paired_obslist(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir):
   on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
   paired_obs_list=[]
   obsid_list=[]
   sister_obsid_list=[]
   if on_off_moon_string=='on_moon':
      obsid_list_filename="%s%s_on_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
      off_moon_dir=os.path.dirname(os.path.dirname(on_off_moon_dir))+'/off_moon/'
      sister_obsid_list_filename="%s%s_off_moon_%s.txt" % (off_moon_dir,epoch_ID,chan) 
      paired_obs_list_filename="%s%s_on_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   elif on_off_moon_string=='off_moon':
      obsid_list_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
      on_moon_dir=os.path.dirname(os.path.dirname(on_off_moon_dir))+'/on_moon/'
      sister_obsid_list_filename="%s%s_on_moon_%s.txt" % (on_moon_dir,epoch_ID,chan)
      paired_obs_list_filename="%s%s_off_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   else:
      print "Bad values for on/off moon"      
   for line in open(obsid_list_filename):
      obsid_list.append(line.strip()) 
   for line in open(sister_obsid_list_filename):
      sister_obsid_list.append(line.strip())    
   for item_index,item in enumerate(obsid_list):
      paired_obs_string="%s%s" % (str(obsid_list[item_index]),str(sister_obsid_list[item_index]))
      paired_obs_list.append(paired_obs_string)
   with open(paired_obs_list_filename,'w') as f:
      f.write('\n'.join(paired_obs_list)) 
   return paired_obs_list_filename
    
         
def write_and_run_default_scripts(epoch_ID,chan,on_off_moon_dir,machine):
   #function to write the default scripts: download, cotter, calibrate, image, pbcorr
   on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
   if on_off_moon_string=='on_moon':
      observations_filename="%s%s_on_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)   
      off_moon_dir=os.path.dirname(os.path.dirname(on_off_moon_dir))+'/off_moon/'
      sister_observations_filename="%s%s_off_moon_%s.txt" % (off_moon_dir,epoch_ID,chan)    
      track_moon_string='--track_moon'
   elif on_off_moon_string=='off_moon':
      observations_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)
      observations_filename_no_path=observations_filename.split(on_off_moon_dir)[1]
      track_off_moon_filename='%strack_off_moon_%s' % (on_off_moon_dir,observations_filename_no_path)
      on_moon_dir=os.path.dirname(os.path.dirname(on_off_moon_dir))+'/on_moon/'
      sister_observations_filename="%s%s_on_moon_%s.txt" % (on_moon_dir,epoch_ID,chan) 
      track_moon_string='--track_off_moon=%s' % track_off_moon_filename
   else:
      print "Bad values for on/off moon of %s" % on_off_moon_string

   database_name='/group/mwaeor/bmckinley/%s_%s_%s.sqlite' % (epoch_ID,chan,on_off_moon_string)


   #things that depend on obs semester
   obs_semester=epoch_ID[0:5]
   if (obs_semester=='2015A' or obs_semester=='2015B'):
       time_averaging='8'
       freq_averaging='80'
       imsize_string='--imsize=2048'
       wsclean_options_string='--wsclean_options=" -niter 0 -datacolumn CORRECTED_DATA  -scale 0.0085 -weight uniform  -smallinversion  -channelsout 24 -make-psf  "'
   elif (obs_semester=='2017B' or obs_semester=='2018A'):
       time_averaging='4'
       freq_averaging='40'  
       imsize_string='--imsize=4096'
       wsclean_options_string='" -niter 0  -datacolumn CORRECTED_DATA  -scale 0.0042 -weight natural  -smallinversion -channelsout 24 -make-psf "'
   else:
       print "observing semester %s not known" % obs_semester
   
   
   if machine=='namorrodor':
      ben_code_base='/data/code/git/'
      srclist_code_base='/data/code/git/srclists/'
   else:
      ben_code_base='/astro/mwaeor/bmckinley/code/'
      srclist_code_base='/group/mwa/software/srclists/master/'
   #1. download
   default_download_script_name="%s1_default_download_%s_%s_%s.sh" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   generate_download_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_obs_download.py --machine=%s --obsid_infile=%s --database=%s' % (ben_code_base,machine,observations_filename,database_name)
   with open(default_download_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_download_string)
   #run script to generate the download script
   cmd = "chmod +x %s" % default_download_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_download_script_name
   print cmd
   os.system(cmd)

   #2. cotter
   default_cotter_script_name="%s2_default_cotter_%s_%s_%s.sh" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   #with cleanup
   if options.cleanup:
      generate_cotter_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_cotter_moon.py --epoch_ID=%s %s --flag_ants="" --cleanup --obsid_infile=%s --sister_obsid_infile=%s' % (ben_code_base,epoch_ID,track_moon_string,observations_filename,sister_observations_filename)
   else:
   #no cleanup
      generate_cotter_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_cotter_moon.py --epoch_ID=%s %s --flag_ants="" --obsid_infile=%s --sister_obsid_infile=%s' % (ben_code_base,epoch_ID,track_moon_string,observations_filename,sister_observations_filename)
   with open(default_cotter_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_cotter_string)
   cmd = "chmod +x %s" % default_cotter_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_cotter_script_name
   print cmd
   os.system(cmd)
      
   #3. selfcal
   sourcelist_filepath='%s%s' % (srclist_code_base,sourcelist)
   default_selfcal_script_name="%s3_default_selfcal_%s_%s_%s.sh" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   generate_selfcal_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_qselfcal_concat_ms.py --epoch_ID=%s --sourcelist=%s --cotter %s  --selfcal=0 --sister_obsid_infile=%s --obsid_infile=%s' % (ben_code_base,epoch_ID,sourcelist_filepath,track_moon_string,sister_observations_filename,observations_filename) 
   with open(default_selfcal_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_selfcal_string)
   cmd = "chmod +x %s" % default_selfcal_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_selfcal_script_name
   print cmd
   os.system(cmd)
   
   #4. image
   default_image_script_name="%s4_default_image_%s_%s_%s.sh" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   generate_image_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_mwac_qimage_concat_ms.py --cotter --no_pbcorr --epoch_ID=%s  %s %s  --pol="xx,xy,yx,yy"  %s  --obsid_infile=%s  ' % (ben_code_base,epoch_ID,track_moon_string,imsize_string, wsclean_options_string,observations_filename) 
   with open(default_image_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_image_string)
   cmd = "chmod +x %s" % default_image_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_image_script_name
   print cmd
   os.system(cmd)
      
   #5. pbcorr
   default_pbcorr_script_name="%s5_default_pbcorr_%s_%s_%s.sh" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   generate_pbcorr_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_qpbcorr_multi.py  --epoch_ID=%s %s  --obsid_infile=%s --dirty  --channelsout=24 %s --machine=%s   ' % (ben_code_base,epoch_ID,track_moon_string,observations_filename, havebeam_string, machine) 
   with open(default_pbcorr_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_pbcorr_string)
   cmd = "chmod +x %s" % default_pbcorr_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_pbcorr_script_name
   print cmd
   os.system(cmd)

   #Don't need this anymore - gator handles the sbatching   
   ##launch the jobs if requested
   #launch_job_filename='%slaunch_jobs.sh'%on_off_moon_dir
   #with open(launch_job_filename,'w') as f:
   #   if options.launch_cotter:
   #      q_cotter_files_list=glob.glob('%s/q_cotter_moon*'%on_off_moon_dir)
   #      number_of_q_cotter_files=len(q_cotter_files_list)
   #      print 'number_of_q_cotter_files is %s ' % number_of_q_cotter_files
   #      if number_of_q_cotter_files==0:
   #         print 'No q_cotter files '
   #      elif number_of_q_cotter_files==1:
   #         cmd1='jobid=`sbatch %sq_cotter_moon_0.sh | cut -d " " -f 4`' % on_off_moon_dir
   #         #print cmd1
   #         #os.system(cmd1)
   #         f.write(cmd1)
   #      else:
   #         for q_cotter_file_index,q_cotter_file in enumerate(q_cotter_files_list):
   #            cmd1='jobid_%s=`sbatch %s | cut -d " " -f 4`' % (str(q_cotter_file_index),q_cotter_file)
   #            f.write(cmd1)
   #      #print 'WARNING! Too many q_cotter files - now you need to write this bit of code!'
   #   if options.launch_selfcal:
   #      if options.launch_cotter:
   #         if number_of_q_cotter_files==0:
   #            print 'No q_cotter files '
   #         elif number_of_q_cotter_files==1:
   #            cmd2='jobid=`sbatch --dependency=afterok:$jobid %sq_selfcal_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
   #            f.write(cmd2)
   #         else:
   #            dependency_string='--dependency=afterok'
   #            for q_cotter_file_index,q_cotter_file in enumerate(q_cotter_files_list):
   #               dependency_string+=':jobid_%s' % (str(q_cotter_file_index))
   #            cmd2='jobid=`sbatch %s %sq_selfcal_moon.sh | cut -d " " -f 4`' % (dependency_string,on_off_moon_dir)
   #            f.write(cmd2)          
   #            #print 'WARNING! Too many q_cotter files - now you need to write this bit of code!'
   #      else:
   #         cmd2='jobid=`sbatch %sq_selfcal_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
   #      #print cmd2
   #      #os.system(cmd2)
   #      
   #   if options.launch_image:
   #      if options.launch_selfcal:
   #         cmd3='jobid=`sbatch --dependency=afterok:$jobid %sq_image_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
   #      else:
   #         cmd3='jobid=`sbatch %sq_image_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
   #      #print cmd3
   #      #os.system(cmd3)
   #      f.write(cmd3)
   #   if options.launch_pbcorr:
   #      if options.launch_image:
   #         cmd4='sbatch --dependency=afterok:$jobid %sq_pbcorr_moon.sh' % on_off_moon_dir
   #      else:
   #         cmd4='sbatch %sq_pbcorr_moon.sh' % on_off_moon_dir
   #      #print cmd4
   #      #os.system(cmd4)
   #      f.write(cmd4)
   #cmd='chmod +x %s' % launch_job_filename
   #print cmd
   #os.system(cmd)
   #cmd='%s' % launch_job_filename
   #print cmd
   #os.system(cmd)

#def launch_gator(epoch_ID,chan,chan_dir):
#   database_name='/group/mwaeor/bmckinley/%s_%s.sqlite' % (epoch_ID,chan)
#   #we only use the on_moon paired observsations filename - gator takes care of the off moon
#   paired_observations_filename="%son_moon/%s_on_moon_paired_%s.txt" % (chan_dir,epoch_ID,chan)
#
#   cmd='gator_add_to_rts_table.rb -d %s %s' % (database_name,paired_observations_filename)
#   print cmd
#   os.system(cmd)
#
#   cmd='nohup gator_rts_daemon.rb -d %s --cotter &' % (database_name)
#   print cmd
#   os.system(cmd)
         
def make_track_off_moon_file(on_moon_obsid_filename,off_moon_obsid_filename,machine,epoch_ID,chan):
   #make the file that has the sister moon obsid moon ra (hh:mm:ss.ss)  moon dec (dd.mm.ss.s) and off_moon obsid
   if machine=='namorrodor':
      mwa_dir = '/md0/moon/data/MWA/'
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'   
   on_moon_directory=os.path.dirname(on_moon_obsid_filename)+'/'
   off_moon_directory=os.path.dirname(off_moon_obsid_filename)+'/'
   off_moon_obsid_filename_no_path=off_moon_obsid_filename.split(off_moon_directory)[1]
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]
   paired_obs_list=[]
   track_off_moon_string_list=[]
   track_off_moon_filename='%strack_off_moon_%s' % (off_moon_directory,off_moon_obsid_filename_no_path)

   paired_obs_list_filename="%s%s_on_moon_paired_%s.txt" % (on_moon_directory,epoch_ID,chan)

   for line in open(on_moon_obsid_filename):
      on_moon_obsid_list.append(line.strip()) 
   for line in open(off_moon_obsid_filename):
      off_moon_obsid_list.append(line.strip()) 
   n_obs = sum(1 for line in open(on_moon_obsid_filename))
   
   for obsid_index,on_moon_obsid in enumerate(on_moon_obsid_list): 
      on_moon_data_dir = "%s%s/" % (mwa_dir,on_moon_obsid)
      off_moon_obsid=off_moon_obsid_list[obsid_index]
      off_moon_data_dir = "%s%s/" % (mwa_dir,off_moon_obsid)
      
      #paired obs list
      paired_obs_string="%s%s" % (on_moon_obsid,off_moon_obsid)
      paired_obs_list.append(paired_obs_string) 

      #Check LSTs
      LST_difference=float(off_moon_obsid)-float(on_moon_obsid)
      LST_remainder = LST_difference % sidereal_day_sec
      LST_difference_in_sidereal_days=abs(LST_remainder/sidereal_day_sec)
      print "LST difference is: %s in sidereal_days" % LST_difference_in_sidereal_days
      
      metafits_filename='%s%s_metafits_ppds.fits' % (on_moon_directory,on_moon_obsid)
      off_moon_metafits_filename='%s%s_metafits_ppds.fits' % (off_moon_directory,off_moon_obsid)
      #download the metafits file
      cmd='wget -O %s http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s' % (metafits_filename,on_moon_obsid)
      print cmd
      os.system(cmd)
      cmd='wget -O %s http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s' % (off_moon_metafits_filename,off_moon_obsid)
      print cmd
      os.system(cmd)
      
      #copy the metafits files to the data dir if galaxy or magnus
      if (machine == 'galaxy' or machine == 'magnus'):
         on_moon_metafits_data_dir_filename = '%s%s_metafits_ppds.fits' % (on_moon_data_dir,on_moon_obsid)
         off_moon_metafits_data_dir_filename = '%s%s_metafits_ppds.fits' % (off_moon_data_dir,off_moon_obsid)
         cmd = "cp %s %s " % (metafits_filename,on_moon_metafits_data_dir_filename)
         print cmd 
         os.system(cmd)
         cmd = "cp %s %s " % (off_moon_metafits_filename,off_moon_metafits_data_dir_filename)
         print cmd 
         os.system(cmd)
                  
      #read the metafits files and see how many tiles are flagged (then decide whether to process this epoch)
      #on moon 
      try:
         print metafits_filename
         HDU_list = fits.open(metafits_filename)
      except IOError, err:
         'Cannot open metadata file %s\n' % metafits_filename
      data = HDU_list[1].data
      tile_flags = data['FLAG']
      tile_flagging_indices = data['ANTENNA']
      tile_name = data['TILENAME']
      #print tile_flags
      
      flag_antenna_indices_list = []
      flag_antenna_tilenames_list = []
      for index,tile_flag in enumerate(tile_flags):
         if (tile_flag==1):
            if index==0:
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
            elif (tile_name[index]!=tile_name[index-1]):
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])

            
      print "On Moon tile flags from metafits, total %s:" % str(len(flag_antenna_tilenames_list))
      print flag_antenna_tilenames_list

      #off moon 
      try:
         print off_moon_metafits_filename
         HDU_list = fits.open(off_moon_metafits_filename)
      except IOError, err:
         'Cannot open metadata file %s\n' % off_moon_metafits_filename
      data = HDU_list[1].data
      tile_flags = data['FLAG']
      tile_flagging_indices = data['ANTENNA']
      tile_name = data['TILENAME']
      #print tile_flags
      
      flag_antenna_indices_list = []
      flag_antenna_tilenames_list = []
      for index,tile_flag in enumerate(tile_flags):
         if (tile_flag==1):
            if index==0:
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
            elif (tile_name[index]!=tile_name[index-1]):
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
            
      print "Off Moon tile flags from metafits, total %s:" % str(len(flag_antenna_tilenames_list))
      print flag_antenna_tilenames_list
           
      
      #get the date and time of the observation from the metafits file

      HDU_list = fits.open(metafits_filename)
      header=HDU_list[0].header 
      date_unix_utc_string = (header)['GOODTIME']
      #datetime.datetime.utcfromtimestamp(posix_time).strftime('%Y-%m-%dT%H:%M:%SZ')
      date_unix_utc=datetime.utcfromtimestamp(date_unix_utc_string)
      print "start time (GOODTIME) of observation is %s" % (date_unix_utc)
      
      #Get the observation length
      obslength=float((header)['EXPOSURE'])
      print "EXPOSURE is %s s" % obslength
      half_obslength=obslength/2.
      half_obslength_timedelta=timedelta(seconds=half_obslength)
      middle_of_observation=date_unix_utc+half_obslength_timedelta
      print "middle_of_observation is %s" % middle_of_observation
      
      middle_of_observation_formatted=middle_of_observation.strftime('%Y-%m-%dT%H:%M:%S')
      print "middle_of_observation_formatted %s " % middle_of_observation_formatted
      
      #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
      new_date=string.replace(string.replace(middle_of_observation_formatted,'-','/'),'T',' ')
      print new_date
      #find position of moon
      src_file_name='%ssrc_file_moon_%s.txt' % (on_moon_directory,on_moon_obsid)
      if machine=="namorrodor":
         cmd="/data/code/git/ben-astronomy/moon/sched/print_src.py --date='"+new_date+"' > "+ src_file_name
      else:
         cmd = "print_src.py --date='"+new_date+"' > "+ src_file_name
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
      track_off_moon_string='%s %s %s %s %s' % (on_moon_obsid,moon_ra,new_moon_dec,off_moon_obsid,LST_difference_in_sidereal_days)
      track_off_moon_string_list.append(track_off_moon_string)
      #write the track off moon string to a file in the on_moon obsid data dir where gator can find it
      track_off_moon_filename_on_moon_data_dir='%strack_off_moon_%s_%s.txt' % (on_moon_data_dir,on_moon_obsid,off_moon_obsid)
      #check if on moon directory exists and if not, make it
      if not os.path.exists(on_moon_data_dir):
         cmd="mkdir %s" % (on_moon_data_dir)
         print cmd
         os.system(cmd)
      with open(track_off_moon_filename_on_moon_data_dir, 'w') as f:
         f.write(track_off_moon_string)
      print "wrote track_off_moon file %s" % (track_off_moon_filename_on_moon_data_dir)
   with open(track_off_moon_filename, 'w') as f:
      f.write("\n".join(track_off_moon_string_list))
   print "wrote track_off_moon file %s" % (track_off_moon_filename)

   with open(paired_obs_list_filename,'w') as f:
      f.write('\n'.join(paired_obs_list))
   print "wrote paired obsid  file %s" % (paired_obs_list_filename)


#Main function:
def setup_moon_process(options):
   
   if options.copy_images_from_magnus:
      filelist = options.copy_images_from_magnus
      copy_images_from_magnus(filelist)
                  
   machine=options.machine
   if machine=='namorrodor':
      mwa_dir = '/md0/moon/data/MWA/'
      from fpdf import FPDF
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'
   base_dir=options.base_dir
   directory_of_epochs="%sepochs/" % base_dir
   moon_exp_filename=options.infile
   moon_exp_filename_base=os.path.basename(moon_exp_filename).split('.')[0]
   database_name="/group/mwaeor/bmckinley/%s.sqlite" % (moon_exp_filename_base)
   download_database_name="/group/mwaeor/bmckinley/%s_download.sqlite" % (moon_exp_filename_base)
   #get the epoch_IDs and on_moon_dat off_moon_date
   #chan_list=[69]
   chan_list=[69,93,121,145,169]
   on_moon_and_off_moon=["on_moon","off_moon"]
   epoch_ID_list=[]
   on_moon_date_list=[]
   off_moon_date_list=[]
   for line in open(moon_exp_filename):
      epoch_ID=line.split()[0]
      on_moon_date=line.split()[1]
      off_moon_date=line.split()[2]
      epoch_ID_list.append(epoch_ID)
      on_moon_date_list.append(on_moon_date)
      off_moon_date_list.append(off_moon_date)

   #check if epochs directory exists:
   if os.path.isdir(directory_of_epochs):
      pass
   else:
      cmd = "mkdir %s " % directory_of_epochs
      print cmd
      os.system(cmd)
   
   for epoch_ID_index,epoch_ID in enumerate(epoch_ID_list):
      on_moon_date=on_moon_date_list[epoch_ID_index]
      off_moon_date=off_moon_date_list[epoch_ID_index]
      print "Setting up epoch ID %s with on-Moon date %s and off-Moon date %s" % (epoch_ID,on_moon_date,off_moon_date)
      epoch_ID_dir="%s%s/" % (directory_of_epochs,epoch_ID)
      if os.path.isdir(epoch_ID_dir):
         pass
      else:
         cmd = "mkdir %s " % epoch_ID_dir
         print cmd
         os.system(cmd)
      for chan in chan_list:
         chan_dir="%s%s/" % (epoch_ID_dir,chan)
         if os.path.isdir(chan_dir):
            pass
         else:
            cmd = "mkdir %s " % chan_dir
            print cmd
            os.system(cmd)
         for on_off_moon_string in on_moon_and_off_moon:
            on_off_moon_dir="%s%s/" % (chan_dir,on_off_moon_string)
            if os.path.isdir(on_off_moon_dir):
               pass
            else:
               cmd = "mkdir %s " % on_off_moon_dir
               print cmd
               os.system(cmd)
            if not options.have_obs_lists:
               write_obs_list(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir)
            
               #if it is an off moon then write the track_off_moon txt file and the 'paired' obsid file
               if (on_off_moon_string=='off_moon'):
                  on_moon_obsid_filename="%s%s%s_on_moon_%s.txt" % (chan_dir,'on_moon/',epoch_ID,chan) 
                  off_moon_obsid_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)                
                  make_track_off_moon_file(on_moon_obsid_filename,off_moon_obsid_filename,machine,epoch_ID,chan)

            if ((machine=='magnus' or machine=='galaxy') and options.setup_gator_download):
               obsid_filename="%s%s_%s_%s.txt" % (on_off_moon_dir,epoch_ID,on_off_moon_string,chan)
               #if on_off_moon_string=='off_moon':
               #cmd='gator_add_to_downloads_table.rb -d %s %s' % (download_database_name,obsid_filename)
               #If you want to do the conversion through ASVO - this is really slow. Took over a week to convert and download one epochID
               #cmd = "gator_add_to_database.rb %s -d --db %s --conversion-options='job_type=c, timeres=8, freqres=80, edgewidth=80, conversion=ms, allowmissing=false, flagdcchannels=true, noantennapruning=true'" % (obsid_filename,download_database_name)
               #go back to just downloading gpubox files
               cmd = "gator_add_to_database.rb %s -d --db %s " % (obsid_filename,download_database_name)
               print cmd
               os.system(cmd)
              
            if options.purge_ms:
               on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
               if on_off_moon_string=='on_moon':
                  obslist_file="%s%s_on_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)
               elif on_off_moon_string=='off_moon':
                  obslist_file="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)
               else:
                  print "bad value for on_off_moon_string"
               with open(obslist_file,'r') as f:
                  lines = f.readlines()
               for line in lines:
                  obsid = line.split()[0]
                  cmd = "rm -rf %s%s/%s.ms" % (mwa_dir,obsid,obsid)
                  print cmd 
                  os.system(cmd)
                  cmd = "rm -rf %s%s/%s_ms.zip" % (mwa_dir,obsid,obsid)
                  print cmd
                  os.system(cmd)
                  
            if options.rename_ms:
               on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
               if on_off_moon_string=='off_moon':
                  continue
               elif on_off_moon_string=='on_moon':
                  obslist_file="%s%s_on_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan)
               else:
                  print "bad value for on_off_moon_string"
               with open(obslist_file,'r') as f:
                  lines = f.readlines()
               for line in lines:
                  on_moon_obsid = line.strip()[0:10]
                  off_moon_obsid = line.strip()[10:20]
                  old_on_moon_ms_name = '%s%s/%s.ms' % (mwa_dir,on_moon_obsid,on_moon_obsid)
                  old_off_moon_ms_name = '%s%s/%s.ms' % (mwa_dir,off_moon_obsid,off_moon_obsid)
                  new_on_moon_ms_name = '%s%s/%s_%s_trackmoon.ms' % (mwa_dir,on_moon_obsid,on_moon_obsid,epoch_ID)
                  new_off_moon_ms_name = '%s%s/%s_%s_track_off_moon_paired_%s.ms' % (mwa_dir,off_moon_obsid,off_moon_obsid,epoch_ID,on_moon_obsid)

                  #cmd = "rm -rf %s" % (new_on_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)
                  #cmd = "rm -rf %s" % (new_off_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)                  
                  cmd = "mv %s %s" % (old_on_moon_ms_name,new_on_moon_ms_name)
                  print cmd 
                  os.system(cmd)
                  cmd = "mv %s %s" % (old_off_moon_ms_name,new_off_moon_ms_name)
                  print cmd 
                  os.system(cmd)
                  #cmd = "ls %s" % (old_on_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)
                  #cmd = "ls %s" % (old_off_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)
                  #cmd = "ls %s" % (new_on_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)
                  #cmd = "ls %s" % (new_off_moon_ms_name)
                  #print cmd 
                  #os.system(cmd)

            if (options.ms_quality_check and machine=="namorrodor"):
               #need to be sshfs'd in:sshfs galaxy:/ /md0/galaxy/ 
               galaxy_mount_dir = "/md0/galaxy"
               galaxy_mwa_dir = "/astro/mwaeor/MWA/data/"
               on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
               if on_off_moon_string=='off_moon':
                  continue
               elif on_off_moon_string=='on_moon':
                  obslist_file="%s%s_on_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan)
               else:
                  print "bad value for on_off_moon_string"
               with open(obslist_file,'r') as f:
                  lines = f.readlines()
               for line_index,line in enumerate(lines):
                  on_moon_obsid = line.strip()[0:10]
                  off_moon_obsid = line.strip()[10:20]
                  on_moon_ms_name = '%s%s%s/%s_%s_trackmoon.ms' % (galaxy_mount_dir,galaxy_mwa_dir,on_moon_obsid,on_moon_obsid,epoch_ID)
                  off_moon_ms_name = '%s%s%s/%s_%s_track_off_moon_paired_%s.ms' % (galaxy_mount_dir,galaxy_mwa_dir,off_moon_obsid,off_moon_obsid,epoch_ID,on_moon_obsid)

                  galaxy_on_moon_data_dir = "%s%s/" % (galaxy_mwa_dir,on_moon_obsid)
                  galaxy_off_moon_data_dir = "%s%s/" % (galaxy_mwa_dir,off_moon_obsid)
                  #apply the same flags to both ms for moon obs, do this in the on moon stage and get flags from metafits 
                  #get flagged antennas from the metafits file and flag them explicitly
                  on_moon_metafits_file_name = "%s%s%s_metafits_ppds.fits" % (galaxy_mount_dir,galaxy_on_moon_data_dir,on_moon_obsid)
                  off_moon_metafits_file_name = "%s%s%s_metafits_ppds.fits" % (galaxy_mount_dir,galaxy_off_moon_data_dir,off_moon_obsid)
                  
                  if options.flag_ants:
                        #Flag additional antennas (do both as you may want to just flag ants without do_flagging)
                        cmd="flagantennae %s %s" % (on_moon_ms_name,options.flag_ants)
                        print cmd
                        os.system(cmd)
                        
                        cmd="flagantennae %s %s" % (off_moon_ms_name,options.flag_ants)
                        print cmd
                        os.system(cmd)
                        
                  if options.do_flagging:
                     #on moon
                     try:
                        print on_moon_metafits_file_name
                        HDU_list = fits.open(on_moon_metafits_file_name)
                     except IOError, err:
                        'Cannot open metadata file %s\n' % on_moon_metafits_file_name
                     data = HDU_list[1].data
                     tile_flags = data['FLAG']
                     tile_flagging_indices = data['ANTENNA']
                     tile_name = data['TILENAME']
                     #print tile_flags
                     
                     flag_antenna_indices_list = []
                     flag_antenna_tilenames_list = []
                     for index,tile_flag in enumerate(tile_flags):
                        if (tile_flag==1):
                           string = "%s with antenna index %s should be flagged." % (tile_name[index],tile_flagging_indices[index])
                           #print string
                           if index==0:
                              flag_antenna_indices_list.append(tile_flagging_indices[index])
                              flag_antenna_tilenames_list.append(tile_name[index])
                           elif (tile_name[index]!=tile_name[index-1]):
                              flag_antenna_indices_list.append(tile_flagging_indices[index])
                              flag_antenna_tilenames_list.append(tile_name[index])
                     
                     print flag_antenna_indices_list
                     print flag_antenna_tilenames_list
                     
                     flag_ant_string = ' '.join(str(x) for x in flag_antenna_indices_list)
                     
                     #flag em
                     cmd="flagantennae %s %s" % (on_moon_ms_name,flag_ant_string)
                     print cmd
                     os.system(cmd)
                     
                     #off moon
                     try:
                        HDU_list = fits.open(off_moon_metafits_file_name)
                     except IOError, err:
                        'Cannot open metadata file %s\n' % off_moon_metafits_file_name
                     header=HDU_list[0].header
                     data = HDU_list[1].data
                     tile_flags = data['FLAG']
                     tile_flagging_indices = data['ANTENNA']
                     tile_name = data['TILENAME']
                     #print tile_flags
                     
                     flag_antenna_indices_list = []
                     flag_antenna_tilenames_list = []
                     for index,tile_flag in enumerate(tile_flags):
                        if (tile_flag==1):
                           string = "%s with antenna index %s should be flagged." % (tile_name[index],tile_flagging_indices[index])
                           #print string
                           if index==0:
                              flag_antenna_indices_list.append(tile_flagging_indices[index])
                              flag_antenna_tilenames_list.append(tile_name[index])
                           elif (tile_name[index]!=tile_name[index-1]):
                              flag_antenna_indices_list.append(tile_flagging_indices[index])
                              flag_antenna_tilenames_list.append(tile_name[index])
                     
                     print flag_antenna_indices_list
                     print flag_antenna_tilenames_list
                     
                     flag_ant_string = ' '.join(str(x) for x in flag_antenna_indices_list)
                     
                     #flag em
                     cmd="flagantennae %s %s" % (off_moon_ms_name,flag_ant_string)
                     print cmd
                     os.system(cmd)
                     
                     
                     with table(on_moon_ms_name,readonly=False) as on_moon_table:
                        with table(off_moon_ms_name,readonly=False) as off_moon_table:
                           with tablecolumn(on_moon_table,'FLAG') as on_moon_table_flag:
                              with tablecolumn(off_moon_table,'FLAG') as off_moon_table_flag:
                                  #new_moon_flag=np.empty_like(on_moon_table_flag)
                                  for row_index,row in enumerate(on_moon_table_flag):
                                     #if (on_moon_table_flag[row_index].all() != off_moon_table_flag[row_index]).all():
                                     #   print on_moon_table_flag[row_index]
                                     #   print off_moon_table_flag[row_index]
                     
                                     on_moon_table_flag[row_index] = np.logical_or(on_moon_table_flag[row_index],off_moon_table_flag[row_index])
                                     off_moon_table_flag[row_index] = np.logical_or(on_moon_table_flag[row_index],off_moon_table_flag[row_index])
                     
                                     on_moon_table.putcell('FLAG',row_index,on_moon_table_flag[row_index])
                                     off_moon_table.putcell('FLAG',row_index,off_moon_table_flag[row_index])
                     
                  
                  #collect new statistics
                  if options.post_cal_quality_check:
                     #only do the aoquality stuff for every third obsid (it takes a long time)
                     if (line_index % 3 == 0):
                        cmd = "aoquality collect -d CORRECTED_DATA %s" % (on_moon_ms_name)
                        print cmd
                        os.system(cmd)
                  
                        cmd = "aoquality collect -d CORRECTED_DATA %s" % (off_moon_ms_name)
                        print cmd
                        os.system(cmd)                     
                  
                     #plot cal solutions:
                     on_moon_solutions_name = '%s%s%s/%s_%s_selfcal_0_concat_solutions_trackmoon.bin' % (galaxy_mount_dir,galaxy_mwa_dir,on_moon_obsid,on_moon_obsid,epoch_ID)
                     off_moon_solutions_name = '%s%s%s/%s_%s_selfcal_0_concat_solutions_track_off_moon_paired_%s.bin' % (galaxy_mount_dir,galaxy_mwa_dir,off_moon_obsid,off_moon_obsid,epoch_ID,on_moon_obsid)
                     
                     cmd = "aocal_plot.py --outdir=%s --title=%s %s" % (base_dir,on_moon_obsid,on_moon_solutions_name)
                     print cmd
                     os.system(cmd)
                     
                     cmd = "aocal_plot.py --outdir=%s --title=%s %s" % (base_dir,off_moon_obsid,off_moon_solutions_name)
                     print cmd
                     os.system(cmd) 
                  
                  #don't need to collect stats as the tables are already written (unaveraged))
                  #else:
                  #   if (line_index % 3 == 0):
                  #      cmd = "aoquality collect %s" % (on_moon_ms_name)
                  #      print cmd
                  #      os.system(cmd)
                  # 
                  #      cmd = "aoquality collect %s" % (off_moon_ms_name)
                  #      print cmd
                  #      os.system(cmd)
                  
                  
                  #save a png of the aoqplot output for StandardDeviation
                  
                  
                  
                  if options.post_cal_quality_check:
                     if (line_index % 3 == 0):
                        cmd = "aoqplot -save aoqplot_stddev_%s_on_moon_%s_corrected_data StandardDeviation %s" % (on_moon_obsid,epoch_ID,on_moon_ms_name)
                        print cmd
                        os.system(cmd)
                  
                        cmd = "aoqplot -save aoqplot_stddev_%s_off_moon_%s_corrected_data StandardDeviation %s" % (off_moon_obsid,epoch_ID,off_moon_ms_name)
                        print cmd
                        os.system(cmd)

                  else:
                     #if (line_index % 3 == 0):
                     cmd = "aoqplot -save aoqplot_stddev_%s_on_moon_%s StandardDeviation %s" % (on_moon_obsid,epoch_ID,on_moon_ms_name)
                     print cmd
                     os.system(cmd)
                  
                     cmd = "aoqplot -save aoqplot_stddev_%s_off_moon_%s StandardDeviation %s" % (off_moon_obsid,epoch_ID,off_moon_ms_name)
                     print cmd
                     os.system(cmd)
                                 
                  
            #can only do this here if you have already made the obslists - redundant now
            #as paired obslist written in make_track_off_moon_file
            #else:
            #   write_paired_obslist(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir)

            if not (machine=='magnus' or machine=='galaxy'):
               write_and_run_default_scripts(epoch_ID,chan,on_off_moon_dir,machine)
            elif (machine=='magnus' or machine=='galaxy') and options.run_gator:
               pass
            else:
               print "Not running gator"

         if (machine=='magnus' or machine=='galaxy'):
            #launch_gator(epoch_ID,chan,chan_dir)   
            #add the paired obs to the database
            paired_observations_filename="%son_moon/%s_on_moon_paired_%s.txt" % (chan_dir,epoch_ID,chan)
            cmd='gator_add_to_rts_table.rb -d %s --epoch_id %s %s' % (database_name,epoch_ID,paired_observations_filename)
            print cmd
            os.system(cmd)

   #combine qa into two pdfs (on and off moon)  here
   if (options.ms_quality_check and machine=="namorrodor"):
      if options.post_cal_quality_check:
         figname = "aoqplot_std_dev_baselines_combined_on_moon_%s_corrected_data.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *on_moon_%s_corrected_data*baselines.pdf -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)

         figname = "aoqplot_std_dev_baselines_combined_off_moon_%s_corrected_data.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *off_moon_%s_corrected_data*baselines.pdf -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)

         #also combine cal solution plots
         figname = "phase_solutions_combined_on_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *%s*trackmoon*_phase.png -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)
         
         figname = "amp_solutions_combined_on_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *%s*trackmoon*_amp.png -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)

         figname = "phase_solutions_combined_off_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *%s*track_off_moon*_phase.png -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)
         
         figname = "amp_solutions_combined_off_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *%s*track_off_moon*_amp.png -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)
           
      else:
         figname = "aoqplot_std_dev_baselines_combined_on_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *on_moon_%s*baselines.pdf -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)

         figname = "aoqplot_std_dev_baselines_combined_off_moon_%s.pdf" % epoch_ID
         pdf = FPDF('P', 'mm', 'A4')
         cmd = '/usr/bin/montage *off_moon_%s*baselines.pdf -mode concatenate -tile 3x5 %s' % (epoch_ID,figname)
         print cmd
         os.system(cmd)
               
   if (machine=='magnus' or machine=='galaxy') and options.ms_download:
      #launch gator
      cmd='nohup gator_rts_daemon.rb -d %s --cotter --ms_download &' % (database_name)
      print cmd
      os.system(cmd)

   elif (machine=='magnus' or machine=='galaxy') and options.run_gator:
      #launch gator
      cmd='nohup gator_rts_daemon.rb -d %s --cotter &' % (database_name)
      print cmd
      os.system(cmd)   
   else:
      print "not launching gator"

from optparse import OptionParser,OptionGroup

usage = 'Usage: setup_moon_process.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--infile',type='string', dest='infile',default='',help='Just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date (dates as YYYY-MM-DD) e.g. --infile="moon_experiment.txt" [default=%default]')
parser.add_option('--base_dir',type='string', dest='base_dir',default='/astro/mwaeor/bmckinley/moon/test_setup/',help='Base directory to set everything up in e.g. --base_dir="/md0/moon/magnus_setup_tests/" [default=%default]')
parser.add_option('--have_obs_lists',action='store_true',dest='have_obs_lists',default=False,help='Set if you already have all the obs lists (dont want to run find_observations.py)[default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help=' e.g. --machine="magnus" [default=%default]')
parser.add_option('--run_gator',action='store_true',dest='run_gator',default=False,help='Set if you want to run gator [default=%default]')
parser.add_option('--ms_download',action='store_true',dest='ms_download',default=False,help='Set if you want to run gator to download ms only [default=%default]')
parser.add_option('--cleanup',action='store_true',dest='cleanup',default=False,help='Set to delete gpubox files after converting to ms in cotter mode [default=%default]')
parser.add_option('--setup_gator_download',action='store_true',dest='setup_gator_download',default=False,help='Set up the gator download on magnus or galaxy. To start download: nohup [default=%default]')
parser.add_option('--copy_images_from_magnus',type='string', dest='copy_images_from_magnus',default=None,help='Give a file containing the paths to the folders where images are to be copied from e.g. --copy_images_from_magnus="trackmoon_file_list.txt" [default=%default]')
parser.add_option('--purge_ms',action='store_true',dest='purge_ms',default=False,help='Set to delete all ms and zip files to start again [default=%default]')
parser.add_option('--rename_ms',action='store_true',dest='rename_ms',default=False,help='Set to rename ms for moon stuff [default=%default]')
parser.add_option('--ms_quality_check',action='store_true',dest='ms_quality_check',default=False,help='Set to run aoquality and aoqplot on each ms and make a single pdf (need to be sshfs in:sshfs galaxy:/ /md0/galaxy/) [default=%default]')
parser.add_option('--post_cal_quality_check',action='store_true',dest='post_cal_quality_check',default=False,help='Must use with ms_qualtiy_check. Run checks on calibrated data and plot solutions and make a single pdf (need to be sshfs in:sshfs galaxy:/ /md0/galaxy/) [default=%default]')
parser.add_option('--do_flagging',action='store_true',dest='do_flagging',default=False,help='Set to do flagging ie flag ants from metafits and combine on/off moon flags (need to be sshfs in:sshfs galaxy:/ /md0/galaxy/) [default=%default]')
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='flag these antennas (zero-indexed andre style numbering), space separated, used with ms_quality_check only e.g. --flag_ants="56 60" [default=%default]')



#parser.add_option('--launch_cotter',action='store_true',dest='launch_cotter',default=False,help='Actually launch the cotter jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
#parser.add_option('--launch_selfcal',action='store_true',dest='launch_selfcal',default=False,help='Actually launch the selfcal jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
#parser.add_option('--launch_image',action='store_true',dest='launch_image',default=False,help='Actually launch the imaging jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
#parser.add_option('--launch_pbcorr',action='store_true',dest='launch_pbcorr',default=False,help='Actually launch the pbcorr jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')


(options, args) = parser.parse_args()

setup_moon_process(options)

