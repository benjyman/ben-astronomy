#script to set up Moon processing in bulk (on magnus)!
#just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date
import sys,os
import cmd
from datetime import datetime, timedelta, date, time
from astropy.io import fits
import string
import glob


sidereal_day_sec = 86164
sourcelist='srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt'
havebeam_string = ''

#Functions
def write_obs_list(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir):
   on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
   print on_off_moon_string
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
      sister_obsid_list_filename="../off_moon/%s_off_moon_%s.txt" % (epoch_ID,chan) 
      paired_obs_list_filename="%s%s_on_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   elif on_off_moon_string=='off_moon':
      obsid_list_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
      sister_obsid_list_filename="../on_moon/%s_on_moon_%s.txt" % (epoch_ID,chan) 
      paired_obs_list_filename="%s%s_off_moon_paired_%s.txt" % (on_off_moon_dir,epoch_ID,chan) 
   else:
      print "Bad values for on/off moon"      
   for line in open(obsid_list_filename):
      obsid_list.append(line.strip()) 
   for line in open(sister_obsid_list_filename):
      sister_obsid_list.append(line.strip())    
   for item_index,item in enumerate(obsid_list):
      paired_obs_string="%s_%s" % (str(obsid_list[item_index]),str(sister_obsid_list[item_index]))
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
   generate_download_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_obs_download.py --machine=%s %s' % (ben_code_base,machine,observations_filename)
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
   #generate_cotter_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_cotter_moon.py --epoch_ID=%s %s --flag_ants="" --cleanup --obsid_infile=%s --sister_obsid_infile=%s' % (ben_code_base,epoch_ID,track_moon_string,observations_filename,sister_observations_filename)
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
   generate_pbcorr_string='python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/generate_qpbcorr_multi.py  --epoch_ID=%s %s  --obsid_infile=%s --dirty  --channelsout=24 %s   ' % (ben_code_base,epoch_ID,track_moon_string,observations_filename, havebeam_string) 
   with open(default_pbcorr_script_name,'w+') as f:
      f.write('#!/bin/bash -l\n')
      f.write(generate_pbcorr_string)
   cmd = "chmod +x %s" % default_pbcorr_script_name
   #print cmd
   os.system(cmd)
   cmd = "%s" % default_pbcorr_script_name
   print cmd
   os.system(cmd)
   
   #launch the jobs if requested
   launch_job_filename='%slaunch_jobs.sh'%on_off_moon_dir
   with open(launch_job_filename,'w') as f:
      if options.launch_cotter:
         q_cotter_files_list=glob.glob('%s/q_cotter_moon*'%on_off_moon_dir)
         number_of_q_cotter_files=len(q_cotter_files_list)
         print 'number_of_q_cotter_files is %s ' % number_of_q_cotter_files
         if number_of_q_cotter_files==0:
            print 'No q_cotter files '
         elif number_of_q_cotter_files==1:
            cmd1='jobid=`sbatch %sq_cotter_moon_0.sh | cut -d " " -f 4`' % on_off_moon_dir
            #print cmd1
            #os.system(cmd1)
            f.write(cmd1)
         else:
            for q_cotter_file_index,q_cotter_file in enumerate(q_cotter_files_list):
               cmd1='jobid_%s=`sbatch %s | cut -d " " -f 4`' % (str(q_cotter_file_index),q_cotter_file)
               f.write(cmd1)
         #print 'WARNING! Too many q_cotter files - now you need to write this bit of code!'
      if options.launch_selfcal:
         if options.launch_cotter:
            if number_of_q_cotter_files==0:
               print 'No q_cotter files '
            elif number_of_q_cotter_files==1:
               cmd2='jobid=`sbatch --dependency=afterok:$jobid %sq_selfcal_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
               f.write(cmd2)
            else:
               dependency_string='--dependency=afterok'
               for q_cotter_file_index,q_cotter_file in enumerate(q_cotter_files_list):
                  dependency_string+=':jobid_%s' % (str(q_cotter_file_index))
               cmd2='jobid=`sbatch %s %sq_selfcal_moon.sh | cut -d " " -f 4`' % (dependency_string,on_off_moon_dir)
               f.write(cmd2)          
               #print 'WARNING! Too many q_cotter files - now you need to write this bit of code!'
         else:
            cmd2='jobid=`sbatch %sq_selfcal_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
         #print cmd2
         #os.system(cmd2)
         
      if options.launch_image:
         if options.launch_selfcal:
            cmd3='jobid=`sbatch --dependency=afterok:$jobid %sq_image_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
         else:
            cmd3='jobid=`sbatch %sq_image_moon.sh | cut -d " " -f 4`' % on_off_moon_dir
         #print cmd3
         #os.system(cmd3)
         f.write(cmd3)
      if options.launch_pbcorr:
         if options.launch_image:
            cmd4='sbatch --dependency=afterok:$jobid %sq_pbcorr_moon.sh' % on_off_moon_dir
         else:
            cmd4='sbatch %sq_pbcorr_moon.sh' % on_off_moon_dir
         #print cmd4
         #os.system(cmd4)
         f.write(cmd4)
   cmd='chmod +x %s' % launch_job_filename
   print cmd
   os.system(cmd)
   cmd='%s' % launch_job_filename
   print cmd
   os.system(cmd)   
   

   
            
              
            
def make_track_off_moon_file(on_moon_obsid_filename,off_moon_obsid_filename):
   #make the file that has the sister moon obsid moon ra (hh:mm:ss.ss)  moon dec (dd.mm.ss.s) and off_moon obsid
   on_moon_directory=os.path.dirname(on_moon_obsid_filename)+'/'
   off_moon_directory=os.path.dirname(off_moon_obsid_filename)+'/'
   off_moon_obsid_filename_no_path=off_moon_obsid_filename.split(off_moon_directory)[1]
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]
   track_off_moon_string_list=[]
   track_off_moon_filename='%strack_off_moon_%s' % (off_moon_directory,off_moon_obsid_filename_no_path)
   for line in open(on_moon_obsid_filename):
      on_moon_obsid_list.append(line.strip()) 
   for line in open(off_moon_obsid_filename):
      off_moon_obsid_list.append(line.strip()) 
   n_obs = sum(1 for line in open(on_moon_obsid_filename))
   
   for obsid_index,on_moon_obsid in enumerate(on_moon_obsid_list): 
      off_moon_obsid=off_moon_obsid_list[obsid_index]
      #Check LSTs
      LST_difference=float(off_moon_obsid)-float(on_moon_obsid)
      LST_remainder = LST_difference % sidereal_day_sec
      LST_difference_in_sidereal_days=abs(LST_remainder/sidereal_day_sec)
      print "LST difference is: %s in sidereal_days" % LST_difference_in_sidereal_days
      
      metafits_filename='%s%s_metafits_ppds.fits' % (on_moon_directory,on_moon_obsid)
      #download the metafits file
      cmd='wget -O %s http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s' % (metafits_filename,on_moon_obsid)
      print cmd
      os.system(cmd)
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
      track_off_moon_string='%s %s %s %s %s' % (on_moon_obsid,moon_ra,new_moon_dec,off_moon_obsid,LST_difference_in_sidereal_days)
      track_off_moon_string_list.append(track_off_moon_string)
      with open(track_off_moon_filename, 'w') as f:
         f.write("\n".join(track_off_moon_string_list))
      
      

#Main function:
def setup_moon_process(options):
   machine=options.machine
   base_dir=options.base_dir
   directory_of_epochs="%sepochs/" % base_dir
   moon_exp_filename=options.infile
   #get the epoch_IDs and on_moon_dat off_moon_date
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
            
               #if it is an off moon then write the track_off_moon txt file
               if (on_off_moon_string=='off_moon'):
                  on_moon_obsid_filename="%s%s%s_on_moon_%s.txt" % (chan_dir,'on_moon/',epoch_ID,chan) 
                  off_moon_obsid_filename="%s%s_off_moon_%s.txt" % (on_off_moon_dir,epoch_ID,chan)                
                  make_track_off_moon_file(on_moon_obsid_filename,off_moon_obsid_filename)
            
            write_paired_obslist(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir)
            
            write_and_run_default_scripts(epoch_ID,chan,on_off_moon_dir,machine)
            if (options.setup_gator_download and (machine=='magnus' or machine=='galaxy')):
               #download_script_directory=os.path.dirname(default_download_script_name)+'/'
               download_script_directory=on_off_moon_dir
               queue_download_script_name='q_obsdownload_wrapper.sh' 
               #cmd = "mv %s %s" % (queue_download_script_name,download_script_directory)
               #print cmd
               #os.system(cmd)
               
               cmd = "chmod +x %s%s" % (download_script_directory,queue_download_script_name)
               #print cmd
               os.system(cmd)
               cmd = "%s%s" % (download_script_directory,queue_download_script_name)
               print cmd
               os.system(cmd)              
              
                  
            

from optparse import OptionParser,OptionGroup

usage = 'Usage: setup_moon_process.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--infile',type='string', dest='infile',default='',help='Just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date (dates as YYYY-MM-DD) e.g. --infile="moon_experiment.txt" [default=%default]')
parser.add_option('--base_dir',type='string', dest='base_dir',default='/astro/mwaeor/bmckinley/moon/test_setup/',help='Base directory to set everything up in e.g. --base_dir="/md0/moon/magnus_setup_tests/" [default=%default]')
parser.add_option('--have_obs_lists',action='store_true',dest='have_obs_lists',default=False,help='Set if you already have all the obs lists (dont want to run find_observations.py)[default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help=' e.g. --machine="magnus" [default=%default]')
parser.add_option('--setup_gator_download',action='store_true',dest='setup_gator_download',default=False,help='Set up the gator download on magnus or galaxy. To start download: nohup [default=%default]')
parser.add_option('--launch_cotter',action='store_true',dest='launch_cotter',default=False,help='Actually launch the cotter jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
parser.add_option('--launch_selfcal',action='store_true',dest='launch_selfcal',default=False,help='Actually launch the selfcal jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
parser.add_option('--launch_image',action='store_true',dest='launch_image',default=False,help='Actually launch the imaging jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')
parser.add_option('--launch_pbcorr',action='store_true',dest='launch_pbcorr',default=False,help='Actually launch the pbcorr jobs (sbatch to queue on HPC) - otherwise just sets everything up [default=%default]')


(options, args) = parser.parse_args()

setup_moon_process(options)

