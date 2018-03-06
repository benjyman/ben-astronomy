#script to set up Moon processing in bulk (on magnus)!
#just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date
import sys,os
import cmd
from datetime import datetime, timedelta, date, time

def write_obs_lists(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir):
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
   
   
def write_default_scripts(epoch_ID,chan,on_off_moon_dir,machine):
   #function to write the default scripts: download, cotter, calibrate, image, pbcorr
   on_off_moon_string=on_off_moon_dir.strip().split('/')[-2]
   default_download_script_name="%s1_default_download_%s_%s_%s.txt" % (on_off_moon_dir,epoch_ID,chan,on_off_moon_string)
   default_download_file = open(default_download_script_name,'w+')
   generate_download_string='python /data/code/git/ben-astronomy/moon/processing_scripts/namorrodor/generate_obs_download.py /md0/moon/epochs/2018A_01/93/on_moon/1199392480.txt'
   default_download_file.write(generate_download_string)
   
   
   
   
   
def setup_moon_process(options):
   machine='magnus'
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
            write_obs_lists(epoch_ID,on_moon_date,off_moon_date,chan,on_off_moon_dir)


from optparse import OptionParser,OptionGroup

usage = 'Usage: setup_moon_process.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--infile',type='string', dest='infile',default='',help='Just give it a txt file with three columns: epoch_ID on_moon_date off_moon_date (dates as YYYY-MM-DD) e.g. --infile="moon_experiment.txt" [default=%default]')
parser.add_option('--base_dir',type='string', dest='base_dir',default='/md0/moon/magnus_setup_tests/',help='Base directory to set everything up in e.g. --base_dir="/md0/moon/magnus_setup_tests/" [default=%default]')

(options, args) = parser.parse_args()

setup_moon_process(options)

