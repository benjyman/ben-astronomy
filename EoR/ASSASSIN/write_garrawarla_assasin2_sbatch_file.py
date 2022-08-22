#!/usr/bin/env python3

import numpy as np
import math
import os
import ephem

fine_chans_per_EDA2_chan = 27
fine_chan_width_Hz =  (400./512.) / 27. * 1.0e6 #(see Marcin email)
fine_chan_width_MHz = fine_chan_width_Hz / 1.0e6

mwa_latitude_ephem = '-26.70331940'
mwa_longitude_ephem = '116.670575'

def get_eda2_lst(eda_time_string):
   year, month, day, hour, minute, second = eda_time_string[0:4], eda_time_string[4:6],eda_time_string[6:8], eda_time_string[9:11],eda_time_string[11:13],eda_time_string[13:15]
   eda2_observer = ephem.Observer()
   eda2_observer.lon, eda2_observer.lat = mwa_longitude_ephem, mwa_latitude_ephem
   eda2_observer.date = '%s/%s/%s %s:%s:%s' % (year,month,day,hour,minute,second)
   eda2_obs_lst = (eda2_observer.sidereal_time()) 
   #print("LST is")
   #print(eda2_obs_lst)
   eda2_obs_lst_hrs = eda2_obs_lst / 2 / np.pi * 24.
   
   return(eda2_obs_lst_hrs)

def make_EDA2_obs_time_list_each_chan(base_dir,eda2_chan_list,ms=False):
   obs_time_list_each_chan = []
   if ms==False:
      temp_txt_filename = 'uvfits_time_list.txt'
      for eda2_chan in eda2_chan_list:
         chan_obs_time_list = []
         chan_dir = "%s%s/" % (base_dir,eda2_chan)
         cmd = "ls -la %schan_%s_*.uvfits  > %s" % (chan_dir,eda2_chan,temp_txt_filename)
         os.system(cmd)
         with open(temp_txt_filename) as f:
            lines=f.readlines()
         for line in lines[1:-1]:
            obs_time = line.split('.uvfits')[0].split()[-1].split('_')[-1]
            chan_obs_time_list.append(obs_time)
         #there might be no obs, if so just put in a zero
         if len(chan_obs_time_list)!=0:
            obs_time_list_each_chan.append(chan_obs_time_list)
         else:
            obs_time_list_each_chan.append([0])
   else:
      temp_txt_filename = 'ms_time_list.txt'
      for eda2_chan in eda2_chan_list:
         chan_obs_time_list = []
         chan_dir = "%s%s/" % (base_dir,eda2_chan)
         cmd = "ls %s*avg8140.ms | grep '.ms' > %s" % (chan_dir,temp_txt_filename)
         os.system(cmd)
         with open(temp_txt_filename) as f:
            lines=f.readlines()
         for line in lines[1:-1]:
            obs_time = line.split('.ms')[0].split()[-1].split('_')[-6]
            obs_date = line.split('.ms')[0].split()[-1].split('_')[-7].split('/')[-1]
            obs_date_time = "%sT%s" % (obs_date,obs_time)
            chan_obs_time_list.append(obs_date_time)
         #there might be no obs, if so just put in a zero
         if len(chan_obs_time_list)!=0:
            obs_time_list_each_chan.append(chan_obs_time_list)
         else:
            obs_time_list_each_chan.append([0])
   return obs_time_list_each_chan
   
def calculate_freq_MHz_fine_chan_subarray(EDA2_chan):
   freq_MHz_centre = (400./512.)*EDA2_chan
   fine_chan_index_array = np.arange(fine_chans_per_EDA2_chan) - 14
   freq_MHz_fine_chan_sub_array = fine_chan_index_array*fine_chan_width_MHz + freq_MHz_centre + (fine_chan_width_MHz/2.)
   return freq_MHz_fine_chan_sub_array
   
def write_garrawarla_assassin2_sbatch_file(EDA2_chan_list,lst_hrs_list,EDA2_obs_time_list,git_repo_dir="/astro/mwaeor/bmckinley/code/"):
   base_freq_MHz = EDA2_chan_list[0] * (400./512.)
   
   launch_sbatch_jobs_filename = "launch_sbatch_jobs.sh"
   sbatch_filename_list = []
   
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):
      
      if (EDA2_chan_index==0 and math.isclose(base_freq_MHz,50.0)):
         array_job_start = 14
         array_job_end = 26 
      else:
         array_job_start = 0
         array_job_end = 26
      
      EDA2_obs_time = EDA2_obs_time_list[EDA2_chan_index]
      lst_hrs = float(lst_hrs_list[EDA2_chan_index])
      
      #year,month,day,hour,min,sec = time_string.split('_')
      #time_formatted = '%d-%02d-%02dT%02d:%02d:%02d' % (float(year),float(month),float(day),float(hour),float(min),float(sec))

      name_base = "assassin2_sbatch_eda2_chan_%03d_%s" % (EDA2_chan,EDA2_obs_time)  
      sbatch_filename = "%s.sh" % name_base
      
      sbatch_filename_list.append(sbatch_filename)
   
      with open('%s' % sbatch_filename,'w') as outfile:
         outfile.write("#!/bin/bash --login\n#SBATCH --nodes=1\n#SBATCH --partition=workq\n")
         outfile.write("#SBATCH --time=00:30:00\n#SBATCH --account=mwaeor\n#SBATCH --nodes=1\n")
         outfile.write("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --array=%s-%s\n\n" % (str(array_job_start),str(array_job_end)))
         

         outfile.write("module use /pawsey/mwa/software/python3/modulefiles\n\n")
         #outfile.write("module load json-c/0.14\nmodule load hdf5/1.10.5\nmodule load cfitsio/3.48\nmodule load cmake/3.15.0\n")
         #outfile.write("module load cuda/10.2\nmodule load pal/0.9.8\nmodule load python/3.8.2\nmodule load astropy/4.0.1.post1\n\n")
         outfile.write("module load module load python/3.8.2\nmodule load astropy/4.0.1.post1\nmodule load healpy/1.13.0\n\n")
         outfile.write("module load scipy/1.4.1\nmodule load h5py/2.10.0\nmodule load healpy/1.13.0\nmodule load numpy/1.18.2\n\n")
         outfile.write("module load ephem/3.7.7.1\nmodule load matplotlib/3.2.1\n\n")
   
         #outfile.write("source /astro/mwaeor/jline/software/WODEN/build/init_WODEN.sh\n")
         #outfile.write("export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH\n\n")
         outfile.write("time python3 %sgarrawarla_sim_with_complex_beams.py --EDA2_chan=%s --EDA2_chan_index=%s --lst_hrs=%s --EDA2_obs_time=%s --freq_MHz_fine_chan_index=${SLURM_ARRAY_TASK_ID} --beam_dir=%s --git_repo_dir=%s\n" % (git_repo_dir,EDA2_chan,EDA2_chan_index,lst_hrs,EDA2_obs_time,beam_dir,git_repo_dir))
        
                     
      print("wrote %s" % sbatch_filename) 
   
   with open('%s' % launch_sbatch_jobs_filename,'w') as outfile:
      outfile.write("#!/bin/bash\n\n")
      for sbatch_filename in sbatch_filename_list:
         outfile.write("sbatch %s\n\n" % sbatch_filename)
      
   
   
 
if __name__ == "__main__":
    import argparse
    
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)
            
            from argparse import RawTextHelpFormatter
            
    parser = argparse.ArgumentParser(description="Run this to make assassin2 sims and \
            run global signal analysis on Garrawarla with array jobs",formatter_class=SmartFormatter)

    parser.add_argument('--EDA2_chan_list_start_and_end', default=None,
        help='first and last (inclusive) EDA2 chans to simulate. Defines list of EDA2 coarse chans. Generate a script for each EDA2 chan that then runs an array job for each fine chan. e.g. --EDA2_chan_list_start_and_end="64,126"')
        
    #parser.add_argument('--EDA2_obs_time_list', default=None,
    #    help='list of EDA2 obs times that correspond to each EDA2 coarse chan. \
    #         must be same length as EDA2_chan_list and for data calibration must correspond to ms filenames e.g. --EDA2_obs_time_list="20200303T000000,20200303T000400"')

    #parser.add_argument('--lst_hrs_list', default=None,
    #    help='list of LSTs that corresponds to the EDA2 obs time list. e.g. --lst_hrs_list="6.0,6.2"')
   
    parser.add_argument('--EDA2_data_dir', default='/astro/mwaeor/bmckinley/EoR/EDA2/20200303_data/',
        help='directory where EDA2_data is kept --EDA2_data_dir="/astro/mwaeor/bmckinley/EoR/EDA2/20200303_data/"')

    parser.add_argument('--git_repo_dir', default='/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/',
        help='directory where my code is kept --git_repo_dir="/astro/mwaeor/bmckinley/code/"')

    parser.add_argument('--beam_dir', default='/astro/mwaeor/bmckinley/EoR/EDA2/EEPs/',
        help='directory where beam files are kept --git_repo_dir="/astro/mwaeor/bmckinley/EoR/EDA2/EEPs/"')
        
    args = parser.parse_args()
    
    if args.EDA2_chan_list_start_and_end:
       EDA2_chan_list_start_and_end = [s.strip() for s in args.EDA2_chan_list_start_and_end.split(",")]
       EDA2_chan_start = int(EDA2_chan_list_start_and_end[0])
       EDA2_chan_end = int(EDA2_chan_list_start_and_end[1])
       EDA2_chan_list = range(EDA2_chan_start,EDA2_chan_end+1)
       
    #if args.EDA2_obs_time_list:
    #   EDA2_obs_time_list = [s.strip() for s in args.EDA2_obs_time_list.split(",")]
       
    #if args.lst_hrs_list:
    #   lst_hrs_list_strings = [s.strip() for s in args.lst_hrs_list.split(",")]
    #   lst_hrs_list = [float(x) for x in lst_hrs_list_strings] 

    if args.EDA2_data_dir:
       EDA2_data_dir = args.EDA2_data_dir
              
    if args.git_repo_dir:
       git_repo_dir = args.git_repo_dir
       
    if args.beam_dir:
       beam_dir = args.beam_dir

    #times
    ms=True
    EDA2_obs_time_list_each_chan = make_EDA2_obs_time_list_each_chan(EDA2_data_dir,EDA2_chan_list,ms=ms)
    EDA2_obs_time_list_each_chan = EDA2_obs_time_list_each_chan[0:]
    EDA2_obs_time_list = [item[0] for item in EDA2_obs_time_list_each_chan] 

    #LSTs
    lst_hrs_list = []
    for EDA2_obs_time_index,EDA2_obs_time in enumerate(EDA2_obs_time_list):
       #there might have been no obs:
       if EDA2_obs_time!=0:
          #print(EDA2_obs_time)
          lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time)
          lst_hrs_list.append(lst_eda2_hrs)
       else:
          lst_eda2_hrs = "%0.5f" % get_eda2_lst(EDA2_obs_time_list[0])
          lst_hrs_list.append(lst_eda2_hrs)
          
        ##this is for LST 60.0 deg, will do rest of lsts later
        #year,month,day,hour,min,sec = 2015,11,29,15,40,29 #LST=60 deg
        #time_string = '%d_%02d_%02d_%02d_%02d_%02d' % (year,month,day,hour,min,sec)
    
    write_garrawarla_assassin2_sbatch_file(EDA2_chan_list=EDA2_chan_list,lst_hrs_list=lst_hrs_list,EDA2_obs_time_list=EDA2_obs_time_list,git_repo_dir=git_repo_dir)  
