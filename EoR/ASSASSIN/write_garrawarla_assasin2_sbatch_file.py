#!/usr/bin/env python3

import numpy as np
import math

fine_chans_per_EDA2_chan = 27
fine_chan_width_Hz =  (400./512.) / 27. * 1.0e6 #(see Marcin email)
fine_chan_width_MHz = fine_chan_width_Hz / 1.0e6

def calculate_freq_MHz_fine_chan_subarray(EDA2_chan):
   freq_MHz_centre = (400./512.)*EDA2_chan
   fine_chan_index_array = np.arange(fine_chans_per_EDA2_chan) - 14
   freq_MHz_fine_chan_sub_array = fine_chan_index_array*fine_chan_width_MHz + freq_MHz_centre + (fine_chan_width_MHz/2.)
   return freq_MHz_fine_chan_sub_array
   
def write_garrawarla_assassin2_sbatch_file(EDA2_chan_list,lst_hrs_list,EDA2_obs_time_list,git_repo_dir="/astro/mwaeor/bmckinley/code/"):
   base_freq_MHz = EDA2_chan_list[0] * (400./512.)
   for EDA2_chan_index,EDA2_chan in enumerate(EDA2_chan_list):
      
      freq_MHz_fine_chan_subarray = calculate_freq_MHz_fine_chan_subarray(EDA2_chan)
      
      if (EDA2_chan_index==0 and math.isclose(base_freq_MHz,50.0)):
         fine_chans_per_EDA2_chan = 13
         freq_MHz_fine_chan_subarray = freq_MHz_fine_chan_subarray[14:]
      else:
         fine_chans_per_EDA2_chan = 27
      
      EDA2_obs_time = EDA2_obs_time_list[EDA2_chan_index]
      lst_hrs = float(lst_hrs_list[EDA2_chan_index])
      freq_MHz_fine_chan_subarray_string = ','.join([str(x) for x in freq_MHz_fine_chan_subarray])
      
      #year,month,day,hour,min,sec = time_string.split('_')
      #time_formatted = '%d-%02d-%02dT%02d:%02d:%02d' % (float(year),float(month),float(day),float(hour),float(min),float(sec))

      name_base = "assassin2_sbatch_eda2_chan_%03d_%s" % (EDA2_chan,EDA2_obs_time)  
      sbatch_filename = "%s.sh" % name_base
   
      with open('%s' % sbatch_filename,'w') as outfile:
         outfile.write("#!/bin/bash --login\n#SBATCH --nodes=1\n#SBATCH --partition=gpuq\n#SBATCH --gres=gpu:1\n")
         outfile.write("#SBATCH --time=00:30:00\n#SBATCH --account=mwaeor\n#SBATCH --nodes=1\n#SBATCH --mem=10gb\n")
         outfile.write("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --array=0-%s\n\n" % str(fine_chans_per_EDA2_chan-1))
         
         
         outfile.write("time python %sgarrawarla_sim_with_complex_beams.py --EDA2_chan=%s --EDA2_chan_index=%s --lst_hrs=%s --EDA2_obs_time=%s --freq_MHz_fine_chan_subarray=%s\n" % (git_repo_dir,EDA2_chan,EDA2_chan_index,lst_hrs,EDA2_obs_time,freq_MHz_fine_chan_subarray_string))
         #for WODEN
         #outfile.write("module use /pawsey/mwa/software/python3/modulefiles\nmodule load erfa/1.7.0\n")
         #outfile.write("module load json-c/0.14\nmodule load hdf5/1.10.5\nmodule load cfitsio/3.48\nmodule load cmake/3.15.0\n")
         #outfile.write("module load cuda/10.2\nmodule load pal/0.9.8\nmodule load python/3.8.2\nmodule load astropy/4.0.1.post1\n\n")
         ##for woden_sourcelists.py
         #outfile.write("module load scipy/1.4.1\nmodule load h5py/2.10.0\nmodule load healpy/1.13.0\nmodule load numpy/1.18.2\n\n")
         #outfile.write("module load ephem/3.7.7.1\nmodule load matplotlib/3.2.1\n\n")
   
         #outfile.write("source /astro/mwaeor/jline/software/WODEN/build/init_WODEN.sh\n")
         #outfile.write("export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH\n\n")

                     
      print("wrote %s" % sbatch_filename) 
 
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

    parser.add_argument('--EDA2_chan_list', default='64,65',
        help='list of EDA2 coarse chans. Generates a scripts for each that runs an array job for each fine chan. e.g. --EDA2_chan_list="64,65"')
        
    parser.add_argument('--EDA2_obs_time_list', default='20200303T000000,20200303T000400',
        help='list of EDA2 obs times that correspond to each EDA2 coarse chan. \
             must be same length as EDA2_chan_list and for data calibration must correspond to ms filenames e.g. --EDA2_obs_time_list="20200303T000000,20200303T000400"')

    parser.add_argument('--lst_hrs_list', default='6.0,6.2',
        help='list of LSTs that corresponds to the EDA2 obs time list. e.g. --lst_hrs_list="6.0,6.2"')

    parser.add_argument('--git_repo_dir', default='/astro/mwaeor/bmckinley/code/',
        help='directory where my code is kept --git_repo_dir="/astro/mwaeor/bmckinley/code/"')

        
    args = parser.parse_args()
    
    if args.EDA2_chan_list:
       EDA2_chan_list_strings = [s.strip() for s in args.EDA2_chan_list.split(",")]
       EDA2_chan_list = [float(x) for x in EDA2_chan_list_strings] 
       
    if args.EDA2_obs_time_list:
       EDA2_obs_time_list = [s.strip() for s in args.EDA2_obs_time_list.split(",")]
       
    if args.lst_hrs_list:
       lst_hrs_list_strings = [s.strip() for s in args.lst_hrs_list.split(",")]
       lst_hrs_list = [float(x) for x in lst_hrs_list_strings] 
       
    if args.git_repo_dir:
       git_repo_dir = args.git_repo_dir
    ##this is for LST 60.0 deg, will do rest of lsts later
    #year,month,day,hour,min,sec = 2015,11,29,15,40,29 #LST=60 deg
    #time_string = '%d_%02d_%02d_%02d_%02d_%02d' % (year,month,day,hour,min,sec)
    
    write_garrawarla_assassin2_sbatch_file(EDA2_chan_list=EDA2_chan_list,lst_hrs_list=lst_hrs_list,EDA2_obs_time_list=EDA2_obs_time_list,git_repo_dir=git_repo_dir)  
