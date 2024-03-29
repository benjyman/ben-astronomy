#script to download a list of obsids
import os.path
import os
import numpy as np
from distutils import command

def generate_download_obs(options): 
   machine=options.machine.strip()
   database=options.database
   infile=options.obsid_infile
   if machine=="namorrodor":
      mwa_dir = "/md0/moon/data/MWA/"
      code_dir="/data/code/git/MWA_Tools/scripts/"
   elif  machine=="magnus":
      mwa_dir= "/astro/mwaeor/MWA/data/"
      code_dir=""
   else:
      print "machine does not match any known"   

   chunk_size=20
   obsid_list=[]
   for line in open(infile):
      obsid_list.append(line.strip())
   n_obs = sum(1 for line in open(infile))

   if not options.gator:
      #split into chunks of 20 obsids
      if (n_obs > chunk_size):
         n_obs_lists=int(np.ceil(n_obs/float(chunk_size)))
         last_list_length=np.remainder(n_obs,float(chunk_size))
      else:
         n_obs_lists=1
         last_list_length=len(obsid_list) 
      
      for obs_list_number in range(n_obs_lists):
      
         imaging_file = open('q_obsdownload_wrapper_chunk_%s.sh' % (str(obs_list_number)),'w+')
         
         start_obs_id_index=int(obs_list_number*chunk_size)
      
         if (obs_list_number<n_obs_lists-1):
             end_obs_id_index=int(obs_list_number*chunk_size+chunk_size)
         else:
             end_obs_id_index=int(obs_list_number*chunk_size+last_list_length)
      
         if machine=="magnus":
            imaging_file.write('#!/bin/bash -l\n')
            imaging_file.write('#SBATCH -o getNGASdata-%A.out\n' )
            imaging_file.write('##SBATCH --ntasks=1\n')
            imaging_file.write('#SBATCH --ntasks-per-node=1\n')
            imaging_file.write('#SBATCH --time=12:00:00\n')
            imaging_file.write('#SBATCH -J download_obs\n')
            #imaging_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
            #imaging_file.write('#SBATCH --clusters=magnus\n')
            imaging_file.write('#SBATCH --partition=workq\n')
            imaging_file.write('#SBATCH --account=mwaeor\n')
            imaging_file.write('#SBATCH --export=NONE\n')
            #imaging_file.write('module swap PrgEnv-cray/6.0.4 PrgEnv-gnu && module swap gcc gcc/5.3.0\n')
            #imaging_file.write('module use /group/mwa/software/modulefiles\n')
            #imaging_file.write('module load MWA_Tools/mwa-sci\n')
            #imaging_file.write('module load setuptools\n')
            imaging_file.write('cd %s\n' % mwa_dir)
      
         elif machine=="namorrodor":
            imaging_file.write('cd %s\n' % mwa_dir)
      
         for obsid in obsid_list[start_obs_id_index:end_obs_id_index]:
            imaging_file.write('%sobsdownload.py -o %s\n' % (code_dir,obsid))
            imaging_file.write('%sobsdownload.py -f -o %s\n' % (code_dir,obsid))
            #imaging_file.write('%sobsdownload.py -f -m -o %s\n' % (code_dir,obsid))
            #imaging_file.write('%sobsdownload.py -f -o %s\n' % (code_dir,obsid))
            #imaging_file.write('cd %s\n' % obsid)
            #imaging_file.write('make_metafits.py -g %s\n' % (obsid))
            #imaging_file.write('cd ..\n')
      
         imaging_file.close()
         print "wrote file q_obsdownload_wrapper_chunk_%s.sh " % (str(obs_list_number)) 
      
         command="chmod +x q_obsdownload_wrapper_chunk_%s.sh " % (str(obs_list_number)) 
         print command
         os.system(command)
   else:
      download_filename_path=os.path.dirname(infile)+'/'
      download_filename='%sq_obsdownload_wrapper.sh' % download_filename_path
      #database='/group/mwaeor/bmckinley/%s_05_93_on_moon.sqlite' % (epoch_)
      #use gator
      with open(download_filename,'w') as f:
         f.write('#!/bin/bash -l\n')
         f.write('gator_add_to_rts_table.rb -d %s %s\n' % (database,infile))
         ##f.write('nohup gator_download_daemon.rb -d %s &\n' % (database))
      
      
      
import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_obs_download.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine being used for processing, current options are namorrodor and magnus [default=%default]')
parser.add_option('--gator',action='store_true',dest='gator',default=True,help='Use gator for downloading on galaxy and magnus [default=%default]')
parser.add_option('--database',type='string', dest='database',default='/group/mwaeor/bmckinley/mwa_downloads.sqlite',help='name of database to be used for downloads if using gator [default=%default]')
parser.add_option('--obsid_infile',type='string', dest='obsid_infile',default='',help='File containing list of obsids to be cottered  e.g. --obsid_infile="20180107_moon_93.txt" [default=%default]')

(options, args) = parser.parse_args()
  
generate_download_obs(options)







