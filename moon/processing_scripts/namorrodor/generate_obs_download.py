#script to download a list of obsids
import os.path
import numpy as np

def generate_namorrodor(infile,options):


   obsid_list=[]
   for line in open(infile):
      obsid_list.append(line.strip())
   n_obs = sum(1 for line in open(infile))

   #split into chunks of 10 obsids
   if (n_obs > chunk_size):
      n_obs_lists=int(np.ceil(n_obs/float(chunk_size)))
      last_list_length=np.remainder(n_obs,float(chunk_size))
   else:
      n_obs_lists=1
      last_list_length=len(obsid_list) 

   for obs_list_number in range(n_obs_lists):

      imaging_file = open('q_obsdownload_wrapper_chunk_%s.sh' % (str(obs_list_number)),'w+')
      imaging_file.write('#!/bin/bash -l\n')
      imaging_file.write('#SBATCH -o getNGASdata-%A.out\n' )
      imaging_file.write('##SBATCH --ntasks=1\n')
      imaging_file.write('#SBATCH --ntasks-per-node=1\n')
      imaging_file.write('#SBATCH --time=12:00:00\n')
      imaging_file.write('#SBATCH -J download_obs\n')
      #imaging_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
      imaging_file.write('#SBATCH --clusters=zeus\n')
      imaging_file.write('#SBATCH --partition=copyq\n')
      imaging_file.write('#SBATCH --account=mwaeor\n')
      imaging_file.write('#SBATCH --export=NONE\n')
      imaging_file.write('module load pyephem\n')
      imaging_file.write('module load setuptools\n')
      imaging_file.write('cd /scratch2/mwaeor/MWA/data\n')
      
      start_obs_id_index=int(obs_list_number*chunk_size)

      if (obs_list_number<n_obs_lists-1):
          end_obs_id_index=int(obs_list_number*chunk_size+chunk_size)
      else:
          end_obs_id_index=int(obs_list_number*chunk_size+last_list_length)

      for obsid in obsid_list[start_obs_id_index:end_obs_id_index]:
         imaging_file.write('obsdownload.py -o %s\n' % obsid)
         imaging_file.write('obsdownload.py -f -o %s\n' % obsid)
         imaging_file.write('cd %s\n' % obsid)
         imaging_file.write('make_metafits.py -g %s\n' % obsid)
         imaging_file.write('cd ..\n')

      imaging_file.close()
      print "wrote file q_obsdownload_wrapper_chunk_%s.sh " % (str(obs_list_number)) 


import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

usage = 'Usage: generate_obs_download.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()

infile = args[0]
  
generate_galaxy(infile,options)







