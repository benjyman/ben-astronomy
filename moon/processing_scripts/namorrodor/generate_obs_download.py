#script to download a list of obsids
import os.path
import os
import numpy as np
from distutils import command

def generate_namorrodor(infile,options):

   chunk_size=20
   obsid_list=[]
   for line in open(infile):
      obsid_list.append(line.strip())
   n_obs = sum(1 for line in open(infile))

   #split into chunks of 20 obsids
   if (n_obs > chunk_size):
      n_obs_lists=int(np.ceil(n_obs/float(chunk_size)))
      last_list_length=np.remainder(n_obs,float(chunk_size))
   else:
      n_obs_lists=1
      last_list_length=len(obsid_list) 

   for obs_list_number in range(n_obs_lists):

      imaging_file = open('q_obsdownload_wrapper_chunk_%s.sh' % (str(obs_list_number)),'w+')
      imaging_file.write('cd /data/MWA\n')
      
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

      command="chmod +x q_obsdownload_wrapper_chunk_%s.sh " % (str(obs_list_number)) 
      print command
      os.system(command)
      
      
import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = "/data/MWA/"

usage = 'Usage: generate_obs_download.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()

infile = args[0]
  
generate_namorrodor(infile,options)







