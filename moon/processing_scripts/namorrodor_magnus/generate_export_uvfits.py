#python script to generate the submission script to cotter a bunch of moon obs


def generate_export_uvfits(options):
    machine=options.machine
    if machine=='namorrodor':
       ben_code_base='/data/code/git/'
    else:
       ben_code_base='/astro/mwaeor/bmckinley/code/'
      
    obsid_infile=options.obsid_infile
   
    chunk_size=20
    obsid_list=[]
    for line in open(obsid_infile):
       obsid_list.append(line.strip()) 
    n_obs = sum(1 for line in open(obsid_infile))
    
    
    #split into chunks of 20 obsids
    if (n_obs > chunk_size):
       n_obs_lists=int(np.ceil(n_obs/float(chunk_size)))
       last_list_length=np.remainder(n_obs,float(chunk_size))
    else:
       n_obs_lists=1
       last_list_length=len(obsid_list) 

    n_obs = sum(1 for line in open(obsid_infile))
    
    #if (options.tagname):
    #   tagname_string=' --tagname=%s ' % (options.tagname)
    #else:
    #   tagname_string=''

    if (options.epoch_ID):
       epoch_ID_string=' --epoch_ID=%s ' % (options.epoch_ID)
    else:
       epoch_ID_string=''

    if (options.machine):
       machine_string=' --machine=%s ' % (options.machine)
    else:
       machine_string=''
 
    for obs_list_number in range(n_obs_lists):
 
       start_obs_id_index=int(obs_list_number*chunk_size)
    
       if (obs_list_number<n_obs_lists-1):
          end_obs_id_index=int(obs_list_number*chunk_size+chunk_size)
       else:
          end_obs_id_index=int(obs_list_number*chunk_size+last_list_length)
       
       q_filename_path=os.path.dirname(obsid_infile)+'/'        
       q_filename='%sq_export_uvfits_%s.sh' % (q_filename_path,str(obs_list_number))
    
       sbatch_file = open(q_filename,'w+')
       sbatch_file.write('#!/bin/bash -l\n')
       if (machine=='magnus' or machine=='galaxy'):          
          #sbatch_file.write('#!/bin/bash -l\n')
          sbatch_file.write('#SBATCH -o cotter-%A.out\n')
          sbatch_file.write('##SBATCH --ntasks=1\n')
          sbatch_file.write('#SBATCH --ntasks-per-node=1\n')
          sbatch_file.write('#SBATCH --time=12:00:00\n')
          sbatch_file.write('#SBATCH -J cotter_%s\n' % (options.epoch_ID))
          #sbatch_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
          #sbatch_file.write('#SBATCH --clusters=magnus\n')
          sbatch_file.write('#SBATCH --partition=workq\n')
          sbatch_file.write('#SBATCH --account=mwaeor\n')
          sbatch_file.write('#SBATCH --export=NONE\n')
          #sbatch_file.write('module swap PrgEnv-cray/6.0.4 PrgEnv-gnu && module swap gcc gcc/5.3.0\n')
          #sbatch_file.write('module use /group/mwa/software/modulefiles\n')
          #sbatch_file.write('module load MWA_Tools/mwa-sci\n')
          #sbatch_file.write('module load setuptools\n')
    
       for obsid_index,obsid in enumerate(obsid_list[start_obs_id_index:end_obs_id_index]):
          obsid_string=' --obsid=%s ' % obsid
          sbatch_file.write('python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/export_uvfits.py  %s %s %s\n' % (ben_code_base,obsid_string,epoch_ID_string,machine_string) )
       sbatch_file.close()
    
    command="chmod +x %s " % (q_filename) 
    print command
    os.system(command)
    
    print "Wrote file %s" % q_filename
    
import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_export_uvfits.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--obsid_infile',type='string', dest='obsid_infile',default='',help='File containing list of obsids to be cottered  e.g. --obsid_infile="20180107_moon_93.txt" [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')

(options, args) = parser.parse_args()


generate_export_uvfits(options)
