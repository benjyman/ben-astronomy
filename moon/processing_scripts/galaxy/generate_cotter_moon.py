#python script to generate the submission script to cotter a bunch of moon obs


def generate_cotter_moon(obsid_infile,options):

    obsid_list=[]
    for line in open(obsid_infile):
       obsid_list.append(line.strip())
    obsid_list_string=",".join(obsid_list)  

    if (options.track_off_moon):
       track_off_moon_list=[]
       track_off_moon_infile=options.track_off_moon
       for line in open(track_off_moon_infile):
          centre_obs_id=line.split()[0]
          new_RA=line.split()[1]
          new_DEC=line.split()[2]
          track_off_moon_list.append(centre_obs_id)
          track_off_moon_list.append(new_RA)
          track_off_moon_list.append(new_DEC)
       track_off_moon_list_string=",".join(track_off_moon_list)
       
       track_off_moon_string=' --track_off_moon '
    else:
       track_off_moon_string=''
       track_off_moon_list_string=''

    n_obs = sum(1 for line in open(obsid_infile))
    
    if (options.track_moon):
       track_moon_string=' --track_moon '
    else:
       track_moon_string=''

    if (options.tagname):
       tagname_string=' --tagname=%s ' % (options.tagname)
    else:
       tagname_string=''

    if (options.cleanup):
       cleanup_string=' --cleanup ' 
    else:
       cleanup_string=''
       
    if (options.flag_ants):
       flag_ants_string=' --flag_ants="%s" ' % (options.flag_ants)
    else:
       flag_ants_string=''   
 
    q_filename='q_cotter_moon.sh'
    
    sbatch_file = open(q_filename,'w+')
    sbatch_file.write('#!/bin/bash -l\n')
    sbatch_file.write('#SBATCH --nodes=1\n')
    sbatch_file.write('#SBATCH --ntasks-per-node=1\n')
    sbatch_file.write('#SBATCH --time=02:00:00\n')
    sbatch_file.write('#SBATCH --partition=workq\n')
    sbatch_file.write('#SBATCH --account=mwaeor\n')
    sbatch_file.write('#SBATCH --export=NONE\n')
    sbatch_file.write('#SBATCH -J cotter_moon\n')
    sbatch_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
    #sbatch_file.write('#SBATCH --mem=32000\n')
    sbatch_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    sbatch_file.write('cd $SLURM_SUBMIT_DIR\n')

    sbatch_file.write('python /group/mwaeor/bmckinley/git_repos/ben-astronomy/moon/processing_scripts/galaxy/cotter_moon.py '+obsid_list_string +' ${SLURM_ARRAY_TASK_ID} '+ track_off_moon_list_string + '  ' + track_moon_string + ' ' + track_off_moon_string + ' '+ tagname_string + flag_ants_string + cleanup_string + ' \n')

    sbatch_file.close()
   
    print "Wrote file %s" % q_filename
    
import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_cotter_moon.py [text file of obsIDs] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--track_off_moon',type='string',dest='track_off_moon',help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files  e.g. --tagname="" [default=%default]')
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='List of antennas (space separated) to flag after cottering (andre indexing as for rfigui etc)  e.g. --flag_ants="56,60" [default=%default]')
parser.add_option('--cleanup',action='store_true',dest='cleanup',default=False,help='Delete the gpubox files after making the ms [default=%default]')


(options, args) = parser.parse_args()

obsid_infile = args[0]

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

generate_cotter_moon(obsid_infile,options)
