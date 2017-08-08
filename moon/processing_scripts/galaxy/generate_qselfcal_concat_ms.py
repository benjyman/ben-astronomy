#! /usr/bin/env python
"""

generate qsbatch script for selfcal with andre's tools


"""
import os.path

def generate_galaxy(infile,options):

    obsid_list=[]
    for line in open(infile):
       obsid_list.append(line.strip())
    obsid_list_string=",".join(obsid_list)
  
    n_obs = sum(1 for line in open(infile))

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
    
    if (options.track_moon):
       track_moon_string=' --track_moon '
    else:
       track_moon_string=''

    if (options.chgcentre):
       chgcentre_string=' --chgcentre '
    else:
       chgcentre_string=''
       
    if (options.minw):
       minw_string=' --minw '
    else:
       minw_string=''
       
    if (options.tagname):
       tagname_string=' --tagname=%s ' % (options.tagname)
    else:
       tagname_string=''

    if (options.ionpeel):
       ionpeel_string=' --ionpeel=%s ' % (options.ionpeel)
    else:
       ionpeel_string=''

    if (options.model):
       model_string=' --model=%s ' % (options.model)
       print 'using model %s' % options.model
    else:
       model_string=''

    if (options.sourcelist):
       sourcelist_string=' --sourcelist=%s ' % (options.sourcelist)
    else:
       sourcelist_string=''
    if (options.cotter):
       cotter_string=' --cotter '
    else:
       cotter_string=''
    
    if (options.selfcal):
       selfcal_string=' --selfcal=%s ' % options.selfcal
    else:
       selfcal_string=''

    if (options.applyonly):
       applyonly_string=' --applyonly=%s ' % options.applyonly
    else:
       applyonly_string=''
       
    imaging_file = open('q_selfcal_concat_ms_wrapper.sh','w+')
    imaging_file.write('#!/bin/bash -l\n')
    imaging_file.write('#SBATCH --nodes=1\n')
    imaging_file.write('#SBATCH --ntasks-per-node=1\n')
    imaging_file.write('#SBATCH --time=12:00:00\n' )
    imaging_file.write('#SBATCH --partition=workq\n')
    imaging_file.write('#SBATCH --account=mwaeor\n')
    imaging_file.write('#SBATCH --export=NONE\n')
    imaging_file.write('#SBATCH -J selfcalconcat_ms\n')
    imaging_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
    imaging_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    imaging_file.write('cd $SLURM_SUBMIT_DIR\n')


    imaging_file.write('/group/mwaeor/bmckinley/git_repos/ben-astronomy/moon/processing_scripts/galaxy/selfcal_concat_ms.py '+obsid_list_string+' ${SLURM_ARRAY_TASK_ID} ' + track_off_moon_list_string + track_off_moon_string + track_moon_string+tagname_string +ionpeel_string+chgcentre_string+minw_string+cotter_string+model_string+selfcal_string+applyonly_string+sourcelist_string+' \n')

    imaging_file.close()
    print "wrote file q_selfcal_concat_ms_wrapper.sh"
###

import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_selfcal_concat_ms.py [text file of obsIDs] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files  e.g. --tagname="" [default=%default]')
parser.add_option('--ionpeel',type='string', dest='ionpeel',default='',help='Use this to specify a clustered model and run ionpeel e.g.ionpeel="clustered_model.txt" [default=%default]')
parser.add_option('--chgcentre',action='store_true',dest='chgcentre',default=False,help='Calibrate the ms that has had a phase centre shift  e.g. --chgcentre   [default=%default]')
parser.add_option('--minw',action='store_true',dest='minw',default=False,help='Calibrate the ms that has had a phase centre shift to minw e.g. --minw   [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=False,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')
parser.add_option('--model',type='string', dest='model',default='',help='Specify a model to calibrate on overrides --sourcelist. e.g. model="clustered_model.txt" [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default='0',help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--applyonly',type='string', dest='applyonly',default=None,help='Just apply the specified calibration solutions to the measurements set e.g. applyonly=new_solutions.bin [default=%default]')
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (specify this instead of model and an ao model will be made from the sourclist), overidden by --model=')
parser.add_option('--track_off_moon',type='string',dest='track_off_moon',help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')

(options, args) = parser.parse_args()

infile = args[0]



mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/scratch2/mwaeor/MWA/'):
    generate_galaxy(infile,options)
