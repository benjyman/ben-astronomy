#! /usr/bin/env python
"""

generate qsbatch script for selfcal with andre's tools


"""
import os.path

def generate_namorrodor(infile,options):

    obsid_list=[]
    for line in open(infile):
       obsid_list.append(line.strip())
  
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
       
    #if (options.tagname):
    #   tagname_string=' --tagname=%s ' % (options.tagname)
    #else:
    #   tagname_string=''

    if (options.epoch_ID):
       epoch_ID_string=' --epoch_ID=%s ' % (options.epoch_ID)
    else:
       epoch_ID_string=''
    
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
   
    q_filename_path=os.path.dirname(infile)+'/'        
    q_filename='%sq_selfcal_moon_%s.sh' % (q_filename_path,str(obs_list_number))
    
       
    imaging_file = open(q_filename,'w+')
    imaging_file.write('#!/bin/bash -l\n')
    if (machine=='magnus' or machine=='galaxy'):          
      #sbatch_file.write('#!/bin/bash -l\n')
      sbatch_file.write('#SBATCH -o selfcal-%A.out\n' )
      sbatch_file.write('##SBATCH --ntasks=1\n')
      sbatch_file.write('#SBATCH --ntasks-per-node=1\n')
      sbatch_file.write('#SBATCH --time=12:00:00\n')
      sbatch_file.write('#SBATCH -J selfcal_%s\n' % (options.epoch_ID))
      #sbatch_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
      #sbatch_file.write('#SBATCH --clusters=magnus\n')
      sbatch_file.write('#SBATCH --partition=workq\n')
      sbatch_file.write('#SBATCH --account=mwaeor\n')
      sbatch_file.write('#SBATCH --export=NONE\n')
      #sbatch_file.write('module swap PrgEnv-cray/6.0.4 PrgEnv-gnu && module swap gcc gcc/5.3.0\n')
      #sbatch_file.write('module use /group/mwa/software/modulefiles\n')
      #sbatch_file.write('module load MWA_Tools/mwa-sci\n')
      #sbatch_file.write('module load setuptools\n')
          
    for obsid_index,obsid in enumerate(obsid_list):
       if (options.track_off_moon):
          track_off_moon_list_string=",".join(track_off_moon_list[int(float(obsid_index)*3):int(float(obsid_index)*3+3)])
       imaging_file.write('/data/code/git/ben-astronomy/moon/processing_scripts/namorrodor_magnus/selfcal_concat_ms.py '+str(obsid) + ' ' + track_off_moon_list_string + track_off_moon_string + track_moon_string+epoch_ID_string +ionpeel_string+chgcentre_string+minw_string+cotter_string+model_string+selfcal_string+applyonly_string+sourcelist_string+' \n')

    imaging_file.close()
    print "wrote file %s" % q_filename
    
    command="chmod +x %s " % (q_filename) 
    print command
    os.system(command)
    
###

import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_selfcal_concat_ms.py [text file of obsIDs] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
#parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files  e.g. --tagname="" [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations  e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--ionpeel',type='string', dest='ionpeel',default='',help='Use this to specify a clustered model and run ionpeel e.g.ionpeel="clustered_model.txt" [default=%default]')
parser.add_option('--chgcentre',action='store_true',dest='chgcentre',default=False,help='Calibrate the ms that has had a phase centre shift  e.g. --chgcentre   [default=%default]')
parser.add_option('--minw',action='store_true',dest='minw',default=False,help='Calibrate the ms that has had a phase centre shift to minw e.g. --minw   [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=False,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')
parser.add_option('--model',type='string', dest='model',default='',help='Specify a model to calibrate on overrides --sourcelist. e.g. model="clustered_model.txt" [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default='0',help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--applyonly',type='string', dest='applyonly',default=None,help='Just apply the specified calibration solutions to the measurements set e.g. applyonly=new_solutions.bin [default=%default]')
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (specify this instead of model and an ao model will be made from the sourclist), overidden by --model=')
parser.add_option('--track_off_moon',type='string',dest='track_off_moon',help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='namorrodor',help=' e.g. --machine="magnus" [default=%default]')


(options, args) = parser.parse_args()

infile = args[0]



mwa_dir = os.getenv('MWA_DIR','/data/MWA/')

if(mwa_dir == '/data/MWA/'):
    generate_namorrodor(infile,options)
