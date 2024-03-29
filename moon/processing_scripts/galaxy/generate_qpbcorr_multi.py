#! /usr/bin/env python
"""
pbcorr a bunch of images at different chans
"""
import os.path

def generate_galaxy(infile,options):

    obsid_list=[]
    for line in open(infile):
       obsid_list.append(line.strip())
    obsid_list_string=",".join(obsid_list)
  
    n_obs = sum(1 for line in open(infile))
    n_images=int(options.channelsout)
    n_array_jobs=int(n_obs*n_images)

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

    if (options.dirty):
       dirty_string=' --dirty '
    else:
       dirty_string=''

    if (options.tagname):
       tagname_string=' --tagname=%s ' % (options.tagname)
    else:
       tagname_string=''

    if (options.image_base_name):
       image_base_name=options.image_base_name
       image_base_name_string=' --image_base_name=%s ' % (image_base_name)
    else:
       image_base_name_string=''

    if (options.havebeam):
       havebeam_string=' --havebeam=%s '  % options.havebeam
    else:
       havebeam_string=''

    if (options.ms_base_name):
       ms_name=options.ms_base_name
       ms_name_string=' --ms_base_name=%s ' % (ms_name)
    else:
       ms_name_string=''
       
    if (options.metafits_name):
       metafits_name=options.metafits_name
       metafits_name_string=' --metafits_name=%s ' % (metafits_name)
    else:
       metafits_name_string=''
   
    if (options.applyion):
       applyion_string=' --applyion=%s ' % (options.applyion)
       clustered_model_string=' --clustered_model=%s ' % (options.clustered_model)
       if (options.render):
          render_string=' --render '
       else:
          render_string=''
       if (options.pbuncorrect):
          pbuncorrect_string=' --pbuncorrect '
       else:
          pbuncorrect_string=''
    else:
       applyion_string=''    
       clustered_model_string=''
       render_string=''
       pbuncorrect_string=''
    

    if (options.selfcal):
       selfcal_string=' --selfcal=%s ' % options.selfcal
    else:
       selfcal_string=''    

    if (options.ionpeeled):
       ionpeeled_string = ' --ionpeeled'
    else:
       ionpeeled_string = ' '

    if (options.chgcentre):
       chgcentre_string=' --chgcentre="%s" ' % (options.chgcentre)
    else:
       chgcentre_string=' '

    if (options.minw):
       minw_string=' --minw="%s" ' % (options.minw)
    else:
       minw_string=' '

    if (options.channelsout):
       channelsout_string=' --channelsout="%s" ' % (options.channelsout)
    else:
       channelsout_string=' --channelsout="1" '
 
    imaging_file = open('q_pbcorr_multi_wrapper.sh','w+')
    imaging_file.write('#!/bin/bash -l\n')
    imaging_file.write('#SBATCH --time=04:00:00\n')
    imaging_file.write('#SBATCH --partition=workq\n')
    imaging_file.write('#SBATCH --account=mwaeor\n')
    imaging_file.write('#SBATCH --export=NONE\n')
    imaging_file.write('#SBATCH -J pbcorr_multi\n')
    imaging_file.write('#SBATCH --array=0-%s\n' % (n_array_jobs-1))
    imaging_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    imaging_file.write('cd $SLURM_SUBMIT_DIR\n')
    
    imaging_file.write('python /group/mwaeor/CODE/MWA_Tools/scripts/pbcorr_multi.py ' +obsid_list_string+' ${SLURM_ARRAY_TASK_ID} '+track_off_moon_list_string+image_base_name_string +havebeam_string+ms_name_string+metafits_name_string+applyion_string+clustered_model_string+render_string+pbuncorrect_string+selfcal_string+track_moon_string+track_off_moon_string+tagname_string+ionpeeled_string+chgcentre_string+minw_string+dirty_string+channelsout_string+' \n')

    imaging_file.close()
    print "wrote file q_pbcorr_multi_wrapper.sh" 

###

import sys,os,glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_qpbcorr_multi.py [list of obsids] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files  e.g. --tagname="" [default=%default]')
parser.add_option('--track_off_moon',type='string',dest='track_off_moon',help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--image_base_name',type='string', dest='image_base_name',default='',help='image_base_name to be pb corrected  e.g. --image_base_name="" [default=%default]')
parser.add_option('--havebeam', type='string', dest='havebeam', default=None,help='Set this if you already have downloaded the right beam model and give the beam prefix  e.g. --havebeam=tagname (without the _beam at the end) [default=%default]')
parser.add_option('--ms_base_name',type='string', dest='ms_base_name',default='',help='name of ms to be used to generate beam (not needed if havebeam=True  e.g. --ms_name="blah.ms" [default=%default]')
parser.add_option('--metafits_name',type='string', dest='metafits_name',default='',help='name of metafits file to be used to generate beam (not needed if havebeam=True, should be same obsid as ms_name  e.g. --metafits_name="blah.metafits" [default=%default]')
parser.add_option('--applyion',type='string', dest='applyion',default=None,help='apply the ionpeel corrections to the pbcorrected, specify the sourcelist from which the ionpeel solutions were created as in selfcal_concat.py i.e. applyion="srclist_IDR2deep_puma-noextended_HydAGuas.txt" [default=%default]')
parser.add_option('--clustered_model',type='string', dest='clustered_model',default=None,help='If applying ionsolutions, or rendering, specify the clustered model e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--render',action='store_true', dest='render',default=False,help='Render the clustered model to the ion-applied image e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--pbuncorrect',action='store_true', dest='pbuncorrect',default=False,help='pbuncorrect the rendered, ion-applied image so you can pbaddimage them together! e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default='0',help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--ionpeeled',dest='ionpeeled',action='store_true',default=False,help='Set if this is a subsequent imaging run after ionpeeling e.g. --ionpeeled [default=%default]')
parser.add_option('--chgcentre',type='string',dest='chgcentre',default=None,help='Use chgcentre to change phase centre to desired phase centre or shift to min w and shiftback e.g. --chgcentre=" --chgcentre=" 00h00m00.0s 00d00m00.0s "  [default=%default]')
parser.add_option('--minw',type='string',dest='minw',default=None,help='Shift to minw position of whatever ms is central in the chunk and then shiftback (must start at same phase centre which eor obs do e.g. --minw="12345678.ms" or minw="self" will use the minw of the observation  [default=%default]')
parser.add_option('--dirty',action='store_true',dest='dirty',default=False,help='perform all operations on the dirty images, not the cleaned images [default=%default]')
parser.add_option('--channelsout',type='string', dest='channelsout',default='1',help='Specify how many channels were output from wsclean e.g.channelsout=24 [default=%default]')

(options, args) = parser.parse_args()

infile = args[0]
generate_galaxy(infile,options)
