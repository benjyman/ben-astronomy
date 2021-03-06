#! /usr/bin/env python
"""
generate_mwac_qimage_auto.py
Generates a qsub scripts to qa, and image a list of obsIDs that have already been processed by the RTS.
04/09/13: Writes a two part 

"""
import os.path

def generate_namorrodor(options):
    infile = options.obsid_infile
    machine=options.machine
    if machine=='namorrodor':
       ben_code_base='/data/code/git/'
    else:
       ben_code_base='/astro/mwaeor/bmckinley/code/'

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
    
    if options.no_pbcorr:
       no_pbcorr_string = ' --no_pbcorr '
    else:
       no_pbcorr_string = ' '

    if (options.track_moon):
       track_moon_string=' --track_moon '
    else:
       track_moon_string=''

    if (options.chgcentre):
       chgcentre_string=' --chgcentre="%s" ' % (options.chgcentre)
    else:
       chgcentre_string=' '
       
    if (options.minw):
       minw_string=' --minw="%s" ' % (options.minw)
    else:
       minw_string=' '
       
    if (options.corrected_data):
       corrected_data_string=' --corrected_data '
    else:
       corrected_data_string=''

    if (options.concat6):
       concat6_string=' --concat6 '
    else:
       concat6_string=''

    if (options.cotter):
       cotter_string=' --cotter '
    else:
       cotter_string=''

    if (options.crop_images):
       crop_images_string=' --crop_images '
    else:
       crop_images_string=''
       
    #if (options.tagname):
    #   tagname=options.tagname
    #   tagname_string=' --tagname=%s ' % (tagname)
    #else:
    #   tagname_string=''

    if (options.epoch_ID):
       epoch_ID=options.epoch_ID
       epoch_ID_string=' --epoch_ID=%s ' % (epoch_ID)
    else:
       epoch_ID_string=''
      
    if (options.machine):
       machine_string=' --machine=%s ' % (options.machine)
    else:
       machine_string=''

    if (options.pol):
       pol_string = ' --pol='+options.pol
    else:
       pol_string = ' '

    if (options.selfcal):
       selfcal_string = ' --selfcal='+options.selfcal
    else:
       selfcal_string = ' '
     
    if (options.ionpeeled):
       ionpeeled_string = ' --ionpeeled'
    else:
       ionpeeled_string = ' '
    
    if (options.multi):
       if (options.cotter):
          vis_suffix=''
       else:
          vis_suffix='_concat_transform'
       #if (options.minw):
       #   vis_suffix+='_minw'
       if (options.ionpeeled):
          vis_suffix+='_peeled'

       
    imsize_string=' --imsize='+options.imsize+' '

    wsclean_options_string=' --wsclean_options="'+options.wsclean_options+'" '

    q_filename_path=os.path.dirname(infile)+'/'   
    if options.track_moon:
       q_filename='%sq_image_on_moon.sh' % (q_filename_path)     
    elif options.track_off_moon:
       q_filename='%sq_image_off_moon.sh' % (q_filename_path)
    else:
       q_filename='%sq_image_moon.sh' % (q_filename_path)

    sbatch_file = open(q_filename,'w+')
    sbatch_file.write('#!/bin/bash -l\n')
    if (machine=='magnus' or machine=='galaxy'):          
       #sbatch_file.write('#!/bin/bash -l\n')
       sbatch_file.write('#SBATCH -o image-%A.out\n' )
       sbatch_file.write('##SBATCH --ntasks=1\n')
       sbatch_file.write('#SBATCH --ntasks-per-node=1\n')
       sbatch_file.write('#SBATCH --time=12:00:00\n')
       sbatch_file.write('#SBATCH -J image_%s\n' % (options.epoch_ID))
       #sbatch_file.write('#SBATCH --array=0-%s\n' % (n_obs-1))
       #sbatch_file.write('#SBATCH --clusters=magnus\n')
       sbatch_file.write('#SBATCH --partition=workq\n')
       sbatch_file.write('#SBATCH --account=mwaeor\n')
       sbatch_file.write('#SBATCH --export=NONE\n')
    for obsid_index,obsid in enumerate(obsid_list):
       if (options.track_off_moon):
          track_off_moon_list_string=",".join(track_off_moon_list[int(float(obsid_index)*3):int(float(obsid_index)*3+3)])
       sbatch_file.write('python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/image_concat_ms.py %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' % (ben_code_base,str(obsid),track_off_moon_list_string,track_off_moon_string, no_pbcorr_string, track_moon_string,chgcentre_string,minw_string,imsize_string,epoch_ID_string, pol_string,concat6_string,wsclean_options_string,cotter_string,selfcal_string,ionpeeled_string,machine_string,crop_images_string) )
    
    sbatch_file.close()
    print "wrote %s" %  q_filename

    command="chmod +x %s " % (q_filename) 
    print command
    os.system(command)

import sys,os,glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_mwac_qimage_concat_ms.py [options]  [text file of obsIDs]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--track_off_moon',type='string',dest='track_off_moon',help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--imsize',type='string', dest='imsize',default='1024',help='Image size string in pixels for wsclean e.g. --imsize="4096" [default=%default]')
parser.add_option('--no_pbcorr',action='store_true',dest='no_pbcorr',default=False,help='Dont generate beam and do pbcorrection at this stage (choose this if doing mosaicking later on). [default=%default]')
#parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files to image  e.g. --tagname="" [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations  e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--corrected_data',action='store_true',dest='corrected_data',default=False,help='If doing self-cal, use the CORRECTED_DATA column instead of the DATA column. [default=%default]')
parser.add_option('--pol',dest='pol',default='xx,yy',help='Polarisations to image with wsclean e.g. pol="xx,yy,xy,yx" [default=%default]')
parser.add_option('--concat6',action='store_true',dest='concat6',default=False,help='Image ms with 6 coarse chans concatenated together (7.68 MHz) [default=%default]')
parser.add_option('--chgcentre',type='string',dest='chgcentre',default=None,help='Use chgcentre to change phase centre to desired phase centre or shift to min w and shiftback e.g. --chgcentre=" --chgcentre=" 00h00m00.0s 00d00m00.0s "  [default=%default]')
parser.add_option('--wsclean_options',dest='wsclean_options',default=' -niter 100000 -mgain 0.80 -scale 0.0085 -weight uniform -threshold 0.05  -smallinversion ',help=' [default=%default]')
parser.add_option('--multi',action='store_true',dest='multi',default=False,help='Image all the obsids in the obsid list at the same time (better SNR and UV coverage) [default=%default]')
parser.add_option('--selfcal',dest='selfcal',type='string',default=None,help='Set if this is a subsequent imaging run after selfcal e.g. --selfcal=1 [default=%default]')
parser.add_option('--ionpeeled',dest='ionpeeled',action='store_true',default=False,help='Set if this is a subsequent imaging run after ionpeeling e.g. --ionpeeled [default=%default]')
parser.add_option('--minw',type='string',dest='minw',default=None,help='Shift to minw position of whatever ms is central in the chunk and then shiftback (must start at same phase centre which eor obs do e.g. --minw="12345678.ms" or minw="self" will use the minw of the observation  [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=False,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help=' e.g. --machine="magnus" [default=%default]')
parser.add_option('--obsid_infile',type='string', dest='obsid_infile',default='',help='File containing list of obsids to be cottered  e.g. --obsid_infile="20180107_moon_93.txt" [default=%default]')
parser.add_option('--crop_images',action='store_true',dest='crop_images',default=False,help='Crop the images to the size required by model_moon.py e.g. --crop_images   [default=%default]')


(options, args) = parser.parse_args()


generate_namorrodor(options)
