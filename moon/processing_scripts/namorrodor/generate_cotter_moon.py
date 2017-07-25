#python script to generate the submission script to cotter a bunch of moon obs


def generate_cotter_moon(obsid_infile,options):

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
 
    for obs_list_number in range(n_obs_lists):
 
       start_obs_id_index=int(obs_list_number*chunk_size)
    
       if (obs_list_number<n_obs_lists-1):
          end_obs_id_index=int(obs_list_number*chunk_size+chunk_size)
       else:
          end_obs_id_index=int(obs_list_number*chunk_size+last_list_length)
          
       q_filename='q_cotter_moon_%s.sh' % (str(obs_list_number))
    
       sbatch_file = open(q_filename,'w+')
       sbatch_file.write('#!/bin/bash -l\n')
    
       for obsid_index,obsid in enumerate(obsid_list[start_obs_id_index:end_obs_id_index]):
          track_off_moon_list_string=",".join(track_off_moon_list[int(float(obsid_index)*3):int(float(obsid_index)*3+3)])
          sbatch_file.write('python /data/code/git/ben-astronomy/moon/processing_scripts/namorrodor/cotter_moon.py '+ str(obsid) + ' ' + track_off_moon_list_string + '  ' + track_moon_string + ' ' + track_off_moon_string + ' '+ tagname_string + flag_ants_string + cleanup_string + ' \n')

       sbatch_file.close()
    
    command="chmod +x %s " % (q_filename) 
    print command
    os.system(command)
    
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
