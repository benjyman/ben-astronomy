#python script to generate the submission script to cotter a bunch of moon obs


def generate_cotter_moon(options):
    machine=options.machine
    if machine=='namorrodor':
       ben_code_base='/data/code/git/'
    else:
       ben_code_base='/astro/mwaeor/bmckinley/code/'
      
    obsid_infile=options.obsid_infile
    if (options.track_moon or options.track_off_moon):
       sister_obsid_infile=options.sister_obsid_infile

    if options.no_dysco:
       dysco_string='--no_dysco'
    else:
       dysco_string=''
    
    chunk_size=20
    obsid_list=[]
    sister_obsid_list=[]
    for line in open(obsid_infile):
       obsid_list.append(line.strip()) 
    if (options.track_moon or options.track_off_moon):
       for line in open(sister_obsid_infile):
          sister_obsid_list.append(line.strip()) 
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
       
    if (options.epoch_ID):
       epoch_semester=options.epoch_ID[0:5]
       print epoch_semester
       if (epoch_semester=="2015A"):
          time_freq_res="8,80"
       elif (epoch_semester=="2015B"):
          time_freq_res="8,80"
       elif (epoch_semester=="2014A"):
          time_freq_res="8,80"
       elif (epoch_semester=="2014B"):
          time_freq_res="8,80" 
       elif (epoch_semester=="2017B"):
          time_freq_res="4,40"         
       elif (epoch_semester=="2018A"):
          time_freq_res="4,40"
       time_freq_res_string=' --time_freq_res=%s ' % (time_freq_res)
    else:
       time_freq_res_string=''
       
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
       
       q_filename_path=os.path.dirname(obsid_infile)+'/'        
       if options.track_moon:
          q_filename='%sq_cotter_on_moon_%s.sh' % (q_filename_path,str(obs_list_number))
       elif: options.track_off_moon:
          q_filename='%sq_cotter_off_moon_%s.sh' % (q_filename_path,str(obs_list_number))
       else:
          q_filename='%sq_cotter_moon_%s.sh' % (q_filename_path,str(obs_list_number))
       sbatch_file = open(q_filename,'w+')
       sbatch_file.write('#!/bin/bash -l\n')
       if (machine=='magnus' or machine=='galaxy'):          
          #sbatch_file.write('#!/bin/bash -l\n')
          sbatch_file.write('#SBATCH -o cotter-%A.out\n')
          sbatch_file.write('##SBATCH --ntasks=1\n')
          sbatch_file.write('#SBATCH --ntasks-per-node=1\n')
          sbatch_file.write('#SBATCH --time=00:30:00\n')
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
          if (options.track_off_moon or options.track_moon):
             sister_obsid=sister_obsid_list[start_obs_id_index+obsid_index]
             sister_obsid_string=' --sister_obsid=%s ' % sister_obsid
          else:
             sister_obsid_string=''
          if (options.track_off_moon):
             track_off_moon_list_string="--track_off_moon_list="+",".join(track_off_moon_list[int(float(obsid_index)*3):int(float(obsid_index)*3+3)])
          sbatch_file.write('python %sben-astronomy/moon/processing_scripts/namorrodor_magnus/cotter_moon.py %s %s %s %s %s %s %s %s %s %s %s\n' % (ben_code_base,obsid_string,sister_obsid_string,track_off_moon_list_string,track_moon_string,track_off_moon_string,epoch_ID_string,flag_ants_string,cleanup_string,time_freq_res_string,machine_string,dysco_string) )

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
#parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files  e.g. --tagname="" [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='List of antennas (space separated) to flag after cottering (andre indexing as for rfigui etc)  e.g. --flag_ants="56,60" [default=%default]')
parser.add_option('--cleanup',action='store_true',dest='cleanup',default=False,help='Delete the gpubox files after making the ms [default=%default]')
parser.add_option('--obsid_infile',type='string', dest='obsid_infile',default='',help='File containing list of obsids to be cottered  e.g. --obsid_infile="20180107_moon_93.txt" [default=%default]')
parser.add_option('--sister_obsid_infile',type='string', dest='sister_obsid_infile',default='',help='File containing list of LST-matched sister observations for flag merging  e.g. --sister_obsid_infile="20180110_off_moon_93.txt" [default=%default]')
#parser.add_option('--time_freq_res',type='string', dest='time_freq_res',default='8,80',help='Time and then frequency resolution, comma separated e.g. --time_freq_res="8,80" [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')
parser.add_option('--no_dysco',action='store_true',dest='no_dysco',default=False,help='Do not use Dysco compression [default=%default]')

(options, args) = parser.parse_args()


generate_cotter_moon(options)
