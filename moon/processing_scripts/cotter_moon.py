#python script to cotter moon obs

from astropy.io import fits
import string
import os.path

def cotter_moon(obsid_string,list_index,track_off_moon_string,options):
   
   obsid_list=obsid_string.split(',')
   obsid=obsid_list[int(list_index)].strip()
   print obsid

   list_index=int(list_index)
   
   if options.track_off_moon:
      track_off_moon_list=track_off_moon_string.split(',')

      track_off_moon_paired_obsid=track_off_moon_list[int(float(list_index)*3)].strip()

      track_off_moon_new_RA=track_off_moon_list[int(list_index*3+1)].strip()

      track_off_moon_new_DEC=track_off_moon_list[int(list_index*3+2)].strip()

      print "obs paired with %s centering at RA:%s DEC:%s" % (track_off_moon_paired_obsid,track_off_moon_new_RA,track_off_moon_new_DEC)
 

   data_dir='%sdata/%s/' % (mwa_dir,obsid)
   
   tagname=options.tagname
   base_name=obsid+'_cotter_'+tagname

   track_moon_string=' '
   if (options.track_moon):
      base_name+='_trackmoon'

   if (options.track_off_moon):
      base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid

   ms_name=data_dir+base_name+'.ms'
   
   metafits_filename=data_dir+obsid+'.metafits'
   flag_file_name=data_dir+obsid+'_flags.zip'

   #check if metafits file exists, if not, download it
   if not os.path.isfile(metafits_filename):
      cmd="make_metafits.py -o %s -g %s" % (metafits_filename,obsid)
      print cmd
      os.system(cmd)
   else:
      pass

   #unpack the flagfile
   cmd="unzip %s -d %s " % (flag_file_name,data_dir)
   print cmd
   os.system(cmd)

   #if tracking moon, find moon position
   if (options.track_moon):
      #get the date and time of the observation from the metafits file
      try:
         HDU_list = fits.open(metafits_filename)
      except IOError, err:
         'Cannot open metadata file %s\n' % str(options.input_file)
      header=HDU_list[0].header 
      date = (header)['DATESTRT']
      print date
      #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
      new_date=string.replace(string.replace(date,'-','/'),'T',' ')
      print new_date
      #find position of moon
      src_file_name='src_file_moon_%s.txt' % obsid 
      cmd="print_src.py --date='"+new_date+"' > "+ src_file_name
      print cmd
      os.system(cmd)
      #read src_file.txt to get Moon ra and dec 
      with open(src_file_name, "r") as infile:
         lines=infile.readlines()
         moon_ra=lines[6].split()[3].split(',')[0]
         moon_dec=lines[6].split()[6].split(',')[0]
         #print moon_ra
         #print moon_dec
      #get ra and dec in right format for cotter
      new_moon_dec=string.replace(moon_dec,":",".")
      track_moon_string=' -centre %s %s ' % (moon_ra,new_moon_dec)
      print track_moon_string

   #if tracking off moon, shift to required moon  position from a previous time
   if (options.track_off_moon):
      off_moon_ra=track_off_moon_new_RA
      off_moon_dec=track_off_moon_new_DEC
      track_moon_string=' -centre %s %s ' % (off_moon_ra,off_moon_dec)
      print track_moon_string

   cmd='aprun -n 1 cotter -flagfiles %s%s_%s.mwaf -norfi -o %s %s -m %s -timeres 8 -freqres 80 %s*gpubox*.fits' % (data_dir,obsid,"%%",ms_name,track_moon_string,metafits_filename,data_dir)
   print cmd
   os.system(cmd)

   if (options.flag_ants):
      flag_ants_cmd_string="flagantennae %s %s " % (ms_name,options.flag_ants)
      print cmd
      os.system(cmd)
   

import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: cotter_moon.py [obsid_string] [list_index] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms file e.g. --tagname="" [default=%default]')
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='List of antennas (space separated) to flag after cottering (andre indexing as for rfigui etc)  e.g. --flag_ants="56,60" [default=%default]')

(options, args) = parser.parse_args()

obsid_string = args[0]
list_index = args[1]
if (options.track_off_moon):
   track_off_moon_string= args[2]
else:
   track_off_moon_string=' '

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

cotter_moon(obsid_string,list_index,track_off_moon_string,options)


