 #!/usr/bin/env python 
#generate the track_off_moon_xxx.txt file
import os

def make_track_off_moon(centre_chan,options):
   centre_chan=options.centre_chan.strip()
   on_moon_dir=options.on_moon_dir.strip()
   on_moon_date=options.on_moon_date.strip()
   off_moon_date=options.off_moon_date.strip()
   
   on_moon_file="%s/%s/%s_moon_%s.txt" % (on_moon_dir,centre_chan,on_moon_date,centre_chan)
   track_off_moon_outfile_name="track_off_moon_%s_%s.txt" % (off_moon_date,centre_chan)

   track_off_moon_list=[]
   obsid_list=[]
   for line in open(moon_file):
      obsid_list.append(line.strip())

   for moon_obsid in obsid_list:
      print moon_obsid
      cmd="cat %s/%s/*.out | grep cotter | grep centre | grep %s | grep -o 'centre.*' | cut -f2- -d' ' > temp_off_moon_%s_%s.txt " % (moon_dir,centre_chan,moon_obsid,moon_obsid,centre_chan)
      print cmd
      os.system(cmd)
      for line in  open("temp_off_moon_%s_%s.txt" % (moon_obsid,centre_chan)):
         centre_coords="%s %s" % (line.split()[0],line.split()[1])
         if "%s %s" % (moon_obsid,centre_coords) not in track_off_moon_list:
            track_off_moon_list.append("%s %s" % (moon_obsid,centre_coords))

   #print track_off_moon_list
   track_off_moon_outfile=open(track_off_moon_outfile_name,"w")
   track_off_moon_outfile.write("\n".join(track_off_moon_list))
   track_off_moon_outfile.close()

import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: make_track_off_moon.py [centre_chan]'

parser = OptionParser(usage=usage)

parser.add_option('--centre_chan',type='string', dest='centre_chan',default='',help='Centre coarse chan of the observations  e.g. --centre_chan="121" [default=%default]')
parser.add_option('--on_moon_dir',type='string', dest='on_moon_dir',default='',help='Directory containing the on_moon files  e.g. --on_moon_dir="/md0/moon/2018A/01_2018A/on_moon" [default=%default]')
parser.add_option('--on_moon_date',type='string', dest='on_moon_date',default='',help='Date of on_moon obs  e.g. --on_moon_date="20180107" [default=%default]')
parser.add_option('--off_moon_date',type='string', dest='off_moon_date',default='',help='Date of off_moon obs  e.g. --off_moon_date="20180110" [default=%default]')


(options, args) = parser.parse_args()

make_track_off_moon(centre_chan,options)
