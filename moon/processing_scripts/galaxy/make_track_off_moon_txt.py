#!/usr/bin/env python 
#generate the track_off_moon_xxx.txt file
import os

def make_track_off_moon(centre_chan,options):
   centre_chan =str(centre_chan).strip()
   moon_dir="/scratch2/mwaeor/bmckinley/Moon/2016/20150926_moon1"
   moon_file="%s/%s/20150926_moon_%s.txt" % (moon_dir,centre_chan,centre_chan)
   track_off_moon_outfile_name="track_off_moon_20150926_%s.txt" % (centre_chan)

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

(options, args) = parser.parse_args()

centre_chan = args[0]

make_track_off_moon(centre_chan,options)
