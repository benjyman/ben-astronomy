#!/usr/bin/env python
#script to take a list of obsids and produce qa plots
#use: sshfs magnus:/ /md0/galaxy/

import os

def qa_check(options):
   obslist_filename = options.obslist
   data_dir = options.data_dir
   with open(obslist_filename) as f:
      lines = f.readlines()
   for line in lines:   
      obsid = line.strip()
      obs_dir = "%s/%s" % (data_dir,obsid)
      cmd = "aoqplot -save %s_data_stddev_baselines StandardDeviation %s/%s.ms" % (obsid,obs_dir,obsid)
      print cmd
      os.system(cmd)
   
   figname = "all_data_stddev_baselines-baselines.pdf"
   cmd = "rm -f %s" % figname
   print cmd
   os.system(cmd) 
   cmd = '/usr/bin/montage  *_data_stddev_baselines-baselines.pdf  -mode concatenate -tile 3x5 %s' % (figname)
   print cmd
   os.system(cmd)  

from optparse import OptionParser,OptionGroup

usage = 'Usage: qa_check.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--obslist',type='string', dest='obslist',default='obslist.txt',help='list of obsids e.g. --obslist="myobs.txt" [default=%default]')
parser.add_option('--data_dir',type='string', dest='data_dir',default='/astro/mwaeor/MWA/data/',help='data directory e.g. --obslist="/astro/mwaeor/bmckinley/data/" [default=%default]')
parser.add_option('--corrected_data',action='store_true',dest='corrected_data',default=False,help='Also collect statistics and run qa on CORRECTED_DATA column [default=%default]')

(options, args) = parser.parse_args()

qa_check(options)

