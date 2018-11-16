#!/usr/bin/env python
#script to make a histogram of the source flux densities in a sourcelist

import matplotlib.pyplot as plt


def flux_density_hist(options):
   #open the sourcelist file and extract data
   with open(options.source_list) as f:
      lines = f.readlines()

   flux_density_list = []
   for line in lines:
      if "fluxdensity" in line:
         flux_density = float(line.split()[2])
         flux_density_list.append(flux_density)

   plt.clf()
   plot_title = "Histogram with 'auto' bins"
   plt.hist(flux_density_list, bins='auto')  # arguments are passed to np.histogram
   plt.title(plot_title)
   fig_name="flux_density_histogram_%s.png" % (options.source_list.split('.')[0])
   figmap = plt.gcf()
   figmap.savefig(fig_name)





import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: flux_density_hist.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--source_list',type='string', dest='source_list',default='',help='Filename of sourcelist in skymodel 1.0 format [default=%default]')
parser.add_option('--plot_only',action='store_true',dest='plot_only',default=False,help='Just make the plots and maps from pre-calculated data e.g. --plot_only')


(options, args) = parser.parse_args()

flux_density_hist(options)