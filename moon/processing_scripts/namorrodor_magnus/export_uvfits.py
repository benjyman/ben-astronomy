#python script to cotter moon obs

from astropy.io import fits
import string
import os.path
import cmd
import numpy as np
from pyuvdata import UVData


def export_uvfits(options):
   machine=options.machine
   obsid=options.obsid
   if machine=='namorrodor':
      mwa_dir = '/md0/moon/data/MWA/'
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'
      
   data_dir='%s%s/' % (mwa_dir,obsid)

   epoch_ID=options.epoch_ID
   base_name= "%s_%s" % (obsid,epoch_ID)

   ms_name=data_dir+base_name+'.ms'
   uvfits_name=data_dir+base_name+'.uvfits'
 
   UV=UVData()
   UV.read_ms(ms_name)
   UV.write_uvfits('uvfits_name',sppof_nonessential=True)
   

import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: export_uvfits.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--obsid',type='string', dest='obsid',default='',help='obsid to be cottered e.g. --obsid="1199394880" [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')


(options, args) = parser.parse_args()


export_uvfits(options)


