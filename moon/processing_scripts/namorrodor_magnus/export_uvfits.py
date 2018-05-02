#python script to cotter moon obs

from astropy.io import fits
import copy
import string
import os.path
import cmd
import numpy as np
from pyuvdata import UVData


def export_uvfits(options):
   machine=options.machine
   obsid=options.obsid
   number_output_chans=int(options.channels_out)
   if machine=='namorrodor':
      mwa_dir = '/md0/moon/data/MWA/'
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'
      
   data_dir='%s%s/' % (mwa_dir,obsid)

   epoch_ID=options.epoch_ID
   base_name= "%s_%s" % (obsid,epoch_ID)

   if options.ionpeeled:
      base_name+='_peeled'

   if (machine=="magnus" or machine=="galaxy"):  
      ms_name=base_name+'.ms'
   else:
      ms_name=data_dir+base_name+'.ms'
   
   uvdata_name_base=ms_name.split('.')[0]
 
   UV=UVData()
   UV.read_ms(ms_name)

   number_input_chans=UV.Nfreqs
   number_input_chans_in_ouput_chan=number_input_chans/number_output_chans
   print "number of frequency channels in ms is %s" % number_input_chans

   for chan in range(1,number_output_chans+1):
      uvfits_name="%s_%02d.uvfits" % (uvdata_name_base,chan)
      #select the data from the ms
      
      channel_start=(chan*number_input_chans_in_ouput_chan)
      channel_end=channel_start+number_input_chans_in_ouput_chan

      print "selecting channel %s to %s from %s for %s" % (channel_start,channel_end,ms_name,uvfits_name)         
      uvout=copy.deepcopy(UV)
      uvout.select(freq_chans=np.arange(channel_start, channel_end))
      uvout.write_uvfits(uvfits_name, spoof_nonessential=True)

import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: export_uvfits.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--obsid',type='string', dest='obsid',default='',help='obsid to be cottered e.g. --obsid="1199394880" [default=%default]')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')
parser.add_option('--ionpeeled',action='store_true',dest='ionpeeled',default=False,help='Set if the data have been ionpeeled e.g. --ionpeeled [default=%default]')
parser.add_option('--channels_out',type='string', dest='channels_out',default='1',help='Number of channels to split ms into and output seperate uvfits files e.g. --channels_out=24 [default=%default]')

(options, args) = parser.parse_args()


export_uvfits(options)


