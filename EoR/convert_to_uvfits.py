#Script using pyuvdata to convert MWA cotter ms to uvfits
from pyuvdata import UVData
import copy
import numpy as np

def convert(options):
   ms_name=options.ms_name
   number_output_chans=int(options.channels_out)
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


###########################
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage convert_to_uvfits.py [options]'

parser = OptionParser(usage=usage)

#parser.add_option('--do_beam_stuff',action='store_true',dest='do_beam_stuff',default=False,help='set if you want to do all the beam weighted averaging stuff e.g. --do_beam_stuff')
parser.add_option('--ms_name',type='string', dest='ms_name',default='',help='Name of ms to convert to uvfits e.g. --ms_name=1102865728.ms [default=%default]')
parser.add_option('--channels_out',type='string', dest='channels_out',default='1',help='Number of channels to split ms into and output seperate uvfits files e.g. --channels_out=24 [default=%default]')


(options, args) = parser.parse_args()

convert(options)



