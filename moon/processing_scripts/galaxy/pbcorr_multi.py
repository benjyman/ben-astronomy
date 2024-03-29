#! /usr/bin/env python

#script to pbcorrect images

import numpy as np

def pbcorr_multi(obsid_string,list_index,track_off_moon_string,options):

   #do all the next steps for each chan
   number_images=int(options.channelsout)

   list_index=int(list_index)
   obsid_list=obsid_string.split(',')
   n_obs=len(obsid_list)
   print n_obs
   obsid_list_index=int(np.floor(list_index/number_images))
   print obsid_list_index
   obsid=obsid_list[obsid_list_index].strip()
   print obsid

   image_number=list_index-(obsid_list_index*number_images)
   print image_number

   if (options.track_off_moon):
      track_off_moon_list=track_off_moon_string.split(',')
      track_off_moon_paired_obsid=track_off_moon_list[int(float(obsid_list_index)*3)].strip()
   
   data_dir='%sdata/%s/' % (mwa_dir,obsid)

   if (options.tagname):
      tagname=options.tagname
   else:
      tagname='uvdump'

   if (options.ms_base_name):
      ms_name='%s%s_%s.ms' % (data_dir,obsid,options.ms_base_name)
   else:
      ms_base_name='%s_cotter_%s' % (obsid,tagname)
      if (options.track_moon):
         ms_base_name+='_trackmoon'
      if (options.track_off_moon):
         ms_base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
      if (options.ionpeeled):
         ms_base_name+='_peeled'
      ms_name=data_dir+ms_base_name+'.ms'

   if (options.image_base_name):
      image_base=options.image_base_name
      image_base_name='%s_%s' % (obsid,image_base)
   else:
      image_base_name='%s_cotter_%s' % (obsid,tagname)

      if (options.chgcentre):
         image_base_name+='_newcentre'
      if (options.minw):
         image_base_name+='_minw'
      if (options.track_moon):
         image_base_name+='_trackmoon'
      if (options.track_off_moon):
         image_base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
      if (options.ionpeeled):
         image_base_name+='_peeled'
      
   if (options.metafits_name):
      metafits_name=options.metafits_name
   else:
      metafits_name=' %s%s.metafits ' % (data_dir,obsid)
   
   if (options.havebeam):
      if (options.track_off_moon):
         beam_base_name=image_base_name+'_beam'
      else:
         beam_base_name='%s_%s_beam' % (obsid,options.havebeam)
   else:
      beam_base_name=image_base_name+'_beam'

   if (options.applyion):
      solutions_base_name='%s_%s_selfcal_%s_concat_solutions' % (obsid,tagname,options.selfcal)
      if (options.track_moon):
         solutions_base_name+='_trackmoon'
      if (options.track_off_moon):
         solutions_base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
      solutions_base_name=solutions_base_name+'_ionpeel'
      ionpeel_sourcelist=options.applyion
      clustered_model="clustered_" + options.applyion.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt"

      ion_solutions=solutions_base_name+'.bin'

      #flag the bad ionsolutions so applyion doesn't crash
      flagged_ion_solutions='flagged_%s' % (ion_solutions)
      cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/ionsol -o %s -threshold-dl 4e-4 -threshold-dm 4e-4 -max-gain 2 -min-gain 0.2 -replace-unset %s %s ' % (flagged_ion_solutions,ion_solutions,clustered_model)
      print cmd
      os.system(cmd)

   #go!
   if number_images>1:
      image_base_name_loop='%s-%04d' % (image_base_name,image_number)
      beam_base_name_loop='%s-%04d' % (beam_base_name,image_number)
      print image_base_name_loop
      print beam_base_name_loop

   if not (options.havebeam):
      cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/beam -2014i -proto '+ image_base_name_loop+'-XX-image.fits -name '+ beam_base_name_loop+ ' -ms '+ms_name+ ' -m ' +metafits_name
      print cmd
      os.system(cmd)
   
   if (options.dirty):
      cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/pbcorrect '+image_base_name_loop+ ' dirty.fits '+beam_base_name_loop+' '+ image_base_name_loop+'_dirty'
   else:
      cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/pbcorrect '+image_base_name_loop+ ' image.fits '+beam_base_name_loop+' '+ image_base_name_loop
   print cmd
   os.system(cmd)

   if (options.applyion):
      if (options.dirty):
         pbcorr_image_name='%s_dirty-I.fits' % (image_base_name_loop)
         applied_image_base='%s_dirty_applied' % (image_base_name_loop)
      else:
         pbcorr_image_name='%s-I.fits' % (image_base_name_loop)
         applied_image_base='%s_applied' % (image_base_name_loop)
      
      applied_image_name=applied_image_base+'-I.fits' 
      cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/applyion %s %s %s %s ' % (pbcorr_image_name,applied_image_name,clustered_model,flagged_ion_solutions)
      print cmd
      os.system(cmd)
      if (options.pbuncorrect and not options.render):
         cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/pbcorrect -uncorrect %s image.fits %s %s ' % (applied_image_base,beam_base_name_loop,applied_image_base)
         print cmd
         os.system(cmd)

      if (options.render):
         rendered_image_base='%s_rendered' % (applied_image_base)
         rendered_image_name='%s-I.fits' % rendered_image_base
         cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/render -a -r -ion %s rendered -t %s -o %s %s' % (ion_solutions,applied_image_name,rendered_image_name,clustered_model)
         print cmd
         os.system(cmd)
         
         if (options.pbuncorrect):
            cmd='aprun -n 1 /group/mwasci/bmckinley/code/anoko/mwa-reduce/build/pbcorrect -uncorrect %s image.fits %s %s ' % (rendered_image_base,beam_base_name_loop,rendered_image_base) 
            print cmd
            os.system(cmd)
         
##
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: pbcorr_multi.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--image_base_name',type='string', dest='image_base_name',default='',help='image_base_name to be pb corrected, not required if using standard names based on tagname (i.e. leave it out unless you are doing something special)  e.g. --image_base_name="" [default=%default]')
parser.add_option('--havebeam', type='string', dest='havebeam', default=None,help='Set this if you already have downloaded the right beam model and give the beam prefix, not inclusing the leading obsid  e.g. --havebeam=tagname (ie without the _beam at the end) [default=%default]')
parser.add_option('--ms_base_name',type='string', dest='ms_base_name',default='',help='name of ms to be used to generate beam, not including the leading obsid or the .ms (not needed if havebeam=True, not required if using standard names based on tagname (i.e. leave it out unless you are doing something special)  e.g. --ms_name="blah.ms" [default=%default]')
parser.add_option('--metafits_name',type='string', dest='metafits_name',default=None,help='name of metafits file to be used to generate beam (not needed if havebeam=True, should be same obsid as ms_name) only specify if correcting an image made from several obsids - for single obs images leave blank  e.g. --metafits_name="blah.metafits" [default=%default]')
parser.add_option('--applyion',type='string', dest='applyion',default=None,help='apply the ionpeel corrections to the pbcorrected, specify the sourcelist from which the ionpeel solutions were created as in selfcal_concat.py i.e. applyion="srclist_IDR2deep_puma-noextended_HydAGuas.txt" [default=%default]')
parser.add_option('--clustered_model',type='string', dest='clustered_model',default=None,help='If applying ionsolutions, or rendering, specify the clustered model e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--render',action='store_true', dest='render',default=False,help='Render the clustered model to the ion-applied image e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--pbuncorrect',action='store_true', dest='pbuncorrect',default=False,help='pbuncorrect the rendered, ion-applied image so you can pbaddimage them together! e.g. clustered_model="555_ionpeel.bin" [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files to calbrate  e.g. --tagname="" [default=%default]')
parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default="0",help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--ionpeeled',action='store_true',dest='ionpeeled',default=False,help='Set if this is a subsequent imaging run after ionpeeling e.g. --ionpeeled [default=%default]')
parser.add_option('--chgcentre',type='string',dest='chgcentre',default=None,help='Use chgcentre to change phase centre to desired phase centre or shift to min w and shiftback e.g. --chgcentre=" -minw -shiftback " or --chgcentre=" 00h00m00.0s 00d00m00.0s "  [default=%default]')
parser.add_option('--minw',type='string',dest='minw',default=None,help='Shift to minw position of whatever ms is central in the chunk and then shiftback (must start at same phase centre which eor obs do e.g. --minw="12345678.ms"  [default=%default]')
parser.add_option('--dirty',action='store_true',dest='dirty',default=False,help='perform all operations on the dirty images, not the cleaned images [default=%default]')
parser.add_option('--channelsout',type='string', dest='channelsout',default='1',help='Specify how many channels were output from wsclean e.g.channelsout=24 [default=%default]')

(options, args) = parser.parse_args()

obsid_string = args[0]
list_index = args[1]
if (options.track_off_moon):
   track_off_moon_string= args[2]
else:
   track_off_moon_string=' '

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

pbcorr_multi(obsid_string,list_index,track_off_moon_string,options)


