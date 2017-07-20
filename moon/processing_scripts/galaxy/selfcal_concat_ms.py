#! /usr/bin/env python

#This is to be run after cleaning so that the model from the clean is in the ms already - andre's calibrate tool will use this model
#can also be run to use a model in Andre format to calibrate, or just to apply solutions already made.
#option to track the moon.
import string 

def selfcal_concat_ms(obsid_string,list_index,track_off_moon_string,options):

   list_index=int(list_index)
   obsid_list=obsid_string.split(',')
   obsid=obsid_list[list_index].strip()
   print obsid
   
   if (options.track_off_moon):
      track_off_moon_list=track_off_moon_string.split(',')
      track_off_moon_paired_obsid=track_off_moon_list[int(float(list_index)*3)].strip()
      track_off_moon_new_RA=track_off_moon_list[int(list_index*3+1)].strip()
      track_off_moon_new_DEC=track_off_moon_list[int(list_index*3+2)].strip()
      print "obs paired with %s centering at RA:%s DEC:%s" % (track_off_moon_paired_obsid,track_off_moon_new_RA,track_off_moon_new_DEC)
   
   data_dir='%sdata/%s/' % (mwa_dir,obsid)
   metafits_file_name = data_dir+obsid+".metafits"

   if (options.tagname):
      tagname=options.tagname
   else:
      tagname='uvdump'

   if (options.model):
      model_string=' -m %s -applybeam ' % options.model
   else:
      model_string=' '

   if (options.sourcelist and not options.model):
      #make an ao model from the the rts sourclist
      ao_model_name=options.sourcelist.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt" 
      print "making %s " % ao_model_name
      cmd='aprun -n 1  python /group/mwaeor/CODE/MWA_Tools/scripts/srclist_by_beam.py -x --aocalibrate -m %s -n 1000 -s %s ' % (metafits_file_name,options.sourcelist)
      print cmd
      os.system(cmd)
      model_string=' -m %s -applybeam ' % ao_model_name
   else:
      model_string=' '  

   solutions_base_name='%s_%s_selfcal_%s_concat_solutions' % (obsid,tagname,options.selfcal) 
   if (options.chgcentre):
      solutions_base_name+='_newcentre'
   if (options.track_moon):
      solutions_base_name+='_trackmoon'
   if (options.track_off_moon):
      solutions_base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid   
   if (options.ionpeel):
      solutions_base_name=solutions_base_name+'_ionpeel'
      ionpeel_sourcelist=options.ionpeel
      clustered_model_name="clustered_" + options.ionpeel.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt"
    
      print "making %s " % clustered_model_name
      cmd='aprun -n 1 cluster %s %s 25 ' % (options.ionpeel.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt",clustered_model_name)
      print cmd
      os.system(cmd)


   solutions_name=solutions_base_name+'.bin'

   if(options.cotter):
      concat_vis_base='%s_cotter_%s' % (obsid,tagname)
   else:
      concat_vis_base='%s_%s_concat_transform' % (obsid,tagname)

   if (options.chgcentre):
      concat_vis_base+='_newcentre'
   if (options.track_moon):
      concat_vis_base=concat_vis_base+'_trackmoon'
   if (options.track_off_moon):
      concat_vis_base+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
   if (options.minw):
      concat_vis_base+='_minw'
    
   if (options.chgcentre  or options.minw):
      concat_vis_name=concat_vis_base+'.ms'
   else:
      concat_vis_name=data_dir+concat_vis_base+'.ms'
   
   if (options.applyonly):
      solutions=options.applyonly
      #apply solutions
      cmd='aprun -n 1 applysolutions %s %s ' % (concat_vis_name,solutions)
      print cmd
      os.system(cmd)
   
   else:
      if (options.ionpeel):
         peeled_ms_name=data_dir+concat_vis_base+'_peeled.ms'
         #remove any previous attempt
         cmd='aprun -n 1 rm -rf %s' % (peeled_ms_name)
         print cmd
         os.system(cmd)
         #make a copy of the unpeeled ms just in case
         cmd='aprun -n 1 cp -r %s %s' % (concat_vis_name,peeled_ms_name)
         print cmd
         os.system(cmd) 

         #run andres ionpeel 
         cmd='aprun -n 1 ionpeel -t 8 %s %s %s ' % (peeled_ms_name,clustered_model_name,solutions_name)
         print cmd
         os.system(cmd) 
   
      else:
         #run andres calibrate 
         cmd='aprun -n 1 calibrate -minuv 60  %s %s %s ' % (model_string,concat_vis_name,solutions_name)
         print cmd
         os.system(cmd)

         #apply solutions
         cmd='aprun -n 1 applysolutions %s %s ' % (concat_vis_name,solutions_name)
         print cmd
         os.system(cmd)


##
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: selfcal_concat_ms.py [obsid_list] [list_index] [options] '

parser = OptionParser(usage=usage)

parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files to calbrate  e.g. --tagname="" [default=%default]')
parser.add_option('--chgcentre',action='store_true',dest='chgcentre',default=False,help='Calibrate the ms that has had a phase centre shift  e.g. --chgcentre   [default=%default]')
parser.add_option('--ionpeel',type='string', dest='ionpeel',default='',help='Use this to specify a base sourcelist from which a clustered model will be made for each obsid and then run ionpeel e.g. ionpeel="sourcelist_PUMA1.txt" [default=%default]')
parser.add_option('--minw',action='store_true',dest='minw',default=False,help='Calibrate the ms that has had a phase centre shift to minw e.g. --minw   [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=False,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')
parser.add_option('--model',type='string', dest='model',default=None,help='Specify a model to calibrate on, model should have absolute, not apparent, flux densities (i.e. no beam appplied) e.g. model="clustered_model.txt" [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default="0",help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--applyonly',type='string', dest='applyonly',default=None,help='Just apply the specified calibration solutions to the measurements set e.g. applyonly=new_solutions.bin [default=%default]')
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (specify this instead of model and an ao model will be made from the sourclist), overidden by --model=')

(options, args) = parser.parse_args()

obsid_string = args[0]
list_index = args[1]
if (options.track_off_moon):
   track_off_moon_string= args[2]
else:
   track_off_moon_string=' '

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

selfcal_concat_ms(obsid_string,list_index,track_off_moon_string,options)

