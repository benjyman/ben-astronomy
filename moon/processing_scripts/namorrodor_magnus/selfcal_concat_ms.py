#! /usr/bin/env python

#This is to be run after cleaning so that the model from the clean is in the ms already - andre's calibrate tool will use this model
#can also be run to use a model in Andre format to calibrate, or just to apply solutions already made.
#option to track the moon.
import string 
from casacore.tables import table,tablecolumn,tablerow
import numpy as np

def selfcal_concat_ms(obsid,track_off_moon_string,options):
   #tagname=options.tagname
   epoch_ID=options.epoch_ID
   
   machine=options.machine
   if machine=='namorrodor':
      if epoch_ID=='2015B_05':
         mwa_dir='/data/MWA/'
      else:
         mwa_dir = '/md0/moon/data/MWA/'
   else:
      mwa_dir = '/astro/mwaeor/MWA/data/'
  
   print obsid

   if options.sister_obsid:   
      sister_obsid=options.sister_obsid
   
      
   if (options.track_off_moon):
      track_off_moon_list=track_off_moon_string.split(',')

      track_off_moon_paired_obsid=track_off_moon_list[0].strip()

      track_off_moon_new_RA=track_off_moon_list[1].strip()

      track_off_moon_new_DEC=track_off_moon_list[2].strip()

      print "obs paired with %s centering at RA:%s DEC:%s" % (track_off_moon_paired_obsid,track_off_moon_new_RA,track_off_moon_new_DEC)
   
   data_dir='%s%s/' % (mwa_dir,obsid)
   metafits_file_name = "%s%s_metafits_ppds.fits" % (data_dir,obsid)

   if (options.model):
      model_string=' -m %s -applybeam ' % options.model
   else:
      model_string=' '

   if (options.sourcelist and not options.model):
      #make an ao model 
      ao_model_name_cwd=options.sourcelist.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt" 
      ao_model_name="%s%s" % (data_dir,ao_model_name_cwd)
      #if (not os.path.exists(ao_model_name)):
      print "making %s " % ao_model_name
      cmd='srclist_by_beam.py -x --aocalibrate -m %s -n 1000 -s %s ' % (metafits_file_name,options.sourcelist)
      os.system(cmd)
      cmd="mv %s %s" % (ao_model_name_cwd,data_dir)
      os.system(cmd)
      #else:
      #   print "using sourcelist %s " % ao_model_name
      model_string=' -m %s -applybeam ' % ao_model_name
   else:
      model_string=' '  

   solutions_base_name='%s%s_%s_selfcal_%s_concat_solutions' % (data_dir,obsid,epoch_ID,options.selfcal) 
   if (options.chgcentre):
      solutions_base_name+='_newcentre'
   if (options.track_moon):
      solutions_base_name+='_trackmoon'
   if (options.track_off_moon):
      solutions_base_name+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid   
   if (options.ionpeel):
      solutions_base_name=solutions_base_name+'_ionpeel'
      ionpeel_sourcelist=options.ionpeel
      clustered_model_name=data_dir+"clustered_10dirs_" + options.ionpeel.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt"
    
      print "making %s " % clustered_model_name
      #cmd='cluster %s %s 10 ' % (options.ionpeel.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt",clustered_model_name)
      ao_model_name_cwd=options.ionpeel.split("/")[-1].split(".")[0] + "_" + obsid + "_aocal1000.txt"
      ao_model_name="%s%s" % (data_dir,ao_model_name_cwd)
      cmd='cluster %s %s 10 ' % (ao_model_name,clustered_model_name)
      print cmd
      os.system(cmd)


   solutions_name=solutions_base_name+'.bin'

   if(options.cotter):
      #concat_vis_base='%s_%s' % (obsid,epoch_ID)
      if epoch_ID=='2015B_05':
         if not options.track_off_moon:
            concat_vis_base='%s_cotter_%s' % (obsid,'20150926_moon_69')
         else:
            concat_vis_base='%s_cotter_%s' % (obsid,'20150929_moon_69')
      else:
         #concat_vis_base='%s_cotter_%s' % (obsid,tagname)
         concat_vis_base='%s_%s' % (obsid,epoch_ID)

   if (options.chgcentre):
      concat_vis_base+='_newcentre'
   if (options.minw):
      concat_vis_base+='_minw'
   
   #do all this stuff in the selfcal stage instead - so you are sure that both ms have already meen created
   #here merge the flag columns for both the on_moon and off_moon ms
   #only do this step for track_off_moon, (always need to make on_moon ms first)
   if (options.track_moon):
      #base_name='%s_%s' % (obsid,epoch_ID)
      on_moon_basename=concat_vis_base + '_trackmoon'
      on_moon_ms_name=data_dir+on_moon_basename+'.ms'
      off_moon_base_name="%s_%s_track_off_moon_paired_%s" % (sister_obsid,epoch_ID,obsid)
      off_moon_ms_name=mwa_dir+sister_obsid+'/'+off_moon_base_name+'.ms'
      concat_vis_name=on_moon_ms_name 
   elif (options.track_off_moon):
      on_moon_basename="%s_%s_trackmoon" % (sister_obsid,epoch_ID) 
      on_moon_ms_name=mwa_dir+sister_obsid+'/'+on_moon_basename+'.ms'
      off_moon_base_name=concat_vis_base+'_track_off_moon_paired_' + sister_obsid
      off_moon_ms_name=data_dir+off_moon_base_name+'.ms'
      concat_vis_name=off_moon_ms_name
   else:
      if (options.chgcentre  or options.minw):
         concat_vis_name=concat_vis_base+'.ms'
      else:
         if (machine == "magnus" or machine == "galaxy"):
            if (options.track_moon or options.track_off_moon):
               concat_vis_name=data_dir+concat_vis_base+'.ms'
            else:
               concat_vis_name=concat_vis_base+'.ms'
         else:
            concat_vis_name=data_dir+concat_vis_base+'.ms'

   if (options.track_off_moon):
      on_moon_table=table(on_moon_ms_name,readonly=False)
      #print on_moon_table
      #on_moon_table_UVW=tablecolumn(on_moon_table,'UVW')
      #on_moon_table_time=tablecolumn(on_moon_table,'TIME')
      on_moon_table_flag=tablecolumn(on_moon_table,'FLAG')
   
      off_moon_table=table(off_moon_ms_name,readonly=False)
      #print off_moon_table
      #off_moon_table_UVW=tablecolumn(off_moon_table,'UVW')
      #off_moon_table_time=tablecolumn(off_moon_table,'TIME')
      off_moon_table_flag=tablecolumn(off_moon_table,'FLAG')
   
      new_off_moon_table_flag=np.logical_not(off_moon_table_flag)

   
      #look in obs_list_generator.py for how to compare LSTs
   
      #Not sure if 2018A obs are well enough LST matched
      #Continue anyway (come back and look at 2015A) 
   
      #OR the two flag columns
      new_flag_column=np.logical_or(on_moon_table_flag,off_moon_table_flag)
      #print new_flag_column[20100:20101]
   
      on_moon_table.putcol('FLAG',new_flag_column)
      off_moon_table.putcol('FLAG',new_flag_column)
   
      on_moon_table.close()
      off_moon_table.close()
   
   
   if (options.applyonly):
      solutions=options.applyonly
      #apply solutions
      cmd='applysolutions %s %s ' % (concat_vis_name,solutions)
      print cmd
      os.system(cmd)
   
   else:
      if (options.ionpeel):
         if (machine == "magnus" or machine == "galaxy"): 
            if (options.track_moon or options.track_off_moon):
               peeled_ms_name=data_dir + concat_vis_base+'_peeled.ms'
            else:
               peeled_ms_name=concat_vis_base+'_peeled.ms'
         else:
            peeled_ms_name=data_dir+concat_vis_base+'_peeled.ms'
         #remove any previous attempt
         cmd='rm -rf %s' % (peeled_ms_name)
         print cmd
         os.system(cmd)
         #make a copy of the unpeeled ms just in case
         cmd='cp -r %s %s' % (concat_vis_name,peeled_ms_name)
         print cmd
         os.system(cmd) 

         #run andres ionpeel 
         cmd='ionpeel -t 8 %s %s %s ' % (peeled_ms_name,clustered_model_name,solutions_name)
         print cmd
         os.system(cmd) 
   
      else:
         #run andres calibrate 
         cmd='calibrate -minuv 60  %s %s %s ' % (model_string,concat_vis_name,solutions_name)
         print cmd
         os.system(cmd)

         #apply solutions
         cmd='applysolutions %s %s ' % (concat_vis_name,solutions_name)
         print cmd
         os.system(cmd)


##
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: selfcal_concat_ms.py [obsid_list] [list_index] [options] '

parser = OptionParser(usage=usage)

parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
#parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files to calbrate  e.g. --tagname="" [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='',help='epoch_ID of observations  e.g. --epoch_ID="2018A_01" [default=%default]')
parser.add_option('--chgcentre',action='store_true',dest='chgcentre',default=False,help='Calibrate the ms that has had a phase centre shift  e.g. --chgcentre   [default=%default]')
parser.add_option('--ionpeel',type='string', dest='ionpeel',default='',help='Use this to specify a base sourcelist from which a clustered model will be made for each obsid and then run ionpeel e.g. ionpeel="sourcelist_PUMA1.txt" [default=%default]')
parser.add_option('--minw',action='store_true',dest='minw',default=False,help='Calibrate the ms that has had a phase centre shift to minw e.g. --minw   [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=True,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')
parser.add_option('--model',type='string', dest='model',default=None,help='Specify a model to calibrate on, model should have absolute, not apparent, flux densities (i.e. no beam appplied) e.g. model="clustered_model.txt" [default=%default]')
parser.add_option('--selfcal',type='string', dest='selfcal',default="0",help='Specify how many times the image has been through selfcal, labels the solutions so that you can revert to an older calibration e.g. selfcal=2 [default=%default]')
parser.add_option('--applyonly',type='string', dest='applyonly',default=None,help='Just apply the specified calibration solutions to the measurements set e.g. applyonly=new_solutions.bin [default=%default]')
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (specify this instead of model and an ao model will be made from the sourclist), overidden by --model=')
parser.add_option('--machine',type='string', dest='machine',default='magnus',help='machine can be galaxy, magnus or namorrodor e.g. --machine="namorrodor" [default=%default]')
parser.add_option('--sister_obsid',type='string', dest='sister_obsid',default='',help='sister_obsid e.g. --sister_obsid="1199396880" [default=%default]')


(options, args) = parser.parse_args()

obsid = args[0]
if (options.track_off_moon):
   track_off_moon_string= args[1]
else:
   track_off_moon_string=' '


selfcal_concat_ms(obsid,track_off_moon_string,options)

