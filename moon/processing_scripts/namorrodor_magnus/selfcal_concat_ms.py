#! /usr/bin/env python

#This is to be run after cleaning so that the model from the clean is in the ms already - andre's calibrate tool will use this model
#can also be run to use a model in Andre format to calibrate, or just to apply solutions already made.
#option to track the moon.
import string 
from casacore.tables import table,tablecolumn,tablerow
import numpy as np
from astropy.io import fits
import os.path
import cmd
import datetime


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
      if (not os.path.exists(data_dir+ao_model_name_cwd)):
         print "making %s " % ao_model_name
         cmd='srclist_by_beam.py -x --aocalibrate -m %s -n 1000 -s %s ' % (metafits_file_name,options.sourcelist)
         os.system(cmd)
         cmd="mv %s %s" % (ao_model_name_cwd,data_dir)
         os.system(cmd)
      else:
         print "using sourcelist %s " % ao_model_name
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
      #cluster in just 5 directions to be consistent with what RTS is doing
      cmd='cluster %s %s 5 ' % (ao_model_name,clustered_model_name)
      print cmd
      os.system(cmd)


   solutions_name=solutions_base_name+'.bin'

   if(options.cotter):
      concat_vis_base='%s_%s' % (obsid,epoch_ID)
      #if epoch_ID=='2015B_05':
      #   if not options.track_off_moon:
      #      concat_vis_base='%s_cotter_%s' % (obsid,'20150926_moon_69')
      #   else:
      #      concat_vis_base='%s_cotter_%s' % (obsid,'20150929_moon_69')
      #else:
      #   #concat_vis_base='%s_cotter_%s' % (obsid,tagname)
      #   concat_vis_base='%s_%s' % (obsid,epoch_ID)

   if (options.chgcentre):
      concat_vis_base+='_newcentre'
   if (options.minw):
      concat_vis_base+='_minw'
   
   #do all this stuff in the selfcal stage instead - so you are sure that both ms have already meen created
   #here merge the flag columns for both the on_moon and off_moon ms
   #only do this step for track_off_moon, (always need to make on_moon ms first)
   if (options.track_moon):
      #base_name='%s_%s' % (obsid,epoch_ID)
      on_moon_obsid = obsid
      off_moon_obsid = sister_obsid
      
      on_moon_basename=concat_vis_base + '_trackmoon'
      on_moon_ms_name=data_dir+on_moon_basename+'.ms'
      off_moon_base_name="%s_%s_track_off_moon_paired_%s" % (sister_obsid,epoch_ID,obsid)
      off_moon_ms_name=mwa_dir+sister_obsid+'/'+off_moon_base_name+'.ms'
      concat_vis_name=on_moon_ms_name 
   elif (options.track_off_moon):
      on_moon_obsid = sister_obsid
      off_moon_obsid = obsid
      
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

   if (options.manta_ray):      
      #if tracking moon, find moon position
      if (options.track_moon):
         #get the date and time of the observation from the metafits file
         try:
            HDU_list = fits.open(metafits_file_name)
         except IOError, err:
            'Cannot open metadata file %s\n' % str(options.input_file)
         header=HDU_list[0].header 
         date_unix_utc_string = (header)['GOODTIME']
         #datetime.datetime.utcfromtimestamp(posix_time).strftime('%Y-%m-%dT%H:%M:%SZ')
         date_unix_utc=datetime.datetime.utcfromtimestamp(date_unix_utc_string)
         print "start time (GOODTIME) of observation is %s" % (date_unix_utc)
         
         #Get the observation length
         obslength=float((header)['EXPOSURE'])
         print "EXPOSURE is %s s" % obslength
         half_obslength=obslength/2.
         half_obslength_timedelta=datetime.timedelta(seconds=half_obslength)
         middle_of_observation=date_unix_utc+half_obslength_timedelta
         print "middle_of_observation is %s" % middle_of_observation
         
         middle_of_observation_formatted=middle_of_observation.strftime('%Y-%m-%dT%H:%M:%S')
         print "middle_of_observation_formatted %s " % middle_of_observation_formatted
         
         #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
         new_date=string.replace(string.replace(middle_of_observation_formatted,'-','/'),'T',' ')
         print new_date
         #find position of moon
         src_file_name='src_file_moon_%s.txt' % obsid 
         cmd="print_src.py --date='"+new_date+"' > "+ src_file_name
         print cmd
         os.system(cmd)
         #read src_file.txt to get Moon ra and dec 
         with open(src_file_name, "r") as infile:
            lines=infile.readlines()
            moon_ra=lines[6].split()[3].split(',')[0]
            moon_dec=lines[6].split()[6].split(',')[0]
            #print moon_ra
            #print moon_dec
         #get ra and dec in right format for cotter
         new_moon_dec=string.replace(moon_dec,":",".")
         track_moon_string=' %s %s ' % (moon_ra,new_moon_dec)
         print track_moon_string
      
      
      
      #if tracking off moon, shift to required moon  position from a previous time
      if (options.track_off_moon):
         off_moon_ra=track_off_moon_new_RA
         off_moon_dec=track_off_moon_new_DEC
         track_moon_string=' %s %s ' % (off_moon_ra,off_moon_dec)
         print track_moon_string
      
      #for on_moon:  
      ##this unzip is done by gator now 
      #cmd = "unzip -o -d %s %s%s_ms.zip " % (data_dir,data_dir,obsid)
      #print cmd
      #os.system(cmd)
      
      ##remove if there is an existing ms:
      ##This is dont using moon_process_setup now
      #cmd = "rm -rf %s" % (concat_vis_name)
      #print cmd
      #os.system(cmd)  
      
      ##rename ms:
      #cmd = "mv %s%s.ms %s" % (data_dir,obsid,concat_vis_name)
      #print cmd
      #os.system(cmd)
      
      #change the phase centre of the ms
      cmd = "chgcentre %s %s" % (concat_vis_name,track_moon_string)
      print cmd
      os.system(cmd)

      ##repeat for off moon
      #off_moon_data_dir=mwa_dir+sister_obsid+'/'
      #cmd = "unzip -o -d %s %s%s_ms.zip " % (off_moon_data_dir,off_moon_data_dir,sister_obsid)
      #print cmd
      #os.system(cmd)
      
      ##remove if there is an existing ms:
      #cmd = "rm -rf %s" % (off_moon_ms_name)
      #print cmd
      #os.system(cmd)  
      
      ##rename ms:
      #cmd = "mv %s%s.ms %s" % (off_moon_data_dir,sister_obsid,off_moon_ms_name)
      #print cmd
      #os.system(cmd)
   
   #if (options.track_off_moon):
   #   #change the phase centre of the ms
   #   cmd = "chgcentre %s %s" % (concat_vis_name,track_moon_string)
   #   print cmd
   #   os.system(cmd)      

   if (options.flag_ants):
      #do both on and off moon
      flag_ants_cmd_string="flagantennae %s %s " % (on_moon_ms_name,options.flag_ants)
      print flag_ants_cmd_string
      os.system(flag_ants_cmd_string)

      flag_ants_cmd_string="flagantennae %s %s " % (off_moon_ms_name,options.flag_ants)
      print flag_ants_cmd_string
      os.system(flag_ants_cmd_string)     



   #apply the same flags to both ms for moon obs, do this in the on moon stage and get flags from metafits 
   if (options.track_moon and not options.ionpeel):      
      #get flagged antennas from the metafits file and flag them explicitly
      on_moon_metafits_file_name = "%s_metafits_ppds.fits" % on_moon_obsid
      off_moon_metafits_file_name = "%s_metafits_ppds.fits" % off_moon_obsid
      
      #on moon
      try:
         HDU_list = fits.open(metafits_file_name)
      except IOError, err:
         'Cannot open metadata file %s\n' % str(options.input_file)
      data = HDU_list[1].data
      tile_flags = data['FLAG']
      tile_flagging_indices = data['ANTENNA']
      tile_name = data['TILENAME']
      #print tile_flags
      
      flag_antenna_indices_list = []
      flag_antenna_tilenames_list = []
      for index,tile_flag in enumerate(tile_flags):
         if (tile_flag==1):
            string = "%s with antenna index %s should be flagged." % (tile_name[index],tile_flagging_indices[index])
            #print string
            if index==0:
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
            elif (tile_name[index]!=tile_name[index-1]):
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
      
      print flag_antenna_indices_list
      print flag_antenna_tilenames_list
      
      flag_ant_string = ' '.join(str(x) for x in flag_antenna_indices_list)
      
      #flag em
      cmd="flagantennae %s %s" % (on_moon_ms_name,flag_ant_string)
      print cmd
      os.system(cmd)

      #off moon
      try:
         HDU_list = fits.open(off_moon_metafits_file_name)
      except IOError, err:
         'Cannot open metadata file %s\n' % str(options.input_file)
      header=HDU_list[0].header
      data = HDU_list[1].data
      tile_flags = data['FLAG']
      tile_flagging_indices = data['ANTENNA']
      tile_name = data['TILENAME']
      #print tile_flags
      
      flag_antenna_indices_list = []
      flag_antenna_tilenames_list = []
      for index,tile_flag in enumerate(tile_flags):
         if (tile_flag==1):
            string = "%s with antenna index %s should be flagged." % (tile_name[index],tile_flagging_indices[index])
            #print string
            if index==0:
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
            elif (tile_name[index]!=tile_name[index-1]):
               flag_antenna_indices_list.append(tile_flagging_indices[index])
               flag_antenna_tilenames_list.append(tile_name[index])
      
      print flag_antenna_indices_list
      print flag_antenna_tilenames_list
      
      flag_ant_string = ' '.join(str(x) for x in flag_antenna_indices_list)
      
      #flag em
      cmd="flagantennae %s %s" % (off_moon_ms_name,flag_ant_string)
      print cmd
      os.system(cmd)

      
      with table(on_moon_ms_name,readonly=False) as on_moon_table:
         with table(off_moon_ms_name,readonly=False) as off_moon_table:
            with tablecolumn(on_moon_table,'FLAG') as on_moon_table_flag:
               with tablecolumn(off_moon_table,'FLAG') as off_moon_table_flag:
                   #new_moon_flag=np.empty_like(on_moon_table_flag)
                   for row_index,row in enumerate(on_moon_table_flag):
                      #if (on_moon_table_flag[row_index].all() != off_moon_table_flag[row_index]).all():
                      #   print on_moon_table_flag[row_index]
                      #   print off_moon_table_flag[row_index]
      
                      on_moon_table_flag[row_index] = np.logical_or(on_moon_table_flag[row_index],off_moon_table_flag[row_index])
                      off_moon_table_flag[row_index] = np.logical_or(on_moon_table_flag[row_index],off_moon_table_flag[row_index])
      
                      on_moon_table.putcell('FLAG',row_index,on_moon_table_flag[row_index])
                      off_moon_table.putcell('FLAG',row_index,off_moon_table_flag[row_index])

      #collect new statistics
      cmd = "aoquality collect %s" % (on_moon_ms_name)
      print cmd
      os.system(cmd)
      
      cmd = "aoquality collect %s" % (off_moon_ms_name)
      print cmd
      os.system(cmd)
      
      
      #save a png of the aoqplot output for StandardDeviation
      cmd = "aoqplot -save aoqplot_stddev_%s_on_moon StandardDeviation %s" % (on_moon_obsid, on_moon_ms_name)
      print cmd
      os.system(cmd)
      
      cmd = "aoqplot -save aoqplot_stddev_%s_off_moon StandardDeviation %s" % (off_moon_obsid, off_moon_ms_name)
      print cmd
      os.system(cmd)
   
   
   
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

         #run andres ionpeel, changed -t 1 from 8, since doing 8s timeres in cotter for EoR2 
         cmd='ionpeel -t 1 -minuv 30 %s %s %s ' % (peeled_ms_name,clustered_model_name,solutions_name)
         print cmd
         os.system(cmd) 
   
      else:
         #run andres calibrate , changed minuv to 30 to be more in line with rts processing (I think)
         cmd='calibrate -t 1 -minuv 30  %s %s %s ' % (model_string,concat_vis_name,solutions_name)
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
parser.add_option('--manta_ray',action='store_true',dest='manta_ray',default=True,help='set if manta_ray used for downloading (need to unzip and chgcentre) [default=%default]') 
parser.add_option('--flag_ants',type='string', dest='flag_ants',default='',help='List of antennas (space separated) to flag after cottering (andre indexing as for rfigui etc)  e.g. --flag_ants="56,60" [default=%default]')


(options, args) = parser.parse_args()

obsid = args[0]
if (options.track_off_moon):
   track_off_moon_string= args[1]
else:
   track_off_moon_string=' '


selfcal_concat_ms(obsid,track_off_moon_string,options)

