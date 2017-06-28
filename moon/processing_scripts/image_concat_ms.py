#! /usr/bin/env python

#script to import rts uvfits files to ms and image with wsclean
#option to track the moon.

def image_concat_ms(obsid_string,list_index,track_off_moon_string,options):

   obsid_list=obsid_string.split(',')
   obsid=obsid_list[int(list_index)].strip()
   print obsid

   list_index=int(list_index)
   if (options.track_off_moon):
      track_off_moon_list=track_off_moon_string.split(',')
      track_off_moon_paired_obsid=track_off_moon_list[int(float(list_index)*3)].strip()
      track_off_moon_new_RA=track_off_moon_list[int(list_index*3+1)].strip()
      track_off_moon_new_DEC=track_off_moon_list[int(list_index*3+2)].strip()
      print "obs paired with %s centering at RA:%s DEC:%s" % (track_off_moon_paired_obsid,track_off_moon_new_RA,track_off_moon_new_DEC)
 
   data_dir='%sdata/%s/' % (mwa_dir,obsid)
   
   imsize=options.imsize
   #wsclean, beam correct
   if (options.tagname):
      tagname=options.tagname
   else:
      tagname='uvdump'

   if (options.corrected_data):
      data_column=' CORRECTED_DATA '
   else:
      data_column=' DATA '
      
   if (options.pol):
      pol_string = ' -pol '+ options.pol
      if (options.pol == 'xx,yy,xy,yx' or options.pol == 'xx,yy,yx,xy' or options.pol == 'yy,xx,xy,yx' or options.pol == 'yy,xx,yx,xy' or options.pol == 'xx,xy,yx,yy' or options.pol == 'yy,xy,yx,xx' or options.pol == 'xx,yx,xy,yy' or options.pol == 'yy,yx,xy,xx'):
         pol_string = ' -pol '+ options.pol+ ' -joinpolarizations '
   else:
      pol_string = '  '

   wsclean_options_string=options.wsclean_options

   if (options.multi):
      multi_string=options.multi
      obsid='multi'
   else:
      multi_string='' 
      
   if (options.chgcentre):
      chgcentre_string=options.chgcentre
   else:
      chgcentre_string=' ' 

   if (options.minw):
      minw_string=options.minw
   else:
      minw_string=' '
      
   metafits_filename=obsid+'.metafits'

   if not (options.concat6):
      #concat_vis_base='%s_%s_concat' % (obsid,tagname)
      if (options.cotter):
         concat_vis_base='%s_cotter_%s' % (obsid,tagname)
      else:
         #need to image the transformed ms
         concat_vis_base='%s_%s_concat_transform' % (obsid,tagname)
      concat_image_base=concat_vis_base
      
      if (options.selfcal):
         concat_image_base+='_selfcal_%s' % (options.selfcal)      
      if (options.chgcentre):
         concat_image_base+='_newcentre' 
      if (options.minw):
         concat_image_base+='_minw' 
      if (options.track_moon):
         concat_vis_base+='_trackmoon'
         concat_image_base+='_trackmoon'
      if (options.track_off_moon):
         concat_vis_base+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid
         concat_image_base+='_track_off_moon_paired_%s' % track_off_moon_paired_obsid    
      if (options.ionpeeled):
         concat_image_base+='_peeled'
         concat_vis_base+='_peeled'
         
      concat_vis_name=data_dir+concat_vis_base+'.ms'

   #This is now done earlier at the cotter stage 
   #if tracking moon, use chgcentre
   #if (options.track_moon and not options.concat6):
   #   #get the date and time of the observation from the metafits file
   #   try:
   #      HDU_list = pyfits.open(metafits_filename)
   #   except IOError, err:
   #      'Cannot open metadata file %s\n' % str(options.input_file)
   #   header=HDU_list[0].header 
   #   date = (header)['DATESTRT']
   #   print date
   #   #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
   #   new_date=string.replace(string.replace(date,'-','/'),'T',' ')
   #   print new_date
   #   #find position of moon
   #   cmd="print_src.py --date='"+new_date+"' > src_file.txt"
   #   print cmd
   #   os.system(cmd)
   #   #read src_file.txt to get Moon ra and dec 
   #   with open('src_file.txt', "r") as infile:
   #      lines=infile.readlines()
   #      moon_ra=lines[6].split()[3].split(',')[0]
   #      moon_dec=lines[6].split()[6].split(',')[0]
   #      print moon_ra
   #      print moon_dec
   #   #get ra and dec in right format for chg_centre 
   #   new_moon_dec=string.replace(moon_dec,":",".")
   #   #make a copy of the ms to change phase centre of
   #   ms_copy_name=concat_vis_base+"_chgcentre.ms"
   #   cmd="rm -rf "+ms_copy_name
   #   print cmd
   #   os.system(cmd)
   #   cmd="cp -r "+concat_vis_name+ " "+ms_copy_name
   #   print cmd
   #   os.system(cmd)
   #
   #   concat_vis_name = ms_copy_name
   #   #change the phase centre of the ms to the Moon position
   #   cmd="aprun -n 1 chgcentre "+concat_vis_name+" "+moon_ra+" "+new_moon_dec
   #   print cmd
   #   os.system(cmd)
   #      
   #   #clean up?

   #   #if not tracking moon, just image the original ms.

   if (options.chgcentre and not options.multi):
      #make a copy of the ms to change phase centre of
      ms_newcentre_base=concat_vis_base+"_newcentre"
      ms_newcentre_name=ms_newcentre_base+".ms"
      if (not options.selfcal):
         cmd="rm -rf "+ms_newcentre_name
         print cmd
         os.system(cmd)
         cmd="cp -r "+concat_vis_name+ " "+ms_newcentre_name
         print cmd
         os.system(cmd)

         #change the phase centre of the ms to the Moon position
         cmd="aprun -n 1 chgcentre %s %s " % (ms_newcentre_name, chgcentre_string)
         print cmd
         os.system(cmd)
         
      concat_vis_base=ms_newcentre_base
         
      if (options.minw):
         ms_copy_base=ms_newcentre_base+"_minw"
         ms_copy_name=ms_copy_base+".ms"
         if (not options.selfcal):
            cmd="rm -rf "+ms_copy_name
            print cmd
            os.system(cmd)
            cmd="cp -r "+ms_newcentre_name+ " "+ms_copy_name
            print cmd
            os.system(cmd)
            #find the minw position for the specified ms
            minw_filename=ms_copy_base+'.txt'
            if (minw_string=='self'):
                cmd="aprun -n 1 chgcentre %s > %s " % (ms_copy_name,minw_filename)    
            else:
                cmd="aprun -n 1 chgcentre %s > %s " % (minw_string,minw_filename)
            print cmd
            os.system(cmd)
            with open(minw_filename) as f:
                minw_value=f.readlines()[8].strip()
                print minw_value
            #change the phase centre of the minw position and shiftback
            cmd="aprun -n 1 chgcentre -shiftback %s %s " % (ms_copy_name,minw_value)
            print cmd
            os.system(cmd)
         
         concat_vis_base=ms_copy_base
         
      concat_vis_name = concat_vis_base+".ms"
   
   if not (options.concat6):
      #remove any old images before re-imaging:
      cmd='rm -f '+concat_image_base+'*.fits'
      print cmd
      os.system(cmd)
      cmd='rm -f '+concat_image_base+'*.tmp'
      print cmd
      os.system(cmd)

      if (options.multi):
         if (options.chgcentre):
            #make a copy of the mss to change phase centre (if this is not a selfcal)
            new_multi_string=""
            for ms_file in multi_string.strip().split(' '):
               ms_copy_base_full=ms_file.split('.ms')[0]
               ms_copy_base=ms_copy_base_full.split('/')[-1]
               ms_copy_name=ms_copy_base+'_newcentre.ms'
               if (not options.selfcal):
                  cmd="rm -rf "+ms_copy_name
                  print cmd
                  os.system(cmd)
                  cmd="cp -r "+ms_file+ " "+ms_copy_name
                  print cmd
                  os.system(cmd)
                  #change the phase centre of the ms to the Moon position
                  cmd="aprun -n 1 chgcentre %s %s " % (ms_copy_name,chgcentre_string)
                  print cmd
                  os.system(cmd)
               new_multi_string+=(ms_copy_name + " ")
               if not (options.minw):
                  multi_string=new_multi_string
         
         if (options.minw):
            new_multi_string=""
            for ms_file in multi_string.strip().split(' '):
               ms_copy_base_full=ms_file.split('.ms')[0]
               ms_copy_base=ms_copy_base_full.split('/')[-1]
               if (options.chgcentre):
                  ms_copy_base+='_newcentre'
               ms_copy_base+='_minw'   
               ms_copy_name=ms_copy_base+'.ms'
               if (not options.selfcal):
                  cmd="rm -rf "+ms_copy_name
                  print cmd
                  os.system(cmd)
                  cmd="cp -r "+ms_file+ " "+ms_copy_name
                  print cmd
                  os.system(cmd)
                  #find the minw position for the specified ms
                  minw_filename=ms_file.split('.ms')[0]+'_minw.txt'
                  if (minw_string=='self'):
                      cmd="aprun -n 1 chgcentre %s > %s " % (ms_file,minw_filename)    
                  else:
                      cmd="aprun -n 1 chgcentre %s > %s " % (minw_string,minw_filename)
                  print cmd
                  os.system(cmd)
                  with open(minw_filename) as f:
                     minw_value=f.readlines()[8].strip()
                     print minw_value
                  #change the phase centre of the minw position and shiftback
                  cmd="aprun -n 1 chgcentre -shiftback %s %s " % (ms_copy_name,minw_value)
                  print cmd
                  os.system(cmd)
               new_multi_string+=(ms_copy_name + " ")
               multi_string=new_multi_string
         
         cmd='aprun -n 1 wsclean -name '+ concat_image_base+' -size '+imsize+' '+imsize+ ' ' +wsclean_options_string+' '+ pol_string + ' ' + multi_string
         print cmd
         os.system(cmd)


      #image with wsclean
      else: 
         cmd='aprun -n 1 wsclean -name '+ concat_image_base+' -size '+imsize+' '+imsize+ ' ' +wsclean_options_string+' '+ pol_string + ' ' + concat_vis_name
         print cmd
         os.system(cmd)
		 #delete all the unwanted files
		 #cmd='rm -f '+concat_image_base+'*-residual.fits'
		 #print cmd
		 #os.system(cmd)
		 #cmd='rm -f '+concat_image_base+'*-dirty.fits'
		 #print cmd
		 #os.system(cmd)
		 #cmd='rm -f '+concat_image_base+'*-model.fits'
		 #print cmd
		 #os.system(cmd)
		 #cmd='rm -f '+concat_image_base+'*-psf.fits'
		 #print cmd
		 #os.system(cmd)
      
      #only generate beam if any one of the beam files doesn't already exist and you are doing pbcorrecting at this stage (not mosaicking later
      if not (options.no_pbcorr):
         if not os.path.exists(concat_image_base+'_beam-xxi.fits') or not os.path.exists(concat_image_base+'_beam-xxr.fits') or not os.path.exists(concat_image_base+'_beam-xyi.fits') or not os.path.exists(concat_image_base+'_beam-xyr.fits') or not os.path.exists(concat_image_base+'_beam-yxi.fits')  or not os.path.exists(concat_image_base+'_beam-yxr.fits')  or not os.path.exists(concat_image_base+'_beam-yyi.fits') or not os.path.exists(concat_image_base+'_beam-yyr.fits'):
            cmd='aprun -n 1 beam -2014i -proto '+ concat_image_base+'-XX-image.fits -name '+concat_image_base+'_beam -ms '+concat_vis_name+ ' -m ' +obsid+ '.metafits'
            os.system(cmd)
         cmd='aprun -n 1 pbcorrect '+concat_image_base+ ' image.fits '+concat_image_base+'_beam '+ concat_image_base
         print cmd
         os.system(cmd)

   if (options.concat6):
      for band in range(1,5):

         concat_vis_base='%s_%s_concat_%sof4' % (obsid,tagname,str(band))
         concat_vis_name=concat_vis_base+'.ms'
         concat_image_base='%s_%s_concat_%sof4' % (obsid,tagname,str(band))

         #remove any old images before re-imaging:
         cmd='rm -f '+concat_image_base+'*.fits'
         print cmd
         os.system(cmd)

         #image with wsclean
         cmd='aprun -n 1 wsclean -name '+ concat_image_base+' -size '+imsize+' '+imsize+ '  -scale .015 -niter 1000 -weight uniform -threshold 0.05 -mgain 0.80 -datacolumn '+data_column+' -absmem 32 '+ pol_string + ' ' + concat_vis_name
         print cmd
         os.system(cmd)
         #delete all the unwanted files
         cmd='rm -f '+concat_image_base+'*-residual.fits'
         print cmd
         os.system(cmd)
         cmd='rm -f '+concat_image_base+'*-dirty.fits'
         print cmd
         os.system(cmd)
         cmd='rm -f '+concat_image_base+'*-model.fits'
         print cmd
         os.system(cmd)
         cmd='rm -f '+concat_image_base+'*-psf.fits'
         print cmd
         os.system(cmd)

         #only generate beam if any one of the beam files doesn't already exist and you are doing pbcorrecting at this stage (not mosaicking later
         if not (options.no_pbcorr):
            if not os.path.exists(concat_image_base+'_beam-xxi.fits') or not os.path.exists(concat_image_base+'_beam-xxr.fits') or not os.path.exists(concat_image_base+'_beam-xyi.fits') or not os.path.exists(concat_image_base+'_beam-xyr.fits') or not os.path.exists(concat_image_base+'_beam-yxi.fits')  or not os.path.exists(concat_image_base+'_beam-yxr.fits')  or not os.path.exists(concat_image_base+'_beam-yyi.fits') or not os.path.exists(concat_image_base+'_beam-yyr.fits'):
               cmd='aprun -n 1 beam -2014i -proto '+ concat_image_base+'-XX-image.fits -name '+concat_image_base+'_beam -ms '+concat_vis_name+ ' -m ' +obsid+ '.metafits'
               os.system(cmd)
            cmd='aprun -n 1 pbcorrect '+concat_image_base+ ' image.fits '+concat_image_base+'_beam '+ concat_image_base
            print cmd
            os.system(cmd)




##
import sys,os
from optparse import OptionParser,OptionGroup

#If --multi is true, then put in chunk_no as obs_id, that forms base of image name
usage = 'Usage: image_concat_ms.py [obsid_string list_index'

parser = OptionParser(usage=usage)

parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--track_off_moon',action='store_true',dest='track_off_moon',default=False,help='Track the Moons position on a previous night. Provide the name of a text file with two columns (obsid_from_previous_night cotter_ms_phase_centre  [default=%default]')
parser.add_option('--imsize',dest='imsize',default='1024',help='Image size string in pixels for wsclean e.g. imsize="4096" [default=%default]')
parser.add_option('--no_pbcorr',action='store_true',dest='no_pbcorr',default=False,help='Dont generate beam and do pbcorrection at this stage (choose this if doing mosaicking later on). [default=%default]')
parser.add_option('--tagname',type='string', dest='tagname',default='',help='Tagname for ms files to image  e.g. --tagname="" [default=%default]')
parser.add_option('--corrected_data',action='store_true',dest='corrected_data',default=False,help='If doing self-cal, use the CORRECTED_DATA column instead of the DATA column. [default=%default]')
parser.add_option('--pol',dest='pol',default='xx,yy',help='Polarisations to image with wsclean e.g. pol="xx,yy,xy,yx" [default=%default]')
parser.add_option('--concat6',action='store_true',dest='concat6',default=False,help='Only image concatenated ms of 7.68 MHz (6 coarse chans) [default=%default]')
parser.add_option('--chgcentre',type='string',dest='chgcentre',default=None,help='Use chgcentre to change phase centre to desired phase centre or shift to min w and shiftback e.g. --chgcentre=" -minw -shiftback " or --chgcentre=" 00h00m00.0s 00d00m00.0s "  [default=%default]')
parser.add_option('--wsclean_options',dest='wsclean_options',default='',help='wsclean options [default=%default]')
parser.add_option('--multi',type='string', dest='multi',default='',help='list of ms to image at the same time  e.g. --multi="" [default=%default]')
parser.add_option('--selfcal',dest='selfcal',default=None,help='Set if this is a subsequent imaging run after selfcal e.g. --selfcal=1 [default=%default]')
parser.add_option('--ionpeeled',action='store_true',dest='ionpeeled',default=False,help='Set if this is a subsequent imaging run after ionpeeling e.g. --ionpeeled [default=%default]')
parser.add_option('--minw',type='string',dest='minw',default=None,help='Shift to minw position of whatever ms is central in the chunk and then shiftback (must start at same phase centre which eor obs do e.g. --minw="12345678.ms"  [default=%default]')
parser.add_option('--cotter',action='store_true',dest='cotter',default=False,help='Use an ms from cotter, not imported from RTS e.g. --cotter   [default=%default]')


(options, args) = parser.parse_args()

obsid_string = args[0]
list_index = args[1]
if (options.track_off_moon):
   track_off_moon_string= args[2]
else:
   track_off_moon_string=' '


mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

image_concat_ms(obsid_string,list_index,track_off_moon_string,options)


