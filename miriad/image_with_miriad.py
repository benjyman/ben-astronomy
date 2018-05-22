#!/usr/bin/env python 
import numpy as np
import os
import cmd


#export the fits files as jpg or png with external grid
def export_image_ds9(options):
   basename=options.basename
   subarray_list=['full','centre','north','south','east','west']
   #subarray_list=['full']
   n_image_time_chunks = 60
   for subarray in subarray_list:
      for i in range(0,n_image_time_chunks):
         image_basename="%s_%s_%02d" % (basename,subarray,i)
         #cmd = "ds9 -fits %s_xx_restor.fits -scale limits 1 100 -cmap invert yes -colorbar yes -grid yes -grid axes type exterior -export jpg %s_xx_restor.jpg -exit " % (image_basename,image_basename)
         cmd = "ds9 %s_xx_restor.fits -invert -colorbar yes -view buttons no -view panner no -view magnifier no -view info no -scale limits 0.0 100 -zoom 1 -width 512 -height 512 -grid  -saveimage png %s_xx_restor.png -quit" % (image_basename,image_basename)
         print cmd
         os.system(cmd)

         cmd = "ds9 %s_yy_restor.fits -invert -colorbar yes -view buttons no -view panner no -view magnifier no -view info no -scale limits 0.0 100 -zoom 1 -width 512 -height 512 -grid  -saveimage png %s_yy_restor.png -quit" % (image_basename,image_basename)
         print cmd
         os.system(cmd)

   
   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_full_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_full_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_full_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_full_yy.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_centre_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_centre_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_centre_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_centre_yy.mp4" 
   print cmd
   os.system(cmd)
   
   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_north_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_north_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_north_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_north_yy.mp4" 
   print cmd
   os.system(cmd)
   
   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_south_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_south_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_south_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_south_yy.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_east_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_east_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_east_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_east_yy.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_west_%02d_xx_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_west_xx.mp4" 
   print cmd
   os.system(cmd)

   cmd = "ffmpeg -framerate 5 -i sun_20180405_aavs1_west_%02d_yy_restor.png -c:v libx264 -r 30 -pix_fmt yuv420p asun_20180405_aavs1_west_yy.mp4" 
   print cmd
   os.system(cmd)      
   #ffmpeg -framerate 5 -i sun_aavs1_${stokes}_%d.png -c:v libx264 -r 30 -pix_fmt yuv420p sun_aavs1_${stokes}.mp4

#script to do imaging with miriad incl subarrays     
def image_with_miriad(options):
   
   basename=options.basename
   #subarray_list=['full','centre','north','south','east','west']
   subarray_list=['full']
   n_flag_chunks = 8
   n_image_time_chunks = 60
   date = "18apr05"
   
   for subarray in subarray_list:
            
      #read in the uvdata
      uvfits_name="/md0/AAVS-1/data/%s.uvfits" % (basename)
      uvdata="/md0/AAVS-1/data/%s.uv" % (basename)
      
      if options.read_uvdata:
         cmd= "fits in=%s out=%s op=uvin" % (uvfits_name,uvdata)
         print cmd
         os.system(cmd)
      
      #sort out the subarrays
      if (subarray=='full'):
         flagged_ant_filename=''
         print "imaging with full array"
      elif (subarray=='centre'):
         flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedCentralArea_FlaggedAntID.txt"
      elif (subarray=='north'):
         flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedNorthArea_FlaggedAntID.txt"
      elif (subarray=='south'):
         flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedSouthArea_FlaggedAntID.txt"
      elif (subarray=='east'):
         flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedEastArea_FlaggedAntID.txt"
      elif (subarray=='west'):
         flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedWestArea_FlaggedAntID.txt"
      else:
         print "no subarray"
      
      flagged_ant_list=[]
      #read the files to get the antenna arrays
      if subarray != 'full':
         
         new_uvdata_name='/md0/AAVS-1/data/%s_%s.uv' % (basename,subarray)
         
         #copy uvdata
         if options.copy_uvdata_subarrays:
            cmd = "cp -r %s %s" % (uvdata,new_uvdata_name)
            print cmd
            os.system(cmd)
            
            for line in open(flagged_ant_filename):
               flagged_ant_list.append(line.strip()) 
         
            #have to split up flagged ant array cause miriad can't handle a big string
            #split into 57 chunks
            split_flagged_ant_array = np.array_split(flagged_ant_list,n_flag_chunks)
            #print flagged_ant_filename
            #print split_flagged_ant_array
      
            for i in range(0,n_flag_chunks):
               list=split_flagged_ant_array[i]
               flagged_ant_string = ','.join(map(str, list)) 
               flagged_ant_string_select="select=antennae\(%s\)" % flagged_ant_string
            
               cmd = "uvflag vis='%s' %s flagval=flag" % (new_uvdata_name,flagged_ant_string_select)
               print cmd
               os.system(cmd)
            
         uvdata = new_uvdata_name
      
      #calibrate the data
      if options.calibrate:
         cmd = "mfcal vis=%s flux=10000,0.150,0 refant=2 interval=1 stokes=xx,yy,xy,yx" % (uvdata)
         print cmd
         os.system(cmd)
         
         
      #make an image for every two minutes
      
      #use varplt to get the times
      #mirFix the images to put correct keywords in header for Slant Orthographic Projection
      cmd = "if [ ! -s lst.txt ] ; then varplt vis=%s yaxis=lst log=lst.txt ; fi " % (uvdata)
      print cmd
      os.system(cmd)
      
      time_list=[]
      lst_list=[]
      for line in open('lst.txt'):
         time_list.append(line.split()[1].strip()) 
         lst_list.append(line.split()[2].strip()) 
         
      #two hours of data - image every 2 mins - 60 images
   
      split_time_list = np.array_split(time_list[4:],n_image_time_chunks)
      split_lst_list = np.array_split(lst_list[4:],n_image_time_chunks)
   
      for i in range(0,n_image_time_chunks):
         time_1 = split_time_list[i][0]
         time_2 = split_time_list[i][-1]
         
         lst_index=int(np.round(len(split_time_list[i])/2))
         lst_hrs = split_lst_list[i][lst_index]
         
         t1 = "%s:%s.0" % (date,time_1)
         t2 = "%s:%s.0" % (date,time_2)
      
         #t1="18apr05:04:00:00"
         #t2="18apr05:04:02:00"
      
         select_time_string="select=time\(%s,%s\)" % (t1,t2)
         
         image_basename="%s_%s_%02d" % (basename,subarray,i)
      
         cmd = "invert vis=%s map=%s_xx.map beam=%s.beam options=double imsize=512 stokes=xx robust=-0.5 cell=1800 %s " % (uvdata,image_basename,image_basename,select_time_string)
         print cmd
         os.system(cmd)
         cmd = "clean map=%s_xx.map beam=%s.beam out=%s_xx.model niters=2000" % (image_basename,image_basename,image_basename)
         print cmd
         os.system(cmd)
         cmd = "restor model=%s_xx.model  beam=%s.beam map=%s_xx.map out=%s_xx.restor " % (image_basename,image_basename,image_basename,image_basename)
         print cmd
         os.system(cmd)
      
         cmd = "invert vis=%s map=%s_yy.map options=double imsize=512 stokes=yy robust=-0.5 cell=1800 %s " % (uvdata,image_basename,select_time_string)
         print cmd
         os.system(cmd)
         cmd = "clean map=%s_yy.map beam=%s.beam out=%s_yy.model niters=2000" % (image_basename,image_basename,image_basename)
         print cmd
         os.system(cmd)
         cmd = "restor model=%s_yy.model  beam=%s.beam map=%s_yy.map out=%s_yy.restor " % (image_basename,image_basename,image_basename,image_basename)
         print cmd
         os.system(cmd)   
      
         #cmd = "lst_hrs=`head -$(( (nlines-4)/2 +5 )) lst.txt | tail -1 | tr -s " "| cut -f 4 -d " "` "
         #print cmd
         #os.system(cmd)
      
         cmd = "/data/code/git/ben-astronomy/miriad/mirFixhdr.sh %s_xx.restor %s" % (image_basename,lst_hrs)
         print cmd
         os.system(cmd)
      
         cmd = "/data/code/git/ben-astronomy/miriad/mirFixhdr.sh %s_yy.restor %s" % (image_basename,lst_hrs)
         print cmd
         os.system(cmd)
      
         #export to fits
         cmd= "fits in=%s_xx.restor out=%s_xx_restor.fits op=xyout" % (image_basename,image_basename)
         print cmd
         os.system(cmd)
         

         cmd= "fits in=%s_yy.restor out=%s_yy_restor.fits op=xyout" % (image_basename,image_basename)
         print cmd
         os.system(cmd)
         
         
import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: image_with_miriad.py [options]'

parser = OptionParser(usage=usage)


parser.add_option('--basename',type='string', dest='basename',default='sun',help=' e.g. --basename="sun_01" [default=%default]')
parser.add_option('--read_uvdata',action='store_true',dest='read_uvdata',default=False,help='Read in the uvdata  [default=%default]')
parser.add_option('--copy_uvdata_subarrays',action='store_true',dest='copy_uvdata_subarrays',default=False,help='Make new copies of data for subarray flagging  [default=%default]')
parser.add_option('--calibrate',action='store_true',dest='calibrate',default=False,help='Do calibration  [default=%default]')
parser.add_option('--export_only',action='store_true',dest='export_only',default=False,help='Just export fits files as images  [default=%default]')



(options, args) = parser.parse_args()

if options.export_only:
   export_image_ds9(options)
else:
   image_with_miriad(options)

