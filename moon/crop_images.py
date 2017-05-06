#script to crop  Moon images to sacve disk space

from astropy.io import fits
from astropy import wcs
import numpy as np
import scipy.ndimage as ndimage
import os
import os.path

#set if processing on Galaxy
galaxy=False

#final_image_size=520
#band_centre_chans=[69,93,121,145,169]
band_centre_chans=[69]

for centre_chan in band_centre_chans:
   #read the info files to find out how many observations there are and what the on and off-moon obsids are:
   if (galaxy):
     on_moon_filename="20150926_moon1/%s/20150926_moon_%s_test.txt" % (str(centre_chan),str(centre_chan))
     off_moon_filename="20150929_off_moon1/%s/20150929_off_moon1_%s_test.txt" % (str(centre_chan),str(centre_chan))   
   else:
     on_moon_filename="20150926_moon_%s_test.txt" % (str(centre_chan))
     off_moon_filename="20150929_off_moon1_%s_test.txt" % (str(centre_chan))
   
   print on_moon_filename
   
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]
   
   with open(on_moon_filename,'r') as on_moon_file:
      on_moon_lines = on_moon_file.readlines()
   for line in on_moon_lines:
      on_moon_obsid_list.append(line.strip())
   
   with open(off_moon_filename,'r') as off_moon_file:
      off_moon_lines = off_moon_file.readlines()
   for line in off_moon_lines:
      off_moon_obsid_list.append(line.strip())

   n_obs=len(on_moon_obsid_list)
   n_chans=24 
   
   #need to do all this stuff for each frequency channel and each obsid
   #sort out filenames 
   for chan in range(n_chans):
      for obsid_index in range(n_obs):
         on_moon_obsid=on_moon_obsid_list[obsid_index]
         off_moon_obsid=off_moon_obsid_list[obsid_index] 
         chan_string='%.04d' % chan
         if (galaxy):
            large_on_moon_image_name="20150926_moon1/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            cropped_on_moon_outname="20150926_moon1/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            large_off_moon_image_name="20150929_off_moon1/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            cropped_off_moon_outname="20150929_off_moon1/%s/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I_cropped.fits" % (str(centre_chan),off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            large_psf_image_name="20150926_moon1/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
            cropped_psf_outname="20150926_moon1/%s/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf_cropped.fits" % (str(centre_chan),on_moon_obsid,str(centre_chan),chan_string)
         else:
            large_on_moon_image_name="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            cropped_on_moon_outname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s_dirty_applied-I_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            large_off_moon_image_name="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            cropped_off_moon_outname="images/%s_cotter_20150929_moon_%s_track_off_moon_paired_%s_peeled-%s_dirty_applied-I_cropped.fits" % (off_moon_obsid,str(centre_chan),on_moon_obsid,chan_string)
            large_psf_image_name="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf.fits" % (on_moon_obsid,str(centre_chan),chan_string)
            cropped_psf_outname="images/%s_cotter_20150926_moon_%s_trackmoon_peeled-%s-psf_cropped.fits" % (on_moon_obsid,str(centre_chan),chan_string)

         #On Moon
         if os.path.isfile(large_on_moon_image_name) and os.access(large_on_moon_image_name, os.R_OK):
             f = fits.open(large_on_moon_image_name)
         else:
            print "Either file %s is missing or is not readable" % large_on_moon_image_name
            continue
         
         cdelt1=f[0].header['CDELT1']
         print cdelt1 
         w = wcs.WCS(f[0].header)
         newf = fits.PrimaryHDU()

         #69, 93 and 121 are all 5120 x 5120 pix
         if (centre_chan <= 121): 
            xstart_moon,xend_moon,ystart_moon,yend_moon=2300,2820,2300,2820
         #145,169 are 3840 x 3840
         else:
            xstart_moon,xend_moon,ystart_moon,yend_moon=1660,2180,1660,2180
         
         image=f[0].data[:,:,xstart_moon:xend_moon,ystart_moon:yend_moon]

         newf.data = image
         newf.header = f[0].header
         newf.header.update(w[:,:,xstart_moon:xend_moon,ystart_moon:yend_moon].to_header())
         #newf.header['CDELT1']=f[0].header['CDELT1']
         #print newf.header['CDELT1'] 
 
         fits.writeto(cropped_on_moon_outname,newf.data,clobber=True)
         fits.update(cropped_on_moon_outname, newf.data, newf.header)
         f.close()
         
         #Off Moon
         if os.path.isfile(large_off_moon_image_name) and os.access(large_off_moon_image_name, os.R_OK):
             f = fits.open(large_off_moon_image_name)
         else:
            print "Either file %s is missing or is not readable" % large_off_moon_image_name
            continue
         w = wcs.WCS(f[0].header)
         newf = fits.PrimaryHDU()

         #69, 93 and 121 are all 5120 x 5120 pix
         if (centre_chan <= 121): 
            xstart_moon,xend_moon,ystart_moon,yend_moon=2300,2820,2300,2820
         #145,169 are 3840 x 3840
         else:
            xstart_moon,xend_moon,ystart_moon,yend_moon=1660,2180,1660,2180
         
         image=f[0].data[:,:,xstart_moon:xend_moon,ystart_moon:yend_moon]

         newf.data = image
         newf.header = f[0].header
         newf.header.update(w[:,:,xstart_moon:xend_moon,ystart_moon:yend_moon].to_header())
 
         fits.writeto(cropped_off_moon_outname,newf.data,clobber=True)
         fits.update(cropped_off_moon_outname, newf.data, newf.header)
         f.close()
         
         #PSF
         if os.path.isfile(large_psf_image_name) and os.access(large_psf_image_name, os.R_OK):
             f = fits.open(large_psf_image_name)
         else:
            print "Either file %s is missing or is not readable" % large_psf_image_name
            continue
         w = wcs.WCS(f[0].header)
         newf = fits.PrimaryHDU()
         if (centre_chan <= 121): 
            xstart_psf,xend_psf,ystart_psf,yend_psf=2299,2819,2301,2821
         #145,169 are 3840 x 3840
         else:
            xstart_psf,xend_psf,ystart_psf,yend_psf=1659,2179,1659,2179

         image=f[0].data[:,:,xstart_psf:xend_psf,ystart_psf:yend_psf]

         newf.data = image
         newf.header = f[0].header
         newf.header.update(w[:,:,xstart_psf:xend_psf,ystart_psf:yend_psf].to_header())
 
         fits.writeto(cropped_psf_outname,newf.data,clobber=True)
         fits.update(cropped_psf_outname, newf.data, newf.header)
         f.close()
            
            
            
            
            
            
