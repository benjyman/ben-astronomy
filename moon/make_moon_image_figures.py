#Script to make moon image figures for papers and presentations
import pyfits
import pylab as py
import numpy as np
 
def make_images(options):
   stokes=options.stokes
   image_dir=options.image_dir
   image_size=[1024,1024]
   #indexes for cropping
   xstart,xend,ystart,yend=437,587,437,587
   centre_chan=options.centre_chan
   #sub_chan=options.sub_chan
   min_val=float(options.min_val)
   max_val=float(options.max_val)
   on_moon_obsid_infile=options.on_moon_obsid_infile
   off_moon_obsid_infile=options.off_moon_obsid_infile
   
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]  
   #n_obs = sum(1 for line in open(on_moon_obsid_infile))
   
   for line in open(on_moon_obsid_infile):
       on_moon_obsid_list.append(line.strip()) 
     
   for line in open(off_moon_obsid_infile):
       off_moon_obsid_list.append(line.strip()) 
 
   #Do all this stuff for each of the 24 sub-bands
   for sub_chan in range(0,24):
       
      chan_string='%.04d' % float(sub_chan)
      
      average_moon_image_fitsname="average_moon_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_off_moon_image_fitsname="average_off_moon_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_difference_image_fitsname="average_difference_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_difference_jyperpix_image_fitsname="average_difference_jyperpix_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_rfi_modelled_image_fitsname="average_rfi_modelled_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_moon_modelled_image_fitsname="average_moon_modelled_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_residual_modelled_image_fitsname="average_residual_modelled_image_chan_%s_%s_%s.fits" % (centre_chan,chan_string,stokes)
      average_moon_image_figname="average_moon_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_off_moon_image_figname="average_off_moon_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_difference_image_figname="average_difference_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_difference_jyperpix_image_figname="average_difference_jyperpix_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_rfi_modelled_image_figname="average_rfi_modelled_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_moon_modelled_image_figname="average_moon_modelled_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_residual_modelled_image_figname="average_residual_modelled_image_chan_%s_%s_%s.png" % (centre_chan,chan_string,stokes)
      average_moon_image_array=np.zeros(image_size)
      average_off_moon_image_array=np.zeros(image_size)
      average_difference_image_array=np.zeros(image_size)
      average_difference_jyperpix_image_array=np.zeros(image_size)
      average_rfi_modelled_image_array=np.zeros(image_size)
      average_moon_modelled_image_array=np.zeros(image_size)
      average_residual_modelled_image_array=np.zeros(image_size)
      
      difference_image_counter=0
      difference_jyperpix_image_counter=0
      moon_image_counter=0
      off_moon_image_counter=0
      modelled_moon_image_counter=0
      modelled_rfi_image_counter=0
      modelled_residual_image_counter=0
      
      for obsid_index,on_moon_obsid in enumerate(on_moon_obsid_list):
         off_moon_obsid=off_moon_obsid_list[obsid_index]
         moon_fitsname="%s/moon_zoom_%s_%s-%s-%s.fits" % (image_dir,on_moon_obsid,str(centre_chan),chan_string,stokes)
         off_moon_fitsname="%s/off_moon_zoom_%s_paired_with_%s_%s-%s-%s.fits" % (image_dir,off_moon_obsid,on_moon_obsid,str(centre_chan),chan_string,stokes)
         moon_difference_fitsname="%s/difference_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
         moon_difference_jyperpix_fitsname="%s/difference_%s_%s_on_off_moon_%s-%s-%s_jyperpix.fits" % (image_dir,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
         modelled_moon_fitsname="%s/moon_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
         modelled_rfi_fitsname="%s/rfi_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
         modelled_residual_fitsname="%s/residual_modelled_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
      
         on_moon_title="Dirty Moon image obsid %s chan %s %s Stokes %s" % (on_moon_obsid,centre_chan,chan_string,stokes)
         off_moon_title="Dirty off-Moon image obsid %s chan %s %s Stokes %s" % (off_moon_obsid,centre_chan,chan_string,stokes)
         diff_moon_title="Difference image on-moon %s off-moon %s chan %s %s Stokes %s" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         diff_moon_jyperpix_title="Difference image Jy per pix on-moon %s off-moon %s chan %s %s Stokes %s" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_moon_title="Modelled Moon image on-moon %s off-moon %s chan %s %s Stokes %s" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_rfi_title="Modelled RFI image on-moon %s off-moon %s chan %s %s Stokes %s" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_residual_title="Modelled Residual image on-moon %s off-moon %s chan %s %s Stokes %s" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         
         on_moon_figname="Dirty_Moon_image_obsid_%s_chan_%s_%s_stokes%s.png" % (on_moon_obsid,centre_chan,chan_string,stokes)      
         off_moon_figname="Dirty_off_Moon_image_obsid_%s_chan_%s_%s_stokes%s.png" % (off_moon_obsid,centre_chan,chan_string,stokes)  
         diff_moon_figname="Difference_image_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         diff_moon_jyperpix_figname="Difference_image_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s_jyperpix.png" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_moon_figname="Modelled_moon_image_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_rfi_figname="Modelled_rfi_image_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         modelled_residual_figname="Modelled_residual_image_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         
         print moon_fitsname  
         if os.path.isfile(moon_fitsname) and os.access(moon_fitsname, os.R_OK):
            moon_hdulist = pyfits.open(moon_fitsname)
         else:
            print "Either file %s is missing or is not readable" % moon_fitsname
            continue       
         moon_data=moon_hdulist[0].data
         average_moon_image_array+=moon_data
         moon_image_counter+=1
      
         print off_moon_fitsname  
         if os.path.isfile(off_moon_fitsname) and os.access(off_moon_fitsname, os.R_OK):
            off_moon_hdulist = pyfits.open(off_moon_fitsname)
         else:
            print "Either file %s is missing or is not readable" % off_moon_fitsname
            continue       
         off_moon_data=off_moon_hdulist[0].data
         average_off_moon_image_array+=off_moon_data
         off_moon_image_counter+=1
                     
         print moon_difference_fitsname  
         if os.path.isfile(moon_difference_fitsname) and os.access(moon_difference_fitsname, os.R_OK):
            moon_difference_hdulist = pyfits.open(moon_difference_fitsname)
         else:
            print "Either file %s is missing or is not readable" % moon_difference_fitsname 
            continue      
         moon_difference_data=moon_difference_hdulist[0].data
         average_difference_image_array+=moon_difference_data
         difference_image_counter+=1
         
         print moon_difference_jyperpix_fitsname  
         if os.path.isfile(moon_difference_jyperpix_fitsname) and os.access(moon_difference_jyperpix_fitsname, os.R_OK):
            moon_difference_jypepix_hdulist = pyfits.open(moon_difference_jyperpix_fitsname)
         else:
            print "Either file %s is missing or is not readable" % moon_difference_jyperpix_fitsname 
            continue      
         moon_difference_jyperpix_data=moon_difference_jypepix_hdulist[0].data
         average_difference_jyperpix_image_array+=moon_difference_jyperpix_data
         difference_jyperpix_image_counter+=1
         
         print modelled_moon_fitsname  
         if os.path.isfile(modelled_moon_fitsname) and os.access(modelled_moon_fitsname, os.R_OK):
            modelled_moon_hdulist = pyfits.open(modelled_moon_fitsname)
         else:
            print "Either file %s is missing or is not readable" % modelled_moon_fitsname 
            continue      
         modelled_moon_data=modelled_moon_hdulist[0].data
         average_moon_modelled_image_array+=modelled_moon_data
         modelled_moon_image_counter+=1
         
         print modelled_rfi_fitsname  
         if os.path.isfile(modelled_rfi_fitsname) and os.access(modelled_rfi_fitsname, os.R_OK):
            modelled_rfi_hdulist = pyfits.open(modelled_rfi_fitsname)
         else:
            print "Either file %s is missing or is not readable" % modelled_rfi_fitsname 
            continue      
         modelled_rfi_data=modelled_rfi_hdulist[0].data
         average_rfi_modelled_image_array+=modelled_rfi_data
         modelled_rfi_image_counter+=1
      
         print modelled_residual_fitsname  
         if os.path.isfile(modelled_residual_fitsname) and os.access(modelled_residual_fitsname, os.R_OK):
            modelled_residual_hdulist = pyfits.open(modelled_residual_fitsname)
         else:
            print "Either file %s is missing or is not readable" % modelled_residual_fitsname 
            continue      
         modelled_residual_data=modelled_residual_hdulist[0].data
         average_residual_modelled_image_array+=modelled_residual_data
         modelled_residual_image_counter+=1
         
         if options.plot_each_obsid:
            #Plot images:
            py.figure(1)
            py.clf()
            py.title(on_moon_title)
            py.imshow( ( moon_data ), cmap=py.cm.Greys,vmin=min_val, vmax=max_val,origin='lower')
            py.colorbar()
            py.savefig(on_moon_figname)
            py.close()
      
            py.figure(2)
            py.clf()
            py.title(off_moon_title)
            py.imshow( ( off_moon_data ), cmap=py.cm.Greys,vmin=min_val, vmax=max_val,origin='lower')
            py.colorbar()
            py.savefig(off_moon_figname)
            py.close()
                  
            py.figure(3)
            py.clf()
            py.title(diff_moon_title)
            py.imshow( ( moon_difference_data ), cmap=py.cm.Greys,origin='lower')
            py.colorbar()
            py.savefig(diff_moon_figname)
            py.close()
      
            py.figure(4)
            py.clf()
            py.title(modelled_moon_title)
            py.imshow( ( modelled_moon_data ), cmap=py.cm.Greys,origin='lower',vmin=min_val, vmax=max_val)
            py.colorbar()
            py.savefig(modelled_moon_figname)
            py.close()
            
            py.figure(5)
            py.clf()
            py.title(modelled_rfi_title)
            py.imshow( ( modelled_rfi_data ), cmap=py.cm.Greys,origin='lower')
            py.colorbar()
            py.savefig(modelled_rfi_figname)
            py.close()
         
            py.figure(6)
            py.clf()
            py.title(modelled_residual_title)
            py.imshow( ( modelled_residual_data ), cmap=py.cm.Greys,origin='lower')
            py.colorbar()
            py.savefig(modelled_residual_figname)
            py.close()   
      
      #crop the images to make smaller
      
      #form the average Moon image and plot
      average_moon_image=average_moon_image_array/moon_image_counter
      average_moon_image_crop = average_moon_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_moon_image_fitsname,average_moon_image,clobber=True)
      py.figure(1)
      py.clf()
      py.title("Average Moon image from %s obsids chan %s %s" % (str(moon_image_counter),centre_chan,chan_string))
      py.imshow( ( average_moon_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_moon_image_figname)
      py.close()
      
      #form the average Off Moon image and plot
      average_off_moon_image=average_off_moon_image_array/off_moon_image_counter
      average_off_moon_image_crop = average_off_moon_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_off_moon_image_fitsname,average_off_moon_image,clobber=True)
      py.figure(2)
      py.clf()
      py.title("Average Off Moon image from %s obsids chan %s %s" % (str(off_moon_image_counter),centre_chan,chan_string))
      py.imshow( ( average_off_moon_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_off_moon_image_figname)
      py.close()
         
      #form the average difference image and plot
      vmin_difference=-2.0
      vmax_difference=2.0
      average_difference_image=average_difference_image_array/difference_image_counter
      average_difference_image_crop = average_difference_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_difference_image_fitsname,average_difference_image,clobber=True)
      py.figure(3)
      py.clf()
      py.title("Average difference image from %s obsids chan %s %s" % (str(difference_image_counter),centre_chan,chan_string))
      py.imshow( ( average_difference_image_crop ), cmap=py.cm.Greys,origin='lower',vmin=vmin_difference, vmax=vmax_difference)
      py.colorbar()
      py.savefig(average_difference_image_figname)
      py.close()

      #form the average difference image in jy per pix and plot
      vmin_difference_jyperpix=-0.002
      vmax_difference_jyperpix=0.002
      average_difference_jyperpix_image=average_difference_jyperpix_image_array/difference_jyperpix_image_counter
      average_difference_jyperpix_image_crop = average_difference_jyperpix_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_difference_jyperpix_image_fitsname,average_difference_jyperpix_image,clobber=True)
      py.figure(3)
      py.clf()
      py.title("Average difference Jy per pix image from %s obsids chan %s %s" % (str(difference_jyperpix_image_counter),centre_chan,chan_string))
      py.imshow( ( average_difference_jyperpix_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_difference_jyperpix_image_figname)
      py.close()
            
      #form the average modelled moon image and plot
      average_moon_modelled_image=average_moon_modelled_image_array/modelled_moon_image_counter
      average_moon_modelled_image_crop = average_moon_modelled_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_moon_modelled_image_fitsname,average_moon_modelled_image,clobber=True)
      py.figure(4)
      py.clf()
      py.title("Average modelled moon image from %s obsids chan %s %s" % (str(modelled_moon_image_counter),centre_chan,chan_string))
      py.imshow( ( average_moon_modelled_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_moon_modelled_image_figname)
      py.close()
      
      #form the average modelled rfi image and plot
      average_rfi_modelled_image=average_rfi_modelled_image_array/modelled_rfi_image_counter
      average_rfi_modelled_image_crop = average_rfi_modelled_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_rfi_modelled_image_fitsname,average_rfi_modelled_image,clobber=True)
      py.figure(5)
      py.clf()
      py.title("Average modelled rfi image from %s obsids chan %s %s" % (str(modelled_rfi_image_counter),centre_chan,chan_string))
      py.imshow( ( average_rfi_modelled_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_rfi_modelled_image_figname)
      py.close()   
      
      #form the average modelled residual image and plot
      average_residual_modelled_image=average_residual_modelled_image_array/modelled_residual_image_counter
      average_residual_modelled_image_crop = average_residual_modelled_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_residual_modelled_image_fitsname,average_residual_modelled_image,clobber=True)
      py.figure(6)
      py.clf()
      py.title("Average modelled residual image from %s obsids chan %s %s" % (str(modelled_rfi_image_counter),centre_chan,chan_string))
      py.imshow( ( average_residual_modelled_image_crop ), cmap=py.cm.Greys,origin='lower')
      py.colorbar()
      py.savefig(average_residual_modelled_image_figname)
      py.close() 
    


  
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: cotter_moon.py [obsid] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--centre_chan',type='string', dest='centre_chan',default=None,help='Centre channel (e.g.coarse chan number 69, 93, 121, 145 or 169)) to make Moon images. e.g. --centre_chan="69" [default=%default]')
#parser.add_option('--sub_chan',type='string', dest='sub_chan',default=None,help='Sub channel (e.g. chan number within centre chan boundaries 0-23)) to make Moon images. e.g. --sub_chan="12" [default=%default]')
parser.add_option('--min_val',type='string', dest='min_val',default="-13.0",help='Minimum data value in Moon image stretch e.g. --min_val="-13.0"')
parser.add_option('--max_val',type='string', dest='max_val',default="22.0",help='Minimum data value in Moon image stretch e.g. --max_val="22.0"')
parser.add_option('--on_moon_obsid_infile',type='string', dest='on_moon_obsid_infile',default=None,help='Name of file containing the on moon obsids to plot e.g. --on_moon_obsid_infile="/data/moon/20150926_69_moon.txt"')
parser.add_option('--off_moon_obsid_infile',type='string', dest='off_moon_obsid_infile',default=None,help='Name of file containing the off moon obsids to plot e.g. --off_moon_obsid_infile="/data/moon/20150926_69_off_moon.txt"')
parser.add_option('--plot_each_obsid',action='store_true',dest='plot_each_obsid',default=False,help='Make individual images for each obsid, rather than just the average [default=%default]')
parser.add_option('--stokes',type='string', dest='stokes',default='I',help='stokes parameter of Moon images. Can be I,Q,U,V or "linear" (sqrt(Q^2+U^2)). e.g. --stokes="Q" [default=%default]')
parser.add_option('--image_dir',type='string', dest='image_dir',default='/data/moon/2017',help='Directory where the Moon images from model_moon.py are stored. e.g. --image_dir="/md0/moon/2015/stokes_I" [default=%default]')




(options, args) = parser.parse_args()

#freq = args[0]

make_images(options)






   
