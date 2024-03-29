#Script to make moon image figures for papers and presentations
import matplotlib
matplotlib.use('Agg')
import pyfits
import pylab as py
import numpy as np
from astropy.wcs import WCS
import astropy 
import datetime
import copy
import aplpy
import matplotlib.pyplot as plt
from matplotlib.dates import  DateFormatter


def set_shared_ylabel(a, ylabel, labelpad = 0.01):
    """Set a y label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    ylabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""

    f = a[0].get_figure()
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    top = a[0].get_position().y1
    bottom = a[-1].get_position().y0

    # get the coordinates of the left side of the tick labels 
    x0 = 1
    for at in a:
        at.set_ylabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.yaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.inverse_transformed(f.transFigure)
        xt = bboxes.x0
        if xt < x0:
            x0 = xt
    tick_label_left = x0

    # set position of label
    a[-1].set_ylabel(ylabel)
    a[-1].yaxis.set_label_coords(tick_label_left - labelpad,(bottom + top)/2, transform=f.transFigure)

 
def make_images(options):
   stokes=options.stokes
   epoch_ID=options.epoch_ID
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

   vmin_new=-0.005
   vmax_new=0.005
   
   on_moon_obsid_list=[]
   off_moon_obsid_list=[]  
   #n_obs = sum(1 for line in open(on_moon_obsid_infile))
   
   for line in open(on_moon_obsid_infile):
       on_moon_obsid_list.append(line.strip()) 
     
   for line in open(off_moon_obsid_infile):
       off_moon_obsid_list.append(line.strip()) 
 
   #Do all this stuff for each of the 24 sub-bands
   for sub_chan in range(23,24):
       
      chan_string='%.04d' % float(sub_chan)
      
      average_moon_image_fitsname="average_moon_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_off_moon_image_fitsname="average_off_moon_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_difference_image_fitsname="average_difference_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_difference_jyperpix_image_fitsname="average_difference_jyperpix_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_rfi_modelled_image_fitsname="average_rfi_modelled_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_moon_modelled_image_fitsname="average_moon_modelled_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      new_cropped_moon_modelled_image_fitsname="new_moon_modelled_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_residual_modelled_image_fitsname="average_residual_modelled_image_%s_chan_%s_%s_%s.fits" % (epoch_ID,centre_chan,chan_string,stokes)
      average_moon_image_figname="average_moon_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_off_moon_image_figname="average_off_moon_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_difference_image_figname="average_difference_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_difference_jyperpix_image_figname="average_difference_jyperpix_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_rfi_modelled_image_figname="average_rfi_modelled_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_moon_modelled_image_figname="average_moon_modelled_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      new_moon_modelled_image_figname="new_moon_modelled_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
      average_residual_modelled_image_figname="average_residual_modelled_image_%s_chan_%s_%s_%s.png" % (epoch_ID,centre_chan,chan_string,stokes)
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
   
      if (options.new_modelled_only):
         #make the new modelled moon image (i.e. after diffuse RFI subtraction)
         new_moon_modelled_fitsname="%s/new_moon_modelled_%s-%s-%s.fits" % (image_dir,str(centre_chan),chan_string,stokes)
         print new_moon_modelled_fitsname
         if os.path.isfile(new_moon_modelled_fitsname) and os.access(new_moon_modelled_fitsname, os.R_OK):
            new_modelled_moon_hdulist = pyfits.open(new_moon_modelled_fitsname)
         else:
            print "Either file %s is missing or is not readable" % new_moon_modelled_fitsname
            continue
         new_moon_modelled_image=new_modelled_moon_hdulist[0].data
         new_moon_modelled_image_crop=new_moon_modelled_image[ystart:yend,xstart:xend]
         pyfits.writeto(new_cropped_moon_modelled_image_fitsname,new_moon_modelled_image_crop,clobber=True)
         py.figure(6)
         py.clf()
         py.title("New modelled moon image chan %s %s" % (centre_chan,chan_string))
         py.imshow( ( new_moon_modelled_image_crop ), cmap=py.cm.Greys,origin='lower',vmin=vmin_new, vmax=vmax_new)
         py.colorbar()
         py.savefig(new_moon_modelled_image_figname)
         py.close()
         continue    

      for obsid_index,on_moon_obsid in enumerate(on_moon_obsid_list):
         if options.old_labelling:
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
         else:
            off_moon_obsid=off_moon_obsid_list[obsid_index]
            moon_fitsname="%s/moon_zoom_%s_%s_%s-%s-%s.fits" % (image_dir,epoch_ID,on_moon_obsid,str(centre_chan),chan_string,stokes)
            off_moon_fitsname="%s/off_moon_zoom_%s_%s_paired_with_%s_%s-%s-%s.fits" % (image_dir,epoch_ID,off_moon_obsid,on_moon_obsid,str(centre_chan),chan_string,stokes)
            moon_difference_fitsname="%s/difference_%s_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,epoch_ID,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            moon_difference_jyperpix_fitsname="%s/difference_%s_%s_%s_on_off_moon_%s-%s-%s_jyperpix.fits" % (image_dir,epoch_ID,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            modelled_moon_fitsname="%s/moon_modelled_%s_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,epoch_ID,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            modelled_rfi_fitsname="%s/rfi_modelled_%s_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,epoch_ID,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
            modelled_residual_fitsname="%s/residual_modelled_%s_%s_%s_on_off_moon_%s-%s-%s.fits" % (image_dir,epoch_ID,on_moon_obsid,off_moon_obsid,str(centre_chan),chan_string,stokes)
      
            on_moon_title="Dirty Moon image %s obsid %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,centre_chan,chan_string,stokes)
            off_moon_title="Dirty off-Moon image %s obsid %s chan %s %s Stokes %s" % (epoch_ID,off_moon_obsid,centre_chan,chan_string,stokes)
            diff_moon_title="Difference image %s on-moon %s off-moon %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            diff_moon_jyperpix_title="Difference image %s Jy per pix on-moon %s off-moon %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_moon_title="Modelled Moon image %s on-moon %s off-moon %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_rfi_title="Modelled RFI image %s on-moon %s off-moon %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_residual_title="Modelled Residual image %s on-moon %s off-moon %s chan %s %s Stokes %s" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         
            on_moon_figname="Dirty_Moon_image_%s_obsid_%s_chan_%s_%s_stokes%s.png" % (epoch_ID,on_moon_obsid,centre_chan,chan_string,stokes)      
            off_moon_figname="Dirty_off_Moon_image_%s_obsid_%s_chan_%s_%s_stokes%s.png" % (epoch_ID,off_moon_obsid,centre_chan,chan_string,stokes)  
            diff_moon_figname="Difference_image_%s_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            diff_moon_jyperpix_figname="Difference_image_%s_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s_jyperpix.png" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_moon_figname="Modelled_moon_image_%s_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_rfi_figname="Modelled_rfi_image_%s_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
            modelled_residual_figname="Modelled_residual_image_%s_on_moon_%s_off_moon_%s chan_%s_%s_stokes%s.png" % (epoch_ID,on_moon_obsid,on_moon_obsid,centre_chan,chan_string,stokes)
         
         print moon_fitsname  
         if os.path.isfile(moon_fitsname) and os.access(moon_fitsname, os.R_OK):
            moon_hdulist = pyfits.open(moon_fitsname)
         else:
            print "Either file %s is missing or is not readable" % moon_fitsname
            continue 
         
         print "opening image %s" % moon_fitsname
         moon_data=moon_hdulist[0].data
         moon_header=moon_hdulist[0].header
         #print moon_header
         #del moon_header[5]
         #del moon_header[5]
         #del moon_header['history']
         #del moon_header['ctype3']
         #del moon_header['crval3']
         #del moon_header['crpix3']
         #del moon_header['cdelt3']
         #del moon_header['ctype4']
         #del moon_header['crval4']
         #del moon_header['crpix4']
         #del moon_header['cdelt4']
         
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
      print "moon_image_counter is %s" % moon_image_counter
      average_moon_image=average_moon_image_array/moon_image_counter
      average_moon_image_crop = average_moon_image[ystart:yend,xstart:xend]

      wcs = WCS(moon_header)
      #print wcs
      wcs=wcs.dropaxis(2)
      wcs=wcs.dropaxis(2)
      wcs_cropped = wcs[ystart:yend,xstart:xend]
      #print wcs_cropped
      moon_header.update(wcs_cropped.to_header())
      
      #print moon_header
      del moon_header['ctype3']
      del moon_header['cdelt3']
      del moon_header['crval3']
      del moon_header['crpix3']
      del moon_header['cunit3']
      del moon_header['ctype4']
      del moon_header['cdelt4']
      del moon_header['crval4']
      del moon_header['crpix4']
      del moon_header['cunit4']
      #print moon_header
      
      pyfits.writeto(average_moon_image_fitsname,average_moon_image_crop,clobber=True)
      pyfits.update(average_moon_image_fitsname,average_moon_image_crop,clobber=True,header=moon_header)
      
      #py.figure(1)
      #py.clf()
      #py.title("Average Moon image from %s obsids chan %s %s" % (str(moon_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_moon_image_crop ), cmap=py.cm.Greys,origin='lower')
      #py.colorbar()
      #py.savefig(average_moon_image_figname)
      #py.close()
      
      matplotlib.rcParams['xtick.direction'] = 'in'
      matplotlib.rcParams['ytick.direction'] = 'in'

      #gc = aplpy.FITSFigure(average_moon_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_moon_image_figname)
      #print "saved %s with aplpy" % average_moon_image_figname
      #gc.close()
      
      #form the average Off Moon image and plot
      average_off_moon_image=average_off_moon_image_array/off_moon_image_counter
      average_off_moon_image_crop = average_off_moon_image[ystart:yend,xstart:xend]

      
      pyfits.writeto(average_off_moon_image_fitsname,average_off_moon_image_crop,clobber=True)
      pyfits.update(average_off_moon_image_fitsname,average_off_moon_image_crop,clobber=True,header=moon_header)
      
      #gc = aplpy.FITSFigure(average_off_moon_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_off_moon_image_figname)
      #print "saved %s with aplpy" % average_off_moon_image_figname
      #gc.close()
      
      #py.figure(2)
      #py.clf()
      #py.title("Average Off Moon image from %s obsids chan %s %s" % (str(off_moon_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_off_moon_image_crop ), cmap=py.cm.Greys,origin='lower')
      #py.colorbar()
      #py.savefig(average_off_moon_image_figname)
      #py.close()
         
      #form the average difference image and plot
      vmin_difference=-2.0
      vmax_difference=2.0
      average_difference_image=average_difference_image_array/difference_image_counter
      average_difference_image_crop = average_difference_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_difference_image_fitsname,average_difference_image_crop,clobber=True)
      pyfits.update(average_difference_image_fitsname,average_difference_image_crop,clobber=True,header=moon_header)

      #gc = aplpy.FITSFigure(average_difference_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_difference_image_figname)
      #print "saved %s with aplpy" % average_difference_image_figname
      #gc.close()
      
      #py.figure(3)
      #py.clf()
      #py.title("Average difference image from %s obsids chan %s %s" % (str(difference_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_difference_image_crop ), cmap=py.cm.Greys,origin='lower',vmin=vmin_difference, vmax=vmax_difference)
      #py.colorbar()
      #py.savefig(average_difference_image_figname)
      #py.close()

      #form the average difference image in jy per pix and plot
      vmin_difference_jyperpix=-0.002
      vmax_difference_jyperpix=0.002
      average_difference_jyperpix_image=average_difference_jyperpix_image_array/difference_jyperpix_image_counter
      average_difference_jyperpix_image_crop = average_difference_jyperpix_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_difference_jyperpix_image_fitsname,average_difference_jyperpix_image_crop,clobber=True)
      pyfits.update(average_difference_jyperpix_image_fitsname,average_difference_jyperpix_image_crop,clobber=True,header=moon_header)
      
     #gc = aplpy.FITSFigure(average_difference_jyperpix_image_fitsname)
     #gc.show_grayscale()
     #gc.add_colorbar()
     #gc.set_theme('publication')
     #gc.save(average_difference_jyperpix_image_figname)
     #print "saved %s with aplpy" % average_difference_jyperpix_image_figname
     #gc.close()
     
      #py.figure(4)
      #py.clf()
      #py.title("Average difference Jy per pix image from %s obsids chan %s %s" % (str(difference_jyperpix_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_difference_jyperpix_image_crop ), cmap=py.cm.Greys,origin='lower')
      #py.colorbar()
      #py.savefig(average_difference_jyperpix_image_figname)
      #py.close()
            
      #form the average modelled moon image and plot
      moon_modelled_vmax=0.01
      moon_modelled_vmin=0.00
      average_moon_modelled_image=average_moon_modelled_image_array/modelled_moon_image_counter
      average_moon_modelled_image_crop = average_moon_modelled_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_moon_modelled_image_fitsname,average_moon_modelled_image_crop,clobber=True)
      pyfits.update(average_moon_modelled_image_fitsname,average_moon_modelled_image_crop,clobber=True,header=moon_header)
      
      #gc = aplpy.FITSFigure(average_moon_modelled_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_moon_modelled_image_figname)
      #print "saved %s with aplpy" % average_moon_modelled_image_figname
      #gc.close()
      
      #py.figure(5)
      #py.clf()
      #py.title("Average modelled moon image from %s obsids chan %s %s" % (str(modelled_moon_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_moon_modelled_image_crop ), cmap=py.cm.Greys,origin='lower',vmax=moon_modelled_vmax, vmin=moon_modelled_vmin)
      #py.colorbar()
      #py.savefig(average_moon_modelled_image_figname)
      #py.close()

      #form the average modelled rfi image and plot
      average_rfi_modelled_image=average_rfi_modelled_image_array/modelled_rfi_image_counter
      average_rfi_modelled_image_crop = average_rfi_modelled_image[ystart:yend,xstart:xend]
      pyfits.writeto(average_rfi_modelled_image_fitsname,average_rfi_modelled_image_crop,clobber=True)
      pyfits.update(average_rfi_modelled_image_fitsname,average_rfi_modelled_image_crop,clobber=True,header=moon_header)
      
      #gc = aplpy.FITSFigure(average_rfi_modelled_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_rfi_modelled_image_figname)
      #print "saved %s with aplpy" % average_rfi_modelled_image_figname
      #gc.close()
            
      #py.figure(7)
      #py.clf()
      #py.title("Average modelled rfi image from %s obsids chan %s %s" % (str(modelled_rfi_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_rfi_modelled_image_crop ), cmap=py.cm.Greys,origin='lower')
      #py.colorbar()
      #py.savefig(average_rfi_modelled_image_figname)
      #py.close()   
      
      #form the average modelled residual image and plot
      average_residual_modelled_image=average_residual_modelled_image_array/modelled_residual_image_counter
      average_residual_modelled_image_crop = average_residual_modelled_image[ystart:yend,xstart:xend]
      #average_residual_modelled_image_crop[75,0:50] = 1
      pyfits.writeto(average_residual_modelled_image_fitsname,average_residual_modelled_image_crop,clobber=True)
      pyfits.update(average_residual_modelled_image_fitsname,average_residual_modelled_image_crop,clobber=True,header=moon_header)

      #gc = aplpy.FITSFigure(average_residual_modelled_image_fitsname)
      #gc.show_grayscale()
      #gc.add_colorbar()
      #gc.set_theme('publication')
      #gc.save(average_residual_modelled_image_figname)
      #print "saved %s with aplpy" % average_residual_modelled_image_figname
      #gc.close()
      
      #py.figure(8)
      #py.clf()
      #py.title("Average modelled residual image from %s obsids chan %s %s" % (str(modelled_rfi_image_counter),centre_chan,chan_string))
      #py.imshow( ( average_residual_modelled_image_crop ), cmap=py.cm.Greys,origin='lower')
      #py.colorbar()
      #py.savefig(average_residual_modelled_image_figname)
      #py.close() 

      #1D profiles:
      # row and column sharing

      average_moon_modelled_hdulist = pyfits.open(average_moon_modelled_image_fitsname)
      average_moon_modelled_header=average_moon_modelled_hdulist[0].header
      #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
      f, (AX) = plt.subplots(2, 2, sharex='col', sharey='row')
      #f.tight_layout()
      ymin = -.015
      ymax = .055
      ystep = .02
      AX[0,0].invert_xaxis()
      AX[0,0].yaxis.set_ticks(np.arange(ymin, ymax, ystep))
      AX[0,0].set_ylim([ymin,ymax])
      #AX[0,0].add_label(0.07,0.93,'A',relative=True,size=25,style='normal',weight=boldness,color='black')
      AX[0,0].annotate("A", xy=(0.05,0.88), xycoords="axes fraction",fontweight='bold')
      AX[0,1].invert_xaxis()
      AX[1,0].yaxis.set_ticks(np.arange(ymin, ymax, ystep))
      AX[0,1].annotate("B", xy=(0.05,0.88), xycoords="axes fraction",fontweight='bold')
      AX[1,0].set_ylim([ymin,ymax])
      AX[1,0].annotate("C", xy=(0.05,0.88), xycoords="axes fraction",fontweight='bold')
      AX[1,1].annotate("D", xy=(0.05,0.88), xycoords="axes fraction",fontweight='bold')
      #plt.xlabel('Right Ascension (hrs:mins)')
      #f.gca().invert_xaxis()
      set_shared_ylabel(AX[:,0], 'Brightness (Jy/pixel)')
      f.text(0.5, 0.09, 'Right Ascension (hrs:mins)', ha='center', va='center')
      #f.text(0.06, 0.5, 'common ylabel', ha='center', va='center', rotation='vertical')

      x_ref_deg = average_moon_modelled_header['crval1']
      print "x_ref_deg is %s" % x_ref_deg
      x_delt_deg= average_moon_modelled_header['cdelt1']
      print "x_delt_deg is %s" % x_delt_deg
      x_length = average_moon_modelled_header['naxis1']
      print "x_length is %s" % x_length
      x_start_deg = x_ref_deg - (x_length/2)*x_delt_deg
      print "x_start_deg is %s" % x_start_deg
      x_end_deg=x_start_deg + (x_length*x_delt_deg)
      #x_start_hrs = (x_start_deg + 360)/15.
      #print "x_start_hrs is %s" % x_start_hrs      
      
      x_deg = np.linspace(x_start_deg,x_end_deg,x_length)
      x_hrs = (x_deg + 360)/15.
      #print x_hrs
      date_time_initial = datetime.datetime(2015, 9, 26, 00, 00)
      formatter = DateFormatter('%H:%M') 
      
      t = [date_time_initial + datetime.timedelta(hours=i) for i in x_hrs]
      
      #This is a horizontal line
      y=average_difference_jyperpix_image_crop[75,:]
      AX[0,0].plot(t, y,color='black')
      AX[0,0].xaxis.set_major_formatter(formatter)  

       
      y=average_moon_modelled_image_crop[75,:]
      AX[0,1].plot(t, y,color='black')
      # [left, bottom, width, height] 
      #AX[0,1].set_position([0.0,0.5,0.45,0.45])
      AX[0,1].xaxis.set_major_formatter(formatter) 
   
      y=average_rfi_modelled_image_crop[75,:]
      AX[1,0].plot(t, y,color='black')
      AX[1,0].xaxis.set_major_formatter(formatter) 

      y=average_residual_modelled_image_crop[75,:]
      AX[1,1].plot(t, y,color='black')
      AX[1,1].xaxis.set_major_formatter(formatter)
      
      f.autofmt_xdate()
      f.savefig('profile_fig_four_panel.png',dpi=900)

      #TODO:
      #figureout what slice you are plotting and if the RA matches!
      #set vmin vmax for yscale
      #set RA (hh:mm) and brightness (Jy/pixel) labels
      #smoosh subplots together more
      
      
      

      #make a four panel:
      font_size=15
      dx,dy = .47,.43
      boldness = 580
      four_panel_figname = 'fig1_four_panel.png'
      fig = py.figure(figsize=(17, 15))

      f1 = aplpy.FITSFigure(average_difference_jyperpix_image_fitsname, figure=fig, subplot=[0.05,0.5,dx,dy])
      f1.set_axis_labels_font(size=font_size)
      f1.add_label(0.07,0.93,'A',relative=True,size=25,style='normal',weight=boldness,color='black')
      f1.show_grayscale()
      f1.add_colorbar()
      f1.colorbar.set_font(size=font_size)
      f1.tick_labels.set_xformat('hh:mm')
      f1.tick_labels.set_yformat('dd:mm')
      f1.tick_labels.set_font(size=font_size)
      f1.set_theme('publication')
      
      f1.hide_xaxis_label()
      f1.hide_xtick_labels()
      
      f2 = aplpy.FITSFigure(average_moon_modelled_image_fitsname, figure=fig, subplot=[0.5,0.5,dx,dy])
      f2.add_label(0.07,0.93,'B',relative=True,size=25,style='normal',weight=boldness,color='black')
      f2.tick_labels.set_xformat('hh:mm')
      f2.tick_labels.set_yformat('dd:mm')
      f2.tick_labels.set_font(size=font_size)
      f2.show_grayscale()
      f2.add_colorbar()
      f2.colorbar.set_font(size=font_size)
      f2.colorbar.set_axis_label_font(size=font_size)
      f2.colorbar.set_axis_label_text('Jy/pixel')
      f2.set_theme('publication')
      
      f2.hide_yaxis_label()
      f2.hide_ytick_labels()
      f2.hide_xaxis_label()
      f2.hide_xtick_labels()
      
      f3 = aplpy.FITSFigure(average_rfi_modelled_image_fitsname, figure=fig, subplot=[0.05,0.05,dx,dy])
      f3.add_label(0.07,0.93,'C',relative=True,size=25,style='normal',weight=boldness,color='black')
      f3.set_axis_labels_font(size=font_size)
      f3.tick_labels.set_xformat('hh:mm')
      f3.tick_labels.set_yformat('dd:mm')
      f3.tick_labels.set_font(size=font_size)
      f3.show_grayscale()
      f3.add_colorbar()
      f3.colorbar.set_font(size=font_size)
      f3.set_theme('publication')
      
      f4 = aplpy.FITSFigure(average_residual_modelled_image_fitsname, figure=fig, subplot=[0.5,0.05,dx,dy])
      f4.add_label(0.07,0.93,'D',relative=True,size=25,style='normal',weight=boldness,color='white')
      f4.set_axis_labels_font(size=font_size)
      f4.tick_labels.set_xformat('hh:mm')
      f4.tick_labels.set_yformat('dd:mm')
      f4.tick_labels.set_font(size=font_size)
      f4.show_grayscale()
      f4.add_colorbar()
      f4.colorbar.set_font(size=font_size)
      f4.colorbar.set_axis_label_font(size=font_size)
      f4.colorbar.set_axis_label_text('Jy/pixel')
      f4.set_theme('publication')
      
      f4.hide_yaxis_label()
      f4.hide_ytick_labels()
      
      fig.canvas.draw()
      py.savefig(four_panel_figname,dpi=500)


  
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
parser.add_option('--new_modelled_only',action='store_true',dest='new_modelled_only',default=False,help='Only make the new modelled moon images [default=%default]')
parser.add_option('--epoch_ID',type='string', dest='epoch_ID',default='2018A_01',help='Epoch ID is a unique identifier for a set of paired of Moon/off Moon observations e.g.epoch)ID="2018A_01" [default=%default]')
parser.add_option('--old_labelling',action='store_true',dest='old_labelling',default=False,help='Use images with old labelling for 201B_05 without epoch ID in /data/moon/2017 [default=%default]')



(options, args) = parser.parse_args()

#freq = args[0]

make_images(options)






   
