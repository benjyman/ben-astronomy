#!/usr/bin/env python
#Direction-dpendent calibration with 'pre-peel'

from astropy.io import fits
import numpy as np
import os

generate_new_beams = False

subtraction_radius_arcmin = 7.0
subtraction_radius_deg = subtraction_radius_arcmin/60.0

#obsids_2015_2018 = ['1117031728','1199663088']
obsids_2015_2018 = ['1117031728']
coarse_chan = 145

def create_circular_mask(h, w, centre=None, radius=None):
    #give centre and array of two (x and y centre)
    if centre is None: # use the middle of the image
        centre = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(centre[0], centre[1], w - centre[0], h - centre[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_centre = np.sqrt((X - centre[0])**2 + (Y - centre[1])**2)
    mask = dist_from_centre <= radius
    return mask*1.


#/go to /md0/ATeam/CenA/DDCal

for obsid in obsids_2015_2018:
   print "obsid is %s" % obsid
   #1. Take the best image of CenA - want the model image
   fits_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-model-pb-I.fits"
   masked_fits_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-model-pb-masked-I.fits"
   pbuncorrect_input_name = masked_fits_image_filename.split('-I.fits')[0]
   beam_name = "%s_%s_beam" % (obsid,coarse_chan)
   uncorrected_image_name = "%s_%s_uncorr" % (obsid,coarse_chan)
   hi_res_image_for_self_cal_01_name = "%s_%s_selfcal_01" % (obsid,coarse_chan)
   self_cal_01_sols_name = "%s_%s_selfcal_01.bin" % (obsid,coarse_chan)
   hi_res_image_for_self_cal_02_name = "%s_%s_selfcal_02" % (obsid,coarse_chan)
   self_cal_02_sols_name = "%s_%s_selfcal_02.bin" % (obsid,coarse_chan)
   hi_res_image_for_self_cal_03_name = "%s_%s_selfcal_03" % (obsid,coarse_chan)
   self_cal_03_sols_name = "%s_%s_selfcal_03.bin" % (obsid,coarse_chan)
   hi_res_image_for_self_cal_04_name = "%s_%s_selfcal_04" % (obsid,coarse_chan)
   self_cal_04_sols_name = "%s_%s_selfcal_04.bin" % (obsid,coarse_chan)
   normal_res_image_subtr_core_name = "%s_%s_subtr_core" % (obsid,coarse_chan)
   
   
   #experiment with imaging only the core and then differencing the full and core-only images!
   #...
   core_only_image_name = "%s_%s_core_only" % (obsid,coarse_chan)
   
   
   
   hdulist = fits.open(fits_image_filename)
   header = hdulist[0].header
   image_data = hdulist[0].data[0,0,:,:]
   
   pixel_scale_deg = float(header['CDELT1'])
   imsize = int(header['NAXIS1'])
   image_centre = imsize/2.
   
   radius_pix = abs(subtraction_radius_deg/pixel_scale_deg)
   
   #2. Set the core to zero (-43 +- 7 arcmin)
   mask = create_circular_mask(imsize,imsize,centre=[image_centre,image_centre],radius=radius_pix)
   
   masked_image_data = mask * image_data
   #print masked_image_data
   #print masked_image_data[imsize/2-20:imsize/2+20,imsize/2-20:imsize/2+20]
   
   fits.writeto(masked_fits_image_filename,masked_image_data,clobber=True)
   fits.update(masked_fits_image_filename,masked_image_data,header=header)
   print "wrote  image %s" %  masked_fits_image_filename
   
   #Do I need to pb-uncorrect to get xx,yy,xy,yx?? I think so... 
   #Need the beams for each obsid 
   #get metafits
   
   #see if its 2015 data or 2018 data
   metafits_filename = "/md0/ATeam/CenA/2015/%s/%s/%s_metafits_ppds.fits" % (coarse_chan,obsid,obsid)
   subtr_ms_name = "/md0/ATeam/CenA/DDCal/%s_subtr.ms " % (obsid)
   core_only_ms_name = "/md0/ATeam/CenA/DDCal/%s_core_only.ms " % (obsid)
   ms_copy_name = "/md0/ATeam/CenA/DDCal/%s_copy.ms " % (obsid)
   if os.path.isfile(metafits_filename) and os.access(metafits_filename, os.R_OK):
      obsid_year = 2015
      orig_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s.ms " % (obsid_year,coarse_chan,obsid,obsid)
   else:
      obsid_year = 2018
      metafits_filename = "/md0/ATeam/CenA/2018/%s/%s/%s_metafits_ppds.fits" % (coarse_chan,obsid,obsid)
      orig_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s.ms " % (obsid_year,coarse_chan,obsid,obsid)
      
   if generate_new_beams:
      cmd = "beam -2016 -proto %s -name %s -m %s " % (fits_image_filename,beam_name,metafits_filename)
      print cmd
      os.system(cmd)
      
   #uncorrect image with beam
   cmd = "pbcorrect -uncorrect %s model.fits %s %s" % (uncorrected_image_name,beam_name,pbuncorrect_input_name)
   #print cmd
   #os.system(cmd)
 
   #3  make a copy of the ms's
   cmd = "rm -rf %s" % (ms_copy_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "rm -rf %s" % (subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "rm -rf %s" % (core_only_ms_name)
   print cmd
   os.system(cmd)
   
   cmd = "cp -r %s %s" % (orig_ms_name,ms_copy_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "cp -r %s %s" % (orig_ms_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "cp -r %s %s" % (orig_ms_name,core_only_ms_name)
   print cmd
   os.system(cmd)
   
   #4. for subtracting ... make a uvmodel with only the core (wsclean -predict)
   cmd = "wsclean -predict -name %s -size %s %s -scale %.4f -pol xx,xy,yx,yy %s " % (uncorrected_image_name,int(imsize),int(imsize),abs(float(pixel_scale_deg)),subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   #5. subtract the core model from the visibilities (use subtrmodel -usemodelcol)
   cmd = "subtrmodel -usemodelcol -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)


   #for core only image ... make a uvmodel with only the core (wsclean -predict)
   cmd = "wsclean -predict -name %s -size %s %s -scale %.4f -pol xx,xy,yx,yy %s " % (uncorrected_image_name,int(imsize),int(imsize),abs(float(pixel_scale_deg)),core_only_ms_name)
   print cmd
   os.system(cmd) 
   
   #image the core-only model
   cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 3 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column MODEL_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 12 -join-channels -fit-spectral-pol 3 %s" % (core_only_image_name,core_only_ms_name)
   print cmd
   os.system(cmd)
   
   sys.exit()
   
   
   
   #5a image subtracted model at smaller imsize and higher res (more w terms?) shallow clean 
   ###Dont do this anymore###############
   #cmd = "wsclean -name %s -size 2048 2048 -auto-threshold 8 -auto-mask 10 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale 8asec -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (hi_res_image_for_self_cal_01_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   #It worked!!
   
   
   #6. Self-cal on the core-only ms
   #cmd = "calibrate -minuv 60 -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_01_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "applysolutions -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_01_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "wsclean -name %s -size 2048 2048 -auto-threshold 4 -auto-mask 6 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale 8asec -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (hi_res_image_for_self_cal_02_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   #again
   #cmd = "calibrate -minuv 60 -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_02_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "applysolutions  -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_02_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "wsclean -name %s -size 2048 2048 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale 8asec -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (hi_res_image_for_self_cal_03_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   #again
   #cmd = "calibrate -minuv 60 -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_03_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "applysolutions  -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_03_sols_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "wsclean -name %s -size 2048 2048 -auto-threshold 0.3 -auto-mask 1 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale 8asec -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (hi_res_image_for_self_cal_04_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   #7. Reimage the 'just core' at higher angular res
   #8. Do another self cal
   
   #take a look at the cal sols
   #cmd = "aocal_plot.py %s" % self_cal_01_sols_name
   #print cmd
   #os.system(cmd)
   
   #cmd = "aocal_plot.py %s" % self_cal_02_sols_name
   #print cmd
   #os.system(cmd)
   
   #cmd = "aocal_plot.py %s" % self_cal_03_sols_name
   #print cmd
   #os.system(cmd)
   ######################################################################
   
   #9. Subtract the model from the unsubtracted ms    subtrmodel <model> <ms>
   #nno this doesnt work
   #cmd = "subtrmodel -usemodelcol -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,ms_copy_name)
   #print cmd
   #os.system(cmd)
   
   
   #10 re-image at full image size 
   cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 3 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 12 -join-channels -fit-spectral-pol 3 %s" % (normal_res_image_subtr_core_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)

   #It seemed to work, try self cal, BUT you are probably going to have to peel Cen B (big source out of image to south)
   cmd = "calibrate -minuv 60 -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_01_sols_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "applysolutions  -datacolumn CORRECTED_DATA %s %s " % (subtr_ms_name,self_cal_01_sols_name)
   #print cmd
   #os.system(cmd)
   
   cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 3 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 12 -join-channels -fit-spectral-pol 3 %s" % (hi_res_image_for_self_cal_01_name,subtr_ms_name)
   #print cmd
   #os.system(cmd)
   
   
   
   
   
   

