#!/usr/bin/env python
#try ionpeeling out the core!

from astropy.io import fits
import numpy as np
import os

obsids_2015_2018 = ['1117031728']
coarse_chan = 145
hi_uv_cutoff_m = 1500
normal_uv_cutoff = 60

for obsid in obsids_2015_2018:
   print "obsid is %s" % obsid
   image_filename = "/md0/ATeam/CenA/CenA_2015_2018_joint_idg_12_obs_145_selfcal_02-MFS-image-pb.fits"
   skymodel_filename = "/md0/ATeam/CenA/CenA_2015_2018_joint_idg_12_obs_145_selfcal_02-sources.txt"
   aoskymodel_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02-aosources.txt"
   aoskymodel_core_only_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_core_only-aosources.txt"
   aoskymodel_nil_only_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_NIL_only-aosources.txt"
   aoskymodel_sil_only_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_SIL_only-aosources.txt"
   rendered_output_core_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_core_only_aomodel.fits"
   rendered_output_nil_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_NIL_only_aomodel.fits"
   rendered_output_sil_image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_SIL_only_aomodel.fits"
   clustered_core_model_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_core_only_clustered-aosources.txt"
   clustered_nil_model_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_NIL_only_clustered-aosources.txt"
   clustered_sil_model_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_02_SIL_only_clustered-aosources.txt"
   high_uv_cutoff_solutions_filename = "%s_core_model_high_uv_cut.bin" % (obsid)
   timestep_solutions_filename = "%s_core_model_timestep_cal.bin" % (obsid)
   solutions_filename = "%s_core_model_cal.bin" % (obsid)
   ionpeel_solution_filename = "%s_ionpeel_core.bin" % (obsid)
   ionpeel_solution_filename_nil = "%s_ionpeel_nil.bin" % (obsid)
   metafits_filename = "/md0/ATeam/CenA/2015/%s/%s/%s_metafits_ppds.fits" % (coarse_chan,obsid,obsid)
   wsclean_timestep_cal_imagename = "%s_%s_core_model_timestep_cal" % (obsid,coarse_chan)
   wsclean_high_uvcutoff_cal_imagename = "%s_%s_core_model_cal_high_uv_cut" % (obsid,coarse_chan)
   wsclean_no_negative_imagename = "%s_%s_core_model_cal_no_negative" % (obsid,coarse_chan)
   wsclean_no_negative_no_si_240_wlayers_imagename = "%s_%s_core_model_cal_no_neg_no_si_240w" % (obsid,coarse_chan)
   wsclean_core_ionpeeled_imagename = "%s_%s_core_ionpeeled" % (obsid,coarse_chan)
   wsclean_nil_ionpeeled_imagename = "%s_%s_nil_core_ionpeeled" % (obsid,coarse_chan)
   wsclean_sil_ionpeeled_imagename = "%s_%s_sil_core_ionpeeled" % (obsid,coarse_chan)
   core_model_name = "/md0/ATeam/CenA/CenA_core_wsclean_model.txt"
   
   if os.path.isfile(metafits_filename) and os.access(metafits_filename, os.R_OK):
      obsid_year = 2015
      orig_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s.ms " % (obsid_year,coarse_chan,obsid,obsid)
      ionpeel_ms_name = "%s_ionpeel_core.ms " % (obsid)
   else:
      obsid_year = 2018
      metafits_filename = "/md0/ATeam/CenA/2018/%s/%s/%s_metafits_ppds.fits" % (coarse_chan,obsid,obsid)
      orig_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s.ms " % (obsid_year,coarse_chan,obsid,obsid)
      ionpeel_ms_name = "%s_ionpeel_core.ms " % (obsid)
   
   #make copy
   #cmd = "rm -rf %s" % (ionpeel_ms_name)
   #print cmd
   #os.system(cmd)
   
   #cmd = "cp -r %s %s" % (orig_ms_name,ionpeel_ms_name)
   #print cmd
   #os.system(cmd)
   


   #cmd = "calibrate -minuv %s -ch 1 -m %s  -applybeam %s %s" % (normal_uv_cutoff,core_model_name,ionpeel_ms_name,solutions_filename)
   #print cmd
   #os.system(cmd)
   #cmd = "applysolutions %s %s" % (ionpeel_ms_name,solutions_filename)
   #print cmd
   #os.system(cmd)

   #image with -no-negative, -no-small-inversion and 240 wlayers!
   cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 3 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale .004 -weight uniform -no-small-inversion -make-psf -pol I -apply-primary-beam -pb-undersampling 4 -channels-out 12 -join-channels -no-negative -nwlayers 240  %s" % (wsclean_no_negative_no_si_240_wlayers_imagename,ionpeel_ms_name)
   print cmd
   os.system(cmd)
   
   

   #cmd = "calibrate -minuv %s -ch 1 -datacolumn DATA -m %s -applybeam %s %s" % (normal_uv_cutoff,core_model_name,ionpeel_ms_name, solutions_filename)
   #print cmd
   #os.system(cmd)
   #cmd = "applysolutions %s %s" % (ionpeel_ms_name,timestep_solutions_filename)
   #print cmd
   #os.system(cmd)

   #image!
   #cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale .004 -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels %s" % (wsclean_timestep_cal_imagename,ionpeel_ms_name)
   #print cmd
   #os.system(cmd)
   
   

   #convert to andre format
   #cmd = "bbs2model %s %s" % (skymodel_filename,aoskymodel_filename)
   #print cmd
   #os.system(cmd)
   #cut out just the core (7 amin / 0.116666666 deg radius):
   #cmd = "editmodel -m %s  -near  13h25m27.6s -43d01m08.8s 0.116666666  %s" % (aoskymodel_core_only_filename,aoskymodel_filename)
   #print cmd
   #os.system(cmd)
   #render the core model to check
   #cmd = "render -t %s  -o %s  -r  %s " % (image_filename,rendered_output_core_image_filename,aoskymodel_core_only_filename)
   #print cmd
   #os.system(cmd)
   
   #So now we have a core model to ionpeel (may need to edit it more to get rid of negative and small components etc)
   #I think I need to run cluster - but with only 1 direction since it is just the core
   #cmd = "cluster %s %s 1 " % (aoskymodel_core_only_filename,clustered_core_model_filename)
   #print cmd
   #os.system(cmd)
  
   #observation 1117031728 is 112 s long so lets do 4 sec intervals -t 28 (or 8 sec? -t 14)
   #cmd = "ionpeel -t 28 %s %s %s " % (ionpeel_ms_name, clustered_core_model_filename, ionpeel_solution_filename)
   #print cmd
   #os.system(cmd)
   
   #image the ionpeeled ms
   #cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale .004 -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (wsclean_core_ionpeeled_imagename,ionpeel_ms_name)
   #print cmd
   #os.system(cmd)
   
   ###
   #cut out just the NIL (3 amin / 0.05 deg radius):
   #cmd = "editmodel -m %s  -near  13h25m43.3s -42d57m32.7s 0.05  %s" % (aoskymodel_nil_only_filename,aoskymodel_filename)
   #print cmd
   #os.system(cmd)
   #render the core model to check
   #cmd = "render -t %s  -o %s  -r  %s " % (image_filename,rendered_output_nil_image_filename,aoskymodel_nil_only_filename)
   #print cmd
   #os.system(cmd)
   #cmd = "cluster %s %s 1 " % (aoskymodel_nil_only_filename,clustered_nil_model_filename)
   #print cmd
   #os.system(cmd)
  
   #observation 1117031728 is 112 s long so lets do 4 sec intervals -t 28 (or 8 sec? -t 14)
   #cmd = "ionpeel -t 1 -groupchannels 768 -datacolumn CORRECTED_DATA %s %s %s " % (ionpeel_ms_name, clustered_nil_model_filename, ionpeel_solution_filename_nil)
   #print cmd
   #os.system(cmd)
   
   #image the ionpeeled ms
   #cmd = "wsclean -name %s -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -data-column CORRECTED_DATA -scale .004 -weight uniform -small-inversion -make-psf -pol xx,xy,yx,yy -join-polarizations -channels-out 12 -join-channels  %s" % (wsclean_nil_ionpeeled_imagename,ionpeel_ms_name)
   #print cmd
   #os.system(cmd)
   
   
   
