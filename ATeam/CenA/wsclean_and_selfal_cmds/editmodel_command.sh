#bbs2model CenA_2015_2018_joint_idg_stokesI_spec_pol_1_selfcal_01-sources.txt CenA_2015_2018_joint_idg_stokesI_spec_pol_1_selfcal_01-aosources.txt
#ra dec same as chgcentre, dist in deg, -t flux in Jy
#editmodel -m CenA_01.txt  -near  13h25m27.6s -43d01m08.8s 4.0  CenA_2015_2018_joint_idg_stokesI_spec_pol_1_selfcal_01-aosources.txt 
#editmodel -m CenA_only_01.txt -collect CenA CenA_01.txt
#editmodel -m CenA_01_spec_index.txt -replace-si -0.7 CenA_only_01.txt
#don't need to scale , total flux looks about right
#editmodel -m Fornax_A_only_01_scaled.txt 
#render -t CenA_2015_2018_joint_idg_stokesI_spec_pol_1_selfcal_01-MFS-image-pb.fits  -o  CenA_01_aomodel.fits -r  CenA_01_spec_index.txt  

#new one 12 obs
#bbs2model CenA_2015_2018_joint_idg_12_obs_self_cal_04-sources.txt CenA_2015_2018_joint_idg_12_obs_self_cal_04-aosources.txt 
#ra dec same as chgcentre, dist in deg, -t flux in Jy
#editmodel -m CenA_core.txt  -near 13h25m27.6s -43d01m08.8s  7amin CenA_2015_2018_joint_idg_12_obs_self_cal_04-aosources.txt
 
#check to see flux scale etc
#render -r  -t  CenA_2015_2018_joint_idg_12_obs_self_cal_04-MFS-image-pb.fits -o CenA_core.fits CenA_core.txt

#collect into one source
#editmodel -m CenA_core_collect.txt -collect CenA_core CenA_core.txt 

#total flux density of inner lobes only 278420.41 mJy (this is because the image used for initial cal was condon et al 1996 CenA from VLA, with
#core 278000 mJy! Alvarez et al 2000 show inner lobes have an SI of -0.7, so just scale to 185 MHz
#At 185MHz should be S_185 = 278*(185/1425)^-0.7 = 1161 Jy (!)
#So scale the model 
#editmodel -m CenA_core_collect_scaled.txt -sc 1161 185 278 1425 CenA_core_collect.txt
#check to see flux scale etc
#render -r  -t  CenA_2015_2018_joint_idg_12_obs_self_cal_04-MFS-image-pb.fits -o CenA_core_collect_scaled.fits CenA_core_collect_scaled.txt
#Looks good
#fix small position error (see initial_offset.png)
#editmodel -m CenA_core_collect_scaled_shift.txt -shift -14asec -10asec CenA_core_collect_scaled.txt
#check
#render -r  -t  CenA_2015_2018_joint_idg_12_obs_self_cal_04-MFS-image-pb.fits -o CenA_core_collect_scaled_shift.fits CenA_core_collect_scaled_shift.txt

#All looks good, model ready for testing
#give it a better name:
mv CenA_core_collect_scaled_shift.txt CenA_core_wsclean_model.txt



