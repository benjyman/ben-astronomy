#bbs2model 2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-sources.txt  2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-ao_sources.txt
#ra dec same as chgcentre, dist in deg, -t flux in Jy
#editmodel -m Fornax_A_01.txt  -near  03h22m41.7s -37d12m30.0s 0.65  2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-ao_sources.txt 
#editmodel -m Fornax_A_only_01.txt -collect Fornax_A Fornax_A_01.txt 
#mv Fornax_A_only_01.txt  Fornax_A_187_MHz.txt
#sed -i -e 's/polynomial/spectral-index/g' Fornax_A_187_MHz.txt
#editmodel -m Fornax_A_187_MHz_spec_index.txt -replace-si -0.8 Fornax_A_187_MHz.txt
#editmodel -rts Fornax_A_187_MHz_spec_index_RTS.txt Fornax_A_187_MHz_spec_index.txt
#don't need to scale , total flux looks about right
#editmodel -m Fornax_A_only_01_scaled.txt 
#render -t 2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-MFS-image-pb.fits  -o  Fornax_A_01_aomodel_gaussians.fits  -gaussians Fornax_A_01.txt  
#render -t 2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-MFS-image-pb.fits  -o  Fornax_A_01_aomodel_points.fits  Fornax_A_01.txt

#Fix position shift of .024deg DEC .016deg ra
#editmodel -m Fornax_A_187_MHz_spec_index_shift.txt -shift -42asec -34asec Fornax_A_187_MHz_spec_index.txt
#render -r  -t 2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-MFS-image-pb.fits  -o  Fornax_A_shift.fits Fornax_A_187_MHz_spec_index_shift.txt
#render  -t 2014B_2018A_joint_idg_stokesI_deeper_self_cal_02_spec_pol_1-MFS-image-pb.fits  -o  Fornax_A_shift_nobeam.fits Fornax_A_187_MHz_spec_index_shift.txt
#editmodel -rts Fornax_A_187_MHz_spec_index_shift_RTS.txt Fornax_A_187_MHz_spec_index_shift.txt


#editmodel -m Fornax_A_187_MHz_spec_index_shift_flux_cut_01.txt -tc 0.05 Fornax_A_187_MHz_spec_index_shift.txt 
editmodel -rts Fornax_A_187_MHz_spec_index_shift_flux_cut_01_RTS.txt Fornax_A_187_MHz_spec_index_shift_flux_cut_01.txt



