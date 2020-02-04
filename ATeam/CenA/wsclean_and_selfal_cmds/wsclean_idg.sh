#Try Chan 169!
#wsclean -name CenA_2015_2018_joint_idg_11_obs_169 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

#wsclean -name CenA_2015_2018_joint_idg_11_obs_169_selcal_01 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

#wsclean -name CenA_2015_2018_joint_idg_11_obs_169_selcal_02 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

#chan 169 robust 0:
#wsclean -name CenA_2015_2018_joint_idg_11_obs_169_selfcal_03_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

#Chan 121
#wsclean -name CenA_2015_2018_joint_idg_11_obs_121 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.0045 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/121/1200777400/1200777400.ms 2018/121/1200345160/1200345160.ms 2018/121/1200604504/1200604504.ms 2018/121/1200518056/1200518056.ms 2018/121/1200431608/1200431608.ms 2015/121/1117031608/1117031608.ms 2015/121/1121334352/1121334352.ms 2015/121/1121325688/1121325688.ms 2015/121/1121247928/1121247928.ms 2015/121/1121321184/1121321184.ms 2015/121/1121322984/1121322984.ms

#wsclean -name CenA_2015_2018_joint_idg_11_obs_121_selfcal_01 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.0045 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/121/1200777400/1200777400.ms 2018/121/1200345160/1200345160.ms 2018/121/1200604504/1200604504.ms 2018/121/1200518056/1200518056.ms 2018/121/1200431608/1200431608.ms 2015/121/1117031608/1117031608.ms 2015/121/1121334352/1121334352.ms 2015/121/1121325688/1121325688.ms 2015/121/1121247928/1121247928.ms 2015/121/1121321184/1121321184.ms 2015/121/1121322984/1121322984.ms

#wsclean -name CenA_2015_2018_joint_idg_11_obs_121_selfcal_02 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.0045 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2018/121/1200777400/1200777400.ms 2018/121/1200345160/1200345160.ms 2018/121/1200604504/1200604504.ms 2018/121/1200518056/1200518056.ms 2018/121/1200431608/1200431608.ms 2015/121/1117031608/1117031608.ms 2015/121/1121334352/1121334352.ms 2015/121/1121325688/1121325688.ms 2015/121/1121247928/1121247928.ms 2015/121/1121321184/1121321184.ms 2015/121/1121322984/1121322984.ms


#wsclean -name CenA_2015_2018_joint_idg_11_obs_121_selfcal_03_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.0045 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2 2018/121/1200777400/1200777400.ms 2018/121/1200345160/1200345160.ms 2018/121/1200604504/1200604504.ms 2018/121/1200518056/1200518056.ms 2018/121/1200431608/1200431608.ms 2015/121/1117031608/1117031608.ms 2015/121/1121334352/1121334352.ms 2015/121/1121325688/1121325688.ms 2015/121/1121247928/1121247928.ms 2015/121/1121321184/1121321184.ms 2015/121/1121322984/1121322984.ms

#No change:
#radial spikes much worse withlower freq. try making finer scale images at low freq.
#wsclean -name CenA_2015_2018_joint_idg_11_obs_121_selfcal_03_robust0_scale_change -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 24 -join-channels -fit-spectral-pol 2 2018/121/1200777400/1200777400.ms 2018/121/1200345160/1200345160.ms 2018/121/1200604504/1200604504.ms 2018/121/1200518056/1200518056.ms 2018/121/1200431608/1200431608.ms 2015/121/1117031608/1117031608.ms 2015/121/1121334352/1121334352.ms 2015/121/1121325688/1121325688.ms 2015/121/1121247928/1121247928.ms 2015/121/1121321184/1121321184.ms 2015/121/1121322984/1121322984.ms

#redo image for chan 145, to get calibration/flux scale same as other chans
#wsclean -name CenA_2015_2018_joint_idg_12_obs_145 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms 

#145 self cal
#wsclean -name CenA_2015_2018_joint_idg_12_obs_145_selfcal_01 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs_145_selfcal_02 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2 2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms

#145 robust 0 (new cal for flux scale
#wsclean -name CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2 2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms

#flag extra side chan (3)
#wsclean -name CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_flag3 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2 2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms


#wsclean -name CenA_2015_2018_joint_idg_stokesI_spec_pol_1 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 1  2015/1117031728/1117031728.ms  2018/1199663088/1199663088.ms

#wsclean -name CenA_2015_2018_joint_idg_stokesI_spec_pol_1 -size 4096 4096 -auto-threshold 3 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 1  2015/1117031728/1117031728.ms  2018/1199663088/1199663088.ms

#wsclean -name CenA_2015_2018_joint_idg_stokesI_spec_pol_1_selfcal_01 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 1  2015/1117031728/1117031728.ms  2018/1199663088/1199663088.ms

#wsclean -name CenA_2015_2018_joint_idg_stokesI_spec_pol_3_selfcal_02 -size 4096 4096 -auto-threshold 1 -auto-mask 5 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 3  2015/1117031728/1117031728.ms  2018/1199663088/1199663088.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms

#done:
#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_01 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms 

#run this next (but keep in mind that self-cal may not fix things - need to correct using Torrance's ATeams to properly take into account ionospheric shifts
# and jointly deconvolve with wscelan+idg (but maybe peeling the core using an individualised model for each snapshot might help...)
#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_02 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_04 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_04_robust0 -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_04_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms

#wsclean -name CenA_2015_2018_joint_idg_12_obs_self_cal_05_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 10 -join-channels -fit-spectral-pol 2  2015/1117031728/1117031728.ms 2015/1121420968/1121420968.ms  2015/1121593824/1121593824.ms 2015/1121334536/1121334536.ms  2015/1121507392/1121507392.ms  2015/1121680256/1121680256.ms 2018/1199663088/1199663088.ms  2018/1200691136/1200691136.ms  2018/1200864032/1200864032.ms 2018/1200604688/1200604688.ms  2018/1200777584/1200777584.ms  2018/1200950480/1200950480.ms


#try joint deconvolution of 145 and 169 robust 0! (ambitious)
#Can't do this on desktop machine - not enough memory - need gpu - or more memory (magnus/galaxy?) 
##August 2019 - retry this by only using a small portion of the time in the mss, SNR is high ... maybe less time and more uvcoverage will help....?
#wsclean -name CenA_2015_2018_joint_idg_23_obs_145_169_self_cal_03_robust0 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.86 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 24 -join-channels -fit-spectral-pol 2  2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

#wsclean -name CenA_2015_2018_joint_idg_23_obs_145_169_self_cal_03_robust0_1timestep -interval 1 3 -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.86 -save-source-list -data-column CORRECTED_DATA -scale 0.003 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 24 -join-channels -fit-spectral-pol 2  2015/145/1117031728/1117031728.ms 2015/145/1121420968/1121420968.ms  2015/145/1121593824/1121593824.ms 2015/145/1121334536/1121334536.ms  2015/145/1121507392/1121507392.ms  2015/145/1121680256/1121680256.ms 2018/145/1199663088/1199663088.ms  2018/145/1200691136/1200691136.ms  2018/145/1200864032/1200864032.ms 2018/145/1200604688/1200604688.ms  2018/145/1200777584/1200777584.ms  2018/145/1200950480/1200950480.ms 2018/169/1200864208/1200864208.ms 2018/169/1200777760/1200777760.ms 2018/169/1200691312/1200691312.ms 2018/169/1200604864/1200604864.ms 2018/169/1200431968/1200431968.ms 2015/169/1117031848/1117031848.ms 2015/169/1121594000/1121594000.ms 2015/169/1121507576/1121507576.ms 2015/169/1121421144/1121421144.ms 2015/169/1121334712/1121334712.ms 2015/169/1121248288/1121248288.ms

wsclean -name CenA_2015_2018_joint_idg_23_obs_145_169_self_cal_03_robust0_16secsplit -size 4096 4096 -niter 200000  -threshold 0.005  -multiscale -mgain 0.86 -save-source-list -data-column DATA -scale 0.003 -weight briggs 0  -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 24 -join-channels -fit-spectral-pol 2  2015/145/1117031728/1117031728_split.ms 2015/145/1121420968/1121420968_split.ms  2015/145/1121593824/1121593824_split.ms 2015/145/1121334536/1121334536_split.ms  2015/145/1121507392/1121507392_split.ms  2015/145/1121680256/1121680256_split.ms 2018/145/1199663088/1199663088_split.ms  2018/145/1200691136/1200691136_split.ms  2018/145/1200864032/1200864032_split.ms 2018/145/1200604688/1200604688_split.ms  2018/145/1200777584/1200777584_split.ms  2018/145/1200950480/1200950480_split.ms 2018/169/1200864208/1200864208_split.ms 2018/169/1200777760/1200777760_split.ms 2018/169/1200691312/1200691312_split.ms 2018/169/1200604864/1200604864_split.ms 2018/169/1200431968/1200431968_split.ms 2015/169/1117031848/1117031848_split.ms 2015/169/1121594000/1121594000_split.ms 2015/169/1121507576/1121507576_split.ms 2015/169/1121421144/1121421144_split.ms 2015/169/1121334712/1121334712_split.ms 2015/169/1121248288/1121248288_split.ms


#What does it look like if you just image one obsid from 2018?
#wsclean -name CenA_2018_1199663088_145_self_cal_05_robust0 -size 4096 4096 -niter 100000  -threshold 0.025  -multiscale -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0  -small-inversion -make-psf -pol I -apply-primary-beam -pb-undersampling 4  -channels-out 1  2018/145/1199663088/1199663088.ms

#weird aliasing on the above. Try uniform weighting and no multiscale?
#wsclean -name CenA_2018_1199663088_145_self_cal_05_uniform_auto -size 4096 4096 -auto-threshold 1 -auto-mask 3 -multiscale -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -pol I -apply-primary-beam -pb-undersampling 4  2018/145/1199663088/1199663088.ms

#still weird aliasing but interestingly, no radial features. Try without multiscale
#wsclean -name CenA_2018_1199663088_145_self_cal_05_uniform_autoi_no_multiscale -size 4096 4096 -auto-threshold 1 -auto-mask 3 -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -pol I -apply-primary-beam -pb-undersampling 4  2018/145/1199663088/1199663088.ms

#removing multiscale gets rid of the weird aliasing effects.
#Try 2015 Ph1 only and see difference
#wsclean -name CenA_2015_1117031728_145_self_cal_05_uniform_auto_no_multiscale -size 4096 4096 -auto-threshold 1 -auto-mask 3 -niter 1000000 -mgain 0.85 -save-source-list -data-column CORRECTED_DATA -scale 0.004 -weight uniform -small-inversion -pol I -apply-primary-beam -pb-undersampling 4 2015/145/1117031728/1117031728.ms

