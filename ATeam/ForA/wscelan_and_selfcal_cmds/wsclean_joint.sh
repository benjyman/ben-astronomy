#joint deconvolution
wsclean -name 2014B_2018A_joint_idg_stokesI -size 4096 4096 -niter 1000 -threshold 2.5 -multiscale -mgain 0.85 -data-column CORRECTED_DATA -scale 0.004 -weight briggs 0 -small-inversion -make-psf -pol I -use-idg -grid-with-beam  -channels-out 8  2014B/1102865728/1102865728.ms 2014B/1102865128/1102865128.ms 2014B/1102864528/1102864528.ms 2018A/1202816352/1202816352.ms 2018A/1202815752/1202815752.ms 2018A/1202815152/1202815152.ms
