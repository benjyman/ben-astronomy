#!/usr/bin/env bash

flagmwa   2018/145/1199663088/1199663088.ms 2 24 1 3
flagmwa   2018/145/1200604688/1200604688.ms 2 24 1 3
flagmwa   2018/145/1200691136/1200691136.ms 2 24 1 3
flagmwa   2018/145/1200777584/1200777584.ms 2 24 1 3
flagmwa   2018/145/1200864032/1200864032.ms 2 24 1 3
flagmwa   2018/145/1200950480/1200950480.ms 2 24 1 3

flagmwa   2015/145/1117031728/1117031728.ms 2 24 1 3
flagmwa   2015/145/1121334536/1121334536.ms 2 24 1 3
flagmwa   2015/145/1121420968/1121420968.ms 2 24 1 3
flagmwa   2015/145/1121507392/1121507392.ms 2 24 1 3
flagmwa   2015/145/1121593824/1121593824.ms 2 24 1 3
flagmwa   2015/145/1121680256/1121680256.ms 2 24 1 3


#applysolutions 2018/1200864032/1200864032.ms 2018/1200777584/1200777584_selfcal_joint_06.bin
#applysolutions 2018/1200950480/1200950480.ms 2018/1200777584/1200777584_selfcal_joint_06.bin


#try to calibrate on model I gave to Gemma:
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1199663088/1199663088.ms 2018/145/1199663088/1199663088_core_model.bin
#applysolutions 2018/145/1199663088/1199663088.ms 2018/145/1199663088/1199663088_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1200604688/1200604688.ms 2018/145/1200604688/1200604688_core_model.bin
#applysolutions 2018/145/1200604688/1200604688.ms 2018/145/1200604688/1200604688_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1200691136/1200691136.ms 2018/145/1200691136/1200691136_core_model.bin
#applysolutions 2018/145/1200691136/1200691136.ms 2018/145/1200691136/1200691136_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1200777584/1200777584.ms 2018/145/1200777584/1200777584_core_model.bin
#applysolutions 2018/145/1200777584/1200777584.ms 2018/145/1200777584/1200777584_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1200864032/1200864032.ms 2018/145/1200864032/1200864032_core_model.bin
#applysolutions 2018/145/1200864032/1200864032.ms 2018/145/1200864032/1200864032_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/145/1200950480/1200950480.ms 2018/145/1200950480/1200950480_core_model.bin
#applysolutions 2018/145/1200950480/1200950480.ms 2018/145/1200950480/1200950480_core_model.bin

#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1117031728/1117031728.ms 2015/145/1117031728/1117031728_core_model.bin
#applysolutions 2015/145/1117031728/1117031728.ms 2015/145/1117031728/1117031728_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1121334536/1121334536.ms 2015/145/1121334536/1121334536_core_model.bin
#applysolutions 2015/145/1121334536/1121334536.ms 2015/145/1121334536/1121334536_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1121420968/1121420968.ms 2015/145/1121420968/1121420968_core_model.bin
#applysolutions 2015/145/1121420968/1121420968.ms 2015/145/1121420968/1121420968_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1121507392/1121507392.ms 2015/145/1121507392/1121507392_core_model.bin
#applysolutions 2015/145/1121507392/1121507392.ms 2015/145/1121507392/1121507392_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1121593824/1121593824.ms 2015/145/1121593824/1121593824_core_model.bin
#applysolutions 2015/145/1121593824/1121593824.ms 2015/145/1121593824/1121593824_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/145/1121680256/1121680256.ms 2015/145/1121680256/1121680256_core_model.bin
#applysolutions 2015/145/1121680256/1121680256.ms 2015/145/1121680256/1121680256_core_model.bin

##selfcal
#calibrate -minuv 60 -ch 1  2018/145/1199663088/1199663088.ms 2018/145/1199663088/1199663088_145_selfcal_joint_03.bin
#applysolutions 2018/145/1199663088/1199663088.ms 2018/145/1199663088/1199663088_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2018/145/1200604688/1200604688.ms 2018/145/1200604688/1200604688_145_selfcal_joint_03.bin
#applysolutions 2018/145/1200604688/1200604688.ms 2018/145/1200604688/1200604688_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2018/145/1200691136/1200691136.ms 2018/145/1200691136/1200691136_145_selfcal_joint_03.bin
#applysolutions 2018/145/1200691136/1200691136.ms 2018/145/1200691136/1200691136_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2018/145/1200777584/1200777584.ms 2018/145/1200777584/1200777584_145_selfcal_joint_03.bin
#applysolutions 2018/145/1200777584/1200777584.ms 2018/145/1200777584/1200777584_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2018/145/1200864032/1200864032.ms 2018/145/1200864032/1200864032_145_selfcal_joint_03.bin
#applysolutions 2018/145/1200864032/1200864032.ms 2018/145/1200864032/1200864032_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2018/145/1200950480/1200950480.ms 2018/145/1200950480/1200950480_145_selfcal_joint_03.bin
#applysolutions 2018/145/1200950480/1200950480.ms 2018/145/1200950480/1200950480_145_selfcal_joint_03.bin

#calibrate -minuv 60 -ch 1  2015/145/1117031728/1117031728.ms 2015/145/1117031728/1117031728_145_selfcal_joint_03.bin
#applysolutions 2015/145/1117031728/1117031728.ms 2015/145/1117031728/1117031728_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2015/145/1121334536/1121334536.ms 2015/145/1121334536/1121334536_145_selfcal_joint_03.bin
#applysolutions 2015/145/1121334536/1121334536.ms 2015/145/1121334536/1121334536_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2015/145/1121420968/1121420968.ms 2015/145/1121420968/1121420968_145_selfcal_joint_03.bin
#applysolutions 2015/145/1121420968/1121420968.ms 2015/145/1121420968/1121420968_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2015/145/1121507392/1121507392.ms 2015/145/1121507392/1121507392_145_selfcal_joint_03.bin
#applysolutions 2015/145/1121507392/1121507392.ms 2015/145/1121507392/1121507392_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2015/145/1121593824/1121593824.ms 2015/145/1121593824/1121593824_145_selfcal_joint_03.bin
#applysolutions 2015/145/1121593824/1121593824.ms 2015/145/1121593824/1121593824_145_selfcal_joint_03.bin
#calibrate -minuv 60 -ch 1  2015/145/1121680256/1121680256.ms 2015/145/1121680256/1121680256_145_selfcal_joint_03.bin
#applysolutions 2015/145/1121680256/1121680256.ms 2015/145/1121680256/1121680256_145_selfcal_joint_03.bin
