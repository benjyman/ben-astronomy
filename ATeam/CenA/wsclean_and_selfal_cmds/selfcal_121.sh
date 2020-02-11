#!/usr/bin/env bash

#flagmwa   2018/121/1200777400/1200777400.ms 0
#flagmwa   2018/121/1200690952/1200690952.ms 0
#flagmwa   2018/121/1200604504/1200604504.ms 0
#flagmwa   2018/121/1200518056/1200518056.ms 0
#flagmwa   2018/121/1200431608/1200431608.ms 0
#flagmwa   2018/121/1200345160/1200345160.ms 0

#flagmwa   2015/121/1117031608/1117031608.ms 0
#flagmwa   2015/121/1121334352/1121334352.ms 0
#flagmwa   2015/121/1121325688/1121325688.ms 0
#flagmwa   2015/121/1121247928/1121247928.ms 0
#flagmwa   2015/121/1121321184/1121321184.ms 0
#flagmwa   2015/121/1121322984/1121322984.ms 0


#applysolutions 2018/1200431608/1200431608.ms 2018/1200518056/1200518056_selfcal_joint_06.bin
#applysolutions 2018/1200345160/1200345160.ms 2018/1200518056/1200518056_selfcal_joint_06.bin


#try to calibrate on model I gave to Gemma:
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200777400/1200777400.ms 2018/121/1200777400/1200777400_core_model.bin
#applysolutions 2018/121/1200777400/1200777400.ms 2018/121/1200777400/1200777400_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200690952/1200690952.ms 2018/121/1200690952/1200690952_core_model.bin
#applysolutions 2018/121/1200690952/1200690952.ms 2018/121/1200690952/1200690952_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200604504/1200604504.ms 2018/121/1200604504/1200604504_core_model.bin
#applysolutions 2018/121/1200604504/1200604504.ms 2018/121/1200604504/1200604504_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200518056/1200518056.ms 2018/121/1200518056/1200518056_core_model.bin
#applysolutions 2018/121/1200518056/1200518056.ms 2018/121/1200518056/1200518056_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200431608/1200431608.ms 2018/121/1200431608/1200431608_core_model.bin
#applysolutions 2018/121/1200431608/1200431608.ms 2018/121/1200431608/1200431608_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/121/1200345160/1200345160.ms 2018/121/1200345160/1200345160_core_model.bin
#applysolutions 2018/121/1200345160/1200345160.ms 2018/121/1200345160/1200345160_core_model.bin

#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1117031608/1117031608.ms 2015/121/1117031608/1117031608_core_model.bin
#applysolutions 2015/121/1117031608/1117031608.ms 2015/121/1117031608/1117031608_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1121334352/1121334352.ms 2015/121/1121334352/1121334352_core_model.bin
#applysolutions 2015/121/1121334352/1121334352.ms 2015/121/1121334352/1121334352_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1121325688/1121325688.ms 2015/121/1121325688/1121325688_core_model.bin
#applysolutions 2015/121/1121325688/1121325688.ms 2015/121/1121325688/1121325688_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1121247928/1121247928.ms 2015/121/1121247928/1121247928_core_model.bin
#applysolutions 2015/121/1121247928/1121247928.ms 2015/121/1121247928/1121247928_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1121321184/1121321184.ms 2015/121/1121321184/1121321184_core_model.bin
#applysolutions 2015/121/1121321184/1121321184.ms 2015/121/1121321184/1121321184_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/121/1121322984/1121322984.ms 2015/121/1121322984/1121322984_core_model.bin
#applysolutions 2015/121/1121322984/1121322984.ms 2015/121/1121322984/1121322984_core_model.bin

##selfcal
calibrate -minuv 60 -ch 1  2018/121/1200777400/1200777400.ms 2018/121/1200777400/1200777400_selfcal_joint_03.bin
applysolutions 2018/121/1200777400/1200777400.ms 2018/121/1200777400/1200777400_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2018/121/1200690952/1200690952.ms 2018/121/1200690952/1200690952_selfcal_joint_03.bin
applysolutions 2018/121/1200690952/1200690952.ms 2018/121/1200690952/1200690952_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2018/121/1200604504/1200604504.ms 2018/121/1200604504/1200604504_selfcal_joint_03.bin
applysolutions 2018/121/1200604504/1200604504.ms 2018/121/1200604504/1200604504_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2018/121/1200518056/1200518056.ms 2018/121/1200518056/1200518056_selfcal_joint_03.bin
applysolutions 2018/121/1200518056/1200518056.ms 2018/121/1200518056/1200518056_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2018/121/1200431608/1200431608.ms 2018/121/1200431608/1200431608_selfcal_joint_03.bin
applysolutions 2018/121/1200431608/1200431608.ms 2018/121/1200431608/1200431608_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2018/121/1200345160/1200345160.ms 2018/121/1200345160/1200345160_selfcal_joint_03.bin
applysolutions 2018/121/1200345160/1200345160.ms 2018/121/1200345160/1200345160_selfcal_joint_03.bin

calibrate -minuv 60 -ch 1  2015/121/1117031608/1117031608.ms 2015/121/1117031608/1117031608_selfcal_joint_03.bin
applysolutions 2015/121/1117031608/1117031608.ms 2015/121/1117031608/1117031608_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2015/121/1121334352/1121334352.ms 2015/121/1121334352/1121334352_selfcal_joint_03.bin
applysolutions 2015/121/1121334352/1121334352.ms 2015/121/1121334352/1121334352_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2015/121/1121325688/1121325688.ms 2015/121/1121325688/1121325688_selfcal_joint_03.bin
applysolutions 2015/121/1121325688/1121325688.ms 2015/121/1121325688/1121325688_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2015/121/1121247928/1121247928.ms 2015/121/1121247928/1121247928_selfcal_joint_03.bin
applysolutions 2015/121/1121247928/1121247928.ms 2015/121/1121247928/1121247928_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2015/121/1121321184/1121321184.ms 2015/121/1121321184/1121321184_selfcal_joint_03.bin
applysolutions 2015/121/1121321184/1121321184.ms 2015/121/1121321184/1121321184_selfcal_joint_03.bin
calibrate -minuv 60 -ch 1  2015/121/1121322984/1121322984.ms 2015/121/1121322984/1121322984_selfcal_joint_03.bin
applysolutions 2015/121/1121322984/1121322984.ms 2015/121/1121322984/1121322984_selfcal_joint_03.bin