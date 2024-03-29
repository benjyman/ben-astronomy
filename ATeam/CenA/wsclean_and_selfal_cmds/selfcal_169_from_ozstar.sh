#!/usr/bin/env bash

#flagmwa   2018/169/1200864208/1200864208.ms 0
#flagmwa   2018/169/1200777760/1200777760.ms 0
#flagmwa   2018/169/1200691312/1200691312.ms 0
#flagmwa   2018/169/1200604864/1200604864.ms 0
#flagmwa   2018/169/1200518416/1200518416.ms 0
#flagmwa   2018/169/1200431968/1200431968.ms 0

#flagmwa   2015/169/1117031848/1117031848.ms 0
#flagmwa   2015/169/1121594000/1121594000.ms 0
#flagmwa   2015/169/1121507576/1121507576.ms 0
#flagmwa   2015/169/1121421144/1121421144.ms 0
#flagmwa   2015/169/1121334712/1121334712.ms 0
#flagmwa   2015/169/1121248288/1121248288.ms 0


#applysolutions 2018/1200518416/1200518416.ms 2018/1200604864/1200604864_selfcal_joint_06.bin
#applysolutions 2018/1200431968/1200431968.ms 2018/1200604864/1200604864_selfcal_joint_06.bin


#try to calibrate on model I gave to Gemma:
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200864208/1200864208.ms 2018/169/1200864208/1200864208_core_model.bin
#applysolutions 2018/169/1200864208/1200864208.ms 2018/169/1200864208/1200864208_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200777760/1200777760.ms 2018/169/1200777760/1200777760_core_model.bin
#applysolutions 2018/169/1200777760/1200777760.ms 2018/169/1200777760/1200777760_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200691312/1200691312.ms 2018/169/1200691312/1200691312_core_model.bin
#applysolutions 2018/169/1200691312/1200691312.ms 2018/169/1200691312/1200691312_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200604864/1200604864.ms 2018/169/1200604864/1200604864_core_model.bin
#applysolutions 2018/169/1200604864/1200604864.ms 2018/169/1200604864/1200604864_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200518416/1200518416.ms 2018/169/1200518416/1200518416_core_model.bin
#applysolutions 2018/169/1200518416/1200518416.ms 2018/169/1200518416/1200518416_core_model.bin
#calibrate -minuv 60 -ch 1 -m CenA_core_wsclean_model.txt  -applybeam 2018/169/1200431968/1200431968.ms 2018/169/1200431968/1200431968_core_model.bin
#applysolutions 2018/169/1200431968/1200431968.ms 2018/169/1200431968/1200431968_core_model.bin

#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1117031848/1117031848.ms 2015/169/1117031848/1117031848_core_model.bin
#applysolutions 2015/169/1117031848/1117031848.ms 2015/169/1117031848/1117031848_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1121594000/1121594000.ms 2015/169/1121594000/1121594000_core_model.bin
#applysolutions 2015/169/1121594000/1121594000.ms 2015/169/1121594000/1121594000_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1121507576/1121507576.ms 2015/169/1121507576/1121507576_core_model.bin
#applysolutions 2015/169/1121507576/1121507576.ms 2015/169/1121507576/1121507576_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1121421144/1121421144.ms 2015/169/1121421144/1121421144_core_model.bin
#applysolutions 2015/169/1121421144/1121421144.ms 2015/169/1121421144/1121421144_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1121334712/1121334712.ms 2015/169/1121334712/1121334712_core_model.bin
#applysolutions 2015/169/1121334712/1121334712.ms 2015/169/1121334712/1121334712_core_model.bin
#calibrate -minuv 60 -ch 1  -m CenA_core_wsclean_model.txt -applybeam 2015/169/1121248288/1121248288.ms 2015/169/1121248288/1121248288_core_model.bin
#applysolutions 2015/169/1121248288/1121248288.ms 2015/169/1121248288/1121248288_core_model.bin

##selfcal
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200864208.ms from_ozstar/2018/169/1200864208_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200864208.ms from_ozstar/2018/169/1200864208_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200777760.ms from_ozstar/2018/169/1200777760_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200777760.ms from_ozstar/2018/169/1200777760_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200691312.ms from_ozstar/2018/169/1200691312_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200691312.ms from_ozstar/2018/169/1200691312_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200604864.ms from_ozstar/2018/169/1200604864_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200604864.ms from_ozstar/2018/169/1200604864_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200518416.ms from_ozstar/2018/169/1200518416_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200518416.ms from_ozstar/2018/169/1200518416_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2018/169/1200431968.ms from_ozstar/2018/169/1200431968_selfcal_joint_04.bin
applysolutions from_ozstar/2018/169/1200431968.ms from_ozstar/2018/169/1200431968_selfcal_joint_04.bin

calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1117031848.ms from_ozstar/2015/169/1117031848_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1117031848.ms from_ozstar/2015/169/1117031848_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1121594000.ms from_ozstar/2015/169/1121594000_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1121594000.ms from_ozstar/2015/169/1121594000_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1121507576.ms from_ozstar/2015/169/1121507576_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1121507576.ms from_ozstar/2015/169/1121507576_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1121421144.ms from_ozstar/2015/169/1121421144_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1121421144.ms from_ozstar/2015/169/1121421144_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1121334712.ms from_ozstar/2015/169/1121334712_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1121334712.ms from_ozstar/2015/169/1121334712_selfcal_joint_04.bin
calibrate -minuv 60 -ch 1  from_ozstar/2015/169/1121248288.ms from_ozstar/2015/169/1121248288_selfcal_joint_04.bin
applysolutions from_ozstar/2015/169/1121248288.ms from_ozstar/2015/169/1121248288_selfcal_joint_04.bin
