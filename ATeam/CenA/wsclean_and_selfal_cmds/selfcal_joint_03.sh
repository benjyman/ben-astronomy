#!/usr/bin/env bash
calibrate -minuv 60 -ch 4  2018/1199663088/1199663088.ms   2018/1199663088/1199663088_selfcal_joint_04.bin
applysolutions 2018/1199663088/1199663088.ms 2018/1199663088/1199663088_selfcal_joint_04.bin
calibrate -minuv 60 -ch 2  2015/1117031728/1117031728.ms 2015/1117031728/1117031728_selfcal_joint_04.bin
applysolutions 2015/1117031728/1117031728.ms 2015/1117031728/1117031728_selfcal_joint_04.bin