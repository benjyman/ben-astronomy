#!/usr/bin/env bash
#script calibrates data based on model already in ms (after joint deconvolution)
calibrate -minuv 60 2014B/1102865728/1102865728.ms 2014B/1102865728/1102865728_selfcal_joint_model_deeper_02.bin 
applysolutions 2014B/1102865728/1102865728.ms 2014B/1102865728/1102865728_selfcal_joint_model_deeper_02.bin
calibrate -minuv 60 2014B/1102865128/1102865128.ms 2014B/1102865128/1102865128_selfcal_joint_model_deeper_02.bin 
applysolutions 2014B/1102865128/1102865128.ms 2014B/1102865128/1102865128_selfcal_joint_model_deeper_02.bin
calibrate -minuv 60 2014B/1102864528/1102864528.ms 2014B/1102864528/1102864528_selfcal_joint_model_deeper_02.bin 
applysolutions 2014B/1102864528/1102864528.ms 2014B/1102864528/1102864528_selfcal_joint_model_deeper_02.bin
calibrate -minuv 60 2018A/1202816352/1202816352.ms 2018A/1202816352/1202816352_selfcal_joint_model_deeper_02.bin 
applysolutions 2018A/1202816352/1202816352.ms 2018A/1202816352/1202816352_selfcal_joint_model_deeper_02.bin
calibrate -minuv 60 2018A/1202815752/1202815752.ms 2018A/1202815752/1202815752_selfcal_joint_model_deeper_02.bin 
applysolutions 2018A/1202815752/1202815752.ms 2018A/1202815752/1202815752_selfcal_joint_model_deeper_02.bin
calibrate -minuv 60 2018A/1202815152/1202815152.ms 2018A/1202815152/1202815152_selfcal_joint_model_deeper_02.bin 
applysolutions 2018A/1202815152/1202815152.ms 2018A/1202815152/1202815152_selfcal_joint_model_deeper_02.bin
