#!/bin/bash
# script to image AAVS-1 sun data with miriad

if [ $# -lt 1 ] ; then
    echo "Usage: $0 basename subarray "
    echo "  basename is the uvfits file name without the .uvfits"
    echo "  subarray can be centre,north,south,east,west, or left blank"
    exit 0
fi

basename=$1
subarray=$2

t1="18apr05:04:00:00"
t2="18apr05:04:02:00"

select_time_string="select=time(${t1},${t2})"

#sort out the subarrays
flagged_ant_filename="/data/code/git/ben-astronomy/miriad/AAVS1_SelectedCentralArea_FlaggedAntID.txt"

#read the files to get the antenna arrays
readarray  flagged_ants < ${flagged_ant_filename}
#IFS=$'\n' flagged_ants=($(cat ${flagged_ant_filename}))

flagged_ant_string=""
for i in "${flagged_ants[@]}"
do
   #echo $i
   flagged_ant_string="$flagged_ant_string,${i//[!0-9]/}"
   # do whatever on $i
done

flagged_ant_string_1="select=antennae(${flagged_ant_string:1:30})"
flagged_ant_string_2="select=antennae(${flagged_ant_string:32:47})"
flagged_ant_string_3="select=antennae(${flagged_ant_string:80:45})"
flagged_ant_string_4="select=antennae(${flagged_ant_string:126:59})"
flagged_ant_string_5="select=antennae(${flagged_ant_string:186:59})"
flagged_ant_string_6="select=antennae(${flagged_ant_string:246:59})"
flagged_ant_string_7="select=antennae(${flagged_ant_string:306:59})"
#flagged_ant_string_8="select=antennae(${flagged_ant_string:218:50})"

#echo ${flagged_ant_string_1}
#echo ${flagged_ant_string_2}
#echo ${flagged_ant_string_3}
#echo ${flagged_ant_string_4}
#echo ${flagged_ant_string_5}
#echo ${flagged_ant_string_6}
#echo ${flagged_ant_string_7}
#echo ${flagged_ant_string_8}

uvfits_name="/md0/AAVS-1/data/${basename}.uvfits"
uvdata="/md0/AAVS-1/data/${basename}.uv"

#read in the data
#fits in=${uvfits_name} out=${uvdata}  op=uvin

if [ ${subarray} == "centre" ] ; then
   echo ${subarray}
   uvdata_name="/md0/AAVS-1/data/${basename}_centre.uv"
   #cp -r ${uvdata} ${uvdata_name}

   #uvflag vis=${uvdata_name} ${flagged_ant_string_1} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_2} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_3} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_4} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_5} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_6} flagval=flag
   #uvflag vis=${uvdata_name} ${flagged_ant_string_7} flagval=flag

   uvdata=${uvdata_name} 
fi


#calibrate the data
#mfcal vis=${uvdata} flux=10000,0.150,0 refant=2 interval=1 stokes=xx,yy,xy,yx

invert vis=${uvdata} map=${basename}_xx.map beam=${basename}.beam options=double imsize=512 stokes=xx robust=-0.5 cell=1800 ${select_time_string}
clean map=${basename}_xx.map beam=${basename}.beam out=${basename}_xx.model niters=2000
restor model=${basename}_xx.model  beam=${basename}.beam map=${basename}_xx.map out=${basename}_xx.restor

invert vis=${uvdata} map=${basename}_yy.map  options=double imsize=512 stokes=yy robust=-0.5 cell=1800 ${select_time_string}
clean map=${basename}_yy.map beam=${basename}.beam out=${basename}_yy.model niters=2000
restor model=${basename}_yy.model  beam=${basename}.beam map=${basename}_yy.map out=${basename}_yy.restor
