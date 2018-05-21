#!/bin/bash
# script to add slant orthographic keyword headers to a miriad image

if [ $# -lt 1 ] ; then
    echo "Usage: $0 imagename LST_hours"
    echo "  imagename is a the miriad image"
    exit 0
fi

img=$1
lst_hrs=$2

lat_radian=-0.46606
dec_radian=`itemize in=${img} | grep crval2 | tr -d " " | cut -f 2 -d "="`
ra_radian=`itemize in=${img} | grep crval1 | tr -d " " | cut -f 2 -d "="`

ha_radian=`echo $lst_hrs $ra_radian | awk ' { print $1/3.81972 - $2 } ' `
# calc the zenith angle (Z) and the parallacitc angle (chi)
cosZ=`echo $lat_radian $dec_radian $ha_radian | awk ' { print sin($1)*sin($2) + cos($1)*cos($2)*cos($3) } ' `
# use trig identity tan(acos(x))=sqrt(1-x^2)/x
tanZA=`echo $cosZ | awk ' { print sqrt(1.0 - $1*$1)/$1 }' ` 
chi_radian=`echo $lat_radian $dec_radian $ha_radian | awk ' { print atan2(cos($1)*sin($3),(sin($1)*cos($2)-cos($1)*sin($2)*cos($3))) } ' `
xi=`echo $tanZA $chi_radian | awk ' { print $1*sin($2) } ' `
eta=`echo $tanZA $chi_radian | awk ' { print $1*cos($2) } ' `
#update the pv items in the image data
puthd in=${img}/pv1 value=$xi > /dev/null
puthd in=${img}/pv2 value=$eta > /dev/null

