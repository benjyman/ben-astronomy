python /data/code/git/ben-astronomy/moon/processing_scripts/namorrodor/generate_mwac_qimage_concat_ms.py --cotter --track_off_moon=/data/moon/2017/track_off_moon_20150926_93.txt --no_pbcorr --tagname=20150929_moon_93  --imsize=2048 --pol="xx,xy,yx,yy"  --wsclean_options=" -niter 0  -datacolumn CORRECTED_DATA  -scale 0.0085 -weight uniform  -smallinversion -channelsout 24 -make-psf " /data/moon/2017/20150929_off_moon1_93.txt    