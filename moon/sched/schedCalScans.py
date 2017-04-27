#!/usr/bin/python
# A schedule calibration scans on a UT day at sunset and the following sunrise
# Randall Wayth. Sep 2014.

import os,sys,math,ephem,logging
import mwasources,track_srcs

logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger(__name__)   # default logger level is WARNING
side_per_sol = 24.0/23.93447    # sidereal days per solar days.
proj_id="D0006"
cal_dwell=180
cent_freqs=['69','93','121','145','169']
picket_freqs="62;63;69;70;76;77;84;85;93;94;103;104;113;114;125;126;139;140;153;154;169;170;187;188"

total_time_hrs=(len(cent_freqs)*cal_dwell+30.0)/3600.0

def usage():
    print >> sys.stderr,"""Scripty to schedule the evening and morning cal scans on a date"""
    print >> sys.stderr,"Usage:"
    print >> sys.stderr,"schedCalScans <date>"
    print >> sys.stderr,"\t<date> is in YYYY/MM/DD"
    sys.exit(1)

if __name__ == "__main__":

    logger.setLevel(logging.DEBUG)
    # parse command-line args
    if (len(sys.argv) < 2): usage();

    allsrcs = mwasources.LoadSources()
    calsrcs = mwasources.CalSources(allsrcs)

    # find the time of the sunset and following sunrise on a day
    site = ephem.Observer()
    site.long = '116:40:14.93'
    site.lat = '-26:42:11.95'
    site.date = sys.argv[1]
    lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
    # note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
    logger.info("LST at 0 hours UT on "+sys.argv[1]+" is "+str(lst_utc0))
    sun = ephem.Sun()
    sun.compute(site)
    settime=site.next_setting(sun)  # returns a Date object
    logger.info("Sun sets at "+str(settime)+" UT")
    risetime=site.next_rising(sun)  # returns a Date object
    logger.info("Sun rise at "+str(risetime)+" UT")

    # calc LST at sunset
    lst_sunset_hrs = (settime - site.date)/ephem.hour*side_per_sol + lst_utc0*12.0/math.pi
    if lst_sunset_hrs > 24.0: lst_sunset_hrs -= 24.0
    logger.info("LST at sunset "+str(lst_sunset_hrs)+" hours")
    # calc LST at sunrise
    lst_sunrise_hrs = (risetime - site.date)/ephem.hour*side_per_sol + lst_utc0*12.0/math.pi
    if lst_sunrise_hrs > 24.0: lst_sunrise_hrs -= 24.0
    logger.info("LST at sunrise "+str(lst_sunrise_hrs)+" hours")

    sunset_cals = mwasources.FindCal(calsrcs,lst_sunset_hrs*15.0)
    sunrise_cals =mwasources.FindCal(calsrcs,lst_sunrise_hrs*15.0)
    
    if len(sunset_cals) == 0:
        logger.warning("WARNING: Unable to find calibrator for sunset")
    else:
        # get the first item in the cals list, without the zenith angle arg
        sunset_cal = sunset_cals[0][0]
        logger.info("Top cal source at sunset: "+str(sunset_cal.name))

    if len(sunrise_cals) == 0:
        logger.warning("WARNING: Unable to find calibrator for sunrise")
    else:
        # get the first item in the cals list, without the zenith angle arg
        sunrise_cal = sunrise_cals[0][0]
        logger.info("Top cal source at sunrise: "+str(sunrise_cal.name))

    sunset_cmd ="./track_a_source.py -m 40_2 -d %d %s 0.01 %.3f %s %s" % (cal_dwell,sys.argv[1],lst_sunset_hrs,sunset_cal.name,proj_id)
    sunrise_cmd="./track_a_source.py -m 40_2 -d %d %s 0.01 %.3f %s %s" % (cal_dwell,sys.argv[1],lst_sunrise_hrs,sunrise_cal.name,proj_id)
    logger.info("Sunset track command: "+sunset_cmd)
    logger.info("Sunrise track command: "+sunrise_cmd)
    lst_sunset_hrs += total_time_hrs*side_per_sol;
    lst_sunrise_hrs += total_time_hrs*side_per_sol;
    sunset_cmd_picket ="./track_a_source.py -m 40_2 -f '%s' -d %d %s 0.01 %.3f %s %s" % (picket_freqs,cal_dwell,sys.argv[1],lst_sunset_hrs,sunset_cal.name,proj_id)
    sunrise_cmd_picket="./track_a_source.py -m 40_2 -f '%s' -d %d %s 0.01 %.3f %s %s" % (picket_freqs,cal_dwell,sys.argv[1],lst_sunrise_hrs,sunrise_cal.name,proj_id)
    logger.info("Sunset picket track command: "+sunset_cmd_picket)
    logger.info("Sunrise picket track command: "+sunrise_cmd_picket)
    os.system(sunset_cmd)
    os.system(sunrise_cmd)
    os.system(sunset_cmd_picket)
    os.system(sunrise_cmd_picket)

