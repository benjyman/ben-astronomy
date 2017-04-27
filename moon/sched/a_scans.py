#!/usr/bin/python

# A scripty to track a list of sources over a nominated set of freqs (can be
# different for each source) to start at a specified time and run for a 
# specified time
# Randall Wayth. April 2011.

import ephem
import math
import time
import sys
import os
import mwasources
    
def e2h(ha,dec,lat):
    """Convert equatorial ha/dec coords (in hours, degrees) to az,el (degs) give an observer latitute (degs). Returns (az,el)"""
    ha_rad = ha*15.0*math.pi/180.0
    dec_rad = dec*math.pi/180.0
    lat_rad = lat*math.pi/180.0
    sh = math.sin(ha_rad)
    ch = math.cos(ha_rad)
    sd = math.sin(dec_rad)
    cd = math.cos(dec_rad)
    sp = math.sin(lat_rad)
    cp = math.cos(lat_rad)
    x = - ch * cd * sp + sd * cp;
    y = - sh * cd;
    z = ch * cd * cp + sd * sp;
    r = math.sqrt(x*x + y*y)
    a = math.atan2( y, x )
    if a < 0.0:
        az = a + 2.0*math.pi
    else:
        az = a
    el = math.atan2( z, r );
    return (az*180.0/math.pi,el*180.0/math.pi)


def usage():
    print """Scripty to run a bunch of single_observations for sources at various freqs and dwell times"""
    print "Usage:"
    print "a_scans <date> <LST> <n_hours>"
    print "\t<date> in YYYY/MM/DD for start day"
    sys.exit(1)

allsrcs = mwasources.LoadSources()
srclist = mwasources.CalSources(allsrcs)

# parse command-line args
if (len(sys.argv) < 2): usage();

lst_starthour = float(sys.argv[2])    # hours
sched_secs = int(3600*float(sys.argv[3]))     # seconds of schedule time to fill
gap_secs = 8             # gap in seconds between scans
int_time=2
freqres=40

side_per_sol = 24.0/23.9344696    # sidereal days per solar days.
timenow = time.gmtime()
site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
site.date = ephem.Date(sys.argv[1])
lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
# note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
utc_starthour = (lst_starthour - lst_utc0*12.0/math.pi)/side_per_sol
if utc_starthour < 0: utc_starthour += 24.0/side_per_sol

datecomps=sys.argv[1].split('/') # split start date into YYYY MM and DD

# write actual observation commands to a text file
outfile = "%d%02d%02d_cmdfile_cals.txt" % (int(datecomps[0]), int(datecomps[1]), int(datecomps[2]))
fpcmd = open(outfile, 'w')

print >> fpcmd, "# Start day: "+str(site.date)
print >> fpcmd, "# LST at UTC 0 on start date: "+str(lst_utc0)
print >> fpcmd, "# Desired start LST: "+str(lst_starthour)
print >> fpcmd, "# Time offset in hours: "+str(utc_starthour)
print >> fpcmd, "# Approx local time at start: "+str(utc_starthour+8.)

# sanity check
target = ephem.FixedBody()
target._ra = sys.argv[2]
target._dec=site.lat
target._epoch=site.date
target_transit = site.next_transit(target)
print >> fpcmd, "# Target RA transits at "+str(target_transit)

starthr =int(utc_starthour)
startmin=int((utc_starthour-starthr)*60)
startsec=int(((utc_starthour-starthr)*60 - startmin)*60)
starttime = time.mktime((int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), starthr, startmin, startsec, 0, 0, 0))
endtime = starttime + sched_secs
delay = starttime - time.mktime(timenow)
if delay > 0:
    print "Will schedule observations starting "+str(delay)+" seconds from now."
    thetime = starttime
else:
    print >> sys.stderr, "Desired start time passed "+str(-delay)+" seconds ago."
    sys.exit(1)

# set the correlator mode 2 mins before start if necessary
corrmode_startime=time.localtime(starttime-60)
str_starttime="%s-%s-%s,%s:%s:%s" % (corrmode_startime[0:6])
cmd="single_observation.py --obsname=CORR_MODE_%s_%s --az=0 --el=90 --starttime=%s --stoptime=++48s --freq=69,24 --mode=CORR_MODE_CHANGE --inttime=%s --freqres=%s" % (freqres,int_time,str_starttime,int_time,freqres)
print >> fpcmd,"# correlator mode setup\n"+cmd

lst_hours = lst_starthour
srcindex = 0
elapsed_time = 0.0
n_done =0
while thetime < endtime:
    curr_src = srclist[srcindex]
    for freq in curr_src.freqs:
        ha_hrs = lst_hours - curr_src.ra_degs/15.0
        (az,el) = e2h(ha_hrs,curr_src.dec_degs,float(site.lat)*180.0/math.pi)
        if (el < curr_src.el_limit):
            print >> fpcmd, "# Skipping source "+curr_src.name+" due to low el: "+str(el)+". Limit: "+str(curr_src.el_limit)
            continue
        else:
            n_done += 1
        # create time string suitable for single_observation. Here we need to
        # use localtime because the time is already in gmtime and we don't want
        # to add/subtract anything from it.
        timestr = time.strftime("%Y-%m-%d,%H:%M:%S",time.localtime(thetime))
        print >> fpcmd, "# Source "+curr_src.name+" for freq "+str(freq)+",at UTC time "+timestr+". Elapsed time: "+str(elapsed_time)
        print >> fpcmd, "#\tThe HA is: "+str(ha_hrs)+" hours. The LST is: "+str(lst_hours)+" hours. az/el: "+str(az)+","+str(el)

        # create single_observation command
        cmd = "single_observation.py --starttime="+timestr+" --stoptime=++"+str(curr_src.dwelltime_sec)+"s"
        cmd += " --freq="+str(freq)+",24 --obsname="+curr_src.name+"_"+str(freq)
        cmd += " --creator=Randall --useazel --usegrid= --project=D0006 --calibration --calibrator="+curr_src.name
        if curr_src.gain is not None:
            cmd += " --gain_control_value="+str(curr_src.gain)
        if curr_src.known == False:
            cmd += " --ra="+str(curr_src.ra_degs)+" --dec="+str(curr_src.dec_degs)
        else:
            cmd += " --source="+curr_src.name
        cmd += " --inttime="+str(int_time)+" --freqres="+str(freqres)

        print >> fpcmd, cmd


        # increment time
        thetime += curr_src.dwelltime_sec + gap_secs
        lst_hours += (curr_src.dwelltime_sec + gap_secs)*side_per_sol/3600.0
        elapsed_time += curr_src.dwelltime_sec + gap_secs
    # go to next source
    srcindex = (srcindex+1) % len(srclist)
    if srcindex == 0:
        if n_done > 0:
            # new cycle
            n_done = 0
        else:
            print >> fpcmd, "# Hmm. Seem to have skipped ALL the sources. Nothing interesting at this LST... Advancing 1 minute"
            thetime += 60.0
            lst_hours += side_per_sol/60.0
            elapsed_time += 60.0

