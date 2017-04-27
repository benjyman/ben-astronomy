#!/usr/bin/python

# A scripty to schedule the fixed dec drift scans for MWA commissioning
# Randall Wayth. July 2012.

import math, time, sys, os, getopt
import ephem
import numpy
import mwasources
import track_srcs

# set some defaults
# overall loop is over pointings (outer loop) then freq (inner loop).
# define centre coarse chans for freq coverage
#cent_freqs = ['71','95','121','145','169']    # can be changed on command-line
cent_freqs = ['69','93','121','145','169']    # can be changed on command-line
#define scan lengths in sec per freq
scan_times = numpy.zeros(len(cent_freqs)) + 116
gap_secs = 4            # gap in seconds between scans
startcal=False
endcal=False
cal_interval=0    # time in seconds between cal scans
cal_dwell=116     # dwell time for calibration scans
creator="Randall"
side_per_sol = 24.0/23.93447    # sidereal days per solar days.
MODES=('40_0.5','40_1','40_2','20_1','20_2','10_2')
mode=MODES[0]
ha_hrs=0.0

def usage():
    print >> sys.stderr,"""Scripty to schedule the fixed dec survey"""
    print >> sys.stderr,"Usage:"
    print >> sys.stderr,"fds [options] <date> <LST_hrs> <dec_degs> <proj_id> <hrs_to_sched>"
    print >> sys.stderr,"\toptions:"
    print >> sys.stderr,"\t-s   do a calibration scan at the start"
    print >> sys.stderr,"\t-e   do a calibration scan at the end"
    print >> sys.stderr,"\t-f   comma separated list of central coarse chan IDs to observse. Default: "+str(cent_freqs)
    print >> sys.stderr,"\t-u   user    specify a creator (no whitespace)"
    print >> sys.stderr,"\t-c   interval    do calibration scans at this interval (hours)"
    print >> sys.stderr,"\t-d   dwelltime    specify dwell time (seconds) including gap seconds. default: "+str(scan_times[0]+gap_secs)
    print >> sys.stderr,"\t-m   mode  specify correlator mode. Supported modes: "+str(MODES)
    print >> sys.stderr,"\t-h   hours  specify HA for off-meridian scans. Default: "+str(ha_hrs)
    print >> sys.stderr,"\t<date>\tin YYYY/MM/DD for start day"
    print >> sys.stderr,"\t<LST_hrs>\tstart LST in hours"
    print >> sys.stderr,"\t<decs>\tin degrees. Comma separated list of DECs"
    print >> sys.stderr,"\t      \t\tSweet spot DECs: 18.6, 9.7, 1.6, -5.9, -13.0, -19.9, -26.7, -33.5, -40.2, -47.5, -55.0, -63.0, -72.0"
    print >> sys.stderr,"\t<proj_id>\tneeds to be an existing MWA project. Use G0008 for GLEAM, G0001 for RSM."
    print >> sys.stderr,"\t<hrs_to_sched>\tnumber of hours to observe for"
    sys.exit(1)


def doCalScans(calsrcs,lst_hours,cal_dwell_sec,cal_gap_secs,timeres,freqres):
    """include a set of calibration scans choosing the best cal for the given LST"""
    print >> fpcmd, "# *** Dedicated calibration scans for LST: %.2f hrs" % (lst_hours)
    cals = mwasources.FindCal(calsrcs,lst_hours*15.0)   # this returns a list of tuples where each tuple is (mwasource,za)
    if len(cals) == 0:
        print >> sys.stderr, "WARNING: Unable to find calibrator for LST "+str(lst_hours)
    else:
        # send the first item in the cals list, without the zenith angle arg
        track_srcs.trackSources([x[0] for x in cals], lst_hours, thetime, 1, cal_gap_secs, proj_id, fpcmd=fpcmd, freqs_override=cent_freqs, name_override=creator, dwell_override_secs=cal_dwell,timeres=timeres,freqres=freqres)


# parse command-line args
if (len(sys.argv) < 5): usage();

opts,args = getopt.getopt(sys.argv[1:], 'esf:u:d:c:m:h:')

for o, a in opts:
        if o == "-s":
            startcal = True
        elif o == "-e":
            endcal = True;
        elif o == "-f":
            cent_freqs = a.strip().split(',')
            print >> sys.stderr, "Overriding default freq list with "+str(cent_freqs)
        elif o == "-u":
            creator = a
            scan_times[0] = int(a)-gap_secs
        elif o == "-c":
            cal_interval = int(float(a)*3600)
        elif o == "-m":
            assert a in MODES, "Unknown mode %r" % a
            mode = a
        elif o == "-h":
            ha_hrs = float(a)
            print >> sys.stderr,"Using HA of "+str(ha_hrs)+" hours."
        else:
            assert False, "unhandled option"

lst_starthour = float(args[1])    # hours
decs = numpy.array(map(float,args[2].split(',')))  # degrees
proj_id = args[3]
sched_secs=3600*float(args[4])
cal_gap = (int(cal_dwell/60) + 1)*60 - cal_dwell    # gap between cal scans to round up to integer minutes

# extract freq and time res from mode
m = mode.split('_')
freqres = m[0]
timeres = m[1]

timenow = time.gmtime()
site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
site.date = args[0]
lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
# note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
utc_starthour = (lst_starthour - lst_utc0*12.0/math.pi)/side_per_sol
if utc_starthour < 0: utc_starthour += 24.0/side_per_sol

datecomps=args[0].split('/') # split start date into YYYY MM and DD

# write actual observation commands to a text file
outfile = "%d%02d%02d_cmdfile_FDS%s_LST%s.txt" % (int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), args[2],args[1])
fpcmd = open(outfile, 'w')

print >> fpcmd, "# Start day: "+str(site.date)
print >> fpcmd, "# LST at UTC 0 on start date: "+str(lst_utc0*12.0/math.pi)
print >> fpcmd, "# Desired start LST: "+str(lst_starthour)
print >> fpcmd, "# UTC start hour: "+str(utc_starthour)
print >> fpcmd, "# Approx local time at start: "+str(utc_starthour+8.)
print >> fpcmd, "# Observing coarse chan IDs: "+str(cent_freqs)

# sanity check
target = ephem.FixedBody()
target._ra = args[1]
target._dec=site.lat
target._epoch=site.date
target_transit = site.next_transit(target)
print >> fpcmd, "# Target RA transits at "+str(target_transit)

# calc how long until schedule starts
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
    print "Desired start time passed "+str(-delay)+" seconds ago."
    sys.exit(1)

# define DEC ranges for meridian scans
zas  = decs*0.0
azs  = decs*0.0
for i in range(len(decs)):
    (az,el) = track_srcs.e2h(ha_hrs,decs[i],site.lat*180/math.pi)
    zas[i] = 90.0-el
    azs[i] = az
    print >> fpcmd, "# Az/El for Dec "+str(decs[i])+", HA "+str(ha_hrs)+" (hours) is "+str(az)+","+str(el)
# gain = 1.0    # use default

# calculate the total time of the scans over all DECs and freqs:
stripe_time_hrs = 0.0
for dec in decs:
  for i in range(len(cent_freqs)):
    stripe_time_hrs += (scan_times[i]+gap_secs)/3600.0
print >> fpcmd, "# Total time for all freqs: "+str(stripe_time_hrs)+" hours"

calsrcs = mwasources.CalSources(mwasources.LoadSources())
print >> fpcmd, "# There are "+str(len(calsrcs))+" cal sources"

# set some counters
stripe_index=0
lst_hours = lst_starthour
elapsed_time = 0
time_since_cal_secs = 0

# set the correlator mode 2 mins before start if necessary
corrmode_startime=time.localtime(starttime-60)
str_starttime="%s-%s-%s,%s:%s:%s" % (corrmode_startime[0:6])
cmd="single_observation.py --creator=%s --obsname=CORR_MODE_%s_%s --az=0 --el=90 --starttime=%s --stoptime=++48s --freq=%s,24 --mode=CORR_MODE_CHANGE --inttime=%s --freqres=%s" % (creator, freqres,timeres,str_starttime,str(cent_freqs[0]),timeres,freqres)
print >> fpcmd,"# correlator mode setup\n"+cmd

# reduce the end time if a cal scan at the end is required
if endcal:
    # time to run a cal is all freqs for a single pointing
    cal_time_secs = len(cent_freqs)*(cal_dwell+cal_gap)
    endtime -= cal_time_secs

# do a start calibration if desired
if startcal:
    doCalScans(calsrcs,lst_hours,cal_dwell,cal_gap,timeres,freqres)
    time_inc_sec=len(cent_freqs)*(cal_dwell+cal_gap)
    thetime += time_inc_sec
    lst_hours += time_inc_sec*side_per_sol/3600.0
    if lst_hours > 24.0: lst_hours -= 24.0
    print >> fpcmd, "# Incremented start time by "+str(time_inc_sec)+" seconds for cal scans"

while thetime < endtime:
  # calc the LST in the middle and end of the stripe
  lst_end_hrs = lst_hours + stripe_time_hrs*side_per_sol
  if lst_end_hrs >= 24.0: lst_end_hrs -= 24.0
  lst_mid_hrs = lst_hours + stripe_time_hrs*side_per_sol/2.0
  if lst_mid_hrs >= 24.0: lst_mid_hrs -= 24.0

  lst_stripe_range = numpy.array([lst_hours,lst_end_hrs])
  print >> fpcmd, "# Stripe %d. LST range: %f %f. Mid: %f" % (stripe_index,lst_stripe_range[0],lst_stripe_range[1],lst_mid_hrs)

  for p in range(len(decs)):
    # set the pointing
    az = azs[p]
    el = 90.0 - zas[p]
    ra_hrs = lst_hours
    dec = decs[p]
    print >> fpcmd, "# Pointing %d. RA: %g, DEC: %g, Az: %g, El: %g" %(p,ra_hrs,dec,az,el)
    beam_srcs = mwasources.InBeam(calsrcs,ra_hrs*15.0,dec,10.0)
    if len(beam_srcs) > 0:
        print >> fpcmd, "# There are "+str(len(beam_srcs))+" cal sources in this pointing."
        for s in beam_srcs:
            print >> fpcmd, "# "+str(s)+" Distance: "+str(s.AngularSep(ra_hrs*15.0,dec)*180.0/math.pi)
    for f in range(len(cent_freqs)):
        # create time string suitable for single_observation. Here we need to
        # use localtime because the time is already in gmtime and we don't want
        # to add/subtract anything from it.
        timestr = time.strftime("%Y-%m-%d,%H:%M:%S",time.localtime(thetime))
        print >> fpcmd, "# The LST is: %g. Coarse chan: %s. UTC: %s. Elapsed time: %g" % (lst_hours,cent_freqs[f],timestr,elapsed_time)

        # create single_observation command
        cmd = "single_observation.py"
        cmd += " --starttime="+timestr+" --stoptime=++"+str(int(scan_times[f]))+"s"
        cmd += " --freq="+str(cent_freqs[f])+",24"
        cmd += " --obsname=FDS_DEC"+str(decs[p])+"_"+str(cent_freqs[f])
        cmd += " --creator="+creator
        cmd += " --useazel --usegrid= "
        cmd += " --project="+proj_id 
#        cmd += " --gain_control_value="+str(gain) # use default
        cmd += " --az="+str(az)+" --el="+str(el)
        cmd += " --inttime="+str(timeres)
        cmd += " --freqres="+str(freqres)
        if len(beam_srcs) > 0:
            cmd += " --calibration --calibrator="
            ind=0
            for s in beam_srcs:
                cmd += s.name
                if ind < len(beam_srcs)-1:
                    cmd += ","
                ind += 1
        print >> fpcmd, cmd

        # increment time
        this_time_inc = scan_times[f] + gap_secs
        thetime += this_time_inc
        lst_hours += (this_time_inc)*side_per_sol/3600.0
        if lst_hours >= 24.0: lst_hours -= 24.0
        elapsed_time += this_time_inc 
        time_since_cal_secs += this_time_inc
        if thetime > endtime: break

        if cal_interval > 0 and time_since_cal_secs >= cal_interval:
            print >> fpcmd, "# It has been %d seconds since a cal scan. Doing one now" %(time_since_cal_secs)
            doCalScans(calsrcs,lst_hours,cal_dwell,cal_gap,timeres,freqres)
            # need to skip as many pointings as there are DECs
            time_inc_sec=(cal_dwell + cal_gap)*len(cent_freqs)
            if len(decs) > 1:
                time_inc_sec *= len(decs)
                print >> fpcmd, "# delaying after cal scan to keep drift scans equally spaced for multiple DECs"
            thetime += time_inc_sec
            lst_hours += time_inc_sec*side_per_sol/3600.0
            time_since_cal_secs = 0
            elapsed_time += time_inc_sec
            print >> fpcmd, "# Cal scans done. New LST: %.2f hrs" % (lst_hours)

if endcal:
    doCalScans(calsrcs,lst_hours,cal_dwell,cal_gap,timeres,freqres)
    time_inc_sec=3600*stripe_time_hrs/len(decs)
    thetime += time_inc_sec
    lst_hours += time_inc_sec*side_per_sol/3600.0


