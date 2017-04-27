#!/usr/bin/python

import math,time,sys,os,getopt
import ephem
import numpy
import mwasources
import track_srcs

# set some defaults
# define centre coarse chans for freq coverage
gap_secs = 4            # gap in seconds between scans
freqs = ['69','93','121','145','169']    # can be changed on command-line
#define scan lengths in sec per freq
scan_times = numpy.zeros(10) + 116
startcal=False
endcal=False
cal_interval=0    # time in seconds between cal scans
cal_dwell=236     # dwell time for calibration scans
cal_dwell=112
creator=None
side_per_sol = 24.0/23.93447    # sidereal days per solar day.
MODES=('40_0.5','40_1','40_2','20_1','20_2','10_2')
mode=MODES[0]
doModeChange=True

def usage():
    print >> sys.stderr,"""Scripty to run a bunch of single_observations at various freqs to track a known source from mwasources"""
    print >> sys.stderr,"Usage:"
    print >> sys.stderr,sys.argv[0]+" [options ] <date> <n_hours> <lst_starthour> <srcname> <projid>"
    print >> sys.stderr,"\t-s   do a calibration scan at the start"
    print >> sys.stderr,"\t-e   do a calibration scan at the end"
    print >> sys.stderr,"\t-f   comma separated list of central coarse chan IDs to observse. Default: "+str(freqs)
    print >> sys.stderr,"\t     OR specify list of semicolon separated coarse chan IDs for non-contiguous (picket fence) mode"
    print >> sys.stderr,"\t-u   user    specify a creator (no whitespace). Default: username"
    print >> sys.stderr,"\t-c   interval    do calibration scans at this interval (hours). No default"
    print >> sys.stderr,"\t-d   dwelltime    specify dwell time (seconds) including gap seconds. default: "+str(scan_times[0]+gap_secs)
    print >> sys.stderr,"\t-m   mode  specify correlator mode. Supported modes: "+str(MODES)
    print >> sys.stderr,"\t-M   do not issue a mode change command (assumes correlator is already in desired mode)"
    print >> sys.stderr,"\t-C   track coordinate specified as RA,DEC (decimal hours,decimal degrees). Will be overridden if srcname is known"
    print >> sys.stderr,"\tdate           in YYYY/MM/DD for start day"
    print >> sys.stderr,"\tn_hours        is how many hours total including cal scans"
    print >> sys.stderr,"\tlst_starthour  is the desired LST [hour] to start obs or 'now' for now"
    print >> sys.stderr,"\tsrcname        is the name of the source as defined in mwasources"
    print >> sys.stderr,"\tprojid         is an existing project id. Use C001 if in doubt"
    sys.exit(1)

def doCalScans(calsrcs,lst_hours,cal_dwell_sec,cal_gap_secs,timeres,freqres):
    """include a set of calibration scans choosing the best cal for the given LST"""
    print >> fpcmd, "# *** Dedicated calibration scans for LST: %.2f hrs" % (lst_hours)
    cals = mwasources.FindCal(calsrcs,lst_hours*15.0)   # this returns a list of tuples where each tuple is (mwasource,za)
    if len(cals) == 0:
        print >> sys.stderr, "WARNING: Unable to find calibrator for start of scan"
    else:
        # send the first item in the cals list, without the zenith angle arg
        track_srcs.trackSources([x[0] for x in cals], lst_hours, thetime, 1, cal_gap_secs, proj_id, fpcmd=fpcmd, freqs_override=freqs, name_override=creator,timeres=timeres,freqres=freqres)

# parse command-line args
if (len(sys.argv) < 5): usage();
opts,args = getopt.getopt(sys.argv[1:], 'esMf:u:d:c:m:C:')

for o, a in opts:
        if o == "-s":
            startcal = True
        elif o == "-e":
            endcal = True;
        elif o == "-M":
            doModeChange = False;
        elif o == "-f":
            freqs = a.strip().split(',')
            print >> sys.stderr, "Overriding default freq list with "+str(freqs)
        elif o == "-u":
            creator = a
        elif o == "-d":
            scan_times[:] = int(a)-gap_secs
        elif o == "-c":
            cal_interval = int(float(a)*3600)
        elif o == "-m":
            assert a in MODES, "Unknown mode %r" % a
            mode = a
        elif o == '-C':
            radec = a.strip().split(',')
        else:
            assert False, "unhandled option"

proj_id=args[4]
srcname=args[3]
sched_secs = 3600*float(args[1])     # seconds of schedule time to fill
cal_gap = (int(cal_dwell/60) + 1)*60 - cal_dwell    # gap between cal scans to round up to integer minutes

timenow = time.gmtime()
site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
lst_now = site.sidereal_time()
site.date = ephem.Date(args[0])
lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
if args[2]=='now':
    lst_starthour = lst_now*12.0/math.pi + 32./3600.
else:
    lst_starthour = float(args[2])

# note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
utc_starthour = (lst_starthour - lst_utc0*12.0/math.pi)/side_per_sol
if utc_starthour < 0: utc_starthour += 24.0/side_per_sol

# load all known sources
allsrcs=mwasources.LoadSources()
# find the source
foundsrc=False
for src in allsrcs:
  if src.name==srcname:
    foundsrc=True
    break
if foundsrc==False:
  print >> sys.stderr, "Unknown source "+srcname
  if radec is not None:
    print >> sys.stderr,"Creating manually specified source with ra "+radec[0]+", dec "+radec[1]
    src=mwasources.Source(srcname,float(radec[0])*15,float(radec[1]))
  else:
    sys.exit(1)

# extract freq and time res from mode
m = mode.split('_')
freqres = m[0]
timeres = m[1]

datecomps=args[0].split('/') # split start date into YYYY MM and DD

# write actual observation commands to a text file
outfile = "%d%02d%02d_cmdfile_%s_LST%s.txt" % (int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), args[3],args[2])
fpcmd = open(outfile, 'w')

print >> fpcmd, "# Start day: "+str(site.date)
print >> fpcmd, "# LST at UTC 0 on start date: "+str(lst_utc0*12.0/math.pi)
print >> fpcmd, "# Desired start LST: "+str(lst_starthour)
print >> fpcmd, "# UTC start hour: "+str(utc_starthour)
localtime_start=utc_starthour+8.
if localtime_start > 24.0: localtime_start -= 24.0
print >> fpcmd, "# WST time at start: "+str(localtime_start)
print >> fpcmd, "# Will track for "+str(sched_secs)+" seconds"

# sanity check
target = ephem.FixedBody()
target._ra = args[2]
target._dec=site.lat
target._epoch=site.date
target_transit = site.next_transit(target)
print >> fpcmd, "# Target RA transits at "+str(target_transit)

calsrcs = mwasources.CalSources(allsrcs)
print >> fpcmd, "# There are "+str(len(calsrcs))+" cal sources"

starthr =int(utc_starthour)
startmin=int((utc_starthour-starthr)*60)
startsec=int(((utc_starthour-starthr)*60 - startmin)*60)
starttime = time.mktime((int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), starthr, startmin, startsec, 0, 0, 0))
delay = starttime - time.mktime(timenow)
if delay > 0:
    print "Will schedule observations starting "+str(delay)+" seconds from now."
    thetime = starttime
else:
    print "Desired start time passed "+str(-delay)+" seconds ago."
    sys.exit(1)

stripe_time_secs=0 # total time for all freqs
for i in range(len(freqs)):
  stripe_time_secs += gap_secs+scan_times[i]
stripe_time_hrs = stripe_time_secs/3600.0
cal_stripe_time_hrs = len(freqs)*(cal_dwell+cal_gap)/3600.0
print >> fpcmd, "# Total time for all freqs: "+str(cal_stripe_time_hrs)+" hours for cals, "+str(stripe_time_hrs)+" hours for sources"

# set some counters
stripe_index=0
lst_hours = lst_starthour
elapsed_time = 0
time_since_cal_secs = 0

# set the correlator mode shortly before start
if doModeChange:
    corrmode_startime=time.localtime(starttime-24)
    str_starttime="%s-%s-%s,%s:%s:%s" % (corrmode_startime[0:6])
    c_str=''
    if creator:
        c_str="--creator="+str(creator)
    cmd="single_observation.py %s --obsname=CORR_MODE_%s_%s --az=0 --el=90 --starttime=%s --stoptime=++16s --freq=69,24 --mode=CORR_MODE_CHANGE --inttime=%s --freqres=%s --project=%s" % (c_str,freqres,timeres,str_starttime,timeres,freqres,proj_id) 
    print >> fpcmd,"# correlator mode setup\n"+cmd

# reduce the end time if a cal scan at the end is required
if endcal:
    # time to run a cal is all freqs for a single pointing
    sched_secs -= int(cal_stripe_time_hrs*3600)
    print >> fpcmd, "# Reducing total time by "+str(cal_stripe_time_hrs)+" hours for end cal scans. New schedule time: "+str(sched_secs)+" seconds"

# do a start calibration if desired
if startcal:
    doCalScans(calsrcs,lst_hours,cal_dwell,cal_gap,timeres,freqres)
    time_inc_sec=3600*cal_stripe_time_hrs
    thetime += time_inc_sec
    lst_hours += time_inc_sec*side_per_sol/3600.0
    if lst_hours > 24.0: lst_hours -= 24.0
    print >> fpcmd, "# Incremented start time by "+str(3600*cal_stripe_time_hrs)+" seconds for cal scans"

track_srcs.trackSources([src], lst_hours, thetime, sched_secs, gap_secs, proj_id, fpcmd, freqs, dwell_override_secs=scan_times[0],name_override=creator,timeres=timeres,freqres=freqres)

# advance the LST and wallclock time. Need to work out how much time trackSources takes.
# it will be an integer number of (sched_secs/dwell_override) rounded up
nsrcscans=int(sched_secs/stripe_time_secs)
nfracscan=(sched_secs/stripe_time_secs)-nsrcscans
if (nfracscan > 0): nsrcscans += 1
srctime=stripe_time_secs*nsrcscans
print >> fpcmd, "# There are "+str(nsrcscans)+" scans of time "+str(stripe_time_secs)+" for the source totalling "+str(srctime)+" seconds"
thetime += srctime
lst_hours += srctime*side_per_sol/3600.0
if lst_hours > 24.0: lst_hours -= 24.0

if endcal:
    doCalScans(calsrcs,lst_hours,cal_dwell,cal_gap,timeres,freqres)

