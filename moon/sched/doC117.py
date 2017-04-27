#!/usr/bin/python

import ephem
import math
import time
import sys
import os
import mwasources
import track_srcs

def usage():
    print """Scripty to run a bunch of single_observations at various freqs for a tracked C117 scan"""
    print "Usage:"
    print "ForA.py <date> <n_hours>"
    print "\t<date> in YYYY/MM/DD for start day"
    print "\t<n_hours> is how many hours total"
    sys.exit(1)

# parse command-line args
if (len(sys.argv) < 3): usage();

src = mwasources.Source('0924-2201',141.0,-22.0,freqs=[175],dwelltime_sec=296)
hyda= mwasources.Source('HydA',139.5236198,-12.0955426,iscal=True,known=True,freqs=[175],dwelltime_sec=56)

sched_secs = 3600*int(sys.argv[2])     # seconds of schedule time to fill
gap_secs = 4            # gap in seconds between scans
lst_starthour = src.ra_degs/15.0 - sched_secs/2.0/3600.0    # hours
# special hack for epsilon
lst_starthour -= 1.0

# write actual observation commands to a text file
outfile = "C117_cmdfile.txt"
fpcmd = open(outfile, 'w')

timenow = time.gmtime()
site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
site.date = ephem.Date(sys.argv[1])
lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
# note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
utc_starthour = (lst_starthour - lst_utc0*12.0/math.pi)/1.00278551532
if utc_starthour < 0: utc_starthour += 24.0

print >> fpcmd, "# Start day: "+str(site.date)
print >> fpcmd, "# LST at UTC 0 on start date: "+str(lst_utc0)
print >> fpcmd, "# Desired start LST: "+str(lst_starthour)
print >> fpcmd, "# Time offset in hours: "+str(utc_starthour)
localtime_start=utc_starthour+8.
if localtime_start > 24.0: localtime_start -= 24.0
print >> fpcmd, "# Approx local time at start: "+str(localtime_start)

datecomps=sys.argv[1].split('/') # split start date into YYYY MM and DD

starttime = time.mktime((int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), int(utc_starthour), int((utc_starthour-int(utc_starthour))*60), 0, 0, 0, 0))
endtime = starttime + sched_secs
delay = starttime - time.mktime(timenow)
if delay > 0:
    print "Will schedule observations starting "+str(delay)+" seconds from now."
    thetime = starttime
else:
    print "Desired start time passed "+str(-delay)+" seconds ago."
    sys.exit(1)

track_srcs.trackSources([src,hyda],lst_starthour,starttime,sched_secs,gap_secs,'C117',fpcmd)

