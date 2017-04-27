#!/usr/bin/python

# A scripty to track a list of sources over a nominated set of freqs (can be
# different for each source) to start at a specified time and run for a 
# specified time
# Randall Wayth. July 2012.

import ephem
import math
import time
import sys
import os
import numpy
import mwasources

def usage():
    print >> sys.stderr,"""Scripty to schedule the meridian sky survey"""
    print >> sys.stderr,"Usage:"
    print >> sys.stderr,"mss <date> <LST>"
    print >> sys.stderr,"\t<date> in YYYY/MM/DD for start day"
    sys.exit(1)

# parse command-line args
if (len(sys.argv) < 2): usage();

lst_starthour = float(sys.argv[2])    # hours
sched_secs = 3600*12     # seconds of schedule time to fill
gap_secs = 16            # gap in seconds between scans
side_per_sol = 24.0/23.93447    # sidereal days per solar days.
timenow = time.gmtime()
site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
site.date = sys.argv[1]
lst_utc0 = site.sidereal_time()      # calc LST at 0 hours UTC on day of observations
# note: lst result is natively in radian, but will print in hour/min/seconds. Don't be caught out.
utc_starthour = (lst_starthour - lst_utc0*12.0/math.pi)/side_per_sol

# write actual observation commands to a text file
outfile = "mss_cmdfile.txt"
fpcmd = open(outfile, 'w')

print >> fpcmd, "# Start day: "+str(site.date)
print >> fpcmd, "# LST at UTC 0 on start date: "+str(lst_utc0)
print >> fpcmd, "# Desired start LST: "+str(lst_starthour)
print >> fpcmd, "# Time offset in hours: "+str(utc_starthour)
print >> fpcmd, "# Approx local time at start: "+str(utc_starthour+8.)

datecomps=sys.argv[1].split('/') # split start date into YYYY MM and DD

starttime = time.mktime((int(datecomps[0]), int(datecomps[1]), int(datecomps[2]), int(utc_starthour), int((utc_starthour-int(utc_starthour))*60), 0, 0, 0, 0))
endtime = starttime + sched_secs
delay = starttime - time.mktime(timenow)
if delay > 0:
    print >> fpcmd, "Will schedule observations starting "+str(delay)+" seconds from now."
    thetime = starttime
else:
    print >> fpcmd, "Desired start time passed "+str(-delay)+" seconds ago."
    sys.exit(1)

# overall loop is over pointings (outer loop) then freq (inner loop).
# define centre coarse chans for freq coverage
cent_freqs = [93,117,141]
#define scan lengths in sec per freq
scan_times = [120,120,120]
# define DEC ranges for meridian scans
decs = numpy.array([18.6, 1.6, -13.0, -26.7, -40.2, -55.0, -72.0])
zas  = decs + 26.7
azs  = (zas < 0)*180.0
zas = abs(zas)
gain=10.0

# calculate the total time of a full meridian stripe:
stripe_time_hrs = 0.0
stripe_sep_hrs = 0.8
for dec in decs:
  for i in range(len(scan_times)):
    stripe_time_hrs += scan_times[i]/3600.0
print >> fpcmd, "# Total time for stripe of meridian: "+str(stripe_time_hrs)+" hours"
print >> fpcmd, "# Time between stripes: "+str(stripe_sep_hrs)+" hours"
slack=stripe_sep_hrs-stripe_time_hrs
print >> fpcmd, "# Slack: "+str(slack)+" hours"
if (slack < 0):
  print >> sys.stderr, "No slack! Adjust stripe_sep_hrs."
  exit(1)

calsrcs = mwasources.CalSources(mwasources.LoadSources())
print >> fpcmd, "# There are "+str(len(calsrcs))+" cal sources"

stripe_index=0
lst_hours = lst_starthour
elapsed_time = 0
while thetime < endtime:
  # calc the LST in the middle and end of the stripe
  lst_end_hrs = lst_hours + stripe_time_hrs*side_per_sol
  if lst_end_hrs >= 24.0: lst_end_hrs -= 24.0
  lst_mid_hrs = lst_hours + stripe_time_hrs*side_per_sol/2.0
  if lst_mid_hrs >= 24.0: lst_mid_hrs -= 24.0

  lst_stripe_range = numpy.array([lst_hours,lst_end_hrs])
  print >> fpcmd, "# Stripe %d. LST range: %f %f. Mid: %f" % (stripe_index,lst_stripe_range[0],lst_stripe_range[1],lst_mid_hrs)

  # find the sources that transit during this stripe
  transit_srcs = mwasources.InRA_Range(calsrcs, lst_mid_hrs,stripe_time_hrs*side_per_sol/2.0)
  print >> fpcmd, "# There are "+str(len(transit_srcs))+" cal sources that transit in this stripe"
  calpointings = abs(decs * 0.0)
  for s in transit_srcs:
    ang_dist = abs(decs - s.dec_degs)
    min_dist = ang_dist.min()
    for i in range(len(ang_dist)):
      if ang_dist[i] == min_dist: calpointings[i] = 1
    print >> fpcmd, "# Src: "+s.name+". RA: "+str(s.ra_degs*24.0/360.0)+" hrs. DEC: "+str(s.dec_degs)+" degs."
    print >> fpcmd, "# "+str(calpointings)

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
            print >> fpcmd, str(s)+"# Distance: "+str(s.AngularSep(ra_hrs*15.0,dec)*180.0/math.pi)
    for f in range(len(cent_freqs)):
        # create time string suitable for single_observation. Here we need to
        # use localtime because the time is already in gmtime and we don't want
        # to add/subtract anything from it.
        timestr = time.strftime("%Y-%m-%d,%H:%M:%S",time.localtime(thetime))
        print >> fpcmd, "# The LST is: %g. Coarse chan: %d. UTC: %s. Elapsed time: %g" % (lst_hours,cent_freqs[f],timestr,elapsed_time)

        # create single_observation command
        cmd = "single_observation.py --starttime="+timestr+" --stoptime=++"+str(scan_times[f])+"s --freq="+str(cent_freqs[f])+",24 --obsname=MSS_DEC"+str(decs[p])+"_"+str(cent_freqs[f])+" --projectid C101 --creator=Randall --useazel --usegrid= "
        cmd += " --ra="+str(ra_hrs*15.0)+" --dec="+str(dec)
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
        thetime += scan_times[f] + gap_secs
        lst_hours += (scan_times[f] + gap_secs)*side_per_sol/3600.0
        if lst_hours >= 24.0: lst_hours -= 24.0
        elapsed_time += scan_times[f] + gap_secs

