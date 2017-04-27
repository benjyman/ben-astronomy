#!/usr/bin/python

# A scripty to track a list of sources over a nominated set of freqs (can be
# different for each source) to start at a specified time and run for a 
# specified time
# Randall Wayth. April 2011.

import math,time,sys,os,logging
import ephem
import mwasources

logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger(__name__)   # default logger level is WARNING

    
def e2h(ha,dec,lat):
    """Convert equatorial ha/dec coords (in hours, degrees) to az,el (degs) give an observer latitute (degs). Returns (az,el) in degs"""
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


def trackSources(srclist,lst_starthour,starttime,total_time_sec,gap_secs,projid='C001',fpcmd=sys.stdout,freqs_override=None,dwell_override_secs=None,name_override=None,freqres=40,timeres=0.5):
  """track a list of mwasources. starttime is a time object of the built-in type like from time.mktime"""
  srcindex = 0
  elapsed_time = 0.0
  n_done =0
  thetime = starttime
  endtime = starttime + total_time_sec
  if lst_starthour > 24.0: lst_starthour -= 24.0
  lst_hours = lst_starthour
  site_lat_degs = -26.7033

  while thetime < endtime:
    curr_src = srclist[srcindex]
    if dwell_override_secs is None:
        dwelltime_sec = curr_src.dwelltime_sec
    else:
        dwelltime_sec = dwell_override_secs 
    if freqs_override is None:
        freqs=curr_src.freqs
    else:
        freqs=freqs_override
    for freq in freqs:
        ha_hrs = lst_hours - curr_src.ra_degs/15.0
        (az,el) = e2h(ha_hrs,curr_src.dec_degs,site_lat_degs)
#        print >> fpcmd,"# Src: "+curr_src.name+". LST: "+str(lst_hours)
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
        cmd = "single_observation.py --starttime="+timestr
        cmd += " --stoptime=++"+str(int(dwelltime_sec))+"s"
        if freq.find(';') != -1:
            cmd += " --freq='"+str(freq)+"'"
            cmd += " --obsname="+curr_src.name
        else:
            cmd += " --freq="+str(freq)+",24"
            cmd += " --obsname="+curr_src.name+"_"+str(freq)
        cmd += " --inttime="+str(timeres)
        cmd += " --freqres="+str(freqres)
        if name_override is not None:
            cmd += " --creator="+name_override
        cmd += " --useazel --usegrid= --project="+projid
        if curr_src.gain is not None:
            cmd += " --gain_control_value="+str(curr_src.gain)
        if curr_src.known == False:
            cmd += " --ra="+str(curr_src.ra_degs)+" --dec="+str(curr_src.dec_degs)
        else:
            cmd += " --source="+curr_src.name
        if curr_src.iscal:
            cmd += " --calibration --calibrator="+curr_src.name

        print >> fpcmd, cmd


        # increment time
        thetime += dwelltime_sec + gap_secs
        lst_hours += (dwelltime_sec + gap_secs)*1.00278551532/3600.0
        if lst_hours > 24.0: lst_hours -= 24.0
        elapsed_time += dwelltime_sec + gap_secs
    # go to next source
    srcindex = (srcindex+1) % len(srclist)
    if srcindex == 0:
        if n_done > 0:
            # new cycle
            n_done = 0
        else:
            print >> fpcmd, "# Hmm. Seem to have skipped ALL the sources. Nothing interesting at this LST... Advancing 1 minute"
            thetime += 60.0
            lst_hours += 1.00278551532/60.0
            elapsed_time += 60.0

if __name__ == "__main__":
  # simple test of this code
  src = mwasources.Source('0924-2201',141.0,-22.0,freqs=['175'],dwelltime_sec=296)
  hyda= mwasources.Source('HydA',139.5236198,-12.0955426,iscal=True,known=True,freqs=['175'],dwelltime_sec=56)
  pupa= mwasources.Source('PupA',126.0307,-42.9963,known=True)
  sgra= mwasources.Source('SgrA',266.417,-29.0078,gain=10,known=True)
  starttime=time.mktime((2013, 05, 02, 20,0, 0, 0, 0, 0))
  trackSources([pupa],8.0,starttime,1800.0/2,4,'C001',sys.stdout)
  sys.exit(0)
