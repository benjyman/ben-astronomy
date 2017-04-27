"""Utilities for creating an MWA schedule with various observation modes
and correlator modes. This is meant to be used as a library, with top-level
programs making use of this code.
Randall Wayth. 2012-2014
"""
import math, time, sys, os, getopt
import ephem
import numpy
import mwasources
import track_srcs

class ObsMode:
    """Object to capture MWA observing mode (freq, duration etc) and
    correlator mode info
    """
    def __init__(self,dwell=120,freq="121,24"):
        self.dwell=dwell
        self.freq=freq
        
def doCalScans(calsrcs,lst_hours,cal_dwell_sec,cal_gap_secs):
    """include a set of calibration scans choosing the best cal for the given LST"""
    print >> fpcmd, "# *** Dedicated calibration scans for LST: %.2f hrs" % (lst_hours)
    cals = mwasources.FindCal(calsrcs,lst_hours*15.0)   # this returns a list of tuples where each tuple is (mwasource,za)
    if len(cals) == 0:
        print >> sys.stderr, "WARNING: Unable to find calibrator for start of scan"
    else:
        # send the first item in the cals list, without the zenith angle arg
        track_srcs.trackSources([x[0] for x in cals], lst_hours, thetime, 1, cal_gap_secs, proj_id, fpcmd=fpcmd, freqs_override=cent_freqs, name_override=creator, dwell_override_secs=cal_dwell,timeres=timeres,freqres=freqres)


