#!/usr/bin/python
import math, time, sys, os, getopt
import ephem

def usage():
    print >> sys.stderr,"Usage: %s <date>" % (sys.argv[0])
    print >> sys.stderr,"\tdate: YYYY/MM/DD"
    sys.exit(1)

site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
if len(sys.argv) > 1:
    site.date=sys.argv[1]
sun = ephem.Sun()
sun.compute(site)
print str(sun.ra*12/math.pi) + " hours"
