#!/usr/bin/python
import ephem
import math

site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
lst = site.sidereal_time() # returns an Angle object
now = ephem.now()   # returns a Date object
sun = ephem.Sun()
sun.compute(now)
settime=site.next_setting(sun)  # returns a Date object
risetime=site.next_rising(sun)  # returns a Date object
print "LST at MRO: "+str(lst) + " ("+str(lst/math.pi*180.0/15.0)+")"
print "Sunset at UTC: "+str(settime)
print "Sunset WA time: "+str(ephem.date(settime+8.0*ephem.hour))
print "Sunrise at UTC: "+str(risetime)
print "Sunrise WA time: "+str(ephem.date(risetime+8.0*ephem.hour))
print "RA of Sun: "+str(sun.ra)
time_to_set = (settime-now)/ephem.hour  # convert time diff to hours
lst_settime=lst/math.pi*180.0/15.0+time_to_set/0.99726958
while lst_settime > 24.0: lst_settime -= 24.0
print "LST at sunset: "+str(lst_settime)
