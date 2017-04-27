#!/usr/bin/python
# script to generate the LSTs and sunset times for Randall's google spreadsheet shedules.
import ephem
import math

dfile=open('/tmp/dates.txt','w')
lfile=open('/tmp/lsts.txt','w')
sfile=open('/tmp/sunsets.txt','w')
rfile=open('/tmp/sunrises.txt','w')

site = ephem.Observer()
site.long = '116:40:14.93'
site.lat = '-26:42:11.95'
site.date= '2015/01/05'
sun = ephem.Sun()
for d in range(365/2):
  settime=site.next_setting(sun)
  risetime=site.next_rising(sun)
  lst = site.sidereal_time() # returns an Angle object
  print >> dfile, site.date
  print >> lfile, str(lst*12/math.pi)
  print >> sfile, ephem.Date(settime + 8./24)
  print >> rfile, ephem.Date(risetime + 8./24)
  print site.date, str(lst*12/math.pi), ephem.Date(settime + 8./24),ephem.Date(risetime + 8./24)
  site.date += 1
