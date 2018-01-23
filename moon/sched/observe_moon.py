#!/opt/local/bin/python
###!/usr/bin/python
#
#B. McKinley
#Search through a range of times and find where (as seen from the MWA):
#1. The Moon is above a certain threshold elevation, and
#2. The Moon is below a certain threshold galactic latitude, and
#3. The sun is down
#4. Other bright sources not too close.
#
# Used print_src by as a starting point (D. Mitchell)
#
# Requires:
#   - PyEphem (install "setuptools" and then type "sudo easy_install pyephem")
#

import sys, getopt, string, os
from ephem import *
from time import *
from dateutil import rrule
from datetime import datetime, timedelta
from astropy import units as u
from astropy.coordinates import SkyCoord

#function
def time_string_to_decimals(time_string):
   fields = time_string.split(":")
   hours = fields[0] if len(fields) > 0 else 0.0
   minutes = fields[1] if len(fields) > 1 else 0.0
   seconds = fields[2] if len(fields) > 2 else 0.0
   return float(hours) + (float(minutes) / 60.0) + (float(seconds) / pow(60.0, 2)) 


#----------------------------------------------------------------------------------------------------------------------#
# check for command-line options

def help_message():
    print '''
Options: --startdate=str  -- Optional. Date string (UTC) for start of search, e.g., --startdate='2008/2/28 16:22:56'.
                                  Default is the current time.
         --duration=str  -- Optional. Duration of search in days, default is 1, e.g. --duration='100'
         --frequency=str    ---Optional. Frequency that the program computes the Moon position within the time range specified. Default is hourly, e.g. --frequency='minutely'
         --moon_elev=str -- Optional. Minimum elevation of the Moon in degrees. Default is 60, e.g. --moon_elev='70'
         --gal_sep=str  --Optional. Minimum angular distance of Moon from the Galactic plane in degrees. Default is 50, e.g. --gal_sep='60'
         --sun_elev=str  --Optional. Maximum elevation of the Sun in degrees. Defualt is -5, e.g. --sun_elev='-10'  
         --epoch=str -- Optional. Epoch string, e.g., --epoch='2008/2/28 16:22:56'.,
                                  Default is J2000 (--epoch='2000/1/1 12').,
         --radec=str -- Optional. Add an additional location to track (in J2000),
                                  e.g., --radec='05:19:49.7,-45:46:44'.,
         --schedule=str  -- Optional. Number of hours to observe at transit. Prints the transit RA and the start LST in decimal hours (sets duration to 1), Default is 6 hours e.g. --schedule='6'
         --mode=str --Optional. Correlator mode e.g. --mode='40_1' or --mode='10_0.5', default is '40_1'
    '''
    sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------#

# initialise command-line variables.

r2d = 180.0/pi
d2r = pi/180.0

r2h = 12.0/pi
h2r = pi/12.0

obs_date  = '';
obs_epoch = '';
duration = 1
frequency='hourly'
moon_elev=60
gal_sep=50
sun_elev=-5
sched_hours=6
mode='40_1'

add_radec='';
print_type='l';

#----------------------------------------------------------------------------------------------------------------------#
# read command-line arguments

try:
    options, xarguments = getopt.getopt(sys.argv[1:], '', ['startdate=','duration=','frequency=','moon_elev=','gal_sep=','sun_elev=','epoch=','radec=','schedule=','mode='])
except getopt.error:
    help_message()

for a in options[:]:
    if a[0] == '--startdate' and a[1] != '':
        obs_date = a[1]
        options.remove(a)
        break
    elif a[0] == '--startdate' and a[1] == '':
        print '--startdate expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--duration' and a[1] != '':
        duration = float(a[1])
        options.remove(a)
        break
    elif a[0] == '--duration' and a[1] == '':
        print '--duration expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--frequency' and a[1] != '':
        frequency = a[1]
        options.remove(a)
        break
    elif a[0] == '--frequency' and a[1] == '':
        print '--frequency expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--moon_elev' and a[1] != '':
        moon_elev = float(a[1])
        options.remove(a)
        break
    elif a[0] == '--moon_elev' and a[1] == '':
        print '--moon_elev expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--gal_sep' and a[1] != '':
        gal_sep = float(a[1])
        options.remove(a)
        break
    elif a[0] == '--gal_sep' and a[1] == '':
        print '--gal_sep expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--sun_elev' and a[1] != '':
        sun_elev = float(a[1])
        options.remove(a)
        break
    elif a[0] == '--sun_elev' and a[1] == '':
        print '--sun_elev expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--epoch' and a[1] != '':
        obs_epoch = a[1]
        options.remove(a)
        break
    elif a[0] == '--epoch' and a[1] == '':
        print '--epoch expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--radec' and a[1] != '':
        add_radec = "User-defined,f|J,%s,0.0" % a[1]
        options.remove(a)
        break
    elif a[0] == '--radec' and a[1] == '':
        print '--radec expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--schedule' and a[1] != '':
        print_type = 't'
        sched_hours=float(a[1])
        options.remove(a)
        break
    elif a[0] == '--schedule' and a[1] == '':
        print '--schedule expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--mode' and a[1] != '':
        print_type = 't'
        mode=a[1]
        options.remove(a)
        break
    elif a[0] == '--mode' and a[1] == '':
        print '--mode expects an argument'
        sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------#

# MAIN

#setup file to write for observations

#outfile="G0017_2015_schedule_obs_test.txt"
outfile="G0017_2018A_schedule_obs.txt"
fileprint = open(outfile, 'w')
#print >> fileprint, "# set of commands to schedule G0017 (moon EoR global signal) for 2015-A"
print >> fileprint, "# set of commands to schedule G0017 (moon EoR global signal)"

# Create an observer
MRO = Observer()

# Set the observer at Boolardy
MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 377.83

if obs_date  != '': 
   start=datetime.strptime(obs_date, '%Y/%m/%d %H:%M:%S')
else:
   start=datetime.now()

end=start+timedelta(days=duration)

if frequency=='hourly':
   frequency=rrule.HOURLY
if frequency=='minutely':
   frequency=rrule.MINUTELY

good_night=False

#to check whether we already have observations on this date
previous_date=''

for dt in rrule.rrule(frequency,dtstart=start,until=end):
   #work out if this night is a good night to observe
   current_date=str(dt)[0:10]
   #calculate the Moon's position for each time
   #set the time at the MRO (UT)
   MRO.date=dt
   moon = Moon(MRO);
   az = moon.az  * r2d
   el = moon.alt * r2d
   ra0, dec0 = MRO.radec_of(moon.az, moon.alt)
   #galactic coords:
   gal=Galactic(moon)
   long0 = gal.long
   lat0 = gal.lat * r2d
   LST_time=MRO.sidereal_time() 
   LST_time_r2h=MRO.sidereal_time()*r2h

   #check sun position
   sun = Sun(MRO);
   sun_el = sun.alt * r2d
   
   #work out if its a good night to observe 
   if sun_el < sun_elev and el > moon_elev and (lat0 < -1.0*gal_sep or lat0 > gal_sep):
      good_night=True
   else:
      good_night=False
   #if its a good night and we haven't already written this night to the file and we want to then do it:   
   if (good_night==True) and (current_date != previous_date) and print_type=='t':
      #work out Moon transit
      MRO.date=dt
      moon_transit_date=str(MRO.next_transit(moon)).split()[0]
      print moon_transit_date
      print current_date
      #set MRO date to transit date
      MRO.date=moon_transit_date   
      moon_transit = Moon(MRO);
      transit_az = moon_transit.az  * r2d
      transit_el = moon_transit.alt * r2d
      ra_t, dec_t = MRO.radec_of(moon_transit.az, moon_transit.alt)
      transit_radec=SkyCoord(ra=str(ra_t),dec=str(dec_t), frame='icrs',unit=(u.hour,u.deg))
      dec_t_decimal=transit_radec.dec.deg 
      ra_t_string=str(dt)[0:10] + ' ' +  str(ra_t).split(':')[0]+':'+ str(ra_t).split(':')[1]+':'+str(ra_t).split(':')[2].split('.')[0]
      ra_t_time=datetime.strptime(ra_t_string, '%Y-%m-%d %H:%M:%S')
      start_LST=ra_t_time - timedelta(hours=sched_hours/2)
      start_LST_string=str(start_LST)
      start_LST_short=start_LST_string[11:20]
      start_LST_decimal=time_string_to_decimals(start_LST_short)
      ra_t_decimal=float(time_string_to_decimals(str(ra_t)))
    
      #if --schedule is set and we haven't already written this night to the file
      #print "UTC Date: %1s, start_LST: %8.3f, transit_ra = %8.3f, transit_dec = %8.3f, el = %8.3f, obs length = %8.1f hrs " % (dt, start_LST_decimal, ra_t_decimal, dec_t_decimal, el, sched_hours)
      #print >> fileprint, "UTC Date: %1s, start_LST: %8.3f, transit_ra = %8.3f, transit_dec = %8.3f, el = %8.3f, obs length = %8.1f hrs " % (dt, start_LST_decimal, ra_t_decimal, dec_t_decimal, el, sched_hours)
      print >> fileprint, '# On Moon for transit at %1s' %  (moon_transit_date) 
      print >> fileprint, 'radec="%1.2f,%1.2f"' % (ra_t_decimal, dec_t_decimal)
<<<<<<< Updated upstream
      print >> fileprint, './track_a_source.py -s -e -d 120 -m 40_1 -C $radec -u BenMcKinley %1s %2.1f %2.1f EoRMoon G0017' % (str(moon_transit_date)[0:10], sched_hours,start_LST_decimal)
=======
      print >> fileprint, './track_a_source.py -s -e -d 120 -m %s -C $radec -u BenMcKinley %1s %2.1f %2.1f EoRMoon G0017' % (mode,str(moon_transit_date)[0:10], sched_hours,start_LST_decimal)
>>>>>>> Stashed changes
      print >> fileprint, "#Off Moon\n"
      previous_date = current_date



   #only print all this stuff if --schedule not set  if the criteria above are met:
   if good_night and print_type=='l':
      print "UTC Date: %1s, LST: %s=%f, %1s: ra = %11s, dec = %11s, gal_lat = %11s, el = %8.4f, %1s: el = %8.4f " % (dt,LST_time, LST_time_r2h,'Moon', ra0, dec0, lat0, el, 'Sun',sun_el)

print ""

