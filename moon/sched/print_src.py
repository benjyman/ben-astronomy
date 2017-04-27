#!/usr/bin/python
#
# D. Mitchell: Print out the location of various sources.
#
# Requires:
#   - PyEphem (install "setuptools" and then type "sudo easy_install pyephem")
#

import sys, getopt, string, os
from ephem import *
from time import *

#----------------------------------------------------------------------------------------------------------------------#
# check for command-line options

def help_message():
    print '''
Options: --date=str  -- Optional. Date string (UTC), e.g., --date='2008/2/28 16:22:56'.
                                  Default is the current time.
         --epoch=str -- Optional. Epoch string, e.g., --epoch='2008/2/28 16:22:56'.
                                  Default is J2000 (--epoch='2000/1/1 12').
         --radec=str -- Optional. Add an additional location to track (in J2000),
                                  e.g., --radec='05:19:49.7,-45:46:44'.
         --short     -- Optional. Only print RA & Declination for each source.
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

add_radec = '';

print_type = 'l';

#----------------------------------------------------------------------------------------------------------------------#
# read command-line arguments

try:
    options, xarguments = getopt.getopt(sys.argv[1:], '', ['date=','epoch=','radec=','short'])
except getopt.error:
    help_message()

for a in options[:]:
    if a[0] == '--date' and a[1] != '':
        obs_date = a[1]
        options.remove(a)
        break
    elif a[0] == '--date' and a[1] == '':
        print '--date expects an argument'
        sys.exit(0)
for a in options[:]:
    if a[0] == '--epoch' and a[1] != '':
        obs_epoch = a[1]
        options.remove(a)
        break
for a in options[:]:
    if a[0] == '--radec' and a[1] != '':
        add_radec = "User-defined,f|J,%s,0.0" % a[1]
        options.remove(a)
        break
for a in options[:]:
    if a[0] == '--short':
        print_type = 's'
        options.remove(a)
        break
    elif a[0] == '--epoch' and a[1] == '':
        print '--epoch expects an argument'
        sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------#
# MAIN

# Create an observer
MRO = Observer()

# Set the observer at Boolardy
MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 377.83

# Set the time.
if obs_date  != '': MRO.date  = obs_date
if obs_epoch != '': MRO.epoch = obs_epoch

sname = [];
slist = [];
SJy   = [];

sname.append("Sun");              SJy.append(" 1e4 Jy - 1e8 Jy "); slist.append("");
sname.append("Moon");              SJy.append(" 60 Jy in FM band! "); slist.append("");

sname.append("SCP");              SJy.append("-----------------"); slist.append("");

sname.append("Fornax A");         SJy.append(" 259 Jy @ 408 MHz"); slist.append(readdb("ForA,f|J,03:22:41.7,-37:12:30,0.0"));
sname.append("EoR field 1");      SJy.append("-----------------"); slist.append(readdb("EoR,f|J,04:00:00.0,-30:00:00,0.0"));
sname.append("J0444-2809");       SJy.append("  65 Jy @  80 MHz"); slist.append(readdb("EoR,f|J,04:44:36.5,-28:09:44,0.0"));
sname.append("Pictor A");         SJy.append(" 410 Jy @  80 MHz"); slist.append(readdb("PicA,f|J,05:19:49.7,-45:46:44,0.0"));
sname.append("Crab SNR (Tau A)"); SJy.append("1530 Jy @ 178 MHz"); slist.append(readdb("Crab,f|J,05:34:31.97,+22:00:52.1,0.0"));
sname.append("3C 161");           SJy.append("  96 Jy @  80 MHz"); slist.append(readdb("3C161,f|J,06:27:10,-05:53:10,0.0"));
sname.append("PSR B0628-28");     SJy.append("  ~3 Jy @  80 MHz"); slist.append(readdb("J0630-2834,f|J,06:30:49.4,-28:34:43,0.0")); 
sname.append("Puppis A");         SJy.append(" 198 Jy @ 408 MHz"); slist.append(readdb("PupA,f|J,08:24:10,-43:00:00,0.0"));
sname.append("PSR B0834+06");     SJy.append("  ~1 Jy @  74 MHz"); slist.append(readdb("J0837+0610,f|J,08:37:05.9,+06:10:11,0.0")); 
sname.append("Hydra A");          SJy.append(" 441 Jy @  80 MHz"); slist.append(readdb("HydA,f|J,09:18:05.7,-12:05:44,0.0"));
sname.append("PSR B0950+08");     SJy.append("  ~8 Jy @  80 MHz"); slist.append(readdb("J0953+0755,f|J,09:53:09.3,+07:55:37,0.0")); 
sname.append("EoR field 2");      SJy.append("-----------------"); slist.append(readdb("EoR,f|J,10:20:00.0,-10:00:00,0.0"));
sname.append("Virgo A");          SJy.append("1414 Jy @  80 MHz"); slist.append(readdb("VirA,f|J,12:30:49.4,+12:23:28,0.0"));
sname.append("Centaurus A");      SJy.append("2746 Jy @ 408 MHz"); slist.append(readdb("CenA,f|J,13:25:27.6,-43:01:09,0.0"));
sname.append("Hercules A");       SJy.append(" 751 Jy @  80 MHz"); slist.append(readdb("HerA,f|J,16:51:08.1,+04:59:33,0.0"));
sname.append("3C 353");           SJy.append(" 349 Jy @  80 MHz"); slist.append(readdb("3C353,f|J,17:20:28.1,-00:58:47,0.0"));
sname.append("Sgr A*");           SJy.append("-----------------"); slist.append(readdb("SgrA,f|J,17:45:40.0,-29:00:28,0.0"));
sname.append("J2358-6054");       SJy.append("  57 Jy @ 408 MHz"); slist.append(readdb("PKS2152,f|J,21:57:05.8,-69:41:22,0.0"));
sname.append("J2157-6941");       SJy.append("  80 Jy @ 408 MHz"); slist.append(readdb("PKS2356,f|J,23:58:59.2,-60:54:57,0.0"));

if add_radec != '':
    sname.append("User-defined");
    SJy.append("-----------------");
    slist.append(readdb(add_radec));

print "\nUTC: %s\n" % ( MRO.date )
print "At %s N, %s E (LST %s = %f):\n" % ( MRO.lat, MRO.long, MRO.sidereal_time(), MRO.sidereal_time()*r2h )

for snum in range(0,len(slist)):

    src = slist[snum]

    if sname[snum] == "Sun":
        sun = Sun(MRO);
        az = sun.az  * r2d
        el = sun.alt * r2d
        ra0, dec0 = MRO.radec_of(sun.az, sun.alt)
    elif sname[snum] == "Moon":
        moon = Moon(MRO);
        az = moon.az  * r2d
        el = moon.alt * r2d
        ra0, dec0 = MRO.radec_of(moon.az, moon.alt)
    elif sname[snum] == "SCP":
        scp = Equatorial('0', '-90', epoch=J2000)
        az = 180
        el = -MRO.lat * r2d
        ra0, dec0 = scp.ra, scp.dec
    else:
        src.compute(MRO);
        ra0, dec0 = MRO.radec_of(src.az, src.alt)
        az = src.az  * r2d
        el = src.alt * r2d

    if print_type == 'l':
        print "%17s: ra = %11s, dec = %11s, az = %8.4f, el = %8.4f: %s" % (sname[snum], ra0, dec0, az, el, SJy[snum])

    if print_type == 's':
        print "%17s: ra = %11s, dec = %11s" % (sname[snum], ra0, dec0)

print ""
