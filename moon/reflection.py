#code to compute the reflected emission received from the Moon at the MWA
#algorithm from Appendix B Vedantham et al 2015
#first obsid is 11273217455 '2015/09/26 16:55:28'

from ephem import *
from time import *
from astropy import units as u
from astropy.coordinates import SkyCoord


r2d = 180.0/pi
d2r = pi/180.0

r2h = 12.0/pi
h2r = pi/12.0


#time
obs_date='2015/09/26 16:55:28'
obs_epoch='2000/1/1 12'


# Create an observer
MRO = Observer()

# Set the observer at the MWA
MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 377.83

# Set the time.
if obs_date  != '': MRO.date  = obs_date
if obs_epoch != '': MRO.epoch = obs_epoch


#RA Dec of the Moon
moon = Moon(MRO);
az = moon.az  * r2d
el = moon.alt * r2d
ra0, dec0 = MRO.radec_of(moon.az, moon.alt)
ra_deg=r2d*(hours(ra0))
dec_deg=r2d*(degrees(dec0))

print "Moon is at RA %s degrees Dec %s " % (ra_deg,dec_deg)
print "Alt %s, Az %s" % (el,az)

#define HEALPIX (N=32) grid for Lunar surface

#extract the pixels visible to the MWA

#For each pixel:
#calculate position vector i 
pixel_position_ra_dec=SkyCoord(ra_deg*u.degree, dec_deg*u.degree,distance=1*u.m)
print pixel_position_ra_dec
pixel_position_cartesian_x=pixel_position_ra_dec.cartesian.x
print pixel_position_cartesian_x
pixel_position_cartesian_y=pixel_position_ra_dec.cartesian.y
print pixel_position_cartesian_y
pixel_position_cartesian_z=pixel_position_ra_dec.cartesian.z
print pixel_position_cartesian_z

#make the vector

#do the same of rthe MWA position with EarthLocation


#calculate normal vector n

#define the plane of incidence and plane of reflection 

