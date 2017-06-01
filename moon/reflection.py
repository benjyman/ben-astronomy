#code to compute the reflected emission received from the Moon at the MWA
#algorithm from Appendix B Vedantham et al 2015
#first obsid is 11273217455 '2015/09/26 16:55:28'

#from ephem import *
from time import *
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from skyfield.api import Topos, load
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt  

moon_radius_km=1738.1
moon_radius_deg=0.25
NSIDE=32

r2d = 180.0/np.pi
d2r = np.pi/180.0

r2h = 12.0/np.pi
h2r = np.pi/12.0

#Functions for converting healpix pixel index to RA, DEC and vice versa
def IndexToRaDec(index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(np.pi*2.-phi),-np.degrees(theta-np.pi/2.)

def RaDecToIndex(RA,Dec):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-Dec+90.),np.radians(360.-RA))


#set time
#ephem
##obs_date='2015/09/26 16:55:28'
##obs_epoch='2000/1/1 12'
#skyfield:
ts = load.timescale()
t=ts.utc(2015, 9, 26, 16, 55, 28.0)


##With Pyephem:
## Create an observer
##MRO = Observer()

## Set the observer at the MWA
##MRO.lat, MRO.long, MRO.elevation = '-26:42:11.95', '116:40:14.93', 377.83

## Set the time.
##if obs_date  != '': MRO.date  = obs_date
##if obs_epoch != '': MRO.epoch = obs_epoch

planets = load('de421.bsp')
earth=planets['earth']

MWA = earth + Topos(latitude_degrees=-26.70331940, longitude_degrees=116.67081524 ,elevation_m=377.83)
moon = planets['moon']
astrometric_moon_centre_ephem=MWA.at(t).observe(moon)
moon_centre_ra, moon_centre_dec, moon_centre_distance = astrometric_moon_centre_ephem.radec()

#define a healpix grid (N=32) on the lunar surface
healpix_grid_full = np.zeros(hp.nside2npix(NSIDE))
hp.mollview(healpix_grid_full, title="Mollview image RING")
fig_name="mollview_image_blank"
figmap = plt.gcf()
figmap.savefig(fig_name,dpi=100)
plt.clf()

#now extract all the moon pixels
for index,pix in enumerate(healpix_grid_full):
  #define moon pixels as a disk with radius 15 arcmin (0.25 deg) and centre moon_centre_ra, moon_centre_dec
  pixel_RA,pixel_Dec=IndexToRaDec(index)
  #print pixel_RA
  #print pixel_Dec

  min_moon_ra=(moon_centre_ra.to(u.deg)/u.deg)-moon_radius_deg
  #print min_moon_ra
  max_moon_ra=(moon_centre_ra.to(u.deg)/u.deg)+moon_radius_deg
  #print max_moon_ra
  min_moon_dec=(moon_centre_dec.to(u.deg)/u.deg)-moon_radius_deg
  #print min_moon_dec
  max_moon_dec=(moon_centre_dec.to(u.deg)/u.deg)+moon_radius_deg
  #print max_moon_dec
  #if (1.8 < hp.pix2ang(NSIDE,index)[0] < 2.2 and 1.8 < hp.pix2ang(NSIDE,index)[1] < 2.2):
  if (min_moon_ra < pixel_RA < max_moon_ra and min_moon_dec < pixel_Dec < max_moon_dec):
     print min_moon_ra
     healpix_grid_full[index]=1

print healpix_grid_full

hp.mollview(healpix_grid_full, title="Mollview Moon Pix")

fig_name="mollview_moon_pix"
figmap = plt.gcf()
figmap.savefig(fig_name,dpi=100)
plt.clf()

#print(moon_centre_ra.hstr())
#print(moon_centre_dec.dstr())
#print(moon_centre_distance.km)

dec_offset=0.0*u.degree
ra_offset=0.0*u.degree

moon_pix_dec = moon_centre_dec.to(u.deg)+dec_offset
moon_pix_ra=moon_centre_ra.to(u.deg)+ra_offset
moon_pix_distance_km=moon_centre_distance.km - moon_radius_km



moon_centre_skycoord=SkyCoord(ra=moon_centre_ra.to(u.deg), dec=moon_centre_dec.to(u.deg), distance=moon_centre_distance.km*u.km)
#print float(moon_centre_skycoord.cartesian.x/u.km)

moon_centre_vector=np.array([moon_centre_skycoord.cartesian.x/u.km,moon_centre_skycoord.cartesian.y/u.km,moon_centre_skycoord.cartesian.z/u.km])
#print moon_centre_vector

moon_pixel_skycoord=SkyCoord(ra=moon_pix_ra, dec=moon_pix_dec, distance=moon_pix_distance_km*u.km)
moon_pixel_vector=np.array([moon_pixel_skycoord.cartesian.x/u.km,moon_pixel_skycoord.cartesian.y/u.km,moon_pixel_skycoord.cartesian.z/u.km])
#print moon_pixel_vector

normal_vector=moon_centre_vector-moon_pixel_vector
#print normal_vector

#calculate tangential vector
tangent=np.cross(np.cross(normal_vector,moon_centre_vector),normal_vector)
#print tangent











##RA Dec of the Moon centre
##moon_centre = Moon(MRO);
##moon_centre_az = moon_centre.az  * r2d
##moon_centre_el = moon_centre.alt * r2d
##moon_centre_ra0, moon_centre_dec0 = MRO.radec_of(moon_centre.az, moon_centre.alt)
##moon_centre_ra_deg=r2d*(hours(moon_centre_ra0))
##moon_centre_dec_deg=r2d*(degrees(moon_centre_dec0))

##print "Moon centre is at RA %s degrees Dec %s " % (moon_centre_ra_deg,moon_centre_dec_deg)
##print "Alt %s, Az %s" % (moon_centre_el,moon_centre_az)

#define HEALPIX (N=32) grid for Lunar surface

#extract the pixels visible to the MWA



#For Moon centre

##calculate position vector i 
##moon_centre_pixel_position_ra_dec=SkyCoord(moon_centre_ra_deg*u.degree, moon_centre_dec_deg*u.degree,distance=1*u.m)
##print moon_centre_pixel_position_ra_dec
##moon_centre_pixel_position_cartesian_x=moon_centre_pixel_position_ra_dec.cartesian.x
##print moon_centre_pixel_position_cartesian_x
##moon_centre_pixel_position_cartesian_y=moon_centre_pixel_position_ra_dec.cartesian.y
##print moon_centre_pixel_position_cartesian_y
##moon_centre_pixel_position_cartesian_z=moon_centre_pixel_position_ra_dec.cartesian.z
##print moon_centre_pixel_position_cartesian_z

##define the moon centre position vector
##moon_centre_pixel_position_vector=[moon_centre_pixel_position_cartesian_x,moon_centre_pixel_position_cartesian_y,moon_centre_pixel_position_cartesian_z]


#then for each Moon pixel

#Get RA Dec from healpix function
#define position vector as above (for now use centre)
##pixel_position_vector=moon_centre_pixel_position_vector

#calculate vector normal to Moon surface 
#normal_vector=

#do the same of rthe MWA position with EarthLocation.from_geodetic(lon, lat, height=0.0, ellipsoid=None)


#calculate normal vector n

#define the plane of incidence and plane of reflection 

