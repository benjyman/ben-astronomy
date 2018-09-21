#code to compute the reflected emission received from the Moon at the MWA
#algorithm from Appendix B Vedantham et al 2015
#first obsid is 11273217455 '2015/09/26 16:55:28'

#from ephem import *
import math
from time import *
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance, ICRS, cartesian_to_spherical
from skyfield.api import Topos, load
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
#from pygsm import GlobalSkyModel2016

from pygsm import GSMObserver2016
from datetime import datetime

#set if data is already saved - just need to plot
plot_only = False



#moon_radius_km=1738.1
#moon_radius_deg=0.25

#Don't degrade NSIDE
#USE 512 FOR 2008 gsm AND 1024 FOR 2016 gsm
NSIDE_interim=1024
#NSIDE_final=32

r2d = 180.0/np.pi
d2r = np.pi/180.0

r2h = 12.0/np.pi
h2r = np.pi/12.0

#define frequencies
#number_coarse_chans=int(5*24+4)
#big_chan_array=np.arange(number_coarse_chans)+57
#big_freq_array=big_chan_array*1.28
#freq_array=np.arange(70,235,5)
freq_array=[150]
#freq_array_filename="freq_array_%s.npy" % date_time_string
#np.save(freq_array_filename,freq_array)

#normalise vectors to unit vectors:
def magnitude(v):
    return math.sqrt(sum(v[i]*v[i] for i in range(len(v))))

def normalise(v):
    vmag = magnitude(v)
    return [ v[i]/vmag  for i in range(len(v)) ]

#Convert HMS to degrees and vice versa
def HMS2deg(ra='', dec=''):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    D, M, S = [float(i) for i in dec.split()]
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    H, M, S = [float(i) for i in ra.split()]
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC



#Need to find where the Moon is in relation to the MWA
#set time
#ephem
##obs_date='2015/09/26 16:55:28'
##obs_epoch='2000/1/1 12'
#skyfield:


#Set the time
#This is is equivalent to 2015/09/26 16:55:28

#start: "2015_09_26_11_15_00", end:    "2015_09_26_17_15_00"
date_time_string_list = ["2015_09_26_11_15_00","2015_09_26_11_45_00","2015_09_26_12_15_00","2015_09_26_12_45_00","2015_09_26_13_15_00","2015_09_26_13_45_00","2015_09_26_14_15_00","2015_09_26_14_45_00","2015_09_26_15_15_00","2015_09_26_15_45_00","2015_09_26_16_15_00","2015_09_26_16_45_00","2015_09_26_17_15_00"]

#date_time_string_list = ["2015_09_26_16_55_28"]

if not plot_only:
   #date_time_string="2015_09_26_16_55_28"
   for date_time_string in date_time_string_list:
      print date_time_string
      year=float(date_time_string.split("_")[0])
      month=float(date_time_string.split("_")[1])
      day=float(date_time_string.split("_")[2])
      hour=float(date_time_string.split("_")[3])
      minute=float(date_time_string.split("_")[4])
      second=float(date_time_string.split("_")[5])
      
      
      ts = load.timescale()
      t=ts.utc(year, month, day, hour, minute, second)
      t_HMS=[hour, minute, second]
      
      planets = load('de421.bsp')
      earth=planets['earth']
      
      MWA = earth + Topos(latitude_degrees=-26.70331940, longitude_degrees=116.67081524 ,elevation_m=377.83)
      moon = planets['moon']
      astrometric_moon_centre_ephem=MWA.at(t).observe(moon)
      moon_centre_ra, moon_centre_dec, moon_centre_distance = astrometric_moon_centre_ephem.radec()
      
      #Define the RA and DEC of the MWA site on Earth as viewed from the Moon 
      earth_mwa_RA_deg=((moon_centre_ra.to(u.deg)/u.deg)-180.0)
      earth_mwa_dec_deg=(-1.0*(moon_centre_dec.to(u.deg)/u.deg))
      
      #So we want to imagine there is an observatory on Earth with zenith pointing at this RA DEC
      #and define a 'moon observatory with zenith facing the MWA
      moon_observatory_lat=earth_mwa_dec_deg
      
      #work out the longitude from HA = UTC - observer_lon - RA
      #for HA=0, 0=UTC-observer_lon - RA_zenith
      #observer_lon=UTC-RA_zenith
      
      t_to_ra_string="%s %s %s " % (t_HMS[0],t_HMS[1],t_HMS[2])
      t_degrees=HMS2deg(ra=t_to_ra_string, dec="0 0 0 ")
      #print t_degrees
      #print earth_mwa_RA_deg
      
      #Not sure why we need to use the time again here - hasn't the time been factored in above when computing the moon position?
      #Naybe it is just moon_observatory_lon = earth_mwa_RA_deg ?                          
      #moon_observatory_lon=float(t_degrees[0])-float(earth_mwa_RA_deg)
      moon_observatory_lon = float(earth_mwa_RA_deg)
      print "moon_observatory_lat %s" %moon_observatory_lat
      print "moon_observatory_lon %s" %moon_observatory_lon
      
      
      
      # Setup observatory location - in this case, MWA, Australia
      #latitude_degrees=-26.70331940, longitude_degrees=116.67081524 ,elevation_m=377.83
      (latitude, longitude, elevation) = (moon_observatory_lat, moon_observatory_lon, 0)
      ov = GSMObserver2016()
      ov.lon = longitude
      ov.lat = latitude
      ov.elev = elevation
      ov.date = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
      
      #initialise array of average brightness temp
      disk_averaged_temp_array=np.zeros(len(freq_array))
      disk_average_temp_filename="disk_av_temp_array_%s.npy" % date_time_string
      
      #Do for each freq
      for freq_index,freq_MHz in enumerate(freq_array):
         print str(freq_MHz)
         ##new celestial from Moon
         gsm_map_from_moon=ov.generate(freq_MHz)
      
         #plot this
         plt.clf()
         map_title="GSM Map from moon"
         hp.orthview(map=gsm_map_from_moon,coord='C',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
         fig_name="gsm_from_moon_map_%s_%sMHz.png" % (date_time_string,str(freq_MHz))
         figmap = plt.gcf()
         figmap.savefig(fig_name,dpi=100)
         plt.close()
      
         #This messes it all up! Don't do it =(
         ##want a smaller map to speed things up but still higher res than the final nside = 32
         #hp.ud_grade(gsm_map_from_moon, NSIDE_interim)
      
         moon_map=np.zeros(hp.nside2npix(NSIDE_interim))
         #zenith_pixel=hp.ang2pix(32,np.pi/2,0)
         #test_map[zenith_pixel]=1
         
         #make some vectors.
         
         #The zenith vector is important, it points from the Moon to the MWA. We need it to get the incident moon surface vector
         zenith_theta=np.pi/2
         zenith_phi=0
         zenith_vector=hp.ang2vec(zenith_theta,zenith_phi)
         
         #The vector of the moon pixel is a normal to the surface of the Moon
         print "moon map has %s pixels" % len(moon_map)
         for moon_pixel_index,moon_pixel in enumerate(moon_map):
            pixel_theta,pixel_phi=hp.pix2ang(NSIDE_interim,moon_pixel_index)
            
            #print pixel_theta
            #print pixel_phi
            
            moon_normal_vector=hp.ang2vec(pixel_theta,pixel_phi)
         
            earth_moon_vector=zenith_vector*float(moon_centre_distance.to(u.km)/u.km)
            #print earth_moon_vector
            #incident vector is:  -earth_moon_vector+moon_normal_vector
            incident_vector=-earth_moon_vector+moon_normal_vector
            #incident vector is actually the reflected vector so:
            #incident_vector=-1.0*incident_vector
         
            #print incident_vector
         
            #Need to normalise the incident vector
            incident_vector_unit=normalise(incident_vector)
            #print incident_vector_unit
            
            #dot product of these two vectors
            #dot_product = np.dot(incident_vector_unit,moon_normal_vector)
            #print dot_product
            
            #now we have the normal and incident unit vectors, we can get the reflected vector (this is all assuming a ray coming from Earth/MWA!)
            #Snells law:
            reflected_vector=2.0*np.dot(moon_normal_vector,incident_vector_unit)*moon_normal_vector - incident_vector_unit
            #print reflected_vector
            #change the sign of the reflected vector:
            reflected_vector=reflected_vector*-1.0
            #print reflected_vector
            
            #now find the pixel number  corresponding to the reflected vector
            gsm_pixel_mapped=hp.vec2pix(NSIDE_interim,reflected_vector[0],reflected_vector[1],reflected_vector[2])
            gsm_pixel_temp=gsm_map_from_moon[gsm_pixel_mapped]
            
            
            #the actual reflected temp will be scaled by the dot product between the normal and incident unit vectors
            #and the Moon albedo of 7%
            dot_product=np.dot(-moon_normal_vector,incident_vector_unit)
            gsm_pixel_temp_reflected=gsm_pixel_temp*dot_product*0.07
            #print dot_product
            #print gsm_pixel_temp_reflected
            
            #print "moon pixel index %s has value %s " % (moon_pixel_index,gsm_pixel_temp_reflected)
            percentage_done = (float(moon_pixel_index)/float(len(moon_map)))*100.
            #print "%0.9f per cent done" % percentage_done
            moon_map[moon_pixel_index]=gsm_pixel_temp_reflected
   
   
         #hp.ud_grade(moon_map, NSIDE_final)
         
         #mask the negative values
         moon_map[moon_map < 0] = np.nan
         
         #Calculate the disk-averaged temperature
         disk_averaged_temp=np.nanmean(moon_map)
         disk_averaged_temp_array[freq_index]=disk_averaged_temp
         print "disk_averaged_temp is %s K for datetime %s and freq %s MHz" % (disk_averaged_temp,date_time_string,freq_MHz)
         
         #save the moon map data
         moon_map_array_filename="moon_map_%s_%sMHz.npy" % (date_time_string,freq_MHz)
         np.save(moon_map_array_filename,moon_map)
         
         plt.clf()
         map_title="Moon Map"
         hp.orthview(map=moon_map,coord='C',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
         fig_name="moon_map_%s_%sMHz.png" % (date_time_string,str(freq_MHz))
         figmap = plt.gcf()
         figmap.savefig(fig_name,dpi=100)
         plt.close()
         
         #also print the generated gsm map for comparison
         #plt.clf()
         #map_title="gsm map"
         #hp.orthview(map=gsm_map_from_moon,coord='C',half_sky=False,xsize=400,title=map_title,rot=(0,0,0))
         #fig_name="gsm_map_%s_%sMHz.png" % (date_time_string,str(freq_MHz))
         #figmap = plt.gcf()
         #figmap.savefig(fig_name,dpi=100)
         
      #save the disk averaged temp array
      np.save(disk_average_temp_filename,disk_averaged_temp_array)

else:
   #(plot_only)
   freq_array = [150]
   for date_time_string in date_time_string_list:
      year=float(date_time_string.split("_")[0])
      month=float(date_time_string.split("_")[1])
      day=float(date_time_string.split("_")[2])
      hour=float(date_time_string.split("_")[3])
      minute=float(date_time_string.split("_")[4])
      second=float(date_time_string.split("_")[5])
      
      for freq_index,freq_MHz in enumerate(freq_array):
         print str(freq_MHz)
         #moon_map_array_filename="/data/moon/2017/GSM_2016_reflection/moon_map_%s_%sMHz.npy" % (date_time_string,freq_MHz)
         moon_map_array_filename="moon_map_%s_%sMHz.npy" % (date_time_string,freq_MHz)
         moon_map=np.load(moon_map_array_filename)
         #mask the negative values
         moon_map[moon_map < 0] = np.nan
         
         plt.clf()
         map_title="Moon Reflection UTC %.0f:%.0f:%.0f" % (hour,minute,second)
         hp.orthview(map=moon_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
         fig_name="paper_moon_map_%s_%sMHz.png" % (date_time_string,str(freq_MHz))
         figmap = plt.gcf()
         figmap.savefig(fig_name,dpi=500)


