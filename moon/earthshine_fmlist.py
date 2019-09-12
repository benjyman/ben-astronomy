#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_moon
import pandas as pd
import math
import sys
from mpl_toolkits.basemap import Basemap
from itertools import chain
from datetime import datetime, timedelta, date, time
import matplotlib.dates as mdates

def draw_map(m, scale=0.2):
    # draw a shaded-relief image
    m.shadedrelief(scale=scale)
    
    # lats and longs are returned as a dictionary
    lats = m.drawparallels(np.linspace(-90, 90, 13))
    lons = m.drawmeridians(np.linspace(-180, 180, 13))

    # keys contain the plt.Line2D instances
    lat_lines = chain(*(tup[1][0] for tup in lats.items()))
    lon_lines = chain(*(tup[1][0] for tup in lons.items()))
    all_lines = chain(lat_lines, lon_lines)
    
    # cycle through these lines and set the desired style
    for line in all_lines:
        line.set(linestyle='-', alpha=0.3, color='w')
        
plot_only = True
fmlist_data_filename = '/md0/moon/lauren/moonraker.csv'
#These are numpy arrays. When you open them with python numpy you will see they have a shape of (124,50,3). This corresponds to (n_channels, n_observations, 3), the 3 is because there is specular, diffuse, and observation ID.
mwa_data_filename = '/md0/moon/lauren/big_RFI_vs_time_freq_2015B_05.npy'

time_delta=np.linspace(0, 24, 96)*u.hour
start_time_utc = Time('2015-09-26 00:00:00')
utc_times=start_time_utc+time_delta
MWA=EarthLocation(lat=-26.7*u.deg, lon=116.67*u.deg)

powers_sum_timestep_array_janskys_filename = "powers_sum_timestep_array_janskys.npy"

moon_radius = 1737.4 * 1000.0 #metres, mean radius from NASA website https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
pi=math.pi
Gain=1.0
#Radar_Cross_Section_old = 1.04278e10+12

#from Evans 1961 and WInter 1962
Radar_Cross_Section = 0.081 * pi * moon_radius**2


##Distance_to_Moon=3844000000
#below from https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
Distance_to_Moon = 378000000
#Carson's rule: 2 * (deviation + modulating frequency) = 2 * (75 + 15) = 180kHz 
transmitter_BW = 180000.0
MWA_bandwidth = 30720000

#Since the observations were made over the full 30 MHz BW of the MWA 
#each transmitter is only 180kHz in BW so the power is diluted by this factor
dilution_factor = transmitter_BW/MWA_bandwidth

#df=pd.read_csv('moonraker (1).csv')
#df=pd.read_csv('moonraker_test.csv')

#newdf=df[["id","lat", "long", "mhz","erpkw"]]

#newdf = pd.read_csv('moonraker.txt',sep=';',header=2,usecols=["id","lat", "long", "mhz","erpkw"])
newdf = pd.read_csv(fmlist_data_filename,header=2,sep=',',index_col=False,engine='python',usecols=["id","lat","long", "mhz","erpkw"],error_bad_lines=True)

lats = newdf["lat"]
longs = newdf["long"]
ids = newdf["id"]
frequencies = newdf["mhz"]
powers_transmitted_kW = newdf["erpkw"]
powers_transmitted_kW = np.nan_to_num(powers_transmitted_kW)
#randomly, id 2802303 is huge: 2100000000000.000000, remove it
powers_transmitted_kW[powers_transmitted_kW > 1e11] = 0

powers_transmitted_W = powers_transmitted_kW * 1000.
Powers_Recieved=dilution_factor * ((powers_transmitted_W*Gain)/(4*pi*Distance_to_Moon**2))*(Radar_Cross_Section/(4*pi*Distance_to_Moon**2)) 
Powers_Recieved_per_Hz = Powers_Recieved / transmitter_BW
        
#plot frequency distribution
figname="frequency_histogram.png"
plt.clf()
plt.hist(frequencies)
plt.xlabel("frequencies (MHz)")
#plt.ylabel("Flux Density of Moon (Jy)")
plt.savefig(figname)
print("saved %s" % figname)

figname="transmit_powers_histogram.png"
plt.clf()
plt.hist(powers_transmitted_kW)
plt.xlabel("Power Transmitted (kW)")
plt.savefig(figname)
print("saved %s" % figname)

# create new figure, axes instances.

nylat = 40.78
nylon = -73.98

figname='FM_transmitters_on_Earth.png'
fig = plt.figure(figsize=(8, 6), edgecolor='w')
m = Basemap(projection='cyl', resolution=None,
            llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180, )
m.scatter(longs, lats, latlon=True)
#m.scatter(nylon, nylat, latlon=True)
draw_map(m)

# nylat, nylon are lat/lon of New York

fig.savefig(figname,dpdi=500)
print("saved %s" % figname)


if plot_only:
   powers_sum_timestep_array_janskys = np.load(powers_sum_timestep_array_janskys_filename)
   mwa_data = np.load(mwa_data_filename)
   
   freq_100_MHz_index = 22
   mwa_data_100_Mz_specular = mwa_data[freq_100_MHz_index,:,0]
   mwa_obsid_array_100_MHz=mwa_data[freq_100_MHz_index,:,2]
   #work out the UT time of the gps time stamps
   mwa_time_list=[]
   for gps_time in mwa_obsid_array_100_MHz:
      #utc = 1980-01-06UTC + (gps - (leap_count(2014) - leap_count(1980)))
      #from http://maia.usno.navy.mil/ser7/tai-utc.dat
      leap_count_1980=19
      leap_count_2015=36
      if not np.isnan(gps_time):
         utc = datetime(1980, 1, 6) + timedelta(seconds=gps_time - (leap_count_2015 - leap_count_1980))
         utc_for_print_src=utc.strftime('%d/%m/%Y %H:%M:%S')
         mwa_time_list.append(utc)

   
   print mwa_time_list
   print time_delta[45:64]


   figname="moon_flux_density_vs_time_for_earthsine_mwa.png"
   plt.clf()
   plt.plot(mwa_time_list, mwa_data_100_Mz_specular[0:len(mwa_time_list)])
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
   plt.gcf().autofmt_xdate()
   plt.savefig(figname)
   print("saved %s" % figname)
   
   figname="moon_flux_density_vs_time_for_earthsine_sim.png"
   plt.clf()
   plt.plot(time_delta[45:64], powers_sum_timestep_array_janskys[45:64])
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.savefig(figname)
   print("saved %s" % figname)


   sys.exit()
   
#print(newdf["erpkw"])
#print(newdf["erpkw"]*2)





powers_per_timestep_array = np.zeros(len(utc_times))

for timestep_index, timestep in enumerate(utc_times):
     print("timestep %0.2f of %0.2f" % (timestep_index,len(utc_times)))
     print("Time is %s" % timestep)
    
     #If Moon is below the horizon in this tiemstep (at the MRO) then don't bother calculating the rest ....
     frame_MWA=AltAz(obstime=timestep, location=MWA)
     Moon_MWA=get_moon(timestep).transform_to(frame_MWA)
     moon_alt = Moon_MWA.alt
    
     print moon_alt
     
     if moon_alt < 0.:
        continue
     else:   
        city_Locations=EarthLocation(lat=(newdf["lat"])*u.deg, lon=(newdf["long"])*u.deg)
        #print(city_Locations)
   
        city_frames=AltAz(obstime=timestep, location=city_Locations)
        city_moons=get_moon(timestep).transform_to(city_frames)
        
        cities_alt=city_moons.alt
        
        print(len(cities_alt))
        
        powers_received_city_above_horiz = Powers_Recieved_per_Hz[cities_alt > 0]
        
        print(len(powers_received_city_above_horiz))
        
        #print(powers_received_city_above_horiz)
     
        timestep_powers_sum = np.sum(powers_received_city_above_horiz)

        powers_per_timestep_array[timestep_index] = timestep_powers_sum
        




#        #print(powers_sum)
#     powers_sum_timestep_array[timestep_index] = powers_sum
      
#print(powers_sum_timestep_array)  
powers_sum_timestep_array_janskys = powers_per_timestep_array/1e-26

np.save(powers_sum_timestep_array_janskys_filename,powers_sum_timestep_array_janskys)
#print(powers_sum_timestep_array_janskys)

    
#Janskys=[]
#time_delta2=time_delta.value
#for i in powers_sum:
#    a=i/(1e-26)
#    Janskys.append(a)

figname="flux_density_vs_time_for_earthsine.png"
plt.clf()
plt.plot(time_delta, powers_sum_timestep_array_janskys)
plt.xlabel("24 Hour UTC time on 2015-09-26")
plt.ylabel("Flux Density of Moon (Jy)")
plt.savefig(figname)
print("saved %s" % figname)

#plot just where Moon is above horizon



