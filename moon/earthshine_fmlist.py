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
from datetime import datetime, timedelta, date
import matplotlib.dates as mdates
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
from astropy.io import fits
import os

tile_only = False

def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)


def draw_map(m,scale=1):
    # draw a shaded-relief image
    #m.shadedrelief(scale=scale)
    m.bluemarble(scale=scale)
    #m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    
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

def extract_moon_from_mwa_freq_images(image_dir,mwa_image_name_base_list,nchans,moon_pix,start_freq,freq_step):
   freq_list = []
   moon_flux_peak_list = []
   for image_name_base_index,image_name_base in enumerate(mwa_image_name_base_list):
      freq_base = start_freq + float(image_name_base_index) * (nchans*freq_step)
      for image_number in range(0,nchans):
         freq = freq_base + (float(image_number) * freq_step)
         freq_list.append(freq)
         image_name = "%s%s-%04d-image.fits" % (image_dir,image_name_base,image_number)
         #print image_name
         with fits.open(image_name) as hdulist:
            data = hdulist[0].data[0,0,:,:]
            #print data.shape
            moon_flux_peak = data[int(moon_pix[1]),int(moon_pix[0])]
            #print moon_flux_peak
            moon_flux_peak_list.append(moon_flux_peak)
   
   freq_array = np.asarray(freq_list)
   moon_flux_peak_array = np.asarray(moon_flux_peak_list)  
   
   figname="mwa_moon_flux_peak_vs_freq.png"
   plt.clf()
   plt.plot(freq_array,moon_flux_peak_array)
   plt.xlabel("frequency (MHz)")
   plt.ylabel("Peak flux density of Moon (Jy)")
   plt.savefig(figname)
   print("saved %s" % figname)
   plt.close()
         

start_freq = 72.355
freq_step = 0.200
moon_pix = (525,1107)
image_dir = "/md0/lspringer/mwa_data/new_download/"
nchans = 154
mwa_image_name_base_list = ['1127301352_moon_69','1127301592_moon_93']
extract_moon_from_mwa_freq_images(image_dir,mwa_image_name_base_list,nchans,moon_pix,start_freq,freq_step)
#sys.exit()
        
plot_mwa_data = False

fmlist_data_filename = '/md0/moon/lauren/moonraker.csv'
#These are numpy arrays. When you open them with python numpy you will see they have a shape of (124,50,3). This corresponds to (n_channels, n_observations, 3), the 3 is because there is specular, diffuse, and observation ID.
mwa_data_filename = '/md0/moon/lauren/big_RFI_vs_time_freq_2015B_05.npy'

time_delta=np.linspace(0, 24, 96)*u.hour
start_time_utc = Time('2015-09-26 00:00:00')
utc_times=start_time_utc+time_delta
MWA=EarthLocation(lat=-26.7*u.deg, lon=116.67*u.deg)

if tile_only:
   for timestep_index, timestep in enumerate(utc_times):
     
      timestep_string = timestep.datetime.strftime("%Y_%m_%d_%H_%M_%S")

      print("timestep %0.2f of %0.2f" % (timestep_index,len(utc_times)))
      
      flux_figname="flux_density_vs_time_for_earthsine_%s.png" % timestep_string
      transmitters_figname='FM_transmitters_on_Earth_timestep_%s.png' % (timestep_string)
      moon_figname="moon_alt_at_MWA_%s.png"  % (timestep_string)
      top_figname = "top_%s.png" % (timestep_string)

      #Tile images together at each timestep
      top_list_im = [flux_figname, moon_figname]
      
      imgs    = [ Image.open(i) for i in top_list_im ]
      # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
      min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
      
      #horiz:
      imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
      imgs_comb = Image.fromarray( imgs_comb)
      imgs_comb.save(top_figname)   
      
      #then combine top and bottom vertically
      # for a vertical stacking it is simple: use vstack
      bottom_list_im = [top_figname,transmitters_figname]
      imgs    = [ Image.open(i) for i in bottom_list_im ]
      min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
      imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
      imgs_comb = Image.fromarray( imgs_comb)
      trifecta_name = 'trifecta_%03d.png' % (timestep_index) 
      imgs_comb.save(trifecta_name)

      
      ## save that beautiful picture
      #imgs_comb = Image.fromarray( imgs_comb)
      #imgs_comb.save( 'trifecta.png' )    
      

      
   sys.exit() 


powers_sum_timestep_array_janskys_filename = "powers_sum_timestep_array_janskys.npy"
moon_alt_array_filename = "moon_alt_array.npy"

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
plt.close()

figname="transmit_powers_histogram.png"
plt.clf()
plt.hist(powers_transmitted_kW)
plt.xlabel("Power Transmitted (kW)")
plt.savefig(figname)
print("saved %s" % figname)
plt.close()

figname='FM_transmitters_on_Earth.png'
fig = plt.figure(figsize=(100, 75), edgecolor='w')
m = Basemap(projection='cyl', resolution=None,
            llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180, )
m.scatter(longs, lats,s=2, marker='o',color='orange',latlon=True)
#m.scatter(nylon, nylat, latlon=True)
draw_map(m)

fig.savefig(figname,dpdi=500,bbox_inches='tight')
print("saved %s" % figname)
plt.close()

  

if plot_mwa_data:
   powers_sum_timestep_array_janskys = np.load(powers_sum_timestep_array_janskys_filename)
   mwa_data = np.load(mwa_data_filename)
   
   freq_100_MHz_index = 22
   mwa_data_100_Mz_specular = mwa_data[freq_100_MHz_index,:,0]
   mwa_data_100_Mz_diffuse = mwa_data[freq_100_MHz_index,:,1]
   mwa_data_100_Mz_total = mwa_data_100_Mz_specular + mwa_data_100_Mz_diffuse
   mwa_obsid_array_100_MHz = mwa_data[freq_100_MHz_index,:,2]
   #work out the UT time of the gps time stamps
   mwa_time_list=[]
   sim_time_list = []
   powers_sum_timestep_array_janskys_overplot = []
   for gps_time in mwa_obsid_array_100_MHz:
      #utc = 1980-01-06UTC + (gps - (leap_count(2014) - leap_count(1980)))
      #from http://maia.usno.navy.mil/ser7/tai-utc.dat
      leap_count_1980=19
      leap_count_2015=36
      if not np.isnan(gps_time):
         utc = datetime(1980, 1, 6) + timedelta(seconds=gps_time - (leap_count_2015 - leap_count_1980))
         utc_for_print_src=utc.strftime('%d/%m/%Y %H:%M:%S')
         mwa_time_list.append(utc)
         print utc
         sim_index = nearest_ind(utc_times,Time(utc))
         powers_sum_timestep_array_janskys_overplot.append(powers_sum_timestep_array_janskys[sim_index])
         sim_time_min = (sim_index + 1) * 15
         sim_utc =  datetime(2015, 9, 26) + timedelta(minutes=sim_time_min)
         sim_time_list.append(sim_utc)

   print mwa_time_list
   print sim_time_list


   figname="moon_flux_density_vs_time_for_earthsine_mwa.png"
   plt.clf()
   plt.plot(mwa_time_list, mwa_data_100_Mz_total[0:len(mwa_time_list)])
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
   plt.gcf().autofmt_xdate()
   plt.savefig(figname)
   print("saved %s" % figname)
   plt.close()
   
   figname="moon_flux_density_vs_time_for_earthsine_sim.png"
   plt.clf()
   plt.plot(sim_time_list, powers_sum_timestep_array_janskys_overplot)
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
   plt.gcf().autofmt_xdate()
   plt.savefig(figname)
   print("saved %s" % figname)
   plt.close()
   
   figname="moon_flux_density_vs_time_for_earthsine_sim_and_data.png"
   plt.clf()
   plt.plot(sim_time_list, powers_sum_timestep_array_janskys_overplot,label='sims')
   plt.plot(mwa_time_list, mwa_data_100_Mz_total[0:len(mwa_time_list)],label='data')
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.legend(loc='upper left')
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
   plt.gcf().autofmt_xdate()
   plt.savefig(figname)
   print("saved %s" % figname)
   plt.close()

   figname="flux_density_vs_time_for_earthsine_24hrs.png"
   plt.clf()
   plt.plot(time_delta, powers_sum_timestep_array_janskys)
   plt.xlabel("24 Hour UTC time on 2015-09-26")
   plt.ylabel("Flux Density of Moon (Jy)")
   plt.savefig(figname)
   print("saved %s" % figname)
   plt.close()
   
   sys.exit() 
   
#print(newdf["erpkw"])
#print(newdf["erpkw"]*2)





#powers_per_timestep_array = np.zeros(len(utc_times))
#moon_alt_per_timestep = np.zeros(len(utc_times))
powers_per_timestep_array = np.full(len(utc_times), np.nan)
moon_alt_per_timestep = np.full(len(utc_times), np.nan)
zero_alt_timestep = np.zeros(len(utc_times))

#for waterfall plot
n_bins = 100
waterfall = np.full((len(utc_times),n_bins),np.nan)


for timestep_index, timestep in enumerate(utc_times):
     
     timestep_string = timestep.datetime.strftime("%Y_%m_%d_%H_%M_%S")

     print("timestep %0.2f of %0.2f" % (timestep_index,len(utc_times)))
     print("Time is %s" % timestep)
    
     #If Moon is below the horizon in this tiemstep (at the MRO) then don't bother calculating the rest ....
     frame_MWA=AltAz(obstime=timestep, location=MWA)
     Moon_MWA=get_moon(timestep).transform_to(frame_MWA)
     moon_alt = Moon_MWA.alt
     
     moon_alt_per_timestep[timestep_index] = moon_alt/u.deg
     
     print moon_alt
     
     
     
     if moon_alt < 0.:
        
        powers_per_timestep_array[timestep_index] = 0.
        
        #make progressive plots (zero)
        powers_sum_timestep_array_janskys = powers_per_timestep_array/1e-26
        figname="flux_density_vs_time_for_earthsine_%s.png" % timestep_string
        plt.clf()
        fig = plt.figure(figsize=(8, 6))
        plt.plot(time_delta, powers_sum_timestep_array_janskys)
        plt.xlabel("24 Hour UTC time on 2015-09-26")
        plt.ylabel("Flux Density of Moon (Jy)")
        ax = plt.gca()
        ax.set_ylim(0, 460)
        ax.set_xlim(0, 24)
        #plt.ylim([0, 460])
        plt.savefig(figname)
        print("saved %s" % figname)
        plt.close()
        
        plt.clf()
        figname='FM_transmitters_on_Earth_timestep_%s.png' % (timestep_string)
        fig = plt.figure(figsize=(100, 75), edgecolor='w')
        m = Basemap(projection='cyl', resolution=None,
                    llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=-180, urcrnrlon=180, )
        m.scatter(longs, lats,s=2, marker='o',color='orange',latlon=True)
        draw_map(m)
        
        fig.savefig(figname,dpdi=500,bbox_inches='tight')
        print("saved %s" % figname)
        
        #plot Moon alt
        figname="moon_alt_at_MWA_%s.png"  % (timestep_string)
        plt.clf()
        fig = plt.figure(figsize=(8, 6))
        plt.plot(time_delta, moon_alt_per_timestep)
        plt.plot(time_delta, zero_alt_timestep)
        plt.xlabel("24 Hour UTC time on 2015-09-26")
        plt.ylabel("Moon altitude (deg)")
        ax = plt.gca()
        ax.set_ylim(-80, 80)
        ax.set_xlim(0, 24)
        plt.savefig(figname)
        print("saved %s" % figname)
        plt.close()
        
        powers_received_city_above_horiz = Powers_Recieved_per_Hz * 0.0
        powers_received_city_above_horiz_freqs = frequencies
        
     else:   
        city_Locations=EarthLocation(lat=(newdf["lat"])*u.deg, lon=(newdf["long"])*u.deg)
        #print(city_Locations)
   
        city_frames=AltAz(obstime=timestep, location=city_Locations)
        city_moons=get_moon(timestep).transform_to(city_frames)
        
        cities_alt=city_moons.alt
        
        print(len(cities_alt))
        
        powers_received_city_above_horiz = Powers_Recieved_per_Hz[np.nan_to_num(cities_alt) > 0]
        print(len(powers_received_city_above_horiz))
        
        powers_received_city_above_horiz_freqs = frequencies[np.nan_to_num(cities_alt) > 0]
        #print(powers_received_city_above_horiz)
     
        timestep_powers_sum = np.sum(powers_received_city_above_horiz)

        powers_per_timestep_array[timestep_index] = timestep_powers_sum
        
        powers_sum_timestep_array_janskys = powers_per_timestep_array/1e-26
        
        #make progressive plots
        figname="flux_density_vs_time_for_earthsine_%s.png" % timestep_string
        plt.clf()
        fig = plt.figure(figsize=(8, 6))
        plt.plot(time_delta, powers_sum_timestep_array_janskys)
        plt.xlabel("24 Hour UTC time on 2015-09-26")
        plt.ylabel("Flux Density of Moon (Jy)")
        ax = plt.gca()
        ax.set_ylim(0, 460)
        ax.set_xlim(0, 24)
        plt.savefig(figname)
        print("saved %s" % figname)
        plt.close()


        ##plot cities on the map that have Moon above horiz.
        longs_above = longs[cities_alt > 0]
        lats_above = lats[cities_alt > 0]
        longs_below = longs[cities_alt <= 0]
        lats_below = lats[cities_alt <= 0]
        figname='FM_transmitters_on_Earth_timestep_%s.png' % (timestep_string)
        plt.clf()
        fig = plt.figure(figsize=(100, 75), edgecolor='w')
        m = Basemap(projection='cyl', resolution=None,
                    llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=-180, urcrnrlon=180, )
        longs_above_m, lats_above_m = m(longs_above, lats_above)
        longs_below_m, lats_below_m = m(longs_below, lats_below)
        m.scatter(longs_above_m, lats_above_m,s=2, marker='o',color='red')
        m.scatter(longs_below_m, lats_below_m,s=2, marker='o',color='orange')
        draw_map(m)
        
        fig.savefig(figname,dpdi=500,bbox_inches='tight')
        print("saved %s" % figname)
        plt.close()

        #plot Moon alt
        figname="moon_alt_at_MWA_%s.png"  % (timestep_string)
        plt.clf()
        fig = plt.figure(figsize=(8, 6))
        plt.plot(time_delta, moon_alt_per_timestep)
        plt.plot(time_delta, zero_alt_timestep)
        plt.xlabel("24 Hour UTC time on 2015-09-26")
        plt.ylabel("Moon altitude (deg)")
        ax = plt.gca()
        ax.set_ylim(-80, 80)
        ax.set_xlim(0, 24)
        plt.savefig(figname)
        print("saved %s" % figname)
        plt.close()
     
     #Tile images together at each timestep
     
     flux_figname="flux_density_vs_time_for_earthsine_%s.png" % timestep_string
     transmitters_figname='FM_transmitters_on_Earth_timestep_%s.png' % (timestep_string)
     moon_figname="moon_alt_at_MWA_%s.png"  % (timestep_string)
     top_figname = "top_%s.png" % (timestep_string)
      
     top_list_im = [flux_figname, moon_figname]
     
     imgs    = [ Image.open(i) for i in top_list_im ]
     # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
     min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
     
     #horiz:
     imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
     imgs_comb = Image.fromarray( imgs_comb)
     imgs_comb.save(top_figname)   
     
     #then combine top and bottom vertically
     # for a vertical stacking it is simple: use vstack
     bottom_list_im = [top_figname,transmitters_figname]
     imgs    = [ Image.open(i) for i in bottom_list_im ]
     min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
     imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
     imgs_comb = Image.fromarray( imgs_comb)
     trifecta_name = 'trifecta_%03d.png' % (timestep_index) 
     imgs_comb.save(trifecta_name) 
     print "saved %s" % (trifecta_name)
     
     #waterfall stuff
     n,bins,patches = plt.hist(frequencies,bins=n_bins)
     power_at_freqs_array=np.full(n_bins,0.0)
     power_at_freq =0.0
     for i in range(0, n_bins):
        minfreq=bins[i]
        #print 'min freq %s' % (minfreq)
        maxfreq=bins[i+1]
        #print 'max freq %s' % (maxfreq)
        freq_indices = np.where(np.logical_and(powers_received_city_above_horiz_freqs>=minfreq, powers_received_city_above_horiz_freqs<maxfreq))
        #print freq_indices
        power_at_freq = np.sum(np.nan_to_num(np.array(powers_received_city_above_horiz)[freq_indices]))
        #print power_at_freq
        power_at_freqs_array[i]=power_at_freq
     
     waterfall[timestep_index,:]=power_at_freqs_array
     waterfall_jy = waterfall/10**(-26)
     
     #save a progressive waterfall figure
     filename = 'waterfall_%s.fits' % timestep_string
     fits.writeto(filename,waterfall_jy,overwrite=True)
     print "wrote %s" % filename
     
     plot_figname = 'waterfall_%s.png' % timestep_string
     fig = plt.figure()
     cur_axes = plt.gca()
     img = cur_axes.matshow(waterfall_jy,vmin=0,vmax=10)
     plt.grid(None) 
     #get rid of tick labels
     
     cur_axes.axes.get_xaxis().set_ticklabels([])
     cur_axes.axes.get_yaxis().set_ticklabels([])
     #get rid of ticks
     cur_axes = plt.gca()
     cur_axes.axes.get_xaxis().set_ticks([])
     cur_axes.axes.get_yaxis().set_ticks([])
     colorbar = plt.colorbar(img, ax=cur_axes)
     colorbar.set_label('Flux density (Jy)')
     
     plt.xlabel('Frequency')
     plt.ylabel('Time')
     
     plt.savefig('%s' % plot_figname, overwrite=True)
     print "saved %s" % plot_figname
     plt.close()

     #tile the waterfalls with the earth map and moon alt
     
     waterfall_figname='waterfall_%s.png' % timestep_string
     transmitters_figname='FM_transmitters_on_Earth_timestep_%s.png' % (timestep_string)
     moon_figname="moon_alt_at_MWA_%s.png"  % (timestep_string)
     top_figname = "top_%s.png" % (timestep_string)
      
     top_list_im = [waterfall_figname, moon_figname]
     
     imgs    = [ Image.open(i) for i in top_list_im ]
     # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
     min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
     
     #horiz:
     imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
     imgs_comb = Image.fromarray( imgs_comb)
     imgs_comb.save(top_figname)   
     
     #then combine top and bottom vertically
     # for a vertical stacking it is simple: use vstack
     bottom_list_im = [top_figname,transmitters_figname]
     imgs    = [ Image.open(i) for i in bottom_list_im ]
     min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
     imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
     imgs_comb = Image.fromarray( imgs_comb)
     trifecta_name = 'waterfall_trifecta_%03d.png' % (timestep_index) 
     imgs_comb.save(trifecta_name) 
     print "saved %s" % (trifecta_name)
     

#print(powers_sum_timestep_array)  
powers_sum_timestep_array_janskys = powers_per_timestep_array/1e-26

np.save(powers_sum_timestep_array_janskys_filename,powers_sum_timestep_array_janskys)
np.save(moon_alt_array_filename,moon_alt_per_timestep)
#print(powers_sum_timestep_array_janskys)

    
#Janskys=[]
#time_delta2=time_delta.value
#for i in powers_sum:
#    a=i/(1e-26)
#    Janskys.append(a)

figname="flux_density_vs_time_for_earthsine.png"
plt.clf()
fig = plt.figure(figsize=(8, 6))
plt.plot(time_delta, powers_sum_timestep_array_janskys)
plt.xlabel("24 Hour UTC time on 2015-09-26")
plt.ylabel("Flux Density of Moon (Jy)")
plt.ylim=[0,460]
plt.savefig(figname)
print("saved %s" % figname)

#plot Moon alt
figname="moon_alt_at_MWA.png"
plt.clf()
fig = plt.figure(figsize=(8, 6))
plt.plot(time_delta, moon_alt_per_timestep)
plt.xlabel("24 Hour UTC time on 2015-09-26")
plt.ylabel("Moon altitude (deg)")
ax = plt.gca()
ax.set_ylim(-80, 80)
ax.set_xlim(0, 24)
plt.savefig(figname)
print("saved %s" % figname)
#plot just where Moon is above horizon

#stitch together images into movie

filename_base = 'trifecta_%03d.png'
movie_name = 'earthshine.mp4'

cmd = "ffmpeg -framerate 6 -i %s -c:v libx264 -r 30 -pix_fmt yuv420p %s" % (filename_base,movie_name)
print cmd
os.system(cmd)
print('made movie %s' % movie_name)

#stitch together images into movie waterfall version

filename_base = 'waterfall_trifecta_%03d.png'

#draw labels on waterfall trifecta
#cmd = " convert waterfall_trifecta_000.png   -background Orange  label:'Faerie Dragon' +swap  -gravity Center -append    anno_label2.png"

movie_name = 'earthshine_waterfall.mp4'

cmd = "ffmpeg -framerate 6 -i %s -c:v libx264 -r 30 -pix_fmt yuv420p %s" % (filename_base,movie_name)
print cmd
os.system(cmd)
print('made movie %s' % movie_name)

