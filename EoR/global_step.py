#!/usr/bin/env python
#script to calculate sensitivity of interferometers to a global signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cmath
from astropy.wcs.docstrings import phi0
from datetime import datetime, date
from pygsm import GSMObserver
import sys

plot_only = False

#This must be a 'feature' of gsmobserver() - lat and long need to be put in as strings, elevation can be a float (I think)...this is terrible (is my previous gsm_reflection work wrong? rerun?)
mwa_latitude_pyephem = "-26:42.199"
mwa_longitude_pyephem = "116:40.2488"



#in m
#mwa_elevation = 377.83
mwa_elevation = 0

freq_MHz = 100.
#start with an interferometer array layout. E.g. text file with x,y,z coordinates of antennas

#definitions: 
#origin is position of antenna index 0

#read interferometer layout text file (x,y,z in metres)
#interferometer_layout_dir = "/data/code/git/ben-astronomy/EoR/"
#interferometer_layout_filename = "%spanda_locations.txt" % interferometer_layout_dir
#interferometer_layout_filename = "%seight_dipole_locations.txt" % interferometer_layout_dir
#interferometer_layout_filename = "%stwo_dipole_locations.txt" % interferometer_layout_dir


#with open(interferometer_layout_filename,'r') as f:
#   lines = f.readlines()
#locations_list_x = []
#locations_list_y = []
#locations_list_z = []
#for ant_index,line in enumerate(lines):
#   #print "For ant index %s" % ant_index
#   location_x = line.split()[0]
#   location_y = line.split()[1]
#   location_z = line.split()[2]
#   locations_list_x.append(location_x)
#   locations_list_y.append(location_y)
#   locations_list_z.append(location_z)  

#CASE 1: two short dipoles in parallel configuration separated by baseline b
#healpix resolution
NSIDE=512 #(to match pygsm)

baseline_length_m = 1.66
linear_array_ant_separation = 1.66

freq_MHz_array = np.arange(0,300,10)

wavelength_array = 300.0/freq_MHz_array

baseline_length_lambda_array = baseline_length_m/wavelength_array
#print baseline_length_lambda_array[34]
#print baseline_length_lambda_array



#####################
#Now to assess sensitivity to angular structure

#sky map - GSM (use instead of GMOSS as you can use pyGSM with Observer() to generate a sky for different times, allowing you to rotate the sky)
year,month,day,hour,minute,second = 2019,1,1,0,0,0

latitude, longitude, elevation = mwa_latitude_pyephem, mwa_longitude_pyephem, mwa_elevation
#Parkes (pygsm example)
#latitude, longitude, elevation = '-32.998370', '148.263659', 100
ov = GSMObserver()
ov.lon = longitude
ov.lat = latitude
ov.elev = elevation

visibility_real_array_short_para_norm_sum = np.full(len(baseline_length_lambda_array),0.)
visibility_imag_array_short_para_norm_sum = np.full(len(baseline_length_lambda_array),0.)
visibility_real_array_short_para_norm_angular_sum = np.full(len(baseline_length_lambda_array),0.)
visibility_imag_array_short_para_norm_angular_sum = np.full(len(baseline_length_lambda_array),0.)
      
      
# do all this stuff for each hour over 12 hours
for hour_index, hour in enumerate(np.arange(12,24)):
   time_string = "%02d_%02d_%02d" % (hour,minute,second)
   print time_string
   date_obs = datetime(year, month, day, hour, minute, second)
   ov.date = date_obs
   #hack pygsm to do full sky (no masking below horizom)
   gsm_map_full_sky=ov.generate(freq_MHz)
   
   plt.clf()
   map_title="Sky from MWA at h:m:s %02d:%02d:%02d" % (hour,minute,second)
   hp.orthview(map=gsm_map_full_sky,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,0,0))
   #ov.view()
   fig_name="sky_from_mwa_at_h_m_s_%02d_%02d_%02d_full_sky.png" % (hour,minute,second)
   figmap = plt.gcf()
   figmap.savefig(fig_name,dpi=500)
   print "saved %s" % fig_name
   plt.close()
   
   #Want to just seee the response to angular variations on the sky (subtract mean)
   gsm_map_full_sky_mean = np.mean(gsm_map_full_sky)
   gsm_map_full_sky_mean_subtr = gsm_map_full_sky - gsm_map_full_sky_mean
   
   plt.clf()
   map_title="Sky from MWA at h:m:s %02d:%02d:%02d" % (hour,minute,second)
   hp.orthview(map=gsm_map_full_sky_mean_subtr,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,0,0))
   #ov.view()
   fig_name="sky_from_mwa_at_h_m_s_%s_full_sky_mean_subtr.png" % (time_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name,dpi=500)
   print "saved %s" % fig_name
   plt.close()
   
   #global
   #sky map - global - this needs to be the average of the sky map so we can compare with the angular response
   #sky_map =  np.ones(hp.nside2npix(NSIDE))
   sky_map =  np.full(hp.nside2npix(NSIDE),gsm_map_full_sky_mean)
   
   #average over a night  
   #for hour in range(0,24):
   #   date_obs = datetime(year, month, day, hour, minute, second)
   #   ov.date = date_obs
   #   print date_obs
   #   print ov.date
   #   gsm_map=ov.generate(freq_MHz)
   #   
   #   
   #   #d = ov.view(logged=True)
   #   
   #   plt.clf()
   #   map_title="Sky from MWA at h:m:s %02d:%02d:%02d" % (hour,minute,second)
   #   #hp.orthview(map=gsm_map,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,0,0))
   #   ov.view()
   #   fig_name="sky_from_mwa_at_h_m_s_%02d_%02d_%02d.png" % (hour,minute,second)
   #   figmap = plt.gcf()
   #   figmap.savefig(fig_name,dpi=500)
   #   print "saved %s" % fig_name
   #   plt.close()

   visibility_amp_array_iso_filename = "visibility_amp_list_iso_%s.npy" % (time_string)
   visibility_real_array_short_para_filename = "visibility_amp_array_short_para_%s.npy" % (time_string)
   visibility_imag_array_short_para_filename = "visibility_imag_array_short_para_%s.npy" % (time_string)
   visibility_real_array_short_para_angular_filename = "visibility_real_array_short_para_angular_%s.npy" % (time_string)
   visibility_imag_array_short_para_angular_filename = "visibility_imag_array_short_para_angular_%s.npy" % (time_string)
      
   if not plot_only:
      #beam maps
      #isotropic case:
      iso_beam_map = np.ones(hp.nside2npix(NSIDE))
   
   
      #short dipole beam map for parallel case
      short_dipole_parallel_beam_map=np.empty_like(iso_beam_map)
      for hpx_index,beam_value in enumerate(short_dipole_parallel_beam_map):
         theta,phi=hp.pix2ang(NSIDE,hpx_index)
         #populate the dipole model assuming a short dipole (see "New comparison of MWA tile beams" by Benjamin McKinley on Twiki)
         #parallel
         theta_parallel=np.arccos(np.sin(theta)*np.sin(phi))
         voltage_parallel=np.sin(theta_parallel)
         power_parallel=voltage_parallel**2
         short_dipole_parallel_beam_map[hpx_index] = power_parallel
         
         #set sky map to zero below the horizon
         #(don't actually want to do this to compare to Sing as they put it in space with no horizon. But makes no difference here as this is just for global 
         #signal so will be identical above and below horiz (not true for angular structure)))
         #if (theta > np.pi/2.):
         #   sky_map[hpx_index]=0.
         
         #make the arrays vertical by changing where the sky is!
         #if (np.pi/2. < phi < 3*np.pi/2.):
         #   print phi
         #   sky_map[hpx_index]=5.
         #   sky_map[hpx_index]=0.
   
      #gaussian beam for log periodic antennas (21CMA)
      fwhm_deg_x =  82.39
      fwhm_deg_y =  58.71
      sigma_x=fwhm_deg_x/2.355
      sigma_y=fwhm_deg_y/2.355
      gaussian_parallel_beam_map=np.empty_like(iso_beam_map)
      for hpx_index,beam_value in enumerate(gaussian_parallel_beam_map):
         theta,phi=hp.pix2ang(NSIDE,hpx_index)
      
         theta_deg=theta/np.pi*180.
         theta_deg_x = theta_deg*np.cos(phi)
         theta_deg_y = theta_deg*np.sin(phi)
      
         #populate the dipole model assuming a short dipole (see "New comparison of MWA tile beams" by Benjamin McKinley on Twiki)
         #parallel
         #beam_value=np.exp(-4*np.log(2) * ((theta_x-theta_x_pointing)**2 + (theta_y-theta_y_pointing)**2) / fwhm_deg**2)  
         beam_value=np.exp(-( theta_deg_x**2./(2.*sigma_x**2.) + theta_deg_y**2./(2.*sigma_y**2.)))
         gaussian_parallel_beam_map[hpx_index] = beam_value
        
      
      #plot the beam maps:
      beam_map_fitsname="iso_beam_map.fits"
      beam_map_figname="iso_beam_map.png"
      beam_map_title = "iso beam map"
      #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,0,0),return_projected_map=True)
      #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
      #iso_beam_map = hp.read_map(beam_map_fitsname)
      #iso_beam_map_projected=hp.orthview(map=iso_beam_map,return_projected_map=True,coord='C')
      
      #plt.imshow(iso_beam_map_projected)
      #plt.title(beam_map_title)
      #plt.colorbar()
      #figmap = plt.gcf()
      #figmap.savefig(beam_map_figname) 
      #print "saved figure: %s" % beam_map_figname
      #plt.clf()
   
      #Short dipole
      beam_map_fitsname="short_dipole_para_beam_map.fits"
      beam_map_figname="short_dipole_para_beam_map.png"
      beam_map_title = "short dipole para beam map"
      #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,0,0),return_projected_map=True)
      #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
      #iso_beam_map = hp.read_map(beam_map_fitsname)
      short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=beam_map_title,rot=(0,90,90))
      
      plt.imshow(short_dipole_parallel_beam_map_projected)
      plt.title(beam_map_title)
      plt.colorbar()
      figmap = plt.gcf()
      figmap.savefig(beam_map_figname) 
      print "saved figure: %s" % beam_map_figname
      plt.clf()
      
      #log periodic dipole
      beam_map_fitsname="log_periodic_para_beam_map.fits"
      beam_map_figname="log_periodic_para_beam_map.png"
      beam_map_title = "log periodic para beam map"
      #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,0,0),return_projected_map=True)
      #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
      #iso_beam_map = hp.read_map(beam_map_fitsname)
      gaussian_parallel_beam_map_projected=hp.orthview(map=gaussian_parallel_beam_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=beam_map_title,rot=(0,90,90))
      
      plt.imshow(gaussian_parallel_beam_map_projected)
      plt.title(beam_map_title)
      plt.colorbar()
      figmap = plt.gcf()
      figmap.savefig(beam_map_figname) 
      print "saved figure: %s" % beam_map_figname
      plt.clf()
   
      visibility_element_array_iso = np.empty_like(iso_beam_map,dtype=complex)
      visibility_amp_list_iso = []
      visibility_element_array_short_para = np.empty_like(iso_beam_map,dtype=complex)
      visibility_real_list_short_para = []
      visibility_imag_list_short_para = []
      visibility_element_array_log= np.empty_like(iso_beam_map,dtype=complex)
      visibility_amp_list_log = []
      visibility_element_array_short_para_four = np.empty_like(iso_beam_map,dtype=complex)
      visibility_amp_list_short_para_four = []
      visibility_element_array_log_four = np.empty_like(iso_beam_map,dtype=complex)
      visibility_amp_list_log_four = []
      
      visibility_element_array_short_para_angular = np.empty_like(iso_beam_map,dtype=complex)
      visibility_real_list_short_para_angular = []
      visibility_imag_list_short_para_angular = []
      
      for baseline_length_lambda_index,baseline_length_lambda in enumerate(baseline_length_lambda_array):
         print baseline_length_lambda
      
         baseline_theta = np.pi/2.
      
         baseline_phi = 0.
         baseline_vector=baseline_length_m*hp.ang2vec(baseline_theta,baseline_phi)
         #print baseline_vector
         #baseline_vector=np.array([baseline_length_m,0,0])
         #print np.linalg.norm(baseline_vector)
         #print baseline_length_lambda
      
         #for each sky position pixel:
         for hpx_index,beam_value in enumerate(iso_beam_map):
            sky_vector=hp.pix2vec(NSIDE,hpx_index)
            #print np.linalg.norm(sky_vector)
            b_dot_r = np.dot(baseline_vector,sky_vector)
            #print b_dot_r
            #print b_dot_r
            phase_angle = 2.*np.pi*b_dot_r/wavelength_array[baseline_length_lambda_index]
            #print phase_angle
            
          
            element_iso = sky_map[hpx_index]*iso_beam_map[hpx_index]*np.exp(-1j*phase_angle)
            visibility_element_array_iso[hpx_index] = element_iso
            element_short_para = sky_map[hpx_index]*short_dipole_parallel_beam_map[hpx_index]*np.exp(-1j*phase_angle)
            visibility_element_array_short_para[hpx_index] = element_short_para
            element_log = sky_map[hpx_index]*gaussian_parallel_beam_map[hpx_index]*np.exp(-1j*phase_angle)
            visibility_element_array_log[hpx_index] = element_log  
            
            #Response to angular structure:
            element_short_para_angular = gsm_map_full_sky_mean_subtr[hpx_index]*short_dipole_parallel_beam_map[hpx_index]*np.exp(-1j*phase_angle)
            visibility_element_array_short_para_angular[hpx_index] = element_short_para_angular
         
            #forget arrays for now
            ##now do array:
            #theta,phi=hp.pix2ang(NSIDE,hpx_index)
            ##psi is the angle between the sky and the long axis of the 1-D array
            ##what happens if you don't fix the antennas at lambda/D spacing ... ie new array factor for each frequency:
            #k=2.*np.pi/wavelength_array[baseline_length_lambda_index]
            ##print k
            #psi_parallel=k*linear_array_ant_separation*np.cos(np.arccos(np.sin(theta)*np.sin(phi)))
            ##psi_parallel=np.pi*np.cos(np.arccos(np.sin(theta)*np.sin(phi)))
            ##number of elements in array
            #N=4.
            #array_factor = (1./N)*((np.sin(N*psi_parallel/2.))/(np.sin(psi_parallel/2.)))
            #element_short_para_four = sky_map[hpx_index]*short_dipole_parallel_beam_map[hpx_index]*array_factor*np.exp(-1j*phase_angle)
            #visibility_element_array_short_para_four[hpx_index] = element_short_para_four      
            #element_log_four = sky_map[hpx_index]*gaussian_parallel_beam_map[hpx_index]*array_factor*np.exp(-1j*phase_angle)
            #visibility_element_array_log_four[hpx_index] = element_log_four
            
         
         #integrate across the whole sky (using healpix so it is just a sum)
         visibility_iso = np.sum(visibility_element_array_iso)
         visibility_real_iso = visibility_iso.real
         visibility_amp_list_iso.append(visibility_real_iso)
      
         visibility_log = np.sum(visibility_element_array_log)
         visibility_real_log = visibility_log.real
         visibility_amp_list_log.append(visibility_real_log)
            
         visibility_short_para = np.sum(visibility_element_array_short_para)
         visibility_real_short_para = visibility_short_para.real
         visibility_imag_short_para = visibility_short_para.imag
         visibility_real_list_short_para.append(visibility_real_short_para)
         visibility_imag_list_short_para.append(visibility_imag_short_para)
      
         visibility_short_para_angular = np.sum(visibility_element_array_short_para_angular)
         visibility_real_short_para_angular = visibility_short_para_angular.real
         visibility_imag_short_para_angular = visibility_short_para_angular.imag
         visibility_real_list_short_para_angular.append(visibility_real_short_para_angular) 
         visibility_imag_list_short_para_angular.append(visibility_imag_short_para_angular) 
         
         
         #visibility_short_para_four = np.sum(visibility_element_array_short_para_four)
         #visibility_real_short_para_four = visibility_short_para_four.real
         #visibility_amp_list_short_para_four.append(visibility_real_short_para_four)
      
         #visibility_log_four = np.sum(visibility_element_array_log_four)
         #visibility_real_log_four = visibility_log_four.real
         #visibility_amp_list_log_four.append(visibility_real_log_four)

      #save as numpy arrays!
      visibility_amp_array_iso = np.asarray(visibility_amp_list_iso)
      visibility_real_array_short_para = np.asarray(visibility_real_list_short_para)
      visibility_imag_array_short_para = np.asarray(visibility_imag_list_short_para)
      visibility_real_array_short_para_angular = np.asarray(visibility_real_list_short_para_angular)
      visibility_imag_array_short_para_angular = np.asarray(visibility_imag_list_short_para_angular)
      
      np.save(visibility_amp_array_iso_filename,visibility_amp_array_iso)
      np.save(visibility_real_array_short_para_filename,visibility_real_array_short_para)
      np.save(visibility_imag_array_short_para_filename,visibility_imag_array_short_para)
      np.save(visibility_real_array_short_para_angular_filename,visibility_real_array_short_para_angular)
      np.save(visibility_imag_array_short_para_angular_filename,visibility_imag_array_short_para_angular)
      
   if plot_only:
      visibility_amp_array_iso = np.load(visibility_amp_array_iso_filename)
      visibility_real_array_short_para = np.load(visibility_real_array_short_para_filename)
      visibility_imag_array_short_para = np.load(visibility_imag_array_short_para_filename)
      visibility_real_array_short_para_angular = np.load(visibility_real_array_short_para_angular_filename)
      visibility_imag_array_short_para_angular = np.load(visibility_imag_array_short_para_angular_filename)
   
   #normalise to global value at b=0   
   visibility_amp_array_iso_norm = visibility_amp_array_iso/visibility_amp_array_iso[0]
   visibility_real_array_short_para_norm = visibility_real_array_short_para/visibility_real_array_short_para[0]
   visibility_real_array_short_para_norm_abs = np.abs(visibility_real_array_short_para_norm)
   visibility_imag_array_short_para_norm = visibility_imag_array_short_para/visibility_real_array_short_para[0]
   visibility_imag_array_short_para_norm_abs = np.abs(visibility_imag_array_short_para_norm)
   
   
   #normalise to max of real for angular
   visibility_real_array_short_para_angular_max = np.max(visibility_real_array_short_para_angular)
   visibility_real_array_short_para_angular_norm = visibility_real_array_short_para_angular/visibility_real_array_short_para[0]
   visibility_real_array_short_para_angular_norm_abs = np.abs(visibility_real_array_short_para_angular_norm)
   
   #same for imag
   visibility_imag_array_short_para_angular_max = np.max(visibility_imag_array_short_para_angular)
   visibility_imag_array_short_para_angular_norm = visibility_imag_array_short_para_angular/visibility_real_array_short_para[0]
   visibility_imag_array_short_para_angular_norm_abs = np.abs(visibility_imag_array_short_para_angular_norm)
   

   #Add to sums for averaging
   visibility_real_array_short_para_norm_sum += visibility_real_array_short_para_norm
   visibility_imag_array_short_para_norm_sum += visibility_imag_array_short_para_norm
   visibility_real_array_short_para_norm_angular_sum += visibility_real_array_short_para_angular_norm
   visibility_imag_array_short_para_norm_angular_sum += visibility_imag_array_short_para_angular_norm   
   
   #calc average
   visibility_real_array_short_para_norm_average = visibility_real_array_short_para_norm_sum / (hour_index+1)
   visibility_imag_array_short_para_norm_average = visibility_imag_array_short_para_norm_sum / (hour_index+1)
   visibility_real_array_short_para_angular_norm_average = visibility_real_array_short_para_norm_angular_sum / (hour_index+1)
   visibility_imag_array_short_para_angular_norm_average = visibility_imag_array_short_para_norm_angular_sum / (hour_index+1)

   #abs
   visibility_real_array_short_para_norm_average_abs = np.abs(visibility_real_array_short_para_norm_average)
   visibility_imag_array_short_para_norm_average_abs = np.abs(visibility_imag_array_short_para_norm_average)
   visibility_real_array_short_para_angular_norm_average_abs = np.abs(visibility_real_array_short_para_angular_norm_average)
   visibility_imag_array_short_para_angular_norm_average_abs = np.abs(visibility_imag_array_short_para_angular_norm_average)   
   
   #visibility_amp_list_short_para_four = visibility_amp_list_short_para_four/np.nanmax(visibility_amp_list_short_para_four)
   #visibility_amp_array_log = visibility_amp_array_log/np.nanmax(visibility_amp_array_log)
   #visibility_amp_list_log_four = visibility_amp_list_log_four/np.nanmax(visibility_amp_list_log_four)
   
   #for 21CMA, antenna spacing is 1.66 m, so for this to be lambda/2 the lambda = 3.32, frequency = 90.36 MHz, and the baseline length of 3.7 m = 1.114 lambda
   #baseline_length_21CMA_90_MHz_lambda = 1.114
   #vis_amp_21CMA_90_MHz = visibility_amp_list_short_para_four[34]
   #print "vis_amp_21CMA_90_MHz %s" % vis_amp_21CMA_90_MHz
   
   plt.clf()
   plot=plt.figure()
   plot_title="Visibility amplitude vs baseline length"
   #plot_title="Vis ampl vs freq baseline:%0.1fm and sep %0.1fm" % (baseline_length_m,linear_array_ant_separation)
   plot_figname="visibility_amp_vs_baseline_length_%s.png" % (time_string)
   #plot_figname="visibility_amp_vs_frequency_baseline_%0.1fm_ant_sep_%0.1fm.png" % (baseline_length_m,linear_array_ant_separation)
   #plt.errorbar(freq_array_band,Tb_band_1,yerr=Tb_error_band_1,label="Measured")
   #plt.plot(baseline_length_lambda_array,visibility_amp_array_iso_norm,label="single isotropic")
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_norm,label="global real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_norm,label="global imag")
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_angular_norm,label="angular real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_angular_norm,label="angular imag")
   
   #plt.plot(baseline_length_lambda_array,visibility_amp_list_short_para_four,label="four short para")
   #plt.plot(freq_MHz_array,visibility_amp_list_iso,label="single isotropic")
   #plt.plot(freq_MHz_array,visibility_amp_list_log,label="single log-periodic")
   #plt.plot(freq_MHz_array,visibility_amp_list_short_para,label="single short para")
   #plt.plot(freq_MHz_array,visibility_amp_list_short_para_four,label="four short para")
   #plt.plot(freq_MHz_array,visibility_amp_list_log_four,label="four log-periodic")
   
   #plt.scatter(baseline_length_21CMA_90_MHz_lambda,vis_amp_21CMA_90_MHz,label="21 CMA short para 90 MHz",marker='+',c='r')
   
   plt.title(plot_title)
   plt.ylabel('Visibility amplitude')
   plt.xlabel('Baseline length (lambda)')
   #plt.xlabel('Frequency (MHz)')
   #plt.xticks(np.arange(min(baseline_length_lambda_array), max(baseline_length_lambda_array)+1, 1.0))
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(-0.8, 1.1)
   plot.savefig('%s' % plot_figname)
   print "saved %s" % plot_figname
   plt.close()   

   #abs
   plt.clf()
   plot=plt.figure()
   plot_title="Visibility amplitude abs vs baseline length"
   #plot_title="Vis ampl vs freq baseline:%0.1fm and sep %0.1fm" % (baseline_length_m,linear_array_ant_separation)
   plot_figname="visibility_amp_abs_vs_baseline_length_%s.png" % (time_string)
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_norm_abs,label="global real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_norm_abs,label="global imag")
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_angular_norm_abs,label="angular real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_angular_norm_abs,label="angular imag")
   
   plt.title(plot_title)
   plt.ylabel('Visibility amplitude (abs)')
   plt.xlabel('Baseline length (lambda)')
   #plt.xlabel('Frequency (MHz)')
   #plt.xticks(np.arange(min(baseline_length_lambda_array), max(baseline_length_lambda_array)+1, 1.0))
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(0, 1.1)
   plot.savefig('%s' % plot_figname)
   print "saved %s" % plot_figname
   plt.close() 
   
   #average
   plt.clf()
   plot=plt.figure()
   plot_title="Visibility amplitude vs baseline length"
   plot_figname="visibility_amp_vs_baseline_length_average_to_%s.png" % (time_string)
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_norm_average,label="global real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_norm_average,label="global imag")
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_angular_norm_average,label="angular real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_angular_norm_average,label="angular imag")
   
   plt.title(plot_title)
   plt.ylabel('Visibility amplitude')
   plt.xlabel('Baseline length (lambda)')
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(-0.8, 1.1)
   plot.savefig('%s' % plot_figname)
   print "saved %s" % plot_figname
   plt.close()   

   #average abs
   plt.clf()
   plot=plt.figure()
   plot_title="Visibility amplitude abs vs baseline length"
   plot_figname="visibility_amp_abs_vs_baseline_length_average_to_%s.png" % (time_string)
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_norm_average_abs,label="global real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_norm_average_abs,label="global real")
   plt.plot(baseline_length_lambda_array,visibility_real_array_short_para_angular_norm_average_abs,label="angular real")
   plt.plot(baseline_length_lambda_array,visibility_imag_array_short_para_angular_norm_average_abs,label="angular imag")
   
   plt.title(plot_title)
   plt.ylabel('Visibility amplitude (abs)')
   plt.xlabel('Baseline length (lambda)')
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(0, 1.1)
   plot.savefig('%s' % plot_figname)
   print "saved %s" % plot_figname
   plt.close() 
   
   
   
   
   
   
   

