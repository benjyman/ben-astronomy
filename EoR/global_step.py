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
import math
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
import os


def plot_different_baselines(antenna_locations_filename):
   #print antenna_locations_filename
   
   #1. calc u,v  from antenna positions
   with open(antenna_locations_filename) as f:
      lines = f.readlines()
   n_ants = len(lines)
   #print n_ants
   ant_i_array = np.arange(0,n_ants)
   ant_j_array = np.arange(0,n_ants)
   
   n_baselines = int(n_ants*(n_ants-1.) / 2.)
   #print n_baselines
   
   #visibility array n_baselines * u,v, complex vis
   visibility_array = np.full([n_baselines,3],(np.nan))
   
   visibility_index = 0
   for ant_i in ant_i_array:
      for ant_j in ant_j_array:
         if (ant_j > ant_i):
            ant_i_x = float(lines[int(ant_i)].split()[0])
            ant_j_x = float(lines[int(ant_j)].split()[0])
            u = ant_i_x - ant_j_x
            
            ant_i_y = float(lines[int(ant_i)].split()[1])
            ant_j_y = float(lines[int(ant_j)].split()[1])
            v = ant_i_y - ant_j_y
            
            visibility_index += 1
            print visibility_index
            #visibility_array[]
            
   print "n_baselines %s" % n_baselines
   
   #for each u,v point calc 1/4PI * INTEGRAL( A(r,nu) * e^(-i2PI B dot r / lambda), dOmega) = X
   # Have system of linear equations V = X Tsky + N
   
   #Maximum likelihood (least squares solution):
   # Tsky = (X^t X)^-1 X^t V
   
   #residuals:
   #covariance matrix:
   
   
   

antenna_locations_filename = '/md0/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant'

plot_different_baselines(antenna_locations_filename)
sys.exit()


plot_only = False
in_orbit = False

#This must be a 'feature' of gsmobserver() - lat and long need to be put in as strings, elevation can be a float (I think)...this is terrible (is my previous gsm_reflection work wrong? rerun?)
mwa_latitude_pyephem = "-26:42.199"
mwa_longitude_pyephem = "116:40.2488"
orbit_latitude_pyephem = "00:00.00"
orbit_longitude_pyephem = "00:00.00"

#horizon?
horizon=True
#ground plane ?
use_ground_plane=True

#if using ground plane - dipole height (normally use 0.3 m for MWA tile dipoles)
dipole_height_m = 0.3

#in m
#mwa_elevation = 377.83
mwa_elevation = 0
orbit_elevation = 35786000

sky_map_orthview_min = -400
sky_map_orthview_max = 2000

freq_MHz = 100.
wavelength = 300./freq_MHz
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
#Use October when Galactic plane sets
year,month,day,hour,minute,second = 2019,10,1,0,0,0

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

#beam is always at zenith, only need to do these once
#beam maps
#isotropic case:
iso_beam_map = np.ones(hp.nside2npix(NSIDE))


#short dipole beam map for parallel case (no ground plane)
#and also with ground plane
if (horizon==False):
   use_ground_plane=False
short_dipole_parallel_beam_map=np.empty_like(iso_beam_map)
for hpx_index,beam_value in enumerate(short_dipole_parallel_beam_map):
   theta,phi=hp.pix2ang(NSIDE,hpx_index)
      
   #populate the dipole model assuming a short dipole (see "New comparison of MWA tile beams" by Benjamin McKinley on Twiki)
   #parallel dipoles, they are placed in an EW baseline config, so the dipoles must be yy NOT XX!)
   #theta_parallel_xx=np.arccos(np.sin(theta)*np.sin(phi)) 
   theta_parallel=np.arccos(np.sin(theta)*np.cos(phi)) 
   if use_ground_plane==True:
      #with ground plane see randall's idl code
      #; calculate the E field response of a 2-element end-fire antenna with
      #; isotropic receivers separated by d, and opposite phase
      #; which is the same as a single element with 
      #; a groundplane and distance d/2 above ground.
      #; separation of antennas: d (wavelengths)
      #; angle of incoming radiation = th (radian, from line joining
      #; antennas)
      #; see Krauss equation 4-10, section 4-2.
      #function end_fire_2element,d,th
      #  return,sin(!pi*d*cos(th))*2.0
      #end
      #d is twice the dipole height in wavelengths
      d_in_lambda = (2. * dipole_height_m)/wavelength
      gp_effect = 2.*math.sin(math.pi*d_in_lambda*math.cos(theta))
      voltage_parallel=np.sin(theta_parallel) * gp_effect
      
      ##randall's code uses this (they are the same):
      #voltage_parallel = sqrt(1-(sin(az[p])*sin(za[p]))^2) * gp_effect
      #sin_theta_parallel = np.sin(theta_parallel)
      #randall_proj_x = np.sqrt(1-(math.sin(phi)*math.sin(theta))**2)
      #proj_y = sqrt(1-(cos(az[p])*sin(za[p]))^2)
      #print 'sin_theta_parallel %s' % sin_theta_parallel
      #print 'randall_proj_x %s' % randall_proj_x
      
   else:
      voltage_parallel=np.sin(theta_parallel)
   
   #TEST for comparing to Randalls, add in 1/cos(za) factor
   #as I thought, this makes the beam wider and therefore its Fourier footprint narrower - not what we want!
   #factor = 1./math.cos(theta)
      
   power_parallel=voltage_parallel**2       #* factor
   short_dipole_parallel_beam_map[hpx_index] = power_parallel
   
   if horizon==True:
      if (theta > np.pi/2.):
         short_dipole_parallel_beam_map[hpx_index]=0.

#normalise
beam_max = np.max(short_dipole_parallel_beam_map)
short_dipole_parallel_beam_map = short_dipole_parallel_beam_map/beam_max
   
#set sky map to zero below the horizon
#(don't actually want to do this to compare to Sing as they put it in space with no horizon. But makes no difference here as this is just for global 
#signal so will be identical above and below horiz (not true for angular structure)))
#Instead, just set the beam to zero below the horizon
#if horizon==True:
#   if (theta > np.pi/2.):
#      sky_map[hpx_index]=0.
#      gsm_map_full_sky_mean_subtr[hpx_index]=0.
      
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
#plt.clf()
#short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
#hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
#iso_beam_map = hp.read_map(beam_map_fitsname)
#iso_beam_map_projected=hp.orthview(map=iso_beam_map,return_projected_map=True,coord='C')


#figmap = plt.gcf()
#figmap.savefig(beam_map_figname,bbox_inches='tight') 
#print "saved figure: %s" % beam_map_figname
#plt.close()

plt.clf()    
#Short dipole
beam_map_fitsname="short_dipole_para_beam_map.fits"
beam_map_figname="short_dipole_para_beam_map.png"
beam_map_title = "short dipole para beam map"
#short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
#hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
#iso_beam_map = hp.read_map(beam_map_fitsname)
short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=beam_map_title,rot=(0,90,0))

figmap = plt.gcf()
figmap.savefig(beam_map_figname,bbox_inches='tight') 
print "saved figure: %s" % beam_map_figname
plt.close()

plt.clf()
#log periodic dipole
beam_map_fitsname="log_periodic_para_beam_map.fits"
beam_map_figname="log_periodic_para_beam_map.png"
beam_map_title = "log periodic para beam map"
#short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
#hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
#iso_beam_map = hp.read_map(beam_map_fitsname)
gaussian_parallel_beam_map_projected=hp.orthview(map=gaussian_parallel_beam_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=beam_map_title,rot=(0,90,0))

figmap = plt.gcf()
figmap.savefig(beam_map_figname,bbox_inches='tight') 
print "saved figure: %s" % beam_map_figname
plt.close()


if in_orbit==True:
   latitude, longitude, elevation = orbit_latitude_pyephem, orbit_longitude_pyephem, orbit_elevation
   #Parkes (pygsm example)
   #latitude, longitude, elevation = '-32.998370', '148.263659', 100
   ov_orbit = GSMObserver()
   ov_orbit.lon = longitude
   ov_orbit.lat = latitude
   ov_orbit.elev = elevation

   
   
# do all this stuff for each hour over 12 hours
for hour_index, hour in enumerate(np.arange(0,24)):
   time_string = "%02d_%02d_%02d" % (hour,minute,second)
   print time_string
   date_obs = datetime(year, month, day, hour, minute, second)
   ov.date = date_obs
   if in_orbit:
      ov_orbit.date = date_obs
      #hack pygsm to do full sky (no masking below horizom)
 
      gsm_map_full_sky=ov_orbit.generate(freq_MHz)
      #Want to just seee the response to angular variations on the sky (subtract mean)
      gsm_map_full_sky_mean = np.mean(gsm_map_full_sky)
      gsm_map_full_sky_mean_subtr = gsm_map_full_sky - gsm_map_full_sky_mean
      sky_map =  np.full(hp.nside2npix(NSIDE),gsm_map_full_sky_mean)
   
   
   else:
      gsm_map_full_sky=ov.generate(freq_MHz)
      
      #Want to just seee the response to angular variations on the sky (subtract mean)
      gsm_map_full_sky_mean = np.mean(gsm_map_full_sky)
      gsm_map_full_sky_mean_subtr = gsm_map_full_sky - gsm_map_full_sky_mean
      
      #global
      #sky map - global - this needs to be the average of the sky map so we can compare with the angular response
      #sky_map =  np.ones(hp.nside2npix(NSIDE))
      sky_map =  np.full(hp.nside2npix(NSIDE),gsm_map_full_sky_mean)
      
   plt.clf()
   map_title="Sky from MWA at h:m:s %02d:%02d:%02d" % (hour,minute,second)
   hp.orthview(map=gsm_map_full_sky,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,90,0),min=0,max=sky_map_orthview_max)
   #ov.view()
   fig_name="sky_from_mwa_at_h_m_s_%s_full_sky.png" % (time_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name,dpi=500,bbox_inches='tight')
   print "saved %s" % fig_name
   plt.close()
   

   
   plt.clf()
   map_title="Sky from MWA at h:m:s %02d:%02d:%02d" % (hour,minute,second)
   hp.orthview(map=gsm_map_full_sky_mean_subtr,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,90,0),min=sky_map_orthview_min,max=sky_map_orthview_max)
   #ov.view()
   fig_name="sky_from_mwa_at_h_m_s_%s_full_sky_mean_subtr.png" % (time_string)
   figmap = plt.gcf()
   figmap.savefig(fig_name,dpi=500,bbox_inches='tight')
   print "saved %s" % fig_name
   plt.close()
   

   
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
   #   #hp.orthview(map=gsm_map,half_sky=False,xsize=2000,title=map_title,coord='E',rot=(0,90,0))
   #   ov.view()
   #   fig_name="sky_from_mwa_at_h_m_s_%02d_%02d_%02d.png" % (hour,minute,second)
   #   figmap = plt.gcf()
   #   figmap.savefig(fig_name,dpi=500,bbox_inches='tight')
   #   print "saved %s" % fig_name
   #   plt.close()

   visibility_amp_array_iso_filename = "visibility_amp_list_iso_%s.npy" % (time_string)
   visibility_real_array_short_para_filename = "visibility_amp_array_short_para_%s.npy" % (time_string)
   visibility_imag_array_short_para_filename = "visibility_imag_array_short_para_%s.npy" % (time_string)
   visibility_real_array_short_para_angular_filename = "visibility_real_array_short_para_angular_%s.npy" % (time_string)
   visibility_imag_array_short_para_angular_filename = "visibility_imag_array_short_para_angular_%s.npy" % (time_string)
      

   #save sky map images and data:
   #global
   fitsname="global_sky_map_h_m_s_%s.fits" % (time_string)
   hp.write_map(fitsname,sky_map,dtype=np.float32, overwrite=True)
   
   plt.clf()
   figname="global_sky_map_h_m_s_%s.png" % (time_string)
   title = "global sky map"
   #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
   #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
   #iso_beam_map = hp.read_map(beam_map_fitsname)
   map_projected=hp.orthview(map=sky_map,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=title,rot=(0,90,0))

   figmap = plt.gcf()
   figmap.savefig(figname,dpi=500,bbox_inches='tight') 
   print "saved figure: %s" % figname
   plt.close()
      
   #show the sky maps multiplied by the beam
   #global
   global_sky_with_beam = sky_map * short_dipole_parallel_beam_map
   
   fitsname="global_sky_with_beam_h_m_s_%s.fits" % (time_string)
   hp.write_map(fitsname,global_sky_with_beam,dtype=np.float32, overwrite=True)
   
   plt.clf()
   figname="global_sky_with_beam_h_m_s_%s.png" % (time_string)
   title = "global sky map"
   #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
   #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
   #iso_beam_map = hp.read_map(beam_map_fitsname)
   map_projected=hp.orthview(map=global_sky_with_beam,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=title,rot=(0,90,0))

   figmap = plt.gcf()
   figmap.savefig(figname,dpi=500,bbox_inches='tight') 
   print "saved figure: %s" % figname
   plt.close()   
   
   
   #angular
   angular_sky_with_beam = gsm_map_full_sky_mean_subtr * short_dipole_parallel_beam_map
   
   fitsname="angular_sky_with_beam_h_m_s_%s.fits" % (time_string)
   hp.write_map(fitsname,angular_sky_with_beam,dtype=np.float32, overwrite=True)
   
   plt.clf()
   figname="angular_sky_with_beam_h_m_s_%s.png" % (time_string)
   title = "angular sky map"
   #short_dipole_parallel_beam_map_projected=hp.orthview(map=short_dipole_parallel_beam_map,coord='C',half_sky=True,xsize=400,title=beam_map_title,rot=(0,90,0),return_projected_map=True)
   #hp.write_map(beam_map_fitsname, iso_beam_map,dtype=np.float32, overwrite=True)
   #iso_beam_map = hp.read_map(beam_map_fitsname)
   map_projected=hp.orthview(map=angular_sky_with_beam,return_projected_map=True,coord='E',half_sky=False,xsize=400,title=title,rot=(0,90,0),min=sky_map_orthview_min,max=sky_map_orthview_max)

   figmap = plt.gcf()
   figmap.savefig(figname,dpi=500,bbox_inches='tight') 
   print "saved figure: %s" % figname
   plt.close()
      
   if not plot_only:   
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

   #abs log
   plt.clf()
   plot=plt.figure()
   plot_title="Visibility amplitude log abs vs baseline length"
   #plot_title="Vis ampl vs freq baseline:%0.1fm and sep %0.1fm" % (baseline_length_m,linear_array_ant_separation)
   plot_figname="visibility_amp_abs_vs_baseline_length_%s_log.png" % (time_string)
   plt.plot(baseline_length_lambda_array,np.log10(visibility_real_array_short_para_norm_abs),label="global real")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_imag_array_short_para_norm_abs),label="global imag")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_real_array_short_para_angular_norm_abs),label="angular real")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_imag_array_short_para_angular_norm_abs),label="angular imag")
   
   plt.title(plot_title)
   plt.ylabel('Log10 visibility amplitude (abs)')
   plt.xlabel('Baseline length (lambda)')
   #plt.xlabel('Frequency (MHz)')
   #plt.xticks(np.arange(min(baseline_length_lambda_array), max(baseline_length_lambda_array)+1, 1.0))
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(-4, 0.1)
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

   #average abs log
   plt.clf()
   plot=plt.figure()
   plot_title="Log visibility amplitude abs vs baseline length"
   plot_figname="visibility_amp_abs_vs_baseline_length_average_to_%s_log.png" % (time_string)
   plt.plot(baseline_length_lambda_array,np.log10(visibility_real_array_short_para_norm_average_abs),label="global real")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_imag_array_short_para_norm_average_abs),label="global real")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_real_array_short_para_angular_norm_average_abs),label="angular real")
   plt.plot(baseline_length_lambda_array,np.log10(visibility_imag_array_short_para_angular_norm_average_abs),label="angular imag")
   
   plt.title(plot_title)
   plt.ylabel('Log10 visibility amplitude (abs)')
   plt.xlabel('Baseline length (lambda)')
   plt.legend(loc=1)
   ax = plt.gca()
   ax.set_ylim(-4, 0.1)
   plot.savefig('%s' % plot_figname)
   print "saved %s" % plot_figname
   plt.close() 
    
   #Add some text to the image(s)
   #img = Image.open("%s_yy_restor.png" % out_image_basename)
   #draw = ImageDraw.Draw(img)
   ## font = ImageFont.truetype(<font-file>, <font-size>)
   ##font = ImageFont.truetype("sans-serif.ttf", 16)
   #font = ImageFont.truetype('FreeSans.ttf',30)
   ## draw.text((x, y),"Sample Text",(r,g,b))
   #draw.text((10, 10),"YY Channel %s (%.0f MHz)" % (chan,freq_MHz),(0,0,0),font=font)
   ##draw.text((256, 256),"Channel %s" % chan,(0,0,0))
   #img.save("%s_yy_restor.png" % out_image_basename)
            
            
   #Tile images together at each timestep
   
   #Tile images together at each timestep

   sky_figname="angular_sky_with_beam_h_m_s_%s.png" % (time_string)
   response_figname="visibility_amp_vs_baseline_length_average_to_%s.png" % (time_string)
   av_response_abs_log_figname="visibility_amp_abs_vs_baseline_length_average_to_%s_log.png" % (time_string)
   bottom_figname = "bottom_%s.png" % (time_string)

   bottom_list_im = [response_figname, av_response_abs_log_figname]

   imgs    = [ Image.open(i) for i in bottom_list_im ]
   ## pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
   min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]

   ##horiz:
   imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
   imgs_comb = Image.fromarray( imgs_comb)
   imgs_comb.save(bottom_figname)
   #
   ##then combine top and bottom vertically
   ## for a vertical stacking it is simple: use vstack
   top_list_im = [sky_figname,bottom_figname]
   imgs    = [ Image.open(i) for i in top_list_im ]
   min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
   imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
   imgs_comb = Image.fromarray( imgs_comb)
   trifecta_name = 'trifecta_%03d.png' % (hour_index)
   imgs_comb.save(trifecta_name)
   print "saved %s" % (trifecta_name)

#stitch together images into movie

filename_base = 'trifecta_%03d.png'
movie_name = 'global_step.mp4'

cmd = "ffmpeg -framerate 6 -i %s -c:v libx264 -r 30 -pix_fmt yuv420p %s" % (filename_base,movie_name)
print cmd
os.system(cmd)
print('made movie %s' % movie_name)

   
   
   

