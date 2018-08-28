#!/usr/bin/env python
#script to calculate sensitivity of interferometers to a global signal
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cmath
from astropy.wcs.docstrings import phi0

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
NSIDE=32

baseline_length_m = 1.66
linear_array_ant_separation = 1.66

freq_MHz_array = np.arange(0,300)

wavelength_array = 300.0/freq_MHz_array

baseline_length_lambda_array = baseline_length_m/wavelength_array
print baseline_length_lambda_array[34]
#print baseline_length_lambda_array

#sky map
sky_map =  np.ones(hp.nside2npix(NSIDE))

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
   if (theta > np.pi/2.):
      sky_map[hpx_index]=0.
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
visibility_amp_list_short_para = []
visibility_element_array_log= np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_log = []
visibility_element_array_short_para_four = np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_short_para_four = []
visibility_element_array_log_four = np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_log_four = []

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
      #now do array:
      theta,phi=hp.pix2ang(NSIDE,hpx_index)
      #psi is the angle between the sky and the long axis of the 1-D array
      #what happens if you don't fix the antennas at lambda/D spacing ... ie new array factor for each frequency:
      k=2.*np.pi/wavelength_array[baseline_length_lambda_index]
      #print k
      psi_parallel=k*linear_array_ant_separation*np.cos(np.arccos(np.sin(theta)*np.sin(phi)))
      #psi_parallel=np.pi*np.cos(np.arccos(np.sin(theta)*np.sin(phi)))
      #number of elements in array
      N=4.
      array_factor = (1./N)*((np.sin(N*psi_parallel/2.))/(np.sin(psi_parallel/2.)))
      element_short_para_four = sky_map[hpx_index]*short_dipole_parallel_beam_map[hpx_index]*array_factor*np.exp(-1j*phase_angle)
      visibility_element_array_short_para_four[hpx_index] = element_short_para_four      

   
   #integrate across the whole sky (using healpix so it is just a sum)
   visibility_iso = np.sum(visibility_element_array_iso)
   visibility_real_iso = visibility_iso.real
   visibility_amp_list_iso.append(visibility_real_iso)

   visibility_log = np.sum(visibility_element_array_log)
   visibility_real_log = visibility_log.real
   visibility_amp_list_log.append(visibility_real_log)
      
   visibility_short_para = np.sum(visibility_element_array_short_para)
   visibility_real_short_para = visibility_short_para.real
   visibility_amp_list_short_para.append(visibility_real_short_para)

   visibility_short_para_four = np.sum(visibility_element_array_short_para_four)
   visibility_real_short_para_four = visibility_short_para_four.real
   visibility_amp_list_short_para_four.append(visibility_real_short_para_four)
      
#normalise to value at b=0   
visibility_amp_list_iso = visibility_amp_list_iso/visibility_amp_list_iso[0]
visibility_amp_list_short_para = visibility_amp_list_short_para/visibility_amp_list_short_para[0]
visibility_amp_list_short_para_four = visibility_amp_list_short_para_four/np.nanmax(visibility_amp_list_short_para_four)
visibility_amp_list_log = visibility_amp_list_log/np.nanmax(visibility_amp_list_log)

#for 21CMA, antenna spacing is 1.66 m, so for this to be lambda/2 the lambda = 3.32, frequency = 90.36 MHz, and the baseline length of 3.7 m = 1.114 lambda
#baseline_length_21CMA_90_MHz_lambda = 1.114
#vis_amp_21CMA_90_MHz = visibility_amp_list_short_para_four[34]
#print "vis_amp_21CMA_90_MHz %s" % vis_amp_21CMA_90_MHz

plt.clf()
plot=plt.figure()
#plot_title="Visibility amplitude vs baseline length"
plot_title="Vis ampl vs freq baseline:%0.1fm and sep %0.1fm" % (baseline_length_m,linear_array_ant_separation)
#plot_figname="visibility_amp_vs_baseline_length.png"
plot_figname="visibility_amp_vs_frequency_baseline_%0.1fm_ant_sep_%0.1fm.png" % (baseline_length_m,linear_array_ant_separation)
#plt.errorbar(freq_array_band,Tb_band_1,yerr=Tb_error_band_1,label="Measured")
#plt.plot(baseline_length_lambda_array,visibility_amp_list_iso,label="single isotropic")
#plt.plot(baseline_length_lambda_array,visibility_amp_list_short_para,label="single short para")
#plt.plot(baseline_length_lambda_array,visibility_amp_list_short_para_four,label="four short para")
plt.plot(freq_MHz_array,visibility_amp_list_iso,label="single isotropic")
plt.plot(freq_MHz_array,visibility_amp_list_log,label="single log-periodic")
plt.plot(freq_MHz_array,visibility_amp_list_short_para,label="single short para")
plt.plot(freq_MHz_array,visibility_amp_list_short_para_four,label="four short para")

#plt.scatter(baseline_length_21CMA_90_MHz_lambda,vis_amp_21CMA_90_MHz,label="21 CMA short para 90 MHz",marker='+',c='r')
plt.title(plot_title)
plt.ylabel('Visibility amplitude')
#plt.xlabel('Baseline length (lambda)')
plt.xlabel('Frequency (MHz)')
#plt.xticks(np.arange(min(baseline_length_lambda_array), max(baseline_length_lambda_array)+1, 1.0))
plt.legend(loc=1)
plot.savefig('%s' % plot_figname)
print "saved %s" % plot_figname
plt.close()   



