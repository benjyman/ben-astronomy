#!/usr/bin/env python
#script to calculate sensitivity of interferometers to a global signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cmath

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

baseline_length_m = 10.0

freq_MHz_array = np.arange(0,300)

wavelength_array = 300.0/freq_MHz_array

baseline_length_lambda_array = baseline_length_m/wavelength_array
print baseline_length_lambda_array[34]
#print baseline_length_lambda_array


#beam map
#isotropic case:
iso_beam_map = np.ones(hp.nside2npix(NSIDE))

#make short dipole beam map for parallel case
short_dipole_parallel_beam_map=np.empty_like(iso_beam_map)
for hpx_index,beam_value in enumerate(short_dipole_parallel_beam_map):
   theta,phi=hp.pix2ang(NSIDE,hpx_index)
   #populate the dipole model assuming a short dipole (see "New comparison of MWA tile beams" by Benjamin McKinley on Twiki)
   #parallel 
   theta_parallel=np.arccos(np.sin(theta)*np.sin(phi))
   voltage_parallel=np.sin(theta_parallel)
   power_parallel=voltage_parallel**2
   short_dipole_parallel_beam_map[hpx_index] = power_parallel

visibility_element_array_iso = np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_iso = []
visibility_element_array_short_para = np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_short_para = []
visibility_element_array_short_para_four = np.empty_like(iso_beam_map,dtype=complex)
visibility_amp_list_short_para_four = []


for baseline_length_lambda_index,baseline_length_lambda in enumerate(baseline_length_lambda_array):
   print baseline_length_lambda
   baseline_theta = np.pi/2.
   baseline_phi = 0
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
      element_iso = iso_beam_map[hpx_index]*np.exp(-1j*phase_angle)
      visibility_element_array_iso[hpx_index] = element_iso
      element_short_para = short_dipole_parallel_beam_map[hpx_index]*np.exp(-1j*phase_angle)
      visibility_element_array_short_para[hpx_index] = element_short_para
   
      #now do array:
      theta,phi=hp.pix2ang(NSIDE,hpx_index)
      psi_parallel=np.pi*np.cos(np.arccos(np.sin(theta)*np.sin(phi)))
      #number of elements in array
      N=4.
      array_factor = (1./N)*((np.sin(N*psi_parallel/2.))/(np.sin(psi_parallel/2.)))
      element_short_para_four = short_dipole_parallel_beam_map[hpx_index]*array_factor*np.exp(-1j*phase_angle)
      visibility_element_array_short_para_four[hpx_index] = element_short_para_four      
   
   
   #integrate across the whole sky (using healpix so it is just a sum)
   visibility_iso = np.sum(visibility_element_array_iso)
   visibility_real_iso = visibility_iso.real
   visibility_amp_list_iso.append(visibility_real_iso)

   visibility_short_para = np.sum(visibility_element_array_short_para)
   visibility_real_short_para = visibility_short_para.real
   visibility_amp_list_short_para.append(visibility_real_short_para)

   visibility_short_para_four = np.sum(visibility_element_array_short_para_four)
   visibility_real_short_para_four = visibility_short_para_four.real
   visibility_amp_list_short_para_four.append(visibility_real_short_para_four)
      
#normalise to value at b=0   
visibility_amp_list_iso = visibility_amp_list_iso/visibility_amp_list_iso[0]
visibility_amp_list_short_para = visibility_amp_list_short_para/visibility_amp_list_short_para[0]
visibility_amp_list_short_para_four = visibility_amp_list_short_para_four/visibility_amp_list_short_para_four[0]


#for 21CMA, antenna spacing is 1.66 m, so for this to be lambda/2 the lambda = 3.32, frequency = 90.36 MHz, and the baseline length of 3.7 m = 1.114 lambda
baseline_length_21CMA_90_MHz_lambda = 1.114
vis_amp_21CMA_90_MHz = visibility_amp_list_short_para_four[34]
print "vis_amp_21CMA_90_MHz %s" % vis_amp_21CMA_90_MHz

plt.clf()
plot=plt.figure()
plot_title="Visibility amplitude vs baseline length"
plot_figname="visibility_amp_vs_baseline_length.png"
#plt.errorbar(freq_array_band,Tb_band_1,yerr=Tb_error_band_1,label="Measured")
plt.plot(baseline_length_lambda_array,visibility_amp_list_iso,label="single isotropic")
plt.plot(baseline_length_lambda_array,visibility_amp_list_short_para,label="single short para")
plt.plot(baseline_length_lambda_array,visibility_amp_list_short_para_four,label="four short para")
plt.scatter(baseline_length_21CMA_90_MHz_lambda,vis_amp_21CMA_90_MHz,label="21 CMA short para 90 MHz",marker='+',c='r')
plt.title(plot_title)
plt.ylabel('Visibility amplitude')
plt.xlabel('Baseline length (lambda)')
plt.xticks(np.arange(min(baseline_length_lambda_array), max(baseline_length_lambda_array)+1, 1.0))
plt.legend(loc=1)
plot.savefig('%s' % plot_figname)
print "saved %s" % plot_figname
plt.close()   


