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

#print baseline_length_lambda_array


#beam map
#isotropic case:
beam_map = np.ones(hp.nside2npix(NSIDE))
visibility_element_array = np.empty_like(beam_map,dtype=complex)

#repeat for each frequency/wavelength:

visibility_amp_list = []
for baseline_length_lambda_index,baseline_length_lambda in enumerate(baseline_length_lambda_array):
   baseline_vector=np.array([baseline_length_m,0,0])
   #print np.linalg.norm(baseline_vector)
   #print baseline_length_lambda

   #for each sky position pixel:
   for hpx_index,beam_value in enumerate(beam_map):
      sky_vector=hp.pix2vec(NSIDE,hpx_index)
      #print np.linalg.norm(sky_vector)
      b_dot_r = np.dot(baseline_vector,sky_vector)
      #print b_dot_r
      #print b_dot_r
      phase_angle = 2.*np.pi*b_dot_r/wavelength_array[baseline_length_lambda_index]
      #print phase_angle
      element = beam_map[hpx_index]*np.exp(-1j*phase_angle)
      #print element
      visibility_element_array[hpx_index] = element
      #print element
      #print visibility_element_array[hpx_index]
   
   #integrate across the whole sky (using healpix so it is just a sum)
   visibility = np.sum(visibility_element_array)
   print visibility
   visibility_real = visibility.real
   print visibility_real
   visibility_amp_list.append(visibility_real)
 
#normalise to value at b=0   
visibility_amp_list = visibility_amp_list/visibility_amp_list[0]
   
plt.clf()
plot=plt.figure()
plot_title="Visibility amplitude vs baseline length for parallel config"
plot_figname="visibility_amp_vs_baseline_length_single.png"
#plt.errorbar(freq_array_band,Tb_band_1,yerr=Tb_error_band_1,label="Measured")
plt.plot(baseline_length_lambda_array,visibility_amp_list,label="single iso element")
plt.title(plot_title)
plt.ylabel('Visibility amplitude')
plt.xlabel('Baseline length (lambda)')
plt.legend(loc=1)
plot.savefig('%s' % plot_figname)
print "saved %s" % plot_figname
plt.close()   


