#!/usr/bin/env python
#ASSASSIN: (All Sky SigAl Short SPacing INterferometer)
#Script to replicate Cath's global step simulations

#Daniel:
# The format of the data is in standard spherical coordinates (theta, phi).
# Where theta is the zenith angle and phi is anti-clockwise going from east to north looking down at the array. 
# The 'origin' will be on the top left corner. Going from top to bottom is increasing in phi and going left to right is increasing in theta. 
# The data is in 1 deg steps.

import numpy as np
import healpy as hp
from pygsm import GSMObserver2016
from datetime import datetime, date
import matplotlib.pyplot as plt
import os,sys
from reproject import reproject_from_healpix
import pyfits
from astropy.wcs import WCS
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from pyuvdata import UVBeam
import math


#mu_0 = 4.*np.pi*10**(-7)
k = 1.38065e-23
sq_deg_in_1_sr = (180./math.pi)**2
#Need to specify the location of the observatory (MWA ... put in exact EDA later)
# Setup observatory location - in this case, MWA, Australia
latitude_degrees=-26.70331940
longitude_degrees=116.67081524
elevation_m=377.83


use_analytic_beam = True
generate_new_hpx = True
plot_gsm_map_hpx = True
generate_new_vis = True

#This is for Daniels FEKO model beams
generate_new_average_beam = False
#apply_normalisation = True
plot_all_beams = False
plot_average_beam = False

#specify date
#date_time_string = '2018_09_25_00_00_00'
#lst_list = ['00', '03', '06', '09', '12', '15', '18', '21']
lst_list = ['12']
#pol_list = ['X','Y']
pol_list = ['X']
#freq_MHz_list = np.range(50,200,1)
freq_MHz_list = [160.]

n_ants = 256
template_imsize = 512
template_cell_size_asec = 180./template_imsize*60.*60.

def write_beam_fits_sin(cart_image,fitsname,lst=None):
   if lst==None:
      lst=60
   else:
      lst=lst
   hdu = fits.PrimaryHDU(cart_image)
   hdul = fits.HDUList([hdu])
   hdul.writeto(fitsname,clobber=True)
   header = hdul[0].header
   header['CTYPE1']  = 'RA---SIN'
   header['CRPIX1']  =  91
   header['CRVAL1']  = 60
   header['CDELT1']  = -0.636620
   header['CUNIT1']  = 'deg     '
   header['CTYPE2']  = 'DEC--SIN'
   header['CRPIX2']  = 91
   header['CRVAL2']  = -27
   header['CDELT2']  = 0.636620
   header['CUNIT2']  = 'deg  '
   #target_header = fits.Header.fromstring("NAXIS   =                    2\nNAXIS1  =                  180\nNAXIS2  =                  180\nCTYPE1  = 'RA---SIN'\nCRPIX1  =                   90\nCRVAL1  =                   -4\nCDELT1  =                    1\nCUNIT1  =' deg     '\nCTYPE2  = 'DEC--SIN'\nCRPIX2  =                   90\nCRVAL2  =                  -27\nCDELT2  =                    1\nCUNIT2  ='deg     '\nCOORDSYS= 'icrs    '\n", sep='\n')
   #target_header = fits.Header
   #pyfits.writeto(Power_pattern_average_interp_fitsname,cart_image,clobber=True)
   fits.update(fitsname,cart_image,header=header)
   print "wrote image %s" %  fitsname
   
   
def polar2cartesian(r, t, grid, x, y, order=3):

    X, Y = np.meshgrid(x, y)

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(X, Y)

    ir = interp1d(r, np.arange(len(r)), bounds_error=False)
    it = interp1d(t, np.arange(len(t)))

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    return map_coordinates(grid, np.array([new_ir, new_it]),
                            order=order).reshape(new_r.shape)
                            
#Do this stuff for each lst and each freq:                            
for lst in lst_list:
   print(lst)
   for freq_MHz in freq_MHz_list:
      print(freq_MHz)
      #for date just use 00JAN1$fakedayfrac as in Randall's sims (can't put in a four digit year anyway)
      date_time_string = '1900_01_01_%02d_00_00' % (float(lst))
      print date_time_string
      wavelength = 300./freq_MHz
      freq_GHz = freq_MHz/1000.
      
      gsm_hpx_fits_name = "gsm_map_LST_%s_%s_MHz_hpx.fits" % (lst,int(freq_MHz))
      reprojected_gsm_fitsname = "gsm_map_LST_%s_%s_MHz_reprojected.fits" % (lst,int(freq_MHz))
      reprojected_gsm_im_name = "gsm_map_LST_%s_%s_MHz_reprojected.im" % (lst,int(freq_MHz))
      #lna_impedance_aavs1_filename = "/data/code/git/ben-astronomy/AAVS-1/AAVS1_LNA_impedance_180718.txt"
      
      eda_model_vis_name = "eda_model_LST_%s_%s_MHz.vis" % (lst,int(freq_MHz))
      eda_model_no_source_image_name = "eda_no_source_LST_%s_%s_MHz.map" % (lst,int(freq_MHz))
      eda_model_no_source_beam_name = "eda_no_source_LST_%s_%s_MHz.beam" % (lst,int(freq_MHz))
      eda_model_no_source_image_fits_name = "eda_no_source_LST_%s_%s_MHz.fits" % (lst,int(freq_MHz))
      
      
      
      #if apply_normalisation:
      #   #get the lna impedance for this frequenccy
      #   with open(lna_impedance_aavs1_filename) as f:
      #      lines = f.readlines()
      #   for line in lines:
      #      freq = line.split()[0]
      #      if float(freq) == float(freq_MHz):
      #         lna_z_real = float(line.split()[1])
      #         lna_z_imag = float(line.split()[2])
      #
      #   #print "%s %s" % (lna_z_real,lna_z_imag)
      #   lna_z_complex=complex(lna_z_real,lna_z_imag)
      #   print "LNA impedance is %s" % lna_z_complex
      #
      #   #w = 2pi*freq
      #   w = 2.*np.pi*float(freq_MHz)*10**(6)
      #
      #   #N_factor = (-j4pi/(w u0))*Z_lna
      #   normalisation_factor = (4.*np.pi*1j / (mu_0 * w)) * lna_z_complex
      
      #Generate model visibilities
      #1. Get a sky model (healpix map of GSM2017) - easy
      #2. Get a beam model in healpix (not entirely straightforward, but I think I have done this before) just use short dipole model.
      #3. Generate a point source model (I don't understand how this is done)
      #4. Generate a global EoR signal (I think I understand how to do this from Caths eqn 3 etc
       
      #Get EDA antenna positions (you have these)
      
      #Use zenith pointing centre and phase centre of RA = 4h, Dec= -27 degrees
      
      #Generate simulated visibilities ... with miriad to start?
      #uvgen to make a blank visibility set (no sources, but with correct observing parameters)
      #uvmodel with the 'add' option to take the pb attenuated gsm image and add it to the visibility file from uvgen
      
      #2. Sky model
      #Need a pb-attenuated gsm image in slant-orthographic projection
      # - obtain gsm image at desired frequency and datetime (healpix)

      ov = GSMObserver2016()
      ov.lon = longitude_degrees
      ov.lat = latitude_degrees
      ov.elev = elevation_m
      
      #and datetime of observation (eventually do this for many dates)
      year=int(date_time_string.split("_")[0])
      month=int(date_time_string.split("_")[1])
      day=int(date_time_string.split("_")[2])
      hour=int(date_time_string.split("_")[3])
      minute=int(date_time_string.split("_")[4])
      second=int(date_time_string.split("_")[5])
      
      ov.date = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
      
      
      #1. get 'blank' visibilities in miriad and make an image (1 Jy pt source) to use as a template for the output projection
      if generate_new_vis:
      
         cmd = "rm -rf eda_model*"
         print(cmd)
         os.system(cmd)
      
         cmd = "uvgen source=$MIRCAT/no.source ant='/data/code/git/ben-astronomy/AAVS-1/AAVS1_loc_uvgen.ant' baseunit=-3.33564 corr='1,1,0,0.781' time=19SEP26:00:00:00.0 freq=%.2f,0.0 radec='4,-27' harange=0,0.0005556,0.0005556 lat=-26.70331940 out=%s stokes=xx,yy,xy,yx" % (freq_GHz,eda_model_vis_name)
         print(cmd)
         os.system(cmd)
         
         cmd = "rm -rf %s %s %s " % (eda_model_no_source_image_name,eda_model_no_source_beam_name,eda_model_no_source_image_fits_name)
         print(cmd)
         os.system(cmd)
         
         cmd = "invert vis=%s map=%s beam=%s imsize=%s cell=%s stokes=i options=mfs" % (eda_model_vis_name,eda_model_no_source_image_name,eda_model_no_source_beam_name,template_imsize,template_cell_size_asec)
         print(cmd)
         os.system(cmd)
      
         cmd = "fits in=%s out=%s op=xyout" % (eda_model_no_source_image_name,eda_model_no_source_image_fits_name)
         print(cmd)
         os.system(cmd)
         
      
             
      if generate_new_hpx:
         cmd = "rm -rf %s %s" % (gsm_hpx_fits_name,reprojected_gsm_im_name)
         print(cmd)
         os.system(cmd)
         gsm_map = ov.generate(freq_MHz)
         hp.write_map(gsm_hpx_fits_name,gsm_map,coord='C')
      else:
         gsm_map = hp.read_map(gsm_hpx_fits_name)
      
      #plot?
      if plot_gsm_map_hpx:
         #plot
         plt.clf()
         map_title="GSM from MWA at %.0f:%.0f:%.0f" % (hour,minute,second)
         #hp.orthview(map=gsm_map,half_sky=True,xsize=2000,title=map_title,rot=(0,0,0),min=0, max=100)
         hp.orthview(map=gsm_map,half_sky=False,title=map_title)
         #hp.mollview(map=gsm_map_from_moon,coord='C',xsize=400,title=map_title)
         fig_name="gsm_map_%s_%sMHz.png" % (date_time_string,int(freq_MHz))
         figmap = plt.gcf()
         figmap.savefig(fig_name,dpi=500)
         print "saved %s" % fig_name
      
      #Miriad doesn't seem to be able to import the hpx file
      #Try using reproject_from_healpix
      #output_projection is the header of the pt source made above
      if os.path.isfile(eda_model_no_source_image_fits_name) and os.access(eda_model_no_source_image_fits_name, os.R_OK):
         hdulist = pyfits.open(eda_model_no_source_image_fits_name)
      else:
         print "Either file %s is missing or is not readable" % eda_model_no_source_image_fits_name
         #continue        
      
      data=hdulist[0].data[0,0,:,:]
      no_source_header=hdulist[0].header
      pix_size_deg = float(no_source_header['CDELT1'])
      pix_area_deg_sq = pix_size_deg*pix_size_deg
      pix_area_sr = pix_area_deg_sq / sq_deg_in_1_sr
      
      
      
      #needs cooordsys keyword
      #pt_source_header['COORDSYS'] = 'icrs'
      
      #print(pt_source_header)
      
      #del pt_source_header[8]
      #del pt_source_header[8]
      del no_source_header['history']
                      
      #print(pt_source_header)
      
      target_wcs = WCS(no_source_header)
      
      target_wcs=target_wcs.dropaxis(2)
      target_wcs=target_wcs.dropaxis(2)
                      
      #hdu_hpx = pyfits.open(gsm_hpx_fits_name)[1]
      ##hdu_hpx.info()
      #hpx_header = hdu_hpx.header
      #print(hpx_header)
      
      reprojected_gsm_map,footprint = reproject_from_healpix(gsm_hpx_fits_name, target_wcs,shape_out=(template_imsize,template_imsize), order='bilinear')
      
      #write the map to fits
      pyfits.writeto(reprojected_gsm_fitsname,reprojected_gsm_map,clobber=True)
      pyfits.update(reprojected_gsm_fitsname,reprojected_gsm_map,header=no_source_header)
      print "wrote image %s" %  reprojected_gsm_fitsname
            
      #Maybe don't need this at all:          
      # - use miriad regrid to regrid to slant orthographic? (or 'unhealpix' https://github.com/AlecThomson/unhealpix)
      #cmd = "regrid in=%s" % (gsm_hpx_fits_name)
      #print(cmd)
      #os.system(cmd)
      
      
      
      #Do GSM map stuff here as it doesn't depend on pol
      cmd = "fits in=%s out=%s op=xyin" % (reprojected_gsm_fitsname,reprojected_gsm_im_name)
      print(cmd)
      os.system(cmd)
      
      #uvmodel requires the model to be in Jy/pix
      #This scaling doesn't take into account the changing pixel area across the image - need too account for this somewhere with a 1/cos(za) term (can do it in the beam...)
      
      scale = (2. * k * 1.0e26 * pix_area_sr) / (wavelength**2)
      print "scale map by %s to get to Jy/pix" % scale
      
      reprojected_gsm_im_Jy_per_pix_name =  "gsm_map_%s_%s_MHz_reprojected_Jy_pix.im" % (date_time_string,int(freq_MHz))
      
      cmd = "rm -rf %s" % reprojected_gsm_im_Jy_per_pix_name
      
      cmd = "maths exp=%s*%s out=%s " % (scale,reprojected_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name)
      print(cmd)
      os.system(cmd)
      
      
      #Now for the beam (EDA)!
      if generate_new_average_beam:
         if freq_MHz==160.:
            EDA_chan=125
         else:
            print("no EDA beam data for this frequency (%s MHz)" % freq_MHz)
         for freq_MHz in freq_MHz_list:
            for pol in pol_list:   
               Etheta_mag_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Etheta_mag_ch%s.fits" % (EDA_chan,pol,EDA_chan)
               Etheta_phase_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Etheta_phase_ch%s.fits" % (EDA_chan,pol,EDA_chan)
               Ephi_mag_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Ephi_mag_ch%s.fits" % (EDA_chan,pol,EDA_chan)
               Ephi_phase_filename = "/md0/AAVS-1/beam_pattern_tests/new_data/EDA/Ch_%s/EDA_%spol_Ephi_phase_ch%s.fits" % (EDA_chan,pol,EDA_chan)
         
               Ephi_mag_HDUlist = pyfits.open(Ephi_mag_filename)
               Ephi_mag_header = Ephi_mag_HDUlist[0].header
            
               Ephi_mag_HDUlist.verify('fix')
               Ephi_mag_data = Ephi_mag_HDUlist[0].data
               
               Ephi_phase_HDUlist = pyfits.open(Ephi_phase_filename)
               Ephi_phase_header = Ephi_phase_HDUlist[0].header
               Ephi_phase_HDUlist.verify('fix')
               Ephi_phase_data = Ephi_phase_HDUlist[0].data   
            
               Etheta_mag_HDUlist = pyfits.open(Etheta_mag_filename)
               Etheta_mag_header = Etheta_mag_HDUlist[0].header
               Etheta_mag_HDUlist.verify('fix')
               Etheta_mag_data = Etheta_mag_HDUlist[0].data
             
               Etheta_phase_HDUlist = pyfits.open(Etheta_phase_filename)
               Etheta_phase_header = Etheta_phase_HDUlist[0].header
               Etheta_phase_HDUlist.verify('fix')
               Etheta_phase_data = Etheta_phase_HDUlist[0].data  
               
               #Want an average beam pattern
               average_beam_array = np.zeros([361,91])
               
               for ant_index in range(0,n_ants-1):
                  #ant_name_randall = line.split()[0][3:6]
                  ant_name_randall = "%03d" % (ant_index+1)
                  #print("Ant %s" % ant_name_randall)
                 
                  Power_pattern_ant_interp_title = 'Antenna %s Power Pattern %s Pol %s MHz' % (ant_name_randall,pol,int(freq_MHz))
                  Power_pattern_ant_interp_figname = 'Power_pattern_Ant%s_%s_%s_MHz_interp.png' % (ant_name_randall,pol,int(freq_MHz))
                  power_pattern_ant_filename = '/mnt/md0/AAVS-1/beam_pattern_tests/new_data/%s_MHz/AAVS1_power_pattern_linear_Ant%s_%s_theta_phi.npy' % (int(freq_MHz),ant_name_randall,pol)
                  
                  Ephi_phase_data_ant = Ephi_phase_data[ant_index]
                  Ephi_mag_data_ant = Ephi_mag_data[ant_index]
                  Etheta_phase_data_ant = Etheta_phase_data[ant_index]
                  Etheta_mag_data_ant = Etheta_mag_data[ant_index]
                  
                  power_pattern_ant = Etheta_mag_data_ant**2 + Ephi_mag_data_ant**2
                  
                  #print(power_pattern_ant.shape)
                  #print(power_pattern_ant.max())
                  #need to normalise (can do at end?)
                  #include a 1/cos(za) term to account for pixel area changing in sine projection
                  #need to check I have the order of ZA right...
                  for za_deg in range(0,91):
                     za_rad = (za_deg/180.)*np.pi
                     #print("zenith angle is %s deg, divide by factor cos(za)=%s" % (za_deg,np.cos(za_rad)))
                     power_pattern_ant[:,za_deg] = power_pattern_ant[:,za_deg]/np.cos(za_rad)
                  
                  average_beam_array += power_pattern_ant
                  
                  #Ephi_real_ant = Ephi_mag_data_ant * np.cos(Ephi_phase_data_ant)
                  #Ephi_imag_ant = Ephi_mag_data_ant * np.sin(Ephi_phase_data_ant)
                  #Ephi_complex_ant = Ephi_real_ant+Ephi_imag_ant*1j
                  #Etheta_real_ant = Etheta_mag_data_ant * np.cos(Etheta_phase_data_ant)
                  #Etheta_imag_ant = Etheta_mag_data_ant * np.sin(Etheta_phase_data_ant)
                  #Etheta_complex_ant = Etheta_real_ant+Etheta_imag_ant*1j
               
                  #if apply_normalisation:
                  #            
                  #   Ephi_complex_ant = Ephi_complex_ant * normalisation_factor
                  #   Etheta_complex_ant = Etheta_complex_ant * normalisation_factor
                  
                  
                  if plot_all_beams:
                     #Cartesian grid:
                     power_patter_ant_log = 10*np.log10(power_pattern_ant)
                     
                     # Define original polar grid
                     
                     nr = 91
                     nt = 361
                     r = np.linspace(1, 91, nr)
                     t = np.linspace(-np.pi, np.pi, nt)
                     z = np.swapaxes(power_patter_ant_log,0,1)
                     
                     nx = 180
                     ny = 180
                     x = np.linspace(-90, 90., nx)
                     y = np.linspace(-90., 90., ny)
                     
                     # Interpolate polar grid to cartesian grid (nearest neighbor)
                     
                     #blank out below horizon:
                     centre_x = nx/2
                     centre_y = ny/2
                     #centre_x = 3
                     #centre_y = 3
                     y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
                     mask = x**2+y**2 <= centre_x**2
                     mask = 1*mask.astype(float)
                     
                     
                     fig = plt.figure()
                     ax = fig.add_subplot(111)
                     cart_image = polar2cartesian(r, t, z, x, y, order=3)
                     cart_image = cart_image*mask
                     cart_image[cart_image==-0.] = np.nan
                     ###get east and west right
                     ##cart_image = np.flip(cart_image,axis=0)
                     ###sort out N-S (same as rotating 180
                     ##cart_image = np.flip(cart_image,axis=1)
                     cart_image = np.rot90(cart_image, k=1, axes=(0,1))
                     
                     #max_ant001_db = 10*np.log10(max_ant001)
                     #min_ant001_db = 10*np.log10(min_ant001)
                     #
         
                     img = ax.imshow(cart_image, interpolation='nearest', vmax=0,vmin=-25)
                     #fig.savefig('test4.png')
                  
                     cbar = plt.colorbar(img, ax=ax)
                     cbar.set_label('Power', rotation=270, labelpad=10)
                     #plt.colorbar(img, ax=ax)
                     plt.title(Power_pattern_ant_interp_title)
                     plt.savefig(Power_pattern_ant_interp_figname)
                     plt.close()
                     print("saved %s" % Power_pattern_ant_interp_figname)
            
               #calculate average
               average_beam_array = average_beam_array/n_ants
               #normalise the average beam
               # max of average beam array:
               #print average_beam_array
               average_array_max = average_beam_array.max()
               #print("%s" % average_array_max)
               average_beam_array_norm = average_beam_array/average_beam_array.max()
            
               #sin projection sampling 
               #Need to do a better job than this simple nearest neighbour interpolation - but fine for now.
               #print("av beam array shape is")
               #print average_beam_array.shape
               sin_projected_beam_array=average_beam_array*0
               step=1./91.
               sin_za_array=np.arange(0,1,step)
               #print sin_za_array
               sin_proj_sampling_rad = np.arcsin(sin_za_array)
               sin_proj_sampling_deg = sin_proj_sampling_rad/np.pi*180.
               sin_proj_sampling_deg_round = np.around(sin_proj_sampling_deg)
               sin_proj_sampling_deg_indexes = sin_proj_sampling_deg_round.astype(int)
               #print sin_proj_sampling_deg_indexes.shape
               
               for az in range(0,361):
                  power_vs_za_array_proj = sin_proj_sampling_deg_indexes*0.
                  for orig_za_index in range(0,91):
                     new_za_index = sin_proj_sampling_deg_indexes[orig_za_index]
                     power_vs_za_array_proj[orig_za_index] = average_beam_array[az,new_za_index]
                  sin_projected_beam_array[az,:] = power_vs_za_array_proj
                  
               sin_projected_beam_array_norm = sin_projected_beam_array/sin_projected_beam_array.max()
                  
               #   sin_za_array = np.sin(za_array_rad)
               #   print za_array_rad
               #   print sin_za_array
            
               if plot_average_beam:
                  power_pattern_average_interp_title = 'Average Antenna Power Pattern %s Pol %s MHz' % (pol,int(freq_MHz))
                  power_pattern_average_interp_figname = 'power_pattern_average_%s_%s_MHz_interp.png' % (pol,int(freq_MHz))
                  power_pattern_average_interp_fitsname =  'power_pattern_average_%s_%s_MHz_interp.fits' % (pol,int(freq_MHz))
                  
                  
                  #Cartesian grid:
                  #power_pattern_ant_log = 10*np.log10(average_beam_array_norm)
                  
                  # Define original polar grid
                  
                  nr = 91
                  nt = 361
                  r = np.linspace(1, 91, nr)
                  t = np.linspace(-np.pi, np.pi, nt)
                  z = np.swapaxes(average_beam_array_norm,0,1)
                  
                  nx = 180
                  ny = 180
                  x = np.linspace(-90, 90., nx)
                  y = np.linspace(-90., 90., ny)
                  
                  # Interpolate polar grid to cartesian grid (nearest neighbor)
                  
                  #blank out below horizon:
                  centre_x = nx/2
                  centre_y = ny/2
                  #centre_x = 3
                  #centre_y = 3
                  y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
                  mask = x**2+y**2 <= centre_x**2
                  mask = 1*mask.astype(float)
                  
                  
                  fig = plt.figure()
                  ax = fig.add_subplot(111)
                  cart_image = polar2cartesian(r, t, z, x, y, order=3)
                  cart_image = cart_image*mask
                  cart_image[cart_image==-0.] = np.nan
                  ###get east and west right
                  ##cart_image = np.flip(cart_image,axis=0)
                  ###sort out N-S (same as rotating 180
                  ##cart_image = np.flip(cart_image,axis=1)
                  cart_image = np.rot90(cart_image, k=1, axes=(0,1))
                  
                  #max_ant001_db = 10*np.log10(max_ant001)
                  #min_ant001_db = 10*np.log10(min_ant001)
                  #
              
                  img = ax.imshow(cart_image, interpolation='nearest')
                  #fig.savefig('test4.png')
                  
                  cbar = plt.colorbar(img, ax=ax)
                  cbar.set_label('Power', rotation=270, labelpad=10)
                  #plt.colorbar(img, ax=ax)
                  plt.title(power_pattern_average_interp_title)
                  plt.savefig(power_pattern_average_interp_figname)
                  plt.close()
                  print("saved %s" % power_pattern_average_interp_figname)
                  
                  pyfits.writeto(power_pattern_average_interp_fitsname,cart_image,clobber=True)
                  print "wrote fits image %s" %  power_pattern_average_interp_fitsname
         
         
                  ##SIN projected!
                  power_pattern_average_interp_sin_title = 'Average Antenna Power Pattern SIN Proj %s Pol %s MHz' % (pol,int(freq_MHz))
                  power_pattern_average_interp_sin_figname = 'power_pattern_average_%s_%s_MHz_interp_sin.png' % (pol,int(freq_MHz))
                  power_pattern_average_interp_sin_fitsname =  'power_pattern_average_%s_%s_MHz_interp_sin.fits' % (pol,int(freq_MHz))
                  
                  
                  #Cartesian grid:
                  #power_pattern_ant_log = 10*np.log10(average_beam_array_norm)
                  
                  # Define original polar grid
                  
                  nr = 91
                  nt = 361
                  r = np.linspace(1, 91, nr)
                  t = np.linspace(-np.pi, np.pi, nt)
                  z = np.swapaxes(sin_projected_beam_array_norm,0,1)
                  
                  nx = 180
                  ny = 180
                  x = np.linspace(-90, 90., nx)
                  y = np.linspace(-90., 90., ny)
                  
                  # Interpolate polar grid to cartesian grid (nearest neighbor)
                  
                  #blank out below horizon:
                  centre_x = nx/2
                  centre_y = ny/2
                  #centre_x = 3
                  #centre_y = 3
                  y,x = np.ogrid[-centre_x:centre_x, -centre_y: centre_y]
                  mask = x**2+y**2 <= centre_x**2
                  mask = 1*mask.astype(float)
                  
                  
                  fig = plt.figure()
                  ax = fig.add_subplot(111)
                  cart_image_sin = polar2cartesian(r, t, z, x, y, order=3)
                  cart_image_sin = cart_image_sin*mask
                  cart_image_sin[cart_image_sin==-0.] = np.nan
                  ###get east and west right
                  ##cart_image = np.flip(cart_image,axis=0)
                  ###sort out N-S (same as rotating 180
                  ##cart_image = np.flip(cart_image,axis=1)
                  cart_image_sin = np.rot90(cart_image_sin, k=1, axes=(0,1))
                  
                  #max_ant001_db = 10*np.log10(max_ant001)
                  #min_ant001_db = 10*np.log10(min_ant001)
                  #
              
                  img = ax.imshow(cart_image_sin, interpolation='nearest')
                  #fig.savefig('test4.png')
                  
                  cbar = plt.colorbar(img, ax=ax)
                  cbar.set_label('Power', rotation=270, labelpad=10)
                  #plt.colorbar(img, ax=ax)
                  plt.title(power_pattern_average_interp_sin_title)
                  plt.savefig(power_pattern_average_interp_sin_figname)
                  plt.close()
                  print("saved %s" % power_pattern_average_interp_sin_figname)
                  
              
               ##need to write as a fits file with the correct header
         
               write_beam_fits_sin(cart_image_sin,power_pattern_average_interp_sin_fitsname)
               
               #Okay, I think this looks alright. Have a SIN projected beam and a SIN projected gsm
            
       
      #This stuff outside of beam generation loop
      #do for all freqs and all pols:
      for freq_MHz in freq_MHz_list:
         for pol in pol_list:
            #in miriad
            power_pattern_average_interp_sin_im_name = 'power_pattern_average_%s_%s_MHz_interp_sin.im' % (pol,int(freq_MHz))
            power_pattern_average_interp_sin_regrid_gsm_im_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.im' % (pol,int(freq_MHz))
            power_pattern_average_interp_sin_regrid_gsm_fits_name =  'power_pattern_average_%s_%s_MHz_sin_regrid.fits' % (pol,int(freq_MHz))
         
            cmd = "rm -rf %s %s " % (power_pattern_average_interp_sin_im_name,power_pattern_average_interp_sin_regrid_gsm_im_name)
            print(cmd)
            os.system(cmd)
         
            cmd = "fits in=%s out=%s op=xyin" % (power_pattern_average_interp_sin_fitsname,power_pattern_average_interp_sin_im_name)
            print(cmd)
            os.system(cmd)
         
            #regrid beam to gsm (or should I do other way round? beam has already been regridded twice!?)
            cmd = "regrid in=%s out=%s tin=%s tol=0" % (power_pattern_average_interp_sin_im_name,power_pattern_average_interp_sin_regrid_gsm_im_name,reprojected_gsm_im_name)
            print(cmd)
            os.system(cmd)   
            
         
            #Sweet! finally have a gsm and a beam.
            #now multiply them to get the apparent sky and put that into uvmodel, 
         
            apparent_sky_im_name = "apparent_sky_%s_%s_pol_%s_MHz.im" % (date_time_string,pol,int(freq_MHz))
         
            cmd = "rm -rf %s " % apparent_sky_im_name
            print(cmd)
            os.system(cmd)
         
            cmd = "maths exp=%s*%s out=%s " % (power_pattern_average_interp_sin_regrid_gsm_im_name,reprojected_gsm_im_Jy_per_pix_name,apparent_sky_im_name)
            print(cmd)
            os.system(cmd)
         
         
            #then put into the vis 
            model_sky_vis = "eda_model_plus_sky.vis"
         
            cmd = "rm -rf %s" % model_sky_vis
            print(cmd)
            os.system(cmd)
         
            cmd = "uvmodel vis=%s model=%s options=add out=%s" % (eda_model_vis_name,apparent_sky_im_name,model_sky_vis)
            print(cmd)
            os.system(cmd)
            #then can finally try to repeat Caths stuff (but need multiple freqs!)
         
            #Put in a global signal!
         
            #do for 50 - 200 MHz at 1 MHz res
   
      
      
      
   




