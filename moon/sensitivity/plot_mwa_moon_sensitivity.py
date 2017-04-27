# This script plots the lunar-flux-sensitivty  for the MWA 128T array
# Check the observational parameters below
# Note that the sqrt(2) increase in noise due to inter-night differencing has not been taken into account

import numpy as np
import scipy.special
import matplotlib.pyplot as plt

def make_baselines(pos):
   pos=np.array(pos)
   u = np.meshgrid(pos[0,]) - np.transpose(np.meshgrid(pos[0,])) # Build u,v,w values
   v = np.meshgrid(pos[1,]) - np.transpose(np.meshgrid(pos[1,]))
   w = np.meshgrid(pos[2,]) - np.transpose(np.meshgrid(pos[2,]))
   return [u,v,w]



def read_mwa_tile_pos():
   fname = 'ant_MWA_ph2_extended_miriad.txt'
   #fname = 'ant128_miriad.txt'
   array = np.loadtxt(fname, delimiter=' ', usecols=(0,1,2), unpack=True,skiprows=1)
   return array



#fid=fopen(fname,'r');
#C=textscan(fid,'%f %f %f','CommentStyle','#');
#fclose(fid);
#%
#pos = [C{1} C{2} C{3}];


#Check these parameters (eg SEFD and moon solid angle)
freq=np.arange(73,233) * 1e6 # Freq in MHz --> Hz
wavelength = 3e8 / freq # Wavelength in meters
Aeff = 14.5 # Eff tile area in m2
Tsky = 3000 * (freq/60e6)**(-2.5) # Tsky is 3000K at 60 MHz with index of -2.5
SEFD = 2.0* 1380 * Tsky / Aeff #; % Sys eq flux density in Jy
df = 1.28e6 # Channel width in Hz
tint = 8 * 3600 # Exposure time in Hrs --> sec
Tmoon = 230.0 # Moon temp in Kelvin
Smoon = (2.0 * 1380 * (Tmoon-Tsky) * (0.5*np.pi/180.0)**2) / (wavelength**2)  # Moon flux density in Jy

#%
#%	No need to change parameters below this line
#%
pos = read_mwa_tile_pos() #; % Tile position in meter (on a local tangent plane)
baselines = make_baselines(pos) # Baselines in meter (snapshot, zenith target)
[um, vm, wm] = baselines


bl = np.arange(5,605,10) #[5:10:605]'; % Baseline bin centers
flux_std=np.zeros(len(freq))
temp_std=np.zeros(len(freq))

for ifreq in np.arange(0,len(freq)): # For each freq channel
   u = um/wavelength[ifreq]
   v = vm/wavelength[ifreq]
   w = wm/wavelength[ifreq] # Baselines in wavelengths
   D = np.sqrt(u**2+v**2) # Baseline distance in wavelengths
   Di = np.round(D / 10.0) # Bin index for sensitivity calc (bin width = 10 wavelengths)

   ##signal =  np.sin(np.pi*bl*0.5*np.pi/180) / (np.pi*bl*0.5*np.pi/180) #; % Moon signal per vis in Jy
   signal = 2.0 * scipy.special.jv(1,(np.pi*bl*0.5*np.pi/180)) / (np.pi*bl*0.5*np.pi/180) # Moon signal per vis in Jy
   #noise = repmat(SEFD(ifreq) / sqrt(2*df*tint),length(bl),1) #Noise per vis in Jy
   noise = np.tile(SEFD[ifreq] / np.sqrt(2*df*tint),[1,len(bl)]) #Noise per vis in Jy
 
   Nvis=np.zeros(len(bl))
   for ibin in np.arange(0,len(bl)): # For each baseline bin (calc # vis)
      I = np.argwhere(Di==ibin)  # I has indicies of all baselines that fall in this bin
      Nvis[ibin] = 2 * len(I) #; % Number of visibilities in that bin, 2 for conjugate vis

   snr = signal / noise * np.sqrt(Nvis) #  SNR per baseline bin
   flux_std[ifreq] = (np.sum(snr**2))**(-0.5)
   temp_std[ifreq] = flux_std[ifreq] * wavelength[ifreq]**2 /2/1380 / ((0.5*np.pi/180)**2)

#print flux_std
#print temp_std

plt.clf()
T_sensitivity_plot=plt.figure(1)
plt.semilogy(freq/1e6,temp_std)
plt.grid(True,which="both")
#plt.title("Global EoR Temperature Uncertainy for MWA")
plt.ylabel('Temperature Uncertainty (K)')
plt.xlabel('Frequency (MHz)')
T_sensitivity_plot.savefig('temp_uncertainty.png')

#figure()
#[hAx,hLine1,hLine2]=plotyy(freq/1e6,flux_std*1e3,freq/1e6,temp_std);
#set(hLine1,'linewidth',2); set(hLine2,'linewidth',2);
#set(hAx(1),'fontsize',16); set(hAx(2),'fontsize',16);
#xlabel('Frequency in MHz','fontsize',16);
#ylabel(hAx(1),'Flux uncertainty in mJy','fontsize',16); ylabel(hAx(2),'Temperature uncertainty in K','fontsize',16);
#print('mwa_moon_sensitivity.eps','-depsc');
