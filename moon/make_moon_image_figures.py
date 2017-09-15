#Script to make moon image figures for papers and presentations


 
def make_images():
    
   frequency_chan=options.frequency
    
   #Plot images:
   py.figure(1)
   py.clf()
   py.title("Dirty Moon Image")
   py.imshow( ( moon_zoom ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("moon_image.png")
   
   py.figure(2)
   py.clf()
   py.title("Moon-Sky Difference Image")
   py.imshow( ( moon_minus_sky ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("moon_minus_sky.png")
   
   py.figure(3)
   py.clf()
   py.title("Moon Mask")
   py.imshow( ( moon_mask ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_mask.png")
   
   py.figure(4)
   py.clf()
   py.title("Moon Mask Convolved with PSF")
   py.imshow( ( G_image_shift ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_mask_psf_convol.png")
   
   
   py.figure(5)
   py.clf()
   py.title("Reconstructed Moon Dirty Image")
   py.imshow( ( reconstructed_moon ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_reconstructed.png")
   
   
   py.figure(6)
   py.clf()
   py.title("Reconstructed RFI-only Dirty Image")
   py.imshow( ( RFI_alone ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_only.png")
   
   py.figure(7)
   py.clf()
   py.title("Dirty Off-Moon Image")
   py.imshow( ( off_moon_zoom ), cmap=py.cm.Greys,vmin=-13.0, vmax=22.0,origin='lower')
   py.colorbar()
   py.savefig("off_moon_image.png")
   
   py.figure(8)
   py.clf()
   py.title("PSF Image")
   py.imshow( ( psf_zoom ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("PSF_image.png")
   
   
   py.figure(9)
   py.clf()
   py.title("Residual Image")
   py.imshow( ( R ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("residual_image.png")
   
   #broadened RFI
   py.figure(10)
   py.clf()
   py.title("RFI Broadening Image")
   py.imshow( ( RFI_broadening ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_broadening_image.png")
   
   py.figure(11)
   py.clf()
   py.title("RFI Broadened Convolved Image")
   py.imshow( ( RFI_convolved_shift ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_broadened_convol_image.png")
   
   py.figure(12)
   py.clf()
   py.title("Broadened Reconstructed RFI-only Dirty Image")
   py.imshow( ( RFI_alone2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("RFI_only_broadened.png")
   
   py.figure(13)
   py.clf()
   py.title("Residual Image with Broadened RFI")
   py.imshow( ( R2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("residual_image_broadened_rfi.png")
   
   py.figure(14)
   py.clf()
   py.title("Reconstructed Moon Dirty RFI Broadened Image")
   py.imshow( ( reconstructed_moon2 ), cmap=py.cm.Greys,origin='lower')
   py.colorbar()
   py.savefig("moon_reconstructed_rfi_broadened.png")

  
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: cotter_moon.py [obsid] [options]'

parser = OptionParser(usage=usage)

parser.add_option('--frequency',type='float', dest='frequency',default=None,help='Frequency channel to make Moon images for. e.g. --frequency=150 [default=%default]')


(options, args) = parser.parse_args()

#freq = args[0]

make_images(options)
   