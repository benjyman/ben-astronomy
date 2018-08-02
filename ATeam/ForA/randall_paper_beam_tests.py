#Test size of synthesised beam for wsclean images of ForA for Ph1 and Ph2
import os
import pyfits
import numpy as np

data_dir='/md0/ATeam/ForA/combine_ph1_ph2_ForA_145/'

mwa_phase_list=['phase1','phase2']
#mwa_phase_list=['phase2']
weight_list=['natural','briggs_0','uniform']

D_Ph1 = 3000.
D_Ph2 = 6000.
freq_MHz = 185.
wavelength = 300./freq_MHz
lambda_on_D_rad_Ph1 = wavelength/D_Ph1
lambda_on_D_rad_Ph2 = wavelength/D_Ph2
lambda_on_D_deg_ph1 = lambda_on_D_rad_Ph1/np.pi*180
lambda_on_D_deg_ph2 = lambda_on_D_rad_Ph2/np.pi*180
lambda_on_D_arcmin_ph1 = lambda_on_D_deg_ph1*60.
lambda_on_D_arcmin_ph2 = lambda_on_D_deg_ph2*60.
lambda_on_D_arcsec_ph1 = lambda_on_D_arcmin_ph1*60.
lambda_on_D_arcsec_ph2 = lambda_on_D_arcmin_ph2*60.

print "lambda_on_D_rad_Ph1 is %s deg (%s arcmin or %s arcsec)" % (lambda_on_D_rad_Ph1,lambda_on_D_arcmin_ph1,lambda_on_D_arcsec_ph1)
print "lambda_on_D_rad_Ph2 is %s deg (%s arcmin or %s arcsec)" % (lambda_on_D_rad_Ph2,lambda_on_D_arcmin_ph2,lambda_on_D_arcsec_ph2)


for mwa_phase in mwa_phase_list:
   if mwa_phase=='phase1':
      obsid='1102865728'
      image_prefix='2014B_Phase1_%s' % obsid
      ms_name='%s2014B/%s/%s.ms' % (data_dir,obsid,obsid)
   else:
      obsid='1202816352'
      image_prefix='2018A_Phase2_%s' % obsid
      ms_name='%s2018A/%s/%s.ms' % (data_dir,obsid,obsid)

   for weight in weight_list:
      if weight=='briggs_0':
         weight_string='briggs 0'
         image_name='%s_briggs_0' % (image_prefix)
      else:
         weight_string=weight
         image_name='%s_%s' % (image_prefix,weight)
      cmd='wsclean -name %s -nofitbeam  -theoreticbeam -size 100 100 -niter 10 -mgain 0.95 -datacolumn CORRECTED_DATA -scale 0.004 -weight %s -smallinversion -pol xx  %s' % (image_name,weight_string,ms_name)
      print cmd
      #os.system(cmd)

      XX_image_name="%s-image.fits" % (image_name)
      hdulist = pyfits.open(XX_image_name)
      header=hdulist[0].header
      bmaj_deg=header['bmaj']
      bmin_deg=header['bmin']
      bmaj_arcmin=float(bmaj_deg)*60.
      bmin_arcmin=float(bmin_deg)*60.
      print "%s, %s: BMAJ %s BMIN %s deg, BMAJ %s BMIN %s arcmin " % (mwa_phase,weight,bmaj_deg,bmin_deg,bmaj_arcmin,bmin_arcmin)





#2018A_hdulist = pyfits.open(2018A_XX_image_name)


