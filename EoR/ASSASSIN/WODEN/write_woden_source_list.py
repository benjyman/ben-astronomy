def write_woden_sourcelists(hpx_fits_file_list,nside):
    for hpx_fits_filename in hpx_fits_file_list:
       name_base = hpx_fits_filename.split('.fits')[0]
       data = hp.read_map(hpx_fits_filename,nest=False)
       pix_inds = arange(hp.nside2npix(nside))
       l, b = hp.pix2ang(nside,pix_inds,lonlat=True)
       gal_coords = SkyCoord(l*u.deg, b*u.deg, frame='galactic')
       ra = gal_coords.icrs.ra.value
       dec = gal_coords.icrs.dec.value
       fluxes = data*hp.nside2pixarea(nside,degrees=False)*1e+6
       fig = plt.figure(figsize=(10,10))
       hp.mollview(log10(data), sub=(2,1,1), fig=fig,title='Galactic')
       hp.mollview(log10(fluxes), sub=(2,1,2), fig=fig,title='Equatorial')
       fig.savefig('%s_woden_map.png' % name_base, bbox_inches='tight')
       plt.close()
       source_ind = 0
       with open('%s_woden.xt' % name_base,'w') as outfile:
           for ind,flux in enumerate(fluxes):
               if source_ind == 0:
                   outfile.write('SOURCE pygsm P %d G 0 S 0 0\n' %len(fluxes))
               outfile.write('COMPONENT POINT %.7f %.6f\n' %(ra[ind]/15.0,dec[ind]))
               outfile.write('LINEAR 150e+6 %.10f 0 0 0 0.0\n' %flux)
               outfile.write('ENDCOMPONENT\n')
               source_ind += 1
           outfile.write('ENDSOURCE')
