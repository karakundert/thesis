#!/usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
from astropy.io import fits
from wcsaxes import WCS

files = ['point', 'extended', 'M51', 'psf']
effects_files = ['no_perturbation_point', 'illum_offset', 'pointing_cor_point', \
                 'pointing_point', 'size_diff', 'all_effects']
uv_files = ['aperture_real', 'aperture_imag']
fig_shape = 0.15, 0.1, 0.8, 0.8

pl.rc('text', usetex=True)
pl.rc('font', family='serif')

i = 0
pl.clf()

pb = fits.open('pb.fits')[0].data[0,0,:,:]
psf = fits.open('pb.fits')[0].data[0,0,:,:]

for file in files:
    filename = file + '.fits'
    hdu = fits.open(filename)[0]
    data = hdu.data[0,0,:,:]
    max = np.amax(data)

    # get axis information from header
    wcs = WCS(hdu.header)
    wcs = wcs.dropaxis(3)
    wcs = wcs.dropaxis(2)

    # create figure
    fig = pl.figure(i)
    synth_beam = pl.Circle((50,50), 5, color='white', fill=False, alpha=0.5)
    ax = fig.add_axes(fig_shape,projection=wcs)
    ax.contour(pb, levels=[0.2, 0.4, 0.6, 0.8], colors='white', alpha=0.5)
    ra = ax.coords['ra']
    ra.set_axislabel(r'J2000 Right Ascension', fontsize=18)
    ra.set_major_formatter('hh:mm:ss')
    ra.set_ticks(color='white', exclude_overlapping=True, size=18)
    ra.display_minor_ticks('True')
    dec = ax.coords['dec']
    dec.set_axislabel(r'J2000 Declination', fontsize=18)
    dec.set_major_formatter('dd:mm:ss')
    dec.set_ticks(color='white', exclude_overlapping=True, size=18)
    dec.display_minor_ticks('True')
    fig.gca().add_artist(synth_beam)

    pl.imshow(data/max, origin='lower')
    cbar = pl.colorbar()
    cbar.set_label(r'Jy', fontsize=16)
    pl.savefig(files[i]+'.png')

    i += 1

j = 0

for file in effects_files:
    filename = file + '.fits'
    hdu = fits.open(filename)[0]
    data = hdu.data[0,0,:,:]
    max = np.amax(data)

    # get axis information from header
    wcs = WCS(hdu.header)
    wcs = wcs.dropaxis(3)
    wcs = wcs.dropaxis(2)

    # create figure
    fig = pl.figure(i)
    synth_beam = pl.Circle((50,50), 5, color='black', fill=False, alpha=0.5)
    ax = fig.add_axes(fig_shape,projection=wcs)
    ax.contour(pb, levels=[0.2, 0.4, 0.6, 0.8], colors='black', alpha=0.5)
    ra = ax.coords['ra']
    ra.set_axislabel(r'J2000 Right Ascension', fontsize=18)
    ra.set_major_formatter('hh:mm:ss')
    ra.set_ticks(color='black', exclude_overlapping=True, size=18)
    ra.display_minor_ticks('True')
    dec = ax.coords['dec']
    dec.set_axislabel(r'J2000 Declination', fontsize=18)
    dec.set_major_formatter('dd:mm:ss')
    dec.set_ticks(color='black', exclude_overlapping=True, size=18)
    dec.display_minor_ticks('True')
    fig.gca().add_artist(synth_beam)

    imgplt = pl.imshow(data/max, origin='lower', norm=mpl.colors.LogNorm())
    cbar = pl.colorbar()
    cbar.set_label(r'log(Jy)', fontsize=16)
    pl.savefig(effects_files[j]+'.png')

    i += 1
    j += 1
    
j = 0

for file in uv_files:
    filename = file + '.fits'
    hdu = fits.open(filename)[0]
    data = hdu.data[0,0,:,:]
    shape = data.shape
    hdr = hdu.header
    axes = []

    for k in range(1,5):
        axis = {
                'type': hdr['CTYPE%d' % k], # projection type
                'rval': hdr['CRVAL%d' % k], # value of reference pixel
                'delt': hdr['CDELT%d' % k], # physical increment at reference pixel
                'rpix': hdr['CRPIX%d' % k], # reference pixel
                'unit': hdr['CUNIT%d' % k], # units
        }
        axis['start'] = axis['rval'] - axis['delt'] * (axis['rpix'] - 1)
        axis['end'] = axis['rval'] + axis['delt'] * (axis['rpix'] - 1)
        if axis['unit']:
            axis['label'] = '%s (%s)' % (axis['type'], axis['unit'])
            axis['value'] = '%s = %s %s' % (axis['type'], axis['rval'], axis['unit'])
        else:
            axis['label'] = axis['type']
            axis['value'] = '%s = %s' % (axis['type'], axis['rval'])
        axes.append(axis)

    fig = pl.figure(i)
    extent = [
            axes[0]['start']/2, axes[0]['end']/2,
            axes[1]['start']/2, axes[1]['end']/2
    ]
    ax = fig.add_axes(fig_shape)

    pl.xlabel(axes[0]['label'], fontsize=18)
    pl.ylabel(axes[1]['label'], fontsize=18)
    if j == 0:
        pl.tick_params(axis='both', direction = 'in', color='white', labelsize=18)
    else:
        pl.tick_params(axis='both', direction = 'in', color='black', labelsize=18)

    pl.imshow(data[383:639, 383:639], extent=extent, origin='lower')
    pl.savefig(uv_files[j]+'.png')

    i += 1
    j += 1
    
pl.show()
