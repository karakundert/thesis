#!/usr/bin/env python

import pyfits

extended = pyfits.open('aperture_real.fits')
print extended.info()
header = extended[0].header

for each in header:
    if each != 'HISTORY':
        print each, '=', header[each]
