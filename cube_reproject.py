import numpy as np
from astropy import wcs
from astropy.io import fits
from reproject import reproject_exact
#import matplotlib.pyplot as plt

filename = ['ngc4826_DR5_LL1_cube.fits','ngc4826_DR5_LL2_cube.fits','ngc4826_DR5_SL1_cube.fits','ngc4826_DR5_SL2_cube.fits','ngc4826_DR5_SH_cube.fits','ngc4826_DR5_LH_cube.fits']
#filename = ['ngc4826_DR5_LL2_cube.fits']

for nam in filename:
	hdu1 = fits.open(nam)[0]
	hdu2 = fits.open('ngc4826_DR5_LL2_cube.fits')[0]

	CARD1 = [('NAXIS1',hdu1.header['NAXIS1']),('NAXIS2',hdu1.header['NAXIS2']),('BUNIT',hdu1.header['BUNIT']),('CTYPE1',hdu1.header['CTYPE1']),
	         ('CTYPE2',hdu1.header['CTYPE2']),('CRVAL1',hdu1.header['CRVAL1']),('CRVAL2',hdu1.header['CRVAL2']),('CRPIX1',hdu1.header['CRPIX1']),
	         ('CRPIX2',hdu1.header['CRPIX2']),('PC1_1',hdu1.header['PC1_1']),('PC1_2',hdu1.header['PC1_2']),('PC2_1',hdu1.header['PC2_1']),
	         ('PC2_2',hdu1.header['PC2_2']),('CDELT1',hdu1.header['CDELT1']),('CDELT2',hdu1.header['CDELT2'])]
	hdr1 = fits.Header(cards= CARD1)
	wcs1 = wcs.WCS(hdr1)

	CARD2 = [('NAXIS1',hdu2.header['NAXIS1']),('NAXIS2',hdu2.header['NAXIS2']),('BUNIT',hdu2.header['BUNIT']),('CTYPE1',hdu2.header['CTYPE1']),
	         ('CTYPE2',hdu2.header['CTYPE2']),('CRVAL1',hdu2.header['CRVAL1']),('CRVAL2',hdu2.header['CRVAL2']),('CRPIX1',hdu2.header['CRPIX1']),
	         ('CRPIX2',hdu2.header['CRPIX2']),('PC1_1',hdu2.header['PC1_1']),('PC1_2',hdu2.header['PC1_2']),('PC2_1',hdu2.header['PC2_1']),
	         ('PC2_2',hdu2.header['PC2_2']),('CDELT1',hdu2.header['CDELT1']),('CDELT2',hdu2.header['CDELT2'])]
	hdr2 = fits.Header(cards= CARD2)
	wcs2 = wcs.WCS(hdr2)

#print(wcs1,wcs2)
	Spec = []
	for i in range(0,hdu1.header['NAXIS3']):
		array_exact = reproject_exact((hdu1.data[i],wcs1), wcs2, shape_out = hdu2.data[0].shape)
		Spec.append(list(array_exact[0]))
	Spec = np.array(Spec)

	hdr2['NAXIS'] = 3
	hdr2['NAXIS3'] = hdu2.header['NAXIS3']
#fits.writeto('L1_on2_L2.fits', array_exact[0], header = hdr2, overwrite = True)

	primary_hdu = fits.PrimaryHDU(Spec, header=hdr2)
	hdu_wav = fits.open(nam)[1]
	hdul = fits.HDUList([primary_hdu, hdu_wav])
	hdul.writeto('reprojectimg/' + nam.split('_')[2] + '_on2_LL2.fits', overwrite = True)
