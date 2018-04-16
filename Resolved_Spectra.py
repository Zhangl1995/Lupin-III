import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


filename = ['reprojectimg/LH_on2_LL2.fits','reprojectimg/SH_on2_LL2.fits']
#filename = ['reprojectimg/LL1_on2_LL2.fits','reprojectimg/LL2_on2_LL2.fits',
#            'reprojectimg/SL1_on2_LL2.fits','reprojectimg/SL2_on2_LL2.fits',
#            'reprojectimg/LH_on2_LL2.fits','reprojectimg/SH_on2_LL2.fits']

plt.figure(figsize=[20,5])

for nam in filename:
    cube = fits.open(nam)
    cube.info()

# Re-order FLUX, IVAR, and MASK arrays from (wavelength, DEC, RA) to (RA, DEC, wavelength).
    flux = np.transpose(cube['PRIMARY'].data, axes=(2, 1, 0))
    flux_sum = flux.sum(axis=2)
    flux_sum[np.isnan(flux_sum)] = 0
    mm = flux_sum.shape
#    nn = np.argmax(flux_sum)
#    xx = (nn)//mm[1]
#    yy = (nn)%mm[1]
    xx, yy = np.where(flux_sum == np.max(flux_sum))
#    print(flux_sum[xx,yy],np.max(flux_sum))
#    print(xx,yy)
   
    wave = cube['WCS-TAB'].data[0][0]
    wave.resize(len(wave))

    flux_header = cube['PRIMARY'].header
    
#    plt.figure(figsize=[10,5])

    n = np.min([xx[0],yy[0],mm[0]-xx[0],mm[1]-yy[0]])+1
#    print(n)
    for i in range(0,n):
        flux_part = flux[xx[0]-i:xx[0]+1+i,yy[0]-i:yy[0]+1+i]
        flux_part = flux_part.mean(axis=(0,1))
        fmean = np.mean(flux_part[1:6])
        if (fmean > 1):
            plt.plot(wave, np.log10(flux_part/fmean),label=str(2*i+1)+'*'+str(2*i+1)+'$\ pixel^2$')
    m = np.min([xx[0],mm[0]-xx[0]])
    if(m//5>0):
        for j in range(n,m,m//5):
            flux_part = flux[xx[0]-j:xx[0]+1+i,:]
            flux_part = flux_part.mean(axis=(0,1))
            fmean = np.mean(flux_part[1:6])
            if (fmean > 0.1):
                plt.plot(wave, np.log10(flux_part/fmean),label=str(2*j+1)+'$\ columns$')
    else:
        for j in range(n,m,1):
            flux_part = flux[xx[0]-j:xx[0]+1+i,:]
            flux_part = flux_part.mean(axis=(0,1))
            fmean = np.mean(flux_part[1:6])
            if (fmean > 0.1):
                plt.plot(wave, np.log10(flux_part/fmean),label=str(2*j+1)+'$\ columns$')

#plt.xlim(32)
#    plt.ylim(-1)
#plt.title(nam)
plt.xlabel('$\lambda \, [um]$')
plt.ylabel('$log\ $'+flux_header['BUNIT'])
#plt.legend(loc = 2,fontsize='small')
plt.show()
