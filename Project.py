import numpy as np
from astropy.io import ascii
from scipy.interpolate import interp1d
import astropy.units as u
ls_mic = 2.99792458e14 #micron/s

def blackcgs(T, nu):
    h = 6.625*10**-27
    k = 1.38*10**-16
    cs = 2.99792458e10
    firstterm = 2*h*nu**3/cs**2
    x = h*nu/(k*T)
    exppart = 1/(np.exp(x)-1)
    bnu = firstterm*exppart
# calculate blam and return it if keyword lamb is set
#if keyword_set(lamb) then begin
#blam = bnu*nu/(3*10**10/nu) 
#return(blam)
    return(bnu)

def BandFunc_mono(wave0, wavelength, spc, rsr, alpha):
    freq = (ls_mic / wavelength)
    nu0 = ls_mic / wave0
    Sv = freq**alpha
    fluxFltr = nu0**alpha*np.trapz(rsr*spc, freq) / np.trapz(rsr*Sv, freq)
    return (fluxFltr)

def BandFunc_mips(wave0, wavelength, spc, rsr, T):
    freq = (ls_mic / wavelength)
    nu0 = ls_mic / wave0
    fluxFltr = blackcgs(T, nu0)*np.trapz(rsr*spc, freq) / np.trapz(rsr*blackcgs(T, freq), freq)
    return (fluxFltr)

def BandFunc_mean(wavelength, spc, rsr):
    signal = np.trapz(rsr*spc/wavelength, x=wavelength)
    norm   = np.trapz(rsr/wavelength, x=wavelength)
    fluxFltr = signal/norm
    return (fluxFltr)

monoFilters = ["Herschel_PACS_70", "Herschel_PACS_100", "Herschel_PACS_160",
               "Herschel_SPIRE_250", "Herschel_SPIRE_350", "Herschel_SPIRE_500",
               "Herschel_SPIRE_250_e", "Herschel_SPIRE_350_e", "Herschel_SPIRE_500_e",
               "Spitzer_IRAC1", "Spitzer_IRAC2", "Spitzer_IRAC3", "Spitzer_IRAC4",
               "IRAS_12", "IRAS_25", "IRAS_60", "IRAS_100","IRC_S7","IRC_S11","IRC_L15","IRC_L24"]
mipsFilters = ["Spitzer_MIPS_24", "Spitzer_MIPS_70", "Spitzer_MIPS_160"]
meanFilters = ["2MASS_J", "2MASS_H", "2MASS_Ks",
               "WISE_w1", "WISE_w2", "WISE_w3", "WISE_w4"]

filterDict = {
    "2MASS_J": 1.235,
    "2MASS_H": 1.662,
    "2MASS_Ks": 2.159,
    "WISE_w1": 3.353,
    "WISE_w2": 4.603,
    "WISE_w3": 11.561,
    "WISE_w4": 22.088,
    "IRC_S7": 7.0,
    "IRC_S11": 11.0,
    "IRC_L15": 15.0,
    "IRC_L24": 24.0,
    "Herschel_PACS_70": 70.,
    "Herschel_PACS_100": 100.,
    "Herschel_PACS_160": 160.,
    "Herschel_SPIRE_250": 250.,
    "Herschel_SPIRE_350": 350.,
    "Herschel_SPIRE_500": 500.,
    "Spitzer_IRAC1": 3.550,
    "Spitzer_IRAC2": 4.493,
    "Spitzer_IRAC3": 5.731,
    "Spitzer_IRAC4": 7.872,
    "IRAS_12": 12.,
    "IRAS_25": 25.,
    "IRAS_60": 60.,
    "IRAS_100": 100.,
    "Spitzer_MIPS_24": 24.,
    "Spitzer_MIPS_70": 70.,
    "Spitzer_MIPS_160": 160.,
    "JCMT_SCUBA1_450": 450.,
    "JCMT_SCUBA1_850": 850.,
}

def Photometry(wavelength, spectrum, band):

    bandPath = '/Users/zhangl/Applications/Fitter/sedfit/filters/'
    bandFile = "{0}.dat".format(band)
    bandPck = np.genfromtxt(bandPath+bandFile)
    bandWave = bandPck[:, 0]
    bandRsr = bandPck[:, 1]
    ma = bandWave>=5.24281693
    bandWave = bandWave[ma]
    bandRsr = bandRsr[ma]
    all_wave = np.concatenate((wavelength, bandWave))
    sort_all = sorted(all_wave)
    unique_sort_all = np.unique(sort_all)
    ma1 = unique_sort_all > np.min(bandWave)
    ma2 = unique_sort_all < np.max(bandWave)
    new_wavelengths = unique_sort_all[ma1*ma2]
    fs = interp1d(wavelength, spectrum, kind='cubic')
    spc = fs(new_wavelengths)
    fr = interp1d(bandWave, bandRsr, kind='cubic')
    rsr = fr(new_wavelengths)
    if band in monoFilters:
        Phot = BandFunc_mono(filterDict[band], new_wavelengths, spc, rsr, -1)
    elif band in mipsFilters:
        Phot = BandFunc_mips(filterDict[band], new_wavelengths, spc, rsr, 10000)
    elif band in meanFilters:
        Phot = BandFunc_mean(new_wavelengths, spc, rsr)
    else:
        raise ValueError("The input band ({0}) is incorrect!".format(bn))
    return(Phot)

def wcss(hdu, x, y):
    from astropy import wcs
    from astropy.io import fits
    if('CROTA2' in hdu):
        CARD = [('NAXIS1',hdu['NAXIS1']),('NAXIS2',hdu['NAXIS2']),('BUNIT',hdu['BUNIT']),('CTYPE1',hdu['CTYPE1']),
                ('CTYPE2',hdu['CTYPE2']),('CRVAL1',hdu['CRVAL1']),('CRVAL2',hdu['CRVAL2']),('CRPIX1',hdu['CRPIX1']),
                ('CRPIX2',hdu['CRPIX2']),('CROTA2',hdu['CROTA2']),('CDELT1',hdu['CDELT1']),('CDELT2',hdu['CDELT2'])]
    else:
        CARD = [('NAXIS1',hdu['NAXIS1']),('NAXIS2',hdu['NAXIS2']),('BUNIT',hdu['BUNIT']),('CTYPE1',hdu['CTYPE1']),
                ('CTYPE2',hdu['CTYPE2']),('CRVAL1',hdu['CRVAL1']),('CRVAL2',hdu['CRVAL2']),('CRPIX1',hdu['CRPIX1']),
                ('CRPIX2',hdu['CRPIX2']),('PC1_1',hdu['PC1_1']),('PC2_1',hdu['PC2_1']),('PC1_2',hdu['PC1_2']),
                ('PC2_2',hdu['PC2_2']),('CDELT1',hdu['CDELT1']),('CDELT2',hdu['CDELT2'])]
    hdr = fits.Header(cards= CARD)
    w = wcs.WCS(hdr)
    recx, recy = w.wcs_pix2world(x, y ,0)
    return [recx, recy]

def wcss3(hdu, x, y):
    from astropy import wcs
    from astropy.io import fits
    if('CROTA2' in hdu):
        CARD = [('NAXIS1',hdu['NAXIS1']),('NAXIS2',hdu['NAXIS2']),('BUNIT',hdu['BUNIT']),('CTYPE1',hdu['CTYPE1']),
                ('CTYPE2',hdu['CTYPE2']),('CRVAL1',hdu['CRVAL1']),('CRVAL2',hdu['CRVAL2']),('CRPIX1',hdu['CRPIX1']),
                ('CRPIX2',hdu['CRPIX2']),('CROTA2',hdu['CROTA2']),('CDELT1',hdu['CDELT1']),('CDELT2',hdu['CDELT2'])]
    else:
        CARD = [('NAXIS1',hdu['NAXIS1']),('NAXIS2',hdu['NAXIS2']),('BUNIT',hdu['BUNIT']),('CTYPE1',hdu['CTYPE1']),
                ('CTYPE2',hdu['CTYPE2']),('CRVAL1',hdu['CRVAL1']),('CRVAL2',hdu['CRVAL2']),('CRPIX1',hdu['CRPIX1']),
                ('CRPIX2',hdu['CRPIX2']),('PC1_1',hdu['PC1_1']),('PC2_1',hdu['PC2_1']),('PC1_2',hdu['PC1_2']),
                ('PC2_2',hdu['PC2_2']),('CDELT1',hdu['CDELT1']),('CDELT2',hdu['CDELT2'])]
    hdr = fits.Header(cards= CARD)
    w = wcs.WCS(hdr)
    recx, recy = w.wcs_world2pix(x, y ,0)
    return [recx, recy]

def Montage_Reproject(name, mode, O2, dirsc, scale, rp, rpu, rpo):
    import os
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table

    path = '/Users/zhangl/desktop/SINGS/{0}/{1}/'.format(dirsc, name)
    rp_file = path + rp
    rpu_file = path + rpu
    if(os.path.isfile(rp_file)==False):
        print('No such path or file:{0}'.format(rp_file))
        return ()
    rp_fits = fits.open(rp_file)
    data = rp_fits[0].data
    data = data.sum(axis=0)
    data = np.transpose(data, axes=(1,0))
    header = rp_fits[0].header
    indcube = np.full((header['NAXIS1']+20, header['NAXIS2']+20), 0)
    mask = np.isnan(data)
    indcube[10:-10,10:-10][~mask] = 1
    header['NAXIS1'] = header['NAXIS1'] + 20
    header['NAXIS2'] = header['NAXIS2'] + 20
    header['CRPIX1'] = header['CRPIX1'] + 10
    header['CRPIX2'] = header['CRPIX2'] + 10
    indcube = np.transpose(indcube, axes=(1, 0))
    img = fits.PrimaryHDU(indcube, header=header)
    dirs = 'IndC/{0}'.format(name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    img.writeto('{0}/ind_{1}.fits'.format(dirs, mode), overwrite = True)

    p_file = path + rpo
    if(os.path.isfile(p_file)==False):
        print('No such path or file:{0}'.format(p_file))
        return()
    flux_header = fits.open(p_file)[0].header
    dirs = 'HdrText_{0}/{1}_{2}/{3}'.format(O2, dirsc, scale, name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    hdr = open('HdrText_{0}/{1}_{2}/{3}/{3}_{4}_hdr.txt'.format(O2, dirsc, scale, name, mode),'w')
    hdr.write('SIMPLE  = T\nBITPIX  = -64\nNAXIS   = 3\n')
    hdr.write('NAXIS1  = {0}\n'.format(str(flux_header['NAXIS1']+2)))
    hdr.write('NAXIS2  = {0}\n'.format(str(flux_header['NAXIS2']+2)))
    hdr.write('NAXIS3  = {0}\n'.format(header['NAXIS3']))
    hdr.write('CTYPE1  = {0}\n'.format(flux_header['CTYPE1']))
    hdr.write('CTYPE2  = {0}\n'.format(flux_header['CTYPE2']))
    hdr.write('EQUINOX  = {0}\n'.format(flux_header['EQUINOX']))
    hdr.write('CRVAL1  = {0}\n'.format(flux_header['CRVAL1']))
    hdr.write('CRVAL2  = {0}\n'.format(flux_header['CRVAL2']))
    hdr.write('CRVAL3  = {0}\n'.format(flux_header['CRVAL3']))
    hdr.write('CRPIX1  = {0}\n'.format(str(flux_header['CRPIX1']+1)))
    hdr.write('CRPIX2  = {0}\n'.format(str(flux_header['CRPIX2']+1)))
    hdr.write('CRPIX3  = {0}\n'.format(flux_header['CRPIX3']))
    hdr.write('CDELT1  = {0}\n'.format(str(10/3600)))
    hdr.write('CDELT2  = {0}\n'.format(str(10/3600)))
    hdr.write('CDELT3  = {0}\n'.format(str(0)))
#        hdr.write('PC1_1  = {0}\n'.format(p_fits[0].header['PC1_1']))
#        hdr.write('PC2_1  = {0}\n'.format(p_fits[0].header['PC2_1']))
#        hdr.write('PC1_2  = {0}\n'.format(p_fits[0].header['PC1_2']))
#        hdr.write('PC2_2  = {0}\n'.format(p_fits[0].header['PC2_2']))
    Ro2_s = np.arcsin(flux_header['PC2_1']*flux_header['CDELT2']/flux_header['CDELT1'])*180/np.pi
    Ro22_s = [round(Ro2_s-(Ro2_s/np.abs(Ro2_s)-1)*180,7), round(Ro2_s*180/np.abs(Ro2_s)-Ro2_s-(Ro2_s/np.abs(Ro2_s)-1)*180,7)]
    Ro2_c = np.arccos(flux_header['PC1_1'])*180/np.pi
    Ro22_c = [round(Ro2_c,7), round(360-Ro2_c,7)]
    Ro2 = list(set(Ro22_s).intersection(Ro22_c))[0]
    hdr.write('CROTA2  = {0}\n'.format(str(Ro2)))
    hdr.write('BUNIT  = MJy_sr^-1\n')
    hdr.write('END')
    hdr.close()
    dirs = 'Projection_{0}/{1}_{2}/{3}'.format(O2, dirsc, scale, name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    ind_file = 'IndC/{0}/ind_{1}.fits'.format(name, mode)
    os.system('mProjectCube -f {0} Projection_{1}/{2}_{3}/{4}/{5}_on2_{1}.fits HdrText_{1}/{2}_{3}/{4}/{4}_{5}_hdr.txt'.format(rp_file,O2,dirsc,scale,name,mode))
    os.system('mProject -f {0} Projection_{1}/{2}_{3}/{4}/{4}_{5}ind.fits HdrText_{1}/{2}_{3}/{4}/{4}_{5}_hdr.txt'.format(ind_file,O2,dirsc,scale,name,mode))
    Montage = '/Users/zhangl/Software/Montage2/bin/'
    os.system('{0}mProjectCube -f {1} Projection_{2}/{3}_{4}/{5}/{6}_on2_{2}_unc.fits HdrText_{2}/{3}_{4}/{5}/{5}_{6}_hdr.txt'.format(Montage,rpu_file,O2,dirsc,scale,name,mode))
#        sr1 = np.abs(rp_fits[0].header['CDELT1'])**2
#        sr2 = np.abs(p_fits[0].header['CDELT1'])**2

    rpd_file = 'Projection_{0}/{1}_{2}/{3}/{4}_on2_{0}.fits'.format(O2,dirsc,scale,name,mode)
    Spec_hdu = fits.open(rpd_file)[0]
#        Spec_hdu.data = Spec_hdu.data*sr1/sr2
    Spec_hdu.header['BUNIT'] = 'MJy/sr'
    Wave_hdu = rp_fits[1]
    hdul = fits.HDUList([Spec_hdu,Wave_hdu])
    hdul.writeto('Projection_{0}/{1}_{2}/{3}/{4}_on2_{0}.fits'.format(O2,dirsc,scale,name,mode), overwrite = True)
    rpud_file = 'Projection_{0}/{1}_{2}/{3}/{4}_on2_{0}_unc.fits'.format(O2,dirsc,scale,name,mode)
    Spec_hdu = fits.open(rpud_file)[0]
#        Spec_hdu.data = Spec_hdu.data*sr1/sr2
    Spec_hdu.header['BUNIT'] = 'MJy/sr'
    Wave_hdu = rp_fits[1]
    hdul = fits.HDUList([Spec_hdu,Wave_hdu])
    hdul.writeto('Projection_{0}/{1}_{2}/{3}/{4}_on2_{0}_unc.fits'.format(O2,dirsc,scale,name,mode), overwrite = True)
    rp_fits.close()
    return ()

def Montage_Reproject_Img(name, band, O2, dirsc, scale, dirsi):
    import os
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table

    p_file = '/Users/zhangl/desktop/SINGS/Cube/{0}/{0}_DR5_{1}_cube.fits'.format(name,O2)
    if(os.path.isfile(p_file)==False):
        print('No such path or file:{0}'.format(p_file))
        return ()
    else:
        Img_header = fits.open(p_file)[0].header
    dirs = 'HdrText_{0}/{1}_{2}/{3}'.format(O2, dirsc, scale, name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    hdr = open('{0}/{1}_hdr.txt'.format(dirs, band),'w')
    rp_file = '/Users/zhangl/Downloads/SINGS/{0}/{1}_{2}.fits'.format(dirsi, name, band)
    if(os.path.isfile(rp_file)==False):
        print('No such path or file:{0}'.format(rp_file))
        return ()
    hdr.write('SIMPLE  = T\nBITPIX  = -64\nNAXIS   = 2\n')
    hdr.write('NAXIS1  = {0}\n'.format(str(Img_header['NAXIS1']*2+1)))
    hdr.write('NAXIS2  = {0}\n'.format(str(Img_header['NAXIS2']*2+1)))
    hdr.write('CTYPE1  = {0}\n'.format(Img_header['CTYPE1']))
    hdr.write('CTYPE2  = {0}\n'.format(Img_header['CTYPE2']))
    hdr.write('EQUINOX  = {0}\n'.format(Img_header['EQUINOX']))
    hdr.write('CRVAL1  = {0}\n'.format(str(Img_header['CRVAL1'])))
    hdr.write('CRVAL2  = {0}\n'.format(str(Img_header['CRVAL2'])))
    hdr.write('CRPIX1  = {0}\n'.format(str(Img_header['CRPIX1']*2)))
    hdr.write('CRPIX2  = {0}\n'.format(str(Img_header['CRPIX2']*2)))
    hdr.write('CDELT1  = {0}\n'.format(str(10/3600)))
    hdr.write('CDELT2  = {0}\n'.format(str(10/3600)))
#        hdr.write('PC1_1  = {0}\n'.format(p_fits[0].header['PC1_1']))
#        hdr.write('PC2_1  = {0}\n'.format(p_fits[0].header['PC2_1']))
#        hdr.write('PC1_2  = {0}\n'.format(p_fits[0].header['PC1_2']))
#        hdr.write('PC2_2  = {0}\n'.format(p_fits[0].header['PC2_2']))
    Ro2_s = np.arcsin(Img_header['PC2_1']*Img_header['CDELT2']/Img_header['CDELT1'])*180/np.pi
    Ro22_s = [round(Ro2_s-(Ro2_s/np.abs(Ro2_s)-1)*180,7), round(Ro2_s*180/np.abs(Ro2_s)-Ro2_s-(Ro2_s/np.abs(Ro2_s)-1)*180,7)]
    Ro2_c = np.arccos(Img_header['PC1_1'])*180/np.pi
    Ro22_c = [round(Ro2_c,7), round(360-Ro2_c,7)]
    Ro2 = list(set(Ro22_s).intersection(Ro22_c))[0]
    hdr.write('CROTA2  = {0}\n'.format(str(Ro2)))
    hdr.write('BUNIT  = MJy_sr^-1\n')
    hdr.write('END')
    hdr.close()
    path = 'Projection_{0}/{1}_{2}/{3}'.format(O2, dirsi, scale, name)
    if(os.path.exists(path) != True):
        os.makedirs(path)
    os.system('mProject -f {0} {1}/{2}_{3}.fits HdrText_{4}/{5}_{6}/{2}/{3}_hdr.txt'.format(rp_file,path,name,band,O2,dirsc,scale))
    photd = fits.open('Projection_{0}/{1}_{2}/{3}/{3}_{4}.fits'.format(O2, dirsi, scale, name, band), mode = 'update')
    photd[0].header['BUNIT'] = 'MJy/sr'
    photd.close()
    return ()

if __name__ == '__main__':

    gal_list = open('nuc_region_gal.txt','r')
    gal_name = gal_list.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_list.close()
    gal_name = ['ngc5194']
    mode_name = ['SL2', 'SL1', 'LL2', 'LL1']
    Band_list = ['IRAC4', 'MIPS24']
    for name in gal_name:
        for mode in mode_name:
            try:
                Montage_Reproject(name, mode, 'LL2', 'Cube_convl', 10, '{0}_{1}_convolved.fits'.format(name,mode), '{0}_{1}_unc_convolved.fits'.format(name,mode), '{0}_LL2_convolved.fits'.format(name))
                print('Projection of {0} {1} image is OK!'.format(name, mode))
            except:
                print('Something worry happens when reprojecting {0} {1} image!'.format(name, mode))
        for band in Band_list:
            try:
                Montage_Reproject_Img(name, band, 'LL2', 'Cube_convl', 10, '{0}_convl'.format(band[:4]))
                print('Projection of {0} {1} imaging is OK!'.format(name, band))
            except:
                print('Something worry happens when reproject {0} {1} imaging!'.format(name))
import numpy as np
from scipy.optimize import leastsq

def func(x, p):
    a, b, c = p
    y = a*x**2 + b*x + c
    return y

def residual(pp, x, y):
    if(x<0):
        x = -x
        diff = pp[-1]*y - func(x, pp[0:3])
    else:
        diff = y - func(x, pp[0:3])
    return diff

def res_vec(pp, x, y):
    dd = np.zeros(x.shape)
    for ii in range(len(dd)):
        dd[ii] = residual(pp, x[ii], y[ii])
    return dd

residual_func = lambda pp, x, y: res_vec(pp, x, y)

def Cal_k(pp, xx, yy):
    k = leastsq(residual_func, pp, args =(xx, yy))[0]
    return k

def Spectra_redis(name, O2, dirsc, scale):
    import os
    import numpy as np
    from astropy import wcs
    from astropy.io import fits
    from astropy.io import ascii
    from astropy import constants
    from astropy.table import Table
    import matplotlib.pyplot as plt


    data_mode = ['SL2', 'SL1', 'LL2', 'LL1']
    nnmm = 'Projection_{0}/{1}_{2}/{3}/'.format(O2, dirsc, scale, name)
    aabb = nnmm + 'SL2_on2_{0}.fits'.format(O2)
    if(os.path.isfile(aabb)):
        cube = fits.open(aabb)
    else:
        print('No such file:{0}'.format(aabb))
        return ()
    flux_header = cube[0].header
    apple = max(flux_header['NAXIS1'],flux_header['NAXIS2']) + 10
    crpixx, crpixy = flux_header['NAXIS1']//2, flux_header['NAXIS2']//2
    coor_mode = wcss(flux_header, [crpixx], [crpixy])
    crval = [[coor_mode[0][0]],[coor_mode[1][0]]]
    cube.close()
    for mode in data_mode:
        aabb = nnmm + '{0}_on2_{1}.fits'.format(mode, O2)
        aacc = nnmm + '{0}_on2_{1}_unc.fits'.format(mode, O2)
        if(os.path.isfile(aabb)):
            cube = fits.open(aabb)
            cube_u = fits.open(aacc)
        else:
            print('No such file:{0}'.format(aabb))
            continue
        flux = np.transpose(cube['PRIMARY'].data, axes=(2, 1, 0))
        flux_u = np.transpose(cube_u['PRIMARY'].data, axes=(2, 1, 0))
        flux_header = cube['PRIMARY'].header
        orange = flux_header['NAXIS3']
        banana = np.full((apple, apple, orange), np.nan)
        banana_u = np.full((apple, apple, orange), np.nan)
        crpix = wcss3(flux_header, [crval[0][0]], [crval[1][0]])
        xx = int(round(crpix[0][0],0))
        yy = int(round(crpix[1][0],0))
        mm1 = max(0, xx-apple//2)
        mm2 = min(xx+apple-apple//2,flux_header['NAXIS1'])
        nn1 = max(0, yy-apple//2)
        nn2 = min(yy+apple-apple//2,flux_header['NAXIS2'])
        for ii in range(mm1,mm2):
            for jj in range(nn1,nn2):
                banana[apple//2 + ii-xx, apple//2 + jj-yy] = flux[ii, jj]
                banana_u[apple//2 + ii-xx, apple//2 + jj-yy] = flux_u[ii, jj]
        CARD = [('NAXIS1',apple),('NAXIS2',apple),('NAXIS3',orange),('CTYPE1',flux_header['CTYPE1']),('CTYPE2',flux_header['CTYPE2']),
                ('EQUINOX',flux_header['EQUINOX']),('CRVAL1',crval[0][0]),('CRVAL2',crval[1][0]),('CRVAL3',flux_header['CRVAL3']),
                ('CRPIX1',apple//2+1),('CRPIX2',apple//2+1),('CRPIX3',flux_header['CRPIX3']),('CDELT1',flux_header['CDELT1']),
                ('CDELT2',flux_header['CDELT2']),('CDELT3',flux_header['CDELT3']),('CROTA2',flux_header['CROTA2']),('BUNIT',flux_header['BUNIT'])]
        hdr = fits.Header(cards= CARD)
        banana = np.transpose(banana, axes=(2, 1, 0))
        primary_hdu = fits.PrimaryHDU(banana, header=hdr)
        hdu_wav = cube['WCS-TAB']
        hdul = fits.HDUList([primary_hdu, hdu_wav])
        dirs = 'ReProjection/{0}_{1}/{2}'.format(dirsc, scale, name)
        if(os.path.exists(dirs) != True):
            os.makedirs(dirs)
        hdul.writeto('ReProjection/{0}_{1}/{2}/{3}_on2_{4}.fits'.format(dirsc,scale,name,mode,O2), overwrite = True)
        banana_u = np.transpose(banana_u, axes=(2, 1, 0))
        primary_hdu = fits.PrimaryHDU(banana_u, header=hdr)
        hdu_wav = cube['WCS-TAB']
        hdul = fits.HDUList([primary_hdu, hdu_wav])
        hdul.writeto('ReProjection/{0}_{1}/{2}/{3}_on2_{4}_unc.fits'.format(dirsc,scale,name,mode,O2), overwrite = True)
        cube.close()
        cube_u.close()
###---------------------------------------------------------------------------------------------------------------
### here to scale the spectra within different wavelength range
###---------------------------------------------------------------------------------------------------------------
    Spec = []
    Spec_u = []
    Wave = []
    for mode in data_mode:
        aadd = 'ReProjection/{0}_{1}/{2}/{3}_on2_{4}.fits'.format(dirsc,scale,name,mode,O2)
        aaee = 'ReProjection/{0}_{1}/{2}/{3}_on2_{4}_unc.fits'.format(dirsc,scale,name,mode,O2)
        if(os.path.isfile(aadd)):
            cube = fits.open(aadd)
            cube_u = fits.open(aaee)
        else:
            print('No such file:{0}'.format(aadd))
            continue
        Sspec = cube[0].data
        Sspec_u = cube_u[0].data
        Wwave = cube[1].data
        Spec.append(Sspec)
        Spec_u.append(Sspec_u)
        Wave.append(Wwave)
    Spec = np.array(Spec)
    Spec_u = np.array(Spec_u)
    Wave = np.array(Wave)
    shape = Spec[0][0].shape
    S2_S1 = np.ones(Spec[0][0].shape)
    ma1 = Wave[0][0][0].sum(axis=1)>7.4
    ma2 = Wave[1][0][0].sum(axis=1)<7.9
    for ii in range(shape[0]):
        for jj in range(shape[1]):
            ww = np.array(list(-Wave[0][0][0].sum(axis=1)[ma1]) + list(Wave[1][0][0].sum(axis=1)[ma2]))
            ss = np.array(list(Spec[0][ma1,ii,jj]) + list(Spec[1][ma2,ii,jj]))
            if(np.isnan(ss.sum()) != True):
                kk = Cal_k([0,1,1,1], ww, ss)
                S2_S1[ii,jj] = kk[-1]
#                plt.scatter(Wave[0][0][0].sum(axis=1)[ma1], S2_S1[ii,jj]*Spec[0][ma1,ii,jj], marker='+', color = 'b')
#                plt.scatter(Wave[1][0][0].sum(axis=1)[ma2], Spec[1][ma2,ii,jj], marker='+', color = 'r')
#                plt.scatter(np.abs(ww), func(np.abs(ww),kk[0:3]), color = 'k', marker='x')
#                plt.show()
    L1_L2 = np.ones(Spec[0][0].shape)
    ma11 = Wave[2][0][0].sum(axis=1)>19.8
    ma12 = Wave[2][0][0].sum(axis=1)<20.6
    ma1 = ma11*ma12
    ma2 = Wave[3][0][0].sum(axis=1)<21.6
    for ii in range(shape[0]):
        for jj in range(shape[1]):
            ww = np.array(list(Wave[2][0][0].sum(axis=1)[ma1]) + list(-Wave[3][0][0].sum(axis=1)[ma2]))
            ss = np.array(list(Spec[2][ma1,ii,jj]) + list(Spec[3][ma2,ii,jj]))
            if(np.isnan(ss.sum()) != True):
                kk = Cal_k([0,1,1,1], ww, ss)
                L1_L2[ii,jj] = kk[-1]
#                plt.scatter(Wave[2][0][0].sum(axis=1)[ma1], Spec[2][ma1,ii,jj], marker='+', color = 'b')
#                plt.scatter(Wave[3][0][0].sum(axis=1)[ma2], L1_L2[ii,jj]*Spec[3][ma2,ii,jj], marker='+', color = 'r')
#                plt.scatter(np.abs(ww), func(np.abs(ww),kk[0:3]), color = 'k', marker='x')
#                plt.show()
    S_L = np.ones(Spec[0][0].shape)
    ma1 = Wave[1][0][0].sum(axis=1)>14.35
    ma2 = Wave[2][0][0].sum(axis=1)<14.9
    for ii in range(shape[0]):
        for jj in range(shape[1]):
            ww = np.array(list(-Wave[1][0][0].sum(axis=1)[ma1]) + list(Wave[2][0][0].sum(axis=1)[ma2]))
            ss = np.array(list(Spec[1][ma1,ii,jj]) + list(Spec[2][ma2,ii,jj]))
            if(np.isnan(ss.sum()) != True):
                kk = Cal_k([0,1,1,1], ww, ss)
                S_L[ii,jj] = kk[-1]
#                plt.scatter(Wave[1][0][0].sum(axis=1)[ma1], S_L[ii,jj]*Spec[1][ma1,ii,jj], marker='+', color = 'b')
#                plt.scatter(Wave[2][0][0].sum(axis=1)[ma2], Spec[2][ma2,ii,jj], marker='+', color = 'r')
#                plt.scatter(np.abs(ww), func(np.abs(ww),kk[0:3]), color = 'k', marker='x')
#                plt.show()
    Spectra = []
    Sigma = []
    Wavelength = []
    mask = Wave[0][0][0].sum(axis=1) < 7.5337
    Spectra = Spectra + list(Spec[0][mask]*S2_S1)        #*S_L)
    Sigma = Sigma + list(Spec_u[0][mask]*S2_S1)          #*S_L)
    Wavelength = Wavelength + list(Wave[0][0][0][mask])
    ma1 = Wave[1][0][0].sum(axis=1) >= 7.5337
    ma2 = Wave[1][0][0].sum(axis=1) < 14.2666
    mask = ma1*ma2
    Spectra = Spectra + list(Spec[1][mask])              #*S_L)
    Sigma = Sigma + list(Spec_u[1][mask])                #*S_L)
    Wavelength = Wavelength + list(Wave[1][0][0][mask])
    ma1 = Wave[2][0][0].sum(axis=1) >= 14.2666
    ma2 = Wave[2][0][0].sum(axis=1) < 20.5201
    mask = ma1*ma2
    Spectra = Spectra + list(Spec[2][mask])
    Sigma = Sigma + list(Spec_u[2][mask])
    Wavelength = Wavelength + list(Wave[2][0][0][mask])
    mask = Wave[3][0][0].sum(axis=1) >= 20.5201
    Spectra = Spectra + list(Spec[3][mask]*L1_L2)
    Sigma = Sigma + list(Spec_u[3][mask]*L1_L2)
    Wavelength = Wavelength + list(Wave[3][0][0][mask])
    Spec = np.array(Spectra)
    Spec_u = np.array(Sigma)
    Wave = Table([[Wavelength]],names = ['WAVELENGTH'], dtype=['float32'])
    
# Another scaling according to the difference between synthetic photometry and observed photometry
    IRAC = fits.open('Slice/{0}_{1}/{2}/{2}_{3}.fits'.format(dirsc, scale, name, 'IRAC4'))[0].data
#    delt = header['CDELT1']*header['CDELT2']
#    sr = delt*(2*np.pi/360)**2*10**9
#    IRAC = IRAC*sr
    IRACScale = np.ones(Spec[0].shape)
    ma = Wave[0][0].sum(axis=1) <= 14.2666
    photwave = Wave[0][0].sum(axis=1)[ma]
    for ii in range(shape[0]):
        for jj in range(shape[1]):
            if(np.isnan(Spec[ma,ii,jj].sum()) != True):
                SynPhot = Photometry(photwave, Spec[ma,ii,jj], 'Spitzer_IRAC4')
                ObsPhot = IRAC[ii,jj]
                IRACScale[ii,jj] = ObsPhot/SynPhot
    Spec[ma] = Spec[ma]*IRACScale
    Spec_u[ma] = Spec_u[ma]*IRACScale
    MIPS = fits.open('Slice/{0}_{1}/{2}/{2}_{3}.fits'.format(dirsc, scale, name, 'MIPS24'))[0].data
    MIPSScale = np.ones(Spec[0].shape)
    photwave = Wave[0][0].sum(axis=1)[~ma]
    for ii in range(shape[0]):
        for jj in range(shape[1]):
            if(np.isnan(Spec[~ma,ii,jj].sum()) != True):
                SynPhot = Photometry(photwave, Spec[~ma,ii,jj], 'Spitzer_MIPS_24')
                ObsPhot = MIPS[ii,jj]
                MIPSScale[ii,jj] = ObsPhot/SynPhot
    Spec[~ma] = Spec[~ma]*MIPSScale
    Spec_u[~ma] = Spec_u[~ma]*MIPSScale
    
    Spec_hdu = fits.PrimaryHDU(Spec, header = cube[0].header)
    Spec_hdu.header['NAXIS3'] = len(Wave[0][0])
    Wave_hdu = fits.BinTableHDU(Wave, name='WCS-TAB')
    Wave_hdu.header['TUNIT1'] = 'um'
    hdul = fits.HDUList([Spec_hdu,Wave_hdu])
    hdul.writeto('ReProjection/{0}_{1}/{2}/{2}_on2_{3}.fits'.format(dirsc,scale,name,O2), overwrite = True)
    Spec_u_hdu = fits.PrimaryHDU(Spec_u, header = cube[0].header)
    Spec_u_hdu.header['NAXIS3'] = len(Wave[0][0])
    hdul_u = fits.HDUList([Spec_u_hdu,Wave_hdu])
    hdul_u.writeto('ReProjection/{0}_{1}/{2}/{2}_on2_{3}_unc.fits'.format(dirsc,scale,name,O2), overwrite = True)
#################################################################################################################
# here need to take care about the difference of dimension order between numerical array and datacube
#        img = np.transpose(img, axes=(2, 1, 0))
#        img_u = np.transpose(img_u, axes=(2, 1, 0))
#################################################################################################################
    lightspd = (constants.c).to('cm/s').value
    wavelength = 10**-4*Wave[0][0].sum(axis=1)
    Spec_lbd = np.array(list(map(lambda x, y: x * y, Spec, lightspd/wavelength**2)))
    Spec_u_lbd = np.array(list(map(lambda x, y: x * y, Spec_u, lightspd/wavelength**2)))
    img = np.trapz(Spec_lbd, x = wavelength, axis = 0)
    img_u = np.trapz(Spec_u_lbd**2, x = wavelength, axis = 0)
    detwav = wavelength[-1] - wavelength[0]
    img = img/detwav
    img_u = np.sqrt(img_u/detwav)
#        img = np.transpose(img, axes=(1, 0))
#        img_u = np.transpose(img_u, axes=(1, 0))
    hdr = cube[0].header
#        hdr['NAXIS3'] = 1
    Iimg = fits.PrimaryHDU(img, hdr)
    Iimg_u = fits.PrimaryHDU(img_u, hdr)
    Iimg.writeto('ReProjection/{0}_{1}/{2}/{2}_on2_{3}_img.fits'.format(dirsc,scale,name,O2), overwrite = True)
    Iimg_u.writeto('ReProjection/{0}_{1}/{2}/{2}_on2_{3}_img_u.fits'.format(dirsc,scale,name,O2), overwrite = True)
    cube.close()
    cube_u.close()
    return()

def Photo_redis(name, band, dirsi, O2, dirsc, scale):
    import os
    import numpy as np
    from astropy import wcs
    from astropy.io import fits

    aabb = 'Projection_{0}/{1}_{2}/{3}/SL2_on2_{0}.fits'.format(O2, dirsc, scale, name)
    if(os.path.isfile(aabb)==False):
        print('No such file:{0}'.format(aabb))
        return ()
    else:
        flux_header = fits.open(aabb)[0].header
    apple = max(flux_header['NAXIS1'],flux_header['NAXIS2']) + 10
    crpixx, crpixy = flux_header['NAXIS1']//2, flux_header['NAXIS2']//2
    coor_mode = wcss(flux_header, [crpixx], [crpixy])
    crval = [[coor_mode[0][0]],[coor_mode[1][0]]]

    aabb = 'Projection_{0}/{1}_{2}/{3}/{3}_{4}.fits'.format(O2, dirsi, scale, name, band)
    if(os.path.isfile(aabb)):
        cube = fits.open(aabb)
    else:
        print('No such file:{0}'.format(aabb))
        return ()
    flux = np.transpose(cube['PRIMARY'].data, axes=(1, 0))
    flux_header = cube['PRIMARY'].header
    banana = np.full((apple, apple), np.nan)
    crpix = wcss3(flux_header, [crval[0][0]], [crval[1][0]])
    xx = int(round(crpix[0][0],0))
    yy = int(round(crpix[1][0],0))
    mm1 = max(0, xx-apple//2)
    mm2 = min(xx+apple-apple//2,flux_header['NAXIS1'])
    nn1 = max(0, yy-apple//2)
    nn2 = min(yy+apple-apple//2,flux_header['NAXIS2'])
    for ii in range(mm1,mm2):
        for jj in range(nn1,nn2):
            banana[apple//2 + ii-xx, apple//2 + jj-yy] = flux[ii, jj]
    CARD = [('NAXIS1',apple),('NAXIS2',apple),('CTYPE1',flux_header['CTYPE1']),('CTYPE2',flux_header['CTYPE2']),
            ('EQUINOX',flux_header['EQUINOX']),('CRVAL1',crval[0][0]),('CRVAL2',crval[1][0]),('CRPIX1',apple//2+1),
            ('CRPIX2',apple//2+1),('CDELT1',flux_header['CDELT1']),('CDELT2',flux_header['CDELT2']),
            ('CROTA2',flux_header['CROTA2']),('BUNIT','MJy/sr')]
    hdr = fits.Header(cards= CARD)
    banana = np.transpose(banana, axes=(1, 0))
    Img_hdu = fits.PrimaryHDU(banana, header=hdr)
    dirs = 'Slice/{0}_{1}/{2}'.format(dirsc, scale, name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    Img_hdu.writeto('Slice/{0}_{1}/{2}/{2}_{3}.fits'.format(dirsc, scale, name, band), overwrite = True)
    cube.close()
    return()

if __name__ == '__main__':

    gal_list = open('nuc_region_gal.txt','r')
    gal_name = gal_list.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_name = ['ngc5194']
    Band_list = ['IRAC4', 'MIPS24']
    Mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    for name in gal_name:
        for band in Band_list:
            try:
                Photo_redis(name, band, '{0}_convl'.format(band[:4]), 'LL2', 'Cube_convl', 10)
                print('{0} {1} image is OK!'.format(name, band))
            except:
                print('Something worry happens when dealing with {0} {1} image!'.format(name, band))
        try:
            Spectra_redis(name, 'LL2', 'Cube_convl', 10)
            print('{0} is OK!'.format(name))
        except:
            print('Something worry happens when dealing with {0} cube!'.format(name))
        for mode in Mode_list:
            try:
                Photo_redis(name, '{0}ind'.format(mode), 'Cube_convl', 'LL2', 'Cube_convl', 10)
                print('{0} {1} indicative image is OK!'.format(name, band))
            except:
                print('Something worry happens when dealing with {0} {1} indicative image!'.format(name, mode))

def ind(name, dirsc, scale):
    import os
    import numpy as np
    from astropy.io import fits

    mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    shape = fits.open('Slice/{0}_{1}/{2}/{2}_LL2ind.fits'.format(dirsc, scale, name))[0].data.shape
    ind = np.full(shape, 0)
    for mode in mode_list:
        banana = fits.open('Slice/{0}_{1}/{2}/{2}_{3}ind.fits'.format(dirsc, scale, name, mode))[0].data
        mask = banana==0
        banana[mask] = np.nan
        ind = ind + banana
    ind = ind/4
    mask = np.isnan(ind)
    ind[mask] = 0
    header = fits.open('Slice/{0}_{1}/{2}/{2}_{3}ind.fits'.format(dirsc, scale, name, mode))[0].header
    find = fits.PrimaryHDU(ind, header=header)
    dirs = 'Slice/{0}_{1}/{2}'.format(dirsc, scale, name)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    find.writeto('Slice/{0}_{1}/{2}/indice.fits'.format(dirsc, scale, name), overwrite=True)
    return ()

if __name__ == '__main__':

    gal_list = open('nuc_region_gal.txt','r')
    gal_name = gal_list.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_name = ['ngc5194']
    for name in gal_name:
        try:
#            ind(name, 'Cube', 25)
            ind(name, 'Cube_convl', 10)
            print('Indicative image of {0} is OK!'.format(name))
        except:
            print('Something worry happens when dealing with sample {0}!'.format(name))

def Optimal_binning(name, O2, SN, band, dirsc, scale):
    import numpy as np
    from astropy.io import fits
    from astropy.io import ascii
    from astropy import constants
    import matplotlib.pyplot as plt
    from astropy.table import Table

    np.seterr(divide='ignore', invalid='ignore')
    cube = fits.open('ReProjection/{0}_{1}/{2}/{2}_on2_{3}.fits'.format(dirsc,scale,name,O2))
    cube_u = fits.open('ReProjection/{0}_{1}/{2}/{2}_on2_{3}_unc.fits'.format(dirsc,scale,name,O2))
    Wave = cube['WCS-TAB'].data[0][0]
    Wave = Wave.sum(axis=1)
    mask1 = Wave > band[0]
    mask2 = Wave < band[1]
    mask = mask1*mask2
    Spec = cube['Primary'].data[mask]
    Spec_u = cube_u['Primary'].data[mask]
    Wave = 10**-4*Wave[mask]
    lightspd = (constants.c).to('cm/s').value
    Spec_lbd = np.array(list(map(lambda x, y: x * y, Spec, lightspd/Wave**2)))
    Spec_u_lbd = np.array(list(map(lambda x, y: x * y, Spec_u, lightspd/Wave**2)))
    Spec = np.transpose(Spec, axes=(2, 1, 0))
    Spec_u = np.transpose(Spec_u, axes=(2, 1, 0))
    Spec_lbd = np.transpose(Spec_lbd, axes=(2, 1, 0))
    Spec_u_lbd = np.transpose(Spec_u_lbd, axes=(2, 1, 0))
#    f2i = np.transpose(cube['Primary'].data, axes=(2, 1, 0))
#    f2i_u = np.transpose(cube_u['Primary'].data, axes=(2, 1, 0))
    S_sum = np.trapz(Spec_lbd, x = Wave, axis = 2)/(Wave[-1]-Wave[0])
#    hdu = fits.PrimaryHDU(S_sum)
#    hdu.writeto('test.fits', overwrite = True)
    S_sum[np.isnan(S_sum)] = 0
    Su_sum = np.sqrt(np.trapz(Spec_u_lbd**2, x = Wave, axis = 2)/(Wave[-1]-Wave[0]))
    Su_sum[np.isnan(Su_sum)] = 0
    size_xx, size_yy = S_sum.shape
    mask = np.full((size_xx, size_yy), False)
    for size_x in range(0,size_xx):
        for size_y in range(0,size_yy):
            leqzero = (Spec[size_x,size_y] == 0).sum()
            if(leqzero > len(Spec[size_x,size_y])/5):
                S_sum[size_x,size_y] = 0
                Su_sum[size_x,size_y] = 0
            leqzero = (Spec_u[size_x,size_y] == 0).sum()
            if(leqzero > len(Spec_u[size_x,size_y])/5):
                Su_sum[size_x,size_y] = 0
    coor = [[], []]
    bins = 0
    binnum = []
    inds = 0
    while(inds == 0):
        xxx, yyy = np.where(S_sum == np.max(S_sum))
        xx , yy = xxx[0], yyy[0]
        sig = S_sum[xx,yy]
        noi = Su_sum[xx,yy]
        if(sig/noi > SN):
            coor[0].append(xx)
            coor[1].append(yy)
            bins = bins + 1
            binnum.append(bins)
            S_sum[xx, yy] = 0
        else:
            coxy = [[[xx-1,xx+1],[yy,yy+1]],[[xx,xx+2],[yy,yy+1]],[[xx,xx+1],[yy-1,yy+1]],[[xx,xx+1],[yy,yy+2]],
                    [[xx-1,xx+1],[yy,yy+2]],[[xx,xx+2],[yy,yy+2]],[[xx-1,xx+1],[yy-1,yy+1]],[[xx,xx+2],[yy-1,yy+1]],
                    [[xx-1,xx+2],[yy,yy+2]],[[xx-1,xx+2],[yy-1,yy+1]],[[xx-1,xx+1],[yy-1,yy+2]],[[xx,xx+2],[yy-1,yy+2]],
                    [[xx-1,xx+2],[yy-1,yy+2]],
                    [[xx-2,xx+1],[yy,yy+2]],[[xx-1,xx+1],[yy,yy+3]],[[xx,xx+3],[yy,yy+2]],[[xx,xx+2],[yy,yy+3]],
                    [[xx-2,xx+1],[yy-1,yy+1]],[[xx-1,xx+1],[yy-2,yy+1]],[[xx,xx+3],[yy-1,yy+1]],[[xx,xx+2],[yy-2,yy+1]],
                    [[xx-2,xx+1],[yy,yy+3]],[[xx,xx+3],[yy,yy+3]],[[xx-2,xx+1],[yy-2,yy+1]],[[xx,xx+3],[yy-2,yy+1]]]
            chi = []
            nn = []
            s_n = []
            for cooo in coxy:
                binSp = Spec[cooo[0][0]:cooo[0][1],cooo[1][0]:cooo[1][1]]
                binSu = Spec_u[cooo[0][0]:cooo[0][1],cooo[1][0]:cooo[1][1]]
                imgSp = S_sum[cooo[0][0]:cooo[0][1],cooo[1][0]:cooo[1][1]]
                imgSu = Su_sum[cooo[0][0]:cooo[0][1],cooo[1][0]:cooo[1][1]]
                mask = (imgSp != 0)
                sig = imgSp[mask].sum()
                noi = np.sqrt((imgSu[mask]**2).sum())
                s_n.append(sig/noi)
                binSp = binSp[mask]
                binSu = binSu[mask]
                nn.append(len(binSp))
                cchi = 0
                for num in range(0,nn[-1]):
                    Smb = (binSp[num]*Spec[xx,yy]/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()/(Spec[xx,yy]**2/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()
                    cchi = cchi + ((binSp[num]-Smb*Spec[xx,yy])**2/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()
                chi.append(cchi)
            chi = np.array(chi)
            nn = np.array(nn)
            s_n = np.array(s_n)
            sn2bin = np.where(s_n > SN)
            if(len(sn2bin[0]) > 0):
                nn2bin = np.where(nn == np.min(nn[sn2bin]))
                co2bin = np.where(chi == np.min(chi[nn2bin]))
                bins = bins + 1
                for cx in range(coxy[co2bin[0][0]][0][0],coxy[co2bin[0][0]][0][1]):
                    for cy in range(coxy[co2bin[0][0]][1][0],coxy[co2bin[0][0]][1][1]):
                        if(S_sum[cx, cy] != 0):
                            coor[0].append(cx)
                            coor[1].append(cy)
                            binnum.append(bins)
                            S_sum[cx, cy] = 0
            else:
                ind = 1
                ii = 2
                while(ind == 1):
                    coxy = [[[xx-ii,xx+1],[yy,yy+ii]],[[xx-ii+1,xx+1],[yy,yy+ii+1]],[[xx,xx+ii+1],[yy,yy+ii]],[[xx,xx+ii],[yy,yy+ii+1]],
                            [[xx-ii,xx+1],[yy-ii+1,yy+1]],[[xx-ii+1,xx+1],[yy-ii,yy+1]],[[xx,xx+ii+1],[yy-ii+1,yy+1]],[[xx,xx+ii],[yy-ii,yy+1]],
                            [[xx-ii,xx+1],[yy,yy+ii+1]],[[xx,xx+ii+1],[yy,yy+ii+1]],[[xx-ii,xx+1],[yy-ii,yy+1]],[[xx,xx+ii+1],[yy-ii,yy+1]],
                            [[xx-ii,xx+ii+1],[yy,yy+ii+1]],[[xx-ii,xx+ii+1],[yy-ii,yy+1]],[[xx-ii,xx+1],[yy-ii,yy+ii+1]],[[xx,xx+ii+1],[yy-ii,yy+ii+1]],
                            [[xx-ii-1,xx+1],[yy,yy+ii+1]],[[xx-ii,xx+1],[yy,yy+ii+2]],[[xx,xx+ii+2],[yy,yy+ii+1]],[[xx,xx+ii+1],[yy,yy+ii+2]],
                            [[xx-ii-1,xx+1],[yy-ii,yy+1]],[[xx-ii,xx+1],[yy-ii-1,yy+1]],[[xx,xx+ii+2],[yy-ii,yy+1]],[[xx,xx+ii+1],[yy-ii-1,yy+1]]]
                    if(ii%2 == 1):
                        coxy.append([[xx-ii,xx+ii+1],[yy-ii,yy+ii+1]])
                    if(ii%2 == 0):
                        coxy.append([[xx-ii,xx+ii],[yy-ii+1,yy+ii+1]])
                        coxy.append([[xx-ii+1,xx+ii+1],[yy-ii+1,yy+ii+1]])
                        coxy.append([[xx-ii,xx+ii],[yy-ii,yy+ii]])
                        coxy.append([[xx-ii+1,xx+ii+1],[yy-ii,yy+ii]])
                    chi = []
                    nn = []
                    s_n = []
                    for cooo in coxy:
                        width = len(S_sum)
                        xmin = max(0, cooo[0][0])
                        xmax = min(width, cooo[0][1])
                        ymin = max(0, cooo[1][0])
                        ymax = min(width, cooo[1][1])
                        binSp = Spec[xmin:xmax,ymin:ymax]
                        binSu = Spec_u[xmin:xmax,ymin:ymax]
                        imgSp = S_sum[xmin:xmax,ymin:ymax]
                        imgSu = Su_sum[xmin:xmax,ymin:ymax]
                        mask = (imgSp != 0)
                        sig = imgSp[mask].sum()
                        noi = np.sqrt((imgSu[mask]**2).sum())
                        s_n.append(sig/noi)
                        binSp = binSp[mask]
                        binSu = binSu[mask]
                        nn.append(len(binSp))
                        cchi = 0
                        for num in range(0,nn[-1]):
                            Smb = (binSp[num]*Spec[xx,yy]/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()/(Spec[xx,yy]**2/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()
                            cchi = cchi + ((binSp[num]-Smb*Spec[xx,yy])**2/(binSu[num]**2+Spec_u[xx,yy]**2)).sum()
                        chi.append(cchi)
                    chi = np.array(chi)
                    nn = np.array(nn)
                    s_n = np.array(s_n)
                    sn2bin = np.where(s_n > SN)
                    if(len(sn2bin[0]) > 0):
                        nn2bin = np.where(nn == np.min(nn[sn2bin]))
                        co2bin = np.where(chi == np.min(chi[nn2bin]))
                        bins = bins + 1
                        width = len(S_sum)
                        xmin = max(0, coxy[co2bin[0][0]][0][0])
                        xmax = min(width, coxy[co2bin[0][0]][0][1])
                        ymin = max(0, coxy[co2bin[0][0]][1][0])
                        ymax = min(width, coxy[co2bin[0][0]][1][1])
                        for cx in range(xmin, xmax):
                            for cy in range(ymin, ymax):
                                if(S_sum[cx, cy] != 0):
                                    coor[0].append(cx)
                                    coor[1].append(cy)
                                    binnum.append(bins)
                                    S_sum[cx, cy] = 0
                        ind = 0
                    else:
                        S_cri = S_sum*1
                        width = len(S_sum)
                        xmin = max(0, xx-ii)
                        xmax = min(width, xx+ii+1)
                        ymin = max(0, yy-ii)
                        ymax = min(width, yy+ii+1)
                        for cx in range(xmin, xmax):
                            for cy in range(ymin, ymax):
                                S_cri[cx, cy] = 0
                        if(S_cri.sum() == 0):
                            bins = bins + 1
                            for cx in range(xmin, xmax):
                                for cy in range(ymin, ymax):
                                    if(S_sum[cx, cy] != 0):
                                        coor[0].append(cx)
                                        coor[1].append(cy)
                                        binnum.append(bins)
                                        S_sum[cx, cy] = 0
                            ind = 0
                        else:
                            ii = ii+1
        if(S_sum.sum() == 0):
            inds = 1
    coor = np.array(coor)
    binnum = np.array(binnum)
    pixelsize = 1
    plt.figure(figsize=[6,6])
    xmin, xmax = np.min(coor[0]), np.max(coor[0])
    ymin, ymax = np.min(coor[1]), np.max(coor[1])
    width = len(S_sum)
    img = np.full((width, width), np.nan)
    for i in range(0, len(binnum)):
        img[coor[0][i]][coor[1][i]] = binnum[i]
    plt.imshow(np.rot90(img), interpolation='nearest', cmap='Spectral',
               extent=[0-pixelsize/2, width-pixelsize/2, 0-pixelsize/2, width-pixelsize/2])
    plt.xlabel('R (pixel)',fontsize=12)
    plt.ylabel('R (pixel)',fontsize=12)
    plt.title('Map of Optimal bins',fontsize=12)
    plt.show()
    return(coor[0], coor[1], binnum)

def exe_bin(name, O2, SN, band, dirsc, scale):
    path = 'Slice/{0}_{1}/{2}/'.format(dirsc, scale, name)
    if(os.path.exists(path) == False):
        os.makedirs(path)
    print(name)
    x, y, binnum = Optimal_binning(name, O2, SN, band, dirsc, scale)
    ascii.write([x, y, binnum], 'ReProjection/{0}_{1}/{2}/{2}_binning_output.txt'.format(dirsc,scale,name), 
                names=['X', 'Y', 'binNum'], format = 'ipac', overwrite = True)
    ascii.write([x, y, binnum], path+'{1}_binning_output.txt'.format(path, name), 
                names=['X', 'Y', 'binNum'], format = 'ipac', overwrite = True)
    return ()

import os
from astropy.io import ascii

gal_list = open('nuc_region_gal.txt','r')
gal_name = gal_list.readlines()
gal_name = list(map(lambda x: x.strip(), gal_name))
gal_name = ['ngc5194']
for name in gal_name:
    try:
        exe_bin(name, 'LL2', 0, [5,20], 'Cube_convl', 10)
    except:
        print('Something worry happens!')

def Bin_Spec(name, dirsc, scale, O2):

    import os
    import numpy as np
    from astropy.io import fits
    from astropy.io import ascii
    from astropy.table import Table
    import matplotlib.pyplot as plt
    
    binR = ascii.read('ReProjection/{0}_{1}/{2}/{2}_binning_output.txt'.format(dirsc, scale, name), format = 'ipac')
    cube = fits.open('ReProjection/{0}_{1}/{2}/{2}_on2_{3}.fits'.format(dirsc, scale, name, O2))
    cube_u = fits.open('ReProjection/{0}_{1}/{2}/{2}_on2_{3}_unc.fits'.format(dirsc, scale, name, O2))
    flux = np.transpose(cube['PRIMARY'].data, axes=(2, 1, 0))
    flux_u = np.transpose(cube_u['PRIMARY'].data, axes=(2, 1, 0))

    x, y, binN = binR['X'], binR['Y'], binR['binNum']
    Wave = cube['WCS-TAB'].data[0][0]
    Wave = Wave.sum(axis = 1)
    lee = len(Wave)
    num = np.max(binN)
    for ii in range (1,num+1):
        mask = (binN == ii)
        binRm = binR[mask]
        Spec = np.zeros(lee)
        Spec_u = np.zeros(lee)
        nn = len(binRm)
        for each in binRm:
            Spec = Spec + flux[each['X'], each['Y']]
            Spec_u = Spec_u + flux_u[each['X'], each['Y']]**2
        delt = cube[0].header['CDELT1']*cube[0].header['CDELT2']
        sr = delt*(2*np.pi/360)**2*10**9
        Spec = Spec*sr
        band = np.zeros(lee, dtype=int)
        Spectra = Table([Wave, Spec, Spec_u, band], names = ('wavelength','flux','sigma','band'))
        Spectra['wavelength'].format = '.8f'
        Spectra['flux'].format = '.8f'
        Spectra['sigma'].format = '.9f'
#                    Spec['band'].format = '6d'
        dirs = 'Spectra/{0}_{1}/{2}/'.format(dirsc, scale, name)
        if(os.path.exists(dirs) != True):
            os.makedirs(dirs)
        Spectra.write('Spectra/{0}_{1}/{2}/{2}_bin_{3}.tbl'.format(dirsc, scale, name, str(ii)),format='ascii.ipac',overwrite=True)
    return ()

def Bin_Phot(name, band, dirsc, scale):

    import numpy as np
    from astropy.io import fits
    from astropy.io import ascii
    from astropy.table import Table
    
    path = 'Slice/{0}_{1}/{2}/'.format(dirsc, scale, name)
    binR = ascii.read(path + '{0}_binning_output.txt'.format(name), format = 'ipac')
    cube = fits.open(path + '{0}_{1}.fits'.format(name, band))
    flux = np.transpose(cube[0].data, axes=(1, 0))
    cover = fits.open('Slice/{0}_{1}/{2}/indice.fits'.format(dirsc, scale, name))[0].data
    cover = np.transpose(cover, axes=(1, 0))

    x, y, binN = binR['X'], binR['Y'], binR['binNum']
    num = np.max(binN)
    ind = []
    photo = []
    area = []
    delt = cube[0].header['CDELT1']*cube[0].header['CDELT2']
    sr = delt*(2*np.pi/360)**2*10**9
    for ii in range (1,num+1):
        mask = (binN == ii)
        binRm = binR[mask]
        Img = 0
        ca = 0
        for each in binRm:
            Img = Img + flux[each['X'], each['Y']]*sr
            ca = ca + cover[each['X'], each['Y']]
        ind.append(ii)
        photo.append(Img)
        area.append(ca/len(binRm))
    Photom = Table([ind, photo, area], names = ('binN','flux', 'area'))
    Photom['flux'].format = '.8f'
    Photom['area'].format = '.3f'
    Photom.write(path + '{0}_{1}_photometry.tbl'.format(name, band), format='ascii.ipac', overwrite=True)
    return ()

if __name__ == '__main__':

    gal_list = open('nuc_region_gal.txt','r')
    gal_name = gal_list.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_name = ['ngc5194']
    Band_list = ['IRAC4', 'MIPS24']
    for name in gal_name:
        try:
            Bin_Spec(name, 'Cube_convl', 10, 'LL2')
            print('Spectrometry of {0} is OK!'.format(name))
        except:
            print('Something worry happens when dealing with sample {0}!'.format(name))
        for band in Band_list:
            try:
                Bin_Phot(name, band, 'Cube_convl', 10)
                print('Photometry of {0} {1} image is OK!'.format(name, band))
            except:
                print('Something worry happens when dealing with sample {0}!'.format(name))

def Plot_bSpec(name, dirsc, scale):

    import numpy as np
    from astropy.io import fits
    from astropy.io import ascii
    import matplotlib.pyplot as plt
    
    binR = ascii.read('Slice/{0}_{1}/{2}/{2}_binning_output.txt'.format(dirsc, scale, name), format = 'ipac')
    cover = fits.open('Slice/{0}_{1}/{2}/indice.fits'.format(dirsc, scale, name))[0].data
    cover = np.transpose(cover, axes=(1, 0))
    
    x, y, binN = binR['X'], binR['Y'], binR['binNum']
    num = np.max(binN)
    plt.figure(figsize=[15,8])
    colors = plt.cm.Spectral(np.linspace(0,1,num+1))
    for ii in range(1,num+1):
        if(cover[x[ii-1],y[ii-1]] == 1):
            Spec = ascii.read('Spectra/{0}_{1}/{2}/{2}_bin_{3}.tbl'.format(dirsc, scale, name, ii),format='ipac')
            plt.plot(Spec['wavelength'], Spec['flux'], c = colors[ii])
    plt.xscale('log')
    plt.xlabel('$\lambda\ (um)$', fontsize=18)
    plt.ylabel('$Flux\ (mJy)$', fontsize=18)
    plt.title(name, fontsize=18)
    plt.show()
    return

if __name__ == '__main__':

    gal_file = open('nuc_region_gal.txt')
    gal_name = gal_file.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_name = ['ngc5194']
    for name in gal_name:
        try:
            Plot_bSpec(name, 'Cube_convl', 10)
        except:
            print('Something worry happens when dealing with sample {0}!'.format(name))

import os
import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table

def Spc_cal(name, calband, dirsc, scale, dirsi):

    np.seterr(divide='ignore', invalid='ignore')
    path = 'ForSFR/Spc_{0}/{1}/'.format(dirsi,name)
    if(os.path.exists(path) == False):
        os.makedirs(path)
    dirs = 'Spectra/{0}_{1}/{2}'.format(dirsc, scale, name)
    SpcList = []
    for f in os.listdir(dirs):
        if f.endswith('.tbl'):
            SpcList.append(f)
    Phot = ascii.read('Slice/{0}_{1}/{2}/{2}_{3}_photometry.tbl'.format(dirsc, scale, name, calband[0]), format = 'ipac')
    ratio = []
    for spc in SpcList:
        Specc = ascii.read('Spectra/{0}_{1}/{2}/{3}'.format(dirsc, scale, name, spc), format = 'ipac')
        Wave = Specc['wavelength']
        Spec = Specc['flux']
        sigma = Specc['sigma']
        band_ind = Specc['band']
        ii = spc.split('_')[-1].split('.')[0]
        SynPhot = Photometry(Wave, Spec, 'Spitzer_{0}'.format(calband[1]))
#            binN, flux = Phot['binN'], Phot['flux']
        ma = Phot['binN']==int(ii)
        Phof = Phot[ma]['flux']
        Spec = Spec*Phof[0]/SynPhot
        sigma = sigma*Phof[0]/SynPhot
        if(Phot[ma]['area']==1):
            ratio.append(Phof[0]/SynPhot)
            mask = Spec>0
            Spectra = Table([Wave[mask]*u.um, Spec[mask]*u.mJy, sigma[mask]*u.mJy, band_ind[mask]], names = ('wavelength','flux','sigma','band'))
            Spectra['wavelength'].format = '.8f'
            Spectra['flux'].format = '.8f'
            Spectra['sigma'].format = '.9f'
            Spectra.write(path + spc, format='ascii.ipac', overwrite=True)
    ratio = np.array(ratio)
    return('{0} {1} {2} {3}'.format(ratio.max(), ratio.min(), ratio.mean(), ratio.std()))

if __name__ == '__main__':

    gal_list = open('nuc_region_gal.txt','r')
    gal_name = gal_list.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    gal_name = ['ngc5194']
    Band_list = [['IRAC4', 'IRAC4'], ['MIPS24', 'MIPS_24']]
    SynPho = open('SynOPho.txt','w+')
    for name in gal_name:
        for band in Band_list:
            try:
                print('Calibration of {0} spectra with {1} photometry'.format(name, band[0]))
                ratio = Spc_cal(name, band, 'Cube_convl', 10, '{0}'.format(band[0][:4]))
                SynPho.write('{0}\n'.format(ratio))
                print('Calibration is OK!'.format(name, band[0]))
            except:
                print('Something wrong happens when dealing with sample {0}!'.format(name))
    SynPho.close()

def Circularization(path, mode, name, recover_edges, path1, path2):
    import numpy as np 
    from astropy.io import fits
    from scipy.ndimage import rotate

    nm1 = name.split('.')[0]
    nm2 = name.split('.')[1]
    Img  = fits.open(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, path1))
    image = Img[0].data
    header = Img[0].header
    size_im = image.shape

    if (recover_edges == 1):
        sz = max(size_im)
        pixels_added_side = int((sz * np.sqrt(2.0) - sz)/ 2.0) + 1 
        pixels_added      = 2 * pixels_added_side
        new_image_size_x = size_im[0] + pixels_added
        new_image_size_y = size_im[1] + pixels_added

        new_image = np.zeros((new_image_size_x,new_image_size_y))
        image_one = np.zeros((new_image_size_x,new_image_size_y))

        new_image [pixels_added_side:pixels_added_side+size_im[0],pixels_added_side:pixels_added_side+size_im[1]] = image
        wh_data  = new_image != 0
        image_one[wh_data] = 1.0

        image = new_image
        header['CRPIX1']=(new_image_size_x-1)/2
        header['CRPIX2']=(new_image_size_y-1)/2
        header['NAXIS1']=new_image_size_x
        header['NAXIS2']=new_image_size_y

    size_im = len(image)
    center_pixel = int((size_im - 1) / 2)

    original_image = image

#First we mask the pixel outside the major circle contained into the image,
#since they are far away from the center an they do not contain useful information...

    xdist, ydist = np.mgrid[0:size_im, 0:size_im] - center_pixel
    distance_sq = xdist**2 + ydist**2 
############################################################################################################
    not_useful  = distance_sq > center_pixel**2
    image[not_useful] = 0
############################################################################################################

#now we rotate 180 deg te image and add it to itself, we rotate the composite 90 deg 
#and add it to itself, then 45 deg... in 10 ierations we have 16384 rotations...
    image_temp = image
    for iteration in range(14,0,-1):
        angle = float(360.0 / (2.0 ** iteration))
        rotate(image, angle, axes=(1, 0), reshape=False, output=image_temp)
        image = (image + image_temp )/2.0

    if (recover_edges == 1):
        image_temp = image_one
        for iteration in range(14,0,-1):
            angle = float(360.0 / (2.0 ** iteration))
            rotate(image_one, angle, axes=(1, 0), reshape=False, output=image_temp)
            image_one = (image_one + image_temp )/2.0

        wh = image_one > 0.1
        corrector = image_one * 0
        corrector[wh] = 1 / image_one[wh]
        image = image * corrector

    image[not_useful] = 0.0
    image /= image.sum()
    wh = np.abs(image) < 1.0e-10
    image[wh] = 0.0

    hdu = fits.PrimaryHDU(image,header)
    hdu.writeto(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, path2), overwrite = True)   

    difference = original_image - image
    difference[wh] = 0.0
    non_axi = np.abs(difference).sum() / np.abs(image).sum()
#print('------------------------------------------------------------------------------')
#print('The anisotropy parameter g of this circlularization is {0}'.format(str(round(non_axi,4))))
#print('The test was succesfully circularized!')
#print('------------------------------------------------------------------------------')
    return(non_axi)
if __name__ == '__main__':
    mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    #non_axi = []
    path = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/PSFsam/'
    for mode in mode_list:
        wavef = open('/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/Wave/{0}_wavelength.txt'.format(mode))
        wavel = wavef.readlines()
        for one in wavel:
            try:
                n_axi = Circularization(path, mode, one.strip(), 0, 'PSFsamp', 'PSFcirc')
#               non_axi.append(n_axi)
            except:
                print('Something worry happens when dealing with slice at wavelength {0}!'.format(one.strip()))


def PSF_matching(path, mode, name, FWHM, path1, path2):
    import numpy as np
    from astropy.io import fits
    from numpy.fft import fft2, fftshift
    from astropy.modeling.models import Gaussian2D
    from photutils import create_matching_kernel, SplitCosineBellWindow
    import matplotlib.pyplot as plt

    nm1 = name.split('.')[0]
    nm2 = name.split('.')[1]
    Img  = fits.open(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, path1))
    psf = Img[0].data
    header = Img[0].header
    size_im = psf.shape
    scale = round(header['PIXSCALX'], 3)
    sigma = FWHM/(2*np.sqrt(2*np.log(2)))
    y, x = np.mgrid[0:round(size_im[0]*scale,3):scale, 0:round(size_im[0]*scale,3):scale] - round(scale*(size_im[0]-1)/2, 3)
    gm = Gaussian2D(100, 0, 0, sigma, sigma)
    gs = gm(x, y)
    gs /= gs.sum()
#    size_im = len(x)
#    center_pixel = int((size_im - 1) / 2)
#    xdist, ydist = np.mgrid[0:size_im, 0:size_im] - center_pixel
#    distance_sq = xdist**2 + ydist**2 
#    not_useful  = distance_sq > center_pixel**2
#    gs[not_useful] = 0
    
    ftg1 = np.abs(np.real(fftshift(fft2(psf))))
    mask = ftg1 < 5*10**-3*np.max(ftg1)
    dist = np.sqrt(x**2 + y**2)
    kha = np.min(dist[mask])/scale
    alpha = 0.3*kha/((len(psf)-1)/2)
    beta = 0.7*kha/((len(psf)-1)/2)
    window = SplitCosineBellWindow(alpha, beta)
    kernel = create_matching_kernel(psf, gs, window=window)
    hdu = fits.PrimaryHDU(kernel, header)
    hdu.writeto(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, path2), overwrite = True)
    return(kha)

if __name__ == '__main__':
    mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    path = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/PSFsam/'
    #kha = []
    #non_axii = []
    for mode in mode_list:
        wavef = open('/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/Wave/{0}_wavelength.txt'.format(mode))
        wavel = wavef.readlines()
        for one in wavel:
            try:
                khaa = PSF_matching(path, mode, one.strip(), 10, 'PSFcirc', 'ker2ten')
                n_axii = Circularization(path, mode, one.strip(), 0, 'ker2ten', 'kercirc')
#               kha.append(khaa)
#               non_axii.append(n_axii)
            except:
                print('Something worry happens when dealing with slice at wavelength {0}!'.format(one.strip()))


def CircKer(image, recover_edges):
    import numpy as np 
    from astropy.io import fits
    from scipy.ndimage import rotate

    image = image
    size_im = image.shape

    if (recover_edges == 1):
        sz = max(size_im)
        pixels_added_side = int((sz * np.sqrt(2.0) - sz)/ 2.0) + 1 
        pixels_added      = 2 * pixels_added_side
        new_image_size_x = size_im[0] + pixels_added
        new_image_size_y = size_im[1] + pixels_added

        new_image = np.zeros((new_image_size_x,new_image_size_y))
        image_one = np.zeros((new_image_size_x,new_image_size_y))

        new_image [pixels_added_side:pixels_added_side+size_im[0],pixels_added_side:pixels_added_side+size_im[1]] = image
        wh_data  = new_image != 0
        image_one[wh_data] = 1.0

        image = new_image

    size_im = len(image)
    center_pixel = int((size_im - 1) / 2)

    original_image = image

#First we mask the pixel outside the major circle contained into the image,
#since they are far away from the center an they do not contain useful information...

    xdist, ydist = np.mgrid[0:size_im, 0:size_im] - center_pixel
    distance_sq = xdist**2 + ydist**2 
############################################################################################################
    not_useful  = distance_sq > center_pixel**2
    image[not_useful] = 0
############################################################################################################

#now we rotate 180 deg te image and add it to itself, we rotate the composite 90 deg 
#and add it to itself, then 45 deg... in 10 ierations we have 16384 rotations...
    image_temp = image
    for iteration in range(14,0,-1):
        angle = float(360.0 / (2.0 ** iteration))
        rotate(image, angle, axes=(1, 0), reshape=False, output=image_temp)
        image = (image + image_temp )/2.0

    if (recover_edges == 1):
        image_temp = image_one
        for iteration in range(14,0,-1):
            angle = float(360.0 / (2.0 ** iteration))
            rotate(image_one, angle, axes=(1, 0), reshape=False, output=image_temp)
            image_one = (image_one + image_temp )/2.0

        wh = image_one > 0.1
        corrector = image_one * 0
        corrector[wh] = 1 / image_one[wh]
        image = image * corrector

    image[not_useful] = 0.0
    image /= image.sum()
    wh = np.abs(image) < 1.0e-10
    image[wh] = 0.0
    return(image)

def Ker_resize(path, samfile, tagtfits, mode, pathker, pathkerf):
    import numpy as np
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from astropy.convolution import convolve_fft
    from photutils.psf.matching import resize_psf

    wavef = open(samfile)
    wavel = wavef.readlines()
    img_fits = fits.open(tagtfits)
    for ii in range(0, len(wavel)):
        wave = wavel[ii].strip()
        nm1 = wave.split('.')[0]
        nm2 = wave.split('.')[1]
        kerf = fits.open(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, pathker))
        ker = kerf[0].data
        scalei = round(kerf[0].header['PIXSCALX'],3)
        scaleii = round(np.abs(img_fits[0].header['CDELT1'])*3600, 3)
        if(mode[0:2] == 'SL'):
            ker = ker[4:-4,4:-4]
# just make sure not to invoke the warning "the returned array has changed.", UserWarning
        rem = round((len(ker)*(scalei/scaleii))%1,2)
        while(rem >= 0.5):
            ker = ker[1:-1,1:-1]
            rem = round((len(ker)*(scalei/scaleii))%1,2)
        ker = resize_psf(ker, scalei, scaleii, order=3)
        ker = CircKer(ker, 0)
        hdu = fits.PrimaryHDU(ker)
        hdu.writeto(path + '{0}_{3}/{1}_{2}.fits'.format(mode[0:2], nm1, nm2, pathkerf), overwrite = True)
        img_fits.close()
# The dimensions do not have to be odd in all directions, unlike in the non-fft convolve function
#        size = len(ker)
#        if ((size % 2) == 0):
#            size = size + 1
#            new_ker = np.zeros((size,size))
#            new_ker[0:size-1,0:size-1] = ker
#            ker = new_ker
    return ()

if __name__ == '__main__':
    mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    path = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/PSFsam/'
    target = 'ngc5194'
    for mode in mode_list:
        samfile = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/Wave/{0}_wavelength.txt'.format(mode)
        tagtfits = '/Users/zhangl/desktop/SINGS/Cube/{0}/{0}_DR5_{1}_cube.fits'.format(target, mode)
        try:
            Ker_resize(path, samfile, tagtfits, mode, 'kercirc', 'ker4conv')
        except:
            print('Something worry happens when dealing with {0} data!'.format(mode))


def Convolve(path, samfile, tagtfits, target, mode, dtype):
    import os
    import numpy as np
    from astropy.io import fits
    from astropy.convolution import convolve_fft

    wavef = open(samfile)
    wavel = wavef.readlines()
    wavel = list(map(lambda x: x.strip(), wavel))
    if(os.path.isfile(tagtfits)==True):
        img_fits = fits.open(tagtfits)
    else:
        print('No such path or file:{0}'.format(tagtfits))
        return()
    img = img_fits[0].data
    img_conv = []
    for ii in range(0, len(wavel)):
        slice_Img = img[ii]
        wave = wavel[ii]
        nm1 = wave.split('.')[0]
        nm2 = wave.split('.')[1]
        kerf = fits.open(path + '{0}_ker4conv/{1}_{2}.fits'.format(mode[0:2], nm1, nm2))
        ker = kerf[0].data
        img_cov = convolve_fft(slice_Img, ker, boundary = 'fill', allow_huge=True, preserve_nan=True)
        img_conv.append(list(img_cov))
    hdu = fits.PrimaryHDU(img_conv, header=img_fits[0].header)
    dirs = '/Users/zhangl/desktop/SINGS/Cube_convl/{0}/'.format(target)
    if(os.path.exists(dirs) != True):
        os.makedirs(dirs)
    hdu.writeto(dirs + '{0}_{1}_{2}.fits'.format(target, mode, dtype), overwrite = True)
    img_fits.close()
    print('Covolution of {0} {1} image is OK!'.format(target, mode))
    return()

if __name__ == '__main__':
    mode_list = ['SL2', 'SL1', 'LL2', 'LL1']
    path = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/PSFsam/'
    gal_file = open('nuc_region_gal.txt')
    gal_name = gal_file.readlines()
    gal_name = list(map(lambda x: x.strip(), gal_name))
    for target in gal_name:
        for mode in mode_list:
            samfile = '/Users/zhangl/desktop/SINGS/Kernel/G_Aniano_Kernels/Wave/{0}_wavelength.txt'.format(mode)
            tagtfits1 = '/Users/zhangl/desktop/SINGS/Cube/{0}/{0}_DR5_{1}_cube.fits'.format(target, mode)
            tagtfits2 = '/Users/zhangl/desktop/SINGS/Cube/{0}/{0}_DR5_{1}_cube_unc.fits'.format(target, mode)
            try:
                Convolve(path, samfile, tagtfits1, target, mode, 'convolved')
                Convolve(path, samfile, tagtfits2, target, mode, 'unc_convolved')
            except:
                print('Something worry happens when dealing with {0} {1} data!'.format(target, mode))

            dirs = '/Users/zhangl/desktop/SINGS/Cube_convl/{0}/'.format(target)
            tagtfits = '/Users/zhangl/desktop/SINGS/Cube/{0}/{0}_DR5_{1}_cube.fits'.format(target, mode)
            try:
                wave = fits.open(tagtfits)[1]
                ctrg1 = dirs + '{0}_{1}_{2}.fits'.format(target, mode, 'convolved')
                ctrg2 = dirs + '{0}_{1}_{2}.fits'.format(target, mode, 'unc_convolved')
                spec1 = fits.open(ctrg1)[0]
                spec2 = fits.open(ctrg2)[0]
                hdul1 = fits.HDUList([spec1,wave])
                hdul2 = fits.HDUList([spec2,wave])
                hdul1.writeto(dirs + '{0}_{1}_{2}.fits'.format(target, mode, 'convolved'), overwrite=True)
                hdul2.writeto(dirs + '{0}_{1}_{2}.fits'.format(target, mode, 'unc_convolved'), overwrite=True)
            except:
                print('Something worry happens when dealing with {0} {1} data!'.format(target, mode))

#This config file is for the objects without silicates emission feature.
#
import numpy as np
from astropy.table import Table
from collections import OrderedDict

################################################################################
#                                    Data                                      #
################################################################################
targname = "NGC5194_sed_1"
redshift = 0.00154
distance = 7.62
#sedFile = "NGC5194_sp.tbl"  
sedFile  = "test/NGC5194_SED/mock/ngc5194_sed_1.tbl" 
dataDict = {
        "phtName": "Phot",
        "spcName": "IRS",
        "bandList_use": ['WISE_w1','Spitzer_IRAC1','Spitzer_IRAC2','WISE_w2','Spitzer_IRAC3'],
        "bandList_ignore":[],
        "frame": "obs",
}

################################################################################
#                                   Model                                      #
################################################################################
#Spec = Table.read(sedFile, format='ascii.ipac')
#wave = Spec['wavelength']
#wavemod = 10**np.linspace(0.5, 2.0, 1000)
#mask = Spec['band']=='0'
#mask1 = wavemod < np.min(wave[mask])
#mask2 = wavemod > np.max(wave[mask])
#waveModel = np.concatenate((wavemod[mask1], wave, wavemod[mask2]))
waveModel = 10**np.linspace(0, 2.0, 1000)
parAddDict_all = {
    "kappa0": 1.92, #MBB
    "lambda0": 350 #MBB
}
modelDict = OrderedDict(
    (
        ("PAH", {
                "function": "pah",
                "logLpah": {
                    "value": 28.53,
                    "range": [25.5, 31.5],
                    "type": "c",
                    "vary": True,
                    "latex": r"$\mathrm{log}\,L_\mathrm{PAH}$",
                }
            }
        ),
        ("Old Star", {
                "function": "BlackBody",
                "logOmega":{
                    "value": -20,
                    "range": [-25, 15],
                    "type": "c",
                    "vary": True,
                    "latex": r"$\mathrm{log}\,\Omega$",
                },
                "T":{
                    "value": 5000.0,
                    "range": [4500.0, 5500.0],
                    "type": "c",
                    "vary": False,
                    "latex": r"$T$",
                },
            }
        ),
        ("MBB1", {
                "function": "Modified_BlackBody",
                "logM":{
                    "value": 3.296,
                    "range": [-10, 10],
                    "type": "c",
                    "vary": True,
                    "latex": r"$\mathrm{log}\,M_\mathrm{d}$",
                },
                "beta":{
                    "value": 2.00,
                    "range": [1.0, 2.5],
                    "type": "c",
                    "vary": False, #True, #
                    "latex": r"$\beta$",
                },
                "T":{
                    "value": 41.28,
                    "range": [37.0, 53.0],
                    "type": "c",
                    "vary": True,
                    "latex": r"$T$",
                },
            }
        ),
        ("MBB2", {
                "function": "Modified_BlackBody",
                "logM":{
                    "value": -0.046,
                    "range": [-10, 10],
                    "type": "c",
                    "vary": True,
                    "latex": r"$\mathrm{log}\,M_\mathrm{d}$",
                },
                "beta":{
                    "value": 2.00,
                    "range": [1.0, 2.5],
                    "type": "c",
                    "vary": False, #True, #
                    "latex": r"$\beta$",
                },
                "T":{
                    "value": 109,
                    "range": [45.0, 155.0],
                    "type": "c",
                    "vary": True,
                    "latex": r"$T$",
                },
            }
        ),
        ("MBB3", {
                "function": "Modified_BlackBody",
                "logM":{
                    "value": -3.206,
                    "range": [-10, 10],
                    "type": "c",
                    "vary": True,
                    "latex": r"$\mathrm{log}\,M_\mathrm{d}$",
                },
                "beta":{
                    "value": 2.00,
                    "range": [1.0, 2.5],
                    "type": "c",
                    "vary": False, #True, #
                    "latex": r"$\beta$",
                },
                "T":{
                    "value": 381.2,
                    "range": [120.0, 550.0],
                    "type": "c",
                    "vary": True,
                    "latex": r"$T$",
                },
            }
        ),
        ("Extinction", {
                "function": "Smith07",
                "logtau": {
                    "value": -6.495,
                    "range": [-8.0, 3],  # [-4.0, 1.5],
                    "type": "c",
                    "vary": True,#
                    "latex": r"$\mathrm{log}\,\tau_\mathrm{ext}$",
                },
                "multiply": ["PAH", "MBB1", "MBB2", "MBB3"]
            }
         ),
    )
)

parTruth = None  #Whether to provide the truth of the model
#modelUnct = False #False #Whether to consider the model uncertainty in the fitting
unctDict = OrderedDict(
    (
        ("lnf"  , [-10., 2]),
        ("lna"  , [-10., 3]),
        ("lntau", [-10., 0.]),
    )
)
################################################################################
#                                   emcee                                      #
################################################################################

burnin1 = OrderedDict(
    (
        ("sampler"  , "EnsembleSampler"),
        ("nwalkers" , 128), #The number of walkers.
        ("iteration", [4000, 4000]),            #[1000, 1000, 1000]), #The iteration of burn-in run.
        ("thin"     , 1), #To thin the recorded sample.
        ("ball-r"   , 0.1), #The radius of the ball as the fraction of the full range.
    )
)
burnin2 = OrderedDict(
    (
        ("sampler"  , "PTSampler"),
        ("ntemps"   , 16), #The number of temperature ladders only for PTSampler.
        ("nwalkers" , 128), #The number of walkers.
        ("iteration", [1000, 1000]), #The iteration of burn-in run.
        ("thin"     , 1), #To thin the recorded sample.
        ("ball-r"   , 0.1), #The radius of the ball as the fraction of the full range.
    )
)
final = OrderedDict(
    (
        ("sampler"  , "EnsembleSampler"),
        ("ntemps"   , 4),
        ("nwalkers" , 128),
        ("iteration", [1000, 600]),              #[1000, 600]),
        ("thin"     , 1),
        ("ball-r"   , 0.01),
    )
)
setup = OrderedDict(
    (
        ("threads"  , 4), #Not used if MPI is using.
        ("printfrac", 0.1),
        ("pslow"    , 16),
        ("pscenter" , 50),
        ("pshigh"   , 84),
    )
)
emceeDict = OrderedDict(
    (
        ("BurnIn-1", burnin1),
        #("BurnIn-2", burnin1),
        ("Final", final),
        ("Setup", setup),
    )
)

#Postprocess#
#-----------#
ppDict = {
    "burn-in" : 300,
    "low"     : 16,
    "center"  : 50,
    "high"    : 84,
    "nuisance": True, #False, #
    "fraction": 0.1, #The fraction of walkers to be dropped.
    "savepath": "test/NGC5194_SED/mock/Fitted_Result/", # The path to save the results
}
