'''
Beam tools should be a collection of helper function to calculate beam properties
in both real and multipole space.
'''

import sys, os
import matplotlib
matplotlib.use('Agg')
import numpy as np
import gaussian_fit as gf
import pdata as pd
import scipy.interpolate
from scipy import optimize
import time

# numel = [nx, ny]
# cr [min(x), min(y), max(x), max(y)]

def wlng(blm):
    '''
    Returns the non-Gaussian beam window function from a blm arry
    '''

    wlms = np.abs(blm**2)
    wlng = np.sum(wlm[:,1:], axis=1)
    wlg = np.sum(wlm[:,0], axis=1)
    maxwl = np.max(wlg)
    wlr = wlng/wlg

    return wlr

def edge_taper(taper_angle, etaper, angle):

    s = np.sqrt(-taper_angle**2/2/np.log(10**(etaper/10.)))
    print('Beam width {:5.3f} deg'.format(s))

    return 10*np.log10(np.exp(-angle**2/2/s**2))

def ghz2mm(nu):

    return 1000*2.998e8/(nu*1e9)

def pixel2fwhm(d, wavel=3.):

    '''
    Converts pixel size to beam FWHM assuming the Rayleigh criterion. Returns
    an error if the pixel size is smaller than the wavelength.

    d       : Pixel size [mm]
    wavel   : Wavelength [mm]
    Returns : FWHM [deg]
    '''

    if 0.973*wavel/d > 1.0:
        raise ValueError('Pixel size smaller than wavelength')

    fwhm_in_lambda_over_d = 0.97354

    #return np.degrees(np.arcsin(1.22*wavel/d))
    return np.degrees(fwhm_in_lambda_over_d*wavel/d)

def aperture2fwhm(d, wavel=3.):

    return pixel2fwhm(d, wavel=wavel)

def fwhm2sigma(fwhm):

    return fwhm/np.sqrt(8*np.log(2))

def sigma2fwhm(sigma):

    return np.sqrt(8*np.log(2))*sigma

def fwhm2s(fwhm):

    return fwhm / np.sqrt(8*np.log(2))

def s2fwhm(s):

    return s * np.sqrt(8*np.log(2))

def pixel2sigma(d, wavel=3.):

    return fwhm2sigma(pixel2fwhm(d, wavel))

def pixel2taper(d, taper_angle, wavel=3.):
    '''
    Give it a pixel size and an angle and it will calculate the dB taper
    at that angle.
    '''

    return -10*np.log10(np.exp(-taper_angle**2/2. \
        /pixel2sigma(d, wavel)**2))

def taper_diff(d, etaper, taper_angle, wavel):

    return np.abs(etaper-pixel2taper(d, taper_angle, wavel))

def taper2pixel(etaper, taper_angle, wavel=3.):

    if etaper < 0:
        raise ValueError('The edge taper (by definition here) should be > 0')

    d0 = 6
    dopt = optimize.minimize(taper_diff, d0, args=(etaper, taper_angle, wavel),
        bounds=((1.3*wavel, 10*wavel),))

    tdiff = taper_diff(float(dopt['x']), etaper, taper_angle, wavel)

    if tdiff < 0.01:
        return float(dopt['x'])

    else:
        print(dopt['message'])
        raise ValueError('Invalid configuration')

def profile2map(theta, profile, cr, numel):

    X, Y, _ = get_mesh(cr, numel)
    R = np.sqrt(X**2 + Y**2)

    profilei = np.interp(R, theta, profile)
    arr1 = np.reshape(profilei, (numel[0], numel[1]))

    # print(profilei)
    # print(arr1)
    # print(np.min(profilei), np.max(arr1), np.mean(profilei))
    # print(np.min(arr1), np.max(arr1), np.mean(arr1))

    return arr1

def profile_bsa(theta, profile, cr, numel):

    arr1 = profile2map(theta, profile, cr, numel)
    dA = get_da(cr, numel)
    bsa_out = bsa(arr1, dA, normalize=True)

    return bsa(arr1, dA, normalize=True)

def bsa(arr, dA, maxval=None, normalize=True):

    if maxval is None and normalize:
        maxval = np.max(arr.flatten())

        return dA * np.sum(arr/maxval)

    return dA * np.sum(arr)

def bsa_two_field(arr1, arr2, dA, maxval=None, normalize=True):

    arr = np.sqrt(arr1**2 + arr2**2)

    if maxval is None and normalize:
        maxval = np.max(arr.flatten())

        return dA*np.sum(arr/maxval)

    return dA*np.sum(arr)

def bsa_cross(arr1, arr2, dA):

    maxval = np.max(arr1.flatten())

    return dA*np.sum(arr1 / maxval), dA*np.sum(arr2 / maxval),




def bsa_r(x, y, arr, dA, rmax, xtrans=0., ytrans=0.):
    '''
    Returns the beam power inside a radius, total beam power, and the
    fractional difference
    '''

    bsa_total = bsa(arr, dA)
    r = np.sqrt((x-xtrans)**2+(y-ytrans)**2)
    ridx = (r < rmax)
    arr_inside = arr[ridx]
    bsa_inside = bsa(arr_inside, dA)

    return bsa_inside, bsa_total, bsa_total-bsa_inside

def power_r(x, y, arr, rmax, xtrans=0., ytrans=0.):
    '''
    Returns the beam power inside a radius, total beam power, and the
    fractional difference
    '''

    power_total = bsa(arr, 1)
    r = np.sqrt((x-xtrans)**2+(y-ytrans)**2)
    ridx = (r < rmax)
    arr_inside = arr[ridx]
    power_inside = bsa(arr_inside, 1)

    return power_inside, power_total, power_total-power_inside

def bsa_ellipse(x, y, arr, dA, ellipse):
    '''
    Based on Stack overflow answer
    '''

    # Ellipse properties
    ex, ey = ellipse[0], ellipse[1]
    # Half semi major and semi minor axis
    esx, esy = ellipse[2], ellipse[3]
    esth = ellipse[4]

    c = np.cos(np.radians(180. - esth))
    s = np.sin(np.radians(180. - esth))

    # Translated
    xc, yc = x - ex, y - ey

    # Rotated
    xct = xc*c - yc*s
    yct = xc*s + yc*c

    rad_cc = xct**2/esx**2 + yct**2/esy**2

    npix = float(len(x.flatten()))
    inidx = (rad_cc < 1.)
    ninsise = np.sum(inidx)

    # print('npix = {:f}'.format(int(npix)))
    # print('ninside {:d}'.format(ninsise))
    # print('Fraction of pixels inside is {:3.2f} %'.format(100*ninsise/npix))
    # print(np.shape(arr))
    # print(np.shape(inidx))

    bsa_total = bsa(arr, dA)
    arr_inside = arr[inidx]

    bsa_inside = bsa(arr_inside, dA)

    return bsa_inside, bsa_total

def elliptical_profile(x, y, arr, dA, ellipse, nbins=200):

    # Ellipse properties
    ex, ey = ellipse[0], ellipse[1]
    # Half semi major and semi minor axis
    esx, esy = ellipse[2], ellipse[3]
    esth = ellipse[4]

    c = np.cos(np.radians(180. - esth))
    s = np.sin(np.radians(180. - esth))

    # Translated
    xc, yc = x - ex, y - ey

    # Rotated
    xct = xc*c - yc*s
    yct = xc*s + yc*c

    rad_cc = np.sqrt(xct**2/esx**2 + yct**2/esy**2)

    rmax = rad_cc.max()

    rs = np.linspace(0, rmax, nbins+1)
    drs = np.median(np.diff(rs))
    profile = np.nan*np.ones(nbins)

    for i in xrange(nbins-1):
        i1 = (rad_cc > rs[i])
        i2 = (rad_cc < rs[i+1])

        profile[i] = np.mean(arr[np.logical_and(i1, i2)])


    return rs[:-1]+drs/2.0, profile

def get_profile(cr, numel, arr, fg=False):
    '''
    '''

    ssx, ssy = step_size(cr, numel)
    _, _, sx, sy = gf.centerofmass(arr)
    idxs = center_idxs(arr)
    cx, cy = idxs[1], idxs[0]

    if fg:
        dA = get_da(cr,numel)
        bsolid = bsa(arr, dA)
        fgain = 4*np.pi/degsq2srad(bsolid)

        arr = arr/arr.max()
        arr = fgain*arr

    r, profile = radial_profile(arr, [cx, cy], step_size=ssx)

    return r, profile

def get_profile_wcustom_weight():
    pass

def radial_profile(arr, center, step_size=None):
    '''
    Creates a radial profile. Function found on stack overflow.
    '''

    y, x = np.indices((arr.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), arr.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr

    ru = np.unique(r.ravel()).astype(np.float)
    if step_size is not None:
        ru *= step_size

    return ru, radialprofile

def x_profile(x, arr, step_size=None):
    '''
    Creates a radial profile. Function found on stack overflow.
    '''

    # print(type(x))
    # print(type(arr))
    # print(np.sum(np.isnan(arr)))

    bins = np.linspace(np.min(x), np.max(x), 50)

    # print(bins)
    digitized = np.digitize(x, bins)
    bin_means = [np.nanmean(arr[digitized == i]) for i in range(1, len(bins))]

    # print(np.shape(bins))
    # print(np.shape(bin_means))

    bins = bins[:-1] + np.median(np.diff(bins))

    return bins, np.array(bin_means)

def slice_profile(x, y, arr):#, center):

    edge1 =  np.abs(0.1 * x)
    edge2 = -np.abs(0.1 * x)

    idx = (y < edge1) & (y > edge2)
    r, profile = x_profile(x[idx], arr[idx])

    print('Testing slice profile:')
    print(np.sum(idx)/ float(arr.size))


    return r, profile

def gfit(cr, numel, arr1, verbose=False, gfwhm=0.03, gell=0.01, fit_radius=0.7,
        return_model=False, return_angle=False, benchmark=False, center_on_peak=True):


    t0 = time.time()
    if center_on_peak:
        _, cr_out, numel_out, arr1_out = center_map(cr, numel, arr1)
    else:
        cr_out = cr
        numel_out = numel
        arr1_out = arr1

    # return [yl, yh, xl, xh], cr2use, numel2use, arr2use

    t1 = time.time()
    if numel_out[0] < 1 or numel_out[1] < 1:
        cr_out = cr
        numel_out = numel
        arr1_out = arr1

    xx, yy, dA = get_mesh(cr_out, numel_out)
    x, y = get_lin(cr_out, numel_out, mult=1)

    t2 = time.time()

    cx, cy, _, _ = gf.centerofmass(np.flipud(arr1_out))
    cx0 = np.interp(cx, np.arange(len(x)), x)
    cy0 = np.interp(cy, np.arange(len(y)), y)

    t3 = time.time()

    Data = {'data':arr1_out, 'x':xx, 'y':yy}

    # cx, cy, angle, fwhm, ellipticity, peak
    srln = np.sqrt(8*np.log(2))
    pguess1 = [cx0, cy0, 0.0, gfwhm, gell, np.max(arr1_out)]
    po, cov, model_out1, res1 = gf.egfit1(Data, fit_radius=fit_radius,
        fit4dc=False, pguess=pguess1)

    t4 = time.time()

    fwhm_x = po[3]/np.sqrt(po[4])
    sx = fwhm_x /np.sqrt(8*np.log(2))
    sy = sx*po[4]

    sx = np.abs(sx)
    sy = np.abs(sy)

    cx = po[0]
    cy = po[1]

    ellipticity = (np.max([sx, sy])-np.min([sx, sy])) / \
        (sx+sy)

    t5 = time.time()

    if verbose:

        print('Elliptical Gaussian beam fit results:')
        print('  Amplitude      : {:5.2f}'.format(po[5]))
        print('  X centroid     : {:5.3f} deg'.format(po[0]))
        print('  Y centroid     : {:5.3f} deg'.format(po[1]))
        print('  FWHM              : {:5.3f} arcmin'.format(np.abs(po[3])*60))
        print('  Rotation angle : {:5.2f} deg'.format(po[2]))
        print('  Ellipticity    : {:5.4f}'.format(ellipticity))

    if benchmark:
        print('Timing results:')
        for t in [t0, t1, t2, t3, t4, t5]:
            print('  t = {}'.format(t-t0))

    if return_model:
        return cx, cy, sx, sy, np.abs(po[3]), ellipticity, cr_out, numel_out, model_out1

    if return_angle:
        return cx, cy, sx, sy, np.abs(po[3]), ellipticity, po[2]

    return cx, cy, sx, sy, np.abs(po[3]), ellipticity

def gfit_plane(cr, numel, arr1, verbose=True, gfwhm=0.03, gell=0.01, fit_radius=0.7,
        return_model=False, benchmark=False, hard_radius=150):


    t0 = time.time()
    # _, cr_out, numel_out, arr1_out = center_map(cr, numel, arr1)

    cr_out, numel_out, arr1_out = np.copy(cr), np.copy(numel), np.copy(arr1)


    t1 = time.time()


    xx, yy, dA = get_mesh(cr_out, numel_out)
    x, y = get_lin(cr_out, numel_out, mult=1)

    r = np.sqrt(xx**2 + yy**2)
    ifit = (r < hard_radius)
    arr1_out[~ifit] = 0.0

    t2 = time.time()
    # print('cr_out', cr_out)
    # print('numel_out', numel_out)

    cx, cy, _, _ = gf.centerofmass(np.flipud(arr1_out))
    cx0 = np.interp(cx, np.arange(len(x)), x)
    cy0 = np.interp(cy, np.arange(len(y)), y)

    t3 = time.time()

    Data = {'data':arr1_out, 'x':xx, 'y':yy}

    # cx, cy, angle, fwhm, ellipticity, peak
    srln = np.sqrt(8*np.log(2))
    pguess1 = [cx0, cy0, 0.0, gfwhm, gell, np.max(arr1_out)]


    print('pguess1', pguess1)
    print('fit_radius', fit_radius)

    po, cov, model_out1, res1 = gf.egfit1(Data, fit_radius=fit_radius,
        fit4dc=False, pguess=pguess1, force_center=True)

    # pguess2 = [np.max(arr1_out), cx0, cy0, gfwhm/srln, gfwhm/srln, 0.0]
    # po, cov, model_out1, res1 = gf.egfit2(xx, yy, arr1_out, fit_radius=fit_radius,
        # pguess=pguess2)


    t4 = time.time()

    fwhm_x = po[3]/np.sqrt(po[4])
    sx = fwhm_x /np.sqrt(8*np.log(2))
    sy = sx*po[4]

    sx = np.abs(sx)
    sy = np.abs(sy)

    cx = po[0]
    cy = po[1]

    ellipticity = (np.max([sx, sy])-np.min([sx, sy])) / \
        (sx+sy)

    t5 = time.time()

    if verbose:

        print('Elliptical Gaussian beam fit results:')
        print('  Amplitude      : {:5.2f}'.format(po[5]))
        print('  X centroid     : {:5.3f} mm'.format(po[0]))
        print('  Y centroid     : {:5.3f} mm'.format(po[1]))
        print('  FWHM              : {:5.3f} mm'.format(np.abs(po[3])))
        print('  Rotation angle : {:5.2f} deg'.format(po[2]))
        print('  Ellipticity    : {:5.4f}'.format(ellipticity))

    if benchmark:
        print('Timing results:')
        for t in [t0, t1, t2, t3, t4, t5]:
            print('  t = {}'.format(t-t0))

    if return_model:
        return cx, cy, sx, sy, np.abs(po[3]), ellipticity, cr_out, numel_out, model_out1

    return cx, cy, sx, sy, np.abs(po[3]), ellipticity

def gfit1d(theta, profile, verbose=True, gfwhm=None, fit_radius=None,
        return_model=False):

    if fit_radius:
        idx = np.argmin(np.abs(theta - fit_radius))
    else:
        idx = None

    sg = np.sum(theta * profile)/np.sum(profile)
    fwhmg = s2fwhm(sg)
    profile2u = profile / np.max(profile)

    def gprofile(r, fwhm):

        s = fwhm/np.sqrt(8*np.log(2))

        return np.exp(-(r**2)/(2.0 * s**2))

    popt, pcov = optimize.curve_fit(gprofile, theta, profile, p0=fwhmg)

    # print(theta)
    # print(profile)
    # print(gprofile(theta, popt[0]))
    # print(popt)

    return popt


    # # cx, cy, angle, fwhm, ellipticity, peak
    # pguess1 = [cx0, cy0, 0.0, gfwhm, gell, np.max(arr1_out)]
    # po, cov, model_out1, res1 = gf.egfit1(Data, fit_radius=fit_radius,
    #     fit4dc=False, pguess=pguess1)

    # fwhm_x = po[3]/np.sqrt(po[4])
    # sx = fwhm_x /np.sqrt(8*np.log(2))
    # sy = sx*po[4]

    # sx = np.abs(sx)
    # sy = np.abs(sy)

    # cx = po[0]
    # cy = po[1]

    # ellipticity = (np.max([sx, sy])-np.min([sx, sy])) / \
    #     (sx+sy)

    # if verbose:

    #     print('Elliptical Gaussian beam fit results:')
    #     print('\tAmplitude      : {:5.2f}'.format(po[5]))
    #     print('\tX centroid     : {:5.3f} deg'.format(po[0]))
    #     print('\tY centroid     : {:5.3f} deg'.format(po[1]))
    #     print('\tFWHM              : {:5.3f} arcmin'.format(np.abs(po[3])*60))
    #     print('\tRotation angle : {:5.2f} deg'.format(po[2]))
    #     print('\tEllipticity    : {:5.4f}'.format(ellipticity))

    # if return_model:
    #     return cx, cy, sx, sy, np.abs(po[3]), ellipticity, cr_out, numel_out, model_out1

    # return cx, cy, sx, sy, np.abs(po[3]), ellipticity

def power_atr(ax, x, y, arr, cr, numel, r):
    '''
    Calculates the power at radius r
    '''

    dppix = (cr[2]-cr[0])/numel[0]
    maxarr = np.max(arr.flatten())
    xc, yc = int(0.5*numel[0]), int(0.5*numel[1]+r/dppix)

    print('Relative power at radius r: {:5.5f} / {:5.2f} dB'.\
        format(arr[xc,yc]/maxarr,
        10*np.log10(arr[xc,yc]/maxarr)))


    ax.plot(xc, yc, zorder=5, marker='o', color='black', ms=2)

def power_ghost(cr, numel, arr, xg=-2.0, yg=0.0, rmax=0.5, verbose=True):

    X, Y, _ = get_mesh(cr, numel)
    x = np.linspace(cr[0], cr[2], numel[0])
    y = np.linspace(cr[1], cr[3], numel[1])

    xidx = np.argmin(np.abs(x - xg))
    yidx = np.argmin(np.abs(y - yg))
    xi = x[xidx]
    yi = y[yidx]

    Rg = np.sqrt((X - xi)**2 + (Y -yi)**2)

    gidx = (Rg < rmax)
    mean_val = np.mean(arr[gidx])
    max_val = np.max(arr[gidx])

    if verbose:
        print('Mean ghost val: {:.3f}'.format(mean_val))
        print('Ratio to peak (mean): {:.3f} dB'.format(10*np.log10(mean_val/np.max(arr))))
        print('Ratio to peak (max) : {:.3f} dB'.format(10*np.log10(max_val/np.max(arr))))

    masked_arr = np.copy(arr)
    masked_arr[gidx] = np.nan

    return masked_arr

def cr_str(cr):

    return '[{:5.2f}, {:5.2f}, {:5.2f}, {:5.2f}]'.format(cr[0], cr[1], cr[2], cr[3])

def get_lin(cr, numel, mult=10):

    x2 = np.linspace(cr[0], cr[2], mult*numel[0])
    y2 = np.linspace(cr[1], cr[3], mult*numel[1])

    return x2, y2

def get_da(cr, numel, mult=1):

    # return ((cr[2]-cr[0])/(mult*numel[0]))**2
    return ((float(cr[2])-float(cr[0]))/(mult*float(numel[0])))**2

def get_mesh(cr, numel, mult=1.0):

    xx, yy = np.meshgrid(np.linspace(cr[0], cr[2], int(mult*numel[0])),
        np.linspace(cr[1], cr[3], int(mult*numel[1])))

    dA = get_da(cr, numel, mult=1)

    return xx, yy, dA

def get_r(cr, numel, mult=1):

    xx, yy = np.meshgrid(np.linspace(cr[0], cr[2], mult*numel[0]),
        np.linspace(cr[1], cr[3], mult*numel[1]))

    r = np.sqrt(xx**2 + yy**2)

    return r

def get_crnumel(xx, yy):

    xxs = np.shape(xx)
    numel = [xxs[1], xxs[0]]
    cr = [xx.min(), yy.min(), xx.max(), yy.max()]

    return cr, numel

def arcminsq2srad(arcmin2):
    '''
    https://en.wikipedia.org/wiki/Square_degree
    '''

    return arcmin2/(41252.9612494/4/np.pi)/3600

def nu2wavel(nu):

    c = 2.99792458e8
    wavel = c/nu
    return wavel


def degsq2srad(deg2):
    '''
    https://en.wikipedia.org/wiki/Square_degree
    '''

    return 4*np.pi*deg2/(129600/np.pi)

    # (129600/np.pi) = 41252.9612494

def srad2sqdeg(srad):

    #return srad/(4*np.pi)*(129600/np.pi)
    return srad/(4*np.pi)*(129600/np.pi)

def srad2sqarcmin(srad):

    return srad/(4*np.pi)*(129600/np.pi)*3600

def wl_gauss(ell, fwhm=7.0):
    '''
    Calulates the Gaussian beam window function where fwhm is the full
    width at half maximum in degrees.
    return = Wl
    '''
    fwhm = fwhm/180.0*np.pi
    sigma = fwhm/np.sqrt(8*np.log(2))

    wl_out = np.exp(-ell*(ell+1)*sigma**2)
    wl_out = wl_out/np.max(wl_out)
    wl_out[0] = 1.0

    return wl_out

def bl_gauss(ell, fwhm=7.0):
    '''
    Calulates the Gaussian beam window function where fwhm is the full
    width at half maximum in degrees.
    return = Wl
    '''
    fwhm = fwhm/180.0*np.pi
    sigma = fwhm/np.sqrt(8*np.log(2))

    bl_out = np.exp(-0.5*ell*(ell+1)*sigma**2)
    bl_out = bl_out/np.max(bl_out)
    bl_out[0] = 1.0

    return bl_out

def center_idxs(arr):
    '''
    Returns the indexes corresponding to the location of the peak pixel value

    The first dimension returned is y-axis (vertical)
    The second dimension returned is x-axis (horizontal)
    '''

    return np.unravel_index(np.argmax(arr), np.shape(arr))

def dist2edge(cr, numel, arr):

    idxs = center_idxs(arr)
    #dists = np.array([idxs[0], numel[1]-idxs[0], idxs[1], numel[0]-idxs[1]])
    dists = np.array([idxs[0], numel[0]-idxs[0], idxs[1], numel[1]-idxs[1]])

    return idxs, dists

def center_map(cr, numel, arr, ix=None):
    '''
    Centers a map on the peak value
    '''

    X, Y, _ = get_mesh(cr, numel)

    if ix is None:

        idxs, dists = dist2edge(cr, numel, arr)
        sh = np.min(dists)
        xl, xh = (idxs[1]-sh), (idxs[1]+sh-1)
        yl, yh = (idxs[0]-sh), (idxs[0]+sh-1)

    else:

        xl, xh = ix[2], ix[3]
        yl, yh = ix[0], ix[1]

    # print('center map')
    # print(xl, xh)
    # print(yl, yh)

    cr2use = [X[0,xl], Y[yl,0], X[0,xh], Y[yh,0]]
    numel2use = [yh-yl, xh-xl, 0]
    arr2use = arr[yl:yh, xl:xh]

    return [yl, yh, xl, xh], cr2use, numel2use, arr2use

def trunc_map(X, Y, arr, ratio=0.2):
    '''
    Shrinks the size of the map by ratio towards center
    '''

    nx, ny = np.shape(arr)

    midx = (nx+1)/2.0
    midy = (ny+1)/2.0

    lx, rx = int(midx-np.round(ratio/2*nx)), int(midx+np.round(ratio/2*nx))
    ly, ry = int(midy-np.round(ratio/2*ny)), int(midy+np.round(ratio/2*ny))

    numel = np.shape(X[lx:rx+1,ly:ry+1])
    cr = [X[0,lx], Y[ly,0], X[0,rx+1], Y[ry+1,0]]

    return [lx, rx+1, ly, ry+1], cr, numel, \
        X[lx:rx+1,ly:ry+1], Y[lx:rx+1,ly:ry+1], arr[lx:rx+1,ly:ry+1]

def center_cr(cr):

    crx = np.mean([cr[0], cr[2]])
    cry = np.mean([cr[1], cr[3]])
    return [cr[0]-crx, cr[1]-cry, cr[2]-crx, cr[3]-cry]

def interp_map(cr, numel, arr, box=None, cri=None, numeli=None,
    mult=5, center=False, centered=False):
    '''
    Outputs an interpolated map according to prescription
    '''
    # Centering map beforehand if requested
    if center:
        cr, numel, arr = center_map(cr, numel, arr)

    # Identifying the location of the maximum value
    if not centered:
        X, Y, _ = get_mesh(cr, numel)
        idxs, dists = dist2edge(cr, numel, arr)
        xc, yc = X[idxs], np.flipud(Y)[idxs]
    else:
        xc, yc = 0.0, 0.0

    if box is not None:
        cri = [xc-box/2., yc-box/2., xc+box/2., yc+box/2.]

    if cri is None and numeli is None:
        xx2, yy2 = get_lin(cr, numel, mult=mult)
        cri, numeli = get_crnumel(xx2, yy2)
    elif box is not None and numeli is not None:
        xx2, yy2 = get_lin(cri, numeli, mult=1)
    elif cri is not None and numeli is not None:
        xx2, yy2 = get_lin(cri, numeli, mult=1)
    else:
        raise ValueError('Incorrect parameter values')

    xx, yy, dA = get_mesh(cr, numel)

    fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0],
        np.rot90(arr, k=-1))
    arr2use = np.rot90(fint(xx2, yy2), k=+1)

    return cri, numeli, arr2use

def simple_interp(cr, numel, arr, box, numeli):

    xx, yy, dA = get_mesh(cr, numel)

    fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0],
        np.rot90(arr, k=-1))

    crx = np.mean([cr[0], cr[2]])
    cry = np.mean([cr[1], cr[3]])

    cri = [crx-box/2., cry-box/2., crx+box/2., cry+box/2.]

    xx2, yy2 = get_lin(cri, numeli, mult=1)
    arr2use = np.rot90(fint(xx2, yy2), k=+1)

    return cri, numeli, arr2use

def null_outside(cr, numel, arr, radius=1.2, center=None):

    X, Y, _ = get_mesh(cr, numel)
    if center is None:
        center = peak_idx(arr)

    R = np.sqrt((X-X[center])**2 + (Y-Y[center])**2)

    idx = (R > radius)

    arr[idx] = 0.0

    return arr

def peak_idx(arr):

    return np.unravel_index(np.nanargmax(arr.flatten()), np.shape(arr))

def peak_loc(arr, cr, numel):

    pidx = peak_idx(arr)
    print('pidx:')
    print(pidx)

    X, Y, da = get_mesh(cr, numel)

    xpeak = X[pidx]
    ypeak = Y[pidx]

    return xpeak, ypeak

def step_size(cr, numel):

    return (cr[2]-cr[0])/float(numel[0]), (cr[3]-cr[1])/float(numel[1])

def gauss_profile(r, fwhm, A=1):

    s = fwhm/np.sqrt(8*np.log(2))

    return A*np.exp(-(r**2)/(2*s**2))

def cal_profile(fpu, r, fwhm):

    beam07={'x1':[0.9955,28.4537,0.0778,16.5750,0.0490,3.7256],
        'x3':[0.9975,28.2710,0.0007,28.5989,0.0533,4.3278],
        'x5':[0.9987,28.7889,0.1539,36.3706,0.0609,5.0126],
        'x2':[1.0016,41.5220,0.9017,28.6508,0.0714,3.8988],
        'x4':[0.9995,41.5459,0.7724,29.8233,0.0697,3.8579],
        'x6':[1.0031,41.2828,0.8139,25.6374,0.0851,4.9568]}

    beam12={'x1':[0.9977, 28.5591, 0.0626, 16.5663, 0.0494, 3.793],
        'x3':[1.0009, 28.4219, -0.0561, 28.6087, 0.055, 4.8192],
        'x5':[0.9999, 28.8655, 0.2435, 36.3747, 0.0607, 4.9926],
        'x2':[1.0089, 41.7397, 0.0576, 12.7253, 0.0649, 3.9704],
        'x4':[1.0042, 41.7285, 0.8071, 29.8306, 0.0646, 3.8466],
        'x6':[1.0055, 41.3063, 0.1024, 11.3075, 0.0629, 4.2918]}

    # A1,w1,A2,k,A3,k3=beam07['x'+str(fpu)]
    A1,w1,A2,k,A3,k3=beam12['x'+str(fpu)]

    w2=w1*k
    w3=w1*k3
    A=A1+A2+A3
    A1/=A
    A2/=A
    A3/=A

    print('Debug:')
    print(fwhm)
    print(w2/60.)
    print(w3/60.)

    print('A1: {:5.5f}'.format(A1))
    print('A2: {:5.5f}'.format(A2))
    print('A3: {:5.5f}'.format(A3))

    print('fwhm {:5.5f} deg'.format(fwhm))
    print('w2 {:5.5f} deg'.format(w2/60.))
    print('w3 {:5.5f} deg'.format(w3/60.))

    profile1 = gauss_profile(r, fwhm, A=A1)
    profile2 = gauss_profile(r, w2/60., A=A2)
    profile3 = gauss_profile(r, w3/60., A=A3)

    profile = profile1+profile2+profile3

    profile /= profile.max()

    return profile

def cal_beam(fpu, fwhm=None, lmax=700, normalize=True):
    """
    WARNING: This function is used for calibration pipeline.
    Make sure you understand what it does before you use it.

    Return the triple gaussian B_ells for calibration analysis.
    The gaussian beam to be nomalized to 1: B_ell(l=0)=1.
    Key feature: adjustable main fwhm.

    Arguments
    ---------
    fpu : int
        Focal Plane that being analysed. Must be a int number within [1,2,3,4,5,6]
    fwhm : float, optional
        Main gaussian-beam fwhm(w1). Default to be the width in
        beam product. Adjustable. Unit:arcmin.
    beam_product: string, optional
        Name of beam_product values. Default to be 'beam07'.
        Options:'beam02''beam03''beam04''beam05''beam06''beam07'
    lmax: int, optional
        Maximum value of multipoles. Default to be 700.

    Returns
    -------
    B_ells: list
        Triple gaussian beam with the second and third gaussian
        fixed and the main beam to be input fwhm.

        Equations:B_ells=[A1*gauss(fwhm)+A2*gauss(w2)+A3*gauss(w3)]/(A1+A2+A3)

    """
    beam07={'x1':[0.9955,28.4537,0.0778,16.5750,0.0490,3.7256],
        'x3':[0.9975,28.2710,0.0007,28.5989,0.0533,4.3278],
        'x5':[0.9987,28.7889,0.1539,36.3706,0.0609,5.0126],
        'x2':[1.0016,41.5220,0.9017,28.6508,0.0714,3.8988],
        'x4':[0.9995,41.5459,0.7724,29.8233,0.0697,3.8579],
        'x6':[1.0031,41.2828,0.8139,25.6374,0.0851,4.9568]}

    beam08={'x1':[0.9963,28.4535,0.0779,16.5750,0.0491,3.7252],
         'x3':[0.9974,28.2710,0.0071,28.5989,0.0533,4.3278],
         'x5':[0.9990,28.7889,0.1537,36.3706,0.0609,5.0130],
         'x2':[1.0016,41.5220,0.9018,28.6508,0.0714,3.8986],
         'x4':[0.9989,41.5396,0.7970,29.8272,0.0696,3.8416],
         'x6':[1.0031,41.2828,0.8139,25.6374,0.0851,4.9568]}

    fplist=[1,2,3,4,5,6]

    Ef = fpu
    fp='x'+str(Ef)

    A1,w1,A2,k,A3,k3=beam08[fp]

    w2=w1*k
    w3=w1*k3
    A=A1+A2+A3

    if normalize:
        A1/=A
        A2/=A
        A3/=A

    #Bl2=A2*hp.sphtfunc.gauss_beam(np.deg2rad(w2/60),lmax)
    #Bl3=A3*hp.sphtfunc.gauss_beam(np.deg2rad(w3/60),lmax)
    Bl2=A2*gauss_beam(np.deg2rad(w2/60),lmax)
    Bl3=A3*gauss_beam(np.deg2rad(w3/60),lmax)

    if fwhm is not None:
        w1=fwhm

    #Bl1=A1*hp.sphtfunc.gauss_beam(np.deg2rad(w1/60),lmax)
    Bl1=A1*gauss_beam(np.deg2rad(w1/60),lmax)

    #pd.keyboard()

    B_ells=Bl1+Bl2+Bl3

    return B_ells

def physical(p):

    A1,w1,A2,k2,A3,k3=p
    w2=w1*k2
    w3=w1*k3

    w1rad=np.deg2rad(w1/60)
    w2rad=np.deg2rad(w2/60)
    w3rad=np.deg2rad(w3/60)

    sigma1=w1rad/np.sqrt(8*np.log(2))
    sigma2=w2rad/np.sqrt(8*np.log(2))
    sigma3=w3rad/np.sqrt(8*np.log(2))

    theta=np.arange(0,1000*sigma1,0.00001)
    thetaarcmin=np.rad2deg(theta)*60

    TG1=A1*1/sigma1*np.exp(-theta**2/2/(sigma1**2))
    TG2=A2*1/sigma2*np.exp(-theta**2/2/(sigma2**2))
    TG3=A3*1/sigma3*np.exp(-theta**2/2/(sigma3**2))

    TG=(TG1+TG2+TG3)/(A1/sigma1+A2/sigma2+A3/sigma3)
    return thetaarcmin,TG

def triple(p):

    A1,w1,A2,k,A3,k3=p
    w2=w1*k
    w3=w1*k3
    #DGBl1=A1*hp.sphtfunc.gauss_beam(np.deg2rad(w1/60),700)
    #DGBl2=A2*hp.sphtfunc.gauss_beam(np.deg2rad(w2/60),700)
    #DGBl3=A3*hp.sphtfunc.gauss_beam(np.deg2rad(w3/60),700)
    DGBl1=A1*gauss_beam(np.deg2rad(w1/60),700)
    DGBl2=A2*gauss_beam(np.deg2rad(w2/60),700)
    DGBl3=A3*gauss_beam(np.deg2rad(w3/60),700)
    DGBls=DGBl1+DGBl2+DGBl3

    return DGBls

def gauss_beam(w, lmax):

    s = w/np.sqrt(8*np.log(2))
    ell = np.arange(lmax+1)

    return np.exp(-0.5*(ell+1)*ell*s**2)

def beam_params(cr, numel):

    midx, midy = (numel[0]+1)/2, (numel[1]+1)/2
    dx, dy = (cr[2]-cr[0])/float(numel[0]-1), (cr[3]-cr[1])/float(numel[1]-1)

    # Hack to get rid of floating point differences
    if dx != dy and np.abs((dx-dy)/dx) < 1e-15:
        dx = dy

    return midx, midy, np.radians(dx), np.radians(dy)

def litebird_pixels():

    nus = np.array([100, 140, 195, 119, 166, 195, 280, 235, 337, 402])
    pixel_sizes = np.array([11.6, 11.6, 11.6, 11.6, 11.6, 6.6, 6.6, 6.6, 6.6, 5.7])
    taper_angles = np.array([12.3356, 8.878, 6.398, 10.41, 7.504, 12.9449, 9.094, 10.798, 7.576, 7.356])

    wavels = 1000. * 2.998e8/(1e9 * nus)

    taper_angle = 16.0
    for nu, psize, wavel, taper_angle in zip(nus, pixel_sizes, wavels, taper_angles):

        taper_angle = 12.7
        print('{} GHz | Edge taper is -{:.2f} dB ({:.1f} deg FWHM) for {:.1f} mm pixel size'.\
            format(nu, pixel2taper(psize, taper_angle, wavel=wavel), pixel2fwhm(psize, wavel=wavel), psize))

def main():

    litebird_pixels()
    return

    taper_angle = 25.
    pixel_size = [4., 5., 8., 10.]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle), ps))


    # etaper = [1.72, 3.40, 10.15, 16.33]
    # for eta in etaper:
    #     print('Pixel is {:.2f} mm'.format(taper2pixel(eta, taper_angle, wavel=3.0)))

    print('40 GHz')
    taper_angle = 16
    # taper_angle = 8
    pixel_size = [8, 10, 12, 14]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=6.75), ps))

    print('90 GHz')
    pixel_size = [4.3, 5.3, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=3.0), ps))

    print('150 GHz')
    pixel_size = [3., 4.3, 5.3, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=2.0), ps))

    print('220 GHz')
    pixel_size = [3., 4.3, 5.3, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=1.35), ps))

    print('270 GHz')
    pixel_size = [3., 4.3, 5.3, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=1.1), ps))

    print('405 GHz')
    pixel_size = [3., 4.3, 5.8, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=0.74), ps))

    print('860 GHz')
    pixel_size = [1.35, 3., 4.3, 5.8, 6., 7., 7.8, 14.0]
    for ps in pixel_size:
        print('Edge taper is -{:.2f} dB for {:.1f} mm pixel size'.\
            format(pixel2taper(ps, taper_angle, wavel=0.3486), ps))


if __name__ == '__main__':

    main()