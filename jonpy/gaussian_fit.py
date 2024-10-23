import numpy as np
from scipy import optimize
import time

def eg1(p, x, y, fit4dc=False, verbose=False):

    # X offset
    beamx = p[0]
    # Y offset
    beamy = p[1]
    # rotation angle
    gam = p[2]
    # mean FWHM
    fwhm_mean = p[3]
    #ellipticity
    ellip = p[4]
    # amplitude
    A = p[5]
    # baseline offset from zero
    if fit4dc:
        offset = p[6]
    else:
        offset = 0.0

    p[2] = np.fmod(p[2], 180.0)
    sgam = np.sin(gam/180.0*np.pi)
    cgam = np.cos(gam/180.0*np.pi)

    fwhm_x = fwhm_mean/np.sqrt(ellip)
    sig_x = fwhm_x /np.sqrt(8*np.log(2))
    sig_y = sig_x*ellip

    u = x-beamx
    v = y-beamy

    xp = cgam*u+sgam*v
    yp = -sgam*u+cgam*v

    if verbose:
        print('Fit 1 returning an elliptical Gaussian beam:')
        print('Amplitude      : {:5.2f}'.format(p[5]))
        print('X centroid     : {:5.3f} deg'.format(p[0]))
        print('Y centroid     : {:5.3f} deg'.format(p[1]))
        print('FWHM           : {:5.3f} arcmin'.format(60*p[3]))
        print('Rotation angle : {:5.2f} deg'.format(p[2]))
        print('Eccentricity    : {:5.4f}'.format(p[4]))

    return A*np.exp(-(xp/sig_x)**2 /2 - (yp/sig_y)**2 / 2) + offset

def eg1_res(p, x, y, data, sigma, fit4dc=False):

    out = ((data-eg1(p, x, y, fit4dc=fit4dc))/sigma).flatten()

    return out

def eg1_res_curve_fit(xy, p):

    x, y = xy[0], xy[1]
    out = eg1(p, x, y).flatten()

    return out

def egfit1(Data, fit_radius=0.5, fit4dc=False, pguess=None, benchmark=False,
    force_center=False, verbose=False):

    x, y, data = Data['x'], Data['y'], Data['data']

    t0 = time.time()

    if fit_radius is not None:

        if (pguess is None) or force_center:
            r = np.sqrt(x**2 + y**2)
        else:
            r = np.sqrt((x-pguess[0])**2 + (y-pguess[1])**2)

        ifit = (r < fit_radius)
        if len(ifit) == 0:
            return 0.0, None, 0.0, 0.0

        t1 = time.time()
        sigma = np.std(data[np.where(r > fit_radius)])

        # print('Pguess:')
        # print(pguess)
        if fit4dc and pguess is None:
            pguess = [0.0, 0.0, 1.0, 0.5, 1.0, data[ifit].max(), 0.0]
        elif pguess is None:
            pguess = [0.0, 0.0, 1.0, 0.5, 1.0, data[ifit].max()]

        t2 = time.time()
        po, cov, infodict, errmsg, success = optimize.leastsq(eg1_res, pguess,
            args=(x[ifit], y[ifit], data[ifit], sigma, fit4dc),
            full_output=1)

        t3 = time.time()

        model = eg1(po, x, y)
        res = np.sum(np.abs(data[ifit]-model[ifit]))/len(model[ifit]) \
            /np.mean(data[ifit])

    else:

        t1 = time.time()
        sigma = np.std(data)
        if fit4dc and pguess is None:
            pguess = [0.0, -5.0, 1.0, 0.5, 1.0, data.max(), 0.0]
        elif pguess is None:
            pguess = [0.0, -5.0, 1.0, 0.5, 1.0, data.max()]
        t2 = time.time()

        # po, cov, infodict, errmsg, success = optimize.leastsq(eg1_res, pguess,
        #     args=(x, y, data, sigma, fit4dc), full_output=0, ftol=1e-6, xtol=1e-6)

        po, cov, infodict, errmsg, success = optimize.leastsq(eg1_res, pguess,
            args=(x, y, data, sigma, fit4dc), full_output=0)

        t3 = time.time()

        model = eg1(po,x,y)
        res = np.sum(np.abs(data-model))/float(len(model))/np.mean(data)

    if verbose:
        print('Fit 1 returning an elliptical Gaussian beam:')
        print('Amplitude      : {:5.2f}'.format(po[5]))
        print('X centroid     : {:5.3f} deg'.format(po[0]))
        print('Y centroid     : {:5.3f} deg'.format(po[1]))
        print('FWHM           : {:5.3f} arcmin'.format(60*po[3]))
        print('Rotation angle : {:5.2f} deg'.format(po[2]))
        print('Eccentricity    : {:5.4f}'.format(po[4]))

    if benchmark:
        # print(infodict)
        # print(errmsg)
        # print(success)
        print('egfit1 timing results:')
        for t in [t0, t1, t2, t3]:
            print('  t = {}'.format(t-t0))


    return po, cov, model, res

def eg2(xy, amplitude, x0, y0, sigma_x, sigma_y, theta):
    '''
    Create the most general 2D gaussian
    '''

    x = xy[0]
    y = xy[1]

    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0)
                            + c*((y-y0)**2)))
    return np.ravel(g)

def egfit2(X, Y, arr, maxfev=100000, fit_radius=0.5, pguess=None):

    '''
    Fit an elliptical Gaussian to a dataset
    '''

    print(np.shape(X))
    print(np.shape(Y))
    print(np.shape(arr))

    x, y = X[0,:], Y[0,:]
    if pguess is None:
        pguess = [np.abs(np.max(arr)-np.min(arr)), 0., 0., 0.5, 0.5, 0.]

    if fit_radius is not None:
        if pguess is None:
            R = np.sqrt(X**2 + Y**2)
        else:
            R = np.sqrt((X-pguess[1])**2 + (Y-pguess[2])**2)

        ifit = (R < fit_radius)

        popt, pcov = optimize.curve_fit(eg2, (X[ifit], Y[ifit]), arr[ifit], \
            p0=pguess, maxfev=maxfev)

        model = np.reshape(eg2((X,Y), *popt), np.shape(arr))
        res = np.sum(np.abs(arr[ifit]-model[ifit]))/float(len(model[ifit])) \
            /np.mean(arr[ifit])

    else:

        popt, pcov = optimize.curve_fit(eg2, (X.flatten(),Y.flatten()),
            arr.flatten(), p0=pguess, maxfev = maxfev)

        model = np.reshape(eg2((X,Y),*popt), np.shape(arr))
        res = np.sum(np.abs(arr-model))/float(len(model)) \
            /np.mean(arr)

    return popt, pcov, model, res

def egfit3():
    '''
    Code that was originally written in Matlab and I have now ported over to
    here.

    ToDo: Write this code

    '''

    pass

def centerofmass(arr):
    '''
    Estimates the center of mass and the width of the Gaussian distribution
    '''

    nx, ny = np.shape(arr)
    mx, my = np.sum(arr, axis=0), np.sum(arr, axis=1)
    mx = mx*(mx > 0).astype(float)
    my = my*(my > 0).astype(float)

    x, y = np.arange(nx), np.arange(ny)
    cx, cy = np.sum(mx*x)/np.sum(mx), np.sum(my*y)/np.sum(my)
    sx = np.sqrt(np.sum(mx*np.abs(x-cx)**2)/np.sum(mx))
    sy = np.sqrt(np.sum(my*np.abs(y-cy)**2)/np.sum(my))

    return cx, cy, sx, sy

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

    if step_size is not None:
        r *= step_size

    return np.unique(r.ravel()), radialprofile

# def center_map(cr, numel, arr):
#     '''
#     A function that we can use to center a map
#     '''

#     idxs = np.unravel_index(np.argmax(arr), np.shape(arr))
#     dists = np.array([idxs[0], numel[1]-idxs[0], idxs[1], numel[0]-idxs[1]])
#     sh = np.min(dists)

#     arr2use = arr[(idxs[0]-sh):(idxs[0]+sh), (idxs[1]-sh):(idxs[1]+sh)]
#     numel2use = [sh, sh, 0]
#     cr2use = []

#     return cr2use, numel2use, arr2use







