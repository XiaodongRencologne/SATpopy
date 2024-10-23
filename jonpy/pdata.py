import code, sys, os
import numpy as np
from pylab import *

def cldata(t,y):
    yout = interp_zeros(t,y)
    yout = np.ma.masked_where(np.isnan(yout),yout)
    #yout = np.ma.masked_invalid(np.isnan(yout),yout)
    return yout

def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y

def interp_zeros(t,y):
    badidx = np.where(y == 0)[0]
    ygood = np.delete(y,badidx)
    tgood = np.delete(t,badidx)
    yout = extrap(t,tgood,ygood)
    return yout

def interp_nans(t,y):
    badidx = np.where(np.isnan(y))[0]
    ygood = np.delete(y,badidx)
    tgood = np.delete(t,badidx)
    yout = extrap(t,tgood,ygood)
    return yout

def dwnsamp(vec,dsratio):
    idx = arange(len(vec))
    gidx = (mod(idx,dsratio) == 0)
    outvec = vec[gidx]
    return outvec

def minall(x):
    xx = np.hstack(x)
    return np.min(xx)

def maxall(x):
    xx = np.hstack(x)
    return np.max(xx)

def dspk_weight(win,stepsize=1.0):

    nclear = 1000
    wout = win
    dsell = np.diff(win)
    # Removing spike
    for i in np.arange(len(win)-1):
        if np.abs(dsell[i]) > stepsize:
            wout[(i+1):(nclear+i+1)] = np.nan

    wout = interp_nans(np.arange(len(wout)),wout)
    return wout

def shiftnorm(x):
    """
    Shifts tod so no negative values are present
    Normalizes to unity
    """
    x = x-np.min(x)
    x = x/np.max(x)
    return x

def mean_ns(x,n):
     s = np.std(x)
     m = np.median(x)
     goodidx = np.where(np.logical_and(x > (m-n*s),x < (m+n*s)))
     tmp = x[goodidx]
     return np.mean(tmp)

def meanwonan(x):
    badidx = np.where(np.isnan(x))[0]
    x = np.delete(x,badidx)
    return np.mean(x)

def mean_noedge(*args,**kwargs):

    if kwargs['p'] == []:
        p = 0.05
    else:
        p = kwargs['p']

    results = []
    for arg in args:
        idx = np.argsort(arg)
        l = len(idx)
        xs = np.round(p*l)
        tmpidx = idx[xs-1:l-(xs+1)]
        tmp = np.mean(arg[tmpidx])
        results.append(tmp)

    return results

def remnan(x):
    """
    Remove's nan elements from an array
    """
    badidx = np.where(np.isnan(x))[0]
    x = np.delete(x,badidx)
    return x

def minabs(x):
    """
    Shorthand version of matlab's very useful min(abs(x)).
    Nan elements are not included
    """
    badidx = np.where(np.isnan(x))[0]
    y = np.delete(x,badidx)
    maxval = np.max(y)
    x[badidx] = maxval+1
    return np.argmin(np.abs(x))

def remdc(tod,flip=False):
    """
    This function removes a minimum dc offset from a timestream. This
    can be useful when you have a timeline with an arbitrary offset
    and you want to produce a beam map.
    """
    if flip:
        tod = -tod

    mintod = np.min(tod)
    return tod
    #return (tod - mintod)


def flip(tod):
    tod = -tod
    return tod

def edge_cut(tod,de):
    tod_min = np.min(tod)
    tod_max = np.max(tod)
    badidx = np.where(np.logical_or(tod < (tod_min+de),tod > (tod_max-de)))
    return badidx

def interp_znans(t,y):
    yout = interp_zeros(t,y)
    yout = interp_nans(t,yout)
    return yout


def log10fp(x):
    badidx = np.where(x <= 0)[0]
    x[badidx] = np.nan
    return np.log10(x)

def log10abs(x):
    return np.log10(np.abs(x))

def norm_alize(d):
    d = interp_nans(np.arange(len(d)),d)
    goodidx = np.where(~np.isnan(d))
    maxd = np.max(d[goodidx])
    dout = d/maxd

    return dout

def interpmask(mask,v):
    vtmp = v
    vtmp[mask] = np.nan
    vout = interp_nans(np.arange(len(vtmp)),vtmp)
    return vout

def maxwonan(x):
    badidx = np.where(np.isnan(x))[0]
    x = np.delete(x,badidx)
    return np.max(x)

def minwonan(x):
    badidx = np.where(np.isnan(x))[0]
    x = np.delete(x,badidx)
    return np.min(x)

def remjump(t,v,start_times,end_times):
    """
    This function removes unwanted jumps in timestreams.
    This can be useful when you're trying to clean the weight
    timestream.
    """
    vout = v
    vchange = 0
    for i in np.arange(len(start_times)):
        idx1 = minabs(t-start_times[i])
        idx2 = minabs(t-end_times[i])
        vdiff = v[idx2]-v[idx1]
        vout[idx1:idx2] = np.nan
        vout[idx1+1:] = vout[idx1+1:]-vdiff
        print(vdiff)
        vchange -= vdiff

    print("Integrated change " + str(vchange))
    vout = interp_nans(np.arange(len(vout)),vout)
    return vout

def bindishit(bins,x,y):
    """
    Bins some shit
    """
    digi = np.digitize(x,bins)
    y_mean = np.array([y[digi == i].mean() for i in range(1, len(bins))])
    bdiff = np.mean(np.diff(bins),axis=0)
    b_m = bins + bdiff/2
    b_m = b_m[0:-1]
    return b_m, y_mean


def keyboard(banner=None):
    ''' Function that mimics the matlab keyboard command '''
    # use exception trick to pick up the current frame

    try:
        raise None
    except:
        frame = sys.exc_info()[2].tb_frame.f_back
    print("# Use quit() to exit :) Happy debugging!")
    # evaluate commands in current namespace

    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner, local=namespace)
    except SystemExit:
        return
