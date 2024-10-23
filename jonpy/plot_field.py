'''
Various beam-related plotting tools. I need to try to seperate

Some functions are at the borderline between belonging in plot_field and
beam_tools

'''

import sys, os
import matplotlib
import itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import patches
from matplotlib.collections import PatchCollection
from custom_cmap import *
import matplotlib.cm as cm
import numpy as np
import parse_grasp as pg
#from test_cmap import *
import scipy.interpolate
import gaussian_fit as gf
import pdata as pd
import beam_tools as bt

from pylab import rcParams
rcParams['font.size'] = 8


cmap = plt.get_cmap('viridis')
cmap_log = plt.get_cmap('magma')
cmap_grey = plt.get_cmap('Greys')

def configure_axis(ax, mult, numel, pix_step, unit_step1, unit_step2):

    plt.xticks(pix_step)
    plt.yticks(pix_step)
    ax.set_xticklabels(['{:5.2f}'.format(us) for us in unit_step1])
    ax.set_yticklabels(['{:5.2f}'.format(us) for us in unit_step2])

    plt.xlim([0, mult*numel[0]])
    plt.ylim([0, mult*numel[0]])

def configure_labels(ax=None, title=None, xlab=None, ylab=None,
        plot_displacement=False):

    if ax is None:

        if title is None:
            plt.title('Field amplitude')
        else:
            plt.title(title)

        if plot_displacement:
            plt.xlabel('Displacement [m]')
            plt.ylabel('Displacement [m]')

        if xlab is None or ylab is None:
            plt.xlabel('Az [deg]')
            plt.ylabel('El [deg]')
        else:
            plt.xlabel(xlab)
            plt.ylabel(ylab)

    else:

        if title is None:
            ax.set_title('Field amplitude')
        else:
            ax.set_title(title)

        if plot_displacement:
            ax.set_xlabel('Displacement [m]')
            ax.set_ylabel('Displacement [m]')

        if xlab is None or ylab is None:
            ax.set_xlabel('Az [deg]')
            ax.set_ylabel('El [deg]')
        else:
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)

def plot_ellipse_old(ax, numel, cr, angle = 0, xcenter=0., ycenter=42.5,
    width=50, height=30):

    dppix = (cr[2]-cr[0])/numel[0]

    theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
    x, y = 0.5 * width * np.cos(theta), 0.5 * height * np.sin(theta)
    rtheta = np.radians(angle)
    R = np.array([
        [np.cos(rtheta), -np.sin(rtheta)],
        [np.sin(rtheta),  np.cos(rtheta)],
        ])

    x, y = np.dot(R, np.array([x, y]))
    x += xcenter
    y += ycenter

    ax.fill(x, y, alpha=0.5, facecolor='white', zorder=3)

    e1 = patches.Ellipse((xcenter, ycenter), width, height,
        angle=angle, linewidth=1, color='white', fill=False, zorder=3)

    ax.add_patch(e1)

def plot_ellipse(ax, ellipse, color='white', ls='-'):

   # Ellipse properties
    ex, ey = ellipse[0], ellipse[1]
    # Half semi major and semi minor axis
    esx, esy = ellipse[2], ellipse[3]
    esth = ellipse[4]

    # c = np.cos(np.radians(180. - esth))
    # s = np.sin(np.radians(180. - esth))

    # # Translated
    # xc, yc = x - ex, y - ey

    # # Rotated
    # xct = xc*c - yc*s
    # yct = xc*s + yc*c

    e1 = patches.Ellipse((ex, ey), 2*esx, 2*esy,
        angle=esth, linewidth=1, color=color, ls=ls,
        fill=False, zorder=3)

    ax.add_patch(e1)

def plot_circle(ax, crad, cx=0, cy=0, color='white', ls='-', lw=1,
    calpha=0.2):

    from matplotlib.patches import Circle, RegularPolygon

    # circle = plt.Circle((cx, cy), crad,
    #         color=color, alpha=calpha,
    #         zorder=2, fill='white', ls=ls, lw=lw) # transform=ax.transData._b,

    circle = Circle([cx, cy], radius=crad,
            fill=False, edgecolor='white', lw=2, zorder=2,
            alpha=1.0)

    ax.add_artist(circle)

def plot_poly(ax,numel,cr,crad):

    npx = numel[0]
    secondary = np.array([[0.35*npx, 0.72*npx],
        [0.45*npx, 0.75*npx],
        [0.5*npx, 0.85*npx],
        [0.55*npx, 0.75*npx],
        [0.65*npx, 0.72*npx]])
    polygon = Polygon(secondary, True)

    p = PatchCollection([polygon], color='white', facecolor='white',
        zorder=3)

    ax.add_collection(p)

def interp_beam(x, y, cr, numel, arr):

    xx, yy, _ = bt.get_mesh(cr, numel)
    fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
    return fint(x, y)

# def turn_off_seaborn():

#     try:
#         import seaborn as sns
#         sns.reset_orig()
#     except ImportError:
#         print('WARNING: Unable to import seaborn.')

def plot_contour(cr, numel, arr, title=None, fname=None, dpi=150,
    vmin=None, vmax=None, crad=10, interp=True, mult=5,
    truncate=False, truncr=0.2, contlabel=True, levels=None,
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None):

    #turn_off_seaborn()

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    if log:
        cmap2u = cmap_log
    else:
        cmap2u = cmap

    plt.figure()
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:

        cmapn = 'plasma'
        cmap = plt.cm.get_cmap('plasma')
        maxv = 10*np.log10(arr2use).max()
        if levels is None:
            levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            #levels = [maxv-5, maxv-3, maxv-2, maxv-1]

            cmn = 5.0
            cs = plt.contour(10*np.log10(arr2use),
                colors=[cmap(0/cmn), cmap(1/cmn), cmap(2/cmn), cmap(3/cmn)],
                vmin=vmin, vmax=vmax, extent=extent, levels=levels,
                linewidths=1)
        else:

            levels = maxv+levels
            cmn = float(len(levels))+1.0
            cs = plt.contour(10*np.log10(arr2use),
                colors=[cmap(0/cmn), cmap(1/cmn), cmap(2/cmn), cmap(3/cmn)],
                vmin=vmin, vmax=vmax, extent=extent, levels=levels,
                linewidths=1)

    else:
        cs = plt.contour(arr2use, cmap=cmap, vmin=vmin, vmax=vmax,
            extent=extent)

    print('Contlabel = ', contlabel)

    if contlabel:
        plt.clabel(cs)

    if fname is None:
        fname2use = 'co'
    else:
        fname2use = fname

    if log:
        fname2use = fname+'_log'

    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.gca().set_aspect('equal')
    plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')
    plt.close()

def plot_beam(cr, numel, arr,
    ax=None,
    title=None, fname='beam', figsize=None,
    dpi=300,
    xlab=None, ylab=None,
    vmin=None, vmax=None, vminl=None, vmaxl=None, interp=True, mult=5,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False, show_cbar=True,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    aspect='equal', return_ax=False, tight_layout=False,
    plot_circ=False, crad=1, cx1=0, cx2=0, circle_color='white', plot_ellip=False,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None, fs=None,
    plot_line=False, xline=np.array([0., 0.]), yline=np.array([0., 0.]),
    anno_str=None, fmt="%5.2f", fmt_log="%5.0f", hexagon=False, transparent=False):
    '''
    Plots a beam
    '''

    extent = [cr[0], cr[2], cr[1], cr[3]]
    #extent = [cr[1], cr[3], cr[0], cr[2]]

    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        _, cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    ##turn_off_seaborn()

    if fs is not None:
        matplotlib.rcParams.update({'font.size':fs})
        matplotlib.rcParams.update({'xtick.labelsize':fs})
        matplotlib.rcParams.update({'ytick.labelsize':fs})

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    if ax is None:
        #plt.figure(figsize=figsize)
        fig=plt.figure(figsize=figsize)
        ax = plt.subplot(111)

    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            vmin=vminl, vmax=vmaxl, extent=extent)
        if show_cbar:

            # cax = fig.add_axes([ax.get_position().x1+0.01,
            #     ax.get_position().y0,0.02,ax.get_position().height])



            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)

            #cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')

    else:
        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent)
        if show_cbar:
            cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:
            maxv = 10*np.log10(arr2use).max()
            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            # c = plt.gca().contour(y, x, 10*np.log10(np.flipud(arr2use)),
            #     cmap=plt.get_cmap('Greys'), linewidths=lw, levels=levels,
            #     vmin=np.min(levels), vmax=np.max(levels))


            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=plt.get_cmap('Greys'), linewidths=lw, levels=levels,
                vmin=np.min(levels), vmax=np.max(levels))

            if contlabel:
                plt.clabel(c, fmt=fmt_log)
        else:

            if levels is not None:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5, levels=levels)
            else:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5)

            if contlabel:
                plt.clabel(c, fmt=fmt)

    if plot_pix:
        plot_pixels()

    if plot_circ:
        # plot_circle(plt.gca(), mult*numel, cr, crad, color=circle_color)
        #print('Plotting circle')
        #print(cx1, cx2, crad)
        # plot_circle(plt.gca(), cx1, cx2, crad, color=circle_color)
        print('Plotting circle123')
        plot_circle(ax, crad, cx2, cx1, color=circle_color)


    if plot_line:
        plt.plot(xline, yline, color='white', lw=1, ls=':')
        plt.plot(xline, -yline, color='white', lw=1, ls=':')

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), mult*numel, cr, ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    configure_labels(ax=ax, title=title, xlab=xlab, ylab=ylab)

    if hexagon:
        from matplotlib.patches import Circle, RegularPolygon
        hexagon = RegularPolygon([0, 0], 6, radius=0.217, orientation=np.radians(30.),
            fill=False, edgecolor='white', lw=2, zorder=2,
            alpha=1.0)

        ax.add_artist(hexagon)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if log:
        plt.clim([vminl, vmaxl])
    else:
        plt.clim([vmin, vmax])

    if tight_layout:
        fig = plt.gcf()
        fig.tight_layout()

    plt.draw()
    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight',
            transparent=transparent)

    if return_ax:
        ax = plt.gca()
        return ax
    else:
        plt.close()


def plot_beam_spill(cr, numel, arr,
    ax=None,
    title=None, fname='beam', figsize=None,
    dpi=150, fs=None,
    xlab=None, ylab=None,
    vmin=None, vmax=None, vminl=None, vmaxl=None, interp=True, mult=5,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False, show_cbar=True,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    aspect='equal',
    xticks=True, yticks=True,
    plot_circ=False, circle_color='white', lwc=1, plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None,
    anno_str=None, fmt="%5.2f", fmt_log="%5.0f", hexagon=False,
    transparent=False):
    '''
    Plots a beam
    '''

    if fs is not None:
        matplotlib.rcParams.update({'font.size':fs})
        matplotlib.rcParams.update({'xtick.labelsize':fs})
        matplotlib.rcParams.update({'ytick.labelsize':fs})

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        _, cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    ##turn_off_seaborn()

    if fs is not None:
        matplotlib.rcParams.update({'font.size':fs})
        matplotlib.rcParams.update({'xtick.labelsize':fs})
        matplotlib.rcParams.update({'ytick.labelsize':fs})

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.subplot(111)

    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        im = ax.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            vmin=vminl, vmax=vmaxl, extent=extent)
        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')

    else:
        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent)
        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:
            maxv = 10*np.log10(arr2use).max()
            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = ax.contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=plt.get_cmap('Greys'), linewidths=lw, levels=levels,
                vmin=np.min(levels), vmax=np.max(levels))

            if contlabel:
                plt.clabel(c, fmt=fmt_log)
        else:

            if levels is not None:
                c = ax.contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5, levels=levels)
            else:
                c = ax.contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5)

            if contlabel:
                plt.clabel(c, fmt=fmt)

    if plot_pix:
        plot_pixels()

    if plot_circ:
        # plot_circle(plt.gca(), mult*numel, cr, crad, color=circle_color)
        plot_circle(ax, crad, color=circle_color, lw=lwc)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(ax, mult*numel, cr, ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(ax, ellip,
                    color=col2use, ls=ls2use)

    configure_labels(ax=ax, title=title, xlab=xlab, ylab=ylab)

    if hexagon:
        from matplotlib.patches import Circle, RegularPolygon
        hexagon = RegularPolygon([0, 0], 6, radius=0.217,
            fill=False, edgecolor='white', lw=2, zorder=2,
            alpha=1.0)

        ax.add_artist(hexagon)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if add_grids:
        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if not xticks:
        ax.set_xticklabels([])

    if not yticks:
        ax.set_yticklabels([])

    if log:
        im.set_clim([vminl, vmaxl])
    else:
        im.set_clim([vmin, vmax])

    if save:
        plt.draw()
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight',
            transparent=transparent)
        plt.close()

    return im


def plot_ground_beam(cr, numel, arr, title=None, fname='beam', figsize=None,
    dpi=150,
    xlab=None, ylab=None,
    vmin=None, vmax=None, vminl=None, vmaxl=None, interp=True, mult=5,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False, show_cbar=True,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    aspect='equal',
    plot_circ=False, circle_color='white', plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None, fs=None,
    anno_str=None, fmt="%5.2f", fmt_log="%5.0f"):
    '''
    Plots a beam
    '''

    ground_angle = np.degrees(np.arctan(2.0/9.0))
    ground_screen_angle = np.degrees(np.arctan(4.5/9.0))


    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    #turn_off_seaborn()

    if fs is not None:
        matplotlib.rcParams.update({'font.size':fs})
        matplotlib.rcParams.update({'xtick.labelsize':fs})
        matplotlib.rcParams.update({'ytick.labelsize':fs})


    if not os.path.exists(imgd):
        os.makedirs(imgd)

    plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:


        xx, yy, dA = bt.get_mesh(cr,numel)
        R = np.sqrt(xx**2 + yy**2)

        badidx = (R > 90+ground_angle)
        faintidx = (R < 90+ground_angle) &  (R > 90-ground_screen_angle)

        arr2use[badidx] = np.nan
        im2 = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            vmin=vminl, vmax=vmaxl, extent=extent, alpha=0.7)

        arr2use[faintidx] = np.nan
        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            vmin=vminl, vmax=vmaxl, extent=extent, alpha=1.0)


        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')

    else:
        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent)
        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:
            maxv = 10*np.log10(arr2use).max()
            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=cmap_grey, linewidths=lw, levels=levels)

            if contlabel:
                plt.clabel(c, fmt=fmt_log)
        else:

            if levels is not None:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5, levels=levels)
            else:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5)

            if contlabel:
                plt.clabel(c, fmt=fmt)

    if plot_pix:
        plot_pixels()

    if plot_circ:
        # plot_circle(plt.gca(), mult*numel, cr, crad, color=circle_color)
        plot_circle(plt.gca(), crad, color=circle_color)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), mult*numel, cr, ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    configure_labels(title=title, xlab=xlab, ylab=ylab)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if log:
        plt.clim([vminl, vmaxl])
    else:
        plt.clim([vmin, vmax])

    plt.draw()

    plt.yticks([-90, -60, -30, 0, 30, 60, 90], [180, 150, 120, 90, 60, 30, 0])


    plot_circle(plt.gca(), 90+ground_angle, color=circle_color, ls=':')
    plot_circle(plt.gca(), 90-26.6, color=circle_color, ls=':')

    plot_circle(plt.gca(), 50, color=circle_color, ls='--')
    plot_circle(plt.gca(), 40, color=circle_color, ls='--')
    plot_circle(plt.gca(), 30, color=circle_color, ls='--')

    plt.annotate('60', xy=(30, 0), color='white', size=8)
    plt.annotate('50', xy=(40, 5), color='white', size=8)
    plt.annotate('40', xy=(50, 10), color='white', size=8)
    plt.annotate('Ground screen ', xy=(-22, 80), color='white', size=8)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')

    plt.close()


def plot_taurus_ground_beam(cr, numel, arr, title=None, fname='beam', figsize=None,
    dpi=150,
    xlab=None, ylab=None,
    vmin=None, vmax=None, vminl=None, vmaxl=None, interp=True, mult=5,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False, show_cbar=True,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    aspect='equal',
    plot_circ=False, circle_color='white', plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None, fs=None,
    anno_str=None, fmt="%5.2f", fmt_log="%5.0f"):
    '''
    Plots a beam
    '''

    ground_angle = np.degrees(np.arctan(2.0/9.0))
    ground_screen_angle = np.degrees(np.arctan(4.5/9.0))


    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    #turn_off_seaborn()

    if fs is not None:
        matplotlib.rcParams.update({'font.size':fs})
        matplotlib.rcParams.update({'xtick.labelsize':fs})
        matplotlib.rcParams.update({'ytick.labelsize':fs})


    if not os.path.exists(imgd):
        os.makedirs(imgd)

    plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:


        xx, yy, dA = bt.get_mesh(cr,numel)
        R = np.sqrt(xx**2 + yy**2)

        badidx = (R > 90+ground_angle)
        faintidx = (R < 90+ground_angle) &  (R > 90-ground_screen_angle)

        arr2use[badidx] = np.nan
        # im2 = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            # vmin=vminl, vmax=vmaxl, extent=extent, alpha=0.7)

        # arr2use[faintidx] = np.nan
        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            vmin=vminl, vmax=vmaxl, extent=extent, alpha=1.0)

        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')

    else:
        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent)
        if show_cbar:
            cbar = plt.colorbar(im)
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:
            maxv = 10*np.log10(arr2use).max()
            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=cmap_grey, linewidths=lw, levels=levels)

            if contlabel:
                plt.clabel(c, fmt=fmt_log)
        else:

            if levels is not None:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5, levels=levels)
            else:
                c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                     linewidths=0.5)

            if contlabel:
                plt.clabel(c, fmt=fmt)

    if plot_pix:
        plot_pixels()

    if plot_circ:
        # plot_circle(plt.gca(), mult*numel, cr, crad, color=circle_color)
        plot_circle(plt.gca(), crad, color=circle_color)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), mult*numel, cr, ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    configure_labels(title=title, xlab=xlab, ylab=ylab)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if log:
        plt.clim([vminl, vmaxl])
    else:
        plt.clim([vmin, vmax])

    plt.draw()

    plt.yticks([-90, -60, -30, 0, 30, 60, 90], [180, 150, 120, 90, 60, 30, 0])
    # plot_circle(plt.gca(), 90+ground_angle, color=circle_color, ls=':')
    # plot_circle(plt.gca(), 90-26.6, color=circle_color, ls=':')

    plot_circle(plt.gca(), 90-49, color=circle_color, ls='--')
    plot_circle(plt.gca(), 90-35, color=circle_color, ls='--')
    plot_circle(plt.gca(), 90-21, color=circle_color, ls='--')

    plt.annotate('49', xy=(30, 0), color='white', size=8)
    plt.annotate('35', xy=(45, 5), color='white', size=8)
    plt.annotate('21', xy=(60, 10), color='white', size=8)
    # plt.annotate('Ground screen ', xy=(-22, 80), color='white', size=8)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')

    plt.close()


def merge_masks(masks):
    maski = masks[0]
    for j in range(len(masks)):
        maski = np.logical_or(maski, masks[j])

    return maski

def mask_patch(patch, xx, yy):
    '''
    Receives two 2d vectors, xx and yy. Returns a binary mask with the same
    shape. The mask is true when the x,y coordinate is inside the borders,
    False ioc.
    '''
    XY = np.dstack((xx,yy))
    XY_flat = XY.reshape((-1,2))

    # print(patch.get_path())
    # print(patch.get_width())
    # print(patch.get_xy())

    width = patch.get_width()

    x, y = patch.get_xy()

    vertices = [[x, y], [x, y+width], [x+width, y+width], [x+width, y]]
    path = matplotlib.path.Path(vertices)

    # masked = patch.get_path().contains_points(XY_flat)
    masked = path.contains_points(XY_flat)
    masked = masked.reshape(np.shape(xx))

    return masked


def plot_beam2(cr, numel, arr, title=None, fname='beam', figsize=None,
    dpi=150, xlab=None, ylab=None, vmin=None, vmax=None,
    vminl=None, vmaxl=None, interp=True, mult=5, logrange=None, aspect='equal',
    scatter_plot=False, xs=None, ys=None, scatter_arr=None,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, log_label='dB', 
    imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    plot_circ=False, calpha=0.1,
    circle_color='white', circle_ls='-', plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None,
    plot_squares=False, square_width=87.5,
    xsq=[-2, -1, 0, 1, -2, -1, 0, 1],
    ysq = [-1, -1, -1, -1, 0, 0, 0, 0],
    anno_str=None, anno_strx=0.03, anno_stry=0.95):
    '''
    Plots a beam
    '''

    xsq = square_width * np.array(xsq)
    ysq = square_width * np.array(ysq)

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)

    if truncate:
        _, cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    #turn_off_seaborn()

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        maxv = 10*np.log10(arr2use).max()
        if logrange is not None:
            vmaxl = maxv
            vminl = maxv - logrange

        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect=aspect,
            norm=plt.Normalize(vmin=vminl, vmax=vmaxl, clip=True),
            extent=extent)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)   
        cbar = plt.colorbar(im, cax=cax)

        #cbar = plt.colorbar(im)
        cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel(log_label)

        plt.sca(ax)

        # cbar = plt.colorbar(cmap=cmapl)
        # cbar = plt.colorbar(sm, boundaries=[vminl, vmaxl],
        #     ticks=[-70, -60, -50, -40, -30, -20, -10])

    else:

        if vmin is None:
            vmin = np.nanmin(arr2use)

        if vmax is None:
            vmax = np.nanmax(arr2use)

        im = plt.imshow(arr2use, cmap=cmap, aspect=aspect, vmin=vmin, vmax=vmax,
            extent=extent)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)

        plt.sca(ax)

        if clab is not None and len(clab) > 7:
            cbar.ax.set_ylabel(clab)

        else:
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:

            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=cmap_grey, linewidths=lw, levels=levels)

            if contlabel:
                plt.clabel(c, fmt="%5.0f")
        else:
            c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                 linewidths=0.5, origin='lower', extent=extent)

            if contlabel:
                plt.clabel(c, fmt="%5.2f")

    if plot_pix:
        plot_pixels()

    if plot_circ:
        plot_circle(plt.gca(), crad, color=circle_color,
            ls=circle_ls, calpha=calpha)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    if plot_squares:

        squares = []
        for xsqi, ysqi in zip(xsq, ysq):
            # square_verts.append(fp.squareVert([xsqi, ysqi], square_width))
            square = patches.Rectangle([xsqi, ysqi], square_width, square_width,
                edgecolor='white', linewidth=2, facecolor='None', alpha=0.2)

            squares.append(square)
            ax.add_artist(square)

        N = len(xsq)
        square_mask = merge_masks(map(mask_patch,
            squares, itertools.repeat(xx, N), itertools.repeat(yy, N)))

    configure_labels(ax=ax, title=title, xlab=xlab, ylab=ylab)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(anno_strx, anno_stry), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'

    if scatter_plot:

        if log:
            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmapl, edgecolor='black',
                vmin=np.nanmin(10*np.log10(arr2use)),
                vmax=np.nanmax(10*np.log10(arr2use)))
        else:
            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmap, edgecolor='black', vmin=vmin, vmax=vmax)

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')

    plt.close()

def plot_beam2_squares(cr, numel, arr, title=None, fname='beam', figsize=None,
    dpi=150, xlab=None, ylab=None,
    vmin=None, vmax=None, vminl=None, vmaxl=None, interp=True, mult=5,
    logrange=None,
    scatter_plot=False, xs=None, ys=None, scatter_arr=None,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False,
    cmapl=plt.get_cmap('magma'), cmap=plt.get_cmap('viridis'),
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None, relative_levels=True,
    add_contour=True, contlabel=True, levels=None, lw=0.5, add_grids=False,
    plot_circ=False, circle_color='white', circle_ls='-', plot_ellip=False, crad=1,
    ellipse=[0, 0, 1, 1, 0], ellipse_pars=None,
    alpha1=1.0, alpha2=1.0, alpha_squares=0.5,
    plot_squares=False, square_width=95, #87.5
    xsq=[-2, -1, 0, 1, -2, -1, 0, 1],
    ysq = [-1, -1, -1, -1, 0, 0, 0, 0],
    anno_str=None):
    '''
    Plots a beam
    '''

    xsq = square_width * np.array(xsq)
    ysq = square_width * np.array(ysq)

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    #turn_off_seaborn()

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr, numel)
        xx2, yy2 = bt.get_lin(cr, numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        maxv = 10*np.log10(np.nanmax(arr2use))
        print(maxv)
        if logrange is not None:
            vmaxl = maxv
            vminl = maxv - logrange

        im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect='equal',
            norm=plt.Normalize(vmin=vminl, vmax=vmaxl, clip=True),
            extent=extent, alpha=alpha1)

        cbar = plt.colorbar(im)
        cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('dB')
        # cbar = plt.colorbar(cmap=cmapl)
        # cbar = plt.colorbar(sm, boundaries=[vminl, vmaxl],
        #     ticks=[-70, -60, -50, -40, -30, -20, -10])

    else:

        if vmin is None:
            vmin = np.nanmin(arr2use)

        if vmax is None:
            vmax = np.nanmax(arr2use)

        im = plt.imshow(arr2use, cmap=cmap, aspect='equal', vmin=vmin, vmax=vmax,
            extent=extent, alpha=alpha1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)

        plt.sca(ax)

        if len(clab) > 7:
            cbar.ax.set_ylabel(clab)

        else:
            cbar.ax.set_xlabel(clab) if clab is not None else cbar.ax.set_xlabel('Power')

    if add_contour:

        x, y = bt.get_lin(cr, numel, mult=mult)
        if log:

            if levels is None:
                levels = [maxv-20, maxv-13, maxv-10, maxv-3]
            else:
                levels = levels
                if relative_levels:
                    levels += maxv

            c = plt.gca().contour(x, y, 10*np.log10(np.flipud(arr2use)),
                cmap=cmap_grey, linewidths=lw, levels=levels)

            if contlabel:
                plt.clabel(c, fmt="%5.0f")
        else:
            c = plt.gca().contour(x, y, np.flipud(arr2use), cmap=cmap_grey,
                 linewidths=0.5, origin='lower', extent=extent,
                 levels=[0.005, 0.01, 0.02, 0.05, 0.1])

            if contlabel:
                plt.clabel(c, fmt="%5.3f")

    if plot_pix:
        plot_pixels()

    if plot_circ:
        plot_circle(plt.gca(), crad, color=circle_color,
            ls=circle_ls)

    if plot_ellip:
        if len(np.shape(ellipse)) == 1:

            if ellipse_pars is not None:
                col2use = ellipse_pars['color']
                ls2use = ellipse_pars['ls']
            else:
                col2use, ls2use = 'white', '-'

            plot_ellipse(plt.gca(), ellipse,
                color=col2use, ls=ls2use)

        else:
            for ei, ellip in enumerate(ellipse):
                if ellipse_pars is not None:
                    col2use = ellipse_pars['color'][ei]
                    ls2use = ellipse_pars['ls'][ei]
                else:
                    col2use, ls2use = 'white', '-'

                plot_ellipse(plt.gca(), ellip,
                    color=col2use, ls=ls2use)

    if plot_squares:

        # square_verts = []
        squares = []
        for xsqi, ysqi in zip(xsq, ysq):
            # square_verts.append(fp.squareVert([xsqi, ysqi], square_width))
            square = patches.Rectangle([xsqi, ysqi], square_width, square_width,
                edgecolor='white', linewidth=2, ls=':', facecolor='None', alpha=alpha_squares)

            squares.append(square)
            ax.add_artist(square)

        N = len(xsq)
        square_mask = merge_masks(map(mask_patch,
            squares, itertools.repeat(xx, N), itertools.repeat(yy, N)))

        arr2use[~square_mask] = np.nan

        # if log:
        #     im = plt.imshow(10*np.log10(arr2use), cmap=cmapl, aspect='equal',
        #         norm=plt.Normalize(vmin=vminl, vmax=vmaxl, clip=True),
        #         extent=extent, alpha=alpha2)

        # else:

        #     im = plt.imshow(arr2use, cmap=cmap, aspect='equal', vmin=vmin, vmax=vmax,
        #         extent=extent, alpha=alpha2)


    configure_labels(title=title, xlab=xlab, ylab=ylab)

    if anno_str is not None:
        plt.annotate(anno_str, xy=(0.03, 0.95), xycoords='axes fraction',
            color='white')

    fname2use = 'co' if fname is None else fname
    if log:
        fname2use = fname+'_log'


    if scatter_plot:

        if log:
            if logrange is not None:
                vmaxl = maxv
                vminl = maxv - logrange


            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmapl, edgecolor='black',
                vmin=vminl,#np.nanmin(10*np.log10(arr2use)),
                vmax=vmaxl)#np.nanmax(10*np.log10(arr2use)))
        else:
            plt.scatter(xs, ys, c=scatter_arr, s=50, marker='h',
                cmap=cmap, edgecolor='black', vmin=vmin, vmax=vmax)

    if add_grids:
        ax = plt.gca()

        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        ax.grid(color='w', linestyle='-', linewidth=1, alpha=0.2)
        ax.set_xticks(np.arange(np.round(xlim[0], decimals=1),
            np.round(xlim[1], decimals=1), 0.1), minor=True)

        ax.set_xticks(np.arange(np.round(ylim[0], decimals=1),
            np.round(ylim[1], decimals=1), 0.1), minor=True)

        #ax.set_yticks(np.arange(-3, 3.1, 0.1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=0.1, alpha=0.3)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi, bbox_inches='tight')

    plt.close()

    return arr2use, square_mask


def plot_surface(cr, numel, arr, title=None, fname='beam', dpi=150,
    vmin=None, vmax=None, vminl=None, vmaxl=None, crad=10, interp=True, mult=5,
    truncate=False, truncr=0.2, clab=None, save=True, plot_pix=False,
    log=True, imgd='img/', nticks=5, xlim=None, ylim=None,
    add_contour=True, levels=None, lw=0.5, show=False):

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    xx0, yy0 = xx.copy(), yy.copy()

    #turn_off_seaborn()

    if not os.path.exists(imgd):
        os.makedirs(imgd)

    if log:
        cmap2u = cmap_log
    else:
        cmap2u = cmap

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    if not interp:
        mult = 1
        arr2use = arr.copy()
    else:
        xx, yy, dA = bt.get_mesh(cr,numel)
        xx2, yy2 = bt.get_lin(cr,numel, mult=mult)

        fint = scipy.interpolate.RectBivariateSpline(xx[0,:], yy[:,0], arr)
        arr2use = fint(xx2, yy2)
        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

    if log:
        ax.plot_surface(xx, yy, 10*np.log10(arr2use), cmap=cmap_log,
            vmin=vminl, vmax=vmaxl)

    else:
        ax.plot_surface(xx, yy, arr2use, cmap=cmap, vmin=vmin, vmax=vmax)

    configure_labels(title=title)

    fname2use = 'co' if fname is None else fname

    if log:
        fname2use = fname+'_log'

    plt.xlim(xlim)
    plt.ylim(ylim)

    if save:
        plt.savefig(imgd+fname2use+'.png', dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close()

def plot_pixels():

    #pix_step = np.linspace(0, mult*numel[0], nticks)
    #unit_step1 = np.linspace(cr[0], cr[2], nticks)
    #unit_step2 = np.linspace(cr[1], cr[3], nticks)
    #configure_axis(ax, mult, numel, pix_step, unit_step1, unit_step2)

    xo = 0.03
    xf1 = np.array([ 0.068, 0.136, 0.204, 0.272, 0.34, 0.408, 0.442,
        0.476, 0.51, 0.544, 0.578, 0.612])
    yf1 = np.array([ 0.068, 0.136, 0.204, 0.272, 0.34, 0.408, 0.442, 0.476,
        0.51, 0.544, 0.578, 0.612])

    xf2 = np.array([ 0.0748, 0.1496, 0.2244, 0.2992, 0.374, 0.4488, 0.5236,
        0.5984, 0.6732, 0.8228, 0.8976])
    yf2 = np.array([ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

    xf3 = np.array([ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    yf3 = np.array([ 0.0748, 0.1496, 0.2244, 0.2992, 0.374, 0.4488, 0.5236,
        0.5984, 0.6732, 0.748, 0.8228])

    plt.plot(xo+xf1,yf1, lw=0, marker='.',ms=8, label='D', color='#4c72b0')
    plt.plot(xo+xf2,yf2, lw=0, marker='.',ms=8, label='X', color='#55a868')
    plt.plot(xo+xf3,yf3, lw=0, marker='.',ms=8, label='Y', color='#c44e52')
    plt.legend(loc=2, numpoints=1)

def plot_grasp(fpath, nfreq=1, pname='beam', save=True, near_field=False,
    cocx=False, comcx=False, contour=False, verbose=True):

    if verbose:
        print('Parsing data...')

    if nfreq == 1:
        cr, numel, data = pg.read_grasp_grid(fpath)

    else:
        cr, numel, data = pg.read_grasp_grids(fpath,
            nfreq=nfreq, debug=False)

    if near_field:
        arr1, arr2, arr3 = pg.parse_data(data, numel, near_field=near_field,
            nfreq=nfreq)
    else:
        arr1, arr2 = pg.parse_data(data, numel, near_field=near_field,
            nfreq=nfreq)

    if verbose:
        print('Done parsing, now plotting...')

    if contour:
        pfunc = plot_contour
    else:
        pfunc = plot_beam

    if comcx:
        pfunc(cr, numel, np.abs(arr1-arr2), fname=pname+'_comcx', vmin=-30, vmax=0)

    pfunc(cr, numel, arr1, fname=pname+'_co')
    pfunc(cr, numel, arr2, fname=pname+'_cx')

    if near_field:
        pfunc(cr, numel, arr3, fname=pname+'_cr')

    if cocx:
        pfunc(cr, numel, arr1+arr2, fname=pname+'_cocx')

def plot2d(cr, numel, arr, title=None, fname=None, dpi=150,
    imgd='img/', vmin=None, vmax=None, nticks=5, log=True,
    vminl=None, vmaxl=None, truncate=False, truncr=0.2, crad=10, interp=False,
    mult=5, plot_displacement=False, fit_arr1=None, fit_arr2=None,
    plot_shapes=False):
    '''
    Plotting function with a number of optional arguments. Originally designed
    to be able to accommodate data and the corresponding fit. I need to fix this.
    '''

    # Getting rid of seaborn configs if they are there
    #turn_off_seaborn()

    extent = [cr[0], cr[2], cr[1], cr[3]]
    xx, yy, dA = bt.get_mesh(cr, numel)
    xx0, yy0 = xx.copy(), yy.copy()

    if truncate:
        cr, numel, xx, yy, arr = bt.trunc_map(xx, yy, arr, ratio=truncr)
        extent = [cr[0], cr[2], cr[1], cr[3]]

    if log:
        cmap2u = cmap_log
    else:
        cmap2u = cmap

    plt.figure()
    ax = plt.subplot(111)
    if not interp:
        mult = 1
        arr2use = arr.copy()
        if fit_arr1 is not None:
            fit_arr2use = fit_arr1.copy()
            if truncate:
                _, _, _, _, fit_arr2use = bt.trunc_map(xx0, yy0,
                    fit_arr2use, ratio=truncr)
        if fit_arr2 is not None:
            fit_arr2use2 = fit_arr2.copy()
            if truncate:
                _, _, _, _, fit_arr2use2 = bt.trunc_map(xx0, yy0,
                    fit_arr2use2, ratio=truncr)

    else:

        xx2, yy2 = bt.get_lin(cr, numel, mult=mult)
        fint = scipy.interpolate.RectBivariateSpline(xx[0,:],yy[:,0],arr)

        arr2use = fint(xx2,yy2)
        if fit_arr1 is not None:
            if truncate:
                _, _, _, _, fit_arr1 = bt.trunc_map(xx0, yy0, fit_arr1,
                    ratio=truncr)

            fint_fit = scipy.interpolate.RectBivariateSpline(xx[0,:],
                yy[:,0], fit_arr1)
            fit_arr2use = fint_fit(xx2,yy2)

        if fit_arr2 is not None:
            if truncate:
                _, _, _, _, fit_arr2 = bt.trunc_map(xx0, yy0, fit_arr2,
                    ratio=truncr)

            fint_fit2 = scipy.interpolate.RectBivariateSpline(xx[0,:],
                yy[:,0], fit_arr2)

            fit_arr2use2 = fint_fit2(xx2,yy2)

        xx, yy, dA = bt.get_mesh(cr, numel, mult=mult)

        if log:
            i1 = (arr2use <= 0)
            arr2use[i1] = np.min(arr2use[~i1])

            if fit_arr1 is not None:
                i2 = (fit_arr2use <= 0)
                fit_arr2use[i2] = np.min(fit_arr2use[~i2])
            if fit_arr2 is not None:
                i3 = (fit_arr2use2 <= 0)
                fit_arr2use2[i3] = np.min(fit_arr2use2[~i3])

    if log:

        plt.imshow(10*np.log10(arr2use), cmap=cmap_log,aspect='equal',
            vmin=vminl, vmax=vmaxl, extent=extent)
        cbar = plt.colorbar()
        cbar.ax.set_xlabel('dB')
    else:
        plt.imshow(arr2use, cmap=cmap, aspect='equal', vmin=vmin, vmax=vmax,
            extent=extent)
        cbar = plt.colorbar()
        cbar.ax.set_xlabel('Power')

    if plot_shapes:
        #plot_poly(ax,numel,cr,crad)
        plot_circle(ax,mult*numel,cr,crad)
        plot_ellipse(ax,mult*numel,cr)
        power_atr(ax, xx, yy, arr2use, cr, mult*numel, crad+0.7)

    pix_step = np.linspace(0, mult*numel[0], nticks)
    unit_step1 = np.linspace(cr[0], cr[2], nticks)
    unit_step2 = np.linspace(cr[1], cr[3], nticks)

    if plot_displacement:
        R = 7.8
        unit_step1 = R*np.sin(np.radians(unit_step1))
        unit_step2 = R*np.sin(np.radians(unit_step2))

    configure_labels(title=title, plot_displacement=plot_displacement)

    if fname is None:
        fname2use = 'co'
    else:
        fname2use = fname
    if log:
        fname2use = fname+'_log'
    if plot_displacement:
        fname2use = fname+'_displacement'

    plt.axis(extent)
    plt.savefig(imgd+fname2use+'.png', dpi=dpi)
    plt.close()

    ## Plotting the fit array as well
    if fit_arr1 is not None:

        plt.figure()
        ax = plt.subplot(111)

        if log:
            plt.imshow(10*np.log10(fit_arr2use), cmap=cmap_log,
                aspect='equal', vmin=vminl, vmax=vmaxl, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('dB')
        else:
            plt.imshow(fit_arr2use, cmap=cmap, aspect='equal',
                vmin=vmin, vmax=vmax, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('Power')

        if plot_shapes:
            #plot_poly(ax,numel,cr,crad)
            plot_circle(ax,mult*numel,cr,crad)
            plot_ellipse(ax,mult*numel,cr)
            power_atr(ax, xx, yy, arr2use, cr, mult*numel, crad+0.7)

        if title is not None:
            configure_labels(title=title+' - Fit',
                plot_displacement=plot_displacement)

        if fname is None:
            fname2use = 'co'
        fname2use = fname+'_fit'
        if log:
            fname2use += '_log'
        if plot_displacement:
            fname2use += '_displacement'

        plt.savefig(imgd+fname2use+'.png', dpi=dpi)
        plt.close()

        # Plotting the difference map
        if not log:
            plt.figure()
            ax = plt.subplot(111)

            plt.imshow((arr2use-fit_arr2use)/arr2use.max(),
                cmap=cmap2u, aspect='equal',
                vmin=vmin, vmax=vmax, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('Diff')

            if plot_shapes:
                #plot_poly(ax,numel,cr,crad)
                plot_circle(ax,mult*numel,cr,crad)
                plot_ellipse(ax,mult*numel,cr)
                power_atr(ax,xx,yy,arr2use,cr,mult*numel,crad+0.7)

            if title is not None:
                configure_labels(title=title+' - Difference',
                    plot_displacement=plot_displacement)

            if fname is None:
                fname2use = 'co'
            fname2use = fname+'_diff'
            if plot_displacement:
                fname2use += '_displacement'

            plt.savefig(imgd+fname2use+'.png', dpi=dpi)
            plt.close()

    ## Plotting the 2nd fit array as well
    if fit_arr2 is not None:

        plt.figure()
        ax = plt.subplot(111)

        if log:
            plt.imshow(10*np.log10(fit_arr2use2), cmap=cmap_log,
                aspect='equal', vmin=vminl, vmax=vmaxl, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('dB')
        else:
            plt.imshow(fit_arr2use2, cmap=cmap, aspect='equal',
                vmin=vmin, vmax=vmax, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('Power')

        if plot_shapes:
            #plot_poly(ax,numel,cr,crad)
            plot_circle(ax,mult*numel,cr,crad)
            plot_ellipse(ax,mult*numel,cr)
            power_atr(ax, xx, yy, arr2use, cr, mult*numel,crad+0.7)

        configure_labels(title=title+' - Fit2',
            plot_displacement=plot_displacement)

        if fname is None:
            fname = 'co'
        fname2use = fname+'_fit2'
        if log:
            fname2use += '_log'
        if plot_displacement:
            fname2use += '_displacement'

        plt.savefig(imgd+fname2use+'.png', dpi=dpi)
        plt.close()

        # Plotting the difference map
        if not log:
            plt.figure()
            ax = plt.subplot(111)

            plt.imshow(arr2use-fit_arr2use2, cmap=cmap, aspect='equal',
                vmin=vmin, vmax=vmax, extent=extent)
            cbar = plt.colorbar()
            cbar.ax.set_xlabel('Diff')

            if plot_shapes:
                #plot_poly(ax,numel,cr,crad)
                plot_circle(ax,mult*numel,cr,crad)
                plot_ellipse(ax,mult*numel,cr)
                power_atr(ax, xx, yy, arr2use, cr, mult*numel, crad+0.7)

            configure_labels(title=title+' - Difference2',
                plot_displacement=plot_displacement)

            if fname is None:
                fname2use = 'co'
            fname2use = fname+'_diff2'
            if plot_displacement:
                fname2use += '_displacement2'

            plt.savefig(imgd+fname2use+'.png', dpi=dpi)
            plt.close()

def plot_profile(cr, numel, arr,
    title=None, fname=None, estr='', rmax=None, label='Data', color=None,
    interp=False, mult= 5, mark_fwhm=False, dashes=[], units='deg',
    fit_arr1=None, fit_arr2=None, log=True, imgd='img/', alpha=1.0,
    vminl=None, vmaxl=None, dpi=150, finish=True, fg=False, normalize=False,
    parse_only=False, zorder=None):
    '''
    Insert text
    '''

    #turn_off_seaborn()

    r, profile = bt.get_profile(cr, numel, arr, fg=fg)
    if units=='arcmin':
        r *= 60

    if parse_only:
        profile /= profile.max()
        return r, profile

    if interp:

        r2u = np.linspace(np.min(r), np.max(r), mult*len(r))
        f2u = scipy.interpolate.interp1d(r, profile, fill_value='extrapolate')
        pr2u = f2u(r2u)

    else:
        r2u, pr2u = r, profile

    if normalize:
        pr2u /= pr2u.max()

    if log:
        plt.plot(r2u, 10*np.log10(pr2u), label=label, color=color, alpha=alpha,
            dashes=dashes, zorder=zorder)
    else:
        plt.plot(r2u, pr2u, label=label, color=color, alpha=alpha,
            dashes=dashes, zorder=zorder)

    if fit_arr1 is not None:
        rf1, profilef1 = bt.get_profile(cr, numel, fit_arr1, fg=fg)
        if units=='arcmin':
            rf1 *= 60

        if normalize:
            profilef1 /= profilef1.max()
        if log:
            plt.plot(rf1, 10*np.log10(profilef1), color=color,dashes=[4, 4],
                alpha=alpha, label='Fit', zorder=zorder)
        else:
            plt.plot(rf1, profilef1, color=color,dashes=[4, 4],
                alpha=alpha, label='Fit', zorder=zorder)

    if fit_arr2 is not None:
        rf2, profilef2 = bt.get_profile(cr, numel, fit_arr2, fg=fg)
        if units=='arcmin':
            rf2 *= 60

        if normalize:
            profilef2 /= profilef2.max()

        if log:
            plt.plot(rf2, 10*np.log10(profilef2), dashes=[8, 8], label='Fit 2',
                alpha=alpha, zorder=zorder)
        else:
            plt.plot(rf2, profilef2, dashes=[8, 8], label='Fit 2',
                alpha=alpha, zorder=zorder)

    if finish:

        if mark_fwhm:
            idx = np.argsort(np.abs(pr2u-np.max(pr2u)/2.))[0]
            plt.axvline(x=r2u[idx], color='black', alpha=0.5, dashes=[8, 8])

            print('Data FWHM = {:5.2f} arcmin'.format(2*r2u[idx]))

            if fit_arr1 is not None:
                idx = np.argsort(np.abs(profilef1-np.max(profilef1)/2.))[0]
                plt.axvline(x=rf1[idx], color='black', alpha=0.3, dashes=[8, 8])

                print('Fit FWHM = {:5.2f} arcmin'.format(2*rf1[idx]))

        if rmax is not None:
            plt.xlim([0, rmax])

        plt.xlabel('Angle [{:s}]'.format(units))
        if log:
            plt.ylabel('Response [dB]')
            if vminl is not None and vmaxl is not None:
                plt.ylim([vminl, vmaxl])
        else:
            plt.ylabel('Linear response')

        if fit_arr1 is not None or fit_arr2 is not None:
            plt.legend()
        if title is not None:
            plt.title(title)

        if fname is not None:
            fname2use = fname
        else:
            fname2use = 'profile'

        plt.savefig(imgd+fname2use+estr+'.png',dpi=dpi)
        plt.close()



    return r, profile
