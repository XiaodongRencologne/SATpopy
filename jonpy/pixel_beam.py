import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set()
import os
import numpy as np
import beam_tools as bt
opj = os.path.join

def hornBeam(data):
    thetas = np.array(data[0][1:]).astype(float)
    freqs  = []
    Eplane = []
    Hplane = []
    for d in data[1:]:
        head = d[0].split("'")
        if 'rEL3Y' in head[0]:
            continue
        freq = float(head[1].strip('GHz'))
        az   = float(head[3].strip('deg'))
        if   az == 0:
            freqs.append(freq)
            Eplane.append(np.array(d[1:]).astype(float))
        elif az == 90:
            Hplane.append(np.array(d[1:]).astype(float))
        else:
            continue

    bm = (np.array(Eplane)**2 + np.array(Hplane)**2)/2.

    #Assume symmetry in theta and mirror
    thetas = np.concatenate((-1*np.flipud(thetas[1:]), thetas))
    bm = np.array([np.concatenate((np.flipud(be[1:]), be)) for be in bm])

    return thetas, freqs, bm/np.amax(bm)

def sinuBeam(data):
    thetas = np.array(data[0][1:]).astype(float)
    freqs  = []
    Eplane = []
    Hplane = []
    for d in data[1:]:
        head = d[0].split("'")
        if 'rEL3Y' in head[0] or '45deg' in head[5]:
            continue
        freq = float(head[3].strip('GHz'))
        az   = float(head[5].strip('deg'))
        if   az == 0:
            freqs.append(freq)
            Eplane.append(np.array(d[1:]).astype(float))
        elif az == 90:
            Hplane.append(np.array(d[1:]).astype(float))
        else:
            continue

    bm = (np.array(Eplane)**2 + np.array(Hplane)**2)/2.

    return thetas, freqs, bm/np.amax(bm)

def pixel_beam(freq, basedir, hornfile):

    hornfile = hornfile.replace('090','MF').replace('150','MF').replace('.cut','.csv')
    data = np.genfromtxt(hornfile, delimiter=',', unpack=True, dtype=np.str)
    if 'Horn' in hornfile:
        thetas, freqs, beams = hornBeam(data)
    else:
        thetas, freqs, beams = sinuBeam(data)

    if freq < np.min(freqs) or freq > np.max(freqs):
        raise ValueError('freq is outside the range defined by the file')

    freqs = np.array(freqs)
    thetas = np.array(thetas)
    freqi = np.argmin(np.abs(freqs-freq))

    return thetas, beams[freqi, :]


def write_hornfile(freq, hornfile, fname, basedir='../20170921_pixel_beams/',
    hprofile=None, hdrstr=None, th0=0, th1=180, phi=[0, 90, 180, 270], pscale=None,
    gaussian_fit=False, debug=True):
    '''
    Function that takes E and H fields from HFSS simulations and write them to a
    file format that GRASP can parse
    '''

    import scipy

    fid = open(fname, 'w')
    if hdrstr is None:
        hdrstr = 'Tabulated feed data\n'

    _, eprofile = pixel_beam(freq, basedir, hornfile)
    ntheta = len(eprofile)
    eprofile = eprofile[int((ntheta-1)/2):]
    ntheta = len(eprofile)

    theta = np.linspace(th0, th1, ntheta)
    # plt.plot(theta, eprofile, label='original')

    if pscale:
        thi = np.linspace(pscale * np.min(theta), pscale * np.max(theta), len(theta))
        ef = scipy.interpolate.interp1d(thi, eprofile, kind='cubic',
            fill_value='extrapolate')
        eprofile = ef(theta)

    if hprofile is None:
        hprofile = np.zeros(ntheta)
        eprofile /= eprofile.max()

    if gaussian_fit:

        popt = bt.gfit1d(theta, eprofile/np.max(eprofile))
        gprofile = bt.gauss_profile(theta, popt[0])

        cr = [-90, -90, 90, 90]
        numel = [1001, 1001, 0]

        bsa = bt.profile_bsa(theta, eprofile/np.max(eprofile), cr, numel)
        bsag = bt.profile_bsa(theta, gprofile, cr, numel)
        bsar = bsa/bsag

        print('bsar: {}'.format(bsar))

        if debug:
            plt.plot(theta, eprofile/np.max(eprofile), label='input')
            plt.plot(theta, gprofile, label='Gaussian fit')
            plt.savefig('gaussian_fit.png')
            plt.close()

    for ph in phi:
        property_str = ' {:1.5e} {:1.5e} {:d} {:.5f}    3    1    2\n'.\
            format(th0, float(th1)/ntheta, ntheta, float(ph))

        fid.write(hdrstr)
        fid.write(property_str)

        if gaussian_fit:
            for gpr, hpr in zip(gprofile, hprofile):
                fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(gpr, hpr))
        else:
            for epr, hpr in zip(eprofile, hprofile):
                fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(epr, hpr))


    fid.close()

    if gaussian_fit:
        return gprofile, hprofile

    else:
        return eprofile, hprofile

def read_hornfile(hornfile):

    out_str = open(hornfile, 'r').read()

    return out_str

def write_multifreq_tabulated_feed(freqs, hornfile, fname,
    basedir='../20170921_pixel_beams/', hdrstr=None, load_csv=True,
    th0=0, th1=180, phi=[0, 90, 180, 270]):
    '''
    Function that takes E and H fields from HFSS simulations and write them to a
    file format that GRASP can parse
    '''

    fid = open(opj('hornfiles/', fname), 'w')

    combined_str = ''

    for i, freq in enumerate(freqs):
        if hdrstr is None:
            hdrstr = 'Tabulated feed data\n'

        if load_csv:

            theta, eprofile = pixel_beam(freq, basedir, hornfile)

            ntheta = len(eprofile)
            eprofile = eprofile[int((ntheta-1)/2):]
            ntheta = len(eprofile)

            hprofile = np.zeros_like(eprofile)
            eprofile /= eprofile.max()

            for ph in phi:
                property_str = ' {:1.5e} {:1.5e} {:d} {:.5f}    3    1    2\n'.\
                    format(th0, float(th1)/ntheta, ntheta, float(ph))

                fid.write(hdrstr)
                fid.write(property_str)
                for epr, hpr in zip(eprofile, hprofile):
                    fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(epr, hpr))

        else:

            if type(hornfile) is not type(list()):
                raise ValueError('hornfile has to be a list in this case')

            out_str = read_hornfile(hornfile[i])
            combined_str += out_str

    if not load_csv:
        fid.write(combined_str)

    fid.close()


if __name__ == '__main__':

    basedir = '/Users/jon/Dropbox/Optics/20170921_pixel_beams/'

    freqs = [84, 88, 92, 96, 100]

    for freq in freqs:

        theta, beam = pixel_beam(freq, basedir, 'Horn_MF_5p3.csv')
        beam /= np.max(beam)

        plt.plot(theta, beam, label='{} GHz'.format(freq))

    plt.legend()
    plt.xlabel('Theta [deg]')
    plt.ylabel('Response')
    plt.savefig('pixel_beam.png')
    plt.close()






