import os, code, sys, pickle, subprocess
import numpy as np

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

def run_command(command):
    #import shlex
    #process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc

def dump_safely(dict2sav, fname):

    attempts = 0
    failed = True
    while attempts < 5:
        try:

            if attempts > 0:
                time.sleep(1)

            pickle.dump(dict2sav, open(fname, 'wb'))
            failed = False
            break
        except IOError:
            attempts += 1
            print('Failed {:d}'.format(attempts))

    if failed:
        raise e

def in_stringlist(str,stringlist):

    for s in stringlist:
        if str == s:
            return True

    return False

def read_lens_file(fname):

    try:
        data = np.loadtxt(fname, skiprows=2)
    except IOError:
        data = np.loadtxt('run/'+fname, skiprows=2)

    r, z = data[:,0], data[:,1]

    return r, z

def write_lens_file(r, z, fname):
    '''
    Takes the r and z vectors defining the radius and
    height of the surface and writes it in the file fname following
    the grasp specifications.
    '''
    header = ' sub.sfc\n   %i           1           0\n'%len(r)

    fid = open(fname, 'w')
    fid.write(header)
    for j, item in enumerate(r):
        fid.write('   {:1.6e} {:1.6e}\n'.format(r[j], z[j]))

    fid.close()

def write_tabulated_feed(fname, eprofile=None, hprofile=None, hdrstr=None,
    th0=0, th1=180, ntheta=181, phi=[0, 90, 180, 270], norm_profiles=True):
    '''
    Function that takes E and H fields from HFSS simulations and write them to a
    file format that GRASP can parse
    '''

    nphi = len(phi)

    fid = open(fname, 'w')
    if hdrstr is None:
        hdrstr = 'Tabulated feed data\n'

    if hprofile is None and eprofile is None:
        raise ValueError('eprofile and hprofile cannot both be None')


    endim = 0 if eprofile is None else eprofile.ndim
    hndim = 0 if hprofile is None else hprofile.ndim

    if (hndim==1 and endim==2) or (hndim==2 and endim==1):
        raise ValueError('Both have to be 2D')

    if hprofile is None:
        hprofile = np.zeros_like(eprofile)
        if endim == 1:
            if norm_profiles:
                eprofile /= eprofile.max()
            eprofile = np.repeat(np.atleast_2d(eprofile), nphi, axis=0)
            hprofile = np.repeat(np.atleast_2d(hprofile), nphi, axis=0)
        else:
            eprofile_norm = np.atleast_2d(np.max(eprofile, axis=1)).T \
                * np.ones_like(eprofile)
            if norm_profiles:
                eprofile = eprofile / eprofile_norm

    elif eprofile is None:
        eprofile = np.zeros_like(hprofile)
        if hndim == 1:
            if norm_profiles:
                hprofile /= hprofile.max()
            hprofile = np.repeat(np.atleast_2d(hprofile), nphi, axis=0)
            eprofile = np.repeat(np.atleast_2d(eprofile), nphi, axis=0)
        else:
            hprofile_norm = np.atleast_2d(np.max(hprofile, axis=1)).T \
                * np.ones_like(hprofile)
            if norm_profiles:
                hprofile = hprofile / hprofile_norm


    else:
        if hndim == 1:
            if norm_profiles:
                hprofile /= hprofile.max()
                eprofile /= eprofile.max()

            hprofile = np.repeat(np.atleast_2d(hprofile), nphi, axis=0)
            eprofile = np.repeat(np.atleast_2d(eprofile), nphi, axis=0)
        else:
            eprofile_norm = np.atleast_2d(np.max(eprofile, axis=1)).T \
                * np.ones_like(eprofile)
            if norm_profiles:
                eprofile = eprofile / eprofile_norm

            hprofile_norm = np.atleast_2d(np.max(hprofile, axis=1)).T \
                * np.ones_like(hprofile)
            if norm_profiles:
                hprofile = hprofile / hprofile_norm

            # hprofile /= np.max(hprofile, axis=1)
            # eprofile /= np.max(eprofile, axis=1)

    for i, ph in enumerate(phi):
        property_str = ' {:1.10e} {:1.10e} {:d} {:.10f}    3    1    2\n'.\
            format(th0, float(th1-th0)/ntheta, ntheta, float(ph))

        fid.write(hdrstr)
        fid.write(property_str)
        # Old version
        # for epr, hpr in zip(eprofile, hprofile):
            # fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(epr, hpr))

        # New version, confusing because the hpr is used as the imagnary
        # component of the E-field
        for epr, hpr in zip(eprofile[i], hprofile[i]):
            fid.write(' {:1.6e} {:1.6e} 0.0 0.0\n'.format(epr, hpr))

    fid.close()

def backup_write_tabulated_feed(fname, eprofile=None, hprofile=None, hdrstr=None,
    th0=0, th1=180, ntheta=181, phi=[0, 90, 180, 270]):
    '''
    Function that takes E and H fields from HFSS simulations and write them to a
    file format that GRASP can parse
    '''

    fid = open(fname, 'w')
    if hdrstr is None:
        hdrstr = 'Tabulated feed data\n'

    if hprofile is None and eprofile is None:
        raise ValueError('eprofile and hprofile cannot both be None')

    if hprofile is None:
        hprofile = np.zeros_like(eprofile)
        eprofile /= eprofile.max()
    elif eprofile is None:
        eprofile = np.zeros_like(hprofile)
        hprofile /= hprofile.max()

    for ph in phi:
        property_str = ' {:1.10e} {:1.10e} {:d} {:.10f}    3    1    2\n'.\
            format(th0, float(th1)/ntheta, ntheta, float(ph))

        fid.write(hdrstr)
        fid.write(property_str)
        # Old version
        # for epr, hpr in zip(eprofile, hprofile):
            # fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(epr, hpr))

        # New version, confusing because the hpr is used as the imagnary
        # component of the E-field
        for epr, hpr in zip(eprofile, hprofile):
            fid.write(' {:1.6e} {:1.6e} 0.0 0.0\n'.format(epr, hpr))

    fid.close()

def write_multifreq_tabulated_feed(fname, eprofile, freqs, hprofile=None, hdrstr=None,
    th0=0, th1=180, ntheta=181, phi=[0, 90, 180, 270]):
    '''
    Function that takes E and H fields from HFSS simulations and write them to a
    file format that GRASP can parse
    '''

    fid = open(fname, 'w')

    for freq in freqs:
        if hdrstr is None:
            hdrstr = 'Tabulated feed data\n'

        if hprofile is None:
            hprofile = np.zeros_like(eprofile)
            eprofile /= eprofile.max()

        for ph in phi:
            property_str = ' {:1.5e} {:1.5e} {:d} {:.5f}    3    1    2\n'.\
                format(th0, float(th1)/ntheta, ntheta, float(ph))

            fid.write(hdrstr)
            fid.write(property_str)
            for epr, hpr in zip(eprofile, hprofile):
                fid.write(' {:1.6e} 0.0 {:1.6e} 0.0\n'.format(epr, hpr))

    fid.close()

def write_text_file(fname, eprofile, hprofile=None, hdrstr=None,
    th0=0, th1=180, ntheta=181, phi=[0, 90, 180, 270]):

    fid = open(fname, 'w')

    header = '# Tabulated pixel beam\n'
    if hdrstr is not None:
        header += '#  {:s}'.format(hdrstr)
    header += '#  column 1: Angle [theta]\n'
    header += '#  column 2: E field, real [V/m]\n'
    header += '#  column 3: E field, imag [V/m]\n'
    fid.write(header)

    if hprofile is None:
        hprofile = np.zeros_like(eprofile)
        eprofile /= eprofile.max()

    theta = np.linspace(th0, th1, ntheta)

    for i, (th, e, h) in enumerate(zip(theta, eprofile, hprofile)):
        fid.write('{:5.5f} {:.8e} {:.8e}\n'.format(th, e, h))

    fid.close()

def write_gxp_file(fname, torfname, tcifname):

    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    fid = open(fname, 'w')

    fid.write('[TOR file]\n')
    fid.write('{:s}\n'.format(torfname))
    fid.write('[TCI file]\n')
    fid.write('{:s}\n'.format(tcifname))
    fid.close()

def write_tor_file(fname, grasp_objects, grasp_object=None, mom_objects=None,
    base_file=None, additional_string=None):
    '''
    Writes a tor file from a collection of objects
    '''

    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    fid = open(fname, 'w')

    # Reading in and writing base tor file template if it exists
    # Normally these are basic components of the model (coordinate systems)
    if base_file is not None:
        with open(base_file, 'r') as file:
            base_str = file.read()

        fid.write(base_str)

    if grasp_objects is not None:
        for go in grasp_objects:
            fid.write(go.grasp_string())

    if mom_objects is not None:
        for mo in mom_objects:
            fid.write(mo.mom_string())

    if grasp_object is not None:
        fid.write(grasp_object.create_frequency_list())

    if additional_string:
        fid.write(additional_string)

    fid.write('\n\n')
    final_str = '//DO NOT MODIFY OBJECTS BELOW THIS LINE.\n' \
        + '//THESE OBJECTS ARE CREATED AND MANAGED BY THE\n' \
        +'//GRAPHICAL USER INTERFACE AND SHOULD NOT BE\n' \
        +'//MODIFIED MANUALLY!\n\n'

    # if grasp_objects is not None:
    #     for go in grasp_objects:
    #         viewf = getattr(go, 'view_string', None)
    #         if callable(viewf):
    #             final_str += go.view_string()

    final_str += '//$$ Saved at 12:34:56 on 29.07.2017 by GRASP ver. 10.5.0 SN=005897\n'

    fid.write(final_str)
    fid.close()

def sequence_str(ref_list, postfix='.po'):

    str_out = ''
    nref_list = len(ref_list)
    for ri, ritem in enumerate(ref_list):
        str_out += 'ref({:s}{:s})'.format(ritem.name, postfix)
        #print(str_out)
        if ri < nref_list-1:
            str_out += ', '
        # else:
        #     str_out += ''

    return str_out

def write_tci_file(det, fname, source_list,
    object_list=None,
    grid=None,
    cut=None,
    r1=None,
    r2=None,
    add_reflectors=False,
    extra_grids=[],
    extra_cuts=[],
    add_cut=False,
    auto_convergence='on',
    other_grids=None,
    ffsources=None,
    additional_sources=None,
    verbose=False):
    '''
    Creates a TCI file

    grid: The canonical output grid for the calculation

    cut: The canonical output cut for the calculation

    r1: (from SO) reflector1 in the system

    r2: (from SO) reflector2 in the system

    add_reflectors: calls a special function that creates commands appropr. for reflectors

    extra_grids: Extra grids to be treated as final output grids

    extra_cuts: Extra cuts to be treated as final output cuts

    auto_convergence (bool): Whether or not to use automatic convergence as part of calculation,

    other_grids: Grids with a well defined source to also calculate

    ffsources: Extra sources to be used when cycling over out_sys.
    These are added together with the last optical element.

    additional_sources: Additional sources that pop up in the chain that aren't necessarily
    supposed to be added with the last optical element.



    '''

    # Formatting source string
    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    out_sys = []
    if grid is None and cut is None:
        raise ValueError('Have to look at results, configure grid or cut')

    if not isinstance(grid, list) and grid is not None:
        grid = [grid]

    if not isinstance(cut, list) and cut is not None:
        cut = [cut]

    if (grid is not None) and not add_reflectors:
        for gridi in grid:
            out_sys.append(gridi)

    if (cut is not None) and not add_reflectors:
        for cuti in cut:
            out_sys.append(cuti)

    if additional_sources is None and object_list is not None:
        additional_sources = len(object_list) * [None]

    if extra_grids:
        out_sys += extra_grids

    if extra_cuts:
        out_sys += extra_cuts

    source_string = ''
    nsources = len(source_list)
    nobjects = len(object_list) if object_list is not None else 0
    nffsources = len(ffsources) if ffsources is not None else 0

    ncmds = 1+nobjects
    cmdn = 1

    for n, source in enumerate(source_list):
        source_string += 'ref({:s})'.format(source.poname)
        if n+1 < nsources:
            source_string += ', '

    string_out = ''
    if nobjects > 0:

        for n, (obj, add_source) in enumerate(zip(object_list, additional_sources)):

            if verbose:
                print('Object name:')
                print(obj.name)

            last_object = True if (n+1)==nobjects else False
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(obj.poname)

            if n == 0:
                # Getting currents on object number one
                string_out += ' source : sequence({:s}),  &\n'.format(source_string)

            else:
                # An object in a sequence of objects

                if add_source is not None:
                    source_list = [object_list[n-1]] + add_source
                    string_out += ' source : sequence({:s}),  &\n'.\
                        format(sequence_str(source_list, postfix='.po'))

                else:
                    string_out += ' source : sequence(ref({:s})),  &\n'.\
                        format(object_list[n-1].poname)

            if np.abs(obj.field_accuracy+80) > 0.001:
                string_out += 'field_accuracy : {:.5f}, '.\
                    format(obj.field_accuracy)

            string_out += 'auto_convergence_of_po : {:s}'.\
                format(obj.auto_convergence)

            if obj.auto_convergence == 'on' and not last_object:
                # Converging on the next object
                string_out += ', convergence_on_scatterer'
                string_out += ' : sequence(ref({:s}.scr))'.\
                    format(object_list[n+1].name)

                if obj.conv_grid is not None:
                    string_out += ', convergence_on_output_grid'

                    try:
                        ncgrids = len(obj.conv_grid)
                    except AttributeError:
                        ncgrids  = 1
                        obj.conv_grid = [obj.conv_grid]

                    if ncgrids == 1:
                        string_out += ' : sequence(ref({:s})))'.\
                            format(obj.conv_grid[0].name)

                    else:
                        tmp_string =''
                        refstr = sequence_str(obj.conv_grid, postfix='')
                        string_out += ' : sequence({:s}))'.format(refstr)

                    if verbose:
                        print('Convergence grid names')
                        for cg in obj.conv_grid:
                            print(cg.name)

                else:
                    string_out += ')'

            elif obj.auto_convergence == 'on':

                # Converging on grid/cut
                string_out += ', convergence_on_output_grid'
                conv_str = ''
                for i, ogrid in enumerate(out_sys):
                    conv_str += 'ref({:s})'.format(ogrid.name)
                    if (i+1) < len(out_sys):
                        conv_str += ','

                string_out += ' : sequence({:s}))'.format(conv_str)

            else:
                string_out += ')'

            # string_out +=  ' cmd_{:d}\n'.format(n+1)
            string_out +=  ' \n'

        last_source = object_list[-1].poname

    else:

        last_source = source_string

    # Looking at field on grid or cut
    conv_str = ''
    for i, ogrid in enumerate(out_sys):

        string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
        if nffsources == 0:
            # print('last_source', last_source)
            if nobjects > 0:
                string_out += ' source : sequence(ref({:s}))) \n'.\
                    format(last_source)
            else:
                string_out += ' source : sequence({:s})) \n'.\
                    format(last_source)

        else:
            refstr = ''
            for i, ffsource in enumerate(ffsources):
                refstr += 'ref({:s})'.format(ffsource.poname)
                refstr += ','

            refstr += 'ref({:s})'.format(last_source)

            string_out += ' source : sequence({:s})) \n'.\
                format(refstr)

    if other_grids is not None:
        for n, ogrid in enumerate(other_grids):
            if ogrid is not None:

                string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
                print(ogrid.sources)
                string_out += ' source : sequence({:s})) \n'.\
                    format(sequence_str(ogrid.sources, postfix='.po'))

    if add_reflectors:
        tmp_string = append2tci(det, r1, r2, source_string, None,
            last_source='ref({:s})'.format(last_source),
            extra_cuts=[cut[0]] if add_cut else None, extra_grids=[grid[0]], add_cut=add_cut)

        string_out += tmp_string


    string_out += 'QUIT\n'

    fid = open(fname,'w+')
    fid.write(string_out)
    fid.close()

def old_write_tci_file(det, fname, source_list,
    object_list=None, grid=None, cut=None, r1=None, r2=None,
    add_reflectors=False, extra_grids=None, extra_cuts=None,
    auto_convergence='on', other_grids=None):
    '''
    Creates a TCI file based on
    '''

    # Formatting source string
    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    out_sys = []
    if grid is None and cut is None:
        raise ValueError('Have to look at results, configure grid or cut')

    if grid is not None:
        out_sys.append(grid)
    if cut is not None:
        out_sys.append(cut)

    # if extra_grids is not None:
    #     out_sys.append(extra_grid)

    # if extra_cuts is not None:
    #     out_sys.append(extra_cut)

    source_string = ''
    nsources = len(source_list)
    nobjects = len(object_list) if object_list is not None else 0

    ncmds = 1+nobjects

    cmdn = 1

    for n, source in enumerate(source_list):
        source_string += 'ref({:s})'.format(source.poname)
        if n+1 < nsources:
            source_string += ', '

    string_out = ''
    if nobjects > 0:

        for n, obj in enumerate(object_list):

            last_object = True if (n+1)==nobjects else False
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(obj.poname)

            if n == 0:
                # Getting currents on object number one
                string_out += ' source : sequence({:s}),  &\n'.format(source_string)

            else:
                # An object in a sequence of objects
                string_out += ' source : sequence(ref({:s})),  &\n'.\
                    format(object_list[n-1].poname)

            if np.abs(obj.field_accuracy+80) > 0.001:
                string_out += 'field_accuracy : {:.5f}, '.\
                    format(obj.field_accuracy)

            string_out += 'auto_convergence_of_po : {:s}'.\
                format(obj.auto_convergence)

            if obj.auto_convergence == 'on' and not last_object:
                # Converging on the next object
                string_out += ', convergence_on_scatterer'
                string_out += ' : sequence(ref({:s}.scr))'.\
                    format(object_list[n+1].name)

                if obj.conv_grid is not None:

                    if len(obj.conv_grid) == 1:
                        string_out += ', convergence_on_output_grid'
                        string_out += ' : sequence(ref({:s})))'.\
                            format(obj.conv_grid.name)

                    else:
                        tmp_string =''
                        for ocg in obj.conv_grid:
                            tmp_string += 'ref({:s})'
                else:
                    string_out += ')'


            elif obj.auto_convergence == 'on':
                # Converging on grid/cur
                string_out += ', convergence_on_output_grid'

                conv_str = ''
                for i, grid in enumerate(out_sys):
                    conv_str += 'ref({:s})'.format(grid.name)
                    if (i+1) < len(out_sys):
                        conv_str += ','

                string_out += ' : sequence({:s}))'.format(conv_str)
            else:
                string_out += ')'

            string_out +=  ' cmd_{:d}\n'.format(n+1)

        last_source = object_list[-1].poname

    else:

        last_source = source_string

    # Looking at field on grid or cut
    conv_str = ''
    for i, grid in enumerate(out_sys):

        string_out += 'COMMAND OBJECT {:s} get_field ('.format(grid.name)
        string_out += ' source : sequence(ref({:s}))) cmd_{:d}\n'.\
            format(last_source, ncmds+i)


    if other_grids is not None:
        for n, ogrid in enumerate(other_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
            string_out += ' source : sequence(ref({:s}))) cmd_{:d}\n'.\
                format(ogrid.source_name, ncmds+n+1)


    if add_reflectors:
        tmp_string = append2tci(det, r1, r2, source_string, None,
            last_source=last_source,
            extra_cuts=extra_cuts, extra_grids=extra_grids)
        string_out += tmp_string


    string_out += 'QUIT\n'

    fid = open(fname,'w+')
    fid.write(string_out)
    fid.close()

def append2tci(det, r1, r2, source_string, cluster=None, last_source=None, extra_cuts=None,
    extra_grids=None, fa=-40, add_cut=False):

    string_out = ''
    if cluster is not None:

        string_out += 'COMMAND OBJECT M2_mod.po get_currents ('
        string_out += ' source : sequence({:s}, ref({:s}.bor_mom)), & \n'.\
            format(source_string, cluster.name)

        string_out += 'field_accuracy : {:.5f}, '.format(r2.field_accuracy)
        string_out += ' convergence_on_scatterer : sequence(ref({:s}.scr)))\n'.\
            format(r1.name)

        string_out += 'COMMAND OBJECT M1_mod.po get_currents ('
        string_out += ' source : sequence(ref(M2_mod.po)), '
        string_out += 'field_accuracy : {:.5f}, & \n'.format(r1.field_accuracy)
        string_out += 'convergence_on_output_grid : sequence('
        if add_cut:
            string_out += 'ref({:s}_reflector_sc), ref({:s}_reflector_sg)))\n'.\
                format(det, det)
        else:
            string_out += 'ref({:s}_reflector_sg)))\n'.format(det)

    else:

        string_out += 'COMMAND OBJECT M2_mod.po get_currents ('
        string_out += ' auto_convergence_of_po : on,'
        string_out += ' source : sequence({:s}), '.\
            format(last_source)

        string_out += 'field_accuracy : {:.5f}, '.format(r2.field_accuracy)
        string_out += ' convergence_on_scatterer : sequence(ref({:s}.scr)))\n'.\
            format(r1.name)

        string_out += 'COMMAND OBJECT M1_mod.po get_currents ('
        string_out += ' auto_convergence_of_po : on,'

        string_out += ' source : sequence(ref(M2_mod.po)), '
        string_out += 'field_accuracy : {:.5f}, & \n'.format(r1.field_accuracy)
        string_out += ' convergence_on_output_grid : sequence('
        if add_cut:
            string_out += 'ref({:s}_reflector_sc), ref({:s}_reflector_sg)))\n'.format(det, det)
        else:
            string_out += 'ref({:s}_reflector_sg)))\n'.format(det)

    # string_out += 'COMMAND OBJECT {:s}_reflector_sg get_field ( source : '.format(det)
    # string_out += 'sequence(ref(M1_mod.po)))\n'

    # string_out += 'COMMAND OBJECT {:s}_reflector_sc get_field ( source : '.format(det)
    # string_out += 'sequence(ref(M1_mod.po)))\n'


    out_sys = []
    if extra_grids is not None:
        for extra_grid in extra_grids:
            out_sys.append(extra_grid)

    if extra_cuts is not None:
        for extra_cut in extra_cuts:
            out_sys.append(extra_cut)

    for out in out_sys:
        string_out += 'COMMAND OBJECT {:s} get_field ('.format(out.name)
        string_out += ' source : sequence(ref(M1_mod.po))) \n'

    return string_out

    # COMMAND OBJECT M2_mod.po get_currents ( source :  &
    # sequence(ref(cluster10.bor_mom), ref(horn.po)), field_accuracy : -40.0,  &
    # convergence_on_scatterer : sequence(ref(M1_mod.scr)))
    # COMMAND OBJECT M1_mod.po get_currents ( source : sequence(ref(M2_mod.po)),  &
    # field_accuracy : -40.0, convergence_on_output_grid :  &
    # sequence(ref(test_mom_h1_t1r000phi000_10_sc2),  &
    # ref(test_mom_h1_t1r000phi000_10_sg2)))
    # COMMAND OBJECT test_mom_h1_t1r000phi000_10_sc2 get_field ( source :  &
    # sequence(ref(M1_mod.po)))
    # COMMAND OBJECT test_mom_h1_t1r000phi000_10_sg2 get_field ( source :  &
    # sequence(ref(M1_mod.po))) cmd_9

def write_mom_tci_file(fname, source_list, cluster,
    grid=None, cut=None, extra_grids=None, extra_cuts=None, other_grids=None,
    telescope_grids=None, auto_convergence=True,
    add_sources=True,
    fa=-40):

    out_sys = []
    if grid is None and cut is None and extra_cuts is None and other_grids is None:
        raise ValueError('Have to look at results, configure grid or cut')

    if grid is not None:
        out_sys.extend(grid if type(grid) == list else [grid])

    if cut is not None:
        out_sys.extend(cut if type(cut) == list else [cut])

    if extra_grids:
        out_sys += extra_grids

    string_out = ''
    nsources = len(source_list)
    source_string = ''
    for n, source in enumerate(source_list):
        source_string += 'ref({:s})'.format(source.poname)
        if n+1 < nsources:
            source_string += ', '

    string_out += 'COMMAND OBJECT {:s} get_currents ('.format(cluster.moname)
    string_out += ' source : sequence({:s})'.format(source_string)
    # string_out += ', convergence_on_output_grid'

    # conv_str = ''
    # for i, ogrid in enumerate(out_sys):
    #     conv_str += 'ref({:s})'.format(ogrid.name)
    #     if (i+1) < len(out_sys):
    #         conv_str += ', '

    string_out += ')\n'

    for out in out_sys:
        print(out.name)
        string_out += 'COMMAND OBJECT {:s} get_field ('.format(out.name)
        if add_sources:
            string_out += ' source : sequence({:s}, ref({:s}.mom))) \n'.\
                    format(source_string, cluster.name)
        else:
            string_out += ' source : sequence(ref({:s}.mom))) \n'.format(cluster.name)

    if other_grids is not None:
        for n, ogrid in enumerate(other_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
            if add_sources:
                string_out += ' source : sequence({:s}, ref({:s}.mom))) \n'.\
                    format(source_string, cluster.name)
            else:
                string_out += ' source : sequence(ref({:s}.mom))) \n'.format(cluster.name)

    if telescope_grids is not None:
        for n, tgrid in enumerate(telescope_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(tgrid.name)
            string_out += ' source : sequence({:s})) \n'.\
                    format(source_string)

    fid = open(fname,'w+')
    fid.write(string_out)
    fid.close()

def write_bor_mom_tci_file(det, fname, source_list, cluster,
    r1=None, r2=None, add_reflectors=False,
    add_absorbing_potube=False, add_absorbing_stop=False,
    potube=None, postop=None,
    grid=None, cut=None, extra_grids=None, extra_cuts=None, other_grids=None,
    auto_convergence=True, fa=-40, square_pixels=False,
    telescope_grids=None):

    # Formatting source string
    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    out_sys = []
    if grid is None and cut is None:
        raise ValueError('Have to look at results, configure grid or cut')

    if grid is not None:
        out_sys.extend(grid if type(grid) == list else [grid])

    if cut is not None:
        out_sys.extend(cut if type(cut) == list else [cut])
        #out_sys.append(cut)

    # if extra_grids is not None:
    #     for extra_grid in extra_grids:
    #         out_sys.append(extra_grid)

    if extra_grids:
        out_sys += extra_grids

    # if extra_cuts is not None:
    #     for extra_cut in extra_cuts:
    #         out_sys.append(extra_cut)

    string_out = ''

    nsources = len(source_list)
    if not square_pixels:

        source_string = ''
        for n, source in enumerate(source_list):
            source_string += 'ref({:s})'.format(source.poname)
            if n+1 < nsources:
                source_string += ', '

        string_out += 'COMMAND OBJECT {:s} get_currents ('.format(cluster.moname)
        string_out += ' source : sequence({:s})'.format(source_string)


    else:
        source_string += 'ref({:s})'.format(source_list[0].poname)

        string_out += 'COMMAND OBJECT {:s} get_currents ('.format(source_list[1].poname)
        string_out += ' source : sequence({:s})'.format(source_string)
        string_out += ', auto_convergence_of_po : on'
        if np.abs(fa+80) > 0.001:
            string_out += ', field_accuracy : {:.5f}'.format(fa)


    conv_str = ''

    if auto_convergence:
        # Converging on grid/cut

        string_out += ', auto_convergence_of_po : on'
        if np.abs(fa+80) > 0.001:
            string_out += ', field_accuracy : {:.5f}'.format(fa)

        if not add_absorbing_stop:
            string_out += ', convergence_on_output_grid'

        for i, ogrid in enumerate(out_sys):
            conv_str += 'ref({:s})'.format(ogrid.name)
            if (i+1) < len(out_sys):
                conv_str += ', '

        if not add_absorbing_stop:
            string_out += ' : sequence({:s})'.format(conv_str)

    if add_absorbing_potube:
        tube_str = ''
        tube_str_po = ''
        for i, pot in enumerate(potube):
            tube_str += 'ref({:s}.scr)'.format(pot.name)
            tube_str_po += 'ref({:s})'.format(pot.poname)
            if (i+1) < len(potube):
                tube_str += ', '
                tube_str_po += ', '

        string_out += ', convergence_on_scatterer'
        string_out += ' : sequence({:s})'.format(tube_str)

    elif add_absorbing_stop:

        stop_str = ''
        stop_str_po = ''
        for i, pot in enumerate(postop):
            stop_str += 'ref({:s}.scr)'.format(pot.name)
            stop_str_po += 'ref({:s})'.format(pot.poname)
            # if (i+1) < len(potube):
            #     stop_str += ', '
            #     stop_str_po += ', '

        string_out += ', convergence_on_scatterer'
        string_out += ' : sequence({:s})'.format(stop_str)


    string_out += ')\n'

    # Adding an absorbing physical optics tube around the system
    if add_absorbing_potube:
        source_string += ', ref({:s})'.format(cluster.moname)
        for i, pot in enumerate(potube):
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(pot.poname)
            string_out += ' source : sequence({:s})'.format(source_string)
            string_out += ', auto_convergence_of_po : on'
            if auto_convergence:
                string_out += ', convergence_on_output_grid'
                string_out += ' : sequence({:s})'.format(conv_str)

                if np.abs(fa+80) > 0.001:
                    string_out += ', field_accuracy : {:.5f}'.format(fa)

            string_out += ')\n'

    elif add_absorbing_stop:
        source_string += ', ref({:s})'.format(cluster.moname)
        for i, pot in enumerate(postop):
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(pot.poname)
            string_out += ' source : sequence({:s})'.format(source_string)
            string_out += ', auto_convergence_of_po : on'

            print('auto_convergence:', auto_convergence)
            print('conv_str:', conv_str)
            if auto_convergence:
                string_out += ', convergence_on_output_grid'
                string_out += ' : sequence({:s})'.format(conv_str)

                if np.abs(fa+80) > 0.001:
                    string_out += ', field_accuracy : {:.5f}'.format(fa)

            string_out += ')\n'

    ### Going through grids and cuts
    for out in out_sys:
        print(out.name)
        string_out += 'COMMAND OBJECT {:s} get_field ('.format(out.name)
        if add_absorbing_potube:
            string_out += ' source : sequence({:s}, {:s})) \n'.\
                format(source_string, tube_str_po)
        elif add_absorbing_stop:
            string_out += ' source : sequence({:s})) \n'.\
                format(stop_str_po)
        else:
            string_out += ' source : sequence({:s}, ref({:s}.bor_mom))) \n'.\
                format(source_string, cluster.name)

    if other_grids is not None:
        for n, ogrid in enumerate(other_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
            if add_absorbing_potube:
                string_out += ' source : sequence({:s}, {:s})) \n'.\
                    format(source_string, tube_str_po)
            elif add_absorbing_stop:
                string_out += ' source : sequence({:s})) \n'.\
                    format(stop_str_po)
            else:
                string_out += ' source : sequence({:s}, ref({:s}.bor_mom))) \n'.\
                    format(source_string, cluster.name)

    if telescope_grids is not None:
        for n, tgrid in enumerate(telescope_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(tgrid.name)
            string_out += ' source : sequence({:s}, ref({:s}.bor_mom))) \n'.\
                    format(source_string, cluster.name)

    if add_reflectors:

        if add_absorbing_potube:
            tmp_string = append2tci(det, r1, r2, [source_string, tube_str_po])
        elif add_absorbing_stop:
            tmp_string = append2tci(det, r1, r2, stop_str_po, last_source=stop_str_po,
                extra_cuts=[cut] if cut is not None else None, extra_grids=[grid])
        else:
            tmp_string = append2tci(det, r1, r2, [source_string, cluster.name])

        string_out += tmp_string

    string_out += 'QUIT\n'

    fid = open(fname,'w+')
    fid.write(string_out)
    fid.close()

def write_hybrid_mom_tci_file(det, fname, source_list, cluster,
    r1=None, r2=None, add_reflectors=False,
    add_absorbing_potube=False, add_absorbing_stop=False,
    potube=None, postop=None,
    add_reflecting_baffle=False, add_baffle_aperture=False,
    baffle=None,
    grid=None, cut=None, extra_grids=None, extra_cuts=None, other_grids=None,
    auto_convergence=False, fa=-40, square_pixels=False,
    additional_source=None):

    # Formatting source string
    if len(fname.split('/')) > 1:
        dstr = os.path.join(*fname.split('/')[:-1])
        if not os.path.exists(dstr):
            os.makedirs(dstr)

    out_sys = []
    if grid is None and cut is None:
        raise ValueError('Have to look at results, configure grid or cut')

    if grid is not None:
        out_sys.extend(grid if type(grid) == list else [grid])

    if cut is not None:
        out_sys.extend(cut if type(cut) == list else [cut])
        #out_sys.append(cut)

    # if extra_grids is not None:
    #     for extra_grid in extra_grids:
    #         out_sys.append(extra_grid)

    if extra_grids:
        out_sys += extra_grids


    # if extra_cuts is not None:
    #     for extra_cut in extra_cuts:
    #         out_sys.append(extra_cut)

    string_out = ''

    nsources = len(source_list)
    if not square_pixels:
        source_string = ''
        for n, source in enumerate(source_list):
            source_string += 'ref({:s})'.format(source.poname)
            if n+1 < nsources:
                source_string += ', '

        string_out += 'COMMAND OBJECT {:s} get_currents ('.format(cluster.moname)
        string_out += ' source : sequence({:s})'.format(source_string)

    else:
        source_string += 'ref({:s})'.format(source_list[0].poname)

        string_out += 'COMMAND OBJECT {:s} get_currents ('.format(source_list[1].poname)
        string_out += ' source : sequence({:s})'.format(source_string)
        string_out += ', auto_convergence_of_po : on'
        if np.abs(fa+80) > 0.001:
            string_out += ', field_accuracy : {:.5f}'.format(fa)


    if auto_convergence:
        # Converging on grid/cut

        string_out += ', auto_convergence_of_po : on'
        if np.abs(fa+80) > 0.001:
            string_out += ', field_accuracy : {:.5f}'.format(fa)

        string_out += ', convergence_on_output_grid'
        conv_str = ''
        for i, grid in enumerate(out_sys):
            conv_str += 'ref({:s})'.format(grid.name)
            if (i+1) < len(out_sys):
                conv_str += ', '

        string_out += ' : sequence({:s})'.format(conv_str)

    if add_absorbing_potube:
        tube_str = ''
        tube_str_po = ''
        for i, pot in enumerate(potube):
            tube_str += 'ref({:s}.scr)'.format(pot.name)
            tube_str_po += 'ref({:s})'.format(pot.poname)
            if (i+1) < len(potube):
                tube_str += ', '
                tube_str_po += ', '

        string_out += ', convergence_on_scatterer'
        string_out += ' : sequence({:s})'.format(tube_str)
        print('tube_str', tube_str)


    print('add_absorbing_stop:', add_absorbing_stop)
    if add_absorbing_stop:

        stop_str = ''
        stop_str_po = ''
        for i, pot in enumerate(postop):
            stop_str += 'ref({:s}.scr)'.format(pot.name)
            stop_str_po += 'ref({:s})'.format(pot.poname)
            if (i+1) < len(postop):
                stop_str += ', '
                stop_str_po += ', '

        string_out += ', convergence_on_scatterer'
        string_out += ' : sequence({:s})'.format(stop_str)
        print('stop_str', stop_str)

    string_out += ')\n'

    if add_absorbing_stop and add_reflecting_baffle:
        baffle_str = 'ref({:s}.scr)'.format(baffle[0].name)

    if add_absorbing_stop and add_baffle_aperture:
        baffle_str = 'ref({:s}.scr)'.format(baffle[0].name)

    # Adding an absorbing physical optics tube around the system
    if add_absorbing_potube and not add_baffle_aperture:
        print('we must get to here')
        source_string += ', ref({:s})'.format(cluster.moname)

        if additional_source:
            source_string += ', ref({:s})'.format(additional_source[0].poname)

        print('potube:', potube)
        for i, pot in enumerate(potube):
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(pot.poname)
            string_out += ' source : sequence({:s})'.format(source_string)
            string_out += ', auto_convergence_of_po : on'
            if auto_convergence:
                string_out += ', convergence_on_output_grid'
                string_out += ' : sequence({:s})'.format(conv_str)

                if np.abs(fa+80) > 0.001:
                    string_out += ', field_accuracy : {:.5f}'.format(fa)

            string_out += ')\n'

    elif add_baffle_aperture:# or add_reflecting_baffle:
        print('are we here')
        source_string += ', ref({:s})'.format(cluster.moname)

        if additional_source:
            source_string += ', ref({:s})'.format(additional_source[0].poname)

        # string_out += 'COMMAND OBJECT {:s} get_currents ('.format(postop[0].poname)
        # string_out += ' source : sequence({:s})'.format(source_string)
        # string_out += ', auto_convergence_of_po : on'
        # if auto_convergence:
        #     if not add_baffle_aperture:
        #         string_out += ', convergence_on_output_grid'
        #         string_out += ' : sequence({:s})'.format(conv_str)

        #     string_out += ', convergence_on_scatterer'
        #     string_out += ' : sequence({:s})'.format(baffle_str)

        #     if np.abs(fa+80) > 0.001:
        #         string_out += ', field_accuracy : {:.5f}'.format(fa)

        #     string_out += ')\n'

        string_out += 'COMMAND OBJECT {:s} get_currents ('.format(baffle[0].poname)
        string_out += ' source : sequence({:s})'.format(source_string)
        string_out += ', auto_convergence_of_po : on'
        if auto_convergence:
            string_out += ', convergence_on_output_grid'
            string_out += ' : sequence({:s})'.format(conv_str)
            if np.abs(fa+80) > 0.001:
                string_out += ', field_accuracy : {:.5f}'.format(fa)

            string_out += ')\n'

        # All fields have to go through the aperture
        if add_baffle_aperture:
            source_string = 'ref({:s})'.format(baffle[0].poname)
        else:
            source_string += ', ref({:s})'.format(baffle[0].poname)

    elif add_absorbing_stop:
        print('how about here')
        source_string += ', ref({:s})'.format(cluster.moname)
        for i, pot in enumerate(postop):
            string_out += 'COMMAND OBJECT {:s} get_currents ('.format(pot.poname)
            string_out += ' source : sequence({:s})'.format(source_string)
            string_out += ', auto_convergence_of_po : on'
            if auto_convergence:
                string_out += ', convergence_on_output_grid'
                string_out += ' : sequence({:s})'.format(conv_str)

                if np.abs(fa+80) > 0.001:
                    string_out += ', field_accuracy : {:.5f}'.format(fa)

            string_out += ')\n'

    print('out_sys:', out_sys)

    ### Going through grids and cuts
    for out in out_sys:

        print('In out_sys')

        print(add_absorbing_potube)
        print(add_baffle_aperture)
        print(add_absorbing_stop)

        string_out += 'COMMAND OBJECT {:s} get_field ('.format(out.name)
        if add_absorbing_potube:
            string_out += ' source : sequence({:s}, {:s})) \n'.\
                format(source_string, tube_str_po)
        elif add_baffle_aperture:
            string_out += ' source : sequence({:s})) \n'.\
                format(source_string)

            print('Here in out_sys')

        # elif add_absorbing_stop and add_reflecting_baffle:
        #     if additional_source is not None:
        #         string_out += ' source : sequence(ref({:s}), ref({:s}), ref({:s}))) \n'.\
        #             format(postop[0].poname, baffle[0].poname, additional_source[0].poname)
        #     else:

        #         string_out += ' source : sequence(ref({:s}), ref({:s}))) \n'.\
        #             format(postop[0].poname, baffle[0].poname)

        elif add_absorbing_stop:
            if additional_source is not None:
                print(additional_source)
                string_out += ' source : sequence(ref({:s}), ref({:s}))) \n'.\
                    format(postop[0].poname, additional_source[0].poname)
            else:
                string_out += ' source : sequence(ref({:s}))) \n'.\
                    format(postop[0].poname)
        else:
            string_out += ' source : sequence({:s}, ref({:s}.bor_mom))) \n'.\
                format(source_string, cluster.name)

    if other_grids is not None:
        for n, ogrid in enumerate(other_grids):
            string_out += 'COMMAND OBJECT {:s} get_field ('.format(ogrid.name)
            if add_absorbing_potube:
                string_out += ' source : sequence({:s}, {:s})) \n'.\
                    format(source_string, tube_str_po)

            elif add_absorbing_stop and add_reflecting_baffle:
                string_out += ' source : sequence(ref({:s}), ref({:s}))) \n'.\
                    format(postop[0].poname, baffle[0].poname)

            elif add_absorbing_stop and add_baffle_aperture:
                string_out += ' source : sequence(ref({:s}))) \n'.\
                    format(baffle[0].poname)

            elif add_absorbing_stop:
                string_out += ' source : sequence(ref({:s}))) \n'.\
                    format(postop[0].poname)
            else:
                string_out += ' source : sequence({:s}, ref({:s}.bor_mom))) \n'.\
                    format(source_string, cluster.name)

    if add_reflectors:
        tmp_string = append2tci(det, r1, r2, source_string, cluster,
            extra_cuts=extra_cuts, extra_grids=extra_grids)
        string_out += tmp_string

    string_out += 'QUIT\n'

    fid = open(fname,'w+')
    fid.write(string_out)
    fid.close()

def create_frequency_list(self):

    str_out = '\n\n'
    str_out += '{:s}  frequency\n'.format(self.frequency_list)
    str_out += '(\n'

    if len(self.frequency) == 1:
        str_out += '  frequency_list   : sequence({:3.1f} GHz)\n'.\
            format(self.frequency[0])
        str_out += ')'

    else:

        freqstr = ''.join(['{:.5f} GHz, '.format(nu) for nu in self.frequency[:-1]])
        freqstr += '{:.5f} GHz'.format(self.frequency[-1])
        str_out += '  frequency_list   : sequence({:s})\n'.format(freqstr)
        str_out += ')'

    return str_out

class CoordinateSystem():

    def __init__(self, x=0, y=0, z=0, units='mm', name=None, base=None):

        self.name = name if name is not None else 'undefined'
        self.units = units
        self.x = x
        self.y = y
        self.z = z
        self.base = base

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}.csy coor_sys\n'.format(self.name)
        gstring += '(\n'
        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
                format(self.z, self.units)
        if self.base is not None:
            gstring += ', \n'
            gstring += '  base             : ref({:s}),\n'.format(self.base)
        else:
            gstring += '\n'
        gstring += ')\n'

        return gstring

class GraspCoordinateSystem():

    def __init__(self, theta=0, phi=0, psi=0,
        x=0, y=0, z=0, units='mm', name=None, base=None):

        self.name = name if name is not None else 'undefined'
        self.units = units
        self.theta = theta
        self.phi = phi
        self.psi = psi
        self.x = x
        self.y = y
        self.z = z
        self.base = base

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}.csy coor_sys_grasp_angles\n'.format(self.name)
        gstring += '(\n'
        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
                format(self.z, self.units)
        gstring += '  theta           : {:5.5f}, \n'.format(self.theta)
        gstring += '  phi           : {:5.5f}, \n'.format(self.phi)
        gstring += '  psi           : {:5.5f}'.format(self.psi)

        if self.base is not None:
            gstring += ', \n'
            gstring += '  base             : ref({:s}),\n'.format(self.base)
        else:
            gstring += '\n'
        gstring += ')\n'

        return gstring

class ScattererCluster():

    def __init__(self, scatterer_list, name=None, frequency_list=None,
        min_ppm=1.0e-4, use_bor=True, max_mesh_length=1.5,
        allow_mlfmm=True):

        self.scatterer_list = scatterer_list
        self.name = name if name is not None else 'undefined'
        self.frequency_list = frequency_list
        self.min_ppm = min_ppm
        self.use_bor = use_bor
        self.allow_mlfmm = allow_mlfmm
        self.max_mesh_length = max_mesh_length

        self.moname = name+'.bor_mom' if use_bor else name+'.mom'

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s} scatterer_cluster\n'.format(self.name)
        gstring += '(\n'

        seq_str = ''
        # for scat in self.scatterer_list:
        #     seq_str += 'ref({:s}.msh),'.format(scat.name)

        for scat in self.scatterer_list:
            seq_str += 'ref({:s}),'.format(scat.name + '.msh')
            # seq_str += 'ref({:s}),'.format(scat.name)

        seq_str = seq_str[:-1] # Removing the last comma

        gstring += '  scatterers           : sequence({:s})\n'.\
                format(seq_str)
        gstring += ')\n'

        if self.use_bor:

            gstring += '{:s}.bor_mom   bor_mom'.format(self.name)
            gstring += '(\n'
            gstring += '  frequency      : ref({:s}),\n'.format(self.frequency_list)
            gstring += '  scatterer      : ref({:s}),\n'.format(self.name)
            gstring += '  min_power_percentage_per_mode : {:.2e}\n'.format(self.min_ppm)
            gstring += ')\n'

        else:

            gstring += '{:s}.mom   mom'.format(self.name)
            gstring += '(\n'
            gstring += '  frequency      : ref({:s}),\n'.format(self.frequency_list)
            gstring += '  max_mesh_length: {:.3f},\n'.format(self.max_mesh_length)
            if not self.allow_mlfmm:
                gstring += '  iterative_solution: struct(use_mlfmm: off),\n'
            gstring += '  scatterer      : ref({:s})\n'.format(self.name)
            gstring += ')\n'


        return gstring

class CircularPlate():

    def __init__(self, name=None, fancy_name=None, base=None, units='mm', radius=1.0,
        x=0., y=0., z=0.):
        '''

        '''

        self.name = name
        self.fancy_name = fancy_name
        self.units = units
        self.radius = radius
        self.base = base
        self.x = x
        self.y = y
        self.z = z

    def grasp_string(self):
        '''

        '''

        gstring = '{:s}  circular_plate\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '   radius          : {:.5f} {:s}\n'.format(self.radius, self.units)
        gstring += ')\n'

        gstring += '\n{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base         : ref({:s}.csy),\n'.\
                format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        gstring += ')\n\n'

        return gstring

class RectangularPlate():

    def __init__(self, name=None, fancy_name=None, base=None, units='mm',
        xw=0.1, yw=0.1, x=0., y=0., z=0.):
        '''

        '''

        self.name = name
        self.fancy_name = fancy_name
        self.units = units
        self.base = base
        self.xw = xw
        self.yw = yw

        self.x = x
        self.y = y
        self.z = z

    def grasp_string(self):
        '''

        '''

        gstring = '{:s}  rectangular_plate\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '   corner_1          : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(-self.xw, self.units, -self.yw, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(0.0, self.units)
        gstring += '   corner_2          : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(self.xw, self.units, -self.yw, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(0.0, self.units)
        gstring += '   opp_point          : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(self.xw, self.units, self.yw, self.units)+'z: {:5.5f} {:s})\n'.\
            format(0.0, self.units)

        gstring += ')\n'

        gstring += '\n{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base         : ref({:s}.csy),\n'.\
                format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        gstring += ')\n\n'

        return gstring

class Baffle():
    '''
    '''

    def __init__(self, name=None, fancy_name=None,
        r1=200, r2=None, angle=10, height=100, t=1,
        units='mm', base=None, invert_x=False, invert_z=False, x=0, y=0, z=0,
        xv=0., yv=0., zv=0., xa=0., ya=0., za=1.0,
        momcsy=None,
        ztrans=0., xtrans=None, ytrans=None,
        n=1.0, PEC=True, closed=False,
        frequency_list=None, field_accuracy=-80, auto_convergence='on',
        conv_grid=None, po1=50, po2=50, ptd=50, spillover=False,
        permeability=1.0, loss_tangent=0.0, flip=False, absorbing=False,
        Zs_real=0., Zs_imag=0., tol=1e-6, Np=51,
        spline_reflector=False, circular_strut=False,
        elevation=None, add_absorbing_back=False, add_pecabs=True):
        '''
        Baffle is implemented as a PO reflector OR as a MoM object
        (either conducting or absorbing)

        angle = 0 is a cylinder
        angle 0 < angle < 90 is a baffle
        angle = 90 is an annulus
        '''

        if spline_reflector and circular_strut:
            raise ValueError('Cannot be both spline_reflector and circular_strut')

        self.name = name

        self.r1 = r1
        self.angle = angle

        if angle == 0.:
            self.angle_fudge = 0.1

        self.height = height
        self.t = t

        # PO-related
        self.poname = name+'.po'
        self.fancy_name = fancy_name
        self.frequency_list = frequency_list
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence

        self.po1 = po1
        self.po2 = po2
        self.ptd = ptd
        self.spillover = spillover

        if conv_grid is not None and not isinstance(conv_grid, list):
            self.conv_grid = [conv_grid]
        else:
            self.conv_grid = conv_grid


        if r2 is None:
            self.r2 = r1 + height*np.tan(np.radians(angle))
        else:
            self.r2 = r2

        self.units = units

        self.PEC = PEC
        self.closed = closed

        self.ztrans = ztrans
        self.xtrans = xtrans
        self.ytrans = ytrans

        self.momcsy = momcsy

        self.tol = tol

        self.base = base
        self.invert_x = invert_x
        self.invert_z = invert_z
        self.flip = flip

        self.x, self.y, self.z = x, y, z

        # Related to the circular cone scatterer
        self.spline_reflector = spline_reflector
        self.xa, self.ya, self.za = xa, ya, za
        self.xv, self.yv, self.zv = xv, yv, zv
        self.circular_strut = circular_strut
        self.permittivity = n**2
        self.permeability = permeability
        self.loss_tangent = loss_tangent
        self.Zs_real = Zs_real
        self.Zs_imag = Zs_imag

        if not self.PEC and Zs_real >= 0.0:
            raise ValueError('Are you sure you want to have a dielectric with Zs_real >= 0? See GRASP Manual')

        self.absorbing = absorbing

        self.elevation = elevation

        self.add_absorbing_back = add_absorbing_back

        self.add_pecabs = add_pecabs

        self.Np = Np


    def mom_mesh(self):

        if self.r1 == self.r2:
            r = self.r1*np.ones(self.Np)
        else:
            r = np.linspace(self.r1, self.r2, self.Np)

        rinv = r[::-1]

        if self.angle == 90: # Stop
            z1 = np.zeros_like(r)

            if self.closed:
                z2 = z1 + self.t
                r = np.hstack((r, rinv))
                z1 = np.hstack((z1, z2))

        elif self.angle == 0 or self.angle == 0.1: # Cylinder
            # z1 = np.linspace(0, self.height, self.Np)
            z1 = np.linspace(self.z, self.z + self.height, self.Np)
            if self.closed:

                z1inv = z1[::-1]
                r = np.hstack((r, rinv+self.t))
                z1 = np.hstack((z1, z1inv))

        else:
            #print('Baffle')
            z1 = (r-self.r1)/np.tan(np.radians(self.angle))

            if self.closed:
                # z2 = np.copy(z1)
                z1inv = z1[::-1]
                r = np.hstack((r, rinv - self.t))
                z1 = np.hstack((z1, z1inv))

            #print(z1)
            #print(r)

        if self.flip:
            z1 = -z1

        # print(self.name)
        # print(self.height)
        # print(z1)
        # print(r)
        # print(self.r1)
        # print(self.angle)

        self.r = r
        self.zag = self.ztrans + z1

        # self.r = np.hstack((r, rinv))
        # self.zag = self.ztrans + np.hstack((z1, z1[::-1] - 0.1))

    def mom_string(self):

        mstring = '{:s}.msh  bor_mesh\n'.format(self.name)
        mstring += '(\n'
        if self.xtrans is not None and self.ytrans is not None:
            mstring += '  coor_sys          : ref({:s}_trans.csy),\n'.format(self.momcsy)
        else:
            mstring += '  coor_sys          : ref({:s}.csy),\n'.format(self.momcsy)
        mstring += '  regions           : table\n'
        mstring += '  (\n'

        if not self.PEC:
            mstring += '    1    {:5.5f} {:5.5f} {:5.5f} \n'.format(self.permittivity,
                self.permeability, self.loss_tangent)
        mstring += '  ),\n'

        mstring += '  nodes             : table\n'
        mstring += '    (\n'

        for i, (ri, zi) in enumerate(zip(self.r, self.zag)):
            mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(i+1, ri, zi)

        mstring += '    ),\n'

        mstring += '  linear_segments   : table\n'
        mstring += '    (\n'
        if self.PEC and self.closed:
            reg2 = -1
        elif self.PEC:
            reg2 = 0
        else:
            reg2 = 1

        for i in range(len(self.r)-1):
            mstring += '      {:d}    {:d}    {:d}    {:d}    {:d}    {:.5e}    {:.5e}\n'.format(
                i+1, i+1, i+2, 0, reg2, self.Zs_real, self.Zs_imag)

        if self.closed:
            mstring += '      {:d}    {:d}    {:d}    {:d}    {:d}    {:.5e}    {:.5e}\n'.format(
                i+2, i+2, 1, 0, reg2, self.Zs_real, self.Zs_imag)

        mstring += '    ),\n'

        mstring += '  length_unit       : {:s},\n'.format(self.units)
        mstring += '  coor_order        : rho_z,\n'
        mstring += '  advanced_regions  : table\n'
        mstring += '    (\n'
        mstring += '    )\n'

        mstring += ')\n\n'

        if self.xtrans is not None and self.ytrans is not None:

            mstring += '{:s}_trans.csy coor_sys  \n'.format(self.momcsy)
            mstring += '(\n'

            mstring += '  base             : ref({:s}.csy),\n'.format(self.momcsy)
            mstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.xtrans, self.units, self.ytrans, self.units)+'z: {:5.5f} {:s})\n'.\
                format(0.0, self.units)
            mstring += ')\n\n'

        mstring += '{:s}_mom.csy coor_sys  \n'.format(self.name)
        mstring += '(\n'
        if self.base is not None:
            mstring += '  base             : ref({:s}),\n'.format(self.base)
        if self.invert_x:
            mstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z:
            mstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        mstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        mstring += ')\n\n'

        return mstring


    def grasp_string(self):

        gstring = ''
        if self.elevation is not None:

            gstring += '{:s}_cylrot.csy coor_sys_grasp_angles\n'.format(self.name)
            gstring += '(\n'
            gstring += '  theta         : {:5.5f}, \n'.format(self.elevation)
            # gstring += '  phi           : {:5.5f}, \n'.format(self.phi)
            # gstring += '  psi           : {:5.5f},'.format(self.psi)
            if self.base is not None:
                gstring += '  base          : ref({:s})\n'.format(self.base)
            gstring += ')\n'

        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.\
            format((self.base if self.elevation is None else '{:s}_cylrot.csy'.format(self.name)))
        if self.invert_x:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.angle == 0.:
            zoff = self.r1/np.tan(np.radians(self.angle_fudge))

        else:
            zoff = self.r1/np.tan(np.radians(self.angle))

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z - zoff, self.units)
        gstring += ')\n\n'

        if self.spline_reflector:

            gstring += '{:s}.scr  spline_circ_sym_reflector\n'.format(self.name)
            gstring +=   '(\n'
            gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
            gstring += '  length_unit       : {:s},\n'.format(self.units)
            gstring += '  nodes             : table\n'
            gstring += '    (\n'

            for i, (ri, zi) in enumerate(zip(self.r, self.zag)):
                gstring += '         {:.5e}    {:.5e}\n'.format(ri, zi)

            gstring += '    ),\n'
            gstring += ')\n\n'

            gstring += '{:s}.po  po_multi_face_scatterer\n'.format(self.name)
            gstring += '(\n'
            gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
            gstring += '  po_points        : struct(po1: {:d}, po2: {:d}),\n'.\
                format(self.po1, self.po2)
            gstring += '  ptd_points       : sequence\n'
            gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
            gstring += '    ),\n'
            gstring += '  spill_over       : {:s},\n'.format(\
                'on' if self.spillover else 'off')
            gstring += '  scatterer        : ref({:s}.scr)\n'.format(self.name)
            gstring += ')\n\n'

            return gstring

        # Using a circular strut to capture geometry
        elif self.circular_strut:

            pass

        # Using a single face scatterer that combines a circular cone and an elliptical rim [default]
        else:

            half_axis = self.height * np.tan(np.radians(self.angle))+self.r1

            gstring += '{:s}.srf  circular_cone\n'.format(self.name)
            gstring += '(\n'

            if self.angle == 0.:
                gstring += '  half_cone_angle  : {:.5f},\n'.format(self.angle_fudge)
            else:
                gstring += '  half_cone_angle  : {:.5f},\n'.format(self.angle)

            gstring += '  vertex           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.xv, self.units, self.yv, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.zv, self.units)

            if self.add_absorbing_back:
                gstring += '  el_prop        : sequence(ref(perfect_conductivity), ref(perfect_absorption)),\n'

            gstring += '  axis             : struct(x: {:5.5f}, y: {:5.5f}, '.format(
                self.xa, self.ya)+'z: {:5.5f})\n'.format(self.za)
            gstring += ')\n'

            # Rim
            gstring += '\n{:s}.rim  elliptical_rim\n'.format(self.name)
            gstring += '(\n'
            gstring += '  half_axis        : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
                format(half_axis, self.units, half_axis, self.units)
            gstring += ')\n'
            gstring += ' \n'

            # if self.absorbing:
            #     gstring += '\n{:s}_el_prop  perfect_absorption\n'.format(self.name)
            #     gstring += '(\n'
            #     gstring += ')\n'
            #     gstring += ' \n'

            # Scatterer
            gstring += '{:s}.scr  reflector\n'.format(self.name)
            gstring += '(\n'
            gstring += '  coor_sys            : ref({:s}.csy),\n'.format(self.name)
            gstring += '  centre_hole_radius  : {:5.4f} {:s},\n'.format(self.r1, self.units)
            gstring += '  surfaces            : sequence(ref({:s}.srf)),\n'.format(self.name)
            if self.absorbing:
                gstring += '  el_prop            : sequence(ref(perfect_absorption)),\n'

            gstring += '  rim                 : ref({:s}.rim)\n'.format(self.name)
            gstring += ')\n\n'

            # PO object
            gstring += '{:s}.po  po_single_face_scatterer\n'.format(self.name)
            gstring += '(\n'
            gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
            gstring += '  po_points        : struct(po1: {:d}, po2: {:d}),\n'.\
                format(self.po1, self.po2)
            gstring += '  ptd_points       : sequence\n'
            gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
            gstring += '    ),\n'
            gstring += '  spill_over       : {:s},\n'.format(\
                'on' if self.spillover else 'off')
            gstring += '  scatterer        : ref({:s}.scr)\n'.format(self.name)
            gstring += ')\n\n'

            if self.add_pecabs:

                gstring += 'perfect_conductivity  perfect_conductivity\n'
                gstring += '(\n'
                gstring += '  displacement     : 1.0 mm \n'
                gstring += ')\n'

                gstring += 'perfect_absorption  perfect_absorption\n'
                gstring += '(\n'
                gstring += ')\n\n'

            return gstring

            # circular_cone  circular_cone
            # (
            #   half_cone_angle  : 20.0,
            #   vertex           : struct(x: 0.0 m, y: 0.0 m, z: 0.0 m),
            #   axis             : struct(x: 0.0, y: 0.0, z: 2.0)
            # )

            # elliptical_rim  elliptical_rim
            # (
            #   half_axis        : struct(x: 0.3 m, y: 0.3 m)
            # )

            # reflector  reflector
            # (
            #   surfaces         : sequence(ref(circular_cone)),
            #   rim              : ref(elliptical_rim),
            #   centre_hole_radius : 0.2 m,
            #   serration        : struct(inner_rim: ref(), shape: linear)
            # )

            # po_reflector  po_single_face_scatterer
            # (
            #   frequency        : ref(freq),
            #   scatterer        : ref(reflector)
            # )

class CadObject(object):

    def __init__(self, name=None):

        self.name = name
        self.poname = name

class Cut():

    def __init__(self,name=None):

        pass

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}  aperture_in_screen\n'.format(self.name)
        gstring += '(\n'
        gstring += '  theta_range         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  phi_range           : ref({:s}_elliptical_rim)\n'.format(self.name)
        gstring += ')'

class Grid():

    def __init__(self,name=None):
        pass

    def grasp_string(self):

        pass

class Aperture():

    def __init__(self, name=None, fancy_name=None,
            units='mm', frequency_list=None,
            field_accuracy=-80, auto_convergence='on',
            conv_grid=None, po1=50, po2=50, ptd=50, spillover=False,
            d=10, x=0, y=0, z=0, base=None):

        self.name = name
        self.fancy_name = fancy_name
        self.poname = name+'.po'
        self.units = units
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.po1 = po1
        self.po2 = po2
        self.ptd = ptd
        self.spillover = spillover
        self.frequency_list = frequency_list
        self.d = d
        self.x, self.y, self.z = x, y, z
        self.base = base


    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Printing out aperture information:')
        print('Aperture diameter {:5.2f} {:s}'.format(self.d, self.units))
        print('Aperture position x/y/z {:3.2f}/{:3.2f}/{:3.2f} {:s}'.\
            format(self.x, self.y, self.z, self.units))

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}.scr  aperture_in_screen\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  rim              : ref({:s}_elliptical_rim)\n'.format(self.name)
        gstring += ')\n'

        gstring += '\n{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.\
            format(self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        gstring += ')\n'

        gstring += '\n{:s}_elliptical_rim  elliptical_rim\n'.format(self.name)
        gstring += '(\n'

        gstring += '  half_axis        : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
            format(self.d/2., self.units, self.d/2., self.units)
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}  po_aperture_in_screen\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  po_points        : struct(po1: {:d}, po2: {:d}),\n'.\
            format(self.po1, self.po2)
        gstring += '  ptd_points       : sequence\n'
        gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
        gstring += '    ),\n'
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  scatterer        : ref({:s}.scr)\n'.format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        return gstring

class Aperture2():

    def __init__(self, name=None, fancy_name=None,
            units='mm', frequency_list=None,
            field_accuracy=-80, auto_convergence='on',
            conv_grid=None, po1=50, po2=50, ptd=50, spillover=False,
            hexagonal=False, hexang=0, flip_trans=True,
            d=10, x=0, y=0, z=0, rot_matrix=None, base=None):

        self.name = name
        self.fancy_name = fancy_name
        self.poname = name+'.po'
        self.units = units
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.po1 = po1
        self.po2 = po2
        self.ptd = ptd
        self.hexagonal = hexagonal
        self.hexang = hexang
        self.spillover = spillover
        self.frequency_list = frequency_list
        self.d = d
        self.x, self.y, self.z = x, y, z
        self.rot_matrix = rot_matrix
        self.base = base
        self.flip_trans = flip_trans


    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Printing out aperture information:')
        print('Aperture diameter {:5.2f} {:s}'.format(self.d, self.units))
        print('Aperture position x/y/z {:3.2f}/{:3.2f}/{:3.2f} {:s}'.\
            format(self.x, self.y, self.z, self.units))

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}.scr  aperture_in_screen\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        if self.hexagonal:
            gstring += '  rim              : ref({:s}_polygonal_rim)\n'.format(self.name)
        else:
            gstring += '  rim              : ref({:s}_elliptical_rim)\n'.format(self.name)

        gstring += ')\n'


        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z, self.units)

        if self.flip_trans:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
        gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
        gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_rot.csy)\n'.\
            format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        if self.hexagonal:

            r = self.d/2.
            hexa = self.hexang

            gstring += '\n{:s}_polygonal_rim  polygonal_rim\n'.format(self.name)
            gstring += '(\n'
            gstring += '  length_unit  : {:s}, \n'.format(self.units)
            gstring += '  nodes        : table \n'
            gstring += '  ( \n'
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(30+hexa)), r*np.sin(np.deg2rad(30+hexa)))
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(90+hexa)), r*np.sin(np.deg2rad(90+hexa)))
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(150+hexa)), r*np.sin(np.deg2rad(150+hexa)))
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(210+hexa)), r*np.sin(np.deg2rad(210+hexa)))
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(270+hexa)), r*np.sin(np.deg2rad(270+hexa)))
            gstring += '      {:.6e}    {:.6e}  \n'.\
                format(r*np.cos(np.deg2rad(330+hexa)), r*np.sin(np.deg2rad(330+hexa)))
            gstring += '  ) \n'
            gstring += ')\n'
            gstring += ' \n'


        else:
            gstring += '\n{:s}_elliptical_rim  elliptical_rim\n'.format(self.name)
            gstring += '(\n'
            gstring += '  half_axis        : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
                format(self.d/2., self.units, self.d/2., self.units)
            gstring += ')\n'
            gstring += ' \n'

        gstring += '{:s}  po_aperture_in_screen\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  po_points        : struct(po1: {:d}, po2: {:d}),\n'.\
            format(self.po1, self.po2)
        gstring += '  ptd_points       : sequence\n'
        gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
        gstring += '    ),\n'
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  scatterer        : ref({:s}.scr)\n'.format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        return gstring

class RectangularAperture():

    def __init__(self, name=None, fancy_name=None,
            units='mm', frequency_list=None,
            auto_convergence='on', field_accuracy=-80,
            conv_grid=None, ptd=50, spillover=False,
            xw=8, yw=8, x=0, y=0, z=0, base=None):

        self.name = name
        self.fancy_name = fancy_name
        self.poname = name+'.po'
        self.units = units
        self.frequency_list = frequency_list
        self.auto_convergence = auto_convergence
        self.field_accuracy = field_accuracy
        self.conv_grid = conv_grid
        self.ptd = ptd
        self.spillover = spillover
        self.xw = xw
        self.yw = yw
        self.x, self.y, self.z = x, y, z
        self.base = base


    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Printing out aperture information:')
        print('Rectangular Aperture {:5.2f} {:s}'.format(self.xw, self.units))
        print('Rectangular Aperture position x/y/z {:3.2f}/{:3.2f}/{:3.2f} {:s}'.\
            format(self.x, self.y, self.z, self.units))

    def grasp_string(self):

        #gstring = '\n'
        gstring = '{:s}  aperture_in_screen\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  rim              : ref({:s}_rectangular_rim)\n'.format(self.name)
        gstring += ')\n'

        gstring += '\n{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base         : ref({:s}.csy),\n'.\
                format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '\
        .format(self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        gstring += ')\n'

        gstring += '\n{:s}_rectangular_rim  rectangular_rim\n'.format(self.name)
        gstring += '(\n'
        gstring += '  side_lengths     : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
            format(self.xw, self.units, self.yw, self.units)
        gstring += ')\n\n'

        gstring += '{:s}  po_aperture_in_screen\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  ptd_points       : sequence\n'
        gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
        gstring += '    ),\n'
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  scatterer        : ref({:s})\n'.format(self.name)
        gstring += ')\n\n'

        return gstring

class PiecewiseLBR():
    '''
    This is a MoM object only
    '''

    def __init__(self, name=None, fancy_name=None, units='mm',
        frequency_list=None, d=10, momcsy='mom',
        ztrans=0.,
        permittivity=1.0, permeability=1.0, loss_tangent=0.0,
        Zs_real=0, Zs_imag=0.,
        Np=31, PEC=True, lossy_conductor=False, add_pecabs=False,
        x=0, y=0, r1=0, r2=0, z1=0, z2=0,
        nodes=None, radius=None,
        z=0, surf1=None, surf2=None, auto_convergence='on',
        field_accuracy=-80, spillover=False, conv_grid=None,
        base=None, invert_x=False, invert_z=False, rot_matrix=np.identity(3),
        tol=1e-6):

        self.name = name
        self.poname = name+'.po'
        self.fancy_name = fancy_name
        self.units = units
        self.frequency_list = frequency_list
        self.auto_convergence = auto_convergence
        self.field_accuracy = field_accuracy
        self.spillover = spillover
        # self.conv_grid = conv_grid
        if conv_grid is not None and not isinstance(conv_grid, list):
            self.conv_grid = [conv_grid]
        else:
            self.conv_grid = conv_grid

        self.base = base

        self.x = x
        self.y = y
        self.z = z
        self.r1 = r1
        self.r2 = r2
        self.z1 = z1
        self.z2 = z2

        self.nodes = nodes
        self.radius = radius

        self.PEC = PEC
        self.lossy_conductor = lossy_conductor

        self.ztrans = ztrans
        self.momcsy = momcsy
        self.tol = tol

        self.permittivity = permittivity
        self.permeability = permeability
        self.loss_tangent = loss_tangent
        self.Zs_real = Zs_real
        self.Zs_imag = Zs_imag
        # self.absorbing = absorbing

        self.add_pecabs = add_pecabs
        self.Np = Np

    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Piecewise linear body of rotation:')


    def mom_mesh(self):

        if self.nodes is None:

            r = np.linspace(self.r1, self.r2, self.Np)
            z = np.linspace(self.z1, self.z2, self.Np)

        else:

            r = self.radius * np.ones_like(self.nodes[:, 0])
            #np.sqrt((self.nodes[:, 0]-self.x)**2 + (self.nodes[:, 1]-self.y)**2)
            # print('r')
            # print(r)
            # print(self.nodes[:,0])
            # print(self.nodes[:,1])
            # print(self.x)
            # print(self.y)

            z = self.nodes[:, 2]

        self.r = r
        self.zag = self.ztrans + z

        # self.r = np.hstack((r, rinv))
        # self.zag = self.ztrans + np.hstack((z1, z1[::-1] - 0.1))

    def mom_string(self):

        mstring = '{:s}.msh bor_mesh\n'.format(self.name)
        mstring += '(\n'
        mstring += '  coor_sys          : ref({:s}.csy),\n'.format(self.momcsy)
        mstring += '  regions           : table\n'
        mstring += '  (\n'

        if not self.PEC and not self.lossy_conductor:
            mstring += '    1    {:5.5f} {:5.5f} {:5.5f} \n'.format(self.permittivity,
                self.permeability, self.loss_tangent)

        mstring += '  ),\n'

        mstring += '  nodes             : table\n'
        mstring += '    (\n'

        for i, (ri, zi) in enumerate(zip(self.r, self.zag)):
            mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(i+1, ri, zi)

        mstring += '    ),\n'

        mstring += '  linear_segments   : table\n'
        mstring += '    (\n'

        # if self.PEC and self.closed:
        #     reg2 = -1
        if self.PEC or self.lossy_conductor:
            reg2 = 0
        else:
            reg2 = 1

        for i in range(len(self.r)-1):
            mstring += '      {:d}    {:d}    {:d}    {:d}    {:d}    {:.5e}    {:.5e}\n'.format(
                i+1, i+1, i+2, 0, reg2, self.Zs_real, self.Zs_imag)

        # if self.closed:
        #     mstring += '      {:d}    {:d}    {:d}    {:d}    {:d}    {:.5e}    {:.5e}\n'.format(
        #         i+2, i+2, 1, 0, reg2, self.Zs_real, self.Zs_imag)

        mstring += '    ),\n'

        mstring += '  length_unit       : {:s},\n'.format(self.units)
        mstring += '  coor_order        : rho_z,\n'
        mstring += '  advanced_regions  : table\n'
        mstring += '    (\n'
        mstring += '    )\n'

        mstring += ')\n\n'

        mstring += '{:s}_mom.csy coor_sys  \n'.format(self.name)
        mstring += '(\n'
        if self.base is not None:
            mstring += '  base             : ref({:s}),\n'.format(self.base)

        mstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)
        mstring += ')\n\n'

        return mstring

class Reflector():

    def __init__(self, name=None, fancy_name=None,
            surface=None, units='mm', frequency_list=None,
            radius=20,
            x=0, y=0, z=0,
            ha_x=0, ha_y=0, xc=0, yc=0,
            po1=50, po2=50, ptd=50, spillover=False,
            elliptical_rim=True,
            rectangular_rim=False,
            re_x=0, re_y=0,
            xygrid=True,
            paraboloid=False, hyperboloid=False,
            parab_focal_length=1, hyperb_vertex_distance=-1, hyperb_foci_distance=100,
            field_accuracy=-80, auto_convergence='on',
            base=None, add_absorbing_back=False,
            add_pecabs=True,
            conductivity=2.5e7,
            finite_conductivity=False,
            add_fec=False,
            fec_name='finite_conductivity'):

        self.name = name
        self.fancy_name = fancy_name
        self.poname = name+'.po'
        self.units = units
        self.surface = surface
        self.frequency_list = frequency_list
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = None
        self.po1 = po1
        self.po2 = po2
        self.ptd = ptd
        self.spillover = spillover
        self.radius = radius
        self.x, self.y, self.z = x, y, z
        self.ha_x = ha_x
        self.ha_y = ha_y
        self.re_x = re_x
        self.re_y = re_y
        self.d = np.sqrt(self.ha_x * self.ha_y)
        self.xc = xc
        self.yc = yc
        self.base = base
        self.elliptical_rim = elliptical_rim
        self.rectangular_rim = rectangular_rim
        self.add_pecabs = add_pecabs

        if self.rectangular_rim and self.elliptical_rim:
            self.elliptical_rim = False

        self.xygrid = xygrid
        self.paraboloid = paraboloid
        self.hyperboloid = hyperboloid

        if sum([self.xygrid, self.paraboloid, self.hyperboloid]) > 1:
            raise ValueError('Reflector can only have one surface type')

        self.parab_focal_length = parab_focal_length
        self.hyperb_vertex_distance = hyperb_vertex_distance
        self.hyperb_foci_distance = hyperb_foci_distance

        self.add_absorbing_back = add_absorbing_back

        self.finite_conductivity = finite_conductivity
        self.conductivity = conductivity
        self.add_fec = add_fec
        self.fec_name = fec_name

    def grasp_string(self):

        gstring = '{:s}.scr  reflector\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  surfaces         : sequence(ref({:s}.srf)),\n'.format(self.name)
        if self.add_absorbing_back:
            gstring += '  el_prop        : sequence(ref(perfect_conductivity), ref(perfect_absorption)),\n'

        if self.finite_conductivity:
            gstring += '  el_prop        : sequence(ref({:s})),\n'.format(self.fec_name)

        gstring += '  rim              : ref({:s}.rim)\n'.format(self.name)
        gstring += ')\n\n'

        if self.elliptical_rim:

            gstring += '{:s}.rim  elliptical_rim\n'.format(self.name)
            gstring += '(\n'
            gstring += '  centre           : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s}),\n'.\
                format(self.xc, self.units, self.yc, self.units)
            gstring += '  half_axis        : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
                format(self.ha_x, self.units, self.ha_y, self.units)
            gstring += ')\n'
        elif self.rectangular_rim:
            gstring += '{:s}.rim  rectangular_rim\n'.format(self.name)
            gstring += '(\n'
            gstring += '  centre           : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s}),\n'.\
                format(self.xc, self.units, self.yc, self.units)
            gstring += '  side_lengths     : struct(x: {:5.4f} {:s},  y: {:5.5f} {:s})\n'.\
                format(self.re_x, self.units, self.re_y, self.units)
            gstring += ')\n\n'

        if self.xygrid:
            gstring += '{:s}.srf  regular_xy_grid\n'.format(self.name)
            gstring += '(\n'
            gstring += '  file_name        : {:s}.sfc,\n'.format(self.name)
            gstring += '  xy_unit          : {:s},\n'.format(self.units)
            gstring += '  z_unit           : {:s}\n'.format(self.units)
            gstring += ')\n'

        elif self.paraboloid:

            gstring += '{:s}.srf  paraboloid\n'.format(self.name)
            gstring += '(\n'
            gstring += '  focal_length        : {:.5f} {:s} \n'.\
                format(self.parab_focal_length, self.units)
            gstring += ')\n'

        elif self.hyperboloid:

            gstring += '{:s}.srf  hyperboloid\n'.format(self.name)
            gstring += '(\n'
            gstring += '  vertex_distance        : {:.5f} {:s}, \n'.\
                format(self.hyperb_vertex_distance, self.units)
            gstring += '  foci_distance          : {:.5f} {:s}, \n'.\
                format(self.hyperb_foci_distance, self.units)
            gstring += ')\n'

        gstring += '{:s}.po  po_single_face_scatterer\n'.format(self.name)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  po_points        : struct(po1: {:d}, po2: {:d}),\n'.\
            format(self.po1, self.po2)
        gstring += '  ptd_points       : sequence\n'
        gstring += '    (    struct(edge: -1, ptd: {:d})\n'.format(self.ptd)
        gstring += '    ),\n'
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  scatterer        : ref({:s}.scr)\n'.format(self.name)
        gstring += ')\n\n'

        if self.add_pecabs:
            gstring += 'perfect_conductivity perfect_conductivity\n'
            gstring += '(\n'
            gstring += '  displacement     : 1.0 mm \n'
            gstring += ')\n'

            gstring += 'perfect_absorption  perfect_absorption\n'
            gstring += '(\n'
            gstring += ')\n\n'

        if self.add_fec:
            gstring += '{:s} finite_conductivity\n'.format(self.fec_name)
            gstring += '(\n'
            gstring += '  conductivity        : {:.5e} S/m,\n'.format(self.conductivity)
            gstring += ')\n\n'

        return gstring

class PlaneWave():

    def __init__(self, name=None, fancy_name=None, units='mm',
        frequency_list=None, radius=20, x=0, y=0, z=0, rmin=0, rmax=2500,
        pmin=0, pmax=360,
        rstart=0., rend=1., npr=11, npp=11, rot_matrix=None, base=None,
        invert_x=False, invert_z=False):

        self.name = name
        self.fancy_name = fancy_name
        self.poname = name+'.po'
        self.units = units
        self.frequency_list = frequency_list
        self.radius = radius
        self.x, self.y, self.z = x, y, z
        self.rmin = rmin
        self.rmax = rmax
        self.pmin = pmin
        self.pmax = pmax
        self.rstart = rstart
        self.rend = rend
        self.npr = npr
        self.npp = npp
        self.invert_x = invert_x
        self.invert_z = invert_z
        self.rot_matrix = rot_matrix
        self.base = base

    def grasp_string(self):

        gstring = '{:s}.po  plane_wave\n'.format(self.name)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  aperture_radius  : {:5.5f} {:s}\n'.format(self.radius, self.units)
        gstring += ')\n\n'

        if self.rot_matrix is not None:

            gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            if self.base is not None:
                gstring += '  base             : ref({:s}),\n'.format(self.base)

            gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
                format(self.z, self.units)
            gstring += '  x_axis           : struct(x: {:5.5f}, y: 0.0, z: {:5.5f}),\n'.\
                format(self.rot_matrix[0][0], self.rot_matrix[0][2])
            gstring += '  y_axis           : struct(x: 0.0, y: 1.0, z: 0.0)\n'
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_rot2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 0.0, y: {:5.5f}, z: {:5.5f}),\n'.\
                format(self.rot_matrix[1][1], self.rot_matrix[1][2])
            gstring += '  base             : ref({:s}_rot.csy)\n'.\
                format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            if self.invert_x:
                gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            if self.invert_z:
                gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
            gstring += '  base             : ref({:s}_rot2.csy)\n'.\
                format(self.name)

        else:

            gstring += '{:s}.csy coor_sys\n'.format(self.name)
            gstring += '(\n'
            gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.\
                format(self.x, self.units, self.y, self.units) + \
                'z: {:5.5f} {:s})'.format(self.z, self.units)

            if self.base is not None:
                gstring += ', \n  base             : ref({:s})\n'.format(self.base)
            else:
                gstring += '\n'

        gstring += ')\n\n'

        return gstring

    def view_string(self):

        vstring = '{:s}.rfp  rays_from_plane_waves\n'.format(self.name)
        vstring += '(\n'
        vstring += '  objects           : sequence(ref({:s})),\n'.format(self.name)
        vstring += '  rho_range         : struct(start: {:5.5f} {:s}, end: {:5.5f} {:s}, np: {:d}),\n'.\
            format(self.rmin, self.units, self.rmax, self.units, self.npr)
        vstring += '  phi_range         : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.pmin, self.pmax, self.npp)
        vstring += '  ray_path_range    : struct(start: {:5.5f} {:s}, end: {:5.5f} {:s})\n'.\
            format(self.rstart, self.units, self.rend, self.units)
        vstring += ')\n\n'

        return vstring

class SimpleLens():
    '''
    A simple lens object encapsulates all the
    information needed to describe a simple lens. The class comes with functions
    that allow one to map lens properties from one function bases to another.
    '''

    def __init__(self, name=None, fancy_name=None, units='mm',
        frequency_list=None, d=10, t=2, n=1.5,
        r1=0, r2=0, bs1=0, bs2=0, f1po1=50, f1po2=50, f2po1=50, f2po2=50,
        x=0, y=0, z=0, surf1=None, surf2=None, auto_convergence='on',
        field_accuracy=-80, spillover=False, conv_grid=None,
        base=None, invert_x=False, invert_z=False, rot_matrix=None,
        flip_trans=True):
        '''
        Attributes

        d : diameter
        t : thickness
        n : refractive index
        r1 : radius of curvature for entry surface
        r2 : radius of curvature for exit surface
        bs1 : conic constant for entry surface
        bs2 : conic constant for exit surface
        '''

        self.name = name
        self.poname = name+'.po'
        self.fancy_name = fancy_name
        self.units = units
        self.frequency_list = frequency_list
        self.auto_convergence = auto_convergence
        self.field_accuracy = field_accuracy
        self.spillover = spillover
        self.conv_grid = conv_grid
        self.d = d
        self.n = n
        self.t = t
        self.r1 = r1 if not np.isinf(r1) else 0.0
        self.r2 = r2 if not np.isinf(r2) else 0.0
        self.f1po1, self.f1po2 = f1po1, f1po2
        self.f2po1, self.f2po2 = f2po1, f2po2
        self.bs1, self.bs2 = bs1, bs2
        self.x, self.y, self.z = x, y, z
        self.surf1 = surf1
        self.surf2 = surf2
        self.base = base
        self.invert_x = invert_x
        self.invert_z = invert_z
        self.rot_matrix = rot_matrix
        self.flip_trans = flip_trans


    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Printing out lens information:')
        print('Lens diameter {:5.2f} {:s}'.format(self.d, self.units))
        print('Lens thickness {:5.2f} {:s}'.format(self.t, self.units))
        print('Lens radius of curvature  {:5.2f}/{:5.2f} {:s}'.\
            format(self.r1, self.r2, self.units))
        print('Lens conic constant  {:5.2f}/{:5.2f}'.format(self.bs1, self.bs2))
        print('Lens position x/y/z {:3.2f}/{:3.2f}/{:3.2f} {:s}'.\
            format(self.x, self.y, self.z, self.units))

    def grasp_string(self):
        #gstring = '\n'
        gstring = '{:s}.scr  simple_lens\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  diameter         : {:3.5f} {:s},\n'.format(self.d, self.units)
        gstring += '  thickness        : {:3.5f} {:s},\n'.format(self.t, self.units)
        gstring += '  refractive_index : {:3.5f},\n'.format(self.n)
        gstring += '  r1               : {:3.5f} {:s},\n'.format(self.r1,self.units)
        gstring += '  r2               : {:3.5f} {:s},\n'.format(self.r2,self.units)
        gstring += '  bs1              : {:3.5f},\n'.format(self.bs1)
        gstring += '  bs2              : {:3.5f},\n'.format(self.bs2)

        if self.surf1 is not None:
            gstring += '  surface1_file    : {:s},\n'.format(self.surf1)
        if self.surf2 is not None:
            gstring += '  surface2_file    : {:s},\n'.format(self.surf2)

        gstring += '  length_unit_in_files : {:s}\n'.format(self.units)
        gstring += ')\n'
        gstring += ' \n'

        if self.rot_matrix is not None:

            gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            if self.base is not None:
                gstring += '  base             : ref({:s}),\n'.format(self.base)

            gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
                format(self.z, self.units)
            gstring += '  x_axis           : struct(x: {:5.5f}, y: 0.0, z: {:5.5f}),\n'.\
                format(self.rot_matrix[0][0], self.rot_matrix[0][2])
            gstring += '  y_axis           : struct(x: 0.0, y: 1.0, z: 0.0)\n'
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_rot2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 0.0, y: {:5.5f}, z: {:5.5f}),\n'.\
                format(self.rot_matrix[1][1], self.rot_matrix[1][2])
            gstring += '  base             : ref({:s}_rot.csy)\n'.\
                format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            if self.invert_x or self.invert_z:

                gstring += '{:s}_invert.csy coor_sys  \n'.format(self.name)
                gstring += '(\n'

                if self.invert_x:
                    gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
                if self.invert_z:
                    gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
                gstring += '  base             : ref({:s}_rot2.csy)\n'.\
                    format(self.name)
                gstring += ')\n'

            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            # gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            #     self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
            #     format(self.z, self.units)

            if self.invert_x or self.invert_z:
                gstring += '  base             : ref({:s}_invert.csy)\n'.\
                    format(self.name)
            else:
                gstring += '  base             : ref({:s}_rot2.csy)\n'.\
                    format(self.name)

        else:
            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            if self.base is not None:
                gstring += '  base             : ref({:s}),\n'.format(self.base)
            if self.invert_x:
                gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            if self.invert_z:
                gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

            gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
                format(self.z, self.units)

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}  po_lens\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  lens             : ref({:s}.scr),\n'.format(self.name)
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  po_points        : struct(face1_po1: {:d}, face1_po2: {:d},'.\
            format(self.f1po1, self.f1po2) + ' face2_po1: {:d}, face2_po2: {:d})\n'.\
            format(self.f2po1, self.f2po2)

        gstring += ')\n\n'

        return gstring


    def set_t(self,t):
        self.t = t

    def get_t(self):
        return self.t

class SimpleLens2(object):
    '''
    A simple lens object encapsulates all the
    information needed to describe a simple lens. The class comes with functions
    that allow one to map lens properties from one function bases to another.
    '''

    def __init__(self, name=None, fancy_name=None, units='mm',
        frequency_list=None, d=10, t=2, n=1.5, momcsy='mom',
        xtrans=None, ytrans=None, ztrans=0.,
        permeability=1.0, loss_tangent=0.0, Zs_real=-1., Zs_imag=0.,
        Np=31,
        ar1=False, ar2=False, t_ar1=1.0, t_ar2=1.0,
        ar1b=False, ar2b=False, t_ar1b=0.0, t_ar2b=0.0,
        ar1c=False, ar2c=False, t_ar1c=0.0, t_ar2c=0.0,
        permittivity_ar1=1.0, permeability_ar1=1.0, loss_tangent_ar1=0.0,
        permittivity_ar1b=1.0, permeability_ar1b=1.0, loss_tangent_ar1b=0.0,
        permittivity_ar1c=1.0, permeability_ar1c=1.0, loss_tangent_ar1c=0.0,
        permittivity_ar2=1.0, permeability_ar2=1.0, loss_tangent_ar2=0.0,
        permittivity_ar2b=1.0, permeability_ar2b=1.0, loss_tangent_ar2b=0.0,
        permittivity_ar2c=1.0, permeability_ar2c=1.0, loss_tangent_ar2c=0.0,
        r1=0, r2=0, bs1=0, bs2=0, f1po1=50, f1po2=50, f2po1=50, f2po2=50,
        x=0, y=0, z=0, surf1=None, surf2=None, auto_convergence='on',
        field_accuracy=-80, spillover=False, conv_grid=None,
        base=None, invert_x=False, invert_z=False, rot_matrix=np.identity(3),
        theta=0.0, phi=0.0, use_grasp_angles=False,
        tol=1e-6,
        flip_trans=True):
        '''
        Attributes

        d : diameter
        t : thickness
        n : refractive index
        r1 : radius of curvature for entry surface
        r2 : radius of curvature for exit surface
        bs1 : conic constant for entry surface
        bs2 : conic constant for exit surface
        '''

        self.name = name
        self.poname = name+'.po'
        self.fancy_name = fancy_name
        self.units = units
        self.frequency_list = frequency_list
        self.auto_convergence = auto_convergence
        self.field_accuracy = field_accuracy
        self.spillover = spillover
        # self.conv_grid = conv_grid

        if conv_grid is not None and not isinstance(conv_grid, list):
            self.conv_grid = [conv_grid]
        else:
            self.conv_grid = conv_grid

        self.d = float(d)
        self.n = float(n)
        self.Np = Np

        self.momcsy = momcsy

        self.xtrans = xtrans
        self.ytrans = ytrans
        self.ztrans = ztrans

        self.permittivity = n**2
        self.permeability = permeability
        self.loss_tangent = loss_tangent

        self.ar1 = ar1
        self.ar1b = ar1b
        self.ar1c = ar1c
        self.ar2 = ar2
        self.ar2b = ar2b
        self.ar2c = ar2c

        self.t_ar1 = float(t_ar1)
        self.permittivity_ar1 = permittivity_ar1
        self.permeability_ar1 = permeability_ar1
        self.loss_tangent_ar1 = loss_tangent_ar1

        self.t_ar2 = float(t_ar2)
        self.permittivity_ar2 = permittivity_ar2
        self.permeability_ar2 = permeability_ar2
        self.loss_tangent_ar2 = loss_tangent_ar2

        self.t_ar1b = float(t_ar1b)
        self.permittivity_ar1b = permittivity_ar1b
        self.permeability_ar1b = permeability_ar1b
        self.loss_tangent_ar1b = loss_tangent_ar1b

        self.t_ar1c = float(t_ar1c)
        self.permittivity_ar1c = permittivity_ar1c
        self.permeability_ar1c = permeability_ar1c
        self.loss_tangent_ar1c = loss_tangent_ar1c

        self.t_ar2b = float(t_ar2b)
        self.permittivity_ar2b = permittivity_ar2b
        self.permeability_ar2b = permeability_ar2b
        self.loss_tangent_ar2b = loss_tangent_ar2b

        self.t_ar2c = float(t_ar2c)
        self.permittivity_ar2c = permittivity_ar2c
        self.permeability_ar2c = permeability_ar2c
        self.loss_tangent_ar2c = loss_tangent_ar2c

        self.Zs_real = Zs_real
        self.Zs_imag = Zs_imag
        self.t = t
        self.r1 = float(r1) if not np.isinf(r1) else 0.0
        self.r2 = float(r2) if not np.isinf(r2) else 0.0
        self.f1po1, self.f1po2 = f1po1, f1po2
        self.f2po1, self.f2po2 = f2po1, f2po2
        self.bs1, self.bs2 = float(bs1), float(bs2)
        self.x, self.y, self.z = x, y, z
        self.surf1 = surf1
        self.surf2 = surf2
        self.base = base
        self.invert_x = invert_x
        self.invert_z = invert_z
        self.rot_matrix = rot_matrix
        self.tol = tol

        self.use_grasp_angles = use_grasp_angles
        self.theta = 0.0
        self.phi = 0.0

        self.flip_trans = flip_trans

    def info(self):
        '''
        Prints some information about this lens to the screen
        '''

        print('Printing out lens information:')
        print('Lens diameter {:5.2f} {:s}'.format(self.d, self.units))
        print('Lens thickness {:5.2f} {:s}'.format(self.t, self.units))
        print('Lens radius of curvature  {:5.2f}/{:5.2f} {:s}'.\
            format(self.r1, self.r2, self.units))
        print('Lens conic constant  {:5.2f}/{:5.2f}'.format(self.bs1, self.bs2))
        print('Lens position x/y/z {:3.2f}/{:3.2f}/{:3.2f} {:s}'.\
            format(self.x, self.y, self.z, self.units))

    def mom_mesh(self, Np=None, Np_lens=None, verbose=False):

        if Np is not None and Np_lens is not None:
            if Np != Np_lens:
                raise ValueError('Not currently implemented')

        if Np is None:
            Np = self.Np
        else:
            self.Np = Np

        r = np.linspace(0, self.d/2.0, Np)
        rinv = r[::-1]

        # print('r', r)
        # print('r_max', np.max(r))
        # print('r_min', np.min(r))

        # print('rass')
        # print(self.tol)
        # print(self.r1)
        # print(self.r2)

        if self.r1 == 0 or abs(1/self.r1) < self.tol:
            k1 = 0
        else:
            k1 = 1/self.r1

        if self.r2 == 0 or abs(1/self.r2) < self.tol:
            k2 = 0
        else:
            k2 = 1/self.r2

        z1 = k1 * r**2/(1+np.sqrt(1-(1+self.bs1)*k1**2*r**2))
        z2 = self.t - k2 * r**2/(1+np.sqrt(1-(1+self.bs2)*k2**2*r**2))

        # print('here')
        # print(z1)
        # print(z2)
        # print(k1, k2)
        # print(self.bs1, self.bs2)
        # print(r)
        # print(self.surf1)
        # print(self.surf2)
        # print(z1)

        if self.surf1 is not None:
            r1t, z1t = read_lens_file(self.surf1)
            z1 = z1t
            r = r1t

            # print(z1t)

        if self.surf2 is not None:
            r2t, z2t = read_lens_file(self.surf2)
            z2 = self.t - z2t
            rinv = r2t[::-1]

        if len(r) != len(rinv):
            raise ValueError('Length of r is not the same as of rinv')

        if self.surf1 is not None and self.surf2 is not None:
            self.Np = len(r)

        #keyboard()

        print('AR:')
        print(self.ar1)
        print(self.ar1b)
        print(self.ar1c)
        print(self.ar2)
        print(self.ar2b)
        print(self.ar2c)
        # print('self.t:', self.t)
        # print('k1:', k1)
        # print('k2:', k2)
        # print('z1:', z1)
        # print('z2:', z2)
        # if self.ar1:
        #     print('t_ar1:', self.t_ar1)
        # if self.ar2:
        #     print('t_ar2:', self.t_ar2)

        if verbose:
            print(self.name)
            print('Even asphere lens, maximum abs(zag) = {:5.5f} {:s}'.\
                format(np.max(np.abs(z)), self.units))

        if self.ar1:
            z_ar1 = z1[::-1] - self.t_ar1
            r_ar1 = rinv.copy()

            self.z_ar1 = z_ar1 + self.ztrans
            self.r_ar1 = r_ar1

        if self.ar1b:
            z_ar1b = z1[::-1] - self.t_ar1b - self.t_ar1
            r_ar1b = rinv.copy()

            self.z_ar1b = z_ar1b + self.ztrans
            self.r_ar1b = r_ar1b

        if self.ar1c:
            z_ar1c = z1[::-1] - self.t_ar1c - self.t_ar1b - self.t_ar1
            r_ar1c = rinv.copy()

            self.z_ar1c = z_ar1c + self.ztrans
            self.r_ar1c = r_ar1c

        if self.ar2:
            z_ar2 = z2[::-1] + self.t_ar2
            r_ar2 = rinv.copy()

            self.z_ar2 = z_ar2 + self.ztrans
            self.r_ar2 = r_ar2

        if self.ar2b:
            z_ar2b = z2[::-1] + self.t_ar2b + self.t_ar2
            r_ar2b = rinv.copy()

            self.z_ar2b = z_ar2b + self.ztrans
            self.r_ar2b = r_ar2b

        if self.ar2c:
            z_ar2c = z2[::-1] + self.t_ar2c + self.t_ar2b + self.t_ar2
            r_ar2c = rinv.copy()

            self.z_ar2c = z_ar2c + self.ztrans
            self.r_ar2c = r_ar2c


        self.rr = np.hstack((r, rinv))
        self.zag = np.hstack((z1, z2[::-1]))

        self.zag += self.ztrans

        print(np.shape(z1))
        print(np.shape(z2))
        print(np.shape(z2[::-1]))
        print(np.shape(self.zag))
        print(self.zag)

        #keyboard()

        # print('zag:')
        # print(len(self.zag))
        # print(self.zag)

        # print(z1)
        # print(z2)
        # print(self.rr)
        # print(self.zag)

        # if self.ar1:
        #     print('z_ar1:')
        #     print(self.t_ar1)
        #     print(self.z_ar1)

        # if self.ar1b:
        #     print('z_ar1b:')
        #     print(self.t_ar1b)
        #     print(self.z_ar1b)

        # if self.ar1c:
        #     print('z_ar1c:')
        #     print(self.t_ar1c)
        #     print(self.z_ar1c)

        # if self.ar2:
        #     print('z_ar2:')
        #     print(self.t_ar2)
        #     print(self.z_ar2)

        # if self.ar2b:
        #     print('z_ar2b:')
        #     print(self.t_ar2b)
        #     print(self.z_ar2b)

        # if self.ar2c:
        #     print('z_ar2c:')
        #     print(self.t_ar2c)
        #     print(self.z_ar2c)

        # print('self.ztrans', self.ztrans)
        # print('\n')


    def mom_string(self):

        mstring = '{:s}.msh  bor_mesh\n'.format(self.name)
        mstring += '(\n'

        if self.xtrans is not None and self.ytrans is not None:
            mstring += '  coor_sys          : ref({:s}_trans.csy),\n'.format(self.momcsy)
        else:
            mstring += '  coor_sys          : ref({:s}.csy),\n'.format(self.momcsy)

        mstring += '  regions           : table\n'
        mstring += '  (\n'
        mstring += '    1    {:5.5f} {:5.5f} {:5.5f} \n'.format(self.permittivity,
            self.permeability, self.loss_tangent)

        nreg = 1
        if self.ar1:
            mstring += '    2    {:5.5f} {:5.5f} {:5.5f} \n'.format(self.permittivity_ar1,
            self.permeability_ar1, self.loss_tangent_ar1)
            nreg += 1

        if self.ar2:
            mstring += '    {:d}    {:5.5f} {:5.5f} {:5.5f} \n'.format(nreg+1,
                self.permittivity_ar2, self.permeability_ar2, self.loss_tangent_ar2)
            nreg += 1

        if self.ar1b:
            mstring += '    {:d}    {:5.5f} {:5.5f} {:5.5f} \n'.format(nreg+1,
                self.permittivity_ar1b, self.permeability_ar1b, self.loss_tangent_ar1b)
            nreg += 1

        if self.ar1c:
            mstring += '    {:d}    {:5.5f} {:5.5f} {:5.5f} \n'.format(nreg+1,
                self.permittivity_ar1c, self.permeability_ar1c, self.loss_tangent_ar1c)
            nreg += 1

        if self.ar2b:
            mstring += '    {:d}    {:5.5f} {:5.5f} {:5.5f} \n'.format(nreg+1,
                self.permittivity_ar2b, self.permeability_ar2b, self.loss_tangent_ar2b)
            nreg += 1

        if self.ar2c:
            mstring += '    {:d}    {:5.5f} {:5.5f} {:5.5f} \n'.format(nreg+1,
                self.permittivity_ar2c, self.permeability_ar2c, self.loss_tangent_ar2c)
            nreg += 1

        mstring += '  ),\n'

        ### nodes
        mstring += '  nodes             : table\n'
        mstring += '    (\n'

        Np = self.Np

        print(np.shape(self.rr))
        print(np.shape(self.zag))
        for i, (ri, zi) in enumerate(zip(self.rr, self.zag)):
            mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(i+1, ri, zi)
            print(ri, zi)

        #raise

        # Adding AR coating 1
        if self.ar1:
            for i, (ri, zi) in enumerate(zip(self.r_ar1, self.z_ar1)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(2*Np+i+1, ri, zi)

        # Adding AR coating 2
        if self.ar2:
            for i, (ri, zi) in enumerate(zip(self.r_ar2, self.z_ar2)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(3*Np-1 + i + 2, ri, zi)

        # Adding AR coating 1b
        if self.ar1b:
            # nadd = 4*Np-1 if self.ar1 else 3*Np-1
            nadd = (2+int(self.ar1)+int(self.ar2))*Np
            for i, (ri, zi) in enumerate(zip(self.r_ar1b, self.z_ar1b)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(nadd+i+1, ri, zi)

        # Adding AR coating 1c
        if self.ar1c:
            # nadd = 4*Np-1 if self.ar1 else 3*Np-1
            nadd = (2+int(self.ar1)+int(self.ar2)+int(self.ar1b))*Np
            for i, (ri, zi) in enumerate(zip(self.r_ar1c, self.z_ar1c)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(nadd+i+1, ri, zi)

        # Adding AR coating 2b
        if self.ar2b:
            # nadd = 4*Np-1 if self.ar1 else 3*Np-1
            nadd = (2+int(self.ar1) + int(self.ar2) + int(self.ar1b) + int(self.ar1c))*Np
            for i, (ri, zi) in enumerate(zip(self.r_ar2b, self.z_ar2b)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(nadd+i+1, ri, zi)

        # Adding AR coating 1c
        if self.ar2c:
            # nadd = 4*Np-1 if self.ar1 else 3*Np-1
            nadd = (2+int(self.ar1) + int(self.ar2) + int(self.ar1b) + int(self.ar1c) + int(self.ar2b))*Np
            for i, (ri, zi) in enumerate(zip(self.r_ar2c, self.z_ar2c)):
                mstring += '      {:d}   {:.5e}    {:.5e}\n'.format(nadd+i+1, ri, zi)



        mstring += '    ),\n'

        ### linear segments
        mstring += '  linear_segments   : table\n'
        mstring += '    (\n'

        template = '      {:d}    {:d}    {:d}    {:d}    {:d}    {:.5e}    {:.5e}\n'

        # Main dielectric
        reg1 = 1
        for i in range(2*Np-1):

            if self.ar1 and (i < Np - 1):
                reg2 = 2
            elif self.ar2 and (i > Np -1):
                if self.ar1:
                    reg2 = 3
                else:
                    reg2 = 2
            else:
                reg2 = 0

            mstring += template.format(i+1, i+1, i+2, reg1, reg2, self.Zs_real, self.Zs_imag)

        # Adding AR coating 1
        if self.ar1:

            reg1, reg2 = 0, 2
            # reg1l = 0 + 1+int(self.ar1) + int(self.ar2) + int(self.ar1b)
            reg1l = 0 + (1+int(self.ar1) + int(self.ar2) + int(self.ar1b)) if self.ar1b else 0
            mstring += template.format(2*Np+1, Np, 2*Np+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            for i in range(Np-1):
                mstring += template.format(
                    2*Np+i+2, 2*Np+i+1, 2*Np+i+2, reg1l, reg2, self.Zs_real, self.Zs_imag)

        # Adding AR coating 2
        if self.ar2:

            nadd = 3*Np if self.ar1 else 2*Np

            reg1 = 0
            reg1l = 1+int(self.ar1)+int(self.ar2)+int(self.ar1b)+int(self.ar1c)+int(self.ar2b) if self.ar2b else 0
            reg2 = 3 if self.ar1 else 2
            mstring += template.format(nadd+1, Np+1, nadd+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            for i in range(Np-1):
                mstring += template.format(nadd+i+2, nadd+i+1, nadd+i+2, reg1l, reg2, self.Zs_real, self.Zs_imag)

        if self.ar1b:

            ioff = int(self.ar1) + int(self.ar2)
            nadd = (2+ioff)*Np
            nadd1 = (1+int(self.ar1))*Np

            reg1 = 0
            reg2 = (2 + ioff)
            reg1l = 0 if not self.ar1c else  (2 + ioff+1)

            # side
            mstring += template.format(nadd+1, nadd+1, nadd1+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            # layer
            for i in range(Np-1):
                mstring += template.format(nadd+i+2, nadd+i+1, nadd+i+2, reg1l, reg2, self.Zs_real, self.Zs_imag)

        if self.ar1c:

            ioff = int(self.ar1) + int(self.ar2) + int(self.ar1b)
            nadd = (2+ioff)*Np
            nadd1 = (1 + int(self.ar1) + int(self.ar1b) + int(self.ar2))*Np

            reg1 = 0
            reg2 = (2 + ioff)

            # side
            mstring += template.format(nadd+1, nadd+1, nadd1+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            # layer
            for i in range(Np-1):
                mstring += template.format(nadd+i+2, nadd+i+1, nadd+i+2, reg1, reg2, self.Zs_real, self.Zs_imag)

        if self.ar2b:

            ioff = int(self.ar1) + int(self.ar2) + int(self.ar1b) + int(self.ar1c)
            nadd = (2+ioff)*Np
            nadd1 = (1 + int(self.ar1) + int(self.ar2))*Np

            reg1 = 0
            reg2 = (2 + ioff)
            reg1l = 0 if not self.ar1c else (2 + ioff+1)

            # side
            mstring += template.format(nadd+1, nadd+1, nadd1+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            # layer
            for i in range(Np-1):
                mstring += template.format(nadd+i+2, nadd+i+1, nadd+i+2, reg1l, reg2, self.Zs_real, self.Zs_imag)

        if self.ar2c:

            ioff = int(self.ar1) + int(self.ar2) + int(self.ar1b) + int(self.ar1c) + int(self.ar2b)
            nadd = (2+ioff)*Np
            nadd1 = (1 + int(self.ar1) + int(self.ar1b) + int(self.ar1c) + int(self.ar2) + int(self.ar2b))*Np

            reg1 = 0
            reg2 = (2 + ioff)

            # side
            mstring += template.format(nadd+1, nadd+1, nadd1+1, reg1, reg2, self.Zs_real, self.Zs_imag)
            # layer
            for i in range(Np-1):
                mstring += template.format(nadd+i+2, nadd+i+1, nadd+i+2, reg1, reg2, self.Zs_real, self.Zs_imag)


        mstring += '    ),\n'

        mstring += '  length_unit       : {:s},\n'.format(self.units)
        mstring += '  coor_order        : rho_z,\n'
        mstring += '  advanced_regions  : table\n'
        mstring += '    (\n'
        mstring += '    )\n'

        mstring += ')\n\n'

        if self.xtrans is not None and self.ytrans is not None:

            mstring += '{:s}_trans.csy coor_sys  \n'.format(self.momcsy)
            mstring += '(\n'

            mstring += '  base             : ref({:s}.csy),\n'.format(self.momcsy)
            mstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
                self.xtrans, self.units, self.ytrans, self.units)+'z: {:5.5f} {:s})\n'.\
                format(0.0, self.units)
            mstring += ')\n\n'


        # mstring += '{:s}.mcsy coor_sys  \n'.format(self.name)
        # mstring += '(\n'
        # if self.base is not None:
        #     mstring += '  base             : ref({:s}),\n'.format(self.base)
        # if self.invert_x:
        #     mstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # if self.invert_z:
        #     mstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        # mstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
        #     self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
        #     format(self.z, self.units)
        # mstring += ')\n\n'

        return mstring


    def grasp_string(self):
        #gstring = '\n'
        gstring = '{:s}.scr  simple_lens\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  diameter         : {:3.5f} {:s},\n'.format(self.d, self.units)
        gstring += '  thickness        : {:3.5f} {:s},\n'.format(self.t, self.units)
        gstring += '  refractive_index : {:3.5f},\n'.format(self.n)
        gstring += '  r1               : {:3.5f} {:s},\n'.format(self.r1,self.units)
        gstring += '  r2               : {:3.5f} {:s},\n'.format(self.r2,self.units)
        gstring += '  bs1              : {:3.5f},\n'.format(self.bs1)
        gstring += '  bs2              : {:3.5f},\n'.format(self.bs2)

        if self.surf1 is not None:
            gstring += '  surface1_file    : {:s},\n'.format(self.surf1)
        if self.surf2 is not None:
            gstring += '  surface2_file    : {:s},\n'.format(self.surf2)

        gstring += '  length_unit_in_files : {:s}\n'.format(self.units)
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z, self.units)


        if self.flip_trans:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        if self.use_grasp_angles:
            gstring += '{:s}.csy coor_sys_grasp_angles  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}_trans.csy ),\n'.format(self.name)
            gstring += '  theta            : {:.5f},\n'.format(self.theta)
            gstring += '  phi            : {:.5f},\n'.format(self.phi)
            gstring += '  psi            : {:.5f}\n'.format(90)
            gstring += ')\n'
            gstring += ' \n'


        else:

            gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
            gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
                format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
            gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
                format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
            gstring += ')\n'
            gstring += ' \n'

        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_rot.csy)\n'.\
            format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}  po_lens\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  lens             : ref({:s}.scr),\n'.format(self.name)
        gstring += '  spill_over       : {:s},\n'.format(\
            'on' if self.spillover else 'off')
        gstring += '  po_points        : struct(face1_po1: {:d}, face1_po2: {:d},'.\
            format(self.f1po1, self.f1po2) + ' face2_po1: {:d}, face2_po2: {:d})\n'.\
            format(self.f2po1, self.f2po2)

        gstring += ')\n\n'

        return gstring


    def set_t(self,t):
        self.t = t

    def get_t(self):
        return self.t

class FiniteConductivity(object):

    def __init__(self, name='finite_conductivity', conductivity=2.5e7):
        self.name = name
        self.conductivity = conductivity

    def grasp_string():
        gstring = '{:s} finite_conductivity\n'.format(self.name)
        gstring += '(\n'
        gstring += '  conductivity        : {:.5e} S/m,\n'.format(self.conductivity)
        gstring += ')\n\n'

        return gstring

class Lens():

    def __init__(self, R, bs, Rmax, Np, a2=0, a4=0, a6=0, a8=0, a10=0):

        self.R = R
        self.bs = bs
        self.Rmax = Rmax
        self.Np = Np

    #SimpleLens.__init__(self)

class Lens_EvenAsphere(SimpleLens2):
    '''
    See definition in Zemax manual
    '''

    def __init__(self,
        a1=0, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, invert=False,
        verbose=False, **kwargs):

        super(Lens_EvenAsphere, self).__init__(**kwargs)

        self.a1, self.a2, self.a3, self.a4 = a1, a2, a3, a4
        self.a5, self.a6, self.a7, self.a8 = a5, a6, a7, a8

        r = np.linspace(0, self.d/2.0, self.Np)
        k = 1/self.r1

        self.r = r
        self.k = k

        z = k * r**2/(1+np.sqrt(1-(1+self.bs1)*k**2*r**2)) + a1*r**2 + a2*r**4 \
            + a3*r**6 + a4*r**8 + a5*r**10 + a6*r**12 + a7*r**14 + a8*r**16

        self.invert = invert
        if self.invert:
            z = -z

        if verbose:
            print('{:s}, even asphere lens, maximum abs(zag) = {:5.5f} {:s}'.\
                format(self.name, np.max(np.abs(z)),self.units))

        self.z = z
        self.zmax = np.max(np.abs(z))

    def write_surf(self, fname=None):

        if fname is None and self.name is not None:
            fname = self.name

        write_lens_file(self.r, self.z, fname)

class Lens_OddAsphere():
    '''
    See definition in Zemax manual
    '''

    def __init__(self, name=None, units='mm', d=150, r1=0, bs1=0, Np=1001,
        a1=0, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, invert=False,
        verbose=False):

        self.name = name
        self.units = units
        self.d = d
        self.r1 = r1 if not np.isinf(r1) else 0.0

        self.a1, self.a2, self.a3, self.a4 = a1, a2, a3, a4
        self.a5, self.a6, self.a7, self.a8 = a5, a6, a7, a8

        r = np.linspace(0, d/2.0, Np)
        k = 1/r1

        self.r = r
        self.k = k

        z = k * r**2/(1+np.sqrt(1-(1+bs1)*k**2*r**2)) + a1*r**1 + a2*r**2 \
            + a3*r**3 + a4*r**4 + a5*r**5 + a6*r**6 + a7*r**7 + a8*r**8

        if verbose:

            print('{:s}, odd asphere lens, maximum abs(zag) = {:5.5f} {:s}'.\
                format(self.name, np.max(np.abs(z)),self.units))

        self.invert = invert
        if self.invert:
            z = -z

        self.z = z

    def write_surf(self, fname=None):

        if fname is None and self.name is not None:
            fname = self.name

        write_lens_file(self.r, self.z, fname)

class HexagonalRim():
    def __init__(self, name=None, units='mm', r=150):

        self.name = name
        self.units = units
        self.r = r

    def write_rim():
        pass

class Horn():
    '''
    This is a Gaussian beam pattern
    '''

    def __init__(self, name=None, units='mm', frequency_list=None,
        taper_angle=25, taper=-3, field_accuracy=-80, auto_convergence='on',
        conv_grid=None, x=0, y=0, z=0, base=None,
        invert=False, invert_x=False, invert_z=False, factor=None,
        ypol=False, az0=None, el0=None, tol=1e-6):

        self.name = name
        self.poname = name+'.po'
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.taper_angle = taper_angle
        self.taper = taper
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.x = x
        self.y = y
        self.z = z

        self.az0 = az0
        self.el0 = el0
        self.tol = tol


        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z
        self.ypol = ypol

        self.factor = factor

    def grasp_string(self, invert=False):

        gstring = '\n'
        gstring += '{:s}  gaussian_beam_pattern\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  taper_angle      : {:5.5f},\n'.format(self.taper_angle)
        if self.factor is not None:
            # factor           : struct(db: 10.0, deg: 0.0)
            gstring += '  factor           : struct(db: {:.5f}, deg: 0.0),\n'.format(self.factor)

        if self.ypol:
            gstring += '  polarisation     : linear_y,\n'

        gstring += '  taper            : {:5.5f}\n'.format(self.taper)
        gstring += ')\n'
        gstring += ' \n'

        if self.az0 is not None and self.el0 is not None:
            if np.abs(self.az0) < self.tol and np.abs(self.el0) < self.tol:
                phi, theta = 0, 0

            else:
                phi = np.degrees(np.arctan(-self.az0/self.el0))
                theta = -self.el0/np.cos(np.radians(phi))
                phi -= 90

            gstring += '{:s}.csy coor_sys_grasp_angles  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}2.csy ),\n'.format(self.name)
            gstring += '  theta            : {:.5f},\n'.format(theta)
            gstring += '  phi            : {:.5f},\n'.format(phi)
            gstring += '  psi            : {:.5f}\n'.format(90)
            gstring += ')\n'
            gstring += ' \n'
            gstring += '{:s}2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'

            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'

        else:

            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'


        if invert or self.invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)

        gstring += ')\n'

        return gstring

class Horn2():
    '''
    This is a Gaussian beam pattern
    '''

    def __init__(self, name=None, units='mm', frequency_list=None,
        taper_angle=25, taper=-3, field_accuracy=-80, auto_convergence='on',
        conv_grid=None, rot_matrix=np.identity(3),
        x=0, y=0, z=0, x1=0, y1=0, z1=0, base=None,
        invert=False, invert_x=False, invert_z=False, flip_pol=False):


        self.name = name
        self.poname = name+'.po'
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.taper_angle = taper_angle
        self.taper = taper
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.x = x
        self.y = y
        self.z = z
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.rot_matrix = rot_matrix

        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z

        self.flip_pol = flip_pol

    def grasp_string(self, invert=False):

        gstring = '\n'
        gstring += '{:s}  gaussian_beam_pattern\n'.format(self.poname)
        gstring += '(\n'
        if self.flip_pol:
            gstring += '  coor_sys         : ref({:s}_flip1.csy),\n'.format(self.name)
        else:
            gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  taper_angle      : {:5.5f},\n'.format(self.taper_angle)
        gstring += '  taper            : {:5.5f}\n'.format(self.taper)
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x1, self.units, self.y1, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z1, self.units)

        # gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
        gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
        gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
        gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
        gstring += ')\n'
        gstring += ' \n'

        if self.flip_pol:
            gstring += '{:s}_flip1.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_flip2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_flip1.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'


        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: 0.0 {:s}),\n'.\
            format(self.units)

        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.flip_pol:
            gstring += '  base             : ref({:s}_flip2.csy)\n'.format(self.name)
        else:
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)

        gstring += ')\n'
        gstring += ' \n'

        return gstring

class Horn2_nf():
    '''
    This is a Gaussian beam pattern (near field)
    '''

    def __init__(self, name=None, units='mm', frequency_list=None,
        beam_radius=10., field_accuracy=-80, auto_convergence='on',
        conv_grid=None, rot_matrix=np.identity(3),
        x=0, y=0, z=0, x1=0, y1=0, z1=0, base=None,
        invert=False, invert_x=False, invert_z=False, flip_pol=False):

        self.name = name
        self.poname = name+'.po'
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.beam_radius = beam_radius
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.x = x
        self.y = y
        self.z = z
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.rot_matrix = rot_matrix

        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z

        self.flip_pol = flip_pol

    def grasp_string(self, invert=False):

        gstring = '\n'
        gstring += '{:s}  gaussian_beam\n'.format(self.poname)
        gstring += '(\n'
        if self.flip_pol:
            gstring += '  coor_sys         : ref({:s}_flip1.csy),\n'.format(self.name)
        else:
            gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  beam_radius      : {:5.5f} {:s},\n'.format(self.beam_radius, self.units)
        gstring += '  phase_front_radius  :  0.0 m\n'
        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x1, self.units, self.y1, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z1, self.units)

        # gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
        gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
        gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
        gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
        gstring += ')\n'
        gstring += ' \n'

        if self.flip_pol:
            gstring += '{:s}_flip1.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_flip2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_flip1.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'


        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: 0.0 {:s}),\n'.\
            format(self.units)

        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.flip_pol:
            gstring += '  base             : ref({:s}_flip2.csy)\n'.format(self.name)
        else:
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)

        gstring += ')\n'
        gstring += ' \n'

        return gstring

class EllipticalPattern():

    def __init__(self, name=None, frequency_list=None, taper=1.0,
        taper_anglex='', taper_angley='', polarisation_angle=90.0, 
        far_forced='on',
        units='mm',
        field_accuracy=-80, auto_convergence='on',
        conv_grid=None, rot_matrix=np.identity(3),
        x=0, y=0, z=0, x1=0, y1=0, z1=0,
        base=None, invert=False,
        invert_x=False, invert_z=False,
        flip_pol=False,
        power_norm=False, az0=None, el0=None, tol=1e-6,
        factor=None, use_grasp_angles=False,
        theta=None, phi=None, psi=None
        ):

        self.name = name
        self.poname = name+'.po'
        self.frequency_list = frequency_list
        self.taper = taper
        self.taper_anglex = taper_anglex
        self.taper_angley = taper_angley
        self.polarisation_angle = polarisation_angle
        self.far_forced = far_forced

        self.base = base
        self.units = units
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        
        self.x = x
        self.y = y
        self.z = z
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.rot_matrix = rot_matrix

        self.az0 = az0
        self.el0 = el0
        self.tol = tol

        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z

        self.theta = theta
        self.phi = phi
        self.psi = psi
        self.use_grasp_angles = use_grasp_angles

        self.flip_pol = flip_pol
        self.power_norm = power_norm

        self.factor = factor

    def grasp_string(self, invert=False):

        gstring = '\n'
        gstring += '{:s}  elliptical_pattern\n'.format(self.poname)
        gstring += '(\n'

        gstring += '  frequency           : ref({:s}),\n'.format(self.frequency_list)
        #gstring += '  coor_sys            : ref({:s}.csy),\n'.format(self.name)        
        
        if isinstance(self.taper, str):
            gstring += '  taper               : "ref({:s})",\n'.format(self.taper)
        else:    
            gstring += '  taper               : {:5.8f},\n'.format(self.taper)

        gstring += '  taper_angles         : struct(zx: "ref({:s})", zy: "ref({:s})"),\n'.\
            format(self.taper_anglex, self.taper_angley)
        gstring += '  polarisation_angle  : {:5.8f},\n'.format(self.polarisation_angle)
        gstring += '  far_forced          : {:s},\n'.format(self.far_forced)        

        ############################################
        #### Coordinate system related, taken from TabulatedPattern2
        if self.flip_pol and (self.az0 is None or self.el0 is None):
            gstring += '  coor_sys         : ref({:s}_flip1.csy),\n'.format(self.name)
        else:
            gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)

        if self.power_norm:
            gstring += '  power_norm       : on,\n'

        if self.factor is not None:
            gstring += '  factor           : struct(db: {:.5f}, deg: 0.0),\n'.format(self.factor)

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x1, self.units, self.y1, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z1, self.units)

        # if self.flip_pol:
        gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
        gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # else:
        # gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
        gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
        gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
        gstring += ')\n'
        gstring += ' \n'

        if self.flip_pol:
            gstring += '{:s}_flip1.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_flip2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_flip1.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

        if (self.az0 is not None and self.el0 is not None) or self.use_grasp_angles:

            if (np.abs(self.az0) < self.tol and np.abs(self.el0) < self.tol) and not self.use_grasp_angles:
                phi, theta = 0, 0
            elif (np.abs(self.el0) < self.tol) and not self.use_grasp_angles:
                phi, theta = 0, self.az0

            else:
                #print(self.tol)
                #print(self.el0)
                phi = np.degrees(np.arctan(-self.az0/self.el0)) if self.phi is None else self.phi
                theta = -self.el0/np.cos(np.radians(phi)) if self.theta is None else self.theta
                phi -= 90

            gstring += '{:s}.csy coor_sys_grasp_angles  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}2.csy ),\n'.format(self.name)
            gstring += '  theta            : {:.5f},\n'.format(theta)
            gstring += '  phi            : {:.5f},\n'.format(phi)
            gstring += '  psi            : {:.5f}\n'.\
                format((0. if self.flip_pol else 90.) if self.psi is None else self.psi)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'

        else:

            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: 0.0 {:s}),\n'.\
            format(self.units)

        if invert or self.invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            # gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.flip_pol:
            gstring += '  base             : ref({:s}_flip2.csy)\n'.format(self.name)
        else:
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        return gstring


class RealVariable():
    def __init__(self, name=None, value=0.0):

        self.name = name
        self.value = value

    def grasp_string(self):
        gstring = '\n'
        gstring += '{:s}  real_variable\n'.format(self.name)
        gstring += '(\n'
        gstring += '  value         : {:.7f}\n'.format(self.value)
        gstring += ')'

        return gstring

class TabulatedPattern():

    def __init__(self, name=None, filename=None, units='mm',
        frequency_list=None, field_accuracy=-80, auto_convergence='on',
        conv_grid=None,
        x=0, y=0, z=0, base=None, invert=False, invert_x=False, invert_z=False):

        self.name = name
        self.poname = name+'.po'
        self.filename = filename
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.x = x
        self.y = y
        self.z = z

        self.az0 = az0
        self.el0 = el0

        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z

    def grasp_string(self, invert=False):

        gstring = '\n'
        gstring += '{:s}  tabulated_pattern\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  file_name        : {:s}\n'.format(self.filename)
        gstring += ')\n'
        gstring += ' \n'
        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'


        # if self.az0 is not None and self.el0 is not None:
        #     if np.abs(self.az0) < self.tol and np.abs(self.el0) < self.tol:
        #         phi, theta = 0, 0

        #     else:
        #         phi = np.degrees(np.arctan(-self.az0/self.el0))
        #         theta = -self.el0/np.cos(np.radians(phi))
        #         phi -= 90

        #     gstring += '{:s}_grasp_angles.csy coor_sys_grasp_angles  \n'.format(self.name)
        #     gstring += '(\n'
        #     gstring += '  base             : ref({:s}.csy ),\n'.format(self.name)
        #     gstring += '  theta            : {:.5f},\n'.format(theta)
        #     gstring += '  phi            : {:.5f},\n'.format(phi)
        #     gstring += '  psi            : {:.5f}\n'.format(90)
        #     gstring += ')\n'
        #     gstring += ' \n'

        if invert or self.invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)


        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)


        gstring += ')\n'

        return gstring

class TabulatedPattern2():

    def __init__(self, name=None, filename=None, units='mm',
        frequency_list=None, field_accuracy=-80, auto_convergence='on',
        conv_grid=None, rot_matrix=np.identity(3),
        x=0, y=0, z=0, x1=0, y1=0, z1=0, number_of_cuts=4,
        base=None, invert=False,
        invert_x=False, invert_z=False,
        flip_pol=False,
        power_norm=False, az0=None, el0=None, tol=1e-6,
        factor=None, use_grasp_angles=False,
        theta=None, phi=None, psi=None):

        self.name = name
        self.poname = name+'.po'
        self.filename = filename
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.field_accuracy = field_accuracy
        self.auto_convergence = auto_convergence
        self.conv_grid = conv_grid
        self.number_of_cuts = number_of_cuts
        self.x = x
        self.y = y
        self.z = z
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.rot_matrix = rot_matrix

        self.az0 = az0
        self.el0 = el0
        self.tol = tol

        self.invert = invert
        self.invert_x = invert_x
        self.invert_z = invert_z

        self.theta = theta
        self.phi = phi
        self.psi = psi
        self.use_grasp_angles = use_grasp_angles

        self.flip_pol = flip_pol
        self.power_norm = power_norm

        self.factor = factor

    def grasp_string(self, invert=False):

        gstring = '\n'

        gstring += '{:s}  tabulated_pattern\n'.format(self.poname)
        gstring += '(\n'

        if self.flip_pol and (self.az0 is None or self.el0 is None):
            gstring += '  coor_sys         : ref({:s}_flip1.csy),\n'.format(self.name)
        else:
            gstring += '  coor_sys         : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency        : ref({:s}),\n'.format(self.frequency_list)
        gstring += '  number_of_cuts   : {:d},\n'.format(self.number_of_cuts)
        if self.power_norm:
            gstring += '  power_norm       : on,\n'

        if self.factor is not None:
            gstring += '  factor           : struct(db: {:.5f}, deg: 0.0),\n'.format(self.factor)

        gstring += '  file_name        : {:s}\n'.format(self.filename)

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_trans.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x1, self.units, self.y1, self.units)+'z: {:5.5f} {:s}),\n'.\
            format(self.z1, self.units)

        # if self.flip_pol:
        gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
        gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # else:
        # gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        # gstring += '  y_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'

        gstring += ')\n'
        gstring += ' \n'

        gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'
        gstring += '  base             : ref({:s}_trans.csy),\n'.format(self.name)
        gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
        gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
            format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
        gstring += ')\n'
        gstring += ' \n'

        if self.flip_pol:
            gstring += '{:s}_flip1.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: -1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: 1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}_flip2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  base             : ref({:s}_flip1.csy)\n'.format(self.name)
            gstring += ')\n'
            gstring += ' \n'

        if (self.az0 is not None and self.el0 is not None) or self.use_grasp_angles:

            if (np.abs(self.az0) < self.tol and np.abs(self.el0) < self.tol) and not self.use_grasp_angles:
                phi, theta = 0, 0
            elif (np.abs(self.el0) < self.tol) and not self.use_grasp_angles:
                phi, theta = 0, self.az0

            else:
                #print(self.tol)
                #print(self.el0)
                phi = np.degrees(np.arctan(-self.az0/self.el0)) if self.phi is None else self.phi
                theta = -self.el0/np.cos(np.radians(phi)) if self.theta is None else self.theta
                phi -= 90

            gstring += '{:s}.csy coor_sys_grasp_angles  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}2.csy ),\n'.format(self.name)
            gstring += '  theta            : {:.5f},\n'.format(theta)
            gstring += '  phi            : {:.5f},\n'.format(phi)
            gstring += '  psi            : {:.5f}\n'.\
                format((0. if self.flip_pol else 90.) if self.psi is None else self.psi)
            gstring += ')\n'
            gstring += ' \n'

            gstring += '{:s}2.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'

        else:

            gstring += '{:s}.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: 0.0 {:s}),\n'.\
            format(self.units)

        if invert or self.invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            # gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'
        if self.invert_x and not invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
        if self.invert_z and not invert:
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if self.flip_pol:
            gstring += '  base             : ref({:s}_flip2.csy)\n'.format(self.name)
        else:
            gstring += '  base             : ref({:s}_rot.csy)\n'.format(self.name)
        gstring += ')\n'
        gstring += ' \n'

        return gstring

class RectangularHorn():

    def __init__(self, name=None, units='mm', frequency_list=None,
        aperture_width=6, aperture_height=6,
        flare_length_xz=12, flare_length_yz=12, field_accuracy=-80,
        x=0, y=0, z=0, factor=None, phase=None, base=None):

        self.name = name
        self.poname = name+'.po'
        self.base = base
        self.units = units
        self.frequency_list = frequency_list
        self.aperture_width = aperture_width
        self.aperture_height = aperture_height
        self.flare_length_xz = flare_length_xz
        self.flare_length_yz = flare_length_yz
        self.field_accuracy = field_accuracy
        self.x = x
        self.y = y
        self.z = z
        self.factor = factor
        self.phase = phase

    def grasp_string(self, invert=False, rotate=False):

        gstring = '\n'
        gstring += '{:s}  rectangular_horn\n'.format(self.poname)
        gstring += '(\n'
        gstring += '  coor_sys          : ref({:s}.csy),\n'.format(self.name)
        gstring += '  frequency         : ref({:s}),\n'.format(self.frequency_list)
        if self.factor is not None or self.phase is not None:
            # factor           : struct(db: 10.0, deg: 0.0)
            phase2u = 0.0 if self.phase is None else self.phase
            factor2u = 0.0 if self.factor is None else self.factor
            gstring += '  factor           : struct(db: {:.5f}, deg: {:.3f}),\n'.\
                format(factor2u, phase2u)

        gstring += '  aperture_width    : {:5.5f} {:s},\n'.format(self.aperture_width, self.units)
        gstring += '  aperture_height   : {:5.5f} {:s},\n'.format(self.aperture_height, self.units)
        gstring += '  flare_length_xz    : {:5.5f} {:s},\n'.format(self.flare_length_xz, self.units)
        gstring += '  flare_length_yz   : {:5.5f} {:s}\n'.format(self.flare_length_yz, self.units)
        gstring += ')\n'
        gstring += ' \n'
        gstring += '{:s}.csy coor_sys  \n'.format(self.name)
        gstring += '(\n'

        if invert:
            gstring += '  x_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'
            gstring += '  z_axis           : struct(x: 0.0, y: 0.0, z: -1.0),\n'

        if rotate:
            gstring += '  x_axis           : struct(x: 0.0, y: 1.0, z: 0.0),\n'
            gstring += '  y_axis           : struct(x: -1.0, y: 0.0, z: 0.0),\n'


        if self.base is not None:
            gstring += '  base             : ref({:s}),\n'.format(self.base)

        gstring += '  origin           : struct(x: {:5.5f} {:s}, y: {:5.5f} {:s}, '.format(
            self.x, self.units, self.y, self.units)+'z: {:5.5f} {:s})\n'.\
            format(self.z, self.units)

        gstring += ')\n'

        return gstring

class VoltageGenerator():

    def __init__(self, name=None, fancy_name=None, base=None,
        coord=np.array([0., 0., 0.]),
        amplitude=1.0, phase=0., units='mm'):

        self.name = name
        self.poname = name+'.po'
        self.fancy_name = fancy_name
        self.base = base
        self.coord = coord
        self.units = units
        self.amplitude = amplitude
        self.phase = phase

    def grasp_string(self):

        gstring = '{:s}  voltage_generator\n'.format(self.poname)
        gstring += '(\n'

        if self.base is not None:
            gstring += '  coor_sys             : ref({:s}),\n'.format(self.base)
        else:
            gstring += '  coor_sys             : ref({:s}.csy),\n'.format(self.name)

        gstring += '  generators       : sequence\n    (    '
        gstring += 'struct(x: {:.4e} {:s}, y: {:.4e} {:s}, z: {:.4e} {:s}'.\
            format(self.coord[0], self.units,
                   self.coord[1],  self.units,
                   self.coord[2], self.units)
        gstring += ',\n         amplitude: {:.4e} V, \n         phase: {:.4e})\n'.\
            format(self.amplitude, self.phase)
        gstring += '    )\n'
        gstring += ')\n\n'

        return gstring

class PiecewiseStraightWire():
    '''
    Creates a Piecewise Straight Wire object
    '''

    def __init__(self, name=None, fancy_name=None, base=None,
        radius=0.1, nodes=None, length=1., units='mm'):

        self.name = name
        self.fancy_name = fancy_name
        self.base = base
        self.nodes = nodes
        self.length = length
        self.radius = radius
        self.units = units

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}.msh  piecewise_straight_wire\n'.format(self.name)
        gstring += '(\n'

        if self.base is not None:
            gstring += '  coor_sys             : ref({:s}),\n'.format(self.base)
        else:
            gstring += '  coor_sys             : ref({:s}.csy),\n'.format(self.name)

        ## Writing nodes to file
        if self.nodes is None:
            self.nodes = np.array([[0., 0., 0.], [0., 0., self.length]] )

        (nx, ny) = np.shape(self.nodes)
        if ny != 3:
            raise ValueError('The nodes array has to be a N x 3 matrix')

        nstr = '  nodes           : sequence(\n'
        for i, node in enumerate(self.nodes):

            if i == 0:
                nstr += '      struct(x: {:.5e} {:s}, y: {:.5e} {:s}, z: {:.5e} {:s}),\n'.\
                 format(node[0], self.units, node[1], self.units, node[2], self.units, )

            elif i == (nx-1):
                nstr += '      struct(x: {:.5e} {:s}, y: {:.5e} {:s}, z: {:.5e} {:s})\n'.\
                 format(node[0], self.units, node[1], self.units, node[2], self.units, )

            else:
                nstr += '      struct(x: {:.5e} {:s}, y: {:.5e} {:s}, z: {:.5e} {:s}),\n'.\
                 format(node[0], self.units, node[1], self.units, node[2], self.units, )

        nstr += '      ),\n'


        gstring += nstr
        gstring += '   radius             : {:.5f} {:s}\n'.format(self.radius, self.units)

        gstring += ')\n'
        gstring += ' \n'

        return gstring

class SphericalGrid():

    def __init__(self, name=None, fancy_name=None, units='m', az1=-1, az2=1, el1=-1, el2=1, Np=101,
            ref_coor='base', source_name=None, sources=[], ref_obj=None, near_field=False, near_dist=10,
            gtype='elevation_and_azimuth'):

        self.name = name
        self.fancy_name = fancy_name
        self.units = units
        self.ref_coor = ref_coor
        self.gtype = gtype
        self.az1 = az1
        self.az2 = az2
        self.el1 = el1
        self.el2 = el2
        self.Np = Np
        self.ref_obj = ref_obj
        self.near_field = near_field
        self.near_dist = near_dist
        self.source_name = source_name
        self.sources = sources

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}  spherical_grid\n'.format(self.name)
        gstring += '(\n'
        if self.ref_coor is not None:
            gstring += '  coor_sys         : ref({:s}),\n'.format(self.ref_coor)

        if self.near_field:
            gstring += '  near_far         : near,\n'
            gstring += '  near_dist        : {:.5f} {:s},\n'.\
                format(self.near_dist, self.units)

        if self.gtype != 'uv':
            gstring += '  grid_type        : {},\n'.format(self.gtype)

        gstring += '  x_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.az1, self.az2, self.Np)
        gstring += '  y_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.el1, self.el2, self.Np)
        gstring += ')'

        return gstring

class SphericalCut():

    def __init__(self, name=None, units='m', th1=-90, th2=90, phi1=0, phi2=180,
        ntheta=101, nphi=5, ref_coor='base', near_field=False, near_dist=10,
        rot_matrix=None, az0=None, el0=None, tol=1e-6):

        self.name = name
        self.units = units
        self.ref_coor = ref_coor
        self.th1 = th1
        self.th2 = th2
        self.phi1 = phi1
        self.phi2 = phi2
        self.ntheta = ntheta
        self.nphi = nphi
        self.near_field = near_field
        self.near_dist = near_dist
        self.rot_matrix = rot_matrix
        self.az0 = az0
        self.el0 = el0
        self.tol = tol

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}  spherical_cut\n'.format(self.name)
        gstring += '(\n'

        gstring += '  coor_sys         : ref({:s}_rot.csy),\n'.format(self.name)
        if self.near_field:
            gstring += '  near_far         : near,\n'
            gstring += '  near_dist        : {:.5f} {:s},\n'.\
                format(self.near_dist, self.units)

        gstring += '  theta_range      : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.th1, self.th2, self.ntheta)
        gstring += '  phi_range        : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.phi1, self.phi2, self.nphi)
        gstring += ')\n\n'

        if self.rot_matrix is not None:

            gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}),\n'.format(self.ref_coor)
            gstring += '  x_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
                format(self.rot_matrix[0][0], self.rot_matrix[0][1], self.rot_matrix[0][2])
            gstring += '  y_axis           : struct(x: {:.5e}, y: {:.5e}, z: {:.5e}),\n'.\
                format(self.rot_matrix[1][0], self.rot_matrix[1][1], self.rot_matrix[1][2])
            gstring += ')\n'
            gstring += ' \n'

        elif self.az0 is not None and self.el0 is not None:
            if np.abs(self.az0) < self.tol and np.abs(self.el0) < self.tol:

                print('In this rare exception')
                phi, theta = 0, 0

            else:

                # phi = np.degrees(np.arctan(-self.el0/self.az0))
                # theta = -self.az0/np.cos(np.radians(phi))

                if np.abs(self.el0) < self.tol:
                    phi = np.sign(self.az0) * 90
                else:
                    phi = np.degrees(np.arctan(-self.az0/self.el0))

                theta = -self.el0/np.cos(np.radians(phi))

                phi -= 90

                phi   *= -1
                theta *= -1

            gstring += '{:s}_rot.csy coor_sys_grasp_angles  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}),\n'.format(self.ref_coor)
            gstring += '  theta            : {:.5f},\n'.format(theta)
            gstring += '  phi            : {:.5f},\n'.format(phi)
            gstring += '  psi            : {:.5f}\n'.format(90)
            gstring += ')\n'
            gstring += ' \n'

        else:
            gstring += '{:s}_rot.csy coor_sys  \n'.format(self.name)
            gstring += '(\n'
            gstring += '  base             : ref({:s}),\n'.format(self.ref_coor)
            gstring += ')\n'
            gstring += ' \n'

        return gstring

class PlanarGrid():

    def __init__(self, name=None, fancy_name=None,
        xmin=-1, xmax=1, ymin=-1, ymax=1,
        Np=101, units='m', ref_coor='base', near_dist=0.2, ref_obj=None,
        d=None, sources=[]):

        self.name = name
        self.fancy_name = fancy_name
        self.ref_coor = ref_coor
        self.xmin= xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.units = units
        self.near_dist = near_dist

        self.Np = Np
        self.sources = sources

        if ref_obj is not None:
            ref_units = ref_obj.units
            ref_d = ref_obj.d

            if ref_units not in ['m', 'cm', 'mm']:
                raise ValueError('Only know how to parse m, cm, and mm')

            if ref_units == 'mm':
                self.d = ref_d/1000.
            elif ref_units == 'cm':
                self.d = ref_d/100.
            else:
                self.d = ref_d

        elif d is not None:
            self.d = d

    def grasp_string(self):

        gstring = '\n'

        print('near_dist: {}'.format(self.near_dist))

        gstring += '{:s}  planar_grid\n'.format(self.name)
        gstring += '(\n'
        gstring += '  coor_sys         : ref({:s}),\n'.format(self.ref_coor)
        gstring += '  near_dist        : {:.5f} {:s},\n'.format(self.near_dist, self.units)
        gstring += '  x_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}, unit: {:s}),\n'.\
            format(self.xmin, self.xmax, self.Np, self.units)
        gstring += '  y_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.ymin, self.ymax, self.Np)
        gstring += ')'

        return gstring

class CylindricalGrid():

    def __init__(self, name=None, fancy_name=None,
        r=1.0, phi_min=0.0, phi_max=360.0, zmin=0.0, zmax=1.0,
        nz=101, nphi=101, ref_coor='base',): #, elevation=None

        self.name = name
        self.fancy_name = fancy_name
        self.ref_coor = ref_coor
        self.r = r
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.zmin = zmin
        self.zmax = zmax
        self.nphi = nphi
        self.nz = nz
        self.units = 'm'

    def grasp_string(self):

        gstring = '\n'
        # if self.elevation is not None:

        #     gstring += '{:s}_cylrot.csy coor_sys_grasp_angles\n'.format(self.name)
        #     gstring += '(\n'
        #     gstring += '  theta         : {:5.5f}, \n'.format(self.elevation)
        #     # gstring += '  phi           : {:5.5f}, \n'.format(self.phi)
        #     # gstring += '  psi           : {:5.5f},'.format(self.psi)
        #     gstring += '  base          : ref({:s})\n'.format(self.ref_coor)
        #     gstring += ')\n'

        gstring += '{:s}  cylindrical_grid\n'.format(self.name)
        gstring += '(\n'
        if self.ref_coor is not None:
            gstring += '  coor_sys         : ref({:s}),\n'.format(self.ref_coor)
            # gstring += '  coor_sys         : ref({:s}),\n'.\
            #     format(self.ref_coor if self.elevation is None else '{:s}_cylrot.csy'.format(self.name))

        gstring += '  radius        : {:.5f} {:s},\n'.format(self.r, self.units)
        gstring += '  phi_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.phi_min, self.phi_max, self.nphi)
        gstring += '  z_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.zmin, self.zmax, self.nz)
        gstring += ')'

        return gstring

class PlanarCut():

    def __init__(self, name=None, fancy_name=None,
        rmin=0.0, rmax=1.0, nr=101, phi_min=-180.0, phi_max=180.0, nphi=101,
        z=1.0, ref_coor='base', cut_type='planar'):

        self.name = name
        self.fancy_name = fancy_name
        self.ref_coor = ref_coor
        self.cut_type = cut_type
        self.rmin = rmin
        self.rmax = rmax
        self.nr = nr
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.nphi = nphi
        self.z = z
        self.units = 'm'

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}  planar_cut\n'.format(self.name)
        gstring += '(\n'
        if self.ref_coor is not None:
            gstring += '  coor_sys         : ref({:s}),\n'.format(self.ref_coor)

        if self.cut_type != 'planar':
            gstring += '  cut_type      : {:s},\n'.format(self.cut_type)
        gstring += '  near_dist     : {:.5f} {:s},\n'.format(self.z, self.units)
        gstring += '  phi_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.phi_min, self.phi_max, self.nphi)
        gstring += '  rho_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.rmin, self.rmax, self.nr)
        gstring += ')\n'

        return gstring

class CylindricalCut():

    def __init__(self, name=None, fancy_name=None,
        r=1.0, phi_min=-180.0, phi_max=180.0, zmin=0.0, zmax=1.0,
        nz=101, nphi=101, ref_coor='base', cut_type='circular'): #, elevation=None

        self.name = name
        self.fancy_name = fancy_name
        self.ref_coor = ref_coor
        self.r = r
        self.phi_min = phi_min
        self.phi_max = phi_max
        self.zmin = zmin
        self.zmax = zmax
        self.nphi = nphi
        self.nz = nz
        self.units = 'm'

    def grasp_string(self):

        gstring = '\n'
        gstring += '{:s}  cylindrical_cut\n'.format(self.name)
        gstring += '(\n'
        if self.ref_coor is not None:
            gstring += '  coor_sys         : ref({:s}),\n'.format(self.ref_coor)

        gstring += '  cut_type      : circular,\n'
        gstring += '  radius        : {:.5f} {:s},\n'.format(self.r, self.units)
        gstring += '  phi_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d}),\n'.\
            format(self.phi_min, self.phi_max, self.nphi)
        gstring += '  z_range          : struct(start: {:5.5f}, end: {:5.5f}, np: {:d})\n'.\
            format(self.zmin, self.zmax, self.nz)
        gstring += ')'

        return gstring

class Graspy():
    '''
    '''

    def __init__(self, filename=None, frequency=[10], units=None,
        frequency_list='freq'):

        self.filename = filename
        self.units = units
        self.frequency = frequency # [GHz]
        self.frequency_list = frequency_list


    def create_frequency_list(self):

        str_out = '{:s}  frequency\n'.format(self.frequency_list)
        str_out += '(\n'

        if len(self.frequency) == 1:
            str_out += '  frequency_list   : sequence({:.5f} GHz)\n'.\
                format(self.frequency[0])
            str_out += ')'

        else:

            freqstr = ''.join(['{:.5f} GHz, '.format(nu) for nu in self.frequency[:-1]])
            freqstr += '{:.5f} GHz'.format(self.frequency[-1])
            str_out += '  frequency_list   : sequence({:s})\n'.format(freqstr)
            str_out += ')'

        return str_out

if __name__ == '__main__'        :

    pass
