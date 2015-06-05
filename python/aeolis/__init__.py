'''AeoLiS : sediment transport model Aeolian transport with Limited Supply

Usage:
    aeolis run <path>

Positional arguments:
    path               model run location

Options:
    -h, --help         show this help message and exit
'''

__all__ = ['filesys', 'plot']

import filesys

import os
import re
import sys
import time
import glob
import netCDF4
import itertools
import numpy as np
import pandas as pd
from datetime import datetime
from bmi.wrapper import BMIWrapper

class AeoLiS:
    '''AeoLiS model class'''

    t = 0
    it = 0
    iout = 0
    tstop = None
    output = None
    output_stats = None
    ios = 0

    ncattrs = dict()
    
    params = dict(
        mixtoplayer   = True,
        sweeptoplayer = True,
        th_grainsize  = True,
        th_bedslope   = True,
        th_moisture   = True,
        th_humidity   = True,
        bedupdate     = True,
        evaporation   = True,
        gusts         = True,
    
        VS    = 1.0,
        Tp    = 1.0,
        u_th  = 4.0,
        z0    = 1.0,
        k     = 0.01,
        Cb    = 1.8,
        phi   = 40.0,
        g     = 9.81,
        dt    = 0.05,
        dx    = 1.0,
        dy    = 1.0,
        tstop = 3600.0,
        CFL   = -1.0,
        accfac = 1.0,

        method_moist = 'belly_johnson',
        scheme = 'euler_backward',

        max_iter = 100,
        max_error = 1e-6,

        nfractions = 1,
        nlayers = 3,
        layer_thickness = 0.0005,
        minfrac = 0.0001,
        rhop = 2650.0,
        rhom = 1650.0,
        rhow = 1025.0,
        rhoa = 1.25,
        porosity = 0.4,
        A = 100.0,
        Hs = 1.0,
        gamma = 0.5,
        facDOD = 0.1,
        F = 1e-4,
        Cw = 0.0,
        w = 0.03,
        bi = 1.0,

        grain_size = [300e-6],
        grain_dist = [1.],

        wind_file = '',
        bed_file = '',
        tide_file = '',
    )

    def __init__(self, fpath, outputtimes=None, outputvars=['Ct', 'u', 'zb'], outputtypes=[],
                 configfile='aeolis.txt', outputfile='aeolis.nc', **kwargs):
        self.fpath = fpath
    
        self.set_outputtimes(outputtimes)
        self.set_outputvars(outputvars)
        self.set_outputtypes(outputtypes)

        self.set_configfile(configfile)
        self.set_outputfile(outputfile)

        self.set_params(**kwargs)
        

    def __repr__(self):
        s =  'AeoLiS_Input instance:\n'
        s += '  %10s:  %s\n' % ('configfile', self.configfile)
        s += '  %10s:  %s\n' % ('outputfile', self.outputfile)
        for i, (k, v) in enumerate(sorted(self.params.iteritems())):
            if i == 0:
                s += '  %10s: %15s = %s\n' % ('params', k, v)
            else:
                if self.__isiterable(v):
                    s += '  %10s  %15s = %s\n' % ('', k, v[0])
                    for vi in v[1:]:
                        s += '  %10s  %15s = %s\n' % ('', '', vi)
                else:
                    s += '  %10s  %15s = %s\n' % ('', k, v)
        return s


    def set_outputtimes(self, outputtimes):
        self.outputtimes = outputtimes


    def set_outputvars(self, outputvars):
        self.outputvars = outputvars


    def set_outputtypes(self, outputtypes):
        self.outputtypes = outputtypes


    def set_configfile(self, configfile):
        self.configfile = os.path.abspath(os.path.join(self.fpath, configfile))

        
    def set_outputfile(self, outputfile):
        self.outputfile = os.path.abspath(os.path.join(self.fpath, outputfile))

        
    def set_params(self, **kwargs):
        for k, v in kwargs.iteritems():
            self.params[k] = v


    def set_ncattrs(self, **kwargs):
        for k, v in kwargs.iteritems():
            self.ncattrs[k] = v

            
    def write_params(self, standalone=False):

        print 'Writing input to %s' % self.configfile

        params = self.params
        
        # set tout to large to prevent slow computations
        if not params.has_key('tout'):
            params['tout'] = params['tstop']

        if standalone:
            params['output_dir'] = self.outputfile.replace('.nc','')
            params['outputvars'] = self.outputvars
            params['outputtypes'] = self.outputtypes
        
        with open(self.configfile, 'w') as fp:
            for k, v in sorted(params.iteritems()):
                if isinstance(v, bool):
                    fp.write('%s = %s\n' % (k, 'T' if v else 'F'))
                elif isinstance(v, int):
                    fp.write('%s = %d\n' % (k, v))
                elif isinstance(v, float):
                    fp.write('%s = %f\n' % (k, v))
                elif isinstance(v, str):
                    fp.write('%s = %s\n' % (k, v))
                elif self.__isiterable(v):
                    w = [''] * len(v)
                    for i in range(len(v)):
                        if isinstance(v[i], int):
                            w[i] = '%d' % v[i]
                        elif isinstance(v[i], float):
                            w[i] = '%f' % v[i]
                        elif isinstance(v[i], str):
                            w[i] = v[i]
                        else:
                            raise TypeError('Unknown parameter type [%s, %d]' % (k, i))
                    fp.write('%s = %s\n' % (k, ' '.join(w)))
                else:
                    raise TypeError('Unknown parameter type [%s]' % k)


    def init_output(self, model):

        ny, nx, nl, nf = model.get_var_shape('mass')
        
        with netCDF4.Dataset(self.outputfile, 'w') as nc:

            ## add dimensions
            nc.createDimension('x', nx)
            nc.createDimension('y', ny)
            nc.createDimension('time', 0)
            nc.createDimension('nv', 2)
            nc.createDimension('nv2', 4)
            nc.createDimension('layers', nl)
            nc.createDimension('fractions', nf)
            
            ## add global attributes
            # see http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/formats/DataDiscoveryAttConvention.html
            nc.Conventions = 'CF-1.6'
            nc.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
            nc.featureType = 'grid'
            nc.cdm_data_type = 'grid'
            nc.standard_name_vocabulary = 'CF-1.6'
            nc.title = ''
            nc.summary = ''
            nc.source = 'AeoLiS'
            nc.id = ''
            nc.naming_authority = ''
            nc.time_coverage_start = ''
            nc.time_coverage_end = ''
            nc.time_coverage_resolution = ''
            nc.geospatial_lat_min = 0
            nc.geospatial_lat_max = 0
            nc.geospatial_lat_units = 'degrees_north'
            nc.geospatial_lat_resolution = ''
            nc.geospatial_lon_min = 0
            nc.geospatial_lon_max = 0
            nc.geospatial_lon_units = 'degrees_east'
            nc.geospatial_lon_resolution = ''
            nc.geospatial_vertical_min = 0
            nc.geospatial_vertical_max = 0
            nc.geospatial_vertical_units = ''
            nc.geospatial_vertical_resolution = ''
            nc.geospatial_vertical_positive = ''
            nc.institution = ''
            nc.creator_name = ''
            nc.creator_url = ''
            nc.creator_email = ''
            nc.project = ''
            nc.processing_level = ''
            nc.references = ''
            nc.keywords_vocabulary = 'NASA/GCMD Earth Science Keywords. Version 6.0'
            nc.keywords = ''
            nc.acknowledgment = ''
            nc.comment = ''
            nc.contributor_name = ''
            nc.contributor_role = ''
            nc.date_created = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
            nc.date_modified = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
            nc.date_issued = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
            nc.publisher_name = ''
            nc.publisher_email = ''
            nc.publisher_url = ''
            nc.history = ''
            nc.license = ''
            nc.metadata_link = '0'
            
            ## add variables
            nc.createVariable('x', 'float32', (u'x',))
            nc.variables['x'].long_name = 'x-coordinate'
            nc.variables['x'].standard_name = 'projection_x_coordinate'
            nc.variables['x'].units = 'm'
            nc.variables['x'].axis = 'X'
            nc.variables['x'].valid_min = 0
            nc.variables['x'].valid_max = 0
            nc.variables['x'].bounds = 'x_bounds'
            nc.variables['x'].grid_mapping = 'crs'
            nc.variables['x'].comment = ''
            
            nc.createVariable('y', 'float32', (u'y',))
            nc.variables['y'].long_name = 'y-coordinate'
            nc.variables['y'].standard_name = 'projection_y_coordinate'
            nc.variables['y'].units = 'm'
            nc.variables['y'].axis = 'Y'
            nc.variables['y'].valid_min = 0
            nc.variables['y'].valid_max = 0
            nc.variables['y'].bounds = 'y_bounds'
            nc.variables['y'].grid_mapping = 'crs'
            nc.variables['y'].comment = ''

            nc.createVariable('layers', 'float32', (u'layers',))
            nc.variables['layers'].long_name = 'bed layers'
            nc.variables['layers'].standard_name = ''
            nc.variables['layers'].units = '-'
            nc.variables['layers'].valid_min = 0
            nc.variables['layers'].valid_max = 0
            nc.variables['layers'].comment = ''

            nc.createVariable('fractions', 'float32', (u'fractions',))
            nc.variables['fractions'].long_name = 'sediment fractions'
            nc.variables['fractions'].standard_name = ''
            nc.variables['fractions'].units = 'm'
            nc.variables['fractions'].valid_min = 0
            nc.variables['fractions'].valid_max = 0
            nc.variables['fractions'].comment = ''
            
            nc.createVariable('lat', 'float32', (u'y', u'x'))
            nc.variables['lat'].long_name = 'latitude'
            nc.variables['lat'].standard_name = 'latitude'
            nc.variables['lat'].units = 'degrees_north'
            nc.variables['lat'].valid_min = 0
            nc.variables['lat'].valid_max = 0
            nc.variables['lat'].bounds = 'lat_bounds'
            nc.variables['lat'].ancillary_variables = ''
            nc.variables['lat'].comment = ''
            
            nc.createVariable('lon', 'float32', (u'y', u'x'))
            nc.variables['lon'].long_name = 'longitude'
            nc.variables['lon'].standard_name = 'longitude'
            nc.variables['lon'].units = 'degrees_east'
            nc.variables['lon'].valid_min = 0
            nc.variables['lon'].valid_max = 0
            nc.variables['lon'].bounds = 'lon_bounds'
            nc.variables['lon'].ancillary_variables = ''
            nc.variables['lon'].comment = ''
            
            nc.createVariable('time', 'float64', (u'time',))
            nc.variables['time'].long_name = 'time'
            nc.variables['time'].standard_name = 'time'
            nc.variables['time'].units = 'seconds since 1970-01-01 00:00:00 0:00'
            nc.variables['time'].calendar = 'julian'
            nc.variables['time'].axis = 'T'
            nc.variables['time'].bounds = 'time_bounds'
            nc.variables['time'].ancillary_variables = ''
            nc.variables['time'].comment = ''
            
            nc.createVariable('x_bounds', 'float32', (u'x', u'nv'))
            nc.variables['x_bounds'].units = 'm'
            nc.variables['x_bounds'].comment = 'x-coordinate values at the upper and lower bounds of each pixel.'
            
            nc.createVariable('y_bounds', 'float32', (u'y', u'nv'))
            nc.variables['y_bounds'].units = 'm'
            nc.variables['y_bounds'].comment = 'y-coordinate values at the left and right bounds of each pixel.'
            
            nc.createVariable('lat_bounds', 'float32', (u'y', u'x', u'nv2'))
            nc.variables['lat_bounds'].units = 'degrees_north'
            nc.variables['lat_bounds'].comment = 'latitude values at the north and south bounds of each pixel.'
            
            nc.createVariable('lon_bounds', 'float32', (u'y', u'x', u'nv2'))
            nc.variables['lon_bounds'].units = 'degrees_east'
            nc.variables['lon_bounds'].comment = 'longitude values at the west and east bounds of each pixel.'
            
            nc.createVariable('time_bounds', 'float32', (u'time', u'nv'))
            nc.variables['time_bounds'].units = 'seconds since 1970-01-01 00:00:00 0:00'
            nc.variables['time_bounds'].comment = 'time bounds for each time value'

            for var0 in self.outputvars:

                dims = AeoLiS.get_dims(var0)
                if dims is None:
                    continue

                for typ in [None] + self.outputtypes:

                    if typ:
                        var = '%s.%s' % (var0, typ)
                    else:
                        var = var0
                        
                    nc.createVariable(var, 'float32', dims)
                    nc.variables[var].long_name = var
                    nc.variables[var].standard_name = ''
                    nc.variables[var].units = ''
                    nc.variables[var].scale_factor = 1.0
                    nc.variables[var].add_offset = 0.0
                    nc.variables[var].valid_min = 0
                    nc.variables[var].valid_max = 0
                    nc.variables[var].coordinates = ' '.join(dims)
                    nc.variables[var].grid_mapping = 'crs'
                    nc.variables[var].source = ''
                    nc.variables[var].references = ''
                    nc.variables[var].cell_methods = ''
                    nc.variables[var].ancillary_variables = ''
                    nc.variables[var].comment = ''
            
            nc.createVariable('crs', 'int32', ())
            nc.variables['crs'].grid_mapping_name = 'oblique_stereographic'
            nc.variables['crs'].epsg_code = 'EPSG:28992'
            nc.variables['crs'].semi_major_axis = 6377397.155
            nc.variables['crs'].semi_minor_axis = 6356078.96282
            nc.variables['crs'].inverse_flattening = 299.1528128
            nc.variables['crs'].latitude_of_projection_origin = 52.0922178
            nc.variables['crs'].longitude_of_projection_origin = 5.23155
            nc.variables['crs'].scale_factor_at_projection_origin = 0.9999079
            nc.variables['crs'].false_easting = 155000.0
            nc.variables['crs'].false_northing = 463000.0
            nc.variables['crs'].proj4_params = '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +towgs84=565.4174,50.3319,465.5542,-0.398957388243134,0.343987817378283,-1.87740163998045,4.0725 +no_defs'

            # set netcdf attributes
            for k, v in self.ncattrs.iteritems():
                if isinstance(v, bool):
                    nc.setncattr(k, int(v))
                else:
                    nc.setncattr(k, v)
                
            # store static data
            nc.variables['x'][:] = np.arange(0, nx, 1) * self.params['dx']
            nc.variables['y'][:] = np.arange(0, ny, 1) * self.params['dy']
            nc.variables['layers'][:] = np.arange(0, nl, 1) * self.params['layer_thickness']
            nc.variables['fractions'][:] = self.params['grain_size']
            
            # store model settings
            grp = nc.createGroup('settings')
            for k, v in self.params.iteritems():
                if isinstance(v, bool):
                    grp.setncattr(k, int(v))
                else:
                    grp.setncattr(k, v)

        # initialize output stats
        self.output_stats = {}
        for var in self.outputvars:
            self.output_stats[var] = {}
            for typ in self.outputtypes:
                self.output_stats[var][typ] = None
        self.clear_output(model)
                
        # initialize output property
        self.output = AeoLiS_Output(self.outputfile)

            
    def update_output(self, model):

        if len(self.outputtypes) > 0:
            for var in self.outputvars:
                data = model.get_var(var)
                if 'sum' in self.outputtypes:
                    self.output_stats[var]['sum'] += data
                if 'var' in self.outputtypes:
                    self.output_stats[var]['var'] += data**2
                if 'min' in self.outputtypes:
                    self.output_stats[var]['min'] = np.min(self.output_stats[var]['min'], data)
                if 'max' in self.outputtypes:
                    self.output_stats[var]['max'] = np.max(self.output_stats[var]['max'], data)
                
        self.ios += 1

                
    def clear_output(self, model):

        n = {}
        n['y'], n['x'], n['layers'], n['fractions'] = model.get_var_shape('mass')

        for var, types in self.output_stats.iteritems():
            dims = AeoLiS.get_dims(var)
            if dims is None:
                continue
            if 'sum' in self.outputtypes:
                self.output_stats[var]['sum'] = np.zeros([n[k] for k in dims[1:]])
            if 'var' in self.outputtypes:
                self.output_stats[var]['var'] = np.zeros([n[k] for k in dims[1:]])
            if 'min' in self.outputtypes:
                self.output_stats[var]['min'] = np.zeros([n[k] for k in dims[1:]]) + np.inf
            if 'max' in self.outputtypes:
                self.output_stats[var]['max'] = np.zeros([n[k] for k in dims[1:]]) - np.inf

        self.ios = 0

                
    def write_output(self, model):

        i = self.iout
        with netCDF4.Dataset(self.outputfile, 'a') as nc:
            nc.variables['time'][i] = self.t
            for var in self.outputvars:
                data = model.get_var(var)
                nc.variables[var][i,...] = data
                if 'sum' in self.outputtypes:
                    nc.variables['%s.sum' % var][i,...] = self.output_stats[var]['sum']
                if 'avg' in self.outputtypes:
                    nc.variables['%s.avg' % var][i,...] = self.output_stats[var]['sum'] / self.ios
                if 'var' in self.outputtypes:
                    nc.variables['%s.var' % var][i,...] = (self.output_stats[var]['var'] - \
                                                           self.output_stats[var]['sum']**2 / \
                                                           self.ios) / (self.ios - 1)
                if 'min' in self.outputtypes:
                    nc.variables['%s.min' % var][i,...] = self.output_stats[var]['min']
                if 'max' in self.outputtypes:
                    nc.variables['%s.max' % var][i,...] = self.output_stats[var]['max']

        self.clear_output(model)
        self.iout += 1


    def run(self, overwrite=True):

        if os.path.exists(self.outputfile) and not overwrite:
            return

        print ' '
        print '         d8888                   888      d8b  .d8888b.   ' 
        print '        d88888                   888      Y8P d88P  Y88b  ' 
        print '       d88P888                   888          Y88b.       ' 
        print '      d88P 888  .d88b.   .d88b.  888      888  "Y888b.    ' 
        print '     d88P  888 d8P  Y8b d88""88b 888      888     "Y88b.  ' 
        print '    d88P   888 88888888 888  888 888      888       "888  ' 
        print '   d8888888888 Y8b.     Y88..88P 888      888 Y88b  d88P  ' 
        print '  d88P     888  "Y8888   "Y88P"  88888888 888  "Y8888P"   '
        print ' '

        self.write_params()

        print ' '
            
        with BMIWrapper(engine='aeolis', configfile=self.configfile) as model:

            self.init_output(model)
        
            t = 0
            self.t = 0
            self.it = 0
            self.iout = 0
            self.tstop = model.get_end_time()
            
            tstart = time.time()
            tlog = tstart

            while self.t <= self.tstop+1:
                t = self.t

                model.update(-1)

                self.update_output(model)
                
                self.t = model.get_current_time()
                self.it += 1
                        
                if self.outputtimes is not None:
                    if np.mod(self.t, self.outputtimes) < self.t - t:
                        self.write_output(model)

                if (np.mod(self.t, self.tstop/10.) < self.t - t or \
                   time.time() - tlog > 60.):
                    
                    p = self.t / self.tstop
                    dt1 = time.time() - tstart
                    dt2 = dt1 / p
                    dt3 = dt2 * (1-p)

                    fmt = '[%5.1f%%] %s / %s / %s (avg. dt=%5.3f)'
                    
                    print fmt % (p*100.,
                                 time.strftime('%H:%M:%S', time.gmtime(dt1)),
                                 time.strftime('%H:%M:%S', time.gmtime(dt2)),
                                 time.strftime('%H:%M:%S', time.gmtime(dt3)),
                                 self.t / self.it)

                    tlog = time.time()


    @staticmethod
    def get_dims(var):

        var = var.split('.')[0]
        
        if var in ['mass']:
            dims = (u'time', u'y', u'x', u'layers', u'fractions')
        elif var in ['d10', 'd50', 'd90', 'moist', 'thlyr']:
            dims = (u'time', u'y', u'x', u'layers')
        elif var in ['Cu', 'Ct', 'uth', 'supply', 'p']:
            dims = (u'time', u'y', u'x', u'fractions')
        elif var in ['x', 'z', 'zb']:
            dims = (u'time', u'y', u'x')
        elif var in ['uw']:
            dims = (u'time',)
        else:
            dims = None
        return dims
    

    def __isiterable(self, x):
        if isinstance(x, str):
            return False
        try:
            [xi for xi in x]
        except TypeError:
            return False
        return True

                    
class AeoLiS_Output:
    '''AeoLiS output reader
    
    Reads one or more AeoLiS output directories. Can be initialized
    using a glob structure. It reads multiple directories into a
    single MultiIndex pandas structure.
    
    Examples
    --------
    >>> a = AeoLiS('cases/case_T1.0_Cw*')
    >>> a.variables['Ct'][:,-1,:]

    '''
    
    fpath = []
    fpaths = []
    
    def __init__(self, fpath):
        '''Initialize paths
        
        Parameters
        ----------
        fpath : str or list
            String or list with strings with paths of glob structures

        '''
        
        if type(fpath) is str:
            self.fpath = [fpath]
        else:
            self.fpath = fpath
            
        self.variables = AeoLiS_Variables(self)
        
        self.read_paths()
    
    
    def __repr__(self):
        s =  'AeoLiS_Output instance:\n'
        s += '  paths: %s\n' % self.fpath[0]
        for fpath in self.fpath[1:]:
            s += '         %s\n' % fpath
        s += '\n'
        for fpath in self.fpaths:
            s += '         %s\n' % fpath
        return s
        
        
    def read_paths(self):
        '''Create list of paths from glob'''
        
        self.fpaths = []
        for fpath in self.fpath:
            for fname in glob.glob(fpath):
                self.fpaths.append(fname)

            
    def iterate_paths(self, absolute_paths=True):
        '''Iterate all paths in structure
        
        Parameters
        ----------
        absolute_paths : bool, optional
            Flag to toggle absolute path output

        '''
        
        for fpath in self.fpaths:
            if os.path.exists(fpath):
                if absolute_paths:
                    yield fpath
                else:
                    yield os.path.split(fpath)[1]
            
            
    def parse_path(self, fpath):
        '''Split path into multiple parts that represent the model settings
        
        If paths are constructed such that they identify the model
        settings used, this function splits the path into these
        identification parts. Each part consists of a non-numeric
        keyword directly followed by a numeric value. Parts are
        separated by a non-alphanumeric charachter, like an underscore
        or dash.
        
        Parameters
        ----------
        fpath : str
            Path to be splitted
            
        Returns
        -------
        list of tuples
            List of 2-tuples with keyword/value pairs
            
        Examples
        --------
        >>> parse_path('aeolis_Cw0_T1.0_bi3.0')
            [('Cw', 0.0), ('T', 1.0), ('bi', 3.0)]

        '''
        
        parts = []
        fpath = os.path.split(fpath)
        if len(fpath[1]) > 0:
            fpath = fpath[1]
        else:
            fpath = fpath[0]
        for fpart in re.split('[_-]+', fpath):
            m = re.match('^([^\d]+)([\d\.]+)(\.nc|\/)?$', fpart)
            if m:
                parts.append((m.groups()[0], float(m.groups()[1])))
        if len(parts) == 0:
            parts = [fpath]
        return parts


    def get_settings(self, diff=True):
        ncattrs = {}
        
        for i, fpath in enumerate(self.fpaths):
            if fpath.endswith('.nc'):
                with netCDF4.Dataset(fpath, 'r') as nc:
                    for attr in nc.groups['settings'].ncattrs():
                        if not ncattrs.has_key(attr):
                            ncattrs[attr] = [[]] * len(self.fpaths)
                        ncattrs[attr][i] = nc.groups['settings'].getncattr(attr)
                        
        if diff and len(self.fpaths) > 1:
            for attr in ncattrs.keys():
                n = 0
                for i, val1 in enumerate(ncattrs[attr]):
                    n += 1
                    for val2 in ncattrs[attr][:i]:
                        chk = val1 == val2
                        if not isinstance(chk, bool):
                            chk = chk.all()
                        if chk:
                            n -= 1
                            break
                if n == 1:
                    del ncattrs[attr]
                            
        return ncattrs
        
            
class AeoLiS_Variables:
    '''AeoLiS variables reader
    
    Implements a __getitem__ function that first selects a variable,
    reads it into a Pandas DataFrame and indexes the DataFrame. In
    case multiple output directories are read, these are concatenated
    into a MultiIndex DataFrame. Any singular levels in the MultiIndex
    DataFrame are removed.

    '''
    
    varname = None
    vartype = None
    
    def __init__(self, parent):
        self.parent = parent
    
    
    def __repr__(self):
        s =  'AeoLiS_Variables instance:\n'
        s += '  variable name: %s\n' % self.varname
        s += '  variable type: %s\n' % self.vartype
        return s
    
    
    def __getitem__(self, s):
        '''Implements a __getitem__ function for variable names and indices
        
        Parameters
        ----------
        s : str, slice or tuple of slices
            String to select variable or slices for variable index if
            variable is already selected
            
        Returns
        -------
        self or pandas.DataFrame
            self if type(s) is str and a Pandas DataFrame otherwise
        
        Example
        -------
        >>> a = AeoLiS('cases/case_T1.0_Cw*')
        >>> a.variables
            AeoLiS_Variables instance:
                variable: None
        >>> a.variables['Ct']
            AeoLiS_Variables instance:
                variable: Ct
        >>> a.variables['Ct'][:,-1,:]

        '''
        
        if type(s) == str:
            s = s.split('.')
            if len(s) > 0:
                self.varname = s[0]
            if len(s) > 1:
                self.vartype = s[1]
            return self
        else:
            dfs = []
            for fpath in self.parent.iterate_paths():
                try:
                    df = self.read_variable(fpath, s)
                    p = zip(*self.parent.parse_path(fpath))
                    
                    if len(p) < 2:
                        p.append(None)
                        
                    if type(df) is pd.Panel or type(df) is pd.Panel4D:
                        df = self.to_dataframe(df)

                    df.index = pd.MultiIndex.from_tuples([p[1] + (ix,)
                                                          for ix in df.index], names=p[0] + ('',))

                    dfs.append(df)
                except:
                    raise
                    print 'Skipped %s' % fpath
                    
            df = pd.concat(dfs)
                
            return self.consolidate(self.consolidate(df).T).T
        
        
    def read_variable(self, fpath, s, variable=None):
        '''Reads a variable from file or interprets as alias if not found

        Parameters
        ----------
        fpath : str
            Path to look for variable file
        variable : str, optional
            Variable to look for, if not given self.varname is used

        Returns
        -------
        pandas.DataFrame
            Pandas DataFrame with data from variable file

        '''
        
        if variable is None:
            var = self.varname
        else:
            var = variable

        # check if netcdf
        if fpath.endswith('.nc'):
            var = self.__filename(var, include_ext=False)
            return self.read_variable_netcdf(fpath, s, variable=var)

        # check if fortran
        fname = os.path.join(fpath, self.__filename(var))
        if os.path.exists(fname):
            return filesys.load_dataframe(fname).iloc[s]

        # check if alias
        if variable is None:
            return self.read_alias(fpath, s)

        raise ValueError('Variable not found: %s' % variable)


    def read_variable_netcdf(self, fpath, s, variable=None):

        if not isinstance(s, tuple):
            s = (s,)

        if variable is None:
            var = self.varname
        else:
            var = variable

        df = None
        with netCDF4.Dataset(fpath, 'r') as nc:

            # read axes
            dims = AeoLiS.get_dims(var)
            if dims is None:
                NotImplemented('No support for variable %s' % var)
            axs = []
            for i, dim in enumerate(dims):
                if dim == 'time':
                    ax = netCDF4.num2date(nc.variables['time'][s[i]],
                                          nc.variables['time'].units)
                else:
                    ax = nc.variables[dim][s[i]]
                if isinstance(ax, np.ndarray):
                    axs.append(pd.Index(ax, name=dim))

            # read data
            data = nc.variables[var][s]

            # construct pandas object
            if len(axs) == 1:
                df = pd.Series(data, index=axs[0])
            elif len(axs) == 2:
                df = pd.DataFrame(data, index=axs[0], columns=axs[1])
            elif len(axs) == 3:
                df = pd.Panel(data, items=axs[0], major_axis=axs[1], minor_axis=axs[2])
            elif len(axs) == 4:
                df = pd.Panel4D(data, labels=axs[0], items=axs[1], major_axis=axs[2], minor_axis=axs[3])
            else:
                raise NotImplemented('No pandas structure with more than four dimensions, reduce dimensionality')

        return df
    
    
    def read_alias(self, fpath, s, variable=None):
        '''Interprets variable as alias

        Currently supported aliases are:
        * saturation: Ct / Cu
        * q: Ct * u

        Alias can be added by overloading this method as follows:
        
            class AeoLiS_Variables_ext(AeoLiS_Variables):
                def read_alias(self, fpath, variable=None):
                    try:
                        df = super(AeoLiS_Variables_ext, self).read_alias(fpath, variable=variable)
                    except:
                        if self.varname == '<alias name>':
                            <implement alias here>
                        else:
                            raise ValueError('Unknown alias [%s]' % self.varname)
                    return df

        Parameters
        ----------
        fpath : str
            Path to look for variable file(s)
        variable : str, optional
            Alias to construct, if not given self.varname is used

        Returns
        -------
        pandas.DataFrame
            Composite Pandas DataFrame with data according to alias

        '''
        
        if self.varname == 'saturation':
            Ct = self.read_variable(fpath, s, 'Ct')
            Cu = self.read_variable(fpath, s, 'Cu')
            df = Ct.divide(np.maximum(1e-10, Cu))
        elif self.varname == 'q':
            Ct = self.read_variable(fpath, s, 'Ct')
            u = self.read_variable(fpath, s, 'u')
            df = Ct.multiply(u, axis=0)
        else:
            raise ValueError('Unknown alias [%s]' % self.varname)
        return df
    
        
    def consolidate(self, df):
        '''Consolidates Pandas DataFrame
        
        Consolidates Pandas MultiIndex DataFrame by removing all
        MultiIndex levels with a length equal or smaller than 1
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to be consolidates
            
        Returns
        -------
        pandas.DataFrame
            Consolidated DataFrame

        '''
        
        try:
            for i in range(len(df.index.levels),0,-1):
                if len(df.index.levels[i-1]) <= 1:
                    df.reset_index(level=i-1, drop=True, inplace=True)
        except:
            pass
        return df
    

    @staticmethod
    def to_dataframe(obj):
        '''Convert any higher dimensional pandas object to a MultiIndex DataFrame

        Convert any higher dimensional pandas object to a MultiIndex
        DataFrame by keeping the index of the first dimension and
        stacking all others into a MultiIndex to be used as columns.

        Currently supported objects are:
        * Panel
        * Panel4D
        * PanelND

        Parameters
        ----------
        obj : pandas object
            Pandas object to be converted

        Returns
        -------
        pandas.DataFrame
            Converted pandas DataFrame
        '''

        if hasattr(obj, 'axis_orders'):
            axs = obj.axis_orders
        elif type(obj) is pd.Panel4D:
            axs = ['labels', 'items', 'major_axis', 'minor_axis']
        elif type(obj) is pd.Panel:
            axs = ['items', 'major_axis', 'minor_axis']
        else:
            raise NotImplemented('Unsupported object type')

        dims = []
        names = []
        levels = []
        for ax in axs:
            ax = getattr(obj, ax)
            dims.append(range(len(ax)))
            levels.append(ax.values)
            names.append(ax.name)

        ix1 = pd.Index(data=levels[0], name=names[0])
        ix2 = pd.MultiIndex(levels=levels[1:],
                            labels=zip(*[r for r in itertools.product(*dims[1:])]),
                            names=names[1:])

        data = obj.as_matrix().reshape((len(dims[0]),-1))
        df = pd.DataFrame(data, index=ix1, columns=ix2)

        return df

    
    def __filename(self, varname=None, vartype=None, include_ext=True):

        '''Constructs filename'''
        
        fname = ''
        if varname is None:
            varname = self.varname
        if vartype is None:
            vartype = self.vartype
        if varname is not None:
            fname += varname
            if vartype is not None:
                fname += '.%s' % vartype
            if include_ext:
                fname += '.out'
        return fname

    
def cmd():
    import docopt

    arguments = docopt.docopt(__doc__)

    a = AeoLiS(arguments['<path>'])
    a.run()


if __name__ == '__main__':
    cmd()
