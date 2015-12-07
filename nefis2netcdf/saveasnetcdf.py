#!/usr/bin/env python
# encoding: utf-8

#Libs import
from __future__ import division
import numpy as np
import netCDF4 as nc

def save_as_netcdf(nefdata, filename, exceptions=[], compression=False, debug=False):
    """
    NefData object in a netcdf file

    inputs:
      - nefdata = NefData object
      - filename = file name, string
    options:
      - exceptions = list of variables to exclude from output file
                     , list of strings
      - compression = compresses data with zlib and uses at least 3 significant digits, boolean
        Note: Works only with netcdf format
    """
    if debug: print "Computing bounding box..."
    if compression:
        zlib = True
        least_significant_digit = 3
    else:
        zlib = False
        least_significant_digit = None

    # Check if netcdffilename has an extension
    if not filename[-3:] == '.nc':
        filename = filename + '.nc'
    f = nc.Dataset(filename, 'w') #, format='NETCDF4_CLASSIC') # not sure why it doesn't work

    # Create dimensions dictionnary
    dimlist = []
    dims = {}

    # first spatial dimension
    spacelist = ['NMAX', 'KMAX', 'MMAX']
    spacedims = []
    elmlist = nefdata.elmlist

    for d in spacelist:
            try:
                dim = int(nefdata.load_data(d))
                dims[d[0].lower()+"dim"] = dim
                spacedims.append(dim)
            except KeyError:
                pass

    # then other dims including time dimension
    dims['tdim'] = 0
    for elm in elmlist:
        des, type, bytesize, quantity, units, dimensions, ntimes, size = nefdata.elm_info(elm, print_out=False)
        dimlist += list(dimensions)
        if ntimes > dims['tdim']:
            dims['tdim'] = ntimes

    # clean up dimensions
    spacedims.append(dims['tdim'])
    dimlist = list(np.unique(dimlist))
    for k in spacedims:
        dimlist.remove(k)

    for ii, dim in enumerate(dimlist):
        dims['dim'+str(ii)] = dim
    if debug: print dims
    for key in dims.keys():
        f.createDimension(key, dims[key])

    # Load in netcdf file
    if debug: print "Loading in nc file..."
    for var in nefdata.elmlist:
        if not var in exceptions:  # check if in list of exceptions
            # Load data
            data = nefdata.load_data(var)
            # defines data's dimensions
            dim = []
            try:  # in case of string
                s = data.shape
            except AttributeError:
                s = [1]
            for d in s:
                flag = 1
                count = 0
                while flag:
                    if count == len(dims.keys()):  # when two dimensions are the same by coincidence
                        # dim.append(dim[-1])  # TODO search in the entire list == d rather than last index
                        for k in dim:
                            if dims[k] == d:
                                newkey = k
                        dim.append(newkey)
                        flag = 0
                    else:
                        key = dims.keys()[count]
                        if dims[key] == d:
                            if len(dim) == len(s):  # in case of similar dims
                                count += 1
                                pass
                            else:
                                if key not in dim:  # make sure dimension doesn't get recorded twice
                                    dim.append(key)
                                    flag = 0
                                    count += 1
                                else:
                                    count += 1
                        else:
                            count += 1
            dim = tuple(dim)
            if debug: print "Dim: "+str(dim)
            # Create variable and metadata
            #  data type
            des, format, bytesize, quantity, units, dimensions, ntimes, size = nefdata.elm_info(var, print_out=False)
            if format.lower().replace(' ', '') == 'real':
                format = float
                # format = 'f8'
            elif format.lower().replace(' ', '') == 'integer':
                format = int
                # format = 'i8'
            else:
                format = str
                # format = 'c'
            if debug:
                print "Format: " + str(format)
                print "Dim: " + str(dim)
                print "... saving "+var+"..."

            tmp_var = f.createVariable(var.lower(), format, dim,
                                       zlib=zlib, least_significant_digit=least_significant_digit)
            #  metadata
            tmp_var.long_name = des
            tmp_var.units = units.replace(' ', '').replace('[', '').replace(']', '').lower()
            #  load in data
            try:  # in case of string and VLEN slices
                tmp_var[:] = data[:]
            except IndexError:  # todo: come up with proper fix
                tmp_var[0] = data

    # Clean exit
    f.close()