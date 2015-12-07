#!/usr/bin/env python
# encoding: utf-8

import nefis
import sys
from os.path import exists

# Local import
from nefis2netcdf.nefdata import NefData
from nefis2netcdf.saveasnetcdf import save_as_netcdf


def convert_2_netcdf(fileheader, netcdffilename, exceptions=[], compression=False, debug=False):
    """
    Convert Delft output in Nefis format (*.dat, *.def) into a netcdf file (*.nc)

    :param fileheader: Delft output file header (header_ex.dat, header_ex.def), string
                       Note that the header can also include a path (ex.: path_2_file/header_ex.dat...)
    :param netcdffilename: Netcdf file name, string
                           Note that the header can also include a path (ex.: path_2_file/netcdffilename)
    :param exceptions: list of variables to exclude from output file , list of strings
    :param compression: compresses data with zlib and uses at least 3 significant digits, boolean
    :param debug: debug flag, boolean
    """
    # Check if netcdffilename has an extension
    if not netcdffilename[-3:] == '.nc':
        netcdffilename = netcdffilename + '.nc'
    # Check if file already exist
    if exists(netcdffilename):
        input1 = raw_input(netcdffilename+" already exists. Overwrite? (y/n): ")
        if input1 in ['n', 'N', 'no', 'No', 'NO']:
            netcdffilename = raw_input("New netcdf file name: ")

    ndata = NefData(fileheader, debug=debug)
    save_as_netcdf(ndata, netcdffilename, exceptions=exceptions, compression=compression, debug=debug)