#!/usr/bin/env python
# encoding: utf-8
import numpy as np
import struct
import os
import sys
import nefis


class NefData(object):
    def __init__(self, fileheader, ac_type='r', coding=' ', debug=False):
        """
        NefData class is a python wrapper for the Nefis c-lib.
        It permits to read or write Nefis type files.

        :param fileheader: Nefis file name header, string
        :param ac_type: access type, character.
                        'w’ create a NEFIS file set, existing NEFIS file will be deleted
                        ’r’ read only access of the NEFIS file set
                        ’r+’ or 'w+' update of the NEFIS file set, write and read access of the NEFIS file set.
        :param coding: Coding type, character.
                       ’B’ big-endian representation of the data
                       ’L’ little-endian representation of the data
                       ' ' Endianess-representation is taken from the machine
        :param debug: debug flag, boolean
        """
        # ---Attributes definition---
        # groups
        self.grplist = []
        self.grp = {}
        # elements
        self.elmlist=[]
        self.elm={}
        # file names
        self.dat_file = fileheader+'.dat'
        self.def_file = fileheader+'.def'
        # access type
        if ac_type == 'w':
            self.ac_type = 'c'
        elif ac_type in ['r+', 'w+']:
            self.ac_type = 'u'
        else:
            self.ac_type = 'r'
        # Coding aka endianess-representation
        if coding in ['L', 'l']:
            self.coding = 'L'
        elif coding in ['B', 'b']:
            self.coding = 'B'
        else:
            self.coding = 'automatic'
        # debug flag
        if debug:
            self.debug=True
        else:
            self.debug=False

        # read/write/update delft output files
        self.fp = -1
        error, self.fp = nefis.crenef(self.dat_file, self.def_file, coding, ac_type)
        self._print_error(error)

        # populate a list of all groups and elements
        self._grp_list()
        self._elm_list()

        # debug message
        if self.debug:
            print("---------")
            print(self.dat_file)
            print(self.def_file)
            print(self.coding)
            print(self.ac_type)
            print("---------")

    # Special methods
    def __del__(self):
        """
        Close the data and definition file. Depending on the access type the output buffers will be
        written to the files before closing the files
        """
        error = nefis.clsnef(self.fp)
        self._print_error(error)

    def _print_error(self, error):
        """
        Prints error message.

        :param error: error code, int
        """
        if not error == 0:
            print("---------")
            print nefis.neferr()[1]
            print("---------")
            sys.exit(1)

    def _grp_list(self):
        """
        Reads groups from definition file
        """
        ngrp = 0
        while True:
            if ngrp == 0:
                error, grp_name, cel_name, grp_count_dimensions, grp_dimensions, grp_order = nefis.inqfgr(self.fp)
            else:
                error, grp_name, cel_name, grp_count_dimensions, grp_dimensions, grp_order = nefis.inqngr(self.fp)

            if error == -6028: # no more groups to read
                return self.grplist
            elif not error == 0:  # some other error
                self._print_error(error)
                break
            else:
                ngrp += 1
            self.grplist.append(grp_name)
            # Give the maximum used index of the free dimension of a data group.
            error, ntimes = nefis.inqmxi(self.fp, grp_name)
            self._print_error(error)
            self.grp[grp_name]={}
            self.grp[grp_name]['ntimes'] = ntimes

    def _elm_list(self):
        """
        Read elements from the definition file.
        """
        nelm = 0
        while True:
            if nelm == 0:
                error, elm_name, elm_type, elm_quantity, elm_unit, elm_description, elm_single_bytes, elm_bytes, elm_count, elm_dimensions = nefis.inqfel(self.fp)
            else:
                error, elm_name, elm_type, elm_quantity, elm_unit, elm_description, elm_single_bytes, elm_bytes, elm_count, elm_dimensions = nefis.inqnel(self.fp)

            if error == -6024: # no more elements to read
                return self.elmlist
            elif not error == 0:  # some other error
                self._print_error(error)
                break
            else:
                nelm+=1
            self.elmlist.append(elm_name)

            # define dictionaries for each element name
            self.elm[elm_name]={}
            self.elm[elm_name]['info']={}
            info = self.elm[elm_name]['info']
            info['type'] = elm_type
            info['unit'] = elm_unit
            info['desc'] = elm_description
            info['bytesize'] = elm_single_bytes
            info['bytes'] = elm_bytes
            info['ndim'] = elm_count
            info['dims'] = elm_dimensions[::-1]  # reverse dimensions (C vs Fortran ??)

            elm_grp = self._find_elm(elm_name)[0]
            info['group'] = elm_grp

            # Give the maximum used index of the free dimension of a data group
            error, ntimes = nefis.inqmxi(self.fp, elm_grp)
            self._print_error(error)
            info['ntimes'] = ntimes

    def _find_elm(self, search_elm_name):
        """
        Finds group name and associated element's number

        :param search_elm_name: element's name, string
        :return: group name (string) and element's number (int)
        """
        ngrp=0
        while True:
            if ngrp == 0:
                error, grp_name, cel_name, grp_count_dimensions, grp_dimensions, grp_order = nefis.inqfgr(self.fp)
            else:
                error, grp_name, cel_name, grp_count_dimensions, grp_dimensions, grp_order = nefis.inqngr(self.fp)

            ngrp += 1
            if error == -6028: # no more groups to read
                return '', -1   # element not found

            # Read a cell definition from the definition file
            error, count, elm_names = nefis.inqcel(self.fp, cel_name)
            self._print_error(error)

            for i in range(count):
                if search_elm_name.strip() == elm_names[17*i:17*(i+1)].strip():
                    return grp_name, i    # return group and element number

    # Instant methods
    def grp_info(self, grp_name, print_out=True):
        """
        Print out and return information on any given group

        :param grp_name: group name, string
        :param print_out: print-out flag, boolean
        :return: count dimensions, dimensions, order, cell's name
        """
        error, cell_name, grp_count_dimensions, grp_dimensions, grp_order = nefis.inqgrp(self.fp, grp_name)
        self._print_error(error)
        if print_out:
            print('  Group count  : "%d"' % grp_count_dimensions)
            print('  Dimensions   : "%s"' % grp_dimensions)
            print('  Group order  : "%s"' % grp_order)
            print('  Cell name    : "%s"' % cell_name)
            # maximum used index of the free dimension of a data group
            error, ntimes = nefis.inqmxi (self.fp, grp_name)
            self._print_error(error)
            print('  Nb time step : "%d"' % ntimes)

        return grp_count_dimensions, grp_dimensions, grp_order, cell_name

    def cell_info(self, cel_name, print_out=True):
        """
        Read a cell definition from the definition file

        :param cel_name: cell's name, string
        :param print_out: print-out flag, boolean
        :return: element's count, elements' names
        """
        error, count, elm_names = nefis.inqcel(self.fp, cel_name)
        self._print_error(error)
        if print_out:
            print('  Element Count: "%d"' % count)
            print('  Element names: "%s"' % elm_names)

        return count, elm_names

    def elm_info(self, elm_name, print_out=True):
        """
        Print out and return an element definitions from the definition file
        :param elm_name:
        :param print_out: print-out flag, boolean
        :return: description, type, byte size, quantity, units, dimensions, nb. time step, size
        """

        error, elm_type, elm_single_byte, elm_quantity, elm_unit, elm_description, elm_count, elm_dimensions = nefis.inqelm(self.fp, elm_name)
        self._print_error(error)
        elm = self.elm[elm_name]
        info = elm['info']
        ntimes = info['ntimes']

        if print_out:
            print '  Element      : ', elm_name
            print '  Description  : ', elm_description
            print '  Type         : ', elm_type
            print '  Byte size    : ', elm_single_byte
            print '  Quantity     : ', elm_quantity
            print '  Units        : ', elm_unit
            print '  Dimensions   : ', elm_dimensions
            print '  Nb time steps: ', ntimes
            print '  Size         : ', elm_count

        return elm_description, elm_type, elm_single_byte, elm_quantity, elm_unit, elm_dimensions, ntimes, elm_count

    def file_info(self):
        """
        Print out all the file info
        """
        print "====Group Information===="
        for grp in self.grplist:
            print " - "+grp+": "
            out1, out2, out3, out4 = self.grp_info(grp)
            del out1, out2, out3, out4

        print "====Element Information===="
        for elm in self.elmlist:
            print " - "+elm+": "
            out1, out2, out3, out4, out5, out6, out7 = self.elm_info(elm)
            del out1, out2, out3, out4, out5, out6, out7

    def load_data(self, elm_name):
        """
        Load element's data

        :param elm_name: element's name, string
        :return: element's data, numpy array
        """
        if not elm_name in self.elm.keys():
            print 'ERROR: ', elm_name, ' does not exist'
            return None

        elm = self.elm[elm_name]
        info = elm['info']
        grp_name = info['group']

        usr_index = np.zeros(15, dtype=np.int32).reshape(5,3)

        # format inputs for nefis.getelt
        ntimes = info['ntimes']
        usr_index[0,0] = 1  # first timestep
        usr_index[0,1] = ntimes  # last time step
        usr_index[0,2] = 1   # increment
        np.ascontiguousarray(usr_index, dtype=np.int32)

        length = info['bytes']*ntimes

        usr_order = np.arange(1,6,dtype=np.int32).reshape(5)
        np.ascontiguousarray(usr_order, dtype=np.int32)

        error, buffer_res = nefis.getelt(self.fp, grp_name, elm_name, usr_index, usr_order, length)
        self._print_error(error)

        if ntimes > 1:
            dims = np.append(ntimes,info['dims']) # add time as a dimension
        else:
            dims = info['dims']

        if info['type'].strip() == 'REAL':
            if info['bytesize'] == 4:
                fmt = "%df" % (length/4)
            if info['bytesize'] == 8:
                fmt = "%dd" % (length/8)
            elm_data = np.asarray(struct.unpack(fmt, buffer_res)).reshape( dims )  # unpack data into array
        elif info['type'].strip() == 'INTEGER':
            if info['bytesize'] == 4:
                fmt = "%di" % (length/4)
            if info['bytesize'] == 8:
                fmt = "%dl" % (length/8)
            elm_data = np.asarray(struct.unpack(fmt, buffer_res)).reshape( dims )  # unpack data into array
        elif info['type'].strip() == 'CHARACTE':
            elm_data = buffer_res

        ## TODO -- determine if it makes sense to reorganize data into t,i,j,k ordera -- could cause problems
        #data = elm['data']
        #
        #if info['ndim'] == 3:
        #    if ntimes > 1:
        #        # dimensions are (t, k, i, j)
        #        data = np.transpose(data,(0,2,3,1))   # now it's t,i,j,k
        #        info['dims'] = info['dims'][[0,2,3,1]]  # update dims to match
        #    else:
        #        # dimensions are (k, i, j)
        #        data = np.transpose(data,(1,2,0))  # now it's i,j,k
        #        info['dims'] = info['dims'][[1,2,0]]  # update dims to match

        return elm_data
       


