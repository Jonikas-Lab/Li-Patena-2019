#!/usr/bin/env python2
"""
Generate Biomek transfer files to combine N samples into M pools, using a provided binary code.

Given a number of samples, a number of pools, a binary code, and information about the physical plates used for the samples and pools, assign one codeword from the code to each sample.  Write a varying number of Biomek robot transfer files (depends on -m/-o and -M options) that when given the appropriate plates will execute the combinatorial pooling.  Also write a general output file containing all the script information, the plate/well positions and codewords for each sample, and the plate/well positions for each pool. 

Assumes the samples will be in sequential order over the provided number of same-sized input plates: [plate1 well A1, plate1 well A2, ..., plate1 well B2, .., plate1 well A1, ...].  The pools will be similarly distributed over the provided number of output plates.  The number and size of output plates does not need to have any particular relationship to the number and size of input plates.

Arbitrarily assigns codewords from the provided binary code to each sample; the pools correspond to each bit of the codeword.  Whenever a bit of the codeword for the sample is 1, that sample should be added to the corresponding pool, so a line of the form "sample_plate_name,sample_plate_well,pool_plate_name,pool_plate_well,volume" is added to the Biomek command list. 

For more details on inputs and available options, run script with -h, or see help strings in define_option_parser function.

 -- Weronika Patena, Jonikas lab, Carnegie Institution, July 2011

USAGE:  robotic_plate_transfer.py [options] outfile_base_name
        robotic_plate_transfer.py [-h] [-t] [-T]
"""

# standard libraries
import sys, os
import unittest
from collections import defaultdict, Counter
import itertools
from math import ceil
from string import ascii_uppercase, ascii_lowercase     # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' and 'abcdef...'
# my modules
import binary_code_utilities
from general_utilities import invert_list_to_dict, save_line_list_as_file, write_header_data
from testing_utilities import run_functional_tests

class PlateTransferError(Exception):
    """ Exception for this file (does nothing interesting)."""
    pass

class Plate_type:
    """ A plate type (384-, 96-, 24-, 6-well, or custom). Converts between sequential numbers and well positions."""

    ### defined plate types: standard
    _standard_plate_shapes = {6: (2,3), 24: (4,6), 96: (8,12), 384: (16,24)}
    ### defined plate types: custom
    _custom_plate_types = {}
    ## "fake 6-well" plate types: really using a 6-well plate in the Biomek but pretending it's a 96-well (by giving
    #   the 96-well positions that correspond to 6-wells), since the Biomek won't accept a 6-well plate description.
    # (two versions - one with a single deep-96-well per 6-well, one with four to get around volume constraints)
    _custom_plate_types['fake6_simple'] = (6, ['B2','B7','B11','G2','G7','G11'])
    _custom_plate_types['fake6_complex'] = (6, [('B2','B3','C2','C3'), ('B6','B7','C6','C7'), ('B10','B11','C10','C11'), 
                                                ('F2','F3','G2','G3'), ('F6','F7','G6','G7'), ('F10','F11','G10','G11')])

    ### all defined plate types as a name:(size,well_ID_list) dict - standard and custom together, for easier usage
    _plate_type_definitions = dict([(str(size), (size,None)) for size in _standard_plate_shapes])
    # add custom sizes
    _plate_type_definitions.update(_custom_plate_types)
    # making a nice string-only no-duplicates prettily-sorted version of the allowed plate type names, for printing
    defined_plate_types = sorted( sorted( set([str(x) for x in _plate_type_definitions.keys()])), key = lambda x: len(x) )
    del x   # just to get it off the class variable list
    defined_plate_types_str = ', '.join(defined_plate_types)

    def __init__(self, standard_type=None, size=None, well_ID_list=None):
        """ Set up the plate: preexisting standard_type, standard based on size, or custom based on well_ID_list. 

        If standard_type is given, look it up in Plate_type._plate_type_definitions and use the size/well_ID_list 
         from there; if size is also given, check it.  Cannot give well_ID_list together with standard_type.
        If well_ID_list is given, generate a custom plate based on that (and check that its size matches size, if given).
        If only size is given, look it up in Plate_type._standard_plate_shapes and generate a standard plate 
         based on the number of rows and columns from there.
        """
        # MAYBE-TODO could give the plate types names?
        # MAYBE-TODO could make Plate_type class collect a dictionary of all the initialized plate types by name?
        # MAYBE-TODO could give each plate type a max volume
        size = None if (size is None) else int(size)
        # look up standard_type in _plate_type_definitions, set size and well_ID_list based on that, proceed normally.
        if standard_type is not None:
            if well_ID_list is not None:
                raise PlateTransferError("Only one of standard_type and well_ID_list should be specified!")
            expected_size = size
            try:
                size, well_ID_list = self.__class__._plate_type_definitions[str(standard_type)]
            except KeyError:
                raise PlateTransferError("standard_type %s is not defined! Defined types are %s."%(standard_type, 
                                                                                                   defined_plate_types))
        # if well_ID_list is specified, just use it (and check that the size matches, if not None)
        if well_ID_list is not None:    self._make_custom_well_ID_list_from_list(well_ID_list, size)
        # if well_ID_list wasn't passed, infer the plate shape from the size (it must be one of the standard sizes)
        else:                           self._make_standard_well_ID_list_from_size(size)
        # if we used a standard_type and a size was passed, check that the size is as expected
        if standard_type is not None and expected_size is not None:
            if not expected_size == self.size:
                raise PlateTransferError("Standard plate type %s size is %s, not the %s passed by user!"%(standard_type, 
                                                                                              self.size, expected_size))
        self._make_well_generators_and_everything()
    
    def _make_standard_well_ID_list_from_size(self, size):
        """ Generate a well ID list (such that list[number]==ID), based on self.rows/columns, like ['A1','A2', ...]."""
        self.size = size
        # figure out rows/columns based on standard plate shapes
        try:                
           rows, columns = self.__class__._standard_plate_shapes[size]
        except KeyError:    
            raise PlateTransferError("Plate size must be standard (one of %s)"%self.__class__._standard_plate_shapes.keys()
                                     + "when well_ID_list isn't given!  Size %s given is unacceptable."%size)
        assert rows*columns == size, "Malformed plate %s in _standard_plate_shapes!"
        # make the well_ID_list based on rows/columns
        self.full_well_ID_list = []
        for row in range(rows):
            self.full_well_ID_list.extend(['%s%d'%(ascii_uppercase[row],col+1) for col in range(columns)])
        assert len(self.full_well_ID_list) == self.size,\
                "Plate generated by _make_standard_well_ID_list_from_size is the wrong size!"

    def _make_custom_well_ID_list_from_list(self, well_ID_list, size=None):
        """ Set up custom plate: copy well_ID_list to self.full_well_ID_list, check size if given.
        """
        if (size is not None) and size != len(well_ID_list):
            raise PlateTransferError("If both size and well_ID_list are passed when creating new custom Plate_type, "
                                     +"size must match list length! (currently %s and %s)"%(size, well_ID_list))
        self.full_well_ID_list = well_ID_list
        self.size = len(self.full_well_ID_list)

    def _make_well_generators_and_everything(self):
        """ Set up the rest

        We have a full well ID list, the elements of which can be either single wells or tuples of wells.  
        For example, ['A', ('B','C')]
        We'll make two things:
          - a well ID generator list - for the example, it'll be [genA, genBC], where genA.next() is always A, 
            and genBC.next() alternates between B and C
          - a dictionary to get the well number from the ID - for the example, it'll be {A:1, B:2, C:2}
            (note how B and C have the same position, since they were in the same tuple)
        """
        self.well_ID_generator_list = []
        self.well_ID_dict = {}
        for N, wells in enumerate(self.full_well_ID_list):
            if type(wells) != tuple:    wells = (wells,)
            self.well_ID_generator_list.append(itertools.cycle(wells))
            for well in wells:
                if well in self.well_ID_dict:
                    raise PlateTransferError("Well %s shows up multiple times in %s! Each well should only show up once"%(
                                                well, self.full_well_ID_list))
                self.well_ID_dict[well] = N

    def get_well_ID_from_number(self,number):
        """ Given a 0-based well number (4), return the ID (B2 for 6-well plate, A5 for 96-well plate)."""
        try:                return self.well_ID_generator_list[number].next()
        except IndexError:  raise PlateTransferError("Can't get well %s from a %s-well plate!"%(number,self.size))

    # TODO this doesn't actually ever get used!
    def get_well_number_from_ID(self,ID):
        """ Given a well ID (B1), return the sequential 0-based well number (3 for 6-well plate, 12 for 96-well plate)."""
        try:                return self.well_ID_dict[ID]
        except IndexError:  raise PlateTransferError("Can't get well %s from a %s-well plate!"%(ID,self.size))

defined_plate_types = Plate_type.defined_plate_types
defined_plate_types_str = Plate_type.defined_plate_types_str


### Unit-tests for the classes/setup above and "general functions" below
# MAYBE-TODO do I want separate test classes for each function, or should I merge them? What's the standard? Unittest output doesn't give the class name/docstring when it fails, only the function name/docstring...

class Testing__Plate_type(unittest.TestCase):
    """ Unit-tests for the Plate_type class and methods. """

    ### testing the plate types that have already been generated on module import

    def test__get_well_ID_from_number__known_wells(self):
        # try the first, second, first-in-second-row, and last well of each plate size
        assert Plate_type('6').get_well_ID_from_number(0) == 'A1'
        assert Plate_type('6').get_well_ID_from_number(1) == 'A2'
        assert Plate_type('6').get_well_ID_from_number(3) == 'B1'
        assert Plate_type('6').get_well_ID_from_number(5) == 'B3'
        assert Plate_type('24').get_well_ID_from_number(0) == 'A1'
        assert Plate_type('24').get_well_ID_from_number(1) == 'A2'
        assert Plate_type('24').get_well_ID_from_number(6) == 'B1'
        assert Plate_type('24').get_well_ID_from_number(23) == 'D6'
        assert Plate_type(96).get_well_ID_from_number(0) == 'A1'
        assert Plate_type(96).get_well_ID_from_number(1) == 'A2'
        assert Plate_type(96).get_well_ID_from_number(12) == 'B1'
        assert Plate_type(96).get_well_ID_from_number(95) == 'H12'
        assert Plate_type(384).get_well_ID_from_number(0) == 'A1'
        assert Plate_type(384).get_well_ID_from_number(1) == 'A2'
        assert Plate_type(384).get_well_ID_from_number(24) == 'B1'
        assert Plate_type(384).get_well_ID_from_number(383) == 'P24'
        assert Plate_type('fake6_simple').get_well_ID_from_number(0) == 'B2'
        assert Plate_type('fake6_simple').get_well_ID_from_number(1) == 'B7'
        assert Plate_type('fake6_simple').get_well_ID_from_number(3) == 'G2'
        assert Plate_type('fake6_simple').get_well_ID_from_number(5) == 'G11'
        # for the fake6_complex type, the same well number should alternate between four IDs!
        P = Plate_type('fake6_complex')
        assert P.get_well_ID_from_number(0) == 'B2'
        assert P.get_well_ID_from_number(0) == 'B3'
        assert P.get_well_ID_from_number(0) == 'C2'
        assert P.get_well_ID_from_number(0) == 'C3'
        assert P.get_well_ID_from_number(0) == 'B2'
        P = Plate_type('fake6_complex')
        assert P.get_well_ID_from_number(0) == 'B2'
        assert P.get_well_ID_from_number(5) == 'F10'
        assert P.get_well_ID_from_number(5) == 'F11'
        assert P.get_well_ID_from_number(5) == 'G10'
        assert P.get_well_ID_from_number(5) == 'G11'
        assert P.get_well_ID_from_number(5) == 'F10'
        assert P.get_well_ID_from_number(0) == 'B3'
        assert P.get_well_ID_from_number(5) == 'F11'
        assert P.get_well_ID_from_number(0) == 'C2'

    def test__get_well_ID_from_number__last_wells(self):
        # you can also use a negative list index to get wells from the end, why not
        assert Plate_type(6).get_well_ID_from_number(-1) == 'B3'
        assert Plate_type(24).get_well_ID_from_number(-1) == 'D6'
        assert Plate_type('96').get_well_ID_from_number(-1) == 'H12'
        assert Plate_type('384').get_well_ID_from_number(-1) == 'P24'
        assert Plate_type('fake6_simple').get_well_ID_from_number(-1) == 'G11'
        assert Plate_type('fake6_complex').get_well_ID_from_number(-1) == 'F10'

    def test__get_well_ID_from_number__bad_numbers(self):
        # basic tests with obviously wrong values
        self.assertRaises(PlateTransferError, Plate_type('6').get_well_ID_from_number, 100)
        self.assertRaises(TypeError, Plate_type('6').get_well_ID_from_number, 0.5)
        self.assertRaises(TypeError, Plate_type('6').get_well_ID_from_number, 'A')
        self.assertRaises(TypeError, Plate_type('6').get_well_ID_from_number, [1])
        # well numbers are 0-based, so there should be no well N in an N-well plate (the wells are 0..N-1)
        for plate_size in Plate_type._standard_plate_shapes.keys():
            self.assertRaises(PlateTransferError, Plate_type(plate_size).get_well_ID_from_number, plate_size)

    def test__get_well_number_from_ID__known_wells(self):
        # try the first, second, first-in-second-row, and last well of each plate size
        assert Plate_type('6').get_well_number_from_ID('A1') == 0
        assert Plate_type('6').get_well_number_from_ID('A2') == 1
        assert Plate_type('6').get_well_number_from_ID('B1') == 3
        assert Plate_type('6').get_well_number_from_ID('B3') == 5
        assert Plate_type(24).get_well_number_from_ID('A1') == 0
        assert Plate_type(24).get_well_number_from_ID('A2') == 1
        assert Plate_type(24).get_well_number_from_ID('B1') == 6
        assert Plate_type(24).get_well_number_from_ID('D6') == 23
        assert Plate_type(96).get_well_number_from_ID('A1') == 0
        assert Plate_type(96).get_well_number_from_ID('A2') == 1
        assert Plate_type(96).get_well_number_from_ID('B1') == 12
        assert Plate_type(96).get_well_number_from_ID('H12') == 95
        assert Plate_type('384').get_well_number_from_ID('A1') == 0
        assert Plate_type('384').get_well_number_from_ID('A2') == 1
        assert Plate_type('384').get_well_number_from_ID('B1') == 24
        assert Plate_type('384').get_well_number_from_ID('P24') == 383
        assert Plate_type('fake6_simple').get_well_number_from_ID('B2') == 0
        assert Plate_type('fake6_simple').get_well_number_from_ID('B7') == 1
        assert Plate_type('fake6_simple').get_well_number_from_ID('G2') == 3
        assert Plate_type('fake6_simple').get_well_number_from_ID('G11') == 5
        # for the fake6_complex type, four well IDs should have the same well number!
        assert Plate_type('fake6_complex').get_well_number_from_ID('B2') == 0
        assert Plate_type('fake6_complex').get_well_number_from_ID('B3') == 0
        assert Plate_type('fake6_complex').get_well_number_from_ID('C2') == 0
        assert Plate_type('fake6_complex').get_well_number_from_ID('C3') == 0
        assert Plate_type('fake6_complex').get_well_number_from_ID('B7') == 1
        assert Plate_type('fake6_complex').get_well_number_from_ID('G2') == 3
        assert Plate_type('fake6_complex').get_well_number_from_ID('F10') == 5
        assert Plate_type('fake6_complex').get_well_number_from_ID('G11') == 5
        
    def test__get_well_number__fromID__bad_numbers(self):
        """ For each size, test a row and column that shouldn't exist). """
        # MAYBE-TODO make.get_well_number_from_ID check for this explicitly, raise PlateTransferError instead of KeyError?
        self.assertRaises(KeyError, Plate_type('6').get_well_number_from_ID, 'A4')
        self.assertRaises(KeyError, Plate_type('6').get_well_number_from_ID, 'C1')
        self.assertRaises(KeyError, Plate_type('24').get_well_number_from_ID, 'A7')
        self.assertRaises(KeyError, Plate_type('24').get_well_number_from_ID, 'E1')
        self.assertRaises(KeyError, Plate_type(96).get_well_number_from_ID, 'A13')
        self.assertRaises(KeyError, Plate_type(96).get_well_number_from_ID, 'I1')
        self.assertRaises(KeyError, Plate_type(384).get_well_number_from_ID, 'A25')
        self.assertRaises(KeyError, Plate_type(384).get_well_number_from_ID, 'Q1')
        for P in (Plate_type('fake6_simple'), Plate_type('fake6_complex')):
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'B1')
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'B4')
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'B12')
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'A2')
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'D2')
            self.assertRaises(KeyError, P.get_well_number_from_ID, 'H2')

    ### Testing creating new plate types (legal and not)

    def test__creating_bad_plate_types_from_size(self):
        # Note: creating GOOD plate types from size was already tested in the section above
        # Shouldn't ever be able to create a 0-well plate
        self.assertRaises(PlateTransferError,Plate_type,size=0)
        self.assertRaises(PlateTransferError,Plate_type,size='0')
        # Shouldn't be able to create a 10-well plate, it wasn't defined
        self.assertRaises(PlateTransferError,Plate_type,size=10)
        self.assertRaises(PlateTransferError,Plate_type,size='10')
        # Non-integer values won't work
        self.assertRaises(ValueError,Plate_type,size='T')
        self.assertRaises(PlateTransferError,Plate_type,size=None)

    def test__creating_bad_plate_types_from_standard(self):
        # Note: creating GOOD plate types from standard_type was already tested in the section above
        # Shouldn't ever be able to create a 0-well plate
        for bad_type in [0, '0', 10, 'T', 'fake', [0], None]:
            self.assertRaises(PlateTransferError,Plate_type,standard_type=bad_type)

    def test__creating_good_plate_types_from_list(self):
        # Note: this is the only method of creating non-predefined types, so it gets its own "good" creation test
        new_plate_type = Plate_type(well_ID_list=['a','b','c'])
        assert new_plate_type.size == 3
        assert new_plate_type.get_well_ID_from_number(0) == 'a'
        assert new_plate_type.get_well_ID_from_number(2) == 'c'
        assert new_plate_type.get_well_number_from_ID('a') == 0
        assert new_plate_type.get_well_number_from_ID('c') == 2
        new_plate_type = Plate_type(well_ID_list=['A29','D1'])
        assert new_plate_type.size == 2
        assert new_plate_type.get_well_ID_from_number(0) == 'A29'
        assert new_plate_type.get_well_ID_from_number(1) == 'D1'
        assert new_plate_type.get_well_number_from_ID('A29') == 0
        assert new_plate_type.get_well_number_from_ID('D1') == 1

    def test__creating_bad_plate_types_from_list(self):
        # if size is given, it must match the lenght of well_ID_list
        self.assertRaises(PlateTransferError,Plate_type, size=1, well_ID_list=['a','b'])
        self.assertRaises(PlateTransferError,Plate_type, size=5, well_ID_list=['a','b'])
        # well_ID_list must be a list and must not contain duplicates
        for bad_list in [True, 4, 1232]:
            self.assertRaises(TypeError,Plate_type, well_ID_list=bad_list)
        for bad_list in [['a','a'], [0,0], [('a','b'), 'a'], [('a','a'), 'b']]:
            self.assertRaises(PlateTransferError,Plate_type, well_ID_list=bad_list)


class Testing__generate_outfile_names(unittest.TestCase):
    """ Unit-tests for the generate_outfile_names function. """

    def test__single_Biomek_file(self):
        assert generate_outfile_names('X',0,0) == ('X.txt',['X_Biomek.csv'],[])

    def test__different_input_formats(self):
        # file_plate_names can be specified in multiple ways
        assert generate_outfile_names('X',1,0,1,['A']) == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,1,'A') == ('X.txt',['X_Biomek_A.csv'],[])
        assert generate_outfile_names('X',1,0,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],[])
        assert generate_outfile_names('X',1,0,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],[])

    def test__mirroring(self):
        # mirroring
        assert generate_outfile_names('X',0,1) == ('X.txt',['X_Biomek.csv'],['X_Biomek_mirror.csv'])
        assert generate_outfile_names('X',1,1,1,['A']) == ('X.txt',['X_Biomek_A.csv'],['X_Biomek_mirror_A.csv'])
        assert generate_outfile_names('X',1,1,2,['A','B']) == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                               ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A,B') == ('X.txt',['X_Biomek_A.csv','X_Biomek_B.csv'],
                                                           ['X_Biomek_mirror_A.csv','X_Biomek_mirror_B.csv'])
        assert generate_outfile_names('X',1,1,2,'A') == ('X.txt',['X_Biomek_A1.csv','X_Biomek_A2.csv'],
                                                         ['X_Biomek_mirror_A1.csv','X_Biomek_mirror_A2.csv'])

    def test__second_arg_true_requires_last_two_args(self):
        # the last two args must be given if second arg is True
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,0,['A'])
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1,1)
        self.assertRaises(PlateTransferError, generate_outfile_names, 'X',1,1,['A'])


class Testing__get_plate_name_list_from_input(unittest.TestCase):
    """ Unit-tests for the get_plate_name_list_from_input function. """

    def test__input_types(self):
        # input can be a list of appropriate length (return unchanged) or a string (split on ,)
        assert get_plate_name_list_from_input(1,['A']) == ['A']
        assert get_plate_name_list_from_input(1,'A') == ['A']
        assert get_plate_name_list_from_input(4,['A','B','C','D']) == ['A','B','C','D']
        assert get_plate_name_list_from_input(4,'A,B,C,D') == ['A','B','C','D']
        # passing a one-element string is allowed for any N - plates are numbered automatically
        assert get_plate_name_list_from_input(4,'A') == ['A1','A2','A3','A4']
        # if a number >=10 is given, plate numbers should be zero-padded so that they sort correctly
        assert get_plate_name_list_from_input(10,'A') == ['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10']
        # if the second arg is a LIST, not string, with one element, fail 
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,['A'])

    def test__correct_number_of_plate_names(self):
        # first arg should match length of second arg (unless the latter is a one-element string)
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 1,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 1,['A','B','C'])
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 2,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 2,['A','B','C'])
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,'A,B,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 4,['A','B','C'])

    def test__duplicates_not_allowed(self):
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 3,'A,A,C')
        self.assertRaises(PlateTransferError, get_plate_name_list_from_input, 3,['A','A','C'])


class Testing__numbers_to_plate_and_well_IDs(unittest.TestCase):
    """ Unit-tests for the numbers_to_plate_and_well_IDs function. """

    def test__correct_cases(self):
        assert numbers_to_plate_and_well_IDs(10, 6, 2, ['plate1','plate2']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate2,A1', 'plate2,A2', 'plate2,A3', 'plate2,B1']
        assert numbers_to_plate_and_well_IDs(10, '6', 2, ['plate1','plate2']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate2,A1', 'plate2,A2', 'plate2,A3', 'plate2,B1']
        assert numbers_to_plate_and_well_IDs(10, 24, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,B1', 'plate1,B2', 'plate1,B3', 'plate1,B4']
        assert numbers_to_plate_and_well_IDs(10, 96, 1, ['plate1']) == ['plate1,A1', 'plate1,A2', 'plate1,A3', 'plate1,A4', 'plate1,A5', 'plate1,A6', 'plate1,A7', 'plate1,A8', 'plate1,A9', 'plate1,A10']
        assert numbers_to_plate_and_well_IDs(10, 'fake6_simple', 2, ['plate1','plate2']) == ['plate1,B2', 'plate1,B7', 'plate1,B11', 'plate1,G2', 'plate1,G7', 'plate1,G11', 'plate2,B2', 'plate2,B7', 'plate2,B11', 'plate2,G2']
        assert numbers_to_plate_and_well_IDs(10, 'fake6_complex', 2, ['plate1','plate2']) == ['plate1,B2', 'plate1,B6', 'plate1,B10', 'plate1,F2', 'plate1,F6', 'plate1,F10', 'plate2,B2', 'plate2,B6', 'plate2,B10', 'plate2,F2']

    def test__bad_input(self):
        # Shouldn't work, N_plates doesn't match the plate ID list
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 10, 6, 2, ['plate1'])
        # Shouldn't work, not enough plates
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 20, 6, 2, ['plate1','plate2'])
        # Shouldn't work, too many plates
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 2, 6, 2, ['plate1','plate2'])
        # Shouldn't work, 10 isn't a valid plate size
        self.assertRaises(PlateTransferError, numbers_to_plate_and_well_IDs, 2, 10, 2, ['plate1','plate2'])


class Testing__assign_codewords(unittest.TestCase):
    """ Unit-tests for the assign_codewords function. """

    def setUp(self):
        # make some test binary codes
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        [self.b01,self.b10,self.b11,self.b00] = [b01,b10,b11,b00]
        self.b_list_with_00 = [b01,b10,b11,b00]
        self.b_list_without_00 = [b01,b10,b11]
        self.B_with_00 = binary_code_utilities.Binary_code(2,self.b_list_with_00)
        self.B_without_00 = binary_code_utilities.Binary_code(2,self.b_list_without_00)
        [b01111,b10000,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['01111','10000','10001']]
        [self.b01111,self.b10000,self.b10001] = [b01111,b10000,b10001]
        self.b_list_longer = [b01111,b10000,b10001]
        self.B_longer = binary_code_utilities.Binary_code(5,self.b_list_longer)

    def test__result_is_subset_of_right_length(self):
        """ Result should be a subset of the expected set and of the expected length (don't check the ordering yet). """
        # (should work the same regardless of take_high value)
        for h in [True,False]:
            s = set(assign_codewords(3,2,self.B_without_00,take_high=h,quiet=True))
            assert s == set(self.b_list_without_00) and len(s) == 3
            s = set(assign_codewords(2,2,self.B_without_00,take_high=h,quiet=True))
            assert s.issubset(set(self.b_list_without_00)) and len(s) == 2
            s = set(assign_codewords(1,2,self.B_without_00,take_high=h,quiet=True))
            assert s.issubset(set(self.b_list_without_00)) and len(s) == 1
            s = set(assign_codewords(3,5,self.B_longer,take_high=h,quiet=True))
            assert s == set(self.b_list_longer) and len(s) == 3
            s = set(assign_codewords(2,5,self.B_longer,take_high=h,quiet=True))
            assert s.issubset(set(self.b_list_longer)) and len(s) == 2
            s = set(assign_codewords(1,5,self.B_longer,take_high=h,quiet=True))
            assert s.issubset(set(self.b_list_longer)) and len(s) == 1

    def test__all_zero_codeword_always_removed(self):
        """ The all-zero codeword should always be thrown away. """
        # (should work the same regardless of take_high value)
        for h in [True,False]:
            assert self.b00 not in set(assign_codewords(3,2,self.B_with_00,take_high=h,quiet=True))
            assert self.b00 not in set(assign_codewords(2,2,self.B_with_00,take_high=h,quiet=True))
            assert self.b00 not in set(assign_codewords(1,2,self.B_with_00,take_high=h,quiet=True))

    def test__correct_codewords(self):
        """ The codewords chosen should be in the expected bit-sum subset (low or high depending on take_high). 
        If the expected bit-sum subset contains exactly N codewords, the result is exact; otherwise a random subset.
        The codewords should be sorted lexicographically."""
        # do test multiple times, since randomness is involved
        for i in range(10):
            assert set(assign_codewords(1,2,self.B_without_00,quiet=True)).issubset([self.b01, self.b10])
            assert assign_codewords(2,2,self.B_without_00,quiet=True) == [self.b01,self.b10]
            assert assign_codewords(3,2,self.B_without_00,quiet=True) == [self.b01,self.b10,self.b11]
            assert assign_codewords(1,5,self.B_longer,quiet=True) == [self.b10000]
            assert assign_codewords(2,5,self.B_longer,quiet=True) == [self.b10000,self.b10001]
            assert assign_codewords(3,5,self.B_longer,quiet=True) == [self.b01111,self.b10000,self.b10001]
            # if take_high is set to True, the sorting by weight should be reversed.
            assert assign_codewords(1,2,self.B_without_00,take_high=True,quiet=True) == [self.b11]
            assert set(assign_codewords(2,2,self.B_without_00,take_high=True,quiet=True)).issubset([self.b01,
                                                                                                    self.b10,self.b11])
            assert assign_codewords(3,2,self.B_without_00,take_high=True,quiet=True) == [self.b01,self.b10,self.b11]
            assert assign_codewords(1,5,self.B_longer,take_high=True,quiet=True) == [self.b01111]
            assert assign_codewords(2,5,self.B_longer,take_high=True,quiet=True) == [self.b01111,self.b10001]
            assert assign_codewords(3,5,self.B_longer,take_high=True,quiet=True) == [self.b01111,self.b10000,self.b10001]

    def test__error_when_wrong_sample_or_pool_number(self):
        """ Function should fail when the number of samples or pools doesn't match. """
        # (should work the same regardless of take_high value)
        for h in [True,False]:
            # Should be too many samples for given code
            self.assertRaises(PlateTransferError, assign_codewords, 4,2,self.B_without_00,take_high=h,quiet=True)
            # Should be too many samples after '00' removal
            self.assertRaises(PlateTransferError, assign_codewords, 4,2,self.B_with_00,take_high=h,quiet=True)
            # Number of pools doesn't match code length
            self.assertRaises(PlateTransferError, assign_codewords, 4,1,self.B_without_00,take_high=h,quiet=True)
            self.assertRaises(PlateTransferError, assign_codewords, 4,3,self.B_without_00,take_high=h,quiet=True)


class Testing__make_Biomek_file_commands(unittest.TestCase):
    """ Unit-tests for the make_Biomek_file_commands function. """

    def test__basic_functionality(self):
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        # basic functionality for combinatorial pooling (note that 'x', 'A' etc here would really be 'plate1,A1' or such)
        assert make_Biomek_file_commands([b10],['x'],['A','B'],5) == ['x,A,5']
        assert make_Biomek_file_commands([b01],['x'],['A','B'],5) == ['x,B,5']
        assert make_Biomek_file_commands([b11],['x'],['A','B'],5) == ['x,A,5','x,B,5']
        assert make_Biomek_file_commands([b01,b10,b11],['x','y','z'],['A','B'],5) == ['x,B,5','y,A,5','z,A,5','z,B,5']
        assert make_Biomek_file_commands([b11,b10,b01],['x','y','z'],['A','B'],5) == ['x,A,5','x,B,5','y,A,5','z,B,5']
        assert make_Biomek_file_commands([b01,b01,b01],['x','y','z'],['A','B'],5) == ['x,B,5','y,B,5','z,B,5']

    def test__fail_for_length_mismatches(self):
        [b01,b10,b11,b00] = [binary_code_utilities.Binary_codeword(x) for x in ['01','10','11','00']]
        [b1,b10001] = [binary_code_utilities.Binary_codeword(x) for x in ['1','10001']]
        # Lengths of codewords and sample_positions are mismatched
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b11],['x','y','z'],['A','B'],5)
        # Lengths of codewords and sample_positions are mismatched
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b11],['x','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b1],['x','y','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b10001],['x','y','z'],['A','B'],5)
        # Number of pools and codeword length mismatch
        self.assertRaises(PlateTransferError, make_Biomek_file_commands, [b01,b10,b11],['x','y','z'],['B'],5)


class Testing__split_command_list_by_source(unittest.TestCase):
    """ Unit-tests for the split_command_list_by_source function. """

    def test__empty_list(self):
        """ empty list -> empty dictionary """
        assert split_command_list_by_source([]) == {}

    def test__single_plate(self):
        """ if there's only one plate, the result should be a one-item dictionary. """
        assert split_command_list_by_source(['p1,A1,x,5','p1,A2,y,5']) == {'p1':['p1,A1,x,5','p1,A2,y,5']}

    def test__multiple_plates(self):
        """ if there are multiple plates, return one dict per plate. """
        assert split_command_list_by_source(['p1,A1,x,5','p2,A1,y,5']) == {'p1':['p1,A1,x,5'], 'p2': ['p2,A1,y,5']}
        assert split_command_list_by_source(['p1,A1,x,5','p2,A1,y,5','p2,A2,y,5']) == {'p1':['p1,A1,x,5'], 
                                                                                       'p2': ['p2,A1,y,5','p2,A2,y,5']}


class Testing__split_command_list_to_max_commands(unittest.TestCase):

    def test__bad_inputs(self):
        """ N must be a positive integer (yes, I'm not checking for other wrong types, sue me). """
        for N in [-1,0]:
            self.assertRaises(PlateTransferError, split_command_list_to_max_commands, [], N)
            self.assertRaises(PlateTransferError, split_command_list_to_max_commands, ['a','b','c'], N)

    def test__max_greater_than_len(self):
        """ if N<len(list), return the original list. """
        for N in range(1,10):
            assert split_command_list_to_max_commands([], N) == []
            assert split_command_list_to_max_commands(['a','b','c'], N+3) == [['a','b','c']]

    def test__various(self):
        """ Testing specific cases. """
        assert split_command_list_to_max_commands(['a','b','c'], 1) == [['a'],['b'],['c']]
        assert split_command_list_to_max_commands(['a','b','c'], 2) == [['a','b'],['c']]
        assert split_command_list_to_max_commands(['a','b','c'], 3) == [['a','b','c']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 1) == [['a'],['b'],['c'],['d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 2) == [['a','b'],['c','d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 3) == [['a','b'],['c','d']]
        assert split_command_list_to_max_commands(['a','b','c','d'], 4) == [['a','b','c','d']]


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    parser = define_option_parser()
    if not os.access("./error-correcting_codes",os.F_OK):
        print "Error: there is not error-correcting_codes folder in this directory - can't run tests."
        return 1

    ### 1) TESTS WITH ACTUAL OUTPUT REFERENCE FILES - THE OUTPUT IS CHECKED FOR CORRECTNESS
    # test setup - giving the name/description/arguments for each test
    test_folder = "test_data"
    tests = [
             ("test_basic", "Basic test: one 96-well source plate, one 6-well destination plate, [3,2,1] code", 
              "-n7 -N3  -p1 -s96 -P1 -S6   -o          -i Source -c error-correcting_codes/3-3-1_list -q"), 
             ("test_multi-source", "Multiple source plates (two 6-well), single Biomek file", 
              "-n7 -N3  -p2 -s6  -P1 -S6   -o          -i Source -c error-correcting_codes/3-3-1_list -q"),
             ("test_multi-outfile", "Multiple Biomek files: one per source plate - two source plates (two 6-well)", 
              "-n7 -N3  -p2 -s6  -P1 -S6   -m          -i Source -c error-correcting_codes/3-3-1_list -q"),
             ("test_split-outfile", "Biomek outfile split by size (max 4 lines/file) (one 6-well source plate)", 
              "-n7 -N3  -p1 -s96 -P1 -S6   -o -x4      -i Source -c error-correcting_codes/3-3-1_list -q"),
             ("test_mirror-outfile", "Mirror Biomek file: same as test_basic but with extra mirror Biomek file", 
              "-n7 -N3  -p1 -s96 -P1 -S6   -o     -M   -i Source -c error-correcting_codes/3-3-1_list -q"),
             ("test_other-code-432", "Using a [4,3,2] code; multiple source plates, single outfile, basic", 
              "-n7 -N4  -p2 -s6  -P1 -S6   -o          -i Source -c error-correcting_codes/4-3-2_list -q"),
             ("test_432-all-outfiles", "Same as test_other-code-432 but with multiple/split/mirror outfiles", 
              "-n7 -N4  -p2 -s6  -P1 -S6   -m -x4 -M   -i Source -c error-correcting_codes/4-3-2_list -q"),
             ("test_fake6-simple-plate", "Same as test_basic but using fake6_simple destination plate", 
              "-n7 -N3  -p1 -s96 -P1 -Sfake6_simple   -o      -i Source -c error-correcting_codes/3-3-1_list -q"),
             ("test_fake6-complex-plate", "Same as test_basic but using fake6_complex destination plate", 
              "-n7 -N3  -p1 -s96 -P1 -Sfake6_complex   -o      -i Source -c error-correcting_codes/3-3-1_list -q"),
            ] 
    return1 = run_functional_tests(tests, parser, run_main_function, test_folder)
    # MAYBE-TODO right now I'm using regular expressions and compare_files_with_regex to avoid having the tests fail due to different date or some such. The right way to do this is probably with Mock library - read up on that and change to it that method some point? (See my stackoverflow question http://stackoverflow.com/questions/9726214/testing-full-program-by-comparing-output-file-to-reference-file-whats-it-calle)

    ### 2) "SMOKE TESTS" WITH NO OUTPUT REFERENCE FILES - just make sure they work and don't give errors
    # MAYBE-TODO may want to remove those after I have enough real tests above... Or at least make them optional.
    smoketest_folder = "test_data/smoke-test_outputs"
    test_runs = [("test1", "-n63  -N15 -P3     -i Source1 -o    -C error-correcting_codes/15-6-6_generator -q"),
                 ("test2", "-n63  -N15 -P3     -i Source1 -o -M -C error-correcting_codes/15-6-6_generator -q"),
                 ("test3", "-n384 -N18 -p4 -P3 -i Source  -m    -C error-correcting_codes/18-9-6_generator -q")]
    # MAYBE-TODO add name/description strings to the test cases?
    return2 = run_functional_tests(test_runs, parser, run_main_function, smoketest_folder, smoke_tests=True)

    return (return1 or return2)

### General functions (not input/output or optparse-related or testing or main), in approximate order of use

def generate_outfile_names(outfile_basename,if_multiple_files,if_mirror_files,number_of_files=None,file_plate_names=None):
    """ Given the base outfile name, generate full outfile names: (general_outfile, Biomek, Biomek_mirror).
    If if_multiple_files is false:  X -> (X.txt, [X_Biomek.csv], M)  (the last two arguments are ignored in this case)
    Otherwise, something like this:  X -> (X.txt, [X_Biomek_A.csv,X_Biomek_B.csv,...], M)  
      (assuming get_plate_name_list_from_input(number_of_files,file_plate_names) returns something like [A,B,...]) 
    If if_mirror_files is true, Biomek_mirror (M above) is []; otherwise it's the same as Biomek with a _mirror suffix."""
    ### 1) generate the basic outfile names
    main_outfile = outfile_basename+'.txt'
    ### 2) generate the basic outfile names
    if not if_multiple_files:
        # note: just ignore the last two arguments, they may be inconsistent with a single file, that's FINE. 
        outfiles_Biomek = [outfile_basename+'_Biomek.csv']
    else:
        if not (number_of_files and file_plate_names):
            raise PlateTransferError("If outputting multiple Biomek files, must specify number and names!")
        # file_plate_names here can still be a single string or any number of other things - get a list
        file_plate_names = get_plate_name_list_from_input(number_of_files, file_plate_names)
        outfiles_Biomek = [outfile_basename+'_Biomek_'+plate+'.csv' for plate in file_plate_names]
    ### 3) generate the mirror outfile names (if requested)
    if not if_mirror_files:
        outfiles_Biomek_mirror = []
    else:
        outfiles_Biomek_mirror = [name.replace('_Biomek','_Biomek_mirror') for name in outfiles_Biomek]
    return (main_outfile,outfiles_Biomek,outfiles_Biomek_mirror)


def get_plate_name_list_from_input(N_plates,ID_input):
    """ Return a list of plate names of length N_plates, generated using ID_input. 
    If ID_input is already a correct list, return it.  Otherwise assume it's a string: if splitting it 
    on commas yields N results, return those; if it yields a single result X, return [X1,X2,...,XN].
    Note that if N_plates is 1, a count won't be appended: (1,'X') returns ['X'], not ['X1']. """
    # if we already have a list with the right number of arguments, just return it after checking it makes sense
    if isinstance(ID_input,list): 
        if not len(ID_input)==len(set(ID_input)):
            raise PlateTransferError("You specified a list of plate names with duplicate values! Don't do that.")
        if not len(ID_input)==N_plates:
            raise PlateTransferError("Passed a plate name list of wrong length to get_plate_name_list_from_input!")
        return ID_input
    # otherwise generate list from string:
    split_IDs = ID_input.split(',')
    # either it's a comma-separated 'list' of plate names, of correct length, and should be treated as above
    #  (testing for len==N_plates first: (1,'X') returns ['X'], not ['X1'] - easier, since (1,['X']) returns ['X'] too.)
    if len(split_IDs) == N_plates:  
        if not len(split_IDs)==len(set(split_IDs)):
            raise PlateTransferError("You specified a list of plate names with duplicate values! Don't do that.")
        return split_IDs
    # or it's a single plate name and sequential numbers should be appended to it to get the list.
    elif len(split_IDs) == 1:
        N_digits = len(str(N_plates))
        return ['%s%0*d'%(ID_input, N_digits, n+1) for n in range(N_plates)]
    else:
        raise PlateTransferError("Can't figure out how to name %s plates using input \"%s\"!"%(N_plates,ID_input))


def numbers_to_plate_and_well_IDs(N_samples, plate_type_name, N_plates, plate_IDs, zero_padding=False):
    """ Given the number of samples and plate information, return a list of with plate/well positions for each sample.

    The returned will be of length N_samples, about like this: ['plate1,A1','plate1,A2',...,'plate1,H12','plate2,A1',...]. 

    N_samples is the total number of samples; N_plates is the number of plates, which should match the length of plate_IDs; 
     plate_IDs is a list of arbitrary plate ID strings (like ['plate1', 'plate2']);
     plate_type_name should be one of defined_plate_types_str.
    If zero_padding is True, the wells will be A01 instead of A1.
    """

    if not len(plate_IDs)==N_plates:
        raise PlateTransferError("The number of plates must match the number of plate IDs provided!")

    # Make a separate Plate_type instance for each destination plate, just in case they have special behavior
    plate_instances = [Plate_type(standard_type=plate_type_name) for _ in range(N_plates)]
    plate_size = plate_instances[0].size

    if N_samples > plate_size*N_plates:
        raise PlateTransferError("Can't fit %s samples in %s %s-well plates!"%(N_samples,N_plates,plate_type_name))
    if N_samples <= plate_size*(N_plates-1):
        raise PlateTransferError("Why use %s %s-well plates "%(N_plates,plate_type_name)
                                 + "when you can fit %s samples in %s plates?"%(N_samples,N_plates-1))

    position_list = []
    for i in range(N_samples):
        plate_number = i // plate_size
        well_number = i % plate_size
        plate_ID = plate_IDs[plate_number]
        well_ID = plate_instances[plate_number].get_well_ID_from_number(well_number)
        if zero_padding:
            well_ID = well_ID[0] + '%02d'%int(well_ID[1:])
        position_list.append("%s,%s"%(plate_ID,well_ID))

    return position_list

# TODO TODO TODO this isn't good! If I pass position lists around like that, each transfer will always use the same position, even when it's a fake6_complex plate type where different transfers should use different positions!


def assign_codewords(N_samples, N_pools, binary_code, take_high=False, quiet=False):
    """ Return a list of Binary_codeword objects, of length N_samples - one codeword per sample, ordered. """

    if not N_pools == binary_code.length:
        raise PlateTransferError("The codeword length in the provided binary code must be equal to N_pools!")
    if N_samples > binary_code.size():
        raise PlateTransferError("N_samples cannot exceed the size of the provided binary code!")
    if binary_code.remove_extreme_codeword(bit=0):
        if not quiet:
            print("Warning: the binary code provided contained the all-zero codeword, which cannot be used, "
                  + "as it would result in a sample being added to none of the pools and thus not sequenced. "
                  + "The all-zero codeword was removed.")
        if N_samples > binary_code.size():
            raise PlateTransferError("After removing the all-zero codeword, there aren't enough codewords for all "
                                     + "samples! Aborting. You could consider inverting the code to get around this.")
    if N_samples < binary_code.size() and not quiet:
        print("Warning: N_samples is lower than the size of the provided binary code - an arbitrary subset of codewords "
              + "will be used.  You may want to reduce your code size manually for improved attributes.")

    # get the desired number of samples from the binary code (returns a sorted list, original code is unchanged)
    #  (taking either the low-weight or high-weight codewords depending on the take_high argument value)
    codeword_list = sorted(list(binary_code.give_N_codewords_by_bit_sum(N_samples,take_high=take_high)))
    return codeword_list


def make_Biomek_file_commands(sample_codewords, sample_positions, pool_positions, volume):
    """ Return a list of Biomek transfer commands to perform combinatorial pooling based on sample_codewords.

    Inputs:
     - sample_codewords - sequential list of Binary_codeword object corresponding to each sample
     - sample_positions - list of plate/well position strings (like "Source1,A4") for each sample, in the same order
     - pool_positions - same-format list of pool position strings
     - volume - integer giving the volume of all the transfers 

    Combinatorial pooling: the pools correspond to each bit of the codeword.  Sample A should be added to pool X
    whenever bit X of the codeword for sample A is 1. 
    Biomek command list format:  a list of strings of the form "plateA,wellA,plateX,wellX,volume" 
    where "plateA,wellA" is the value of sample_positions[A], and "plateX,wellX" is pool_positions[X]. """

    # make sure the inputs make sense
    if not len(sample_codewords)==len(sample_positions):
        raise PlateTransferError("The number of sample positions doesn't match the number of codewords!")
    if not set([len(x) for x in sample_codewords]) == set([len(pool_positions)]):
        raise PlateTransferError("Not all codeword lentgths match the number of pools *%s)!"%len(pool_positions))
        # note that the second set is always of size 1, so this implicitly makes sure all codewords are the same length

    Biomek_file_commands = []
    for (sample_number, (sample_codeword, sample_position)) in enumerate(zip(sample_codewords,sample_positions)):
        pools_to_add_sample_to = [pool_number for (pool_number,if_add) in enumerate(sample_codeword.list()) if if_add==1]
        for pool_number in pools_to_add_sample_to:
            Biomek_file_commands.append("%s,%s,%s"%(sample_position,pool_positions[pool_number],volume))
    return Biomek_file_commands


def split_command_list_by_source(Biomek_file_commands):
    """ Split list of "x,_,_,_,_" strings into multiple lists by x value, return a (x_val: line_list) dictionary."""
    data_dict = defaultdict(lambda: [])
    for line in Biomek_file_commands:
        source_plate = line.split(',')[0]
        data_dict[source_plate].append(line)
    return data_dict


def split_command_list_to_max_commands(Biomek_file_commands, max_lines=300):
    """ Split list of strings into multiple lists no longer than max_lines; return list of those lists."""
    if max_lines<=0: raise PlateTransferError("max_lines must be a positive integer!")
    if len(Biomek_file_commands)==0:    return []
    N_lines = len(Biomek_file_commands)
    N_lists = int(ceil(float(N_lines)/max_lines))
    # rather than just take max_lines until nothing is left, split the original list up evenly: 
    #  even if max_lines is 3, len4 should become [len2,len2] rather than [len3,len1]
    N_lines_per_list = int(ceil(float(N_lines)/N_lists))
    #print "%s lines, %s max lines per list -> %s lists with <=%s lines"%(N_lines, max_lines, N_lists, N_lines_per_list)
    new_lists = []
    for i in range(N_lists):
        new_lists.append(Biomek_file_commands[i*N_lines_per_list : (i+1)*N_lines_per_list])
    return new_lists


### Input/output functions - no need/ability to unit-test, all the complicated functionality should be elsewhere.

def get_binary_code(length,listfile=None,matrixfile=None):
    """ Given a listfile or matrixfile name (but not both!), return the generated binary code."""
    if not (listfile or matrixfile):
        raise PlateTransferError("You must provide either a listfile or a matrixfile to generate a binary code!")
    if listfile and matrixfile:
        raise PlateTransferError("You cannot provide both a listfile and a matrixfile to generate a single binary code!")
    if listfile:
        method = 'listfile'
        infile = listfile
    elif matrixfile:
        method = 'matrixfile'
        infile = matrixfile
    binary_code = binary_code_utilities.Binary_code(length=length, val=infile, method=method)
    return binary_code


def write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, max_commands_per_file=0, 
                               Biomek_header="", quiet=False):
    """ Write Biomek_file_commands to outfiles_Biomek, optionally splitting; return filename:real_filename(s) dict.  

    Biomek_file_commands and outfiles_Biomek must be lists of matching length (giving file contents and filenames).
    Each output file will start with the header line (Biomek_header argument), then all the command lines.

    The outfiles_Biomek argument must be a list: containing a single element if there will be one outfile 
      (in which case Biomek_file_commands should be a single list of command strings), 
     or more for multiple outfiles (in which case Biomek_file_commands should be a dict of matching length, 
      with each key being the plate name which should match the Biomek file name, and each value being 
      a list of command strings to be written to the corresponding Biomek file). 

    If max_commands_per_file is 0, each command list is simply written to the corresponding file; if it's N>0, 
     each command list is split into sublists with at most N lines each, and written to files with _a/_b/_c/... suffixes.
    
    The return value is a outfile_name:real_outfile_name(s) dictionary: the keys will be the elements of outfiles_Biomek,
     and the values will be either tuples of *_partXofY partial files due to splitting into max_commands_per_file, 
     or identical to the keys if no splitting was one (i.e. max_commands_per_file was 0).
    """
    ### First make a list of which commands should go with which outfile: (the actual printing is done at the end)
    data_filename_list = []
    # if there's just one outfile, simply print the whole command list to it
    if len(outfiles_Biomek)==1:
        data_filename_list.append((Biomek_file_commands, outfiles_Biomek[0]))
    # if there are multiple outfiles, split the commands by source and print subsets to corresponding outfiles
    else:
        # split the commands by source
        transfer_file_command_sets = split_command_list_by_source(Biomek_file_commands)
        # Sort both the command set list (make list from dict first) and the Biomek outfile set.  
        #  (They're based on the same underlying plate names passed to the function, plus invariant prefixes/suffixes, 
        #  so they should always match once they're sorted, even if the plate names weren't sorted sensibly themselves.)
        transfer_file_command_sets = sorted(list(transfer_file_command_sets.items()))
        outfiles_Biomek.sort()
        # make sure the resulting lists match by length (if there are more files than command sets, it may be all right)
        if len(transfer_file_command_sets) > len(outfiles_Biomek):
            raise PlateTransferError("ERROR: More Biomek command sets than outfile names were provided - can't write all!"
                                     +"\n%s command sets, %s outfiles (%s)"%(len(transfer_file_command_sets), 
                                                                             len(outfiles_Biomek), outfiles_Biomek))
        elif len(transfer_file_command_sets) < len(outfiles_Biomek) and not quiet:
            print("WARNING: Fewer Biomek command sets than outfile names were provided - some outfiles will be empty! "
                  +"(may RARELY be expected, for mirror files if the all-ones keyword was present in original file)"
                  +"\n%s command sets, %s outfiles (%s)"%(len(transfer_file_command_sets), 
                                                          len(outfiles_Biomek), outfiles_Biomek))
        # for each command set, find a single matching outfile by name (raise exception if found none/multiple), 
        #  and queue the data to be written to that file.
        outfiles_matched = set()
        for (set_name,command_list) in transfer_file_command_sets:
            matching_outfiles = [outfile for outfile in outfiles_Biomek if set_name in outfile]
            if not matching_outfiles:
                raise PlateTransferError("Can't match the command set %s to a Biomek outfile! (outfiles: %s)"%(set_name, 
                                                                                                         outfiles_Biomek))
            if len(matching_outfiles)>1:
                raise PlateTransferError("The command set %s matched to multiple Biomek outfiles! (%s)"%(set_name, 
                                                                                                       matching_outfiles)
                                         +" - if you generated your source plate names by hand, try changing them.")
            # MAYBE-TODO multiple matches may be an issue if the user gives a list of source plate names and they're of different lengths - for instance Plate1 will match both Plate1 and Plate10 (when my program auto-generates numbered plates, they're Plate01, so this won't happen).  How to correct for that?  Eh, let the user see the warning and rename her plates.
            outfile = matching_outfiles[0]
            outfiles_matched.add(outfile)
            data_filename_list.append((command_list, outfile))
        # for outfiles with no matching contents, just queue them to be written empty
        for outfile in set(outfiles_Biomek) - outfiles_matched:
            data_filename_list.append(([], outfile))

    ### now for each (data,filename) tuple, actually write it to a file (or to multiple files if splitting is necessary)
    # keep track of the final output filenames in a outfile_Biomek:final_outfile(s) dictionary
    #  (the names can change if file needs to be split due to max_commands_per_file)
    final_output_filenames = {}
    for command_list,filename in data_filename_list:
        # if there's no line-count max per file, just write lines to file
        if max_commands_per_file==0:
            save_line_list_as_file(command_list, filename, header=Biomek_header)
            final_output_filenames[filename] = filename
        # otherwise split lines into multiple files with <X lines each, with _n/N filename suffixes.
        #  (if the file has few enough lines not to be split, just give it a _1/1 suffix to make that clear)
        else:
            command_sublists = split_command_list_to_max_commands(command_list, max_commands_per_file)
            # force printing empty files - make command_sublists contain an empty list rather than being empty itself
            if not command_sublists:    command_sublists = [[]]
            N_total_files = len(command_sublists)
            N_digits = len(str(N_total_files))
            final_output_filenames[filename] = []
            for n,command_sublist in enumerate(command_sublists):
                basename, ext = os.path.splitext(filename)
                file_suffix = "part%0*dof%d"%(N_digits, n+1, N_total_files)
                filename_with_suffix = basename + '_' + file_suffix + ext
                save_line_list_as_file(command_sublist, filename_with_suffix, header=Biomek_header)
                final_output_filenames[filename].append(filename_with_suffix)
    return final_output_filenames


def write_data_to_outfile(main_outfile, sample_codewords, sample_positions, pool_positions, outfiles_Biomek, 
                          mirror_sample_codewords=[], mirror_pool_positions=[], transfer_volume=0, options=None):
    """ Write data to main_outfile: header, detailed sample/pool data, info on Biomek outfiles and overall counts/volumes.

    Header information: command, path, date/time, options - all as #-start comments. 
    Sample/pool data: tab-separated tables (with headers) containing numbers, plate/well positions, codewords or pooling 
     schemes, transfer counts and total transfer volumes - one table for samples, one for pools, one for mirror pools. 
    Footer: list of corresponding Biomek command files, info on total sample/pool/mirrorpool numbers, 
     info on the min/max transfer/count/volume for samples/pools.
    """
    ### print all the usual header information (command, path, date/time, options)
    OUTFILE = open(main_outfile,'w')
    write_header_data(OUTFILE,options)
    ### write the detailed sample/pool number/position/codeword/volume data
    # sample data
    sample_transfers = {}
    OUTFILE.write("\n")
    for setname, curr_sample_codewords in [('', sample_codewords), ('mirror', mirror_sample_codewords)]:
        if curr_sample_codewords:     # only print the header if there's any content
            OUTFILE.write("sample_number\tplate_and_well_position\t%scodeword\ttransfers\tvolume (ul)\n" 
                          %('' if setname=='' else setname+'_'))
            sample_transfers[setname] = []
        for (number,(codeword,position)) in enumerate(zip(curr_sample_codewords,sample_positions)):
            total_transfers = codeword.weight()
            total_volume = total_transfers * transfer_volume
            OUTFILE.write("%s\t%s\t%s\t%s\t%s\n"%(number, position, codeword.string(), total_transfers, total_volume))
            sample_transfers[setname].append(total_transfers)
    # pool data (first normal, then mirror, with a header for each)
    pool_transfers = {}
    OUTFILE.write("\n")
    for setname, curr_pool_positions, curr_sample_codewords in [('', pool_positions, sample_codewords), 
                                                            ('mirror', mirror_pool_positions, mirror_sample_codewords)]:
        if curr_pool_positions:     # only print the header if there's any content
            OUTFILE.write("%spool_number\tplate_and_well_position\tpooling_scheme\ttransfers\tvolume (ul)\n" 
                          %('' if setname=='' else setname+'_'))
            pool_transfers[setname] = []
        for (number,position) in enumerate(curr_pool_positions):
            pooling_scheme = ''.join([codeword.string()[number] for codeword in curr_sample_codewords])
            # MAYBE-TODO implement this pooling_scheme thing as a binary_code_utilities.py Binary_code method?
            total_transfers = pooling_scheme.count('1')
            total_volume = total_transfers * transfer_volume
            OUTFILE.write("%s\t%s\t%s\t%s\t%s\n"%(number, position, pooling_scheme, total_transfers, total_volume))
            pool_transfers[setname].append(total_transfers)
    ### print footer info: corresponding Biomek outfile list, min/max transfers/volume per sample/pool
    # make nice outfile list for printing: strip outermost [] pair (with [1:-1]), get rid of quotes, 
    #  and remove the folder name, since they're in the same folder as the main_outfile
    nice_outfile_list = str(outfiles_Biomek)[1:-1].replace("'",'').replace(os.path.dirname(main_outfile)+os.path.sep,'')
    OUTFILE.write("\n# Corresponding Biomek command file(s): %s\n"%nice_outfile_list)
    OUTFILE.write("# Total %s samples into %s pools (and %s mirror pools)\n"%(len(sample_positions), len(pool_positions), 
                                                        len(mirror_pool_positions)))
    for setname, curr_sample_transfers in sample_transfers.items():
        min_transfers, max_transfers = min(curr_sample_transfers), max(curr_sample_transfers)
        OUTFILE.write("%stransfers from samples: "%('' if setname=='' else setname+' '))
        OUTFILE.write("%s-%s per sample (%s-%s ul), "%(min_transfers, max_transfers, 
                                                       min_transfers*transfer_volume, max_transfers*transfer_volume))
        OUTFILE.write("total %s transfers\n"%sum(curr_sample_transfers))
    for setname, curr_pool_transfers in pool_transfers.items():
        min_transfers, max_transfers = min(curr_pool_transfers), max(curr_pool_transfers)
        OUTFILE.write("transfers into %spools: "%('' if setname=='' else setname+' '))
        OUTFILE.write("%s-%s per pool (%s-%s ul), "%(min_transfers, max_transfers, 
                                                    min_transfers*transfer_volume, max_transfers*transfer_volume))
        OUTFILE.write("total %s transfers\n"%sum(curr_pool_transfers))
    OUTFILE.close()
    # MAYBE-TODO print some info about the minimum Hamming distance?


### Option parser functions - no real need to unit-test, and it'd be complicated and weird.

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    parser.add_option('-q','--quiet', action='store_true', default=False, help="Don't print warnings (default False).")
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default False).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run with a set of predetermined realistic options to make sure it more or less works "
                      + "(default False). Ignores all other options/arguments; output files will start with 'test'.")

    parser.add_option('-n','--number_of_samples', type='int', metavar='N', help="Number of samples to pool (required).")
    parser.add_option('-N','--number_of_pools', type='int', metavar='M', help="Number of resulting pools (required).")
    parser.add_option('-c','--binary_code_list_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a list of codewords).")
    parser.add_option('-C','--binary_code_generator_file', metavar='FILE', 
                      help="File containing the binary code to use for the pooling (as a generator matrix).")
    parser.add_option('-M','--add_mirror_pooling_files', action='store_true', default=False, 
                      help="In addition to the normal Biomek file, also make files with commands for a 'mirrored' set: " 
                      + "if sample A is in pool B in the normal set it isn't in the mirrored set, and vice versa.")
    parser.add_option('-u','--mirror_pool_plate_suffix', default='_mirror', metavar='S', 
                      help="Append S to the pool plate names for the mirror pooling Biomek files (default %default).")

    parser.add_option('-s','--size_of_sample_plates', type='choice', choices=defined_plate_types, default='96',
                      help="Sample (source) plate size (allowed values: %s) "%defined_plate_types_str 
                          +"(default %default)", metavar='M')
    parser.add_option('-S','--size_of_pool_plates',type='choice', choices=defined_plate_types, default='fake6_simple',
                      help="Pool (destination) plate size (allowed values: %s) "%defined_plate_types_str 
                          +"(default %default)", metavar='M')
    parser.add_option('-p','--number_of_sample_plates', type='int', default=1, metavar='N', 
                      help="Total number of sample (source) plates to use (default %default).")
    parser.add_option('-P','--number_of_pool_plates', type='int', default=4, metavar='N', 
                      help="Total number of pool (destination) plates to use (default %default).")
    parser.add_option('-i','--sample_plate_IDs', default='Source', metavar='S', 
                      help="Sample plate IDs (must match the IDs in the Biomek deck setup): either a comma-separated list "
                      + "with the number of values matching -n ('A,B,C,D' - plates will be named A, B, C and D), "
                      + "or a single string with no commas ('x' - plates will be named x1, x2 etc.) (default %default).")
    parser.add_option('-I','--pool_plate_IDs', default='Destination', metavar='S', 
                      help="Pool plate IDs - see -i for details on how it works. Default %default.")

    parser.add_option('-m','--multiple_Biomek_files', action='store_true', default=True,
                      help="Generate multiple Biomek files, one per source plate (on by default).")
    parser.add_option('-o','--one_Biomek_file', action='store_false', dest='multiple_Biomek_files',
                      help="Generate a single Biomek file regardless of number of plates involved (off by default).")
    parser.add_option('-x','--max_commands_per_file', type='int', default=0, metavar='N',
                      help="Split Biomek command files so they contain to more than N commands "
                          +"(use a/b/c/... suffixes for the split files) (0 means no maximum; default %default).")

    parser.add_option('-v','--volume_per_transfer', type='int', default=20, metavar='V', 
                      help="Liquid volume to use for each sample-to-pool transfer (default %default).")
    parser.add_option('-H','--Biomek_file_header', default="SourcePlt,SourceWell,DestPlt,DestWell,Volume", metavar='S', 
                      help="Header line for the Biomek transfer file (won't be printed if set to ''; default %default).")

    # MAYBE-TODO add max sample and pool volume options to make the program auto-check that they aren't exceeded

    # MAYBE-TODO add a minimum Hamming distance option, make the program check that all pairs in the final codeword set satisfy that?
    # MAYBE-TODO same for bit-sum limit for when I want to specify that

    # MAYBE-TODO implement more ways of dealing with the clonality issue?  
    #   1. Simple check - check that the bitwise OR of no two codewords generates a codeword that's in the set
    #   2. More complicated check - minimum Hamming distance we require between any codeword and any bitwise OR of two codewords (would have to be specified as an option).  Or just get a histogram of that printed. 
    #   3. Actually try to get a set that doesn't have clonality issues? (based on some minimum distance given as an option).  (The actual mechanism for this would have to be implemented in binary_code_utilities.py - reduce_by_Hamming_distance or something else, see possible methods in the notes above)

    return parser


def check_options_and_args(parser,options,args):
    """ Take optparse parser/options/args, check number of args, check required options and option conflicts. """
    # MAYBE-TODO rewrite this so it doesn't require the parser as an argument, by moving the first try/except to __main__?  Why is it a problem, really? 

    try:
        [outfile_basename] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: There must be exactly one output file base name (shown as X in the examples). "
                 + "The general output file will be X.txt, the Biomek output file will be X_Biomek.csv "
                 + "(or X_Biomek_Source1.csv, X_Biomek_Source2.csv, etc, with the -m option).")

    for curr_plate_size in [options.size_of_sample_plates, options.size_of_pool_plates]:
        if not curr_plate_size in defined_plate_types:
            sys.exit("Plate sizes (-s and -S) must be one of the defined sizes (%s)!"%defined_plate_types_str)
    if not (options.number_of_samples>0 and options.number_of_pools>0):
        sys.exit("Positive -n and -N values required!")
    if not bool(options.binary_code_list_file) ^ bool(options.binary_code_generator_file):  # is xor
        sys.exit("Exactly one of -c and -C must be provided!")

    # MAYBE-TODO could allow -p/-P to be automatically calculated from -n/-N and -s/-S?

    return options,outfile_basename


### Main functionality - difficult to unit-test due to optparse options as an argument, can just be tested with full runs.

def run_main_function(parser,options,args):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    Takes an optparse parser object, and the options object and arg list generated by the parser."""
    # MAYBE-TODO may be more convenient for interactive use if this just took an input string, and generated/defined the parser itself...
    options,outfile_basename = check_options_and_args(parser,options,args)
    outfiles = generate_outfile_names(outfile_basename, options.multiple_Biomek_files, options.add_mirror_pooling_files, 
                                      options.number_of_sample_plates, options.sample_plate_IDs)
    (main_outfile, outfiles_Biomek, outfiles_Biomek_mirror) = outfiles
    # assign codewords to samples
    binary_code = get_binary_code(options.number_of_pools, 
                                  options.binary_code_list_file, options.binary_code_generator_file)
    sample_codewords = assign_codewords(options.number_of_samples, options.number_of_pools, binary_code, 
                                        quiet=options.quiet)
    # generate plate names from strings if they weren't given as lists
    input_plate_names = get_plate_name_list_from_input(options.number_of_sample_plates, options.sample_plate_IDs)
    output_plate_names = get_plate_name_list_from_input(options.number_of_pool_plates, options.pool_plate_IDs)
    # generate the plate+well position strings for each input sample and each output pool
    sample_positions = numbers_to_plate_and_well_IDs(options.number_of_samples, options.size_of_sample_plates, 
                                                     options.number_of_sample_plates, input_plate_names)
    pool_positions = numbers_to_plate_and_well_IDs(options.number_of_pools, options.size_of_pool_plates, 
                                                   options.number_of_pool_plates, output_plate_names)
    # generate the Biomek transfer command list based on sample codewords and sample/pool positions
    Biomek_file_commands = make_Biomek_file_commands(sample_codewords, sample_positions, pool_positions, 
                                                     options.volume_per_transfer)
    # optionally generate mirror Biomek files: invert the codewords, add suffix to pool plate names, run same functions.
    if options.add_mirror_pooling_files:
        mirror_sample_codewords = [~codeword for codeword in sample_codewords]
        mirror_output_plate_names = [plate_name+options.mirror_pool_plate_suffix for plate_name in output_plate_names]
        mirror_pool_positions = numbers_to_plate_and_well_IDs(options.number_of_pools, options.size_of_pool_plates, 
                                                              options.number_of_pool_plates, mirror_output_plate_names)
        mirror_Biomek_file_commands = make_Biomek_file_commands(mirror_sample_codewords, sample_positions, 
                                                                mirror_pool_positions, options.volume_per_transfer)
    else:
        mirror_sample_codewords, mirror_pool_positions = [], []

    ### write data to outfiles, keeping track of real outfile names as returned by data-writing functions
    Biomek_real_outfile_dict = {main_outfile: main_outfile}
    # write commands to Biomek outfiles (normal, and optionally mirror)
    Biomek_normalfile_dict = write_data_to_Biomek_files(outfiles_Biomek, Biomek_file_commands, 
                                        options.max_commands_per_file, options.Biomek_file_header, quiet=options.quiet)
    Biomek_real_outfile_dict.update(Biomek_normalfile_dict)
    if options.add_mirror_pooling_files:
        Biomek_mirrorfile_dict = write_data_to_Biomek_files(outfiles_Biomek_mirror, mirror_Biomek_file_commands, 
                                        options.max_commands_per_file, options.Biomek_file_header, quiet=options.quiet)
        Biomek_real_outfile_dict.update(Biomek_mirrorfile_dict)
    # make nice sorted list of real Biomek outfiles (normal and mirror) to write to main_outfile and return
    outfiles_Biomek = [Biomek_real_outfile_dict[f] for f in outfiles_Biomek]
    outfiles_Biomek_mirror = [Biomek_real_outfile_dict[f] for f in outfiles_Biomek_mirror]
    # write full data (including a header, all Biomek outfile names, samples/destinations/codewords) to main_outfile
    write_data_to_outfile(main_outfile, sample_codewords, sample_positions, pool_positions, 
                          outfiles_Biomek+outfiles_Biomek_mirror, mirror_sample_codewords, mirror_pool_positions, 
                          options.volume_per_transfer, options)
    # return (and optionally print) list of all the outfiles generated
    if not options.quiet:  
        print("Overview output file: %s"%main_outfile)
        print("Biomek output files: %s, %s"%(outfiles_Biomek, outfiles_Biomek_mirror))
    return (main_outfile, outfiles_Biomek, outfiles_Biomek_mirror)


if __name__=='__main__':

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # Unit-testing - don't even look for more options/arguments, just run the test suite
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments (including -T), "
              + "running the built-in simple test suite.")
        print("Defined plate sizes: %s"%defined_plate_types_str)
        # tun unittest.main, passing it no arguments (by default it takes sys.argv and complains about the -t)
        unittest.main(argv=[sys.argv[0]])
        # MAYBE-TODO make it run unit-tests for everything it imports too? Especially binary_code_utilities...
        # MAYBE-TODO unittest.main automatically quits - there's no way of doing -t and -T at once.  Do I care?
        #   may be fixed in a future version, and there is a patch: http://bugs.python.org/issue3379


    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs.")
        test_result = do_test_run()
        sys.exit(test_result)

    # If it's not a test run, just run the main functionality
    run_main_function(parser,options,args)
