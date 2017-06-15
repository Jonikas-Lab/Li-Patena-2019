#! /usr/bin/env python2.7

"""
Various general programming/math utilities I wrote - see docstring for each function for what it does.
Weronika Patena, 2010-2011
"""

# standard library
from __future__ import division 
import sys, os
import unittest
from collections import defaultdict
import pickle as Pickle     # rename so I can use pickle as a function name
import subprocess
import time
# other packages
# my modules

######################################## GENERAL CONVENIENCE STUFF ########################################

def timenow():
    """ Print the current time as a string. """
    print time.ctime(time.time())


######################################## STRUCTURES / CLASSES / ETC ########################################

### list/dictionary/etc utilities

def compare_lists_unordered(list1,list2):
    """ Given two lists, return True if they contain the same elements, False otherwise.
    The same could be done with bool(set(list1)==set(list2)), except when list elements are unhashable. """
    if len(list1)!=len(list2):  return False
    for x in list1:
        if x not in list2:      return False
    return True

def add_to_dict_no_replacement(dictionary, key, new_val, key_name='key', val_name='value', raise_exception=True, print_error=True):
    """ If key isn't in dictionary yet, add it with val; if already has a different value, raise exception or print error message."""
    if key in dictionary and dictionary[key] != new_val:
        error_msg = "Attempting to change existing %s for %s %s! (old %s, new %s)"%(val_name,key_name, key, dictionary[key], new_val)
        if raise_exception:   raise ValueError(error_msg)
        if print_error:       print "ERROR: %s. NOT CHANGING the dictionary."%error_msg
    else:
        dictionary[key] = new_val

def reduce_dicts_to_overlaps(dict_list, exception_on_different_values=False):
    """ Given a list of dictionaries, return a similar new list of dicts with only the keys shared by all the dictionaries.
    Caution: returns a LIST OF DICTS, NOT a new dictionary that contains only the overlap! Easy to get confused."""
    # get the set of keys present in all dictionaries (set intersection)
    if len(dict_list)==0:   return []
    # generate list of keys that are in ALL dictionaries (starting from the list of the first one and reducing)
    overlap_keys = set(dict_list[0].keys())
    for d in dict_list:     overlap_keys &= set(d.keys())
    # make a new dictionary list
    new_dict_list = []
    for d in dict_list:     
        new_reduced_dict = dict([(k,v) for (k,v) in d.iteritems() if k in overlap_keys])
        new_dict_list.append(new_reduced_dict)
    return new_dict_list

def filter_dict_by_keys(input_dict, good_keys):
    """ Return copy of input_dict with all keys not in good_keys removed. """
    # MAYBE-TODO add a bad_keys arg, to do either positive or negative filtering?
    # make good_keys into a set for speed, if it wasn't one already
    good_keys = set(good_keys)
    new_dict = {key:val for (key,val) in input_dict.iteritems() if key in good_keys}
    return new_dict

def count_list_values(input_list):
    """ Given a list, return a value:number_of_times_value_occurred_in_list dictionary. """
    value_counts = defaultdict(lambda: 0)
    for value in input_list:
        value_counts[value] += 1
    return dict(value_counts)
    # Note: this is basically a very simplified version of collections.Counter - which was only added in python 2.7, and I'm still on 2.6, so I can't use it for now.  Should switch to that once I'm on 2.7 though.. 

def add_dicts_of_ints(dict1, dict2, recursive=False, 
                      allow_nonconflicting_values=False, remove_nonallowed_values=False):
    """ Add int-valued dicts together by adding the values: given {a:1, b:1} and {a:2, c:2}, return {a:3, b:1, c:2}. 

    If recursive is True, if both dicts have another dict as a value for a given key, recursively apply 
     add_dicts_of_ints to those values.
    Normally only integer values are allowed (and possibly dictionary values, if recursive is True), and other values
     cause a ValueError. If remove_nonallowed_values is True, other values are silently ignored instead; 
     if allow_nonconflictiong_values is True, other values are kept if the key is only present in one dict
      or if the value is the same in both dicts.
    This uses explicit type checks for ints, not duck typing, sorry - don't want to add strings/lists/etc! 
    """
    # MAYBE-TODO add an argument that'd be a list of other types that should be treated like ints (i.e. added)?
    keys_1_only = set(dict1.keys()) - set(dict2.keys())
    keys_2_only = set(dict2.keys()) - set(dict1.keys())
    keys_both = set(dict2.keys()) & set(dict1.keys())
    new_dict = {}
    ### for keys in both dicts:
    for key in keys_both:
        # if both values are ints, use their sum as new value
        if type(dict1[key]) == type(dict2[key]) == int:
            new_dict[key] = dict1[key] + dict2[key]
        # if recursive is True and both values are dicts, run add_dicts_of_ints on them with same options for new value
        elif recursive and type(dict1[key]) == type(dict2[key]) == dict:
            new_dict[key] = add_dicts_of_ints(dict1[key], dict2[key], recursive=recursive, 
                                              allow_nonconflicting_values=allow_nonconflicting_values, 
                                              remove_nonallowed_values=remove_nonallowed_values)
        # if allow_nonconflicting_values is True and both values are the same, use that for new value
        elif allow_nonconflicting_values and dict1[key] == dict2[key]:
            new_dict[key] = dict1[key]
        # if the values aren't any of the above cases and remove_nonallowed_values is True:
        elif remove_nonallowed_values:
            # if one is an allowed type and the other isn't, keep the allowed one and ignore the other
            #  (really the allowed one should just be moved to the keys_*_only set, in case further checking is needed)
            if recursive:       allowed_types = [int,dict]
            else:               allowed_types = [int]
            if type(dict1[key]) in allowed_types and type(dict2[key]) not in allowed_types:
                keys_1_only.add(key)
            elif type(dict2[key]) in allowed_types and type(dict1[key]) not in allowed_types:
                keys_2_only.add(key)
            # otherwise (if neither are an allowed type, or both are but they're different, like int/dict) just ignore both
            else:
                continue
        # otherwise raise exception
        else:
            raise ValueError("Encountered non-allowed value in one of the dictionaries! See docstring for options.")
    ### for keys that are only in one of the dicts:
    for key_set, curr_dict in [(keys_1_only, dict1), (keys_2_only, dict2)]:
        for key in key_set:     
            # copy ints to new_dict; if allow_nonconflicting_values is True, copy any values without checking
            if allow_nonconflicting_values or type(curr_dict[key])==int:
                new_dict[key] = curr_dict[key]
            # if recursive is True and the value is a dict, apply add_dicts_of_ints to it and {} to get new value
            elif recursive and type(curr_dict[key])==dict:
                new_dict[key] = add_dicts_of_ints(curr_dict[key], {}, recursive=recursive, 
                                                  allow_nonconflicting_values=allow_nonconflicting_values, 
                                                  remove_nonallowed_values=remove_nonallowed_values)
            # if the value isn't any of the above and remove_nonallowed_values is True, just ignore it
            elif remove_nonallowed_values:
                continue
            # otherwise raise exception
            else:
                raise ValueError("Encountered non-allowed value in one of the dictionaries! See docstring for options.")
    return new_dict

def invert_list_to_dict(input_list):
    """ Given a list with no duplicates, return a dict mapping the values to list positions: [a,b,c] -> {a:1,b:2,c:3}."""
    if not len(set(input_list)) == len(input_list):
        raise ValueError("Can't reliably invert a list with duplicate elements!")
    return dict([(value,index) for (index,value) in enumerate(input_list)])

def invert_dict_nodups(input_dict):
    """ Given a dict with no duplicate values, return value:key dict: {a:1,b:2]} -> {1:a,2:b}."""
    if not len(set(input_dict.values())) == len(input_dict.values()):
        raise ValueError("This method can't invert a dictionary with duplicate values! Use invert_dict_tolists.")
    return dict([(value,key) for (key,value) in input_dict.iteritems()])

def invert_dict_tolists(input_dict):
    """ Given a dict (duplicate values allowed), return value:key_list dict: {a:1,b:2,c:1]} -> {1:[a,c],2:b}."""
    inverted_dict = defaultdict(lambda: set())
    for key,value in input_dict.iteritems():
        inverted_dict[value].add(key)
    return dict(inverted_dict)      # changing defaultdict to plain dict to avoid surprises

def invert_listdict_nodups(input_dict):
    """ Given a dict with non-overlapping list/set values, return single_value:key dict: {a:[1,2],b:[3]} -> {1:a,2:a,3:b}.
    """
    inverted_dict = {}
    for key,value_list in input_dict.iteritems():
        try:
            for value in value_list:
                if value in inverted_dict:
                    raise ValueError("This method can't invert a dictionary with duplicate values! "
                                     +"Use invert_listdict_tolists.")
                inverted_dict[value] = key
        except TypeError:
            raise ValueError("invert_listdict_nodups expects all input_dict values to be lists/sets/etc!")
    return inverted_dict

def invert_listdict_tolists(input_dict):
    """ Given a dict with list/set values, return single_value:key_list dict: {a:[1,2],b:[2]]} -> {1:[a],2:[a,b]}."""
    inverted_dict = defaultdict(lambda: set())
    for key,value_list in input_dict.iteritems():
        try:
            for value in value_list:
                inverted_dict[value].add(key)
        except TypeError:
            raise ValueError("invert_listdict_tolists expects all input_dict values to be lists/sets/etc!")
    return dict(inverted_dict)      # changing defaultdict to plain dict to avoid surprises

# MAYBE-TODO refactor to avoid code duplication between the *_tolists and *_nodups pairs above?

def sort_lists_inside_dict(input_dict, reverse=False, key=None):
    """ Given a key:value_list dict, return same dict but with sorted value_lists (using key and reverse args for sort)."""
    new_dict = {}
    for k,l in input_dict.iteritems():
        if key is None:     new_dict[k] = sorted(l, reverse=reverse)
        else:               new_dict[k] = sorted(l, reverse=reverse, key=key)
    return new_dict

def flatten_lists(input_val, unique_only=False):
    """ Given a list or dict of lists/arrays/sets/etc, return list combining all the value lists (or a set if unique_only). """
    if isinstance(input_val, dict):
        input_val = input_val.values()
    all_vals = sum([list(x) for x in input_val], [])
    if unique_only: return set(all_vals)
    else:           return all_vals
    # TODO unit-test!

class keybased_defaultdict(defaultdict):
    """ A defaultdict equivalent that passes the key as an argument to the default-value factory, instead of no argumnets.
    Normal defaultdict(int)[9] is 0, because no argument is passed to the int function, and int() is 0.  
    On the other hand keybased_defaultdict(int)[9] would be 9, keybased_defaultdict(bool)[9] would be True, etc.  """
    # from http://stackoverflow.com/questions/2912231/is-there-a-clever-way-to-pass-the-key-to-defaultdicts-default-factory
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            value = self[key] = self.default_factory(key)
            return value

### Useful class mix-ins

class FrozenClass(object):
    """ Class that allows prevention of adding new attributes at some point after creation - NOT full immutability.

    Use by inheriting from this, and then running self._freeze() at the end of __init__, like this:
        class Test_freezing(FrozenClass):
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self._freeze() # no new attributes after this point.
        a = Test_freezing(2,3)
        a.x = 10    # can modify existing attributes after creation
        a.z = 10    # fails - cannot add NEW atributes after creation
    """
    # by Jochen Ritzel on Stackoverflow 
    #   (http://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    # This has to be a new-style class (i.e. inheriting from object), old-style classes don't seem to have __setattr__?
    # MAYBE-TODO could/should this be done as a decorators instead?

    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)
    
    def _freeze(self):
        self.__isfrozen = True


def reversibly_immutable(wrapped_class):
    """ NOT FINISHED!!!
    
    DECORATOR: gives a class a make_immutable() and make_mutable() method; also makes hashable if _hash() is defined.

    _____
    """
    # MAYBE-TODO not finished!!!  see ~/experiments/mutant_pool_screens/mutant_deepseq_analysis/code/mutant_analysis_classes.py Insertion_position class for one implementation of this - I think it's too class-specific for a general decorator to be easy.

    if not isinstance(object,wrapped_class):
        raise ValueError("Can only apply the reversibly_immutable decorator to new-style classes based on object!")

    def changes_not_allowed(self, *args, **kwargs):
        raise AttributeError("'%s' object is currently immutable, can't make changes!"%type(self))

    # TODO what methods besides __setitem__/__delitem__ do we need to override?  Probably depends on class... COMPLICATED!!
    # MAYBE-TODO this overrides the general methods like __setitem__ etc, but what if the class has other custom methods that modify things? 
    def make_immutable(self):
        # TODO implement removing mutability methods!
        raise NotImplementedError()
        self.__setitem__ = changes_not_allowed
        self.__delitem__ = changes_not_allowed
        # make hashable by setting __hash__ to the private _hash method
        if hasattr(self,'_hash'):
            self.__hash__ = self._hash

    def make_mutable(self):
        # TODO implement re-adding mutability methods!
        raise NotImplementedError()
        # since class is mutable now, it can't be hashable - remove __hash__, 
        #  but keep private _hash method in case we want to make it immutable/hashable again
        if hasattr(self,'__hash__'):
            del self.__hash__

    wrapped_class.make_immutable = make_immutable
    wrapped_class.make_mutable = make_mutable
    return wrapped_class

    
######################################## FILE READING/WRITING, COMMAND-LINE STUFF  ########################################

### Read various file types
# MAYBE-TODO Or I could just use numpy.genfromtext or such for all this...
#    data = numpy.genfromtxt("test.csv",names=True,delimiter=",",dtype=None)

def read_two_column_file(filename,numerical_values=True):
    """ Read in a two-column (name,value) tab-separated file (ignore #-comment lines), return data:float(value) dict. """
    data_dict = {}
    for line in open(filename):
        if line[0]=='#':    continue
        name,value = line.strip().split('\t')
        if numerical_values:    data_dict[name] = float(value)
        else:                   data_dict[name] = value
    return data_dict

def read_tab_separated_file(filename, ignore_comments=True, separator='\t'):
    """ Read in a tab-separated file (ignore #-comment lines), return a list for each column. """
    for line in open(filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip('\n').split(separator)
        try:
            for i in range(len(fields)): data_list[i].append(fields[i])
        except NameError:
            data_list = [[x] for x in fields]
        if not len(fields)==len(data_list):
            sys.exit("Error: Found line with %s columns, but first line had %s columns! (line %s)"%(len(fields), 
                                                                                                    len(data_list), line))
    return data_list
    # TODO add to unit-tests? Or some kind of test.

def read_tab_separated_file_with_headers(filename, ID_column=0, ignore_comments=True, all_values_are_numbers=False):
    """ Read a tab-separated file with column headers in the first line; ignore comment lines (starting with #) unless specified otherwise. Return a list of IDs from ID_column (in order), and a column_header:column_data dictionary, with one entry per non-ID column, where column_data is an ID:value dictionary with values taken from the column with that header."""
    list_of_column_data = read_tab_separated_file(filename,ignore_comments)
    ID_list = list_of_column_data.pop(ID_column)    # grab the ID column and remove it from data column list
    ID_list.pop(0)                                  # discard header from ID column
    data_dict_by_header = {}
    for column_data in list_of_column_data:
        # grab header, put the rest in a dictionary by ID, in order
        colheader = column_data.pop(0)
        if all_values_are_numbers:
            data_dict_by_header[colheader] = dict([(ID_list[i],float(column_data[i])) for i in range(len(column_data))])
        else:
            data_dict_by_header[colheader] = dict([(ID_list[i],column_data[i]) for i in range(len(column_data))])
    return ID_list, data_dict_by_header
    # TODO add to unit-tests? Or some kind of test.


def pickle(data, outfile_name, protocol=0):
    """ Run pickle.dump to save data to the outfile - small convenience function. 

    Pickle protocols: 0 is the backward-compatible version, 1 and 2 are faster but generate binary files.
    Protocol -1 will choose the highest available.
    """
    with open(outfile_name, 'w' if protocol==0 else 'wb') as PICKLEFILE:
        Pickle.dump(data, PICKLEFILE, protocol)


def unpickle(infile_name):
    """ Just run pickle.load on the infile and return the result - small convenience function. """
    with open(infile_name) as PICKLEFILE:
        return Pickle.load(PICKLEFILE)


### Writing to files

class Fake_outfile(object):
    """ Fake file opened for writing - equivalent to writing to /dev/null without touching the OS. All methods do nothing. """
    # see http://stackoverflow.com/questions/2929899/cross-platform-dev-null-in-python
    def write(self, *_):        pass
    def writelines(self, *_):   pass
    def flush(self, *_):        pass
    def seek(self, *_):         pass
    def truncate(self, *_):     pass
    def close(self, *_):        pass
FAKE_OUTFILE = Fake_outfile()
# This is an ALREADY OPEN file.  To get a filename that will yield the equivalent of /dev/null when opened, use os.devnull; not sure if there's a way of doing that without involving the OS, at least without redefining the open builtin.

def replaces_infile_with_outfile(function_to_be_decorated):
    """ DECORATOR, takes a function that takes an infile and outfile and makes it replace infile with outfile instead. 
    The decorated function must still have an argument named outfile, but its value will never be used."""
    def wrapped_function(infile, *args, **kwargs):
        outfile = '.__tmp__'+infile
        kwargs['outfile'] = outfile
        return_val = function_to_be_decorated(infile, *args, **kwargs)
        if os.path.exists(infile):
            os.remove(infile)
        os.rename(outfile,infile)
        return return_val
    return wrapped_function
    # MAYBE-TODO is the way I'm doing this decorator really the best?  Not bad, but requires adding the fake "outfile" argument to all the functions... If I was somehow the function_to_be_decorated's variable list directly instead (with f.__dict__ or f.__setattr__), that wouldn't be an issue... Ask on stackoverflow?
    # TODO add to unit-tests? Or some kind of test.

def save_line_list_as_file(line_list, filename, header="", add_newlines=True):
    """ Given a list of lines, a filename, and an optional header, open file, write header and all lines, close. """
    line_end = "\n" if add_newlines else ""
    # using the with-as syntax automatically closes the file afterward
    with open(filename,'w') as OUTFILE:
        if header:               OUTFILE.write(header+line_end)
        for line in line_list:   OUTFILE.write(line+line_end)

def write_header_data(OUTFILE,options=None):
    """ Print general script run data (command, path, date/time, full optparse options, etc) to given open file object."""
    import pwd,time,socket
    OUTFILE.write("# Command line this file was generated with: %s\n"%(' '.join(sys.argv)))
    OUTFILE.write("# Path: %s\n"%(os.getcwd()))
    OUTFILE.write("# Date: %s,\t\tUser: %s,\t\tSystem: %s\n"%(time.ctime(), pwd.getpwuid(os.getuid())[0], 
                                                              socket.gethostname()))
    if options:     OUTFILE.write("# Full options: %s\n"%options)

def print_text_from_file(infile, OUTFILE=None, printing=True, add_newlines=0):
    """ Write all text from infile to OUTFILE (if not None), also print to stdout if printing is set. 
    Return line counts.  Infile should be a filename; OUTFILE should be an open file object. """
    line_count = 0
    for line in open(infile):
        if OUTFILE is not None:     OUTFILE.write(line)
        if printing:                print line,
        line_count += 1
    if add_newlines:
        if OUTFILE is not None:     OUTFILE.write('\n'*add_newlines)
        if printing:                print '\n'*add_newlines,
    return line_count


### Command/line utilities (running processes, getting output, etc)

def run_command_get_output(command, shell=True):
    """ Run command using subprocess.Popen (with shell arg as given); return (stdout, stderr, returncode). """
    p = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = '' if not p.stdout else p.stdout.read()
    stderr = '' if not p.stderr else p.stderr.read()
    return stdout, stderr, p.returncode


def run_command_and_print_info(command, LOGFILE=None, printing=True, shell=True, program_name=None):
    """ Run command using subprocess.call; first print a line describing that to LOGFILE and/or stdout.
    The shell arg to subprocess.call is given by shell; LOGFILE should be an open file object; 
    program_name is only used for printing, and the first word of the command will be used by default. 
    Note - to print stdout/stderr as well, use run_command_print_info_output. """
    if program_name is None:
        program_name = command.split(' ')[0]
    output = "### Running %s: %s"%(program_name, command)
    if LOGFILE is not None:     LOGFILE.write(output+'\n')
    if printing:                print output
    return run_command_get_output(command, shell)


def run_command_print_info_output(command, LOGFILE=None, printing_level=2, shell=True, program_name=None, add_newlines=0):
    """ Run shell command; write command info line and stdout+stderr to LOGFILE and/or stdout.

    The shell arg to subprocess.Popen is given by shell; LOGFILE should be an open file object; 
    program_name is only used for printing, and the first word of the command will be used by default.
    printing_level 1 means print only the command to stdout; 2 means print command and output.
    The command's return code is returned.
    """
    header_line = "### Running %s: %s\n"%((program_name or command.split(' ')[0]), command) 
    if printing_level>0:        print header_line
    stdout, stderr, retcode = run_command_get_output(command, shell)
    command_output = stdout + stderr
    if not command_output.startswith('\n'):  header_line += '\n'
    full_output = command_output + '\n'*add_newlines
    if printing_level>1:        print full_output
    if LOGFILE is not None:     LOGFILE.write(header_line + full_output)
    return retcode


######################################## NUMERIC DATA MANIPULATION ########################################

NaN = float('NaN')

def nan_func(): 
    """ Simply return the float NaN value (for use in defaultdict). """
    return NaN


def int_or_float(x):
    """ Return int(x) if it's equal to x, else return x (useful to convert whole-number float to int). """
    if int(x) == x:     return int(x)
    else:               return x


def value_and_percentages(val, totals, fractions_not_percentages=False, percentage_format_str='%.2g', value_format_str='%s',
                          NA_for_zero_division=True, exception_for_100='default', insert_word=None, words_for_percentages=None):
    """ Return a string containing val and val as percentage of each total: 1,[2,3] -> "1 (50%, 33%)".

    If fractions_not_percentages=True, 1,[2,3] -> "1 (0.5, 0.33)" instead.
    percentage_format_str is the string used to format each percentage/fraction (in .X, X is the precision in digits).
    If NA_for_zero_division is True, just print N/A for dividing-by-zero cases rather than raising an exception.
    If exception_for_100 is True, 100 is formatted as '100' rather than '1e+02' even if precision is 2.
    If exception_for_100 is 'default', it's True for percentages with '%.2g' format but False otherwise and for fractions.
    If insert_word is given, return "x WORD (y%)" instead of just "x (y%)"
    Words_for_percentages must be a list of strings same length as totals (or None):
     if given, return "x (y% WORD1, z% WORD2)" instead of just "x (y%, z%)"
    """ 
    if exception_for_100=='default': 
        exception_for_100 = True if (not fractions_not_percentages and percentage_format_str in ['%.2g','%.2f']) else False
    def _format_number(number):
        if float(percentage_format_str%number)==100 and exception_for_100:  return '100'
        else:                                                               return percentage_format_str%number
    if fractions_not_percentages:   percentage_getter = lambda x,total: _format_number(x/total)
    else:                           percentage_getter = lambda x,total: _format_number(100*x/total) + '%'
    if NA_for_zero_division:    full_percentage_getter = lambda x,total: 'N/A' if total==0 else percentage_getter(x,total)
    else:                       full_percentage_getter = percentage_getter
    string_for_total = value_format_str%val + ('' if insert_word is None else ' '+insert_word)
    if words_for_percentages is None:   words_for_percentages = [None for _ in totals]
    words_for_percentages = ['' if x is None else ' '+x for x in words_for_percentages]
    strings_for_percentages = ["%s%s"%(full_percentage_getter(val,total), word) 
                               for (total, word) in zip(totals, words_for_percentages)]
    return "%s (%s)"%(string_for_total, ', '.join(strings_for_percentages))


### Get rid of NaN/inf numbers singly or in lists/dicts, replace by input

def clean_number(val,replace_NaN,replace_Inf,replace_NegInf,make_positive=False):
    """ Replace NaN/Inf/NegInf value with arguments provided, otherwise return input. Also make val positive if asked. """
    from numpy import isnan,isinf,isneginf
    if isnan(val):            return replace_NaN
    elif isinf(val):          return replace_Inf
    elif isneginf(val):       return replace_NegInf
    elif make_positive and val<=0:  return replace_NegInf
    else:                           return val
    # TODO add to unit-tests

def clean_data(input_data,replace_NaN,replace_Inf=None,replace_NegInf=None,make_positive=False):
    """ Take a list/tuple/set/dict, return same with NaN/Inf/negative values replaced with the respective arguments. """
    from numpy import isnan,isinf,isneginf
    # if list/tuple/set, just treat as a list and transform back at the end
    if type(input_data) in (list,tuple,set):
        input_list = list(input_data)
    # if dictionary, take the key and value lists, clean the value list, then transform back to dict at the end
    elif type(input_data)==dict:
        key_list = input_data.keys()
        input_list = [input_data[x] for x in key_list]
    else:
        raise ValueError("input_data argument must be a list, tuple, set or dictionary")
    # if replace_NegInf/replace_Inf weren't set, set them to lower/higher than current min/max value
    if replace_NegInf==None:    
        if make_positive:   replace_NegInf = min([x for x in input_list if x>0 and not isnan(x)])
        else:               replace_NegInf = min([x for x in input_list if not isneginf(x) and not isnan(x)])
        if replace_NegInf<0:    replace_NegInf *=1.5
        else:           replace_NegInf *= 0.5
    if replace_Inf==None:    replace_Inf = 1.5 * max([x for x in input_list if not isinf(x) and not isnan(x)])
    # use clean_number to make a cleaned value list
    new_list = [clean_number(x,replace_NaN,replace_Inf,replace_NegInf,make_positive) for x in input_list]
    # output - switch back to original input type
    if type(input_data)==list:      return new_list
    elif type(input_data)==tuple:   return tuple(new_list)
    elif type(input_data)==set:     return set(new_list)
    elif type(input_data)==dict:    return dict(zip(key_list,new_list))
    # TODO add to unit-tests

def clean_data_remove(input_data,remove_NaN=True,remove_Inf=True,remove_NegInf=True,remove_negative=False,remove_zero=False):
    """ Take a list/tuple/set/dict, return same with NaN/Inf/NegInfo/negative values removed accordint to arguments. """
    from numpy import isnan,isinf,isneginf
    def bad_value(x):
        if remove_NaN and isnan(x):         return True
        elif remove_Inf and isinf(x):       return True
        elif remove_NegInf and isneginf(x): return True
        elif remove_negative and x<0:       return True
        elif remove_zero and x==0:          return True
        return False
    # if list/tuple/set, just treat as a list and transform back at the end
    if type(input_data) in (list,tuple,set):
        new_data = [x for x in input_data if not bad_value(x)]
        if type(input_data)==list:      return new_data
        elif type(input_data)==tuple:   return tuple(new_data)
        elif type(input_data)==set:     return set(new_data)
    # if dictionary, have to use a different removal method
    elif type(input_data)==dict:
        return dict([(k,x) for (k,x) in input_data.iteritems() if not bad_value(x)])
    else:
        raise ValueError("input_data argument must be a list, tuple, set or dictionary")
    # TODO add to unit-tests


### cutoffs, ranges, sets, etc 

def parse_cutoffs_into_ranges(cutoff_list, if_overlapping):
    """ Given a cutoff_list [0,2,5], return [(0,inf),(2,inf),(5,inf)] if if_overlapping, else [(0,2),(2,5),(5,inf)]. """ 
    if not if_overlapping:
        return [(x,float('inf')) for x in cutoff_list]
    ranges = []
    for i in range(len(cutoff_list)-1):
        ranges.append((cutoff_list[i],cutoff_list[i+1]))
    ranges.append((cutoff_list[-1],float('inf')))
    return ranges
    # TODO add to unit-tests

def get_sets_from_cutoffs(value_dict, value_ranges):
    """ Given a (name:val) dictionary and a list of (min,max) value ranges, return a (range:names_in_range) dict. """
    ranges_to_sets = {}
    for (vmin,vmax) in value_ranges:
        ranges_to_sets[(vmin,vmax)] = [name for (name,val) in value_dict.items() if val>=vmin and val<vmax]
    return ranges_to_sets
    # TODO add to unit-tests

### local min/max finding

def find_local_maxima_by_width(data, N_surrounding_points=1, include_start_end=True):
    """ Return list of local maxima (value,index) tuples; deals with noise by requiring N lower values on each side.

    Go over the data list searching for a point that's higher than N_surrounding_points on both sides;
     return the list of (value,index) tuples for all such points.
    The point of the N_surrounding_points is to filter out noise by excluding local maxima that are just due to small 
     random changes in what without the noise would be a monotonic line.  There are other approaches to noise problems.

    Dealing with runs of identical local maximum values: the program requires N_surrounding_points on both sides 
      of the entire run of identical values, and only returns the first value/index of the identical run.

    The include_start_end argument governs how to deal with the first and last N_surrounding_points of data:
      - if False, these points will be ignore, since they cannot have N_surrounding_points on both sides
      - if True, these points will be used - the N_surrounding_points requirement will be relaxed to allow fewer
                                                surrounding points if they go up against the data start/end.
    """
    # I also have a full command-line program with options that's a wrapper around this function:
    #    find_local_maxima.py in ~/experiments/other_projects/local_maxima_detector_for_Ute

    # I may be padding the data, so make a copy first! And convert to list - this fails on numpy arrays etc.
    data = list(data)
    if include_start_end:
        padding = [min(data)] * N_surrounding_points
        data = padding + data + padding
    else:   
        padding = []
    local_maxima = []
    for i in range(N_surrounding_points, len(data)-N_surrounding_points):
        curr_value = data[i]
        # if there are multiple identical values, ignore any but the first
        if data[i-1] == data[i]:
            continue
        # grab N_surrounding_points values before the current value
        values_before_curr = data[i-N_surrounding_points : i]
        # grab N_surrounding_points values after the current value, IGNORING any identical ones
        N_plateau_values = 0
        for x in data[i+1:]:
            if x != curr_value:
                break
            N_plateau_values += 1
        values_after_curr = data[i+1+N_plateau_values : i+1+N_plateau_values+N_surrounding_points]
        # if the current value is higher than N_surrounding_points on both sides, add it to local_maxima list
        if curr_value > max(values_before_curr) and curr_value > max(values_after_curr):
            # (if the data was padded, subtract the padding size from the index to get index in original data)
            curr_index = i - len(padding)
            local_maxima.append((curr_value, curr_index))
    return local_maxima

# MAYBE-TODO an alternative implementation would be find_local_maxima_by_height, where we'd take any point that was at least Kx higher than the surrounding M or fewer points, even if M was just 2.  That would also filter out noise reasonably well, assuming the noise is small in amplitude.
# MAYBE-TODO for more flexible dealing with noise, there may be something in scipy... See http://stackoverflow.com/questions/1713335/peak-finding-algorithm-for-python-scipy

### moving average and moving median, plus some help functions (anything starting with _ won't be imported)

def _check_input(data,window_size):
    if window_size<1 or window_size>len(data):   raise ValueError, "window size must be between 1 and data length."

def _make_window_indices(data,window_size):
    half_window = int(window_size/2)
    index_values = [i+half_window for i in range(len(data)-window_size+1)]
    return index_values


def moving_average(data,window_size=10,return_indices=False):
    """ Return the moving average of the input data with the given window size. 

    Optionally also returns the window center indices for plotting. 
    Reasonably efficient implementation, doesn't calculate the whole average for each window.
    """
    _check_input(data,window_size)
    first_average = sum(data[:window_size])/window_size
    average_values = [first_average]
    for i in range(len(data)-window_size):
        # for each new average, just take the previous one, subtract the oldest data value and add the new data value
        new_average = average_values[-1] + (data[i+window_size]-data[i])/window_size
        average_values.append(new_average)
    if return_indices:  return average_values, _make_window_indices(data,window_size)
    else:               return average_values 
    # TODO add to unit-tests


def moving_median(data,window_size=10,return_indices=False):
    """ Return the moving median of the input data with the given window size. 

    Optionally also returns the window center indices for plotting. 
    """
    _check_input(data,window_size)
    from numpy import median
    median_values = []
    for i in range(len(data)-window_size+1):
        median_values.append(median(data[i:i+window_size]))
    if return_indices:  return median_values, _make_window_indices(data,window_size)
    else:               return median_values 
    # TODO add to unit-tests


### Other specialized functions 

def split_into_N_sets_by_counts(ID_counts, N):
    """ Given an ID:count dictionary, return a list of sets of IDs with total counts balanced between the sets. """
    # make a sorted (high to low) list of (count,ID) tuples
    counts_IDs = sorted([(count,ID) for (ID,count) in ID_counts.iteritems()], reverse=True)
    output_counts_sets = [[0,set()] for i in range(N)]
    # now go over all IDs, adding an ID (and the corresponding count) to the smallest set each time
    for (count,ID) in counts_IDs:
        output_counts_sets[0][1].add(ID)
        output_counts_sets[0][0] += count
        output_counts_sets.sort()
    return [ID_set for [count,ID_set] in output_counts_sets]

def merge_values_to_unique(value_list, blank_value='NO_BLANK_VALUES', convert_for_set=(lambda x: x), value_name='', 
                           context='values in value_list passed to merge_multi_values'):
    """ Merge value_list to a single value, ignoring blank_value unless all are blank; raise Exception if impossible.

    If value_list is all blank_values, return blank_value;
    if there's one unique non-blank value in value_list (can show up multiple times), return that;
    if there are multiple unique non-blank values, raise ValueError,
     using value_name and context arguments to provide description (they're only used for that).

    If the elements of value_list are mutable and can't be put in a set, 
     provide a convert_for_set function (like tuple for lists).
    """
    value_set = set([convert_for_set(val) for val in value_list])
    converted_blank_value = convert_for_set(blank_value)
    blank_single_set = set([converted_blank_value])
    # we're using convert_for_set to make the values hashable, but we want to return the ORIGINAL values,
    #  not the hashable versions, so keep track with hashable:original dictionary
    unconverted_values = dict([(convert_for_set(val), val) for val in value_list])
    if value_set == blank_single_set:
        return unconverted_values[converted_blank_value]
    else:
        value_set -= blank_single_set
        if len(value_set) == 1:
            return unconverted_values[value_set.pop()]
        else:
            raise ValueError("Different %s have different %s values! "%(context, value_name)
                             +", ".join([str(val) for val in value_set]))



######################################## TESTS FOR THIS FILE ########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__compare_lists_unordered(self):
        from itertools import permutations
        for input_list in [[1,2,3], [True,False,True], ['a','bb',''], [1,True,'str',321314,None]]:
            assert compare_lists_unordered(input_list, input_list) == True
            assert compare_lists_unordered(input_list, input_list*2) == False
            assert compare_lists_unordered(input_list, []) == False
            for permuted_list in permutations(input_list, len(input_list)):
                assert compare_lists_unordered(input_list, permuted_list) == True
                assert compare_lists_unordered(input_list, permuted_list[:-1]) == False
                assert compare_lists_unordered(input_list[:-1], permuted_list) == False
            for permuted_list in permutations(input_list, len(input_list)-1):
                assert compare_lists_unordered(input_list, permuted_list) == False

    def test__add_to_dict_no_replacement(self):
        # you can add new key:value pairs to dictionary, or re-add previous ones if you don't change the value (no-op)
        d = {}
        add_to_dict_no_replacement(d, 1, 2)
        assert d=={1:2}
        add_to_dict_no_replacement(d, 1, 2)
        assert d=={1:2}
        add_to_dict_no_replacement(d, 2, 2)
        assert d=={1:2, 2:2}
        # if you try changing the value for an existing key, you either get an exception, or it just doesn't work 
        #  (with an optional printed error message), depending on raise_exception argument.
        self.assertRaises(ValueError, add_to_dict_no_replacement, d, 1, 300, raise_exception=True)
        assert d=={1:2, 2:2}
        add_to_dict_no_replacement(d, 1, 300, raise_exception=False, print_error=False)
        assert d=={1:2, 2:2}
        # MAYBE-TODO check that the printed error message is correct?


    def test__reduce_dicts_to_overlaps(self):
        d1 = {1:1}
        d2 = {1:2, 2:2}
        d3 = {2:3, 3:3}
        assert reduce_dicts_to_overlaps([d1,d2]) == [{1:1},{1:2}] 
        assert reduce_dicts_to_overlaps([d2,d3]) == [{2:2},{2:3}] 
        assert reduce_dicts_to_overlaps([d1,d3]) == [{},{}] 
        assert reduce_dicts_to_overlaps([d1,d2,d3]) == [{},{},{}] 
        assert reduce_dicts_to_overlaps([]) == [] 

    def test__filter_dict_by_keys(self):
        # edge cases (good_keys empty or containing all the keys)
        for D in ({}, {1:2}, {x:'a' for x in range(100)}):
            assert filter_dict_by_keys(D, set()) == {}
            assert filter_dict_by_keys(D, []) == {}
            assert filter_dict_by_keys(D, D.keys()) == D
            assert filter_dict_by_keys(D, range(200)) == D 
        # more specific tests for a single dict
        assert filter_dict_by_keys({1:2, 3:4}, [1]) == {1:2}
        assert filter_dict_by_keys({1:2, 3:4}, [2]) == {}
        assert filter_dict_by_keys({1:2, 3:4}, [0,1,2]) == {1:2}
        assert filter_dict_by_keys({1:2, 3:4}, [3]) == {3:4}
        for keys in ([1,3], [0,1,2,3,4,5], [1,2,3,'a','b','c',(1,2)], range(100)):
            assert filter_dict_by_keys({1:2, 3:4}, keys) == {1:2, 3:4}
            assert filter_dict_by_keys({}, keys) == {}

    def test__count_list_values(self):
        assert count_list_values([]) == {}
        assert count_list_values([10,12,11]) == {10:1, 12:1, 11:1}
        assert count_list_values([10,12,10]) == {10:2, 12:1}
        assert count_list_values([1]*100) == {1:100}
        assert count_list_values(['a',None,11]) == {'a':1, None:1, 11:1}
        assert count_list_values([None,None]) == {None:2}

    def test__add_dicts_of_ints(self):
        ### Basic int-value-only cases should be the same regardless of options
        for recursive in True,False:
            for allow_nonconflicting in True,False:
                for remove_nonallowed in True,False:
                    kwargs = {'recursive':recursive, 'allow_nonconflicting_values':allow_nonconflicting, 
                              'remove_nonallowed_values':remove_nonallowed}
                    for dict2 in [{}, {1:1}, {'a':1, 'b':100, 'c':-200}]:
                        assert add_dicts_of_ints({}, dict2, **kwargs) == dict2
                    assert add_dicts_of_ints({1:1}, {1:2}, **kwargs) == {1:3}
                    assert add_dicts_of_ints({1.00:1}, {1.00:2}, **kwargs) == {1.00:3}
                    assert add_dicts_of_ints({1:1}, {2:2}, **kwargs) == {1:1, 2:2}
                    assert add_dicts_of_ints({'one':1}, {'two':2}, **kwargs) == {'one':1, 'two':2}
                    assert add_dicts_of_ints({1:1, 2:1}, {2:2, 3:2}, **kwargs) == {1:1, 2:3, 3:2}
                    assert add_dicts_of_ints({(1,):1, (2,):1}, {(2,):2, (3,):2}, **kwargs) == {(1,):1, (2,):3, (3,):2}
        ### Simplest version - no recursion, no non-int values allowed: all "weird" cases raise exceptions 
        kwargs = {'recursive':False, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':False}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:bad_value}, **kwargs)
        # if nonallowed values are removed, exceptions aren't raised, the bad values are just ignored
        kwargs = {'recursive':False, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':True}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            # if there's only a bad value, or two of them, the result is empty
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {}
            # same if there's also a good value for a different key - only that is kept
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {2:1}
            # if there's a bad and good value for the same key, only the good one is kept
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        ### Allowing recursion: 
        kwargs = {'recursive':True, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':False}
        # still fails for non-dict "bad" values
        for bad_value in ['a', [], {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:bad_value}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:1}, **kwargs)
        # if there's a dictionary-type value in one dictionary and nothing in the other, it's kept
        for dict_value in [{}, {1:1}, {'a':1}]:
            assert add_dicts_of_ints({1:dict_value}, {}, **kwargs) == {1:dict_value}
        # if there's a dictionary-type value in both dictionaries, recursively add them up
        assert add_dicts_of_ints({1:{}}, {1:{}}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:{1:1}}, {1:{1:1}}, **kwargs) == {1:{1:2}}
        assert add_dicts_of_ints({1:{'a':1}}, {1:{'a':1}}, **kwargs) == {1:{'a':2}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {}, **kwargs) == {1:1, 2:{2:2, 3:{3:3}}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:{}}, **kwargs) == {1:2, 2:{2:2, 3:{3:3}}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:{4:4}}, **kwargs) == {1:2, 2:{2:2, 3:{3:3}, 4:4}}
        # if there's a dictionary-type value in one dictionary and an int in the other, raise an exception
        for dict_value in [{}, {1:1}, {'a':1}]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:dict_value}, {1:1}, **kwargs)
        self.assertRaises(ValueError, add_dicts_of_ints, {1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:1}, **kwargs)
        ## if nonallowed values are removed, exceptions aren't raised: the bad values are just ignored
        kwargs = {'recursive':True, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':True}
        # for bad values that don't involve dictionaries, the outcomes are the same as in the non-recursive version
        for bad_value in ['a', [], None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {2:1}
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        # if the two values are mismatched but they're both "good" (like an int and a dict), both are ignored
        for dict_value in [{}, {1:1}, {'a':1}]:
            assert add_dicts_of_ints({1:dict_value}, {1:1}, **kwargs) == {}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:1}, **kwargs) == {1:2}
        ## another special case: if a value is a dictionary but it's weird, it gets recursively passed down:
        #   if remove_nonallowed_values is False, that leads to a ValueError a level later (which was tested above), 
        #   but if it's True, it leads to the dictionary being kept in an empty form (since it contained a bad value)
        weird_value = {1:'a'}
        assert add_dicts_of_ints({1:weird_value}, {}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:weird_value}, {1:weird_value}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:weird_value}, {2:1}, **kwargs) == {1:{}, 2:1}
        assert add_dicts_of_ints({1:weird_value}, {1:1}, **kwargs) == {}
        ### With allow_nonconflicting_values turned on:
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        # this time if the "bad" value is only in one of the dictionaries and nothing in the other, 
        #  or is identical in both dictionaries, it's kept. 
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {1:bad_value, 2:1}
        # but if there's a "good" value in one dictionary and a "bad" one in the other, an exception is still raised.
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:1}, **kwargs)
            # except unfortunately this doesn't work right for True, because 1==True... TODO try to fix that somehow?
        ## if nonallowed values are removed, the first cases are the same, in the last case the good value is kept
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':True}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {1:bad_value, 2:1}
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        ### Recursive takes priority over allow_nonconflicting_values if both apply:
        # if recursive is False, if both values are identical dictionaries, they always get passed on unchanged.
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        assert add_dicts_of_ints({1:{1:2}}, {1:{1:2}}, **kwargs) == {1:{1:2}}
        assert add_dicts_of_ints({1:{1:'a'}}, {1:{1:'a'}}, **kwargs) == {1:{1:'a'}}
        # if recursive is True, the same still happens with bad-valued dictionaries, 
        #  but for int-valued dictionaries the new value is the sum of the originals instead of a copy
        kwargs = {'recursive':True, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        assert add_dicts_of_ints({1:{1:2}}, {1:{1:2}}, **kwargs) == {1:{1:4}}
        assert add_dicts_of_ints({1:{1:'a'}}, {1:{1:'a'}}, **kwargs) == {1:{1:'a'}}

    def test__invert_list_to_dict(self):
        assert invert_list_to_dict([]) == {}
        assert invert_list_to_dict([10,12,11]) == {10:0, 12:1, 11:2}
        assert invert_list_to_dict(['a',None,11]) == {'a':0, None:1, 11:2}
        self.assertRaises(ValueError, invert_list_to_dict, [10,11,10])
        self.assertRaises(ValueError, invert_list_to_dict, [None,None])

    def test__invert_dict_nodups(self):
        assert invert_dict_nodups({}) == {}
        assert invert_dict_nodups({1:2,3:4}) == {2:1,4:3}
        self.assertRaises(ValueError, invert_dict_nodups, {1:2,3:2})

    def test__invert_dict_tolists(self):
        assert invert_dict_tolists({}) == {}
        assert invert_dict_tolists({1:2,3:4}) == {2:set([1]),4:set([3])}
        assert invert_dict_tolists({1:2,3:2}) == {2:set([1,3])}

    def test__invert_listdict_nodups(self):
        assert invert_listdict_nodups({}) == {}
        assert invert_listdict_nodups({1:[2],3:[4]}) == {2:1,4:3}
        assert invert_listdict_nodups({1:[2],3:[4,6]}) == {2:1,4:3,6:3}
        assert invert_listdict_nodups({1:[2,8],3:[4,6]}) == {2:1,8:1,4:3,6:3}
        # no duplicates
        self.assertRaises(ValueError, invert_listdict_nodups, {1:[2,2],3:[4]})
        self.assertRaises(ValueError, invert_listdict_nodups, {1:[2,4],3:[2]})
        # values must be lists/sets
        self.assertRaises(ValueError, invert_listdict_nodups, {1:2,3:2})

    def test__invert_listdict_tolists(self):
        assert invert_listdict_tolists({}) == {}
        assert invert_listdict_tolists({1:[2],3:[4]}) == {2:set([1]),4:set([3])}
        assert invert_listdict_tolists({1:[2],3:[2]}) == {2:set([1,3])}
        assert invert_listdict_tolists({1:[2,4],3:[2]}) == {2:set([1,3]),4:set([1])}
        self.assertRaises(ValueError, invert_listdict_tolists, {1:2,3:2})

    def test__sort_lists_inside_dict(self):
        assert sort_lists_inside_dict({}) == {}
        self.assertRaises(TypeError, sort_lists_inside_dict, {1:1, 2:2})
        self.assertRaises(TypeError, sort_lists_inside_dict, {1:[1,2,3], 2:2})
        assert sort_lists_inside_dict({1:[1,3,2], 2:[4,5,4]}) == {1:[1,2,3], 2:[4,4,5]}
        assert sort_lists_inside_dict({1:[1,3,2], 2:[4,5,4]}, reverse=True) == {1:[3,2,1], 2:[5,4,4]}
        assert sort_lists_inside_dict({1:[0,-1,2]}) == {1:[-1,0,2]}
        assert sort_lists_inside_dict({1:[-1,0,2]}, reverse=True) ==  {1:[2,0,-1]}
        assert sort_lists_inside_dict({1:[-1,0,2]}, key=lambda x: abs(x)) ==  {1:[0,-1,2]}
        assert sort_lists_inside_dict({1:[-1,0,2]}, key=lambda x: abs(x), reverse=True) ==  {1:[2,-1,0]}

    def test__keybased_defaultdict(self):
        D_nodefault = keybased_defaultdict(None)
        self.assertRaises(KeyError, lambda: D_nodefault[1])
        D_constantdefault = keybased_defaultdict(lambda x: 0)
        assert D_constantdefault[1] == 0
        assert D_constantdefault[2] == 0
        D_variabledefault = keybased_defaultdict(lambda x: 2*x)
        assert D_variabledefault[1] == 2
        assert D_variabledefault[2] == 4

    def test__FrozenClass(self):
        class Test_freezing(FrozenClass):
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self._freeze() # no new attributes after this point.
        a = Test_freezing(2,3)
        assert a.x == 2
        a.x = 10    # can modify existing attributes after creation - shouldn't raise an exception
        assert a.x == 10
        # doing a.z = 10 should fail - cannot add NEW atributes after creation
        #  testing this two ways: with __setattr__ as a function, and with writing a test function to explicitly 
        #  test the "a.z = 1" statement (which of course should use __setattr__ anyway, but may as well check)
        self.assertRaises(TypeError, a.__setattr__, 'z', 10)
        def test_function(obj,val):  
            obj.z = val
        self.assertRaises(TypeError, test_function, a, 10)

    def test__FAKE_OUTFILE(self):
        assert FAKE_OUTFILE.write() is None
        assert FAKE_OUTFILE.writelines() is None
        assert FAKE_OUTFILE.flush() is None
        assert FAKE_OUTFILE.seek() is None
        assert FAKE_OUTFILE.truncate() is None
        assert FAKE_OUTFILE.close() is None

    def test__int_or_float(self):
        for bad_input in ['a', [], [1,2,3]]:
            self.assertRaises((ValueError,TypeError), int_or_float, bad_input)
        assert int_or_float(3) == 3
        assert int_or_float(3.0) == 3
        assert int_or_float(3.5) == 3.5

    def test__value_and_percentages(self):
        assert value_and_percentages(1, [2], False) == "1 (50%)"
        assert value_and_percentages(1, [2], True) == "1 (0.5)"
        assert value_and_percentages(1, [2, 3, 100, 10000], False) == "1 (50%, 33%, 1%, 0.01%)"
        assert value_and_percentages(1, [2, 3, 100, 10000], True) == "1 (0.5, 0.33, 0.01, 0.0001)"
        assert value_and_percentages(1, [3, 7000], False, percentage_format_str='%.2g') == "1 (33%, 0.014%)"
        assert value_and_percentages(1, [3, 7000], False, percentage_format_str='%.4g') == "1 (33.33%, 0.01429%)"
        assert value_and_percentages(1, [3, 7000], True, percentage_format_str='%.2g') == "1 (0.33, 0.00014)"
        assert value_and_percentages(1, [3, 7000], True, percentage_format_str='%.4g') == "1 (0.3333, 0.0001429)"
        # NA_for_zero_division - if True, just print N/A for division-by-zero rather than raising an exception
        for zero in (0, 0.0):
            self.assertRaises(ZeroDivisionError, value_and_percentages, 1, [zero], NA_for_zero_division=False)
            assert value_and_percentages(1, [zero], NA_for_zero_division=True) == "1 (N/A)"
        # exception_for_100 (default True for percentages with '%.2g' and False for fractions/otherwise)
        assert value_and_percentages(1, [1], False, '%.2g', exception_for_100='default') == "1 (100%)"
        assert value_and_percentages(1, [1], False, '%.2g', exception_for_100=True) == "1 (100%)"
        assert value_and_percentages(1, [1], False, '%.2g', exception_for_100=False) == "1 (1e+02%)"
        assert value_and_percentages(1, [1], False, '%.1g', exception_for_100='default') == "1 (1e+02%)"
        assert value_and_percentages(1, [1], True) == "1 (1)"
        assert value_and_percentages(1, [0.01], True, '%.2g', exception_for_100='default') == "1 (1e+02)"
        assert value_and_percentages(1, [0.01], True, '%.2g', exception_for_100=False) == "1 (1e+02)"
        assert value_and_percentages(1, [0.01], True, '%.2g', exception_for_100=True) == "1 (100)"
        # testing insert_word option
        assert value_and_percentages(1, [2], False, insert_word='A') == "1 A (50%)"
        assert value_and_percentages(1, [2], True, insert_word='turtle(s)') == "1 turtle(s) (0.5)"
        # testing words_for_percentages option
        assert value_and_percentages(1, [2,4], False, words_for_percentages='A B'.split()) == "1 (50% A, 25% B)"
        assert value_and_percentages(1, [2,4], True, words_for_percentages=['of this', 'of that']) == "1 (0.5 of this, 0.25 of that)"
        # testing value_format_str option
        assert value_and_percentages(3.3, [6.6], False, value_format_str='%s') == "3.3 (50%)"
        assert value_and_percentages(3.3, [6.6], False, value_format_str='%.0f') == "3 (50%)"
        assert value_and_percentages(3.3, [6.6], False, value_format_str='%.4f') == "3.3000 (50%)"

    def test__find_local_maxima_by_width(self):
        # basic functionality - find the local maximum
        for N_surrounding_points in [1,2,3]:
            assert find_local_maxima_by_width([1,2,3,4,3,2,1], N_surrounding_points) == [(4,3)]
            assert find_local_maxima_by_width([1,2,3,4,3,2,1,2,3,4,3,2,1], N_surrounding_points) == [(4,3),(4,9)]
        # require at least N_surrounding_points lower points on both sides of the maximum
        assert find_local_maxima_by_width([1,2,1,2,1,2,1], N_surrounding_points=1) == [(2,1),(2,3),(2,5)]
        for N_surrounding_points in [2,3]:
            assert find_local_maxima_by_width([1,2,1,2,1,2,1], N_surrounding_points) == []
        # require N_surrounding_points on either end 
        assert find_local_maxima_by_width([1,1,2,1,2], N_surrounding_points=2) == []
        assert find_local_maxima_by_width([1,1,2,1,2], N_surrounding_points=1) == [(2,2),(2,4)]
        assert find_local_maxima_by_width([2,1,2,1,1], N_surrounding_points=2) == []
        assert find_local_maxima_by_width([2,1,2,1,1], N_surrounding_points=1) == [(2,0),(2,2)]
        # how to deal with points in the first/last N_surrounding_points depends on include_start_end:
        assert find_local_maxima_by_width([2,1,1,2], N_surrounding_points=1, include_start_end=False) == []
        assert find_local_maxima_by_width([2,1,1,2], N_surrounding_points=1, include_start_end=True) == [(2,0),(2,3)]
        # when there are adjacent points with same value, return either one
        assert find_local_maxima_by_width([1,2,2,1],1) in ( [(2,1)], [(2,2)] )

    def test__split_into_N_sets_by_counts(self):
        input1 = {'a':1000}
        for N in range(1,10):
            assert compare_lists_unordered(split_into_N_sets_by_counts(input1,N), 
                                           [set(['a'])] + [set() for i in range(N-1)])

        input2 = {'a':1002, 'b':1001, 'c':1000}
        assert compare_lists_unordered(split_into_N_sets_by_counts(input2,1), [set(['a','b','c'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input2,2), [set(['a']),set(['b','c'])])
        for N in range(3,10):
            assert compare_lists_unordered(split_into_N_sets_by_counts(input2,N), 
                                           [set(['a']), set(['b']), set(['c'])] + [set() for i in range(N-3)])

        input3 = {'a':5, 'b':4, 'c':3, 'd':2, 'e':1}
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,1), [set(['a','b','c','d','e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,2), [set(['a','d','e']), set(['b','c'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,3), [set(['a']), set(['b','e']), set(['c','d'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,4), 
                                       [set(['a']), set(['b']), set(['c']), set(['d','e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,5), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,6), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set()])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,7), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set(), set()])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,8), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set(), set(), set()])

    def test__merge_values_to_unique(self):
        # with no blank, just return a value if input is multiple copies of it, otherwise raise ValueError
        assert merge_values_to_unique('aaaa') == 'a'
        assert merge_values_to_unique([1,1,1]) == 1
        self.assertRaises(ValueError, merge_values_to_unique, [1,2,1])
        # with a blank, ignore blank values UNLESS they're the only ones present
        assert merge_values_to_unique([1,2,1], blank_value=1) == 2
        assert merge_values_to_unique([1,1,1], blank_value=1) == 1
        assert merge_values_to_unique([1,2,1], blank_value=2) == 1
        assert merge_values_to_unique([2,2,2], blank_value=2) == 2
        self.assertRaises(ValueError, merge_values_to_unique, [1,2,1], blank_value=3)
        self.assertRaises(ValueError, merge_values_to_unique, [1,2,3], blank_value=1)
        self.assertRaises(ValueError, merge_values_to_unique, [1,2,3], blank_value=2)
        # testing the convert_for_set function - unhashable value types should return TypeError without it
        self.assertRaises(TypeError, merge_values_to_unique, [[1,1], [1,1]])
        assert merge_values_to_unique([[1,1], [1,1]], convert_for_set=tuple) == [1,1]
        self.assertRaises(ValueError, merge_values_to_unique, [[1,1], [2,2]], convert_for_set=tuple)
        assert merge_values_to_unique([[1,1], [2,2]], blank_value=[1,1], convert_for_set=tuple) == [2,2]
        # MAYBE-TODO add tests for the value_name and context args, which are only used for the error message?

    # TODO add tests for everything else


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
