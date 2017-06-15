#! /usr/bin/env python2
"""
Various utilities for dealing with binary codes: representation, reading from a file, calculating the bit sum, Hamming distances, comparing values within sets, etc.  See individual function docstrings for details.
"""

import sys
from collections import defaultdict
import itertools
from numpy import array, dot
import random
import unittest
import bitstring

# test bitstring module BitArray.bin behavior!  It changed at some point between early 2012 and early 2013...
_b = bitstring.BitArray(bin='111').bin
if _b == '111':         _bitstring_behavior = 'NEW'
elif _b == '0b111':     _bitstring_behavior = 'OLD'
else:                     
    sys.exit("Unexpected bitstring module behavior!  bitstring.BitArray(bin='111') should be '111' or '0b111', is %s!"%_bin)

# two functions to get binary string from BitArray, depending on bitstring module version
def _get_Bc_string_old(binary_codeword):    return binary_codeword.codeword.bin[2:]
def _get_Bc_string_new(binary_codeword):    return binary_codeword.codeword.bin

# my modules
from general_utilities import invert_dict_tolists, invert_listdict_tolists

class BinaryCodeError(Exception):
    """ Exceptions in the binary_code_utilities module."""
    pass


######### Binary string representations

class Binary_codeword:
    """ A binary string representation (like '01101'), with a defined length. Supports |, &, ^, ~ bitwise operators.

    Can be made from string, int, list.  
    Not just a binary representation of an integer: '001' and '01' are distinct. 

    Not actually represented as strings internally - that would be insanely slow.
    Implemented using the bistring package: 
        http://pypi.python.org/pypi/bitstring/2.2.0 , http://code.google.com/p/python-bitstring/ """
    #  see http://stackoverflow.com/questions/142812/does-python-have-a-bitfield-type for more implementation options

    def __init__(self,val,length=None,check_length=False):
        """ Generate the self.codeword binary word based on val; pad with 0s on the left to desired length if specified.
        If check_length is True, instead of padding make sure the length is as specified, raise BinaryCodeError if not.
        If val is a 0/1 string, strip spaces/newlines and convert straight to a binary string.
        If val is an int, act as if the builtin bin() function was used to convert to a string (stripping initial 0b).
        If val is a list of 0/1 or True/False values, treat it as the corresponding 0/1 string.
        If val is a Binary_codeword instance, just copy it.
        How other argument types will behave is not guaranteed - depends on the package used for the representation."""
        # for binary strings just specify it's bin: '110' and '0b110' and '  110\n' and '11 0' all give the same result
        if isinstance(val,str):   
            self.codeword = bitstring.BitArray(bin=val)
        # for ints, use the uint method if we know the length, otherwise have to convert to string
        elif isinstance(val,int):
            if length is not None:  self.codeword = bitstring.BitArray(uint=val,length=length)
            else:                   self.codeword = bitstring.BitArray(bin(val))
        # for Binary_codeword instances, make a new one with the same bitstring
        elif isinstance(val,Binary_codeword):
            self.codeword = bitstring.BitArray(val.codeword)
        # lists of 0/1 or True/False values natively work as expected; I don't know or care what happens with other types.
        else:                       self.codeword = bitstring.BitArray(val)   
        # pad to given length or check the length if necessary
        if length is not None and not length==len(self):
            # if only checking the length, raise an error
            if check_length:
                raise BinaryCodeError("The created binary codeword didn't match the expected length!")
            # otherwise pad with 0s on the left to the given length (does nothing if length < len(self)
            self.codeword.prepend((length-len(self))*[0])
            # unless it's already too long, then raise an error
            if length-len(self) < 0:
                raise BinaryCodeError("Can't pad the codeword %s to length %s, "%(self.codeword,length)
                                      + "since that's lower than its current length!")

    # TODO this should be a property, really!  Same for a lot of the below.
    def weight(self):
        """ Return the number of 1's in the codeword (i.e. the Hamming weight or bitwise sum)."""
        return self.codeword.count(1)

    # string function depends on bitstring module version/behavior (it changed on me at some point)
    #  (I'm doing it this way instead of putting a test inside the string function for speed)
    if _bitstring_behavior == 'NEW':        string = _get_Bc_string_new
    elif _bitstring_behavior == 'OLD':      string = _get_Bc_string_old
    string.__doc__ = """ Return a plain 0/1 string representation. """

    def list(self):
        """ Return a representation as a list of ints (0 or 1). """
        return [int(x) for x in list(self.codeword)]

    def __len__(self):
        """ Return the length of the codeword."""
        return self.codeword.length

    # Bitwise and, or, xor, not operators - just passed on to the underlying bitstring objects.
    def __and__(self,other):    return Binary_codeword(self.codeword & other.codeword)
    def __or__(self,other):     return Binary_codeword(self.codeword | other.codeword)
    def __xor__(self,other):    return Binary_codeword(self.codeword ^ other.codeword)
    def __invert__(self):       return Binary_codeword(~self.codeword)

    def __cmp__(self,other):
        """ Comparison/sorting based on string representation: Two instances with the same bitstring should be equal."""
        # NOTE: comparison and hashing are related and need to match!
        #   By default, they're based on identity: hash(x) is id(x), and x==y iff id(x)==id(y), i.e. if 
        #     x and y are the SAME object, not just different objects with identical contents. 
        #   So unless I overload the operators, a set can contain two distinct Binary_code('111') instances! BAD.
        #   If I implement __cmp__ but not __hash__, the objects are considered unhashable, because otherwise
        #     there can be cases where x==y but hash(x)!=hash(y), which is BAD for hashing.  
        #   See http://docs.python.org/reference/datamodel.html#object.__hash__ for more on this.
        # MAYBE-TODO in order to make this absolutely right, I should make sure the bitstrings are immutable... How?
        #   Apparently newer version of bistring has an immutable Bits type - could use that instead of BitArray?
        # MAYBE-TODO is string comparison really what I want here?  How about when the lengths are different? Should 
        #   bitstrings with different lengths even be comparable?  I suppose they should just so I can sort stuff and 
        #   get a consistent result. Possibly just sorting by length first would be better, but it doesn't matter much.
        #   As long as identity works correctly and there's SOME reproducible sorting, we're fine.
        # MAYBE-TODO should I also implement comparison to strings and such, by doing binary_codeword(other).string()?
        #   That conversion should probably be taken care of in the function performing the comparison.
        return cmp(self.string(),other.string())

    def __hash__(self):
        """ Hashing based on the string representation."""
        # NOTE: comparison and hashing are related and need to match!  See notes in __cmp__
        return hash(self.string())

    # This is exactly as it should be: __repr__ gives the full evaluatable definition, __str__ gives something readable.
    def __repr__(self):     return "Binary_codeword('%s')"%self.string()
    def __str__(self):      return self.string()


######### basic functions on multiple Binary_codeword objects: distance, mutations, etc

def Hamming_distance(val1,val2):
    """ Given two binary strings, return the number of bits by which their binary representations differ. """
    bitwise_xor = val1^val2
    return bitwise_xor.weight()


def bit_change_count(val1,val2):
    """ Given two binary strings, return A) the count of bits that were 0 in val1 and 1 in val2, and B) the reverse. """
    count_0_to_1 = (~val1)&val2
    count_1_to_0 = val1&(~val2)
    return count_0_to_1.weight(), count_1_to_0.weight()


def _change_all_position_combinations(val, max_N_changes, allowed_change_positions=None):
    """ Just a help function for expand_by_all_mutations.
    Return a set of all new values generated by inverting (0/1) at most max_N_changes of the bits 
     of the original val argument, with changes only allowed in allowed_change_positions (all by default).
    val must be a list of 0/1 ints. """
    try:                    val[0]
    except TypeError:       raise ValueError("val has to be a list of 0/1 ints!")
    if not val[0] in [0,1]: raise ValueError("val has to be a list of 0/1 ints!")
    # if no allowed_change_positions was given, assume all positions are allowed
    if allowed_change_positions is None:
        allowed_change_positions = range(len(val))
    new_value_set = set()
    for N_changes in range(max_N_changes+1):
        if N_changes>len(allowed_change_positions): break   # this is necessary in 2.6.1 but not 2.6.5...
        for change_positions in itertools.combinations(allowed_change_positions,N_changes):
            new_val = list(val)
            for pos in change_positions:
                new_val[pos] = int(not val[pos])
            new_value_set.add(tuple(new_val))   # need to convert to tuple because lists aren't hashable
    return new_value_set


def expand_by_all_mutations(codeword_set, N_permitted_changes):
    """ Return set of all codewords that are too close to codeword_set (distant by at most the given number of changes).
    N_permitted_changes can be given either as a single value, or as a (N_1_to_0_changes, N_0_to_1_changes) tuple."""
    codewords_as_lists = [codeword.list() for codeword in codeword_set]
    expanded_codeword_set = set()

    if isinstance(N_permitted_changes,int):
        for codeword in codewords_as_lists:
            # change all possible combinations of positions (0s or 1s, all positions are allowed)
            new_val_full_set = _change_all_position_combinations(codeword,N_permitted_changes)
            expanded_codeword_set.update([Binary_codeword(new_val) for new_val in new_val_full_set])

    elif isinstance(N_permitted_changes,tuple) and len(N_permitted_changes)==2:
        (permitted_1_to_0_changes,permitted_0_to_1_changes) = N_permitted_changes
        for codeword in codewords_as_lists:
            # change all possible combinations of 1s to 0s up to the allowed limit
            list_of_1_positions = [x for x in range(len(codeword)) if codeword[x]==1]
            new_val_halfway_set = _change_all_position_combinations(codeword,permitted_1_to_0_changes,list_of_1_positions)
            # now take each of the resulting values and change all possible combinations of 0s up to the allowed limit
            for new_val in new_val_halfway_set:
                list_of_0_positions = [x for x in range(len(new_val)) if new_val[x]==0]
                new_val_full_set = _change_all_position_combinations(new_val, permitted_0_to_1_changes,list_of_0_positions)
                expanded_codeword_set.update([Binary_codeword(new_val) for new_val in new_val_full_set])
    else:
        raise ValueError("N_permitted_changes must be an int or a tuple of two ints!")

    return expanded_codeword_set


def expand_by_all_mutations_dict(codeword_set, N_permitted_changes):
    """ Return a dictionary with the keys being all the codeword too close to codeword_set, 
     and the values being the list of base codewords that that particular result codeword was too close to.
    (See expand_by_all_mutations docstring for details on what "too close" means.) """
    codeword_to_expanded_set = dict()
    for C in codeword_set:
        codeword_to_expanded_set[C] = expand_by_all_mutations([C], N_permitted_changes)
    return invert_listdict_tolists(codeword_to_expanded_set)


######### Binary code (set of binary strings) representation

# MAYBE-TODO figure out a naming that won't confuse people!  (Or me!)  Mathematically a "code" is a set of codewords (or a method of encoding things), but IRL a "code" can be either a method of encoding things or just a string, so code/codeword is confusing and even I'm using them wrong!

class Binary_code:
    """ Essentially a set of Binary_codeword objects, all of the same length."""

    def __init__(self,length,val=[],method='list',expected_count=0):
        """ Initialize with given codeword length; add all elements of values to the set of codewords (default empty)."""
        try: 
            self.length = int(length)
        except (ValueError,TypeError):  
            raise BinaryCodeError('Binary_code length argument "%s" is not an int or possible to cast to an int!'%length)
        self.method = method
        self.codewords = set()
        if method=='list':
            for x in val: self.add(x)
        elif method=='listfile':    
            self.read_code_from_file(val)
        elif method=='matrix':    
            self.get_code_from_generator_matrix(generator_matrix=val)
        elif method=='matrixfile':    
            self.get_code_from_generator_matrix(generator_file=val)
        else:
            raise BinaryCodeError("method %s passed to the Binary_code initializer not recognized! "%method 
                                  + "(allowed methods are list, listfile, matrix, matrixfile)")
        if expected_count and not self.size()==expected_count:
            raise BinaryCodeError("Binary_code initializer gave %s codewords, "%self.size() 
                                  + "not %s as expected!"%expected_count)

    def add(self,val):
        """ Add Binary_code(val) codeword to the code, checking for correct length."""
        # it's all right if val is a Binary_codeword already, that works too - similar to sets
        self.codewords.add(Binary_codeword(val,length=self.length,check_length=True))

    def remove(self,val):
        """ Remove Binary_code(val) codeword from the code; fail if val wasn't in the code, or is the wrong length."""
        try:
            self.codewords.remove(Binary_codeword(val,length=self.length,check_length=True))
        # trying to remove an element from a set where it wasn't present raises KeyError - we want a similar behavior.
        except KeyError:
            raise BinaryCodeError("Codeword %s cannot be removed from code because it wasn't present!"%val)

    def remove_extreme_codeword(self,bit=0):
        """ Remove the all-zero codeword (if bit==0; default) or the all-one codeword (if bit==1) from the code.  
        Return 1 if the codeword was present, otherwise return 0; only raise an exception if bit argument isn't 0/1. """
        # MAYBE-TODO or should I make this return a new code instead? Do I want the codes to be immutable or not?
        if bit not in [0,1,'0','1']:
            raise BinaryCodeError("bit argument to remove_extreme_codeword must be 0 or 1!")
        try:
            self.codewords.remove(Binary_codeword(str(bit)*self.length,length=self.length))
            return 1
        except KeyError:
            return 0

    def size(self):
        """ Return the number of codewords currently in the code."""
        return len(self.codewords)

    def read_code_from_file(self,infile,expected_count=0):
        """ Populate the code with codewords read from a plaintext file of 0/1 strings (one per line).
        Skip comment lines (starting with #).  Optionally make sure the codeword count is as expected. """
        read_count = 0
        for line in open(infile):
            if not line.startswith('#'):
                self.add(line)
                read_count += 1
        if expected_count and not read_count==expected_count:   
            raise BinaryCodeError("File %s contained %s codewords, "%(infile,read_count)
                                  + "not %s as expected!"%expected_count)

    def write_code_to_file(self,outfile):
        """ Write all the codewords (in arbitrary order) to a plaintext file of 0/1 strings (one per line). """
        OUTFILE = open(outfile,'w')
        for codeword in self.codewords:
            OUTFILE.write(codeword.string()+'\n')
        OUTFILE.close()

    def get_code_from_generator_matrix(self,generator_file=None,generator_matrix=None,expected_count=0):
        """ Given either a generator matrix (as a numpy array) or a plaintext file containing the matrix, 
        add all codewords generated by that matrix to the current code.  Check if the count is as expected, if given.
        More information on how a generator matrix works: http://en.wikipedia.org/wiki/Generator_matrix."""
        # make sure we didn't get a file AND an arra; read the matrix from the file if it was provided as a file
        if generator_file and generator_matrix:
            raise BinaryCodeError("Provide either a generator_file or a generator_matrix, not both!")
        if generator_file:
            generator_matrix = array([[int(x) for x in line.strip()] for line in open(generator_file)])
        # the encoded word length and the codeword length (i.e. the n and k from the [n,k,d] description)
        #  can be inferred from the generator matrix size; check if they match self.length and expected_count (if given)
        input_code_length, output_code_length = generator_matrix.shape
        codeword_count = 2**input_code_length
        if not self.length==output_code_length:
            raise BinaryCodeError("Trying to use a generator matrix of shape (%s,%s) "%generator_matrix.shape
                                  + "to generate codewords of length %s - the sizes don't match up!"%self.length)
        if expected_count and not expected_count==codeword_count:
            raise BinaryCodeError("The number of codewords generated by the given matrix will be %s, "%codeword_count
                                  + "not %s as expected!"%expected_count)
        # generate all possible input codewords of length k, and convert them to codewords using matrix dot-multiplication
        #   (need to do a %2 at the end because binary addition in coding theory is basically xor - I think...)
        for x in range(codeword_count):
            x = Binary_codeword(x,length=input_code_length).list()
            self.add(dot(x,generator_matrix)%2)

    # MAYBE-TODO I could make one or both of these be the initialization signature instead, but who cares
    # MAYBE-TODO could add the minimum Hamming distance to this?
    def __str__(self):      return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())
    def __repr__(self):     return "<Binary_code instance of length %s and size %s>"%(self.length,self.size())

    # Implementing eq, ne and hashing based on codeword sets; no ge/le comparison (you can't sort sets)
    # NOTE: comparison and hashing are related and need to match!  See notes in Binary_codeword.
    def __eq__(self,other):     return self.codewords == other.codewords
    def __ne__(self,other):     return self.codewords != other.codewords
    # MAYBE-TODO this hashing solution is dangerous, since technically codeword sets ARE mutable... Hmmmm...
    #def __hash__(self):         return hash(frozenset(self.codewords))      

    def find_Hamming_distance_range(self):
        """ Return a tuple containing the lowest and highest Hamming distance between all codeword pairs."""
        # set start lov/high values to impossible extremes to make sure they get overwritten
        if len(self.codewords)==0:
            return None,None
        low = self.length+1
        high = 0
        for x,y in itertools.combinations(self.codewords,2):
            dist = Hamming_distance(x,y)
            low = min(low,dist)
            high = max(high,dist)
        return low,high
        # more on Hamming distance comparisons/implementations: http://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string

    def find_bit_sum_counts(self):
        """ Return the number of codewords with each possible bit-sum value (weight), as a list of (bit-sum, count) tuples.
        The return value is a sorted list of tuples, for readability, but convertion into a dictionary is trivial."""
        bit_sum_counts = defaultdict(lambda: 0)
        for codeword in self.codewords:
            bit_sum_counts[codeword.weight()] += 1
        return sorted(bit_sum_counts.items())

    def total_bit_sum(self):
        """ Return the total sum of bits in all the codewords."""
        return sum([codeword.weight() for codeword in self.codewords])

    def bit_sums_across_digits(self):
        """ Return a list giving the total number of codewords with a 1 at each digit, over codeword length. """
        return [sum([codeword.list()[digit] for codeword in self.codewords]) for digit in range(self.length)]

    def add_parity_bit(self):
        """ Return a new Binary_code object generated by adding a parity bit to the current codewords.
        The new Binary_code will have the same number of codewords, a length higher by 1, and a minimum Hamming distance
        equal to the current one if the current one is even, or higher by one if the current one is odd."""
        new_code = Binary_code(self.length+1)
        for codeword in self.codewords:
            parity_bit = codeword.weight() % 2
            new_code.add(codeword.string()+str(parity_bit))
        assert new_code.size()==self.size()
        return new_code

    def invert(self):
        """ Return a new Binary_code object containing the bitwise inverses of all the codewords in this code."""
        new_code = Binary_code(self.length)
        for codeword in self.codewords:
            new_code.add((~codeword).string())
        return new_code

    def choose_codewords_by_bit_sum(self, low, high, replace_self=False):
        """ Take all codewords with bit sums in low-high range; either return the set, or replace self.codewords with it.
        If high is -1, don't apply an upper bound. 
        """
        if high==-1:    high = self.find_bit_sum_counts()[-1][0]
        new_codewords = set()
        for codeword in self.codewords:
            if low <= codeword.weight() <= high:  new_codewords.add(codeword)
        if replace_self:
            self.codewords = new_codewords
            return
        else:
            return new_codewords

    def give_N_codewords_random(self,N):
        """ Return a set of N randomly chosen codewords.  Raise an error if N is higher than current code size. """
        if N>self.size():
            raise BinaryCodeError("Cannot reduce the code to %s elements, it's already only %s!"%(N,self.size()))
        # Grab the codewords as a list; Sort normally (i.e. by string) just to make sure the result is reproducible.
        new_codeword_set = set(random.sample(self.codewords, N))
        return new_codeword_set

    def give_N_codewords_by_bit_sum(self, N, take_high=False):
        """ Return set of N codewords with the lowest possible bit-sums (or highest if take_high); random within that.
        First determines the bit-sum range (starting at lowest and going up, or the opposite if take_high) that contains
         at least N codewords; then chooses N codewords randomly from all the codewords in that bit-sum range.
        Raise an error if N is higher than current code size. """
        if N>self.size():
            raise BinaryCodeError("Cannot reduce the code to %s elements, it's already only %s!"%(N,self.size()))
        # get a list of (bit-sum,codeword-count) tuples, sorted low-to-high or high-to-low
        bit_sum_counts = sorted(self.find_bit_sum_counts(), reverse=take_high)
        # go over the list in order, keeping the bit-sums as a list, until there are enough codewords
        bit_sums_to_use, codeword_total = [], 0
        for bit_sum, codeword_count in bit_sum_counts:
            bit_sums_to_use.append(bit_sum)
            codeword_total += codeword_count
            if codeword_total>=N:   break
        # now that we know what bit-sum range to use, grab all the keywords from that range, and randomly choose N.
        bit_sum_min, bit_sum_max = min(bit_sums_to_use), max(bit_sums_to_use)
        codewords_by_bit_sum = self.choose_codewords_by_bit_sum(bit_sum_min, bit_sum_max, replace_self=False)
        new_codeword_set = set(random.sample(codewords_by_bit_sum, N))
        return new_codeword_set
        # MAYBE-TODO add an option to make it take all the codewords from all the bit-sum sets except the last one (and random from the last one to get up to N), instead of just taking random N ones from the whole range?  More complicated, not sure if useful.

    def give_N_codewords_even_distribution(self, N, N_tries, return_repeat_summary=False):
        """ Run give_N_codewords_random N_tries times, return result with most even bit_sums_across_digits distribution.
        If return_repeat_summary is True, also return a list containing the max-min range for each try.
        """
        best_code, best_BSAD_range, all_BSAD_ranges = None, self.size(), []
        for i in range(N_tries):
            curr_code = Binary_code(self.length, self.give_N_codewords_random(N))
            bit_sums_across_digits = curr_code.bit_sums_across_digits()
            BSAD_range = max(bit_sums_across_digits) - min(bit_sums_across_digits)
            if BSAD_range < best_BSAD_range:
                best_code = curr_code
                best_BSAD_range = BSAD_range
            all_BSAD_ranges.append(BSAD_range)
        if return_repeat_summary:     return best_code.codewords, all_BSAD_ranges
        else:                         return best_code.codewords

    def add_mirrored_bits(self, bit_position_list):
        """ Return new Binary_code with all codewords extended by mirroring the given bits.

        For example a (0011,1010,1100) code with bit_position_list [0,1] will become (001111,101001,110000).
        The new bits will be added in the order provided; a position may be given multiple times - [0,0,1,1,0] is valid. 
        """
        # MAYBE-TODO add input checking? Make sure that bit_position_list is a sequence/iterator/something, and that new_codeword_01_list[bit_position] isn't an IndexError...
        new_codewords = set()
        new_codeword_length = self.length + len(bit_position_list)
        for codeword in self.codewords:
            new_codeword_01_list = codeword.list()
            for bit_position in bit_position_list:
                new_codeword_01_list.append(not new_codeword_01_list[bit_position])
            new_codewords.add(Binary_codeword(val=new_codeword_01_list))
        return Binary_code(length=new_codeword_length, val=new_codewords, method='list', expected_count=self.size())


    def clonality_count_conflicts(self, N_allowed_changes=(0,0), count_self_conflicts=False,remove_all_zero_codeword=False,
                                  print_conflict_details=False, return_conflict_details=False, quiet=False):
        """ Simple clonality conflict count.  Return a (conflict_count: codeword_set) dictionary.
        Go over all combinations of codewords A,B,C in the code, and whenever the clonality sum A+B is close enough to C  
         according to N_allowed changes (which can be either a single number or a (1_to_0_changes, 0_to_1_changes) tuple)
         add one conflict count to all three codewords. 
        The count_self_conflicts argument specifies whether A+B=A is considered as a problem the same way A+B=C is; 
         if count_self_conflicts is True, you may want to also set remove_all_zero_codeword to True, 
         otherwise it'll conflict with everything (a warning will be printed).
        Codewords with 0 conflicts are guaranteed not to generate problems, but nothing very useful can be said about
         how many and which of the ones with low counts could be added to that set without generating clonality issues.
        The last two arguments allow printing and/or returning of detailed conflict data - a set of (set(A,B),A|B,C,s,n) 
         tuples, where A and B are two codewords in the code, A|B is their clonality result, C is the codeword or set
          of codewords that A|B conflicts with (i.e. is too close to, based on N_allowed_changes), s is 'self' if it was 
          a self-conflict and '' otherwise, and n is N_allowed_changes (since I don't have the real N_changes available).
        """

        # deal with all-zero codeword: remove if requested, print warning if it's False and count_self_conflicts is True
        if remove_all_zero_codeword:    self.remove_extreme_codeword(bit=0)
        elif count_self_conflicts and (Binary_codeword('0'*self.length) in self.codewords) and not quiet:
            print("Warning: you're running a clonality conflict check with count_self_conflicts turned on, and your code "
                  +"contains the all-zero codeword - be aware that it's going to generate a clonality conflict with "
                  +"EVERYTHING. Set the remove_all_zero_codeword argument to True if you'd like to prevent that; "
                  +"set the quiet argument to True to silence this message.")
        # MAYBE-TODO add a remove_all_one_codeword option too?  It's frequently bad to have it in there...

        # set up conflict-count dictionary, with a 0 for each codeword
        codeword_to_conflict_count = dict([(codeword,0) for codeword in self.codewords])
        if return_conflict_details:     all_conflict_details = set()

        ### Special case just for 0 allowed_changes, because the normal way is SLOW
        # MAYBE-TODO it would be better code if the special case wasn't here... But it is faster than the general case.
        # MAYBE-TODO add option to force using the general case even for 0 changes, to make totally sure it's the same?
        if N_allowed_changes in [0, (0,0)]:
            for A,B in itertools.combinations(self.codewords,2):
                clonality_result = A|B
                conflict_details = None
                conflict_info = (frozenset([A,B]),clonality_result,frozenset([clonality_result]))
                # if clonality_result is A or B, what to do depends on if count_self_conflicts is True
                if clonality_result in [A,B]:
                    if count_self_conflicts:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
                        # MAYBE-TODO use short string representations for A,B,clonality_result?  Or only when printing?...
                        conflict_details = conflict_info+('self',N_allowed_changes)
                # if clonality_result is an existing codeword (that's not A or B), add a conflict count
                elif clonality_result in self.codewords: 
                    for codeword in [A,B,clonality_result]:
                        codeword_to_conflict_count[codeword] += 1
                    conflict_details = conflict_info+('',N_allowed_changes)
                if conflict_details and print_conflict_details:     print conflict_details
                if conflict_details and return_conflict_details:    all_conflict_details.add(conflict_details)

        ### Standard case, for when the allowed change count arguments aren't 0: 
        #   pre-calculate a pool of all illegal clonality results based on all the codewords and check against that.
        # MAYBE-TODO why is this still so much slower than the special case version, even with 0,0 arguments?  
        #   Is it the set operations? Ir is it that much slower at all, really?... Do I care enough to fix it?
        else:
            expanded_conflict_values = expand_by_all_mutations_dict(self.codewords, N_allowed_changes)
            for A,B in itertools.combinations(self.codewords,2):
                clonality_result = A|B
                conflict_details = None
                if clonality_result in expanded_conflict_values: 
                    # if the list of base codewords for the conflict contains codewords other than A and B, 
                    #   register a conflict for A, B and all the base codewords (doing a set union 
                    #    ensures that even if A/B were among the base codewords, they only get one conflict)
                    if len(expanded_conflict_values[clonality_result]-set([A,B]))>0:
                        for codeword in expanded_conflict_values[clonality_result].union(set([A,B])):
                            codeword_to_conflict_count[codeword] += 1
                        conflict_set = frozenset(expanded_conflict_values[clonality_result] - set([A,B]))
                        conflict_details = (frozenset([A,B]),clonality_result,conflict_set,'',N_allowed_changes)
                    # if the list of base codewords for the conflict was only A/B, 
                    #   only register a conflict is count_self_conflicts is True
                    elif count_self_conflicts:
                        for codeword in A,B:
                            codeword_to_conflict_count[codeword] += 1
                        conflict_set = frozenset(expanded_conflict_values[clonality_result] & set([A,B]))
                        conflict_details = (frozenset([A,B]),clonality_result,conflict_set,'self',N_allowed_changes)
                if conflict_details and print_conflict_details:     print conflict_details
                if conflict_details and return_conflict_details:    all_conflict_details.add(conflict_details)

        ### generate the final conflict_count:codeword_set dictionary from the codeword:conflict_count one
        conflict_count_to_codeword_set = invert_dict_tolists(codeword_to_conflict_count)
        if return_conflict_details:     return conflict_count_to_codeword_set, all_conflict_details
        else:                           return conflict_count_to_codeword_set

    def clonality_conflict_check(self, N_allowed_changes=(0,0), count_self_conflicts=False, remove_all_zero_codeword=False,
                                  print_conflict_details=False, quiet=False):
        """ Return False if the code contains no clonality conflicts based on arguments, True otherwise.
        Passes all its arguments to clonality_count_conflicts - see that function's docstring for details."""
        # MAYBE-TODO could rewrite this to be much faster by writing it separately and having it just go on until
        #   the first conflict and then return True, instead of using clonality_count_conflicts and thus needing to
        #   go through all the combinations, but I'm not sure it's worth it.
        conflict_count_to_codeword_set = self.clonality_count_conflicts(N_allowed_changes, count_self_conflicts, 
                            remove_all_zero_codeword, print_conflict_details, return_conflict_details=False, quiet=quiet)
        if conflict_count_to_codeword_set.keys() in ([0],[]):   return False
        else:                                                   return True

    # MAYBE-TODO write a function to find the highest N_allowed_changes that results in no clonality conflict?  Maybe just for single-digit N_allowed_changes - otherwise we run into the issue of whether (2,0) or (0,3) is higher, although there could be more options to deal with that.

    def clonality_obvious_no_conflict_subset(self, N_allowed_changes=(0,0), count_self_conflicts=False, 
                                         remove_all_zero_codeword=False, print_conflict_details=False, quiet=False):
        """ Really naive solution of the clonality problem: return a set of codewords with 0 conflicts, 
         as given by clonality_count_conflicts with the same arguments (see that function's docstring for more)."""
        conflict_count_to_codeword_set = self.clonality_count_conflicts(N_allowed_changes, count_self_conflicts, 
                            remove_all_zero_codeword, print_conflict_details, return_conflict_details=False, quiet=quiet)
        try:                return conflict_count_to_codeword_set[0]
        except KeyError:    return set()

    def clonality_grow_no_conflict_subset(self, N_allowed_changes=(0,0), starting_subset=None, more_random=False, 
           N_repeats=1, return_repeat_summary=False, 
           count_self_conflicts=False, remove_all_zero_codeword=False, quiet=False):
        """ Imperfect iterative partially-random codeword addition solution to the clonality problem, close to Goodman2009.

        First use self.clonality_count_conflicts to get the full conflict graph, using the N_allowed_changes, 
         count_self_conflicts and remove_all_zero_codeword argument values provided.
        Then repeat the following N_repeats times and return the best result:
         Starting from starting_subset (empty by default), try adding one codeword at a time to the subset: 
          only add the codeword if the result is conflict-free, and keep going until all codewords have been tried.
          If more_random is True, the codewords will be tried in random order; otherwise they will be sorted by the total
           number of conflicts they participate in, and only random within that.
        Return the resulting no-conflict codeword set; if N_repeats>1 and return_repeat_summary is True, 
         also return a list of the lengths of all the N_repeats codeword sets found.
        """
        ### First get all the detailed conflict-count info
        conflict_count_to_codeword_set, conflict_detail_tuples = self.clonality_count_conflicts(N_allowed_changes, 
                            count_self_conflicts, remove_all_zero_codeword, return_conflict_details=True, quiet=quiet)
        # conflict_count_to_codeword_set is a conflict_count:codeword_set dictionary.
        # conflict_detail_tuples is a set of (set(A,B), A|B, conflicting_codeword_set, if_self, N_allowed_changes) tuples

        # conflict_triples is a set of (A,B,C) frozensets - except they may just be (A,B) if self-conflicts are counted...
        #  I could use tuples instead of frozensets here, but then we'd have the issue of (A,B,C) and (A,C,B) being 
        #    different and counted twice, which seems pointless.
        conflict_triples = set()
        for (AB_set,_,C_set,_,_) in conflict_detail_tuples:
            A,B = tuple(AB_set)
            for C in C_set:     conflict_triples.add(frozenset([A,B,C]))

        # make a codeword:set_of_conflict_pairs dictionary
        codeword_to_conflict_pairs = defaultdict(lambda: set())
        for triple in conflict_triples:
            # not that this works whether the triple is in fact a triple, duple (self-conflicts), or whatever
            for curr_codeword in triple:
                codeword_to_conflict_pairs[curr_codeword].add(triple - frozenset([curr_codeword]))

        # three ways of counting conflicts - ARE they actually the same? Probably not exactly, 
        #  especially when we include self-conflicts and possibilities like a (A,B,C) conflict and an (A,C,B) one... 
        #  MAYBE-TODO add some asserts using approximate relations between those numbers?
        #total_conflicts_1 = sum([ccount*len(cwset) for (ccount,cwset) in conflict_count_to_codeword_set.iteritems()])/3
        #total_conflicts_2 = len(conflict_triples) 
        #total_conflicts_3 = sum([len(pairs) for pairs in codeword_to_conflict_pairs.values()])

        def count_conflicts_codeword_and_set(codeword_conflict_pairs, codeword_set):
            return len([pair for pair in codeword_conflict_pairs if pair.issubset(codeword_set)])

        ### set the starting subset (empty unless provided); make sure it's conflict-free and part of the current code
        starting_subset = set() if starting_subset is None else starting_subset
        # if a code instead of a set of codewords was given as starting_subset, just take its codewords, that's fine
        if isinstance(starting_subset, Binary_code):    starting_subset = starting_subset.codewords
        if remove_all_zero_codeword:
            if '0'*self.length in starting_subset:   starting_subset.remove('0'*self.length)
        for codeword in starting_subset:
            if count_conflicts_codeword_and_set(codeword_to_conflict_pairs[codeword], starting_subset) > 0:
               raise BinaryCodeError("starting_subset provided to clonality_grow_no_conflict_subset is not conflict-free!")
        if not all([codeword in self.codewords for codeword in starting_subset]):
            raise BinaryCodeError("starting_subset provided to clonality_grow_no_conflict_subset is not part of the code!")

        ### repeat randomly generating an order and making a subset multiple times, return the best (and a trial summary)
        best_subset = set()
        all_subset_lengths = []
        multiple_subsets, multiple_codeword_addition_orders, prev_codewords_addition_order = False, False, []
        for i in range(N_repeats):
            ### what order all the codewords should be attempted-added in: (ignore codewords already in starting_subset)
            # default: based on conflict-count per codeword (lowest first), random within that
            if not more_random:
                codewords_to_add = []
                for ccount, cwset in sorted(conflict_count_to_codeword_set.iteritems()):
                    cwlist = [c for c in cwset if c not in starting_subset]
                    random.shuffle(cwlist)
                    codewords_to_add += cwlist
            # if more_random: completely random without regard to conflict-count
            else:
                codewords_to_add = [c for c in self.codewords if c not in starting_subset]
                random.shuffle(codewords_to_add)
            assert len(codewords_to_add) + len(starting_subset) == len(self.codewords)

            # go over the codewords_to_add list: if the current codeword has no conflicts with current_subset, add it, 
            #  otherwise skip and go on to the next one.
            # using set() here so that current_subset is a real copy of starting_subset, not two labels for one object!
            current_subset = set(starting_subset)   
            for codeword in codewords_to_add:
                assert codeword not in current_subset, "Error: shouldn't be adding an already present codeword!"
                if count_conflicts_codeword_and_set(codeword_to_conflict_pairs[codeword], current_subset) == 0:
                    current_subset.add(codeword)

            # check if the codeword addition order and resulting subset is always the same or not
            if not multiple_codeword_addition_orders and prev_codewords_addition_order!=[]:
                if prev_codewords_addition_order != codewords_to_add:
                    multiple_codeword_addition_orders = True
            prev_codewords_addition_order = codewords_to_add
            if not multiple_subsets and best_subset!=set() and current_subset != best_subset:
                multiple_subsets = True

            all_subset_lengths.append(len(current_subset))
            if len(current_subset) > len(best_subset):
                best_subset = current_subset

        # check that the results are sane; print warnings if there's something suspicions
        assert len(best_subset) == max(all_subset_lengths)
        if not quiet and N_repeats>1 and (len(self.codewords)-len(starting_subset) > 1):
            if not multiple_codeword_addition_orders:
                print("WARNING: Only 1 random order of %s elements in %s repeats - RANDOMNESS PROBABLY FAILING!"
                      %(len(codewords_to_add), N_repeats))
        if not quiet and N_repeats>1 and not multiple_subsets:
            print("Warning: The same subset always found in %s repeats - something may be wrong!"%N_repeats)

        if return_repeat_summary:   return best_subset, all_subset_lengths
        else:                       return best_subset

    # MAYBE-TODO implement some other options for reducing the set to one without clonality issues?
    #  * Any other sensible algorithms for doing this?  See Clonality section of ../notes_combinatorial_pooling_theory.txt - I had some new ideas...


class Testing__Binary_codeword(unittest.TestCase):
    """ Testing Binary_codeword functionality. """

    def test__creation_value_types(self):
        """ Test Binary_codeword instance initiation with different value types (string, int, list, Binary_codeword) """
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword(7).string(), '111')
        self.assertEqual(Binary_codeword('0b 111\n').string(), '111')
        self.assertEqual(Binary_codeword([True,True,True]).string(), '111')
        self.assertEqual(Binary_codeword(Binary_codeword('111')).string(), '111')

    def test__initiation_length_check_and_padding(self):
        """ Test if Binary_codeword initiation deals with codeword length correctly (length check and padding). """
        #   '111' has length 3, should work (no assert, just making sure it doesn't raise an exception)
        Binary_codeword('111',3,check_length=True)
        #   '111' doesn't have length 5, should raise an error
        self.assertRaises(BinaryCodeError, Binary_codeword, '111', 5, check_length=True)
        #   padding a length-3 string to the same length shouldn't change anything
        self.assertEqual(Binary_codeword('111',3), Binary_codeword('111'))
        #   padding a length-3 string to a higher length should work
        self.assertEqual(Binary_codeword(7,4).string(), '0111')
        self.assertEqual(Binary_codeword('111',5).string(), '00111')
        #   padding a length-3 string to length 2 should raise an error
        self.assertRaises(BinaryCodeError, Binary_codeword, '111', 2)

    def test__equality_and_comparison_operators(self):
        self.assertTrue(Binary_codeword('111') == Binary_codeword('111'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('101'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('1110'))
        self.assertTrue(Binary_codeword('111') != Binary_codeword('0111'))
        self.assertTrue(Binary_codeword('101') < Binary_codeword('111'))

    def test__bitwise_operations(self):
        self.assertEqual(~Binary_codeword('000'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') | Binary_codeword('011'), Binary_codeword('111'))
        self.assertEqual(Binary_codeword('110') & Binary_codeword('011'), Binary_codeword('010'))
        self.assertEqual(Binary_codeword('110') ^ Binary_codeword('011'), Binary_codeword('101'))
        # comparing bitstrings of different lengths should fail
        self.assertRaises(ValueError, Hamming_distance, Binary_codeword('000'), Binary_codeword('0000'))

    def test__length_calculation(self):
        self.assertEqual(len(Binary_codeword('111')), 3)
        self.assertEqual(len(Binary_codeword('00000')), 5)

    def test__weight_calculation(self):
        self.assertEqual(Binary_codeword('111').weight(), 3)
        self.assertEqual(Binary_codeword('001').weight(), 1)
        self.assertEqual(Binary_codeword('00000').weight(), 0)

    def test__list_string_representations(self):
        self.assertEqual(Binary_codeword('111').string(), '111')
        self.assertEqual(Binary_codeword('111').list(), [1,1,1])


class Testing__other_functions(unittest.TestCase):
    """ Testing functions that aren't part of either of the main classes."""
    
    def test__Hamming_distance_calculation(self):
        self.assertEqual(Hamming_distance(Binary_codeword('000'),Binary_codeword('000')), 0)
        self.assertEqual(Hamming_distance(Binary_codeword('111'),Binary_codeword('000')), 3)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('000')), 2)
        self.assertEqual(Hamming_distance(Binary_codeword('101'),Binary_codeword('010')), 3)

    def test__bit_change_count(self):
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('000')) == (0,0)
        assert bit_change_count(Binary_codeword('111'),Binary_codeword('000')) == (0,3)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('111')) == (3,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('000')) == (0,2)
        assert bit_change_count(Binary_codeword('000'),Binary_codeword('101')) == (2,0)
        assert bit_change_count(Binary_codeword('101'),Binary_codeword('010')) == (1,2)
        assert bit_change_count(Binary_codeword('010'),Binary_codeword('101')) == (2,1)

    def test__change_all_position_combinations(self):
        assert _change_all_position_combinations([1,1], 0, None) == set([(1,1)])
        assert _change_all_position_combinations([1,1], 1, None) == set([(1,1),(0,1),(1,0)])
        assert _change_all_position_combinations([1,1], 2, None) == set([(1,1),(0,1),(1,0),(0,0)])
        assert _change_all_position_combinations([1,1], 10, None) == set([(1,1),(0,1),(1,0),(0,0)])
        assert _change_all_position_combinations([1,1], 1, [0]) == set([(1,1),(0,1)])
        assert _change_all_position_combinations([1,1], 1, [1]) == set([(1,1),(1,0)])
        assert _change_all_position_combinations([1,1], 1, [0,1]) == set([(1,1),(0,1),(1,0)])
        assert _change_all_position_combinations([0,0], 0, None) == set([(0,0)])
        assert _change_all_position_combinations([0,0], 1, None) == set([(0,0),(0,1),(1,0)])

    def test__expand_by_all_mutations(self):
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        all_test_codewords = set([b11,b10,b01,b00])
        single_value_results = dict()
        ### setup of values
        # N_permitted changes as a single number
        single_value_results[(b11,0)] = set([b11])
        single_value_results[(b11,1)] = set([b11,b10,b01])
        single_value_results[(b11,2)] = set([b11,b10,b01,b00])
        single_value_results[(b00,0)] = set([b00])
        single_value_results[(b00,1)] = set([b00,b10,b01])
        single_value_results[(b00,2)] = set([b00,b10,b01,b11])
        single_value_results[(b01,0)] = set([b01])
        single_value_results[(b01,1)] = set([b01,b00,b11])
        single_value_results[(b01,2)] = set([b01,b00,b11,b10])
        single_value_results[(b10,0)] = set([b10])
        single_value_results[(b10,1)] = set([b10,b00,b11])
        single_value_results[(b10,2)] = set([b10,b00,b11,b01])
        # N_permitted changes as a (N_permitted_1_to_0,N_permitted_0_to_1) tuple
        for N_changes in [(0,0),(0,1),(0,2)]:   single_value_results[(b11,N_changes)] = set([b11])
        for N_changes in [(1,0),(1,1),(1,2)]:   single_value_results[(b11,N_changes)] = set([b11,b01,b10])
        for N_changes in [(2,0),(2,1),(2,2)]:   single_value_results[(b11,N_changes)] = set([b11,b01,b10,b00])
        for N_changes in [(0,0),(1,0),(2,0)]:       single_value_results[(b00,N_changes)] = set([b00])
        for N_changes in [(0,1),(1,1),(2,1)]:       single_value_results[(b00,N_changes)] = set([b00,b01,b10])
        for N_changes in [(0,2),(1,2),(2,2)]:       single_value_results[(b00,N_changes)] = set([b00,b01,b10,b11])
        for N_changes in [(0,0)]:                       single_value_results[(b01,N_changes)] = set([b01])
        for N_changes in [(0,1),(0,2)]:                 single_value_results[(b01,N_changes)] = set([b01,b11])
        for N_changes in [(1,0),(2,0)]:                 single_value_results[(b01,N_changes)] = set([b01,b00])
        for N_changes in [(1,1),(1,2),(2,1),(2,2)]:     single_value_results[(b01,N_changes)] = set([b01,b00,b11,b10])
        for N_changes in [(0,0)]:                           single_value_results[(b10,N_changes)] = set([b10])
        for N_changes in [(0,1),(0,2)]:                     single_value_results[(b10,N_changes)] = set([b10,b11])
        for N_changes in [(1,0),(2,0)]:                     single_value_results[(b10,N_changes)] = set([b10,b00])
        for N_changes in [(1,1),(1,2),(2,1),(2,2)]:         single_value_results[(b10,N_changes)] = set([b10,b00,b11,b01])
        # now do the actual testing based on the defined values
        N_changes_values = [0,1,2]
        for N_changes in N_changes_values+list(itertools.product(N_changes_values,repeat=2)):
            for test_set_size in range(5):
                for test_codewords in itertools.combinations(all_test_codewords,test_set_size):
                    # test expand_by_all_mutations - simply a union of the result sets
                    all_result_sets = [single_value_results[(c,N_changes)] for c in test_codewords]
                    assert expand_by_all_mutations(test_codewords, N_changes) == set().union(*all_result_sets)
                    # text expand_by_all_mutations_dict - an inverted dictionary of the above, pretty much
                    all_result_dict = dict([(c,single_value_results[(c,N_changes)]) for c in test_codewords])
                    result_to_base_dict = invert_listdict_tolists(all_result_dict)
                    assert expand_by_all_mutations_dict(test_codewords, N_changes) == result_to_base_dict


class Testing__Binary_code__most_functions(unittest.TestCase):
    """ Testing Binary_code functionality, except for clonality-conflict functions, which have their own test suite."""

    # MAYBE-TODO convert all the assert statements to self.assertEqual or self.assertTrue or such? 
    #   That's the way unittest functions should be written, but the current version works too...
    #   I could probably just use nosetest if I wanted - that catches normal assert statements too. 

    def test__creation_from_list_and_properties(self):
        for l in [0,1,5,100]:
            B = Binary_code(l,[])
            assert B.length == l
            assert B.size() == 0
            assert B.find_Hamming_distance_range() == (None,None)
            assert B.find_bit_sum_counts() == []
            assert B.total_bit_sum() == 0
        B = Binary_code(3,['110','101','011','000'])
        assert B.length == 3
        assert B.size() == 4
        assert B.find_Hamming_distance_range() == (2,2)
        assert B.find_bit_sum_counts() == [(0,1), (2,3)]
        assert B.total_bit_sum() == 6
        # check that creation fails if the length is wrong
        self.assertRaises(BinaryCodeError, Binary_code, 4, ['110','101','011','000'])
        # check that creation fails if the expected count is wrong
        self.assertRaises(BinaryCodeError, Binary_code, 3, ['110','101','011','000'], expected_count=5)
        # check that creation fails with inexistent method keyword
        self.assertRaises(BinaryCodeError, Binary_code, 4, ['110','101','011','000'], method='random')

    def test__creation_from_matrix_generator_file(self):
        # (also implicitly checks generation from a matrix object)
        infile1 = 'error-correcting_codes/19-10-5_generator'
        try:            B19 = Binary_code(19,val=infile1,method='matrixfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run matrix file test."%infile1)
        assert B19.find_bit_sum_counts() == [(0,1), (5,30), (6,64), (7,90), (8,150), (9,180), (10,168), (11,156), (12,104), (13,46), (14,24), (15,10), (16,1)]
        B20 = B19.add_parity_bit()
        assert B20.size() == 2**10
        assert B20.length == 20
        assert B20.find_bit_sum_counts() == [(0, 1), (6, 94), (8, 240), (10, 348), (12, 260), (14, 70), (16, 11)]

    def test__creation_from_code_list_file(self):
        B19 = Binary_code(19,val='error-correcting_codes/19-10-5_generator',method='matrixfile',expected_count=2**10)
        B20 = B19.add_parity_bit()
        infile2 = 'error-correcting_codes/20-10-6_list'
        try:            B20_new = Binary_code(20,val=infile2,method='listfile',expected_count=2**10)
        except IOError: sys.exit("Couldn't find input file %s to run list file test."%infile2)
        assert B20_new == B20

    def test__codeword_add_remove(self):
        B = Binary_code(3,['110','101','011','000'])
        C = Binary_code(3,B.codewords)
        # add an element, verify that the Binary_code properties are correct
        C.add('111')
        assert C.length == B.length
        assert C.size() == B.size() + 1
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == B.find_bit_sum_counts() + [(3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3
        # remove an element, verify that the Binary_code properties are correct
        C.remove('110')
        assert C.length == B.length
        assert C.size() == B.size()
        assert C.find_Hamming_distance_range() == (1,3)
        assert C.find_bit_sum_counts() == [(0,1), (2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum() + 3 - 2
        # remove should fail if the element wasn't there
        self.assertRaises(BinaryCodeError, B.remove, '111')
        # add/remove should fail if the length is wrong
        self.assertRaises(BinaryCodeError, B.add, '1111')
        self.assertRaises(BinaryCodeError, B.remove, '1111')

    def test__remove_extreme_codeword(self):
        B = Binary_code(3,['111','101','011','000'])
        C = Binary_code(3,B.codewords)
        ### removing all-zero codeword
        assert C.remove_extreme_codeword(bit=0) == 1
        assert C.length == B.length
        assert C.size() == B.size() - 1
        assert C.find_Hamming_distance_range() == (1,2)
        assert C.find_bit_sum_counts() == [(2,2), (3,1)]
        assert C.total_bit_sum() == B.total_bit_sum()
        # try removing it again - should have no effect
        assert C.remove_extreme_codeword(bit=0) == 0
        assert C.length == B.length
        assert C.size() == B.size() - 1
        ### removing all-one codeword
        assert C.remove_extreme_codeword(bit=1) == 1
        assert C.length == B.length
        assert C.size() == B.size() - 2
        assert C.find_Hamming_distance_range() == (2,2)
        assert C.find_bit_sum_counts() == [(2,2)]
        assert C.total_bit_sum() == B.total_bit_sum() - 3
        # try removing it again - should have no effect
        assert C.remove_extreme_codeword(bit=1) == 0
        assert C.length == B.length
        assert C.size() == B.size() - 2

    def test__various_bit_sums(self):
        """ Tests find_bit_sum_counts, total_bit_sum, bit_sums_across_digits. """
        A = Binary_code(1,[])
        assert A.find_bit_sum_counts() == []
        assert A.total_bit_sum() == 0
        assert A.bit_sums_across_digits() == [0]
        A3 = Binary_code(3,[])
        assert A3.find_bit_sum_counts() == []
        assert A3.total_bit_sum() == 0
        assert A3.bit_sums_across_digits() == [0,0,0]
        C = Binary_code(1,['1','0'])
        assert C.find_bit_sum_counts() == [(0,1),(1,1)]
        assert C.total_bit_sum() == 1
        assert C.bit_sums_across_digits() == [1]
        D = Binary_code(2,['11','10','01','00'])
        assert D.find_bit_sum_counts() == [(0,1),(1,2),(2,1)]
        assert D.total_bit_sum() == 4
        assert D.bit_sums_across_digits() == [2,2]
        B = Binary_code(3,['110','101','011','000'])
        assert B.find_bit_sum_counts() == [(0,1),(2,3)]
        assert B.total_bit_sum() == 6
        assert B.bit_sums_across_digits() == [2,2,2]

    def test__invert(self):
        B = Binary_code(3,['110','101','011','000'])
        B_ = B.invert()
        assert B_.length == B.length
        assert B_.size() == B_.size()
        assert B_.find_Hamming_distance_range() == B.find_Hamming_distance_range()
        assert B_.find_bit_sum_counts() == sorted([(B.length-w,n) for (w,n) in B.find_bit_sum_counts()])
        assert B_.total_bit_sum() == B.size() * B.length - B.total_bit_sum()

    def test__adding_parity_bit(self):
        D = Binary_code(2,['11','10','01','00'])
        assert D.length == 2
        assert D.size() == 4
        assert D.find_Hamming_distance_range() == (1,2)
        assert D.find_bit_sum_counts() == [(0,1), (1,2), (2,1)]
        assert D.total_bit_sum() == 4
        E = D.add_parity_bit()
        assert E.length == D.length + 1
        assert E.size() == D.size()
        assert E.find_Hamming_distance_range() == (2,2)
        assert E.find_bit_sum_counts() == [(0,1), (2,3)]
        assert E.total_bit_sum() == D.total_bit_sum() + 2

    def test__choose_codewords_by_bit_sum(self):
        D = Binary_code(2,['11','10','01','00'])
        E = Binary_code(2,['11','00'])
        # test normal cases with replace_self False
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(0,-1)]) == set(['00','01','10','11'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(0,2)]) == set(['00','01','10','11'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(0,1)]) == set(['00','01','10'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(0,0)]) == set(['00'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(1,2)]) == set(['01','10','11'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(1,-1)]) == set(['01','10','11'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(1,1)]) == set(['01','10'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(2,2)]) == set(['11'])
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(2,-1)]) == set(['11'])
        assert set([c.string() for c in E.choose_codewords_by_bit_sum(0,2)]) == set(['00','11'])
        assert set([c.string() for c in E.choose_codewords_by_bit_sum(0,-1)]) == set(['00','11'])
        assert set([c.string() for c in E.choose_codewords_by_bit_sum(0,0)]) == set(['00'])
        assert set([c.string() for c in E.choose_codewords_by_bit_sum(2,2)]) == set(['11'])
        assert set([c.string() for c in E.choose_codewords_by_bit_sum(1,1)]) == set()
        # test that empty sets are returned if the low-high range is empty
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(3,-1)]) == set()
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(4,100)]) == set()
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(4,1)]) == set()
        assert set([c.string() for c in D.choose_codewords_by_bit_sum(1,0)]) == set()
        # test the replace_self True version
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(0,2,replace_self=True)
        assert D_ == D
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(0,0,replace_self=True)
        assert D_.size() == 1
        assert D_.find_bit_sum_counts() == [(0,1)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(0,1,replace_self=True)
        assert D_.size() == 3
        assert D_.find_bit_sum_counts() == [(0,1),(1,2)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(1,2,replace_self=True)
        assert D_.size() == 3
        assert D_.find_bit_sum_counts() == [(1,2),(2,1)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(1,1,replace_self=True)
        assert D_.size() == 2
        assert D_.find_bit_sum_counts() == [(1,2)]
        D_ = Binary_code(2,['11','10','01','00'])
        D_.choose_codewords_by_bit_sum(0,-1,replace_self=True)
        assert D_ == D
        E = Binary_code(2,['11','00'])
        E.choose_codewords_by_bit_sum(1,1,replace_self=True)
        assert E.size() == 0

    def test__give_N_codewords_random(self):
        D = Binary_code(2,['11','10','01','00'])
        # do the test several times, since randomness is involved
        for i in range(10):
            # the actual codewords chosen are random, so just check that there's the right number of them,
            #  that they're all unique, and that they all come from the original set
            for N in range(5):
                codewords = D.give_N_codewords_random(N)
                assert len(codewords) == N
                assert len(set(codewords)) == len(codewords)
                assert all([c in D.codewords for c in codewords])
            # check that it's impossible to get 5 elements from a 4-element binary code
            self.assertRaises(BinaryCodeError, D.give_N_codewords_random, 5)

    def test__give_N_codewords_by_bit_sum(self):
        D = Binary_code(2,['11','10','01','00'])
        # do the test several times, since randomness is involved
        for i in range(10):
            # the actual codewords chosen are random, so just check that there's the right number of them,
            #  that they're all unique, and that they all come from the original set
            for take_high in True, False:
                for N in range(5):
                    codewords = D.give_N_codewords_by_bit_sum(N, take_high=take_high)
                    assert len(codewords) == N
                    assert len(set(codewords)) == len(codewords)
                    assert all([c in D.codewords for c in codewords])
            # check that the bit-sums of the codewords are low or high as expected: first for take_high False, then True
            #  (note that when picking 2 codewords, there's some variability in the bit-sum range!)
            for N in 1,3,4:
                assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(N, False)]) == 0
            assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(2, False)]) in (0,1)
            assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(1, False)]) == 0
            assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(2, False)]) == 1
            assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(3, False)]) == 1
            assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(4, False)]) == 2
            for N in 1,3,4:
                assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(N, True)]) == 2
            assert max([c.weight() for c in D.give_N_codewords_by_bit_sum(2, True)]) in (1,2)
            assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(1, True)]) == 2
            assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(2, True)]) == 1
            assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(3, True)]) == 1
            assert min([c.weight() for c in D.give_N_codewords_by_bit_sum(4, True)]) == 0
            # check that it's impossible to get 5 elements from a 4-element binary code
            for take_high in True, False:
                self.assertRaises(BinaryCodeError, D.give_N_codewords_by_bit_sum, 5, take_high=take_high)

    def test__give_N_codewords_even_distribution(self):
        D = Binary_code(2,['11','10','01','00'])
        # just testing that the basic properties are the same as of give_N_codewords_random, 
        #  and that the second return value (added due to return_repeat_summary=True) is sane.
        for N in range(5):
            codewords, all_ranges = D.give_N_codewords_even_distribution(N, N_tries=50, return_repeat_summary=True)
            assert len(codewords) == N
            assert len(set(codewords)) == len(codewords)
            assert all([c in D.codewords for c in codewords])
            assert len(all_ranges)==50
            assert min(all_ranges)>=0 and max(all_ranges)<=2
            result_pool_counts = Binary_code(2,codewords).bit_sums_across_digits()
            assert max(result_pool_counts)-min(result_pool_counts) == min(all_ranges)
            # check that it's impossible to get 5 elements from a 4-element binary code
            for N_tries in [1,5,30]:
                self.assertRaises(BinaryCodeError, D.give_N_codewords_even_distribution, 5, N_tries)
        ### a few specific tests
        D = Binary_code(2,['11','10','01','00'])
        B = Binary_code(3,['110','101','011','000'])
        # when choosing 0 or 4 codewords out of 4, the outcome will always be the same
        for code in D,B:
            codewords, all_ranges = code.give_N_codewords_even_distribution(0, N_tries=50, return_repeat_summary=True)
            assert codewords == set() and set(all_ranges) == set([0])
            codewords, all_ranges = code.give_N_codewords_even_distribution(4, N_tries=50, return_repeat_summary=True)
            assert codewords == code.codewords and set(all_ranges) == set([0])
        ### when choosing 2 or 3 codewords out of 4, the outcome should vary:
        # for code D, it's possible to pick 2 or 3 perfectly distributed codewords (when 3, the set will include 00)
        for N in [2,3]:
            codewords, all_ranges = D.give_N_codewords_even_distribution(N, N_tries=50, return_repeat_summary=True)
            result_pool_counts = Binary_code(2,codewords).bit_sums_across_digits()
            assert set(all_ranges) == set([0,1])
            assert max(result_pool_counts)-min(result_pool_counts) == 0
        # for code B, it's possible to pick 3 but not 2 perfectly distributed codewords
        codewords, all_ranges = B.give_N_codewords_even_distribution(2, N_tries=50, return_repeat_summary=True)
        result_pool_counts = Binary_code(3,codewords).bit_sums_across_digits()
        assert set(all_ranges) == set([1])
        assert max(result_pool_counts)-min(result_pool_counts) == 1
        codewords, all_ranges = B.give_N_codewords_even_distribution(3, N_tries=50, return_repeat_summary=True)
        result_pool_counts = Binary_code(3,codewords).bit_sums_across_digits()
        assert set(all_ranges) == set([0,1])
        assert max(result_pool_counts)-min(result_pool_counts) == 0

    def test__add_mirrored_bits(self):
        D = Binary_code(2,['11','10','01','00'])
        # mirroring no bits should yield identical code
        assert D.add_mirrored_bits([]) == D
        # check a few cases of mirroring some bits
        assert D.add_mirrored_bits([0]) == Binary_code(3,['110','100','011','001'])
        assert D.add_mirrored_bits([1]) == Binary_code(3,['110','101','010','001'])
        assert D.add_mirrored_bits([0,1]) == Binary_code(4,['1100','1001','0110','0011'])
        assert D.add_mirrored_bits([0,1,0,1]) == Binary_code(6,['110000','100101','011010','001111'])
        assert D.add_mirrored_bits([0,0,0,1,1,1]) == Binary_code(8,['11000000','10000111','01111000','00111111'])
        # anything that works as an index is legal, so negative bit positions should just take from the end
        assert D.add_mirrored_bits([-1]) == D.add_mirrored_bits([1])
        assert D.add_mirrored_bits([-2]) == D.add_mirrored_bits([0])
        ### same thing, different code
        B = Binary_code(3,['110','101','011','000'])
        assert B.add_mirrored_bits([]) == B
        assert B.add_mirrored_bits([0]) == Binary_code(4,['1100','1010','0111','0001'])
        assert B.add_mirrored_bits([2,2,2]) == Binary_code(6,['110111','101000','011000','000111'])
        assert B.add_mirrored_bits([-2]) == B.add_mirrored_bits([1])
        ### errors:
        # the argument must be a sequence/iterator/something of ints
        self.assertRaises(TypeError, D.add_mirrored_bits, 5)
        self.assertRaises(TypeError, D.add_mirrored_bits, 'abc')
        self.assertRaises(TypeError, D.add_mirrored_bits, ['a'])
        # the index must be within the range of the length of the code's codewords
        self.assertRaises(IndexError, D.add_mirrored_bits, [100])
        self.assertRaises(IndexError, D.add_mirrored_bits, [2])
        self.assertRaises(IndexError, D.add_mirrored_bits, [-3])
        self.assertRaises(IndexError, B.add_mirrored_bits, [3])
        self.assertRaises(IndexError, B.add_mirrored_bits, [-4])


class Testing__Binary_code__clonality_conflict_functions(unittest.TestCase):
    """ Tests clonality_count_conflicts, clonality_conflict_check, clonality_obvious_no_conflict_subset, 
    clonality_grow_no_conflict_subset (all those functions are related in more or less trivial ways).
    Each test function here tests all the functions listed in parallel - easiest that way."""

    def test__empty_code_gives_no_conflicts(self):
        """ Empty codes should return no conflicts with all option combinations."""
        A = Binary_code(3,[])
        for SC,RZ in [(True,True,),(True,False),(False,True),(False,False)]:
            for N_changes in [0,2,(0,0),(1,0),(0,1),(3,3),(1,10)]:
                assert A.clonality_count_conflicts(N_changes,SC,RZ,return_conflict_details=True,quiet=True) == ({},set())
                assert A.clonality_conflict_check(N_changes,SC,RZ,quiet=True) == False 
                assert A.clonality_obvious_no_conflict_subset(N_changes,SC,RZ,quiet=True) == set()
                for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                    assert A.clonality_grow_no_conflict_subset(N_changes, more_random=more_random,
                           N_repeats=N_repeats, count_self_conflicts=SC,remove_all_zero_codeword=RZ,quiet=True) == set()
                    subset,sizes = A.clonality_grow_no_conflict_subset(N_changes, more_random=more_random, 
                                                       N_repeats=N_repeats, return_repeat_summary=True, 
                                                       count_self_conflicts=SC,remove_all_zero_codeword=RZ,quiet=True)
                    assert subset == set() and set(sizes) == set([0])

    def test__count_self_conflicts__and__remove_all_zero_codeword(self):
        """ testing that count_self_conflicts and remove_all_zero_codeword work, with N_allowed_changes 0."""
        # defining the binary codewords so I can use them to check that the results are right
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        # with no self-conflict counted, whether or not all-zero codeword is removed - no conflicts
        B = Binary_code(3,[b110,b101,b011,b000])
        full_set = frozenset([b110,b101,b011,b000])
        set_no_zero = frozenset([b110,b101,b011])
        for RZ,out_set in [(False,full_set), (True,set_no_zero)]:
            assert B.clonality_count_conflicts(0,False,RZ,return_conflict_details=True,quiet=True) == ({0:out_set}, set())
            assert B.clonality_conflict_check(0,False,RZ,quiet=True) == False 
            assert B.clonality_obvious_no_conflict_subset(0,False,RZ,quiet=True) == out_set
            for starting_set in [None, out_set, set([b110]), set([b101])]:
                for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                    assert B.clonality_grow_no_conflict_subset(0,starting_subset=starting_set, more_random=more_random, 
                                                   N_repeats=N_repeats, remove_all_zero_codeword=RZ,quiet=True) == out_set
                    subset,sizes = B.clonality_grow_no_conflict_subset(0, starting_subset=None, more_random=more_random, 
                                   N_repeats=N_repeats, return_repeat_summary=True, remove_all_zero_codeword=RZ,quiet=True)
                    assert subset == out_set and set(sizes) == set([len(out_set)])
        # with self-conflict counted, and with all-zero keyword - clonality conflicts
        B = Binary_code(3,[b110,b101,b011,b000])    # need to re-make it because the all-zero codeword was removed above
        conflicts = set([(frozenset([b000,x]),x,frozenset([x]),'self',0) for x in set_no_zero])
        result = {1:set_no_zero, 3:set([b000])}
        assert B.clonality_count_conflicts(0,True,False,return_conflict_details=True,quiet=True) == (result,conflicts)
        assert B.clonality_conflict_check(0,True,False,quiet=True) == True 
        assert B.clonality_obvious_no_conflict_subset(0,True,False,quiet=True) == set()
        for N_repeats in [1,5,30]:
            assert B.clonality_grow_no_conflict_subset(0,count_self_conflicts=True, more_random=False, 
                                       N_repeats=N_repeats, remove_all_zero_codeword=False, quiet=True) == set_no_zero
            result_with_more_random = B.clonality_grow_no_conflict_subset(0,count_self_conflicts=True, more_random=False, 
                                       N_repeats=N_repeats, remove_all_zero_codeword=False, quiet=True)
            assert result_with_more_random==set_no_zero or len(result_with_more_random)==1
        subset,sizes = B.clonality_grow_no_conflict_subset(0, more_random=False, N_repeats=50, 
                           return_repeat_summary=True, count_self_conflicts=True,remove_all_zero_codeword=False,quiet=True)
        assert subset == set_no_zero and set(sizes) == set([len(set_no_zero)])
        subset,sizes = B.clonality_grow_no_conflict_subset(0, more_random=True, N_repeats=50, 
                           return_repeat_summary=True, count_self_conflicts=True,remove_all_zero_codeword=False,quiet=True)
        assert subset == set_no_zero and set(sizes) == set([1,len(set_no_zero)])

        # with self-conflict counted, but removing all-zero keyword - no clonality conflicts
        assert B.clonality_count_conflicts(0,True,True,return_conflict_details=True,quiet=True) == ({0:set_no_zero}, set())
        assert B.clonality_conflict_check(0,True,True,quiet=True) == False 
        assert B.clonality_obvious_no_conflict_subset(0,True,True,quiet=True) == set_no_zero
        for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
            assert B.clonality_grow_no_conflict_subset(0,count_self_conflicts=True, more_random=more_random,
                                           N_repeats=N_repeats, remove_all_zero_codeword=True,quiet=True) == set_no_zero
        subset,sizes = B.clonality_grow_no_conflict_subset(0, more_random=True, N_repeats=50, 
                           return_repeat_summary=True, count_self_conflicts=True,remove_all_zero_codeword=False,quiet=True)
        assert subset == set_no_zero and set(sizes) == set([len(set_no_zero)])

    def test__count_self_conflicts__allowed_changes_zero(self):
        """ testing count_self_conflicts with N_allowed_changes 0 but without involving the all-zero codeword."""
        # defining the binary codewords so I can use them to check that the results are right
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        C = Binary_code(2,[b01,b11])
        full_set = frozenset([b01,b11])
        # with no self-conflict - no conflicts
        assert C.clonality_count_conflicts(0,False,return_conflict_details=True,quiet=True) == ({0:full_set}, set())
        assert C.clonality_conflict_check(0,False,quiet=True) == False 
        assert C.clonality_obvious_no_conflict_subset(0,False,quiet=True) == full_set
        for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
            assert C.clonality_grow_no_conflict_subset(0,more_random=more_random,N_repeats=N_repeats,
                                                       count_self_conflicts=False,quiet=True) == full_set
        subset,sizes = C.clonality_grow_no_conflict_subset(0, more_random=True, N_repeats=50, 
                                                       return_repeat_summary=True, count_self_conflicts=False,quiet=True)
        assert subset == full_set and set(sizes) == set([len(full_set)])
        # with self-conflict on - conflicts!
        conflicts = set([(full_set,b11,frozenset([b11]),'self',0)])
        assert C.clonality_count_conflicts(0,True,return_conflict_details=True,quiet=True) == ({1:full_set},conflicts)
        assert C.clonality_conflict_check(0,True,quiet=True) == True 
        assert C.clonality_obvious_no_conflict_subset(0,True,quiet=True) == set()
        for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
            assert len(C.clonality_grow_no_conflict_subset(0,count_self_conflicts=True,
                                                           more_random=more_random,N_repeats=N_repeats,quiet=True)) == 1
        subset,sizes = C.clonality_grow_no_conflict_subset(0, more_random=True, N_repeats=50, 
                                                       return_repeat_summary=True, count_self_conflicts=True,quiet=True)
        assert set(sizes) == set([1])

    def test__count_self_conflicts__allowed_changes_nonzero(self):
        """ testing count_self_conflicts with N_allowed_changes other than 0; remove_all_zero_codeword not involved."""
        # defining the binary codewords so I can use them to check that the results are right
        [b11,b10,b01,b00] = [Binary_codeword(x) for x in ['11','10','01','00']]
        C = Binary_code(2,[b01,b10])
        full_set = frozenset([b01,b10])
        # with 0 allowed 0-to-1 changes - no conflicts even with self-conflict on
        for NC in [0,(0,0),(1,0),(2,0),(9,0)]:
            assert C.clonality_count_conflicts(NC,True,return_conflict_details=True,quiet=True) == ({0:full_set}, set())
            assert C.clonality_conflict_check(NC,True,quiet=True) == False 
            assert C.clonality_obvious_no_conflict_subset(NC,True,quiet=True) == full_set
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert C.clonality_grow_no_conflict_subset(NC,count_self_conflicts=True,
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True) == full_set
        # with at least one allowed 0-to-1 change but no self-conflict - no conflicts
        for NC in [1,2,10,(0,1),(0,2),(0,9),(1,1),(2,2),(9,9)]:
            assert C.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) == ({0:full_set}, set())
            assert C.clonality_conflict_check(NC,False,quiet=True) == False 
            assert C.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == full_set
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert C.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True) == full_set
        # with at least one allowed 0-to-1 change and self-conflict on - conflicts!
        for NC in [1,2,10,(0,1),(0,2),(0,9),(1,1),(2,2),(9,9)]:
            conflicts = set([(full_set,b11,full_set,'self',NC)])
            assert C.clonality_count_conflicts(NC,True,return_conflict_details=True,quiet=True) == ({1:full_set},conflicts)
            assert C.clonality_conflict_check(NC,True,quiet=True) == True 
            assert C.clonality_obvious_no_conflict_subset(NC,True,quiet=True) == set()
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert len(C.clonality_grow_no_conflict_subset(NC,count_self_conflicts=True,
                                                           more_random=more_random,N_repeats=N_repeats, quiet=True)) == 1

    def test__no_self_conflicts__allowed_changes_nonzero__1(self):
        """ Test cases with 4-5-bit numbers that show conflicts with nonzero allowed changes only; self-conflicts ignored.
        I've already tested count_self_conflicts and remove_all_zero_codeword, so using False for both here."""
        # defining the binary codewords so I can use them to check that the results are right:
        #   since we're not looking at self-conflicts here, need at least 3 different codewords;
        #   since we're not interested in cases where 0 changes already give a conflict, we can't use 2-bit numbers. 
        [b0001,b0010,b0011,b0111,b1111] = [Binary_codeword(x) for x in ['0001','0010','0011','0111','1111']]
        ### 0001, 0010 -> 0011; expect a clonality conflict with 0011 always (0111 with 1 allowed change, 1111 with 2, ...)
        #   (since self-conflicts are excluded, there's no need to worry about 0010|0111=0111 etc.)
        D = Binary_code(4,[b0001,b0010,b0011])
        full_set = frozenset([b0001,b0010,b0011])
        # with 0 allowed 0->1 changes, one conflict
        for NC in [0,(0,0),(1,0),(2,0),(9,0)]:
            conflicts = set([(frozenset([b0001,b0010]),b0011,frozenset([b0011]),'',NC)])
            assert D.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) ==({1:full_set},conflicts)
            assert D.clonality_conflict_check(NC,False,quiet=True) == True 
        # with multiple allowed 0->1 changes, three conflicts, because 0011+0010=0011 would conflict with 0001, etc
        for NC in [1,(1,1),4,(3,3),5,10]:
            assert D.clonality_count_conflicts(NC,False,return_conflict_details=False,quiet=True) == {3:full_set}
            assert D.clonality_conflict_check(NC,False,quiet=True) == True 
        # regardless of whether there's one or three conflicts, just grabbing the no-conflict subset gives an empty set, 
        #  and growing the subset until a first conflict is met gives a set of any two codewords.
        for NC in [0,(0,0),(1,0),(2,0),(9,0),1,(1,1),4,(3,3),5,10]:
            assert D.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == set()
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert len(D.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                           more_random=more_random,N_repeats=N_repeats, quiet=True)) == 2
            # giving D.clonality_grow_no_conflict_subset a starting_subset of two codewords gives same subset as result
            for subset in itertools.combinations(D.codewords,2):
                for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                    assert D.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False, starting_subset=subset, 
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True) == set(subset)
            # giving it a starting_subset of all three codewords gives error due to conflicts in starting_subset
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                self.assertRaises(BinaryCodeError,D.clonality_grow_no_conflict_subset,NC,count_self_conflicts=False, 
                              starting_subset=D.codewords, more_random=more_random,N_repeats=N_repeats, quiet=True)
        E = Binary_code(4,[b0001,b0010,b0011,b1111])
        # note that if we add b1111 to the set, with no allowed 1->0 changes that doesn't conflict with anything
        for NC in [0,(0,0),(1,0)]:
            assert E.clonality_conflict_check(NC,False,quiet=True) == True 
            assert E.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == set([b1111])
            # and the grown subset can contain b1111 and any of the other two.
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                grown_subset = E.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                               more_random=more_random,N_repeats=N_repeats, quiet=True)
            assert (b1111 in grown_subset) and len(grown_subset)==3
        ### 0001, 0010 -> 0011; expect a clonality conflict with 0111 with 1 allowed 1->0 change
        F = Binary_code(4,[b0001,b0010,b0111])
        full_set = frozenset([b0001,b0010,b0111])
        # with no 1->0 changes allowed (and up to one 0->1 change), there are no conflicts
        for NC in [0,(0,0),(0,1)]:
            assert F.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) == ({0:full_set},set())
            assert F.clonality_conflict_check(NC,False,quiet=True) == False 
            assert F.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == full_set
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert F.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True) == full_set
        # with one or more 1->0 changes allowed, there's one conflict (0001+0010=0011 with 0111)
        for NC in [1,(1,0),(2,0),(1,1)]:
            conflicts = set([(frozenset([b0001,b0010]),b0011,frozenset([b0111]),'',NC)])
            assert F.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) ==({1:full_set},conflicts)
            assert F.clonality_conflict_check(NC,False,quiet=True) == True 
        # with 2+ 0->1 changes AND 1+ 1->0 changes, all three possible conflicts (0111+0001=0111 conflict with 0010 etc)
        for NC in [2,(1,2),(2,2),9,(9,9)]:
            assert F.clonality_count_conflicts(NC,False,return_conflict_details=False,quiet=True) =={3:full_set}
            assert F.clonality_conflict_check(NC,False,quiet=True) == True 
        # whether there's one or three conflicts, the obvious subset is empty 
        #  and the grown subset can contain any two of the three conflict-codewords
        for NC in [1,(1,0),(2,0),(1,1),2,(1,2),(2,2),9,(9,9)]:
            assert F.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == set()
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert len(F.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                           more_random=more_random,N_repeats=N_repeats, quiet=True)) == 2
        ### 0001, 0010 -> 0011; expect a clonality conflict with 1111 with 2 allowed 1->0 changes
        G = Binary_code(4,[b0001,b0010,b1111])
        full_set = frozenset([b0001,b0010,b1111])
        # with 0 or 1 1->0 changes allowed (and up to one 0->1 change), there are no conflicts
        for NC in [0,(0,0),(0,1),1,(1,0),(1,1)]:
            assert G.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) == ({0:full_set},set())
            assert G.clonality_conflict_check(NC,False,quiet=True) == False 
            assert G.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == full_set
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert G.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True) == full_set
        # with two or more 1->0 changes allowed, there's one conflict (0001+0010=0011 with 1111)
        for NC in [2,(2,0),(2,1),(2,2),(9,0),(9,2)]:
            conflicts = set([(frozenset([b0001,b0010]),b0011,frozenset([b1111]),'',NC)])
            assert G.clonality_count_conflicts(NC,False,return_conflict_details=True,quiet=True) ==({1:full_set},conflicts)
            assert G.clonality_conflict_check(NC,False,quiet=True) == True 
        # with 3+ 0->1 changes AND 2+ 1->0 changes, all three possible conflicts (0111+0001=1111 conflict with 0010 etc)
        for NC in [3,(2,3),9,(9,9)]:
            assert G.clonality_count_conflicts(NC,False,return_conflict_details=False,quiet=True) =={3:full_set}
            assert G.clonality_conflict_check(NC,False,quiet=True) == True 
        # whether there's one or three conflicts, the obvious subset is empty 
        #  and the grown subset can contain any two of the three conflict-codewords
        for NC in [2,(2,0),(2,1),(2,2),(9,0),(9,2),3,(2,3),9,(9,9)]:
            assert G.clonality_obvious_no_conflict_subset(NC,False,quiet=True) == set()
            for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                assert len(G.clonality_grow_no_conflict_subset(NC,count_self_conflicts=False,
                                                   more_random=more_random,N_repeats=N_repeats, quiet=True)) == 2
        # see experiments/generating_library/1110_clonality_check_troubleshooting/notes.txt for more tests/notes

    def test__no_self_conflicts__allowed_changes_nonzero__2(self):
        """ Test cases with 3-bit numbers that show conflicts with nonzero allowed changes only; self-conflicts ignored.
        I've already tested count_self_conflicts and remove_all_zero_codeword, so using False for both here."""
        # derived from old version, with return_conflict_details=False 
        # defining the binary codewords so I can use them to check that the results are right:
        #   since we're not looking at self-conflicts here, need at least 3 different codewords;
        #   since we're not interested in cases where 0 changes already give a conflict, we can't use 2-bit numbers. 
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        [b001,b010,b100,b111] = [Binary_codeword(x) for x in ['001','010','100','111']]
        data_and_outputs = []
        D = Binary_code(3,[b110,b101,b011])
        data_and_outputs.append((D,0, {0:set([b110,b101,b011])}, (3,(3,3)) ))
        data_and_outputs.append((D,1, {3:set([b110,b101,b011])}, (2,(2,2)) ))
        data_and_outputs.append((D,(0,0), {0:set([b110,b101,b011])}, (3,(3,3)) ))
        data_and_outputs.append((D,(1,0), {0:set([b110,b101,b011])}, (3,(3,3)) ))
        data_and_outputs.append((D,(0,1), {3:set([b110,b101,b011])}, (2,(2,2)) ))
        F = Binary_code(3,[b001,b010,b100,b111])
        data_and_outputs.append((F,0, {0:set([b001,b010,b100,b111])}, (4,(4,4)) ))
        data_and_outputs.append((F,1, {2:set([b001,b010,b100]),3:set([b111])}, (3,(2,3)) ))
        data_and_outputs.append((F,(0,0), {0:set([b001,b010,b100,b111])}, (4,(4,4)) ))
        data_and_outputs.append((F,(1,0), {2:set([b001,b010,b100]),3:set([b111])}, (3,(2,3)) ))
        data_and_outputs.append((F,(0,1), {0:set([b001,b010,b100,b111])}, (4,(4,4)) ))
        for (code,N_changes,result,no_conflict_subset_len) in data_and_outputs:
            assert code.clonality_count_conflicts(N_changes,False,return_conflict_details=False,quiet=True) == result 
            expected_check_outcome = (False if result.keys()==[0] else True)
            assert code.clonality_conflict_check(N_changes,False,quiet=True) == expected_check_outcome
            if 0 in result:     expected_obvious_subset = result[0]
            else:               expected_obvious_subset = set()
            assert code.clonality_obvious_no_conflict_subset(N_changes,False,quiet=True) == expected_obvious_subset
            # compare the clonality_grow_no_conflict_subset sizes with the precomputed values in the original dataset
            #  (with more_random=False the sizes are always the same; with more_random=True the outcome sometimes depends
            #   on whether b111 is tried first or not, so I gave the min and max) 
            for N_repeats in [1,5,30]:
                assert len(code.clonality_grow_no_conflict_subset(N_changes,count_self_conflicts=False,
                                      more_random=False,N_repeats=N_repeats, quiet=True)) == no_conflict_subset_len[0]
                subset_with_more_random = code.clonality_grow_no_conflict_subset(N_changes,count_self_conflicts=False,
                                                                      more_random=True,N_repeats=N_repeats, quiet=True)
                assert no_conflict_subset_len[1][0] <= len(subset_with_more_random) <= no_conflict_subset_len[1][1]

            # and if enough repeats are run, we should get both the min and the max (THIS MAKES SURE RANDOMNESS WORKS)
            subset_more_random,all_sizes = code.clonality_grow_no_conflict_subset(N_changes,count_self_conflicts=False,
                                  more_random=True,N_repeats=50, return_repeat_summary=True, quiet=True)                
            assert len(all_sizes)==50
            assert min(all_sizes) == no_conflict_subset_len[1][0] and max(all_sizes) == no_conflict_subset_len[1][1]
            assert len(subset_more_random) == max(all_sizes)

            # in the two cases with b111 and a non-obvious solution, if we start with b111 in starting_subset, 
            #  the result will have only 2 codewords regardless of more_random
            if b111 in code.codewords and no_conflict_subset_len[0]==3:
                for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
                    assert len(code.clonality_grow_no_conflict_subset(N_changes,count_self_conflicts=False,
                               more_random=more_random,N_repeats=N_repeats, quiet=True, starting_subset=set([b111]))) == 2

    def test__clonality_grow_no_conflict_subset(self):
        """ Extra clonality_grow_no_conflict_subset checks with larger real codes - doesn't check the results compared
        to ones calculated by hand, just makes sure they're internally consistent."""
        try:
            B3 = Binary_code(3,val='error-correcting_codes/3-3-1_list',method='listfile',expected_count=2**3)
            B4 = Binary_code(4,val='error-correcting_codes/4-3-2_list',method='listfile',expected_count=2**3)
            B11 = Binary_code(11,val='error-correcting_codes/11-7-3_generator',method='matrixfile',expected_count=2**7)
        except IOError: sys.exit("Couldn't find input file in error-correcting_codes/ folder to run clonality test.")
        NC_full_list = [0,(0,0),(0,1),1,(1,0),(2,0),(1,1),2,(1,2),3]
        NC_short_list = [0,(0,1),1]
        for code, slow_test in [(B3,False), (B4,False), (B11,True)]:
            # how many NC values and other things to try depends on how slow the test is - B11 takes several seconds
            if slow_test:
                NC_list, addition_tries, random_and_repeats = NC_short_list, 2, [(False,30),(True,1)]
            else:
                NC_list, addition_tries, random_and_repeats = NC_full_list, 10, itertools.product([False,True],[1,5,30])
            # actually do the test with the given settings
            for N_changes in NC_list:
                for more_random, N_repeats in random_and_repeats:
                    no_conflict_subset = code.clonality_grow_no_conflict_subset(N_changes, count_self_conflicts=False, 
                                                                                remove_all_zero_codeword=True, quiet=True)
                    no_conflict_subcode = Binary_code(code.length, no_conflict_subset, method='list')
                    # make sure the generated subset is in fact no-conflict
                    assert no_conflict_subcode.clonality_conflict_check(N_changes,False,quiet=True) == False
                    # make sure no codeword could be added to the generated subset without introducing conflict
                    #  (only try out a few, it'd take forever otherwise!)
                    for i in range(addition_tries):
                        additional_codeword = random.choice(list(code.codewords - no_conflict_subset))
                        new_conflict_subset = no_conflict_subset | set([additional_codeword])
                        new_conflict_subcode = Binary_code(code.length, new_conflict_subset, method='list')
                        assert new_conflict_subcode.clonality_conflict_check(N_changes,False,quiet=True) == True

    def test__clonality_grow_no_conflict_subset__bad_starting_subset(self):
        """ Error should be raised if starting subset isn't conflict-free or isn't part of the full code. """
        [b110,b101,b011,b000] = [Binary_codeword(x) for x in ['110','101','011','000']]
        [b001,b010,b100,b111] = [Binary_codeword(x) for x in ['001','010','100','111']]
        for more_random, N_repeats in itertools.product([False,True], [1,5,30]):
            # no error if the starting subset is part of the code; error otherwise
            B = Binary_code(3,[b110,b101,b011,b000])
            B.clonality_grow_no_conflict_subset(0, count_self_conflicts=False, remove_all_zero_codeword=True, 
                                                starting_subset=set([b011]), quiet=True)
            self.assertRaises(BinaryCodeError, B.clonality_grow_no_conflict_subset, 0, count_self_conflicts=False, 
                              remove_all_zero_codeword=True, starting_subset=set([b001]), quiet=True)
            # no error if the starting subset is conflict-free with 1 change; error otherwise
            D = Binary_code(3,[b110,b100,b111,b010])
            B.clonality_grow_no_conflict_subset(1, count_self_conflicts=False, remove_all_zero_codeword=True, 
                                                starting_subset=set([b110,b101]), quiet=True)
            self.assertRaises(BinaryCodeError, D.clonality_grow_no_conflict_subset, 1, count_self_conflicts=False, 
                              remove_all_zero_codeword=True, starting_subset=set([b110,b101]), quiet=True)
            # including and not removing the all-zero codeword: error if counting self-conflicts, no error otherwise.
            B = Binary_code(3,[b110,b101,b011,b000])
            B.clonality_grow_no_conflict_subset(1, count_self_conflicts=False, remove_all_zero_codeword=False, 
                                                starting_subset=set([b110,b000]), quiet=True)
            self.assertRaises(BinaryCodeError, D.clonality_grow_no_conflict_subset, 1, count_self_conflicts=True, 
                              remove_all_zero_codeword=False, starting_subset=set([b110,b000]), quiet=True)


if __name__=='__main__':
    """ If module is ran directly, run tests. """

    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
