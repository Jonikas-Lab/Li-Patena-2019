#!/usr/bin/env python2.7
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011
"""

from __future__ import division
# basic libraries
import sys, os, re
import unittest
from collections import defaultdict
import itertools
import copy
import random
from math import isnan
# other libraries
from numpy import median, round, isnan
import HTSeq
from BCBio import GFF
# my modules
from general_utilities import split_into_N_sets_by_counts, add_dicts_of_ints, sort_lists_inside_dict, invert_listdict_nodups, keybased_defaultdict, value_and_percentages, FAKE_OUTFILE, NaN, nan_func, merge_values_to_unique, pickle, unpickle, write_header_data
from basic_seq_utilities import SEQ_ENDS, SEQ_STRANDS, SEQ_DIRECTIONS, SEQ_ORIENTATIONS, name_seq_generator_from_fasta_fastq, position_test_contains, position_test_overlap, chromosome_sort_key, get_seq_count_from_collapsed_header
from deepseq_utilities import check_mutation_count_by_optional_NM_field, get_HTSeq_optional_field, Fake_deepseq_objects
from parse_annotation_file import get_all_gene_annotation

# TODO TODO TODO add option to get info about adjacent genes for intergenic mutants!  Need to figure out how to structure that...

# MAYBE-TODO it might be good to split this file into multiple files at some point?  At least Insertion_position/etc.

############################ Constants and simple help functions ###########################

### Constants

RELATIVE_READ_DIRECTIONS = 'inward outward'.split()
MAX_POSITION_DISTANCE = 3000
MULTIPLE_GENE_JOIN = ' & '
MULTIPLE_mRNA_JOIN = '|'

class SPECIAL_POSITIONS(object):
    unaligned = "UNALIGNED"
    multi_aligned = "MULTIPLE"
    unknown = "UNKNOWN"
# it seems like I have to set this afterward because I can't access __dict__ from inside the class
SPECIAL_POSITIONS.all_undefined = [value for (name,value) in SPECIAL_POSITIONS.__dict__.items() if not name.startswith('__')]

class SPECIAL_GENE_CODES(object):
    not_determined = "gene_unknown"
    chromosome_not_in_reference = "unknown_chrom"
    not_found = "no_gene_found"
SPECIAL_GENE_CODES.all_codes = [value for (name,value) in SPECIAL_GENE_CODES.__dict__.items() if not name.startswith('__')]

GENE_FEATURE_NAMES = { 'five_prime_UTR' : "5'UTR", 
                       'three_prime_UTR' : "3'UTR",
                     }

GENE_FEATURE_ORDER = defaultdict(lambda: 5, { 'CDS':0, 'intron': 1, 
                                              'five_prime_UTR': 2, "5'UTR": 2, 
                                              'three_prime_UTR': 3, "3'UTR": 3, 
                                              "5'UTR_intron": 4, "3'UTR_intron": 4, 
                                              'boundary': 6, 
                                             })

### Simple help functions

class MutantError(Exception):
    """ Exceptions in the mutant_analysis_classes module; no special behavior."""
    pass

def is_cassette_chromosome(chromosome_name):
    """ Returns True if chromosome_name sounds like one of our insertion cassettes, False otherwise. """
    return ("insertion_cassette" in chromosome_name)

def is_other_chromosome(chromosome_name):
    """ Returns True if chromosome_name is neither cassette nor chromosome/scaffold, False otherwise. """
    if is_cassette_chromosome(chromosome_name):                                     return False
    if chromosome_name.startswith('chr') or chromosome_name.startswith('scaffold'): return False
    else:                                                                           return True

# MAYBE-TODO may want to put the is_*_chromosome functionality under command-line user control someday?  Provide option to give list of cassette chromosome names and "other" ones (or cassette ones and prefixes for genomic ones?)

def check_valid_position_tuple(pos):
    """ Takes a (chrom, start_pos, end_pos, strand) tuple - raises MutantError if it's wrong.

    Start_pos and end_pos should be int, 1-based, inclusive (so in AATTGG, the position of AA is 1-2) - unlike in HTSeq!
    Strand should be +/-.  No checks are done on chrom.
    """
    try:                                chrom, start_pos, end_pos, strand = pos
    except (TypeError, ValueError):     raise MutantError("Didn't get a correct position tuple! %s"%pos)
    if strand not in SEQ_STRANDS:       raise MutantError("Invalid strand %s!"%strand)
    if start_pos < 1:                   raise MutantError("Sequence positions must be positive!")
    if start_pos > end_pos:             raise MutantError("Sequence start can't be after end!")

def HTSeq_pos_to_tuple(HTSeq_pos):
    """ Convert an HTSeq.GenomicPosition instance to a (chrom,start_pos,end_pos,strand) tuple. 
    
    Start_pos and end_pos are 1-based, inclusive (so in AATTGG, the position of AA is 1-2) - unlike in HTSeq!
    """
    try:
        chrom = HTSeq_pos.chrom
    except AttributeError:
        raise MutantError("Invalid position %s! Need an HTSeq iv object. (If empty, maybe read wasn't aligned?)"%(HTSeq_pos,))
    strand = HTSeq_pos.strand
    # HTSeq is 0-based and I want 1-based, thus the +1; end has no +1 because in HTSeq end is the base AFTER the alignment.
    start_pos = HTSeq_pos.start+1
    end_pos = HTSeq_pos.end
    output_pos = (chrom, start_pos, end_pos, strand)
    check_valid_position_tuple(output_pos)
    return output_pos

def parse_flanking_region_aln_or_pos(flanking_region_aln_or_pos):
    """ Return (chrom, start_pos, end_pos, strand) tuple. 

    Input can be same tuple, then just pass through; 
     or an HTSeq alignment object, then if aligned, figure out the tuple and return it, 
      or if unaligned, return SPECIAL_POSITIONS.unaligned or SPECIAL_POSITIONS.multi_aligned depending on optional XM field value.
    """
    try:                                
        check_valid_position_tuple(flanking_region_aln_or_pos)
        return flanking_region_aln_or_pos
    except MutantError:
        try:    
            pos = flanking_region_aln_or_pos.iv
        except AttributeError:
            raise MutantError("parse_flanking_region_aln_or_pos input should be HTSeq aln or position tuple! "
                             +"Got %s"%(flanking_region_aln_or_pos,))
        if pos:                     return HTSeq_pos_to_tuple(pos) 
        # if unaligned, figure out if unaligned or multi-aligned, and just return the appropriate special position code
        else:   
            try:                    XM_val = get_HTSeq_optional_field(flanking_region_aln_or_pos, 'XM')
            except KeyError:        return SPECIAL_POSITIONS.unaligned
            if int(XM_val) > 1:     return SPECIAL_POSITIONS.multi_aligned
            else:                   return SPECIAL_POSITIONS.unaligned

def check_valid_end_info(cassette_end, relative_read_direction):
    if cassette_end not in SEQ_ENDS:
        raise MutantError("cassette_end argument must be one of %s."%SEQ_ENDS)
    if relative_read_direction not in RELATIVE_READ_DIRECTIONS:
        raise MutantError("relative_read_direction argument must be one of %s."%RELATIVE_READ_DIRECTIONS)

###################################### Insertion position functions/classes #####################################

# has to be a new-style object-based class due to the immutability/hashability thing
# MAYBE-TODO implement fuzzy positions?
class Insertion_position(object):
    """ A descriptor of the position of a genomic insertion, with separate before/after sides; optionally immutable.

    Attributes: 
        - chromosome, strand - the chromosome the insertion is on, and the strand it's in sense orientation to.
        - position_before, position_after - positions before and after the insertion site (1-based): 
                                            integers, or None if unknown.
        - min_position and max_position - lowest/highest possible position values as plain numbers, no ambiguity
        - full_position - string describing the full position: 3-4 for exact positions, 3-? or 4-? if one side is unknown. 
        Note that the three *_position attributes are really property-decorated methods, and cannot be assigned to.

    Methods: 
        - comparison/sorting: __cmp__ is based on chromosome name/number, min_/max_position, strand,
                                   and position_before/_after, in that order - see __cmp__ docstring for more detail
        - copy method returns an identical but separate copy of the object (not just another reference to the same object)
        - printing: __str__ returns a string of the chromosome,strand,full_position values
                    __repr__ returns a string of the Insertion_position() call to create a new identical object.
        Mutability and hashability: by default instances are mutable, and thus unhashable, since they implement __cmp__. 
         There are make_immutable and make_mutable_REMEMBER_CLEANUP_FIRST methods to reversibly toggle the state 
          to immutable/hashable and back. 
         This works by defining __hash__ to use the _make_key() value if immutable and raise an exception otherwise, 
          and decorating __setattr__ and __delattr__ to raise an exception if immutable and work normally otherwise.
         It's not perfect immutability, it can be gotten around by using object.__setitem__ etc, but it's good enough.
    """

    # TODO really, do I ever need those to be mutable?  I could just make a new one whenever I need to modify one, instead...  I don't think it happens often, does it?

    # TODO use __slots__  is to optimize memory usage!  This would be good to make work, but it gives me trouble with pickle!  Also it requires changes in make_mutable_REMEMBER_CLEANUP_FIRST - see TODO comments there.
    #__slots__ = ['chromosome', 'strand', 'position_before', 'position_after', 'immutable']

    # NOTE: originally I used biopython SeqFeature objects for position_before and position_after (specifically SeqFeature.ExactPosition(min_val) for exactly positions and SeqFeature.WithinPosition(min_val, min_val-max_val) for ambiguous ones), but then I realized I'm not using that and it's over-complicated and hard to make immutable and may not be what I need anyway even if I do need immutable positions, so I switched to just integers. The function to generate those was called make_position_range, and I removed it on 2012-04-23, if I ever want to look it up again later.

    def __init__(self, chromosome, strand, full_position=None, position_before=None, position_after=None, immutable=False):
        """ Initialize all values - chromosome/strand are just copied from arguments; positions are more complicated. 
        
        You must provide either full_position, OR one or both of position_before/position_after. 
        The two position_* arguments must be castable to ints, or None.
        The full_position argument must be a string of the form '100-200', '?-200' or '100-?', such as would be generated 
         by self.full_position() - self.position_before and _after are set based on the two parts of the string.
        Self.min_/max_position are calculated based on self.position_before/_after - both, or whichever one isn't None.
        If immutable is True, the object is made immutable (by calling self.make_immutable() right after initiation.
        """
        # need to make instance mutable to be able to set anything, due to how __setattr__ is decorated
        self.make_mutable_REMEMBER_CLEANUP_FIRST()  
        # now start setting attributes
        self.chromosome = chromosome
        self.strand = strand
        # parse full_position if provided
        if full_position is not None:
            if (position_before is not None) or (position_after is not None):
                raise ValueError("If providing full_position, cannot also provide position_before/position_after!")
            self.position_before, self.position_after = self._parse_full_position(full_position)
        # otherwise use position_before and/or position_after
        else:
            if position_before is None and position_after is None:
                raise ValueError("Can't create an Insertion_position object with no known position values!")
            try:
                self.position_before = None if position_before is None else int(position_before)
                self.position_after = None if position_after is None else int(position_after)
            except TypeError:  
                raise ValueError("position_before/position_after must be int-castable or None!")
        if immutable:   self.make_immutable()

    @property   # this is a builtin decorator to make an attribute out of a method
    def min_position(self):
        if self.position_after is None:     return self.position_before
        elif self.position_before is None:  return self.position_after-1
        else:                               return min(self.position_before, self.position_after-1)

    @property   # this is a builtin decorator to make an attribute out of a method
    def max_position(self):
        if self.position_after is None:     return self.position_before+1
        elif self.position_before is None:  return self.position_after
        else:                               return min(self.position_before+1, self.position_after)

    @property   # this is a builtin decorator to make an attribute out of a method
    def full_position(self):
        info_before = str(self.position_before) if self.position_before is not None else '?'
        info_after = str(self.position_after) if self.position_after is not None else '?'
        return info_before + '-' + info_after

    @classmethod
    def _parse_full_position(cls, full_position_string):
        """ Parse a full_position string to proper (position_before, position_after) value. """
        try:
            before,after = [cls._parse_single_position(s) for s in full_position_string.split('-')]
        except (ValueError,AttributeError):
            raise ValueError("The full_position argument must be a string of the form '100-200', '?-200' or '100-?'!"
                             "Got '%s'"%(full_position_string,))
        if before is None and after is None:
            raise ValueError("At least one section of the full_position argument must be a number!")
        return before,after

    @staticmethod
    def _parse_single_position(pos_string):
        """ Make a proper position value: cast to int, or return None if '?' or ''. """
        if pos_string in ['?','']:  return None
        else:                       return int(pos_string)

    def __str__(self):
        """ Return short summary of important info. """
        return ' '.join([self.chromosome, self.strand, self.full_position])

    def __repr__(self):
        """ Return full object-construction call (not very complicated) as a string. """
        return "Insertion_position('%s', '%s', full_position='%s', immutable=%s)"%(self.chromosome, self.strand, 
                                                                                   self.full_position, self.immutable)

    def copy(self):
        """ Return a deep-copy of self - NOT just a reference to the same object. """
        return Insertion_position(self.chromosome, self.strand, position_before=self.position_before, 
                                  position_after=self.position_after)

    def _make_key(self):
        """ Make key for sorting/comparison - based on chromosome/position/strand, with improved chromosome-number sorting.

        First two fields are chromosome data - splits chromosome into name/number (both optional), 
         so that 'chr2' sorts before 'chr12' (but 'chr' before 'chr1', and 'other_chr1' after 'chr4'), and also so that chromosomes
         sort first, then other names (cassette, chloroplast/mitochondrial, anything else), then scaffolds.
        Next two fields are min_/max_position - these are always numerically defined, so ?-101 and 100-? will sort together
         (as opposed to if we used position_before/_after, which can be None).
        Next field is strand - we want the sorting on position BEFORE strand, it's more readable/sensible that way.
        Final two fields are position_before/after, to ensure ?-101 isn't considered equal to 100-101.
        """
        all_position_values = (chromosome_sort_key(self.chromosome), self.min_position, self.max_position, 
                               self.strand, self.position_before, self.position_after)
        return all_position_values

    def __cmp__(self,other):
        """ Based on tuple-comparison of (chr_name, chr_number, min_pos, max_pos, strand) - in that order for sane sort.
        
        If other isn't an Insertion_position instance, assume self < other, I suppose (all we really needs is self != other.
        """
        try:
            other_key = other._make_key()
        except AttributeError:
            return -1
        return cmp(self._make_key(), other_key)

    # MAYBE-TODO do I also want a rough_comparison method, which would return True for ?-101 and 100-101 etc?  How about 100-101 and 100-102, then?

    # MAYBE-TODO add some kind of merge function to merge two positions into one?  Do I actually need that?  When I'm merging two position-based read groups together because they're actually probably the same mutant with some sequencing indels, I actually just keep the more common position, since that's presumably the correct one. So are there any cases where I'd need to merge positions?

    ### MUTABILITY AND HASHABILITY SWITCHES
    # Sometimes I want to use positions as dictionary keys or put them in sets - so they need to be hashable, 
    #   and since I also need a sane comparison operator, I can't use the default object id-based hashing, 
    #   so I have to make the objects immutable for them to be hashable. 
    # More info on how/why that is: http://docs.python.org/reference/datamodel.html#object.__hash__ 
    # This implementation is not perfectly immutable: you can get around the "immutability" by tricks like in 
    #   make_mutable_REMEMBER_CLEANUP_FIRST, but it's enough to make it clear to the user that it shouldn't be changed.
    # Some links on implementation: http://stackoverflow.com/questions/9997176/immutable-dictionary-only-use-as-a-key-for-another-dictionary, http://stackoverflow.com/questions/1151658/python-hashable-dicts, http://stackoverflow.com/questions/4996815/ways-to-make-a-class-immutable-in-python, http://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python

    def make_immutable(self):
        """ Reversibly make object immutable (reasonably) and hashable. """
        # just set the flag to make object immutable and hashable
        self.immutable = True

    def make_mutable_REMEMBER_CLEANUP_FIRST(self):
        """ Reversibly make object mutable and non-hashable. REMEMBER TO REMOVE SELF FROM SETS/DICTS BEFORE CALLING! """
        # UNSET the flag to make object immutable and hashable - need to do it in a roundabout way,
        #  because the immutability prevents simply "self.immutable = False" from working!
        self.__dict__['immutable'] = False
        # but if I put __slots__ in, self.__dict__ won't exist any more... TODO Options for then:
        # setattr(self, 'immutable', False)  -  doesn't seem to work?
        # object.__setattr__(self, 'immutable', False)  -   does that work?

    def __hash__(self):
        """ If self.hashable is True, use private _hash method, otherwise raise exception. """
        if  self.immutable:
            return hash(self._make_key())
        else:
            raise MutantError("This %s is currently mutable, and therefore unhashable! "%repr(self)
                              +"Run self.make_immutable() to change this.")

    def exception_if_immutable(function_to_wrap):
        """ Decorator: raise MutantError if self.immutable, else call function as normal. """
        def wrapped_function(self, *args, **kwargs):
            if self.immutable:
                raise MutantError("This %s is currently immutable, cannot change values! "%repr(self)
                                +"You can run self.make_mutable_REMEMBER_CLEANUP_FIRST() to change this - first make SURE "
                                +"to remove it from any sets or dictionary keys that are relying on it being hashable!")
            else:
                return function_to_wrap(self, *args, **kwargs)
        return wrapped_function

    # apply the exception_if_immutable decorator to all methods incompatible with immutability
    __setattr__ = exception_if_immutable(object.__setattr__)
    __delattr__ = exception_if_immutable(object.__delattr__)

def get_position_distance(pos1, pos2, ignore_strand=False):
    """ Return distance between positions, or NaN if on different chromosomes or on different strands and ignore_strand is False.
    """
    NaN = float('nan')
    if pos1 in SPECIAL_POSITIONS.all_undefined:             return NaN
    elif pos2 in SPECIAL_POSITIONS.all_undefined:           return NaN
    elif pos1.chromosome != pos2.chromosome:                return NaN
    elif not ignore_strand and pos1.strand != pos2.strand:  return NaN
    else:                                                   return abs(pos1.min_position - pos2.min_position)

def get_insertion_pos_from_flanking_region_pos(flanking_region_aln_or_pos, cassette_end, relative_read_direction, 
                                               immutable_position=True):
    """ Return a Insertion_position instance giving the cassette insertion position based on HTSeq read position. 

    Flanking_region_aln_or_pos should be a HTSeq.Alignment instance or a (chrom,start_pos,end_pos,strand) tuple
      (the tuple should have 1-based end-inclusive positions, so AA is 1-2 in AATT; HTSeq positions are 0-based end-exclusive); 
     cassette_end gives the side of the insertion the read is on; relative_read_direction should be inward/outward from the cassette.

    The cassette chromosome will be the same as read chromosome; the cassette strand will be either the same as read strand, 
     or opposite of the read strand. It's the same in two relative_read_direction+cassette_end combinations: 
      if the read is inward to the cassette on the 5' end, or outward from the cassette on the 3' end.  Otherwise it's opposite.

    The cassette position depends on read position, cassette strand (not read strand) and cassette_end in a complex way:
     Data I have:  which end of the insertion cassette the read is on, and which orientation the cassette is in. 
     Data I want:  the position of the base before and after the insertion, regardless of cassette orientation.
                       (the read orientation in regard to cassette doesn't matter at all here)
     If read is 5' of cassette and cassette is +, or read is 3' of cassette and cassette is -, read is BEFORE cassette
     If read is 3' of cassette and cassette is +, or read is 5' of cassette and cassette is -, read is AFTER cassette
      If read is before cassette, I care about the end of the read; if it's after, I care about the start)
      SAM alignment position is leftmost/rightmost, i.e. end is always the "later" position in the genome, 
       regardless of the read orientation, which gives me what I want, i.e. the insertion position.
    Insertion_position uses a 1-based position system (as opposed to HTSeq, which is 0-based).

    If immutable_position is True, the position will be made immutable after creation (this is reversible).
    """
    # check that basic values aren't weird
    check_valid_end_info(cassette_end, relative_read_direction)
    # parse flanking_region_aln_or_pos arg - it'll either return a tuple with the basics, or a special position code
    parsed_position = parse_flanking_region_aln_or_pos(flanking_region_aln_or_pos)
    try:                                chrom, start_pos, end_pos, strand = parsed_position
    except (TypeError, ValueError):     return parsed_position
    check_valid_position_tuple(parsed_position)
    ### chromosome is always the same as read, so just leave it as is
    ### cassette strand is the same as read strand, OR the opposite if the read is opposite to cassette (happens in two cases)
    if (cassette_end=='5prime' and relative_read_direction=='inward'):      pass
    elif (cassette_end=='3prime' and relative_read_direction=='outward'):    pass
    else:                                                                   strand = ('+' if strand=='-' else '-')
    ### cassette position depends on the read position and cassette_end in a somewhat complex way (see docstring)
    if (cassette_end=='5prime' and strand=='+') or (cassette_end=='3prime' and strand=='-'): pos_before, pos_after = end_pos, None
    else:                                                                                    pos_before, pos_after = None, start_pos
    return Insertion_position(chrom, strand, position_before=pos_before, position_after=pos_after, immutable=immutable_position)

def get_RISCC_pos_from_read_pos(read_aln_or_pos, cassette_end, relative_read_direction='inward', immutable_position=True):
    """ Return an Insertion_position instance giving the position of the far end of RISCC genome-side read.

    The output strand should be the same as that of the cassette, given relative_read_direction inward/outward from the cassette.

    See get_insertion_pos_from_flanking_region_pos for how the inputs work, and how this is done for insertion positions.  
    This is essentially calculating an "insertion position" from the OTHER SIDE of this read (so if the read goes inward
     toward the cassette, it calculates the position as if it went away from the cassette, and vice versa), 
     and then reversing the strand to make it match the strand of the real cassette.
    """
    check_valid_end_info(cassette_end, relative_read_direction)
    imaginary_relative_direction= ('outward' if relative_read_direction=='inward' else 'inward')
    imaginary_cassette_position = get_insertion_pos_from_flanking_region_pos(read_aln_or_pos, cassette_end, 
                                                                             imaginary_relative_direction)
    if imaginary_cassette_position in SPECIAL_POSITIONS.all_undefined:
        return imaginary_cassette_position
    real_strand = ('-' if imaginary_cassette_position.strand=='+' else '+')
    return Insertion_position(imaginary_cassette_position.chromosome, real_strand, 
                              full_position=imaginary_cassette_position.full_position, immutable=immutable_position)


def _get_insertion_info(insertion_pos, allowed_strand_vals=SEQ_STRANDS):
    """ Transform Insertion_position instance into (strand, ins_start, ins_end) tuple if it's not one already.
    """
    try:
        strand, ins_start, ins_end = insertion_pos.strand, insertion_pos.min_position, insertion_pos.max_position
    except AttributeError:
        strand, ins_start, ins_end = insertion_pos
    if allowed_strand_vals is not None:
        assert strand in allowed_strand_vals, "Strand should be %s, and is %s!"%(' or '.join(allowed_strand_vals), strand)
    return strand, ins_start, ins_end

def find_gene_by_pos_gff3(insertion_pos, chromosome_GFF_record, detailed_features=False, nearest_genes_for_intergenic=False, 
                          quiet=False):
    """ Look up insertion_pos in chromosome_GFF_record; return (gene_ID,orientation,subfeature,dist_to_edges).

    Insertion_pos is an Insertion_position instance, or a (strand, start_pos, end_pos) tuple; 
     insertion_pos is assumed to be in the chromosome given by chromosome_GFF_record (the caller should check that).

    If insertion_pos overlaps a gene in chromosome_GFF_record, 
     return a 4-tuple of geneID,
      orientation ('sense' if insertion_pos and gene are in the same direction, 'antisense' otherwise), 
      the name of the subfeature (exon/intron/UTR) Insertion_pos is in (or '?' if detailed_features is False),
      and the distance from the 5' and 3' edge of the gene.
     If Insertion_pos is on the edge of two features, subfeature will be 'X/Y'; if it's something more 
      unlikely/complicated, subfeature will be 'X/Y/Z??' and a warning will be printed to STDOUT.
     If gene has multiple mRNAs (splice variants) and the insertion is in different features depending on splice variant,
      the features will be separated with MULTIPLE_mRNA_JOIN (like "intron|exon").
    If multiple genes found, return data like this: "gene1 & gene2", "sense & antisense", "intron & CDS/4'UTR" and print warning.

    If insertion_pos is not in a gene, return ('no_gene_found', '-', '-', '-'), UNLESS nearest_genes_for_intergenic is True,
     then return a 4-tuple that's not quite the same as above: IDs of nearest upstream and downstream genes, 
      orientation of the insertion site vs the gene (upstream/downstream), 'intergenic' as the feature name, 
      and the distance from the nearest end of each gene - all fields will be ' & '-separated.

    Chromosome_GFF_record is a record generated by BCBio.GFF parser from a gff3 file (usually one file = multiple records). 
    """
    # MAYBE-TODO add gene name as well as ID?  But gff parsing doesn't seem to grab those at all - leave it for a separate function.
    # (see experiments/arrayed_library/1311_small-lib_mutant-distribution-file/notes.txt)
    # MAYBE-TODO add gene lengths?  Right now they can be inferred from the distances to gene start/end for mutants in genes,
    #  which seems sufficient.

    # get the needed information from either input format
    strand, ins_start, ins_end = _get_insertion_info(insertion_pos)

    ### Go over all the genes in the chromosome record, 
    #   and calculate the distance from the insertion (or 0 if insertion overlaps the gene)
    gene_distances = []
    for gene in chromosome_GFF_record.features:
        # for GFF positions, always add 1 to the gene/feature start, because BCBio uses 0-based and I use 1-based, 
        #  but don't add 1 to the end, because BCBio uses end-exclusive and I use end-inclusive.
        gene_start, gene_end = gene.location.start.position+1, gene.location.end.position
        if position_test_overlap(gene_start, gene_end, ins_start, ins_end):     gene_distances.append((gene, 0))
        elif nearest_genes_for_intergenic:
            if gene_end < ins_start:                                            gene_distances.append((gene, gene_end-ins_start))
            elif gene_start > ins_end:                                          gene_distances.append((gene, gene_start-ins_end))
            else:                                                               raise MutantError("Gene position confusion!")

    ### Pick genes to return - either all the genes that overlap the insertion, OR the closest gene on each side if none overlap
    if not gene_distances:
        nearest_genes = []
    elif min(abs(dist) for (gene,dist) in gene_distances) == 0:
        nearest_genes = [(gene,dist) for (gene,dist) in gene_distances if dist==0]
    elif nearest_genes_for_intergenic:
        # note that sometimes there ARE no genes on one/both sides!
        nearest_genes = []
        genes_upstream =   [(gene,dist) for (gene,dist) in gene_distances if dist < 0]
        genes_downstream = [(gene,dist) for (gene,dist) in gene_distances if dist > 0]
        if genes_upstream:      nearest_genes.append(max(genes_upstream, key = lambda (gene,dist): dist))
        if genes_downstream:    nearest_genes.append(min(genes_downstream, key = lambda (gene,dist): dist))
    else:
        nearest_genes = []

    ### Get all the data for each gene
    #   (see notes_on_GFF_parsing.txt for what a GFF3 record (chromosome_GFF_record) will be like)
    gene_data_list = []
    for (gene,dist) in nearest_genes:
        gene_start, gene_end = gene.location.start.position+1, gene.location.end.position
        gene_ID = gene.id
        # for mutants in genes, calculate orientation of insertion_pos vs gene
        if dist == 0:
            if strand=='both':          orientation = 'both'
            elif gene.strand==1:        orientation = 'sense' if strand=='+' else 'antisense'
            elif gene.strand==-1:       orientation = 'sense' if strand=='-' else 'antisense'
            else:                       orientation = '?'
        # for intergenic mutants, check whether the mutant is upstream or downstream of the gene, instead of orientation
        else:
            if dist * gene.strand < 0:  orientation = 'downstream'
            else:                       orientation = 'upstream'
        # calculate distances:
        #  - if insertion is inside gene, give distance to 5' end and 3' end, separated by a comma
        #  - if it's intergenic, just give the distance to the closest side of the gene (already calculated)
        if dist == 0:
            dist1 = ins_end - gene_start
            dist2 = gene_end - ins_start
            if gene.strand==1:  distances = '%s,%s'%(dist1,dist2)
            else:               distances = '%s,%s'%(dist2,dist1)
        else:                   distances = str(abs(dist))
        # basic features: intergenic, gene, gene edge
        if dist == 0:
            # if it overlaps an edge, note that in features_basic by adding 'gene_edge'
            # (MAYBE-TODO if I look at things like inner/outer flanking regions, this is where that would go as well)
            if position_test_contains(gene_start, gene_end, ins_start,ins_end): features_basic = []
            else:                                                               features_basic = ['gene_edge']
        else:                                                                   features_basic = ['intergenic']
        # figure out which feature of the gene the insertion is in, IF we're looking for detailed ones (it's a lot of code)
        if dist != 0:                           inner_feature = ''
        elif not detailed_features:             inner_feature = '?'
        else:
            if len(gene.sub_features)==0:       inner_features = ['no_mRNA']
            else:                               inner_features = []
            for mRNA in gene.sub_features:
                if gene.sub_features[0].type != 'mRNA':
                    if not quiet:
                        print("Warning: gene %s in gff file has unexpected non-mRNA sub-features! "%gene_ID
                              +"Returning '??' feature.")
                    inner_features.append('??')
                else:
                    mRNA_start, mRNA_end = mRNA.location.start.position+1,mRNA.location.end.position
                    # if insertion_pos is outside the mRNA, use 'outside_mRNA' as inner_feature
                    if not position_test_overlap(mRNA_start, mRNA_end, ins_start, ins_end):
                        inner_features.append('outside_mRNA')
                    # if insertion_pos is inside the mRNA and mRNA has no subfeatures, use 'mRNA_no_exons' as inner_feature
                    elif len(mRNA.sub_features)==0:   
                        if position_test_contains(mRNA_start, mRNA_end, ins_start, ins_end):  inner_features.append('mRNA_no_exons')
                        else:                                                                 inner_features.append('mRNA_edge')
                    else: 
                        # otherwise go over all subfeatures, see which ones contain/overlap insertion_pos
                        #  (check for overlap only if the contains test failed)
                        features_inside = []
                        if position_test_contains(mRNA_start, mRNA_end, ins_start, ins_end): features_edge = []
                        else:                                                                features_edge = ['mRNA_edge']
                        for feature in mRNA.sub_features:
                            feature_start, feature_end = feature.location.start.position+1, feature.location.end.position
                            try:                feature_type = GENE_FEATURE_NAMES[feature.type]
                            except KeyError:    feature_type = feature.type
                            if position_test_contains(feature_start, feature_end, ins_start, ins_end):
                                features_inside.append(feature_type)
                            elif position_test_overlap(feature_start, feature_end, ins_start, ins_end):
                                features_edge.append(feature_type)
                        # MAYBE-TODO may want to treat exons before 5'UTR or after 3'UTR specially? 
                        #   Not worth it, none in current file.
                        # if insertion_pos is inside a single mRNA subfeature, use the type of the subfeature as inner_feature
                        if len(features_inside)==1 and len(features_edge)==0:
                            inner_features.append(features_inside[0])
                        # if insertion_pos is on the edge of two mRNA subfeatures, use 'subfeature1/subfeature2'
                        elif len(features_inside)==0 and len(features_edge)==2:
                            inner_features.append('/'.join(features_edge))
                            # MAYBE-TODO treat insertions CLOSE to an edge specially too? How large is a splice junction?
                        # if insertion_pos is on the edge of ONE mRNA subfeature, or not touching any subfeatures at all, 
                        #  the implied subfeature is either an intron (if between features) or mRNA_before/after_exons, 
                        #   (which shouldn't happen in normal files).
                        elif len(features_inside)==0 and len(features_edge)<=1:
                            # figure out what the implied feature is - outside intron in CDS (normal) or UTR, or outside all exons
                            # note that before/after and 5'/3' are swapped if gene is on minus strand!
                            CDS_features = [feature for feature in mRNA.sub_features if feature.type=='CDS']
                            if ins_start < min([feature.location.start.position+1 for feature in mRNA.sub_features]):
                                if gene.strand==1:    implied_feature = 'mRNA_before_exons'
                                elif gene.strand==-1: implied_feature = 'mRNA_after_exons'
                            elif ins_end > max([feature.location.end.position for feature in mRNA.sub_features]):
                                if gene.strand==1:    implied_feature = 'mRNA_after_exons'
                                elif gene.strand==-1: implied_feature = 'mRNA_before_exons'
                            elif ins_start < min([feature.location.start.position+1 for feature in CDS_features]):
                                if gene.strand==1:    implied_feature = "5'UTR_intron"
                                elif gene.strand==-1: implied_feature = "3'UTR_intron"
                            elif ins_end > max([feature.location.end.position for feature in CDS_features]):
                                if gene.strand==1:    implied_feature = "3'UTR_intron"
                                elif gene.strand==-1: implied_feature = "5'UTR_intron"
                            else:
                                implied_feature = 'intron'
                            # set inner_feature based on whether insertion_pos is on a real/implied feature edge 
                            #  or completely inside an implied feature
                            if len(features_edge)==1:
                                inner_features.append(features_edge[0] + '/' + implied_feature)
                            elif len(features_edge)==0:
                                inner_features.append(implied_feature)
                        # if insertion_pos is inside two features, or inside one and on the edge of another, 
                        #  print a warning, and use all the feature names, with a ?? at the end to mark strangeness
                        else:
                            inner_features.append('/'.join(features_inside+features_edge) + '??')
                            if not quiet:
                                print(("Warning: Location (%s,%s) matched multiple features (%s) "
                                      +"in gene %s!")%(ins_start, ins_end, inner_features[-1], gene_ID)) 
            inner_feature = MULTIPLE_mRNA_JOIN.join(sorted(set(inner_features)))
        
        # prepend whatever gene-level features (edge etc, or []) were found at the start to the full value
        full_feature = '/'.join(features_basic + ([inner_feature] if inner_feature else []))
        gene_data_list.append([gene_ID, orientation, full_feature, distances])

    ### Return appropriate value
    # if no gene matching insertion_pos was found, return special value
    if not gene_data_list:
        return [SPECIAL_GENE_CODES.not_found, '-', '-', '-']
    # if single gene found, return its info
    elif len(gene_data_list) == 1:
        return gene_data_list[0]
    # if multiple genes found, return data like this: "gene1 & gene2", "sense & antisense", "intron & CDS/4'UTR",
    #  except that the feature for intergenic insertions should be "intergenic" instead of "intergenic & intergenic".
    else: 
        full_data = [MULTIPLE_GENE_JOIN.join(multiple_vals) for multiple_vals in zip(*gene_data_list)]
        if full_data[2] == MULTIPLE_GENE_JOIN.join(['intergenic']*2):   full_data[2] = 'intergenic'
        return full_data


def find_gene_by_pos_simple(insertion_pos, chromosome_gene_pos_dict, allow_multiple_genes=False):
    """ Look up insertion_pos in a chromosome:geneID:(start,end) dictionary; return gene_ID only.

    If allow_multiple_genes is False:
     if insertion_pos overlaps a gene in chromosome_gene_pos_dict, return geneID; otherwise return 'no_gene_found';
     raise an exception if the position overlaps multiple genes.
    Otherwise, return a list of geneID values, or ['no_gene_found'] if there were none.

    Insertion_pos is an Insertion_position instance, or a (chromosome, start_pos, end_pos) tuple.
    Chromosome_gene_pos_dict can be generated by gff_examine_file.gene_position_dict, or something else similar, 
     used to parse gene position data files that don't match the GFF3 format (use find_gene_by_pos_gff3 for that)
    """
    # get the needed information from either input format
    if isinstance(insertion_pos, Insertion_position):
        chromosome, ins_start, ins_end = insertion_pos.chromosome, insertion_pos.min_position, insertion_pos.max_position
    else:
        chromosome, ins_start, ins_end = insertion_pos
    # go over all the genes in the chromosome and look for one that matches the position in insertion_pos
    genes_matched = []
    for gene,(gene_start,gene_end) in chromosome_gene_pos_dict[chromosome].items():
        if position_test_overlap(gene_start, gene_end, ins_start,ins_end):
            genes_matched.append(gene)
    # deal with multiple/single/no matching gene cases depending on allow_multiple_genes value
    if allow_multiple_genes:
        if genes_matched:               return genes_matched
        else:                           return [SPECIAL_GENE_CODES.not_found]
    else:
        if len(genes_matched) > 1:      raise Exception("Multiple genes matched insertion position (%s %s-%s)! (%s)"%(
                                                                    chromosome, ins_start, ins_end, ', '.join(genes_matched)))
        elif len(genes_matched) == 1:   return genes_matched[0]
        else:                           return SPECIAL_GENE_CODES.not_found


######################################### Single mutant functions/classes ###############################################

# help functions returning "blank" mutants for defaultdicts (can't use lambdas because pickle doesn't like them)
def blank_readcount_only_mutant():
    return Insertional_mutant_readcount_only()
def blank_full_single_mutant():
    return Insertional_mutant()
def blank_single_mutant_with_IB(IB):
    return Insertional_mutant(IB=IB)
def blank_multi_mutant_with_IB(IB):
    return Insertional_mutant_multi_dataset(IB=IB)


class Insertional_mutant():
    """ Data regarding a particular insertional mutant: insertion position, gene, read numbers/sequences, etc.

    Mutants have the following attributes (data):
     1) Position/gene attributes (except readcount-related-only mutants, which don't have those):
       - IB - a (hopefully unique) internal barcode sequence for that mutant
       - position - an Insertion_position instance giving the insertion chromosome/strand/position
       - gene, orientation, gene feature - what gene the insertion is in (or one of the SPECIAL_GENE_CODES if unknown), 
                   whether it's in the sense or antisense orientation vs the gene, what feature (exon/intron/UTR) it's in.
     2) Readcount-related attributes (multi-dataset mutants don't have those on the top-level):
       - total_read_count, perfect_read_count - number of all and perfectly aligned deepseq reads
       - sequences_counts_positions_errors - a seq:(count, aligned_cassette_position, N_alignment_errors) dictionary
       - RISCC_genome_side_aligned_reads - contains "genome-side" reads (not immediately next to the cassette) that were 
                                    uniquely aligned to the genome/cassette. They can help determine the real insertion position 
                                    (vs. junk fragments inserted with cassette).
            It's a dictionary - position:(position, read_count, seq_read_error_dict, gene, orientation, feature, 
                distances_from_gene_ends, annotation...), where seq_read_error_dict is a seq:(read_count, N_errors) dict 
                (since multiple seqs can have the same position)
       - RISCC_genome_side_unaligned_reads - same, but contains reads that were unaligned or multiply-aligned, 
                                    so the positions are all unspecified, and sequences are used as keys instead. 

    This class also has subclasses for mutant types with additional functionality or purposefully limited functionality.

    Some notes on methods (functions) of Insertional_mutant objects:
     For detailed information on mutant methods, see method docstrings.
     Methods with names starting with _ are private and shouldn't be used from outside the object itself.
     Readcount-related methods take a dataset_name argument for compatibility with the multi-dataset subclass
      - on basic mutants it should always be None (which is the default).
    """
    # MAYBE-TODO add a __slots__ to this to optimize memory usage for when there's tens of thousands of mutants!  But that may cause issues with pickling, etc...
    # TODO make RISCC_genome_side_aligned_reads values namedtuples or something, if not objects...

    def __init__(self, IB=None, insertion_position=SPECIAL_POSITIONS.unknown):
        """ Set self.position based on argument; initialize read/sequence counts to 0 and gene-info to unknown. 

        insertion_position argument should be an Insertion_position instance. 
        """
        # "standard" mutants need general attributes and readcount-related attributes
        self._set_general_attributes(IB, insertion_position)
        self._set_readcount_related_data_to_zero()

    # MAYBE-TODO give each mutant some kind of unique ID at some point in the process?  Or is genomic location sufficient?  If we end up using per-mutant IBs (in addition to the flanking sequences), we could use that, probably, or that plus genomic location.

    def _set_general_attributes(self, IB, insertion_position=SPECIAL_POSITIONS.unknown):
        self.IB = IB
        self.position = insertion_position
        # MAYBE-TODO should I have a class for the gene data? Especially if I add more of it (annotation info etc)
        self.gene = SPECIAL_GENE_CODES.not_determined
        self.orientation = '?'
        self.gene_feature = '?'
        self.gene_distances = '?'

    def _set_readcount_related_data_to_zero(self):
        """ Set all readcount-related data to 0/empty."""
        self.total_read_count      = 0
        self.perfect_read_count    = 0
        self.RISCC_genome_side_aligned_reads = {}
        self.RISCC_genome_side_unaligned_reads = {}
        self.sequences_counts_positions_errors = {}
        # TODO should all this really be readcount-related? Well, it IS, but when I have a multi-dataset mutant, do I really want to keep the seq/position/count details and the genome-side RISCC read data per dataset rather than total?  Hard to tell, really.  In a perfect world I wouldn't be doing multiple RISCC datasets anyway!

    def read_info(self, dataset_name=None, strict=False):
        """ Help function to get read-info-containing object for both multi-dataset and single mutants.

        Here this just returns self; more complicated version is for multi-dataset mutants.
        Strict is ignored and only present to make the implementation consistent with the multi-dataset version.
        """
        if dataset_name is None:    return self
        else:                       raise MutantError("This is NOT a multi-dataset mutant - cannot provide dataset_name arg!")

    @staticmethod
    def _ensure_dataset_None(dataset_name):
        """ Raise Error if dataset_name isn't None; give explanation. 

        Many Insertional_mutant methods take a dataset_name argument for consistency with the multi-dataset subclass, 
         but for non-subclass objects this shouldn't be provided.  The only reason for having them at all is so that you get
         a more useful error message when you try to provide a dataset_name, rather than just getting an ArgumentError or such.
        """
        if dataset_name is not None:
            raise MutantError("Don't try to provide a dataset_name on a single mutant (rather than the multi-dataset subclass)!")
        # MAYBE-TODO this could be accomplished with a decorator instead, right?

    def add_read(self, HTSeq_alignment, position=SPECIAL_POSITIONS.unknown, read_count=1, dataset_name=None):
        """ Add a read to the data (or multiple identical reads, if read_count>1); return True if perfect alignment.

        Specifically: increment total_read_count, increment perfect_read_count if read is a perfect 
         alignment, increment the appropriate field of sequences_counts_positions_errors based on read sequence.

        This method does NOT check that position matches the current position of the mutant - use decide_and_check_position method 
         for that.  This also does NOT check that position is consistent with the HTSeq_alignment.iv position.

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        # TODO instead of taking HTSeq_alignment, this could just take the seq and N_errors, like add_RISCC_read does?
        self._ensure_dataset_None(dataset_name)
        # increment total_read_count, and add read ID to the ID set
        self.total_read_count += read_count
        # figure out if the read is perfect and increment perfect_read_count if yes; return True if perfect else False.
        # TODO may want to come up with a better option than 10 for the "errors" of unaligned seqs
        if position in SPECIAL_POSITIONS.all_undefined:
            N_errors = 10
        else:
            N_errors = check_mutation_count_by_optional_NM_field(HTSeq_alignment, negative_if_absent=False)
        # add sequence position/readcount data the detailed dictionary.
        seq = HTSeq_alignment.read.seq
        try:                self.sequences_counts_positions_errors[seq][0] += read_count
        except KeyError:    self.sequences_counts_positions_errors[seq] = [read_count, position, N_errors]
        if N_errors==0:  
            self.perfect_read_count += read_count
            return True
        else:
            return False

    def decide_and_check_position(self, max_allowed_dist=0, ratio_to_ignore=100, OUTPUT=None):
        """ Set self.position to be the highest-count DEFINED position; check that all positions are within max_allowed_dist of it. 

        If a defined position isn't within max_allowed_dist of the main position:
         - if its readcount is ratio_to_ignore times lower than that of the main position, 
            REMOVE IT and decrease total/perfect readcounts appropriately
         - otherwise raise an exception.
        """
        if not self.sequences_counts_positions_errors:
            self.position = SPECIAL_POSITIONS.unknown
            return
        main_seq, (main_count, main_pos, main_Nerr) = min(self.sequences_counts_positions_errors.items(), 
                                             key = lambda (s, (c, p, e)): (p in SPECIAL_POSITIONS.all_undefined, -c, e, s))
        self.position = main_pos
        for seq, (count, pos, N_err) in self.sequences_counts_positions_errors.items():
            if pos not in SPECIAL_POSITIONS.all_undefined:
                if not get_position_distance(main_pos, pos, ignore_strand=False) <= max_allowed_dist:
                    if count*ratio_to_ignore <= main_count:
                        # TODO removing these reads is a problem, because we don't remove the genome-side reads!
                        del self.sequences_counts_positions_errors[seq]
                        self.total_read_count -= count
                        if N_err==0:    self.perfect_read_count -= count
                    else:
                        if OUTPUT is not None:
                            OUTPUT.write("Warning: Different cassette-side position in same mutant! REMOVING MUTANT. IB %s,"%self.IB 
                                   +" %s %s %serr %s reads, %s %s %serr %s reads\n"%(main_pos, main_seq, main_Nerr, main_count, 
                                                                                     pos, seq, N_err, count))
                        return True

    def update_gene_info(self, gene, orientation, gene_feature, gene_distances):
        """ Update gene/orientation/feature: if both are the same or one is unknown, keep known; if different, raise error.
        """
        # grab the set of non-special values from own and new data
        gene_both = set([self.gene, gene]) - set([SPECIAL_GENE_CODES.not_determined])
        orientation_both = set([self.orientation, orientation]) - set(['?'])
        feature_both = set([self.gene_feature, gene_feature]) - set(['?'])
        # if there are two different non-special values for any of the three attributes, raise MutantError
        if len(gene_both)>1 or len(orientation_both)>1 or len(feature_both)>1:
            raise MutantError("Can't join the two mutants: the gene/orientation/feature data differs!")
        # otherwise set own data to the better one of own/new data (unless neither are present)
        if gene_both:           self.gene = gene_both.pop()
        if orientation_both:    self.orientation = orientation_both.pop()
        if feature_both:        self.gene_feature = feature_both.pop()
        if self.gene == gene:   self.gene_distances = gene_distances

    # MAYBE-TODO implement mutant-merging?  If yes, copy join_mutant method from mutant_analysis_classes.py and edit.
    #  If we merge mutants that correspond to two sides of one insertion, they'll have to have multiple IBs! 
    #   5'+3', or 5'+5' or 3'+3' in cases with tandem tail-to-tail or head-to-head cassettes, and maybe even more complicated...

    @staticmethod
    def _get_main_sequence_from_data(seqs_to_counts_and_data, N=1, aligned_only=False):
        """ Return the most common sequence in the given data, and its count (or Nth most common sequence if N is provided). """
        # use key to sort reverse by count but non-reverse by seq
        if aligned_only:
            filtered_data = [(seq, data) for (seq, data) in seqs_to_counts_and_data.items() 
                             if isinstance(data[1], Insertion_position)]
        else:
            filtered_data = seqs_to_counts_and_data.items()
        sequences_by_count = sorted([(seq,data[0]) for (seq,data) in filtered_data], 
                                    key = lambda (s,c): (-c, s))
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return tuple(sequences_by_count[N-1])
        except IndexError:  return ('',0)
        # MAYBE-TODO should probably make that '-' or something instead of '', empty strings are hard to see. 
        #  On the other hand '-' isn't a valid sequence, and '' is...

    def get_main_sequence(self, N=1, dataset_name=None, aligned_only=False):
        """ Return the most common sequence in this mutant and its count (or Nth most common sequence if N is provided).

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        self._ensure_dataset_None(dataset_name)
        return self._get_main_sequence_from_data(self.sequences_counts_positions_errors, N, aligned_only)

    def _copy_non_readcount_data(self, source_mutant):
        """ Copy non-readcount-related data from source_mutant to self (making new copies of all objects). """
        # COPY the position, not just make another name for the same value - I wrote a copy() function for positions
        self.position       = source_mutant.position.copy() 
        # strings are immutable and thus safe to "copy" by adding another name to the same value
        self.gene           = source_mutant.gene
        self.orientation    = source_mutant.orientation
        self.gene_feature   = source_mutant.gene_feature
        self.gene_distances = source_mutant.gene_distances

    def _copy_readcount_related_data(self, source_mutant):
        """ Copy readcount-related data from source_mutant to self (making new copies of all objects). """
        # integers are immutable and thus safe to "copy" by adding another name to the same value
        self.total_read_count      = source_mutant.total_read_count
        self.perfect_read_count    = source_mutant.perfect_read_count
        # using dict to make a COPY of the dict instead of just creating another name for the same value
        self.sequences_counts_positions_errors = dict(source_mutant.sequences_counts_positions_errors)

    # MAYBE-TODO should there be a copy_mutant function to make a deepcopy?

    ### RISCC-related functions
    # TODO unit-test all this!!
    # TODO add these as appropriate into multi-dataset or readcount-only versions, if desired; give them dataset_name=None args etc.

    def add_RISCC_read(self, seq, new_position, N_errors=None, read_count=1):
        """ Add new genome-side read to RISCC genome-side aligned or unaligned reads, or increment existing one if present."""
        # TODO why are we even using Insertion_position objects here?? Those aren't insertion positions with a start-end, just single positions...  But still need to be able to deal with unaligned/multi as well as proper positions.
        if not isinstance(new_position, Insertion_position) and new_position not in SPECIAL_POSITIONS.all_undefined:
            raise MutantError("RISCC read position %s is unacceptable - must be Insertion_position object or one of %s!"%(
                new_position, ', '.join(SPECIAL_POSITIONS.all_undefined)))
        # self.RISCC_genome_side_aligned_reads is a position:data dict
        if new_position not in SPECIAL_POSITIONS.all_undefined:
            try:
                # MAYBE-TODO check that the same seq isn't present in a different position?
                self.RISCC_genome_side_aligned_reads[new_position][1] += read_count
                try:                self.RISCC_genome_side_aligned_reads[new_position][2][seq][0] += read_count
                except KeyError:    self.RISCC_genome_side_aligned_reads[new_position][2][seq] = [read_count, N_errors]
            except KeyError:
                seq_count_error_dict = {seq: [read_count, N_errors]}
                self.RISCC_genome_side_aligned_reads[new_position] = [new_position, read_count, seq_count_error_dict, 
                                                                      SPECIAL_GENE_CODES.not_determined, '?', '?', '?']
        # self.RISCC_genome_side_unaligned_reads is a seq:data dict, since the positions aren't usable as keys
        else:
            try:
                self.RISCC_genome_side_unaligned_reads[seq][1] += read_count
                self.RISCC_genome_side_aligned_reads[seq][2][seq][0] += read_count
            except KeyError:
                self.RISCC_genome_side_unaligned_reads[seq] = [new_position, read_count, {seq: [read_count, N_errors]}, 
                                                               SPECIAL_GENE_CODES.not_determined, '?', '?', '?']
        # Note: adding gene/annotation info for those is implemented in the dataset methods.

    def _if_confirming_read(self, position, max_distance):
        if get_position_distance(position, self.position, ignore_strand=False) <= max_distance:  return True
        # TODO what about WEIRD cases, when it's the right distance but in the wrong direction? Those are rare enough to ignore for now, but should do something about them!

    def _decide_if_replace_read(self, new_position, max_distance):
        """ Decide whether new position is "better" than old one.
        """
        # if there are no current reads, add new one
        if not len(self.RISCC_genome_side_aligned_reads):               return True
        # if new read isn't "confirming", it can't be better
        if not self._if_confirming_read(new_position, max_distance):    return False
        # if the new one is "confirming" and the old one isn't, new one has to be better
        old_position = self.RISCC_genome_side_aligned_reads.values()[0][0]
        if not self._if_confirming_read(old_position, max_distance):    return True
        # if both the old and new position meet the basic conditions, pick the highest-distance one
        # TODO what about directionality and weird cases?
        new_dist = abs(new_position.min_position - self.position.min_position)
        old_dist = abs(old_position.min_position - self.position.min_position)
        if new_dist > old_dist:                                 return True
        else:                                                   return False

    def improve_best_RISCC_read(self, seq, new_position, N_errors=None, read_count=1, max_distance=MAX_POSITION_DISTANCE):
        """ Compare the read to the current best one - replace the best one if this is better.

        "Better" meaning furthest away from the cassette-side read, while still remaining on the same chromosome, strand,
            and within max_distance of it.
         If both reads are on a different chromosome or strand than the cassette-side position (both are bad, really), 
            the first one is kept.
         If there is no current read, the new one is always used.
        """
        # if there are more than one current reads, you're not using improve_best_RISCC_read consistently!
        if len(self.RISCC_genome_side_aligned_reads) > 1:
            raise MutantError("Don't try using the improve_best_RISCC_read when keeping more than one read!")
        # if decided to replace, discard old genome-side read dict and make new one from just the current read data.
        if self._decide_if_replace_read(new_position, max_distance):
            self.RISCC_genome_side_aligned_reads, self.RISCC_genome_side_unaligned_reads = {}, {}
            self.add_RISCC_read(seq, new_position, N_errors, read_count)
        # TODO make this count unaligned/confirming/non-confirming reads, too, instead of keeping all these counts as functions that read the actual mutant data, which will be missing in this case?  I did something like that in mutant_Carette.py.

    def RISCC_N_confirming_seqs(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of unique aligned genome-side positions that confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        N = 0
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            # skip non-aligned reads; check aligned reads for confirming.
            try:                    chrom = read_data[0].chromosome
            except AttributeError:  continue
            if self._if_confirming_read(read_data[0], max_distance):  
                N += 1
        return N

    def RISCC_N_non_confirming_seqs(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of unique aligned genome-side positions that DON'T confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        return len(self.RISCC_genome_side_aligned_reads) - self.RISCC_N_confirming_seqs(max_distance)

    def RISCC_N_confirming_reads(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of aligned genome-side READS that confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        N = 0
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            # skip non-aligned reads; check aligned reads for confirming.
            try:                    chrom = read_data[0].chromosome
            except AttributeError:  continue
            if self._if_confirming_read(read_data[0], max_distance):  
                N += read_data[1]
        return N

    def RISCC_N_non_confirming_reads(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of aligned genome-side READS that DON'T confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        return sum(x[1] for x in self.RISCC_genome_side_aligned_reads.values()) - self.RISCC_N_confirming_reads(max_distance)

    @property
    def RISCC_N_unaligned_seqs(self):
        return len(self.RISCC_genome_side_unaligned_reads)

    @property
    def RISCC_N_aligned_seqs(self):
        return len(self.RISCC_genome_side_aligned_reads)

    @property
    def RISCC_N_genomic_seqs(self):
        N = 0
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            try:                    chrom = read_data[0].chromosome
            except AttributeError:  continue
            if not is_cassette_chromosome(chrom):
                N += 1
        return N

    @property
    def RISCC_N_cassette_seqs(self):
        return self.RISCC_N_aligned_seqs - self.RISCC_N_genomic_seqs

    @property
    def RISCC_N_unaligned_reads(self):
        return sum(x[1] for x in self.RISCC_genome_side_unaligned_reads.values())

    @property
    def RISCC_N_aligned_reads(self):
        return sum(x[1] for x in self.RISCC_genome_side_aligned_reads.values())

    @property
    def RISCC_N_genomic_reads(self):
        N = 0
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            try:                    chrom = read_data[0].chromosome
            except AttributeError:  continue
            if not is_cassette_chromosome(chrom):
                N += read_data[1]
        return N

    @property
    def RISCC_N_cassette_reads(self):
        return self.RISCC_N_aligned_reads - self.RISCC_N_genomic_reads

    def RISCC_percent_confirming_seqs(self, max_distance=MAX_POSITION_DISTANCE, round_to_int=False):
        """ % of unique genome-side sequences that confirm the cassette-side position (same chrom/strand, within max_distance).
        """
        if not self.RISCC_N_aligned_seqs:   return float('nan')
        else:                               
            percent = self.RISCC_N_confirming_seqs(max_distance) / self.RISCC_N_aligned_seqs * 100
            if round_to_int:    percent = int(round(percent))
            return percent

    def RISCC_percent_confirming_reads(self, max_distance=MAX_POSITION_DISTANCE, round_to_int=False):
        """ % of genome-side READS that confirm the cassette-side position (same chrom/strand, within max_distance).
        """
        if not self.RISCC_N_aligned_reads:  return float('nan')
        else:
            percent = self.RISCC_N_confirming_reads(max_distance) / self.RISCC_N_aligned_reads * 100
            if round_to_int:    percent = int(round(percent))
            return percent

    @property
    def RISCC_N_genomic_chromosomes(self):
        chroms = set()
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            # grab chromosome - unless read is unaligned, then ignore and go on to next one
            try:                    chrom = read_data[0].chromosome
            except AttributeError:  continue
            chroms.add(chrom)
        return sum(1 for chrom in chroms if not is_cassette_chromosome(chrom))
            
    def RISCC_N_distinct_regions(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return number of distinct insertion regions implied by RISCC data (genomic, cassette, and #chromosomes).

        The output is a 3-tuple giving the number of distinct genome and cassette regions, 
         and the number of distinct non-cassette chromosomes the regions are in.
        Positions on different chromosomes/strands are always counted as distinct; positions on same chromosome/strand are 
         counted as distinct regions if the distance between them is >=max_distance (THIS IS SLIGHTLY ROUGH).
        
        Data used to generate the #regions includes the cassette-side position (single) and all the genome-side RISCC positions. 
        Unaligned reads are ignored.
        """
        # TODO add options for minimum #seqs and #reads to count a region as valid!
        positions_by_chrom_strand = defaultdict(list)
        # add the cassette-side position (single)
        try:
            positions_by_chrom_strand[(self.position.chromosome, self.position.strand)].append(self.position.min_position)
        except AttributeError:
            pass
        # add all the genome-side read positions; skip unaligned ones.
        for read_data in self.RISCC_genome_side_aligned_reads.values():
            pos = read_data[0]
            try:                    positions_by_chrom_strand[(pos.chromosome, pos.strand)].append(pos.min_position)
            except AttributeError:  continue
        # count total number of dictinct regions - different chromosomes or strands, or distance > max_distance
        total_distinct_regions_genome = 0
        total_distinct_regions_cassette = 0
        # for each chromosome, go over all positions and only count ones every MAX_POSITION_DISTANCE as distinct
        for chrom_strand, positions in positions_by_chrom_strand.items():
            positions.sort()
            distinct_regions = [positions[0]]
            for pos in positions[1:]:
                if (pos-distinct_regions[-1]) > max_distance:
                    distinct_regions.append(pos)
            if is_cassette_chromosome(chrom_strand[0]):     total_distinct_regions_cassette += len(distinct_regions)
            else:                                           total_distinct_regions_genome += len(distinct_regions)
        return total_distinct_regions_genome, total_distinct_regions_cassette

    def RISCC_max_confirmed_distance(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the distance between the cassette-side read and the furthest-away same-area genome-side read.

        Return 0 if no same-area genome-side reads, and NaN if there are no uniquely aligned genome-side reads at all.
        """
        distances = []
        if self.RISCC_N_confirming_seqs(max_distance) + self.RISCC_N_non_confirming_seqs(max_distance) == 0:
            return float('NaN')
        for RISCC_read_data in self.RISCC_genome_side_aligned_reads.values():
            # Only look at the genome-side reads that match the cassette-side read position!
            # There's a try/except because unaligned reads don't have proper positions.
            try:                    
                if (RISCC_read_data[0].chromosome == self.position.chromosome 
                    and RISCC_read_data[0].strand == self.position.strand):
                    pos_difference = abs(RISCC_read_data[0].min_position - self.position.min_position)
                else:
                    continue
            except AttributeError:  
                continue
            if pos_difference <= max_distance:
                distances.append(pos_difference)
        try:
            return max(distances)
        except ValueError:
            return 0
        # this is basically unit-tested by the tests for add_RISCC_read and improve_best_RISCC_read

    # TODO add new method that infers the approximate real insertion location in cases of pretty obvious junk fragments!

    def RISCC_print_detail(self, OUTPUT=sys.stdout, max_distance=MAX_POSITION_DISTANCE):
        """ Print RISCC detail: overall distinct position counts, cassette-side position, all genome-side reads/positions/counts.
        """
        ### print summary line
        # TODO is that a useful summary?  Should we add more?  Maybe add % confirmed etc?
        N_distinct_genome, N_distinct_cassette = self.RISCC_N_distinct_regions(max_distance)
        OUTPUT.write(" * IB %s: %s distinct genomic regions, %s cassette regions, plus %s unaligned/multi-aligned seqs."%(
            self.IB, N_distinct_genome, N_distinct_cassette, self.RISCC_N_unaligned_seqs)
                     +" Cassette-side position is RISCC-confirmed to %s bp; %s%% of the genome-side reads are in the same region"%(
                     self.RISCC_max_confirmed_distance(), self.RISCC_percent_confirming_reads(round_to_int=True))
                     +" (%s are, %s are not).\n"%(self.RISCC_N_confirming_reads(), self.RISCC_N_non_confirming_reads()))
        data_header = "(chrom strand pos gene orientation feature distances readcount perfect main_seq gene_name)".replace(' ','\t')
        ### print cassette-side position
        OUTPUT.write("Cassette-side_position::: %s\n"%data_header)
        try:                    main_pos_fields = [self.position.chromosome, self.position.strand, self.position.full_position] 
        except AttributeError:  main_pos_fields = [self.position, '?', '?'] 
        main_pos_fields += [self.gene, self.orientation, self.gene_feature, self.gene_distances, 
                            self.total_read_count, self.perfect_read_count, self.get_main_sequence()[0]]
        try:                                    main_pos_fields.append(self.gene_annotation[0])
        except (AttributeError, IndexError):    main_pos_fields.append('-')
        OUTPUT.write('\t'.join([str(x) for x in main_pos_fields]) + '\n')
        ### print lines for each of the RISCC genome-side reads
        # sort the RISCC reads by alignment position - happens automatically, since position is first on the list
        OUTPUT.write("RISCC_genome-side_reads::: %s\n"%data_header)
        for (position, read_data) in sorted(self.RISCC_genome_side_aligned_reads.items()):
            fields = [read_data[0].chromosome, read_data[0].strand, read_data[0].full_position.strip('?-')]
            readcount = read_data[1]
            main_seq = min(read_data[2].items(), key = lambda (s, (r,e)): (-r, e))[0]
            perfect_read_count = sum(read_count for (read_count, N_errors) in read_data[2].values() if N_errors==0)
            fields += read_data[3:7] + [readcount, perfect_read_count, main_seq] 
            try:                    fields.append(read_data[7])
            except IndexError:      fields.append('-')
            OUTPUT.write('\t'.join([str(x) for x in fields]) + '\n')
        for (seq, read_data) in sorted(self.RISCC_genome_side_unaligned_reads.items()):
            fields = [read_data[0], '-', '-'] + read_data[3:7] + [read_data[1], 0, seq] 
            try:                    fields.append(read_data[7])
            except IndexError:      fields.append('-')
            OUTPUT.write('\t'.join([str(x) for x in fields]) + '\n')


class Insertional_mutant_readcount_only(Insertional_mutant):
    """ Simplified readcount-only version of Insertional_mutant, to use as part of Insertional_mutant_multi_dataset.

    Essentially an Insertional_mutant with only readcount-related data, and no position/gene data.
    Thus some methods are overwritten with ones that just raise errors due to not being relevant to this type of object.
    """

    def __init__(self):
        # readcount-only mutants get only readcount data, no position data
        Insertional_mutant._set_readcount_related_data_to_zero(self)

    def _set_general_attributes(self, insertion_position):
        raise MutantError("It makes no sense to run _set_general_attributes on readcount-related-only mutant object"
                          +" - it has no non-readcount-related info!")

    def _copy_non_readcount_data(self, source_mutant):
        raise MutantError("It makes no sense to run _copy_non_readcount_data on readcount-related-only mutant object"
                          +" - it has no non-readcount-related info!")

    def update_gene_info(self, *args, **kwargs):
        raise MutantError("It makes no sense to run update_gene_info on readcount-related-only mutant object - it has no gene info!")


class Insertional_mutant_multi_dataset(Insertional_mutant):
    """ Multi-dataset version of Insertional_mutant: single position/gene, with read numbers/seqs from multiple datasets. 

    It has a single copy of the usual position/gene data. 
    It has a by_dataset attribute - a dataset_name:Insertional_mutant_readcount_only dictionary used to hold readcont-related
     data from multiple datasets.  This can also be accessed by the read_info method, providing the dataset name.
    Readcount-related methods take a dataset_name argument to specify which dataset they should be applied to.

    Some methods may not have multi-dataset functionality implemented, if I didn't think it would be useful.
    """

    def __init__(self, IB=None, insertion_position=SPECIAL_POSITIONS.unknown):
        # multi-dataset mutants get general attributes
        Insertional_mutant._set_general_attributes(self, IB, insertion_position)
        # instead of single readcount-related attributes, multi-dataset mutants get a by_dataset dictionary, 
        #  with dataset names as keys and readcount-related-only mutants as values.
        self.by_dataset = defaultdict(blank_readcount_only_mutant)

    def _set_readcount_related_data_to_zero(self):
        raise MutantError("_set_readcount_related_data_to_zero NOT IMPLEMENTED on multi-dataset mutant object"
                          +" - it doesn't have a simple single set of readcount-related info!")

    def _check_dataset_presence(self, dataset_name):
        if dataset_name not in self.by_dataset:
            raise MutantError("No dataset %s in this multi-dataset mutant! Present datasets are %s"%(dataset_name, 
                                                                                                 self.by_dataset.keys()))

    def _check_dataset_name_return_data(self, dataset_name, strict=False):
        """ Return readcount-related-only object for dataset_name; raise MutantError if None, or if strict=True and not present. """
        if strict:
            _check_dataset_presence(self, dataset_name)
        elif dataset_name is None:
            raise MutantError("Cannot use None as dataset name!")
        return self.by_dataset[dataset_name]

    def read_info(self, dataset_name=None, strict=False):
        """ Help function to get read-info-containing object for both multi-dataset and single mutants.

        For multi-dataset, return self.by_dataset[dataset_name] if present - if not present, raises an exception if strict, 
         otherwise returns an empty read-info object.
        """
        if dataset_name is None:  
            raise MutantError("This is a multi-dataset mutant - must provide dataset_name arg!")
        if strict:
            self._check_dataset_presence(dataset_name)
            return self.by_dataset[dataset_name]
        else:
            try:                return self.by_dataset[dataset_name]
            except KeyError:    return blank_readcount_only_mutant()
        # TODO unit-tests?

    def add_read(self, HTSeq_alignment, position=SPECIAL_POSITIONS.unknown, read_count=1, dataset_name=None):
        """ Add read to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_read(readcount_data_container, HTSeq_alignment, position, read_count)
    
    def add_counts(self, total_count, perfect_count, sequence_variant_count, assume_new_sequences=False, dataset_name=None):
        """ Add counts to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_counts(readcount_data_container, total_count, perfect_count, sequence_variant_count, 
                                      assume_new_sequences, dataset_name=None)

    def add_sequence_and_counts(self, seq, seq_count, add_to_uniqseqcount=True, dataset_name=None):
        """ Add seqs/counts to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_sequence_and_counts(readcount_data_container, seq, seq_count, add_to_uniqseqcount, dataset_name=None)
    
    @property
    def total_read_count(self):
        return sum(m.total_read_count for m in self.by_dataset.values())

    @property
    def perfect_read_count(self):
        return sum(m.perfect_read_count for m in self.by_dataset.values())

    @property
    def sequences_counts_positions_errors(self):
        joint_dict = {}
        for dataset in self.by_dataset.values():
            for seq, (count, pos, err) in dataset.sequences_counts_positions_errors.items():
                try:                joint_dict[seq][0] += count
                except KeyError:    joint_dict[seq] = [count, pos, err]
        return joint_dict

    def get_main_sequence(self, N=1, dataset_name=None, aligned_only=False):
        """ Return the most common sequence in this mutant and its count (or Nth most common sequence if N is provided).

        If dataset_name is given, return the most common sequence for just that dataset; 
         or if dataset_name is None, return most common sequence by total count over all the datasets.
        """
        if dataset_name is not None:
            seqs_to_counts_and_data = self.by_dataset[dataset_name].sequences_counts_positions_errors
        else:
            seqs_to_counts_and_data = self.sequences_counts_positions_errors
        # MAYBE-TODO print a warning if different dataset mutants have different main sequences?
        return Insertional_mutant._get_main_sequence_from_data(seqs_to_counts_and_data, N, aligned_only)
        # MAYBE-TODO should there be a warning/failure/something if it's a multi-dataset mutant and the user wants
        #  an overall main sequence and only some of the mutants have any sequence data?

    def _copy_readcount_related_data(self, source_mutant):
        raise MutantError("It makes no sense to run _copy_readcount_related_data on multi-dataset mutant object - would need to provide dataset_name, NOT IMPLEMENTED!")

    def add_other_mutant_as_dataset(self, other_mutant, other_mutant_dataset_name, 
                                    overwrite=False, check_constant_data=False):
        """ Copy all readcount-related data from other_mutant to self.by_dataset dictionary[other_mutant_dataset_name].

        If self isn't a multi-dataset mutant, raise an exception.
        If check_constant_data is True, check that the position/gene data of self and other_mutant matches.
        If self already has a other_mutant_dataset_name dataset, raise MutantError, unless overwrite=True, then overwrite.
        """
        if other_mutant_dataset_name in self.by_dataset and not overwrite:
            raise MutantError("This mutant already has a %s dataset! Can't overwrite it with "%other_mutant_dataset_name
                              +"new one.  Choose a different name for new dataset, or use overwrite=True argument.")

        # if desired, check that the position/gene data matches (and update if own gene data is unknown)
        #  (probably should be using ifs rather than asserts, but I think since they're wrapped in a try/except it's fine)
        if check_constant_data:
            if not self.position == other_mutant.position:
                raise MutantError("Can't add mutant2 as dataset to mutant1: the mutant position differs! %s and %s"%(
                                    self.position, other_mutant.position))
            try:
                self.update_gene_info(other_mutant.gene, other_mutant.orientation, 
                                      other_mutant.gene_feature, other_mutant.gene_distances)
            except MutantError:
                raise MutantError("Can't add mutant2 as dataset to mutant1: the mutant gene data differs!"
                                  +" %s, %s, %s and"%(self.gene, self.orientation, self.gene_feature)
                                  +" %s, %s, %s."%(other_mutant.gene, other_mutant.orientation, other_mutant.gene_feature))

        # make a new empty Insertional_mutant object to hold the readcount-related data from other_mutant, 
        #  and put it in the self.by_dataset dictionary under other_mutant_dataset_name
        self.by_dataset[other_mutant_dataset_name] = Insertional_mutant_readcount_only()
        # now fill this new object with readcount-related data from other_mutant
        self.by_dataset[other_mutant_dataset_name]._copy_readcount_related_data(other_mutant)

    def give_single_dataset_mutant(self, single_dataset_name, force=False):
        """ Return a single-dataset mutant based on single_dataset_name; don't modify current mutant.

        If there is no single_dataset_name in current mutant's by_dataset dictionary, raise exception, 
         unless force is True, then return new mutant with zero read-count.
        """
        if single_dataset_name not in self.by_dataset.keys() and not force:
            raise MutantError("This mutant doesn't have a %s dataset! "%single_dataset_name
                              +"Use force=True argument if you want a zero-readcount mutant returned anyway.")
        # generate new mutant, fill it with readcount-related data from self.by_dataset[single_dataset_name] 
        #  and general data from self
        new_mutant = Insertional_mutant()
        new_mutant._copy_non_readcount_data(self)
        new_mutant._copy_readcount_related_data(self.by_dataset[single_dataset_name])
        return new_mutant

    def give_all_single_dataset_mutants(self):
        """ Split multi-dataset mutant into a dataset_name:single-dataset_mutant dictionary and return it; 
        don't modify original mutant.
        """
        mutant_dictionary = defaultdict(blank_full_single_mutant)
        for dataset_name in self.by_dataset.iterkeys():
            mutant_dictionary[dataset_name] = self.give_single_dataset_mutant(dataset_name)
        return mutant_dictionary


######################################### Mutant set functions/classes ###############################################

class Dataset_summary_data():
    """ Summary data for a Insertional_mutant_pool_dataset object.  Lots of obvious attributes; no non-default methods.
    """
    # LATER-TODO update docstring

    def __init__(self, dataset, cassette_end, relative_read_direction, dataset_name=None):
        """ Initialize everything to 0/empty/unknown. """
         # make sure the arguments are valid values
        if not cassette_end in SEQ_ENDS+['?']: 
            raise ValueError("The cassette_end variable must be one of %s or '?'!"%SEQ_ENDS)
        if relative_read_direction not in RELATIVE_READ_DIRECTIONS+['?']: 
            raise ValueError("The relative_read_direction variable must be %s, or '?'!"%(', '.join(RELATIVE_READ_DIRECTIONS)))
        # reference to the containing dataset (for read-counting purposes etc), 
        #  and the dataset name (None if it's a single dataset, string for multi-datasets)
        self.dataset_name = dataset_name
        self.dataset = dataset
        # information on reads that aren't included in the dataset mutants - None or 0 by default
        # TODO I should really go over this and figure out what should be None and what should be 0 and why!!
        self.discarded_read_count, self.discarded_wrong_start, self.discarded_no_cassette = None, None, None
        self.discarded_other_end = 0
        self.non_aligned_read_count, self.unaligned, self.multiple_aligned = 0, 0, 0
        self.ignored_region_read_counts = defaultdict(int)
        # MAYBE-TODO should cassette_end and relative_read_direction be specified for the whole dataset, or just for each set of data added, in add_RISCC_alignment_files_to_data? The only real issue with this would be that then I wouldn't be able to print this information in the summary - or I'd have to keep track of what the value was for each alignment reader added and print that in the summary if it's a single value, or 'varied' if it's different values. Might also want to keep track of how many alignment readers were involved, and print THAT in the summary!  Or even print each (infile_name, cassette_end, relative_read_direction) tuple as a separate line in the header.
        self.cassette_end = cassette_end
        self.relative_read_direction = relative_read_direction

    # TODO unit-test all the methods below!

    def add_discarded_reads(self, N_all_discarded, N_wrong_start, N_no_cassette, N_other_end, replace=False):
        """ Add not-None arg values to appropriate self.discarded_* attributes (or replace them). 
        
        If the original values are None, or replace is True, replace instead of adding.
        If either the original or new value is 'unknown', the result is 'unknown' as well. 
        If any of the args is None, don't modify the original value, unless replace is True, then set to 'unknown'.
        """
        # if self doesn't have an N_other_end attribute (some older datasets don't, this is for those), set it to 0
        try:                    self.discarded_other_end
        except AttributeError:  self.discarded_other_end = 0
        # set everything
        if N_all_discarded is not None:
            if 'unknown' in (N_all_discarded, self.discarded_read_count):   self.discarded_read_count = 'unknown'
            elif replace or self.discarded_read_count is None:              self.discarded_read_count = int(N_all_discarded)
            else:                                                           self.discarded_read_count += int(N_all_discarded)
        elif replace:                                                       self.discarded_read_count = 'unknown'
        if N_wrong_start is not None:
            if 'unknown' in (N_wrong_start, self.discarded_wrong_start):    self.discarded_wrong_start = 'unknown'
            elif replace or self.discarded_wrong_start is None:             self.discarded_wrong_start = int(N_wrong_start)
            else:                                                           self.discarded_wrong_start += int(N_wrong_start)
        elif replace:                                                       self.discarded_wrong_start = 'unknown'
        if N_no_cassette is not None:
            if 'unknown' in (N_no_cassette, self.discarded_no_cassette):    self.discarded_no_cassette = 'unknown'
            elif replace or self.discarded_no_cassette is None:             self.discarded_no_cassette = int(N_no_cassette)
            else:                                                           self.discarded_no_cassette += int(N_no_cassette)
        elif replace:                                                       self.discarded_no_cassette = 'unknown'
        if N_other_end is not None:
            if 'unknown' in (N_other_end, self.discarded_other_end):        self.discarded_other_end = 'unknown'
            elif replace or self.discarded_other_end is None:               self.discarded_other_end = int(N_other_end)
            else:                                                           self.discarded_other_end += int(N_other_end)
        elif replace:                                                       self.discarded_other_end = 'unknown'
        # special case for when we don't know the specific discarded categories, but we know total discarded is 0, 
        #  so the specific categories must be 0 too:
        if self.discarded_read_count == 0:   
            self.discarded_wrong_start, self.discarded_no_cassette, self.discarded_other_end = 0, 0, 0

    def add_nonaligned_reads(self, N_all_non_aligned, N_unaligned, N_multiple_aligned, replace=False):
        """ Add not-None arg values to non_aligned_read_count, unaligned and multiple_aligned (or replace them).
        
        If the original values are None, or replace is True, replace instead of adding.
        If either the original or new value is 'unknown', the result is 'unknown' as well. 
        If any of the args is None, don't modify the original value, unless replace is True, then set to 'unknown'.
        """
        if N_all_non_aligned is not None:
            if 'unknown' in (N_all_non_aligned, self.non_aligned_read_count):   self.non_aligned_read_count = 'unknown'
            elif replace or self.non_aligned_read_count is None:                self.non_aligned_read_count = int(N_all_non_aligned)
            else:                                                               self.non_aligned_read_count += int(N_all_non_aligned)
        elif replace:                                                           self.non_aligned_read_count = 'unknown'
        if N_unaligned is not None:
            if 'unknown' in (N_unaligned, self.unaligned):          self.unaligned = 'unknown'
            elif replace or self.unaligned is None:                 self.unaligned = int(N_unaligned)
            else:                                                   self.unaligned += int(N_unaligned)
        elif replace:                                               self.unaligned = 'unknown'
        if N_multiple_aligned is not None:
            if 'unknown' in (N_multiple_aligned, self.multiple_aligned):    self.multiple_aligned = 'unknown'
            elif replace or self.multiple_aligned is None:                  self.multiple_aligned = int(N_multiple_aligned)
            else:                                                           self.multiple_aligned += int(N_multiple_aligned)
        elif replace:                                                       self.multiple_aligned = 'unknown'
        # Note: NO special case for when we don't know the specific categories, but we know total non_aligned is 0, 
        #  because for old-format files non_aligned is initially 0 but gets increased when reading the actual *.sam file, 
        #  which contains lines for unaligned reads (which are unaligned or multiple, both output the same with bowtie -m option)

    @property
    def aligned_read_count(self):
        return sum([m.read_info(self.dataset_name).total_read_count for m in self.dataset])
    @property
    def perfect_read_count(self):
        return sum([m.read_info(self.dataset_name).perfect_read_count for m in self.dataset])
    @property
    def aligned_incl_removed(self):
        return self.aligned_read_count + sum(self.ignored_region_read_counts.values())

    @property
    def processed_read_count(self):
        """ Total processed readcount (integer): aligned + unaligned (ignored if unknown) + removed due to region. """
        known_values = self.aligned_read_count + sum(self.ignored_region_read_counts.values())
        try:                return known_values + self.non_aligned_read_count
        except TypeError:   return known_values
    @property
    def processed_read_count_str(self):
        """ Total processed readcount (string): aligned + unaligned ('unknown' if unknown) + removed due to region. """
        known_values = self.aligned_read_count + sum(self.ignored_region_read_counts.values())
        try:                return str(known_values + self.non_aligned_read_count)
        except TypeError:   return str(known_values) + '+unknown'

    @property
    def full_read_count(self):
        """ Full read count (integer): processed+discarded, or just processed if discarded is unknown. """
        try:                return self.processed_read_count + self.discarded_read_count
        except TypeError:   return self.processed_read_count
    @property
    def full_read_count_str(self):
        """ Full read count as a string: processed+discarded, or processed+'unknown' if discarded is unknown. """
        if 'unknown' not in self.processed_read_count_str:
            try:                return "%s"%(self.processed_read_count + self.discarded_read_count)
            except TypeError:   return "%s+unknown"%self.processed_read_count
        else:
            try:                return "%s+unknown"%(self.processed_read_count + self.discarded_read_count)
            except TypeError:   return "%s+unknown"%self.processed_read_count

    @property
    def strand_read_counts(self):
        strand_dict = {'+': 0, '-': 0}
        for m in self.dataset:
            try:                    strand_dict[m.position.strand] += m.read_info(self.dataset_name).total_read_count
            except AttributeError:  pass
        return strand_dict

    def reads_in_chromosome(self, chromosome):
        """ Return total number of reads in given chromosome."""
        return sum(m.read_info(self.dataset_name).total_read_count 
                   for m in self.dataset if m.position not in SPECIAL_POSITIONS.all_undefined and m.position.chromosome==chromosome)

    @property
    def all_chromosomes(self):
        chromosome_set = set()
        for m in self.dataset:
            if m.read_info(self.dataset_name).total_read_count:
                try:                    chromosome_set.add(m.position.chromosome)
                except AttributeError:  pass
        return chromosome_set

    @property
    def cassette_chromosomes(self):
        return set(chrom for chrom in self.all_chromosomes if is_cassette_chromosome(chrom))
    @property
    def other_chromosomes(self):
        return set(chrom for chrom in self.all_chromosomes if is_other_chromosome(chrom))
    @property
    def non_genome_chromosomes(self):
        return self.cassette_chromosomes | self.other_chromosomes

    @property
    def N_mutants(self):
        return sum(1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count) 
    def N_mutants_over_readcount(self, readcount_min):
        return sum(1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count >= readcount_min) 
    @property
    def median_readcount(self):
        # only including mutants with non-zero reads
        return median([m.read_info(self.dataset_name).total_read_count for m in self.dataset 
                       if m.read_info(self.dataset_name).total_read_count])
    @property
    def mutants_in_genes(self):
        return len([1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count 
                    and m.gene not in SPECIAL_GENE_CODES.all_codes])
    @property
    def mutants_not_in_genes(self):
        return len([1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count 
                    and m.gene==SPECIAL_GENE_CODES.not_found])
    @property
    def mutants_undetermined(self):
        return len([1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count 
                    and m.gene in (SPECIAL_GENE_CODES.chromosome_not_in_reference, SPECIAL_GENE_CODES.not_determined)])
    @property
    def mutant_counts_by_orientation(self):
        orientation_dict = defaultdict(int)
        for m in self.dataset:
            if m.read_info(self.dataset_name).total_read_count:
                orientation_dict[m.orientation] += 1
        if '?' in orientation_dict:     del orientation_dict['?']
        if '-' in orientation_dict:     del orientation_dict['-']
        return orientation_dict
    @property
    def mutant_counts_by_feature(self):
        feature_dict = defaultdict(int)
        for m in self.dataset:
            if m.read_info(self.dataset_name).total_read_count:
                feature_dict[m.gene_feature] += 1
        if '?' in feature_dict:  del feature_dict['?']
        if '-' in feature_dict:  del feature_dict['-']
        return feature_dict

    def mutants_in_chromosome(self, chromosome):
        """ Return total number of mutants in given chromosome."""
        return sum(1 for m in self.dataset if m.read_info(self.dataset_name).total_read_count 
                   and m.position not in SPECIAL_POSITIONS.all_undefined and m.position.chromosome==chromosome)

    def merged_gene_feature_counts(self, merge_multi_splice_variants=True, merge_boundary_features=True, 
                                   merge_confusing_features=False):
        """ Return (gene_feature,count) list, biologically sorted, optionally with all "boundary" features counted as one.

        The source gene feature counts are based on the self.mutant_counts_by_feature dict.
        If merge_confusing_features==True, any locations containing '??' will be listed as '??'.
        If merge_boundary_features==True, any locations containing '|' and no '??' will be listed as 'multiple_splice_variants'.
        If merge_boundary_features==True, any locations containing '/' and no '??' will be listed as 'boundary'.
        The custom sort order (based on what seems sensible biologically) is: CDS, intron, UTR, other, boundary.
        """
        merged_feature_count_dict = defaultdict(int)
        for feature, count in self.mutant_counts_by_feature.items():
            # note that anything containing '??' AND '/' never gets merged as boundary
            if '??' in feature:
                if merge_confusing_features:                         merged_feature_count_dict['??'] += count
                else:                                                merged_feature_count_dict[feature] += count
            elif '|' in feature and merge_multiple_splice_variants:  merged_feature_count_dict['multiple_splice_variants'] += count
            elif '/' in feature and merge_boundary_features:         merged_feature_count_dict['boundary'] += count
            else:                                                    merged_feature_count_dict[feature] += count
        return merged_feature_count_dict

    @property
    def most_common_mutants(self):
        """ Return list of mutants with the most total reads (in dataset if multi-dataset)."""
        highest_readcount = max([mutant.read_info(self.dataset_name).total_read_count for mutant in self.dataset])
        highest_readcount_mutants = [mutant for mutant in self.dataset 
                                     if mutant.read_info(self.dataset_name).total_read_count==highest_readcount]
        return highest_readcount_mutants

    @property
    def all_genes_in_dataset(self):
        """ The set of all genes with at least one mutant in the dataset. """
        # the empty-set argument is needed in case there are no mutants in the dataset - set.union() with empty args is an error.
        return set.union(set(), *[set(genes) for N_mutants,genes 
                                  in self.dataset.get_gene_dict_by_mutant_number(self.dataset_name).items() if N_mutants>0])
    @property
    def N_genes_in_dataset(self):
        """ The number of genes with at least one mutant in the dataset. """
        return len(self.all_genes_in_dataset)
    @property
    def genes_with_multiple_mutants(self):
        """ The set of all genes with at TWO OR MORE mutants in the dataset. """
        # the empty-set argument is needed in case there are no mutants in the dataset - set.union() with empty args is an error.
        return set.union(set(), *[set(genes) for N_mutants,genes 
                                  in self.dataset.get_gene_dict_by_mutant_number(self.dataset_name).items() if N_mutants>1])
    @property
    def N_genes_with_multiple_mutants(self):
        """ The number of genes with TWO OR MORE mutants in the dataset. """
        return len(self.genes_with_multiple_mutants)


class Insertional_mutant_pool_dataset():
    """ A dataset of insertional mutants - contains an iterable of Insertional_mutant objects, and a lot of extra data.

    May be a multi-dataset (i.e. a single dataset that actually contains readcounts for multiple samples) - in that case
     the Insertional_mutant objects will be multi-dataset mutants.

    WILL ADD MORE INFORMATION HERE ONCE I STOP CHANGING THIS SO OFTEN.
    
    Attributes - THIS MAY BE OUT OF DATE
     - cassette_end - specifies which end of the insertion cassette the reads are on
     - relative_read_direction - specifies wiether the reads are directed inward or outward compared to the cassette
     - discarded_read_count - number of reads discarded in preprocessing before alignment (not counted in processed_read_count)
     - ignored_region_read_counts - region_name:read_count dictionary (not counted in processed_read_count) ____ REALLY??
     - processed_read_count, unaligned_read_count - various read counts, obvious
     - strand_read_counts - name:count dictionaries to keep track of reads per strand
     - mutants_in_genes, mutants_not_in_genes, mutants_undetermined - counts of mutants in genes, not in genes, unknown
     - mutant_counts_by_orientation, mutant_count_by_feature - name:count dictionaries for mutant gene location properties

    For methods see method docstrings.
    """
    # TODO update docstring to contain up-to-date sensible info on everything!
    # - make sure to mention that all mutant positions are immutable by default, and how to deal with changing them
    # - explain about normal datasets and multi-datasets
    # - ____

    # Implement new functionality for datasets:
    # - MAYBE-TODO splitting joint datasets into separate ones?  Do we need that? PROBABLY NOT.
    # - MAYBE-TODO adding and subtracting datasets (by readcount) - do we currently need that?
    # - TODO-NEXT calculating growth rates!  Right now it's in mutant_growth_rates.py separately, but that may be temporary

    # MAYBE-TODO should I make multi-datasets a subclass of normal ones instead of having the same class implement both?

    def __init__(self, cassette_end='?', relative_read_direction='?', multi_dataset=False):
        """ Initializes empty dataset; saves properties as provided. """
        # _mutants_by_IB is the main data structure here, but it's also private:
        #      see "METHODS FOR EMULATING A CONTAINER TYPE" section below for proper ways to interact with mutants.
        #   the single mutants should be single-dataset by default, or multi-dataset if the containing dataset is.
        if multi_dataset:   blank_mutant_function = blank_multi_mutant_with_IB
        else:               blank_mutant_function = blank_single_mutant_with_IB
        self._mutants_by_IB = keybased_defaultdict(blank_mutant_function)
        # various dataset summary data - single for a single dataset, a dictionary for a multi-dataset object.
        self.multi_dataset = multi_dataset
        if not multi_dataset:   self.summary = Dataset_summary_data(self, cassette_end, relative_read_direction, None)
        else:                   self.summary = {}
        # data that's NOT related to a particular dataset
        # gene/annotation-related information - LATER-TODO should this even be here, or somewhere else?
        self.gene_annotation_header = []
        self.total_genes_in_genome = 0

    ######### METHODS FOR EMULATING A CONTAINER TYPE (a dataset is essentially a container of mutants)
    # everything is currently based on the private self._mutants_by_IB dictionary, but this may change.
    # Note that all of this will work just fine for multi-datasets without changes, since multi-datasets simply contain
    #  the same dictionary but with multi-dataset mutants instead of single-dataset ones.

    def __len__(self):      
        # Return the number of mutants with some non-zero readcounts (DON'T count zero-read mutants, even if present!).
        #  A bit tricky for multi-datasets. 
        if not self.multi_dataset:  return self.summary.N_mutants
        else:                       return sum(1 for m in self if sum(x.total_read_count for x in m.by_dataset.values()))

    def __iter__(self):     return self._mutants_by_IB.itervalues()
    # for multi-datasets, for iterating over only mutants with 1+ reads in a particular dataset, see mutants_in_dataset()

    @property
    def size(self):         return len(self)

    # instead of __setitem__
    def add_mutant(self, mutant, overwrite=False):
        """ Add mutant to dataset. """
        if mutant.IB in self._mutants_by_IB.keys() and not overwrite:
            raise MutantError("Can't add mutant that would overwrite previous mutant with same IB! "
                              +"Pass overwrite=True argument if you want to overwrite.")
        self._mutants_by_IB[mutant.IB] = mutant

    # instead of __delitem__
    def remove_mutant(self, mutant_or_IB):
        """ Remove mutant (by IB) - can take a mutant or a IB. """
        try:                    IB = mutant_or_IB.IB
        except AttributeError:  IB = mutant_or_IB
        del self._mutants_by_IB[IB]

    # instead of __getitem__
    def get_mutant(self, IB):
        """ Return the mutant with given IB. If mutant doesn't exist, create a new one with no position/reads/sequences. """
        return self._mutants_by_IB[IB]

    def __contains__(self, IB):
        """ Check if dataset contains mutant with given IB.

        You can only check "IB in dataset" at present, not "mutant in dataset", since the latter would probably
         also just check by IB, so the syntax would be misleading.
        """
        return IB in self._mutants_by_IB
    
    # MAYBE-TODO implement mutant search by position?


    ######### READING BASIC DATA INTO DATASET

    @staticmethod
    def _next_until_name_match(generator, ref_name):
        """ grab the next alignment object from the generator until the name matches ref_name. """
        curr_name = ''
        while not curr_name == ref_name:
            curr_aln = generator.next()
            curr_name = curr_aln.read.name.split()[0]
        return curr_aln

    @classmethod
    def _parse_3files_parallel(cls, file1_fastx, file2_sam, file3_sam):
        """ Parse fastq+sam+sam files in parallel - generator yielding (name, seq1, aln2, aln3) tuples.

        It checks that the read names match (except for the paired-end side or anything else after a space).
        It assumes that file1 is the reference, and the other two files may have some extra names (will be ignored).
        """
        generator1 = name_seq_generator_from_fasta_fastq(file1_fastx)
        generator2 = iter(HTSeq.SAM_Reader(file2_sam))
        generator3 = iter(HTSeq.SAM_Reader(file3_sam))
        if_finished_1, if_finished_2, if_finished_3 = False, False, False
        while True:
            try:                    name1, seq1 = generator1.next()
            except StopIteration:   if_finished_1, name1 = True, 'NOTHING_HERE'
            name1 = name1.split()[0]
            try:                    aln2 = cls._next_until_name_match(generator2, name1)
            except StopIteration:   if_finished_2 = True
            try:                    aln3 = cls._next_until_name_match(generator3, name1)
            except StopIteration:   if_finished_3 = True
            # if all the files still contained data, yield it
            if not any([if_finished_1, if_finished_2, if_finished_3]):
                yield (name1, seq1, aln2, aln3)
            # if file1 was finished, we're done - it's okay if the other files had some extra reads
            elif if_finished_1:
                raise StopIteration
            # if file1 WASN'T finished but one of the others was, that's a problem!
            else:
                raise MutantError("Parsing seq/aln files in parallel - inconsistent finished states! "
                                  +"(If finished: %s %s, %s %s, %s %s)"%(file1_fastx, if_finished_1, 
                                                                         file2_sam, if_finished_2, file3_sam, if_finished_3))
        # TODO unit-tests! There are some in experiments/arrayed_library/internal_barcode_processing/code/clustering_tools.py for a similar function - test__parse_3seq_parallel

    # TODO TODO TODO need function to deal with IB-only deepseq reads!  How do we want to handle those?  As multi-dataset mutants?  Or add some new functionality?  How should it interact with the RISCC data?  Does it matter which of the datasets is read in first?
    # TODO how should the IB clustering for IB-only datasets be done?  1) Align them to the RISCC cluster centroids, 2) cluster RISCC and IB-only IB seqs all together, 3) cluster each dataset separately and look for matches. TODO

    def add_RISCC_alignment_files_to_data(self, cassette_side_flank_aligned_file, genome_side_aligned_file, IB_fastq_file, 
                                          allowed_IBs=None, IB_cluster_file=None, 
                                          best_genome_side_only=False, ignore_unaligned=False, 
                                          max_allowed_cassette_side_dist=1, max_cassette_side_ratio_to_ignore=100, 
                                          skip_checks=False, removed_mutant_file='/dev/null', quiet=False):
        """ Add paired-end RISCC reads to dataset mutants, based on IB clustering.

        Input files:
            - IB_cluster_file must be a python pickle file (.pickle) or a python code file (.py) containing a 
                centroid_seq:set_of_IB_seqs dictionary; with the centroids already adjusted as desired.  IB = Internal Barcode.
            - IB_fastq_file should be a fastq file containing the read IDs and sequences of all the clustered IBs.
            - the two aligned files should be SAM files generated from aligning the cassette-side and genome-side
                paired reads separately against the chlamy genome and insertion cassette sequence, probably using 
                deepseq_alignment_wrapper.py.  
            The fastq file and the two aligned files should have matching read IDs in the same order; the IB and 
                cassette-side files can be missing some read IDs, since some reads may not have contained the expected 
                cassette sequence, leading to their removal while trying to extract the IB and cassette-side sequence.

        If best_genome_side_only, only keep the single "best" genome-side read per mutant, rather than all of them - best meaning
          highest confirmed distance.  If ignore_unaligned, don't add unaligned/multi-aligned genome-side reads to mutants at all.
        Make sure that for each mutant, all the cassette-side positions are within max_allowed_cassette_side_dist bp of the
          main position (highest-readcount non-undefined one) - raise exception if not.

        There may be some cases where one mutant has cassette-side reads with different positions (despite having the same IB). 
        In those cases, there are three options:
         - variant positions within max_allowed_cassette_side_dist bp of the most common position are allowed to remain
         - variant positions that don't meet that criterion, but have at least max_cassette_side_ratio_to_ignore times fewer reads
            than the most common position, are removed, and the total/perfect readcounts are decreased to match
            (but the matching genome-side reads from those reads remain, since there's no easy way of removing them)
         - if there's a variant position that doesn't meet either of these criteria, the entire mutant is removed, 
            and a warning is printed to removed_mutant_file.  A total of removed mutants is printed to stdout if quiet is False.
    
        Function goes over the fastq and two aligned files in parallel. 
        """
        # TODO finish docstring
        # MAYBE-TODO add option for not including IBs at all, and making the mutant dict by cassette-side alignment position like before?
        # MAYBE-TODO at some point, maybe add full parsing of multiple alignments, to compare their positions to those of 
        #  unique-aligned cases, rather than just marking them as multiple but treating them as unaligned?
        #  Might not be worth the effort, since if we have unique-aligned cases anyway, we can use those.
        # MAYBE-TODO add ignore_cassette, cassette_only options?
        # MAYBE-TODO add collapsed_readcounts option?  That doesn't make much sense for paired-end reads.

        if self.multi_dataset:  raise MutantError("add_RISCC_alignment_files_to_data not implemented for multi-datasets!")
        if self.summary.cassette_end not in SEQ_ENDS:
            raise MutantError("Cannot add data from an alignment reader if cassette_end isn't specified! Please set the "
                  +"summary.cassette_end attribute of this Insertional_mutant_pool_dataset instance to one of %s first."%SEQ_ENDS)
        if self.summary.relative_read_direction not in RELATIVE_READ_DIRECTIONS:
            raise MutantError("Cannot add data from an alignment reader if relative_read_direction isn't set! "
                              +"Please set the relative_read_direction attribute of this Insertional_mutant_pool_dataset instance "
                              +"to one of %s first."%RELATIVE_READ_DIRECTIONS)

        # read the IB cluster file; make a read_seq:centroid_seq dictionary for fast lookup.
        if IB_cluster_file is not None:
            if type(IB_cluster_file) == dict:
                IB_centroid_to_seqs = IB_cluster_file
            elif IB_cluster_file.endswith('.pickle'):
                IB_centroid_to_seqs = unpickle(IB_cluster_file)
            else:
                raise MutantError("Unknown IB_cluster_file format in add_RISCC_alignment_files_to_data - must be .pickle filename "
                                  +"or a dictionary. Value is %s"%IB_cluster_file)
            IB_seq_to_centroid = invert_listdict_nodups(IB_centroid_to_seqs)

        # set up IB checks - return True if IB is in allowed_IBs or if no allowed_IBs was given.
        if allowed_IBs is None:     _IB_check = lambda IB: True
        else:                       _IB_check = lambda IB: IB in allowed_IBs

        for (readname, IB_seq, cassette_side_aln, genome_side_aln) in self._parse_3files_parallel(
                                            IB_fastq_file, cassette_side_flank_aligned_file, genome_side_aligned_file):
            # get the cassette insertion position (as an Insertion_position object)
            # MAYBE-TODO instead of generating cassette_side_position all the time, even with multiple identical reads, 
            #  check if seq is already present in mutant, or something?  To save time.

            # if the IB isn't in the allowed set, skip this read (unless there is no allowed set, then just keep going)
            if not _IB_check(IB_seq):   continue

            cassette_side_position = get_insertion_pos_from_flanking_region_pos(cassette_side_aln, self.summary.cassette_end, 
                                                                    self.summary.relative_read_direction, immutable_position=True)
            if ignore_unaligned and cassette_side_position in SPECIAL_POSITIONS.all_undefined:
                continue
                # TODO should probably still count it
            # grab mutant based on IB (clustered or not)
            try:                IB_centroid_seq = IB_seq_to_centroid[IB_seq]
            except NameError:   IB_centroid_seq = IB_seq
            except KeyError:    raise MutantError("IB seq %s not found in cluster dict!"%IB_seq)
            mutant = self.get_mutant(IB_centroid_seq)
            mutant.add_read(cassette_side_aln, cassette_side_position, read_count=1, dataset_name=None)
            # Parse the genome-side alignment result to figure out position; add that to the mutant
            # MAYBE-TODO make an option for the genome-side reads to be outward from the cassette? Unlikely to be needed.
            genome_side_position = get_RISCC_pos_from_read_pos(genome_side_aln, self.summary.cassette_end, 'inward')
            if genome_side_position in SPECIAL_POSITIONS.all_undefined:
                N_errors = 10
            else:
                N_errors = check_mutation_count_by_optional_NM_field(genome_side_aln, negative_if_absent=False)
            if best_genome_side_only:
                mutant.improve_best_RISCC_read(genome_side_aln.read.seq, genome_side_position, N_errors, read_count=1, 
                                               max_distance=MAX_POSITION_DISTANCE)
            else:
                if not ignore_unaligned or genome_side_position not in SPECIAL_POSITIONS.all_undefined:
                    mutant.add_RISCC_read(genome_side_aln.read.seq, genome_side_position, N_errors, read_count=1)
                # MAYBE-TODO if ignore_unaligned is True, do we still want to keep a count of unaligned seqs somehow?

        # check that all mutants have consistent cassette positions; remove ones that don't.
        if not skip_checks:
            IBs_to_remove = []
            with open(removed_mutant_file, 'w')  as REMOVED_MUTANT_FILE:
                for mutant in self:
                    if_remove = mutant.decide_and_check_position(max_allowed_cassette_side_dist, 
                                                     ratio_to_ignore=max_cassette_side_ratio_to_ignore, OUTPUT=REMOVED_MUTANT_FILE)
                    if if_remove:   IBs_to_remove.append(mutant.IB)
                summary_text = ("Removed %s/%s mutants due to different flanking seq positions in one mutant "
                                +"(if distance >%s and some are within %sx reads of each other).")%(
                                  len(IBs_to_remove), len(self), max_allowed_cassette_side_dist, max_cassette_side_ratio_to_ignore)
                REMOVED_MUTANT_FILE.write("SUMMARY: " + summary_text + '\n')
            if not quiet: print summary_text + " - see %s for details."%removed_mutant_file
            for IB in IBs_to_remove:
                self.remove_mutant(IB)

        # TODO do we want to add different read category counts to the summary, or make that stuff properties?
        # MAYBE-TODO it might be good to just generate two separate mutant-sets, normal and cassette, with an option called separate_cassette or something, and print them to separate files - but that's more complicated, and right now I don't have the setup for a single dataset having multiple mutant-sets (although I guess I will have to eventually, for removed mutants etc). Right now I do it in mutant_count_alignments.py, which works but there's a lot of code repetition...


    ######### SUMMARY INFORMATION

    @property
    def mutants_by_gene(self):
        """ Return gene_name:mutant_list dict based on full list of dataset mutants; ignore SPECIAL_GENE_CODES) """
        mutants_by_gene = defaultdict(list)
        for mutant in self:
            if mutant.gene not in SPECIAL_GENE_CODES.all_codes:
                mutants_by_gene[mutant.gene].append(mutant)
        return mutants_by_gene
        # LATER-TODO add unit-test

    def get_gene_dict_by_mutant_number(self, dataset=None):
        """ Return mutant_count:gene_ID_set dict (genes with 1/2/etc mutants - in a particular dataset, if multi-dataset). 
        
        If the object is multi-dataset, dataset name must be provided, otherwise it cannot.
        """
        gene_by_mutantN = defaultdict(set)
        for (gene,mutants) in self.mutants_by_gene.iteritems():
            gene_by_mutantN[len([m for m in mutants if m.read_info(dataset).total_read_count])].add(gene)
        # check that the numbers add up to the total number of genes
        all_genes = set([mutant.gene for mutant in self]) - set(SPECIAL_GENE_CODES.all_codes)
        assert sum([len(geneset) for geneset in gene_by_mutantN.itervalues()]) == len(all_genes)
        return gene_by_mutantN
        #LATER-TODO add unit-test


    ######### MULTI-DATASET METHODS

    def _set_joint_genome_info(self, gene_annotation_header_values, total_genes_in_genome_values):
        """ Set self.gene_annotation_header and self.total_genes_in_genome based on inputs. 

        Set to blank if all inputs are blank; otherwise to the single unique non-blank value on the list; 
          if there are multiple distinct non-blank values, raise exception.  
        """
        # Merge any pieces of global information that's not per-dataset
        self.gene_annotation_header = merge_values_to_unique(gene_annotation_header_values, blank_value=[], convert_for_set=tuple, 
                                                value_name='gene_annotation_header', context='datasets in multi-dataset')
        self.total_genes_in_genome = merge_values_to_unique(total_genes_in_genome_values, blank_value=0, 
                                                value_name='total_genes_in_genome', context='datasets in multi-dataset')

    def populate_multi_dataset(self, source_dataset_dict, overwrite=False, check_gene_data=True):
        """ Given a dataset_name:single_dataset_object dictionary, populate current multi-dataset with the data. 

        If dataset already has a dataset with the same name, raise MutantError unless overwrite=True is passed.
        When merging mutants, don't double-check that their positions/genes match, unless check_gene_data=True is passed.
        Separate copies of mutant data are made, but the summary data object is added by reference.
        """
        if not self.multi_dataset:  raise MutantError("populate_multi_dataset can only be run for multi-datasets!")

        for dataset_name,dataset_object in source_dataset_dict.iteritems():

            if any([s in dataset_name for s in ' \t\n']):
                raise MutantError("Dataset name '%s' contains spaces/tabs/newlines - not allowed!"%dataset_name)
            if dataset_name in self.summary and not overwrite:
                raise MutantError("Dataset %s is already present in this multi-dataset! "%dataset_name
                                  +"Pass overwrite=True argument if you want to overwrite the previous data.")

            # copy source dataset summaries to own summary dict:
            #  (just doing a shallow-copy to avoid copying the whole dataset... copying ignored_region_read_counts
            summary_copy = copy.copy(dataset_object.summary)
            summary_copy.ignored_region_read_counts = defaultdict(int, dataset_object.summary.ignored_region_read_counts)
            summary_copy.dataset = self
            summary_copy.dataset_name = dataset_name
            self.summary[dataset_name] = summary_copy

            # Join source dataset mutant data into own multi-dataset mutants
            #  (get_mutant will create new empty mutant if one doesn't exist)
            for mutant in dataset_object:
                # TODO should a both-strand mutant and a + or -strand mutant at the same position be considered the same mutant during joining datasets??  Currently they're not, and this can give pretty odd results!  Make it an option?  Make it depend on readcounts?  I'm not really sure...  (Also see "TODO how should mutant lookup by position deal with both-strand mutants?")
                curr_mutant = self.get_mutant(mutant.IB)
                curr_mutant.add_other_mutant_as_dataset(mutant, dataset_name, overwrite=overwrite, 
                                                        check_constant_data=check_gene_data)
            # MAYBE-TODO add option to only include mutants with non-zero reads (total or perfect) in all datasets?  Or should that only happen during printing?  Or do we even care?  If I ever want to do that, there was code for it in the old version of mutant_join_datasets.py (before 2012-04-26)

        # Join any pieces of global information that's not per-dataset
        #  using getattr instead of just d.total_genes_in_genome because some older datasets don't HAVE the total_genes_in_genome
        #   attribute, and getattr lets me give a default of 0 when the attribute is missing 
        #   (and 0 is used as blank_value in the merge_values_to_unique call in self._set_joint_genome_info).
        self._set_joint_genome_info([d.gene_annotation_header for d in source_dataset_dict.values()], 
                                    [getattr(d,'total_genes_in_genome',0) for d in source_dataset_dict.values()])
        # This has no unit-tests, but has a run-test in mutant_join_datasets.py

    def mutants_in_dataset(self, dataset_name=None):
        """ List of all mutants with non-zero reads in dataset_name (or all mutants if dataset_name=None)."""
        return [mutant for mutant in self if dataset_name is None or mutant.by_dataset[dataset_name].total_read_count>0]

    def _check_dataset_consistency(self):
        """ Raise MutantError if self.summary, mutants, and sef.dataset_order don't all have the same set of datasets! """
        if not self.multi_dataset:  
            raise MutantError("_check_dataset_consistency only makes sense for multi-datasets!")
        def _check_sets_raise_error(set1, set2, set1_name, set2_name):
            if not set1==set2:
                raise MutantError("Multi-dataset mutant pool has different %s and %s dataset sets! %s, %s"%(set1_name, 
                                                                                              set2_name, set1, set2))
        datasets_from_summary = set(self.summary.keys())
        datasets_from_mutants = set.union(*[set(m.by_dataset.keys()) for m in self])
        _check_sets_raise_error(datasets_from_summary, datasets_from_mutants, "from summary", "from mutants")
        try:
            if self._dataset_order is not None:     
                datasets_from_order = set(self._dataset_order)
                _check_sets_raise_error(datasets_from_order, datasets_from_summary, "from dataset_order", "from summary")
        except AttributeError:
            pass

    @property
    def dataset_order(self):
        """ A specific order of datasets, for printing - can be set directly, defaults to alphabetical sort. """
        self._check_dataset_consistency()
        try:                    return self._dataset_order
        except AttributeError:  return sorted(self.summary.keys())

    @dataset_order.setter
    def dataset_order(self, order):
        self._dataset_order = order
        self._check_dataset_consistency()


    ######### PROCESSING/MODIFYING DATASET

    # MAYBE-TODO add function to join two datasets together, merge mutant seqs/counts etc?

    def remove_mutants_in_other_dataset(self, other_dataset, readcount_min=1, perfect_reads=False):
        """ Remove any mutants with at least readcount_min reads in other_dataset (or perfect reads, if perfect_reads=True)

        This is based on EXACT position equality: a ?-101 mutant won't be removed due to a 100-? or 100-101 or 100-102 one.
        """
        # TODO do I want this to be based on non-exact position equality instead?
        if perfect_reads:   get_readcount = lambda m: m.perfect_read_count
        else:               get_readcount = lambda m: m.total_read_count
        # go over all mutants in self; need to convert the iterator to a list to make a separate copy, 
        #  otherwise we'd be modifying the iterator while iterating through it, which isn't allowed.
        for mutant in list(self):
            if get_readcount(other_dataset.get_mutant(mutant.IB)) >= readcount_min:
                self.remove_mutant(mutant.IB)
        # TODO really I shouldn't be removing mutants outright, just noting them as removed or something...  In that case should they or should they not show up in "for m in self"?  Probably not - they should have a separate dictionary?
        # TODO should I keep track of removed reads, and print in summary?  PROBABLY.
        # LATER-TODO unit-test - it does have run-tests though.

    def remove_mutants_not_in_other_dataset(self, other_dataset, readcount_min=1, perfect_reads=False):
        """ Remove any mutants with at least readcount_min reads in other_dataset (or perfect reads, if perfect_reads=True)

        This is based on EXACT position equality: a ?-101 mutant won't be removed due to a 100-? or 100-101 or 100-102 one.
        """
        # TODO do I want this to be based on non-exact position equality instead?
        if perfect_reads:   get_readcount = lambda m: m.perfect_read_count
        else:               get_readcount = lambda m: m.total_read_count
        # go over all mutants in self; need to convert the iterator to a list to make a separate copy, 
        #  otherwise we'd be modifying the iterator while iterating through it, which isn't allowed.
        for mutant in list(self):
            if get_readcount(other_dataset.get_mutant(mutant.IB)) < readcount_min:
                self.remove_mutant(mutant.IB)
        # TODO really I shouldn't be removing mutants outright, just noting them as removed or something...  In that case should they or should they not show up in "for m in self"?  Probably not - they should have a separate dictionary?
        # TODO should I keep track of removed reads, and print in summary?  PROBABLY.
        # LATER-TODO unit-test - it does have run-tests though.

    def remove_mutants_below_readcount(self, min_readcount, perfect_reads=False):
        """ Remove any mutants with below readcount_min reads (or perfect reads, if perfect_reads=True)
        """
        if perfect_reads:   get_readcount = lambda m: m.perfect_read_count
        else:               get_readcount = lambda m: m.total_read_count
        # go over all mutants in self; need to convert dataset to a list to make a separate copy, 
        #  otherwise we'd be modifying the dataset while iterating through it, which isn't allowed.
        for mutant in list(self):
            if get_readcount(mutant) < min_readcount:
                self.remove_mutant(mutant.IB)
        # TODO really I shouldn't be removing mutants outright, just noting them as removed or something...  In that case should they or should they not show up in "for m in self"?  Probably not - they should have a separate dictionary?
        # TODO should I keep track of removed reads, and print in summary?  MAYBE.

    def find_genes_for_mutants(self, genome_version, genefile, detailed_features=True, include_RISCC_reads=False, 
                               nearest_genes_for_intergenic=False, N_run_groups=3, verbosity_level=1):
        """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

        ALSO add gene data to all the RISCC-genome-side-read mutants inside each mutant!

        If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR.
        Read the file in N_run_groups passes to avoid using up too much memory/CPU.
        """ 
        if self.multi_dataset:  raise MutantError("find_genes_for_mutants not implemented for multi-datasets!")
        # MAYBE-TODO implement for multi-datasets?  The actual gene-finding would be easy, since it'd just work on 
        #  multi-dataset mutants instead of single-dataset ones; adding stuff to summary would be harder.

        # Group all the mutants by chromosome, so that I can go over each chromosome in genefile separately
        #   instead of reading in all the data at once (which uses a lot of memory)
        #  Inclue both the main mutants, AND all the RISCC genome-side read sub-mutants if wanted.
        insertion_data_by_chromosome = defaultdict(list)
        for mutant in self:
            if mutant.position not in SPECIAL_POSITIONS.all_undefined:
                insertion_data_by_chromosome[mutant.position.chromosome].append(mutant)
            if include_RISCC_reads:
                for RISCC_read_data in mutant.RISCC_genome_side_aligned_reads.values():
                    insertion_data_by_chromosome[RISCC_read_data[0].chromosome].append(RISCC_read_data)
        self._find_genes_for_list(insertion_data_by_chromosome, genome_version, genefile, 
                                  detailed_features, nearest_genes_for_intergenic, N_run_groups, verbosity_level)

    def _find_genes_for_list(self, insertion_data_by_chromosome, genome_version, genefile, detailed_features=False, 
                             nearest_genes_for_intergenic=False, N_run_groups=3, verbosity_level=1):
        # First get the list of all chromosomes in the file, WITHOUT reading it all into memory
        with open(genefile) as GENEFILE:
            GFF_limit_data = GFF.GFFExaminer().available_limits(GENEFILE)
            chromosomes_and_counts = dict([(c,n) for ((c,),n) in GFF_limit_data['gff_id'].iteritems()])
            all_reference_chromosomes = set(chromosomes_and_counts.keys())

        # Now lump the chromosomes into N_run_groups sets with the feature counts balanced between sets, 
        #  to avoid using too much memory (by reading the whole file at once), 
        #   or using too much time (by reading the whole file for each chromosome/scaffold)
        chromosome_sets = split_into_N_sets_by_counts(chromosomes_and_counts, N_run_groups)

        ### go over all mutants/insertions on each chromosome, figure out which gene they're in (if any), keep track of totals
        # keep track of all the mutant and reference chromosomes to catch chromosomes that are absent in reference
        summ = self.summary
        self.total_genes_in_genome = 0
        for chromosome_set in chromosome_sets:
            genefile_parsing_limits = {'gff_id': list(chromosome_set)}
            if not detailed_features: 
                genefile_parsing_limits['gff_type'] = ['gene']
            with open(genefile) as GENEFILE:
                for chromosome_record in GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
                    if verbosity_level>1:   print "    parsing %s for mutant gene locations..."%chromosome_record.id
                    self.total_genes_in_genome += len(chromosome_record.features)
                    for thing in insertion_data_by_chromosome[chromosome_record.id]:
                        if isinstance(thing, Insertional_mutant):  position = thing.position
                        else:                                      position = thing[0]
                        gene_data = find_gene_by_pos_gff3(position, chromosome_record, detailed_features, 
                                                          nearest_genes_for_intergenic, quiet=(verbosity_level==0))
                        # for genome v5.5, have to strip .v5.5 suffix from gene IDs from gff file
                        if genome_version == 5.5:
                            gene_data[0] = gene_data[0].replace('.v5.5', '')
                        if isinstance(thing, Insertional_mutant):
                            thing.gene, thing.orientation, thing.gene_feature, thing.gene_distances = gene_data
                        else:
                            thing[3:] = gene_data
                        # TODO gene_data now includes distances from gene ends as the fourth thing - use that!
                    if verbosity_level>1:   print "    ...found total %s genes."%(len(chromosome_record.features))
        if verbosity_level>1:   print "    found total %s genes in full genome."%(self.total_genes_in_genome)

        # for mutants or positions in chromosomes that weren't listed in the genefile, use special values
        for chromosome in set(insertion_data_by_chromosome.keys())-set(all_reference_chromosomes):
            if not is_cassette_chromosome(chromosome):
                print 'Warning: chromosome "%s" not found in genefile data!'%(chromosome)
            for thing in insertion_data_by_chromosome[chromosome]:
                gene_data = (SPECIAL_GENE_CODES.chromosome_not_in_reference,'-','-','-')
                if isinstance(thing, Insertional_mutant):
                    thing.gene, thing.orientation, thing.gene_feature, thing.gene_distances = gene_data
                else:
                    thing[3:] = gene_data

    @staticmethod
    def _get_annotation_for_gene(gene, gene_annotation_dict):
        """ Add gene annotation data to mutant.gene_annotation, including multiple-gene cases; return True if annotations found.

        If mutant.gene is a single, this just gets the annotation list for it from gene_annotation_dict 
         and puts that in mutant.gene_annotation (or [] if the gene isn't in gene_annotation_dict);
        If mutant.gene is a MULTIPLE_GENE_JOIN-separated string with multiple gene IDs, this finds the annotation for all of them, 
         zips them together and joins similarly with MULTIPLE_GENE_JOIN - for instance if gene was "A | B" and the annotations
          for A were [a, 1, 100] and for B were [b, 2, 100], the resulting annotations would be ["a | b", "1 | 2", "100 | 100"].
          No duplicate-removal is done; genes with no annotation are skipped.
        """
        # grab annotations for each gene
        annotations = []
        for gene in gene.split(MULTIPLE_GENE_JOIN):
            try:                annotations.append(gene_annotation_dict[gene])
            except KeyError:    pass
        # make joint annotation (each field for all genes); 
        #  make this look better by dealing with empty data specially - turn " & " into "" and " & x" into "- & x", 
        joint_annotations = []
        for ann in zip(*annotations):
            if any(ann):
                ann = [a if a else '-' for a in ann]
                joint_annotations.append(MULTIPLE_GENE_JOIN.join(ann))
            else:
                joint_annotations.append('')
        # MAYBE-TODO do duplicate-removal etc?  But that might just make things confusing - not obvious what goes with which gene.
        return joint_annotations
        # TODO unit-test!

    def add_gene_annotation(self, genome_version, include_RISCC_reads=False, print_info=False):
        """ Add gene annotation to each mutant, based on multiple annotation_files for that genome version.
        """
        # add the annotation info to each mutant (or nothing, if gene has no annotation)
        # MAYBE-TODO should I even store gene annotation in each mutant (AND in each genome-side LEAPseq read), or just keep a separate per-gene dictionary to save space?
        gene_annotation_dict, gene_annotation_header = get_all_gene_annotation(genome_version, print_info=False)
        if gene_annotation_header:  self.gene_annotation_header = gene_annotation_header
        else:                       self.gene_annotation_header = 'GENE_ANNOTATION_DATA'
        # add the annotation info to each mutant (or nothing, if gene has no annotation) 
        N_annotated = 0
        for mutant in self:
            annotation = self._get_annotation_for_gene(mutant.gene, gene_annotation_dict)
            mutant.gene_annotation = annotation
            if annotation:          N_annotated += 1
            if include_RISCC_reads:
                for RISCC_data in mutant.RISCC_genome_side_aligned_reads.values():
                    annotation = self._get_annotation_for_gene(RISCC_data[3], gene_annotation_dict)
                    RISCC_data[7:] = annotation
                    if annotation:      N_annotated += 1
        if print_info:          print "Added %s annotations"%N_annotated
        elif not N_annotated:   print "Warning: No gene annotations found!"
        # LATER-TODO add this to the gene-info run-test case!  But the get_all_gene_annotation method has tests.

    # MAYBE-TODO implement mutant-merging or counting adjacent mutants?   See code and unit/run-tests in mutant_analysis_classes.py.

    ######### WRITING DATA TO FILES

    def _get_summary(self, dataset=None):
        """ Help function to unify single/multi-datasets: return summary if dataset=None, else summary[dataset]. """
        # TODO should this be in some other section?
        if dataset is None:     return self.summary
        else:                   return self.summary[dataset]

    def _most_common_mutants_info(self, dataset=None):
        """ Return a string containing details for most common mutant, or count of most common mutants if multiple. """
        summ = self._get_summary(dataset)
        most_common_mutants = summ.most_common_mutants
        m = most_common_mutants[0]
        # calculate the fraction of total reads per mutant, assuming each mutant has the same readcount
        assert len(set([m.read_info(dataset).total_read_count for m in most_common_mutants])) == 1
        readcount_info = value_and_percentages(m.read_info(dataset).total_read_count, [summ.aligned_read_count])
        if len(most_common_mutants) == 1:   return "%s (%s)"%(readcount_info, m.position)
        else:                               return "%s (%s mutants)"%(readcount_info, len(most_common_mutants))

    @staticmethod
    def _sort_feature_list(feature_list):
        """ Sort feature list: biologically (CDS, intron, UTR, other, boundary), then alphabetically in each category. """
        # using a dictionary to assign each category a sorting position, with "other" (default) between UTR and boundary.
        feature_count_sorted_list = sorted(feature_list, key=lambda f: (GENE_FEATURE_ORDER[f],f))
        return feature_count_sorted_list

    @staticmethod
    def _make_genelist_str(geneset, N_genes_to_print):
        """ Given a set of names, return a string containing at most N_genes of the names: '(A,B,...)', '(A,B)' or ''. """
        if N_genes_to_print>0:
            genelist_to_print = ', '.join(sorted(geneset)[:N_genes_to_print])
            if len(geneset)<=N_genes_to_print:  return '(%s)'%genelist_to_print
            else:                               return '(%s, ...)'%genelist_to_print
        else:                                   return ''

    def print_summary(self, OUTPUT=sys.stdout, N_genes_to_print=5, line_prefix='    ', header_prefix=' * ', 
                      merge_multiple_splice_variants=True, merge_boundary_features=True, count_cassette=True, count_other=True):
        """ Print basic summary info about the dataset/s: read and mutant counts and categories, gene numbers, etc.

        Prints tab-separated table for multi-datasets-1.
        Prints to stdout by default, can also pass an open file object).
        """
        # TODO need proper extra summary stuff for RISCC data!
        #   #correct/incorrect positions, #cassette fragments, for correct positions distribution of distance checked, ...

        ### define a list of datasets+summaries that we'll be dealing with
        if not self.multi_dataset:  summaries_and_datasets = [(self.summary, None)]
        else:                       summaries_and_datasets = [(self.summary[dataset],dataset) for dataset in self.dataset_order]
        summaries, datasets = zip(*summaries_and_datasets)

        ### First set up a list of line descriptions and value-getter functions
        # Note on lambdas: wherever there's a for loop, I have to use lambda x=x (region,orientation,etc) 
        #  to bind x to the CURRENT value of x - otherwise whenever the lambda got called, 
        #   the current environment value of x would be used instead.
        descriptions_and_value_getters = DVG = []

        DVG.append((header_prefix+"Total reads in sample:", 
                    lambda summ: "%s"%(summ.full_read_count_str) ))
        def _fraction_or_unknown(value, totals):
            try:                return value_and_percentages(value, totals)
            except TypeError:   return "%s (unknown)"%value
        DVG.append((header_prefix+"Reads discarded in preprocessing (% of total):", 
                   lambda summ: _fraction_or_unknown(summ.discarded_read_count, [summ.full_read_count])))
        DVG.append((line_prefix+"discarded due to wrong start (% of total):", 
                   lambda summ: _fraction_or_unknown(summ.discarded_wrong_start, [summ.full_read_count])))
        DVG.append((line_prefix+"discarded due to no cassette (% of total):", 
                   lambda summ: _fraction_or_unknown(summ.discarded_no_cassette, [summ.full_read_count])))
        # using getattr because some older datasets don't HAVE the discarded_other_end attribute, and getattr lets me give 
        #  a default of 0 when the attribute is missing (which is what it should be - old datasets didn't have that functionality) 
        all_other_end = [getattr(summ,'discarded_other_end',0) for summ in summaries]
        all_other_end = sum([0 if x=='unknown' else x for x in all_other_end])
        if all_other_end:
            DVG.append((line_prefix+"separated other-end reads (3'/5') (% of total):", 
                       lambda summ: _fraction_or_unknown(getattr(summ,'discarded_other_end',0), [summ.full_read_count])))

        DVG.append((header_prefix+"Reads without a unique alignment (% of total, % of post-preprocessing):", 
                    lambda summ: _fraction_or_unknown(summ.non_aligned_read_count, 
                                                              [summ.full_read_count, summ.processed_read_count]) ))
        DVG.append((line_prefix+"unaligned reads (% of total, % of post-preprocessing):", 
                    lambda summ: _fraction_or_unknown(summ.unaligned, 
                                                              [summ.full_read_count, summ.processed_read_count]) ))
        DVG.append((line_prefix+"multiply aligned reads (% of total, % of post-preprocessing):", 
                    lambda summ: _fraction_or_unknown(summ.multiple_aligned, 
                                                              [summ.full_read_count, summ.processed_read_count]) ))

        DVG.append((header_prefix+"Uniquely aligned reads (% of total, % of post-preprocessing):",
                    lambda summ: value_and_percentages(summ.aligned_incl_removed, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))

        all_ignored_regions = set.union(*[set(summ.ignored_region_read_counts) for summ in summaries])
        for region in sorted(all_ignored_regions):
            DVG.append((line_prefix+"Removed reads aligned to %s (%% of total, %% of post-preprocessing):"%region, 
                        lambda summ,region=region: value_and_percentages( summ.ignored_region_read_counts[region], 
                                                                    [summ.full_read_count, summ.processed_read_count]) ))
        if all_ignored_regions:
            DVG.append((header_prefix+"Remaining aligned reads (% of total, % of post-preprocessing):", 
                        lambda summ: value_and_percentages(summ.aligned_read_count, 
                                                                   [summ.full_read_count, summ.processed_read_count]) ))

        DVG.append((line_prefix+"Perfectly aligned reads, no mismatches (% of aligned):", 
                    lambda summ: value_and_percentages(summ.perfect_read_count, [summ.aligned_read_count]) ))

        for strand in sorted(set.union(*[set(summ.strand_read_counts) for summ in summaries])):
            DVG.append((line_prefix+"Reads with cassette direction matching chromosome %s strand (%% of aligned):"%strand,
                        lambda summ,strand=strand: value_and_percentages(summ.strand_read_counts[strand], 
                                                                                 [summ.aligned_read_count]) ))

        special_chromosomes = []
        if count_cassette:  special_chromosomes += sorted(set.union(*[set(summ.cassette_chromosomes) for summ in summaries]))
        if count_other:     special_chromosomes += sorted(set.union(*[set(summ.other_chromosomes) for summ in summaries]))

        for chromosome in special_chromosomes:
            DVG.append((line_prefix+"Reads aligned to %s (%% of aligned):"%chromosome, 
                        lambda summ,chromosome=chromosome: value_and_percentages(summ.reads_in_chromosome(chromosome), 
                                                                                 [summ.aligned_read_count]) ))
        
        DVG.append((header_prefix+"Distinct mutants (read groups) by cassette insertion position:", 
                    lambda summ: "%s"%summ.N_mutants ))
        DVG.append((line_prefix+"(mutants with 2+, 10+, 100+, 1000+ reads):",
                    lambda summ: "(%s, %s, %s, %s)"%tuple([summ.N_mutants_over_readcount(X) for X in (2,10,100,1000)]) ))
        DVG.append((line_prefix+"(read location with respect to cassette: which end, which direction):", 
                    lambda summ: "(%s, %s)"%(summ.cassette_end, summ.relative_read_direction) ))
        DVG.append((line_prefix+"(average and median reads per mutant):", 
                    lambda summ: "(%d, %d)"%(round((summ.aligned_read_count)/summ.N_mutants), round(summ.median_readcount))))

        DVG.append((line_prefix+"Most common mutant(s): reads (% of aligned) (position or count):",
                    lambda summ: self._most_common_mutants_info(summ.dataset_name) ))
        # MAYBE-TODO may also be a good idea to keep track of the most common SEQUENCE, not just mutant...

        for chromosome in special_chromosomes:
            DVG.append((line_prefix+"Mutant cassettes in %s (%% of total):"%chromosome, 
                        lambda summ,chromosome=chromosome: value_and_percentages(summ.mutants_in_chromosome(chromosome), 
                                                                                 [summ.N_mutants]) ))

        # print the gene annotation info, but only if there is any
        if any([summ.mutants_in_genes+summ.mutants_not_in_genes for summ in summaries]):
            DVG.append((line_prefix+"Mutant cassettes on chromosomes with no gene data"
                                   +" (cassette, some scaffolds, maybe chloroplast/mito) (% of total):", 
                        lambda summ: value_and_percentages(summ.mutants_undetermined, [summ.N_mutants]) ))
            DVG.append((line_prefix+"Mutant cassettes in intergenic spaces (% of total, % of known):", 
                        lambda summ: value_and_percentages(summ.mutants_not_in_genes, 
                                                   [summ.N_mutants, summ.mutants_not_in_genes+summ.mutants_in_genes]) ))
            DVG.append((header_prefix+"Mutant cassettes inside genes (% of total, % of known):", 
                        lambda summ: value_and_percentages(summ.mutants_in_genes, 
                                                   [summ.N_mutants, summ.mutants_not_in_genes+summ.mutants_in_genes]) ))
            for orientation in sorted(set.union(*[set(summ.mutant_counts_by_orientation) for summ in summaries]), reverse=True):
                DVG.append((line_prefix+"Mutant cassettes in %s orientation to gene (%% of ones in genes):"%orientation, 
                            lambda summ,o=orientation: value_and_percentages(summ.mutant_counts_by_orientation[o], 
                                                                                     [summ.mutants_in_genes]) ))
            # custom order for features to make it easier to read: CDS, intron, UTRs, everything else alphabetically after
            # MAYBE-TODO also give print_summary an option for merge_confusing_features arg to merged_gene_feature_counts?
            for feature in self._sort_feature_list(set.union(
                                *[set(summ.merged_gene_feature_counts(merge_boundary_features, merge_multiple_splice_variants)) 
                                  for summ in summaries])):
                DVG.append((line_prefix+"Mutant cassettes in gene feature %s (%% of ones in genes):"%feature, 
                            lambda summ,feature=feature: value_and_percentages(
                                summ.merged_gene_feature_counts(merge_boundary_features, merge_multiple_splice_variants)[feature], 
                                [summ.mutants_in_genes]) ))

            DVG.append((header_prefix+"Genes containing a mutant (% of all genes):", 
                        lambda summ: value_and_percentages(summ.N_genes_in_dataset, 
                                                                   [self.total_genes_in_genome]) ))
            DVG.append((line_prefix+"Genes containing at least two mutants (% of all genes):", 
                        lambda summ: value_and_percentages(summ.N_genes_with_multiple_mutants,
                                                                   [self.total_genes_in_genome]) ))
            DVG.append((line_prefix+"  (total genes in genome annotation data):", 
                        lambda summ: "(%s)"%self.total_genes_in_genome ))
            # MAYBE-TODO put some kind of maximum on this or collapse into ranges rather than listing all the numbers?
            for mutantN in sorted(set.union(*[set(self.get_gene_dict_by_mutant_number(dataset)) for dataset in datasets])):
                DVG.append((line_prefix+"Genes with %s mutants (%% of all genes):"%mutantN, 
                            lambda summ,N=mutantN: value_and_percentages(
                                                            len(self.get_gene_dict_by_mutant_number(summ.dataset_name)[N]),
                                                            [self.total_genes_in_genome]) ))
                DVG.append((line_prefix+"  (some gene names):",
                            lambda summ,N=mutantN: self._make_genelist_str(
                                            self.get_gene_dict_by_mutant_number(summ.dataset_name)[N], N_genes_to_print) ))
            # LATER-TODO Add some measure of mutations, like how many mutants have <50% perfect reads (or something - the number should probably be a command-line option).  Maybe how many mutants have <20%, 20-80%, and >80% perfect reads (or 10 and 90, or make that a variable...)
        
        ### print the data: line description and tab-separated list of values for each dataset.
        # for multi-dataset mutants, write header line with dataset names in order
        if self.multi_dataset:      
            OUTPUT.write(header_prefix + 'DATASETS\t' + '\t'.join(datasets) + '\n')
        for line_description, line_value_getter in descriptions_and_value_getters:
            OUTPUT.write(line_description)
            for summ,dataset in summaries_and_datasets:
                OUTPUT.write('\t' + line_value_getter(summ))
            OUTPUT.write('\n')

    def _sort_data(self, sort_data_by='position'):
        """ Sort the mutants by position or readcount, or leave unsorted. """
        all_mutants = iter(self)
        if sort_data_by=='position':
            sorted_data = sorted(all_mutants, key = lambda m: (m.position, m.IB))
            # x.position here is an Insertion_position object and has a sensible cmp function
            # TODO do unaligned/multi-aligned/unknown positions sort sensibly here?
        elif sort_data_by=='read_count':
            if self.multi_dataset:  
                raise MutantError("Sorting by readcount in print_data not implemented for multi-datasets!")
            sorted_data = sorted(all_mutants, key = lambda m: (m.total_read_count, m.perfect_read_count, m.position, m.IB), 
                                 reverse=True)
        else:
            raise MutantError("Can't sort mutants by %s - only position or readcount are implemented!"%sort_data_by)
        return sorted_data

    def print_data(self, OUTPUT=sys.stdout, sort_data_by=None, header_line=True, header_prefix="# "):
        """ Print full data, one line per mutant: position data, gene info, read counts, RISCC status/info.
        (see the file header line for exactly what all the output fields are).

        Data is printed to OUTPUT, which should be an open file object (stdout by default).
         Output is tab-separated, with optional header starting with "# ".  
        """
        ### print the header line (different for normal and multi-datasets)
        # TODO should separate the number of RISCC reads and IB-only reads, really...
        if header_line:
            header = ['chromosome', 'strand', 'min_position', 'full_position', 'total_reads', 'perfect_reads',
                      'confirmed-dist', 'N-genome-side-conf-reads', 'N-genome-side-nonconf-reads',
                      'gene', 'orientation', 'feature', 'gene-end-distances', 'IB_seq', 'main_cassette_side_seq']
            if self.multi_dataset:
                for dataset_name in self.dataset_order:
                    header += ['reads_in_%s'%dataset_name, 'perfect_in_%s'%dataset_name]
            header += self.gene_annotation_header
            OUTPUT.write(header_prefix + '\t'.join(header) + "\n")

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        sorted_mutants = self._sort_data(sort_data_by)

        # create "empty" annotation line with the correct number of fields, for genes that weren't in annotation file
        if self.gene_annotation_header:
            missing_gene_annotation_data = ['NO GENE DATA'] + ['' for x in self.gene_annotation_header[:-1]]

        ### for each mutant, print the mutant data line (different for normal and multi-datasets)
        for mutant in sorted_mutants:
            try:
                mutant_data = [mutant.position.chromosome, mutant.position.strand, 
                               mutant.position.min_position, mutant.position.full_position]
            except AttributeError:
                mutant_data = [mutant.position, '?', '?', '?']
            mutant_data += [mutant.total_read_count, mutant.perfect_read_count]
            mutant_data += [mutant.RISCC_max_confirmed_distance(), 
                            mutant.RISCC_N_confirming_seqs(), mutant.RISCC_N_non_confirming_seqs()]
            mutant_data += [mutant.gene, mutant.orientation, mutant.gene_feature, mutant.gene_distances] 
            mutant_data += [mutant.IB, mutant.get_main_sequence()[0]]
            if self.multi_dataset:
                for dataset_name in self.dataset_order:
                    mutant_data += [mutant.by_dataset[dataset_name].total_read_count, 
                                    mutant.by_dataset[dataset_name].perfect_read_count]
            # add gene annotation, or a line with the right number of fields if gene annotation is missing
            #  (or if the mutant has no such attribute at all - possible for older-format datasets)
            if self.gene_annotation_header:
                try:
                    if any(mutant.gene_annotation):     mutant_data += mutant.gene_annotation
                    else:                               mutant_data += missing_gene_annotation_data
                except AttributeError:                  mutant_data += missing_gene_annotation_data
            OUTPUT.write('\t'.join([str(x) for x in mutant_data]))
            OUTPUT.write('\n')

    def print_detailed_RISCC_data(self, OUTPUT=sys.stdout, sort_data_by=None, max_distance=MAX_POSITION_DISTANCE):
        """ Write detailed RISCC data (all reads per mutant) to separate file.
        """
        # TODO docstring!

        # TODO should probably add header
        # TODO add annotation!
        # TODO change this to be a proper tab-separated file?

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        sorted_mutants = self._sort_data(sort_data_by)

        ### Quick summary
        N_total = len(self)
        # using sum because you can't do len on generators
        N_single_genomic = sum(1 for m in self if m.RISCC_N_distinct_regions(max_distance)[0]==1)
        N_single_chrom = sum(1 for m in self if m.RISCC_N_genomic_chromosomes==1)
        N_cassette = sum(1 for m in self if m.RISCC_N_distinct_regions(max_distance)[1]>0)
        OUTPUT.write("# %s mutants total; %s have one genomic location; "%(N_total, 
                                                                           value_and_percentages(N_single_genomic, [N_total]))
                     +"%s have locations on only one genomic chromosome; "%(value_and_percentages(N_single_chrom, [N_total]))
                     +"%s also have one or more cassette locations.\n"%(value_and_percentages(N_cassette, [N_total])) )

        # TODO add info on how many actually have ANY genomic-side reads!  And out of those, how many are confirmed vs not.

        ### Print data for all mutants
        for mutant in sorted_mutants:
            mutant.RISCC_print_detail(OUTPUT, max_distance)


# MAYBE-TODO should this be a Dataset method?  I just moved it from mutant_count_aln_IB+RISCC.py.
def save_dataset_files(dataset, outfile, verbosity_level=0, print_genome_side_details=True, 
                       count_cassette=True, count_other=True, sort_data_by='position', options="custom"):
    """ Print data to file, plus summary, pickle, and optionally detail files; optionally print summary to stdout.
    
    The options argument is only used to be printed in the header to make it clear how the file was generated - 
     it should be the applicable optparse options object if there is one, or a text message otherwise.
    """
    # print summary info to stdout if desired
    if verbosity_level>1: print "\nDATA SUMMARY:"
    if verbosity_level>0: dataset.print_summary(count_cassette=count_cassette, count_other=count_other)
    # print full data to outfile
    if verbosity_level>1: print "printing output - time %s."%time.ctime()
    outfile_basename = os.path.splitext(outfile)[0]
    summary_outfile = outfile_basename + '_summary.txt'
    pickled_outfile = outfile_basename + '.pickle'
    detail_outfile = outfile_basename + '_detail.txt'
    with open(summary_outfile,'w') as OUTFILE:
        write_header_data(OUTFILE,options)      # writes timestamp, generating program/options, folder and computer name, etc
        OUTFILE.write("### SUMMARY:\n")
        dataset.print_summary(OUTFILE, count_cassette=count_cassette, count_other=count_other)
        OUTFILE.write("# Basic mutant data in %s (or pickled in %s); detailed genome-side data in %s.\n"%(outfile, pickled_outfile, 
                                                                                                          detail_outfile))
    with open(outfile,'w') as OUTFILE:
        dataset.print_data(OUTPUT=OUTFILE, sort_data_by=sort_data_by, header_line=True, header_prefix='')
    pickle(dataset, pickled_outfile, protocol=-1)
    if print_genome_side_details:
        with open(detail_outfile,'w') as OUTFILE:
            dataset.print_detailed_RISCC_data(OUTPUT=OUTFILE, sort_data_by=sort_data_by)


def read_mutant_file(infile):
    """ Read mutant input file, return as new dataset (.pickle format only). """
    return unpickle(infile)


################################################### Unit-tests and main ########################################################

Fake_HTSeq_genomic_pos = Fake_deepseq_objects.Fake_HTSeq_genomic_pos
Fake_HTSeq_aln = Fake_deepseq_objects.Fake_HTSeq_alignment

class Testing_position_functionality(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """

    def test__Insertion_position(self):
        # different ways of making the same position or approximately same position, 
        #  by defining position_after, position_before, or both
        a_pos1 = Insertion_position('chr','+', '?-101')
        a_pos2 = Insertion_position('chr','+', full_position='?-101')
        a_pos3 = Insertion_position('chr','+', position_after=101)
        a_pos4 = Insertion_position('chr','+', position_after='101')
        a_pos5 = Insertion_position('chr','+', '-101')
        all_after_positions = [a_pos1, a_pos2, a_pos3, a_pos4, a_pos5]
        b_pos1 = Insertion_position('chr','+', '100-?')
        b_pos2 = Insertion_position('chr','+', full_position='100-?')
        b_pos3 = Insertion_position('chr','+', position_before=100)
        b_pos4 = Insertion_position('chr','+', position_before='100')
        b_pos5 = Insertion_position('chr','+', '100-')
        all_before_positions = [b_pos1, b_pos2, b_pos3, b_pos4, b_pos5]
        c_pos1 = Insertion_position('chr','+', '100-101')
        c_pos2 = Insertion_position('chr','+', full_position='100-101')
        c_pos3 = Insertion_position('chr','+', position_before=100, position_after=101)
        c_pos4 = Insertion_position('chr','+', position_before='100', position_after=101)
        c_pos5 = Insertion_position('chr','+', '100-101')
        all_both_positions = [c_pos1, c_pos2, c_pos3, c_pos4, c_pos5]
        all_positions = all_after_positions + all_before_positions + all_both_positions
        # most things are exactly the same for all these positions
        assert all([pos.chromosome == 'chr' for pos in all_positions])
        assert all([pos.strand == '+' for pos in all_positions])
        assert all([pos.min_position == 100 for pos in all_positions])
        assert all([pos.max_position == 101 for pos in all_positions])
        # position_before and position_after is defined for some and not others
        assert all([pos.full_position == '?-101' for pos in all_after_positions])
        assert all([pos.full_position == '100-?' for pos in all_before_positions])
        assert all([pos.full_position == '100-101' for pos in all_both_positions])
        assert all([pos.position_before is None for pos in all_after_positions])
        assert all([pos.position_before == 100 for pos in all_before_positions+all_both_positions])
        assert all([pos.position_after is None for pos in all_before_positions])
        assert all([pos.position_after == 101 for pos in all_after_positions+all_both_positions])
        # printing - str gives just the basic info, repr gives the function to generate the object
        assert all([str(pos) == 'chr + 100-?' for pos in all_before_positions])
        assert all([repr(pos) == "Insertion_position('chr', '+', full_position='100-?', immutable=False)" 
                    for pos in all_before_positions])
        # comparison - positions are only equal if they're exactly identical, so ?-101 and 100-101 aren't equal
        assert a_pos1 == a_pos2 == a_pos3 == a_pos4 == a_pos5
        assert b_pos1 == b_pos2 == b_pos3 == b_pos4 == b_pos5
        assert c_pos1 == c_pos2 == c_pos3 == c_pos4 == c_pos5
        assert a_pos1 != b_pos1 != c_pos1 != a_pos1     # need a_pos1 twice to get it compared to both of the others
        # sorting - based on chromosome names/numbers, then min_pos, then strand
        pos1 = Insertion_position('chr','+', '10-201')
        pos2 = Insertion_position('chr','+', '?-101')
        pos3 = Insertion_position('chr','+', '200-?')
        pos4 = Insertion_position('chr','-', '200-?')
        pos5 = Insertion_position('chr2','+', '?-101')
        pos6 = Insertion_position('chr12','+', '?-101')
        pos7 = Insertion_position('other_chr','+', '?-101')
        pos8 = Insertion_position('other_chr_4','+', '?-101')
        pos9 = Insertion_position('other_chr_13','+', '?-101')
        assert pos1 < pos2 < pos3 < pos4 < pos5 < pos6 < pos7 < pos8 < pos9
        ### copying position - same contents, different ID
        pos5 = pos1.copy()
        assert pos5 == pos1
        assert pos5 is not pos1
        ### invalid creation options
        # at least one position arg is required, and there must be at least one meaningful number overall
        self.assertRaises(ValueError, Insertion_position, 'chr','+')
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position=None, position_before=None)
        self.assertRaises(ValueError, Insertion_position, 'chr','+', position_before=None, position_after=None)
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position='?-?')
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position='-')
        # if providing full_position, can't provide another position arg
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_before=100)
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_after=100)
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_before=100,position_after=100)
        # full_position arg must be proper format
        for bad_fullpos in ['100', '100?', 100, (100,200), [100,200], True, False, '1-2-3', 'adaf', 'as-1']:
            self.assertRaises(ValueError, Insertion_position, 'chr','+', bad_fullpos)
        ### testing immutability/mutability and hashability
        # object is mutable and unhashable to start with unless specified otherwise
        pos1 = Insertion_position('chr','+', '0-?')
        assert str(pos1) == 'chr + 0-?'
        pos1.position_after = 5
        assert str(pos1) == 'chr + 0-5'
        pos1.new_attribute = 'test'
        assert pos1.new_attribute == 'test'
        self.assertRaises(MutantError, set, [pos1])
        self.assertRaises(MutantError, dict, [(pos1, 1)])
        # if we make it immutable, or initialize it as immutable from the start, it's immutable and hashable
        pos1.make_immutable()
        pos2 = Insertion_position('chr','+', '0-?', immutable=True)
        for pos in (pos1,pos2):
            # have to use execute workaround here to use self.assertRaises on statements instead of expressions;
            #  if this was python 2.7 I could just do "with self.assertRaises(MutantError): <bunch of statements>"
            def execute(S, context):  exec S in context
            self.assertRaises(MutantError, execute, "pos.chromosome = 'chr3'", locals())
            self.assertRaises(MutantError, execute, "pos.strand = '-'", locals())
            self.assertRaises(MutantError, execute, "pos.position_after = '5'", locals())
            self.assertRaises(MutantError, execute, "pos.position_before = '4'", locals())
            self.assertRaises(MutantError, execute, "pos.new_attribute = '4'", locals())
            set([pos1])
            dict([(pos1, 1)])
        # if we make it mutable, it's mutable and unhashable again
        pos1.make_mutable_REMEMBER_CLEANUP_FIRST()
        pos2.make_mutable_REMEMBER_CLEANUP_FIRST()
        for pos in pos1,pos2:
            pos.position_before = None
            pos.position_after = 500
            assert str(pos) == 'chr + ?-500'
            pos.new_attribute = 'test'
            assert pos.new_attribute == 'test'
            self.assertRaises(MutantError, set, [pos])
            self.assertRaises(MutantError, dict, [(pos, 1)])

    def test__HTSeq_pos_to_tuple(self):
        ### should raise exception for argument that's not an HTSeq_pos
        for bad_HTSeq_pos in [None, '', 'aaa', 0, 1, 0.65, [], {}, True, False, [1,2,3,4], ('C',4,6,'-')]:
            self.assertRaises(MutantError, HTSeq_pos_to_tuple, bad_HTSeq_pos)
        ### should raise exception for invalid HTSeq_pos argument (bad strand, or start not before end)
        for strand,start,end in [('+',3,3), ('-',3,3), ('-',4,3), ('+',100,3), ('a',1,2), (1,1,2)]:
            self.assertRaises(MutantError, HTSeq_pos_to_tuple, Fake_HTSeq_genomic_pos('C', strand, start, end))
        ### testing normal functionality: should return a (chrom,start_pos,end_pos,strand) tuple with the same chromosome/strand, 
        #    and start/end positions converted from 0-based end-exclusive to 1-based end-inclusive.
        # (so the HTSeq position of AA in AATT would be 0-2, and converted would be 1-2; of TT, 2-4, and 3-4.)
        for strand in '+-':
            assert HTSeq_pos_to_tuple(Fake_HTSeq_genomic_pos('C', strand, 3, 7)) == ('C', 4, 7, strand)
            for (start,end) in [(0,5), (0,100), (10,12), (5,44)]:
                assert HTSeq_pos_to_tuple(Fake_HTSeq_genomic_pos('C', strand, start, end)) == ('C', start+1, end, strand)

    def _check_unaligned_alns(self, aln_parse_function, *extra_args):
        """ Check that the function returns the right SPECIAL_POSITIONS object when given an unaligned HTSeq aln.  """
        fake_aln_unaligned_1     = Fake_HTSeq_aln('AAA', 'name', unaligned=True, optional_field_data={'XM':1})
        fake_aln_unaligned_2     = Fake_HTSeq_aln('AAA', 'name', unaligned=True, optional_field_data={})
        fake_aln_multi_aligned_1 = Fake_HTSeq_aln('AAA', 'name', unaligned=True, optional_field_data={'XM':2})
        fake_aln_multi_aligned_2 = Fake_HTSeq_aln('AAA', 'name', unaligned=True, optional_field_data={'XM':20})
        assert aln_parse_function(fake_aln_unaligned_1, *extra_args) == SPECIAL_POSITIONS.unaligned
        assert aln_parse_function(fake_aln_unaligned_2, *extra_args) == SPECIAL_POSITIONS.unaligned
        assert aln_parse_function(fake_aln_multi_aligned_1, *extra_args) == SPECIAL_POSITIONS.multi_aligned
        assert aln_parse_function(fake_aln_multi_aligned_2, *extra_args) == SPECIAL_POSITIONS.multi_aligned

    def test__parse_flanking_region_aln_or_pos(self):
        # very basic test with an HTSeq alignment - test__HTSeq_pos_to_tuple tests more details of this
        refpos = ('chr1',1,5,'+')
        pos = parse_flanking_region_aln_or_pos(Fake_HTSeq_aln('AAA', 'name', unaligned=False, pos=('chr1','+',0,5)))
        assert pos == refpos == HTSeq_pos_to_tuple(Fake_HTSeq_genomic_pos('chr1','+',0,5))
        # check straight position tuples
        assert parse_flanking_region_aln_or_pos(refpos) == refpos
        # check unaligned and multi-aligned positions
        self._check_unaligned_alns(parse_flanking_region_aln_or_pos)

    def _check_basic_pos_inputs(self, get_pos_function):
        """ Check basic inputs for function that takes (read_aln_or_pos, cassette_end, relative_read_direction). """
        # should raise exception for invalid argument (valid arguments: HTSeq position object or (chrom,start,end,strand) tuple
        #  (strand must be +/-, and start can't be after end)
        for bad_flanking_region in [None, '', 'aaa', 0, 1, 0.65, [], {}, True, False, ('C',2,3,4),('C',2,3,'x'),('C',3,2,'-')]:
            for cassette_end in SEQ_ENDS:
                for relative_read_direction in RELATIVE_READ_DIRECTIONS:
                    self.assertRaises(MutantError, get_pos_function, bad_flanking_region, cassette_end, relative_read_direction)
        # should raise exception for invalid cassette_end or relative_read_direction
        bad_vals = ['','aaa',0,1,[],{},None,True,False,'start','end','middle','read','leftmost','rightmost']
        for bad_val in bad_vals:
            for relative_read_direction in RELATIVE_READ_DIRECTIONS:
                self.assertRaises(MutantError, get_pos_function, ('C',1,5,'+'), bad_val, relative_read_direction)
            for cassette_end in SEQ_ENDS:
                self.assertRaises(MutantError, get_pos_function, ('C',1,5,'+'), cassette_end, bad_val)

    def test__get_insertion_pos_from_flanking_region_pos(self):
        ### check for invalid position or cassette_end/relative_read_direction inputs, and for unaligned aln inputs
        self._check_basic_pos_inputs(get_insertion_pos_from_flanking_region_pos)
        self._check_unaligned_alns(get_insertion_pos_from_flanking_region_pos, '5prime', 'inward')

        ### Testing normal functionality: should return an Insertion_position instance with the same chromosome, 
        #    and strand/position depending on the arguments in a somewhat complicated way.
        #   Checks based on the example at the end of "Possible read sides/directions" section in ../notes.txt
        #      (using HTSeq objects directly, because that's what the examples used)
        #    remember HTSeq position is 0-based and end-exclusive, and the position I want is 1-based end-inclusive!  
        #     So in the end the two relevant 1-based numbers end up being the same as the 0-based positions,
        #     because in the case of the start, min_position is actually start-1, and in the case of the end, we're adding 1
        #     to switch from 0-based to 1-based but then subtracting 1 to make the end the last base instead of the base after
        def _check_outputs(aln, side, rdir, expect_minpos, expect_strand, expect_chrom):
            pos = get_insertion_pos_from_flanking_region_pos(aln, side, rdir)
            assert (pos.min_position, pos.strand, pos.chromosome) == (expect_minpos, expect_strand, expect_chrom)
        for (start,end) in [(1,5), (1,100), (11,12), (6,44), (1000000001, 9000000000)]:
            fake_aln_plus  = Fake_HTSeq_aln('AA', 'x', pos=('C', '+', start, end))
            fake_aln_minus = Fake_HTSeq_aln('AA', 'x', pos=('C', '-', start, end))
            tuple_plus = ('C', start+1, end, '+')
            tuple_minus = ('C', start+1, end, '-')
            for input_plus in (fake_aln_plus, tuple_plus):
                _check_outputs(input_plus, '5prime', 'inward',  end,   '+', 'C')
                _check_outputs(input_plus, '5prime', 'outward', start, '-', 'C')
                _check_outputs(input_plus, '3prime', 'inward' , end,   '-', 'C')
                _check_outputs(input_plus, '3prime', 'outward', start, '+', 'C')
            for input_minus in (fake_aln_minus, tuple_minus):
                _check_outputs(input_minus,'5prime', 'inward',  start, '-', 'C')
                _check_outputs(input_minus,'5prime', 'outward', end,   '+', 'C')
                _check_outputs(input_minus,'3prime', 'inward' , start, '+', 'C')
                _check_outputs(input_minus,'3prime', 'outward', end,   '-', 'C')

    def test__get_RISCC_pos_from_read_pos(self):
        ### check for invalid position or cassette_end/relative_read_direction inputs, and for unaligned aln inputs
        self._check_basic_pos_inputs(get_RISCC_pos_from_read_pos)
        self._check_unaligned_alns(get_RISCC_pos_from_read_pos, '5prime', 'inward')
        ### check for full positions - similar to before
        def _check_outputs(aln, side, rdir, expect_minpos, expect_strand, expect_chrom):
            pos = get_RISCC_pos_from_read_pos(aln, side, rdir)
            assert (pos.min_position, pos.strand, pos.chromosome) == (expect_minpos, expect_strand, expect_chrom)
        for (start,end) in [(1,5), (1,100), (11,12), (6,44), (1000000001, 9000000000)]:
            fake_aln_plus  = Fake_HTSeq_aln('AA', 'x', pos=('C', '+', start, end))
            fake_aln_minus = Fake_HTSeq_aln('AA', 'x', pos=('C', '-', start, end))
            tuple_plus = ('C', start+1, end, '+')
            tuple_minus = ('C', start+1, end, '-')
            for input_plus in (fake_aln_plus, tuple_plus):
                _check_outputs(input_plus, '5prime', 'inward',  start, '+', 'C')
                _check_outputs(input_plus, '5prime', 'outward', end,   '-', 'C')
                _check_outputs(input_plus, '3prime', 'inward' , start, '-', 'C')
                _check_outputs(input_plus, '3prime', 'outward', end,   '+', 'C')
            for input_minus in (fake_aln_minus, tuple_minus):
                _check_outputs(input_minus,'5prime', 'inward',  end,   '-', 'C')
                _check_outputs(input_minus,'5prime', 'outward', start, '+', 'C')
                _check_outputs(input_minus,'3prime', 'inward' , end,   '+', 'C')
                _check_outputs(input_minus,'3prime', 'outward', start, '-', 'C')

    def test__find_gene_by_pos_gff3(self):
        with open('test_data/INPUT_gene-data-2_simple.gff3') as GENEFILE:
            chromosome_records_simple = list(GFF.parse(GENEFILE))
        assert len(chromosome_records_simple) == 2
        chr_rec = chromosome_records_simple[0]
        assert chr_rec.id == 'chrA'
        # the three last args to find_gene_by_pos_gff3 are detailed_features, nearest_genes_for_intergenic, and quiet.
        ### insertion in a gene, with/without details - should be the same regardless of nearest_genes_for_intergenic
        # convenience function to check insertions in genes:
        # - makes the position (actually multiple ones)
        # - checks all possible values of side and of nearest_genes_for_intergenic
        # - checks the version with and without detailed features (auto-generates the expected output without detailed features)
        def _check_pos_in_gene(min_pos, strand, expected_output, chrom='chrA'):
            chr_rec = [x for x in chromosome_records_simple if x.id == chrom][0]
            basic_features = 'gene_edge/?' if 'gene_edge' in expected_output[2] else '?'
            expected_output_basic_features = expected_output[:2] + [basic_features, expected_output[3]]
            for make_full_pos in (lambda x: '%d-?'%x, lambda x: '?-%d'%(x+1), lambda x: '%s-%s'%(x,x+1)):
                pos = Insertion_position(chrom, strand, full_position=make_full_pos(min_pos))
                for either in (True,False):
                    self.assertEquals(find_gene_by_pos_gff3(pos, chr_rec, True, either, True), expected_output)
                    self.assertEquals(find_gene_by_pos_gff3(pos, chr_rec, False, either, True), expected_output_basic_features)
        for (strand, orientation) in [('+', 'sense'), ('-', 'antisense')]:
            _check_pos_in_gene(100, strand, ['gene1_plus', orientation, "gene_edge/mRNA_edge/5'UTR", '0,600'])
            _check_pos_in_gene(130, strand, ['gene1_plus', orientation, "5'UTR",                    '30,570'])
            _check_pos_in_gene(200, strand, ['gene1_plus', orientation, "5'UTR/CDS",                '100,500'])
            _check_pos_in_gene(230, strand, ['gene1_plus', orientation, "CDS",                      '130,470'])
            _check_pos_in_gene(300, strand, ['gene1_plus', orientation, "CDS/intron",               '200,400'])
            _check_pos_in_gene(330, strand, ['gene1_plus', orientation, "intron",                   '230,370'])
            _check_pos_in_gene(400, strand, ['gene1_plus', orientation, "CDS/intron",               '300,300'])
            _check_pos_in_gene(500, strand, ['gene1_plus', orientation, "CDS",                      '400,200'])
            _check_pos_in_gene(600, strand, ['gene1_plus', orientation, "CDS/3'UTR",                '500,100'])
            _check_pos_in_gene(650, strand, ['gene1_plus', orientation, "3'UTR",                    '550,50'])
            _check_pos_in_gene(700, strand, ['gene1_plus', orientation, "gene_edge/mRNA_edge/3'UTR", '600,0'])
        for (strand, orientation) in [('+', 'antisense'), ('-', 'sense')]:
            _check_pos_in_gene(1100, strand, ['gene2_minus', orientation, "gene_edge/mRNA_edge/CDS", '600,0'])
            _check_pos_in_gene(1150, strand, ['gene2_minus', orientation, "CDS",                    '550,50'])
            _check_pos_in_gene(1200, strand, ['gene2_minus', orientation, "CDS/intron",             '500,100'])
            _check_pos_in_gene(1250, strand, ['gene2_minus', orientation, "intron",                 '450,150'])
            _check_pos_in_gene(1300, strand, ['gene2_minus', orientation, "CDS/intron",             '400,200'])
            _check_pos_in_gene(1350, strand, ['gene2_minus', orientation, "CDS",                    '350,250'])
            _check_pos_in_gene(1400, strand, ['gene2_minus', orientation, "CDS/intron",             '300,300'])
            _check_pos_in_gene(1450, strand, ['gene2_minus', orientation, "intron",                 '250,350'])
            _check_pos_in_gene(1500, strand, ['gene2_minus', orientation, "CDS/intron",             '200,400'])
            _check_pos_in_gene(1600, strand, ['gene2_minus', orientation, "CDS",                    '100,500'])
            _check_pos_in_gene(1700, strand, ['gene2_minus', orientation, "gene_edge/mRNA_edge/CDS", '0,600'])
        ### insertions in genes with multiple splice variants
        _check_pos_in_gene(100, '+', ['gene3_multi', 'sense', "gene_edge/mRNA_edge/CDS", '0,600'], 'chrB')
        _check_pos_in_gene(130, '+', ['gene3_multi', 'sense', "CDS",                     '30,570'], 'chrB')
        _check_pos_in_gene(230, '+', ['gene3_multi', 'sense', "intron",                  '130,470'], 'chrB')
        _check_pos_in_gene(330, '+', ['gene3_multi', 'sense', "CDS|intron",              '230,370'], 'chrB')
        _check_pos_in_gene(650, '+', ['gene3_multi', 'sense', "CDS",                     '550,50'], 'chrB')
        ### intergenic insertions - before a gene, between, after; with and without printing intergenic info
        for orientation in '+-':
          for make_full_pos in (lambda x: '%d-?'%x, lambda x: '?-%d'%(x+1)):
            for either in (True,False):
                [pos1, pos2, pos3, pos4] = [Insertion_position('chrA', orientation, full_position=make_full_pos(x)) 
                                      for x in (50, 99, 1000, 2000)]
                for pos in (pos1, pos2, pos3, pos4):
                    self.assertEquals(find_gene_by_pos_gff3(pos, chr_rec, either, False, True),
                                      [SPECIAL_GENE_CODES.not_found, '-', '-', '-'])
                self.assertEquals(find_gene_by_pos_gff3(pos1, chr_rec, either, True, True),
                                  ['gene1_plus', 'upstream', 'intergenic', '50'])
                self.assertEquals(find_gene_by_pos_gff3(pos2, chr_rec, either, True, True),
                                  ['gene1_plus', 'upstream', 'intergenic', '1'])
                self.assertEquals(find_gene_by_pos_gff3(pos3, chr_rec, either, True, True), 
                                  ['gene1_plus & gene2_minus', 'downstream & downstream', 'intergenic', '300 & 100'])
                self.assertEquals(find_gene_by_pos_gff3(pos4, chr_rec, either, True, True),
                                  ['gene2_minus', 'upstream', 'intergenic', '300'])
        # LATER-TODO add tests for weird cases based on test_data/INPUT_gene-data-1_all-cases.gff3


class Testing_Insertional_mutant(unittest.TestCase):
    """ Unit-tests for the Insertional_mutant class and its methods. """

    def test__init(self):
        """ For single, multi-dataset and readcount-only mutants. """
        for chromosome in ['chr1', 'chromosome_2', 'chrom3', 'a', 'adfads', '100', 'scaffold_88']:
            for strand in ['+','-']:
                for position in [1,2,5,100,10000,4323423]:
                    ins_pos_5prime = Insertion_position(chromosome,strand,position_before=position)
                    ins_pos_3prime = Insertion_position(chromosome,strand,position_after=position)
                    # test "normal" mutants - check all details, including position
                    mutant_5prime = Insertional_mutant(insertion_position=ins_pos_5prime)
                    mutant_3prime = Insertional_mutant(insertion_position=ins_pos_3prime)
                    mutant_readcount_only = Insertional_mutant_readcount_only()
                    mutant_multi_dataset = Insertional_mutant_multi_dataset(insertion_position=ins_pos_5prime)
                    # test position details (only for the two "normal" mutants)
                    assert mutant_5prime.position.min_position == position
                    assert mutant_3prime.position.min_position == position-1
                    assert mutant_5prime.position.max_position == position+1
                    assert mutant_3prime.position.max_position == position
                    assert mutant_5prime.position.full_position == "%s-?"%(position)
                    assert mutant_3prime.position.full_position == "?-%s"%position
                    # test non-readcount-related info for all mutants except mutant_readcount_only
                    for mutant in [mutant_5prime, mutant_3prime, mutant_multi_dataset]:
                        assert mutant.position.chromosome == chromosome
                        assert mutant.position.strand == strand
                        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
                        assert mutant.orientation == '?'
                        assert mutant.gene_feature == '?'
                        assert mutant.gene_distances == '?'
                    # test readcount-related info for all mutants except mutant_multi_dataset
                    for mutant in [mutant_5prime, mutant_3prime, mutant_readcount_only]:
                        assert mutant.total_read_count == 0
                        assert mutant.perfect_read_count == 0
                        assert mutant.sequences_counts_positions_errors == {}
                    # test readcount-related info for mutant_multi_dataset
                    assert all([x.total_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.perfect_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.sequences_counts_positions_errors == {} for x in mutant_multi_dataset.by_dataset.values()])

    def test__add_read(self):
        """ For single, multi-dataset and readcount-only mutants. """
        # using fake HTSeq alignment class from deepseq_utilities; defining one perfect and one imperfect alignment
        # note: the detailed mutation-counting methods are imported from deepseq_utilities and unit-tested there.
        position =  Insertion_position('chr1',  '+', position_before=3)
        perfect_aln = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        imperfect_aln = Fake_HTSeq_aln(seq='GGG', optional_field_data={'NM':1})
        # adding perfect and imperfect to mutant increases all the counts as expected
        mutant = Insertional_mutant(insertion_position=position)
        mutant.add_read(perfect_aln, read_count=3, position=position)
        assert mutant.total_read_count == mutant.perfect_read_count == 3
        assert mutant.sequences_counts_positions_errors == {'AAA': [3, position, 0]}
        mutant.add_read(imperfect_aln, read_count=1, position=position)
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 3
        assert mutant.sequences_counts_positions_errors == {'AAA': [3, position, 0], 'GGG': [1, position, 1]}
        # same for a multi-dataset mutant - this time we need to specify which dataset we're adding to
        mutant = Insertional_mutant_multi_dataset(insertion_position=position)
        assert len(mutant.by_dataset) == 0
        mutant.add_read(perfect_aln, read_count=3, dataset_name='d1', position=position)
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].sequences_counts_positions_errors == {'AAA': [3, position, 0]}
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d1', position=position)
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].sequences_counts_positions_errors == {'AAA': [3, position, 0], 'GGG': [1, position, 1]}
        # now adding a read to another dataset - nothing changes in dataset d1, but we have new dataset d2 numbers
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d2', position=position)
        assert len(mutant.by_dataset) == 2
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d2'].total_read_count == 1
        assert mutant.by_dataset['d2'].perfect_read_count == 0
        assert mutant.by_dataset['d2'].sequences_counts_positions_errors == {'GGG': [1, position, 1]}
        # it should be impossible to add a read to a specific dataset in a single-dataset mutant 
        mutant = Insertional_mutant(insertion_position=position)
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3, dataset_name='d1')
        # it should be impossible to add a read to a multi-dataset mutant without giving a dataset_name
        mutant = Insertional_mutant_multi_dataset(insertion_position=position)
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3)

    def test__decide_and_check_position(self):
        # note that add_read DOESN'T check position - need to run mutant.decide_and_check_position.
        position0            = Insertion_position('chr1',  '+', position_before=3)
        position_diff_chrom  = Insertion_position('chr33', '+', position_before=3)
        position_diff_strand = Insertion_position('chr1',  '-', position_before=3)
        position_10bp_away   = Insertion_position('chr1',  '+', position_before=13)
        position_2bp_away    = Insertion_position('chr1',  '+', position_before=5)
        alnAAA = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        alnAAT = Fake_HTSeq_aln(seq='AAT', optional_field_data={'NM':0})
        mutant = Insertional_mutant()
        assert mutant.position == SPECIAL_POSITIONS.unknown
        assert not mutant.decide_and_check_position()
        assert mutant.position == SPECIAL_POSITIONS.unknown
        mutant.add_read(alnAAA, position0)
        assert not mutant.decide_and_check_position()
        assert mutant.position == position0
        # can always add undefined positions
        for curr_position in SPECIAL_POSITIONS.all_undefined:
            mutant.add_read(alnAAT, curr_position, 10)
            assert not mutant.decide_and_check_position()
            assert mutant.position == position0
        # adding positions that don't match the first one gives an error
        # (re-making mutant in each iteration here, since each bad position gets added BEFORE being checked)
        for curr_position in (position_diff_chrom, position_diff_strand, position_10bp_away, position_2bp_away):
            mutant = Insertional_mutant()
            mutant.add_read(alnAAA, position0)
            mutant.add_read(alnAAT, curr_position)
            assert mutant.decide_and_check_position() == True
        # but if I allow some distance, position_2bp_away can be okay, 
        #  and the main position switches to that and back depending on read counts.
        for max_allowed_dist in (2, 5, 10):
            mutant = Insertional_mutant()
            mutant.add_read(alnAAA, position0, 2)
            mutant.add_read(alnAAT, position_2bp_away, 1)
            assert not mutant.decide_and_check_position(max_allowed_dist)
            assert mutant.position == position0
            mutant.add_read(alnAAT, position_2bp_away, 2)
            assert not mutant.decide_and_check_position(max_allowed_dist)
            assert mutant.position == position_2bp_away
            mutant.add_read(alnAAA, position0, 2)
            assert not mutant.decide_and_check_position(max_allowed_dist)
            assert mutant.position == position0
            assert mutant.decide_and_check_position(0) == True
        for max_allowed_dist in (0, 1):
            mutant = Insertional_mutant()
            mutant.add_read(alnAAA, position0)
            mutant.add_read(alnAAT, position_2bp_away)
            assert mutant.decide_and_check_position(max_allowed_dist) == True
        # you can also make a mutant with a missing position and then add whatever position you want later
        for special_position in SPECIAL_POSITIONS.all_undefined:
            for curr_position in (position0, position_diff_chrom, position_diff_strand, position_10bp_away, position_2bp_away):
                mutant = Insertional_mutant()
                mutant.add_read(alnAAA, special_position, 10)
                mutant.add_read(alnAAT, curr_position, 1)
                assert not mutant.decide_and_check_position()
                assert mutant.position == curr_position
        # TODO add test for ratio_to_ignore=100!

    def test__update_gene_info(self):
        mutant = Insertional_mutant(insertion_position=Insertion_position('chr','+',position_before=3))
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == mutant.gene_distances == '?'
        # updating no-info mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?', '?')
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == mutant.gene_distances == '?'
        # updating no-info mutant with useful info - update goes through
        mutant.update_gene_info('gene1', '+', 'f', (10,10))
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        assert mutant.gene_distances == (10,10)
        # updating info-containing mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?', '?')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        assert mutant.gene_distances == (10,10)
        # updating info-containing mutant with same info - no change
        mutant.update_gene_info('gene1', '+', 'f', (10,10))
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        assert mutant.gene_distances == (10,10)
        # updating info-containig mutant with OTHER info - exception
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene2', '+', 'f', (10,10))
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '-', 'f', (10,10))
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '+', 'g', (10,10))

    def test__get_main_sequence(self):
        # single-dataset mutant
        position = Insertion_position('chr','+',position_before=3)
        perfect_aln = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        imperfect_aln = Fake_HTSeq_aln(seq='GGG', optional_field_data={'NM':1})
        imperfect_aln2 = Fake_HTSeq_aln(seq='CCC', optional_field_data={'NM':1})
        mutant = Insertional_mutant(insertion_position=position)
        assert mutant.get_main_sequence() == ('',0)
        assert mutant.get_main_sequence(1) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_read(perfect_aln, read_count=1, position=position)
        mutant.add_read(imperfect_aln, read_count=2, position=position)
        assert mutant.sequences_counts_positions_errors == {'AAA': [1, position, 0], 'GGG': [2, position, 1]}
        assert mutant.get_main_sequence() == ('GGG',2)
        assert mutant.get_main_sequence(1) == ('GGG',2)
        assert mutant.get_main_sequence(2) == ('AAA',1)
        assert mutant.get_main_sequence(3) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_read(perfect_aln, read_count=2, position=position)
        mutant.add_read(imperfect_aln2, read_count=1, position=position)
        assert mutant.sequences_counts_positions_errors == {'AAA': [3, position, 0], 'GGG': [2, position, 1], 
                                                            'CCC': [1, position, 1]}
        assert mutant.get_main_sequence() == ('AAA',3)
        assert mutant.get_main_sequence(1) == ('AAA',3)
        assert mutant.get_main_sequence(2) == ('GGG',2)
        assert mutant.get_main_sequence(3) == ('CCC',1)
        assert mutant.get_main_sequence(4) == ('',0)
        assert mutant.get_main_sequence(5) == ('',0)
        # multi-dataset mutant - getting the top sequence from single dataset or from all datasets added together
        mutant = Insertional_mutant_multi_dataset(insertion_position=Insertion_position('chr','+',position_before=3))
        mutant.add_read(perfect_aln, read_count=3, position=position, dataset_name='d1')
        mutant.add_read(imperfect_aln, read_count=2, position=position, dataset_name='d1')
        mutant.add_read(imperfect_aln2, read_count=3, position=position, dataset_name='d2')
        mutant.add_read(imperfect_aln, read_count=2, position=position, dataset_name='d2')
        assert mutant.get_main_sequence(1, dataset_name='d1') == ('AAA',3)
        assert mutant.get_main_sequence(1, dataset_name='d2') == ('CCC',3)
        assert mutant.get_main_sequence(1) == ('GGG',4)     # GGG is the most common sequence if we add both datasets
        # unaligned reads - make sure they get skipped correctly with aligned_only=True
        unaligned = Fake_HTSeq_aln(seq='CAC', unaligned=True)
        perfect_aln = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        mutant = Insertional_mutant(insertion_position=position)
        mutant.add_read(perfect_aln, read_count=1, position=position)
        mutant.add_read(unaligned, read_count=10)
        assert mutant.get_main_sequence(aligned_only=False) == ('CAC',10)
        assert mutant.get_main_sequence(aligned_only=True) == ('AAA',1)

    ### RISCC data stuff

    @staticmethod
    def _make_pos(pos):
        """ Convenience function to make the args to Fake_HTSeq_genomic_pos from an Insertion_position """
        return pos.chromosome, pos.strand, pos.min_position, pos.min_position+20

    def test__add_RISCC_read(self):
        """ Also tests RISCC_max_confirmed_distance """
        # make the cassette-side read
        pos0 = Insertion_position('chr1','+',position_before=100, immutable=True)
        aln0 = Fake_HTSeq_aln(seq='AAA', readname='read1', pos=self._make_pos(pos0), optional_field_data={'NM':0})
        mutant = Insertional_mutant(insertion_position=pos0)
        assert mutant.total_read_count == 0
        assert isnan(mutant.RISCC_max_confirmed_distance(10))
        mutant.add_read(aln0)
        assert mutant.total_read_count == 1
        assert isnan(mutant.RISCC_max_confirmed_distance(10))
        # add unaligned or multi-aligned RISCC read - the max confirmed distance should still be NaN
        mutant.add_RISCC_read('AAA', SPECIAL_POSITIONS.unaligned)
        assert isnan(mutant.RISCC_max_confirmed_distance(10))
        mutant.add_RISCC_read('AAC', SPECIAL_POSITIONS.multi_aligned)
        assert isnan(mutant.RISCC_max_confirmed_distance(10))
        # add RISCC reads - a confirming one, a non-confirming one, and a weird one.
        pos1 = Insertion_position('chr1','+',position_before=0, immutable=True)
        mutant.add_RISCC_read('CAA', pos1)
        pos2 = Insertion_position('chr2','+',position_before=0, immutable=True)
        mutant.add_RISCC_read('CAC', pos2)
        pos3 = Insertion_position('chr1','-',position_before=0, immutable=True)
        mutant.add_RISCC_read('CAG', pos3)
        assert (len(mutant.RISCC_genome_side_aligned_reads), len(mutant.RISCC_genome_side_unaligned_reads)) == (3, 2)
        assert sorted(data[0].chromosome for data in mutant.RISCC_genome_side_aligned_reads.values()) == 'chr1 chr1 chr2'.split()
        assert sorted(data[0].strand for data in mutant.RISCC_genome_side_aligned_reads.values()) == '+ + -'.split()
        assert mutant.RISCC_max_confirmed_distance(1000) == 100
        assert mutant.RISCC_max_confirmed_distance(100) == 100
        assert mutant.RISCC_max_confirmed_distance(10) == 0
        # same seq can't be added with different defined position, but can be added with unknown - TODO this doesn't work for now
        #self.assertRaises(MutantError, mutant.add_RISCC_read, 'AAA', pos1)
        #self.assertRaises(MutantError, mutant.add_RISCC_read, 'CAA', pos2)
        #self.assertRaises(MutantError, mutant.add_RISCC_read, 'CAA', pos2)
        mutant.add_RISCC_read('AAA', SPECIAL_POSITIONS.unaligned)
        mutant.add_RISCC_read('CAA', pos1)
        mutant.add_RISCC_read('CAG', pos3)
        mutant.add_RISCC_read('AAA', SPECIAL_POSITIONS.unknown)
        mutant.add_RISCC_read('CAA', SPECIAL_POSITIONS.unknown)
        mutant.add_RISCC_read('CAG', SPECIAL_POSITIONS.unknown)

    def test__improve_best_RISCC_read(self):
        """ Also tests RISCC_max_confirmed_distance """
        # make the cassette-side read
        pos0 = Insertion_position('chr1','+',position_before=100, immutable=True)
        aln0 = Fake_HTSeq_aln(seq='AAA', readname='read1', pos=self._make_pos(pos0), optional_field_data={'NM':0})
        mutant = Insertional_mutant(insertion_position=pos0)
        assert mutant.position.min_position == 100
        assert mutant.total_read_count == 0
        assert isnan(mutant.RISCC_max_confirmed_distance(1000))
        mutant.add_read(aln0, pos0)
        assert mutant.total_read_count == 1
        assert isnan(mutant.RISCC_max_confirmed_distance(1000))
        # add bad read first, then better ones, make sure they're replaced (or not, if both are equally bad)
        #  (max confirmed distance should change from NaN to 0 once there are some non-confirming reads, 
        #   then to a positive value once there are some confirming reads)
        pos1 = Insertion_position('chr2','+',position_before=10, immutable=True)        # bad - wrong chromosome
        mutant.improve_best_RISCC_read('AAA', pos1, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 10
        assert mutant.RISCC_max_confirmed_distance(1000) == 0
        pos2 = Insertion_position('chr1','-',position_before=11, immutable=True)        # bad - wrong strand
        mutant.improve_best_RISCC_read('AAA', pos2, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 10
        assert mutant.RISCC_max_confirmed_distance(1000) == 0
        pos3 = Insertion_position('chr1','+',position_before=12, immutable=True)        # bad - distance larger than max_distance
        mutant.improve_best_RISCC_read('AAA', pos3, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 10
        assert mutant.RISCC_max_confirmed_distance(1000) == 0
        pos4 = Insertion_position('chr1','+',position_before=80, immutable=True)        # good - 20bp distance
        mutant.improve_best_RISCC_read('AAA', pos4, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 80
        assert mutant.RISCC_max_confirmed_distance(1000) == 20
        pos5 = Insertion_position('chr1','+',position_before=70, immutable=True)        # better - 30bp distance
        mutant.improve_best_RISCC_read('AAA', pos5, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30          
        mutant.improve_best_RISCC_read('AAA', pos3, max_distance=500)   # trying pos3 again with higher max_distance
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 12
        assert mutant.RISCC_max_confirmed_distance(1000) == 88
        # add good read first, then worse ones, make sure it's NOT replaced (remember to start with new mutant!)
        mutant = Insertional_mutant(insertion_position=pos0)
        mutant.add_read(aln0)
        mutant.improve_best_RISCC_read('AAA', pos5, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        mutant.improve_best_RISCC_read('AAA', pos1, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        mutant.improve_best_RISCC_read('AAA', pos2, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        mutant.improve_best_RISCC_read('AAA', pos3, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        mutant.improve_best_RISCC_read('AAA', pos4, max_distance=50)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        # make sure wrong-strand closer-distance new read won't replace old right-strand one! (this was a bug once)
        pos0 = Insertion_position('chr1','+',position_before=100, immutable=True)
        aln0 = Fake_HTSeq_aln(seq='AAA', readname='read1', pos=self._make_pos(pos0), optional_field_data={'NM':0})
        mutant = Insertional_mutant(insertion_position=pos0)
        mutant.add_read(aln0)
        pos5 = Insertion_position('chr1','+',position_before=70, immutable=True)
        mutant.improve_best_RISCC_read('AAA', pos5, max_distance=1000)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30
        pos6 = Insertion_position('chr1','-',position_before=11, immutable=True)
        mutant.improve_best_RISCC_read('AAA', pos6, max_distance=1000)
        assert len(mutant.RISCC_genome_side_aligned_reads) == 1
        assert mutant.RISCC_genome_side_aligned_reads.values()[0][0].min_position == 70
        assert mutant.RISCC_max_confirmed_distance(1000) == 30

    ### multi-dataset mutants

    def test__add_other_mutant_as_dataset(self):
        position = Insertion_position('chr','+',position_before=3, immutable=True)
        perfect_aln = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        imperfect_aln = Fake_HTSeq_aln(seq='GGG', optional_field_data={'NM':1})
        mutant1 = Insertional_mutant(insertion_position=position)
        mutant1.add_read(perfect_aln, read_count=2, position=position)
        mutant2 = Insertional_mutant(insertion_position=position)
        mutant2.add_read(imperfect_aln, read_count=1, position=position)
        # adding a mutant to a single-dataset mutant should fail
        if hasattr(mutant1, 'add_other_mutant_as_dataset'):
            self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
        # adding a mutant to a multi-dataset mutant should work
        mutantM = Insertional_mutant_multi_dataset(insertion_position = mutant1.position)
        mutantM.add_other_mutant_as_dataset(mutant1, 'd1')
        mutantM.add_other_mutant_as_dataset(mutant2, 'd2')
        assert mutantM.by_dataset['d1'].total_read_count == 2
        assert mutantM.by_dataset['d1'].perfect_read_count == 2
        assert mutantM.by_dataset['d1'].sequences_counts_positions_errors == {'AAA': [2, position, 0]}
        assert mutantM.by_dataset['d2'].total_read_count == 1
        assert mutantM.by_dataset['d2'].perfect_read_count == 0
        assert mutantM.by_dataset['d2'].sequences_counts_positions_errors == {'GGG': [1, position, 1]}
        # can't add overwrite an existing dataset name, unless overwrite==True
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant2, 'd2')
        mutantM.add_other_mutant_as_dataset(mutant2, 'd2', overwrite=True)
        # if the two mutants have different positions, it should fail, unless check_constant_data=False
        mutant3 = Insertional_mutant(Insertion_position('chr','+',position_before=5, immutable=True))
        mutant4 = Insertional_mutant(Insertion_position('chr2','+',position_before=3, immutable=True))
        mutant5 = Insertional_mutant(Insertion_position('chr','-',position_before=3, immutable=True))
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant3, 'd3', check_constant_data=True)
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant4, 'd4', check_constant_data=True)
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant5, 'd5', check_constant_data=True)
        mutantM.add_other_mutant_as_dataset(mutant3, 'd3', check_constant_data=False)
        mutantM.add_other_mutant_as_dataset(mutant4, 'd4', check_constant_data=False)
        mutantM.add_other_mutant_as_dataset(mutant5, 'd5', check_constant_data=False)

    def _make_multi_mutant(self):
        position = Insertion_position('chr','+',position_before=3, immutable=True)
        perfect_aln = Fake_HTSeq_aln(seq='AAA', optional_field_data={'NM':0})
        imperfect_aln = Fake_HTSeq_aln(seq='GGG', optional_field_data={'NM':1})
        mutant = Insertional_mutant_multi_dataset(insertion_position=position)
        mutant.add_read(perfect_aln, read_count=1, position=position, dataset_name='d1')
        mutant.add_read(imperfect_aln, read_count=1, position=position, dataset_name='d1')
        mutant.add_read(imperfect_aln, read_count=1, position=position, dataset_name='d2')
        return mutant, position

    def _check_single_mutants(self, mutant, position, new_mutant_1, new_mutant_2):
        for new_mutant in (new_mutant_1, new_mutant_2):
            assert new_mutant.position == mutant.position
            assert new_mutant.gene == mutant.gene
        assert new_mutant_1.total_read_count == 2
        assert new_mutant_1.perfect_read_count == 1
        assert new_mutant_1.sequences_counts_positions_errors == {'AAA': [1, position, 0], 'GGG': [1, position, 1]}
        assert new_mutant_2.total_read_count == 1
        assert new_mutant_2.perfect_read_count == 0
        assert new_mutant_2.sequences_counts_positions_errors == {'GGG': [1, position, 1]}

    def test__give_single_dataset_mutant(self):
        mutant, position = self._make_multi_mutant()
        new_mutant_1 = mutant.give_single_dataset_mutant('d1')
        new_mutant_2 = mutant.give_single_dataset_mutant('d2')
        self._check_single_mutants(mutant, position, new_mutant_1, new_mutant_2)
        # trying to extract an inexistent dataset should fail, unless force==True
        self.assertRaises(MutantError, mutant.give_single_dataset_mutant, 'd0', force=False)
        new_mutant_0 = mutant.give_single_dataset_mutant('d0', force=True)
        assert new_mutant_0.position == mutant.position
        assert new_mutant_0.gene == mutant.gene
        assert new_mutant_0.total_read_count == 0
        assert new_mutant_0.perfect_read_count == 0
        assert new_mutant_0.sequences_counts_positions_errors == {}
        # trying to extract a dataset from a single-dataset mutant should fail
        mutant1 = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3, immutable=True))
        self.assertRaises(MutantError, mutant1.give_single_dataset_mutant, 'd2')

    def test__give_all_single_dataset_mutants(self):
        mutant, position = self._make_multi_mutant()
        # extracting two mutants and checking the values
        all_single_mutants = mutant.give_all_single_dataset_mutants()
        assert len(all_single_mutants) == 2
        new_mutant_1 = all_single_mutants['d1']
        new_mutant_2 = all_single_mutants['d2']
        self._check_single_mutants(mutant, position, new_mutant_1, new_mutant_2)


class Testing_Insertional_mutant_pool_dataset(unittest.TestCase):
    """ Unit-tests for the Insertional_mutant_pool_dataset class and its methods. """

    @staticmethod
    def _make_test_mutant_dataset(positions_and_readcounts_string, raw_chrom_names=False):
        """ Help method to quickly make a dataset based on a string of mutant positions/readcounts.

        Comma-separated mutants, first word is position, second is readcount.  "+100 5, A-100 10/10, cassette+400 1"

        Position is chromosome+strand+minpos, with chromosome optional and more complicated.  If raw_chrom_names, 
         just take the chromosome name as given; otherwise prepend 'chromosome_' to what is given, or just use
         'chromosome_1' if chromosome string is empty.

        If readcount is a single number, it's the total read count; if it's two numbers, it's total/perfect.
        """
        dataset = Insertional_mutant_pool_dataset()
        if not positions_and_readcounts_string: 
            return dataset
        for N, string in enumerate(positions_and_readcounts_string.split(', ')):
            raw_pos, readcount = string.split(' ')
            if '/' in readcount:    readcount, perfect = [int(x) for x in readcount.split('/')]
            else:                   readcount = perfect = int(readcount)
            assert readcount >= perfect, "In mutant string %s, perfect readcount is over total - not allowed!"%string
            if '+' in raw_pos:      strand = '+'
            elif '-' in raw_pos:    strand = '-'
            else:                   raise Exception("Short-position %s has no strand!"%raw_pos)
            chrom, pos = raw_pos.split(strand)
            pos = int(pos)
            if not raw_chrom_names:
                if chrom:   chrom = 'chromosome_%s'%chrom
                else:       chrom = 'chromosome_1'
            elif not chrom:
                raise Exception("Short-position %s has no chromosome name - can't use with raw_chrom_names!")
            full_pos = Insertion_position(chrom, strand, position_before=pos, immutable=True)
            mutant = Insertional_mutant(IB=str(N), insertion_position=full_pos)
            mutant.total_read_count = readcount
            mutant.perfect_read_count = perfect
            dataset.add_mutant(mutant)
        return dataset

    def test____init__(self):
        # basic check of empty dataset
        for cassette_end in SEQ_ENDS+['?']:
            for relative_read_direction in RELATIVE_READ_DIRECTIONS+['?']:
                data = Insertional_mutant_pool_dataset(cassette_end, relative_read_direction)
                assert data.summary.cassette_end == cassette_end
                assert data.summary.relative_read_direction == relative_read_direction
                assert len(data) == 0
                assert data.summary.processed_read_count == data.summary.aligned_read_count == data.summary.perfect_read_count == 0
                assert data.summary.non_aligned_read_count == 0
                assert data.summary.discarded_read_count == None
                assert data.summary.ignored_region_read_counts == {}
                assert sum(data.summary.strand_read_counts.values()) == 0
                assert data.summary.mutants_in_genes == data.summary.mutants_not_in_genes\
                        == data.summary.mutants_undetermined == 0
                assert data.summary.mutant_counts_by_orientation == data.summary.mutant_counts_by_feature == {}
        # check bad cassette_end and relative_read_direction values
        for cassette_end in [True, False, None, 0, 1, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, cassette_end, '?')
        for relative_read_direction in [True, False, None, 0, 1, 'reverse', 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, '?', relative_read_direction)

    def _check_RISCC1_outputs(self, dataset):
        assert len(dataset) == 1
        mutant = dataset.get_mutant('CCCC')
        assert mutant.IB == 'CCCC'
        assert str(mutant.position) == 'cassette_3_confirming + ?-101'
        assert mutant.total_read_count == mutant.perfect_read_count == 3

    def test__add_RISCC_alignment_files_to_data(self):
        # testing with IB clustering (.py file and .pickle file)
        infiles = ['test_data/INPUT_RISCC1-alignment-cassette-side.sam', 'test_data/INPUT_RISCC1-alignment-genome-side.sam', 
                   'test_data/INPUT_RISCC1-IBs.fq', None, {'CCCC': set('CCCC CCCG'.split())} ]
        dataset = Insertional_mutant_pool_dataset('3prime', 'outward')
        dataset.add_RISCC_alignment_files_to_data(*infiles, quiet=True)
        self._check_RISCC1_outputs(dataset)
        picklefile = 'test_data/INPUT_RISCC1-IB-clusters.pickle'
        pickle(infiles[-1], picklefile)
        dataset = Insertional_mutant_pool_dataset('3prime', 'outward')
        dataset.add_RISCC_alignment_files_to_data(*infiles[:-1], IB_cluster_file=picklefile, quiet=True)
        self._check_RISCC1_outputs(dataset)
        os.unlink(picklefile)
        # testing without IB clustering
        dataset = Insertional_mutant_pool_dataset('3prime', 'outward')
        dataset.add_RISCC_alignment_files_to_data(*infiles[:-1], quiet=True)
        assert len(dataset) == 2
        mutant = dataset.get_mutant('CCCC')
        assert mutant.IB == 'CCCC'
        assert str(mutant.position) == 'cassette_3_confirming + ?-101'
        assert mutant.total_read_count == mutant.perfect_read_count == 2
        mutant = dataset.get_mutant('CCCG')
        assert mutant.IB == 'CCCG'
        assert str(mutant.position) == 'cassette_3_confirming + ?-101'
        assert mutant.total_read_count == mutant.perfect_read_count == 1
        # testing max_allowed_cassette_side_dist
        infiles = ['test_data/INPUT_RISCC2-alignment-cassette-side.sam'] + infiles[1:]
        dataset = Insertional_mutant_pool_dataset('3prime', 'outward')
        for dist in (1, 2, 5, 10, 100):
            dataset.add_RISCC_alignment_files_to_data(*infiles, max_allowed_cassette_side_dist=dist, quiet=True)
            assert len(dataset) == 1
            dataset.add_RISCC_alignment_files_to_data(*infiles, max_allowed_cassette_side_dist=0, quiet=True)
            assert len(dataset) == 0
        # TODO add more tests - but there's also a run-test for this.
        # TODO add more complicated stuff to input files! see comments in input files.

    def test__remove_mutants_below_readcount(self):
        positions_and_readcounts_raw = "A+100 5/5, A-200 2/2, A+300 1/1, A-400 5/2, A+500 5/1"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        assert set([m.position.min_position for m in dataset]) == set([100,200,300,400,500])
        dataset.remove_mutants_below_readcount(2, perfect_reads=False)
        assert set([m.position.min_position for m in dataset]) == set([100,200,400,500])
        dataset.remove_mutants_below_readcount(4, perfect_reads=False)
        assert set([m.position.min_position for m in dataset]) == set([100,400,500])
        dataset.remove_mutants_below_readcount(2, perfect_reads=True)
        assert set([m.position.min_position for m in dataset]) == set([100,400])
        dataset.remove_mutants_below_readcount(4, perfect_reads=True)
        assert set([m.position.min_position for m in dataset]) == set([100])

    def test__pickle_unpickle(self):
        """ What to make sure that any changes I make still allow the datasets to be pickled/unpickled correctly, 
        since some stuff like lambdas and __slots__ etc interferes with that. """
        pass
    # TODO implement

    # TODO lots of functions without unit-tests here!  Need to decide between unit-tests, run-tests and a mix.


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])

