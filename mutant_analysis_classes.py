#!/usr/bin/env python2.7
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011
"""

from __future__ import division
# basic libraries
import sys, re
import unittest
from collections import defaultdict
import itertools
import copy
import random
# other libraries
from numpy import median, round, isnan
import HTSeq
from BCBio import GFF
# my modules
from general_utilities import split_into_N_sets_by_counts, add_dicts_of_ints, sort_lists_inside_dict, keybased_defaultdict, value_and_percentages, FAKE_OUTFILE, NaN, nan_func, merge_values_to_unique, unpickle
from basic_seq_utilities import SEQ_ENDS, SEQ_STRANDS, SEQ_DIRECTIONS, SEQ_ORIENTATIONS, position_test_contains, position_test_overlap, chromosome_sort_key, get_seq_count_from_collapsed_header
from deepseq_utilities import check_mutation_count_try_all_methods
from parse_annotation_file import parse_gene_annotation_file


### Constants
MULTIPLE_GENE_JOIN = ' & '

class SPECIAL_GENE_CODES(object):
    not_determined = "gene_unknown"
    chromosome_not_in_reference = "unknown_chrom"
    not_found = "no_gene_found"
# it seems like I have to set SPECIAL_GENE_CODES.all_codes afterward because I can't access __dict__ from inside the class
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
    if is_cassette_chromosome(chromosome_name):
        return False
    if chromosome_name.startswith('chr') or chromosome_name.startswith('scaffold'):
        return False
    return True

# MAYBE-TODO may want to put the is_*_chromosome functionality under command-line user control someday?  Provide option to give list of cassette chromosome names and "other" ones (or cassette ones and prefixes for genomic ones?)


# MAYBE-TODO it might be good to split this file into multiple files at some point?  At least Insertion_position/etc.

############################ Functions/classes for dealing with alignment/insertion positions ###########################

# has to be a new-style object-based class due to the immutability/hashability thing
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
            raise ValueError("The full_position argument must be a string of the form '100-200', '?-200' or '100-?'!")
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


def HTSeq_pos_to_tuple(HTSeq_pos):
    """ Convert an HTSeq.GenomicPosition instance to a (chrom,start_pos,end_pos,strand) tuple. 
    
    Start_pos and end_pos are 1-based, inclusive (so in AATTGG, the position of AA is 1-2) - unlike in HTSeq!
    """
    try:
        chrom = HTSeq_pos.chrom
    except AttributeError:
        raise MutantError("Invalid position %s! Need an HTSeq iv object. (If empty, maybe read wasn't aligned?)"%(HTSeq_pos,))
    ### cassette strand is the same as read strand, OR the opposite if reads_are_reverse is True
    strand = HTSeq_pos.strand
    if strand not in SEQ_STRANDS:   raise MutantError("Invalid strand %s!"%strand)
    # HTSeq is 0-based and I want 1-based, thus the +1; end has no +1 because in HTSeq end is the base AFTER the alignment.
    start_pos = HTSeq_pos.start+1
    end_pos = HTSeq_pos.end
    if start_pos < 1:           raise MutantError("Sequence positions must be positive!")
    if start_pos > end_pos:     raise MutantError("Sequence start can't be after end!")
    return (chrom, start_pos, end_pos, strand)


def get_insertion_pos_from_flanking_region_pos(flanking_region_pos, cassette_end, reads_are_reverse=False, immutable_position=True):
    """ Return a Insertion_position instance giving the cassette insertion position based on HTSeq read position. 

    Flanking_region_Pos should be a HTSeq.GenomicPosition instance, or a (chrom,start_pos,end_pos,strand) tuple
      (the tuple should have 1-based end-inclusive positions, so AA is 1-2 in AATT; HTSeq positions are 0-based end-exclusive); 
     cassette_end gives the side of the insertion that read is on; 
     reads_are_reverse is True if the read is in reverse orientation to the cassette, False otherwise. 

    The cassette chromosome will be the same as read chromosome; the cassette strand will be the same as read strand, 
     or opposite of the read strand if reads_are_reverse is True.

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
    try:                                chrom, start_pos, end_pos, strand = flanking_region_pos
    except (TypeError, ValueError):     chrom, start_pos, end_pos, strand = HTSeq_pos_to_tuple(flanking_region_pos) 
    if strand not in SEQ_STRANDS:   raise MutantError("Invalid strand %s!"%strand)
    if start_pos < 1:               raise MutantError("Flanking region positions must be positive!")
    if start_pos > end_pos:         raise MutantError("Flanking region start can't be after end!")
    ### chromosome is always the same as read, so just leave it as is
    ### cassette strand is the same as read strand, OR the opposite if reads_are_reverse is True
    if reads_are_reverse:     strand = '+' if strand=='-' else '-'
    ### cassette position depends on the read position and cassette_end in a somewhat complex way 
    #   (see description in docstring, and ../notes.txt for even more detail)
    if (cassette_end=='5prime' and strand=='+') or (cassette_end=='3prime' and strand=='-'):      
        pos_before, pos_after = end_pos, None
    elif (cassette_end=='3prime' and strand=='+') or (cassette_end=='5prime' and strand=='-'):    
        pos_before, pos_after = None, start_pos
    else:                           
        raise ValueError("cassette_end argument must be one of %s."%SEQ_ENDS)
    return Insertion_position(chrom, strand, position_before=pos_before, position_after=pos_after, 
                              immutable=immutable_position)


def find_gene_by_pos_gff3(insertion_pos, chromosome_GFF_record, detailed_features=False, quiet=False):
    """ Look up insertion_pos in chromosome_GFF_record; return (gene_ID,orientation,subfeature) of insertion in gene.

    Insertion_pos is an Insertion_position instance, or a (strand, start_pos, end_pos) tuple; 
     insertion_pos is assumed to be in the chromosome given by chromosome_GFF_record (the caller should check that).

    If insertion_pos overlaps a gene in chromosome_GFF_record, 
     return geneID (with '(edge)' appended if insertion_pos isn't completely contained inside the gene), 
      orientation ('sense' if insertion_pos and gene are in the same direction, 'antisense' otherwise), 
      and the name of the subfeature (exon/intron/UTR) Insertion_pos is in (or '?' if detailed_features is False); 
      if insertion_pos is not in a gene, return ('no_gene_found', '-', '-').
     If Insertion_pos is on the edge of two features, subfeature will be 'X/Y'; if it's something more 
      unlikely/complicated, subfeature will be 'X/Y/Z??' and a warning will be printed to STDOUT.

    Chromosome_GFF_record is a record generated by BCBio.GFF parser from a gff3 file (usually one file = multiple records). 
    """
    # get the needed information from either input format
    try:
        strand, ins_start, ins_end = insertion_pos.strand, insertion_pos.min_position, insertion_pos.max_position
    except AttributeError:
        strand, ins_start, ins_end = insertion_pos
    # see notes_on_GFF_parsing.txt for what a GFF3 record (chromosome_GFF_record) will be like
    assert strand in ['+','-','both'], "Strand should be %s, and is %s!"%(' or '.join(SEQ_STRANDS), strand)

    # make gene:data dictionary to store the results for each gene in!  In case an insertion does hit two genes.
    gene_data_list = []

    # go over all the genes in the chromosome record and look for any that match the position in insertion_pos
    for gene in chromosome_GFF_record.features:
        # for GFF positions, always add 1 to the gene/feature start, because BCBio uses 0-based and I use 1-based, 
        #  but don't add 1 to the end, because BCBio uses end-exclusive and I use end-inclusive.
        gene_start, gene_end = gene.location.start.position+1, gene.location.end.position
        if position_test_overlap(gene_start, gene_end, ins_start,ins_end):
            # if insertion_pos is inside the gene, save the gene ID
            gene_ID = gene.id
            # if it overlaps an edge, note that in features_gene by adding 'gene_edge'
            # (MAYBE-TODO if I look at things like inner/outer flanking regions, this is where that would go as well)
            if position_test_contains(gene_start, gene_end, ins_start,ins_end): features_gene = []
            else:                                                               features_gene = ['gene_edge']
            # calculate orientation of insertion_pos vs gene
            if strand=='both':      orientation = 'both'
            elif gene.strand==1:    orientation = 'sense' if strand=='+' else 'antisense'
            elif gene.strand==-1:   orientation = 'sense' if strand=='-' else 'antisense'
            else:                   orientation = '?'
            # if we're not looking for detailed features, just return the data now, using '?' as gene_feature
            if not detailed_features:
                full_feature = '/'.join(features_gene+['?'])
                return gene_ID, orientation, full_feature

            ### Find which feature of the gene the insertion_pos should be annotated as:
            # if gene has no features listed, use 'no_mRNA' as gene_feature (different from '?' or '-')
            if len(gene.sub_features)==0:
                inner_feature = 'no_mRNA'
            # if multiple mRNAs, just return 'MULTIPLE_SPLICE_VARIANTS' as feature instead of delving into it
            elif len(gene.sub_features)>1:
                inner_feature = 'MULTIPLE_SPLICE_VARIANTS'
                # TODO implement going over all the splice variants of the gene (even ones the insertion doesn't touch - important
                #  to know if all or some are affected!) and giving the feature info for all of them!  Comma-separated list.
            elif gene.sub_features[0].type != 'mRNA':
                if not quiet:
                    print("Warning: gene %s in gff file has unexpected non-mRNA sub-features! "%gene_ID
                          +"Returning '??' feature.")
                inner_feature = '??'
            else:
                mRNA = gene.sub_features[0]
                mRNA_start, mRNA_end = mRNA.location.start.position+1,mRNA.location.end.position
                # if insertion_pos is outside the mRNA, use 'outside_mRNA' as inner_feature
                if not position_test_overlap(mRNA_start, mRNA_end, ins_start, ins_end):
                    inner_feature = 'outside_mRNA'
                # if insertion_pos is inside the mRNA and mRNA has no subfeatures, use 'mRNA_no_exons' as inner_feature
                elif len(mRNA.sub_features)==0:   
                    if position_test_contains(mRNA_start, mRNA_end, ins_start, ins_end):
                        inner_feature = 'mRNA_no_exons'
                    else:
                        inner_feature = 'mRNA_edge'
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
                    # MAYBE-TODO may want to treat exons before 5'UTR or after 3'UTR specially? (EH, none in current file)
                    # MAYBE-TODO may want to treat cases with multiple UTRs specially?  We do have those in current file!
                    # if insertion_pos is inside a single mRNA subfeature, use the type of the subfeature as inner_feature
                    if len(features_inside)==1 and len(features_edge)==0:
                        inner_feature = features_inside[0]
                    # if insertion_pos is on the edge of two mRNA subfeatures, use 'subfeature1/subfeature2'
                    elif len(features_inside)==0 and len(features_edge)==2:
                        inner_feature = '/'.join(features_edge)
                        # MAYBE-TODO treat insertions CLOSE to an edge specially too? How large is a splice junction?
                    # if insertion_pos is on the edge of ONE mRNA subfeature, or not touching any subfeatures at all, 
                    #  the implied subfeature is either an intron (if between features) or mRNA_before/after_exons, 
                    #   (which shouldn't happen in normal files).
                    elif len(features_inside)==0 and len(features_edge)<=1:
                        # figure out what the implied feature is - outside intron in CDS (normal) or either UTR, or outside all exons
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
                            inner_feature = features_edge[0] + '/' + implied_feature
                        elif len(features_edge)==0:
                            inner_feature = implied_feature
                    # if insertion_pos is inside two features, or inside one and on the edge of another, 
                    #  print a warning, and use all the feature names, with a ?? at the end to mark strangeness
                    else:
                        inner_feature = '/'.join(features_inside+features_edge) + '??'
                        if not quiet:
                            print(("Warning: Location (%s,%s) matched multiple features (%s) "
                                  +"in gene %s!")%(ins_start, ins_end, inner_feature, gene_ID)) 
                # MAYBE-TODO also output distance from start/end of gene/feature/CDS?
            
            # prepend whatever gene-level features (edge etc, or []) were found at the start to the full value
            full_feature = '/'.join(features_gene+[inner_feature])
            gene_data_list.append([gene_ID, orientation, full_feature])
    # if no gene matching insertion_pos was found, return special value
    if not gene_data_list:
        return [SPECIAL_GENE_CODES.not_found, '-', '-']
    # if single gene found, return its info
    elif len(gene_data_list) == 1:
        return gene_data_list[0]
    # if multiple genes found, return data like this: "gene1 | gene2", "sense | antisense", "intron | CDS/4'UTR" and print warning
    else: 
        if not quiet:
            print("Warning: Location \"%s\" matched multiple genes! %s"%(insertion_pos, ', '.join(zip(*gene_data_list)[0])))
        return [MULTIPLE_GENE_JOIN.join(multiple_vals) for multiple_vals in zip(*gene_data_list)]
    # MAYBE-TODO add unit tests?  But this is included in a pretty thorough run-test, so may not be necessary.
    # TODO add unit-test or run-test for overlapping genes!!


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


############################### Main classes describing the mutants and mutant sets #####################################

# help functions returning "blank" mutants for defaultdicts (can't use lambdas because pickle doesn't like them)
def blank_readcount_only_mutant():
    return Insertional_mutant_readcount_only()
def blank_full_single_mutant():
    return Insertional_mutant()
def blank_single_mutant_with_pos(pos):
    return Insertional_mutant(insertion_position=pos)
def blank_multi_mutant_with_pos(pos):
    return Insertional_mutant_multi_dataset(insertion_position=pos)


class Insertional_mutant():
    """ Data regarding a particular insertional mutant: insertion position, gene, read numbers/sequences, etc.

    Mutants have the following attributes (data):
     1) Position/gene attributes (except readcount-related-only mutants, which don't have those):
       - position - an Insertion_position instance giving the insertion chromosome/strand/position
       - gene, orientation, gene feature - what gene the insertion is in (or one of the SPECIAL_GENE_CODES if unknown), 
                   whether it's in the sense or antisense orientation vs the gene, what feature (exon/intron/UTR) it's in.
     2) Readcount-related attributes (multi-dataset mutants don't have those on the top-level):
       - total_read_count, perfect_read_count - number of all and perfectly aligned deepseq reads
       - unique_sequence_count, sequences_and_counts - number of unique read sequences, and a seq:count dictionary

    This class also has subclasses for mutant types with additional functionality or missing functionality.

    Some notes on methods (functions) of Insertional_mutant objects:
     For detailed information on mutant methods, see method docstrings.
     Methods with names starting with _ are private and shouldn't be used from outside the object itself.
     Readcount-related methods take a dataset_name argument for compatibility with the multi-dataset subclass
      - on basic mutants it should always be None (which is the default).
    """
    # TODO add a __slots__ to this to optimize memory usage for when there's tens of thousands of mutants!  But that may cause issues with pickling, etc...

    def __init__(self, insertion_position=None):
        """ Set self.position based on argument; initialize read/sequence counts to 0 and gene-info to unknown. 

        insertion_position argument should be an Insertion_position instance. 
        """
        # "standard" mutants need general attributes and readcount-related attributes
        self._set_general_attributes(insertion_position)
        self._set_readcount_related_data_to_zero()

    def read_info(self, dataset_name=None, strict=False):
        """ Help function to get read-info-containing object for both multi-dataset and single mutants.

        Here this just returns self; more complicated version is for multi-dataset mutants.
        Strict is ignored and only present to make the implementation consistent with the multi-dataset version.
        """
        if dataset_name is not None:
            raise MutantError("This is not a multi-dataset mutant - cannot provide dataset_name arg!")
        else:
            return self

    # MAYBE-TODO give each mutant some kind of unique ID at some point in the process?  Or is genomic location sufficient?  If we end up using per-mutant barcodes (in addition to the flanking sequences), we could use that, probably, or that plus genomic location.

    def _set_general_attributes(self, insertion_position):
        self.position = insertion_position
        # MAYBE-TODO should I have a class for the gene data? Especially if I add more of it (annotation info etc)
        self.gene = SPECIAL_GENE_CODES.not_determined
        self.orientation = '?'
        self.gene_feature = '?'

    def _set_readcount_related_data_to_zero(self):
        """ Set all readcount-related data to 0/empty."""
        self.total_read_count      = 0
        self.perfect_read_count    = 0
        self.unique_sequence_count = 0
        self.sequences_and_counts  = defaultdict(int)
        self.original_strand_readcounts = {}

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

    def add_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False, dataset_name=None):
        """ Add a read to the data (or multiple identical reads, if read_count>1); return True if perfect match.

        Specifically: increment total_read_count, increment perfect_read_count if read is a perfect 
         alignment, increment the appropriate field of sequences_and_counts based on read sequence.
        Note: this does NOT check the read position to make sure it matches that of the object.

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        self._ensure_dataset_None(dataset_name)
        # MAYBE-TODO check HTSeq_alignment chromosome/strand to make sure it matches data in self?  Don't check position, that's more complicated (it can be either start or end) - could maybe check that position is within, idk, 10bp of either alignment start or alignment end?  Or not - I may want to cluster things in a non-position-based way anyway!  Hmmm...
        seq = HTSeq_alignment.read.seq
        # if it's a new sequence, increment unique_sequence_count; add a count to the sequences_and_counts dictionary.
        if seq not in self.sequences_and_counts:
            self.unique_sequence_count += 1
        self.sequences_and_counts[seq] += read_count
        # increment total_read_count
        self.total_read_count += read_count
        # figure out if the read is perfect and increment perfect_read_count if yes; return True if perfect else False.
        treat_unknown_as = 'match' if treat_unknown_as_match else 'mutation'
        mutation_count = check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as=treat_unknown_as)
        if mutation_count==0:  
            self.perfect_read_count += read_count
            return True
        else:
            return False

    def update_gene_info(self, gene, orientation, gene_feature):
        """ Update gene/orientation/feature: if both are the same or one is unknown, keep known; if different, raise error.
        """
        # grab the set of non-special values from own and new data
        gene_both = set([self.gene, gene]) - set([SPECIAL_GENE_CODES.not_determined])
        orientation_both = set([self.orientation, orientation]) - set(['?'])
        feature_both = set([self.gene_feature, gene_feature]) - set(['?'])
        # if there are two different non-special values for any of the three attributes, raise MutantError
        if len(gene_both)>1 or len(orientation_both)>1 or len(feature_both)>1:
            raise MutantError("Can't merge the two mutants: the gene/orientation/feature data differs!")
        # otherwise set own data to the better one of own/new data (unless neither are present)
        if gene_both:           self.gene = gene_both.pop()
        if orientation_both:    self.orientation = orientation_both.pop()
        if feature_both:        self.gene_feature = feature_both.pop()

    def merge_mutant(self, other, check_gene_data=True, opposite_strand_tandem=False, other_is_imperfect=False):
        """ Merge other mutant into this mutant: merge counts, sequences, etc; set other's counts to 0.
        Does NOT check that the positions match - this should be done by the caller. 
        Does NOT merge the positions, except for strand in the opposite_strand_tandem==True case.

        If other_is_imperfect, do NOT count the other mutant's reads as perfect (if it's an off-by-one mutant due to indel, or such).
        If opposite_strand_tandem, expect the two mutants to be +/- strand; set strand and orientation to 'both', 
         and store original +/- strand readcounts in self.original_strand_readcounts.
        """
        if isinstance(other, (Insertional_mutant_multi_dataset, Insertional_mutant_readcount_only)):  
            raise MutantError("Can't merge single mutant with a multi-dataset or readcount-only one!")

        # make sure the two mutants don't have conflicting gene data, if required  (and update if own gene data is unknown)
        # (note: NOT checking position) (this needs to be done first, BEFORE we make changes to the mutants!)
        if check_gene_data:
            try:
                # in the opposite_strand_tandem case, we expect orientations to be different, so pass '?' to the check
                if opposite_strand_tandem:  self.update_gene_info(other.gene, '?', other.gene_feature)
                else:                       self.update_gene_info(other.gene, other.orientation, other.gene_feature)
            except MutantError:
                raise MutantError("Can't merge the two mutants: the gene/orientation/feature data differs!")
        
        ### merge positions
        # in the opposite_strand_tandem case, the new strand should be "both" and the new position should give both ends
        if opposite_strand_tandem:
            assert self.position.strand!=other.position.strand
            # keep track of original by-strand readcounts
            self.original_strand_readcounts = dict([(m.position.strand, m.total_read_count) for m in (self, other)])
            self.position.strand = 'both'
            self.position.position_before = self.position.position_before or other.position.position_before
            self.position.position_after = self.position.position_after or other.position.position_after
            # also set the gene-orientation to "both", unless there is no gene
            if self.gene not in SPECIAL_GENE_CODES.all_codes:
                self.orientation = 'both'
        # otherwise, just keep the self position.
        # MAYBE-TODO implement position merging?  Are we ever going to need it, really?  PROBABLY NOT.

        ### merge read counts
        self.total_read_count += other.total_read_count
        # whether the other mutant's reads should be counted as perfect depends on the situation - if we're adding together mutants
        #  from two datasets, or opposite-tandems, both are equally good, but if we're merging same-strand-adjacent mutants 
        #  (where one has a position off by one due to indel), the other mutant shouldn't be considered perfect. 
        if not other_is_imperfect:
            self.perfect_read_count += other.perfect_read_count
        # merge sequences
        for (seq,count) in other.sequences_and_counts.iteritems():
            self.sequences_and_counts[seq] += count
        self.unique_sequence_count = len(self.sequences_and_counts)
        # LATER-TODO may want to keep more data about sequences! Like exact position and strand and number of mutations - may want to store a list of HTSeq.alignment objects instead of just sequences+counts, really.
        # TODO keep full info about the original mutants, somehow?  I don't know if that should be kept in the mutant, or separately as a list of merged mutants in the dataset.
        other._set_readcount_related_data_to_zero()

    def add_counts(self, total_count, perfect_count, sequence_variant_count, assume_new_sequences=False, dataset_name=None):
        """ Increment self.total_read_count, self.perfect_read_count and self.unique_sequence_count based on inputs.

        Note that if self.unique_sequence_count>0, it's impossible to determine the correct new value: 
         if we had old data with one unique sequence and now we have new data with another one, how do we know
          if that's the same or different sequence?  The correct total could be 1 or 2, so it's an option:
         If assume_new_sequences is True, the total is old+new; if it's False, the total is max(old,new).

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        self._ensure_dataset_None(dataset_name)
        self.total_read_count += total_count
        self.perfect_read_count += perfect_count
        if assume_new_sequences:
            self.unique_sequence_count += sequence_variant_count
        else:
            self.unique_sequence_count = max(self.unique_sequence_count, sequence_variant_count)

    def add_sequence_and_counts(self, seq, seq_count, add_to_uniqseqcount=True, dataset_name=None):
        """ Add seq_count to self.sequences_and_counts[seq] (it's created with count 0 if seq wasn't a key before).

        Note: if add_to_uniqseqcount is False, this will never increment self.unique_sequence_count;
         otherwise it only does so if seq was not already present in the self.sequences_and_counts data
          and if the total number of sequences in self.sequences_and_counts is higher than self.unique_sequence_count. 

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        self._ensure_dataset_None(dataset_name)
        # increment unique_sequence_count if desired (needs to be done first because of the checks it's doing)
        if add_to_uniqseqcount:
            if seq not in self.sequences_and_counts and\
               len(self.sequences_and_counts)>self.unique_sequence_count:
                    self.unique_sequence_count += 1
        # main function: add another seq_count counts of seq
        self.sequences_and_counts[seq] += seq_count

    @staticmethod
    def get_main_sequence_from_data(seqs_to_counts, N=1):
        """ Return the most common sequence in the given data, and its count (or Nth most common sequence if N is provided). """
        sequences_by_count = sorted([(count,seq) for (seq,count) in seqs_to_counts.iteritems()], reverse=True)
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return tuple(reversed(sequences_by_count[N-1]))
        except IndexError:  return ('',0)
        # MAYBE-TODO should probably make that '-' or something instead of '', empty strings are hard to see. 
        #  On the other hand '-' isn't a valid sequence, and '' is...

    def get_main_sequence(self, N=1, dataset_name=None):
        """ Return the most common sequence in this mutant and its count (or Nth most common sequence if N is provided).

        Dataset_name should NEVER be set to non-None, and only present for consistency with the multi-dataset subclass.
        """
        self._ensure_dataset_None(dataset_name)
        return self.get_main_sequence_from_data(self.sequences_and_counts, N)

    def _copy_non_readcount_data(self, source_mutant):
        """ Copy non-readcount-related data from source_mutant to self (making new copies of all objects). """
        # COPY the position, not just make another name for the same value - I wrote a copy() function for positions
        self.position     = source_mutant.position.copy() 
        # strings are immutable and thus safe to "copy" by adding another name to the same value
        self.gene         = source_mutant.gene
        self.orientation  = source_mutant.orientation
        self.gene_feature = source_mutant.gene_feature

    def _copy_readcount_related_data(self, source_mutant):
        """ Copy readcount-related data from source_mutant to self (making new copies of all objects). """
        # integers are immutable and thus safe to "copy" by adding another name to the same value
        self.total_read_count      = source_mutant.total_read_count
        self.perfect_read_count    = source_mutant.perfect_read_count
        self.unique_sequence_count = source_mutant.unique_sequence_count
        # using dict to make a COPY of the dict instead of just creating another name for the same value
        self.sequences_and_counts  = dict(source_mutant.sequences_and_counts)
        # the try/except is for dealing with merging older datasets that didn't have self.original_strand_readcounts
        try:
            self.original_strand_readcounts  = dict(source_mutant.original_strand_readcounts)
        except AttributeError:
            self.original_strand_readcounts = {}

    # MAYBE-TODO should there be a copy_mutant function to make a deepcopy?


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

    def merge_mutant(self, *args, **kwargs):
        raise MutantError("merge_mutant NOT IMPLEMENTED on readcount-related-only mutant object!")
        # MAYBE-TODO implement merging for multi-dataset mutants?  Seems like unnecessary complication for now.

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

    def __init__(self, insertion_position=None):
        # multi-dataset mutants get general attributes
        Insertional_mutant._set_general_attributes(self, insertion_position)
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

    def add_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False, dataset_name=None):
        """ Add read to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_read(readcount_data_container, HTSeq_alignment, read_count, treat_unknown_as_match)
    
    def merge_mutant(self, *args, **kwargs):
        raise MutantError("merge_mutant NOT IMPLEMENTED on multi-dataset mutant object!")
        # MAYBE-TODO implement merging for multi-dataset mutants?  Seems like unnecessary complication for now.

    def add_counts(self, total_count, perfect_count, sequence_variant_count, assume_new_sequences=False, dataset_name=None):
        """ Add counts to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_counts(readcount_data_container, total_count, perfect_count, sequence_variant_count, 
                                      assume_new_sequences, dataset_name=None)

    def add_sequence_and_counts(self, seq, seq_count, add_to_uniqseqcount=True, dataset_name=None):
        """ Add seqs/counts to given dataset (see docstring for Insertional_mutant version) - dataset_name is required. """
        readcount_data_container = self._check_dataset_name_return_data(dataset_name)
        Insertional_mutant.add_sequence_and_counts(readcount_data_container, seq, seq_count, add_to_uniqseqcount, dataset_name=None)
    
    def get_main_sequence(self, N=1, dataset_name=None):
        """ Return the most common sequence in this mutant and its count (or Nth most common sequence if N is provided).

        If dataset_name is given, return the most common sequence for just that dataset; 
         or if dataset_name is None, return most common sequence by total count over all the datasets.
        """
        if dataset_name is not None:
            seqs_to_counts = self.by_dataset[dataset_name].sequences_and_counts
        else:
            seqs_to_counts = reduce(add_dicts_of_ints, 
                                    [dataset.sequences_and_counts for dataset in self.by_dataset.itervalues()])
        # MAYBE-TODO print a warning if different dataset mutants have different main sequences?
        return Insertional_mutant.get_main_sequence_from_data(seqs_to_counts, N)
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
                self.update_gene_info(other_mutant.gene, other_mutant.orientation, other_mutant.gene_feature)
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


class Dataset_summary_data():
    """ Summary data for a Insertional_mutant_pool_dataset object.  Lots of obvious attributes; no non-default methods.
    """
    # LATER-TODO update docstring

    def __init__(self, dataset, cassette_end, reads_are_reverse, dataset_name=None):
        """ Initialize everything to 0/empty/unknown. """
         # make sure the arguments are valid values
        if not cassette_end in SEQ_ENDS+['?']: 
            raise ValueError("The cassette_end variable must be one of %s or '?'!"%SEQ_ENDS)
        if not reads_are_reverse in [True,False,'?']: 
            raise ValueError("The reads_are_reverse variable must be True, False, or '?'!")
        # reference to the containing dataset (for read-counting purposes etc), 
        #  and the dataset name (None if it's a single dataset, string for multi-datasets)
        self.dataset_name = dataset_name
        self.dataset = dataset
        # information on reads that aren't included in the dataset mutants - None or 0 by default
        # TODO I should really go over this and figure out what should be None and what should be 0 and why!!
        self.discarded_read_count, self.discarded_wrong_start, self.discarded_no_cassette = None, None, None
        self.discarded_other_end = 0
        self.non_aligned_read_count, self.unaligned, self.multiple_aligned = 0, None, None
        self.ignored_region_read_counts = defaultdict(int)
        # mutant merging information
        self.blank_adjacent_mutant_info()
        # MAYBE-TODO should cassette_end and reads_are_reverse be specified for the whole dataset, or just for each set of data added, in add_alignment_reader_to_data? The only real issue with this would be that then I wouldn't be able to print this information in the summary - or I'd have to keep track of what the value was for each alignment reader added and print that in the summary if it's a single value, or 'varied' if it's different values. Might also want to keep track of how many alignment readers were involved, and print THAT in the summary!  Or even print each (infile_name, cassette_end, reads_are_reverse) tuple as a separate line in the header.
        self.cassette_end = cassette_end
        self.reads_are_reverse = reads_are_reverse

    # TODO unit-test all the methods below!

    def blank_adjacent_mutant_info(self):
        # general adjacent counting/merging settings
        self.adjacent_max_distance = None
        self.merging_which_chromosomes = (None, None)
        # info on merged mutants
        self.merged_adjacent_same_strand_dict = defaultdict(int)
        self.merged_adjacent_same_strand_readcounts_dict = defaultdict(list)
        self.merged_opposite_tandems = 0
        self.merged_opposite_tandems_readcounts = []
        # info on unmerged close-by mutants
        # MAYBE-TODO should the adjacent-mutant-counts (same_position_opposite, adjacent_opposite_toward, adjacent_opposite_away, adjacent_same_strand) be constants or property values?  May want to keep them as constants, they're expensive to calculate...
        self.same_position_opposite = NaN
        self.adjacent_opposite_toward_dict = defaultdict(nan_func)
        self.adjacent_opposite_away_dict = defaultdict(nan_func)
        self.adjacent_same_strand_dict = defaultdict(nan_func)

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
            try:
                strand_dict[m.position.strand] += m.read_info(self.dataset_name).total_read_count
            except KeyError:
                # special case for both-strand mutants - they have a original_strand_readcounts attribute
                if m.position.strand == 'both':
                    for strand,count in m.read_info(self.dataset_name).original_strand_readcounts.items():
                        strand_dict[strand] += count
                else:
                    raise MutantError("Unknown strand %s! %s"%m.position.strand)
        return strand_dict

    def reads_in_chromosome(self, chromosome):
        """ Return total number of reads in given chromosome."""
        return sum(m.read_info(self.dataset_name).total_read_count 
                   for m in self.dataset if m.position.chromosome==chromosome)

    @property
    def all_chromosomes(self):
        return set(m.position.chromosome for m in self.dataset if m.read_info(self.dataset_name).total_read_count)
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
                   and m.position.chromosome==chromosome)

    def merged_gene_feature_counts(self, merge_boundary_features=True, merge_confusing_features=False):
        """ Return (gene_feature,count) list, biologically sorted, optionally with all "boundary" features counted as one.

        The source gene feature counts are based on the self.mutant_counts_by_feature dict.
        If merge_confusing_features==True, any locations containing '??' will be listed as '??'.
        If merge_boundary_features==True, any locations containing '/' and no '??' will be listed as 'boundary'.
        The custom sort order (based on what seems sensible biologically) is: CDS, intron, UTR, other, boundary.
        """
        merged_feature_count_dict = defaultdict(int)
        for feature, count in self.mutant_counts_by_feature.items():
            # note that anything containing '??' AND '/' never gets merged as boundary
            if '??' in feature:
                if merge_confusing_features:                  merged_feature_count_dict['??'] += count
                else:                                         merged_feature_count_dict[feature] += count
            elif '/' in feature and merge_boundary_features:  merged_feature_count_dict['boundary'] += count
            else:                                             merged_feature_count_dict[feature] += count
        return merged_feature_count_dict

    @property
    def most_common_mutants(self):
        """ Return list of mutants with the most total reads (in dataset if multi-dataset)."""
        highest_readcount = max([mutant.read_info(self.dataset_name).total_read_count for mutant in self.dataset])
        highest_readcount_mutants = [mutant for mutant in self.dataset 
                                     if mutant.read_info(self.dataset_name).total_read_count==highest_readcount]
        return highest_readcount_mutants

    @property
    def adjacent_opposite_toward(self):
        max_distance = self.adjacent_max_distance or 1
        return sum([self.adjacent_opposite_toward_dict[i] for i in range(1, max_distance+1)])
    @property
    def adjacent_opposite_away(self):
        max_distance = self.adjacent_max_distance or 1
        return sum([self.adjacent_opposite_away_dict[i] for i in range(1, max_distance+1)])
    @property
    def adjacent_same_strand(self):
        max_distance = self.adjacent_max_distance or 1
        return sum([self.adjacent_same_strand_dict[i] for i in range(1, max_distance+1)])
    @property
    def merged_adjacent_pairs(self):
        return sum(self.merged_adjacent_same_strand_dict.values())

    def adjacent_mutant_summary(self, long_version=False, add_total_counts=True, max_distance=None):
        """ Return string giving adjacent-mutant counts by category. """
        # To be able to use the adjacent_same_strand etc properties with custom max_distance, 
        #  set self.adjacent_max_distance to the custom value temporarily, then reset afterward
        if max_distance is not None:
            orig_max_distance = self.adjacent_max_distance
            self.adjacent_max_distance = max_distance
        if long_version:
            string = "Adjacent mutant counts (max distance %s): %s adjacent same-strand pairs, %s same-position opposite-strand pairs, %s adjacent opposite-strand away-facing pairs (may be tandems with a deletion), %s adjacent opposite-strand toward-facing pairs (definitely two separate mutants)."
            if add_total_counts:
                string = "Dataset has %s mutants (with %s total reads). "%(len(self.dataset), 
                                                                           self.aligned_read_count) + string
        else:
            string = "(max dist %s) %s adjacent same-strand pairs, %s opposite-tandem, %s adjacent-opposite away-facing, %s toward-facing"
            if add_total_counts:
                string = "%s mutants (%s reads); "%(len(self.dataset), self.aligned_read_count) + string
        string = string%(self.adjacent_max_distance, self.adjacent_same_strand, 
                       self.same_position_opposite, self.adjacent_opposite_away, self.adjacent_opposite_toward)
        if max_distance is not None:
            self.adjacent_max_distance = orig_max_distance 
        return string

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
     - cassette_end - specifies which end of the insertion cassette the reads are on, and 
     - reads_are_reverse - True if the reads are in reverse orientation to the cassette, False otherwise
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

    def __init__(self, cassette_end='?', reads_are_reverse='?', multi_dataset=False, infile=None):
        """ Initializes empty dataset; saves properties as provided; optionally reads data from mutant infile. """
        # _mutants_by_position is the main data structure here, but it's also private:
        #      see "METHODS FOR EMULATING A CONTAINER TYPE" section below for proper ways to interact with mutants.
        #   the single mutants should be single-dataset by default, or multi-dataset if the containing dataset is.
        if multi_dataset:   blank_mutant_function = blank_multi_mutant_with_pos
        else:               blank_mutant_function = blank_single_mutant_with_pos
        self._mutants_by_position = keybased_defaultdict(blank_mutant_function)
        # various dataset summary data - single for a single dataset, a dictionary for a multi-dataset object.
        self.multi_dataset = multi_dataset
        if not multi_dataset:   self.summary = Dataset_summary_data(self, cassette_end, reads_are_reverse, None)
        else:                   self.summary = {}
        # data that's NOT related to a particular dataset
        # gene/annotation-related information - LATER-TODO should this even be here, or somewhere else?
        self.gene_annotation_header = []
        self.total_genes_in_genome = 0
        # optionally read mutant data from infile (TODO maybe should remove this, it's supposedly deprecated...)
        if infile is not None:
            self.read_data_from_file(infile)

    ######### METHODS FOR EMULATING A CONTAINER TYPE (a dataset is essentially a container of mutants)
    # everything is currently based on the private self._mutants_by_position dictionary, but this may change.
    # compared to just using a public self.mutants_by_position dictionary, this approach is better because:
    #   - it's more flexible
    #   - it doesn't allow having self.mutants_by_position[pos1] = mutant(pos2), with mismatched positions
    # Note that all of this will work just fine for multi-datasets without changes, since multi-datasets simply contain
    #  the same dictionary but with multi-dataset mutants instead of single-dataset ones.

    def __len__(self):      
        # Return the number of mutants with some non-zero readcounts (DON'T count zero-read mutants, even if present!).
        #  A bit tricky for multi-datasets. 
        if not self.multi_dataset: 
            return self.summary.N_mutants
        else:
            return len([m for m in self if sum([m.by_dataset[x].total_read_count for x in m.by_dataset.keys()])])

    def __iter__(self):     return self._mutants_by_position.itervalues()
    # for multi-datasets, for iterating over only mutants with 1+ reads in a particular dataset, see mutants_in_dataset()

    @property
    def size(self):         return len(self)

    # instead of __setitem__
    def add_mutant(self, mutant, overwrite=False):
        """ Add mutant to dataset. """
        if mutant.position in self._mutants_by_position.keys() and not overwrite:
            raise MutantError("Can't add mutant that would overwrite previous mutant at same position! "
                              +"Pass overwrite=True argument if you want to overwrite.")
        self._mutants_by_position[mutant.position] = mutant

    # instead of __delitem__
    def remove_mutant(self, *args, **kwargs):
        """ Remove mutant (by position) - can take a mutant, a position, or arguments to create a position. """
        if len(args)==1 and not kwargs and isinstance(args[0], Insertional_mutant):
            position = args[0].position
        else:
            position = self._make_position_object_if_needed(*args, **kwargs)
        del self._mutants_by_position[position]

    # instead of __getitem__
    def get_mutant(self, *args, **kwargs):
        """ Return the mutant with given position (given an Insertion_position instance, or arguments to create one).
        If mutant doesn't exist, create a new one with no reads/sequences. """
        return self._mutants_by_position[self._make_position_object_if_needed(*args,**kwargs)]

    def get_matching_position_mutants(self, *args, **kwargs):
        """ Return a list of mutants "matching" given position - same-position both-strand and +/-strand positions count as matching.
        Input should be an Insertion_position instance, or arguments to create one.

        Currently NOT QUITE FULLY IMPLEMENTED OR TESTED - see TODO comments in the code.
        """
        # TODO in addition to strand, this should return 100-101 when given ?-101 and 100-?, probably?  And vice versa?  TODO think about that!  Currently those are ALL CONSIDERED DISTINCT. 
        # TODO should this be a method of Insertion_position instead?
        position = self._make_position_object_if_needed(*args,**kwargs)
        chromosome = position.chromosome
        pos_before, pos_after = position.position_before, position.position_after
        matching_positions = [position]
        # if strand is both, look for all possible matching single-strand mutants (with known position before, after, or both)
        if position.strand=='both':
            for strand in '+-':
                for (before, after) in [(pos_before, None), (None, pos_after), (pos_before, pos_after)]:
                    matching_positions.append(Insertion_position(chromosome, strand, position_before=before, position_after=after))
        # if strand is +/-, look for a both-strand mutant, based on whichever one of positions before/after is defined
        else:
            if pos_before is None:      matching_positions.append(Insertion_position(chromosome, 'both', 
                                                             position_before=pos_after-1, position_after=pos_after))
            elif pos_before is None:    matching_positions.append(Insertion_position(chromosome, 'both', 
                                                             position_before=pos_before, position_after=pos_before+1))
            else:                       matching_positions.append(Insertion_position(chromosome, 'both', 
                                                             position_before=pos_before, position_after=pos_after))
            # TODO what if our mutant is +strand 100-?, and there's an existing both-strand mutant that's 100-105 or something, NOT the expected 100-101?
        return matching_positions

    # TODO how should mutant lookup by position deal with both-strand mutants?  You can have two mutants in the same position on opposite strands, but you CAN'T have one on both strands and one on + or -...  There should be safeguards to fold the + or - into the both-strand one if it's searched for, and things! THINK ABOUT THAT.

    # MAYBE-TODO implement get_Nth_mutant_by_position and get_Nth_mutant_by_readcount or some such?

    def __contains__(self, *args, **kwargs):
        """ Check if dataset contains mutant with given position (Insertion_position instance, or arguments to create one).

        If there's one (non-keyword) argument and it's an Insertion_position instance, use that for the check, 
         otherwise create a new Insertion_position instance with args/kwargs and check for that.

        You can only check "position in dataset" at present, not "mutant in dataset", since the latter would probably
         also just check by position, so the syntax would be misleading.
        """
        return self._make_position_object_if_needed(*args,**kwargs) in self._mutants_by_position
    
    @staticmethod
    def _make_position_object_if_needed(*args, **kwargs):
        """ Given either an Insertion_position instance or arguments to create one, return Insertion_position instance."""
        # Check if there's a single argument that already is an insertion position
        #  (I used to check for isinstance here instead of hasattr, but that seems to break with pickling?...)
        if len(args)==1 and not kwargs and all([hasattr(args[0], attr) for attr in 'chromosome strand full_position'.split()]):
            return args[0]
        else:
            if 'immutable' not in kwargs:
                kwargs['immutable'] = True
            return Insertion_position(*args,**kwargs)


    ######### READING BASIC DATA INTO DATASET

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, uncollapse_read_counts=False, 
                                     ignore_cassette=False, cassette_only=False, treat_unknown_as_match=False):
        """ Adds all alignments from the reader to the mutant data; currently based only on position, but that may change. 

        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader).

        Set uncollapse_read_counts to True if the original deepseq data was collapsed to unique sequences using
         fastx_uncollapser before alignment, to get the correct original read counts.

        Treat_unknown_as_match governs whether alignments with no detailed information are treated as perfect or not.

        Different ways of treating cassette reads:
         - by default (if all *cassette* args are false) treat cassette reads normally
         - if ignore_cassette==True, ignore cassette reads in the data and list them as removed in the header
         - if cassette_only==True, ignore all OTHER reads and only include cassette reads!
        """
        # LATER-TODO actually instead of cassette_only it might be good to just generate two separate mutant-sets, normal and cassette, with an option called separate_cassette or something, and print them to separate files - but that's more complicated, and right now I don't have the setup for a single dataset having multiple mutant-sets (although I guess I will have to eventually, for removed mutants etc). Right now I do it in mutant_count_alignments.py, which works but there's a lot of code repetition...
        if self.multi_dataset:  raise MutantError("add_alignment_reader_to_data not implemented for multi-datasets!")
        if ignore_cassette and cassette_only:
            raise MutantError("Only one of ignore_cassette and cassette_only arguments can be True - mutually exclusive!")

        summ = self.summary
        if summ.cassette_end == '?':
            raise MutantError("Cannot add data from an alignment reader if cassette_end isn't specified! Please set the "
          +"summary.cassette_end attribute of this Insertional_mutant_pool_dataset instance to one of %s first."%SEQ_ENDS)
        if summ.reads_are_reverse == '?':
            raise MutantError("Cannot add data from an alignment reader if reads_are_reverse isn't set! Please set the "
                  +"reads_are_reverse attribute of this Insertional_mutant_pool_dataset instance to True/False first.")

        # Go over all the reads in the HTSeq_alignment_reader, add them to dataset
        for aln in HTSeq_alignment_reader:
            if uncollapse_read_counts:      read_count = get_seq_count_from_collapsed_header(aln.read.name)
            else:                           read_count = 1
            # if read is unaligned, add to unaligned count and skip to the next read
            if (not aln.aligned) or (aln.iv is None):
                summ.non_aligned_read_count += read_count
                continue
            # get the cassette insertion position (as an Insertion_position object) - USE IMMUTABLE POSITIONS BY DEFAULT
            position = get_insertion_pos_from_flanking_region_pos(aln.iv, summ.cassette_end, summ.reads_are_reverse, 
                                                                  immutable_position=True)
            # if read is aligned to cassette and should be ignored, add to the right count and skip to the next read
            if ignore_cassette and is_cassette_chromosome(position.chromosome):
                summ.ignored_region_read_counts[position.chromosome] += read_count
                continue
            # or if we only want cassette reads, skip non-cassette ones!
            elif cassette_only and not is_cassette_chromosome(position.chromosome):
                summ.ignored_region_read_counts['NON-CASSETTE'] += read_count
                continue
            # grab the right mutant based on the position, and add the reads to it; 
            curr_mutant = self.get_mutant(position)
            curr_mutant.add_read(aln, read_count, treat_unknown_as_match=treat_unknown_as_match)
        # special case for when we don't know the specific unaligned categories, but we know total non-aligned is 0, 
        #  so the specific categories must be 0 too:
        if summ.non_aligned_read_count==0:  summ.unaligned, summ.multiple_aligned = 0, 0

    def read_data_from_file(self, infile, assume_new_sequences=False):
        """ Read data from a file made by self.print_data, add mutants to dataset. Ignores some things. DEPRECATED. 

        Populates most of the dataset total read/mutant count values correctly, but ignores unaligned and discarded reads, 
         anything involving special regions, mutant-merging information, and possibly other things.
        Cannot deal with gene annotation fields.
        Ignores single sequence/count fields; the total number of sequence variants is unreliable if you add to preexisting
        data (if the data originally listed 1 unique sequence, and new data adds another 2 sequences, is the total 2 or 3 
         unique sequences?  If assume_new_sequences is True, the total is old+new; if it's False, it's max(old,new)). 

        DEPRECATED: All datasets should now have pickled versions for easy reading into python - use those instead.  
         This method will be kept around for reading old datasets without pickled versions, but will NOT BE UPDATED.
        """
        if self.multi_dataset:  raise MutantError("read_data_from_file not implemented for multi-datasets!")
        for line in open(infile):
            # try to extract info from the comment lines that can't be reproduced from the data, ignore the rest
            # NOTE: won't be adding any missing information - use pickled files instead, THIS IS DEPRECATED.
            if line.startswith('#'):                                        
                line = line.strip('#').strip()
                if line.startswith("Reads discarded in preprocessing"):
                    data = line.split('\t')[-1].split(' ')[0]
                    try:                data = int(data)
                    except ValueError:  pass
                    self.summary.discarded_read_count = data
                if line.startswith("Unaligned reads") or line.startswith("Reads without a unique alignment"):
                    data = line.split('\t')[-1].split(' ')[0]
                    try:                data = int(data)
                    except ValueError:  pass
                    self.summary.non_aligned_read_count = data
                    # special case for when we don't know the specific unaligned categories, but we know 
                    #  total non-aligned is 0, so the specific categories must be 0 too:
                    if self.summary.non_aligned_read_count==0:  self.summary.unaligned, self.summary.multiple_aligned = 0,0
                if line.startswith("(read location with respect to cassette: which end, which direction)"):
                    data = line.split('\t')[-1].strip('()').split(', ')
                    self.summary.cassette_end = data[0]
                    self.summary.reads_are_reverse = {'reverse':True, 'forward':False, '?':'?'}[data[1]]
                if line.startswith("(total genes in genome annotation data)"):
                    data = line.split('\t')[-1].strip('()')
                    try:                data = int(data)
                    except ValueError:  pass
                    self.total_genes_in_genome = data
                continue
            # ignore special-comment and header lines, parse other tab-separated lines into values
            if line.startswith('<REGEX>#') or line.startswith('<IGNORE>'):  continue       
            if line.startswith('chromosome\tstrand\tmin_position\t'):       continue
            # parse non-comment tab-separated lines into values
            fields = line.split('\t')
            chromosome = fields[0]
            strand = fields[1]
            min_pos = int(fields[2])
            full_pos = fields[3]
            gene, orientation, gene_feature = fields[4:7]
            total_reads,perfect_reads,sequence_variants = [int(x) for x in fields[7:10]]
            # generate new mutant if necessary; add counts and gene info to mutant (USE IMMUTABLE POSITIONS BY DEFAULT)
            position = Insertion_position(chromosome, strand, full_position=full_pos, immutable=True)
            curr_mutant = self.get_mutant(position)
            curr_mutant.add_counts(total_reads,perfect_reads,sequence_variants,assume_new_sequences)
            curr_mutant.update_gene_info(gene, orientation, gene_feature)
            # get however many specific sequences/counts are listed (this is variable)
            sequence_fields = fields[10::2]
            count_fields = fields[11::2]
            for seq, count in zip(sequence_fields, count_fields):
                if int(count)>0:
                    assert seq!=''
                    curr_mutant.sequences_and_counts[seq] += int(count)
            # add to dataset total read/mutant counts
            summ = self.summary
    

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

    def _set_merged_genome_info(self, gene_annotation_header_values, total_genes_in_genome_values):
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

            # Merge source dataset mutant data into own multi-dataset mutants
            #  (get_mutant will create new empty mutant if one doesn't exist)
            for mutant in dataset_object:
                # TODO should a both-strand mutant and a + or -strand mutant at the same position be considered the same mutant during joining datasets??  Currently they're not, and this can give pretty odd results!  Make it an option?  Make it depend on readcounts?  I'm not really sure...  (Also see "TODO how should mutant lookup by position deal with both-strand mutants?")
                curr_mutant = self.get_mutant(mutant.position)
                curr_mutant.add_other_mutant_as_dataset(mutant, dataset_name, overwrite=overwrite, 
                                                        check_constant_data=check_gene_data)
            # MAYBE-TODO add option to only include mutants with non-zero reads (total or perfect) in all datasets?  Or should that only happen during printing?  Or do we even care?  If I ever want to do that, there was code for it in the old version of mutant_join_datasets.py (before 2012-04-26)

        # Merge any pieces of global information that's not per-dataset
        #  using getattr instead of just d.total_genes_in_genome because some older datasets don't HAVE the total_genes_in_genome
        #   attribute, and getattr lets me give a default of 0 when the attribute is missing 
        #   (and 0 is used as blank_value in the merge_values_to_unique call in self._set_merged_genome_info).
        self._set_merged_genome_info([d.gene_annotation_header for d in source_dataset_dict.values()], 
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

    def merge_other_dataset(self, other_dataset):
        """ Move all mutants from other_dataset to self (merge with existing same-position mutants if present).

        Currently NOT IMPLEMENTED if either dataset contains both-stranded mutants (merged opposite-tandems). 
        Only merges mutants with identical positions - 100-?, ?-101 and 100-101 are all considered distinct.
        Mutants are removed from other_dataset.
        """

        ### Merge basic data from summary
        # set name to None, as the safest option - only multi-datasets have names anyway, I think
        self.summary.dataset_name = None
        # for read/direction, keep the old value if old/new match, otherwise change to '?'
        if self.summary.cassette_end != other_dataset.summary.cassette_end:             self.summary.cassette_end = '?'
        if self.summary.reads_are_reverse != other_dataset.summary.reads_are_reverse:   self.summary.reads_are_reverse = '?'
        # MAYBE-TODO if merging a 3' and 5' dataset, maybe should use 'both' or something instead of '?' ?

        ### Merge extra readcount data, i.e. discarded/unaligned/removed/etc read counts from summary
        # using getattr because some older datasets don't have that, should be 0 by default
        other_discarded_other_end = getattr(other_dataset.summary,'discarded_other_end',0)
        self.summary.add_discarded_reads(other_dataset.summary.discarded_read_count, other_dataset.summary.discarded_wrong_start, 
                             other_dataset.summary.discarded_no_cassette, other_discarded_other_end, replace=False)
        self.summary.add_nonaligned_reads(other_dataset.summary.non_aligned_read_count, other_dataset.summary.unaligned, 
                                          other_dataset.summary.multiple_aligned, replace=False)
        self.summary.ignored_region_read_counts = add_dicts_of_ints(self.summary.ignored_region_read_counts, 
                                                                    other_dataset.summary.ignored_region_read_counts)

        ### Deal with the merging info from the summary:
        # Throw an exception if attempting to merge two datasets mutant-merging has been done on - mutant-merging really 
        #  should be done AFTER all the dataset merging, otherwise it may turn out inconsistent
        if self.summary.merged_opposite_tandems or sum(self.summary.merged_adjacent_same_strand_dict.values())\
           or other_dataset.summary.merged_opposite_tandems or sum(other_dataset.summary.merged_adjacent_same_strand_dict.values()):
            raise MutantError("Cannot merge together two datasets that had mutant-merging done - the results may be inconsistent!")
        # Blank the adjacent-counts to force a counting re-run; the merged counts should still be blank, we just checked that.
        self.summary.blank_adjacent_mutant_info()
        
        ### Merge non-mutant info from the dataset itself, non-summary (gene annotation header and total genes in genome)
        # (see comments on self._set_merged_genome_info in populate_multi_dataset for explanation of the getattr)
        self._set_merged_genome_info([d.gene_annotation_header for d in (self, other_dataset)], 
                                     [getattr(d,'total_genes_in_genome',0) for d in (self, other_dataset)])

        ### Merge the mutants
        if any([m.position.strand=='both' for m in self]):
            raise MutantError("Cannot run merge_other_dataset into a dataset with both-stranded mutants!")
        if any([m.position.strand=='both' for m in other_dataset]):
            raise MutantError("Cannot run merge_other_dataset from a dataset with both-stranded mutants!")
        for other_mutant in list(other_dataset):
            this_mutant = self.get_mutant(other_mutant.position)
            this_mutant.merge_mutant(other_mutant)
            other_dataset.remove_mutant(other_mutant)
        assert len(other_dataset) == other_dataset.summary.aligned_read_count == 0

        # LATER-TODO what to do with cases where one dataset has a both-strand mutant and the other has a +strand or -strand one?  Should probably merge them, or if not, give an error (make that an option?).  And what about mutants with similar positions, like ?-101 and 100-? and 100-101?  (This would be relevant when merging 5' and 3' sets, for instance).  See half-implemented get_matching_position_mutants method.

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
            if get_readcount(other_dataset.get_mutant(mutant.position)) >= readcount_min:
                self.remove_mutant(mutant.position)
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
            if get_readcount(other_dataset.get_mutant(mutant.position)) < readcount_min:
                self.remove_mutant(mutant.position)
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
                self.remove_mutant(mutant.position)
        # now redo adjacent-mutant counting, with same parameters as before (most get auto-detected by default)
        if self.summary.adjacent_max_distance is not None:
            cassette,other = self.summary.merging_which_chromosomes
            self.count_adjacent_mutants(count_cassette_chromosomes=cassette, count_other_chromosomes=other)
        # TODO really I shouldn't be removing mutants outright, just noting them as removed or something...  In that case should they or should they not show up in "for m in self"?  Probably not - they should have a separate dictionary?
        # TODO should I keep track of removed reads, and print in summary?  MAYBE.

    def find_genes_for_mutants(self, genefile, detailed_features=False, N_run_groups=3, verbosity_level=1):
        """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

        If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR.
        Read the file in N_run_groups passes to avoid using up too much memory/CPU.
        """ 
        if self.multi_dataset:  raise MutantError("find_genes_for_mutants not implemented for multi-datasets!")
        # MAYBE-TODO implement for multi-datasets?  The actual gene-finding would be easy, since it'd just work on 
        #  multi-dataset mutants instead of single-dataset ones; adding stuff to summary would be harder.

        # group all the mutants by chromosome, so that I can go over each chromosome in genefile separately
        #   instead of reading in all the data at once (which uses a lot of memory)
        mutants_by_chromosome = defaultdict(set)
        for mutant in self:
            mutants_by_chromosome[mutant.position.chromosome].add(mutant)
        self._find_genes_for_mutant_list(mutants_by_chromosome, genefile, detailed_features, N_run_groups, verbosity_level)

    def _find_genes_for_mutant_list(self, mutants_by_chromosome, genefile, detailed_features=False, 
                                    N_run_groups=3, verbosity_level=1):
        # First get the list of all chromosomes in the file, WITHOUT reading it all into memory
        with open(genefile) as GENEFILE:
            GFF_limit_data = GFF.GFFExaminer().available_limits(GENEFILE)
            chromosomes_and_counts = dict([(c,n) for ((c,),n) in GFF_limit_data['gff_id'].iteritems()])
            all_reference_chromosomes = set(chromosomes_and_counts.keys())

        # Now lump the chromosomes into N_run_groups sets with the feature counts balanced between sets, 
        #  to avoid using too much memory (by reading the whole file at once), 
        #   or using too much time (by reading the whole file for each chromosome/scaffold)
        chromosome_sets = split_into_N_sets_by_counts(chromosomes_and_counts, N_run_groups)

        ### go over all mutants on each chromosome, figure out which gene they're in (if any), keep track of totals
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
                    for mutant in mutants_by_chromosome[chromosome_record.id]:
                        gene_ID, orientation, feature = find_gene_by_pos_gff3(mutant.position, chromosome_record, 
                                                                              detailed_features, quiet=(verbosity_level==0))
                        mutant.gene, mutant.orientation, mutant.gene_feature = gene_ID, orientation, feature
                    if verbosity_level>1:   print "    ...found total %s genes."%(len(chromosome_record.features))
        if verbosity_level>1:   print "    found total %s genes in full genome."%(self.total_genes_in_genome)

        # for mutants in chromosomes that weren't listed in the genefile, use special values
        for chromosome in set(mutants_by_chromosome.keys())-set(all_reference_chromosomes):
            if not is_cassette_chromosome(chromosome):
                print 'Warning: chromosome "%s" not found in genefile data!'%(chromosome)
            for mutant in mutants_by_chromosome[chromosome]:
                mutant.gene,mutant.orientation,mutant.gene_feature = SPECIAL_GENE_CODES.chromosome_not_in_reference,'-','-'


    def _get_gene_annotation_dict(self, annotation_file, if_standard_Phytozome_file=None, custom_header=None, print_info=False):
        gene_annotation_dict, gene_annotation_header = parse_gene_annotation_file(annotation_file, 
                                                 standard_Phytozome_file=if_standard_Phytozome_file, header_fields=custom_header, 
                                                 strip_gene_fields_start='.t', verbosity_level=int(print_info))
        # store the annotation header in self.summary, for printing
        if gene_annotation_header:  self.gene_annotation_header = gene_annotation_header
        else:                       self.gene_annotation_header = 'GENE_ANNOTATION_DATA'
        return gene_annotation_dict

    @staticmethod
    def _add_gene_annotation_to_mutant(mutant, gene_annotation_dict):
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
        for gene in mutant.gene.split(MULTIPLE_GENE_JOIN):
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
        mutant.gene_annotation = joint_annotations
        if joint_annotations:   return True
        else:                   return False
        # TODO unit-test!

    def add_gene_annotation(self, annotation_file, if_standard_Phytozome_file=None, custom_header=None, print_info=False):
        """ Add gene annotation to each mutant, based on annotation_file. See parse_gene_annotation_file doc for detail."""
        # add the annotation info to each mutant (or nothing, if gene has no annotation)
        # MAYBE-TODO should I even store gene annotation in each mutant, or just keep a separate per-gene dictionary?
        gene_annotation_dict = self._get_gene_annotation_dict(annotation_file, if_standard_Phytozome_file, custom_header, print_info)
        # add the annotation info to each mutant (or nothing, if gene has no annotation) 
        N_annotated = 0
        for mutant in self:
            if_annotated = self._add_gene_annotation_to_mutant(mutant, gene_annotation_dict)
            if if_annotated:    N_annotated += 1
        if print_info:          print "Added %s annotations"%N_annotated
        elif not N_annotated:   print "Warning: No gene annotations found!"
        # LATER-TODO add this to the gene-info run-test case!


    ######### MUTANT-MERGING AND ADJACENT-COUNTING

    def _possibly_adjacent_positions(self, max_distance, same_strand_only, include_cassette_chromosomes, include_other_chromosomes):
        """ Generates all mutant position pairs that might be within max_distance of each other (same strand, or not).

        May or may not be well optimized (the most naive way would just be combinations(all_positions,2)).

        If same_strand_only is True, both-strand positions are counted as both +strand and -strand.
        """
        # using chromosome_sort_key to sort the chromosomes sensibly! Matters when we're printing merge/count details.
        for chromosome in sorted(self.summary.all_chromosomes, key=chromosome_sort_key):
            if include_cassette_chromosomes==False and is_cassette_chromosome(chromosome): continue
            if include_other_chromosomes==False and is_other_chromosome(chromosome):       continue

            all_positions = sorted([mutant.position for mutant in self if mutant.position.chromosome==chromosome])

            # MAYBE-TODO if same_strand_only, could split up by strand, too!  (both-strand ones should go on both lists?)
            # MAYBE-TODO could optimize this further by only considering pairs of close-by positions, instead of all pairs:
            #  positions are sorted by min_position and then strand, so if you sorted all_positions, 
            #  two mutants on the same strand X bp away can be at most 2X positions away on the sorted position list, 
            #  so use a sliding window of size 2X+1.  #  (or of size X+1 if they're all on the same strand already)

            # go over each pos1,pos2 pair (only once)
            for (pos1, pos2) in itertools.combinations(all_positions,2):
                yield (pos1, pos2)

    def _merge_two_mutants(self, pos1, pos2, mutant2_is_imperfect):
        """ Merge two mutants (given by position), while making sure not to mess up dictionary key immutability. 
        
        Doing the merge may require the position to change, and changing a dictionary key is BAD, 
         so first remove both mutants from the dataset, make the position mutable, 
         merge the mutants (using the Insertional_mutant.merge_mutant method, which does all the real work), 
         make the position immutable again, and finally add the merged mutant back to the dataset.
        """
        assert (pos1 in self and pos2 in self)
        mutant1 = self.get_mutant(pos1)
        mutant2 = self.get_mutant(pos2)
        self.remove_mutant(pos1)
        self.remove_mutant(pos2)
        mutant1.position.make_mutable_REMEMBER_CLEANUP_FIRST()
        mutant1.merge_mutant(mutant2, opposite_strand_tandem=(pos1.strand!=pos2.strand), other_is_imperfect=mutant2_is_imperfect)
        mutant1.position.make_immutable()
        self.add_mutant(mutant1)

    def _readcounts_sorted(self, val1, val2):
        """ Return a tuple of the two readcounts, largest first, for the two args (mutants, readcounts or positions) """
        # val1 and val2 can be either readcounts, or mutants, or positions
        try:
            val1/val2
            readcount1, readcount2 = val1, val2
        except TypeError:
            try:
                readcount1, readcount2 = val1.total_read_count, val2.total_read_count
            except AttributeError:
                readcount1, readcount2 = [self.get_mutant(pos).total_read_count for pos in (val1, val2)]
        # make sure readcount1 is larger
        if readcount1 < readcount2:
            readcount1, readcount2 = readcount2, readcount1
        return (readcount1, readcount2)

    def _readcount_ratio(self, val1, val2):
        """ Return readcount ratio between the two args (mutants, readcounts or positions), ordered so that ratio>=1. """
        readcount1, readcount2 = self._readcounts_sorted(val1, val2)
        return readcount1/readcount2 

    def _grab_mergeable_mutants(self, merge_max_distance, same_strand, max_count_ratio, min_count_ratio, 
                               merge_cassette_chromosomes, merge_other_chromosomes):
        """ Return merge-able mutant dict (to merge key into value) based on the conditions from arguments. 

        Conditions to include a mutant pair in the returned dict:
            - distance, min/max readcount ratio (higher over lower)
            - strand: if same_strand is True, only same-strand pairs; if False, only opposite-strand pairs; if None, both cases.
            - ignore/include cassette/non-nuclear chromosomes

        If one mutant could be merged into multiple different mutants, the target is chosen based on 1) lower distance, 
            2) more reads, 3) earlier position - print info to .
        """
        ### go over all mutant pairs, get a dict of which ones to merge (pos2:pos1 means pos2 should be merged INTO pos1)
        mutants_to_merge_key_into_val = {}
        for pos1,pos2 in self._possibly_adjacent_positions(merge_max_distance, True, merge_cassette_chromosomes, 
                                                           merge_other_chromosomes):
            # if the two positions are on different chromosomes or strands or aren't adjacent, skip
            if pos1.chromosome != pos2.chromosome:              continue
            if same_strand is True and pos1.strand != pos2.strand:      continue
            if same_strand is False and pos1.strand == pos2.strand:  continue
            distance = abs(pos2.min_position-pos1.min_position)
            if distance > merge_max_distance:                   continue
            # all other mutant pairs should be same-strand adjacent - up for merging depending on readcount ratio.
            # grab the two mutants based on position - I'm not doing any merging yet, so no mutants should be missing
            assert (pos1 in self and pos2 in self)
            mutant1 = self.get_mutant(pos1)
            mutant2 = self.get_mutant(pos2)

            # check readcount ratio - if acceptable, save for merging, otherwise go on to next pair
            mutant1_readcount, mutant2_readcount = mutant1.total_read_count, mutant2.total_read_count
            readcount_ratio = self._readcount_ratio(mutant1, mutant2)
            if min_count_ratio is not None and not readcount_ratio>=min_count_ratio:  continue
            if max_count_ratio is not None and not readcount_ratio<=max_count_ratio:  continue
            # want to merge mutant2 into mutant1, so make sure mutant1 is the one with more reads, swap if otherwise
            if mutant1_readcount < mutant2_readcount:
                pos1, pos2 = pos2, pos1
            # or if they have the same readcount numbers, make sure mutant1 is the one with the earlier position
            elif mutant1_readcount == mutant2_readcount:
                pos1, pos2 = sorted([pos2, pos1])
            # normally just save the mutants for merging, except special case below
            if pos2 not in mutants_to_merge_key_into_val:
                mutants_to_merge_key_into_val[pos2] = pos1
            # EXCEPT if pos2 is already set to be merged into pos3 (like three adjacent mutants with readcounts 100,1,100)
            #  see docstring for what to do then.
            #   LATER-TODO is this the right solution, or should I be doing something else with those cases?
            else:
                posA, posB = pos1, mutants_to_merge_key_into_val[pos2]
                distA, distB = [abs(posX.min_position-pos2.min_position) for posX in (posA,posB)]
                if distA<distB:         mutants_to_merge_key_into_val[pos2] = posA
                elif distA>distB:       mutants_to_merge_key_into_val[pos2] = posB
                else:
                    countA, countB = [self.get_mutant(posX).total_read_count for posX in (posA,posB)]
                    if countA>countB:   mutants_to_merge_key_into_val[pos2] = posA
                    elif countA>countB: mutants_to_merge_key_into_val[pos2] = posB
                    else:               mutants_to_merge_key_into_val[pos2] = min(mutants_to_merge_key_into_val[pos2],pos1)
        return mutants_to_merge_key_into_val

    def _pick_mutants_to_merge(self, mutants_to_merge, leave_N_mutants, leave_method, 
                               merge_cassette_chromosomes, merge_other_chromosomes, 
                               leave_auto_multiplier, distance_max, readcount_ratio_sort_reverse):
        """ Help function to decide which mutants to merge based on leave_N_mutants/leave_method/etc - 
        see merge_adjacent_mutants/merge_opposite_tandem_mutants docstrings for details. 
        """
        ### based on leave_N_mutants and leave_method, determine which ones to actually merge
        # if we're just using the ratio directly, or leave_N_mutants is 0 (i.e. merge all mutants, leave none), 
        # just merge all the mutant pairs
        if leave_N_mutants=='use_ratio' or leave_N_mutants==0:  pass
        # for other leave_N_mutants values, decide which mutant pairs to merge and leave
        else:
            # determine the actual leave_N_mutants number if 'auto' (see docstring for how this is done and why)
            #  (only run self.count_adjacent_mutants if it hasn't been run before)
            if leave_N_mutants=='auto':
                if self.summary.adjacent_max_distance is None:
                    self.count_adjacent_mutants(count_cassette_chromosomes=merge_cassette_chromosomes, 
                                                count_other_chromosomes=merge_other_chromosomes) 
                opposite_away =   sum([self.summary.adjacent_opposite_away_dict[i]   for i in range(1, distance_max+1)])
                opposite_toward = sum([self.summary.adjacent_opposite_toward_dict[i] for i in range(1, distance_max+1)])
                leave_N_mutants = int(round(min(opposite_away+opposite_toward, 2*opposite_toward)*leave_auto_multiplier))
            # determine which mutants to merge/leave - method depends on leave_method
            # convert mutants_to_merge_key_into_val from dict to list to make it easier to manipulate
            mutants_to_merge = list(mutants_to_merge.items())
            # sort mutants_to_merge based on the factor specified by leave_method - random, or reverse sort by readcount ratio
            #  (reverse, since we want to first sort the pairs with very different readcounts, i.e. high ratios)
            if leave_method=='random':      
                random.shuffle(mutants_to_merge)
            elif leave_method=='by_ratio':  
                mutants_to_merge.sort(key = lambda (x,y): self._readcount_ratio(x,y), reverse=readcount_ratio_sort_reverse)
            else:                           
                raise MutantError("Unknown leave_method value %s!"%leave_method)
            # take the first N mutant pairs for merging, leaving leave_N_mutants unmerged; convert mutants_to_merge back into a dict
            if leave_N_mutants:     mutants_to_merge = dict(mutants_to_merge[:-leave_N_mutants])
            else:                   mutants_to_merge = dict(mutants_to_merge)
        return mutants_to_merge

    def _merge_mutant_list(self, mutants_to_merge_key_into_val, count_merged_mutant_as_imperfect, 
                           OUTPUT=None, description="", join_string="into", extra_check=None, extra_note=""):
        """ Given a dict of mutants to be merged (key into val), do the merging; return number merged and other info.

        Merging is transitive - if A should be merged into B and B into C, just merge A into C 
         (avoids errors if B has been merged into C already).
        Return total number merged, distance:N_merged dict, and distance:pair_readcounts_list dict.
        Write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no printing).
        """
        N_merged_dict = defaultdict(int)
        merged_readcounts_dict = defaultdict(list)
        for pos2, pos1 in sorted(mutants_to_merge_key_into_val.items(), key = lambda (x,y): min(x,y)):
            # check if pos1 (the merge target) is set to be merged into some pos3 - if so, just merge pos2 into pos3 
            #  directly (especially since maybe pos1 was merged already and is gone now, which would cause an error)
            while pos1 in mutants_to_merge_key_into_val:
                pos1 = mutants_to_merge_key_into_val[pos1]
            # grab the actual mutants (they SHOULD always be there, unless something strange is going on!)
            assert (pos1 in self and pos2 in self)
            mutant1 = self.get_mutant(pos1)
            mutant2 = self.get_mutant(pos2)
            mutant1_readcount, mutant2_readcount = mutant1.total_read_count, mutant2.total_read_count
            distance = abs(pos2.min_position-pos1.min_position)
            readcount_ratio = self._readcount_ratio(mutant1, mutant2)
            # write merging info
            OUTPUT.write(" MERGING %s mutants: %s %s %s, %s and %s reads.\n"%(description, pos2, join_string, pos1,
                                                                                mutant2_readcount, mutant1_readcount))
            # do actual merging, save counts!
            self._merge_two_mutants(pos1, pos2, count_merged_mutant_as_imperfect)
            N_merged_dict[distance] += 1
            merged_readcounts_dict[distance].append((mutant1_readcount, mutant2_readcount))
        N_merged = sum(N_merged_dict.values())
        assert N_merged == len(mutants_to_merge_key_into_val), "Merged a different number than given!"
        return N_merged, N_merged_dict, merged_readcounts_dict

    def merge_adjacent_mutants(self, merge_max_distance=1, leave_N_mutants='auto', min_count_ratio=100, leave_method='by_ratio',
                               merge_cassette_chromosomes=False, merge_other_chromosomes=False, OUTPUT=sys.stdout):
        """ Merge adjacent mutants based on strand, distance, and count ratio or number; save counts.

        The idea is to merge mutants that are likely to actually be a single mutant, but have different positions due to PCR or
         deepseq indels or such - mutants that are close to each other and on the same strand.  
        How likely it is that a pair of nearby same-strand mutants is actually one mutant depends on:
            - the distance between the two mutants - anything larger than 1bp is iffy, really
            - the ratio of the mutant readcounts - indels are assumed to be uncommon, so if one of the mutants has 1000x more
                reads than the other, that seems more likely to be an indel than if both have the same number of reads.
                (readcount ratio always takes the higher readcount over the lower, so 1,100 and 100,1 both have a ratio of 100)
            - the numbers of close-by OPPOSITE-STRAND mutants - if there are a lot more close-by same-strand mutants
                than opposite-strand ones, we can reasonably assume they're not real, since randomly you'd expect both cases 
                to be equally likely.  There's an additional complication in that opposite-strand "toward-facing" mutants really
                do have to always be real (up to a distance of about 20bp, anyway - higher than that and they could be 
                 opposite-tandem cassettes with insertions in between), but "away-facing" ones could be opposite-tandem cassettes 
                 with deletions, so we're balancing the two when this is taken into account 
                (see the "if leave_N_mutants is 'auto'" case explanation below)

        The implementation is a bit complicated.  
         * we're ONLY looking at mutant pairs on the same strand, within merge_max_distance bases of each other.
         * the basic functionality depends on the leave_N_mutants value, and then the min_count_ratio and leave_method values:
            - if leave_N_mutants is 'use_ratio', all mutant pairs with a readcount ratio >=min_count_ratio will be merged, 
                regardless of how many will be merged/left.  leave_method must be "by_ratio" in this case.
            - if leave_N_mutants is a number N>=0, enough mutants will be merged so that N unmerged pairs are left. 
              Which pairs are merged/left depends on leave_method:  
                - if leave_method is 'random', mutants to merge will be chosen randomly
                - if leave_method is 'by_ratio', mutants with highest readcount ratios will be merged
                - other methods may be implemented later!
            - if leave_N_mutants is 'auto', the number of mutants to leave will be determined based on opposite-strand adjacent
                mutant counts, and then the merging will proceed as if leave_N_mutants was that number. 
               The number will be either the total number of opposite-strand adjacent mutants within merge_max_distance, 
                or 2x the number of opposite-strand adjacent toward-facing mutants within merge_max_distance=1,
                whichever is smaller (this makes sense biologically - the toward-facing ones are definitely really separate mutants,
                 and the away-facing ones may not be, so if there are more away- than toward-facing, we don't want to use them, 
                 but if there are fewer, we can assume they're also real, and then it's better to use both numbers for lower error).
                (The counts will be determined by self.count_adjacent_mutants.)
         * if merge_cassette_chromosomes or merge_other_chromosomes is False, don't do merging on cassette chromosomes 
             and other non-genome chromosomes (chloroplast, mitochondrial, etc) respectively.
         * when merging two mutants, merge their read counts, sequences etc, and remove the lower-read-count one from dataset.

        Cannot be run after merge_opposite_tandem_mutants - NOT IMPLEMENTED for both-stranded mutants.

        Write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no writing).
         Save information on merged counts to various self.summary attributes.
        """
        ### check inputs (I just put it in an if block for folding-readability purposes)
        if True:
            allowed_leave_N_mutants_vals = ['auto', 'use_ratio']
            allowed_leave_method_vals = ['by_ratio', 'random']
            if self.multi_dataset:  
                raise MutantError("merge_adjacent_mutants not implemented for multi-datasets!")
            if leave_N_mutants not in allowed_leave_N_mutants_vals and leave_N_mutants<0:
                raise MutantError("Bad leave_N_mutants value!")
            if leave_N_mutants==allowed_leave_N_mutants_vals[1] and leave_method!=allowed_leave_method_vals[0]:
                raise MutantError("Bad leave_N_mutants and leave_method value combination!")
            if min_count_ratio < 1 and min_count_ratio is not None:
                raise MutantError("min_count_ratio must be at least 1, or None!")
            if leave_method not in ['by_ratio', 'random']:  
                raise MutantError("Bad by_ratio value in merge_adjacent_mutants!")
            # if we're not calculating leave_N_mutants based on min_count_ratio, min_count_ratio should be ignored - set to None
            if leave_N_mutants!=allowed_leave_N_mutants_vals[1]:
                min_count_ratio = None
        ### check if values are consistent with previous counting/merging (if block for folding-readability purposes)
        if True:
            if self.summary.adjacent_max_distance not in (None, merge_max_distance):
                raise MutantError("Cannot change max adjacent distance once it's been set! (by merging or counting) "
                                  "Currently %s, can't change to %s."%(self.summary.adjacent_max_distance, merge_max_distance))
                # MAYBE-TODO or allow this and just re-run count_adjacent_mutants afterward with new value?
            if not self.summary.merging_which_chromosomes in [(merge_cassette_chromosomes,merge_other_chromosomes), (None,None)]:
                raise MutantError("Cannot change which chromosomes are subject to mutant-merging! "
                                 +"Currently %s for cassette, %s for non-nuclear. "%self.summary.merging_which_chromosomes
                                 +"Trying to change to %s and %s."%(merge_cassette_chromosomes, merge_other_chromosomes))
            if any([mutant.position.strand=='both' for mutant in self]):
                raise MutantError("merge_adjacent_mutants should be run BEFORE merge_opposite_tandem_mutants "
                                 +"- NOT IMPLEMENTED for both-stranded mutants!")
                # LATER-TODO implement for both-stranded mutants too? In that case all strand-comparisons should be
                #   rewritten so 'both' and either + or - would compare as the same, or something, 
                #  and then for the readcount ratios, for the both-strand mutant, only the same-strand readcount 
                #   should be used! (from mutant.original_strand_readcounts)
                # BUT what would I do if I have a both-stranded mutant with fewer reads than an adjacent single-stranded mutant??
                #  Can't really merge a both-stranded mutant into a single-stranded one, that would make NO SENSE,
                #   since the other-strand reads cannot be a result of an indel.  
                #  COULD read a single-strand mutant into a both-strand one, though...  
                #  But in that case why not just do adjacent-merging before tandem-merging?  
                #  So maybe this should continue not being allowed, and I should rewrite the error message to make it clear that
                #   it's actually a bad idea, not just not implemented. 

        ### write basic info to OUTPUT (if block for readability)
        if True:
            if OUTPUT is None: OUTPUT = FAKE_OUTFILE
            OUTPUT.write("# Merging adjacent mutants: max distance %s, leave_N_mutants %s, "%(merge_max_distance, leave_N_mutants) 
                         +"min_count_ratio %s, leave_method %s\n"%(min_count_ratio, leave_method))
            OUTPUT.write(" (The merged mutant will have the position of whichever original mutant had more reads)\n")
            # MAYBE-TODO come up with a better solution for both low min_count_ratio cases below, or is this all right?
            if min_count_ratio == 1:  
                OUTPUT.write(" (Warning: minimum ratio is 1, so sometimes both mutants will have the same read count "
                             "- the earlier position is used in that case; and if many mutants in a row have the same, "
                             "read count, all will be merged into the first one.)\n")
            elif min_count_ratio < 2:  
                OUTPUT.write(" (Warning: minimum ratio is below 2, so there may be odd behavior with multiple adjacent mutants"
                             " in a row - for instance mutants with 9,8,7 reads would get merged into the first one, not the "
                             "middle one - not clear what would be better. [example numbers may not match real min ratio])\n")
            OUTPUT.write(" (%sincluding cassette chromosomes)"%('' if merge_cassette_chromosomes else 'not ')
                        +" (%sincluding non-nuclear chromosomes)\n"%('' if merge_other_chromosomes else 'not '))

        ### go over all mutant pairs, get a dict of which ones to possibly merge (adjacent position, same strand, min ratio, etc)
        mutants_to_merge_key_into_val = self._grab_mergeable_mutants(merge_max_distance, same_strand=True, 
                                                                    max_count_ratio=None, min_count_ratio=min_count_ratio, 
                                                                    merge_cassette_chromosomes=merge_cassette_chromosomes, 
                                                                    merge_other_chromosomes=merge_other_chromosomes)

        ### based on leave_N_mutants and leave_method, determine which ones to actually merge
        mutants_to_merge_key_into_val = self._pick_mutants_to_merge(mutants_to_merge_key_into_val, leave_N_mutants, leave_method, 
                                                    merge_cassette_chromosomes, merge_other_chromosomes, leave_auto_multiplier=1, 
                                                    distance_max=merge_max_distance, readcount_ratio_sort_reverse=True)
        # MAYBE-TODO write more info to OUTPUT about how the number of mutants to merge is determined, etc?

        ### now that we have all the pairs-to-merge, do the actual merging!
        N_merged, N_merged_dict, merged_readcounts_dict = self._merge_mutant_list(mutants_to_merge_key_into_val, 
                                                                              count_merged_mutant_as_imperfect=True, OUTPUT=OUTPUT, 
                                                                              description="same-strand adjacent", join_string="into")
        ### add overall counts to dataset summary
        self.summary.adjacent_max_distance = merge_max_distance
        self.summary.merging_which_chromosomes = (merge_cassette_chromosomes, merge_other_chromosomes)
        for distance in set(N_merged_dict.keys() + merged_readcounts_dict.keys()):
            self.summary.merged_adjacent_same_strand_dict[distance] += N_merged_dict[distance]
            self.summary.merged_adjacent_same_strand_readcounts_dict[distance].extend(merged_readcounts_dict[distance])
        # sort the lists inside the readcounts dict, since the order doesn't matter
        self.summary.merged_adjacent_same_strand_readcounts_dict = sort_lists_inside_dict(
                                                                self.summary.merged_adjacent_same_strand_readcounts_dict)
        # Always do an adjacent-mutant re-count after merging!
        self.count_adjacent_mutants(count_cassette_chromosomes=merge_cassette_chromosomes, 
                                    count_other_chromosomes=merge_other_chromosomes, OUTPUT=None)
        OUTPUT.write("# Finished merging adjacent mutants: %s pairs merged\n"%self.summary.merged_adjacent_pairs) 

    # TODO was there some other bug when doing tandem-merging and adjacent-merging at once?  I think something weird came up in actual data analysis - see ../../1206_Ru-screen1_deepseq-data-early/notes.txt  "Mutants" section.

    # MAYBE-TODO implement mutant-merging for multi-datasets?  Either independent merging for each component dataset (but that's a bit silly, they could just be done on the singles), OR actually potentially useful different functionality of merging the SAME mutants in all the datasets, using either one particular dataset or the sum of all datasets as reference to decide which pairs to merge.

    # MAYBE-TODO add checking for more complicated cases like three mutants next to each other?  Apparently this does happen, see mutant_pool_screens/1211_positions_Ru-screen1-for-paper

    # TODO keep the merging info in the dataset in addition to printing it!  Maybe just make a separate merged_mutant dictionary with copies of the original before-merging mutants?  And a list of the merged mutant position pairs and the resulting mutant position pairs, so they can all be looked up if desired.

    def merge_opposite_tandem_mutants(self, leave_N_mutants='auto', max_count_ratio=None, leave_method='by_ratio',
                                      merge_cassette_chromosomes=False, merge_other_chromosomes=False, OUTPUT=sys.stdout):
        """ Merge opposite-strand tandem mutants (in same position but opposite strands) depending on ratio; set strand to 'both'. 

        The idea is to merge mutants that are likely to actually be a single mutant that has two opposite-direction cassette copies
         that were inserted together, so that it generates flanking regions on both sides, mapping to the same position
         but on opposite strands.  
        How likely it is that a pair of same-position opposite-strand mutants is actually one mutant depends on:
            - the ratio of the mutant readcounts - the 
                reads than the other, that seems more likely to be an indel than if both have the same number of reads.
                (readcount ratio always takes the higher readcount over the lower, so 1,100 and 100,1 both have a ratio of 100)
            - the numbers of close-by OPPOSITE-STRAND mutants - if there are a lot more same-position opposite-strand mutants
                than 1bp-distant opposite-strand ones, we can reasonably assume they're not real, since randomly you'd expect 
                both cases to be equally likely. Additional complication 1bp-distant opposite-strand mutants can be "toward-facing"
                or "away-facing", and these two cases are treated differently - see merge_adjacent_mutants docstring for detail.

        The implementation is a bit complicated.  
         * we're ONLY looking at mutant pairs on the opposite strands, in the same position
         * for how leave_N_mutants, max_count_ratio, leave_method and merge_*_chromosomes work, see merge_adjacent_mutants docstring
            (except that max_count_ratio is the opposite of merge_adjacent_mutants's min_count_ratio, and we preferentially merge
             mutants with LOWER readcount ratios, i.e. more similar readcounts)
         * when merging two mutants, merge their read counts, sequences etc; remove both from dataset, and add a new one to the 
            dataset, with 'both' for strand and gene orientation, and a  original_strand_readcounts attribute containing the 
            original +/- strand readcounts. If the full positions of the two mutants were 100-? and ?-101, the new one is 100-101.

        Write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no writing).
         Save information on merged counts to various self.summary attributes.
        """
        if self.multi_dataset:  raise MutantError("merge_opposite_tandem_mutants not implemented for multi-datasets!")
        if not self.summary.merging_which_chromosomes in [(merge_cassette_chromosomes,merge_other_chromosomes), 
                                                          (None,None)]:
            raise MutantError("Cannot change which chromosomes are subject to mutant-merging! "
                              "Currently %s for cassette, %s for non-nuclear."%self.summary.merging_which_chromosomes)
        if max_count_ratio < 1 and max_count_ratio is not None:
            raise MutantError("max_count_ratio must be at least 1, or None!")
        # if we're not calculating leave_N_mutants based on ratio, ratio should be ignored - set to None
        if leave_N_mutants!='use_ratio':
            max_count_ratio = None
        if OUTPUT is None: OUTPUT = FAKE_OUTFILE

        OUTPUT.write("# Merging opposite-strand same-position mutants (presumably tail-to-tail tandems): "
                    +"leave_N_mutants %s, max_count_ratio %s, leave_method %s\n"%(leave_N_mutants, max_count_ratio, leave_method)
                    +" (%sincluding cassette chromosomes)"%('' if merge_cassette_chromosomes else 'not ')
                    +" (%sincluding non-nuclear chromosomes)\n"%('' if merge_other_chromosomes else 'not '))

        ### go over all mutant pairs, get a dict of which ones to possibly merge (same position, opposite strands, max ratio, etc)
        mutants_to_merge = self._grab_mergeable_mutants(merge_max_distance=0, same_strand=False, 
                                                        max_count_ratio=max_count_ratio, min_count_ratio=None, 
                                                        merge_cassette_chromosomes=merge_cassette_chromosomes, 
                                                        merge_other_chromosomes=merge_other_chromosomes)

        ### based on leave_N_mutants and leave_method, determine which ones to actually merge
        mutants_to_merge = self._pick_mutants_to_merge(mutants_to_merge, leave_N_mutants, leave_method, merge_cassette_chromosomes, 
                           merge_other_chromosomes, leave_auto_multiplier=1/2, distance_max=1, readcount_ratio_sort_reverse=False)
        # MAYBE-TODO write more info to OUTPUT about how the number of mutants to merge is determined, etc?

        ### now that we have all the pairs-to-merge, do the actual merging!
        N_merged, N_merged_dict, merged_readcounts_dict = self._merge_mutant_list(mutants_to_merge, 
                                                              count_merged_mutant_as_imperfect=False, OUTPUT=OUTPUT, 
                                                              description="opposite-strand same-position tandem", join_string="and")
        assert (N_merged_dict == {0: N_merged} or (N_merged==0 and N_merged_dict=={}))

        ### add overall counts to dataset summary
        self.summary.merged_opposite_tandems += N_merged
        self.summary.merged_opposite_tandems_readcounts.extend(merged_readcounts_dict[0])
        # need to sort the ratio list, since the order doesn't matter
        self.summary.merged_opposite_tandems_readcounts.sort()
        self.summary.merging_which_chromosomes = (merge_cassette_chromosomes, merge_other_chromosomes)
        # Always do an adjacent-mutant re-count after merging!
        self.count_adjacent_mutants(count_cassette_chromosomes=merge_cassette_chromosomes, 
                                    count_other_chromosomes=merge_other_chromosomes, OUTPUT=None)
        OUTPUT.write("# Finished merging opposite-strand same-position mutants: %s pairs merged\n"%\
                     self.summary.merged_opposite_tandems)

    def count_adjacent_mutants(self, max_distance_to_print=None, max_distance_to_count=10000, 
                               count_cassette_chromosomes=False, count_other_chromosomes=False, 
                               different_parameters=False, OUTPUT=None):
        """ Count various categories of adjacent mutants (as dist:count dictionaries); save data to self.summary. 
        
        All mutant pairs with distance<=max_distance_to_count are counted, in a distance-keyed dictionary.
        There are four categories, counted separately:
            - adjacent mutants on the same strand
            - same-position mutants on opposite strands (100-? and ?-101 - both have min_position 100)
            - adjacent mutants on opposite strands, facing away from each other (100-? and ?-103)
            - adjacent mutants on opposite strands, facing toward each other (?-101 and 102-?)
            (the reason for differentiating between the "away" and "toward" cases is that 
             the "away" case could be a tail-to-tail tandem with a genomic deletion in the middle, 
             but the "toward" case really HAS TO be random unrelated mutants (right?), so it can be used as a reference.)
        For each category, two pieces of information are generated, and stored in self.summary:
            - a distance:adjacent_pair_count dictionary
            - a distance:list_of_pair_readcount_readcounts dictionary
            (except for the same-position opposite-strand mutants, which always have distance 0, 
             so there's simply a single adjacent_pair_count and single list_of_pair_readcount_readcounts.) 
          self.summary provides properties/methods to deal with the data.

        If count_cassette_chromosomes or count_other_chromosomes is False, don't include in the count cassette chromosomes 
         and other non-nuclear chromosomes (chloroplast, mitochondrial, etc) respectively.
        If different_parameters is False, check that adjacent_max_distance, count_cassette_chromosomes, 
         and count_other_chromosomes match previously used merging/counting settings as saved in self.summary.

        Write details/totals to OUTPUT (should be an open filehandle; pass sys.stdout if desired, or None for no printing),
         but ONLY for pairs with distance<=max_distance_to_print (for both single lines and total counts).
        """
        if self.multi_dataset:  raise MutantError("count_adjacent_mutants not implemented for multi-datasets!")

        # set the two max_distance arg values, then make sure count < print
        if max_distance_to_count is None:   max_distance_to_count = float('inf')
        if max_distance_to_print is None:   
            if self.summary.adjacent_max_distance is None:  max_distance_to_print = 1
            else:                                           max_distance_to_print = self.summary.adjacent_max_distance
        if max_distance_to_count < max_distance_to_print:
            raise MutantError("max_distance_to_count must be at least as high as max_distance_to_print!. "
                              "Current values: %s and %s."%(max_distance_to_count, max_distance_to_print))

        # make sure the parameters are consistent with the ones already stored; if none are stored, save current ones.
        # MAYBE-TODO set it so if count_cassette_chromosomes or count_other_chromosomes are None, they get set to the self.summary.merging_which_chromosomes values?
        if not different_parameters:
            if self.summary.merging_which_chromosomes == (None,None):
                self.summary.merging_which_chromosomes = (count_cassette_chromosomes, count_other_chromosomes)
            elif not self.summary.merging_which_chromosomes == (count_cassette_chromosomes,count_other_chromosomes):
                raise MutantError("If trying to change which chromosomes are included when doing adjacent-mutant-count "
                                  "compared to mutant-merging, use different_parameters argument."
                                  "(Merging: %s for cassette, %s for non-nuclear.)"%self.summary.merging_which_chromosomes)
            self.summary.merging_which_chromosomes = (count_cassette_chromosomes, count_other_chromosomes)
            if self.summary.adjacent_max_distance is None:
                self.summary.adjacent_max_distance = max_distance_to_print
            elif self.summary.adjacent_max_distance > max_distance_to_count:
                raise MutantError("Using a lower adjacent-mutant max distance for adjacent-mutant-count than used for" 
                                  " mutant-merging!  Use different_parameters argument if this is on purpose."
                                  "(Merging: max distance %s.)"%self.summary.adjacent_max_distance)
            
        if OUTPUT is None: OUTPUT = FAKE_OUTFILE
        max_distance_info = ("actually counting mutants up to distance %s, but only printing info for up to distance %s "
                         "- you can get the full data from the pickle file")%(max_distance_to_count, max_distance_to_print)
        OUTPUT.write("# Counting adjacent mutants: %s\n"%max_distance_info
                    +" (%sincluding cassette chromosomes)"%('' if count_cassette_chromosomes else 'not ')
                    +" (%sincluding non-nuclear chromosomes)\n"%('' if count_other_chromosomes else 'not '))
        adjacent_same_strand_dict = defaultdict(int)
        adjacent_same_strand_readcounts_dict = defaultdict(list)
        same_position_opposite_strands = 0
        same_position_opposite_strands_readcounts = []
        adjacent_opposite_strands_away_dict = defaultdict(int)
        adjacent_opposite_strands_away_readcounts_dict = defaultdict(list)
        adjacent_opposite_strands_toward_dict = defaultdict(int)
        adjacent_opposite_strands_toward_readcounts_dict = defaultdict(list)

        for pos1,pos2 in self._possibly_adjacent_positions(max_distance_to_count, False, count_cassette_chromosomes, 
                                                           count_other_chromosomes):
            # if the two positions are on different chromosomes or aren't adjacent, skip
            if pos1.chromosome != pos2.chromosome:  continue
            distance = abs(pos2.min_position-pos1.min_position)
            if distance > max_distance_to_count:    continue
            readcount1 = self.get_mutant(pos1).total_read_count
            readcount2 = self.get_mutant(pos2).total_read_count
            sorted_readcounts = tuple(sorted([readcount1,readcount2], reverse=True))
            # same position, opposite strands
            if pos1.min_position==pos2.min_position:
                assert pos1.strand != pos2.strand, "Two mutants with same position and strand shouldn't happen!"
                assert 'both' not in (pos1.strand,pos2.strand), "A both-strand mutant can't be same-position to another!"
                same_position_opposite_strands += 1
                same_position_opposite_strands_readcounts.append(sorted_readcounts)
                OUTPUT.write("  opposite-strand same-position tandem mutants: %s and %s, %s and %s reads.\n"%(pos1,pos2,
                                                                                                         readcount1,readcount2))
            # adjacent positions, same strand - remember that both-strand mutants are same-strand with everything!
            elif pos1.strand == pos2.strand or 'both' in (pos1.strand,pos2.strand): 
                assert pos1.min_position!=pos2.min_position, "Two mutants with same position and strand shouldn't happen!"
                adjacent_same_strand_dict[distance] += 1
                adjacent_same_strand_readcounts_dict[distance].append(sorted_readcounts)
                if distance <= max_distance_to_print:
                    OUTPUT.write("  same-strand adjacent mutants: %s and %s, %s and %s reads.\n"%(pos1,pos2, readcount1,readcount2))
            # adjacent positions, opposite strands
            else:
                assert pos1.min_position != pos2.min_position and pos1.strand != pos2.strand, "Not adjacent-opposite!"
                assert 'both' not in (pos1.strand,pos2.strand), "A both-strand mutant can't be opposite-strand to another!"
                #  ?-X is -strand, X-? is +strand.  
                # away  =  100-? and ?-103  =  +strand has a lower position than -strand;  toward  =  other way around
                first_pos = min([pos1, pos2])
                if first_pos.strand == '+':
                    adjacent_opposite_strands_away_dict[distance] += 1
                    adjacent_opposite_strands_away_readcounts_dict[distance].append(sorted_readcounts)
                    if distance <= max_distance_to_print:
                        OUTPUT.write('  adjacent opposite-strand away-facing mutants: %s and %s, %s and %s reads.\n'%(pos1,pos2,
                                                                                                         readcount1,readcount2))
                else:
                    adjacent_opposite_strands_toward_dict[distance] += 1
                    adjacent_opposite_strands_toward_readcounts_dict[distance].append(sorted_readcounts)
                    if distance <= max_distance_to_print:
                        OUTPUT.write('  adjacent opposite-strand toward-facing mutants: %s and %s, %s and %s reads.\n'%(pos1,pos2,
                                                                                                         readcount1,readcount2))

        # add results to dataset summary; sort all the ratio-lists, since the order doesn't matter
        self.summary.adjacent_same_strand_dict = adjacent_same_strand_dict
        self.summary.adjacent_same_strand_readcounts_dict = sort_lists_inside_dict(adjacent_same_strand_readcounts_dict)
        self.summary.same_position_opposite = same_position_opposite_strands
        self.summary.same_position_opposite_readcounts = sorted(same_position_opposite_strands_readcounts)
        self.summary.adjacent_opposite_away_dict = adjacent_opposite_strands_away_dict
        self.summary.adjacent_opposite_away_readcounts_dict = sort_lists_inside_dict(
                                                                        adjacent_opposite_strands_away_readcounts_dict)
        self.summary.adjacent_opposite_toward_dict = adjacent_opposite_strands_toward_dict
        self.summary.adjacent_opposite_toward_readcounts_dict = sort_lists_inside_dict(
                                                                        adjacent_opposite_strands_toward_readcounts_dict)
        OUTPUT.write("# Finished counting (%s). %s\n"%(max_distance_info.replace('counting','counted'), 
                                                       self.summary.adjacent_mutant_summary(long_version=True, 
                                                            add_total_counts=False, max_distance=max_distance_to_print)))

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
                      merge_boundary_features=True, count_cassette=True, count_other=True):
        """ Print basic summary info about the dataset/s: read and mutant counts and categories, gene numbers, etc.

        Prints tab-separated table for multi-datasets-1.
        Prints to stdout by default, can also pass an open file object).
        """
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
        
        DVG.append((header_prefix+"Mutant merging/counts (deciding when different-position reads should be one mutant)", 
                    lambda summ: '' ))
        DVG.append((line_prefix+" (adjacent-merging/counting max distance):", 
                    lambda summ: "(%s)"%summ.adjacent_max_distance ))
        DVG.append((line_prefix+" (if we're including mutants in cassette and in non-nuclear chromosomes):", 
                    lambda summ: "(%s, %s)"%summ.merging_which_chromosomes )) 
        DVG.append((line_prefix+"merged same-strand adjacent mutant pairs and opposite-strand tandem pairs:", 
                    lambda summ: "%s, %s"%(summ.merged_adjacent_pairs, summ.merged_opposite_tandems) ))
        DVG.append((line_prefix+"remaining same-position opposite-strand pairs (if not merged as tandems):", 
                    lambda summ: "%s"%summ.same_position_opposite )) 
        DVG.append((line_prefix+'remaining adjacent opposite-strand "toward-facing" pairs (those are definitely real):', 
                    lambda summ: "%s"%summ.adjacent_opposite_toward )) 
        DVG.append((line_prefix+'remaining adjacent opposite-strand "away-facing" pairs (% of toward-facing):', 
                    lambda summ: value_and_percentages(summ.adjacent_opposite_away, 
                                                       [summ.adjacent_opposite_toward], percentage_format_str='%.0f') )) 
        DVG.append((line_prefix+'remaining adjacent same-strand unmerged pairs (% of 2*toward-facing):', 
                    lambda summ: value_and_percentages(summ.adjacent_same_strand, 
                                                       [2*summ.adjacent_opposite_toward], percentage_format_str='%.0f') )) 

        DVG.append((header_prefix+"Distinct mutants (read groups) by cassette insertion position:", 
                    lambda summ: "%s"%summ.N_mutants ))
        DVG.append((line_prefix+"(mutants with 2+, 10+, 100+, 1000+ reads):",
                    lambda summ: "(%s, %s, %s, %s)"%tuple([summ.N_mutants_over_readcount(X) for X in (2,10,100,1000)]) ))
        DVG.append((line_prefix+"(read location with respect to cassette: which end, which direction):", 
                    lambda summ: "(%s, %s)"%(summ.cassette_end, 
                                             {'?': '?', True: 'reverse', False: 'forward'}[summ.reads_are_reverse]) ))
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
                                *[set(summ.merged_gene_feature_counts(merge_boundary_features)) for summ in summaries])):
                DVG.append((line_prefix+"Mutant cassettes in gene feature %s (%% of ones in genes):"%feature, 
                            lambda summ,feature=feature: value_and_percentages(
                                                        summ.merged_gene_feature_counts(merge_boundary_features)[feature], 
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


    def _sort_data(self, sort_data_by=None):
        """ Sort the mutants by position or readcount, or leave unsorted. """
        all_mutants = iter(self)
        if sort_data_by=='position':
            sorted_data = sorted(all_mutants, key = lambda m: m.position)
            # x.position here is an Insertion_position object and has a sensible cmp function
        elif sort_data_by=='read_count':
            if self.multi_dataset:  
                raise MutantError("Sorting by readcount in print_data not implemented for multi-datasets!")
            sorted_data = sorted(all_mutants, key = lambda m: (m.total_read_count,m.perfect_read_count,m.position), reverse=True)
        else:
            sorted_data = all_mutants
        return sorted_data

    def print_data(self, OUTPUT=sys.stdout, sort_data_by=None, N_sequences=None, header_line=True, header_prefix="# "):
        """ Print full data, one line per mutant: position data, gene info, read counts, optionally sequences.
        (see the file header line for exactly what all the output fields are).

        For a normal dataset, print position/gene info (7 fields), total_reads, perfect_reads, N_sequence_variants, 
         the N_sequences most common sequences and counts (alternating, taking 2*N_sequences fields), 
         and gene annotation if it exists (arbitrary number of tab-separated fields).

        For a multi-dataset, print position/gene info (7 fields), the main sequence (taken over all the datasets),
         then reads_in_<x> and perfect_in_<x> for each dataset (total and perfect reads - total 2*N_datasets fields), 
         then gene annotation if it exists.

        Data is printed to OUTPUT, which should be an open file object (stdout by default).
         Output is tab-separated, with optional header starting with "# ".  
        """
        if self.multi_dataset and N_sequences is not None and N_sequences!=1:
            raise MutantError("Only one sequence can currently be printed in print_data for multi-datasets!")
        if N_sequences is None and not self.multi_dataset:     N_sequences = 2

        ### print the header line (different for normal and multi-datasets)
        # MAYBE-TODO should printing the gene info be optional?  Maybe... 
        #   Would save space when there isn't any meaningful gene info to print anyway.
        if header_line:
            header = ['chromosome','strand','min_position','full_position', 'gene','orientation','feature']
            if not self.multi_dataset:
                header += ['total_reads','perfect_reads', 'N_sequence_variants']
                for N in range(1,N_sequences+1):    
                    header += ['read_sequence_%s'%N, 'seq_%s_count'%N]
            else:
                header += ['main_sequence']
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
            mutant_data = [mutant.position.chromosome, mutant.position.strand, mutant.position.min_position, 
                           mutant.position.full_position, mutant.gene, mutant.orientation, mutant.gene_feature] 
            if not self.multi_dataset:
                mutant_data += [mutant.total_read_count, mutant.perfect_read_count, mutant.unique_sequence_count]
                for N in range(1,N_sequences+1):
                    mutant_data += list(mutant.get_main_sequence(N))
                    # MAYBE-TODO also give the length and number of mutations for each sequence? Optionally?  
                    #   Length is easy, but do I even keep track of mutation number?  I probably should...
            else:
                mutant_data += [mutant.get_main_sequence(1)[0]]
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


def read_mutant_file(infile):
    """ Read mutant input file, return as new dataset (old .txt format, or new .pickle format). """
    if infile.endswith('.pickle'):
        current_dataset = unpickle(infile)
    else:
        current_dataset = Insertional_mutant_pool_dataset(infile=infile)
        current_dataset.count_adjacent_mutants(OUTPUT=None)
    return current_dataset


############################################### Unit-tests ##############################################################

class Testing_position_functionality(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """

    from deepseq_utilities import Fake_deepseq_objects
    Fake_HTSeq_genomic_pos = Fake_deepseq_objects.Fake_HTSeq_genomic_pos

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
            self.assertRaises(MutantError, HTSeq_pos_to_tuple, self.Fake_HTSeq_genomic_pos('C', strand, start, end))
        ### testing normal functionality: should return a (chrom,start_pos,end_pos,strand) tuple with the same chromosome/strand, 
        #    and start/end positions converted from 0-based end-exclusive to 1-based end-inclusive.
        # (so the HTSeq position of AA in AATT would be 0-2, and converted would be 1-2; of TT, 2-4, and 3-4.)
        for strand in '+-':
            assert HTSeq_pos_to_tuple(self.Fake_HTSeq_genomic_pos('C', strand, 3, 7)) == ('C', 4, 7, strand)
            for (start,end) in [(0,5), (0,100), (10,12), (5,44)]:
                assert HTSeq_pos_to_tuple(self.Fake_HTSeq_genomic_pos('C', strand, start, end)) == ('C', start+1, end, strand)


    def test__get_insertion_pos_from_flanking_region_pos(self):
        # should raise exception for invalid argument (valid arguments: HTSeq position object or (chrom,start,end,strand) tuple
        #  (strand must be +/-, and start can't be after end)
        for bad_flanking_region in [None, '', 'aaa', 0, 1, 0.65, [], {}, True, False, ('C',2,3,4),('C',2,3,'x'),('C',3,2,'-')]:
            for cassette_end in SEQ_ENDS:
                self.assertRaises(MutantError, get_insertion_pos_from_flanking_region_pos, bad_flanking_region, cassette_end)
        # should raise exception for invalid cassette_end
        for bad_cassette_end in ['','aaa',0,1,[],{},None,True,False,'start','end','middle','read','leftmost','rightmost']:
            self.assertRaises(ValueError, get_insertion_pos_from_flanking_region_pos, ('C',1,5,'+'), bad_cassette_end)

        ### testing normal functionality: should return an Insertion_position instance with the same chromosome, 
        #    and strand/position depending on the arguments in a somewhat complicated way.
        # Mostly I should be testing position-tuple inputs here, since HTSeq pos object conversion to position tuples 
        #  is dealt with by HTSeq_pos_to_tuple and tested separately, but I can test some HTSeq inputs too.

        ## 1) A few spot-checks, based on the example at the end of "Possible read sides/directions" section in ../notes.txt
        #      (using HTSeq objects directly, because that's what the examples used)
        #    remember HTSeq position is 0-based and end-exclusive, and the position I want is 1-based end-inclusive!  
        #     So in the end the two relevant 1-based numbers end up being the same as the 0-based positions,
        #     because in the case of the start, min_position is actually start-1, and in the case of the end, we're adding 1
        #     to switch from 0-based to 1-based but then subtracting 1 to make the end the last base instead of the base after
        fake_HTSeq_pos_plus = self.Fake_HTSeq_genomic_pos('C', '+', 3, 7)
        fake_HTSeq_pos_minus = self.Fake_HTSeq_genomic_pos('C', '-', 3, 7)
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_plus, '5prime',reads_are_reverse=False).min_position == 7
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_plus, '3prime',reads_are_reverse=False).min_position == 3
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_minus,'5prime',reads_are_reverse=True ).min_position == 7
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_minus,'3prime',reads_are_reverse=True ).min_position == 3
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_plus, '5prime',reads_are_reverse=True ).min_position == 3
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_plus, '3prime',reads_are_reverse=True ).min_position == 7
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_minus,'5prime',reads_are_reverse=False).min_position == 3
        assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos_minus,'3prime',reads_are_reverse=False).min_position == 7

        ## 2) More checks in a loop, using position tuples and HTSeq_pos objects
        for (strand_in, if_reverse, strand_out) in [('+',False,'+'), ('-',False,'-'), ('+',True,'-'), ('-',True,'+')]:
            for (start,end) in [(1,5), (1,100), (11,12), (6,44)]:
                pos_tuple = ('C', start, end, strand_in)
                fake_HTSeq_pos = self.Fake_HTSeq_genomic_pos('C', strand_in, start-1, end)
                result_5prime = get_insertion_pos_from_flanking_region_pos(pos_tuple, '5prime', if_reverse)
                result_3prime = get_insertion_pos_from_flanking_region_pos(pos_tuple, '3prime', if_reverse)
                assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos, '5prime', if_reverse) == result_5prime
                assert get_insertion_pos_from_flanking_region_pos(fake_HTSeq_pos, '3prime', if_reverse) == result_3prime
                assert result_5prime.chromosome == result_3prime.chromosome == 'C'
                assert result_5prime.strand == result_3prime.strand == strand_out
                if strand_out=='+':
                    assert result_5prime.min_position == end
                    assert result_3prime.min_position == start-1
                elif strand_out=='-':
                    assert result_5prime.min_position == start-1
                    assert result_3prime.min_position == end


class Testing_Insertional_mutant(unittest.TestCase):
    """ Unit-tests for the Insertional_mutant class and its methods. """

    def test__init(self):
        for chromosome in ['chr1', 'chromosome_2', 'chrom3', 'a', 'adfads', '100', 'scaffold_88']:
            for strand in ['+','-']:
                for position in [1,2,5,100,10000,4323423]:
                    ins_pos_5prime = Insertion_position(chromosome,strand,position_before=position)
                    ins_pos_3prime = Insertion_position(chromosome,strand,position_after=position)
                    # test "normal" mutants - check all details, including position
                    mutant_5prime = Insertional_mutant(ins_pos_5prime)
                    mutant_3prime = Insertional_mutant(ins_pos_3prime)
                    mutant_readcount_only = Insertional_mutant_readcount_only()
                    mutant_multi_dataset = Insertional_mutant_multi_dataset(ins_pos_5prime)
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
                    # test readcount-related info for all mutants except mutant_multi_dataset
                    for mutant in [mutant_5prime, mutant_3prime, mutant_readcount_only]:
                        assert mutant.total_read_count == 0
                        assert mutant.perfect_read_count == 0
                        assert mutant.unique_sequence_count == 0
                        assert mutant.sequences_and_counts == {}
                    # test readcount-related info for mutant_multi_dataset
                    assert all([x.total_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.perfect_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.unique_sequence_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.sequences_and_counts == {} for x in mutant_multi_dataset.by_dataset.values()])

    def test__add_read(self):
        # using fake HTSeq alignment class from deepseq_utilities; defining one perfect and one imperfect alignment
        # note: the detailed mutation-counting methods are imported from deepseq_utilities and unit-tested there.
        from deepseq_utilities import Fake_deepseq_objects
        Fake_HTSeq_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment
        perfect_aln = Fake_HTSeq_alignment(seq='AAA', cigar_string='===')  # CIGAR = means match
        imperfect_aln = Fake_HTSeq_alignment(seq='GGG', cigar_string='=X=')  # CIGAR X means mismatch
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        # adding perfect and imperfect to mutant increases all the counts as expected
        mutant.add_read(perfect_aln, read_count=3)
        assert mutant.total_read_count == mutant.perfect_read_count == 3
        assert mutant.unique_sequence_count == 1
        assert mutant.sequences_and_counts == {'AAA':3}
        mutant.add_read(imperfect_aln, read_count=1)
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 3
        assert mutant.unique_sequence_count == 2
        assert mutant.sequences_and_counts == {'AAA':3, 'GGG':1}
        # it should be impossible to add a read to a specific dataset in a single-dataset mutant 
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3, dataset_name='d1')
        # same for a multi-dataset mutant - this time we need to specify which dataset we're adding to
        mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        assert len(mutant.by_dataset) == 0
        mutant.add_read(perfect_aln, read_count=3, dataset_name='d1')
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].unique_sequence_count == 1
        assert mutant.by_dataset['d1'].sequences_and_counts == {'AAA':3}
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d1')
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].unique_sequence_count == 2
        assert mutant.by_dataset['d1'].sequences_and_counts == {'AAA':3, 'GGG':1}
        # now adding a read to another dataset - nothing changes in dataset d1, but we have new dataset d2 numbers
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d2')
        assert len(mutant.by_dataset) == 2
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d2'].total_read_count == 1
        assert mutant.by_dataset['d2'].perfect_read_count == 0
        assert mutant.by_dataset['d2'].unique_sequence_count == 1
        assert mutant.by_dataset['d2'].sequences_and_counts == {'GGG':1}
        # it should be impossible to add a read to a multi-dataset mutant without giving a dataset_name
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3)

    def test__update_gene_info(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == '?'
        # updating no-info mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?')
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == '?'
        # updating no-info mutant with useful info - update goes through
        mutant.update_gene_info('gene1', '+', 'f')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containing mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containing mutant with same info - no change
        mutant.update_gene_info('gene1', '+', 'f')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containig mutant with OTHER info - exception
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene2', '+', 'f')
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '-', 'f')
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '+', 'g')


    def test__merge_mutant(self):
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant1.add_counts(2,2,1)
        mutant1.add_sequence_and_counts('AAA',2)
        mutant2 = Insertional_mutant(Insertion_position('chr','+',position_before=5))
        mutant2.add_counts(1,0,1)
        mutant2.add_sequence_and_counts('AA',1)
        # merging two mutants - read counts act as expected
        assert (mutant1.total_read_count, mutant1.perfect_read_count) == (2, 2)
        assert (mutant2.total_read_count, mutant2.perfect_read_count) == (1, 0)
        assert (mutant1.sequences_and_counts, mutant2.sequences_and_counts) == ({'AAA':2}, {'AA':1})
        mutant1.merge_mutant(mutant2)
        assert (mutant1.total_read_count, mutant1.perfect_read_count) == (3, 2)
        assert (mutant2.total_read_count, mutant2.perfect_read_count) == (0, 0)
        assert (mutant1.sequences_and_counts, mutant2.sequences_and_counts) == ({'AAA':2, 'AA':1}, {})
        # if we define the gene for one of the mutants but not the other, still works
        #  (note that starting here mutant2 has 0 reads, but we can still merge it all we want just to check for errors)
        mutant2.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True)
        # if we define the gene for both mutants and it's the same gene, also works
        mutant1.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True)
        # if we define the gene for both mutants and it's a DIFFERENT gene, we get an exception unless we don't check
        mutant2.gene = 'Y'
        self.assertRaises(MutantError, mutant1.merge_mutant, mutant2, check_gene_data=True)
        mutant1.merge_mutant(mutant2, check_gene_data=False)
        # mutant-merging for multi-dataset mutants currently NOT IMPLEMENTED

    def test__add_counts(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0)
        assert mutant.total_read_count == 0
        assert mutant.perfect_read_count == 0
        assert mutant.sequences_and_counts == {}
        mutant.add_counts(2,2,1)
        assert mutant.total_read_count == 2
        assert mutant.perfect_read_count == 2
        assert mutant.sequences_and_counts == {}
        mutant.add_counts(2,1,1)
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 3
        assert mutant.sequences_and_counts == {}
        # how mutant.unique_sequence_count changes depends on assume_new_sequences:
        #  - if False, mutant.unique_sequence_count is max(new_seq_count, mutant.unique_sequence_count)
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 0
        mutant.add_counts(1,1,1,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 1
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        mutant.add_counts(1,1,1,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        #  - if True and mutant.unique_sequence_count is new_seq_count + mutant.unique_sequence_count
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 0
        mutant.add_counts(1,1,1,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 1
        mutant.add_counts(2,2,2,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 3
        mutant.add_counts(1,1,1,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 4
        mutant.add_counts(2,2,2,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 6
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 6
        ### works for multi-dataset mutants (not looking at too much detail - should work the same as for single)
        multi_mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        multi_mutant.add_counts(2,2,1, dataset_name='d')
        assert multi_mutant.by_dataset['d'].total_read_count == 2
        assert multi_mutant.by_dataset['d'].perfect_read_count == 2
        assert multi_mutant.by_dataset['d'].unique_sequence_count == 1
        # doesn't work for multi-dataset mutants if dataset_name not given, or for single if given
        self.assertRaises(MutantError, multi_mutant.add_counts, 1,1,1)
        self.assertRaises(MutantError, mutant.add_counts, 1,1,1, dataset_name='d1')

    def test__add_sequence_and_counts(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        # adding sequence/count to mutant.sequences_and_counts, WITHOUT touching mutant.unique_sequence_count
        mutant.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=False)
        assert mutant.sequences_and_counts == {'AAA':2}
        assert mutant.unique_sequence_count == 0
        # adding sequence/count to mutant.sequences_and_counts, and incrementing mutant.unique_sequence_count if warranted:
        #  - if adding a sequence that was already there, don't increment
        mutant.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4}
        assert mutant.unique_sequence_count == 0
        #  - if adding a new sequence, increment
        mutant.add_sequence_and_counts('GGG',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4, 'GGG':2}
        assert mutant.unique_sequence_count == 1
        #  - if adding a new sequence but mutant.unique_sequence_count is already higher than expected, don't increment
        mutant.unique_sequence_count = 5
        mutant.add_sequence_and_counts('CCC',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4, 'GGG':2, 'CCC':2}
        assert mutant.unique_sequence_count == 5
        # make sure it raises an error if given a non-numeric argument
        for not_a_number in [None,'','a','GGG']:
            self.assertRaises(TypeError,mutant.add_sequence_and_counts,'CCC',not_a_number)
        ### works for multi-dataset mutants (not looking at too much detail - should work the same as for single)
        multi_mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        assert multi_mutant.by_dataset['d'].sequences_and_counts == {}
        multi_mutant.add_sequence_and_counts('GGG',2, dataset_name='d')
        assert multi_mutant.by_dataset['d'].sequences_and_counts == {'GGG':2}
        # doesn't work for multi-dataset mutants if dataset_name not given, or for single if given
        self.assertRaises(MutantError, multi_mutant.add_sequence_and_counts, 'GGG',1)
        self.assertRaises(MutantError, mutant.add_sequence_and_counts, 'GGG',1, dataset_name='d1')

    def test__get_main_sequence(self):
        # single-dataset mutant
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        assert mutant.get_main_sequence() == ('',0)
        assert mutant.get_main_sequence(1) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_sequence_and_counts('AAA',1)
        mutant.add_sequence_and_counts('GGG',2)
        assert mutant.sequences_and_counts == {'AAA':1, 'GGG':2}
        assert mutant.get_main_sequence() == ('GGG',2)
        assert mutant.get_main_sequence(1) == ('GGG',2)
        assert mutant.get_main_sequence(2) == ('AAA',1)
        assert mutant.get_main_sequence(3) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_sequence_and_counts('CCC',1)
        mutant.add_sequence_and_counts('AAA',2)
        assert mutant.sequences_and_counts == {'AAA':3, 'GGG':2, 'CCC':1}
        assert mutant.get_main_sequence() == ('AAA',3)
        assert mutant.get_main_sequence(1) == ('AAA',3)
        assert mutant.get_main_sequence(2) == ('GGG',2)
        assert mutant.get_main_sequence(3) == ('CCC',1)
        assert mutant.get_main_sequence(4) == ('',0)
        assert mutant.get_main_sequence(5) == ('',0)
        # multi-dataset mutant - getting the top sequence from single dataset or from all datasets added together
        mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        mutant.add_sequence_and_counts('AAA',3, dataset_name='d1')
        mutant.add_sequence_and_counts('GGG',2, dataset_name='d1')
        mutant.add_sequence_and_counts('TTT',3, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',2, dataset_name='d2')
        assert mutant.get_main_sequence(1, dataset_name='d1') == ('AAA',3)
        assert mutant.get_main_sequence(1, dataset_name='d2') == ('TTT',3)
        assert mutant.get_main_sequence(1) == ('GGG',4)     # GGG is the most common sequence if we add both datasets

    def test__add_other_mutant_as_dataset(self):
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant1.add_counts(2,1,1)
        mutant1.add_sequence_and_counts('AAA',2)
        mutant2 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant2.add_counts(1,0,1)
        mutant2.add_sequence_and_counts('GGG',1)
        # adding a mutant to a single-dataset mutant should fail
        if hasattr(mutant1, 'add_other_mutant_as_dataset'):
            self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
        # adding a mutant to a multi-dataset mutant should work
        mutantM = Insertional_mutant_multi_dataset(mutant1.position)
        mutantM.add_other_mutant_as_dataset(mutant1, 'd1')
        mutantM.add_other_mutant_as_dataset(mutant2, 'd2')
        assert mutantM.by_dataset['d1'].total_read_count == 2
        assert mutantM.by_dataset['d1'].perfect_read_count == 1
        assert mutantM.by_dataset['d1'].unique_sequence_count == 1
        assert mutantM.by_dataset['d1'].sequences_and_counts == {'AAA':2}
        assert mutantM.by_dataset['d2'].total_read_count == 1
        assert mutantM.by_dataset['d2'].perfect_read_count == 0
        assert mutantM.by_dataset['d2'].unique_sequence_count == 1
        assert mutantM.by_dataset['d2'].sequences_and_counts == {'GGG':1}
        # can't add overwrite an existing dataset name, unless overwrite==True
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant2, 'd2')
        mutantM.add_other_mutant_as_dataset(mutant2, 'd2', overwrite=True)
        # if the two mutants have different positions, it should fail, unless check_constant_data=False
        mutant3 = Insertional_mutant(Insertion_position('chr','+',position_before=5))
        mutant4 = Insertional_mutant(Insertion_position('chr2','+',position_before=3))
        mutant5 = Insertional_mutant(Insertion_position('chr','-',position_before=3))
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant3, 'd3', check_constant_data=True)
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant4, 'd4', check_constant_data=True)
        self.assertRaises(MutantError, mutantM.add_other_mutant_as_dataset, mutant5, 'd5', check_constant_data=True)
        mutantM.add_other_mutant_as_dataset(mutant3, 'd3', check_constant_data=False)
        mutantM.add_other_mutant_as_dataset(mutant4, 'd4', check_constant_data=False)
        mutantM.add_other_mutant_as_dataset(mutant5, 'd5', check_constant_data=False)

    def test__give_single_dataset_mutant(self):
        mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(2,1,1, dataset_name='d1')
        mutant.add_sequence_and_counts('AAA',2, dataset_name='d1')
        mutant.add_counts(1,0,1, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',1, dataset_name='d2')
        # extracting two mutants and checking the values
        new_mutant_1 = mutant.give_single_dataset_mutant('d1')
        new_mutant_2 = mutant.give_single_dataset_mutant('d2')
        for new_mutant in (new_mutant_1, new_mutant_2):
            assert new_mutant.position == mutant.position
            assert new_mutant.gene == mutant.gene
        assert new_mutant_1.total_read_count == 2
        assert new_mutant_1.perfect_read_count == 1
        assert new_mutant_1.unique_sequence_count == 1
        assert new_mutant_1.sequences_and_counts == {'AAA':2}
        assert new_mutant_2.total_read_count == 1
        assert new_mutant_2.perfect_read_count == 0
        assert new_mutant_2.unique_sequence_count == 1
        assert new_mutant_2.sequences_and_counts == {'GGG':1}
        # trying to extract an inexistent dataset should fail, unless force==True
        self.assertRaises(MutantError, mutant.give_single_dataset_mutant, 'd0', force=False)
        new_mutant_0 = mutant.give_single_dataset_mutant('d0', force=True)
        assert new_mutant_0.position == mutant.position
        assert new_mutant_0.gene == mutant.gene
        assert new_mutant_0.total_read_count == 0
        assert new_mutant_0.perfect_read_count == 0
        assert new_mutant_0.unique_sequence_count == 0
        assert new_mutant_0.sequences_and_counts == {}
        # trying to extract a dataset from a single-dataset mutant should fail
        mutant1 = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        self.assertRaises(MutantError, mutant1.give_single_dataset_mutant, 'd2')

    def test__give_all_single_dataset_mutants(self):
        mutant = Insertional_mutant_multi_dataset(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(2,1,1, dataset_name='d1')
        mutant.add_sequence_and_counts('AAA',2, dataset_name='d1')
        mutant.add_counts(1,0,1, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',1, dataset_name='d2')
        # extracting two mutants and checking the values
        all_single_mutants = mutant.give_all_single_dataset_mutants()
        assert len(all_single_mutants) == 2
        new_mutant_1 = all_single_mutants['d1']
        new_mutant_2 = all_single_mutants['d2']
        for new_mutant in (new_mutant_1, new_mutant_2):
            assert new_mutant.position == mutant.position
            assert new_mutant.gene == mutant.gene
        assert new_mutant_1.total_read_count == 2
        assert new_mutant_1.perfect_read_count == 1
        assert new_mutant_1.unique_sequence_count == 1
        assert new_mutant_1.sequences_and_counts == {'AAA':2}
        assert new_mutant_2.total_read_count == 1
        assert new_mutant_2.perfect_read_count == 0
        assert new_mutant_2.unique_sequence_count == 1
        assert new_mutant_2.sequences_and_counts == {'GGG':1}


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
        for string in positions_and_readcounts_string.split(', '):
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
            mutant = Insertional_mutant(insertion_position=full_pos)
            mutant.total_read_count = readcount
            mutant.perfect_read_count = perfect
            dataset.add_mutant(mutant)
        return dataset

    def test__init(self):
        for cassette_end in SEQ_ENDS+['?']:
            for reads_are_reverse in [True,False,'?']:
                data = Insertional_mutant_pool_dataset(cassette_end=cassette_end, reads_are_reverse=reads_are_reverse)
                assert data.summary.cassette_end == cassette_end
                assert data.summary.reads_are_reverse == reads_are_reverse
                assert len(data) == 0
                assert data.summary.processed_read_count == data.summary.aligned_read_count\
                        == data.summary.perfect_read_count == 0
                assert data.summary.non_aligned_read_count == 0
                assert data.summary.discarded_read_count == None
                assert data.summary.ignored_region_read_counts == {}
                assert sum(data.summary.strand_read_counts.values()) == 0
                assert data.summary.mutants_in_genes == data.summary.mutants_not_in_genes\
                        == data.summary.mutants_undetermined == 0
                assert data.summary.mutant_counts_by_orientation == data.summary.mutant_counts_by_feature == {}
        for cassette_end in [True, False, None, 0, 1, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, cassette_end, '?')
        for reads_are_reverse in ['forward', 'reverse', None, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            # note that this list doesn't include 0/1 because 0==False and 1==True, and that's what "in" tests 
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, '?', reads_are_reverse)

    # LATER-TODO add unit-test for add_discarded_reads, find_genes_for_mutants, most_common_mutants, 

    def test__add_alignment_reader_to_data(self):
        pass
        # MAYBE-TODO implement using a mock-up of HTSeq_alignment?  (see Testing_single_functions for how I did that)
        #   make sure it fails if self.cassette_end isn't defined...

    def test__merge_other_dataset(self):
        ### testing that the various overall information gets merged correctly
        dataset = Insertional_mutant_pool_dataset()
        other_dataset = Insertional_mutant_pool_dataset()
        for name1,name2 in [('TEST1', 'TEST2'), ('TEST', 'TEST'), ('TEST', None), (None, 'TEST'), (None, None)]: 
            dataset.summary.dataset_name = name1
            other_dataset.summary.dataset_name = name2
            dataset.merge_other_dataset(other_dataset)
            assert dataset.summary.dataset_name is None
        for end_info in ['3prime 3prime 3prime', '5prime 5prime 5prime', '3prime 5prime ?', '3prime ? ?', '5prime ? ?', '? ? ?']:
            end_a, end_b, end_both = end_info.split()
            for end1,end2 in [(end_a,end_b),(end_b,end_a)]:
                dataset.summary.cassette_end = end1
                other_dataset.summary.cassette_end = end2
                dataset.merge_other_dataset(other_dataset)
                assert dataset.summary.cassette_end == end_both
        for reverse_a,reverse_b,reverse_both in [(True,True,True),(False,False,False), (True,False,'?'), 
                                                 (True,'?','?'), (False,'?','?'), ('?','?','?')]:
            for reverse1,reverse2 in [(reverse_a,reverse_b),(reverse_b,reverse_a)]:
                dataset.summary.reads_are_reverse = reverse1
                other_dataset.summary.reads_are_reverse = reverse2
                dataset.merge_other_dataset(other_dataset)
                assert dataset.summary.reads_are_reverse == reverse_both
        ### testing that the extra readcount information gets merged correctly
        # if both datasets have numbers, they get added together
        dataset.summary.add_discarded_reads(60, 30, 20, 10, replace=True)
        other_dataset.summary.add_discarded_reads(6, 3, 2, 1, replace=True)
        dataset.summary.add_nonaligned_reads(30, 20, 10, replace=True)
        other_dataset.summary.add_nonaligned_reads(3, 2, 1, replace=True)
        dataset.summary.ignored_region_read_counts = {'test1':10, 'test2':20}
        other_dataset.summary.ignored_region_read_counts = {'test1':1, 'test3':3}
        dataset.merge_other_dataset(other_dataset)
        assert (dataset.summary.discarded_read_count, dataset.summary.discarded_wrong_start, 
                dataset.summary.discarded_no_cassette, dataset.summary.discarded_other_end) == (66, 33, 22, 11)
        assert (dataset.summary.non_aligned_read_count, dataset.summary.unaligned, 
                dataset.summary.multiple_aligned) == (33, 22, 11)
        assert dataset.summary.ignored_region_read_counts == {'test1':11, 'test2':20, 'test3':3}
        # if one dataset has numbers and the other has unknowns, the result should be unknowns either way, 
        #  except summary.ignored_region_read_counts, which should still just get added)
        for version in (1,2):
            dataset.summary.add_discarded_reads(60, 30, 20, 10, replace=True)
            other_dataset.summary.add_discarded_reads(None, None, None, None, replace=True)
            dataset.summary.add_nonaligned_reads(30, 20, 10, replace=True)
            other_dataset.summary.add_nonaligned_reads(None, None, None, replace=True)
            dataset.summary.ignored_region_read_counts = {'test1':10, 'test2':20}
            other_dataset.summary.ignored_region_read_counts = {}
            if version==1:
                dataset.merge_other_dataset(other_dataset)
            else:
                other_dataset.merge_other_dataset(dataset)
                dataset = other_dataset
            assert (dataset.summary.discarded_read_count, dataset.summary.discarded_wrong_start, 
                    dataset.summary.discarded_no_cassette, dataset.summary.discarded_other_end) == ('unknown',) * 4
            assert (dataset.summary.non_aligned_read_count, dataset.summary.unaligned, 
                    dataset.summary.multiple_aligned) == ('unknown', 'unknown', 'unknown')
            assert dataset.summary.ignored_region_read_counts == {'test1':10, 'test2':20}
        ### testing that the mutant merging/counting gets dealt with correctly
        # checking that all that data is made blank, if it's non-blank in either dataset
        # these involve no merging
        positions_and_readcounts_raw_1 = "A+100 1, B-201 1, C-300 1, E+500 1, D-400 1"
        positions_and_readcounts_raw_2 = "A+100 1, A-101 1, A+301 1, A-300 1, A-400 1"
        positions_and_readcounts_raw_3 = ""
        all_raws = [positions_and_readcounts_raw_1, positions_and_readcounts_raw_2, positions_and_readcounts_raw_3]
        for positions_and_readcounts_raw_a,positions_and_readcounts_raw_b in itertools.product(all_raws, all_raws):
            # do merging-and-counting for one or both datasets
            for merge_1,merge_2 in [(True,True), (True,False), (False,True)]:
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_a)
                other_dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_b)
                if merge_1:
                    dataset.merge_adjacent_mutants(OUTPUT=None)
                    dataset.merge_opposite_tandem_mutants(OUTPUT=None)
                    dataset.count_adjacent_mutants()
                if merge_2:
                    other_dataset.merge_adjacent_mutants(OUTPUT=None)
                    other_dataset.merge_opposite_tandem_mutants(OUTPUT=None)
                    other_dataset.count_adjacent_mutants()
                dataset.merge_other_dataset(other_dataset)
                assert dataset.summary.adjacent_max_distance == None
                assert dataset.summary.merging_which_chromosomes == (None, None)
                assert dataset.summary.merged_adjacent_same_strand_dict == defaultdict(int)
                assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == defaultdict(list)
                assert dataset.summary.merged_opposite_tandems == 0
                assert dataset.summary.merged_opposite_tandems_readcounts == []
                assert isnan(dataset.summary.same_position_opposite)
                assert dataset.summary.adjacent_opposite_toward_dict == defaultdict(nan_func)
                assert dataset.summary.adjacent_opposite_away_dict == defaultdict(nan_func)
                assert dataset.summary.adjacent_same_strand_dict == defaultdict(nan_func)
        # checking that an exception is raised if either dataset had mutant-merging done
        # 1 and 2 have merging, 3 doesn't
        positions_and_readcounts_raw_1 = "A+100 5, A+101 5, B+100 5, B-100 5, C+100 5, C-101 1, D-100 1, D+101 1"
        positions_and_readcounts_raw_2 = "A+100 1, A-100 1, A+301 1, A-300 1, A-400 1"
        positions_and_readcounts_raw_3 = "A+100 1, A-101 1, A+301 1, A-300 1, A-400 1"
        all_raws = [positions_and_readcounts_raw_1, positions_and_readcounts_raw_2, positions_and_readcounts_raw_3]
        for positions_and_readcounts_raw_a,positions_and_readcounts_raw_b in itertools.permutations(all_raws, 2):
            # do merging-and-counting for one or both datasets
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_a)
                other_dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_b)
                dataset.merge_adjacent_mutants(leave_N_mutants=0, OUTPUT=None)
                dataset.merge_opposite_tandem_mutants(leave_N_mutants=0, OUTPUT=None)
                other_dataset.merge_adjacent_mutants(leave_N_mutants=0, OUTPUT=None)
                other_dataset.merge_opposite_tandem_mutants(leave_N_mutants=0, OUTPUT=None)
                self.assertRaises(MutantError, dataset.merge_other_dataset, other_dataset)
        ### basic test that the correct mutants are being merged:  merge identical positions (first pair), 
        #  but not one-off positions, opposite strands, different positions or different chromosomes (remaining pairs)
        positions_and_readcounts_raw_1 = "A+100 5, B-200 5, C+300 5, D+500 5, E-400 5"
        positions_and_readcounts_raw_2 = "A+100 1, B-201 1, C-300 1, E+500 1, D-400 1"
        chromosomes = set([s[0] for s in positions_and_readcounts_raw_1.split(', ')+positions_and_readcounts_raw_2.split(', ')])
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_1, raw_chrom_names=True)
        other_dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw_2, raw_chrom_names=True)
        for chrom in chromosomes: 
            assert dataset.summary.mutants_in_chromosome(chrom) == other_dataset.summary.mutants_in_chromosome(chrom) == 1
            assert (dataset.summary.reads_in_chromosome(chrom), other_dataset.summary.reads_in_chromosome(chrom)) == (5, 1)
        dataset.merge_other_dataset(other_dataset)
        for chrom,N_mutants in zip(chromosomes,(1,2,2,2,2)):
            assert (dataset.summary.mutants_in_chromosome(chrom),other_dataset.summary.mutants_in_chromosome(chrom)) == (N_mutants,0)
            assert (dataset.summary.reads_in_chromosome(chrom), other_dataset.summary.reads_in_chromosome(chrom)) == (6, 0)
        ### checking that an exception is raised if attempting to merge a both-strand mutant (in either dataset)
        mutant_both = Insertional_mutant(Insertion_position('A','both',position_before=100, immutable=True))
        mutant_both.add_counts(2,2,1)
        mutant_both.add_sequence_and_counts('AAA',2)
        other_dataset.add_mutant(mutant_both)
        self.assertRaises(MutantError, dataset.merge_other_dataset, other_dataset)
        self.assertRaises(MutantError, other_dataset.merge_other_dataset, dataset)
        # MAYBE-TODO add more detailed checks to see if the mutant sequences are merged correctly etc?  Or should that just be done in the mutant.merge_mutant unit-test?  It looks good enough.

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

    def test___readcounts_sorted(self):
        # if the two values are readcounts
        dataset = Insertional_mutant_pool_dataset()
        assert dataset._readcounts_sorted(1,1) == (1,1)
        assert dataset._readcounts_sorted(1,2) == dataset._readcounts_sorted(2,1) == (2,1)
        assert dataset._readcounts_sorted(5,15) == dataset._readcounts_sorted(15,5) == (15,5)
        assert dataset._readcounts_sorted(4,6) == dataset._readcounts_sorted(6,4) == (6,4)
        # if the two values are mutants (these are fake mutants with no attributes except total_read_count)
        mutant1, mutant2 = Insertional_mutant(), Insertional_mutant()
        mutant1.total_read_count, mutant2.total_read_count = 1, 2
        assert dataset._readcounts_sorted(mutant1, mutant1) == (1,1)
        assert dataset._readcounts_sorted(mutant2, mutant2) == (2,2)
        assert dataset._readcounts_sorted(mutant1, mutant2) == dataset._readcounts_sorted(mutant2, mutant1) == (2,1)
        # if the two values are positions (now the fake mutants must have positions too, and must be in a dataset)
        pos1 = Insertion_position('chr1', '+', position_before=100, immutable=True)
        pos2 = Insertion_position('chr1', '+', position_before=200, immutable=True)
        mutant1.position, mutant2.position = pos1, pos2
        dataset.add_mutant(mutant1)
        dataset.add_mutant(mutant2)
        assert dataset._readcounts_sorted(pos1, pos1) == (1,1)
        dataset._readcounts_sorted(pos2, pos2) == (2,2)
        assert dataset._readcounts_sorted(pos1, pos2) == dataset._readcounts_sorted(pos2, pos1) == (2,1)

    def test___readcount_ratio(self):
        # if the two values are readcounts
        dataset = Insertional_mutant_pool_dataset()
        assert dataset._readcount_ratio(1,1) == 1
        assert dataset._readcount_ratio(1,2) == dataset._readcount_ratio(2,1) == 2
        assert dataset._readcount_ratio(5,15) == dataset._readcount_ratio(15,5) == 3
        assert dataset._readcount_ratio(4,6) == dataset._readcount_ratio(6,4) == 1.5
        # if the two values are mutants (these are fake mutants with no attributes except total_read_count)
        mutant1, mutant2 = Insertional_mutant(), Insertional_mutant()
        mutant1.total_read_count, mutant2.total_read_count = 1, 2
        assert dataset._readcount_ratio(mutant1, mutant1) == dataset._readcount_ratio(mutant2, mutant2) == 1
        assert dataset._readcount_ratio(mutant1, mutant2) == dataset._readcount_ratio(mutant2, mutant1) == 2
        # if the two values are positions (now the fake mutants must have positions too, and must be in a dataset)
        pos1 = Insertion_position('chr1', '+', position_before=100, immutable=True)
        pos2 = Insertion_position('chr1', '+', position_before=200, immutable=True)
        mutant1.position, mutant2.position = pos1, pos2
        dataset.add_mutant(mutant1)
        dataset.add_mutant(mutant2)
        assert dataset._readcount_ratio(pos1, pos1) == dataset._readcount_ratio(pos2, pos2) == 1
        assert dataset._readcount_ratio(pos1, pos2) == dataset._readcount_ratio(pos2, pos1) == 2

    def test__merge_adjacent_mutants(self):
        ### basic tests
        # making a dataset with several different same-strand adjacent pairs (alternating strand between pairs):
        #   1*100 - dist 1 ratio 5
        #   1*200 - dist 1 ratio 5
        #   1*300 - dist 1 ratio 2
        #   1*400 - dist 1 ratio 1
        #   2*100 - dist 2 ratio 5
        #   5*100 - dist 5 ratio 5
        #   5*200 - dist 5 ratio 2
        positions_and_readcounts_raw = ("1+100 5, 1+101 1, 1-200 1, 1-201 5, 1-300 4, 1-301 2, 1+400 1, 1+401 1, "
                                        "2-100 5, 2-102 1, "
                                        "5+100 5, 5+105 1, 5+200 4, 5+205 2")
        # adding negative cases - NO same-strand adjacent pairs, only different-chromosome or different-strand ones:
        positions_and_readcounts_raw += ", A+100 5, B+101 1, C+200 5, C-200 1, C+300 5, C-301 1"
        # test1 - max distance 1, min ratio 5
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dataset.summary.merged_adjacent_same_strand_dict == {1:2}
        assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(5,1),(5,1)]}
        # test2 - max distance 2, min ratio 5
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=2, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dataset.summary.merged_adjacent_same_strand_dict == {1:2, 2:1}
        assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(5,1),(5,1)], 2:[(5,1)]}
        # test3 - max distance 10, min ratio 5
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=10, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dataset.summary.merged_adjacent_same_strand_dict == {1:2, 2:1, 5:1}
        assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(5,1),(5,1)], 2:[(5,1)], 5:[(5,1)]}
        # test4 - max distance 10, min ratio 2
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=10, min_count_ratio=2, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dataset.summary.merged_adjacent_same_strand_dict == {1:3, 2:1, 5:2}
        assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(4,2),(5,1),(5,1)], 
                                                                               2:[(5,1)], 5:[(4,2),(5,1)]}
        # test5 - max distance 10, min ratio 1
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=10, min_count_ratio=1, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dataset.summary.merged_adjacent_same_strand_dict == {1:4, 2:1, 5:2}
        assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(1,1),(4,2),(5,1),(5,1)], 
                                                                               2:[(5,1)], 5:[(4,2),(5,1)]}
        ### tests for the new leave_N_mutants and leave_method arguments
        # making a dataset with several different same-strand adjacent pairs (alternating strand between pairs):
        #   1*100 - dist 1 ratio 5
        #   1*200 - dist 1 ratio 5
        #   1*300 - dist 1 ratio 2
        #   1*400 - dist 1 ratio 1
        positions_and_readcounts_raw = ("1+100 5, 1+101 1, 1-200 1, 1-201 5, 1-300 4, 1-301 2, 1+400 1, 1+401 1")
        # when leave_N_mutants isn't 'use_ratio', the ratio is ignored, so repeat the test for a few different numbers
        for R in (1,2,10,100, None):
            # merge everything, regardless of the ratio given
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=0, min_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_adjacent_same_strand_dict == {1:4}
            # don't merge anything
            for N in (4,5,10,100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=N, min_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_adjacent_same_strand_dict == {}
            # leave 1 (merge 3), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=1, leave_method='by_ratio', 
                                           min_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_adjacent_same_strand_dict == {1:3}
            assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(4,2),(5,1),(5,1)]}
            # leave 2 (merge 2), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=2, leave_method='by_ratio', 
                                           min_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_adjacent_same_strand_dict == {1:2}
            assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(5,1),(5,1)]}
            # leave 3 (merge 1), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=3, leave_method='by_ratio', 
                                           min_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_adjacent_same_strand_dict == {1:1}
            assert dataset.summary.merged_adjacent_same_strand_readcounts_dict == {1:[(5,1)]}
            # leave 1 (merge 3), random
            readcount_pairs = set()
            for _ in range(100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=1, leave_method='random', 
                                               min_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_adjacent_same_strand_dict == {1:3}
                for pair in dataset.summary.merged_adjacent_same_strand_readcounts_dict[1]: readcount_pairs.add(pair)
            assert readcount_pairs == set([(1,1),(4,2),(5,1),(5,1)])
            # leave 3 (merge 1), random
            readcount_pairs = set()
            for _ in range(100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants=3, leave_method='random', 
                                               min_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_adjacent_same_strand_dict == {1:1}
                for pair in dataset.summary.merged_adjacent_same_strand_readcounts_dict[1]: readcount_pairs.add(pair)
            assert readcount_pairs == set([(1,1),(4,2),(5,1),(5,1)])

        ### leave_N_mutants='auto' tests
        # the resulting leave_N_mutants value should be min(away+toward, toward*2)
        # (+ before - is away, - before + is toward)
        extra_pairs_and_expected_merged = [
            # no opposite-strand adjacent cases - leave none, merge all
            (''                                     + '',                                       4),
            # 0 away case, 1 toward - leave 1, merge 3
            (''                                     + ', 1-550 1, 1+551 1',                     3),
            # 1 away case, 0 toward - leave 0, merge all
            (', 1+500 1, 1-501 1'                   + '',                                       4),
            # 1 away case, 1 toward - leave 2, merge 2
            (', 1+500 1, 1-501 1'                   + ', 1-550 1, 1+551 1',                     2),
            # 2 away case, 1 toward - leave 2, merge 2
            (', 1+500 1, 1-501 1, 1+600 1, 1-601 1' + ', 1-550 1, 1+551 1',                     2),
            # 1 away case, 2 toward - leave 3, merge 1
            (', 1+500 1, 1-501 1'                   + ', 1-550 1, 1+551 1, 1-650 1, 1+651 1',   1),
            # 2 away case, 2 toward - leave 4, merge 0
            (', 1+500 1, 1-501 1, 1+600 1, 1-601 1' + ', 1-550 1, 1+551 1, 1-650 1, 1+651 1',   0),
        ]
        for (extra_pairs, expected_merged) in extra_pairs_and_expected_merged:
            for R in (1,2,10,100, None):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw + extra_pairs)
                dataset.merge_adjacent_mutants(merge_max_distance=1, leave_N_mutants='auto', leave_method='by_ratio', 
                                               min_count_ratio=R, OUTPUT=None)
                if expected_merged:     assert dataset.summary.merged_adjacent_same_strand_dict == {1:expected_merged}
                else:                   assert dataset.summary.merged_adjacent_same_strand_dict == {}

        ### tests with more complicated and rare cases like three adjacent mutants
        # transitive merging: make sure that A merges into C if A should merge into B and B into C, 
        positions_and_readcounts_raw = "1+100 1, 1+101 2, 1+102 4"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=2, leave_N_mutants='use_ratio', OUTPUT=None)
        assert [m.position.min_position for m in dataset] == [102]
        # transitive merging addendum - check that in the same case, A wouldn't merge into C without B there
        positions_and_readcounts_raw = "1+100 1, 1+102 4"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=2, leave_N_mutants='use_ratio', OUTPUT=None)
        assert len(dataset) == 2
        # another transitive merging case, with opposite merge direction and higher distances
        positions_and_readcounts_raw = "1+100 4, 1+102 2, 1+105 1"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=4, min_count_ratio=2, leave_N_mutants='use_ratio', OUTPUT=None)
        assert [m.position.min_position for m in dataset] == [100]
        # another transitive merging case, with opposite strand and a longer chain
        positions_and_readcounts_raw = "1-100 8, 1-101 4, 1-102 2, 1-103 1"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=2, leave_N_mutants='use_ratio', OUTPUT=None)
        assert [m.position.min_position for m in dataset] == [100]
        # split merging 1 - if B could be merged into either A or C, pick the closer one
        positions_and_readcounts_raw = "1+99 5, 1+100 1, 1+102 10"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=2, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dict([(m.position.min_position,m.total_read_count) for m in dataset]) == {99:6, 102:10}
        # split merging 2 - if B could be merged into either A or C and distances are the same, pick the higher-read one
        positions_and_readcounts_raw = "1+99 5, 1+100 1, 1+101 10"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=2, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dict([(m.position.min_position,m.total_read_count) for m in dataset]) == {99:5, 101:11}
        # split merging 3 - if B could be merged into either A or C and distances and readcounts are the same, 
        #                       pick the earlier-position one (arbitrary but nonrandom for ease of testing)
        positions_and_readcounts_raw = "1+99 5, 1+100 1, 1+101 5"
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_adjacent_mutants(merge_max_distance=2, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None)
        assert dict([(m.position.min_position,m.total_read_count) for m in dataset]) == {99:6, 101:5}

        ### tests for tests for merge_cassette_chromosomes and merge_other_chromosomes args
        positions_and_readcounts_raw = ("chromosome_1+100 5, chromosome_1+101 1, "
                                        "insertion_cassette+100 5, insertion_cassette+101 1, "
                                        "something_else+100 5, something_else+101 1")
        # merge all
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None, 
                                       merge_cassette_chromosomes=True, merge_other_chromosomes=True)
        assert len(dataset) == 3
        # don't merge cassette
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None, 
                                       merge_cassette_chromosomes=False, merge_other_chromosomes=True)
        assert len(dataset) == 4
        assert dataset.summary.mutants_in_chromosome('insertion_cassette') == 2
        # don't merge other
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None, 
                                       merge_cassette_chromosomes=True, merge_other_chromosomes=False)
        assert len(dataset) == 4
        assert dataset.summary.mutants_in_chromosome('something_else') == 2
        # don't merge cassette or other
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_adjacent_mutants(merge_max_distance=1, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None, 
                                       merge_cassette_chromosomes=False, merge_other_chromosomes=False)
        assert len(dataset) == 5
        assert dataset.summary.mutants_in_chromosome('chromosome_1') == 1
        ### tests for bad arguments
        # multi-dataset
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.multi_dataset = True
        self.assertRaises(MutantError, dataset.merge_adjacent_mutants, OUTPUT=None)
        # ratio<1
        for ratio in (0, 0.1, 0.5, 0.99):
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
            self.assertRaises(MutantError, dataset.merge_adjacent_mutants, min_count_ratio=ratio, OUTPUT=None)
        # max_distance inconsistent with previous setting in dataset.summary
        distances = set(range(10))
        for dist1 in distances:
            for dist2 in distances-set([dist1]):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
                dataset.summary.adjacent_max_distance = dist1
                self.assertRaises(MutantError, dataset.merge_adjacent_mutants, merge_max_distance=dist2, OUTPUT=None)
        # any both-strand mutants present
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.add_mutant(Insertional_mutant(insertion_position=Insertion_position('chromosome_1', 'both', 
                                                                                    position_before=100, immutable=True)))
        self.assertRaises(MutantError, dataset.merge_adjacent_mutants, OUTPUT=None)
        # settings for whether to merge cassette/other chromosomes inconsistent with dataset.summary settings
        for other1,cassette1 in [(True,True), (True,False), (False,True), (False,False)]:
            dataset.summary.merging_which_chromosomes = (cassette1,other1)
            for other2,cassette2 in [(other1,not cassette1), (not other1, cassette1), (not other1, not cassette1)]:
                self.assertRaises(MutantError, dataset.merge_adjacent_mutants, OUTPUT=None, 
                                 merge_cassette_chromosomes=cassette1, merge_other_chromosomes=other1)

    def test__merge_opposite_tandem_mutants(self):
        ### setup
        # four same-position opposite-strand cases, with ratios 5, 2, 1, 1 (5:1, 5:5, 4:2, 1:1)
        positions_and_readcounts_raw = "A+100 5, A-100 1, A+200 5, A-200 5, A+300 4, A-300 2, A+400 1, A-400 1"
        # three cases that AREN'T same-position opposite-strand
        positions_and_readcounts_raw += ", B+100 1, B+101 5, B+200 5, B-201 5, B-300 5, B+301 5"

        ### merging by ratio
        # merge everything (with ratio None, or anything 5 or higher)
        for R in (None, 5,10,100):
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=R,leave_method='by_ratio', OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 4
            assert dataset.summary.merged_opposite_tandems_readcounts == [(1,1),(4,2),(5,1),(5,5)]
            assert dataset.summary.mutants_in_chromosome('chromosome_A') == 4
            assert dataset.summary.mutants_in_chromosome('chromosome_B') == 6
        # merge with ratio 1 only
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio', max_count_ratio=1,leave_method='by_ratio', OUTPUT=None)
        assert dataset.summary.merged_opposite_tandems == 2
        assert dataset.summary.merged_opposite_tandems_readcounts == [(1,1),(5,5)]
        assert dataset.summary.mutants_in_chromosome('chromosome_A') == 6
        assert dataset.summary.mutants_in_chromosome('chromosome_B') == 6
        # merge with ratio 2 and less (but not 5)
        for R in (2,3,4):
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=R,leave_method='by_ratio', OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 3
            assert dataset.summary.merged_opposite_tandems_readcounts == [(1,1),(4,2),(5,5)]
            assert dataset.summary.mutants_in_chromosome('chromosome_A') == 5
            assert dataset.summary.mutants_in_chromosome('chromosome_B') == 6

        ### tests for the new leave_N_mutants and leave_method arguments
        # when leave_N_mutants isn't 'use_ratio', the ratio is ignored, so repeat the test for a few different numbers
        for R in (1,2,10,100, None):
            # merge everything, regardless of the ratio given
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants=0, max_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 4
            # don't merge anything
            for N in (4,5,10,100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_opposite_tandem_mutants(leave_N_mutants=N, max_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_opposite_tandems == 0
            # leave 1 (merge 3), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants=1, leave_method='by_ratio', 
                                           max_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 3
            assert dataset.summary.merged_opposite_tandems_readcounts == [(1,1),(4,2),(5,5)]
            # leave 2 (merge 2), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants=2, leave_method='by_ratio', 
                                           max_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 2
            assert dataset.summary.merged_opposite_tandems_readcounts == [(1,1),(5,5)]
            # leave 3 (merge 1), by ratio
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
            dataset.merge_opposite_tandem_mutants(leave_N_mutants=3, leave_method='by_ratio', 
                                           max_count_ratio=R, OUTPUT=None)
            assert dataset.summary.merged_opposite_tandems == 1
            assert dataset.summary.merged_opposite_tandems_readcounts in [ [(1,1)], [(5,5)] ]
            # leave 1 (merge 3), random
            readcount_pairs = set()
            for _ in range(100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_opposite_tandem_mutants(leave_N_mutants=1, leave_method='random', 
                                               max_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_opposite_tandems == 3
                for pair in dataset.summary.merged_opposite_tandems_readcounts: readcount_pairs.add(pair)
            assert readcount_pairs == set([(1,1),(4,2),(5,1),(5,5)])
            # leave 3 (merge 1), random
            readcount_pairs = set()
            for _ in range(100):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
                dataset.merge_opposite_tandem_mutants(leave_N_mutants=3, leave_method='random', 
                                               max_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_opposite_tandems == 1
                for pair in dataset.summary.merged_opposite_tandems_readcounts: readcount_pairs.add(pair)
            assert readcount_pairs == set([(1,1),(4,2),(5,1),(5,5)])

        ### leave_N_mutants='auto' tests
        positions_and_readcounts_raw = "A+100 5, A-100 1, A+200 5, A-200 5, A+300 4, A-300 2, A+400 1, A-400 1"
        # the resulting leave_N_mutants value should be min((away+toward/2, toward), 0.5 rounded to the even number
        # (+ before - is away, - before + is toward)
        extra_pairs_and_expected_merged = [
            # no opposite-strand adjacent cases - leave none, merge all
            (''                                     + '',                                       4),
            # 0 away case, 1 toward - leave 1, merge 3
            (''                                     + ', 1-550 1, 1+551 1',                     4),
            # 1 away case, 0 toward - leave 0, merge all
            (', 1+500 1, 1-501 1'                   + '',                                       4),
            # 1 away case, 1 toward - leave 2, merge 2
            (', 1+500 1, 1-501 1'                   + ', 1-550 1, 1+551 1',                     3),
            # 2 away case, 1 toward - leave 2, merge 2
            (', 1+500 1, 1-501 1, 1+600 1, 1-601 1' + ', 1-550 1, 1+551 1',                     3),
            # 1 away case, 2 toward - leave 3, merge 1
            (', 1+500 1, 1-501 1'                   + ', 1-550 1, 1+551 1, 1-650 1, 1+651 1',   2),
            # 2 away case, 2 toward - leave 4, merge 0
            (', 1+500 1, 1-501 1, 1+600 1, 1-601 1' + ', 1-550 1, 1+551 1, 1-650 1, 1+651 1',   2),
        ]
        for (extra_pairs, expected_merged) in extra_pairs_and_expected_merged:
            for R in (1,2,10,100, None):
                dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw + extra_pairs)
                dataset.merge_opposite_tandem_mutants(leave_N_mutants='auto', leave_method='by_ratio', 
                                               max_count_ratio=R, OUTPUT=None)
                assert dataset.summary.merged_opposite_tandems == expected_merged

        ### test merge_cassette_chromosomes and merge_other_chromosomes args
        positions_and_readcounts_raw = ("chromosome_1+100 5, chromosome_1-100 1, "
                                        "insertion_cassette+100 5, insertion_cassette-100 1, "
                                        "something_else+100 5, something_else-100 1")
        # merge all
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=None, 
                                              OUTPUT=None, merge_cassette_chromosomes=True, merge_other_chromosomes=True)
        assert len(dataset) == 3
        # don't merge cassette
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=None, 
                                              OUTPUT=None, merge_cassette_chromosomes=False, merge_other_chromosomes=True)
        assert dataset.summary.mutants_in_chromosome('insertion_cassette') == 2
        assert len(dataset) == 4
        # don't merge other
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=None, 
                                              OUTPUT=None, merge_cassette_chromosomes=True, merge_other_chromosomes=False)
        assert len(dataset) == 4
        assert dataset.summary.mutants_in_chromosome('something_else') == 2
        # don't merge cassette or other
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_opposite_tandem_mutants(leave_N_mutants='use_ratio',max_count_ratio=None, 
                                              OUTPUT=None, merge_cassette_chromosomes=False, merge_other_chromosomes=False)
        assert len(dataset) == 5
        assert dataset.summary.mutants_in_chromosome('chromosome_1') == 1
        ### tests for bad arguments
        # multi-dataset
        dataset.multi_dataset = True
        self.assertRaises(MutantError, dataset.merge_opposite_tandem_mutants, 
                          leave_N_mutants='use_ratio',max_count_ratio=None, OUTPUT=None)
        # settings for whether to merge cassette/other chromosomes inconsistent with dataset.summary settings
        for other1,cassette1 in [(True,True), (True,False), (False,True), (False,False)]:
            dataset.summary.merging_which_chromosomes = (cassette1,other1)
            for other2,cassette2 in [(other1,not cassette1), (not other1, cassette1), (not other1, not cassette1)]:
                self.assertRaises(MutantError, dataset.merge_opposite_tandem_mutants, leave_N_mutants='use_ratio',
                                  max_count_ratio=None, OUTPUT=None, 
                                  merge_cassette_chromosomes=cassette1, merge_other_chromosomes=other1)

    def test__count_adjacent_mutants(self):
        positions_and_readcounts_raw = '+100 10, -100 5, +101 5, -99 5, -101 5, +105 1'
        # so we have 6*5/2 = 15 pairs
        #  (remember, for 5' cases, + before - is away-facing, - before + is toward-facing
        #   +100 and -100 -- same-pos opposite, ratio 10:5
        #   +100 and +101 -- same-strand adj, dist 1, ratio 10:5
        #   +100 and -99  -- toward-facing, dist 1, ratio 10:5
        #   +100 and -101 -- away-facing, dist 1, ratio 10:5
        #   +100 and +105 -- same-strand adj, dist 5, ratio 10:1
        #   -100 and +101 -- toward-facing, dist 1, ratio 5:5
        #   -100 and -99  -- same-strand adj, dist 1, ratio 5:5
        #   -100 and -101 -- same-strand adj, dist 1, ratio 5:5
        #   -100 and +105 -- toward-facing, dist 5, ratio 5:1
        #   +101 and -99  -- toward-facing, dist 2, ratio 5:5
        #   +101 and -101 -- same-pos opposite, ratio 5:5
        #   +101 and +105 -- same-strand adj, dist 4, ratio 5:1
        #   -99 and -101  -- same-strand adj, dist 2, ratio 5:5
        #   -99 and +105  -- toward-facing, dist 6, ratio 5:1
        #   -101 and +105 -- toward-facing, dist 4, ratio 5:1
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw)
        dataset.count_adjacent_mutants(max_distance_to_count=None)
        # sorted version of the above list, by category, distance, ratio:
        #   +100 and +101 -- same-strand adj, dist 1, ratio 10:5
        #   -100 and -99  -- same-strand adj, dist 1, ratio 5:5
        #   -100 and -101 -- same-strand adj, dist 1, ratio 5:5
        #   -99 and -101  -- same-strand adj, dist 2, ratio 5:5
        #   +101 and +105 -- same-strand adj, dist 4, ratio 5:1
        #   +100 and +105 -- same-strand adj, dist 5, ratio 10:1
        #   +100 and -100 -- same-pos opposite, ratio 10:5
        #   +101 and -101 -- same-pos opposite, ratio 5:5
        #   +100 and -101 -- away-facing, dist 1, ratio 10:5
        #   +100 and -99  -- toward-facing, dist 1, ratio 10:5
        #   -100 and +101 -- toward-facing, dist 1, ratio 5:5
        #   +101 and -99  -- toward-facing, dist 2, ratio 5:5
        #   -101 and +105 -- toward-facing, dist 4, ratio 5:1
        #   -100 and +105 -- toward-facing, dist 5, ratio 5:1
        #   -99 and +105  -- toward-facing, dist 6, ratio 5:1
        assert dataset.summary.adjacent_same_strand_dict == {1:3, 2:1, 4:1, 5:1}
        assert dataset.summary.adjacent_same_strand_readcounts_dict == {1:[(5,5),(5,5),(10,5)], 2:[(5,5)], 4:[(5,1)], 5:[(10,1)]}
        assert dataset.summary.same_position_opposite == 2
        assert dataset.summary.same_position_opposite_readcounts == [(5,5),(10,5)]
        assert dataset.summary.adjacent_opposite_away_dict == {1:1}
        assert dataset.summary.adjacent_opposite_away_readcounts_dict == {1:[(10,5)]}
        assert dataset.summary.adjacent_opposite_toward_dict == {1:2, 2:1, 4:1, 5:1, 6:1}
        assert dataset.summary.adjacent_opposite_toward_readcounts_dict == {1:[(5,5),(10,5)], 2:[(5,5)], 4:[(5,1)], 
                                                                            5:[(5,1)], 6:[(5,1)]}
        # another case, this time with only dist=1 cases counted
        dataset.count_adjacent_mutants(max_distance_to_count=1)
        assert dataset.summary.adjacent_same_strand_dict == {1:3}
        assert dataset.summary.adjacent_same_strand_readcounts_dict == {1:[(5,5),(5,5),(10,5)]}
        assert dataset.summary.same_position_opposite == 2
        assert dataset.summary.same_position_opposite_readcounts == [(5,5),(10,5)]
        assert dataset.summary.adjacent_opposite_away_dict == {1:1}
        assert dataset.summary.adjacent_opposite_away_readcounts_dict == {1:[(10,5)]}
        assert dataset.summary.adjacent_opposite_toward_dict == {1:2}
        assert dataset.summary.adjacent_opposite_toward_readcounts_dict == {1:[(5,5),(10,5)]}
        ### test merge_cassette_chromosomes and merge_other_chromosomes args
        positions_and_readcounts_raw = ("insertion_cassette+100 5, insertion_cassette-100 1, "
                                        "something_else+100 5, something_else+101 1")
        for cassette,other,tandem,adj in [(True,True,1,1), (True,False,1,0), (False,True,0,1), (False,False,0,0)]:
            dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
            dataset.count_adjacent_mutants(max_distance_to_count=1,OUTPUT=None, 
                                           count_cassette_chromosomes=cassette, count_other_chromosomes=other)
            assert (dataset.summary.same_position_opposite, dataset.summary.adjacent_same_strand_dict) == (tandem, {1:adj})
        ### test different_parameters arg
        dataset = self._make_test_mutant_dataset(positions_and_readcounts_raw, raw_chrom_names=True)
        dataset.merge_adjacent_mutants(merge_max_distance=2, min_count_ratio=5, leave_N_mutants='use_ratio', OUTPUT=None, 
                                       merge_cassette_chromosomes=True, merge_other_chromosomes=False)
        # using same count_* args as for merging, and max_distance >= that for merging - works
        for dist in (2,5):
            dataset.count_adjacent_mutants(max_distance_to_count=dist,OUTPUT=None, 
                                           count_cassette_chromosomes=True, count_other_chromosomes=False)
        # using different count_* args, or max_distance < that for merging - Exception
        self.assertRaises(MutantError, dataset.count_adjacent_mutants, max_distance_to_count=0,OUTPUT=None, 
                                           count_cassette_chromosomes=True, count_other_chromosomes=False)
        self.assertRaises(MutantError, dataset.count_adjacent_mutants, max_distance_to_count=1,OUTPUT=None, 
                                           count_cassette_chromosomes=True, count_other_chromosomes=False)
        self.assertRaises(MutantError, dataset.count_adjacent_mutants, max_distance_to_count=2,OUTPUT=None, 
                                           count_cassette_chromosomes=False, count_other_chromosomes=False)
        self.assertRaises(MutantError, dataset.count_adjacent_mutants, max_distance_to_count=2,OUTPUT=None, 
                                           count_cassette_chromosomes=True, count_other_chromosomes=True)


    def _check_infile1(self, data):
        """ Check that data is as expected in test_data_v0/count-aln__cassette-end-5prime.txt. """
        # summary
        assert data.summary.processed_read_count == 40
        assert data.summary.aligned_read_count == 30
        assert data.summary.non_aligned_read_count == 10
        assert data.summary.perfect_read_count == 22
        assert data.summary.strand_read_counts == {'+':27, '-':3}
        assert len(data) == 17
        assert data.summary.mutants_in_genes == data.summary.mutants_not_in_genes == 0
        assert data.summary.mutant_counts_by_orientation == {}
        assert data.summary.mutant_counts_by_feature == {}
        for mutant in data:
            assert mutant.gene == SPECIAL_GENE_CODES.not_determined
            assert mutant.orientation == mutant.gene_feature == '?'
        # just spot-checking some of the outputs
        mutant = data.get_mutant('reads_2_seqs_1','+',position_before=204)
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 2
        assert mutant.perfect_read_count == 2 
        assert mutant.unique_sequence_count == 1
        mutant = data.get_mutant('mutation_yes','+',position_before=204)
        assert mutant.position.chromosome == 'mutation_yes'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 6
        assert mutant.perfect_read_count == 0
        assert mutant.unique_sequence_count == 1
        mutant = data.get_mutant('strandedness_-_normal_+_reverse','-',position_after=101)
        assert mutant.position.chromosome == 'strandedness_-_normal_+_reverse'
        assert mutant.position.min_position == 100
        assert mutant.position.strand == '-'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1
        assert mutant.unique_sequence_count == 1

    def _check_infile2(self, data2):
        """ Check that data is as expected in test_data_v0/count-aln__with-gene-info_merged.txt. """
        # summary
        assert data2.summary.processed_read_count == data2.summary.aligned_read_count\
                == data2.summary.perfect_read_count == 40
        assert data2.summary.non_aligned_read_count == 'unknown'
        assert data2.summary.strand_read_counts == {'+':38, '-':2}
        assert len(data2) == 40
        assert data2.summary.mutants_in_genes == 39
        assert data2.summary.mutants_not_in_genes == 1
        assert data2.summary.mutants_undetermined == 0
        assert data2.summary.mutant_counts_by_orientation == {'sense':37, 'antisense':2}
        assert data2.summary.mutant_counts_by_feature['CDS'] == 6
        assert data2.summary.mutant_counts_by_feature['??'] == 1
        assert data2.summary.mutant_counts_by_feature["CDS/3'UTR"] == 1
        # mutant spot-checks
        mutant = data2.get_mutant('chromosome_A','+',position_before=20)
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 20
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == SPECIAL_GENE_CODES.not_found
        assert mutant.orientation == mutant.gene_feature == '-'
        mutant = data2.get_mutant('chromosome_A','+',position_before=150)
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 150
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == "test.geneA0_proper_plus"
        assert mutant.orientation == 'sense'
        assert mutant.gene_feature == "5'UTR"

    def test__read_data_from_file(self):
        """ uses _check_infile1 and _check_infile2 functions for detail. """
        ## 1. input file with no gene information but more variation in other features
        input_file = 'test_data_v0/count-aln__cassette-end-5prime.txt'
        data = Insertional_mutant_pool_dataset()
        data.read_data_from_file(input_file)
        self._check_infile1(data)
        ## 2. input file with gene information
        input_file2 = 'test_data_v0/count-aln__with-gene-info_merged.txt'
        data2 = Insertional_mutant_pool_dataset()
        data2.read_data_from_file(input_file2)
        self._check_infile2(data2)
        ## 3. adding more data to a file that already has some...
        data.read_data_from_file(input_file, assume_new_sequences=False)
        mutant = data.get_mutant('reads_2_seqs_1','+',position_before=204)
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 4
        # how mutant.unique_sequence_count should act in this case depends on the value of assume_new_sequences
        assert mutant.unique_sequence_count == 1
        data.read_data_from_file(input_file, assume_new_sequences=True)
        assert mutant.unique_sequence_count == 2

    def test__pickle_unpickle(self):
        """ uses _check_infile1 and _check_infile2 functions for detail. """
        import pickle, os
        picklefile = 'test.pickle'
        ## 1. input file with no gene information but more variation in other features
        input_file = 'test_data_v0/count-aln__cassette-end-5prime.txt'
        data = Insertional_mutant_pool_dataset(infile=input_file)
        self._check_infile1(data)
        # MAYBE-TODO the pickling/unpickling should probably use the general_utilities convenience functions
        with open(picklefile, 'w') as PICKLE:    pickle.dump(data, PICKLE)
        with open(picklefile) as PICKLE:         data_new = pickle.load(PICKLE)
        self._check_infile1(data_new)
        os.unlink(picklefile)
        ## 2. input file with gene information
        input_file2 = 'test_data_v0/count-aln__with-gene-info_merged.txt'
        data2 = Insertional_mutant_pool_dataset(infile=input_file2)
        self._check_infile2(data2)
        with open(picklefile, 'w') as PICKLE:    pickle.dump(data2, PICKLE)
        with open(picklefile) as PICKLE:         data2_new = pickle.load(PICKLE)
        self._check_infile2(data2_new)
        os.unlink(picklefile)


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])

