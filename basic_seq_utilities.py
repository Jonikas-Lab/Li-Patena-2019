#! /usr/bin/env python2.7

"""
Various basic DNA utilities I wrote (parse_fasta, complement, reverse_complement, find_seq_between, generate_seq_slices, ...) - most of them just imported from their own files, since I need them as separate programs too.
Weronika Patena, 2008-2010
"""

### basic library
from __future__ import division 
import unittest
import sys
import re
import random
import collections
import itertools
import math
### other packages
### my modules
# help functions
from parse_fasta import parse_fasta,print_seq
from read_input import read_input
from transform_sequence_input import transform_sequence_input
# independently useful functions
from complement import complement
from reverse_complement import reverse_complement
from reverse import reverse

################## global constants ################## 

NORMAL_DNA_BASES = 'ACTG'

SEQ_ENDS = ['5prime','3prime']
SEQ_ENDS_SHORT = ["5'","3'"]
SEQ_STRANDS = ['+','-']
SEQ_DIRECTIONS = ['forward','reverse']
SEQ_ORIENTATIONS = ['sense','antisense']

# the standard extensions for fasta/fastq files (basic form, processed later)
FASTA_EXTENSIONS = "fa fas fasta"
FASTQ_EXTENSIONS = "fq fastq"

# Biopython fastq quality encodings (standard Phred+33 one, new Illumina Phred+64 one, old Illumina nonPhred+64)
# Can be used like this:  "for read in SeqIO.parse(INFILE, "fastq-illumina"):"
FASTQ_QUALITY_ENCODINGS = ["fastq-sanger", "fastq-illumina", "fastq-solexa"]

# which option to pass to Fastx Toolkit for which encoding:
FASTQ_ENCODINGS_FASTX_TOOLKIT = {'auto': '', 'sanger': '-Q33', 'illumina': '-Q64'}

# colors to use for different chromosome types, by default
CHROMOSOME_TYPE_COLORS = { 'chromosome': 'black', 'scaffold': 'blue', 'chloroplast': 'green', 'mitochondrial': 'red', 
                          'cassette': '0.6'}

################## end of global constants ################## 

### processing global constants
FASTA_EXTENSIONS = (FASTA_EXTENSIONS+" "+FASTA_EXTENSIONS.upper()).split(' ')
FASTQ_EXTENSIONS = (FASTQ_EXTENSIONS+" "+FASTQ_EXTENSIONS.upper()).split(' ')
assert not (set(FASTA_EXTENSIONS) & set(FASTQ_EXTENSIONS))

### formatting of biological data

def format_base_distance(x, approximate=True):
    """ Format base distance as a string (bp, kb, Mb), with or without approximations. 
        
    With approximations:    1bp for 1, 2kb for 1800, 12kb for 12345, 2kb for 2000, 1Mb for 1000000, 2Mb for 2111000. 
    Without approximations - only use the kb/Mb suffixes if the number is evenly divisible by 1k or 1M: 
        1bp for 1, 1800bp for 1800, 12345bp for 12345, 2kb for 2000, 2Mb for 2000000, 2111kb for 2111000.. 
    """
    x = int(x)
    k, M = 10**3, 10**6
    if approximate:
        if x >= M:      return "%sMb"%(int(round(x/M)))
        elif x >= k:    return "%skb"%(int(round(x/k)))
        else:           return "%sbp"%x
    else:
        if x==0:        return "0bp"
        elif not x % M: return "%sMb"%(int(x/M))
        elif not x % k: return "%skb"%(int(x/k))
        else:           return "%sbp"%x

### basic fasta/fastq functions

def write_fasta_line(seqname, seq, OUTFILE=sys.stdout):
    """ Given a name and sequence, print in one-line fasta format to OUTFILE (default STDOUT). """
    OUTFILE.write(">%s\n%s\n"%(seqname, seq))


def write_fastq_line(seqname, seq, qual, OUTFILE=sys.stdout):
    """ Given a name, sequence and quality, print in fastq format to OUTFILE (default STDOUT). """
    OUTFILE.write("@%s\n%s\n+\n%s\n"%(seqname, seq, qual))


def name_seq_generator_from_fasta(fasta_infile):
    """ Yield successive (name,seq) pairs read from fasta file. """
    from Bio import SeqIO   # Bio is the biopython package
    with open(fasta_infile) as INFILE:
        for record in SeqIO.parse(INFILE, 'fasta'):
            # record.seq is a biopython Seq object - we want a string
            yield record.name, str(record.seq)


def name_seq_generator_from_fastq(fastq_infile):
    """ Yield successive (name,seq) pairs read from fastq file (ignores qualities!). """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator    # Bio is the biopython package
    with open(fastq_infile) as INFILE:
        for name,seq,_ in FastqGeneralIterator(INFILE):
            yield name, seq


def check_fasta_fastq_format(infilename, verbose=False):
    """ Check fasta/fastq format based on infilename extension; return "fasta" or "fastq", raise ValueError if neither. 

    Note: the result can be used as the format argument to Bio.SeqIO, but for fastq ONLY if you don't care 
     about the quality encoding (otherwise you should use one of FASTQ_QUALITY_ENCODINGS, not just "fastq").
    """

    import os
    extension = os.path.splitext(infilename)[1][1:]
    if extension in FASTA_EXTENSIONS:       seq_format = "fasta"
    elif extension in FASTQ_EXTENSIONS:     seq_format = "fastq"
    else:
        raise ValueError("File %s has an unknown extension %s! Allowed extensions are fasta (%s) and fastq (%s)."%(
                            infilename, extension, ', '.join(FASTA_EXTENSIONS), ', '.join(FASTQ_EXTENSIONS)))
    if verbose:     
        formatted_output.append("File %s recognized as %s.\n"%(infilename, seq_format))
    return seq_format


def name_seq_generator_from_fasta_fastq(infile, verbose_filetype_detection=False):
    """ Yield successive (name,seq) pairs read from fasta or fastq file (filetype detected by extension). """
    seq_format = check_fasta_fastq_format(infile, verbose=verbose_filetype_detection)
    if seq_format=='fasta':     return name_seq_generator_from_fasta(infile)
    elif seq_format=='fastq':   return name_seq_generator_from_fastq(infile)
    else:                       raise ValueError("Unknown input file format %s (file %s)."%(seq_format, infile))
    

# Note: normally it's probably better to parse fastq using biopython or HTSeq or such!  But it can be convenient to have a simple locally defined function with no dependencies requiring installation, so I'll keep this.
def parse_fastq(infile):
    """ Given a fastq file, yield successive (header,sequence,quality) tuples. """
    with open(infile) as INFILE:
        while True:
            header = INFILE.next().strip()
            try:
                seq = INFILE.next().strip()
                header2 = INFILE.next().strip()
                qual = INFILE.next().strip()
            except (StopIteration):
                raise Exception("Input FastQ file is malformed! Last record didn't have all four lines!")

            if not header[0]=='@':  
                raise Exception("Malformed input fastq file! Expected seq-header line (@ start), found %s"%header)
            if not header2[0]=='+':  
                raise Exception("Malformed input fastq file! Expected qual-header line (+ start), found %s"%header2)
            header,header2 = header[1:],header2[1:]
            if not (len(header2)==0 or header==header2):   
                raise Exception("Malformed input fastq file! Qual-header %s doesn't match seq-header %s"%(header2,header))
            if not len(seq)==len(qual):             
                raise Exception("Malformed input fastq file! Seq length doesn't match qual length (%s,%s)"%(seq, qual))

            yield (header, seq, qual)


def get_seq_count_from_collapsed_header(header, return_1_on_failure=False):
    """ Given a sequence header from fastx_collapser, return the original sequence count ('>1-243' means 243 sequences).
    If cannot parse the header, exits with an error message, unless return_1_on_failure is True (then returns 1). """
    try:                        header_fields = header.split('-')
    except AttributeError:      header_fields = [header]
    if len(header_fields) > 1:
        try:                    return int(header_fields[-1])
        except ValueError:      pass
    # if didn't find a '-' to split on, or couldn't get an int from the string, either return 1 or fail
    if return_1_on_failure: 
        return 1
    else:                   
        raise ValueError("Can't parse header %s to get original pre-fastx_collapser sequence count!"%header)


### analyzing sequences

def GC_content(seq, ignore_other_bases=False):
    """ Return the fraction of the sequence bases that are GC. 

    If sequence contains non-%s bases, raise an exception, or ignore them (not counting them in the total) if ignore_other_bases is True.
    """%NORMAL_DNA_BASES
    seq = seq.upper()
    base_counts = collections.Counter(seq)
    total_normal_bases = sum([base_counts[base] for base in NORMAL_DNA_BASES])
    if not ignore_other_bases and not total_normal_bases==len(seq):
        raise ValueError("sequence passed to GC_content contains non-%s bases! Pass ignore_other_bases=True to ignore them."%NORMAL_DNA_BASES)
    if total_normal_bases==0:
        raise ValueError("sequence passed to GC_content has no %s bases - can't calculate GC content!"%NORMAL_DNA_BASES)
    GC_bases = sum([base_counts[base] for base in 'GC'])
    return GC_bases/total_normal_bases


def check_seq_against_pattern(seq, pattern):
    """ Check if sequence string matches pattern string (ACTGN only); return True/False. Inputs must be same length. """
    seq = seq.upper()
    pattern = pattern.upper()
    if not len(seq)==len(pattern):
        raise ValueError("In check_seq_against_pattern, sequence and pattern must be same length! %s, %s"%(seq,pattern))
        # MAYBE-TODO add an option to pad one or the other (left,right,centered) if different lengths?
    if not all([b in NORMAL_DNA_BASES+'N' for b in seq+pattern]):
        raise ValueError("In check_seq_against_pattern, sequence and pattern must contain ACTGN only! %s, %s"%(seq,pattern))
        # MAYBE-TODO expand to the other ambiguous bases, like W, R, etc?
        # MAYBE-TODO expand to RNA?
    return all([(p==s or 'N' in (s,p)) for (s,p) in zip(seq,pattern)])


def get_all_seq_length(seq_list):
    """ Given a list of sequence strings, return N if all seqs are length N; raise ValueError if lengths differ or list was empty."""
    if not seq_list:
        raise ValueError("Empty sequence list passed - cannot give length!")
    seq_lengths = set([len(seq) for seq in seq_list])
    if len(seq_lengths) > 1:
        raise ValueError("Sequences have different lengths - cannot give single length! %s"%seq_lengths)
    return seq_lengths.pop()


def base_count_dict(seq_count_list, convert_counts=lambda x: x, skip_seq_length_check=False):
    """ Given a list of (seq,count) tuples, return a base:base_count_list_for_each_position dict.

    (For example input [(ACG,1),(AAC,1)] will yield base counts of {A:[2,1,0], C:[0,1,1], G:[0,0,1], T:[0,0,0]}.)
    Convert_counts is a function that takes the readcount for a sequence and outputs how many times it should be counted:
     the reasonable output values are always 1, the readcount itself, or some rounding of readcount/average_readcount_per_mutant.

    Normally program checks that all sequences are the same length first - to disable this (for large inputs etc), 
     pass skip_seq_length_check=True.

    Ignore bases other than %s.
    """%NORMAL_DNA_BASES
    # if input list is empty, return empty lists for each base
    if not seq_count_list:
        return {base: [] for base in NORMAL_DNA_BASES}
    # get seq length, make sure they're all the same (or skip if desired
    if skip_seq_length_check:   seq_length = len(seq_count_list[0][0])
    else:                       seq_length = get_all_seq_length(zip(*seq_count_list)[0])
    # initialize the base-count lists to the right length, and fill it out by going over all the seqs
    base_count_dict = {base: [0 for _ in range(seq_length)] for base in NORMAL_DNA_BASES}
    for (seq,count) in seq_count_list:
        for position, base in enumerate(seq.upper()):
            try:
                base_count_dict[base][position] += convert_counts(count)
            except KeyError:
                pass
                # MAYBE-TODO add an option to NOT ignore bases that aren't in NORMAL_DNA_BASES?
    return base_count_dict


def base_fraction_dict_from_count_dict(base_count_list_dict):
    """ Take output of base_count_dict, return normalized fractions instead of counts. 

    Ignore bases other than %s.  If there are no valid bases at some position, return all 'NaN' fractions for that position.
    """%NORMAL_DNA_BASES
    base_fraction_list_dict = {}
    base_totals = [sum(single_pos_counts) for single_pos_counts in zip(*base_count_list_dict.values())]
    for base in NORMAL_DNA_BASES:
        base_fraction_list_dict[base] = [(count/total if total else float('nan')) 
                                         for (count,total) in zip(base_count_list_dict[base], base_totals)]
        # MAYBE-TODO add an option to NOT ignore bases that aren't in NORMAL_DNA_BASES?
    return base_fraction_list_dict


def base_fraction_dict(seq_count_list, convert_counts=lambda x: x):
    """ Same as base_count_dict, but returns normalized fractions instead of counts. 

    Ignore bases other than %s.  If there are no valid bases at some position, return all 'NaN' fractions for that position.
    """%NORMAL_DNA_BASES
    return base_fraction_dict_from_count_dict(base_count_dict(seq_count_list, convert_counts))


def base_fractions_from_GC_content(overall_GC_content):
    """ Given a GC content, return base:fraction dict. Example: 0.6 gives {0.3 for G/C, 0.2 for A/T}. """
    if not 0 <= overall_GC_content <= 1:
        raise ValueError("overall_GC_content must be a number between 0 and 1 (inclusive)!")
    return { 'G': overall_GC_content/2, 'A': (1-overall_GC_content)/2, 
             'C': overall_GC_content/2, 'T': (1-overall_GC_content)/2 }


### utilities to deal with standard genomic chromosome names (sort them correctly etc)

def chromosome_type(chromosome, separate_other=False):
    """ Given the chromosome name, return type: chromosome/scaffold/chloroplast/mitochondrial/cassette/other.

    Check for chromosome/scaffold/chloroplast/mitochondrial/cassette (by substrings or startswith);
     if a name didn't match any, or matched more than one, return 'other', 
     or if separate_other is True, return the first _-separated word of the original name.
    """
    # first see which words are found in the name (check for all of them, just in case more than one is found)
    types = set()
    if chromosome.startswith('chr'):    types.add('chromosome')
    if 'chromosome' in chromosome:      types.add('chromosome')
    if 'scaffold' in chromosome:        types.add('scaffold')
    if 'chloro' in chromosome:          types.add('chloroplast')
    if 'mito' in chromosome:            types.add('mitochondrial')
    if 'cassette' in chromosome:        types.add('cassette')
    # now if the name matched exactly one type, return that, otherwise return 'other' or the first part of the name.
    if len(types)==1:       return types.pop()
    elif separate_other:    return chromosome.split('_')[0]
    else:                   return 'other'


def chromosome_sort_key(chromosome_name, special_search_strings=[]):
    """ Sort key: 1) sort chromosomes, then other, then scaffolds, then special; 2) in each category sort numerically/alpha.

    Numerical sort sorts 'chr1' before 'chr12' correctly.  Only numbers at the end of the chromosome name are recognized.
    The full original chromosome_name is used as the last part of the key, in case of two names resulting in the same key
     (like chr1 and chr_1 and chr01).
    """
    chromosome_data = re.search('^(?P<name>.*[^\d])?(?P<number>\d*)', chromosome_name)
    # chromosome_base is the original chromosome name without the last digits; if there were no non-digit characters, it's None
    chromosome_base = chromosome_data.group('name').strip('_')
    # chromosome_number is based on the number at the end of the chromosome name, or 0 if none was present
    chromosome_number = int(chromosome_data.group('number')) if chromosome_data.group('number') else 0
    if chromosome_base in ('chromosome', 'chr'):    return (1, 'chromosome', chromosome_number, chromosome_name)
    elif chromosome_base=='scaffold':               return (3, 'scaffold', chromosome_number, chromosome_name)
    elif any(x in chromosome_base.lower() for x in special_search_strings):    
                                                    return (4, chromosome_base, chromosome_number, chromosome_name)
    else:                                           return (2, chromosome_base, chromosome_number, chromosome_name) 


def chromosome_color(chromosome, color_for_other=None):
    """ Return a color by chromosome type: black=chromosome, blue=scaffold, green=chloroplast, red=mitochondrial, grey=cassette.

    Uses chromosome_type to determine type. 
    For other types, raise exception if color_for_other is None, otherwise return color_for_other.
    """
    chr_type = chromosome_type(chromosome)
    try:
        return CHROMOSOME_TYPE_COLORS[chr_type]
    except KeyError:
        if color_for_other is not None:     return color_for_other
        else:                               raise ValueError("Can't parse chromosome name %s!"%chromosome)


### testing whether two sequences contain one another, or overlap (based purely on position, no sequence involved)

def position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):
    """ Return True if seq2 is completely contained within seq1, False otherwise. All positions should be numbers."""
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    # the test is simple
    return bool(seq1_start <= seq2_start and seq2_end <= seq1_end)

def position_test_overlap(seq1_start, seq1_end, seq2_start, seq2_end, return_overlap=False):
    """ Returns True if seq1/seq2 overlap, else False; or the number of bases of overlap if return_overlap is True (min 0).
    Assumes the positions are inclusive (i.e. length = end-start+1). """
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    ### if we don't want to measure the overlap, the test is simple
    if return_overlap is False:     return bool(seq2_start <= seq1_end and seq2_end >= seq1_start)
    ### if we do want to measure the overlap, need to do more work:
    # if seq1 contains seq2, or the opposite, return the length of the shorter sequence (see docstring for why +1)
    if position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):  return seq2_end - seq2_start + 1
    if position_test_contains(seq2_start, seq2_end, seq1_start, seq1_end):  return seq1_end - seq1_start + 1
    # if seq2 starts before seq1, or ends after seq1, the overlap is the distance between the "inner" start and end
    #   (there's a +1 in both cases, because if seq1_end==seq2_start, there's 1bp overlap)
    if seq2_start < seq1_start:                                             return max(seq2_end - seq1_start + 1, 0)  
    if seq2_end > seq1_end:                                                 return max(seq1_end - seq2_start + 1, 0)  


### other position-based functions

def find_seq_between(seq, flankA, flankB, exclusive=True):
    """ Return the fragment of seq between flankA and flankB (first occurence). """
    posA = seq.find(flankA)
    posB = seq[posA:].find(flankB)  # only look for flankB in the region AFTER flankA
    if not posB==-1:    posB += posA        # if found, add the length of the omitted region to result 
    lenA = len(flankA)
    lenB = len(flankB)
    if posA==-1 or posB==-1:    return ''
    elif exclusive:             return seq[posA+lenA:posB]
    else:                       return seq[posA:posB+lenB]

def generate_seq_slices(seq, slice_len, step=1, extra_last_slice=True):
    """ Generate (start_pos,slice_seq) pairs with slices covering seq with the given slice length and step size. 

    Start_pos is one-based.
    
    When seq doesn't cleanly divide into slices with given step/slice_len (i.e. len(seq)-slice_len isn't divisible by step)
     (for example for 'aaccttgg' with step=3 and slice_len=4):
     - if extra_last_slice is False, ignore the remainder of seq (at most slice_len-1) 
        (just generate 'aacc', 'cttg' - the last 'g' is ignored)
     - if it's True, yield an extra slice from the last slice_len bases of seq, that's less than step distant from the previous slice
        (after 'aacc' and 'cttg' yield 'ttgg' to cover the last 'g', even though the offset between 'cttg' and 'ttgg' is 1, not 3)
    """
    # note: convert all positions to 1-based by adding 1!
    if slice_len <= 0:   raise ValueError("slice_len argument to generate_seq_slices must be a positive integer!")
    if step <= 0:        raise ValueError("step argument to generate_seq_slices must be a positive integer!")
    seqlen = len(seq)
    # if the seq will fit in one slice, just return that
    if slice_len >= seqlen: 
        yield 1, seq
        return
    # otherwise go over all slices and yield them
    for pos in range(0, seqlen-slice_len, step):
        yield pos+1, seq[pos:pos+slice_len]
    # Special case after the last normal slice, if there's "leftover" sequence:
    #  if the sequence isn't evenly divided into correct-length slices, make one last slice covering the last slice_len 
    #   of the seq (the difference between this slice and the previous one will be <step)
    #  But only if the next slice would be inside seq at all!  
    #   if slice<step, for example 'actg' with len 1 and step 2 should be 'a t', never 'a t g'.
    if (pos+slice_len < seqlen) and (pos+step < seqlen) and extra_last_slice:
        yield seqlen-slice_len+1, seq[-slice_len:]


####################################### Unit-tests #########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__format_base_distance(self):
        for approx in (True,False):
            assert format_base_distance(0, approx) == "0bp"
        # with no approximation
        assert format_base_distance(1, False) == '1bp'
        assert format_base_distance(101, False) == '101bp'
        assert format_base_distance(999, False) == '999bp'
        assert format_base_distance(1001, False) == '1001bp'
        assert format_base_distance(1000001, False) == '1000001bp'
        for x in [1,10,11,100,999, 1001, 10001, 1234, 13987291876]:
            assert format_base_distance(x, False) == "%sbp"%x
            assert format_base_distance(x*1000, False) == "%skb"%x
            assert format_base_distance(x*1000000, False) == "%sMb"%x
            assert format_base_distance(x*1000+1, False) == "%sbp"%(x*1000+1)
            assert format_base_distance(x*1000000+1, False) == "%sbp"%(x*1000000+1)
        # with approximation
        assert format_base_distance(1, True) == '1bp'
        assert format_base_distance(101, True) == '101bp'
        assert format_base_distance(999, True) == '999bp'
        assert format_base_distance(1001, True) == '1kb'
        assert format_base_distance(1000001, True) == '1Mb'
        for base in (1000, 2000, 50000, 999000):
            for plus in (0,1, 100, 499):
                assert format_base_distance(base+plus, True) == '%skb'%(int(base/1000))
                assert format_base_distance(base*1000+plus, True) == '%sMb'%(int(base/1000))
                assert format_base_distance((base+plus)*1000, True) == '%sMb'%(int(base/1000))
        for x in ['a', '1.23', [], [12]]:
            self.assertRaises((ValueError,TypeError), format_base_distance, x)


    def test__parse_fastq(self):
        # need to actually run through the whole iterator to test it - defining it isn't enough
        def parse_fastq_get_first_last(infile):
            seq_iter = parse_fastq(infile)
            seq1 = seq_iter.next()
            for seq in seq_iter:    
                seqN = seq
            return seq1, seqN
        # checking first and last record of _test_inputs/test.fq
        seq1, seqN = parse_fastq_get_first_last("_test_inputs/test.fq")
        assert seq1 == ("ROCKFORD:4:1:1680:975#0/1", "NCTAATACGCGGCCTGGAGCTGGACGTTGGAACCAA", 
                        "BRRRQWVWVW__b_____^___bbb___b_______")
        assert seqN == ("ROCKFORD:4:1:3367:975#0/1", "NCTAAGGCAGATGGACTCCACTGAGGTTGGAACCAA", 
                        "BQQQNWUWUUbbb_bbbbbbbbb__b_bb_____b_") 
        # non-fastq input files
        self.assertRaises(Exception, parse_fastq_get_first_last, "_test_inputs/test.fa")
        self.assertRaises(Exception, parse_fastq_get_first_last, "_test_inputs/textcmp_file1.txt")


    def test__get_seq_count_from_collapsed_header(self):
        for bad_header in ['aaa','aaa-aa', 'a-3-a', 'a-3a', '3-a','a+3','a:3','a 3','3',3,0,100,None,True,False,[],{}]:
            assert get_seq_count_from_collapsed_header(bad_header, return_1_on_failure=True) == 1
            self.assertRaises(ValueError, get_seq_count_from_collapsed_header, bad_header, return_1_on_failure=False) 
        for header_prefix in ['a','a ','>a','> a','a b c','a-b-c','a-3-100']:
            for count in [0,1,3,5,123214]:
                assert get_seq_count_from_collapsed_header(header_prefix+'-'+str(count)) == count


    def test__GC_content(self):
        for no_GC_seq in 'aaaa AAA a TT ATTTta'.split():
            assert GC_content(no_GC_seq) == 0
        for all_GC_seq in 'GGGG gggg C cc GCcg'.split():
            assert GC_content(all_GC_seq) == 1
        for half_GC_seq in 'ATGC tagc AAAGGC ttGGGt'.split():
            assert GC_content(half_GC_seq) == 0.5
        seq_with_Ns = 'ATGCNNNNNNNN'
        self.assertRaises(ValueError, GC_content, seq_with_Ns, ignore_other_bases=False)
        assert GC_content(seq_with_Ns, ignore_other_bases=True) == 0.5
        assert GC_content('ATTG') == GC_content('ATTGnnnnnnnnnn', True) == 0.25

    def test__check_seq_against_pattern(self):
        # fails with illegal characters
        for bad_seq in 'fas 123 hjakdsh ATGW'.split():
            same_length_seq = ''.join('G' for _ in bad_seq)
            self.assertRaises(ValueError, check_seq_against_pattern, bad_seq, same_length_seq)
            self.assertRaises(ValueError, check_seq_against_pattern, same_length_seq, bad_seq)
        # fails with different lengths
        for seq in 'ACT AAAAAA G'.split():
            self.assertRaises(ValueError, check_seq_against_pattern, seq, '')
            self.assertRaises(ValueError, check_seq_against_pattern, '', seq)
            self.assertRaises(ValueError, check_seq_against_pattern, seq, 'NNNNNNNNNNN')
            self.assertRaises(ValueError, check_seq_against_pattern, 'NNNNNNNNNNN', seq)
            self.assertRaises(ValueError, check_seq_against_pattern, seq, seq+'A')
            self.assertRaises(ValueError, check_seq_against_pattern, seq+'A', seq)
        # True if seq==pattern
        for seq in 'ACT AAAAAA G'.split():
            assert check_seq_against_pattern(seq, seq)
        assert check_seq_against_pattern('', '')
        # True if either seq or pattern is all Ns
        for seq in 'ACT AAAAAA G'.split():
            same_length_Ns = ''.join('N' for _ in seq)
            assert check_seq_against_pattern(seq, same_length_Ns)
            assert check_seq_against_pattern(same_length_Ns, seq)
        # more complex cases by hand
        assert check_seq_against_pattern('AAA', 'ANN')
        assert check_seq_against_pattern('AAT', 'ANN')
        assert check_seq_against_pattern('AAA', 'ANA')
        assert not check_seq_against_pattern('AAT', 'ANA')
        assert not check_seq_against_pattern('AAA', 'GNN')
        assert check_seq_against_pattern('GGG', 'GNN')

    def test__get_all_seq_length(self):
        self.assertRaises(ValueError, get_all_seq_length, [])
        self.assertRaises(ValueError, get_all_seq_length, ['AA', 'C'])
        self.assertRaises(ValueError, get_all_seq_length, ['AA', 'CCCCCCCC'])
        for length in (0,1,10,1000):
            assert get_all_seq_length(['A'*length]) == length
            assert get_all_seq_length(['N'*length]) == length
            assert get_all_seq_length(['A'*length, 'A'*length]) == length
            assert get_all_seq_length(['N'*length, 'C'*length, 'G'*length]) == length
            assert get_all_seq_length(['AC'*length, 'GT'*length]) == length*2
            self.assertRaises(ValueError, get_all_seq_length, ['A'*length, 'A'*(length+1)])

    def test__base_count_dict(self):
        # basic functionality
        assert base_count_dict([]) == {base:[] for base in NORMAL_DNA_BASES}
        assert base_count_dict([('ACG',1),('AAC',1)]) == {'A':[2,1,0], 'C':[0,1,1], 'G':[0,0,1], 'T':[0,0,0]}
        # ignores Ns
        assert base_count_dict([('ACG',1),('AAC',1),('NNN',10)]) == {'A':[2,1,0], 'C':[0,1,1], 'G':[0,0,1], 'T':[0,0,0]}
        # fails with different-length sequences
        self.assertRaises(ValueError, base_count_dict, [('ACG',1),('AA',1)])
        # testing the convert_counts function
        assert base_count_dict([('ACG',1),('AAC',1)], lambda x:0) ==  {base:[0,0,0] for base in NORMAL_DNA_BASES}
        assert base_count_dict([('ACG',1),('AAC',1)], lambda x:10*x) ==  {'A':[20,10,0], 'C':[0,10,10], 'G':[0,0,10], 'T':[0,0,0]}
        for N1, N2 in itertools.product((1,2,10,100), (1,2,10,100)):
            assert base_count_dict([('ACG',N1),('AAC',N2)], lambda x:1) == {'A':[2,1,0], 'C':[0,1,1], 'G':[0,0,1], 'T':[0,0,0]}
        assert base_count_dict([('ACG',1),('AAC',10)]) ==  {'A':[11,10,0], 'C':[0,1,10], 'G':[0,0,1], 'T':[0,0,0]}

    def test__base_fraction_dict(self):
        # this mostly just uses base_count_dict, so I'm not going to test it too carefully
        assert base_fraction_dict([]) == {base:[] for base in NORMAL_DNA_BASES}
        assert base_fraction_dict([('ACG',1),('AAC',1)]) == {'A':[1,0.5,0], 'C':[0,0.5,0.5], 'G':[0,0,0.5], 'T':[0,0,0]}
        assert base_fraction_dict([('ACG',1),('AAC',1),('NNN',10)]) == {'A':[1,0.5,0], 'C':[0,0.5,0.5], 'G':[0,0,0.5], 'T':[0,0,0]}
        # checking a case where some bases are always N, so the total is 0 and all the fractions are NaN - more complicated
        #  due to NaN comparison (NaN!=NaN), so need to split it up and use isnan function.
        NaN = 'nan'
        output = base_fraction_dict([('CCN',1),('ANN',1)]) 
        transformed_output = {key:[NaN if math.isnan(x) else x for x in val] for key,val in output.items()}
        assert transformed_output == {'A':[0.5,0,NaN], 'C':[0.5,1,NaN], 'G':[0,0,NaN], 'T':[0,0,NaN]}

    def test__base_fractions_from_GC_content(self):
        for bad_GC_content in (-1, -0.000000000001, 2, 100, 1.000000001, 'ACTA', []):
            self.assertRaises(ValueError, base_fractions_from_GC_content, bad_GC_content)
        for GC_content in [x/100 for x in range(100)]:
            # fuzzy comparison of floats
            assert 0.999999 <= sum(base_fractions_from_GC_content(GC_content).values()) <= 1.000001
            assert set(base_fractions_from_GC_content(GC_content).keys()) == set(NORMAL_DNA_BASES)
        assert base_fractions_from_GC_content(0.5) == {base:0.25 for base in NORMAL_DNA_BASES}
        assert base_fractions_from_GC_content(0.2) == {'C':0.1, 'G':0.1, 'A':0.4, 'T':0.4}

    def test__chromosome_type(self):
        for chrom in ('chromosome_1 chromosome_12 chromosome_FOO chromosome_1_2_3 chromosome1 chr_3 chr4'.split()):
            assert chromosome_type(chrom) == 'chromosome'
        for chrom in ('scaffold_1 scaffold_12 scaffold_FOO scaffold_1_2_3 scaffold2 scaffold'.split()):
            assert chromosome_type(chrom) == 'scaffold'
        for chrom in ('mitochondrial mitochondrion mito_3 mito3'.split()):
            assert chromosome_type(chrom) == 'mitochondrial'
        for chrom in ('chloroplast chloroplast_4 chloroplastABC chloroplast123 chloro_3 chloro3'.split()):
            assert chromosome_type(chrom) == 'chloroplast'
        for chrom in ('cassette insertion_cassette cassette_pMJ0013 insertion_cassette_A'.split()):
            assert chromosome_type(chrom) == 'cassette'
        for chrom in ('123 something random mitochondrial_chromosome chr_chloro mito_and_chloro mito_cassette'.split()):
            assert chromosome_type(chrom, separate_other=False) == 'other'
        assert chromosome_type('123', separate_other=True) == '123'
        assert chromosome_type('12_3', separate_other=True) == '12'
        assert chromosome_type('something', separate_other=True) == 'something'
        assert chromosome_type('mito_and_chloro', separate_other=True) == 'mito'
        assert chromosome_type('foo_bar_baz', separate_other=True) == 'foo'

    def test__chromosome_sort_key(self):
        # testing standard sorting
        chroms_sorted = (
            'chr chromosome_1 chr_2 chr03 chr3 chr_3 chromosome_3 chromosome_12 chromosome_101 chromosome_300 chr301'.split()
            +'31a AAA AAA3 cassette chloroplast insertion_cassettes_31_a some_thing something31a'.split()
            +'scaffold02 scaffold2 scaffold_3 scaffold_3 scaffold_21'.split()
        )
        for _ in range(100):
            chroms = list(chroms_sorted)
            random.shuffle(chroms)
            assert sorted(chroms, key=chromosome_sort_key) == chroms_sorted
        # testing putting custom special cases at the end 
        chroms_sort_1, chroms_sort_2 = 'chromosome_1 X_special scaffold_3'.split(), 'chromosome_1 scaffold_3 X_special'.split()
        assert sorted(chroms_sort_1, key = chromosome_sort_key) == chroms_sort_1
        assert sorted(chroms_sort_1, key = lambda c: chromosome_sort_key(c, special_search_strings=['special'])) == chroms_sort_2

    def test__chromosome_color(self):
        assert all([chromosome_color(c) == 'black' for c in 'chromosome_1 chromosome chromosome_A'.split()])
        assert all([chromosome_color(c) == 'blue' for c in 'scaffold_1 scaffold scaffold_A'.split()])
        assert all([chromosome_color(c) == 'green' for c in 'chloroplast chloroplast_A'.split()])
        assert all([chromosome_color(c) == 'red' for c in 'mitochondrial mitochondrial_A'.split()])
        assert all([chromosome_color(c) == '0.6' for c in 'cassette insertion_cassette cassette_pMJ013b'.split()])
        for other_chr in 'some_thing foo_bar mito_and_chloro'.split():
            assert chromosome_color(other_chr, color_for_other='cyan') == 'cyan'
            self.assertRaises(ValueError, chromosome_color, other_chr, color_for_other=None)


    def test__position_test_contains(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_contains, 11,10,1,1)
        self.assertRaises(ValueError, position_test_contains, 1,1,11,10)
        # seq1 contains seq2 - True
        assert position_test_contains(1,10,2,3)==True
        assert position_test_contains(1,10,9,10)==True
        assert position_test_contains(1,10,10,10)==True
        assert position_test_contains(10,10,10,10)==True
        # seq2 contains seq1 - False
        assert position_test_contains(2,3,1,10)==False
        assert position_test_contains(9,10,1,10)==False
        assert position_test_contains(10,10,1,10)==False
        # seq2 and seq1 overlap - False
        assert position_test_contains(2,3,3,10)==False
        assert position_test_contains(2,4,3,10)==False
        assert position_test_contains(3,10,2,3)==False
        assert position_test_contains(3,10,2,4)==False
        # seq2 and seq1 don't overlap - False
        assert position_test_contains(2,3,4,10)==False
        assert position_test_contains(2,3,6,10)==False
        assert position_test_contains(4,10,2,3)==False
        assert position_test_contains(6,10,2,3)==False


    def test__position_test_overlap__returns_bool(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, False)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, False)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,False)==True
        assert position_test_overlap(1,10,9,10,False)==True
        assert position_test_overlap(1,10,10,10,False)==True
        assert position_test_overlap(10,10,10,10,False)==True
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,False)==True
        assert position_test_overlap(9,10,1,10,False)==True
        assert position_test_overlap(10,10,1,10,False)==True
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,False)==True
        assert position_test_overlap(2,4,3,10,False)==True
        assert position_test_overlap(3,10,2,3,False)==True
        assert position_test_overlap(3,10,2,4,False)==True
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,False)==False
        assert position_test_overlap(2,3,6,10,False)==False
        assert position_test_overlap(4,10,2,3,False)==False
        assert position_test_overlap(6,10,2,3,False)==False

    def test__position_test_overlap__returns_number(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, True)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, True)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,True)==2
        assert position_test_overlap(1,10,9,10,True)==2
        assert position_test_overlap(1,10,10,10,True)==1
        assert position_test_overlap(10,10,10,10,True)==1
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,True)==2
        assert position_test_overlap(9,10,1,10,True)==2
        assert position_test_overlap(10,10,1,10,True)==1
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,True)==1
        assert position_test_overlap(2,4,3,10,True)==2
        assert position_test_overlap(3,10,2,3,True)==1
        assert position_test_overlap(3,10,2,4,True)==2
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,True)==0
        assert position_test_overlap(2,3,6,10,True)==0
        assert position_test_overlap(4,10,2,3,True)==0
        assert position_test_overlap(6,10,2,3,True)==0


    def test__find_seq_between(self):
        # basic
        for input_seq in ['TTTAGAGCCC', 'CCTTTAGAGCCCTT', 'AAAGAGAGATTTAGAGCCCAGAGAGAGGGGNNNNNN']:
            assert find_seq_between(input_seq, 'TTT', 'CCC', exclusive=True)  == 'AGAG'
            assert find_seq_between(input_seq, 'TTT', 'CCC', exclusive=False) == 'TTTAGAGCCC'
        # if the edge seqs qre not in input_seq, or are in the wrong order, return ''
        for input_seq in ['TTAGAGCCC', 'TTTAGAGCC', 'TTAGAGCC', 'AGAG', 'CCCAGAGTTT', 'CCCTTT']:
            assert find_seq_between(input_seq, 'TTT', 'CCC', exclusive=True)  == ''
            assert find_seq_between(input_seq, 'TTT', 'CCC', exclusive=False) == ''
        # if there are multiple occurrences of the edge seqs, use the first one
        assert find_seq_between('TTTAATTTAAAACCC', 'TTT', 'CCC', exclusive=True)  == 'AATTTAAAA'
        assert find_seq_between('TTTAATTTAAAACCCAACCC', 'TTT', 'CCC', exclusive=True)  == 'AATTTAAAA'
        assert find_seq_between('TTTAAAACCCAACCC', 'TTT', 'CCC', exclusive=True)  == 'AAAA'

    def test__generate_seq_slices(self):
        # arguments to generate_seq_slices are (seq, slice_len, step);
        #  the output is a list of (pos,slice_seq) tuples, 
        #   which I'm mostly switching to [pos_tuple, seq_tuple] with zip for easier typing.
        # testing step==1
        assert zip(*generate_seq_slices('actg',1,1)) == [(1,2,3,4), tuple('a c t g'.split())]
        assert zip(*generate_seq_slices('actg',2,1)) == [(1,2,3), tuple('ac ct tg'.split())]
        assert zip(*generate_seq_slices('actg',3,1)) == [(1,2), tuple('act ctg'.split())]
        for slice_len in (4,5,10,100,12345):
            assert list(generate_seq_slices('actg',slice_len,1)) == [(1, 'actg')]
        # testing step==2, and the extra_last_slice argument
        assert list(generate_seq_slices('actg',1,2)) == [(1,'a'), (3,'t')]
        assert list(generate_seq_slices('actg',2,2)) == [(1,'ac'), (3,'tg')]
        assert list(generate_seq_slices('actg',3,2,False)) == [(1,'act')]
        assert list(generate_seq_slices('actg',3,2,True)) == [(1,'act'),(2,'ctg')]
        for slice_len in (4,5,10,100,12345):
            assert list(generate_seq_slices('actg',slice_len,1)) == [(1,'actg')]
        ### bad inputs - note that just making the generator doesn't raise an exception, only iterating over it with G.next!
        # step can't be 0, negative, or a float (except that it can be a float if seqlen<=slice_len, since it's not used at all then)
        for slice_len in (1,2,3,4,5,10,100,12345):
            for step in (0, -1, -2, -100, -1235):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises(ValueError, G.next)
        for slice_len in (1,2,3):
            for step in (.5, 1.5, 2.1, 9.0001, 4.999, 1235453534.1):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises(TypeError, G.next)
        # slice_len can't be 0 or negative, or a float (but it can be a float if it's >=seqlen, since then it never gets used)
        for step in (1,2,3,4,5,10,100,12345):
            for slice_len in (0, -1, -2, -100, -1314, 0.5, 2.1, -100.1, -0.01):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises((ValueError,TypeError), G.next)


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()

