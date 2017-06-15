#! /usr/bin/env python2.7

"""
Various help functions for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import sys
import unittest
import os
from collections import defaultdict
import glob
# other packages
import numpy
# my modules
import general_utilities
import basic_seq_utilities
import mutant_IB_RISCC_classes

# various defaults/constants
STRAND_VAR_VALUES = ('+', '-', 'both', None)
DEFAULT_BIN_SIZE = 20000

DEFAULT_NUCLEAR_GENOME_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Creinhardtii_281_v5.0.fa')
DEFAULT_ALL_GENOME_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre5.5-nm_chl-mit.fa')
DEFAULT_GENOME_CASSETTE_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre5.5-nm_chl-mit_cassette-CIB1.fa')
DEFAULT_GENE_POS_FILE = os.path.expanduser('~/experiments/reference_data/chlamy_annotation/Creinhardtii_281_v5.5.gene.gff3')

def get_chromosome_lengths(genome_file=None):
    """ Return chromosome:length dictionary based on reading a genome fasta file. """
    original_input = genome_file
    if genome_file is None:
        genome_file = DEFAULT_GENOME_CASSETTE_FILE
    chromosome_lengths = defaultdict(int)
    try:
        for header,seq in basic_seq_utilities.parse_fasta(genome_file):
            chromosome_lengths[header] = len(seq)
        return dict(chromosome_lengths)
    except IOError:
        file_info = "default " if original_input is None else ""
        raise ValueError("%sgenome fasta file %s not found! Provide filename."%(file_info, genome_file))
    # MAYBE-TODO should this be in basic_seq_utilities or somewhere?  Except for the specific default value...


def get_mutant_positions_from_dataset(dataset, strand=None):
    """  Return chromosome_name:mutant_position_list for dataset.

    Dataset must be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list/set/something of mutant_analysis_classes.Insertional_mutant instances.
    Use the (known or assumed) position of the base before the insertion (min_position).

    If strand is None, take mutants regardless of strand; 
     if it's '+', '-' or 'both', take only mutants on that strand (both-stranded mutants are opposite-tandems); 
     all other strand values are illegal.
    """
    if not strand in STRAND_VAR_VALUES:
        raise ValueError("Illegal strand value %s! Must be one of %s"%(strand, STRAND_VAR_VALUES))
    chromosome_position_dict = defaultdict(list)
    for mutant in dataset:
        position = mutant.position
        if strand is None or position.strand==strand:
            chromosome_position_dict[position.chromosome].append(position.min_position)
    return chromosome_position_dict


def merge_dataset_files(file_list=[], file_glob_pattern=None):
    """ Return single mutant dataset from adding all the input ones together (input can be filename list or glob pattern).
    """
    if (file_list and file_glob_pattern) or not (file_list or file_glob_pattern):
        raise Exception("Must provide exactly one of file_list and file_glob_pattern!")
    if file_glob_pattern:
        file_list = glob.glob(file_glob_pattern)
    # make empty dataset
    full_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset()
    # read each listed dataset file, merge into the full dataset, then delete single dataset
    for infile in file_list:
        curr_dataset = mutant_analysis_classes.read_mutant_file(infile)
        full_dataset.merge_other_dataset(curr_dataset)
        del curr_dataset
    return full_dataset


def get_histogram_data_from_positions(position_dict, bin_size=DEFAULT_BIN_SIZE, chromosome_lengths=None, chromosomes=None, 
                                      first_bin_offset=0, special_last_bin=True, merge_last_bin_cutoff=0.5, normalize_last_bin=True):
    """ Given a chromosome:position_list dict, return a chromosome:counts_per_bin dict, with counts_per_bin a numpy array. 
    
    The input can actually be either a chromosome:position_list or a (chromosome,strand):position_list dict - 
     the results will be the same for either, the different-strand lists will just be added together. 

    The positions in position_list will be binned into bin_size-sized bins over the length of each chromosome, 
     giving counts_per_bin numpy array, with length matching the number of bins in the chromosome. 
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used) - use the lengths to decide how many bins there will be in the chromosome.

    Skip and ignore the first first_bin_offset of each chromosome.
    If special_last_bin is true, when the chromosome length isn't evenly divided into bin_size-sized bins, 
     if the leftover is less than merge_last_bin_cutoff, merge it into the last bin, otherwise add it as an extra bin 
      (this also ensures that if a chromosome is shorter than a bin, it'll be treated as single bin anyway).
     If normalize_last_bin is also True, normalize the number of positions in the special last bin by its size, 
      so that its value reflects the position density rather than raw count, to match all the other bins for heatmap display.

    If chromosomes is not None, only include those chromosomes; otherwise include all chromosomes in position_dict only.
    """
    # get chromosome list from position_dict - also take care of the case if position_dict keys are (chrom,strand) tuples
    if isinstance(position_dict.keys()[0], str): 
        chromosomes_in_position_dict = position_dict.keys()
    else:                                        
        chromosomes_in_position_dict = set(chrom for (chrom,strand) in position_dict.keys())
    if chromosomes is None:   
        chromosomes = chromosomes_in_position_dict
    # get chromosome lengths, make sure we have them for all chromosomes
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        chromosomes_no_lengths = set(chromosomes) - set(chromosome_lengths)
        if chromosomes_no_lengths:
            raise Exception("some chromosomes have no length data! %s"%chromosomes_no_lengths)
    chromosome_bin_count_lists = {}
    for chromosome in chromosomes:
        chromosome_length = chromosome_lengths[chromosome]
        # position_dict can be either just chromosome or (chromosome,strand)
        #  - get all the positions on the chromosome in either case, by just adding together the two strand lists if needed.
        #       (convert them to lists first, in case they're numpy arrays, for which + works differently!)
        #       (MAYBE-TODO this list conversion is probably inefficient - could do something with numpy array concatenation...)
        try:                position_list = list(position_dict[chromosome])
        except KeyError:    position_list = list(position_dict[(chromosome,'+')]) + list(position_dict[(chromosome,'-')])
        # divide the chromosome into bin_size-sized ranges, using an x.5 cutoff for clarity
        bin_edges = [x-.5 for x in range(first_bin_offset+1, chromosome_length+1, bin_size)]
        # Make sure there's at least one bin, even if the total length is smaller than a bin!  (for scaffolds/etc)
        # There'll be a smaller-than-bin_size chunk left over at the end (if length doesn't divide evenly into bin_size), 
        #  so add an extra bin for that if it's at least half a bin_size, otherwise make the last bin bigger to compensate.
        if special_last_bin:
            last_bin_edge = chromosome_length-.5
            if len(bin_edges)==1:                               bin_edges.append(last_bin_edge)
            if (last_bin_edge - bin_edges[-1])/bin_size < merge_last_bin_cutoff:  bin_edges[-1] = last_bin_edge
            else:                                                                 bin_edges.append(last_bin_edge)
        # use numpy.histogram to get the actual bin counts, IF we have any bins - otherwise the counts are an empty array
        #  (using numpy.array instead of list just to make sure the types match)
        if len(bin_edges)>1:    bin_count_list, _  = numpy.histogram(position_list, bin_edges)
        else:                   bin_count_list     = numpy.array([])
        # for the last bin in each chromosome, the value should be scaled by the smaller/bigger bin-size...
        if special_last_bin and normalize_last_bin and len(bin_count_list):
            bin_count_list[-1] = bin_count_list[-1] / ((bin_edges[-1] - bin_edges[-2])/bin_size)
        chromosome_bin_count_lists[chromosome] = bin_count_list
    return chromosome_bin_count_lists
    # TODO should probably unit-test this!


def get_20bp_fraction(dataset):
    """ Return the fraction of mutants in dataset that only have 20bp flanking regions (no 21bp ones).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).
    """
    max_lengths = [max([len(s) for s in mutant.sequences_and_counts.keys()]) for mutant in dataset]
    assert max_lengths.count(20) + max_lengths.count(21) == len(dataset), "Some sequences are outside the 20-21bp range!"
    return max_lengths.count(20)/len(dataset)
    # LATER-TODO should this be here, or a mutant_analysis_classes.Insertional_mutant_pool_dataset method or something?


def filter_data_by_min_reference(reference_list, other_lists, min_reference_value):
    other_data_tuples = zip(*other_lists)
    return [other_data for (x,other_data) in zip(reference_list, other_data_tuples) if x>=min_reference_value]
    # MAYBE-TODO should this be in general_utilities or something?


def replicate_reproducibility_info(dataset1, dataset2, readcount_cutoffs=[1,10,100,1000], ratio_cutoffs=[1.2,1.5,2,5,10], 
                                   both_ways=True, higher_only=False, real_min_reads_1=None, real_min_reads_2=None, quiet=False):
    """ Return/print reproducibility information between two replicates: % of readcounts within 2x of each other, etc. 

    Answers two basic questions, for each N in readcount_cutoffs:
        1) what percent of flanking sequences with N+ reads in one replicate were also observed in the other replicate?  
            (observed at all, not with N+ reads)
        2) what % of the mutants with N+ reads in one replicate have the two replicate readcounts within 2x/etc of each other? 
            Using each value in ratio_cutoffs for the ratio; using real or normalized readcount 
            (normalized is so that if one replicate has 2x more reads total than the other, the 2x ratio is counted as 1x etc; 
             still using raw readcounts for the cutoffs, because using normalized ones would probably be too hard to explain...)
           If higher_only is True, this only excludes readcounts that are 2x/etc HIGHER in dataset1 than in dataset1 (normalized), 
             rather than either higher or lower. 
        Also reports the total number of mutants with N+ reads: total, as a percentage of all mutants, and optionally 
         as a percentage of all "real" mutants, i.e. ones with at least real_min_reads_1 in dataset1 or real_min_reads_2 in dataset2 
          (if real_min_reads_* are both not None).
        Note that we're looking at mutants with >=X reads in one dataset, and putting no cutoff on the other; 
        If both_ways is False, we just do it one way (looking at mutants with N+ counts in A and checking presence/ratio in B); 
         if True, we do it both ways (A to B and B to A) and add the results together.
        Both_ways and higher_only cannot both be True, since the first one implies we're treating the two samples the same 
         and the second that we're treating them differently.

    For example, if we have two replicates A and B, with readcounts of (0,1,3,4) and (1,0,2,5) for the same 4 mutants:
        * looking at mutants with 1+ read (3 in A, 3 in B):
          - out of mutants with 1+ read in A, 2 are present in B; same the other way around - so total 67% present. 
          - out of mutants with 1+ read in A, the A:B ratios are 1:0, 3:2, 4:5 - so 1/3 within 1.3x, 2/3 within 2x (and 5x etc) 
             (and the reverse for mutants with 1+ read in B, leading to the same final counts)
            if we're looking higher_only, then 1/3 ratios are below 1.3x (and even 1x), 2/3 below 2x (and 5x etc)
            (or if we were looking at mutants with 1+ read in B, the B:A ratios are 1:0, 2:3, 5:4 - same except 2 below 1.3x.)
        * looking at mutants with 3+ reads (2 in A, 1 in B):
          - out of mutants with 3+ reads in A, 2 are present in B; out of 1 in B, 1 present in A - so total 100% present.
          - the ratios are 2:3 and 5:4 for the two in A, and 4:5 for the one in B, so 2/3 within 1.3x, 3/3 within 2x.
        (the ratios are the same normalized and unnormalized here; if A was (0,2,6,8) instead, then the normalized ratios from 
         the last section would be the same as above, but the raw ratios would be 2:6, 5:8, 8:5 - 0/3 within 1.3x, 2/3 within 2x.)
    
    The first two inputs should be either lists if integer readcounts (must be for the same list of mutants, in same order!), 
     or mutant_analysis_classes.Insertional_mutant_pool_dataset instances.

    Outputs, five total, in order:
        1-4: four dictionaries with all the readcount cutoffs as keys (readcount cutoff is N), and various values:
          1: number of mutants with N+ reads in one dataset
          2: number of mutants with N+ reads in one dataset AND present in the other dataset
          3-4: two dictionaries with ratio cutoffs as keys, and the values being the number of mutants 
             with N+ reads in one dataset that had a readcount within that ratio cutoff to the other dataset: 
              first dictionary uses raw readcounts to check ratios, the second uses readcounts normalized to total per dataset.
        5: a mutant_totals list: first element is the total number of mutants with non-zero readcounts, 
            second element is the total number of "real" mutants (i.e. ones with at least real_min_reads_1 in dataset1 
             or real_min_reads_2 in dataset2) - second element only present if real_min_reads_* are both not None.

    If quiet is True, don't print the data, just return it.  Note: printing is nicely formatted and has percentages.
    """
    # MAYBE-TODO add example with a lot more mutants in A than B, or such, to illustrate difference between A/B, B/A, and both_ways?
    # check if the arguments are consistent
    if both_ways and higher_only:
            raise ValueError("Both_ways and higher_only arguments cannot both be True!")
    # if the inputs are datasets instead of readcount lists, convert to readcount lists 
    #   (join the two datasets first so that both readcount lists are for the same mutant list)
    if isinstance(dataset1[0], int):
        readcounts1, readcounts2 = dataset1, dataset2
        if not len(readcounts1) == len(readcounts2):
            raise ValueError("First two arguments to replicate_reproducibility_info must be same length!")
    else:
        joint_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset=True)
        joint_dataset.populate_multi_dataset({'1': dataset1, '2': dataset2})
        readcounts1, readcounts2 = zip(*[[mutant.by_dataset[D].total_read_count for D in '12'] for mutant in joint_dataset])
    # make normalized data for the normalized ratios
    total_1, total_2 = sum(readcounts1), sum(readcounts2)
    norm_readcounts1, norm_readcounts2 = [x/total_1 for x in readcounts1], [x/total_2 for x in readcounts2]
    # we're adding together the comparisons both ways, then redefine the two readcount lists appropriately
    if both_ways:   readcountsA, readcountsB = readcounts1+readcounts2, readcounts2+readcounts1
    else:           readcountsA, readcountsB = readcounts1, readcounts2
    if both_ways:   norm_readcountsA, norm_readcountsB = norm_readcounts1+norm_readcounts2, norm_readcounts2+norm_readcounts1
    else:           norm_readcountsA, norm_readcountsB = norm_readcounts1, norm_readcounts2
    # calculate the total number of mutants that are present in dataset A
    total_N_mutants = sum(x>0 for x in readcountsA)
    # if we got real_min_reads_* data, use that to calculate a separate "real mutant" total and use it
    if real_min_reads_1 is not None and real_min_reads_2 is not None:
        if both_ways:   total_N_real_mutants = sum(x>=real_min_reads_1 for x in readcounts1)\
                                              +sum(x>=real_min_reads_2 for x in readcounts2)
        else:           total_N_real_mutants = sum(x>=real_min_reads_1 for x in readcounts1)
        mutant_count_totals = [total_N_mutants, total_N_real_mutants]
        total_descriptions = ['of all', 'of after-standard-cutoffs']
    else:
        mutant_count_totals = [total_N_mutants]
        total_descriptions = ['of all']
    # make a few dictionaries to return the data in
    N_filtered_dict, N_present_in_B_dict, N_within_ratio_raw_dict, N_within_ratio_norm_dict = {}, {}, {}, {}
    # for each readcount cutoff, filter the data by the cutoff in datasetA, then look at presence and ratio in datasetB:
    for readcount_cutoff in readcount_cutoffs:
        # filter data by minimum readcount in A
        filtered_readcount_pairs = filter_data_by_min_reference(readcountsA, [readcountsA,readcountsB], readcount_cutoff)
        filtered_norm_readcount_pairs = filter_data_by_min_reference(readcountsA, [norm_readcountsA,norm_readcountsB], 
                                                                     readcount_cutoff)
        N_filtered_in_A = len(filtered_readcount_pairs)
        N_filtered_dict[readcount_cutoff] = N_filtered_in_A
        if not quiet:
            print " * looking at mutants with at least %s reads in at least one dataset - %s"%(readcount_cutoff, 
                   general_utilities.value_and_percentages(N_filtered_in_A, mutant_count_totals, 
                                                           percentage_format_str='%.0f', words_for_percentages=total_descriptions))
        # see how many are present in B
        N_present_in_B = sum(b>0 for (a,b) in filtered_readcount_pairs)
        if not quiet:
            print "   - out of those, %s are present in the other."%(general_utilities.value_and_percentages(N_present_in_B, 
                                                                                                             [N_filtered_in_A]))
        N_present_in_B_dict[readcount_cutoff] = N_present_in_B
        # look at how many mutants are within given B:A readcount ratios (raw and normalized readcounts)
        if not quiet:
            if higher_only:     print '   - how many mutants have a readcount replicate1:2 ratio exceeding the given value:'
            else:               print '   - how many mutants are within a given readcount ratio between the replicates:'
        info_strings_raw, info_strings_norm = [], []
        N_within_ratio_raw_dict[readcount_cutoff], N_within_ratio_norm_dict[readcount_cutoff] = {}, {}
        for ratio_cutoff in ratio_cutoffs:
            # note that we're always taking a b/a ratio, not a/b, because we know a is non-zero, but we don't know that about b!
            if higher_only:
                # for the higher-only option we want to test a/b <= ratio_cutoff, i.e. 1/ratio_cutoff <= b/a
                N_within_ratio_raw = sum((1/ratio_cutoff <= b/a) for a,b in filtered_readcount_pairs)
                N_within_ratio_norm = sum((1/ratio_cutoff <= b/a) for a,b in filtered_norm_readcount_pairs)
            else:
                N_within_ratio_raw = sum((1/ratio_cutoff <= b/a <= ratio_cutoff) for a,b in filtered_readcount_pairs)
                N_within_ratio_norm = sum((1/ratio_cutoff <= b/a <= ratio_cutoff) for a,b in filtered_norm_readcount_pairs)
            N_within_ratio_raw_dict[readcount_cutoff][ratio_cutoff] = N_within_ratio_raw
            N_within_ratio_norm_dict[readcount_cutoff][ratio_cutoff] = N_within_ratio_norm
            info_strings_raw.append('%sx - %s'%(ratio_cutoff, 
                                                'N/A' if N_filtered_in_A==0 else '%.0f%%'%(N_within_ratio_raw/N_filtered_in_A*100)))
            info_strings_norm.append('%sx - %s'%(ratio_cutoff, 
                                                'N/A' if N_filtered_in_A==0 else '%.0f%%'%(N_within_ratio_norm/N_filtered_in_A*100)))
        if not quiet:
            print '       raw readcount ratios       : '+', '.join(info_strings_raw)
            print '       normalized readcount ratios: '+', '.join(info_strings_norm)
    return N_filtered_dict, N_present_in_B_dict, N_within_ratio_raw_dict, N_within_ratio_norm_dict, mutant_count_totals


def insertion_pair_type_distance(pos1, side1, pos2, side2):
    """ Given two insertion positions and sides (5'/3'), return (type, distance), where type is inner/outer-cassette or same-dir.

    pos1 and pos2 should be Insertion_position instances (from mutant_IB_RISCC_classes or mutant_analysis_classes),
      or they can be (chromosome, strand, min_position) tuples. (Strand must be + or -, chromosome is a string, min_pos is an int.)
    side1 and side2 should be 5' or 3' (they can be different or the same).

    If the positions are on different chromosomes, the type will be 'diff-chrom' and the distance will be NaN.
    Otherwise, the distance will be the distance between the minimum_position property of the two positions, and 
    the type shows how the two flanking sequences would relate to each other, as follows:

    Inner-cassette would be two flanking seqs in any arrangement like this: 
      (| is the cassette position, ----> is the flanking sequence; the distance between the two |s is the distance.)
        <-------|               <-------|                 <-------|           
                |------>                  |------>                      |------>  
      The first case would be distance 0, like two sides of a perfect cassette insertion (both +strand but different ends, 
        or same end but different strands if it was an opposite-tandem cassette insertion).  The other cases could be
        two sides of a cassette insertion with a genomic deletion.
    Outer-cassette pairs would be any of those (at different distances, high to low):
                     <-------|             <-------|          <-------|         <-------|            
        |------>                      |------>               |------>                |------>              
    The minimum distance is 1 - a distance of 0 would be considered inner-cassette. If the distance is low, those could be 
      two sides of an insertion cassette with a short duplication; if it's long, it could be two sides of a fragment between
      two insertion cassettes.
    Same-direction pairs would be any of these (different distances, low to high):
        |------>               |------>                 |------>              
        |------>                    |------>                         |------>
      These are likely to be independent rather than two sides of the same insertion. They could also possibly be two sides
      of an insertion with a genomic inversion on one side - we've seen one like this in Ru's first screen.
    """
    try:                    chrom1, strand1, minpos1 = pos1.chromosome, pos1.strand, pos1.min_position
    except AttributeError:  chrom1, strand1, minpos1 = pos1
    try:                    chrom2, strand2, minpos2 = pos2.chromosome, pos2.strand, pos2.min_position
    except AttributeError:  chrom2, strand2, minpos2 = pos2
    # check that sides are valid, and truncate them to first character so that 5' and 5prime are the same
    allowed_sides = basic_seq_utilities.SEQ_ENDS + basic_seq_utilities.SEQ_ENDS_SHORT
    if side1 not in allowed_sides:     raise Exception("Invalid side %s - must be one of %s!"%(side1, ', '.join(allowed_sides)))
    if side2 not in allowed_sides:     raise Exception("Invalid side %s - must be one of %s!"%(side2, ', '.join(allowed_sides)))
    side1, side2 = side1[0], side2[0]
    # figure out the case
    if chrom1 != chrom2:                                    return ('diff-chrom', float('NaN'))
    dist = abs(minpos1 - minpos2)
    if side1==side2 and strand1==strand2:                   return ('same-direction', dist)
    if side1!=side2 and strand1!=strand2:                   return ('same-direction', dist)
    if '5' in side1 and strand1=='+' and minpos1 <= minpos2: return ('inner-cassette', dist)
    if '3' in side1 and strand1=='-' and minpos1 <= minpos2: return ('inner-cassette', dist)
    if '5' in side2 and strand2=='+' and minpos2 <= minpos1: return ('inner-cassette', dist)
    if '3' in side2 and strand2=='-' and minpos2 <= minpos1: return ('inner-cassette', dist)
    else:                                                   return ('outer-cassette', dist)

######### unit-tests

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__filter_data_by_min_reference(self):
        for min_val in range(-10,1):
            assert filter_data_by_min_reference([1,3,2], [[3,2,1],[4,5,4]], min_val) == [(3,4), (2,5), (1,4)]
        if True:    # this if statement is only here to make things indent and line up nicely
            assert filter_data_by_min_reference([1,3,2], [[3,2,1],[4,5,4]], 1)       == [(3,4), (2,5), (1,4)]
            assert filter_data_by_min_reference([1,3,2], [[3,2,1],[4,5,4]], 2)       == [       (2,5), (1,4)]
            assert filter_data_by_min_reference([1,3,2], [[3,2,1],[4,5,4]], 3)       == [       (2,5)       ]
        for min_val in range(4,10):
            assert filter_data_by_min_reference([1,3,2], [[3,2,1],[4,5,4]], min_val) == [                   ]

    def test__replicate_reproducibility_info(self):
        # the outputs are: N_filtered_dict, N_present_in_B_dict, N_within_ratio_raw_dict, N_within_ratio_norm_dict
        ### BASIC TEST, based on docstring example.
        # checking both one-way options, and the both-way option (which should just be the sum of the two one-way ones)
        # In all of those, the two datasets have the same readcount sum, so the raw and normalized data should be the same
        #   (last two outputs), so we're only checking output 0-2, and checking that output3 is the same as output2.
        dataA, dataB = [0,1,3,4], [1,0,2,5]
        data1 = replicate_reproducibility_info(dataA, dataB, readcount_cutoffs=[1, 3], ratio_cutoffs=[1.3, 2], 
                                               both_ways=False, quiet=True)
        assert data1[:3] == ({1:3, 3:2}, {1:2, 3:2}, {1:{1.3:1, 2:2}, 3:{1.3:1, 2:2}})
        assert data1[3] == data1[2]
        data2 = replicate_reproducibility_info(dataB, dataA, readcount_cutoffs=[1, 3], ratio_cutoffs=[1.3, 2], 
                                               both_ways=False, quiet=True)
        assert data2[3] == data2[2]
        assert data2[:3] == ({1:3, 3:1}, {1:2, 3:1}, {1:{1.3:1, 2:2}, 3:{1.3:1, 2:1}})
        data3 = replicate_reproducibility_info(dataA, dataB, readcount_cutoffs=[1, 3], ratio_cutoffs=[1.3, 2], 
                                               both_ways=True, quiet=True)
        assert data3[:3] == ({1:6, 3:3}, {1:4, 3:3}, {1:{1.3:2, 2:4}, 3:{1.3:2, 2:3}})
        assert data3[3] == data3[2]
        # now make sure that doing the both-way option by hand (in both directions) gives the same result
        data4 = replicate_reproducibility_info(dataB+dataA, dataA+dataB, readcount_cutoffs=[1, 3], ratio_cutoffs=[1.3, 2], 
                                               both_ways=False, quiet=True)
        data5 = replicate_reproducibility_info(dataA+dataB, dataB+dataA, readcount_cutoffs=[1, 3], ratio_cutoffs=[1.3, 2], 
                                               both_ways=False, quiet=True)
        assert data3 == data4 == data5
        ### NORMALIZATION TEST, again based on docstring example.
        # the normalized ratio counts should be the same as above, but the raw ones should be different.
        data = replicate_reproducibility_info([x*2 for x in dataA], dataB, readcount_cutoffs=[3], ratio_cutoffs=[1.3, 2], 
                                              both_ways=False, quiet=True)
        assert data[:4] == ({3:2}, {3:2}, {3:{1.3:0, 2:1}}, {3:{1.3:1, 2:2}})
        data = replicate_reproducibility_info(dataB, [x*2 for x in dataA], readcount_cutoffs=[3], ratio_cutoffs=[1.3, 2], 
                                              both_ways=False, quiet=True)
        assert data[:4] == ({3:1}, {3:1}, {3:{1.3:0, 2:1}}, {3:{1.3:1, 2:1}})
        data = replicate_reproducibility_info([x*2 for x in dataA], dataB, readcount_cutoffs=[3], ratio_cutoffs=[1.3, 2], 
                                              both_ways=True, quiet=True)
        assert data[:4] == ({3:3}, {3:3}, {3:{1.3:0, 2:2}}, {3:{1.3:2, 2:3}})
        ### HIGHER_ONLY TEST, again based on docstring example, trying A>B and B>A
        # only using one readcount cutoff for simplicity; normalized data should be same as raw, since totals are same
        data1 = replicate_reproducibility_info(dataA, dataB, readcount_cutoffs=[1], ratio_cutoffs=[1, 1.3, 2, 5], 
                                               both_ways=False, higher_only=True, quiet=True)
        assert data1[:3] == ({1:3}, {1:2}, {1:{1:1, 1.3:1, 2:2, 5:2}})
        assert data1[3] == data1[2]
        data2 = replicate_reproducibility_info(dataB, dataA, readcount_cutoffs=[1], ratio_cutoffs=[1, 1.3, 2, 5], 
                                               both_ways=False, higher_only=True, quiet=True)
        assert data2[:3] == ({1:3}, {1:2}, {1:{1:1, 1.3:2, 2:2, 5:2}})
        assert data2[3] == data2[2]
        ### MUTANT-TOTAL AND REAL-MUTANT-TOTAL TEST
        # Just checking the last output (which I haven't checked before) - the mutant total and "real" total list.
        # Ignoring everything else (leaving defaults for readcount/ratio cutoffs)
        dataA, dataB = [0,1,3,4], [1,0,2,5]
        real_min_reads_and_expected_output = [(None, None, [6]), 
                                              (1,    1,    [6,6]), 
                                              (1,    2,    [6,5]), 
                                              (2,    1,    [6,5]), 
                                              (2,    2,    [6,4]), 
                                              (2,    3,    [6,3]), 
                                              (3,    2,    [6,4]), 
                                              (4,    5,    [6,2]), 
                                              (5,    4,    [6,1]), 
                                              (6,    6,    [6,0]), 
                                             ]
        for real_min_1, real_min_2, expected_output in real_min_reads_and_expected_output:
            data = replicate_reproducibility_info(dataA, dataB, both_ways=True, 
                                                  real_min_reads_1=real_min_1, real_min_reads_2=real_min_2, quiet=True)
            assert data[4] == expected_output
        ### NO-CUTOFFS TEST
        # if cutoff lists are empty, we get empty dictionaries, no errors
        data = replicate_reproducibility_info(dataA, dataB, readcount_cutoffs=[], ratio_cutoffs=[], both_ways=True, quiet=True)
        assert data == ({}, {}, {}, {}, [6])
        ### DIVISION BY ZERO TEST
        # checks that we don't get an error raised when readcount cutoff is so high that nothing's left
        data = replicate_reproducibility_info(dataA, dataB, readcount_cutoffs=[100], ratio_cutoffs=[1,2], both_ways=True, quiet=True)
        ### WRONG-ARGUMENTS TEST
        # higher_only and both_ways cannot both be True
        self.assertRaises(ValueError, replicate_reproducibility_info, dataA, dataB, both_ways=True, higher_only=True)
        # the two datasets must be the same length, if given as readcount lists
        self.assertRaises(ValueError, replicate_reproducibility_info, dataA+[3], dataB, both_ways=True, higher_only=True)
        # MAYBE-TODO test the actual printed data?  It's pretty straightforwardly derived from the outputted/tested data, though.

    def test__insertion_pair_type_distance(self):
        side0, pos0 = "5'", mutant_IB_RISCC_classes.Insertion_position('chr1', '+', full_position='100-?')
        assert insertion_pair_type_distance(pos0, side0, pos0, side0) == ('same-direction', 0)
        # other-side cases, same chromosome
        for (side1, strand1) in [("3'", "+"), ("5'", "-")]:
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='?-101')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('inner-cassette', 0) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='?-102')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('inner-cassette', 1) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='?-301')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('inner-cassette', 200) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='?-100')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('outer-cassette', 1) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='?-50')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('outer-cassette', 51) 
        for (side1, strand1) in [("3'", "-"), ("5'", "+")]:
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='100-?')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('same-direction', 0) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='101-?')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('same-direction', 1) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='50-?')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('same-direction', 50) 
            pos1 = mutant_IB_RISCC_classes.Insertion_position('chr1', strand1, full_position='500-?')
            assert insertion_pair_type_distance(pos0, side0, pos1, side1) == ('same-direction', 400) 
        for side1 in "5' 3'".split():
            for chrom1 in 'chr2 chromosome_5 insertion_cassette scaffold_88 chloroplast'.split():
                for strand1 in "+=":
                    for fullpos1 in '?-1 ?-100 ?-10000 1-? 100-? 10000-?'.split():
                        pos1 = mutant_IB_RISCC_classes.Insertion_position(chrom1, strand1, full_position=fullpos1)
                        category, dist = insertion_pair_type_distance(pos0, side0, pos1, side1)
                        assert category == 'diff-chrom'
                        assert numpy.isnan(dist)

    # LATER-TODO add unit-tests for the remaining functions!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
