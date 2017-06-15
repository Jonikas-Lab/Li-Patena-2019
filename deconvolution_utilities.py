#! /usr/bin/env python

"""
Basic functions for combinatorial deconvolution of mutant pools, i.e. matching observed deepseq readcounts for each insertion in each pool to the expected presence/absence codewords of each original sample in each pool, to determine which sequenced insertion corresponds to which original sample.

Basic experimental pipeline (and nomenclature):
 1) Biological SAMPLES are pooled into POOLS, based on predetermined CODEWORDS, which are then used as "EXPECTED" CODEWORDS.
    (sometimes the "samples" are already pools themselves, and the "pools" are called superpools, but that doesn't matter here)
    This is all done in the separate combinatorial_pooling package.
 2) After each POOL is sequenced, each INSERTION found in any pool will have an "OBSERVED" CODEWORD based on 
    which pools it was present in.  Those "observed" codewords are then matched to the closest "expected" codeword, 
    to determine which INSERTION corresponds to which original SAMPLE - this is what this module is for. 

 -- Weronika Patena, 2013
"""

# standard library
from __future__ import division
import os
import collections
import unittest
import random
import math
from collections import defaultdict
from string import ascii_uppercase # this is 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', for well numbers
# other packages
import numpy
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import general_utilities, plotting_utilities
import binary_code_utilities
import mutant_analysis_classes
import mutant_utilities

class DeconvolutionError(Exception):
    """ Exceptions in the deconvolution_utilities module."""
    pass


################################################## Various utilities #######################################################

def get_original_384_well_numbers(sample_number, zero_padding=True):
    """ Get a well name (A1-P24) form a sample number (1-384) assuming the 384-well plate was split into 4 staggered 96-well plates.
    
    It's a complicated setup!  The assumption is that the original 384-well plate was split into four staggered 96-well plates, 
     with every other row/column being in the same plate, like this:
                col1:       col2:       col3:       col4:       col5:       ...     col23:      col24:
        rowA:   plate1-A1   plate2-A1   plate1-A2   plate2-A2   plate1-A3   ...     plate1-A12  plate2-A12
        rowB:   plate3-A1   plate4-A1   plate3-A2   plate4-A2   plate3-A3   ...     plate3-A12  plate4-A12
        rowC:   plate1-B1   plate2-B1   plate1-B2   plate2-B2   plate1-B3   ...     plate1-B12  plate2-B12
                ...         ...         ...         ...         ...         ...     ...         ...
        rowO:   plate1-H1   plate2-H1   plate1-H2   plate2-H2   plate1-H3   ...     plate1-H12  plate2-H12
        rowP:   plate3-H1   plate4-H1   plate3-H2   plate4-H2   plate3-H3   ...     plate3-H12  plate4-H12

     And then that the resulting four 96-well plates were taken in order, with each well in order, 
      so the final sample order would be like this: plate1-A1, plate1-A2, ..., plate2-A1, ..., plate4-H12.

     So the full conversion from final number to original 384-well number is like this:
         - sample 1 would be plate1-A1, i.e. original well A1
         - sample 2 would be plate1-A2, i.e. original well A3
         - sample 3 would be plate1-A3, i.e. original well A5
            ...
         - sample 13 would be plate1-B1, i.e. original well C1
         - sample 14 would be plate1-B2, i.e. original well C3
            ...
         - sample 95 would be plate1-H11, i.e. original well O21
         - sample 96 would be plate1-H12, i.e. original well O23
            ...
         - sample 97 would be plate2-A1, i.e. original well A2
         - sample 98 would be plate2-A2, i.e. original well A4
            ...
         - sample 383 would be plate4-H11, i.e. original well P22
         - sample 384 would be plate4-H12, i.e. original well P24
    """
    # first convert everything to 0-based - it makes all the math MUCH simpler
    sample_number -= 1
    # first 96 are plate1, then plate2, plate3, plate4.
    plate_number_96plate = sample_number // 96
    well_number_96plate =  sample_number %  96
    # split well number into row/column number (96-well plates are 8x12)
    row_number_96plate = well_number_96plate // 12
    col_number_96plate = well_number_96plate %  12
    # now get the original 384-well rows/columns, based on the staggering
    if_second_row = plate_number_96plate > 1
    if_second_col = plate_number_96plate % 2
    row_number_384plate = row_number_96plate*2 + if_second_row
    col_number_384plate = col_number_96plate*2 + if_second_col
    assert 0 <= row_number_384plate < 16
    assert 0 <= col_number_384plate < 24
    return '%s%02d'%(ascii_uppercase[row_number_384plate], col_number_384plate+1)


############################################ Conversion from readcounts to codewords ############################################

# MAYBE-TODO should this be in mutant_utilities or something? Not really deconvolution-specific...
def datasets_to_readcount_table(insertion_pool_joint_dataset, dataset_names_in_order, use_perfect_reads=False):
    """ Given an insertional pool multi-dataset, return a insertion_position:readcount_list dict and a list of datasets in order.

    Insertion_pool_joint_dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, multi-dataset, 
     with each member dataset corresponding to a pool.
    dataset_names_in_order should be a list of the dataset names in the order in which they should be listed in the outputs.

    The two outputs are:
     1) a insertion_position:readcount_list, with keys being the positions of all the insertions in the original multi-dataset, 
        and the values being lists of readcounts for that insertion position in each single dataset, 
        in the order given by dataset_names_in_order.

    If use_perfect_reads, the values in the final table will be perfect read counts; otherwise total read counts.
    """
    # MAYBE-TODO implement option for taking a list of datasets rather than a single joint dataset, too?  (would then also have to take a list of names, which I guess is the same as dataset_names_in_order)
    insertion_readcount_table = {}
    for insertion in insertion_pool_joint_dataset:
        if use_perfect_reads:   readcounts = [insertion.by_dataset[dataset].perfect_read_count for dataset in dataset_names_in_order]
        else:                   readcounts = [insertion.by_dataset[dataset].total_read_count for dataset in dataset_names_in_order]
        insertion_readcount_table[insertion.position] = readcounts
    return insertion_readcount_table


# MAYBE-TODO should this be in mutant_utilities or something? Not really deconvolution-specific...
def normalize_readcount_table(insertion_readcount_table, multiplier=10**6):
    """ Given an insertion:dataset_readcount_list dict, normalize the readcounts to the totals for each dataset, IN PLACE.

    Input should be same as the output from datasets_to_readcount_table.  Input will be modified in-place.
    """
    dataset_readcount_totals = [sum(all_readcounts) for all_readcounts in zip(*insertion_readcount_table.values())]
    for insertion_pos, dataset_readcounts in list(insertion_readcount_table.items()):
        insertion_readcount_table[insertion_pos] = [readcount/total * multiplier
                                                   for readcount,total in zip(dataset_readcounts,dataset_readcount_totals)]
    # doesn't return anything - modifies the input in place.
    # TODO unit-test!


def readcounts_to_presence__cutoffs(readcount_list, cutoffs):
    """ Given a readcount list and a cutoff list (or one cutoff), return True/False list for whether readcount>=cutoff.
    """
    # if cutoffs is a single value, make it a list repeating that single value
    try:                len(cutoffs)
    except TypeError:   cutoffs = [cutoffs for _ in range(len(readcount_list))]
    if not len(cutoffs)==len(readcount_list):
        raise DeconvolutionError("readcount_list and cutoffs must be the same length!")
    return [readcount>=cutoff for readcount,cutoff in zip(readcount_list, cutoffs)]


def readcounts_to_presence__mutant_1(readcount_list, N_always_present, N_always_absent, overall_min=2, 
                                     present_level_function=numpy.median, absent_level_function=numpy.median, 
                                     cutoff_position=0.5, min_cutoff=1):
    """ Given a readcount list for a mutant, use a complex method to decide which readcounts indicate presence/absence.

    ___
    """
    # TODO give more details in docstring!
    # TODO should cutoff_positon be a function instead?
    if N_always_present+N_always_absent > len(readcount_list):
        raise DeconvolutionError("N_always_present+N_always_absent should never be more than len(readcount_list)!")
    if N_always_present<1 or N_always_absent<0:
        raise DeconvolutionError("N_always_present must be >0, and N_always_absent >=0!")
    sorted_readcounts = sorted(readcount_list)
    if N_always_absent == 0:
        absent_level = 0
    else:
        absent_readcounts = sorted_readcounts[:N_always_absent]
        absent_level = absent_level_function(absent_readcounts)
    present_readcounts = sorted_readcounts[-N_always_present:]
    present_level = present_level_function(present_readcounts)
    # TODO how exactly should the overall_min be applied?  
    #  to throw away entire codewords as all-0, should we be comparing it to present_level, or to max?
    if present_level < overall_min:
        return [False for readcount in readcount_list]
    if absent_level >= present_level:
        raise DeconvolutionError("Absent level >= present level - should never happen! (%s and %s, based on %s and %s)"%(
            absent_level, present_level, ([] if N_always_absent==0 else absent_readcounts), present_readcounts))
    cutoff = absent_level + cutoff_position*(present_level-absent_level)
    # adjust cutoff
    cutoff = max(cutoff, min_cutoff)
    # MAYBE-TODO what if cutoff ends up below max(absent_readcounts), or above min(present_readcounts)?  Fix?
    # MAYBE-TODO should there be some minimum difference between absent_level and present_level?
    return [readcount>=cutoff for readcount in readcount_list]


def readcounts_to_codewords(insertion_readcount_table, conversion_function):
    """ Given insertion:reads-per-pool data, return insertion:codeword dict showing which insertion was present in which pool. 

    Insertion_readcount_table should be in the same format as the output of datasets_to_readcount_table.
    Conversion_function should be a function that takes a readcount list and returns a True/False list indicating 
     which of those readcounts are considered present/absent (choose one of the readcounts_to_presence__* functions above, 
      and do a lambda to set the other arguments).  

    The output is an insertion_position:presence_codeword dictionary with the same shape as the first input, 
     with presence_codeword a binary_code_utilities.Binary_codeword 
     representing whether a given insertion is considered present in each pool.
    """
    insertion_codeword_dict = {}
    for insertion_pos, readcount_list in insertion_readcount_table.items():
        presence_list = conversion_function(readcount_list)
        insertion_codeword_dict[insertion_pos] = binary_code_utilities.Binary_codeword(presence_list)
    return insertion_codeword_dict


########################################## Matching observed to expected codewords ##############################################

def read_codewords_from_file(infile_name, new_sample_names=None):
    """ Read sample codewords used for combinatorial pooling, from file generated by robotic_plate_transfer.write_data_to_outfile.

    (robotic_plate_transfer is in the ../../combinatorial_pooling/code folder.)
    Only reads the first "sample_number plate_and_well_position codeword transfers volume" table in the infile, ignores the rest.

    Returns sample_number:codeword dictionary, where codewords are binary_code_utilities.Binary_codeword objects.

    If new_sample_names isn't None, it should be an old_number:new_name dictionary - the new names will then be used in the output.
    """
    if new_sample_names is not None:
        if len(set(new_sample_names.values())) != len(new_sample_names):
            raise DeconvolutionError("All values in the new_sample_names dict must be unique!")
    sample_codewords = {}
    inside_relevant_section = False
    for line in open(infile_name):
        # skip lines until inside the section we want to parse
        if not inside_relevant_section:
            if line.startswith('sample_number\tplate_and_well_position\tcodeword'):
                inside_relevant_section=True
                continue
        # inside the section we want to parse:
        if inside_relevant_section:
            # an empty line or another header means the end of the section - just finish the loop.
            if (not line.strip()) or line.startswith('sample_number'):
                break
            # for each line in the section, just grab sample name and codeword from each line and put in the data dict
            fields = line.strip('\n').split('\t')
            try:                sample_name, codeword_string = fields[0], fields[2]
            except IndexError:  raise DeconvolutionError("Cannot parse line! \"%s\""%line)
            # optionally convert original to new sample names
            if new_sample_names is not None:
                sample_name = new_sample_names[sample_name]
            sample_codewords[sample_name] = binary_code_utilities.Binary_codeword(codeword_string)
    return sample_codewords


def find_closest_sample_codeword(ref_codeword, sample_codewords, min_distance_difference=1, ignore_diff_lengths=False):
    """ Returns the sample with the closest codeword match to the ref_codeword, and the Hamming distance.

    Inputs: a 0/1 string codeword (a binary_code_utilities.Binary_codeword instance, or just a string), 
    and a sample:codeword dictionary containing more such codewords (values must be unique).

    Finds the sample(s) with codeword(s) with the lowest and second-lowest Hamming distance to the ref_codeword:
        - if there is a single sample with the lowest distance, and the difference between the lowest and second-lowest distances 
            is <=min_distance_difference, return the (closest_sample_name, min_Hamming_distance) tuple for the lowest.
        - otherwise, return a (None, min_Hamming_distance) tuple.
    """
    # TODO add docstring and unit-test for ignore_diff_lengths!  Ignore codeword with length mismatches.
    # useful stuff in binary_code_utilities: Binary_codeword object, Hamming_distance, bit_change_count (separate 0->1 and 1->0)
    # make sure codewords are unique
    if sample_codewords is not None:
        if len(set(sample_codewords.values())) != len(sample_codewords):
            raise DeconvolutionError("All values in the sample_codewords dict must be unique!")
    if min_distance_difference<1:
        raise DeconvolutionError("min_distance_difference must be an integer 1 or higher!")
    # LATER-TODO should probably rename Hamming_distance to be lowercase, since it's a function...
    if not isinstance(ref_codeword, binary_code_utilities.Binary_codeword):
        ref_codeword = binary_code_utilities.Binary_codeword(ref_codeword)
    # Just calculate the Hamming distance to all expected codewords
    #  MAYBE-TODO could probably optimize this a lot!  If only with caching...
    if ignore_diff_lengths:
        sample_to_distance = {sample:binary_code_utilities.Hamming_distance(ref_codeword, codeword) 
                                   for sample,codeword in sample_codewords.items() if len(codeword)==len(ref_codeword)}
    else:
        sample_to_distance = {sample:binary_code_utilities.Hamming_distance(ref_codeword, codeword) 
                                   for sample,codeword in sample_codewords.items()}
    min_distance = min(sample_to_distance.values())
    low_dist_samples = [sample for sample,distance in sample_to_distance.items() 
                        if distance < min_distance+min_distance_difference]
    if len(low_dist_samples)==1:    return (low_dist_samples[0], min_distance) 
    else:                           return (None, min_distance) 
    # TODO is this what I actually want to return, or something else?...  Could optionally return the full sorted distance_to_N_samples, or the top 2 of that, or something...


def match_insertions_to_samples(insertion_codewords, sample_codewords, min_distance_difference=1):
    """ Given observed insertion codewords and expected sample codewords, find best sample codeword match for each insertion.

    First two inputs should be insertion_pos:codeword and sample_name:codeword dictionaries, 
     with codewords beint binary_code_utilities.Binary_codeword instances.
    The find_closest_sample_codeword function is used to match insertion to sample codewords; 
     the min_distance_difference arg should be a number - see find_closest_sample_codeword for how it works. 

    The outputs are an insertion_pos:sample_name (or None if no match) and an insertion_pos:codeword_distance dictionary.
    """
    # TODO add option for how many errors should be allowed at all, rather than just always returning the closest unique match, even if it has 4 errors!  Maybe allow splitting that into 0->1 and 1->0 errors separately.
    insertion_samples = {}
    insertion_codeword_distances = {}
    for insertion_pos, insertion_codeword in insertion_codewords.items():
        best_match_sample, min_distance = find_closest_sample_codeword(insertion_codeword, sample_codewords, min_distance_difference)
        insertion_samples[insertion_pos] = best_match_sample
        insertion_codeword_distances[insertion_pos] = min_distance
    return insertion_samples, insertion_codeword_distances
    # TODO unit-test!


############################################### General deconvolution functions ##################################################

def combinatorial_deconvolution(insertion_pool_joint_dataset, sample_codeword_filename, pool_names_in_order, conversion_function, 
                                min_distance_difference=1, new_sample_names=None, 
                                normalize_pool_readcounts=True, normalization_multiplier=10**6, use_perfect_reads=False, 
                                return_readcounts=False):
    """ ___
    """
    # TODO write docstring! Avoid repetition though.
    # TODO wait, what's the relationship between pool_names_in_order and new_sample_names??
    sample_codewords = read_codewords_from_file(sample_codeword_filename, new_sample_names)
    insertion_readcount_table = datasets_to_readcount_table(insertion_pool_joint_dataset, pool_names_in_order, use_perfect_reads)
    if normalize_pool_readcounts: normalize_readcount_table(insertion_readcount_table, normalization_multiplier)
    insertion_codewords = readcounts_to_codewords(insertion_readcount_table, conversion_function)
    insertion_samples, insertion_codeword_distances = match_insertions_to_samples(insertion_codewords, sample_codewords, 
                                                                                  min_distance_difference)
    # MAYBE-TODO add more options for what is and isn't returned?  Or are those a bad idea?
    return_data = (insertion_codewords, insertion_samples, insertion_codeword_distances)
    if return_readcounts:   return_data = (insertion_readcount_table,) + return_data
    return return_data
    # TODO unit-test!


def get_deconvolution_summary(insertion_samples, insertion_codeword_distances):
    """ Given the outputs of match_insertions_to_samples or combinatorial_deconvolution, generate some summary numbers.

    Inputs are insertion_pos:sample_name and insertion_pos:min_codeword_distance dictionaries; 
        sample_name should be None if multiple samples matched.

    Outputs are the total numbers of matched and unmatched insertions, 
     and min_distance:N_insertions dicts for matched and unmatched insertions.
    """
    matched_insertion_distances = {i:d for i,d in insertion_codeword_distances.items() if insertion_samples[i] is not None}
    unmatched_insertion_distances = {i:d for i,d in insertion_codeword_distances.items() if insertion_samples[i] is None}
    N_total_insertions = len(insertion_samples)
    N_matched_insertions = len(matched_insertion_distances)
    N_unmatched_insertions = len(unmatched_insertion_distances)
    matched_insertion_counts_by_distance = collections.Counter(matched_insertion_distances.values())
    unmatched_insertion_counts_by_distance = collections.Counter(unmatched_insertion_distances.values())
    return N_matched_insertions, N_unmatched_insertions, matched_insertion_counts_by_distance, unmatched_insertion_counts_by_distance
    # TODO unit-test


def print_deconvolution_summary(description, N_matched_insertions, N_unmatched_insertions, matched_insertion_counts_by_distance, 
                                unmatched_insertion_counts_by_distance):
    """ Nicely print the data generated by get_deconvolution_summary.
    """
    # MAYBE-TODO add a collapse_higher_than arg, so that "All insertion counts with min distances >collapse_higher_than will be counted together instead of separately."? (put that bit in docstring)
    total_insertions = N_matched_insertions + N_unmatched_insertions
    print "%s: %s total insertions, %s uniquely matched to a sample, %s matched to 2+ samples (unmatched)."%(description, 
                total_insertions, general_utilities.value_and_percentages(N_matched_insertions, [total_insertions]), 
                general_utilities.value_and_percentages(N_unmatched_insertions, [total_insertions]))
    print " * matched insertions by codeword distance (% of matched, % of all): "
    print "    %s"%(', '.join(['%s: %s'%(d,general_utilities.value_and_percentages(n, [N_matched_insertions, total_insertions]))
                                         for d,n in sorted(matched_insertion_counts_by_distance.items())]))
    print " * unmatched insertions by codeword distance (% of unmatched, % of all): "
    print "    %s"%(', '.join(['%s: %s'%(d,general_utilities.value_and_percentages(n, [N_unmatched_insertions, total_insertions]))
                                         for d,n in sorted(unmatched_insertion_counts_by_distance.items())]))


def per_mutant_deconvolution(joint_dataset, dataset_name, pooling_scheme_file, sample_and_pool_names, N_present, N_absent, 
                             cutoff_position, overall_min=2, short=False):
    """ Convenience function for deconvolution and summary/printing. 
    Uses per-mutant cutoffs WITH normalization, assumes we're always taking the max absent and min present value 
     out of however many we're given; can print full info or a short tab-separated line. 
     """
    pool_names, sample_numbers_to_names = sample_and_pool_names
    conversion_function = lambda x: readcounts_to_presence__mutant_1(x, N_present, N_absent, overall_min, 
                                                                     min, max, cutoff_position, overall_min)
    output = combinatorial_deconvolution(joint_dataset, pooling_scheme_file, pool_names, conversion_function, 
                                         new_sample_names=sample_numbers_to_names, normalize_pool_readcounts=True, 
                                         return_readcounts=True)
    insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances = output
    summary = get_deconvolution_summary(insertion_samples, codeword_distances)
    if short:
        matched, unmatched, matched_by_dist, _ = summary
        total = matched + unmatched
        perfect_matched = matched_by_dist[0]
        good_matched = matched_by_dist[0] + matched_by_dist[1]
        okay_matched = matched_by_dist[0] + matched_by_dist[1] + matched_by_dist[2]
        print '\t'.join([str(N_present), str(N_absent), str(cutoff_position), str(overall_min), dataset_name, str(matched), 
                         general_utilities.value_and_percentages(perfect_matched, [matched]), 
                         general_utilities.value_and_percentages(good_matched, [matched]), 
                         general_utilities.value_and_percentages(okay_matched, [matched])])
    else:
        print_deconvolution_summary("%s data, per mutant (max top %s, min bottom %s, cutoff pos %s, min %s)"%(dataset_name, 
                                                                  N_present, N_absent, cutoff_position, overall_min), *summary)
    return insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances, summary


def pick_best_parameters(deconv_data_by_parameters, N_errors_allowed, N_good_factor, percent_good_factor, max_N_high=None):
    """ Pick the best deconvolution parameters, depending on what aspect of the result we want to optimize.

    The deconv_data_by_parameters arg should be a dict with (N_high, N_low, cutoff_pos, overall_min) tuples as keys and 
     (insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances, summary) tuples as values.
        (The first four should be outputs from combinatorial_deconvolution with return_readcounts=True, or the equivalent;
         summary should be the output from get_deconvolution_summary.)

    Look at the parameters/results in deconv_data_by_parameters.  
    Pick the best set of parameters, optimizing a combination of two factors:
        1) the total number of insertion positions mapped to expected codewords with N_errors_allowed or fewer errors
        2) the percentage of those out of all uniquely mapped insertion positions.
    These two numbers will be multiplied by N_good_factor and percent_good_factor respectively, 
     and the resulting combination will be  maximized.

    The output will be the deconv_data_by_parameters key that gives the optimal result.
    """
    parameters_to_result = {}
    for parameters, deconv_data in deconv_data_by_parameters.items():
        N_matched, N_unmatched, N_matched_by_dist, _ = summary = deconv_data[-1]
        N_good_matched = sum(N_matched_by_dist[x] for x in range(N_errors_allowed+1))
        percent_good_matched = N_good_matched/N_matched*100
        final_result = N_good_matched*N_good_factor + percent_good_matched*percent_good_factor
        if max_N_high is None or parameters[0]<=max_N_high:
            parameters_to_result[parameters] = final_result
    best_parameters = max(parameters_to_result.items(), key=lambda (p,r): r)[0]
    return best_parameters
        

def merge_deconvolution_levels(best_data, good_data, decent_data, include_non_unique=False):
    """ Merge best/good/decent deconvolution data: take best set w/ 0 errors, remaining good set w/ 0-1 errors, remaining decent set.

    All three args should be (insertion_readcounts, insertion_codewords, insertion_samples, insertion_codeword_distances, *extras) 
     tuples: the first four elements should be dictionaries with insertion_position as keys (and values being readcount list, 
      observed codeword, best-matched sample, and #errors between observed codeword and that of the sample) - outputs from 
      combinatorial_deconvolution with return_readcounts=True, or the equivalent; additional elements will be ignored.

    Output will be a dictionary with insertion positions as keys - the full joint set of keys from the three inputs, 
     but excluding ones with no unique match - and 5-tuples with the following elements as values:
        1) readcount list
        2) observed codeword
        3) name of sample with the best-matched expected codeword (or None if no unique match)
        4) number of differences between observed codeword and that of the best-match sample(s)
        5) group name:  "best" if that position had a 0-error match in best_data,
                        else "good" if it had a 0-1 error match in good_data,
                        else "decent" if it had a 0-2 error match in decent_data, 
                        else "poor" if it has a unique match with 3+ errors in decent_data.
    """
    all_data = {}
    def _use_insertion_data(insertion_pos, source, groupname):
        all_data[insertion_pos] = (source[0][insertion_pos], source[1][insertion_pos], source[2][insertion_pos], 
                                   source[3][insertion_pos], groupname)
    all_insertions = set(best_data[0].keys()) | set(good_data[0].keys()) | set(decent_data[0].keys())
    for insertion_pos in all_insertions:
        if best_data[3][insertion_pos] == 0:            _use_insertion_data(insertion_pos, best_data, "best")
        elif good_data[3][insertion_pos] <=1:           _use_insertion_data(insertion_pos, good_data, "good")
        elif decent_data[3][insertion_pos] <=2:         _use_insertion_data(insertion_pos, decent_data, "decent")
        elif decent_data[2][insertion_pos] is not None: _use_insertion_data(insertion_pos, decent_data, "poor")
        elif include_non_unique:                        _use_insertion_data(insertion_pos, decent_data, "not_found")
    return all_data
    # TODO unit-test!  I had an error in here!


def get_genes_and_annotation(insertion_position_set, genefile, annotation_file=None, standard_Phytozome_file_version=None):
    """ Make mutant dataset that gives gene/feature/annotation for each insertion position (all readcounts are 0). 

    Make mutant_analysis_classes.Insertional_mutant_pool_dataset containing zero-readcount mutants for each position, 
     then look up gene/feature for the position and annotation for the gene based on genefile and annotation_file.
    """
    # Note: this is separate from get_flanking_seqs because it should be still possible to get the gene annotation when we 
    #           DON'T have flanking seq data!
    dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset()
    # make dataset containing zero-readcount mutants for each position
    #  (using get_mutant rather than add_mutant because it's faster)
    for insertion_pos in insertion_position_set:
        dataset.get_mutant(insertion_pos)
    dataset.find_genes_for_mutants(genefile, detailed_features=True)
    if annotation_file is not None:
        dataset.add_gene_annotation(annotation_file, standard_Phytozome_file_version)
    return dataset


def full_plate_well_gene_results(dict_full_results_plate, dict_full_results_well, outfilename='', genefile=None, 
                                 annotation_file=None, standard_Phytozome_file_version=None, flanking_sequence_data=None):
    """ Generate single table with plate/well deconvolution results, and including gene/feature/annotation if available.

    Dict_full_results_* should be a name:results dict, where name is a description like 5'/3' and results is the output of 
     merge_deconvolution_levels; the key set must be the same for both.  
    Genefile is the path/name of the gff file giving gene locations; annotation_file is the path/name of the tab-separated
     gene annotation file; standard_Phytozome_file_version should be 4 or 5 if it's in normal Phytozome annotation format 
     for v4.3 or v5 chlamy genome, None otherwise.
    Flanking_seq_data should be a mutant dataset that will be used to get the most common observed flanking sequence for each 
     insertion, or None if flanking sequences are not needed in the file.

    Return (header, full_data_table) tuple: full_data_table is a list of mapped insertion position data lists, with the information 
     for each insertion position arranged in fields described by the header; full_data_table is sorted by plate/well.

    If outfilename is given, also print the data to a tab-separated file with that name. 
    """
    if not dict_full_results_plate.keys() == dict_full_results_well.keys(): 
        raise DeconvolutionError("Plate and well data have different subsets! %s, %s"%(dict_full_results_plate.keys(),
                                                                                       dict_full_results_well.keys()))
    # convenience function to get plate/well data if it's present, otherwise defaults
    def _get_deconv_data_for_ins(insertion_pos, deconv_dict):
        try:                ins_data = deconv_dict[insertion_pos]
        except KeyError:    ins_data = None
        if ins_data is not None:
            readcount_list, _, sample_name, N_errors, category = ins_data
            return (sample_name, category, N_errors, numpy.average(readcount_list))
        else:
            return ('-', 'unmapped', '-', None)

    # grab all insertion position and side combinations - we can't just go by insertion positions, because there may be
    #  5' and 3' results with the same position!  (okay, actually not entirely, but who knows)
    insertions_and_sides = []
    for side in dict_full_results_plate.keys():
        for insertion_pos in set(dict_full_results_plate[side].keys() + dict_full_results_well[side].keys()):
            insertions_and_sides.append((insertion_pos, side))
    # now grab all the deconvolution information
    insertion_basics = []
    for insertion_pos, side in insertions_and_sides:
        plate,plate_category,plate_errors,plate_avg_readcount =_get_deconv_data_for_ins(insertion_pos, dict_full_results_plate[side])
        well, well_category, well_errors, well_avg_readcount = _get_deconv_data_for_ins(insertion_pos, dict_full_results_well[side])
        insertion_basics.append((plate, well, side, insertion_pos, plate_category, well_category, plate_errors, well_errors, 
                                 plate_avg_readcount, well_avg_readcount))
    # grab gene/feature/annotation info through making a proper mutant-position dataset 
    #  (merging 5'/3' - we already have the correct insertion positions as +/-, so we can treat them the same now)
    all_insertion_positions = set.union(*[set(d.keys()) for d in dict_full_results_plate.values() + dict_full_results_well.values()])
    if genefile:
        dataset_with_gene_info = get_genes_and_annotation(all_insertion_positions, genefile, 
                                                          annotation_file, standard_Phytozome_file_version)
    # make header
    # MAYBE-TODO is that the best field order?
    header = 'plate well plate_category well_category plate_errors well_errors side'.split()
    header += 'chromosome strand min_position full_position'.split()
    if genefile:
        header += 'gene orientation feature'.split()
    header += 'plate_avg_readcount well_avg_readcount'.split()
    if flanking_sequence_data:
        header += 'flanking_seq'.split()
    if annotation_file:
        header += dataset_with_gene_info.gene_annotation_header
        missing_gene_annotation_data = tuple(['NO GENE DATA'] + ['' for x in dataset_with_gene_info.gene_annotation_header[:-1]])
    # make full data table
    full_data_table = []
    _sort_dash_last = lambda x: 'zzz' if x=='-' else x
    _sort_function = lambda (p,w,s,i,c,d,e,f,a,b): ((_sort_dash_last(p), _sort_dash_last(w), 1 if s=="5'" else 2, i))
    for insertion_data in sorted(insertion_basics, key = _sort_function):
        side = insertion_data[2]
        line_fields = insertion_data[:2] + insertion_data[4:8] + (side,)
        pos = insertion_data[3]
        line_fields += (pos.chromosome, pos.strand, pos.min_position, pos.full_position)
        if genefile:
            mutant = dataset_with_gene_info.get_mutant(pos)
            line_fields += (mutant.gene, mutant.orientation, mutant.gene_feature)
        line_fields += insertion_data[8:10]
        if flanking_sequence_data:
            line_fields += (flanking_sequence_data[side].get_mutant(pos).get_main_sequence()[0], )
        if annotation_file:
            if mutant.gene_annotation:  line_fields += tuple(mutant.gene_annotation)
            else:                       line_fields += missing_gene_annotation_data
        full_data_table.append(line_fields)
    # print to file if desired
    if outfilename:
        with open(outfilename, 'w') as OUTFILE:
            OUTFILE.write('# ' + '\t'.join(header) + '\n')
            for insertion_data in full_data_table:
                OUTFILE.write('\t'.join([str(x) for x in insertion_data]) + '\n')
    return header, full_data_table
    # TODO unit-test?


################################################## Plotting the data #######################################################

def absent_present_readcounts(readcounts, codeword):
    """ Given a list of readcounts and a binary codeword, return separate lists of readcounts with 0s and with 1s in the codeword.
    """
    readcounts_and_presence = zip(readcounts, codeword.list())
    absent_readcounts = [readcount for (readcount,presence) in readcounts_and_presence if not presence]
    present_readcounts = [readcount for (readcount,presence) in readcounts_and_presence if presence]
    return absent_readcounts, present_readcounts


def plot_mutants_and_cutoffs(insertion_readcount_table, insertion_codewords, order='by_average_readcount', min_readcount=20, 
                             N_rows=60, N_columns=13, N_insertions=None, x_scale='symlog', markersize=4, marker='d', 
                             filename=None, title='insertion readcounts (white - absent, black - present)'):
    """ ____

    If filename is not given, plots will be interactive; otherwise they won't, and will be saved as multiple files.
    """
    # TODO write docstring!
    # MAYBE-TODO add option to only plot insertions present in a specific dataset/list?
    if filename is not None:
        mplt.ioff()
    # TODO add give option for insertion_codewords to be None - all readcounts treated the same, no present/absent distinction
    # make a joint [(readcounts, codeword)] list, discarding insertion position data
    # originally the insertion order is by position - reorder it if desired
    # MAYBE-TODO add different order options?
    readcounts_and_codewords = [(readcounts, insertion_codewords[insertion_pos])
                                 for (insertion_pos,readcounts) in sorted(insertion_readcount_table.items())]
    if order=='by_position':            pass
    elif order=='random':               random.shuffle(readcounts_and_codewords)
    elif order=='by_average_readcount':  readcounts_and_codewords.sort(key=lambda (r,c): numpy.mean(r))
    elif order=='by_median_readcount':  readcounts_and_codewords.sort(key=lambda (r,c): numpy.median(r))
    else:                               raise DeconvolutionError("Unknown order value %s!"%order)
    # if desired, remove any data points that have no readcounts above min_readcount AND have empty codewords
    if min_readcount is not None:       
        readcounts_and_codewords = filter(lambda (r,c): max(r)>=min_readcount or c.weight(), 
                                          readcounts_and_codewords)
    # if desired, take some number of insertions, from the start
    if N_insertions is not None:        readcounts_and_codewords = readcounts_and_codewords[:N_insertions]
    N_per_page = N_rows*N_columns
    N_pages = int(math.ceil(len(readcounts_and_codewords)/N_per_page))
    max_readcount = max(max(readcounts) for (readcounts,_) in readcounts_and_codewords)
    dot_plot_kwargs = dict(linestyle='None', markeredgecolor='k', markersize=markersize)
    for page_N in range(N_pages):
        mplt.figure(figsize=(16,12))
        mplt.suptitle(title + '; page %s/%s'%(page_N+1, N_pages), y=0.92)      # y position to avoid wasting space
        page_first_ins_N = page_N*N_per_page
        for column_N in range(N_columns):
            col_first_ins_N = page_first_ins_N + column_N*N_rows
            if col_first_ins_N >= len(readcounts_and_codewords):
                break
            mplt.subplot(1, N_columns, column_N+1)
            # plot a grey line for each insertion, with dots showing readcounts (filled if "present", unfilled if "absent")
            # TODO currently we don't explicitly know the cutoffs, just the codewords!  Should I calculate (or save) the real cutoffs used, instead?  Or infer them based on the highest absent and lowest present readcount and plot those?
            mplt.hlines(range(min(N_rows, len(readcounts_and_codewords)-col_first_ins_N)), 0, max_readcount, colors='0.7')
            # accumulate all dot positions before plotting them all together for the whole column
            absent_dot_positions_x, absent_dot_positions_y = [], []
            present_dot_positions_x, present_dot_positions_y = [], []
            for row_N in range(N_rows):
                curr_ins_N = col_first_ins_N+row_N
                if curr_ins_N >= len(readcounts_and_codewords):
                    break
                absent_readcounts, present_readcounts = absent_present_readcounts(*readcounts_and_codewords[curr_ins_N])
                absent_dot_positions_x.extend(absent_readcounts)
                absent_dot_positions_y.extend([row_N for _ in absent_readcounts])
                present_dot_positions_x.extend(present_readcounts)
                present_dot_positions_y.extend([row_N for _ in present_readcounts])
            mplt.plot(absent_dot_positions_x, absent_dot_positions_y, marker, markerfacecolor='w', **dot_plot_kwargs)
            mplt.plot(present_dot_positions_x, present_dot_positions_y, marker, markerfacecolor='k', **dot_plot_kwargs)
            mplt.xscale(x_scale)
            # remove the axes/etc we don't want, make things pretty
            mplt.xlim(-1, max_readcount*2)
            mplt.ylim(-1, N_rows)
            mplt.yticks([], [])
            mplt.xticks(mplt.xticks()[0], [])
            ax = mplt.gca()
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.yaxis.set_ticks_position('none')
            mplt.draw()
            # TODO tighten layout
            mplt.subplots_adjust(wspace=0.1)
            #mplt.tight_layout() - this isn't great, messes up suptitle
        if filename is not None:
            basename, ext = os.path.splitext(filename)
            pagenum_len = len(str(N_pages+1))
            format_str = '_%%0%sd'%pagenum_len
            plotting_utilities.savefig(basename + format_str%(page_N+1) + ext)
            mplt.close()
    if filename is not None:
        mplt.ion()


def plot_parameter_space_1(deconv_data_by_parameters, N_errors_allowed, chosen_parameters, info='', 
                           marker='x', alpha=1, colors_by_N_high='orangered black blue'.split()):
    """ Given deconvolution data for various parameter sets, plot # and % of matched with <= N errors, and highlight chosen set.

    The deconv_data_by_parameters arg should be a dict with (N_high, N_low, cutoff_pos, overall_min) tuples as keys and 
     (insertion_readcounts, insertion_codewords, insertion_samples, codeword_distances, summary) tuples as values.
        (The first four should be outputs from combinatorial_deconvolution with return_readcounts=True, or the equivalent;
         summary should be the output from get_deconvolution_summary.)
        
    Each item of deconv_data_by_parameters will be plotted as a dot. 
    The chosen_parameters should be one of the deconv_data_by_parameters keys - that dot will be highlighted.

    N_errors_allowed should be an integer 0 or higher, and will be used to define what the axes are:
        - x axis will the the number of insertion positions that were mapped to expected codewords with N or fewer errors
        - y axis will be that number divided by the total number of uniquely mapped insertion positions (with any #errors)

    Info will be used in the plot title.

    The marker, alpha, and colors_by_N_high args set the physical appearance of the markers.

    """
    if not chosen_parameters in deconv_data_by_parameters.keys():
        raise DeconvolutionError("Chosen parameters not present in result set! %s, %s"%(chosen_parameters, 
                                                                                        deconv_data_by_parameters.keys()))
    N_good_matched_val_dict, percent_good_matched_val_dict = defaultdict(list), defaultdict(list)
    for parameters, deconv_data in deconv_data_by_parameters.items():
        N_matched, N_unmatched, N_matched_by_dist, _ = summary = deconv_data[-1]
        N_good_matched = sum(N_matched_by_dist[x] for x in range(N_errors_allowed+1))
        percent_good_matched = N_good_matched/N_matched*100
        N_good_matched_val_dict[parameters[0]].append(N_good_matched)
        percent_good_matched_val_dict[parameters[0]].append(percent_good_matched)
        if parameters==chosen_parameters:
            chosen_N_good_matched, chosen_percent_good_matched = N_good_matched, percent_good_matched
    min_N_high = min(N_good_matched_val_dict.keys())
    for N_high,N_good_matched_vals in N_good_matched_val_dict.items():
        mplt.plot(N_good_matched_vals, percent_good_matched_val_dict[N_high], 
                  marker=marker, color=colors_by_N_high[N_high-min_N_high], markersize=5, alpha=alpha, linestyle='None')
    mplt.plot(chosen_N_good_matched, chosen_percent_good_matched,  
              marker='o', markerfacecolor='None', markeredgecolor='r', markersize=10, markeredgewidth=2, linestyle='None')
    mplt.title("Deconvolution results for different parameters%s.\n Allowing %s errors; red circle marks chosen parameters."%(
        ': %s'%info if info else '', N_errors_allowed))
    mplt.xlabel("number of insertions mapped with at most %s errors"%N_errors_allowed)
    mplt.ylabel("% of those out of all uniquely mapped insertions")
    # MAYBE-TODO could try changing the marker shape/size/color to reflect the other parameters?


def codeword_position_errors(observed_and_expected_codewords, name='?'):
    """ Return and plot how often which codeword position is different between observed and closest matchd expected codewords. 
    
    Input should be a list of (observed, closest_matched_expected) codeword tuples.
    """
    N_errors_per_position = [0 for x in str(observed_and_expected_codewords[0][0])]
    for observed,expected in observed_and_expected_codewords:
        for pos,(obs,exp) in enumerate(zip(str(observed), str(expected))):
            if not obs==exp:
                N_errors_per_position[pos] += 1
    return N_errors_per_position


def plot_plate_well_positions(well_list, marker='o', markersize=5, markeredgecolor='None', markerfacecolor='black', offset=0):
    """ Plot the given well_list on a 384-well plate. """
    columns = [int(well[1:]) for well in well_list]
    rows = [ord(well[0]) for well in well_list]
    mplt.plot([x+offset for x in columns], rows, marker, 
              markersize=markersize, markeredgecolor=markeredgecolor, markerfacecolor=markerfacecolor)

################################################## Extra analysis/utilities #######################################################

def grab_matching_5prime_3prime(DECONV_data, max_distance, position_index=6):
    """ Grab all the cases with perfectly or near-perfectly matching 5' and 3' ends.

    Skip "colonies" where the plate or well is unknown (i.e. '-')! (MAYBE-TODO
    look at those too? If the distance is 0, they're almost certainly really
    the 5' and 3' end of the same thing...)

    When there are multiple 3' or 5' cases, look at all combinations, and take
    the best one that's within 1kb (must be same chromosome and strand)

    Return a (plate,well):(distance_between_sides, position_5prime, position_3prime) dict, 
     where positions are (chrom, strand, min_pos, full_pos) tuples.
    """
    positions_by_colony = defaultdict(lambda: defaultdict(list))
    for line in DECONV_data:
        plate, well = line[:2]
        side, chrom, strand, min_pos, full_pos = line[position_index:position_index+5]    
        if '-' not in (plate, well):
            positions_by_colony[(plate, well)][side].append((chrom, strand, min_pos, full_pos))
    # When there are multiple 3' or 5' cases, look at all combinations, and take the best one that's within 1kb (must be same chromosome and strand)
    distance_and_positions_by_colony = {}
    for colony, data in positions_by_colony.iteritems():
        distances_and_positions = []
        for pos_5prime in data["5'"]:
          for pos_3prime in data["3'"]:
            if pos_5prime[:2] == pos_3prime[:2]:
                dist = abs(int(pos_5prime[2]) - int(pos_3prime[2]))
                if dist < max_distance:  
                    distances_and_positions.append((dist, pos_5prime, pos_3prime))
        if distances_and_positions:
            distance_and_positions_by_colony[colony] = sorted(distances_and_positions)[0]
    return distance_and_positions_by_colony

def get_matching_sides_from_table(DECONV_data_table, ID_cols=[0,1], side_chrom_strand_minpos_cols=[2,3,4,5], 
                                  LEAPseq_cols=[10,12,13], min_conf_dist=500, max_percent_wrong=30, 
                                  distance_cutoffs=[0,10,100,1000,10000], warn_if_3pos_within=100000, only_2ins_cases=False,
                                  print_data=True):
    """ ___
    """
    # TODO write docstring based on comments
    # Grab all the insertions per colony - for real colonies only, no unknown-plate or unknown-well cases, no cassette positions (since those are obviously not real insertion positions anyway, and the cassette is small enough that we could get approximately matching ), no unaligned positions (since they can't match anything anyway).
    # Output: positions_by_colony is a (plate,well):data dict, where data is a (side,chr,strand,minpos):extra_data_list dict, where extra_data_list matches extra_info_header.
    positions_by_colony = collections.defaultdict(set)
    for line in DECONV_data_table:
        ID = tuple(line[x] for x in ID_cols)
        chrom = line[side_chrom_strand_minpos_cols[1]]
        if any(x in ID for x in ('-', None, '?', '')):              continue
        if mutant_analysis_classes.is_cassette_chromosome(chrom):   continue
        if mutant_analysis_classes.is_other_chromosome(chrom):      continue
        side_chrom_strand_minpos = tuple(line[x] for x in side_chrom_strand_minpos_cols)
        positions_by_colony[ID].add(side_chrom_strand_minpos)
    # When there are multiple insertions per colony, sort them by position (in order: chrom, minpos, strand, side), 
    #  and take each adjacent pair - seems like the most sensible way.  Could take all pairs instead, but that would be too many.
    insertion_pairs = []
    if only_2ins_cases:
        for colony, data in positions_by_colony.iteritems():
            if len(data) == 2:
                pos1, pos2 = sorted(data, key = lambda (side,chrom,strand,pos):(chrom,pos,strand,side))
                insertion_pairs.append((pos1, pos2, colony))
    else:
        for colony, data in positions_by_colony.iteritems():
            sorted_data = sorted(data, key = lambda (side,chrom,strand,pos):(chrom,pos,strand,side))
            for (pos1, pos2) in zip(sorted_data, sorted_data[1:]):
                insertion_pairs.append((pos1, pos2, colony))
            if warn_if_3pos_within is not None:
                for (pos1, pos2, pos3) in zip(sorted_data, sorted_data[1:], sorted_data[2:]):
                    if pos1[1]==pos3[1] and pos3[3]-pos1[3] <= warn_if_3pos_within:
                        print "Warning: colony %s has 3 positions within %sbp!\n %s, %s, %s"%(colony, warn_if_3pos_within, 
                                                                                              pos1, pos2, pos3)
    # Now figure out the relative direction for each position pair - depends on strand and side!  And how many options are there? outer-cassette, inner-cassette, same-strand, different-chromosome.  I have a function to do this.
    # Categories:
    # - same-direction: both flanking seqs were read in the same direction
    # - inner-cassette: the flanking seqs face "away" from each other, as if there was an insertion between them, possibly with a deletion
    # - outer-cassette: the flanking seqs face "toward" each other, so that if they're within 20bp of each other, they actually overlap, as if there was a duplication with the insertion or we have two cassettes with a junk fragment in between.
    # TODO change those terms to something more intuitive?
    pair_categories_and_distances = []
    for pos1, pos2, colony in insertion_pairs:
        side1, chrom1, strand1, minpos1 = pos1
        side2, chrom2, strand2, minpos2 = pos2
        category, dist = mutant_utilities.insertion_pair_type_distance((chrom1, strand1, minpos1), side1, 
                                                                       (chrom2, strand2, minpos2), side2)
        pair_categories_and_distances.append((dist, category, side1, side2, pos1, pos2, colony))
    # Bin results into distance categories
    data_by_ID_pos = {tuple([line[x] for x in ID_cols+side_chrom_strand_minpos_cols]): line for line in DECONV_data_table}
    def _LEAPseq_data(line):
        conf_dist, N_conf, N_wrong = [line[x] for x in LEAPseq_cols]
        percent_wrong = N_wrong/(N_conf+N_wrong)*100 if N_conf+N_wrong else float('nan')
        return conf_dist, percent_wrong
    def _is_high_conf(pos, colony):
        line = data_by_ID_pos[tuple(list(colony) + list(pos))]
        conf_dist, percent_wrong = _LEAPseq_data(line)
        return (conf_dist >= min_conf_dist and percent_wrong <= max_percent_wrong)

    category_dist_count_data = {}
    for curr_category in 'inner-cassette outer-cassette same-direction diff-chrom'.split():
        dist_counts = [[0, 0, 0] for _ in range(len(distance_cutoffs)+1)]
        category_dist_count_data[curr_category] = dist_counts
        for (dist, category, side1, side2, pos1, pos2, colony) in pair_categories_and_distances:
            if category != curr_category:   continue
            cat_N = len(distance_cutoffs)
            for N,dist_cutoff in enumerate(distance_cutoffs):
                if dist <= dist_cutoff: 
                    cat_N = N
                    break
            dist_counts[cat_N][0] += 1
            if side1!=side2:    dist_counts[cat_N][1] += 1
            dist_counts[cat_N][2] += (_is_high_conf(pos1, colony) and _is_high_conf(pos2, colony))

        if print_data:
            # MAYBE-TODO make the bin names nice (use kb etc)
            bin_names = ["0bp"] + ["%s-%sbp"%(start+1, end) for (start,end) in zip(distance_cutoffs, distance_cutoffs[1:])]
            bin_names.append("%s+bp"%distance_cutoffs[-1])
            info = "%s (%% 5'+3', %% both-high-conf)::: "%curr_category
            if curr_category == 'diff-chrom':
                counts = [sum(x) for x in zip(*dist_counts)]
                info += "%s all (%.0f%%, %.0f%%), "%(counts[0], 0 if not counts[0] else counts[1]/counts[0]*100, 
                                                     0 if not counts[0] else counts[2]/counts[0]*100)
            else:
                for size_range, counts in zip(bin_names, dist_counts):
                    info += "%s %s (%.2g%%, %.2g%%), "%(counts[0], size_range, 0 if not counts[0] else counts[1]/counts[0]*100, 
                                                        0 if not counts[0] else counts[2]/counts[0]*100)
            print info[:-2]

    return positions_by_colony, insertion_pairs, pair_categories_and_distances, category_dist_count_data




def add_RISCC_to_deconvolution_data(DECONV_header, DECONV_data, RISCC_5prime, RISCC_3prime, 
                                      max_allowed_distance, min_weird_distance):
    """ Add RISCC data (two extra fields) to deconvolution output mutant list.
    """
    # TODO more detail in docstring!
    NEW_deconv_data = []
    NEW_deconv_header = DECONV_header[:2] + 'RISCC_status RISCC_distance'.split() + DECONV_header[2:]
    status_counts = defaultdict(int)
    for fields in DECONV_data:
        position = Insertion_position(fields[7], fields[8], full_position=fields[10], immutable=True)
        side = fields[6]
        if side=="5'":      mutant = RISCC_5prime.get_mutant(position)
        else:               mutant = RISCC_3prime.get_mutant(position)
        status = mutant.RISCC_confirmation_status(side, max_allowed_distance, min_weird_distance=min_weird_distance)
        status_distance = mutant.RISCC_max_confirmed_distance(max_allowed_distance)
        if isnan(status_distance): status_distance = '-'
        else:                            status_distance = str(status_distance)
        status_counts[status] += 1
        NEW_deconv_data.append(fields[:2] + tuple([status, status_distance]) + fields[2:])
    for status,count in sorted(status_counts.items()):
       print "%s: %s"%(status, value_and_percentages(count, [len(NEW_deconv_data)]))
    return NEW_deconv_header, NEW_deconv_data, status_counts
     # TODO unit-test!

def get_all_distances(deconv_data):
    """ Return lists of max confirmed distances for all mutants and separated into various categories.
    """
    # TODO more detail in docstring!
    distances_all, distances_genomic, distances_cassette = [], [], []
    distances_by_side_genomic, distances_by_quality_genomic = defaultdict(list), defaultdict(list)
    distances_by_side_cassette, distances_by_quality_cassette = defaultdict(list), defaultdict(list)
    for line in deconv_data:
        # ignore cases with no max confirmed distance
        if line[3] == '-':    continue
        dist = int(line[3])
        distances_all.append(dist)
        if is_cassette_chromosome(line[9]):
            distances_cassette.append(dist)
            distances_by_side, distances_by_quality = distances_by_side_cassette, distances_by_quality_cassette
        else:
            distances_genomic.append(dist)
            distances_by_side, distances_by_quality = distances_by_side_genomic, distances_by_quality_genomic
        distances_by_side[line[8]].append(dist)
        categories = line[4:6]
        if 'unmapped' in categories:    distances_by_quality['unmapped'].append(dist)
        elif 'poor' in categories:      distances_by_quality['poor'].append(dist)
        elif 'decent' in categories:    distances_by_quality['decent'].append(dist)
        elif 'good' in categories:      distances_by_quality['good'].append(dist)
        elif 'best' in categories:      distances_by_quality['best'].append(dist)
        else:                           print "weird categories??? %s"%' '.join(line)
    return distances_all, distances_genomic, distances_cassette, distances_by_side_genomic, distances_by_quality_genomic, distances_by_side_cassette, distances_by_quality_cassette
    # TODO unit-test!

def seqs_per_colony(DECONV_data, side_index, quiet=False):
    """ Return the number of colonies with 1, 2 etc sequences, total and separated by 5'/3'.  Also print.
    
    DECONV_data should be a list of lines, with line[0] being the mutant ID and line[side_index] being 5' or 3'.
    Returns N_insertions:(N_colonies, side_details) dict, where side_details is a (N_5prime, N_3prime):N_colonies dict.
    """
    colony_Nins = collections.defaultdict(lambda: {"5'": 0, "3'": 0})
    for line in DECONV_data:
        colony = line[0]
        side = line[side_index]
        colony_Nins[colony][side] += 1
    Nins_counts = collections.defaultdict(lambda: [0, collections.defaultdict(int)])
    for Nins_data in colony_Nins.values():
        N_5prime, N_3prime = Nins_data["5'"], Nins_data["3'"]
        total = N_5prime+N_3prime
        Nins_counts[total][0] += 1
        Nins_counts[total][1][(N_5prime, N_3prime)] += 1
    if not quiet:
        for Nins, (Nins_count, Nins_detail) in sorted(Nins_counts.items()):
            print "%s insertions - %s total. By number of 5' and 3': %s."%(Nins, Nins_count, 
                           ', '.join("%s+%s - %s"%(N5, N3, count) for ((N5, N3), count) in sorted(Nins_detail.items())))
    return Nins_counts

def plot_seqs_per_colony(Nins_counts, collapse_from=4, separate_3p_5p=True):
    """ Plot histogram of how many colonies have how many flanking sequences, based on output of seqs_per_colony.

    Flanking sequence numbers >= collapse_from will be collapsed together.
    """
    import copy
    Nins_counts_collapsed = copy.deepcopy(Nins_counts)
    max_Nins_count = max(Nins_counts.keys())
    if collapse_from < max_Nins_count:
        for x in range(collapse_from+1, max_Nins_count+1):
            del Nins_counts_collapsed[x]
            Nins_counts_collapsed[collapse_from][0] += Nins_counts[x][0]
            Nins_counts_collapsed[collapse_from][1].update(Nins_counts[x][1])
    data_range = sorted(Nins_counts_collapsed.keys())
    data_names = data_range if collapse_from >= max_Nins_count else data_range[:-1] + ['%s-%s'%(collapse_from, max_Nins_count)]
    if not separate_3p_5p:
        mplt.bar(sorted(Nins_counts_collapsed.keys()), 
                 [Nins_counts_collapsed[x][0] for x in sorted(Nins_counts_collapsed.keys())], align='center')
        mplt.xticks(data_range, data_names)
        mplt.xlim(0.4, collapse_from + .6)
    else:
        # make new details dict that's just separated into "all 5'", "all 3'" and "some of each" instead of by number of 5'/3'
        Nins_counts_details = collections.defaultdict(lambda: {"5'": 0, "3'": 0, "both": 0})
        for (N, (_, details)) in Nins_counts_collapsed.items():
            for ((N5, N3), count) in details.items():
                if N5==0:   Nins_counts_details[N]["3'"] += count
                elif N3==0: Nins_counts_details[N]["5'"] += count
                else:       Nins_counts_details[N]["both"] += count
        sides = "3' both 5'".split()
        data = [[Nins_counts_details[x][side] for side in sides] for x in data_range]
        bar_list = plotting_utilities.stacked_bar_plot(data, sample_names='1 2 3 4-6'.split(), colors='bcg')
        mplt.legend(bar_list, sides)
    xticks = list(range(1,collapse_from+1))
    mplt.ylabel('number of colonies')
    mplt.xlabel('number of insertion junctions')


###################################################### Testing ###########################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__get_original_384_well_numbers(self):
        # assuming the 384-well plate was split into four 96-well plates, every other row/column - SEE DOCSTRING FOR INFO
        #   all the examples here are straight from the get_original_384_well_numbers docstring.
        assert get_original_384_well_numbers(1)   == 'A01'
        assert get_original_384_well_numbers(2)   == 'A03'
        assert get_original_384_well_numbers(3)   == 'A05'
        assert get_original_384_well_numbers(13)  == 'C01'
        assert get_original_384_well_numbers(14)  == 'C03'
        assert get_original_384_well_numbers(95)  == 'O21'
        assert get_original_384_well_numbers(96)  == 'O23'
        assert get_original_384_well_numbers(97)  == 'A02'
        assert get_original_384_well_numbers(98)  == 'A04'
        assert get_original_384_well_numbers(383) == 'P22'
        assert get_original_384_well_numbers(384) == 'P24'

    def _make_positions(self, N=3):
        """ Help function - make three insertion positions, on chromosomes 1, 2, 3, varying strands/etc, immutable. """
        positions = []
        _I = mutant_analysis_classes.Insertion_position
        positions.append(_I('chr1', '+', position_before=100, immutable=True))
        positions.append(_I('chr2', '-', position_after=501, immutable=True))
        positions.append(_I('chr3', '+', position_before=200, position_after=201, immutable=True))
        positions.append(_I('chr4', '-', position_after=20, immutable=True))
        positions.append(_I('chr5', '_', position_before=444, immutable=True))
        return positions[:N]

    def test__datasets_to_readcount_table(self):
        # make a test case with 3 pools and 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        pools = ['A', 'B', 'C']
        pos1, pos2, pos3 = self._make_positions()
        insertion1 = mutant_analysis_classes.Insertional_mutant_multi_dataset(pos1)
        insertion2 = mutant_analysis_classes.Insertional_mutant_multi_dataset(pos2)
        insertion3 = mutant_analysis_classes.Insertional_mutant_multi_dataset(pos3)
        insertions = [insertion1, insertion2, insertion3]
        # the three numerical arguments to add_counts are total_reads, perfect_reads, sequence_variants (not used).
        insertion1.add_counts(0,  0,  0, dataset_name=pools[0])
        insertion1.add_counts(1,  0,  1, dataset_name=pools[1])
        insertion1.add_counts(10, 9,  1, dataset_name=pools[2])
        insertion2.add_counts(9,  7,  1, dataset_name=pools[0])
        insertion2.add_counts(20, 20, 1, dataset_name=pools[1])
        insertion2.add_counts(0,  0,  0, dataset_name=pools[2])
        insertion3.add_counts(1,  0,  1, dataset_name=pools[0])
        insertion3.add_counts(0,  0,  0, dataset_name=pools[1])
        insertion3.add_counts(3,  3,  1, dataset_name=pools[2])
        dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset=True)
        for insertion in insertions:  dataset.add_mutant(insertion)
        assert datasets_to_readcount_table(dataset, pools, use_perfect_reads=False) == {pos1:[0,1,10], pos2:[9,20,0], pos3:[1,0,3]}
        assert datasets_to_readcount_table(dataset, pools, use_perfect_reads=True) ==  {pos1:[0,0,9],  pos2:[7,20,0], pos3:[0,0,3]}

    def _make_basic_readcount_table(self):
        # make a test case with 3 pools and 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        pos1, pos2, pos3 = self._make_positions()
        return {pos1:[0,1,10], pos2:[9,20,0], pos3:[1,0,3]}

    def _simplify_codeword_output(self, insertion_codeword_dict, ins_positions):
        """ Convenience function for easier processing of codeword tables - return them as a string in given order. """
        return ' '.join([str(insertion_codeword_dict[ins_position]) for ins_position in ins_positions])

    # Note: all the readcounts_to_presence__* tests below ALSO test readcounts_to_codewords, 
    #  since that's just a convenience function to do readcounts_to_presence on a larger scale and change outputs to codewords.

    def test__readcounts_to_presence__cutoffs(self):
        # 3 insertions: readcounts 0,1,10; 9,20,0; 1,0,3
        insertions = self._make_positions()
        readcounts = self._make_basic_readcount_table()
        _S = self._simplify_codeword_output
        function_with_cutoff = lambda cutoff: (lambda x: readcounts_to_presence__cutoffs(x, cutoff))
        # single cutoff (checking all relevant values)
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(0)), insertions)  == '111 111 111'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(1)), insertions)  == '011 110 101'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(2)), insertions)  == '001 110 001'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(3)), insertions)  == '001 110 001'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(4)), insertions)  == '001 110 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(9)), insertions)  == '001 110 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(10)), insertions) == '001 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(11)), insertions) == '000 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(20)), insertions) == '000 010 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(21)), insertions) == '000 000 000'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff(99)), insertions) == '000 000 000'
        # per-pool cutoffs
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([1,1, 1])), insertions) == '011 110 101'
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([3,10,1])), insertions) == '001 110 001' 
        assert _S(readcounts_to_codewords(readcounts, function_with_cutoff([2,10,4])), insertions) == '001 110 000' 

    def test__readcounts_to_presence__mutant_1(self):
        # 3 insertions with various readcount numbers: perfect, noisy, absent, weird
        # MAYBE-TODO try changing the order of the readcounts and making sure the results are the same (but in different order)?
        insertions = self._make_positions(4)
        readcounts = dict(zip(insertions, ([0,0,0,0,10,10,10], [0,1,1,3,8,10,15], [0,0,0,0,1,1,1], [0,1,2,3,4,5,6])))
        _S = self._simplify_codeword_output
        # arguments to readcounts_to_presence__mutant_1 are: readcount_list, N_always_present, N_always_absent, overall_min, 
        #                                               present_level_function, absent_level_function, cutoff_position, min_cutoff
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 1, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000111 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 2, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # upping the overall_min makes more all-0 results
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 6, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000000'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 10, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000000'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 11, numpy.median, numpy.median, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000000 0000000 0000000 0000000'
        # upping min_cutoff makes more single 0s
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 5)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0000011'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.5, 10)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000011 0000000 0000000'
        # changing cutoff_position up/down from 0.5 - more/fewer 0s
        #  (readcounts from above: [0,0,0,0,10,10,10], [0,1,1,3,8,10,15], [0,0,0,0,1,1,1], [0,1,2,3,4,5,6])
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.8, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000011 0000000 0000011'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.median, numpy.median, 0.2, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0001111 0000000 0011111'
        ### changing the level functions from median to max/min:
        # mean of 2 numbers is the same as median; not true for 3 numbers, but in this case it ends up the same
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 2, 2, 2, numpy.mean, numpy.mean, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, numpy.mean, numpy.mean, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # min/max - actually ends up pretty boring...
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, max, min, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        conv_function = lambda x: readcounts_to_presence__mutant_1(x, 3, 3, 2, min, max, 0.5, 1)
        assert _S(readcounts_to_codewords(readcounts, conv_function), insertions)  == '0000111 0000111 0000000 0001111'
        # TODO make more interesting min/max tests!  Or other functions
        # TODO test N_always_absent == 0

    def test__read_codewords_from_file(self):
        # convenience function: compare real output (includes Binary_codeword objects) to simple string representation of dict
        def _compare(output, expected_string):
            assert len(output) == len(expected_string.split(' '))
            for single_string in expected_string.split(' '):
                key,val = single_string.split(':')
                assert str(output[key]) == val
        # basic file
        infile1 = 'test_data/INPUT_codewords_1.txt'
        output1 = read_codewords_from_file(infile1)
        _compare(output1, '0:001 1:010 2:011 3:100 4:101 5:110 6:111')
        # a file with another code, and with a mirror-codeword section too
        infile2 = 'test_data/INPUT_codewords_2.txt'
        output2 = read_codewords_from_file(infile2)
        _compare(output2, '0:0011 1:0101 2:0110 3:1001 4:1010 5:1100 6:1111')
        # trying out the new_sample_names optional dict
        output1b = read_codewords_from_file(infile1, {str(x):str(x+10) for x in range(7)})
        _compare(output1b, '10:001 11:010 12:011 13:100 14:101 15:110 16:111')
        output2b = read_codewords_from_file(infile2, {'0':'dA', '1':'dB', '2':'dC', '3':'dD', '4':'dE', '5':'dF', '6':'dG'})
        _compare(output2b, 'dA:0011 dB:0101 dC:0110 dD:1001 dE:1010 dF:1100 dG:1111')
        # fail if new_sample_names values are non-unique
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):'1' for x in range(7)})
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):min(x,5) for x in range(7)})
        self.assertRaises(DeconvolutionError, read_codewords_from_file, infile1, {str(x):('A' if x<3 else 'B') for x in range(7)})

    def test__find_closest_sample_codeword(self):
        sample_codewords = {x[0]:binary_code_utilities.Binary_codeword(x[2:]) for x in 'A:1100 B:1010 C:0011 D:1000'.split()}
        # diff 1 - if the top and second Hamming distance differ by 1, take the top one
        inputs_outputs_diff1 = ('1100 A 0, 1010 B 0, 0011 C 0, 1000 D 0, '
                                +'1111 None 2, 0000 D 1, 0110 None 2, 0111 C 1, 1110 None 1, 0001 C 1, 0100 A 1')
        # diff 1 - if the top and second Hamming distance differ by 1, count that as None - they must differ by at least 2
        inputs_outputs_diff2 = ('1100 None 0, 1010 None 0, 0011 C 0, 1000 None 0, '
                                +'1111 None 2, 0000 None 1, 0110 None 2, 0111 C 1, 1110 None 1, 0001 None 1, 0100 None 1')
        for diff, inputs_outputs in [(1, inputs_outputs_diff1), (2, inputs_outputs_diff2)]:
            for input_output_str in inputs_outputs.split(', '):
                input_str, sample, distance = input_output_str.split(' ')
                distance = int(distance)
                for input_val in (input_str, binary_code_utilities.Binary_codeword(input_str)):
                    out_sample, out_dist = find_closest_sample_codeword(input_val, sample_codewords, diff)
                    if sample == 'None':    assert (out_sample is None and out_dist == distance)
                    else:                   assert (out_sample, out_dist) == (sample, distance)

    def test__seqs_per_colony(self):
        assert seqs_per_colony([], 1) == {}
        # doing comparisons indirectly to get around all the defaultdicts in the output
        results = seqs_per_colony([[1, "3'"], [1, "5'"], [2, "3'"], [2, "3'"], [3, "5'"]], 1, quiet=True)
        assert len(results) == 2
        assert len(results[1]) == len(results[2]) == 2
        assert results[1][0] == 1
        assert results[2][0] == 2
        assert dict(results[1][1]) == {(1,0):1}
        assert dict(results[2][1]) == {(1,1):1, (0,2):1}
    # LATER-TODO add more unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
