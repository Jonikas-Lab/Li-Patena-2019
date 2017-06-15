#! /usr/bin/env python2.7
"""
Simulation functions for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import time
import os
import sys
import math
import random
import collections
# other packages
import numpy
import scipy
import scipy.stats
# my modules
import general_utilities
import statistics_utilities
import basic_seq_utilities
import mutant_analysis_classes
import mutant_utilities
from mutant_utilities import get_histogram_data_from_positions, get_mutant_positions_from_dataset, get_chromosome_lengths, get_20bp_fraction
import gff_examine_file
from gff_examine_file import MultipleRNAError, NoRNAError, SeqIO

# constants, help functions, etc
format_bp = lambda x: basic_seq_utilities.format_base_distance(x, False)    # the False is to not approximate

######################################## Calculating mappability over genome #######################################

### MEMORY USAGE ISSUES
    # Keeping the full set of mappable/unmappable genome positions, in whatever form, takes a lot of memory!
    # Use sys.getsizeof to get the sizes of basic objects:  
    #  - an int or bool takes 24 bytes; 
    #  - a 21-string takes 58, a 1-string takes 38
    #  - sets and dicts take 1-2x as much memory as the objects themselves.
    #  - lists actually take a lot less!  About 9 bits per int or char (TODO is that really true? Geoffrey says not. Did I test it with big ints and non-identical chars?).
    #  - numpy arrays take a lot less, but I think sys.getsizeof lies about them because they're objects... Still, in practice on real mappability data they seem to take about 2x less than lists, which gets it down to something manageable.  They're also about 2x bigger than the non-numpy version when pickled (1.3G vs 800M), but that's not as big an issue, I'd say.
    # 
    # The genome's about 113M bases, and we already know that most of it's mappable, i.e. unique sequences.  
    #  - that many ints in a set (say 1.5x the cost) will be 4G!  
    #  - storing that many 21-strings in a set (say 1.5x the cost) will be 9G!  
    #  - plus even more for the seq:pos dict, depending on details of pos format...

    ### MAYBE-TODO Improvement options for memory issues:
    # 
    # For the second pass and returned value: Instead of a list of genome positions, keep a genome-sized list filled with 0s and 1s?  With 0/1 as normal ints that would still be about 1G, but I'm sure switching to a more efficient representation (bitstring.BitArray?) would make it a LOT smaller.  Actually it might be smaller as a numpy array too, those may account for element size...
    # 
    # Could I do everything as numpy arrays instead of making lists (with append/extend) and then converting to numpy arrays?  NOT SURE.  I could write the list-creation as a generator, but you can't make numpy arrays from generators without converting to a list first, I checked...
    # 
    # For the first pass - NOT SURE!  How can I more efficiently store sets of strings?  One thing I could do would be use ints or bitarrays or something instead of strings, since all my strings are actually just made of ACTG, with 2 bits per character, not the full range of characters (well, there may be Ns, but I can skip those...). 
    # MAYBE-TODO think about it more - it's currently the part that takes the most memory.  But I only run it once, so that's not really a big issue - it's mostly the returned values that are important.

    ## Issues with python not releasing memory, or something
    # 
    # Technically, on the whole genome, genome_mappable_slices should take 5-8GB, but RELEASE it once done (see "del" calls) - the output is only <1G.  But python seems to keep the memory, even if the output is about the size I expected!  So if I do two runs, say 20bp and 21bp flanks, I end up exceeding my 16G and going into swap... 
    # genome_mappable_insertion_sites* doesn't take as much memory to run, though still some; the output is actually 2x bigger than genome_mappable_slices for genome_mappable_insertion_sites, since there are two mappable insertion positions per slice, and 2Nx bigger for genome_mappable_insertion_sites_multi (4x for two slice lengths, etc - comes out to about 3G). 
    #
    # So it's probably a good idea to pickle the output (after running genome_mappable_slices or genome_mappable_insertion_sites*), close python to release the memory (even though it'll take a LONG TIME, especially if there was swap space used), and open a new instance to unpickle the data and work on it.  
    # BUT if I pickle the output of any of those, and unpickle it in a new interactive python shell, again the output itself doesn't end up too huge (see sizes above), BUT the unpickling seems to take more memory and not release it again! Ugh. So maybe I really should convert it to bitarrays or something.  Although numpy arrays helped some. 

### genome mappability by position

def genome_mappable_slices(slice_len, genome_seq, print_info=True):
    """ Return a chrom:position_array dict giving the positions of all the unique slice_len seq slices. 
    
    Look at all slices of length slice_len in each chromosome (overlapping - for AATGG with len 2, look at AA, AT, TG, GG)
     return the positions of only the ones that showed up exactly once (in forward or reverse-complement - 
      sequences that have a reverse-complement version later in the genome, or that are palindromes, don't show up).
    Positions are the 1-based position of the first base of each slice (i.e. in AATTGGCC, the position of AAT is 1).

    Genome_seq can be a chr_name:chr_seq dict, or the path to a fasta file containing the genome; default file will be used if None. 

    CAVEAT: if you run this on the whole Chlamy genome, around 100Mb, it can take 5-8G of RAM, and returns a 1G data structure!
    """
    # keep original genome_seq, because we'll want it for the second pass, and generators only work once
    original_genome_seq = genome_seq
    if isinstance(genome_seq, str):  genome_seq = basic_seq_utilities.parse_fasta(genome_seq)
    else:                            genome_seq = genome_seq.iteritems()
    ### This is done in two passes over the whole genome, to improve memory usage. The original version was done in one pass.
    ### First pass - go over all chromosome and get a set of unique and non-unique sequences. 
    sequences_seen_once = set()
    sequences_seen_twice = set()
    chrom_names = set()
    N_total_slices = 0
    if print_info:  print "PASS 1 - finding unique sequences (printing only chromosomes over 1Mb):"
    for chrom_name,chrom_seq in genome_seq:
        if print_info and len(chrom_seq)>1000000: 
            print "  %s (%s)...  %s"%(chrom_name, basic_seq_utilities.format_base_distance(len(chrom_seq)), time.ctime())
        if chrom_name in chrom_names:
            raise Exception("%s seen twice in genome_seq!"%(chrom_name))
        chrom_names.add(chrom_name)
        # go over each slice_len-sized sequence slice and check if it's unique
        for slice_start,slice_seq in basic_seq_utilities.generate_seq_slices(chrom_seq, slice_len, step=1):
            N_total_slices += 1
            # make sure everything's uppercase
            slice_seq = slice_seq.upper()
            # if there are Ns in the sequence, just consider it unmappable
            #  MAYBE-TODO do something more complicated?  Technically sequences with one N might be unique,
            #  line ANA if there are no AAA, AGA, ACA or ATA in the whole genome...  Ignoring that for now.  
            #  I'm not even sure if bowtie ever considers genome sequences with Ns mappable...
            if 'N' in slice_seq:
                continue
            # look at both the forward and the reverse-complement sequence of the current slice, 
            #  but only store the forward (+strand) versions (no point in storing both if we're RC-ing each current sequence)
            slice_seq_RC = basic_seq_utilities.reverse_complement(slice_seq)
            # if a sequence is its own reverse-complement, it cannot be uniquely mappable - ignore
            #  (no point in even adding it to any dictionaries, since any future copies of that sequence will be caught here too)
            if slice_seq == slice_seq_RC:
                continue
            # if this is the first time slice_seq shows up, save it, and save its position
            #  (note - originally I did "{slice_seq,RC} & sequences_seen_once | sequences_seen_twice" here, HORRIBLY SLOW)
            if not (slice_seq in sequences_seen_once or slice_seq_RC in sequences_seen_once 
                    or slice_seq in sequences_seen_twice or slice_seq_RC in sequences_seen_twice):
                sequences_seen_once.add(slice_seq)
            # if this is the second time slice_seq shows up, note it as non-unique, 
            #  and remove both it AND its RC from the unique set and the unique position dict.
            #   (use set.discard(x) to avoid raising an error when one of them isn't there)
            elif (slice_seq in sequences_seen_once or slice_seq_RC in sequences_seen_once):
                sequences_seen_twice.add(slice_seq)
                sequences_seen_once.discard(slice_seq)
                sequences_seen_once.discard(slice_seq_RC)
            # if slice_seq is already known to be non-unique (in seen-twice but not seen-once), ignore
    # this is a HUGE data structure, release it as soon as possible - MAYBE-TODO not sure if that works...
    del sequences_seen_twice    

    ### Second pass - go over all chromosomes again, and save the positions of known unique sequences.
    # restart genome_seq generator
    genome_seq = original_genome_seq
    if isinstance(genome_seq, str):     genome_seq = basic_seq_utilities.parse_fasta(genome_seq)
    else:                               genome_seq = genome_seq.iteritems()
    unique_seq_positions_by_chrom = {}
    if print_info:  print "PASS 2 - getting uniquely mappable positions (printing only chromosomes over 1Mb):"
    for chrom_name,chrom_seq in genome_seq:
        if print_info and len(chrom_seq)>1000000: 
            print "  %s (%s)...  %s"%(chrom_name, basic_seq_utilities.format_base_distance(len(chrom_seq)), time.ctime())
        unique_seq_positions_by_chrom[chrom_name] = []
        # go over each slice_len-sized sequence slice and check if it's unique
        for slice_start,slice_seq in basic_seq_utilities.generate_seq_slices(chrom_seq, slice_len, step=1):
            # we already have the full set of unique sequences - now just store the position in unique_seq_positions_by_chrom
            #  if it's a unique sequence.
            if slice_seq in sequences_seen_once:
                unique_seq_positions_by_chrom[chrom_name].append(slice_start)
    # this is a HUGE data structure, release it as soon as possible - MAYBE-TODO not sure if that works...
    del sequences_seen_once

    if print_info:
        N_mappable_slices = sum([len(pos_list) for pos_list in unique_seq_positions_by_chrom.values()])
        fraction_mappable = N_mappable_slices/N_total_slices if N_total_slices else float('NaN')
        print("%.0f%% of %sbp slices are mappable (%s out of %s total on %s chromosomes)"
              %(fraction_mappable*100, slice_len, N_mappable_slices, N_total_slices, len(chrom_names)))

    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.
    unique_seq_positions_by_chrom = {chrom: numpy.array(pos_list) for (chrom, pos_list) in unique_seq_positions_by_chrom.items()}
    return unique_seq_positions_by_chrom


def genome_mappable_insertion_sites(flanking_region_length=21, mappable_slice_pos_dict=None, genome_seq=None, 
                                    include_strand=True, end_sequenced='5prime', print_info=True):
    """ Return all uniquely mappable genomic insertion sites, as dict with chrom or (chrom,strand) keys and pos_array values.

    The mappability is based on the mappability of the flanking region, using flanking_region_length, 
     on either side of the insertion (both sides are counted). 

    The actual positions are the positions of the base BEFORE the insertion, 1-based; 
     the positions are given as numpy arrays (less memory usage than lists), and are sorted. 
    For example if the position 100-120 flanking region is mappable, that means an insertion at 99-100 is mappable, 
     and one at 120-121 is mappable, depending on the orientation: 
      if the flanking region is 5prime, then a +strand 120-121 and a -strand 99-100 insertion is mappable; 
      if 3prime, the strands are inverted. 
    Each cassette orientation gets a separate entry in the position array.
    If include_strand is True, all positions are separated by strand, so the keys in the returned dictionary are
     (chromosome,strand) tuples - so if both a +strand and a -strand insertion at position X of chr1 is mappable, 
      position X will show up once in the ('chr1','+') array and once in ('chr1','-').
     If include_strand is False, the positions are merged, and the return dict keys are just chromosome names,
      so in the same example, position X will just show up twice on the 'chr1' array; 
      if position X is mappable only in one orientation, it'll show up once in the array.
     (The False option is meant for plotting, e.g. with mutant_plotting_utilities.mutant_positions_and_data)

    If mappable_slice_pos_dict is not None, it'll be assumed to be the output of genome_mappable_slices (with slice_len 
     the same as flanking_region_length), and genome_seq will be ignored; 
     if it's None, genome_mappable_slices will be used to generate that data (takes a while, and a lot of memory!).

    Genome_seq can be a list of (chr_name,chr_seq) tuples, or the path to a fasta file containing the genome.
    """
    if end_sequenced not in basic_seq_utilities.SEQ_ENDS:
        raise Exception('Invalid end_sequenced %s! Must be in %s.'%(end_sequenced, basic_seq_utilities.SEQ_ENDS))
    reads_are_reverse = False if end_sequenced=='5prime' else True
    # get a chromosome:mappable_slice_start_pos_list dict, if not provided already (this is the only use of genome_seq)
    if mappable_slice_pos_dict is None:
        if genome_seq is not None:
            mappable_slice_pos_dict = genome_mappable_slices(flanking_region_length, genome_seq, print_info)
        else:
            raise Exception('Either mappable_slice_pos_dict or genome_seq must be given!')
    # convert the mappable slice chromosome/position values into mappable insertion locations, 
    #  treating them as end_sequenced flanking regions - each flanking region should give TWO mappable positions, 
    #   one on each side, with opposite strands (with the insertion position strand depending on end_sequenced)
    mappable_position_data = collections.defaultdict(list)
    for chrom,pos_list in mappable_slice_pos_dict.iteritems():
        if print_info and len(pos_list)>500000: 
            print "  %s (%s mappable slices)...  %s"%(chrom, len(pos_list), time.ctime())
        for pos in pos_list:
            for strand in basic_seq_utilities.SEQ_STRANDS:
                flanking_region_pos = (chrom, pos, pos+flanking_region_length-1, strand)
                # the easiest way of getting the right insertion position is to just use the same function I normally use
                #  for making actual Insertion_position objects from sequenced flanking region position data
                insertion_pos_object = mutant_analysis_classes.get_insertion_pos_from_flanking_region_pos(flanking_region_pos, 
                                                                                              end_sequenced, reads_are_reverse)
                if include_strand:
                    mappable_position_data[chrom,insertion_pos_object.strand].append((insertion_pos_object.min_position))
                else:
                    mappable_position_data[chrom].append(insertion_pos_object.min_position)

    # make sure the data is sorted (MAYBE-TODO does it really matter?)
    mappable_position_data = general_utilities.sort_lists_inside_dict(mappable_position_data)
    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.  MAYBE-TODO or could I do it in-place somehow?  
    mappable_position_data = {key: numpy.array(pos_list) for (key, pos_list) in mappable_position_data.items()}
    return mappable_position_data


def genome_mappable_insertion_sites_multi(flanking_region_lengths=[20,21], mappable_slice_pos_dicts=None, genome_seq=None, 
                                          include_strand=True, end_sequenced='5prime', print_info=True):
    """ Run genome_mappable_insertion_sites multiple times, with different flank length and optionally mappable dict values.

    Flanking_region_lengths should be a list of flanking_region_length values; 
     mappable_slice_pos_dicts can be None, or a same-length list of values. 
    For each pair of those values, genome_mappable_insertion_sites will be run with the matching arguments. 
    The final output will include the sum of position lists from all runs. 

    See genome_mappable_insertion_sites docstring for more info.
    """
    if mappable_slice_pos_dicts == None:    mappable_slice_pos_dicts = [None for _ in flanking_region_lengths]
    # need to start with lists, not arrays, because arrays can't be appended/extended to
    full_mappable_position_data = collections.defaultdict(list)
    # for each set of inputs, run genome_mappable_insertion_sites to generate new data, 
    #  and then merge the new data into full_mappable_position_data
    for flanking_region_length,mappable_slice_pos_dict in zip(flanking_region_lengths, mappable_slice_pos_dicts):
        # get data for 
        curr_mappable_position_data = genome_mappable_insertion_sites(flanking_region_length, mappable_slice_pos_dict, genome_seq, 
                                                                      include_strand, end_sequenced, print_info)
        for chrom,new_pos_list in curr_mappable_position_data.iteritems():
            full_mappable_position_data[chrom].extend(new_pos_list)
    # make sure the data is sorted (MAYBE-TODO does it really matter?) (also converts from defaultdict to dict)
    full_mappable_position_data = general_utilities.sort_lists_inside_dict(full_mappable_position_data)
    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.  
    full_mappable_position_data = {key: numpy.array(pos_list) for (key, pos_list) in full_mappable_position_data.items()}
    return full_mappable_position_data


### mappability data by gene

Gene = collections.namedtuple('Gene', 'record name start end strand bin_start bin_end feature_iter')
Feature = collections.namedtuple('Feature', 'type start end')

def _next_feature(feature_iterator, curr_pos):
    """ Return the first feature in feature_iterator that ends after curr_pos, or None if there aren't any. 

    This assumes the feature does need to be advanced at least once (i.e. that curr_pos is over the current feature end).
    """
    feature = Feature('blank', 0, 0)
    while curr_pos is None or curr_pos > feature.end:
        try:
            feature = feature_iterator.next()
        except StopIteration:
            feature = None
            break
    return feature

def _next_gene_for_mapp(gene_iterator, curr_mapp_pos, gene_mappable_lengths, gene_mappable_lengths_binned, 
                        N_bins, exclude_UTRs_in_bins, skip_multiple_splice_variants=True):
    """ Advance the gene until curr_mapp_pos <= gene_end, and return (gene_record,ID,start,end); don't catch StopIteration. 

    The point of this is that curr_mapp_pos should be the next mappable position, and the current position of gene_iterator
     is around the previous mappable position - so advance gene_iterator to the first gene that has any mappable positions, 
     and for all genes in between, just add them to gene_mappable_lengths and gene_mappable_lengths_binned 
      with all mappable lengths set to zero.
    If curr_mapp_pos is None, assume there are no more mappable positions until the end of gene_iterator 
     (which is probably the chromosome), so go to the end.

    This also advances the feature_iterator of the gene, likewise going on to the first feature that overlaps curr_mapp_pos.

    If skip_multiple_splice_variants is True, genes with multiple splice variants are treated as follows:
        - STILL COUNTED for total gene mappable lengths, AND for binned mappable lengths IF exclude_UTRs_in_bins is False, 
            since splice variants aren't relevant to that!
        - SKIPPED for the purposes of binned mappable lengths if exclude_UTRs_in_bins is True: 
            bin_start and bin_end are set to gene_start so that the binned length is 0
        - SKIPPED for feature mappable length - feature_iterator is set to contain a single gene-length feature 
            named MULTIPLE_SPLICE_VARIANTS, which can be thrown away later if desired.
      If there's a gene with multiple splice variants and skip_multiple_splice_variants is False, an exception will be raised.

    This assumes the gene does need to be advanced at least once (i.e. that curr_mapp_pos > current gene end).
    """
    # if no curr_mapp_pos given, just go over all the genes
    gene_end = 0
    while curr_mapp_pos is None or curr_mapp_pos > gene_end:
        gene_record = gene_iterator.next()
        gene_name = gene_record.id
        gene_start, gene_end = gff_examine_file.get_feature_start_end(gene_record)
        # get the start/end of the length we want binned mappabilies of
        if exclude_UTRs_in_bins is False:   bin_start, bin_end = gene_start, gene_end
        else:                               
            try:                            bin_start, bin_end = gff_examine_file.get_gene_start_end_excluding_UTRs(gene_record)
            except MultipleRNAError:        
                # if "skipping" a gene due to multiple splice variants, just make a fake zero-length region for binning,
                #   since we don't want to skip it entirely (we want to add it to gene_mappable_lengths!), and at the end 
                #   to make a full Gene object, bin_start and bin_end are required.
                #   doing features, to make a 
                if skip_multiple_splice_variants:   bin_start, bin_end = gene_start, gene_start
                else:                               raise
        # when starting a new gene, set its mappable total and binned lengths to 0
        gene_mappable_lengths[gene_name] = 0
        gene_mappable_lengths_binned[gene_name] = [0 for _ in range(N_bins)]
    # now create the final gene record, including a feature iterator
    if len(gene_record.sub_features) == 0:  
        raise NoRNAError("Gene %s has no RNA - can't determine feature mappability!"%gene_name)
    # if we have a gene with one splice variant, make a feature iterator for it. 
    if len(gene_record.sub_features) == 1:   
        feature_iterator = iter(sorted([Feature(feature.type, *gff_examine_file.get_feature_start_end(feature)) 
                                        for feature in gene_record.sub_features[0].sub_features], key=lambda feature: feature.start))
    # if multiple splice variants, raise exception if desired, otherwise make the feature iterator be a fake single feature.
    else:
        if skip_multiple_splice_variants:
            feature_iterator = iter([Feature('MULTIPLE_SPLICE_VARIANTS', gene_start, gene_end)])
        else:
            raise MultipleRNAError("Gene %s has multiple RNAs - can't determine single CDS start/end!"%gene_name)
    gene_strand = gff_examine_file.GFF_strands[gene_record.strand]
    gene = Gene(gene_record, gene_name, gene_start, gene_end, gene_strand, bin_start, bin_end, feature_iterator)
    feature = _next_feature(feature_iterator, curr_mapp_pos)
    return gene, feature


def _save_gene_bin_mappability(mappable_positions, gene, gene_mappable_lengths_binned, N_bins):
    """ Divide mappable_positions list into N bins, save bin counts in gene_mappable_lengths_binned[gene_name]. """
    mappable_lengths_binned = list(numpy.histogram(mappable_positions, bins=N_bins, range=(gene.bin_start,gene.bin_end))[0])
    # if the gene is -strand, the order of the bins should just be reversed
    if gene.strand == '-': 
        mappable_lengths_binned.reverse()
    gene_mappable_lengths_binned[gene.name] = mappable_lengths_binned


def gene_mappability_bins_features(mappability_by_flanking_region, genefile=mutant_utilities.DEFAULT_GENE_POS_FILE, 
                                   split_into_bins=20, exclude_UTRs_in_bins=False):
    """ Calculate mappable gene lengths (total and binned by length) and overall mappable lengths of different gene features.

    Three outputs: 
        1) a gene_name:total_mappable_length dict
        2) a gene_name:list_of_mappable_lengths_by_bin dict - divide the gene length into N equal bins, give their mappable lengths.
            Normally the full mRNA length is binned; if exclude_UTRs_in_bins is True, the region binned will be only
             from the start of the first CDS feature to the end of the last CDS feature.
        3) a seq_type:total_mappable_length dict, with seq_type each of 
                CDS/intron/5'UTR/3'UTR/gene/intergenic (same keys as in gff_examine_file.feature_total_lengths output)
    (all these things are done in parallel instead of by separate functions because it's faster that way, if not very clean)
    If you just want the first output, run gene_mappability instead.

    The mappability_by_flanking_region should be a dict generated by genome_mappable_insertion_sites/_multi.

    The mappable lengths correspond to real lengths (will be equal to real lengths if all positions are mappable), 
     but can be fractions if some position is mappable only for some flanking region lengths or strands.
    """
    if not split_into_bins > 1:  raise Exception("split_into_bins must be >1!")
    gene_mappable_lengths = {}
    gene_mappable_lengths_binned = {}
    feature_total_mappable_lengths = {'CDS':0, 'five_prime_UTR':0, 'three_prime_UTR':0, 'MULTIPLE_SPLICE_VARIANTS':0}
    # define some convenience functions to avoid re-typing things - note that these don't use the given variables as defaults, 
    #  so they're NOT closures, they don't store the at-definition value of the given variables, 
    #   but just using their at-call-time values as if they were globals (which they are), as far as I can tell.
    next_gene = lambda (curr_pos): _next_gene_for_mapp(genes_by_start_pos, curr_pos, gene_mappable_lengths, 
                                                       gene_mappable_lengths_binned, split_into_bins, exclude_UTRs_in_bins)
    save_gene_binned_mappability = lambda: _save_gene_bin_mappability(curr_gene_mappable_positions, curr_gene, 
                                                                      gene_mappable_lengths_binned, split_into_bins)
    # go over the gff genefile by chromosome
    # MAYBE-TODO this currently mixes GFF parsing and getting the mappability data - would it be better to separate them?
    #  I could write a gff_examine_file.py function to get the positions of all the features, 
    #   but dealing with that would probably be more complicated than just doing it the current way...
    # chromosomes that aren't in GENEFILE presumably have no genes - no need to care about them.
    with open(os.path.expanduser(genefile)) as GENEFILE:
        for chromosome_record in gff_examine_file.GFF.parse(GENEFILE):
            chromosome = chromosome_record.id
            # if this chromosome has no genes, go on to the next one - we only care about genes here.
            if not len(chromosome_record.features):
                continue
            genes_by_start_pos = iter(sorted(chromosome_record.features, key=lambda gene_record:gene_record.location.start.position))
            # make a single list of mappable positions, collapsing all flanking region lengths and strands together
            chromosome_mappable_positions = []
            for mappability_dict in mappability_by_flanking_region.values():
                # convert the values to lists instead of numpy arrays (MAYBE-TODO this takes more memory, but it's easier...)
                try:    
                    chromosome_mappable_positions += list(mappability_dict[chromosome])
                # using dict.get to return an empty list if that chromosome/strand isn't in the dict at all
                except KeyError:
                    chromosome_mappable_positions += list(mappability_dict.get((chromosome, '+'), []))
                    chromosome_mappable_positions += list(mappability_dict.get((chromosome, '-'), []))
            # MAYBE-TODO get separate sense/antisense mappable lengths? But it can't be much different, 
            #  since only the 20-21bp on gene/feature edges could make a difference - probably not worth it.
            ### go over mappable positions concurrently, iterating over all genes/features by hand in parallel, 
            #    adding mappable length to all the output values
            # start the first gene
            curr_gene, curr_feature = next_gene(curr_pos=1)
            curr_gene_mappable_positions = []
            for curr_mapp_pos in sorted(chromosome_mappable_positions):
                # if the position is past the end of the gene, go on to the next gene (and the next etc if needed)
                #  (if there are no more genes, break out of the "while True" loop and go on to the next chromosome)
                if curr_mapp_pos > curr_gene.end:
                    # save the current gene binned mappability data, since we're done with it
                    save_gene_binned_mappability()
                    # go on to next gene
                    try:                    curr_gene, curr_feature = next_gene(curr_mapp_pos)
                    except StopIteration:   break
                    curr_gene_mappable_positions = []
                # if the position is BEFORE the start of the current gene, just go on to the next position
                if curr_mapp_pos < curr_gene.start:   continue
                # if the position is inside the gene, add it to the mappability counts by gene/bin/feature
                if curr_gene.start <= curr_mapp_pos <= curr_gene.end:
                    # add it to the gene total mappability
                    gene_mappable_lengths[curr_gene.name] += 1
                    # add it to the binned mappability data if it's inside the range we want to bin
                    if curr_gene.bin_start <= curr_mapp_pos <= curr_gene.bin_end:
                        curr_gene_mappable_positions.append(curr_mapp_pos)
                    # iterate over features and positions in parallel - similar to genes:
                    if curr_mapp_pos > curr_feature.end:
                        curr_feature = _next_feature(curr_gene.feature_iter, curr_mapp_pos)
                    if curr_feature is None or curr_mapp_pos < curr_feature.start:
                        continue
                    if curr_feature.start <= curr_mapp_pos <= curr_feature.end:
                        feature_total_mappable_lengths[curr_feature.type] += 1
            ### if there are no more mappable positions, all remaining genes should get a 0 mappable length (full and binned)
            # save the current gene binned mappability data, since we're done with it
            save_gene_binned_mappability()
            # skip all the next genes (by using None as the current mappable position)
            try:                    curr_gene, curr_feature = next_gene(curr_pos=None)
            except StopIteration:   pass

    # add 'gene', 'intron', 'intergenic' and 'all' categories to feature_total_mappable_lengths, calculating from existing data
    # LATER-TODO this is still ignoring UTR introns, i.e. counting them as introns!  Should fix that, but it's complicated...  
    # Well, UTR introns aren't that difficult from other introns - and anyway hopefully there isn't a lot of them.
    feature_total_mappable_lengths['all'] = sum( sum( len(chrom_mappable_positions) 
                                                      for chrom_mappable_positions in mappability_dict.values() ) 
                                                 for mappability_dict in mappability_by_flanking_region.values() )
    feature_total_mappable_lengths['gene'] = sum(gene_mappable_lengths.values())
    feature_total_mappable_lengths['intergenic'] = feature_total_mappable_lengths['all'] - feature_total_mappable_lengths['gene']
    feature_total_mappable_lengths['intron'] = feature_total_mappable_lengths['gene'] - sum(feature_total_mappable_lengths[x] 
                                                    for x in "CDS five_prime_UTR three_prime_UTR MULTIPLE_SPLICE_VARIANTS".split())
    # MAYBE-TODO remove the throwaway MULTIPLE_SPLICE_VARIANTS category?
    # Convert all the results from counts of mappable positions to actual mappable lengths (<= real lengths).
    #  Note that all the mappable lengths are actually just the counts of mappable positions right now: 
    #   so if mappability_by_flanking_region was generated with one flanking region length, 
    #    the maximum is 2 times the gene length in bp (one potential mappable position for each strand), 
    #   and if it was generated with N flanking region lengths, 
    #    the maximum is 2N times the gene length (one per strand and flanking region length). 
    #  So divide all the values by 2N to get the actual mappability in the same units as real length.
    max_mappable_per_base = 2 * len(mappability_by_flanking_region)
    if max_mappable_per_base>0:
        gene_mappable_lengths = {gene: mapp_len/max_mappable_per_base for (gene, mapp_len) in gene_mappable_lengths.iteritems()}
        gene_mappable_lengths_binned = {gene: [mapp_len/max_mappable_per_base for mapp_len in mapp_len_list] 
                                        for (gene, mapp_len_list) in gene_mappable_lengths_binned.iteritems()}
        feature_total_mappable_lengths = {feature: mapp_len/max_mappable_per_base 
                                          for (feature, mapp_len) in feature_total_mappable_lengths.iteritems()}
    return gene_mappable_lengths, gene_mappable_lengths_binned, feature_total_mappable_lengths
    # TODO really all this mappability isn't quite right, because I'm giving equal weight to 20bp an 21bp flanking regions,
    #  but in reality more mutants have 21bp ones... Fix that!


def gene_mappability(mappability_by_flanking_region, exclude_UTRs=False, genefile=mutant_utilities.DEFAULT_GENE_POS_FILE, 
                        strip_from_genename_ends=[], verbose=False, DEBUG_gene='', DEBUG_chrom=''):
    """ Calculate mappable gene lengths based on genome mappability data and a gene position file.

    Outputs a gene_name:total_mappable_length dictionary.  The mappable lengths correspond to real lengths (will be equal to real 
    lengths if all positions are mappable), but can be fractions if some position is mappable only for some versions or strands.

    If exclude_UTRs, the mappable lengths are calculated excluding UTRs; for genes with multiple splice variants, it.
        calculates it based the start of the first exon and the end of the last one.

    The mappability_by_flanking_region should be a dict with flanking region lengths as keys, 
        and values as generated by genome_mappable_insertion_sites/_multi.

    Optionally strips strings from strip_from_genename_ends from filename ends.
    If verbose, prints a summary of what it's doing (per chromosome)
    """
    # grab gene positions (including or excluding UTRs; if excluding, take longest mRNA), and make into a chrom:gene:(start,end) dict
    if verbose: print "getting gene positions..."
    gene_positions = gff_examine_file.gene_positions(genefile, include_strand=False, 
                                                     coding_only=exclude_UTRs, return_longest_if_multiple=True)
    gene_positions_by_chromosome = collections.defaultdict(dict)
    for (gene, (chrom, start, end)) in gene_positions.items():
        gene_positions_by_chromosome[chrom][gene] = (start, end)
    gene_mappable_lengths = {}
    # now go over each chromosome (to save memory)
    for chromosome, curr_gene_positions in gene_positions_by_chromosome.items():
        if verbose: print "getting mappable positions for %s..."%chromosome
        # make a counter of mappable positions on this chromosome 
        chrom_mappos_counts = collections.Counter()
        for mappability_dict in mappability_by_flanking_region.values():
            try:    
                chrom_mappos_counts.update(mappability_dict[chromosome])
            # using dict.get to return an empty list if that chromosome/strand isn't in the dict at all
            except KeyError:
                chrom_mappos_counts.update(mappability_dict.get((chromosome, '+'), []))
                chrom_mappos_counts.update(mappability_dict.get((chromosome, '-'), []))
        if chromosome==DEBUG_chrom:   print "mappability counter", chromosome, chrom_mappos_counts
        if verbose: print "getting mappable gene lengths for %s genes on %s..."%(len(curr_gene_positions), chromosome)
        for (gene, (start, end)) in curr_gene_positions.items():
            # fix gene name if needed
            for string in strip_from_genename_ends:
                if gene.endswith(string):   gene = gene[:-len(string)]
            # find the number of mappable positions inside the gene by checking how many times they show up in the mappable counter, 
            #  then divide by the number of times they COULD show up: 2 strands, N versions.
            gene_positions = set(range(start, end+1))
            try:
                gene_mappable_lengths[gene] = sum(chrom_mappos_counts[p] for p in gene_positions)\
                                                / len(mappability_by_flanking_region) / 2
                if gene==DEBUG_gene:  
                    print {p: chrom_mappos_counts[p] for p in gene_positions if chrom_mappos_counts[p]},\
                          len(mappability_by_flanking_region)*2
            except ZeroDivisionError:
                gene_mappable_lengths[gene] = 0.0
    return gene_mappable_lengths
    # TODO add option to split the gene length into N bins and get mappable lengths for each bin separately?


def feature_mappability(mappability_by_flanking_region, genefile=mutant_utilities.DEFAULT_GENE_POS_FILE, 
                        genome_fasta_file=gff_examine_file.DEFAULT_NUCLEAR_GENOME_FILE, verbose=False):
    """ Calculate mappable length of all gene features - return feature:mappable_length dictionary.

    Position in multiple features (either because of multiple splice variants or overlapping genes) are categorized as 'multiple'.
    Positions in the SAME feature in multiple splice variants or overlapping genes are categorized as that feature and counted once.
    Introns, UTR introns, etc are calculated based on which positions aren't in a UTR/exon, for each mRNA.
    All features are mutually exclusive, except for 'gene' (which includes all features).
    """
    # get the total chromosome lengths from genome_fasta_file if given - needed to calculate total intergenic length 
    feature_mappable_lengths = collections.defaultdict(int)
    if verbose: print "parsing fasta file to get chromosome lengths..."
    with open(genome_fasta_file) as FASTAFILE:  fasta_seq_dict = SeqIO.to_dict(SeqIO.parse(FASTAFILE, "fasta"))
    with open(os.path.expanduser(genefile)) as GENEFILE:
        for chromosome_record in gff_examine_file.GFF.parse(GENEFILE, base_dict=fasta_seq_dict):
            chromosome = chromosome_record.id
            ######### Make a counter of mappable positions on this chromosome 
            if verbose: print "getting mappable positions for %s..."%chromosome
            chrom_mappos_counts = collections.Counter()
            for mappability_dict in mappability_by_flanking_region.values():
                try:    
                    chrom_mappos_counts.update(mappability_dict[chromosome])
                # using dict.get to return an empty list if that chromosome/strand isn't in the dict at all
                except KeyError:
                    chrom_mappos_counts.update(mappability_dict.get((chromosome, '+'), []))
                    chrom_mappos_counts.update(mappability_dict.get((chromosome, '-'), []))
            # if chromosome is unmappable (or wasn't included in mapping for test purposes), no point in going into it
            if not sum(chrom_mappos_counts.values()):   continue
            ######### get the joint positions of all the features
            if verbose:     print "getting overall feature positions for %s..."%chromosome
            # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3)
            #   my mappability data uses 1-based positions, so convert to that.
            feature_positions = collections.defaultdict(set)
            # for each gene AND each mRNA/feature set, add positions to the right set - I'll sort them out later. 
            for gene_record in chromosome_record.features:
                if len(gene_record.sub_features) == 0:
                    raise NoRNAError("Gene %s has no mRNA - can't determine CDS start/end!"%gene_record.id)
                for mRNA in gene_record.sub_features:
                    if mRNA.location.start.position != gene_record.location.start.position:
                        feature_positions['outside_mRNA'] |= set(range(gene_record.location.start.position+1, 
                                                                       mRNA.location.start.position+1))
                    if mRNA.location.end.position != gene_record.location.end.position:
                        feature_positions['outside_mRNA'] |= set(range(mRNA.location.end.position+1, 
                                                                       gene_record.location.end.position+1))
                    curr_mRNA_positions = set(range(mRNA.location.start.position+1, mRNA.location.end.position+1))
                    curr_all_feature_positions = set()
                    for feature in mRNA.sub_features:
                        curr_feature_positions = set(range(feature.location.start.position+1, feature.location.end.position+1))
                        curr_all_feature_positions |= curr_feature_positions
                        feature_positions[feature.type] |= curr_feature_positions
                    # Now I need to deal with features that are implied rather than directly annotated: introns, UTR introns, other
                    #   This method is basically copied from mutant_IB_RISCC_classes.find_gene_by_pos_gff3
                    intron_positions = curr_mRNA_positions - curr_all_feature_positions
                    CDS_features = [feature for feature in mRNA.sub_features if feature.type=='CDS']
                    for position in intron_positions:
                        if position < min([feature.location.start.position+1 for feature in mRNA.sub_features]):
                            if gene_record.strand==1:    implied_feature = 'mRNA_before_exons'
                            elif gene_record.strand==-1: implied_feature = 'mRNA_after_exons'
                        elif position > max([feature.location.end.position for feature in mRNA.sub_features]):
                            if gene_record.strand==1:    implied_feature = 'mRNA_after_exons'
                            elif gene_record.strand==-1: implied_feature = 'mRNA_before_exons'
                        elif position < min([feature.location.start.position+1 for feature in CDS_features]):
                            if gene_record.strand==1:    implied_feature = "5'UTR_intron"
                            elif gene_record.strand==-1: implied_feature = "3'UTR_intron"
                        elif position > max([feature.location.end.position for feature in CDS_features]):
                            if gene_record.strand==1:    implied_feature = "3'UTR_intron"
                            elif gene_record.strand==-1: implied_feature = "5'UTR_intron"
                        else:
                            implied_feature = 'intron'
                        feature_positions[implied_feature].add(position)
            if verbose:     print {feature: len(positions) for (feature, positions) in feature_positions.items()}
            # At the end any position that's only in one feature set stays there, 
            #   and any that's in multiple sets becomes multiple_splice_variants (includes overlapping gene cases)
            multi_positions = set()
            for (feature1, positions1) in feature_positions.items():
                for (feature2, positions2) in feature_positions.items():
                    if feature1 != feature2:
                        overlap = positions1 & positions2
                        multi_positions |= overlap
            for (feature, positions) in feature_positions.items():  positions -= multi_positions
            feature_positions['multiple'] = multi_positions
            if verbose:     print {feature: len(positions) for (feature, positions) in feature_positions.items()}
            ######### Compare the feature positions with mappable position sets to calculate mappable feature lengths.
            if verbose: print "getting mappable feature lengths for %s..."%(chromosome)
            for (feature, position_set) in feature_positions.items():
                # find the number of mappable feature positions by checking how many times they show up in the mappable counter, 
                #  then divide by the number of times they COULD show up: 2 strands, N versions.
                feature_mappable_lengths[feature] += sum(chrom_mappos_counts[p] for p in position_set)\
                                                        / len(mappability_by_flanking_region) / 2
    ######### To save memory, instead of keeping full gene/intergenic/all position lists, just calculate the lengths now
    feature_mappable_lengths['gene'] = sum(feature_mappable_lengths.values())
    genome_mappable_length = sum(sum(len(x) for x in chrom_positions.values()) for chrom_positions 
                                 in mappability_by_flanking_region.values()) / len(mappability_by_flanking_region) / 2
    feature_mappable_lengths['intergenic'] = genome_mappable_length - feature_mappable_lengths['gene']
    return feature_mappable_lengths
    # TODO unit-tests!



############################# Using mappability to finding hotspots/coldspots in real data #############################

### STATISTICAL NOTES (originally from ~/experiments/mutant_pool_screens/1211_positions_Ru-screen1-for-paper/notes.txt)
    # Basic idea - trying to figure out the specific locations of statistically significant hotspots and coldspots:
    # Go over the genome in a sliding window of some size, see how many mutants are there and whether that's statistically significantly different than expected from a random distribution of the same density, INCUDING MAPPABILITY DATA (easiest to do by just comparing to a large number of random datasets, probably).  Figure out if we have any clear hot or cold spots - then we can see if they have anything in common.
    #
    # What statistical test to use for each window?  What we're comparing is the mutant count in the window, X, against one of two things:
    #  A) the mappability % of that window
    #  B) the mutant count in that window for each of the 1000 simulated datasets
    # Which of these things would be the better choice?  Probably A, really, since B is just a random simulation that's entirely based on A.  In either case, REMEMBER TO DO FALSE DISCOVERY RATE CORRECTION after getting the p-values!
    #
    # A) How would I go about getting the p-value for A?  Basically I want to know the probability of getting X mutants in that window randomly, when inserting N mutants in the genome.  That should be mathematically pretty straightforward, right?  All I need is the total mappable size of the genome, and the total mappable size of that window, and I should be able to calculate the actual distribution... That's the binomial distribution - " the discrete probability distribution of the number of successes in a sequence of n independent yes/no experiments, each of which yields success with probability p."  
    # Trying the binomial distribution: the corresponding statistical significance test is the binomial test (https://en.wikipedia.org/wiki/Binomial_test), scipy.stats.binom_test in python.  Say we want the p-value for getting 25 mutants in a 20kb window, out of a total of 12061 mutants, with the probability around 20kb divided by the genome size (really it should be the mappable lengths of that window and the genome, not the total lengths, but let's go with the approximation for now).  The genome size is 113MB, so:
    #     >>> import scipy.stats
    #     >>> scipy.stats.binom_test(25,12061,20000/113000000)
    #     1.3921234115131231e-18
    # (That's pretty significant!  But we'll see how it looks with FDR correction, of course.)
    # How long does this take?  For 50k it takes <1min, and the whole genome in 20bk windows, 113Mb/20kb, is about 5k - not bad.
    #     >>> time.ctime(); _ = [ scipy.stats.binom_test(random.randrange(30),12061,20000/113000000) for x in range(50000)]; time.ctime()
    #     'Sun Dec  2 19:03:25 2012'
    #     'Sun Dec  2 19:04:12 2012'
    # With a reasonably high N (around 1000), the binomial distribution can be modeled by the Poisson distribution (https://en.wikipedia.org/wiki/Poisson_distribution#Derivation_of_Poisson_distribution_.E2.80.94_The_law_of_rare_events).  My N is the number of mutants in the dataset, which is 12k, so I could do that to speed things up - MAYBE-TODO. 
    #
    # B) Alternatively, for option B - I'm trying to test whether the mean of two distributions (real and simulated mutant counts over a window) are the same or different...  So T-test (scipy.stats.ttest_ind, I think), or Welch's t-test (doesn't assume equal variances - apparently present in a newer scipy version), or Mann-Whitney U test (doesn't assume the data is normally distributed).  I'd say it's pretty likely to be randomly distributed, idk about the variances...  Well, I'm only giving a single number for the real dataset, so there's no variance in that - and if anything, the variance of the real dataset may be smaller than of the simulated ones, not larger, so I think we're good with the Student's T-test.
    #
    # How to do FDR correction?  See "OLD NOTES ON FDR CORRECTION" section in ~/programs/basic_work_mine/statistics_utilities.py


def _get_chrom_either_strand(input_val):
    """ input_val can be either just a chrom string, or a (chrom,strand) tuple - return chrom. """
    try:                            
        if len(input_val)==2:   input_val = input_val[0]
    except TypeError:               raise Exception("Input should be either (chrom,strand) or chrom, not %s!"%input_val)
    if not type(input_val)==str:    raise Exception("Input should be either (chrom,strand) or chrom, not %s!"%input_val)
    return input_val


def N_mutants_in_dataset(dataset, all_chromosomes=None):
    """ Number of mutants in dataset, counting only ones in all_chromosomes if given (otherwise counting all mutants).

    Dataset can be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a chrom:position_list or (chrom,strand):position_list dictionary.
    """
    if all_chromosomes is None:     chromosome_test = lambda chrom: True
    else:                           chromosome_test = lambda chrom: (chrom in all_chromosomes)
    try:
        return len([1 for m in dataset if chromosome_test(m.position.chromosome)])
    except AttributeError:
        return sum([len(position_list) for key,position_list in dataset.items() if chromosome_test(_get_chrom_either_strand(key))])
    # TODO unit-test!


def find_average_mutant_count(dataset, window_size, chromosome_lengths=None):
    """ Average mutant count per window (N_mutants / genome_size * window_size) for dataset. 

    Only counts mutants in chromosomes given in chromosome_lengths; if None, takes the full chlamy genome.

    Dataset can be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a chrom:position_list or (chrom,strand):position_list dictionary.
    """
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
    all_chromosomes = set(chromosome_lengths.keys())
    N_mutants = N_mutants_in_dataset(dataset, all_chromosomes)
    total_genome_length = sum(chromosome_lengths.values())
    return N_mutants / total_genome_length * window_size
    # TODO unit-test!


def find_hot_cold_spots(dataset, window_size, mappable_positions_20bp, mappable_positions_21bp, window_offset=0, fraction_20bp=None, 
                        N_samples=None, chromosome_lengths=None, quiet=False):
    """ Find statistically significant hotspots and coldspots in dataset mutant positions, based on genome mappability. 

    Divides the genome into window_size-sized windows, and for each of them compares 
     the number of mutants in that window in dataset to the number expected using the BINOMIAL DISTRIBUTION, 
      given the total number of mutants in dataset and the genome mappability information
      (the probability of a mutant landing in the window is the number of mappable positions in the window 
       divided by the number of mappable positions in the entire genome - these numbers will be different for 20bp 
       and 21bp mappability, so calculate for both, and get a weighted average based on how many mutants in dataset
       have only 20bp flanking region sequences.)
      FALSE DISCOVERY RATE CORRECTION using the Benjamini-Hochberg method is used on all p-values, 
       with the total number of windows tested used as N_samples, unless another N_samples is provided
       (which you may or may not want to do if doing multiple find_hot_cold_spots runs 
        with different window sizes/offsets over the same dataset: you'll be under-correcting if you don't, 
        but since the results for overlapping windows are positively correlated, you'll be over-correcting if you do!)

    Window_offset defines at which position of each chromosome to start: if a chromosome is 500 long and window_size is 200, 
     - with window_offset 0, the windows checked will be 1-200 and 201-400;
     - with window_offset 100, the windows checked will be 101-300 and 301-500.
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default chlamy file will be used) - need the lengths because the dataset positions and mappability data
     don't give the end position of each chromosome.

    Return a 5-tuple of items: FDR-adjusted p-values, raw p-values, if_coldspot (True if fewer mutants than expected, False if more),
      raw mutant counts, and expected mutant counts (based on mappability of given window and total genome mutant count), 
      with each of these items being a chromosome:list_of_values_per_window dictionary.
     The caller/recipient must take care of keeping track of the window size and offset in order for this data to be meaningful.

    Fraction_20bp gives the fraction of mutants that only have 20bp flanking regions (no 21bp ones) - this matters for mappability
     calculations, because the mappability of 20bp and 21bp flanking regions is different.

    Dataset can be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, or a chrom:position_list dictionary; 
     if it's the former, fraction_20bp can calculated from the dataset itself; if the latter, it must be given explicitly. 

    The two mappable_positions_* arguments should be chrom:position_list or (chrom,strand):position_list dictionaries 
     giving all the mappable positions, as generated by genome_mappable_insertion_sites with either include_strand value. 
    """
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
    if fraction_20bp is None:
        fraction_20bp = get_20bp_fraction(dataset)
    total_N_mutants = N_mutants_in_dataset(dataset)
    mappable_position_data = {20: mappable_positions_20bp, 21: mappable_positions_21bp}
    # get total mappable genome lengths (for 20 and 21bp cases separately) 
    genome_mappable_lengths = {flank_len: sum(len(pos_list) for chrom_strand,pos_list in mappable_positions.items() 
                                              if _get_chrom_either_strand(chrom_strand) in chromosome_lengths)
                               for (flank_len,mappable_positions) in mappable_position_data.items()}
    # print real/mappable lengths as a sanity check - divide the mappable lengths by two because both strands are counted
    if not quiet:
        print "Genome total length %s; mappable length %s for 20bp flank and %s for 21bp"%(sum(chromosome_lengths.values()), 
                                                                   genome_mappable_lengths[20]/2, genome_mappable_lengths[21]/2)
    # get mappable lengths per window, for 20 and 21bp cases separately, with the right window size and offset; 
    #  don't include special_last_bin to avoid complications (so the last part of each chromosome that doesn't 
    #   fit the window size will be ignored - same as the first part before offset)
    #   (chromosomes smaller than window_size+offset will be ignored entirely!)
    window_kwargs = {'bin_size':window_size, 'chromosome_lengths':chromosome_lengths, 'chromosomes':chromosome_lengths.keys(), 
                     'first_bin_offset': window_offset, 'special_last_bin':False}
    window_mappable_length_lists = {flank_len: get_histogram_data_from_positions(mappable_positions, **window_kwargs) 
                                    for (flank_len,mappable_positions) in mappable_position_data.items()}
    # similarly, get mutant counts per window (ignoring the part before the offset and the small part after the last full window)
    #  (if dataset isn't already a chrom:position_list dict, converit it to that first)
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        dataset_positions = get_mutant_positions_from_dataset(dataset)
    else:
        dataset_positions = dataset
    window_mutant_counts = get_histogram_data_from_positions(dataset_positions, **window_kwargs)
    total_N_windows = sum(len(len_list) for len_list in window_mutant_counts.values())
    if N_samples is None:   N_samples = total_N_windows
    # average mutants per window should only count the mutants covered by any of the windows 
    #  (exclude ones skipped due to non-zero offset at chromosome start or a non-full window at chromosome end)
    average_mutants_per_window = find_average_mutant_count(dataset, window_size, chromosome_lengths)
    average_mutants_per_window_2 = sum(sum(counts) for counts in window_mutant_counts.values()) / total_N_windows
    if not quiet:
        print "%s window (offset %s) - average %.3g or %.3g mutants per window "%(format_bp(window_size), format_bp(window_offset), 
                                                                          average_mutants_per_window, average_mutants_per_window_2)
    window_raw_pvalues = {}
    window_adjusted_pvalues = {}
    window_if_coldspot = {}
    window_expected_mutant_counts = {}
    for chrom, window_mutant_count_list in window_mutant_counts.items():
        # get the probability of a mutant landing in the window, 
        #  given the mappable lengths of the window and the total genome for 20 and 21bp cases, 
        #  and the fraction of mutants that is 20bp-only.
        # note that all the values in the histogram dicts are numpy arrays, not lists, so they can be operated on directly
        window_probabilities_20 = window_mappable_length_lists[20][chrom] / genome_mappable_lengths[20]
        window_probabilities_21 = window_mappable_length_lists[21][chrom] / genome_mappable_lengths[21]
        window_probabilities = window_probabilities_20 * fraction_20bp + window_probabilities_21 * (1-fraction_20bp)
        window_expected_mutant_count_list = window_probabilities * total_N_mutants
        window_expected_mutant_counts[chrom] = window_expected_mutant_count_list
        window_if_coldspot[chrom] = window_mutant_count_list < window_expected_mutant_count_list
        # calculate p-value for each window, using the binomial distribution 
        #  based on the window probability and window mutant number; convert to numpy array again
        raw_pvalues = [scipy.stats.binom_test(x=N_mutants_in_window, n=total_N_mutants, p=window_probability) 
                       for (N_mutants_in_window, window_probability) in zip(window_mutant_count_list, window_probabilities)]
        window_raw_pvalues[chrom] = numpy.array(raw_pvalues)
        # adjust the p-values for FDR, using the total N_samples (rather than just the number of samples in this chromosome)
        #  (if there are no windows on a given chromosome (it was too short for the size+offset), 
        #   just set the lists to empty and skip rather than running statistics on empty lists, which can give warnings)
        if len(window_mutant_count_list):
            adjusted_pvalues = statistics_utilities.FDR_adjust_pvalues(window_raw_pvalues[chrom], method='BH', N=N_samples)
        else:
            adjusted_pvalues = []
        window_adjusted_pvalues[chrom] = numpy.array(adjusted_pvalues)
    return window_adjusted_pvalues, window_raw_pvalues, window_if_coldspot, window_mutant_counts, window_expected_mutant_counts
    # TODO unit-test this? How?


def get_hot_cold_spot_data(hc_windowsizes_offsets, DATASET, chromosome_lengths, m20, m21, fraction_20bp=None, quiet=False):
    """ Convenience function to run find_hot_cold_spots for many windowsizes and offsets and collate the data. """
    average_mutants_per_window = {}
    hc_spot_qvals, hc_spot_pvals, hc_spot_ifcold, hc_spot_mcounts, hc_spot_expected_mcounts = {}, {}, {}, {}, {}
    for window_size,offset_list in sorted(hc_windowsizes_offsets.items()):
        average_mutants_per_window[window_size] = find_average_mutant_count(DATASET, window_size, chromosome_lengths)
        hc_spot_qvals[window_size], hc_spot_pvals[window_size], hc_spot_ifcold[window_size], hc_spot_mcounts[window_size], hc_spot_expected_mcounts[window_size] = {}, {}, {}, {}, {}
        for offset in sorted(offset_list):
            hc_data = find_hot_cold_spots(DATASET, window_size, m20, m21, offset, fraction_20bp, chromosome_lengths=chromosome_lengths, quiet=quiet)
            hc_spot_qvals[window_size][offset], hc_spot_pvals[window_size][offset], hc_spot_ifcold[window_size][offset], hc_spot_mcounts[window_size][offset], hc_spot_expected_mcounts[window_size][offset] = hc_data
    return average_mutants_per_window, hc_spot_qvals, hc_spot_pvals, hc_spot_ifcold, hc_spot_mcounts, hc_spot_expected_mcounts


def get_hot_cold_spot_list(pvalue_data_window_size_offset_dict, ifcold_data_window_size_offset_dict, 
                           mcount_data_window_size_offset_dict=None, expected_mcount_data_window_size_offset_dict=None,
                           average_mutants_per_window=None, pval_cutoff=0.05, print_result_counts=True, print_info=True, 
                           OUTFILE=sys.stdout):
    """ Return list of (chrom, start, end, pvalue, kind, window_offset, found_mutants, expected_mutants) tuples with pvalue<=cutoff.

    First three arguments should be window_size:(window_offset:(chromosome:list_of_window_values))) triple nested dictionaries.
      - the p-values should be 0-1, and you should probably get FDR correction done on them first; 
      - the ifcold values should be -1 for coldspots and 1 for hotspots. 
    Pval_cutoff should be a number between 0 and 1. 

    Mcount_data_window_size_offset_dict should contain the number of mutants found in each window, 
     and expected_mcount_data_window_size_offset_dict the expected number of mutants. 
     Both are optional - if None is given, mutant numbers will be shown as 'unknown'. 
    Average_mutants_per_window should be a window_size:average_mutants dict, optional - the value is just printed next to the real 
     mutant count for each hot/coldspot reported, to help see the effect size.  If None, average values are printed as 'unknown'

    In the output, kind is 'hotspot' or 'coldspot'. 
    Sort the output data by chromosome/position;  If print_info, print the output data.
    """
    defaultdict = collections.defaultdict
    if mcount_data_window_size_offset_dict is None:
        mcount_data_window_size_offset_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'unknown')))
    if expected_mcount_data_window_size_offset_dict is None:
        expected_mcount_data_window_size_offset_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'unknown')))
    if average_mutants_per_window is None:
        average_mutants_per_window = defaultdict(lambda: 'unknown')
    hc_data_list = []
    for window_size, pvalue_data_window_offset_dict in pvalue_data_window_size_offset_dict.items():
        for window_offset, pvalue_data in pvalue_data_window_offset_dict.items():
            N_spots_for_window_offset = 0
            ifcold_data = ifcold_data_window_size_offset_dict[window_size][window_offset]
            mcount_data = mcount_data_window_size_offset_dict[window_size][window_offset]
            expected_mcount_data = expected_mcount_data_window_size_offset_dict[window_size][window_offset]
            for chrom, pvalues in pvalue_data.items():
                ifcold_vals = ifcold_data[chrom]
                mcounts = mcount_data[chrom]
                expected_mcounts = expected_mcount_data[chrom]
                color_values = []
                for N,(pvalue,ifcold,mcount,expected) in enumerate(zip(pvalues,ifcold_vals,mcounts,expected_mcounts)):
                    if pvalue <= pval_cutoff:
                        N_spots_for_window_offset += 1
                        kind = 'coldspot' if ifcold else 'hotspot'
                        start_pos = N*window_size + window_offset
                        end_pos = (N+1)*window_size + window_offset
                        hc_data_list.append((chrom, start_pos, end_pos, pvalue, kind, window_offset, mcount, expected))
            if print_result_counts:
                OUTFILE.write(
                    "%s results below adjusted p-value %s for data with %s window and offset %s\n"%(N_spots_for_window_offset, 
                                                                  pval_cutoff, format_bp(window_size), format_bp(window_offset)))
    hc_data_list.sort(key=lambda x: (basic_seq_utilities.chromosome_sort_key(x[0]), x[1:]))
    if print_info:
        for N,(chrom, start_pos, end_pos, pvalue, kind, offset, mcount, expected) in enumerate(hc_data_list):
            window_size = end_pos-start_pos
            OUTFILE.write("%s %s-%s (window size %s) - %s, %.3g adj. p-value (%s mutants, expected %.3g, avg %.3g) (#%s)\n"%(
                                        chrom, format_bp(start_pos), format_bp(end_pos), format_bp(window_size), kind, 
                                        pvalue, mcount, expected, average_mutants_per_window[window_size], N))
    return hc_data_list
    # TODO should unit-test this!


def filter_hot_cold_spot_list(hc_significant_list, filter_by='effect_size', outfile=None):
    """ Filter overlapping sets of hot/coldspots.

    Because we're doing hot/coldspot analysis for multiple window sizes, we get multiple overlapping spots in many areas.
    Those need to be filtered out to the "best" one, or else a 10bp hotspot with a ton of insertions will also show up as 
        a 1Mb hotspot and all sizes in between, etc.
    The idea of this filtering is: if spot A completely contains spot B or vice versa, only keep one of them. 

    The two methods currently implemented (in the filter_by parameter) are:
        - pvalue - pick the one with the lower pvalue
        - effect_size - pick the one with the higher effect size, i.e. the higher observed:expected ratio (or lower for coldspots).

    The hc_significant_list is the output of mutant_simulations.get_hot_cold_spot_list.
    """
    if filter_by == 'pvalue':
        hc_tmp_list = [tuple(x) for x in hc_significant_list]
        hc_tmp_list.sort(key = lambda x: x[3])
    elif filter_by == 'effect_size':
        hc_tmp_list_hot = sorted([tuple(x) for x in hc_significant_list if x[4]=='hotspot'], key = lambda x: x[7]/x[6])
        hc_tmp_list_cold = sorted([tuple(x) for x in hc_significant_list if x[4]=='coldspot'], key = lambda x: x[6]/x[7])
        hc_tmp_list = hc_tmp_list_hot + hc_tmp_list_cold
    else:
        raise Exception("Unknown filtering method %s!"%filter_by)
    hc_filtered_list = []
    while hc_tmp_list:
        curr_best = hc_tmp_list[0]
        del hc_tmp_list[0]
        hc_filtered_list.append(curr_best)
        new_hc_tmp_list = []
        for x in hc_tmp_list:
            if x[0] == curr_best[0] and x[4] == curr_best[4]:
                if x[1] <= curr_best[1] and x[2] >= curr_best[2]:     continue
                if x[1] >= curr_best[1] and x[2] <= curr_best[2]:     continue
            new_hc_tmp_list.append(x)
        hc_tmp_list = new_hc_tmp_list
    print "%s spots filtered down to %s (%s hot, %s cold)"%(len(hc_significant_list), len(hc_filtered_list), 
        sum(1 for x in hc_filtered_list if x[4]=='hotspot'), sum(1 for x in hc_filtered_list if x[4]=='coldspot'))
    print "number by max pvalue: " + ", ".join("%s: %s"%(p, sum(1 for x in hc_filtered_list if x[3]<p)) 
        for p in (0.01, 1e-4, 1e-6, 1e-10, 1e-20, 1e-50, 1e-100, 1e-200, 1e-400))
    print "number by size: " + ", ".join("%s: %s"%(basic_seq_utilities.format_base_distance(s), n) 
        for (s, n) in sorted(collections.Counter(x[2]-x[1] for x in hc_filtered_list).items()))
    if outfile:
        format_bp = lambda x: basic_seq_utilities.format_base_distance(x, False)
        with open(outfile, 'w') as OUTFILE:
          for N,(chrom, start_pos, end_pos, pvalue, kind, offset, mcount, expected) in enumerate(hc_filtered_list):
            window_size = end_pos-start_pos
            OUTFILE.write("%s %s-%s (window size %s) - %s, %.3g adj. p-value (%s mutants, expected %.3g) (#%s)\n"%(
                                        chrom, format_bp(start_pos), format_bp(end_pos), format_bp(window_size), kind, 
                                        pvalue, mcount, expected, N))
    return hc_filtered_list


def observed_mutants_in_window(chromosome, start_pos, end_pos, dataset):
    """ Return the number of mutants in dataset that are on chromosome between start and end (inclusive). 

    Dataset can be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a chrom:position_list or (chrom,strand):position_list dictionary (and position_list may be a list or numpy array!).
    """
    # If dataset is an object, get a generator of min_pos
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        positions_in_chromosome = (m.position.min_position for m in dataset if m.position.chromosome==chromosome)
    # Otherwise assume dataset is a chrom:position_list or (chrom,strand):position_list dict, 
    #  and just grab the position list for that chromosome, or [] if the chromosome isn't in the keys
    elif isinstance(dataset.keys()[0], basestring):
        positions_in_chromosome = dataset.get(chromosome, [])
    else:
        # CAREFUL: need to convert to lists before taking a sum, normally the values may be numpy arrays, which add funny!!
        positions_in_chromosome = sum((list(dataset[key]) for key in dataset.keys() if key[0]==chromosome), [])
    # get the number of positions that are between the start/end given
    return sum(int(start_pos <= pos <= end_pos) for pos in positions_in_chromosome)
    # TODO unit-test


def expected_mutants_in_window(chromosome, start_pos, end_pos, mappable_positions_20bp, mappable_positions_21bp, 
                               fraction_20bp, N_mutants, all_chromosomes=None):
    """ Return the expected number of mutants in chromosome between start and end (inclusive), given genome mappability etc.

    For descriptions of mappable_positions_*, fraction_20bp arguments, see help for find_hot_cold_spots function.
    N_mutants should be the total number of mutants in the dataset (we're calculating how many of those are expected in the window);
    all_chromosomes should be a list of chromosomes in which we're looking at mutants; if None, all mappable chromosomes are allowed.

    The probability of a mutant falling in the given window is defined by the number of mappable positions in that window
     vs the whole genome, separately for 20bp and 21bp flanking regions (mappable_positions_* arguments).
    Given the total number of mutants (N_mutants), and how many of them are 20 vs 21bp (fraction_20bp),
     the expected number of mutants in the given window can be calculated (again, separately for 20 and 21bp, then added together).
    """
    mappable_pos_in_window_20bp = observed_mutants_in_window(chromosome, start_pos, end_pos, mappable_positions_20bp)
    mappable_pos_in_window_21bp = observed_mutants_in_window(chromosome, start_pos, end_pos, mappable_positions_21bp)
    mappable_pos_total_20bp = sum(len(chrom_mapp_positions) for chrom,chrom_mapp_positions in mappable_positions_20bp.items() 
                                  if (all_chromosomes is None or _get_chrom_either_strand(chrom) in all_chromosomes))
    mappable_pos_total_21bp = sum(len(chrom_mapp_positions) for chrom,chrom_mapp_positions in mappable_positions_21bp.items() 
                                  if (all_chromosomes is None or _get_chrom_either_strand(chrom) in all_chromosomes))
    total_N_20bp = N_mutants * fraction_20bp
    total_N_21bp = N_mutants - total_N_20bp
    expected_N_20bp = total_N_20bp * (mappable_pos_in_window_20bp/mappable_pos_total_20bp)
    expected_N_21bp = total_N_21bp * (mappable_pos_in_window_21bp/mappable_pos_total_21bp)
    return expected_N_20bp+expected_N_21bp
    # TODO unit-test!  I had some worries about this in the past, but I don't remember what they were.


################################# Simulating random size-N dataset, various options #######################################

###### help functions


### Two weighted_random_choice_* fuctions instead of one with an option, so that they're both as fast as possible!
#  see http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/ - using weighted_choice_sub
#  MAYBE-TODO make this faster by using the WeightedRandomGenerator class from the same site?
#  also see http://docs.python.org/3/library/random.html - similar, but for python3 only, (2.7 is missing itertools.accumulate)

def weighted_random_choice_single(values, weights, sum_of_weights=None):
    """ Given a list of values and weights, pick a random value, with probability=weight for each value. """
    if sum_of_weights is None:
        sum_of_weights = sum(weights)
    rnd = random.random() * sum_of_weights 
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return values[i]
    # LATER-TODO move this to general_utilities?

def weighted_random_choice_multi(values, weights, N=10):
    """ Given a list of values and weights, pick N random values, with probability=weight for each value.  """
    sum_of_weights = sum(weights)
    results = []
    for _ in range(N):
        rnd = random.random() * sum_of_weights 
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                results.append(values[i])
                continue
    return results
    # LATER-TODO move this to general_utilities?


###### dataset simulation functions (different methods)

def simulate_dataset_from_mappability(N_mutants, fraction_20bp, mappable_positions_20bp, mappable_positions_21bp, 
                                      all_chromosomes=None, include_strand=False, fraction_plus_strand=0.5):
    """ Return position info from simulating inserting N mutants into the genome randomly, taking mappability into account. 
    
    Return a (chrom,strand):position_list dictionary if include_strand, or else just a chrom:position_list one...
     (the latter is useful for plotting or gap size analysis, since both of those ignore strand.

    Fraction_20bp should be the fraction of mutants that have 20bp vs 21bp flanking regions 
     (since the mappability data is slightly different for those two cases).
    Fraction_plus_strand is the fraction of mutants on the +strand (vs the -strand); 
     note that this is STILL RELEVANT even if include_strand is False, because mappability differs slightly between strands!

    All_chromosomes should give a list of chromosomes to include (e.g. if trying to simulate a dataset without chloroplast/mito
     genomes, don't put those on a list!); if None, all chromosomes in mappable_positions_20bp will be used.

    The two mappable_positions_* arguments should be (chrom,strand):position_list dictionaries given all the mappable positions,
     as generated by genome_mappable_insertion_sites with include_strand==True.
    """
    # MAYBE-TODO add a way of copying N_mutants and fraction_* from an existing dataset instead of putting them in separately?  Could do that by putting it here as an argument, or with a help function or something
    if all_chromosomes is None:
        all_chromosomes = set(chrom for chrom,s in mappable_positions_20bp) | set(chrom for chrom,s in mappable_positions_21bp)
    mappable_position_data = {20: mappable_positions_20bp, 21: mappable_positions_21bp}
    # chromosome mappable lengths - a list of (chrom,mappable_lengths) tuples for each (flank_length,chrom,strand) combination, 
    #  to use as a list of values and weights when randomly choosing a chromosome - only include ones in all_chromosomes!
    chrom_mappable_len = collections.defaultdict(list)
    for flank_len,mappable_positions in mappable_position_data.iteritems():
        for (chrom,strand),pos_list in mappable_positions.iteritems():
            if chrom in all_chromosomes:
                chrom_mappable_len[flank_len,strand].append((chrom, len(pos_list)))
    simulated_positions = collections.defaultdict(list)
    for _ in range(N_mutants):
        # first choose a length (20 or 21bp) based on Fraction_20bp, and a strand based on fraction_plus_strand.
        #   random.random() gives a value x in the interval [0, 1) (including 0, excluding 1).
        #   if fraction_* is 0, we want to never pick that option, so the right test is random<fraction.
        if random.random() < fraction_20bp: flank_length = 20
        else:                               flank_length = 21
        if random.random() < fraction_plus_strand: strand = '+'
        else:                                      strand = '-'
        # next choose a chromosome, with a probability corresponding to the mappable length of each chromosome
        #  (given the flank_length and the strand)
        chrom = weighted_random_choice_single(*zip(*chrom_mappable_len[(flank_length,strand)]))
        # next choose a position from the list of mappable positions on the chromosome
        pos = random.choice(mappable_position_data[flank_length][(chrom,strand)])
        # save the data
        if include_strand:  simulated_positions[(chrom,strand)].append(pos)
        else:               simulated_positions[chrom].append(pos)
    # convert output to numpy arrays for speed and memory efficiency
    for key,val_list in simulated_positions.iteritems():
        simulated_positions[key] = numpy.array(val_list)
    return simulated_positions
    # TODO how would I test this sensibly??  Complicated... But I did run it, visualize the results, and compare to a real dataset, and it looked all right.  (../1211_positions_Ru-screen1-for-paper/positions_over_chromosomes/*simulated*)


def simulate_genelist_from_mappability(N_mutants, gene_mappable_lengths, intergenic_mappable_length):
    """ Return gene list from simulating inserting N mutants into the genome randomly, based on gene/intergenic mappable lengths.
    
    Gene_mappable_lengths should be a gene:mappable_length dict; intergenic_mappable_length should be a number.
    They will be used as weights in a weighted random choice for each new mutant insertion (so the probability of each mutant
     ending up in a gene or intergenic region will be proportional to the gene/intergenic mappable lengths).

    Output will be just a list of the gene names for each simulated mutant.
     output will include duplicates, i.e. if there are two mutants in a gene, the gene will show up twice on the list; 
     for intergenic mutants, mutant_analysis_classes.SPECIAL_GENE_CODES.not_found will be used; 
     output will be in no particular order;  output can be used as input to gene_counts_for_mutant_subsets.
    Note: the output should be the same as from running simulate_dataset_from_mappability and then find_genes_for_simulated, 
     but should be a lot faster.
    """
    # make gene ID and length lists to use as weighted_random_choice_multi arguments; also pre-calculate sum_of_lengths
    #  (it'll be faster to do this rather than redo the zip/sum/etc for each N)
    #  (and I think prepending instead of appending the intergenic regions will make it faster, since more cases will stop early)
    gene_IDs, mapp_lengths = zip(*gene_mappable_lengths.items())
    gene_IDs = (mutant_analysis_classes.SPECIAL_GENE_CODES.not_found,) + gene_IDs
    mapp_lengths = (intergenic_mappable_length,) + mapp_lengths
    # use simple weighed random choice to choose a gene N times
    simulated_genelist = weighted_random_choice_multi(gene_IDs, mapp_lengths, N_mutants)
    return simulated_genelist
    # TODO unit-test!
    # TODO this gives different results than simulate_dataset_from_mappability - search for it in 1211_positions_Ru-screen1-for-paper/notes.txt for more info. I ended up not using it - if I want to use it, I should figure out what's wrong and fix it!          Is it because I'm treating 20 and 21bp flanking regions equally when really I shouldn't be?


### simulate dataset with N randomly positioned mutants, matching the gap-size distribution of the real dataset, and taking into account mappable/unmappable positions
# MAYBE-TODO - is this even a useful idea?  It won't match the real hot/coldspot locations...

# MAYBE-TODO would it be possible to also match the real hotspots and coldspots of the real dataset?


###### getting more information for simulated datasets

def dataset_object_from_simulated(simulated_dataset, get_gene_info=True, genefile=mutant_utilities.DEFAULT_GENE_POS_FILE):
    """ Convert a simulated dataset position dictionary into an Insertional_mutant_pool_dataset object, with optional gene info.
    
    Output will be a mutant_analysis_classes.Insertional_mutant_pool_dataset object, with all readcounts set to 1.

    If get_gene_info is True, get the gene information for the positions (from genefile if given, otherwise the default).

    Input should be a (chrom,strand):position_list dictionary, like from simulate_dataset_from_mappability with include_strand True.

    Warning: can take a LOT of time/memory if the simulated dataset is large!  May be best to avoid if possible.
    """
    new_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset()
    # Go over all the simulated read positions, add them to dataset
    for (chrom,strand),pos_list in simulated_dataset.items():
        for pos in pos_list:
            read_count = 1
            # make a simple position - just always specify position_before, regardless of strand
            position = mutant_analysis_classes.Insertion_position(chrom, strand, position_before=pos, immutable=True)
            # grab the right mutant based on the position, and add the reads to it; 
            curr_mutant = new_dataset.get_mutant(position)
            # just add 1 readcount, and don't bother with the sequences etc
            curr_mutant.add_counts(1,0,1)
    # Add gene info to dataset if desired
    if get_gene_info:
        new_dataset.find_genes_for_mutants(genefile, detailed_features=True)
    return new_dataset


def find_genes_for_simulated(sim_dataset, genefile=mutant_utilities.DEFAULT_GENE_POS_FILE):
    """ Get a list of gene names for all the mutants in the dataset (unordered, with duplicates)

    Input should be a (chrom,strand):position_list dictionary, like from simulate_dataset_from_mappability.

    Output will be just a list of the gene names for all the positions in the input, EXCEPT chromosomes not present in genefile;
     output will include duplicates, i.e. if there are two mutants in a gene, the gene will show up twice on the list; 
     when a position isn't in a gene, one of the names in mutant_analysis_classes.SPECIAL_GENE_CODES will be used; 
     output will be in no particular order;  output can be used as input to gene_counts_for_mutant_subsets.

    This function is based on mutant_analysis_classes.Insertional_mutant_pool_dataset.find_genes_for_mutants,
     but acts on simple position data instead of mutant objects; it uses mutant_analysis_classes.find_gene_by_pos_gff3.
    """ 
    gene_list = []
    # only look at genes, not sub-features
    genefile_parsing_limits = {'gff_type': ['gene']}
    with open(os.path.expanduser(genefile)) as GENEFILE:
        # go over all chromosomes in genefile
        #  (just ignore positions in chromosomes that weren't listed in the genefile - MAYBE-TODO do something else about those?)
        for chromosome_record in gff_examine_file.GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
            chromosome = chromosome_record.id
            # for each position in sim_dataset on that chromosome, get the gene info based on the position
            #   (position input to mutant_analysis_classes.find_gene_by_pos_gff3 is a (strand, pos_before, pos_after tuple), 
            #  and add gene name to output list (ignore the other two outputs - orientation and feature).
            for strand in '+-':
                for position in sim_dataset[(chromosome, strand)]:
                    gene, _, _ = mutant_analysis_classes.find_gene_by_pos_gff3((strand, position, position+1), chromosome_record, 
                                                                     detailed_features=False, quiet=True)
                    gene_list.append(gene)
    return gene_list


################################# Randomly chosen subsets of real dataset #######################################

### number of genes with 1+/2+/etc mutants vs number of mutants (randomly chosen mutant subsets)
   
def gene_counts_for_mutant_subsets(dataset, subset_sizes, max_N_mutants=3):
    """ Return numbers of genes with N mutants for different-sized random subsets of dataset, or a single subset size.

    Return an N:list_of_gene_numbers dictionary, where N is each value between 1 and max_N_mutants, 
     and list_of_gene_numbers contains the number of genes with at least N mutants 
      in randomly chosen subsets of dataset, of sizes given by subset_sizes (a sorted list).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
      or a list of mutants (mutant_analysis_classes.Insertional_mutant instances), 
     OR just a list of genes for each mutant, since that's the only information that's actually used.
    """
    # if dataset is a list of mutant objects or a full mutant dataset object, extract just the gene list
    try:                        genes = [mutant.gene for mutant in dataset]
    except AttributeError:      genes = dataset
    subset_sizes.sort()
    subset_sizes = [x for x in subset_sizes if x <= len(dataset)]
    # only shuffle the dataset if we're not going to be using the full one anyway
    if not subset_sizes == [len(dataset)]:
        random.shuffle(genes)
    mutantN_to_gene_count = collections.defaultdict(list)
    for subset_size in subset_sizes:
        # get the gene counts for each mutant number (ignoring not-a-gene markers) 
        #  (dict.pop with a second argument doesn't give an error if that key is absent)
        gene_to_mutantN = collections.Counter(genes[:subset_size])
        for not_a_real_gene in mutant_analysis_classes.SPECIAL_GENE_CODES.all_codes:
            gene_to_mutantN.pop(not_a_real_gene, None)
        # now append the number of genes with N+ mutants to the appropriate mutantN_to_gene_count values
        for N_mutants in range(1,max_N_mutants+1):
            mutantN_to_gene_count[N_mutants].append(len([1 for count in gene_to_mutantN.values() if count>=N_mutants]))
    return mutantN_to_gene_count


################################################# Testing etc ##################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__genome_mappable_everything(self):
        # testing three related functions in parallel, since they use the same cases
        #  - arguments to genome_mappable_slicesare (slice_len, genome_seq, print_info=True)
        #  - arguments to genome_mappable_insertion_sites are (flanking_region_length=21, mappable_slice_pos_dict=None, 
        #           genome_seq=None, include_strand=True, end_sequenced='5prime', print_info=True):

        ### only testing whether there are 0 unique sequences, or a non-zero number
        def _test_all_empty(slice_len, genome):
            """ Test whether all the genome_mappable_* functions return empty sets/lists/dicts/whatever. """
            def _test_empty(D): return sum([len(l) for l in D.values()])==0
            outcomes = set()
            # test genome_mappable_slices
            slice_data = genome_mappable_slices(slice_len, genome_seq=genome, print_info=False)
            outcomes.add(_test_empty(slice_data))
            # test genome_mappable_insertion_sites and genome_mappable_insertion_sites_multi in all variants
            for include_strand in (True,False):
                for end in basic_seq_utilities.SEQ_ENDS:
                    args = (include_strand,end, False)
                    outcomes.add(_test_empty(genome_mappable_insertion_sites(slice_len, slice_data, None, *args)))
                    outcomes.add(_test_empty(genome_mappable_insertion_sites(slice_len, None, genome, *args)))
                    for N in (1,2,5):
                        outcomes.add(_test_empty(genome_mappable_insertion_sites_multi([slice_len]*N, [slice_data]*N, None, *args)))
                        outcomes.add(_test_empty(genome_mappable_insertion_sites_multi([slice_len]*N, None, genome, *args)))
            if len(outcomes) != 1:
                raise Exception("Inconsistent outcomes in _test_all_empty in test__genome_mappable_everything!")
            return outcomes.pop()
        # no unique sequences regardless of slice_len: 
        #  empty genome, two identical or reverse-complement chromosomes, one palindromic chromosome, 
        for slice_len in (1,2,3,4):
            assert _test_all_empty(slice_len, {})
            assert _test_all_empty(slice_len, {'a':'AAA', 'b':'AAA'})
            assert _test_all_empty(slice_len, {'a':'AAA', 'b':'TTT'})
            assert _test_all_empty(slice_len, {'a':'AATT'})
        # one chromosome with 3+ tandem repeats - no unique sequences as long as slice_len <= repeat_len
        #  (this doesn't apply for 2 tandem repeats, because in ACTACT, CTA is still unique! But in ACTACTACT, it's present twice.)
        for repeat_seq in ('ATC', 'AATCCG', 'ACATGACGAGACGGG'):
            for N_repeats in (3,4,5,10):
                chrom_seq = repeat_seq * N_repeats
                for slice_len in range(1, len(repeat_seq)+2):
                    assert _test_all_empty(slice_len, {'a':chrom_seq})
        # no unique sequences if slice_len<chromosome_len:
        #  one chromosome that has repeating substrings or all matching substrings with own reverse-complement
        for slice_len in (1,2):
            assert _test_all_empty(slice_len, {'a':'AAA'})
            assert _test_all_empty(slice_len, {'a':'ATA'})
        # no unique sequences ONLY IF slice_len==1:  any case with more than one AT or GC base.
        for slice_len,if_empty in ((1,True), (2,False), (3,False), (10,False)):
            assert _test_all_empty(slice_len, {'a':'AA'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'ATT'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'GCC'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'G', 'b':'GC'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'AG', 'b':'TC'}) == if_empty
        # cases that should have non-zero unique sequences
        assert not _test_all_empty(3, {'a':'AAA'})
        assert not _test_all_empty(3, {'a':'AAA', 'b':'GGG'})
        assert not _test_all_empty(1, {'a':'AGA'})
        assert not _test_all_empty(2, {'a':'AGA'})
        assert not _test_all_empty(2, {'a':'AAT'})

        ### more detailed testing of actual non-empty outputs
        # help functions for testing all the functions in parallel
        def _read_raw_data(raw_data):
            """ get a dict of lists of ints or (int.str) tuples from a raw string. """
            formatted_data = collections.defaultdict(list)
            for chr_data in raw_data.split(', '):
                chrom, pos_data = chr_data.strip().split(': ')
                # pos_data can be ints (1 25 301) or ints with strand info (1- 25+ 301+)
                for x in pos_data.split(' '):
                    if x[-1] in '+-':   formatted_data[chrom,x[-1]].append(int(x[:-1]))
                    else:               formatted_data[chrom].append(int(x))
            for val in formatted_data.values():
                val.sort()
            return dict(formatted_data)
        def _compare_dicts(numpy_array_dict, list_dict):
            """ Compare a key:val_numpy_array dict to a key:val_list dict, make sure they match. """
            assert numpy_array_dict.keys() == list_dict.keys()
            for key,val in numpy_array_dict.iteritems():
                assert list(val) == list_dict[key]
        def _test_all(slice_len, genome, raw_slice_data, raw_pos_data_5prime):
            """ Test genome_mappable_slices and all variants of genome_mappable_insertion_sites* against expected output. """
            # get the expected output data from the simplified string formats
            slice_data = _read_raw_data(raw_slice_data)
            pos_data_5prime = _read_raw_data(raw_pos_data_5prime)
            # from pos_data_5prime, make pos_data_3prime (just switch all strands) and pos_data_no_strand (remove strand info)
            pos_data_3prime, pos_data_no_strand = {}, collections.defaultdict(list)
            for (chrom,strand),pos_list in pos_data_5prime.items():
                pos_data_3prime[chrom, '+' if strand=='-' else '-'] = pos_list
                pos_data_no_strand[chrom] += pos_list
                pos_data_no_strand[chrom].sort()
            # check genome_mappable_slices output (and save it for later)
            new_slice_data = genome_mappable_slices(slice_len, genome, False)
            _compare_dicts(new_slice_data, slice_data)
            # now try running genome_mappable_insertion_sites with both the raw slice_len/genome data, and the new_slice_data; 
            for include_strand,end,pos_data in ((True,'5prime',pos_data_5prime), (True,'3prime',pos_data_3prime), 
                                                (False,'5prime',pos_data_no_strand), (False,'3prime',pos_data_no_strand)):
                    args = (include_strand,end, False)
                    _compare_dicts(genome_mappable_insertion_sites(slice_len, slice_data, None, *args), pos_data)
                    _compare_dicts(genome_mappable_insertion_sites(slice_len, None, genome, *args), pos_data)
                    _compare_dicts(genome_mappable_insertion_sites_multi([slice_len], None, genome, *args), pos_data)
                    for extra in ([], [{}], [{}, {}, {}]):
                        _compare_dicts(genome_mappable_insertion_sites_multi([slice_len]*(len(extra)+1), [slice_data]+extra, 
                                                                                    None, *args), pos_data)
        # how genome_mappable_insertion_sites works if end is 5prime: 
        #   a mappable 1-2 flanking region means a +strand insertion at 2-? and a -strand insertion at ?-1 will be mappable;
        #  if end is 3prime, it's the same positions but opposite strands.

        # a single-base-repeat chromosome with slice size equal to chromosome size has all unique sequences;
        #  same for two such chromosomes that aren't reverse-complement
        _test_all(3, {'a':'AAA'},                  'a: 1',        'a: 0- 3+')
        _test_all(3, {'a':'AAA', 'b':'GGG'},       'a: 1, b: 1',  'a: 0- 3+, b: 0- 3+')
        # same test, but with the genome read from a fasta file
        _test_all(3, 'test_data/INPUT_genome0.fa', 'a: 1, b: 1',  'a: 0- 3+, b: 0- 3+')

        # a single chromosome that's not a palindrome or a base-repeat has all unique sequences
        #  (or just some if slice_len is 1 and some bases show up twice)
        # (note - the first case here has two different-strand mappable insertions in one position, 
        #  so the genome_mappable_insertion_sites has one position repeated twice - important to check that!)
        _test_all(1, {'a':'AG'},   'a: 1 2',  'a: 0- 1+ 1- 2+') 
        _test_all(2, {'a':'AAT'},  'a: 1',    'a: 0- 2+')
        curr_genome = {'a':'AGA'}
        _test_all(1, curr_genome,  'a: 2',    'a: 1- 2+')
        _test_all(2, curr_genome,  'a: 1 2',  'a: 0- 2+ 1- 3+')

        # test genome_mappable_insertion_sites_multi on the two curr_genome cases added together, different flank sizes
        #  (getting the data directly from genome, or from two genome_mappable_slices results)
        #  (only testing the no-strand version for simplicity - MAYBE-TODO test all versions?)
        _compare_dicts(genome_mappable_insertion_sites_multi([1,2], None, curr_genome, False, '5prime', False), {'a':[0,1,1,2,2,3]})
        slices_1 = genome_mappable_slices(1, curr_genome, False)
        slices_2 = genome_mappable_slices(2, curr_genome, False)
        for fl_both,slices_both in (([1,2],[slices_1,slices_2]), ([2,1],[slices_2,slices_1])):
            _compare_dicts(genome_mappable_insertion_sites_multi(fl_both, slices_both, None, False, '5prime', False), 
                           {'a':[0,1,1,2,2,3]})

 #  def test__gene_mappability(self):
 #      test_gene_infile = 'test_data/INPUT_gene-data-2_simple.gff3'
 #      ### test 1: no mappable positions - zero mappability
 #      for exclude_UTRs_in_bins in (False, True):
 #          for no_mappability_val in [{}, {'chrA':[]}, {'chrB':[]}, {('chrA','+'):[], ('chrA','-'):[]}]:
 #              for test_mappability_data in [{}, {10:no_mappability_val}, {20: no_mappability_val, 21: no_mappability_val}]:
 #                  gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(test_mappability_data, 
 #                                                                                        test_gene_infile, 2, exclude_UTRs_in_bins)
 #                  assert gene_mapp == {'gene1_plus':0, 'gene2_minus':0}
 #                  assert gene_mapp_2bins == {'gene1_plus':[0,0], 'gene2_minus':[0,0]}
 #                  assert feature_mapp == {'CDS':0, 'five_prime_UTR':0, 'three_prime_UTR':0, 
 #                                          'intron':0, 'gene':0, 'intergenic':0, 'all':0, 'MULTIPLE_SPLICE_VARIANTS':0}
 #      ### test 2: intergenic mappable positions only (for simplicity); testing mappability multiplier and strandedness.
 #      #  If we have mappability on a DIFFERENT chromosome that has no genes on it, 
 #      #   we get non-zero intergenic/all feature mappability values.
 #      #  Tests that the mappability multiplier is correct: with one flanking region length and 4 mappable positions, 
 #      #   the mappability is 2; with two lengths with 4 positions each, still 2; with one length with none and one with 4, 1.
 #      #   Tests whether it works for both stranded and non-stranded data, too.
 #      for (mapp_count, test_mappability_data ) in [(2,   {10: {'chrB':[1,2,3,4]} }), 
 #                                                   (2,   {10: {'chrB':[1,1,2,2]} }), 
 #                                                   (2,   {10: {('chrB','+'):[1,2], ('chrB','-'):[1,2]} }), 
 #                                                   (1,   {10: {('chrB','+'):[1,2], ('chrB','-'):[]} }), 
 #                                                   (2,   {10: {'chrB':[1,2,3,4]},   11: {'chrB':[1,2,3,4]} }), 
 #                                                   (1,   {10: {'chrB':[1,2,3,4]},   11: {'chrB':[]} }),
 #                                                   (1,   {10: {'chrB':[1,2]},   11: {'chrB':[1,2]} }),
 #                                                   (1,   {10: {'chrB':[1,2]},   11: {'chrB':[100,200]} }),
 #                                                   (0.5, {10: {'chrB':[1,2]},       11: {'chrB':[]} })]:
 #          for exclude_UTRs_in_bins in (False, True):
 #              gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(test_mappability_data, test_gene_infile, 
 #                                                                                        2, exclude_UTRs_in_bins)
 #              assert (gene_mapp, gene_mapp_2bins) == ({'gene1_plus':0, 'gene2_minus':0}, {'gene1_plus':[0,0], 'gene2_minus':[0,0]})
 #              assert feature_mapp == {'CDS':0, 'five_prime_UTR':0, 'three_prime_UTR':0, 'intron':0, 'gene':0, 
 #                                      'intergenic':mapp_count, 'all':mapp_count, 'MULTIPLE_SPLICE_VARIANTS':0}
 #      ### test 3: some realistic mappable positions (some gene, some intergenic; remember gene2 is -strand!)
 #      # help function to auto-make mappability_data: assume flank-length and chromosome; use every position twice, once per strand
 #      _make_mapp_data = lambda mappability_val: {10: {'chrA': sorted(mappability_val*2)}}
 #      ## test 3a: mappable positions in the middles of various features and intergenic
 #      gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(_make_mapp_data(
 #                                                                                      [150,250,251, 1000, 1150,1250, 2000]), 
 #                                                                              test_gene_infile, 2, exclude_UTRs_in_bins=False)
 #      assert (gene_mapp,gene_mapp_2bins) == ({'gene1_plus':3, 'gene2_minus':2}, {'gene1_plus':[3,0], 'gene2_minus':[0,2]})
 #      assert feature_mapp == {'CDS':3, 'five_prime_UTR':1, 'three_prime_UTR':0, 'intron':1, 'gene':5, 'intergenic':2, 'all':7, 
 #                              'MULTIPLE_SPLICE_VARIANTS':0}
 #      # version with exclude_UTRs_in_bins=True - only gene_mapp_2bins is changed, since one gene1 position was UTR
 #      gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(_make_mapp_data(
 #                                                                                          [150,250,251, 1000, 1150,1250, 2000]), 
 #                                                                              test_gene_infile, 2, exclude_UTRs_in_bins=True)
 #      assert (gene_mapp,gene_mapp_2bins) == ({'gene1_plus':3, 'gene2_minus':2}, {'gene1_plus':[2,0], 'gene2_minus':[0,2]})
 #      assert feature_mapp == {'CDS':3, 'five_prime_UTR':1, 'three_prime_UTR':0, 'intron':1, 'gene':5, 'intergenic':2, 'all':7, 
 #                              'MULTIPLE_SPLICE_VARIANTS':0}
 #      ## test 3b: checking edges - mappable positions right at the edge of intron/feature/intergenic
 #      gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(_make_mapp_data([400,401, 600,601, 700,701, 
 #                                                                                   1100,1101, 1200,1201]), 
 #                                                                  test_gene_infile, 2, exclude_UTRs_in_bins=False)
 #      assert (gene_mapp,gene_mapp_2bins) == ({'gene1_plus':5, 'gene2_minus':3}, {'gene1_plus':[1,4], 'gene2_minus':[0,3]})
 #      assert feature_mapp == {'CDS':4, 'five_prime_UTR':0, 'three_prime_UTR':2, 'intron':2, 'gene':8, 'intergenic':2, 'all':10, 
 #                              'MULTIPLE_SPLICE_VARIANTS':0}
 #      # version with exclude_UTRs_in_bins=True - only gene_mapp_2bins is changed, since two gene1 positions were UTR
 #      gene_mapp, gene_mapp_2bins, feature_mapp = gene_mappability_bins_features(_make_mapp_data([400,401, 600,601, 700,701, 
 #                                                                                   1100,1101, 1200,1201]), 
 #                                                                  test_gene_infile, 2, exclude_UTRs_in_bins=True)
 #      assert (gene_mapp,gene_mapp_2bins) == ({'gene1_plus':5, 'gene2_minus':3}, {'gene1_plus':[1,2], 'gene2_minus':[0,3]})
 #      assert feature_mapp == {'CDS':4, 'five_prime_UTR':0, 'three_prime_UTR':2, 'intron':2, 'gene':8, 'intergenic':2, 'all':10, 
 #                              'MULTIPLE_SPLICE_VARIANTS':0}
 #      # MAYBE-TODO add separate tests of the help functions?
 #      # MAYBE-TODO add tests for error/weird cases? overlapping features, 0-length stuff, etc...
 #      # TODO add multiple splice variant cases to unit-tests!
        
    def test__gene_mappability(self):
        test_gene_infile = 'test_data/INPUT_gene-data-2_simple.gff3'
        ### test 1: no mappable positions - zero mappability
        for no_mappability_val in [{'chrA':[], 'chrB':[]}, {('chrA','+'):[], ('chrA','-'):[], 'chrB':[]}]:
            for exclude_UTRs in (False, True):
                for test_mappability_data in [{10:no_mappability_val}, {20: no_mappability_val, 21: no_mappability_val}]:
                    self.assertEqual(gene_mappability(test_mappability_data, exclude_UTRs, test_gene_infile), 
                                     {'gene1_plus':0.0, 'gene2_minus':0.0, 'gene3_multi':0.0})
        ### test 2: testing mappability multiplier and strandedness.
        #  Tests that the mappability multiplier is correct: with one flanking region length and 4 mappable positions, 
        #   the mappability is 2; with two lengths with 4 positions each, still 2; with one length with none and one with 4, 1.
        #   Tests whether it works for both stranded and non-stranded data, too.
        for (mapp_count, test_mappability_data ) in [(2,   {10: {'chrA':[201,201,202,202]} }), 
                                                     (2,   {10: {'chrA':[201,202,203,204]} }), 
                                                     (4,   {10: {'chrA':[201,201,202,202,203,203,204,204]} }), 
                                                     (2,   {10: {('chrA','+'):[201,202], ('chrA','-'):[201,202]} }), 
                                                     (1,   {10: {('chrA','+'):[201,202], ('chrA','-'):[]} }), 
                                                     (2,   {10: {'chrA':[201,201,202,202]},     11: {'chrA':[201,201,202,202]} }), 
                                                     (2,   {10: {'chrA':[201,202,203,204]},     11: {'chrA':[201,202,203,204]} }), 
                                                     (1,   {10: {'chrA':[201,202,203,204]},     11: {'chrA':[]} }),
                                                     (1,   {10: {'chrA':[201,202]},             11: {'chrA':[201,202]} }),
                                                     (1,   {10: {'chrA':[201,202]},             11: {'chrA':[203,204]} }),
                                                     (0.5, {10: {'chrA':[201,202]},             11: {'chrA':[]} })
                                                    ]:
            for exclude_UTRs in (False, True):
                #print "test data", test_mappability_data, "expected count", mapp_count
                self.assertEqual(gene_mappability(test_mappability_data, exclude_UTRs, test_gene_infile), 
                                 {'gene1_plus': float(mapp_count), 'gene2_minus':0.0, 'gene3_multi':0.0})
        ### test 3: some realistic mappable positions (some gene, some intergenic; remember gene2 is -strand!)
        # help function to auto-make mappability_data: assume flank-length and chromosome; use every position twice, once per strand
        _make_mapp_data = lambda mappability_val: {10: {'chrA': sorted(mappability_val*2)}}
        ## test 3a: mappable positions in the middles of various features and intergenic
        gene_mapp = gene_mappability(_make_mapp_data([150,250,251, 1000, 1150,1250, 2000]), 
                                        exclude_UTRs=False, genefile=test_gene_infile)
        self.assertEqual(gene_mapp, {'gene1_plus':3.0, 'gene2_minus':2.0, 'gene3_multi': 0.0})
        # version with exclude_UTRs=True - one gene1 position was UTR
        gene_mapp = gene_mappability(_make_mapp_data([150,250,251, 1000, 1150,1250, 2000]), 
                                        exclude_UTRs=True, genefile=test_gene_infile)
        self.assertEqual(gene_mapp, {'gene1_plus':2.0, 'gene2_minus':2.0, 'gene3_multi': 0.0})
        ## test 3b: checking edges - mappable positions right at the edge of intron/feature/intergenic
        gene_mapp = gene_mappability(_make_mapp_data([400,401, 600,601, 700,701, 1100,1101, 1200,1201]), 
                                        exclude_UTRs=False, genefile=test_gene_infile)
        self.assertEqual(gene_mapp, {'gene1_plus':5.0, 'gene2_minus':3.0, 'gene3_multi': 0.0})
        # version with exclude_UTRs=True - two gene1 positions were UTR
        gene_mapp = gene_mappability(_make_mapp_data([400,401, 600,601, 700,701, 1100,1101, 1200,1201]), 
                                        exclude_UTRs=True, genefile=test_gene_infile)
        self.assertEqual(gene_mapp, {'gene1_plus':3.0, 'gene2_minus':3.0, 'gene3_multi': 0.0})
        # MAYBE-TODO add separate tests of the help functions?
        # MAYBE-TODO add tests for error/weird cases? overlapping features, 0-length stuff, etc...
        # TODO add multiple splice variant cases to unit-tests!
        
    # LATER-TODO add unit-tests for other stuff!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
