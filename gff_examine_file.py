#!/usr/bin/env python2.7
""" Look at a gff file, print various overview info (to stdout, or outfile if given), check that everything looks sane.

See option help messages (run with -h) for details.. 

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: gff_examine_file.py [options] gff_infile [outfile]"""

# basic libraries
import sys, os, time
import unittest
import itertools
from collections import defaultdict, Counter
import pprint
# other libraries
# MAYBE-TODO try https://github.com/daler/gffutils instead of bcbio-gff?
from BCBio import GFF
from Bio import SeqIO
# my modules
from general_utilities import count_list_values, split_into_N_sets_by_counts, value_and_percentages, add_to_dict_no_replacement
from testing_utilities import run_functional_tests
from mutant_utilities import DEFAULT_NUCLEAR_GENOME_FILE

# At once point BCBio.GFF stopped working after upgrading my computer from Xubuntu 12.10 to 14.04, see https://github.com/chapmanb/bcbb/issues/89 for report/fix.

### Constants

GFF_strands = {1:'+', -1:'-'}

### Help functions used in main - the inputs are GFF objects from an already parsed file

class NoRNAError(Exception):        pass
class MultipleRNAError(Exception):  pass

def check_for_overlapping_genes(sequence_record):
    """ Given a GFF sequence record (from BCBio.GFF parser), return list of tuples of IDs of overlapping genes.  
    Not perfect: may miss some of the overlaps if some genes are completely contained within other genes (raises warning), 
     but if it returns no overlaps, there definitely aren't any."""
    overlapping_gene_pairs = []
    all_gene_positions = []
    for gene in sequence_record.features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_gene_positions.append((gene.location.start.position+1, gene.location.end.position, gene.id))
    all_gene_positions.sort()
    for gene1_data,gene2_data in itertools.izip(all_gene_positions,all_gene_positions[1:]):
        (gene1_start,gene1_end,gene1_name), (gene2_start,gene2_end,gene2_name) = gene1_data, gene2_data
        if gene1_end>=gene2_start:
            overlapping_gene_pairs.append((gene1_name,gene2_name))
        # check for "gene1 contains gene2", print a warning, since it can make other things not work right
        if gene1_end>=gene2_end:
            print("WARNING: gene %s is completely inside gene %s! "%(gene1_name, gene2_name)
                      +"Various gene-position-related results may be inaccurate.")
    return overlapping_gene_pairs
    # MAYBE-TODO rewrite it so it actually detects ALL overlaps?  Right now if gene A contains nonoverlapping genes B and C, it'll sort them as (A,B,C) since A starts first, so it'll detect the (A,B) overlap, but it won't detect the (A,C) overlap because it doesn't CHECK (A,C), only (A,B) and (B,C).  This could be fixed either by just brute-force checking all gene pairs (and then using DNA_basic_utilities.position_test_overlap), or by writing something prettier.  In any case, not a priority, since generally genes DON'T OVERLAP AT ALL.


def check_for_overlapping_features(mRNA, gene_name):
    """ Given a GFF gene record (from BCBio.GFF parser), return list of tuples of IDs of overlapping features.  
    Not perfect: may miss some overlaps if some features are completely contained within other features, 
     but if it returns no overlaps, there definitely aren't any."""
    overlap_found = False
    all_feature_positions = []
    for feature in mRNA.sub_features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_feature_positions.append((feature.location.start.position+1, feature.location.end.position, feature.id))
    all_feature_positions.sort()
    for feature1_data,feature2_data in itertools.izip(all_feature_positions,all_feature_positions[1:]):
        (_,feature1_end,feature1_name), (feature2_start,feature2_end,feature2_name) = feature1_data, feature2_data
        if feature1_end>=feature2_start:
            overlap_found = True
        # check for "feature1 contains feature2", print a warning, since it can make other things not work right
        if feature1_end>=feature2_end:
            print("WARNING: feature %s is completely inside feature %s in gene %s!"%(feature1_name,feature2_name,gene_name)
                      +" There may be additional undetected overlaps downstream.")
    return overlap_found


def find_min_gene_distance(sequence_record, starting_values=None):
    """ Given a GFF record from BCBio.GFF parser, return the lowest distance between adjacent genes, and the two gene IDs. 
    May not be accurate if some genes are completely contained inside other genes. """
    min_distance = len(sequence_record.seq) if starting_values is None else starting_values[0]
    min_gene1 = 'none' if starting_values is None else starting_values[1]
    min_gene2 = 'none' if starting_values is None else starting_values[2]
    all_gene_positions = []
    for gene in sequence_record.features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_gene_positions.append((gene.location.start.position, gene.location.end.position-1, gene.id))
    all_gene_positions.sort()
    for (_,gene1_end,gene1_name), (gene2_start,_,gene2_name) in itertools.izip(all_gene_positions,all_gene_positions[1:]):
        # subtract 1 from distance, so if gene1 is 1-4 and gene2 is 5-9 the distance is 0
        gene_distance = gene2_start - gene1_end - 1
        if gene_distance < min_distance:
            min_distance = gene_distance 
            min_gene1, min_gene2 = gene1_name, gene2_name
    return min_distance, min_gene1, min_gene2


def check_gene_coverage(sequence_records, check_for_overlap=True):
    """ Given a list of GFF sequence records (from BCBio.GFF parser), return the fraction of their length covered by genes.
    Also returns a feature_coverage_fractions dictionary (feature_name: fraction_of_gene_size), 
      but note that if the records were parsed with a gene-only limit, genes have no features and the dict will be empty.
    (Caution: length may be approximate! If gff parsing was done based on a fasta file, length is accurate; 
      otherwise it's an underestimate based on last feature position).
    """
    length_total = 0
    gene_length_total = 0
    total_length_by_feature = defaultdict(lambda: 0)
    for sequence_record in sequence_records:
        length_total += len(sequence_record.seq)
        for gene in sequence_record.features:
            gene_length_total += gene.location.end.position - gene.location.start.position
            # this section tries to keep track of subfeature types
            for feature in gene.sub_features:
                total_length_by_feature[feature.type] += len(feature)
                for subfeature in feature.sub_features:
                    total_length_by_feature[subfeature.type] += len(subfeature)
    gene_coverage_fraction = float(gene_length_total)/length_total
    feature_coverage_fractions = [(feature,float(length)/gene_length_total) for feature,length 
                                                                                  in total_length_by_feature.items()]

    # TODO the by-feature coverage doesn't work because I'm only parsing the file for genes, not features!!!  If I want to parse for features, I need to split things up into multiple passes etc again...
    #print total_length_by_feature

    # Check for overlapping genes and print a warning, since overlapping genes will make the measurement inaccurate
    if check_for_overlap:
        if check_for_overlapping_genes(sequence_record):
            print "WARNING: There are overlapping genes! %% of length covered by genes may not be accurate."
    # MAYBE-TODO actually adjust the measurement for overlapping genes?  Nah, too much work, not enough need for now.

    return gene_coverage_fraction, feature_coverage_fractions


### Independent functions not used in main - take gff file as input, parse it, get some outputs

def get_feature_start_end(feature_record):
    """ Get the start and end feature positions as numbers, 1-based, end-inclusive (so the first two bases are 1-2).

    BCBio uses 0-based end-exclusive positions, so the first two bases are 0-2.
    """
    return (feature_record.location.start.position+1, feature_record.location.end.position)


def get_gene_start_end_excluding_UTRs(gene_record, return_longest_if_multiple=False):
    """ Get the start of the first CDS feature and the end of the last one, (as numbers, 1-based, end-inclusive).

    If gene has no mRNAs or more than one, raise exception.

    BCBio uses 0-based end-exclusive positions, so the first two bases are 0-2; this uses 1-based end-inclusive, so they're 1-2.
    """
    if len(gene_record.sub_features) == 0:  
        raise NoRNAError("Gene %s has no RNA - can't determine CDS start/end!"%gene_record.id)
    if len(gene_record.sub_features) > 1 and not return_longest_if_multiple:   
        raise MultipleRNAError("Gene %s has multiple RNAs - can't determine single CDS start/end!"%gene_record.id)
    CDS_starts, CDS_ends = [], []
    for mRNA in gene_record.sub_features:
        features = mRNA.sub_features
        CDS_positions = [get_feature_start_end(feature) for feature in features if feature.type=='CDS']
        curr_CDS_starts, curr_CDS_ends = zip(*CDS_positions)
        CDS_starts += curr_CDS_starts
        CDS_ends += curr_CDS_ends
    return min(CDS_starts), max(CDS_ends)
    # TODO unit-test!


def gene_positions(genefile, include_chromosome=True, include_strand=True, coding_only=False, 
                   return_longest_if_multiple=False, ignore_strange_cases=False):
    """ Return a gene_ID:(chromosome, strand, start_pos, end_pos) dictionary based on GFF input file. 
    
    The positions are 1-based, end-inclusive. 
    If include_chromosome and/or include_strand is False, the corresponding values are missing from the output tuples.

    If coding_only is True, the start/end positions are the start and end of the first and last exon (i.e. excluding the UTRs). 
     In that case, if  a gene doesn't have an mRNA with exons, or has multiple mRNAs, raise an Exception, 
      unless ignore_strange_cases is True, then just don't include it in the output.
    """
    gene_positions = {}
    with open(os.path.expanduser(genefile)) as GENEFILE:
        # if coding_only is False, only look at genes, not sub-features
        genefile_parsing_limits = {'gff_type': ['gene']} if not coding_only else {}
        for chromosome_record in GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
            for gene_record in chromosome_record.features:
                # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3) - 
                #  convert to 1-based end-inclusive (so first-third base is bases 1,2,3, i.e. range 1-3)
                if include_chromosome:      full_pos_info = (chromosome_record.id,)
                else:                       full_pos_info = ()
                if include_strand:          full_pos_info += (GFF_strands[gene_record.strand],)
                if not coding_only:
                    full_pos_info += get_feature_start_end(gene_record)
                else:
                    try:    start_end = get_gene_start_end_excluding_UTRs(gene_record, return_longest_if_multiple)
                    except (NoRNAError, MultipleRNAError):
                        if ignore_strange_cases:    continue
                        else:                       raise
                    full_pos_info += start_end
                gene_positions[gene_record.id] = full_pos_info
    return gene_positions


def gene_lengths(genefile, exclude_UTRs=False, ignore_strange_cases=False):
    """ Return a gene_ID:gene_length dictionary based on GFF input file (full or coding length). 

    If exclude_UTRs is True, the length is counted from first to last coding exon (excluding UTRs). 
     In that case, if  a gene doesn't have an mRNA with exons, or has multiple mRNAs, raise an Exception, 
      unless ignore_strange_cases is True, then just don't include it in the output.
    """
    gene_lengths = {}
    for gene_ID, (start_pos, end_pos) in gene_positions(genefile, include_chromosome=False, include_strand=False, 
                                                coding_only=exclude_UTRs, ignore_strange_cases=ignore_strange_cases).iteritems():
        # gene_positions returns 1-based end-inclusive positions, so add 1 to get length (since a 1bp gene is 1-1)
        gene_lengths[gene_ID] = end_pos - start_pos + 1
    return gene_lengths


def gene_coding_exon_positions(genefile, ignore_strange_cases=False):
    """ Find list of positions of all coding exons for each gene.
    
    Return a gene_ID:coding_exons dict, where coding_exons is a list of (start_pos, end_pos) tuples. Based on GFF input file. 
    The positions are 1-based, end-inclusive. 

    If a gene doesn't have an mRNA or has multiple RNAs, raise an exception, 
        unless ignore_strange_cases is True, then just don't include those cases in the output.
    """
    # MAYBE-TODO add options for all features?  If needed.
    # MAYBE-TODO add option to deal with multiple-splice-variant cases?  But they might just make the data too messy.
    gene_feature_positions = {}
    with open(os.path.expanduser(genefile)) as GENEFILE:
        for chromosome_record in GFF.parse(GENEFILE, limit_info={}):
            for gene_record in chromosome_record.features:
                if len(gene_record.sub_features) != 1:  
                    if ignore_strange_cases:    continue
                    else:                       raise NoRNAError("Gene %s has no mRNA or multiple mRNAs!"%gene_record.id)
                features = gene_record.sub_features[0].sub_features
                CDS_positions = sorted(get_feature_start_end(feature) for feature in features if feature.type=='CDS')
                gene_feature_positions[gene_record.id] = CDS_positions
    return gene_feature_positions


def feature_total_lengths(genefile, ignore_multiple_splice_variants=False, genome_fasta_file=DEFAULT_NUCLEAR_GENOME_FILE):
    """ Return a dict containing the total lengths of all genes, exons, introns, 5' and 3' UTRs, and intergenic spaces. 

    Using the feature names from the gff file - usually "CDS" is used instead of "exon", 
     since technically "exon" would include the UTRs, and here we want to separate them.

    If  a gene doesn't have an mRNA, raise exception. 
    If it has multiple RNAs, raise exception, unless ignore_multiple_splice_variants is True, then just ignore those cases 
     for the purpose of total gene feature lengths (but not for total gene length).

    To get the correct total intergenic space size, input the path to a genome_fasta_file containing the chromosome sequences, 
     to get chromosome lengths, since the gff genefile file doesn't have those.  
     (The default value is the chlamy v4.3 genome fasta file.)
    """
    # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), so length is straightforward
    feature_total_lengths = defaultdict(int)
    total_genome_length = 0
    # get the total chromosome lengths from genome_fasta_file if given - needed to calculate total intergenic length 
    if genome_fasta_file:
        with open(genome_fasta_file) as FASTAFILE:  fasta_seq_dict = SeqIO.to_dict(SeqIO.parse(FASTAFILE, "fasta"))
    else:                                           fasta_seq_dict = {}
    with open(os.path.expanduser(genefile)) as GENEFILE:
        for chromosome_record in GFF.parse(GENEFILE, base_dict=fasta_seq_dict):
            total_genome_length += len(chromosome_record.seq)
            for gene_record in chromosome_record.features:
                gene_length = gene_record.location.end.position - gene_record.location.start.position
                feature_total_lengths['gene'] += gene_length
                if len(gene_record.sub_features) == 0:
                    raise NoRNAError("Gene %s has no RNA - can't determine CDS start/end!"%gene_record.id)
                if len(gene_record.sub_features) > 1:
                    if ignore_multiple_splice_variants:    
                        feature_total_lengths['MULTIPLE_SPLICE_VARIANTS'] += gene_length
                        continue
                    else:                       
                        raise MultipleRNAError("Gene %s has multiple RNAs - can't determine single CDS start/end!"%gene_record.id)
                else:
                    features = gene_record.sub_features[0].sub_features
                    for feature in features:
                        feature_total_lengths[feature.type] += (feature.location.end.position - feature.location.start.position)
    # calculate total intron length from total gene and exon/UTR length
    feature_total_lengths['intron'] = feature_total_lengths['gene'] - sum(length for feature,length 
                                                                          in feature_total_lengths.items() if feature != 'gene')
    # calculate total intergenic length from total genome and gene length
    feature_total_lengths['all'] = total_genome_length
    feature_total_lengths['intergenic'] = total_genome_length - feature_total_lengths['gene']
    return dict(feature_total_lengths)


### Functions for parsing JGI GFF2-format file 

def parse_JGI_GFF2_file(infile_GFF2, stop_on_error=False):
    """ Given a JGI-style GFF2 file, yield data_fields,annotation_dict tuples for each line. 

    The input should be the filename.
    Data_fields is just a list of the tab-separated fields (including the annotation field as a raw string).
    Annotation_dict is a tag:value dictionary of the annotations from the last field.

    If stop_on_error is True, raise a ValueError when can't parse an annotation field; 
     otherwise print an error line and keep going, leaving the annotation parsing for that line unfinished.
    """ 
    for line in open(infile_GFF2):
        fields = line.strip().split('\t')
        annotations = {}
        for ann in fields[8].split('; '):
            try:
                if '"' in ann:
                    key, val = ann.split(' "')
                    val = val.strip('"')
                else:
                    key, val = ann.split(' ')
                annotations[key] = val
            except ValueError:
                error_msg = "Can't parse this annotation into a key:value pair! %s"%ann
                if stop_on_error:   raise ValueError(error_msg)
                else:               print "ERROR: %s"%error_msg
                continue
        yield fields[:8], annotations
    # TODO unit-test?


def name_to_ID_dicts(infile_GFF2, convert_to_singles=True):
    """ Given a JGI-style GFF2 file, return a set of all gene names, and name:transcriptID and name:proteinID dictionaries.

    The input should be the filename.

    Assumes the file has a 'name' annotation on each line (prints an error if not found), 
     and 'transcriptId' and 'proteinId' on some lines; ignores other annotations.

    If sometimes a single name has more than one transcript/protein, both returned dictionaries have sets as values; 
     if there's always exactly one, they'll be converted to name:ID dictionaries if convert_to_singles is True.
    """
    # go over the file, parse each line and make name:proteinID and name:transcriptID dictionaries
    all_names = set()
    name_to_transcriptIDs = defaultdict(set)
    name_to_proteinIDs = defaultdict(set)
    for fields, annotations in parse_JGI_GFF2_file(infile_GFF2):
        if 'name' not in annotations:
            print "ERROR: line has no 'name' tag! %s"%line
            continue
        name = annotations['name'].strip('"')
        all_names.add(name)
        if 'transcriptId' in annotations:
            name_to_transcriptIDs[name].add(annotations['transcriptId'])
        if 'proteinId' in annotations:
            name_to_proteinIDs[name].add(annotations['proteinId'])
    # print some summaries
    print "Total %s gene names in file; %s have transcript IDs, %s have protein IDs."%(len(all_names), 
                           value_and_percentages(len(name_to_transcriptIDs), [len(all_names)]), 
                           value_and_percentages(len(name_to_proteinIDs), [len(all_names)]))
    genes_with_both_IDs = set(name_to_transcriptIDs) & set(name_to_proteinIDs)
    N_same = sum(name_to_proteinIDs[name]==name_to_transcriptIDs[name] for name in genes_with_both_IDs)
    print "Out of genes that have both (%s), the protein ID and transcript ID are the same for %s."%(len(genes_with_both_IDs), 
                     value_and_percentages(N_same, [len(genes_with_both_IDs)]))
    # See how many unique protein/transcript IDs each name has
    N_transcript_IDs = Counter(len(x) for x in name_to_transcriptIDs.values()) 
    print "Number of unique transcript IDs per gene name, and how many genes have that number: %s"%(
        ', '.join('%s - %s'%(N_IDs,N_genes) for N_IDs,N_genes in sorted(N_transcript_IDs.items())))
    N_protein_IDs = Counter(len(x) for x in name_to_proteinIDs.values())    
    print "Number of unique protein IDs per gene name, and how many genes have that number: %s"%(
        ', '.join('%s - %s'%(N_IDs,N_genes) for N_IDs,N_genes in sorted(N_protein_IDs.items())))
    # Make new dictionaries that just have a single protein/transcript ID per name, instead of a set
    if list(N_transcript_IDs.keys()) == [1]:
        name_to_transcriptIDs = {name: transcriptIDs.pop() for name, transcriptIDs in name_to_transcriptIDs.items()}
    if list(N_protein_IDs.keys()) == [1]:
        name_to_proteinIDs = {name: proteinIDs.pop() for name, proteinIDs in name_to_proteinIDs.items()}
    return all_names, name_to_transcriptIDs, name_to_proteinIDs
    # TODO unit-test?


def gene_position_dict(infile_GFF2):
    """ Given a JGI-style GFF2 file, return chromosome:proteinID:(start,end) dictionary; ignores strand and features.

    The input should be the filename.
    """
    # go over the file, parse each line and grab chromosome/start/end for each 
    name_to_positions = defaultdict(set)
    name_to_proteinID = {}
    name_to_chromosome = {}
    for fields, annotations in parse_JGI_GFF2_file(infile_GFF2):
        try:
            name = annotations['name'].strip('"')
        except KeyError:
            print "ERROR: line has no 'name' tag! %s"%line
            continue
        name_to_positions[name].update([int(x) for x in fields[3:5]])
        add_to_dict_no_replacement(name_to_chromosome, name, fields[0], 'gene', 'chromosome', raise_exception=False)
        if 'proteinId' in annotations:
            proteinID = annotations['proteinId']
            add_to_dict_no_replacement(name_to_proteinID, name, proteinID, 'gene', 'proteinID', raise_exception=False)
    # Now make the final dictionary, with the smallest and largest position for each gene, and using proteinIDs
    #  (note: GFF positions are 1-based inclusive, so no need to adjust the positions, see
    #   https://www.sanger.ac.uk/resources/software/gff/spec.html)
    chromosomes_to_proteinIDs_to_positions = defaultdict(dict)
    for gene, positions in name_to_positions.items():
        chromosome = name_to_chromosome[gene]
        proteinID = name_to_proteinID[gene]
        start = min(positions)
        end = max(positions)
        chromosomes_to_proteinIDs_to_positions[chromosome][proteinID] = (start, end)
    return chromosomes_to_proteinIDs_to_positions
    # TODO unit-test?


######### Test code #########

class Testing_(unittest.TestCase):
    """ Unit-tests for some of the additional functions not used in main (main has a run-test, that's enough). """

    def test__gene_lengths(self):
        output = gene_lengths('test_data/INPUT_gene-data-1_all-cases.gff3')
        assert len(output) == 28
        assert output['test.geneA0_proper_plus'] == output['test.geneA1_proper_minus'] == 700
        assert output['test.geneB5_only_exon'] == 100
        assert output['test.geneD8a_bad_overlapping_genes'] == 400

    # TODO unit-test the other non-main functions!


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_folder = "test_data"

    infile = "test_data/INPUT_gene-data-1_all-cases.gff3"
    tests = [("gff-ex__gene-data-analysis", "-E -n -1 %s"%infile)]

    parser = define_option_parser()
    argument_converter = lambda parser,options,args: (args[0], options, args[1])
    return run_functional_tests(tests, parser, main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-r', '--record_structure', action="store_true", default=False, 
                      help='Show the record structures (for example gene->mRNA->CDS/UTR). Default %default') 
    parser.add_option('-R', '--no_record_structure', action="store_false", dest='record_structure')

    parser.add_option('-c', '--feature_type_counts', action="store_true", default=True, 
                      help='Count the number of feature types in file (gene, mRNA, exon, etc). Default %default') 
    parser.add_option('-C', '--no_feature_type_counts', action="store_false", dest='feature_type_counts')

    parser.add_option('-g', '--gene_counts', action="store_true", default=False, 
                      help="Count genes per chromosome, and the approximate fraction of each chromosome covered by genes. "
                          +"Default %default") 
    parser.add_option('-G', '--no_gene_counts', action="store_false", dest='gene_counts')
    parser.add_option('-a', '--fasta_sequence_file', default='', metavar='FILE', 
                      help="Fasta file containing the sequences listed in gff infile (default %default).")
    parser.add_option('-d', '--print_seq_details', action="store_true", default=False, 
                      help='Print full GFF details for each chromosome (only if -g). Default %default') 
    parser.add_option('-D', '--no_print_seq_details', action="store_false", dest='print_seq_details')

    parser.add_option('-o', '--check_gene_overlaps', action="store_true", default=True, 
                      help='Check for overlapping genes, distances, ID uniqueness, etc; count genes. Default %default') 
    parser.add_option('-O', '--no_check_gene_overlaps', action="store_false", dest='check_gene_overlaps')

    parser.add_option('-f', '--gene_feature_structure_counts', action="store_true", default=True, 
                      help='Give gene counts by UTR/exon count/order; check feature distances/overlaps. Default %default') 
    parser.add_option('-F','--no_gene_feature_structure_counts', action="store_false",dest='gene_feature_structure_counts')
    parser.add_option('-u', '--full_feature_structures', action="store_true", default=False, 
                      help='With -f option, show full as well as simplified feature structures. Default %default') 
    parser.add_option('-U','--no_full_feature_structures', action="store_false",dest='full_feature_structures')
    parser.add_option('-n', '--genes_to_display', type="int", default=5, metavar='N', 
                      help="When showing gene counts per group (-f), show N example genes (-1: all) (default %default).")
    parser.add_option('-e', '--exon_number_cutoff', type="int", default=30, metavar='N', 
                      help="When categorizing genes by exon number, lump together all above N (default %default).")
    parser.add_option('-Y', '--N_detail_run_groups', type="int", default=5, metavar='N', 
                      help="How many passes to split reading the file into with -f option (default %default) "
                          +"- may take a lot of memory (and CPU) if read in a single pass; too many passes waste CPU.")

    parser.add_option('-s', '--source_counts', action="store_true", default=False, 
                      help='Count the features by source (not very useful unless file is mixed-source). Default %default') 
    parser.add_option('-S', '--no_source_counts', action="store_false", dest='source_counts')

    parser.add_option('-l', '--all_gff_limits', action="store_true", default=False, 
                      help='Output all feature counts: by type, source (-cs), chromosome, maybe other? Default %default')
    parser.add_option('-L', '--no_all_gff_limits', action="store_false", dest='all_gff_limits')

    parser.add_option('-E', '--everything', action='store_true', default=False, 
                      help="Examine the infile in ALL implemented ways (turn on all the True/False options).")

    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on test input file, check output against reference. Ignores all other options/arguments.")

    return parser


# Help function for making sure the eval/exec-based -E implementation is right
_sort_dict_string = lambda dict_repr: sorted(str(dict_repr).strip('{}').split(', '))


def main(infile, options, outfile=None):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument is generated by an optparse parser.  Prints to outfile if given, otherwise stdout.
    """
    examiner = GFF.GFFExaminer()

    # turn all True/False options to True (requires eval magic to detect all options)
    if options.everything:
        option_dict = eval(str(options))
        error_text = "The -E option isn't working right, turn everything on by hand!"
        assert _sort_dict_string(option_dict) == _sort_dict_string(options), error_text
        for (key,value) in option_dict.items():
            if value is False:      # must use 'is', not '==': '0==False' is true, but '0 is False' isn't
                exec('options.%s = True'%key)
                assert eval('options.%s'%key) is True, error_text

    # writing either to sys.stdout or to outfile
    with (sys.stdout if outfile is None else open(outfile, 'w')) as OUT:

        # I'm not sure why I need to open the file separately for each operation, but it doesn't work otherwise...

        if options.record_structure:
            with open(infile) as INFILE:
                OUT.write("\n *** Record structures ***\n")
                pprint.pprint(examiner.parent_child_map(INFILE),stream=OUT)

        if options.feature_type_counts:
            with open(infile) as INFILE:
                OUT.write("\n *** Type counts ***\n")
                pprint.pprint({'gff_type': examiner.available_limits(INFILE)['gff_type']},stream=OUT)

        if options.source_counts:
            OUT.write("\n *** Source and source/type counts ***\n")
            with open(infile) as INFILE:
                pprint.pprint({'gff_source': examiner.available_limits(INFILE)['gff_source']},stream=OUT)
            with open(infile) as INFILE:
                pprint.pprint({'gff_source_type': examiner.available_limits(INFILE)['gff_source_type']},stream=OUT)

        if options.all_gff_limits:
            with open(infile) as INFILE:
                OUT.write("\n *** All GFF file limit field values ***\n")
                pprint.pprint(examiner.available_limits(INFILE),stream=OUT)

        if options.gene_counts or options.print_seq_details or options.check_gene_overlaps:
            if options.gene_counts or options.print_seq_details:
                OUT.write("\n *** Gene and other data per chromosome ***\n")
            if options.gene_counts and not options.fasta_sequence_file:
                OUT.write("       (Caution: approximate sequence length is calculated to last gene only!)\n")
            if options.gene_counts or options.print_seq_details:
                OUT.write("\n")
            if options.check_gene_overlaps:
                total_chromosomes = 0
                total_genes = 0
                overlapping_gene_pairs = []
                min_gene_distance_data = None
                gene_IDs = []
            # get the sequences from the fasta file if present, otherwise leave an empty dict to overlay the gff data on
            if options.fasta_sequence_file:
                with open(options.fasta_sequence_file) as FASTAFILE:
                    fasta_seq_dict = SeqIO.to_dict(SeqIO.parse(FASTAFILE, "fasta"))
            else:
                fasta_seq_dict = {}
            with open(infile) as INFILE:
                all_records = GFF.parse(INFILE, limit_info={'gff_type': ['gene']}, base_dict=fasta_seq_dict)
                if options.check_gene_overlaps:
                    total_length = 0
                for record in all_records:
                    if options.gene_counts or options.print_seq_details:
                        OUT.write(" * sequence %s: %s genes\n"%(record.id, len(record.features)))
                    if options.gene_counts:
                        gene_coverage_fraction, feature_coverage_fractions = check_gene_coverage([record],False)
                        gene_coverage_percent = "%.0f%%"%(100*gene_coverage_fraction)
                        OUT.write("length %s bp, with %s covered by genes\n"%(len(record.seq), gene_coverage_percent))
                        for feature,fraction in feature_coverage_fractions:
                            OUT.write(" * %.1f%% of genes covered by %s\n"%(fraction*100, feature))
                    if options.print_seq_details:               
                        OUT.write("GFF parser details:\n")
                        OUT.write("%s\n"%record)
                    if options.gene_counts or options.print_seq_details:
                        OUT.write("\n")
                    if options.check_gene_overlaps:
                        total_chromosomes += 1
                        total_length += len(record.seq)
                        total_genes += len(record.features)
                        overlapping_gene_pairs += check_for_overlapping_genes(record)
                        min_gene_distance_data = find_min_gene_distance(record, min_gene_distance_data)
                        gene_IDs += [gene.id for gene in record.features]

        if options.check_gene_overlaps:
            OUT.write("\n *** Gene overlaps ***\n")
            with open(infile) as INFILE:
                all_records = GFF.parse(INFILE, limit_info={'gff_type': ['gene']}, base_dict=fasta_seq_dict)
                gene_coverage_fraction, feature_coverage_fractions = check_gene_coverage(all_records,False)
                total_gene_coverage_percent = "%.0f%%"%(100*gene_coverage_fraction)
            OUT.write("Total %s genes on %s chromosomes.\n"%(total_genes, total_chromosomes))
            OUT.write("Total genome length %sbp, with %s covered by genes.\n"%(total_length, total_gene_coverage_percent))
            for feature,fraction in feature_coverage_fractions:
                OUT.write(" * %.1f%% of genes covered by %s\n"%(fraction*100, feature))

            # TODO calculate the number of ambiguous bases, inside and outside genes

            for geneA,geneB in overlapping_gene_pairs:  OUT.write("Overlapping gene pair!  IDs: %s, %s.\n"%(geneA,geneB))
            if not overlapping_gene_pairs:              OUT.write("No overlapping genes.\n")
            OUT.write("Minimum distance between two genes is %s (genes %s, %s).\n"%min_gene_distance_data)
            IDs_are_unique = True
            for gene_ID,count in count_list_values(gene_IDs).iteritems():
                if count>1:    
                    OUT.write("Non-unique gene ID! Gene %s occurs %s times.\n"%(gene_ID,count))
                    IDs_are_unique = False
            if IDs_are_unique is True: 
                OUT.write("All gene IDs are unique.\n")


        if options.gene_feature_structure_counts:
            OUT.write("\n *** Gene counts by feature structure ***\n")
            genes_by_feature_structure = defaultdict(set)       # set() returns an empty set
            genes_with_gene_mRNA_gap = set()
            genes_with_mRNA_feature_gap = set()
            overlapping_feature_genes = []

            # TODO calculate total gene area covered by exons/UTRs/introns?  (and print a warning that it only applies 
            #     to genes with a normal structure, i.e. a single mRNA covering the whole gene.)

            ## Go over subsets of chromosomes at once, to avoid reading the whole file into memory at once
            # First get the list of all chromosomes in the file, WITHOUT reading it all into memory
            with open(infile) as INFILE:
                GFF_limit_data = examiner.available_limits(INFILE)
                chromosomes_and_counts = dict([(c,n) for ((c,),n) in GFF_limit_data['gff_id'].items()])

            # Now lump the chromosomes into N_run_groups sets with the feature counts balanced between sets, 
            #  to avoid using too much memory (by reading the whole file at once), 
            #   or using too much time (by reading the whole file for each chromosome/scaffold)
            chromosome_sets = split_into_N_sets_by_counts(chromosomes_and_counts, options.N_detail_run_groups)

            ### go over all mutants on each chromosome, figure out which gene they're in (if any), keep track of totals
            # keep track of all the mutant and reference chromosomes to catch chromosomes that are absent in reference
            for chromosome_set in chromosome_sets:
                genefile_parsing_limits = {'gff_id': list(chromosome_set)}
                with open(infile) as INFILE:
                    for chromosome_record in GFF.parse(INFILE, limit_info=genefile_parsing_limits):
                        for gene in chromosome_record.features:
                            # figure out subfeature structure (use special cases if it's unexpected)
                            if len(gene.sub_features)==0:
                                genes_by_feature_structure[('NO_mRNA',)].add(gene.id)
                            elif len(gene.sub_features)>1:
                                genes_by_feature_structure[('MULTIPLE_mRNAs',)].add(gene.id)
                            else:
                                [mRNA] = gene.sub_features
                                if check_for_overlapping_features(mRNA, gene.id):
                                    overlapping_feature_genes.append(gene.id)
                                if not mRNA.type=='mRNA':
                                    genes_by_feature_structure[('NON_mRNA_PRIMARY_FEATURE',)].add(gene.id)
                                elif len(mRNA.sub_features)==0:
                                    genes_by_feature_structure[('NO_mRNA_SUBFEATURES',)].add(gene.id)
                                else:
                                    # The features are NOT SORTED in a normal gff file!!  Need to sort them.
                                    features_by_pos = sorted([(f.location.start.position,f.type) 
                                                              for f in mRNA.sub_features])
                                    # reverse structure if the gene's on the minus strand, since it's sorted by position
                                    if mRNA.strand in [-1, '-']:    features_by_pos.reverse()
                                    feature_structure = tuple([t for (p,t) in features_by_pos])
                                    genes_by_feature_structure[feature_structure].add(gene.id)
                                # check for gaps between gene, mRNA, and feature starts/ends
                                # MAYBE-TODO the check isn't done if there are multiple mRNAs! but that never happens.
                                mRNA_start, mRNA_end = mRNA.location.start.position, mRNA.location.end.position
                                if mRNA_start != gene.location.start.position:
                                    genes_with_gene_mRNA_gap.add(gene.id)
                                if mRNA_end != gene.location.end.position:
                                    genes_with_gene_mRNA_gap.add(gene.id)
                                if len(mRNA.sub_features)>0:
                                    if mRNA_start != min([f.location.start.position for f in mRNA.sub_features]):
                                        genes_with_mRNA_feature_gap.add(gene.id)
                                    if mRNA_end != max([f.location.end.position for f in mRNA.sub_features]):
                                        genes_with_mRNA_feature_gap.add(gene.id)


            for gene in overlapping_feature_genes:      
                OUT.write("Overlapping features in gene!  Gene ID: %s.\n"%gene)
            if not overlapping_feature_genes:           
                OUT.write("No genes have overlapping features.\n")

            for gene in genes_with_gene_mRNA_gap:       
                OUT.write("Gene has a gap between gene and mRNA start or end!  %s.\n"%gene)
            if not genes_with_gene_mRNA_gap:            
                OUT.write("The mRNA position is the same as gene position for all genes.\n")
            for gene in genes_with_mRNA_feature_gap:    
                OUT.write("Gene has a gap between mRNA and feature start or end!  %s.\n"%gene)
            if not genes_with_mRNA_feature_gap:         
                OUT.write("No gap between mRNA edges and first/last features in any genes.\n")

            # make new dict with all rows of 'CDS' changed to a single 'exon/s' to end up with fewer structure variants, 
            #   also count genes with different exon numbers
            exon_number_gene_counts = defaultdict(set)              # passing the set function: equivalent to lambda: set()
            UTR_5prime_number_gene_counts = defaultdict(set)
            UTR_3prime_number_gene_counts = defaultdict(set)
            genes_by_simpler_feature_structure = defaultdict(set)
            for feature_structure,gene_set in genes_by_feature_structure.iteritems():
                simple_feature_structure = []
                exon_count, UTR_5prime_count, UTR_3prime_count = 0,0,0
                for feature in feature_structure:
                    if feature=='CDS':
                        exon_count += 1
                        if simple_feature_structure==[] or not simple_feature_structure[-1]=='exon/s':
                            simple_feature_structure.append('exon/s')
                    elif feature=='five_prime_UTR':
                        UTR_5prime_count += 1
                        if simple_feature_structure==[] or not simple_feature_structure[-1]=="5'UTR/s":
                            simple_feature_structure.append("5'UTR/s")
                    elif feature=='three_prime_UTR':
                        UTR_3prime_count += 1
                        if simple_feature_structure==[] or not simple_feature_structure[-1]=="3'UTR/s":
                            simple_feature_structure.append("3'UTR/s")
                    else:
                        simple_feature_structure.append(feature)
                genes_by_simpler_feature_structure[tuple(simple_feature_structure)] |= gene_set
                UTR_5prime_number_gene_counts[UTR_5prime_count] |= gene_set
                UTR_3prime_number_gene_counts[UTR_3prime_count] |= gene_set
                if exon_count >= options.exon_number_cutoff:
                    exon_number_gene_counts['%s+'%options.exon_number_cutoff] |= gene_set
                else:
                    exon_number_gene_counts[exon_count] |= gene_set

            N = options.genes_to_display if options.genes_to_display>=0 else None

            OUT.write( "\n * Gene counts by simplified feature structure (adjacent exons/UTRs combined) "
                  +"(total %s structures)\n"%len(genes_by_simpler_feature_structure))
            genes_by_simpler_feature_structure = [(len(s),n,s) for (n,s) in genes_by_simpler_feature_structure.items()]
            genes_by_simpler_feature_structure.sort(reverse=True)
            for count,feature_structure,gene_set in genes_by_simpler_feature_structure:
                OUT.write("%s genes [%s]:  %s\n"%(count, ', '.join(feature_structure), ', '.join(list(gene_set)[:N])))

            OUT.write("\n * Gene counts by exon number (total %s values)\n"%len(exon_number_gene_counts))
            for exon_number,gene_set in sorted(exon_number_gene_counts.items()):
                OUT.write("%s genes with %s exons:  %s\n"%(len(gene_set), exon_number, ', '.join(list(gene_set)[:N])))

            OUT.write("\n * Gene counts by 5'UTR number (total %s values)\n"%len(UTR_5prime_number_gene_counts))
            for UTR_number,gene_set in sorted(UTR_5prime_number_gene_counts.items()):
                OUT.write("%s genes with %s 5'UTRs:  %s\n"%(len(gene_set), UTR_number, ', '.join(list(gene_set)[:N])))

            OUT.write("\n * Gene counts by 3'UTR number (total %s values)\n"%len(UTR_3prime_number_gene_counts))
            for UTR_number,gene_set in sorted(UTR_3prime_number_gene_counts.items()):
                OUT.write("%s genes with %s 3'UTRs:  %s\n"%(len(gene_set), UTR_number, ', '.join(list(gene_set)[:N])))

            if options.full_feature_structures:
                genes_by_feature_structure = [(len(s),x,s) for (x,s) in genes_by_feature_structure.items()]
                genes_by_feature_structure.sort(reverse=True)
                OUT.write("\n * Gene counts by full feature structure (total %s structures)\n" 
                          %len(genes_by_feature_structure))
                for count,feature_structure,gene_set in genes_by_feature_structure:
                    OUT.write("%s genes [%s]:  %s\n"%(count, ', '.join(feature_structure), ', '.join(list(gene_set)[:N])))



if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    try:                
        [infile] = args
    except ValueError:  
        parser.print_help()
        sys.exit("\nError: exactly one input gff file required!")

    main(infile, options)


