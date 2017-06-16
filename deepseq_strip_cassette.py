#! /usr/bin/env python2.7

"""
Take LEAP-Seq cassette-side read files (fastq), which should be <IB><cassette><genomic-flanking-seq>, and look for the full or truncated cassette sequence after the expected IB length (using bowtie2); output the IB and flanking seqs to separate fasta output files; output truncated cassette lengths and error info to another output file, and a summary to stdout.
Output files will be <outfile_base>_IB.fa, <outfile_base>_flank.fa, <outfile_base>_no-cassette.fa and <outfile_base>_cassette-info.txt (tab-separated).
 -- Weronika Patena, 2015
USAGE: deepseq_strip_cassette [options] infile outfile_base
"""

# standard library
from __future__ import division
import sys, os
import unittest
import collections
# other packages
# my modules
import basic_seq_utilities
import general_utilities
import deepseq_utilities

# to see how I figured out the default bowtie options, see ~/experiments/arrayed_library/1504_trying_new_cassette-trimming/notes.txt
DEFAULT_BOWTIE_OPTIONS = "--local --all --ma 3 --mp 5,5 --np 1 --rdg 5,3 --rfg 4,3 --score-min C,20,0 -N0 -L5 -i C,1,0 -R5 -D30 --norc --reorder"
# another possibly good score set based on that file is match=12, mismatch=23, read_gap=20, read_ext=10, ref_gap=17, ref_ext=10, score_min=80?  This may be better at avoiding multiple alignments with same scores?
# this is REALLY WEIRD, but for some reason if I put "-N0 -L5" in the options it works, but in the other order "-L5 -N0" it doesn't work on the same test sequences... (truncated_20bp_cassette in test_data/INPUT_strip_cassette.fa)


DEBUG = 0

class CassetteStrippingError(Exception): pass

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-C', '--cassette_index', default=None, metavar='index',
                      help="Bowtie2 index containing the expected cassette sequence (5' or 3')")
    parser.add_option('-l', '--allowed_IB_lengths', default='21,23', metavar='N,M', 
                      help="How many IB bases to expect before the cassette (format: 'min,max') (default %default).")
    parser.add_option('-b', '--extra_bowtie_options', default="", metavar='"options"', 
                      help="Additional bowtie options, which will be appended to the default options (%s) "%DEFAULT_BOWTIE_OPTIONS
                          +"- the later options will supersede the earlier ones.")
    parser.add_option('-B', '--existing_bowtie_file', default=None, metavar='samfile',
                      help="Existing bowtie2 output file to use instead of running bowtie2 (default %default). "
                          +"Make sure it was run with the same parameters as this! (-x = -C, -3 = -i, =5 = -l[0] minus one)")
    parser.add_option('-m', '--min_cassette_length', type='int', default=10, metavar='N', 
                      help="minimum length of detected cassette fragment (shorter = untrimmed) (default %default).")
    parser.add_option('-i', '--N_ignore_end', type='int', default=30, metavar='N', 
                      help="don't try to align N last bases to the cassette (default %default).")
    parser.add_option('-e', '--max_percent_errors', type='int', default=30, metavar='N', 
                      help="Maximum percent errors (default %default) - "
                          +"NOTE that quality should mainly be determined through bowtie min-score, this is just a secondary check.")
    parser.add_option('-s', '--max_allowed_cassette_start', type='int', default=5, metavar='N', 
                      help="Highest cassette position the alignment can start at (default %default).")
    parser.add_option('-n', '--n_threads', type='int', default=3, metavar='N', 
                      help="number of threads to run bowtie as (default %default).")

    ### cosmetic options
    parser.add_option('-q','--quiet', action='store_true', default=False, 
                      help="Don't print anything to stdout (default %default).")
    return parser


def run_bowtie(infile, bowtie_outfile, cassette_index, N_trim_start=0, N_trim_end=0, n_threads=1, extra_bowtie_options="", 
               quiet=True):
    """ Run bowtie2 to align infile against cassette_index, which should only contain the cassette seq expected in the reads.
    """
    bowtie_command = "bowtie2 %s %s -5 %s -3 %s -p %s -x %s -U %s -S %s"%(DEFAULT_BOWTIE_OPTIONS, extra_bowtie_options, 
                                                  N_trim_start, N_trim_end, n_threads, cassette_index, infile, bowtie_outfile)
    if not quiet:
        print "Bowtie command: %s"%bowtie_command
    bowtie_output = general_utilities.run_command_get_output(bowtie_command, shell=True)
    return bowtie_output


def get_aln_startpos_from_CIGAR(CIGAR):
    """ Return alignment start position (of the read, not the reference), 1-based - based on softclipping in CIGAR string.
    """
    return CIGAR[0].size+1 if CIGAR[0].type == 'S' else 1


def get_aln_len_from_CIGAR(seqlen, CIGAR):
    """ Return alignment length (of the read, not the reference) - based on readlen minus softclipping from CIGAR string.
    """
    return seqlen - sum(c.size for c in CIGAR if c.type=='S')


def choose_best_alignment(alignment_list, allowed_start_positions, min_cassette_length, max_percent_errors, 
                          max_allowed_cassette_start):
    """ Return best alignment from a list of HTSeq-parsed SAM alignments, or None if unaligned or if none meet the criteria.

    First, all valid alignments must meet the following criteria:
        - the alignment start position in the READ must be x+1, where x is one of allowed_start_positions
        - the alignment start position in the REFERENCE must be at most max_allowed_cassette_start
        - the alignment length is at least min_cassette_length
        - the edit distance is at most max_percent_errors of the alignment length
    If there are no valid alignments (or no alignments at all), return None. 

    Otherwise choose and return the best valid alignment:
        - highest alignment score
        - if multiple have the same score, pick the one with a longer read sequence being part of the alignment
        - if there are still multiples, raise an exception.
    """
    # filter alignments to only get valid ones; if there aren't any, return None.
    filtered_alignments = []
    for aln in alignment_list:
        if not aln.aligned:                             continue
        startpos = get_aln_startpos_from_CIGAR(aln.cigar)
        if startpos not in allowed_start_positions:     continue
        ref_startpos = [c for c in aln.cigar if c.type=='M'][0].ref_iv.start + 1
        if ref_startpos > max_allowed_cassette_start:     continue
        aln_len = get_aln_len_from_CIGAR(len(aln.read), aln.cigar)
        if aln_len < min_cassette_length:               continue
        edit_dist = aln.optional_field('NM')
        if edit_dist/aln_len*100 > max_percent_errors:  continue
        filtered_alignments.append(aln)
    if not filtered_alignments:
        return None
    # now try to find a single best alignment based on a ranked list of criteria; if there are multiple best, error.
    # NOTE that this should actually almost never be necessary, because bowtie2 only returns all DISTINCT alignments, 
    #   meaning ones that don't align the same read base to the same ref base - so to get multiple post-filter alns here, 
    #   the read would have to have alignments to two different cassette areas, 
    #   but if both the read alignment start and the cassette alignment start are constrained to a few bp, this is impossible.
    max_score = max(aln.optional_field('AS') for aln in filtered_alignments)
    filtered_alignments = [aln for aln in filtered_alignments if aln.optional_field('AS') == max_score]
    if len(filtered_alignments) == 1:
        return filtered_alignments[0]
    max_aln_len = max(get_aln_len_from_CIGAR(len(aln.read), aln.cigar) for aln in filtered_alignments)
    filtered_alignments = [aln for aln in filtered_alignments if get_aln_len_from_CIGAR(len(aln.read), aln.cigar) == max_aln_len]
    if len(filtered_alignments) == 1:
        return filtered_alignments[0]
    else:
        raise CassetteStrippingError("Multiple alignments with same quality and length! %s %s - "%(aln.read.name, aln.read.seq) 
                                     +"showing CIGAR, NM and MD strings: " 
                                     +', '.join("%s %s %s"%(aln.original_sam_line.split('\t')[5], aln.optional_field('NM'),  
                                                            aln.optional_field('MD')) for aln in filtered_alignments))


def parse_cassette_alignment(infile, bowtie_outfile, outfile_base, allowed_IB_lengths, min_cassette_length, max_percent_errors, 
                             max_allowed_cassette_start, N_trim_start=0, N_trim_end=0, quiet=False):
    """ 
    Parse bowtie2 output file to separate infile seqs into IB,cassette,flank triples (or untrimmed cases)
    
    Keep track of stripped cassette lengths/errors, and print summary to stdout.
    """
    total_seqs = 0
    results_N_unstripped = 0
    results_IB_lengths = collections.Counter()
    # use length ranges for results_cassette_read_lengths and results_reflengths, to avoid huge result lists
    # MAYBE-TODO add a check for any unusually common lengths in each range?  Or include BOTH ranges and singles?
    length_ranges = {}
    for i in range(21):     length_ranges[i] = '0-20'
    for i in range(21, 31): length_ranges[i] = '21-30'
    for i in range(31, 41): length_ranges[i] = '31-40'
    for i in range(41, 46): length_ranges[i] = i
    for i in range(46, 80): length_ranges[i] = '46+'
    results_cassette_read_lengths = collections.Counter()
    results_reflengths = collections.Counter()
    results_refstartpos = collections.Counter()
    results_total_readlen, results_total_errors = 0, 0
    # MAYBE-TODO count seqs with 0,1,2,3+ errors?

    with open(outfile_base+'_cassette-info.txt', 'w') as OUTFILE_cassette:
        OUTFILE_cassette.write("seq_header\tcassette_read_length\tIB_length\tcassette_seq"
                               +"\tN_errors\tref_startpos\tref_aln_length\tCIGAR_field\tMD_field\talignment_score\n")
        # MAYBE-TODO change field order to put ref_aln_length next to cassette_read_length or something?
        # MAYBE-TODO add more detailed error fields? #mismatches, #insertions, #deletions, total length of insertions and deletions
        with open(outfile_base+'_IB.fa', 'w') as OUTFILE_IB:
          with open(outfile_base+'_flank.fa', 'w') as OUTFILE_flank:
            with open(outfile_base+'_no-cassette.fa', 'w') as OUTFILE_no_cassette:
                # Need to parse the fastq infile in parallel with sam, because some first/last bases were trimmed before alignment
                #  and we want them back in the IB/flank output.
                for (name, seq, alignment_list) in deepseq_utilities.parse_fastx_sam_parallel(infile, bowtie_outfile):
                    total_seqs += 1
                    if DEBUG: print name
                    # bowtie output will have multiple alignments per sequence - need to pick the best one.
                    allowed_start_positions = [x-N_trim_start+1 for x in allowed_IB_lengths]    # +1 to make them 1-based
                    aln = choose_best_alignment(alignment_list, allowed_start_positions, min_cassette_length, max_percent_errors, 
                                                max_allowed_cassette_start)
                    if DEBUG: print aln
                    if aln is None:
                        results_N_unstripped += 1
                        OUTFILE_no_cassette.write('>%s\n%s\n'%(name, seq))
                    else:
                        aln_startpos = get_aln_startpos_from_CIGAR(aln.cigar)
                        aln_len = get_aln_len_from_CIGAR(len(aln.read), aln.cigar)
                        cassette_start = N_trim_start + aln_startpos
                        cassette_end = cassette_start + aln_len - 1
                        IB_seq = seq[:cassette_start-1]     # -1 because cassette_start is one-based
                        cassette_seq = seq[cassette_start-1:cassette_end]
                        flank_seq = seq[cassette_end:]
                        CIGAR_string = aln.original_sam_line.split('\t')[5]
                        edit_dist = aln.optional_field('NM')
                        ref_startpos = [c for c in aln.cigar if c.type=='M'][0].ref_iv.start + 1
                        ref_endpos = [c for c in aln.cigar if c.type=='M'][-1].ref_iv.end
                        ref_len = ref_endpos - ref_startpos + 1
                        OUTFILE_IB.write('>%s\n%s\n'%(name, IB_seq))
                        OUTFILE_flank.write('>%s\n%s\n'%(name, flank_seq))
                        OUTFILE_cassette.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(name, aln_len, len(IB_seq), cassette_seq, 
                              edit_dist, ref_startpos, ref_len, CIGAR_string, aln.optional_field('MD'), aln.optional_field('AS')))
                        results_total_readlen += aln_len
                        results_total_errors += edit_dist
                        results_IB_lengths[len(IB_seq)] += 1
                        results_cassette_read_lengths[length_ranges[aln_len]] += 1
                        results_reflengths[length_ranges[ref_len]] += 1
                        results_refstartpos[ref_startpos] += 1
    if not quiet:
        stripped_total = total_seqs - results_N_unstripped
        valp = general_utilities.value_and_percentages
        print "Total sequences:    %s"%total_seqs
        print "Cassette found in:  %s"%valp(stripped_total, [total_seqs])
        print "Not found in:       %s"%valp(results_N_unstripped, [total_seqs])
        print "IB length distribution (% of cassette-found seqs):"
        print '\n'.join("   %-15s  %s"%(length, valp(N, [stripped_total])) for (length, N) in sorted(results_IB_lengths.items()))
        print "Alignment length distribution - read, cassette (% of cassette-found seqs):"
        all_lengths = sorted(set(results_cassette_read_lengths.keys() + results_reflengths.keys()), key=str)
        print '\n'.join("   %-15s  %-15s %s"%(L, valp(results_cassette_read_lengths[L], [stripped_total]), 
                                              valp(results_reflengths[L], [stripped_total])) for L in all_lengths)
        print "First aligned cassette position distribution (% of cassette-found seqs):"
        print '\n'.join("   %-15s  %s"%(length, valp(N, [stripped_total])) for (length, N) in sorted(results_refstartpos.items()))
        print "Overall error rate: %.2g%%"%(results_total_errors/results_total_readlen*100)
    return (total_seqs, results_N_unstripped, results_IB_lengths, results_cassette_read_lengths, results_reflengths, 
            results_refstartpos, results_total_readlen, results_total_errors)


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    try:
        infile, outfile_base = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one infile and one outfile basename required!")
    allowed_min, allowed_max = [int(x) for x in options.allowed_IB_lengths.split(',')]
    allowed_IB_lengths = range(allowed_min, allowed_max+1)
    # the number of ignored bases has to be 1 lower than the shortest allowed IB, to distinguish it from too-short IBs
    N_ignore_start = allowed_min - 1
    if options.existing_bowtie_file:
        bowtie_outfile = options.existing_bowtie_file
        if not os.path.exists(bowtie_outfile):
            raise CassetteStrippingError("ERROR: provided bowtie file %s doesn't exit!"%bowtie_outfile)
    else:
        bowtie_outfile = outfile_base + "_aligned.sam"
        bowtie_output = run_bowtie(infile, bowtie_outfile, options.cassette_index, N_ignore_start, options.N_ignore_end, 
                                   options.n_threads, options.extra_bowtie_options, options.quiet)
        if not os.path.exists(bowtie_outfile):
            raise CassetteStrippingError("\n\nERROR: bowtie failed!\n%s  %s  %s"%bowtie_output)
    output_summary = parse_cassette_alignment(infile, bowtie_outfile, outfile_base, allowed_IB_lengths, options.min_cassette_length, 
                                              options.max_percent_errors, options.max_allowed_cassette_start, 
                                              N_ignore_start, options.N_ignore_end, options.quiet)


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    # tests in (testname, [test_description,] arg_and_infile_string) format
    # running the last test without -q to see the command-line output formatting
    test_runs = [
        ('strip__basic', 'basic cassette-stripping', '-C expected-cassette-end_CIB1-3p -b "-f" -q %s'%Testing.infile1), 
        ('strip__mism', 'cassette-stripping +mism', '-C expected-cassette-end_CIB1-3p -b "-f" -q %s'%Testing.infile2), 
        ('strip__indel', 'cassette-stripping +indel', '-C expected-cassette-end_CIB1-3p -b "-f" -e 100 %s'%Testing.infile3), 
                ]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, Testing.test_folder, 
                                argument_converter=argument_converter) 


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    test_folder = "test_data"
    infile1 = test_folder + '/INPUT_strip_cassette_basic.fa'
    infile2 = test_folder + '/INPUT_strip_cassette_mism.fa'
    infile3 = test_folder + '/INPUT_strip_cassette_indel.fa'
    infile4 = test_folder + '/INPUT_strip_cassette_extras.fa'
    # Note: INPUT_strip_cassette_indel.fa has multiple alignment cases (where there's a long indel in the middle, 
    #           so either perfect half of the match has a higher score than the whole with the big indel); 
    #       INPUT_strip_cassette_mism.fa DOESN'T, even though it has similar cases with a long mismatch - not sure why.

    def _run_test(self, infile, N_total, N_unstripped, IB_lengths, cassette_readlens, reflens, refstartpos, 
                  total_readlen, total_errors=0, allowed_IB_lengths=[21,22,23], min_cass_len=10, max_percent_error=30, max_start=5):
        bowtie_output = run_bowtie(infile, 'test_data/tmp.sam', 'expected-cassette-end_CIB1-3p', 20, 30, 1, "-f")
        output_summary = parse_cassette_alignment(infile, 'test_data/tmp.sam', 'test_data/tmp', 
                                                  allowed_IB_lengths, min_cass_len, max_percent_error, max_start, 20, 30, True)
        (total_seqs, results_N_unstripped, results_IB_lengths, results_cassette_read_lengths, results_reflengths,
         results_refstartpos, results_total_readlen, results_total_errors) = output_summary
        if reflens is None: 
            reflens = cassette_readlens
        if refstartpos is None:
            refstartpos = {1: N_total}
        self.assertEqual(total_seqs, N_total)
        self.assertEqual(results_N_unstripped, N_unstripped)
        self.assertEqual(results_IB_lengths, IB_lengths)
        self.assertEqual(results_cassette_read_lengths, cassette_readlens)
        self.assertEqual(results_reflengths, reflens)
        self.assertEqual(results_refstartpos, refstartpos)
        self.assertEqual(results_total_readlen, total_readlen)
        self.assertEqual(results_total_errors, total_errors)
        # only remove files if all the results were correct
        os.remove('test_data/tmp.sam')
        os.remove('test_data/tmp_flank.fa')
        os.remove('test_data/tmp_IB.fa')
        os.remove('test_data/tmp_no-cassette.fa')
        os.remove('test_data/tmp_cassette-info.txt')

    def test__everything_summary(self):
        # basic/truncated cases
        self._run_test(Testing.infile1, 
                       N_total = 11, N_unstripped = 4, 
                       IB_lengths = {22: 5, 21:1, 23:1}, 
                       cassette_readlens = {45: 3, '31-40':1, '21-30':1, '0-20':2}, reflens = None, 
                       refstartpos = {1: 7}, 
                       total_readlen = 235, total_errors = 0)
        # mismatches
        self._run_test(Testing.infile2, 
                       N_total = 9, N_unstripped = 0, 
                       IB_lengths = {22: 8, 23: 1}, 
                       cassette_readlens = {45: 5, 44: 1, 43: 1, '21-30': 1, '0-20': 1}, reflens = None, 
                       refstartpos = {1: 8, 2: 1}, 
                       total_readlen = 350, total_errors = 15)
        # try different IB lengths to make sure that's working right
        self._run_test(Testing.infile1, N_total = 11, N_unstripped = 9, IB_lengths = {23:1, 24:1}, 
                       cassette_readlens = {45: 2}, reflens = None, refstartpos = {1: 2}, total_readlen = 90, total_errors = 0, 
                       allowed_IB_lengths=[23,24])
        self._run_test(Testing.infile1, N_total = 11, N_unstripped = 11, IB_lengths = {}, 
                       cassette_readlens = {}, reflens = None, refstartpos = {}, total_readlen = 0, total_errors = 0, 
                       allowed_IB_lengths=[25,30])
        # try different min cassette lengths
        self._run_test(Testing.infile1, N_total = 11, N_unstripped = 8, IB_lengths = {21:1, 22:1, 23:1}, 
                       cassette_readlens = {45: 3}, reflens = None, refstartpos = {1: 3}, total_readlen = 135, total_errors = 0, 
                       min_cass_len=45)
        self._run_test(Testing.infile1, N_total = 11, N_unstripped = 7, IB_lengths = {21:1, 22:2, 23:1}, 
                       cassette_readlens = {45: 3, '31-40':1}, reflens = None, 
                       refstartpos = {1: 4}, total_readlen = 175, total_errors = 0, 
                       min_cass_len=40)
        # try different max_percent_error values
        self._run_test(Testing.infile2, N_total = 9, N_unstripped = 5, IB_lengths = {22: 3, 23: 1}, 
                       cassette_readlens = {44:1, 43:1, '21-30':1, '0-20':1}, reflens = None, 
                       refstartpos = {1: 3, 2: 1}, total_readlen = 125, total_errors = 0, 
                       max_percent_error=0)
        self._run_test(Testing.infile2, N_total = 9, N_unstripped = 2, IB_lengths = {22: 6, 23: 1}, 
                       cassette_readlens = {45: 3, 44: 1, 43: 1, '21-30': 1, '0-20': 1}, reflens = None, 
                       refstartpos = {1: 6, 2: 1}, total_readlen = 260, total_errors = 3,
                       max_percent_error=3)
        # try different max_allowed_cassette_start values
        self._run_test(Testing.infile2, 
                       N_total = 9, N_unstripped = 1, IB_lengths = {22: 8}, 
                       cassette_readlens = {45: 5, 43: 1, '21-30': 1, '0-20': 1}, reflens = None, 
                       refstartpos = {1: 8}, total_readlen = 306, total_errors = 15, 
                       max_start=1)
        # not testing infile3 and infile4 - they have detailed tests, and all the summary stuff has already been tested here.


if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        #print("\n * unit-tests for the ______ module")
        #test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(______
        #unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise pass the arguments to the main function
    main(args, options)
