#! /usr/bin/env python2.7

""" Various basic deepseq-related utilities I wrote - see separate function docstrings.  For importing only.
 --Weronika Patena, 2011
"""

# basic libraries
import unittest
# other packages
import HTSeq
# my modules
import basic_seq_utilities


## NOTE: for fasta/fastq (raw data) utilities see basic_seq_utilities.py


################## aligned data utilities ################## 

### NOTE: I'm currently using HTSeq for parsing SAM-format data, but I could try something else - pysam seems good.

# SAM and GFF are 1-based, with the end being the last base; HTSeq is 0-based, with the end being the base AFTER last.

######### NOTES ON THE SAM FORMAT
### Header:
# MAYBE-TODO do something useful with the SAM header?  Or at least copy it to outfile?
### Alignment line fields:
# * query template name
# * bitwise flag (relevant bits: 4 = unmapped, 16 = reverse-complement, 512 = failed quality control)
# * reference sequence name (* = unmapped)
# * leftmost mapping position (1-based) (0 = unmapped)
# * mapping quality ("-10 * log10(probability that position is wrong)"; 255 = unknown)
# * CIGAR string - descriptions of alignment matches/mismatches/etc (M/= match, I/D ins/del, X mismatch, S/H clipping)
# * (PE only - reference name of the mate fragment)
# * (PE only - position of the mate fragment)
# * template length
# * fragment sequence
# * ASCII of Phred-scaled base quality + 33   (original deepseq read quality)
# * OPTIONAL FIELDS, lots of different possibilities:   MD is mismatch info string, NM is edit distance to reference
#       (for info on MD field format see SAM manual footnote, and 
#        sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info)
#########

CIGAR_TYPES_MATCH = ['=']
CIGAR_TYPES_NOOP = ['S','H','P']
CIGAR_TYPES_MUTATION = ['X','I','D']
CIGAR_TYPES_INTRON = ['N']     # 'N' is for introns, but we shouldn't be paying attention to those for genomic DNA seq
CIGAR_TYPES_UNKNOWN = ['M']
# MAYBE-TODO HTSeq doesn't appear aware of the = and X operations...  http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html#HTSeq.CigarOperation  - I emailed the author about it, no response

class DeepseqError(Exception):
    """ Exception in this module; no special behavior."""
    pass

### Parsing two fastq files in parallel (for paired-end deepseq data)

def parse_2fastq_parallel(file1, file2):
    """ Parse two fastq files in parallel - generator yielding (name, seq1, seq2, qual1, qual2) tuples.

    Doesn't check that the readnames match.
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator    # Bio is the biopython package
    with open(file1) as INFILE1:
        with open(file2) as INFILE2:
            generator1 = FastqGeneralIterator(INFILE1)
            generator2 = FastqGeneralIterator(INFILE2)
            if_finished_1, if_finished_2 = False, False
            while True:
                try:                    name1, seq1, qual1 = generator1.next()
                except StopIteration:   if_finished_1 = True
                try:                    name2, seq2, qual2 = generator2.next()
                except StopIteration:   if_finished_2 = True
                name = name1.split()[0]
                if not if_finished_1 and not if_finished_2:
                    yield (name, seq1, seq2, qual1, qual2)
                elif if_finished_1 and if_finished_2:
                    raise StopIteration
                else:
                    raise DeepseqError("One file finished but the other one didn't! Read name %s"%(
                                                                        name if if_finished_2 else name2.split()[0]))
    # TODO unit-test!

def parse_fastx_sam_parallel(fastx_infile, sam_infile):
    """ Parse fastx and resulting sam file in parallel - generator yielding (name, seq, alignment_list) tuples.

    The sam file may contain multiple alignments per read.  Program checks that the readnames match.
    """
    fastx_generator = basic_seq_utilities.name_seq_generator_from_fasta_fastq(fastx_infile)
    sam_generator = iter(HTSeq.bundle_multiple_alignments(HTSeq.SAM_Reader(sam_infile)))
    if_finished_fastx, if_finished_sam = False, False
    while True:
        try:                    name, seq = fastx_generator.next()
        except StopIteration:   if_finished_fastx = True
        try:                    alns = sam_generator.next()
        except StopIteration:   if_finished_sam = True
        # if both finished, good, we're doine
        if if_finished_fastx and if_finished_sam:
            raise StopIteration
        # if one file was finished but the other wasn't, error!
        elif if_finished_fastx or if_finished_sam:
            raise DeepseqError("Parsing seq/aln files in parallel - inconsistent finished states! "
                              +"(If finished: %s %s, %s %s)"%(fastx_infile, if_finished_fastx, sam_infile, if_finished_sam))
        # if all the files still contained data, yield it
        else:
            name = name.split()[0]
            name2 = alns[0].read.name.split()[0]
            if not name2 == name:
                raise DeepseqError("Non-matching readnames between files! %s in %s, %s in %s"%(fastx_infile, name, 
                                                                                               sam_infile, name2))
            yield (name, seq, alns)


### Getting mutation counts from various SAM alignment format fields, as read by HTSeq

def _get_HTSeq_optional_field_either_version(val_or_tuple):
    """ Different HTSeq versions return either val or (name,val) from aln.optional_field(name) - convert either to val. """
    if isinstance(val_or_tuple, tuple): return val_or_tuple[1]
    else:                               return val_or_tuple

def get_HTSeq_optional_field(HTSeq_alignment, field_name):
    """ Return value of optional field (like NM, XM, etc). """
    return _get_HTSeq_optional_field_either_version(HTSeq_alignment.optional_field(field_name))

def check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment, based on CIGAR string; -1 if unknown ('M') by default.
    If treat_unknown_as is 'unknown', return -1 whenever an unknown (M, may be match or mismatch) operation is found; 
     if treat_unknown_as is 'mutation' or 'match', count unknowns accordingly.  Return -1 if read is unaligned.
    If ignore_introns is False, count introns (N) as mutations; otherwise don't."""
    global CIGAR_TYPES_MUTATION, CIGAR_TYPES_INTRON, CIGAR_TYPES_UNKNOWN
    # just return -1 for unaligned reads
    if HTSeq_alignment.cigar is None:
        return -1
    # figure out whether to consider intron-skipping ('N') as a mutation or not, based on argument
    if ignore_introns:
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_mutation = [] + CIGAR_TYPES_MUTATION
    else:
        cigar_types_mutation = CIGAR_TYPES_MUTATION + CIGAR_TYPES_INTRON   
    # figure out how to treat unknown matches ('M'), based on argument
    if treat_unknown_as=='unknown':
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_unknown = [] + CIGAR_TYPES_UNKNOWN
    elif treat_unknown_as=='mutation':
        cigar_types_mutation += CIGAR_TYPES_UNKNOWN
        cigar_types_unknown = []
    elif treat_unknown_as=='match':
        cigar_types_unknown = []
    else:
        raise ValueError("treat_unknown_as argument value must be 'mutation', 'match' or 'unknown'")
    # count the mutations, return total count (or instantly return -1 on finding an unknonw)
    mutations = 0
    for cigar_op in HTSeq_alignment.cigar:
        if cigar_op.type in cigar_types_mutation:
            mutations += cigar_op.size
        # if there's an unknown, just return -1, no need to count
        elif cigar_op.type in cigar_types_unknown:
            return -1
    return mutations


def check_mutation_count_by_optional_NM_field(HTSeq_alignment, negative_if_absent=True):
    """ Return #errors in HTSeq_alignment, based on optional NM field; -1 or exception if field missing."""
    # for unalign reads NM field is missing - returns -1
    try:                        return get_HTSeq_optional_field(HTSeq_alignment, 'NM')
    except KeyError:    
        if negative_if_absent:  return -1
        else:                   raise DeepseqError("Optional NM field missing in read %s - can't determine #errors!"%HTSeq_alignment.read.name)


def check_mutation_count_by_optional_MD_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional MD field; -1 if unknown (MD field missing)."""
    # for info on MD field format see SAM manual footnote, 
    #   and sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info
    #       basically a number means matches, a letter means a mismatch to reference (or insertion? is that different?), 
    #       letters preceded by ^ mean deletion from the reference
    try:                mutation_string = get_HTSeq_optional_field(HTSeq_alignment, 'MD')
    except KeyError:    return -1
    # for unalign reads MD field is missing - returns -1
    mutation_letters = [c for c in mutation_string if not (c.isdigit() or c=='^')]
    #   (^ is used in describing a mutation but it shouldn't be counted as a separate mutation - only letters count.)
    return len(mutation_letters)


def check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment (look at CIGAR string and NM and MD optional fields); -1 if unknown.
    First check the CIGAR string but only accept the answer if there are no unknown ('M') characters; 
     then check the NM and MD fields and return the result if those fields exist.
    If the CIGAR string is ambiguous and neither of the optional fields exist:
     - if treat_unknown_as is 'unknown', return -1
     - if treat_unknown_as is 'mutation' or 'match', return the CIGAR string result with unknowns counted accordingly.
    If ignore_introns is False, count introns (N) in CIGAR string as mutations; otherwise don't .
    Does NOT guarantee returning a sensible value if the CIGAR, NM and MD fields contain inconsistent information.
    """
    mutation_count = check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', 
                                                          ignore_introns=ignore_introns)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_NM_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_MD_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    if treat_unknown_as=='unknown':     
        return -1
    return check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as=treat_unknown_as, 
                                                          ignore_introns=ignore_introns)


################## unit tests ################## 

class Fake_deepseq_objects:
    """ Fake deepseq data objects for testing. """
    # NOTE: not all of those are used in the unit-tests for this module, but they're also imported elsewhere!

    class Fake_HTSeq_cigar_op:
        """ Fake CIGAR operation, mimicking HTSeq cigar object."""
        size = 1
        def __init__(self,string):  
            self.type = string

    class Fake_HTSeq_genomic_pos:
        """ Fake HTSeq.GenomicPosition. """
        def __init__(self, chrom, strand, start, end):
            self.chrom = chrom
            self.strand = strand
            self.start = start
            self.end = end

    class Fake_HTSeq_read:
        """ Fake read, as in HTSeq_alignment.read. """
        def __init__(self,seq='AAA',name='test'):
            self.seq = seq
            self.name = name

    class Fake_HTSeq_alignment:
        """ Fake HTSeq.Alignment object."""

        def __init__(self, seq='AAA', readname='test', unaligned=False, pos=('chr_','+',0,0), 
                     cigar_string=None, optional_field_data={}):    
            self.read = Fake_deepseq_objects.Fake_HTSeq_read(seq,readname)
            if unaligned:
                self.aligned = False
                self.iv = None
            else:
                self.aligned = True
                self.iv = Fake_deepseq_objects.Fake_HTSeq_genomic_pos(*pos)
            self.optional_field_data = optional_field_data
            if cigar_string is None:    self.cigar = None
            else:                       self.cigar = [Fake_deepseq_objects.Fake_HTSeq_cigar_op(c) for c in cigar_string]

        def optional_field(self,field):             
            return self.optional_field_data[field]


class Testing(unittest.TestCase):
    """ Unit-tests for all the functions/classes in this module. """

    def test__parse_fastx_sam_parallel(self):
        output = list(parse_fastx_sam_parallel('_test_inputs/test_parallel2.fq', '_test_inputs/test_parallel2.sam'))
        assert len(output) == 3
        assert [len(x) for x in output] == [3, 3, 3]
        assert [len(x[2]) for x in output] == [2, 1, 1]
        assert output[0][2][0].read.name == output[0][0] == 'ROCKFORD:4:1:1680:975#0/1'
        assert output[0][2][0].read.seq == output[0][1] ==  'ACTAATACGCGGCCTGGAGCTGGACGTTGGAACCAA'
        # the generator isn't really run until you ask for its results, so I have to run list on it to get the error
        self.assertRaises(DeepseqError, list, parse_fastx_sam_parallel('_test_inputs/test.fq', '_test_inputs/test_parallel2.sam'))

    def test__check_mutation_count_by_CIGAR_string(self):
        # no alignment (CIGAR is None)
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment()
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == -1
        # CIGAR is unambiguous, no MD or NM given (or needed)
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='==')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='XX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='DD')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='S=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='SX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 1
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='N=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=True) == 0
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=False) == 1
        # CIGAR is ambiguous (contains M's) - return -1, 2 or 0 depending on what treat_unknown_as is set to
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='MM')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='MX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 1

    def test__check_mutation_count_by_optional_NM_field(self):
        """ the tested function should return -1 if no NM field, otherwise return value of NM field. """
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment()
        assert check_mutation_count_by_optional_NM_field(fake_alignment, negative_if_absent=True) == -1
        self.assertRaises(DeepseqError, check_mutation_count_by_optional_NM_field, fake_alignment, negative_if_absent=False)
        for x in range(10):
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'NM':x})
            assert check_mutation_count_by_optional_NM_field(fake_alignment) == x

    def test__check_mutation_count_by_optional_MD_field(self):
        """ see ~/experiments/reference_data/aligner_format_info/* files for MD field examples."""
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment({})
        assert check_mutation_count_by_optional_MD_field(fake_alignment) == -1
        for s in [str(x) for x in range(30)]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A0G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A2G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A2G2T2C2N'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 5
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'^A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'^AGC'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 3

    def test__check_mutation_count_try_all_methods(self):
        """ The order of check is CIGAR, NM, MD; CIGAR is skipped if ambiguous; NM and MD skipped if inexistent. 
        Not attempting to deal with inconsistent states sensibly."""
        # all measures agree there are no mutations (with 0-2 of NM/MD fields present)
        for opt_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}, {}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='='*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # all measures agree there is a mutation (with 0-2 of NM/MD fields present)
        for opt_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}, {}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='X'+'='*9,optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 1
        # CIGAR is ambiguous, there are no mutations according to NM/MD (NM, MD or both are present)
        for opt_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M'*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # CIGAR is ambiguous, there is a  mutation according to NM/MD (NM, MD or both are present)
        for opt_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M'*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 1


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main()

