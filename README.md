This is a snapshot of all the programs written in the Jonikas lab that were used in the paper **_A Genome-Wide, Mapped Algal Mutant Library Enables High-Throughput Genetic Studies in a Photosynthetic Eukaryote_ by Xiaobo Li, Weronika Patena, Friedrich Fauser, Robert E. Jinkerson, Shai Saroussi, Nina Ivanova, Jacob M. Robertson, Rebecca Yue, Ru Zhang, Josep Vilarrasa-Blasi, Tyler M. Wittkopp, Silvia Ramundo, Sean R. Blum, Audrey Goh, Matthew Laudon, Paul A. Lefebvre, Arthur R. Grossman and Martin C. Jonikas (journal, year)**, as of the final paper submission date - the code in this repository will never be updated.

Updated versions of the code, along with full git logs, functional test cases, and other programs not included here, can be found in our non-snapshot repositories: http://github.com/Jonikas-Lab/deepseq-processing, http://github.com/Jonikas-Lab/mutant-pools, http://github.com/Jonikas-Lab/combinatorial-pooling, http://github.com/Jonikas-Lab/combinatorial-deconvolution, and http://github.com/Jonikas-Lab/basic-bioinf-utilities.

All code was written by Weronika Patena.  

## File information

Command-line utilities used in the main data processing pipeline.  Run them on the command-line with a single "--help" argument to see detailed description of functionality, input, output and options.
 *   `robotic_plate_transfer.py` - generate the robot instruction files for combinatorial pooling.
 *   `deepseq_strip_cassette.py` - remove the expected cassette sequence from the LEAP-Seq reads, allowing for truncations, mismatches and indels (`bowtie2` is used for the sequence alignment); output barcode and flanking sequences to separate files.
 *   `deepseq_alignment_wrapper.py` - align the processed reads against the genome and cassette, parse both results, categorize all reads as unaligned, cassette, genomic-unique or genomic-multiple.

Modules used for combinatorial pooling and deconvolution - imported by the command-line programs above, and/or used for custom analysis in the interactive python shell (should not be run directly).  See class/function docstrings for details.
 *   `binary_code_utilities.py` - utilities for reading, modifying and analyzing binary error-correcting codes used for combinatorial pooling.
 *   `deconvolution_utilities.py` - functions for combinatorial deconvolution (i.e. analyzing the deep-sequencing datasets from each super-pool to find the original plate/colony location corresponding to each detected flanking sequence), as well as some additional utilities for analyzing and annotating the resulting deconvolution tables.

Modules used for insertion position analysis - used for analysis in the interactive python shell (should not be run directly).  See class/function docstrings for details.
 *   `mutant_IB_RISCC_classes.py` - the main module describing the classes used for insertions, insertion positions, and insertion datasets derived from LEAP-Seq data, along with methods used for creating them based on sequence files, annotating and manipulating them.
 *   `mutant_simulations.py` - dealing with genome mappability, mutant simulations, and finding hot/cold spots.
 *   `mutant_utilities.py` - includes a function for categorizing pairs of insertions in one mutant by distance/orientation, and minor utilities for insertion analysis.
 *   `gff_examine_file.py` - functions for parsing and examining GFF files describing gene positions (mostly GFF3 files from Phytozome, with some very basic functions for GFF2 files from JGI) - used to determine gene positions for insertions.
 *   `parse_annotation_file.py` - functions for parsing batch genome annotation files from Phytozome, annotation term definitions, and ID conversions between genome versions.
 *   `mutant_plotting_utilities.py` - used to generate some of the plots and statistics used in the paper.
 *   `mutant_analysis_classes.py` - a legacy version of insertion analysis code, largely replaced by `mutant_IB_RISCC_classes.py`, but still imported by some files

General utility modules that were imported by the command-line programs above, and/or used for custom analysis in the interactive python shell (should not be run directly).  See class/function docstrings for details.  These were mostly not written specifically for this paper, and not all the functions were used in the analysis for the paper, but I didn't want to modify the files by removing unused functions.
 *   `deepseq_utilities.py` - parsing SAM-format file fields to get detailed alignment data; parsing multiple sequence files in parallel to join the barcodes, flanking sequences and distal genomic sequences derived from LEAP-Seq data.
 *   `basic_seq_utilities.py` - basic functions for dealing with fasta/fastq files, sorting chromosome names correctly, etc
 *   `general_utilities.py` - convenience functions for manipulating numbers and standard data structures, file reading/writing, running command-line processes, etc. 
 *   `statistics_utilities.py` - convenience functions for quick running of chi-square goodness-of-fit and independence tests (using `scipy.stats.chisquare`), and for false-discovery-rate adjustment (using the R `p.adjust` function through `rpy2`)
 *   `plotting_utilities.py` - matplotlib convenience functions, defining custom colormaps, etc
 *   `testing_utilities.py` - convenience functions used for running functional tests of the three main command-line programs at the beginnin of this section (this folder doesn't include the input/output files for running functional tests - see the full repository links at the top for that).

Minor modules that were not actually used in the analysis, but are imported by some of the modules that were, and thus are included here for completeness and to make it easier to get everything to run.
 *   `parse_fasta.py`
 *   `reverse_complement.py`
 *   `reverse.py`
 *   `complement.py`
 *   `transform_sequence_input.py`
 *   `read_input.py`

Note that most of this is continued in-progress work.  I tried to keep the files reasonably documented, cleaned up and tested, but some will be better than others.  Please let me know if you have any questions or find any errors.

## Dependencies

#### python 2.7 (http://www.python.org)

#### Python packages (imported in one or more of the files):

 *   `numpy` (http://www.numpy.org/)
 *   `matplotlib` (http://matplotlib.org/)
 *   `scipy` (and `scipy.stats`) (http://www.scipy.org/)
 *   `rpy2` (for a few statistical functions) (http://rpy.sourceforge.net)
 *   `Bio` (Biopython; mostly only minor functions) (http://biopython.org)
 *   `HTSeq` (for parsing sam files) (http://pypi.python.org/pypi/HTSeq)
 *   `BCBio.GFF` (for parsing gff3 files) (http://github.com/chapmanb/bcbb/tree/master/gff)
 *   `bitstring` (for combinatorial pooling codes) (http://pypi.python.org/pypi/bitstring)

#### Standalone command-line programs (called via subprocess or another method in one or more of the files):

 *   `bowtie` (used in `deepseq_alignment_wrapper.py`) (http://bowtie-bio.sourceforge.net/index.shtml)
 *   `bowtie2` (used in `deepseq_strip_cassette.py`) (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
 *   `R` (used via `rpy2` for a few statistical functions) (http://www.r-project.org/)

## License

All of this code is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License, version 3, can be found here: http://www.gnu.org/licenses/gpl-3.0.html
