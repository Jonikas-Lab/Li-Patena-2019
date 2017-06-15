#!/usr/bin/env python2.7

# standard library
from __future__ import division
import sys
import unittest
import os
from collections import defaultdict

# Gene annotation file data for v5.5 genome: 
#   list of (filename, content_headers, ID_column, content_columns, field_splitter, if_join_all_later_fields) tuples:
DEFAULT_GENE_ANNOTATION_FILES_v5p5 = [
    ('Creinhardtii_281_v5.5.geneName.txt', ['gene_name'], 0, [1], '\t', False),
    ('Creinhardtii_281_v5.5.defline.txt', ['defline'], 0, [1], '\t', False),
    ('Creinhardtii_281_v5.5.description.txt', ['description'], 0, [1], '\t', False),
    ('Creinhardtii_281_v5.5.synonym.txt', ['synonyms'], 0, [1], '\t', True),
    ('Creinhardtii_281_v5.5.annotation_info.txt', 'PFAM Panther KOG KEGG_ec KEGG_Orthology Gene_Ontology_terms best_arabidopsis_TAIR10_hit_name best_arabidopsis_TAIR10_hit_symbol best_arabidopsis_TAIR10_hit_defline'.split(), 
     1, [4, 5, 6, 7, 8, 9, 10, 11, 12], '\t', False),
]
DEFAULT_GENE_ANNOTATION_FOLDER = os.path.expanduser('~/experiments/reference_data/chlamy_annotation')
DEFAULT_GENE_ANNOTATION_FILES_v5p5 = [(os.path.join(DEFAULT_GENE_ANNOTATION_FOLDER, f), h, i, c, s, j) 
                                      for (f, h, i, c, s, j) in DEFAULT_GENE_ANNOTATION_FILES_v5p5]
DEFAULT_ID_CONVERSION_FILE_v5p5 = os.path.join(DEFAULT_GENE_ANNOTATION_FOLDER, 
                                               'ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt')

DEFAULT_ANNOTATION_DEFINITION_FILES_v5p5 = {
    'PFAM': ['Pfam-A.clans.tsv'], 
    'Panther': ['PANTHER9.0_HMM_classifications', 'PANTHER6.1_HMM_classifications', 'PANTHER10.0_HMM_classifications'], 
    'KOG': ['KOG.tsv'], 
    'KEGG_ec': ['KEGG_ec.tsv'], 
    'KEGG_Orthology': ['KEGG_Orthology.tsv'], 
    'Gene_Ontology_terms': ['GO.txt']
}
DEFAULT_ANNOTATION_DEF_FOLDER = os.path.join(DEFAULT_GENE_ANNOTATION_FOLDER, 'annotation_definitions')
DEFAULT_ANNOTATION_DEFINITION_FILES_v5p5 = {term: [os.path.join(DEFAULT_ANNOTATION_DEF_FOLDER, f) for f in files] 
                                            for (term, files) in DEFAULT_ANNOTATION_DEFINITION_FILES_v5p5.items()}

# LATER-TODO for synonyms files, strip the "t.*" part from the synonym ID columns as well as the gene ID column?

# Older deprecated annotation info
HEADER_FIELDS_v4 = ['Phytozome transcript name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']

HEADER_FIELDS_v5 = ['Phytozome internal transcript ID', 
                    'Phytozome gene locus name', 'Phytozome transcript name', 'Phytozome protein name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 'Gene Ontology terms', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']


def parse_gene_annotation_file(gene_annotation_filename, content_header_fields, gene_ID_column=0, content_columns=[1], 
                               if_join_all_later_fields=False, pad_with_empty_fields=True, field_splitter=None, 
                               strip_gene_fields_start=".t", genes_start_with=None, ignore_comments=False, verbosity_level=1):
    """ Parse tab-separated gene annotation file; return gene:annotation_list dictionary.

    Use column gene_ID_column to determine gene IDs; optionally shorten the gene name by truncating it starting with the 
     strip_gene_fields_start value if found (if not None) - e.g. if value is '.t', Cre01.g123450.t2.1 would become Cre01.g123450.
    If genes_start_with is not None, make sure all gene IDs start with it.

    If field_splitter is a string (probably '\t'), split fields by that character; if None, split by any number of whitespace chars.
        (Warning: that will treat multiple tabs as one splitter and thus skip empty fields!)

    If pad_with_empty_fields is True, pad shorter lines to the max length. 
    If ignore_comments is True, skip lines starting with #.

    Use content_columns for the values, or if if_join_all_later_fields is True, then use a comma-separated list of ALL
     fields starting with the content_columns one.

    Print some info/warnings to stdout depending on verbosity_level (0 - nothing, 1 - some, 2 - max).
    """
    if not os.path.lexists(gene_annotation_filename):
        raise Exception("Couldn't find the %s gene annotation file!"%gene_annotation_filename)
    if verbosity_level>0:
        print "  Parsing file %s for gene annotation info..."%os.path.basename(gene_annotation_filename)

    ### Parse the whole file into lists of tab-separated fields
    #    (could change to a generator, but most of the data has to stay in memory anyway in a different format, so probably no point)
    #    (and we need special treatment for the first line which may or may not be a header line...)
    data_by_row = []
    for line in open(gene_annotation_filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip().split(field_splitter)
        data_by_row.append(fields)
    if verbosity_level>0:
        print "  Parsed %s lines"%len(data_by_row)
        
    # if any of the other lines doesn't start with a Cre* gene ID, fail!
    if genes_start_with is not None:
        for row in data_by_row:
            if not row[gene_ID_column].startswith(genes_start_with):
                raise Exception("Can't parse file %s - found line that doesn't start "%gene_annotation_filename
                                +"with a %s gene ID!\n  \"%s\""%(genes_start_with, '\t'.join(row)))

    # check that all the data lengths line up (if they don't, don't throw an error
    data_lengths = set([len(row) for row in data_by_row])
    if len(data_lengths)>1 and not if_join_all_later_fields:
        mismatched_lengths = True
        if verbosity_level:     print "Not all data rows have the same length! Lengths found: %s"%list(data_lengths)
        if pad_with_empty_fields:
            max_length = max(data_lengths)
            if verbosity_level>0:
                print "Data field numbers vary between rows - padding all lower-length data rows to length %s"%max_length
            for row in data_by_row:
                if len(row)<max_length:
                    row += ['' for x in range(max_length-len(row))]
            mismatched_lengths = False
    else:   mismatched_lengths = False

    # LATER-TODO figure out the total #fields over the whole file, and then: 
    #   - if content_header_fields is None, just use the filename and N empties
    #   - if content_columns is None, use all except the ID column, I guess?
    if len(content_header_fields) != len(content_columns):
        raise Exception("Error: content_header_fields has a different number of fields than content_columns!")
    if max(data_lengths) < max(content_columns):
        raise Exception("content_columns specifies a column that doesn't exist! File %s; "%gene_annotation_filename
            +"lines have %s-%s fields, but content_columns wants %s"%(min(data_lengths), max(data_lengths), content_columns))

    # MAYBE-TODO remove empty columns (only if all the data lengths match!)

    ### Convert the list-format data into a by-gene dictionary, grabbing only the fields we want
    # we frequently get multiple lines, for different transcripts!  Just concatenate all of them.
    data_by_gene = defaultdict(lambda: [set() for _ in content_columns])
    for data in data_by_row:
        gene = data[gene_ID_column]
        if strip_gene_fields_start is not None:
            gene = gene.split(strip_gene_fields_start)[0]
        if if_join_all_later_fields:
            if len(content_columns) > 1:
                raise Exception("if_join_all_later_fields not implemented with multple content_columns!")
            data_by_gene[gene][0].add(','.join(str(x) for x in data[content_columns[0]:] if x.strip()))
        else:
            for (new_col, old_col) in enumerate(content_columns):
                data_by_gene[gene][new_col].add(data[old_col])

    # At the end, change the sets to strings, and also remove empty strings
    for gene,data in data_by_gene.items():
        data = [', '.join([f for f in fields if f.strip()]) for fields in data]
        data_by_gene[gene] = [x if x else '-' for x in data]

    if verbosity_level>0:
        print "  DONE Parsing gene annotation file - found %s genes"%len(data_by_gene)
    return defaultdict(lambda: ['-' for _ in content_columns], data_by_gene)


def parse_ID_conversion_file(infile=DEFAULT_ID_CONVERSION_FILE_v5p5, key_header='5.5', val_header=None, 
                             strip_fields_after = ['.t', '_t'], strip_fields_before = ['|', 'au5.', 'Au9.'], 
                             blank_val = '--', warn_on_multiples=True):
    """ Return key_ID:other_ID_list dict and other_ID_header based on infile, or just key_ID:val_ID if val_header is given.

    key_ID is taken from key_header; all the other columns go in the value list, 
        UNLESS val_header is given, in which case the output will just be key:val with the other columns ignored.
    From all field values, all characters after each strip_fields_after val and before each strip_fields_before val are removed:
        so, Cre01.g000550.t1.2 becomes Cre01.g000550, jgi|Chlre4|391833 becomes 391833, au5.g987_t1 becomes g987.
    This means that sometimes there are multiple lines per stripped ID - those lines are merged, excluding repeats and blank_val,
        and the results are comma-separated if there are multiple values; 
        if warn_on_multiples is True, a warning is printed for all multiple-value cases.
    """
    # Note that the file is variable-space-separated, not tab-separated!
    conversion_dict = {}
    header = []
    for line in open(infile):
        if line.startswith('#') and '  ' not in line:   continue        # skip title line
        if line.startswith('#'):
            if header:  raise Exception("Found multiple header lines! (start with # and contain tabs)")
            column_headers = line[1:].strip().split()
            try:                key_column = column_headers.index(key_header)
            except ValueError:  raise Exception("key_header %s didn't appear in header!"%key_header)
            header = column_headers
            del header[key_column]
        if not line.startswith('#'):
            if not header:  raise Exception("Found data line before header line!")
            fields = line.strip().split()
            for sep in strip_fields_after:  fields = [x.split(sep)[0] for x in fields]
            for sep in strip_fields_before: fields = [x.split(sep)[-1] for x in fields]
            key = fields[key_column]
            if key == blank_val:    continue
            del fields[key_column]
            if key not in conversion_dict:  
                conversion_dict[key] = fields
            else:
                old_fields = [set(x.split(',')) for x in conversion_dict[key]]
                joint_fields = [','.join((set_o | set([n])) - set([blank_val])) for (set_o, n) in zip(old_fields, fields)]
                conversion_dict[key] = joint_fields
                if warn_on_multiples:
                    for field, curr_header in zip(joint_fields, header):
                        if ',' in field:
                            print "Warning: %s %s has multiple %s values! %s"%(key_header, key, curr_header, field)
    if val_header:
        try:                val_column = column_headers.index(val_header)
        except ValueError:  raise Exception("val_header %s didn't appear in header!"%val_header)
        return {key: vals[val_column] for (key, vals) in conversion_dict.items()}
    else:
        return conversion_dict, header
    # TODO add unit-tests!


def get_all_gene_annotation(genome_version=None, gene_annotation_files=None, ignore_comments=False, print_info=False):
    """ Grab all the annotation (depends on genome version); return gene:annotation_list dict and header list.

    Can provide a gene_annotation_files dict instead - in the same format as DEFAULT_GENE_ANNOTATION_FILES_v5p5 here.
    """
    if genome_version is None and gene_annotation_files is None:
        raise Exception("User has to provide genome_version or gene_annotation_files!")
    if genome_version is not None and gene_annotation_files is not None:
        raise Exception("User should not provide both genome_version and gene_annotation_files, just one!")
    elif gene_annotation_files is None:
        if genome_version == 5.5:
            gene_annotation_files = DEFAULT_GENE_ANNOTATION_FILES_v5p5
        else:
            raise Exception("Genome version %s not implemented right now!"%genome_version)
    # MAYBE-TODO add options for strip_gene_fields_start etc
    gene_annotation_dicts = [parse_gene_annotation_file(filename, content_header_fields, ID_column, content_columns, 
                                                        if_join_all_later_fields, pad_with_empty_fields=True, 
                                                        field_splitter=splitter, strip_gene_fields_start=".t", genes_start_with=None,
                                                        ignore_comments=ignore_comments, verbosity_level=print_info) 
                             for (filename, content_header_fields, ID_column, content_columns, splitter, if_join_all_later_fields) 
                             in gene_annotation_files]
    full_header = sum([content_headers for (_, content_headers, _, _, _, _) in gene_annotation_files], [])

    all_gene_IDs = set.union(*[set(d.keys()) for d in gene_annotation_dicts])
    full_annotation_dict = defaultdict(lambda: ['-' for _ in full_header])
    for gene in all_gene_IDs:
        full_annotation_dict[gene] = sum([d[gene] for d in gene_annotation_dicts], [])
    return full_annotation_dict, full_header


### Bit of old code to get gene names from gff3 file (no longer necessary with v5.5 genome, which has a tab-sep geneName file):
#   new_fields.append('transcript_names')
#   genename_dict = defaultdict(set)
#   for line in open(gff_file_for_gene_names):
#       if line.startswith('#'):    continue
#       all_fields = line.strip().split('\t')
#       if all_fields[2] != 'mRNA': continue
#       data = all_fields[8]
#       fields = dict(x.split('=') for x in data.split(';'))
#       gene = fields['Name'].split('.t')[0]
#       try:                genename_dict[gene].add(fields['geneName'])
#       except KeyError:    pass
#   genename_dict = defaultdict(set, {key:','.join(vals) for (key,vals) in genename_dict.items()})


def get_term_definitions_from_obo(obo_infile):
    """ Parse OBO file; for each [Term] section, yield name:val dict with names being id, name, namespace, def, extras.
    """
    curr_record = {}
    in_term = False
    for line in open(obo_infile):
        if line.startswith('[Term]'):   in_term = True
        if in_term: 
            for field in 'id name namespace def'.split():
                if line.startswith(field + ': '):   curr_record[field] = line.strip()[len(field + ': '):]
            if not line.strip():
                if 'id' not in curr_record: raise Exception("Record lacks ID! %s"%curr_record)
                if len(curr_record)<2:      raise Exception("Record only has ID! %s"%curr_record['id'])
                raw_def = curr_record['def']
                try:                raw_def, raw_extras = raw_def.rsplit('" ', 1)
                except ValueError:  raise Exception("Can't parse this into def and extras! %s"%raw_def)
                curr_record['def'] = raw_def.strip().strip('"').replace('\\"', '"')
                curr_record['extras'] = raw_extras.strip().strip('[]')
                yield curr_record
                curr_record = {}
                in_term = False


def simple_format_GO_file(infile, outfile):
    """ Read obo format GO file and write term definitions to outfile.

    Was used on the files in ~/experiments/reference_data/chlamy_annotation/annotation_definitions
    """
    columns = 'id name namespace def extras'.split()
    with open(outfile, 'w') as OUTFILE:
        OUTFILE.write('\t'.join(columns) + '\n')
        for record in get_term_definitions_from_obo(infile):
            OUTFILE.write('\t'.join(record.get(x,'-') for x in columns) + '\n')


def GO_relationship_tree(infile):
    """ Read obo format GO file, grab the relationships (X is a subset of Y) to get the full tree.

    Was used on go-basic.obo file in ~/experiments/reference_data/chlamy_annotation/annotation_definitions.

    Generated using http://pythonhosted.org/Orange-Bioinformatics/reference/ontology.html
    Other resources on GO parsing: http://blog.nextgenetics.net/?e=6, https://pypi.python.org/pypi/goatools
    """
    import Orange.bio.ontology
    return Orange.bio.ontology.OBOOntology(infile)


def get_all_annotation_definitions(annotation_types_files=DEFAULT_ANNOTATION_DEFINITION_FILES_v5p5):
    """ Return annotation_type:(term:definition) double dictionary based on types and files in argument. """
    parsing_functions = {'PFAM': get_PFAM_definitions, 'Panther': get_Panther_definitions,  'KOG': get_KOG_definitions, 
                         'KEGG_ec': get_KEGGec_definitions, 'KEGG_Orthology': get_KEGGorthology_definitions, 
                         'Gene_Ontology_terms': get_GO_definitions }
    all_annotation = {}
    for annotation_type, files in annotation_types_files.items():
        # parse all files
        try:
            parsing_function = parsing_functions[annotation_type]
        except KeyError:
            raise Exception("No defined parsing function for %s! There are functions for %s"%(annotation_type, 
                                                                                              ', '.join(parsing_functions)))
        dictionaries = [parsing_function(f) for f in files]
        # if there are multiple files, use the later ones to fill in only what was absent in the first
        final_dict = dictionaries[0]
        if len(files) > 1:
            print "For %s data, merging multiple files - starting with %s (%s terms)"%(annotation_type, os.path.split(files[0])[1], 
                                                                                       len(final_dict))
        for (extra_dict, extra_file) in zip(dictionaries[1:], files[1:]):
            N_added = 0
            for term, definition in extra_dict.items():
                if term not in final_dict:
                    N_added += 1
                    final_dict[term] = definition
            print "  - added %s/%s terms from %s - total %s terms."%(N_added, len(extra_dict), 
                                                                     os.path.split(extra_file)[1], len(final_dict))
        all_annotation[annotation_type] = final_dict
    return all_annotation


def get_PFAM_definitions(infile):
    """ Return term:definition dict based on PFAM term definition file. """
    definitions = {}
    for line in open(infile):
        fields = line[:-1].split('\t')
        ID = fields[0]
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        definitions[ID] = fields[4] + (" (%s)"%fields[2] if fields[2] else '')
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


def get_Panther_definitions(infile):
    """ Return term:definition dict based on Panther term definition file. """
    definitions = {}
    for line in open(infile):
        fields = line[:-1].split('\t')
        ID = fields[0]
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        definitions[ID] = "%s (%s)"%(fields[1], ', '.join(fields[2:]))
        # MAYBE-TODO include definitions for all the GO terms listed in the fields, too?  But there's LOTS
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


def get_KOG_definitions(infile):
    """ Return term:definition dict based on KOG term definition file. """
    definitions = {}
    for line in open(infile):
        if line.startswith('Identifier'):   continue
        fields = line[:-1].split('\t')
        ID = fields[0]
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        definitions[ID] = "%s (namespace %s)"%(fields[1], fields[2])
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


def get_KEGGec_definitions(infile):
    """ Return term:definition dict based on KEGG-ec term definition file. """
    definitions = {}
    for line in open(infile):
        if line.startswith('Identifier'):   continue
        fields = line[:-1].split('\t')
        ID = fields[0]
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        if fields[2] in ['', '""']:     definitions[ID] = fields[1]
        else:                           definitions[ID] = "%s (%s)"%(fields[2], fields[1])
    # Some definitions are just "Transferred entry: 1.21.3.4" - in those cases grab the real ones!
    for ID, definition in definitions.items():
        if definition.startswith('Transferred entry:'):
            new_IDs = definition.split(': ')[1].replace(' and ', ', ').split(', ')
            new_definition = "Transferred entries (%s): %s"%(len(new_IDs), ' & '.join(definitions[x] for x in new_IDs))
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


def get_KEGGorthology_definitions(infile):
    """ Return term:definition dict based on KEGG-Orthology term definition file. """
    definitions = {}
    for line in open(infile):
        if line.startswith('Namespace'):                continue
        fields = line[:-1].split('\t')
        ID = fields[1]
        if ID[0] != 'K' or ID[1] not in '1234567890':   continue
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        definitions[ID] = "%s (%s)"%(fields[2], fields[3])
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


def get_GO_definitions(infile):
    """ Return term:definition dict based on GO term definition file. """
    definitions = {}
    for line in open(infile):
        if line.startswith('id\tname'):                continue
        fields = line[:-1].split('\t')
        ID = fields[0]
        if ID in definitions:
            raise Exception("ID %s shows up twice in file %s!"%(ID, infile))
        definitions[ID] = "%s (%s: %s) [%s]"%(fields[1], fields[2], fields[3], fields[4])
    print "Parsed %s file - %s term definitions."%(os.path.split(infile)[1], len(definitions))
    return definitions


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__all(self):
        # do it twice just to make sure the results are the same - earlier I had a bug where they weren't!
        for i in (1,2,3,4):
            if i < 3:
                data, header = get_all_gene_annotation(5.5, print_info=False)
            else:
                data, header = get_all_gene_annotation(gene_annotation_files=DEFAULT_GENE_ANNOTATION_FILES_v5p5, print_info=False)
            self.assertEquals(header, 'gene_name defline description synonyms PFAM Panther KOG KEGG_ec KEGG_Orthology Gene_Ontology_terms best_arabidopsis_TAIR10_hit_name best_arabidopsis_TAIR10_hit_symbol best_arabidopsis_TAIR10_hit_defline'.split())
            self.assertEquals(len(data), 17741)
            # a gene that doesn't exist
            self.assertEquals(data['Cre01.g000000'], ['-' for x in header])
            # a gene that exists
            self.assertEquals(data['Cre01.g000150'], ['ZRT2', 'Zinc-nutrition responsive permease transporter', 'Zinc/iron permease; related to plant homologs; ZIP family, subfamily I; appears as CrZIP2 in PMID: 15710683; PMID: 16766055', 'g6.t1,Cre01.g000150.t1.1', 'PF02535', 'PTHR11040,PTHR11040:SF30', 'KOG1558', '-', '-', 'GO:0016020,GO:0030001,GO:0046873,GO:0055085', 'AT2G04032.1', 'ZIP7', 'zinc transporter 7 precursor'])
            # a gene that exists and has multiple annotation/synonym lines
            self.assertEquals(data['Cre17.g733850'], 
                              ['-', '-', '-', 'g17740.t1,Cre17.g733850.t1.1, g17740.t2', 'PF01391'] + ['-' for x in range(8)])
    # MAYBE-TODO add a v4 unit-test too, and for other formats?  Once I have them implemented.


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
