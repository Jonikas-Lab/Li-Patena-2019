#! /usr/bin/env python2.7
"""
This program takes one DNA sequence as an argument, prints the complement.
Weronika Patena, nov2008
"""

import string

debug = 0

# The complement table, including all IUPAC codes (S,W and N map to themselves):
#        A = adenine            R = G A (purine)            B = G T C
#        C = cytosine           Y = T C (pyrimidine)        D = G A T
#        G = guanine            K = G T (keto)              H = A C T
#        T = thymine            M = A C (amino)             V = G C A
#        U = uracil             S = G C                     
#                               W = A T                     N = A G C T (any)

# define list indices for DNA/RNA
DNA,RNA = 0,1

### Setup - defining complement tables
def complement_table_setup():
    # make empty tables for RNA and DNA complement info, with DNA/RNA 0/1 indices
    original_bases,complement_bases,complement_table = ['',''],['',''],['','']
    # define complement strings
    original_bases[DNA] =   'ATGCRYMKBVDHWSN'
    complement_bases[DNA] = 'TACGYRKMVBHDWSN'
    original_bases[RNA] =   'AUGCRYMKBVDHWSN'
    complement_bases[RNA] = 'UACGYRKMVBHDWSN'
    for base_table in original_bases,complement_bases:
        for seq_type in DNA,RNA:
            # add to complement strings: lowercase versions of bases
            base_table[seq_type] += base_table[seq_type].lower()
            # all other characters just get preserved, which is good. 
    # define actual complement translation tables
    for seq_type in DNA,RNA:
        complement_table[seq_type] = string.maketrans(original_bases[seq_type],
                                                      complement_bases[seq_type])
    return complement_table

# When called/imported, always set up the complement table
complement_table = complement_table_setup()

### The actual complement function
def complement(input_sequence,input=''):
    # if input type isn't given, detect RNA/DNA sequence type (default to DNA)
    if input=='DNA':                    input_type = DNA
    elif input=='RNA':                  input_type = RNA
    else:
        if input_sequence.count('U')>0:     input_type = RNA
        elif input_sequence.count('u')>0:   input_type = RNA
        else:                               input_type = DNA
    complement_seq = input_sequence.translate(complement_table[input_type])
    return complement_seq


### If called directly, read, parse and complement the input.
if __name__ == '__main__':
    import read_input,transform_sequence_input,parse_fasta
    input = read_input.read_input()
    for line in transform_sequence_input.transform_sequence_input(input,complement):
        parse_fasta.print_seq(line)
