#! /usr/bin/env python2.7
"""
This program takes one DNA sequence as an argument, prints the reverse-complement.
Weronika Patena, nov2008
"""

import complement

### The actual reverse-complement function:
def reverse_complement(input_sequence,input_type=''):
    # use the complement function from the complement module
    complement_seq = complement.complement(input_sequence,input_type)
    # easy reverse method: full list slice with a step of -1
    reverse_complement_seq = complement_seq[::-1]
    return reverse_complement_seq

### If called directly:
# try reading the argument; if none given, read from stdin (which is what makes it work with pipes (echo ctcgag | script.py) and when called on a selection in vi and such). 
if __name__ == '__main__':
    import read_input,transform_sequence_input,parse_fasta
    # check if input is a list of files (i.e. if first arg is a valid filename - good enough)
    input = read_input.read_input()
    # transform_sequence_input takes an input AND the function to apply to it - in this case rev-compl
    for line in transform_sequence_input.transform_sequence_input(input,reverse_complement):
        parse_fasta.print_seq(line)
