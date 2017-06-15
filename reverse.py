#! /usr/bin/env python2.7
"""
This program takes one DNA sequence as an argument, prints the reverse.
Weronika Patena, nov2008
"""

### The actual reverse function:
def reverse(input_seq):
    reverse_seq = input_seq[::-1]
    return reverse_seq

### If called directly:
# try reading the argument; if none given, read from stdin (which is what makes it work with pipes (echo ctcgag | script.py) and when called on a selection in vi and such). 
if __name__ == '__main__':
    import read_input,transform_sequence_input,parse_fasta
    # check if input is a list of files (i.e. if first arg is a valid filename - good enough))
    input = read_input.read_input()
    # transform_sequence_input takes an input AND the function to apply to it - in this case reverse()
    for line in transform_sequence_input.transform_sequence_input(input,reverse):
        parse_fasta.print_seq(line)
