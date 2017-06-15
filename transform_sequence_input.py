#! /usr/bin/env python2.7

### Various help functions for sequence and text processing
### Weronika Patena, nov2008

import parse_fasta

debug = 0

### Transform input with a given function
# Input can be fasta or just sequence. 
# Output is a list of (header,seq) tuples for fasta, or lines for text.
# For fasta, find all the full sequences that need to be processed (by a function passed here 
# as an argument), process them, and keep the remaining content.
# Otherwise just transform each line.
def transform_sequence_input(full_input, transform_function):
    output = []
    # check whether the input is fasta-formatted or just sequence
    # full_input is a list of lines, so convert it to string or this won't work
    if str(full_input).count('>')>0:    fasta = True
    else:                               fasta = False
    if debug: print '\n\t### INPUT:\n%s\t### END_INPUT\n'%full_input
    if fasta:
        for (header,seq) in parse_fasta.parse_fasta(full_input):
            output.append((header,transform_function(seq)))
    else:
        for line in full_input:
            output.append(transform_function(line))
    if debug:   print '\n\t######### FINAL OUTPUT: #########'%output
    return output

### Just a basic swapcase function to test the transform_sequence_input function
def swap_case(line):
    return line.swapcase()

### If called directly, test everything
if __name__ == '__main__':
    input = read_input()
    # use nice printing for fasta (header,seq) tuples, or normal for text lines
    try:                
        for (header,seq) in transform_sequence_input(input,swap_case):  parse_fasta.print_seq(header,seq)
    except ValueError:   
        for line in transform_sequence_input(input,swap_case):          print line
