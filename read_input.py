#! /usr/bin/env python2.7
"""
Read input using various methods (from raw or file arguments or stdin)
"""

import sys,os

debug = 0

### Read the input:
# try reading the argument; if none given, read from stdin (which is what makes it work 
# with pipes (echo ctcgag | script.py) and when called on a selection in vi and such). 
def read_input():
    # if arguments were given:
    try:
        arguments = sys.argv[1:]
        # check if they're a list of files (i.e. if first arg is a valid filename - good enough)
        # if yes, print the outputs one by one, with file names.
        # TODO this is an issue with large files, since it tries to read the whole thing into memory.
        # TODO maybe I could just return the open file (or a fileinput thing) instead of the list of lines?
        if os.path.lexists(arguments[0]):
            if debug:   print "Processing input as a list of files..."
            if debug:   print "\t%s"%'\n\t'.join(arguments)
            input = ''
            for infile in arguments:
                INFILE = open(infile,'r')
                input += INFILE.read()
        # otherwise assume it's straight command-line input: 
        # join all into one string, separated by newlines
        else:
            input = '\n'.join(arguments)
    # if there are no arguments, assume it's going through a pipe and read stdin
    # TODO: if stdin is after a newline, on the command line, it never exits! fix? who cares
    except IndexError:
        input = sys.stdin.read()
    # returns the input as a list of lines instead of a single string
    return input.split('\n')

### If called directly, test everything
if __name__ == '__main__':
    input = read_input()
    # TODO what's '\b'?
    print '\b'.join(input)
