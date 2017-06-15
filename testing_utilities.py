#! /usr/bin/env python2.7

"""
Various testing utilities I wrote (for unit tests, functional tests etc) - see docstring for each function/class for what it does.
 --Weronika Patena, 2012
"""

import itertools
import unittest
import re
import os, shutil

debug=0

######### Output file comparison for functional testing
# The reason for this function is to be able to compare an output file and a reference file but ignore things that are expected to change (like the date in the header) - the reference file can contain lines starting with <REGEX> and then giving a regular expression (.* is fine), in which case the output file lines are compared against the regex instead, which allows wildcards for dates/etc.  See my stackoverflow question http://stackoverflow.com/questions/9726214/testing-full-program-by-comparing-output-file-to-reference-file-whats-it-calle


def _advance_iter_keep_state(iterator, skip_IGNORE_lines=False, skip_empty_lines=False):
    """ Try to advance the iterator; return ('', True) if StopIteration was raised, (next_item, False) otherwise. 
    Skip empty (or whitespace-only) lines and lines starting with <IGNORE> if the appropriate flags are True. 
    """
    stop_iteration = False
    # try to advance the iterator: at least once (that's what first_pass is for), and more times if line should be skipped 
    first_pass = True
    while first_pass or (skip_IGNORE_lines and next_item.startswith('<IGNORE>'))\
                     or (skip_empty_lines and not next_item.strip(' \t\n\r')):
        first_pass = False
        try:
            next_item = iterator.next()
        except StopIteration:
            next_item = ''
            stop_iteration = True
            break
    return next_item, stop_iteration
    # MAYBE-TODO should probably write unit-tests for this too, especially if I end up using it for multiple purposes


def clean_line(curr_line, clean_whitespace=True, make_lowercase=False):
    """ make line lowercase; change between-word whitespace to single spaces, remove start/end whitespace."""
    new_line = curr_line
    if make_lowercase: 
        new_line = new_line.lower()
    # replace all multi-space-tab strings with single space; remove space from start/end
    if clean_whitespace:
        new_line = re.sub('\s+', ' ', new_line)
        new_line = new_line.strip()
    return new_line
    # MAYBE-TODO should probably write unit-tests for this too, especially if I end up using it for multiple purposes


def compare_files_with_regex(iter_1, iter_2, 
                             ignore_empty_lines=False, ignore_whitespace=False, ignore_case=False):
    """ return True if the two string iterators are identical (except for special cases), else description of difference.

    The inputs can be lists of strings, generators, sets, open files, etc.  NO INPUT CHECKING IS DONE.
    The program goes over both iterators in parallel (always looking at element1 and element2), as follows:
     - if neither element starts with either special string (<REGEX> or <IGNORE>), compare them directly:
         if the elements are identical, advance both iterators, otherwise return the two elements
     - if one element starts with <REGEX>, parse the rest of it as a regular expression: 
         (removing whitespace from start/end, and adding ^ and $ to start and end)
         if the other element matches the regex, advance both iterators, otherwise return the two elements
     - if both elements start with <REGEX>, raise an exception - can't compare two regular expressions!
     - if either element starts with <IGNORE>, ignore it and advance its iterator to the next element.
    If the end of both iterators has been reached without finding a mismatch, return True; 
     if the end of one but not the other has been reached, return a string stating that and the first extra line. 

    If ignore_empty_lines is True, all lines containing nothing or just whitespace are ignored (treated like <IGNORE>).
    If ignore_whitespace is True, change all sequences of multiple spaces/tabs to a single space, 
     and remove all space at the beginning and end, before doing comparisons.
    If ignore_case is True, change all elements to lowercase before doing comparisons.
    """
    # MAYBE-TODO add some kind of input-checking? (If I do that, add it to unit-tests too)
    # if the arguments are lists/sets/something, convert them to iterators; this leaves iterators unchanged.
    iter1, iter2 = iter(iter_1), iter(iter_2)

    if debug:   print " *** NEW COMPARISON ***"
    while True:

        # advance both iterators to compare the next line pair (skipping IGNORE lines, and empty lines if requested)
        # if either iterator is stopped, exit immediately
        (line1, iter1_stopped) = _advance_iter_keep_state(iter1,skip_IGNORE_lines=True,skip_empty_lines=ignore_empty_lines)
        (line2, iter2_stopped) = _advance_iter_keep_state(iter2,skip_IGNORE_lines=True,skip_empty_lines=ignore_empty_lines)
        if iter1_stopped or iter2_stopped:  break


        if debug:   print 'raw lines:\n\t1) "%s"\n\t2) "%s"'%(line1, line2)

        # raise exception if both elements are regexes - can't match two regexes!
        if line1.startswith('<REGEX>') and line2.startswith('<REGEX>'):
            raise ValueError("Both elements start with <REGEX>, can't match! %s, %s."%(line1,line2))

        # the comparison can be stopped at any point by a stop marker - if the files matched up to that point, return True
        if any([line.startswith("<STOP COMPARING FILE>") for line in (line1,line2)]):
            return True

        # clean up non-regex lines: clean up whitespace, make lowercase
        for curr_line in line1,line2:
            if not line1.startswith('<REGEX>'):     line1 = clean_line(line1, ignore_whitespace, ignore_case)
            if not line2.startswith('<REGEX>'):     line2 = clean_line(line2, ignore_whitespace, ignore_case)

        if debug:   print 'cleaned lines:\n\t1) "%s"\n\t2) "%s"'%(line1.strip(), line2.strip())

        # if one of the lines is a regex, apply it to the other line; return both lines if they don't match
        flags = re.IGNORECASE if ignore_case else 0
        if line1.startswith('<REGEX>'):
            if not re.match('^%s$'%(line1[7:].strip()), line2, flags=flags):  return (line1.strip(),line2.strip())
        elif line2.startswith('<REGEX>'):
            if not re.match('^%s$'%(line2[7:].strip()), line1, flags=flags):  return (line1.strip(),line2.strip())

        # if neither line is a regex, compare them: return both lines if they don't match
        else:
            if not line1==line2:  return (line1.strip(),line2.strip())

        # if there wasn't a return or exception, just repeat the while loop

    # End condition: if all lines matched and both iterators are empty, return True; 
    #  if one iterator still has non-skipped lines left, return an info line and the next line from the other iterator.
    if iter1_stopped and iter2_stopped: return True
    elif iter1_stopped:                 return ("The first iterator ended. Second iterator next line:\n", line2.strip())
    elif iter2_stopped:                 return ("The second iterator ended. First iterator next line:\n", line1.strip())


######### Running functional tests, with or without input/output file comparison


def run_functional_tests(list_of_test_data, option_parser, function_to_run, test_folder, smoke_tests=False, 
                         argument_converter=(lambda *x: x), outfile_option='', append_to_outfilenames=''):
    """ Run tests in list_of_test_data using function_to_run: check outfiles against reference (unless smoke_tests=True).

    The list_of_test_data argument should be a list of (testname,descr,args) or (testname,args) tuples: 
     testname is a short test name (and will be used for outfile names); descr is an optional longer description string; 
     and args is a string of command-line arguments that will be passed to option_parser.  

    For each test, the following is done:
        - add " outfile_option test_folder/testname+append_to_outfilenames" to args to direct the output, 
           to customize how the outfile is specified, and create proper outfile name from test name.
           (just " filename" by default, or " -o filename" or " --output_file filename" or such if desired)
        - parse args using option_parser (which should be an optparse parser object)
        - run "function_to_run(*final_arguments)" to generate the outfiles, where final_arguments is the output of
            "argument_converter(option_parser,options,args)" - by default argument_converter does nothing, so you end up
            running "function_to_run(option_parser,options,args)", but for example if you want to run 
            "function_to_run(args,options)" instead, pass "lambda p,o,a: a,o" as argument_converter. 
           This should result in the creation of one or more output files per test, in test_folder.
        If smoke_tests is True, nothing else is done - we're just making sure the function runs without issues.
        Otherwise, auto-detect output reference files for each test in test_folder (they should start with test_name), 
         compare to actual output files using compare_files_with_regex, print information if the comparison fails.

    The program will print the test name/descr/args to stdout, along with details on the test failure if there is one.
    It stops with the first failure, and returns 1; if all tests passed, returns 0.
    """
    tmp_outfile_basename = "test_tmp"

    if smoke_tests:
        print("\n*** SMOKE TEST RUNS - THE OUTPUT IS NOT CHECKED, WE'RE JUST MAKING SURE THERE ARE NO ERRORS. ***")
        test_type = 'smoke-test'
    else:
        print("\n*** CHECKED TEST RUNS - THE OUTPUT IS CHECKED AGAINST CORRECT REFERENCE FILES. ***")
        test_type = 'checked-test'

    # check for existence of test_folder: for smoke-tests create new one, otherwise just make sure it exists.
    if smoke_tests:
        if os.access(test_folder,os.F_OK):
            print("Test output files will be saved in the %s folder (already present - removing it now)."%test_folder)
            shutil.rmtree(test_folder)
        else:
            print("Test output files will be saved in the %s folder (not present - creating it now)."%test_folder)
        os.mkdir(test_folder)
    else:
        if not os.access(test_folder,os.F_OK):
            print("Error: test folder %s doesn't exist or can't be accessed - can't run tests."%test_folder)
            return 1

    # remove any output test files that may have been left over from previous runs
    #  (if a previous test run failed, the output files aren't deleted, so if the first new test doesn't generate
    #   an outfile with the same name, it'll register an extra outfile with no matching infile and give an error)
    old_test_tmp_outfiles = [f for f in os.listdir(test_folder) if f.startswith(tmp_outfile_basename)]
    for old_outfile in old_test_tmp_outfiles:   
        old_outfile = os.path.join(test_folder,old_outfile)
        os.remove(old_outfile)

    # running each test
    for full_test_data in list_of_test_data:

        # each test tuple must provide test_args, and optionally provides test_name,test_descr as well - check by length
        if len(full_test_data)==3:
            test_name, test_descr, test_args = full_test_data
        elif len(full_test_data)==2:
            test_name, test_args = full_test_data
            test_descr = ''
        else:
            raise Exception("Couldn't parse test_case %s! Needs to be (name,descr,args) or (name,args)."%full_test_data)

        # print info about the test to the command-line
        if test_descr and test_name:    print(" * New %s run: %s (%s)."%(test_type, test_descr, test_name))
        elif test_descr or test_name:   print(" * New %s run: %s."     %(test_type, (test_descr or test_name)))
        else:                           print(" * New %s run."         %(test_type))
        print("   Arguments: %s"%test_args)

        # smoke-tests use test name as outfile name; normal tests use tmp_outfile_basename
        outfilename = (test_name if smoke_tests else tmp_outfile_basename) + append_to_outfilenames
        test_args_and_outfile = test_args + " %s %s"%(outfile_option, os.path.join(test_folder,outfilename))
        (options, args) = option_parser.parse_args(test_args_and_outfile.split())

        # actually run function_to_run with the given options/arguments (run them through argument_converter first)
        processed_function_args = argument_converter(option_parser,options,args)
        function_to_run(*processed_function_args)

        ### if it's a smoke-test, this is all, just go on to the next test; 
        if smoke_tests: continue

        ### if it's not a smoke-test, compare output files to reference files
        # find output reference files in test_folder to compare the output files against
        test_reference_files = [f for f in os.listdir(test_folder) if f.startswith(test_name)]
        # for each reference output files found, compare the tmp output file, 
        #  and fail the test if the output file doesn't exist or differs from the reference file.
        for reffile in test_reference_files:
            reffile = os.path.join(test_folder,reffile)
            expected_outfile = reffile.replace(test_name, tmp_outfile_basename)
            if not os.path.exists(expected_outfile):
                print("TEST FAILED!!  Output file %s (to match reference file %s) doesn't exist."%(expected_outfile,
                                                                                                   reffile))
                return 1
            with open(reffile,'r') as REFFILE:
                with open(expected_outfile,'r') as OUTFILE:
                    file_comparison_result = compare_files_with_regex(REFFILE, OUTFILE)
            if file_comparison_result==True:
                os.remove(expected_outfile)
            else:
                print("TEST FAILED!!  Reference file %s and output file %s differ. MORE INFO "%(reffile,expected_outfile)
                      +"ON DIFFERENCE (the two mismatched lines or an error message): '%s', '%s'"%file_comparison_result)
                return 1
        # make sure there aren't any extra output files without reference files
        test_output_files = [f for f in os.listdir(test_folder) if f.startswith(tmp_outfile_basename)]
        for outfile in test_output_files:
            outfile = os.path.join(test_folder,outfile)
            expected_reffile = outfile.replace(tmp_outfile_basename, test_name)
            if not os.path.exists(expected_reffile):
                print("TEST FAILED!!  Output file %s has no matching reference file (%s)."%(outfile,expected_reffile))
                return 1

    if smoke_tests:
        print("*** Smoke test runs finished. If you didn't get any errors, that's good (warnings are all right). "
              + "You can check the output files to make sure they look reasonable (this is NOT done automatically!). ***")
    else:
        print("*** Checked test runs finished - EVERYTHING IS FINE. ***")
    return 0

    # I don't have a unit-test or functional test for it, but I rewrote do_test_run functions in ~/experiments/generating_library/combinatorial_pooling/code/robotic_plate_transfer.py and ~/experiments/mutant_pool_screens/mutant_deepseq_analysis/code/mutant_count_alignments.py to use it, and it works right.


############################## unit-tests of the functions in this file ##################################

class Testing__everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__compare_files_with_regex(self):
        # MAYBE-TODO in most cases when there shouldn't be a match I'm only checking that the output!=True, not the exact non return value which describes the mismatch - could fix that, but I'm not sure it's worth it.
        # any list should be identical to itself
        for test_list in [[], ['hi','there','thing'], ['1','2','3']]:
            assert compare_files_with_regex(test_list, test_list) == True
        # a list with all IGNORE lines should match the empty list or any other IGNORE-only lists, but not any other list
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            ignore_list_1 = ['<IGNORE> '+l for l in test_list]
            ignore_list_0 = ['<IGNORE> stuff' for l in test_list]
            ignore_list_2 = ['<IGNORE> '+l+l for l in test_list]
            ignore_list_long = ignore_list_0+ignore_list_1+ignore_list_2
            for ignore_list in [ignore_list_1, ignore_list_0, ignore_list_2]:
                assert compare_files_with_regex(ignore_list, ignore_list) == True
                assert compare_files_with_regex(ignore_list, []) == True
                assert compare_files_with_regex([], ignore_list) == True
                # the two below are test cases for when one iterator stopped and the other didn't:
                #  the first line of the return tuple is an info line I don't feel like checking the details of, 
                #   the second is the extra line from whichever iterator was longer (i.e. first line of test_list).
                assert compare_files_with_regex(ignore_list, test_list)[1] == test_list[0]
                assert compare_files_with_regex(test_list, ignore_list)[1] == test_list[0]
            for list1,list2 in itertools.permutations([[],ignore_list_1,ignore_list_0,ignore_list_2,ignore_list_long],2):
                assert compare_files_with_regex(list1, list2) == True
        # a list should match itself no matter how many IGNORE lines are added. 
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            ignore_list = ['<IGNORE> '+l for l in test_list]
            ignore_list_0 = ['<IGNORE> stuff' for l in test_list]
            ignore_list_2 = ['<IGNORE> '+l+l for l in test_list]
            ignore_list_long = ignore_list_0+ignore_list_1+ignore_list_2
            for list1,list2 in itertools.permutations([[], ignore_list, ignore_list_0, ignore_list_2, ignore_list_long],2):
                assert compare_files_with_regex(test_list, list1+test_list) == True
                assert compare_files_with_regex(test_list+list1, test_list+list2) == True
                assert compare_files_with_regex(test_list+list1+test_list, test_list+list2+test_list) == True
        # a list with all REGEX ".*" lines should match any list of the same length, 
        #  but not of a different length (except IGNORE lines); 
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            any_regex_list = ['<REGEX>.*' for l in test_list]
            assert compare_files_with_regex(any_regex_list, test_list) == True
            assert compare_files_with_regex(test_list, any_regex_list) == True
            ignore_list = ['<IGNORE>']*3
            assert compare_files_with_regex(ignore_list+any_regex_list, test_list) == True
            assert compare_files_with_regex(any_regex_list, ignore_list+test_list) == True
            assert compare_files_with_regex(ignore_list+any_regex_list, test_list+ignore_list+ignore_list) == True
        # comparing two REGEX lines should give an exception; 
        #  comparing two lists with REGEX lines in different positions shouldn't
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            any_regex_list = ['<REGEX>.*' for l in test_list]
            self.assertRaises(ValueError, compare_files_with_regex, any_regex_list, any_regex_list)
            first_regex_list = ['<REGEX>.*'] + test_list[1:]
            last_regex_list = test_list[:-1] + ['<REGEX>.*']
            self.assertRaises(ValueError, compare_files_with_regex, first_regex_list, first_regex_list)
            self.assertRaises(ValueError, compare_files_with_regex, last_regex_list, last_regex_list)
            compare_files_with_regex(first_regex_list, last_regex_list) == True
            compare_files_with_regex(last_regex_list, first_regex_list) == True
        # some specific simple ignore and regex testing
        assert compare_files_with_regex(['1','2','3'],['1','2','3']) == True
        assert compare_files_with_regex(['2','1','3'],['1','2','3']) == ('2','1')
        assert compare_files_with_regex(['12','2','3'],['1','2','3']) == ('12','1')
        assert compare_files_with_regex(['<REGEX>[12]','2','3'],['1','2','3']) == True
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['1','2','3']) == ('2','1')
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['2','3']) == True
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['<IGNORE>1','2','3']) == True
        # some specific testing of more complicated regexes
        assert compare_files_with_regex(['<REGEX>.*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>some.*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>.*some'],['some text here']) == ('<REGEX>.*some','some text here')
        assert compare_files_with_regex(['<REGEX>.*here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>some text here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>.*some text here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[sometxhr ]*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[eosmtxhr ]'],['some text here'])==('<REGEX>[eosmtxhr ]','some text here')
        assert compare_files_with_regex(['<REGEX>[omes]*\s[xet]*\s[ehr]*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[oes]*\s[xet]*\s[ehr]*'],['some text here']) != True
        assert compare_files_with_regex(['<REGEX>[omes]\s[xet]*\s[ehr]*'],['some text here']) != True
        assert compare_files_with_regex(['<REGEX>[omes]*[xet]*\s[ehr]*'],['some text here']) != True
        # testing ignore_empty_lines
        for other_list in [['', 'some text'], ['some text', '    \t'], ['', '\t', 'some text', '', '']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_empty_lines=False) != True
            assert compare_files_with_regex(['some text'], other_list, ignore_empty_lines=True) == True
        # testing ignore_whitespace
        for other_list in [['some  text'], ['some\ttext'], ['some \t\t text \t']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_whitespace=False) == ('some text', 
                                                                                                    other_list[0].strip())
            assert compare_files_with_regex(['some text'], other_list, ignore_whitespace=True) == True
            assert compare_files_with_regex(['\t\tsome text'], other_list, ignore_whitespace=False) != True
            assert compare_files_with_regex(['\t\tsome text'], other_list, ignore_whitespace=True) == True
            # removing the whitespace completely doesn't work, though
            assert compare_files_with_regex(['sometext'], other_list, ignore_whitespace=False) == ('sometext',
                                                                                                   other_list[0].strip())
            assert compare_files_with_regex(['sometext'], other_list, ignore_whitespace=True) != True
        # testing ignore_case
        for other_list in [['SOME text'], ['Some Text'], ['sOmE tExT']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_case=False) == ('some text',other_list[0])
            assert compare_files_with_regex(['some text'], other_list, ignore_case=True) == True
            assert compare_files_with_regex(['SOME TEXT'], other_list, ignore_case=False) == ('SOME TEXT',other_list[0])
            assert compare_files_with_regex(['SOME TEXT'], other_list, ignore_case=True) == True
        ### testing on files instead of lists - remember to restart the iterators every time!
        if debug:   print " ************* file tests **************** "
        # non-regex files match themselves
        file1, file2, file3 = ('_test_inputs/textcmp_file1.txt','_test_inputs/textcmp_file2.txt',
                               '_test_inputs/textcmp_file3.txt')
        with open(file1,'r') as F1:
            with open(file1,'r') as F1_:
                assert compare_files_with_regex(F1, F1_) == True
        with open(file3,'r') as F3:
            with open(file3,'r') as F3_:
                assert compare_files_with_regex(F3, F3_) == True
        # comparing regex file to itself causes an exception
        with open(file2,'r') as F2:
            with open(file2,'r') as F2_:
                self.assertRaises(ValueError, compare_files_with_regex, F2, F2_)
        # comparing a non-regex file to a matching file with regex/ignore matches
        with open(file1,'r') as F1:
            with open(file2,'r') as F2:
                assert compare_files_with_regex(F1, F2) == True
        # comparing to the third file with extra empty lines, whitespace and case differences - all three True/False 
        #  arguments must be true to give a match, no other combination works.
        for first_file in (file1, file2):
            if debug:   print " ****** %s vs %s ******"%(first_file, file3)
            with open(first_file,'r') as FILE1:
                with open(file3,'r') as F3:
                    assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=False, 
                                                    ignore_whitespace=False, ignore_case=False) != True
            with open(first_file,'r') as FILE1:
                with open(file3,'r') as F3:
                    assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=True, 
                                                    ignore_whitespace=True, ignore_case=True) == True
            for v1,v2,v3 in itertools.chain(itertools.permutations([True,False,False],3), 
                                            itertools.permutations([True,True, False],3)):
                with open(first_file,'r') as FILE1:
                    with open(file3,'r') as F3:
                        assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=v1, 
                                                        ignore_whitespace=v2, ignore_case=v3) != True
            # MAYBE-TODO note that theoretically compare_files_with_regex is symmetrical 
            #  (order of file1, file2 doesn't matter), but do I actually test that?  Maybe I should.


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
