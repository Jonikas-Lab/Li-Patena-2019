#! /usr/bin/env python2.7

"""
Various statistical convenience functions I wrote - see docstring for each function for what it does.
Weronika Patena, 2010-2013
"""

# standard library
from __future__ import division 
import sys, os
import unittest
import itertools
from collections import defaultdict
import random
# other packages
import numpy
import scipy.stats
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
R_stats = importr('stats')
# my modules


### HELP FUNCTIONS

def array_1D(x):
    """ Convert to 1D numpy array. """
    return numpy.reshape(numpy.array(x), -1)


### STATISTICAL FUNCTIONS

def fisher_exact(contingency_table):
    """ Do a Fisher's exact test using the R version - useful for tables larger than 2x2.

    I'm writing this because scipy.stats.fisher_exact only does 2x2 tables, and sometimes more is useful. 
    Contingency_table should be a list of two lists of arbitrary equal length (like [[2, 3, 1], [10, 100, 1000]])
    Returns a p-value.
    """ 
    # I don't really know how to make a matrix in R, so this may not be the most direct way, but it works
    vector = robjects.IntVector(sum(contingency_table, []))
    matrix = robjects.r.matrix(vector, nrow=2, byrow=True)
    result = R_stats.fisher_test(matrix)
    return result[0][0]


def chisquare_goodness_of_fit(category_counts, expected_frequencies, dof_subtract=0, return_pvalue_only=True, min_count=50):
    """ Gives p-value for whether a list of category counts is different from the expected frequencies, using the chi-square test.

    Simple wrapper around scipy.stats.chisquare, that calculates the expected counts by normalizing expected_frequencies 
     to have the same total as category_counts.

    Dof_subtract is how many degrees of freedom to SUBTRACT from the default (usually no adjustment needed, 0).

    Raises ValueError if any count is below min_count - the chi-square test shouldn't be used for small numbers.

    If return_pvalue_only is True, returns only pvalue, otherwise (chisquare_statistic, pvalue).
    """
    # convert everything to 1D numpy arrays
    category_counts = array_1D(category_counts)
    expected_frequencies = array_1D(expected_frequencies)
    if numpy.min(category_counts) < min_count:
        raise ValueError("Shouldn't use chisquare_goodness_of_fit with small numbers! Adjust min_count if you have to.")
    # using scipy.sum  and numpy.array in case category_counts/expected_frequencies are multi-dimensional matrices
    norm_factor = sum(category_counts) / sum(expected_frequencies)
    expected_counts_same_total = expected_frequencies * norm_factor
    chisq, pval = scipy.stats.chisquare(category_counts, expected_counts_same_total, ddof=dof_subtract)
    if return_pvalue_only:  return pval
    else:                   return chisq, pval


def chisquare_independence(category_counts_a, category_counts_b, dof_subtract=0, return_pvalue_only=True, min_count=50):
    """ Gives p-value for whether two lists of category counts are different, using the chi-square test.

    Dof_subtract is how many degrees of freedom to SUBTRACT from the default (usually no adjustment needed, 0).

    If return_pvalue_only is True, returns only pvalue, otherwise (chisquare_statistic, pvalue).

    Raises ValueError if any count is below min_count - the chi-square test shouldn't be used for small numbers, 
     use Fisher's exact test instead.

    NOTE ON HOW THIS WORKS: to do a chi-square test of independence, you compare the counts of eiter category to expected counts, 
     where expected counts are CALCULATED FROM BOTH CATEGORIES - you DON'T directly compare one to the other,
      the way you do in a goodness-of-fit test!
        EXAMPLE: if counts_A are 110 and 90, and counts_B are 190 and 10, you don't do a chisquare on [110,90], [190,10], 
         but calculate the expected counts for all four categories from the row/column totals: 
          from A+B (110+190=300 and 90+10=100), so since sum(A) and sum(B) are both 200, the expected is [150,50,150,50]. 
          So compare [110,90,190,10] (A and B together) to [150,50,150,50] using the chi-square goodness-of-fit test, 
          BUT make a degree-of-freedom adjustment: the comparison we're REALLY doing is 2*3 (comparing two 3-length datasets), 
           so the degrees of freedom should be (2-1)*(3-1) = 2.  However, we're transforming it in to a comparison between 
            6-length observed and expected datasets, with (2-1)*(6-1)=5 degrees of freedom by default, 
           so we need to subtract 3 from the degrees-of-freedom for the test (in addition to whatever dof_subtract we already have). 

        EXAMPLE 2: what if the totals in A and B are different?  Say A is [110,190] and B is [90,10]. 
            Then the overall counts are [200,200], so we should be comparing [110,190,90,10] to [150,150,50,50]. 
            Again, with DOF adjustment of 3.

    See https://udel.edu/~mcdonald/statchiind.html for source and description of how it should work, 
     and more examples at http://stattrek.com/chi-square-test/independence.aspx 
     and http://omega.albany.edu:8008/mat108dir/chi2independence/chi2in-m2h.html - it was confusing to me at first.  
    """
    # just convert everything to scipy/numpy arrays
    category_counts_a = array_1D(category_counts_a)
    category_counts_b = array_1D(category_counts_b)

    all_observed = numpy.append(category_counts_a, category_counts_b)

    if min(all_observed) < min_count:
        raise ValueError("Shouldn't use chisquare_goodness_of_fit with small numbers! Adjust min_count if you have to.")

    # calculate the expected frequencies for all categories and both datasets
    both_count_totals = category_counts_a + category_counts_b
    sum_a, sum_b = sum(category_counts_a), sum(category_counts_b)
    full_sum = sum_a + sum_b
    both_count_totals_norm_a = both_count_totals*sum_a/full_sum
    both_count_totals_norm_b = both_count_totals*sum_b/full_sum
    all_expected = numpy.append(both_count_totals_norm_a, both_count_totals_norm_b)
    # Degrees-of-freedom adjustment (see docstring example for details)
    proper_dof = len(category_counts_b) - 1
    transformed_dof = len(all_observed) - 1
    extra_dof_adjustment = proper_dof - transformed_dof
    full_dof_subtract = dof_subtract - extra_dof_adjustment
    # do the chi-square goodness-of-fit test of observed vs expected
    return chisquare_goodness_of_fit(all_observed, all_expected, full_dof_subtract, return_pvalue_only, min_count)


def FDR_adjust_pvalues(pvalue_list, N=None, method='BH'):
    """ Adjust a list of p-values for false discovery rate using R's stats::p.adjust function.

    N and method are passed to R_stats.p_adjust: 
     - N is the number of comparisons (if left unspecified, defaults to len(pvalue_list), I think)
     - method is the name of the adjustment method to use (inherited from R)

    Note that this MUST be done after all the p-values are already collected, on the full list of p-values at once:
     trying to do it on single p-values, even with adjusted N, will give different results!
    """
    if not method in R_stats.p_adjust_methods:
        raise ValueError("Unknown method %s - method must be one of (%s)!"%(method, ', '.join(R_stats.p_adjust_methods)))
    if N is None:   return R_stats.p_adjust(FloatVector(pvalue_list), method=method)
    else:           return R_stats.p_adjust(FloatVector(pvalue_list), method=method, n=N)


def binomial_CI(n, N, pct, a=1, b=1, n_pbins=1e3):
    """ Computes binomial confidence interval (Highest Posterior Density Region).

    Function computes the posterior mode along with the upper and lower bounds of the
    Highest Posterior Density Region.
    By mtw729 from Stack Overflow: 
    stackoverflow.com/questions/13059011/is-there-any-python-function-library-for-calculate-binomial-confidence-intervals

    Parameters
    ----------
    n: number of successes 
    N: sample size 
    pct: the size of the confidence interval (between 0 and 1)
    a: the alpha hyper-parameter for the Beta distribution used as a prior (Default=1)
    b: the beta hyper-parameter for the Beta distribution used as a prior (Default=1)
    n_pbins: the number of bins to segment the p_range into (Default=1e3)

    Returns
    -------
    A tuple that contains the mode as well as the lower and upper bounds of the interval
    (mode, lower, upper)

    """
    # fixed random variable object for posterior Beta distribution
    rv = scipy.stats.beta(n+a, N-n+b)
    # determine the mode and standard deviation of the posterior
    stdev = rv.stats('v')**0.5
    mode = (n+a-1.)/(N+a+b-2.)
    # compute the number of sigma that corresponds to this confidence
    # this is used to set the rough range of possible success probabilities
    n_sigma = numpy.ceil(scipy.stats.norm.ppf( (1+pct)/2. ))+1
    # set the min and max values for success probability 
    max_p = mode + n_sigma * stdev
    if max_p > 1:
        max_p = 1.
    min_p = mode - n_sigma * stdev
    if min_p > 1:
        min_p = 1.
    # make the range of success probabilities
    p_range = numpy.linspace(min_p, max_p, n_pbins+1)
    # construct the probability mass function over the given range
    if mode > 0.5:
        sf = rv.sf(p_range)
        pmf = sf[:-1] - sf[1:]
    else:
        cdf = rv.cdf(p_range)
        pmf = cdf[1:] - cdf[:-1]
    # find the upper and lower bounds of the interval 
    sorted_idxs = numpy.argsort( pmf )[::-1]
    cumsum = numpy.cumsum( numpy.sort(pmf)[::-1] )
    j = numpy.argmin( numpy.abs(cumsum - pct) )
    upper = p_range[ (sorted_idxs[:j+1]).max()+1 ]
    lower = p_range[ (sorted_idxs[:j+1]).min() ]    

    return (mode, lower, upper)


def R_clear_environment(garbage_collection_cycles=3):
    """ Attempt to release memory held by R via rpy2.

    Apparently doing a lot of calls to the rpy2-using functions here can cause memory usage to increase until the process is killed.
    This is my attempt at fixing that based on a few StackOverflow answers.  Basically, delete all variables in R, 
    explicitly run python garbage collection and R garbage collection, possibly multiple times.

    Sources:
     - http://stackoverflow.com/questions/5199334/clearing-memory-used-by-rpy2, 
     - http://stackoverflow.com/questions/8144956/python-rpy2-module-refresh-global-r-environment?rq=1
     - http://stackoverflow.com/questions/12277094/memory-leak-with-rpy?noredirect=1&lq=1
    """
    import rpy2.robjects as R
    import gc
    R.r('rm(list = ls(all.names=TRUE))')
    for i in range(garbage_collection_cycles):
        gc.collect()
        R.r('gc()')


# OLD NOTES ON FDR CORRECTION:
    # How to do FDR correction?  According to the Handbook of Biological Statistics (https://udel.edu/~mcdonald/statmultcomp.html), Benjamini-Hochberg correction is probably what I want.  They describe a procedure, but it's slightly odd, because it doesn't give a p-value (q-value?) for each window, just a yes/no significance result based on the p-value and the desired false discovery rate. 
    # I didn't find any obvious way of doing this directly in python, but there's an R function "p.adjust" (http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html), which I can use in python with rpy2 (http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python).  Trying that, with just a few test values:
    #  * get the p_values for a few mutant counts per bin, between 0 and 15:
    #     >>> p_values = [scipy.stats.binom_test(x,12061,20000/113000000) for x in (0,0,1,2,5,10,10,25,25)]
    #     >>> p_values
    #     [0.28620628492789102, 0.28620628492789102, 0.73047985928763548, 1.0, 0.065605526425554839, 7.8933016187778668e-05, 7.8933016187778668e-05, 1.3921234115131231e-18, 1.3921234115131231e-18]
    #  * try the FDR adjustment with default options - the p-values increase a bit:
    #     >>> from rpy2.robjects.packages import importr
    #     >>> R_stats = importr('stats')
    #     >>> from rpy2.robjects.vectors import FloatVector
    #     >>> p_adjusted = R_stats.p_adjust(FloatVector(p_values), method = 'BH')
    #     >>> list(p_adjusted)
    #     [0.36797950919300276, 0.36797950919300276, 0.8217898416985899, 1.0, 0.1180899475659987, 0.000177599286422502, 0.000177599286422502, 6.264555351809054e-18, 6.264555351809054e-18]
    #  * same, but using the correct N - I'm reporting 9 random p-values here, but we actually did 50000 tests (50000 windows), not just 9!  The values went down further - good.
    #     >>> p_adjusted2 = R_stats.p_adjust(FloatVector(p_values), method = 'BH', n=50000)
    #     >>> list(p_adjusted2)
    #     [1.0, 1.0, 1.0, 1.0, 1.0, 0.9866627023472333, 0.9866627023472333, 3.4803085287828076e-14, 3.4803085287828076e-14]


######################################## TESTS FOR THIS FILE ########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__fisher_exact(self):
        # for a 2x2 table, compare to scipy.stats.fisher_exact (using self.assertAlmostEqual for float comparison)
        table_2x2 = [[2, 2], [100, 1000]]
        pval = scipy.stats.fisher_exact(table_2x2)[1]
        pval2 = fisher_exact(table_2x2)
        self.assertAlmostEqual(pval, pval2)
        # for a 2x3 table just compare to a result I got in R
        table_2x3 = [[2, 2, 2], [10, 100, 1000]]
        pval2 = fisher_exact(table_2x3)
        self.assertAlmostEqual(pval2, 0.0001392)

    def test__chisquare_goodness_of_fit(self):
        kwargs = dict(return_pvalue_only=False, min_count=1)
        # test case 1 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_goodness_of_fit([423,133], [3,1], **kwargs)
        assert round(chisq,2) == 0.35
        assert round(pval,3) == 0.557
        # test case 2 from https://udel.edu/~mcdonald/statchigof.html - note that the numbers here are too low for this test, really
        chisq,pval = chisquare_goodness_of_fit([70, 79, 3, 4], [0.54, 0.4, 0.05, 0.01], **kwargs)
        assert round(chisq,3) == 13.593
        assert round(pval,4) == 0.0035
        # test case 3 from https://udel.edu/~mcdonald/statchigof.html - NOTE that this should have 1 degree of freedom
        chisq,pval = chisquare_goodness_of_fit([14, 21, 25], [0.167, 0.483, 0.350], dof_subtract=1, **kwargs)
        assert round(chisq,2) == 4.54
        assert round(pval,3) == 0.033

    def test__chisquare_independence(self):
        # test case 1 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_independence([268,199, 42], [807,759,184], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,2) == 7.26
        assert round(pval,3) == 0.027
        # test case 2 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_independence([127, 99, 264], [116, 67, 161], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,2) == 6.26
        assert round(pval,3) == 0.044
        # test case from http://stattrek.com/chi-square-test/independence.aspx
        chisq,pval = chisquare_independence([200, 150, 50], [250, 300, 50], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,1) == 16.2
        assert round(pval,4) == 0.0003

    def test__FDR_adjust_pvalues(self):
        self.assertRaises(ValueError, FDR_adjust_pvalues, [0,0.1,1], method='FAKE')
        # test based on https://udel.edu/~mcdonald/statmultcomp.html, MODIFIED.
        #  I'm just implementing the math as described on the website, NOT using the output values given on the webpage 
        #   (because they're not q-values but comparison values)
        #  Also, R's p_adjust gives different results than the described math for identical p-values, so I'm not including any. 
        #   (which makes sense - it gives the same adjusted p-value for each identical p-value, rather than different ones)
        #  Actually the results seem different for even somewhat similar values, too... I edited the example to remove them, 
        #   to at least make sure this is APPROXIMATELY right.  I'm pretty sure R isn't wrong, anyway.  MAYBE-TODO better test?
        input_pvalues = [0.010, 0.032, 0.07, 0.20, 0.38, 0.68, 0.97]
        output = FDR_adjust_pvalues(input_pvalues, N=None, method='BH')
        # Calculate the adjusted p-values (which are the largest Q-values for which P<(i/m)Q, i.e. Pm/i), 
        #  and check that they match the output, using self.assertAlmostEqual for approximate float comparison:
        m = len(input_pvalues)
        expected_adjusted_pvalues = [p*m/(i+1) for (i,p) in enumerate(input_pvalues)]
        for obs,exp in zip(output, expected_adjusted_pvalues): self.assertAlmostEqual(obs,exp)
        # check that we get the same adjusted p-values regardless of the order in which they're put in
        for _ in range(10):
            random.shuffle(input_pvalues)
            output = FDR_adjust_pvalues(input_pvalues, N=None, method='BH')
            for obs,exp in zip(sorted(output), expected_adjusted_pvalues): self.assertAlmostEqual(obs,exp)


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
