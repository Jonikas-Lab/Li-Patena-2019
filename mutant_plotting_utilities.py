#! /usr/bin/env python2.7
"""
Plotting utilities specifically for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import glob
import os, sys
import collections
from collections import defaultdict
import math
import random
import itertools
# other packages
import numpy
import scipy
import scipy.stats
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import general_utilities
import statistics_utilities
import basic_seq_utilities
import plotting_utilities
import mutant_analysis_classes
import mutant_utilities
from mutant_utilities import get_histogram_data_from_positions, get_mutant_positions_from_dataset, get_chromosome_lengths
import mutant_simulations


########################################## Plotting mutant positions #############################################

######### Plotting mutant/other positions over chromosomes as heatmaps/lines (and help functions for that)

def _get_chromosomes_by_type(all_chromosomes, include_scaffolds, include_cassette, include_other, output_dict=False):
    """ Filter all_chromosomes based on the include_* arguments; return a list of only ones of the desired types. 

    If output_dict is True, output a type:chromosome_list dictionary instead of a single chromosome_list.
    """
    chosen_chromosomes = defaultdict(list)
    for chromosome in all_chromosomes:
        chr_type = basic_seq_utilities.chromosome_type(chromosome)
        if (chr_type=='chromosome' or (chr_type=='scaffold' and include_scaffolds) or (chr_type=='cassette' and include_cassette) 
            or (chr_type in ('chloroplast', 'mitochondrial', 'other') and include_other)):
            chosen_chromosomes[chr_type].append(chromosome)
    if output_dict:     return dict(chosen_chromosomes)
    else:               return sum(chosen_chromosomes.values(), [])


def _get_plotline_pos(middle_pos, total_width, relative_widths, total_N, N):
    """ Give the left/right edge of the Nth section out of total_N with given relative_widths and total_width.

    Divide a total_width centered at middle_pos into total_N sections (with
    relative widths given by widths, or equal if None), and return the left and
    right edge of the Nth section. 
    """
    if relative_widths is None: 
        relative_widths = [1 for _ in range(total_N)]
    real_width_sum = sum(relative_widths)
    absolute_widths = [x/real_width_sum*total_width for x in relative_widths]
    start_pos = middle_pos - total_width/2
    left_pos = start_pos + sum(absolute_widths[:N])
    right_pos = left_pos + absolute_widths[N]
    return left_pos, right_pos


def mutant_positions_and_data(datasets=[], dataset_formats=0, density_plots=True, colors=None, names=None, maxes_per_base=None, 
                              maxes_per_bin=None, strands=None, widths=None, title='', bin_size=mutant_utilities.DEFAULT_BIN_SIZE, 
                              chromosome_lengths=None, interpolate=False, NaN_color='0.6', total_plotline_width=0.6, 
                              condense_colorbars=True, no_colorbars=False, 
                              include_scaffolds=False, include_cassette=False, include_other=False):
    """ Plot multiple datasets (mutant,position,density) across chromosomes, as position lines or density heatmaps.

    Each element in datasets must match the corresponding value in dataset_formats (or the value if only one is given:
     - format 0 - a mutant_analysis_classes.Insertional_mutant_pool_dataset instance
     - format 1 - a chromosome:position_list dict (with positions as integers)
     - format 2 - a chromosome:bin_density_list dict (with bin densities as numbers - generated with bin sizes matching bin_size.)
    If datasets is a single instance of any of these things instead of a list, it'll be treated as a length-1 list.

    Dataset_formats, density_plots, colors, names, and strands must all be lists of the same length as datasets, or single values
     (which will be converted to lists filled with those values):
     * for each dataset, if the corresponding density_plots value is False, 
        each position will be plotted as a horizontal line (using the corresponding color value, or black if None); 
         (this is impossible if format==2 for this dataset - in that case an error will be returned)
       if True, the positions will be converted to a density map by binning mutants into ranges (with bin_size bases per range), 
        and plotted as a heatmap (using the corresponding color value as the colormap name, or the default colormap if None).
     * the maxes_per_base and maxes_per_bin values can be None or numbers, and are only relevant for density heatmap plots;
        a single sample cannot have both maxes specified, just one or neither. 
        - if both are None, the density heatmaps will be scaled based on the highest value in the dataset 
        - if maxes_per_base is a number X, the density colormap will be scaled so that X positions per BASE is 100%, 
            and the colorbar scale will be 0-100% instead of raw counts/bin 
                (use for values with a defined 0-100% range, like mappability or GC content)
        - if maxes_per_bin is a number X, the density colormap will be scaled so that X positions per BIN is the maximum, 
            (anything higher than X will be the same color as X), and the colorbar scale will remain in raw counts per bin 
                (use for values you wouldn't put on a 0-100% scale, like mutant positions, gene positions, etc)
     * strands is only needed for mutant datasets - each mutant dataset will be converted to a chromosome:position_list dict, 
        taking all mutants if the strands value is None, or only +strand/-strand/opposite-tandem mutants if it's '+'/'-'/'both'.
     * widths are the relative widths of the plots - they will all add up to total_plotline_width. If None, all will be equal.
     * The names values are used as labels for the legend (if density=False) or colorbar (if density=True); 
        if '' is given, that dataset won't be included in the legend or won't have a colorbar
         (use for cases with multiple copies of similar datasets - they should use the same color as a labeled one!)

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).
    Interpolate governs whether color blocks are drawn as rectangles or interpolated on density heatmaps. 
    If condense_colorbars, if there is more than one colorbar to be drawn, they'll be condensed into two rows with reduced spacing.
    Title will be used for the plot title.
    Total_plotline_width gives the fraction of x-axis space given to the chromosome plots vs spaces between chromosomes.
    """
    ### parse the arguments, set defaults, etc
    # for the arguments that should be lists, if a single value is given, change to a correct-length list filled with that value
    if isinstance(datasets, (dict, mutant_analysis_classes.Insertional_mutant_pool_dataset)):    datasets = [datasets]
    if isinstance(dataset_formats, int):                                   dataset_formats = [dataset_formats for _ in datasets]
    if density_plots in (True,False):                                      density_plots = [density_plots for _ in datasets]
    if colors is None or isinstance(colors,str):                           colors = [colors for _ in datasets]
    if names is None or isinstance(names,str):                             names = [names for _ in datasets]
    if maxes_per_base is None or isinstance(maxes_per_base,(int,float)):   maxes_per_base = [maxes_per_base for _ in datasets]
    if maxes_per_bin is None or isinstance(maxes_per_bin,(int,float)):     maxes_per_bin = [maxes_per_bin for _ in datasets]
    if isinstance(strands,str) or strands is None:                         strands = [strands for _ in datasets]
    imshow_kwargs = {} if interpolate else {'interpolation': 'none'}
    # set sensible default values for None colors (different for mutant and other datasets)
    #  - if density is True, use a colorbar (see ~/computers_and_programming/colormaps_all_matplotlib.png; *_r means reverse)
    #  - otherwise use a single color name
    for N,(color,density_plot,d_format) in enumerate(zip(list(colors),density_plots,dataset_formats)):
        if color is None:
            if density_plot:   colors[N]='gist_earth'     
            else:              colors[N]='black'
    # set sensible default values for None names ('mutants' for mutants, more general for other data, since we can't tell)
    for N,(name,d_format) in enumerate(zip(list(names),dataset_formats)):
        if name is None:
            if d_format==0:    names[N]='mutants'
            else:              names[N]='positions'

    ### get the list of chromosomes to plot
    # get the chromosome lengths from a file, if filename or None was given
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
    # grab only the chromosome types we want
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    # sort chromosomes properly by type and name
    all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)
    # MAYBE-TODO for short chromosomes like scaffolds/chlor/mito/cassette, plot in multiple rows? Maybe change the scale to stretch them a bit, since some will probably be smaller than bin_size?  But then I'd have to normalize them separately or something...

    ### get a list of all datasets provided, converted to the appropriate format for plotting, along with their extra args:
    #    - chrom:position_list format if density_plot is False
    #    - chrom:bin_density_list format if density_plot is True
    all_datasets_and_info = []
    for dataset,d_format,density_plot,color,name,max_per_base,max_per_bin,strand in zip(datasets, dataset_formats, density_plots, 
                                                                            colors, names, maxes_per_base, maxes_per_bin, strands):
        # for mutant datasets, use get_mutant_positions_from_dataset to convert each dataset to a chromosome:position_list dict
        if d_format==0:
            dataset = get_mutant_positions_from_dataset(dataset, strand=strand)
            d_format = 1
        # now convert each dataset to the appropriate format based on density_plot and format, or raise exception if impossible
        if (density_plot is False and d_format==1) or (density_plot is True and d_format==2):
            pass
        elif density_plot is False and d_format==2:
            raise ValueError("Dataset %s was provided in density_plot format - cannot plot as lines!"%name)
        elif density_plot is True and d_format==1:
            dataset = get_histogram_data_from_positions(dataset, bin_size, chromosome_lengths, all_chromosomes, 
                                                        special_last_bin=True, merge_last_bin_cutoff=0.5, normalize_last_bin=True)
        all_datasets_and_info.append((dataset, density_plot, color, name, max_per_base, max_per_bin))

    ### plot all the datasets
    all_heatmaps = []
    for N,(dataset,density_plot,color,name,max_per_base,max_per_bin) in enumerate(all_datasets_and_info):
        # for heatmaps, decide what to use as the maximum value to scale the colorbar: 
        #  if max_per_base was provided, use that, otherwise use the actual maximum bin-count
        #  (need to use maximum bin-count over all chromosomes, so that they're all normalized to the same value!)
        if density_plot:
            real_maxcount = max([max(dataset[c]) for c in all_chromosomes if c in dataset])
            if max_per_base is not None and max_per_bin is not None:  
                raise Exception("At least one of max_per_base and max_per_bin for sample name must be None!")
            if max_per_base is not None:    max_count = max_per_base*bin_size
            elif max_per_bin is not None:   max_count = max_per_bin
            else:                           max_count = real_maxcount
            real_mincount = min([min(dataset[c]) for c in all_chromosomes if c in dataset])

        # plot data for each chromosome
        for chr_N,chromosome in enumerate(all_chromosomes):
            chromosome_length = chromosome_lengths[chromosome]
            left_pos, right_pos = _get_plotline_pos(chr_N, total_plotline_width, widths, len(all_datasets_and_info), N)
            # if we're not doing a density plot, just plot the positions as lines
            if not density_plot:
                mplt.vlines(chr_N, left_pos, chromosome_length)
                position_list = dataset[chromosome]
                # only plot the lines if there are any mutants
                if position_list:
                    # only give a label to one chromosome per dataset, so only one shows up in the legend
                    mplt.hlines(position_list, left_pos, right_pos, color=color, 
                                label = ('__nolegend__' if (chr_N!=0 or name=='') else name))
            # if we're doing a density plot, bin the positions into counts, and plot as a heatmap.
            else:
                # inverting the bin value list so that 0 is at the bottom, since otherwise the scale makes no sense;
                #  if no data for this chromosome, no need to plot anything.
                try:                inverted_data = numpy.array(list(reversed(dataset[chromosome])))
                except KeyError:    continue
                # grab the desired colormap, and set it to use NaN_color to display missing data
                cmap = mplt.get_cmap(color)
                cmap.set_bad(NaN_color)
                # actually draw the heatmap image - this is tricky! Matplotlib tends to assume there's only a single heatmap.
                #  - the reshape call is to make the densities from a 1D array into a 2D Nx1 array, for imshow to work properly
                #  - aspect='auto' is to prevent matplotlib from reshaping the entire plot to fit the image shape
                #  - norm=None, vmin, vmax are to prevent matplotlib from normalizing the values to a 0-1 range!
                #    (which would be BAD, because then each chromosome would be normalized separately, and have a different scale)
                #  - DON'T give those a label, so that mplt.legend() will only work on line-plots - colorbars will be separate
                im = mplt.imshow(inverted_data.reshape(-1,1), extent=(left_pos,right_pos,0,chromosome_length), cmap=cmap, 
                                 aspect='auto', norm=None, vmin=min(0,real_mincount), vmax=max_count, **imshow_kwargs)
                # save info image to add colorbars: image, name, and if_normalized (True if max_per_base was given) 
                if chr_N==0 and name!='':  all_heatmaps.append((im, name, (max_per_base is not None), real_maxcount))

    ### set plot limits, ticks, labels, etc
    # it's important to do all this BEFORE colorbars, otherwise it'll apply to the colorbar instead of main axes!
    if title:  mplt.title(title)
    # mplt.imshow has an annoying tendency to reset the plot limits to match a single image, so set them sensibly by hand
    # we want the lower y limit to be slightly below 0, so that the first bin doesn't hide behind the axis line (checked by hand)
    edge_space = 1-total_plotline_width
    mplt.xlim(0 - total_plotline_width/2 - edge_space, len(all_chromosomes)-1 + total_plotline_width/2 + edge_space)
    mplt.ylim(-max(chromosome_lengths.values())*.0015, max(chromosome_lengths.values())*1.01)
    # add xticks with chromosome names as labels!  
    #  - if we only have normal chromosomes, display just the numbers.
    #  - if we have different kinds, need to display the full names - rotate them to keep them from overlapping.
    # MAYBE-TODO add the total number of mutants in each chromosome, and the #mutants/length
    mplt.xlabel('chromosome')
    if not any([include_scaffolds, include_cassette, include_other]):
        mplt.xticks(range(len(all_chromosomes)), [x.split('_')[1] for x in all_chromosomes])
    else:
        mplt.xticks(range(len(all_chromosomes)), all_chromosomes, rotation=90)
    # MAYBE-TODO it'd be nice to get rid of the actual ticks and just keep the labels, the ticks only obscure the data
    # modify yticks to be in Mb - and for some reason this resets ylim, so save it and set it back to previous value
    mplt.ylabel('chromosome position (in Mb) (chromosomes start at the bottom)')
    ylim = mplt.ylim()
    yticks = mplt.yticks()[0]
    mplt.yticks(yticks, [x/1000000 for x in yticks])
    mplt.ylim(ylim)
    plotting_utilities.remove_half_frame()

    ### add legends
    # normal legend for the line-plots (only if there's more than one dataset - if there's just one, presumably the title says it)
    if len(all_datasets_and_info)>1:
        mplt.legend(prop=FontProperties(size='small'))
        # MAYBE-TODO nicer legend formatting?  Get rid of frame, thicker line, line up with colorbars or something?
    # colorbars, if desired
    if not no_colorbars and len(all_heatmaps)>0:
        # if desired, figure out sensible positioning to put smaller colorbars in two rows rather than have the default big ones
        # MAYBE-TODO add option for more than two rows?
        if condense_colorbars:
            # column sizes
            N_columns = numpy.ceil(len(all_heatmaps)/2)
            col_width = 0.1
            col_padding = 0.02
            # bbox is the main figure shape - has attributes like xmin,xmax,ymin,ymax,height,width
            ax = mplt.gca()
            bbox = ax.get_position()
            # set new axes position: decrease width to make room for colorbars, keep left/bottom/height as before (from bbox)
            ax.set_position([bbox.xmin, bbox.ymin, bbox.width - (col_padding+col_width)*(N_columns-1), bbox.height])
            mplt.draw()
            bbox = ax.get_position()
            # row positions/sizes
            row_padding = 0.1*bbox.height
            row_height = (bbox.height - row_padding) / 2
            row_bottoms = [bbox.ymin + bbox.height/2 + row_padding/2, bbox.ymin]
        # add colorbars for each heatmap, labeled with the dataset name
        for N, (plot, name, if_normalized, real_maxcount) in enumerate(all_heatmaps):
            colorbar_kwargs = {}
            # optionally condense the colorbars: put in two rows, reduce spacing
            if condense_colorbars:
                row, column = N % 2, N // 2
                col1_offset = 1 - N_columns * .05
                # add_axes arguments are left,bottom,width,height
                cax = mplt.gcf().add_axes((bbox.xmax + col_padding + col_width*column, row_bottoms[row], 0.012, row_height))
                colorbar_kwargs['cax'] = cax
            c = mplt.colorbar(plot, **colorbar_kwargs)
            c.set_label("%s"%name + ("" if if_normalized else " (per %s)"%(basic_seq_utilities.format_base_distance(bin_size))))
            # MAYBE-TODO put the heatmaps on the plot instead of on the side?  There's room...
            # if if_normalized, scale the ticks to 0-100%
            if if_normalized:
                # c.vmin and c.vmax are the low/high of the colorbar
                assert c.vmin==0, "Colorbar start isn't at 0, can't normalize to 0-100% properly!"
                ticks_and_labels = [(c.vmax*fraction, "%d%%"%(fraction*100)) for fraction in (0, 0.25, 0.5, 0.75, 1)]
            # otherwise, we don't need as many colorbar ticks as matplotlib likes - leave only the full-number ones
            else:
                ticks_and_labels = []
                # c._ticker()[1] gives the tick positions
                for x in c._ticker()[1]:
                    # replace gets around a bug where matplotlib can't render the unicode minus sign in the original x version
                    x = float(x.replace(u'\u2212', '-'))
                    if x==int(x):   ticks_and_labels.append( (x, str(int(x))) )
            # if real_maxcount is higher than the colorbar max, add colorbar ticklabel or change last ticklabel to reflect that
            colorbar_max, highest_tick_pos = c.vmax, ticks_and_labels[-1][0]
            if real_maxcount > colorbar_max:
                if abs(colorbar_max - highest_tick_pos) < 1:
                    ticks_and_labels[-1] = (highest_tick_pos, '%i-%i'%(highest_tick_pos, real_maxcount))
                else:
                    ticks_and_labels.append((int(colorbar_max), '%i-%i'%(int(colorbar_max), real_maxcount)))
            # MAYBE-TODO also add a pointed end to the colorbar if that's the case, or whatever's normally used to signify truncation?  Using the extend arg to mplt.colorbar, so I'd have to do that up above
        c.set_ticks(zip(*ticks_and_labels)[0])
        c.set_ticklabels(zip(*ticks_and_labels)[1])
        mplt.draw()
    # TODO it'd be nice if this actually returned the main axes object... (and maybe a list of colorbar axes too?)  And the numpy.histogram results for all the datasets might be nice too!


def correct_mutant_density_for_mappability(dataset_histogram, mappability_histogram, 
                                           difference_not_ratio=False, mappability_fraction_NaN=0):
    """ Given mutant density and mappability as by-chromosome histograms, return the ratio or difference.

    The first two inputs should be chromosome:bin_density_list dicts, like generated by get_histogram_data_from_positions.
    The keys and value lengths must be the same.  The output will be in the same format, with same keys and value lengths.

    If difference_not_ratio is True, the return dataset will be the difference between the two (can be negative);
    if False, it'll be the ratio of the two inputs, with NaN when mappability is below mappability_fraction_NaN.

    Raise an exception if there are mutants in any zero-mappability range.
    """
    # complain if dataset has chromosomes that mappability doesn't; only warn otherwise
    dataset_chroms = set(dataset_histogram.keys())
    mappability_chroms = set(mappability_histogram.keys())
    if not dataset_chroms.issubset(mappability_chroms):
        raise Exception("the two inputs have different chromosomes!")
    if not dataset_chroms==mappability_chroms:
        print "Warning: some chromosomes present in mappability but not dataset: %s"%(', '.join(mappability_chroms - dataset_chroms))
    corrected_dataset = {}
    NaN = float('NaN')
    # get multiplier to convert the mappability to the expected number of mutants, 
    #  based on the total number of mutants in dataset divided by total mappability
    #   (need to convert everything from arrays to lists to make a sum)
    expected_mutants_multiplier = ( sum(general_utilities.flatten_lists(dataset_histogram)) 
                                   / sum(general_utilities.flatten_lists(mappability_histogram)) )
    # get the number of expected mutants that corresponds to mappability_fraction_NaN of max mappability
    max_mappability = max(general_utilities.flatten_lists(mappability_histogram))
    NaN_cutoff = max_mappability * mappability_fraction_NaN * expected_mutants_multiplier
    for chromosome in dataset_histogram:
        corrected_dataset[chromosome] = []
        expected_mutants_chr_list = mappability_histogram[chromosome] * expected_mutants_multiplier
        if not len(dataset_histogram[chromosome]) == len(expected_mutants_chr_list):
            raise Exception("Mismatched lengths between the two inputs, on %s!"%chromosome)
        for N_mutants, expected_mutants in zip(dataset_histogram[chromosome], expected_mutants_chr_list):
            if expected_mutants==0 and N_mutants>0:
                raise Exception("Mutants found in zero-mappability area on %s!"%chromosome)
            if difference_not_ratio:                corrected_dataset[chromosome].append(N_mutants - expected_mutants)
            else:
                if expected_mutants<=NaN_cutoff:    corrected_dataset[chromosome].append(NaN)
                else:                               corrected_dataset[chromosome].append(N_mutants / expected_mutants)
    return corrected_dataset
    # TODO unit-test!


######### Hotspot/coldspot related plots

def _lowest_cutoff_matched(val, cutoff_list):
    """ Given x and a list of cutoffs, return how many of them x is lower than. """
    return sum(val<cutoff for cutoff in cutoff_list)


def get_hot_cold_spot_colors(pvalue_data, ifcold_data, window_size, window_offset, pval_cutoffs=[0.05, 0.01, 0.005, 0.001], 
                             min_zero=False, print_info=True):
    """ Transform p-value and ifcold (low/high) chrom:value_list dicts into plottable chrom:value lists.

    Both arguments should be chromosome:list_of_window_values, with any number of windows, with given window size and offset. 
      - the p-values should be 0-1, and you should probably get FDR correction done on them first; 
      - the ifcold values should be True for coldspots and False for hotspots. 
    These can be generated by mutant_simulations.find_hot_cold_spots (return values 1 and 3)

    Pval_cutoffs should be a list of 0-1 values of any length, sorted reverse.
    The output will be a chrom:heatmap_value_list, with the value 0 for p-values over pval_cutoffs[1], 
     1 or -1 for p-values under pval_cutoffs[1] but over [2], and so on, with the maximum being len(pval_cutoffs), 
     and the minimum the negative of that.  Negative values indicate coldspots, and positive hotspots. 
    If min_zero is True, len(pval_cutoffs) will be added to each value, to make the lowest possible value 0
     (can be easier for plotting).
    Suitable for plotting with mutant_positions_and_data, format 2, using a blue-white-red colormap or similar, 
     IF the offset and the remainder at the end of the chromosome are relatively small!  Otherwise may end up messy.

    If print_info, print a line giving the number of somewhat significant results, just to give a general idea.
    """
    # TODO modify mutant_positions_and_data to deal with negative values, and to customize colorbar ticklabels for p-values! Add that to docstring both there and here when done
    # TODO implement refactoring the whole list to deal with offsets so it gets plotted correctly, if we see anything worth plotting (right now the positions for different offsets are the same, and the missing chunks at chromosome starts/ends are ignored)
    # TODO is this even used??
    format_bp = lambda x: basic_seq_utilities.format_base_distance(x, False)    # the False is to not approximate
    color_value_dict = {}
    N_significant = 0
    for chrom, pvalues in pvalue_data.items():
        ifcold_vals = ifcold_data[chrom]
        color_values = []
        for N,(pvalue,ifcold) in enumerate(zip(pvalues,ifcold_vals)):
            curr_val = _lowest_cutoff_matched(pvalue, pval_cutoffs)
            if curr_val>0 and print_info:
                N_significant += 1
            if ifcold:      curr_val *= -1
            if min_zero:    curr_val += len(pval_cutoffs)
            color_values.append(curr_val)
        color_value_dict[chrom] = color_values
    if print_info:
        print "%s results below adjusted p-value %s for data with %s window (offset %s)"%(N_significant, pval_cutoffs[0], 
                                                                               format_bp(window_size), format_bp(window_offset))
    return color_value_dict
    # TODO should unit-test this!


def plot_hot_cold_spot_offset(hc_spot_list, pval_cutoffs, all_chromosomes=None, min_offset=0, max_offset=0.4, 
                              linewidth=2, min_height=10000, pval_cutoff_colors=None, chrom_numbers=None):
    """ Plot hot/cold spots as horizontal lines vertically offset by window size/offset, colored by pvalue by chromosome position. 
    
    Hc_spot_list should be a list of (chrom, start, end, pvalue, kind, window_offset, N_mutants, expected) tuples 
     as generated by mutant_simulations.get_hot_cold_spot_list (with kind being 'hotspot' or 'coldspot', 
     and window_offset along with start-end size used to separate potentially overlapping data into lines),

    The offsets will ensure that none of the lines that overlap vertically will completely overlap horizontally, 
     so you can still tell them apart (the smallest window sizes will have the smallest offsets, to make it easier to see).

    Pval_cutoffs should be a reverse-sorted list of cutoffs, 0-1, to be used when coloring the data. 
    Pval_cutoff_colors should be a kind:color_list dict, with each color_list the length of pval_cutoffs, 
     or default if None - so hotspots will be shown in one of the pval_cutoff_colors['hotspot'] colors depending on pvalue, 
      and coldspots likewise in one of the pval_cutoff_colors['coldspot'] colors.

    Min/max offset give the x-axis range relative to the middle of each chromosome over which to spread all the data
     (1 is the distance between two chromosomes on the x axis, so the full range would be about -.4 to .4)
    Min_height is the minimum height 

    All_chromosomes can be a list of all chromosomes to include on the plot, 
     in case you want to leave space for any that DON'T have any hc_spot_list entries.
    Chrom_numbers should be a chromosome_name:number dictionary to give each chromosome a unique position on the x axis; 
     if None, we'll just sort them all sensibly and assign numbers that way. 
    """
    if pval_cutoff_colors is None:
        if len(pval_cutoffs) == 3:
            pval_cutoff_colors = {'hotspot': 'salmon red darkred'.split(), 'coldspot': 'skyblue steelblue darkblue'.split()}
        else:
            raise Exception("If pval_cutoffs isn't length 3, must provide pval_cutoff_colors, no default!")
    if all_chromosomes is None:
        all_chromosomes = list(set(hc_data[0] for hc_data in hc_spot_list))
    else:
        all_chromosomes = list(all_chromosomes)
    if chrom_numbers is None:
        all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)
        chrom_numbers = {chrom:N for (N,chrom) in enumerate(all_chromosomes)}
    # figure out how many potentially overlapping (window_size, window_offset) sets we have, 
    #  give each of them a separate x axis position
    all_lines = sorted(set((end-start, offset) for (chrom, start, end, pvalue, kind, offset, mcount, expected) in hc_spot_list))
    x_offset_range = max_offset-min_offset
    x_offset_per_line = x_offset_range/len(all_lines)
    x_line_offsets = {line: (min_offset)+x_offset_per_line*N for (N,line) in enumerate(all_lines)}
    for (chrom, start, end, pvalue, kind, offset, mcount, expected) in hc_spot_list:
        # plot only the lines with a pvalue that matches at least the highest cutoff
        cutoff_matched = _lowest_cutoff_matched(pvalue, pval_cutoffs)
        if cutoff_matched>0:
            x_offset = x_line_offsets[(end-start,offset)]
            # if min_height is specified, make sure the line height is at least that
            if min_height is not None:
                height = end-start
                missing_height = min_height - height
                if missing_height>0:
                    start, end = start - missing_height/2, end + missing_height/2
            # line color depends on lowest pvalue cutoff matched
            color = pval_cutoff_colors[kind][cutoff_matched-1]
            mplt.vlines(chrom_numbers[chrom]+x_offset, start, end, color=color, linewidth=linewidth)


def plot_hot_cold_spot_overlapping(hc_spot_list, pval_cutoffs, all_chromosomes=None, linewidth=2, x_offset=0.2, 
                                   pval_cutoff_colors=None, sort_by=None, chrom_numbers=None):
    """ Plot hot/cold spots as overlapping horizontal lines with rounded ends, colored by pvalue by chromosome position. 
    
    Hc_spot_list should be a list of (chrom, start, end, pvalue, kind, window_offset, N_mutants, expected) tuples 
     as generated by mutant_simulations.get_hot_cold_spot_list (with kind being 'hotspot' or 'coldspot') 

    Pval_cutoffs should be a reverse-sorted list of cutoffs, 0-1, to be used when coloring the data. 
    Pval_cutoff_colors should be a kind:color_list dict, with each color_list the length of pval_cutoffs, 
     or default if None - so hotspots will be shown in one of the pval_cutoff_colors['hotspot'] colors depending on pvalue, 
      and coldspots likewise in one of the pval_cutoff_colors['coldspot'] colors.

    X_offset gives the x-axis position of the lines relative to the middle of each chromosome 
     (1 is the distance between two chromosomes).

    All the lines will overlap, so if hc_spot_list contains two lines out of which one contains the other, 
     you may not be able to tell them apart - you may want to hand-filter the data passed to this function 
      to only show the "important" lines, if a lot overlap.
    Use the sort_by argument to customize overlapping line issues by plotting lines in a given order (so the last ones are on top):
        - None: don't sort, just use the order as given
        - 'pvalue': lines will be plotted in reverse p-value order, i.e. lowest p-value on top, which usually works sensibly.
        - 'size': lines will be plotted in reverse size order, i.e. smallest on top.

    All_chromosomes can be a list of all chromosomes to include on the plot, 
     in case you want to leave space for any that DON'T have any hc_spot_list entries.
    Chrom_numbers should be a chromosome_name:number dictionary to give each chromosome a unique position on the x axis; 
     if None, we'll just sort them all sensibly and assign numbers that way. 
    """
    if pval_cutoff_colors is None:
        if len(pval_cutoffs) == 3:
            pval_cutoff_colors = {'hotspot': 'salmon red darkred'.split(), 'coldspot': 'skyblue steelblue blue'.split()}
        elif len(pval_cutoffs) == 2:
            pval_cutoff_colors = {'hotspot': 'lightcoral red'.split(), 'coldspot': 'deepskyblue blue'.split()}
        else:
            raise Exception("If pval_cutoffs isn't length 2-3, must provide pval_cutoff_colors, no default!")
    if all_chromosomes is None:
        all_chromosomes = set(hc_data[0] for hc_data in hc_spot_list)
    else:
        all_chromosomes = set(all_chromosomes)
    if chrom_numbers is None:
        all_chromosomes_sorted = list(all_chromosomes)
        all_chromosomes_sorted.sort(key=basic_seq_utilities.chromosome_sort_key)
        chrom_numbers = {chrom:N for (N,chrom) in enumerate(all_chromosomes_sorted)}
    # plot the lines reverse-sorted by p-value, so the lowest p-values are plotted last and thus are on top.
    if sort_by is None:         sorted_hc_spot_list = list(hc_spot_list)
    elif sort_by=='pvalue':     sorted_hc_spot_list = sorted(hc_spot_list, key=lambda x: x[3], reverse=True)
    elif sort_by=='size':       sorted_hc_spot_list = sorted(hc_spot_list, key=lambda x: x[2]-x[1], reverse=True)
    for (chrom, start, end, pvalue, kind, offset, mcount, expected) in sorted_hc_spot_list:
        # plot only ones on all_chromosomes
        if chrom not in all_chromosomes:    continue
        # plot only the lines with a pvalue that matches at least the highest cutoff
        cutoff_matched = _lowest_cutoff_matched(pvalue, pval_cutoffs)
        if cutoff_matched>0:
            # line color depends on lowest pvalue cutoff matched
            color = pval_cutoff_colors[kind][cutoff_matched-1]
            # plot a line, with round caps (so a really short line will be sort of a circle)
            x_pos = chrom_numbers[chrom] + x_offset
            mplt.plot([x_pos, x_pos], [start, end], color=color, linewidth=linewidth, markersize=linewidth, 
                      solid_capstyle='round', marker='o', markeredgewidth=0)
    # MAYBe-TODO refactor this with plot_hot_cold_spot_offset?  They're pretty similar...


######### Gap size plots

def get_gap_sizes(dataset):
    """ Given a dataset, get a list of all the gap sizes between mutant positions, UNSORTED. 
    
    Go over each chromosome, take each pair of adjacent mutants (regardless of strand), get gap size.

    Dataset can be either a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a chrom:position_list dict.
    """
    # convert dataset to chrom:position_list dict if needed
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        dataset = get_mutant_positions_from_dataset(dataset, strand=None)
    all_gap_sizes = []
    # go over each chromosome and give the gaps between each adjacent pair of mutants
    for chrom,pos_list in dataset.items():
        pos_list.sort()
        for pos1,pos2 in zip(pos_list, pos_list[1:]):
            all_gap_sizes.append(pos2-pos1)
    return all_gap_sizes


def gap_size_QQ_plot(reference_dataset, simulated_datasets=None, 
                     N_simulations=None, mappable_positions_20bp=None, mappable_positions_21bp=None, 
                     logscale=True, different_colors=False, markersize=6, plot_padding=(0.1, 0.1)):
    """ Quantile-quantile plot of gap sizes between mutants in real vs simulated reference_dataset (repeated N times). 

    Provide either a list of simulated_datasets (chrom:pos_list dicts), OR information to generate them:
     N_simulations, and the two mappable_positions_* args (output of mutant_simulations.genome_mappable_insertion_sites).
    
    """
    reference_gap_sizes = get_gap_sizes(reference_dataset)
    # Problem: since the genome's separated into chromosomes, simulated datasets with the same number of mutants
    #   can have different numbers of gap sizes!  2 mutants on chr1 and chr2 - 2 gaps; 4 and 0 - 3 gaps. 
    #  This difference can't be more than the number of chromosomes, and likely much less. 
    #  So just remove 10 random gap sizes just in case.  MAYBE-TODO come up with a better solution?
    random.shuffle(reference_gap_sizes)
    reference_gap_sizes = sorted(reference_gap_sizes[:-10])
    max_gap = max(reference_gap_sizes)
    # if simulated datasets were provided, just use them; otherwise generate data to do the simulations
    if simulated_datasets is not None:
        N_simulations = len(simulated_datasets)
    else:
        N_mutants = len(reference_dataset)
        fraction_20bp = mutant_simulations.get_20bp_fraction(reference_dataset)
        # TODO mutant_simulations.get_20bp_fraction only works on a full mutant dataset, not on a chrom:pos_list dict - either require reference_dataset to be real, or add an option to give that separately
    if different_colors:    plot_kwargs = {}
    else:                   plot_kwargs = {'color': 'blue'}
    for N in range(N_simulations):
        # if simulated datasets were provided, use the next one;
        #  otherwise make a new one for each repetition (and discard it afterward, to avoid taking up memory)
        if simulated_datasets is not None:
            simulated_dataset = simulated_datasets[N]
        else:
            simulated_dataset = mutant_simulations.simulate_dataset_mappability(N_mutants, fraction_20bp, 
                                                                                mappable_positions_20bp, mappable_positions_21bp)
        simulated_gap_sizes = get_gap_sizes(simulated_dataset)
        # make sure the number of simulated gap sizes matches the number of reference ones
        #   (see "Problem" comment section up top for why this can happen)
        if len(simulated_gap_sizes) < len(reference_gap_sizes):
            raise Exception("Fewer simulated than real gap sizes! Remove more from the real ones?")
        random.shuffle(simulated_gap_sizes)
        simulated_gap_sizes = sorted(simulated_gap_sizes[:len(reference_gap_sizes)])
        max_gap = max(max_gap, max(simulated_gap_sizes))
        if N==0:    label = "real vs simulated,\nrepeated %s times"%N_simulations
        else:       label = '__nolegend__'
        mplt.plot(reference_gap_sizes, simulated_gap_sizes, '.', markeredgecolor='none', label=label, markersize=markersize, **plot_kwargs)
    mplt.title("Real vs simulated gap sizes between insertions (%s scale),"%('log' if logscale else 'linear')
               +"\nquantile-quantile plot (comparing sorted gap-size lists)"
               +"\n(simulation taking mappability into account, repeated %s times)"%N_simulations)
    mplt.xlabel('Real gap sizes')
    mplt.ylabel('Simulated gap sizes')
    if logscale:
        # use symlog instead of log to deal with gap size 0 - the range between -1 and 1 will be linear
        mplt.xscale('symlog',linthreshx=1)
        mplt.yscale('symlog',linthreshy=1)
        # MAYBE-TODO are those sensible edges?  
        #  Could come up with some way of taking plot_padding[0] into account, but it's tricky on a symlog plot
        plot_edges = (-0.5, 10**((1+plot_padding[1])*math.log10(max_gap)) )
    else:
        plot_edges = (-plot_padding[0]*max_gap, (1+plot_padding[1])*max_gap)
    # make the plot square even if the data isn't; plot a diagonal line to show what identical datasets would look like
    mplt.plot(plot_edges, plot_edges, c='grey', label='if both were identical')
    mplt.xlim(plot_edges)
    mplt.ylim(plot_edges)
    mplt.legend(loc=2)
    # TODO those plots should be square!  How do I square an axes instance?


######### Chromosome mutant density plots

def chromosome_density_scatterplot(mutant_dataset, include_scaffolds=True, include_cassette=True, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a chromosome length vs mutant number scatterplot. 

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.

    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).

    The *_multiplier arguments can be set to N>1 to make the plot reflect that fact that there are multiple copies 
     of the chloroplast/mitochondrial genomes in the cell, so the "functional" length is N times higher.
    """
    # TODO add option to only count effective (mappable) length!

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want;
    #  apply chloro/mito multipliers
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    chromosomes_by_type = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, 
                                                   output_dict=True)

    mutants_in_chromosomes_all = {}
    for chr_type in sorted(chromosomes_by_type, key=basic_seq_utilities.chromosome_sort_key):
        chr_list = chromosomes_by_type[chr_type]
        mutants_in_chromosomes = [mutant_dataset.summary.mutants_in_chromosome(c) for c in chr_list]
        mutants_in_chromosomes_all.update(dict(zip(chr_list, mutants_in_chromosomes)))
        if sum(mutants_in_chromosomes):
            mplt.plot([chromosome_lengths[c] for c in chr_list], mutants_in_chromosomes, 'o', 
                      color=basic_seq_utilities.CHROMOSOME_TYPE_COLORS[chr_type], label=chr_type)

    max_length = max(chromosome_lengths.values())
    max_N_mutants = max(mutants_in_chromosomes_all.values())
    mplt.xlabel('chromosome length (in Mb)')
    mplt.xticks(mplt.xticks()[0], [x/1000000 for x in mplt.xticks()[0]])
    mplt.xlim(-max_length*.05, max_length*1.1)
    mplt.ylabel('number of mutants in chromosome')
    mplt.ylim(-max_N_mutants*.05, max_N_mutants*1.1)
    mplt.legend(loc='lower right', prop=FontProperties(size='medium'))


def chromosome_density_barchart(mutant_dataset, include_scaffolds=True, include_cassette=False, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a simple bar-chart of chromosome mutant densities (per kb).

    See chromosome_density_scatterplot docstring for all the arguments - they're the same.
    """
    # TODO add option to only count effective (mappable) length!

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want; 
    #  apply chloro/mito multipliers; sort
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)

    # calculate and plot mutant densities, colored by chromosome type
    mutant_densities = [mutant_dataset.summary.mutants_in_chromosome(c)/chromosome_lengths[c]*1000 for c in all_chromosomes]

    mplt.bar(range(len(all_chromosomes)), mutant_densities, 
             color=[basic_seq_utilities.chromosome_color(c) for c in all_chromosomes], align='center', linewidth=0)
    # MAYBE-TODO add error bars based on total mutant number?
    # MAYBE-TODO add a line/something showing the max mutants/20kb value for each chromosome?  But that would require more data processing (similar to what was done in mutant_positions_and_data)

    mplt.ylabel('mutants per kb')
    mplt.xticks(range(len(all_chromosomes)), all_chromosomes, rotation=90)
    mplt.xlim(-1, len(all_chromosomes))


######################################### Plotting gene-related data ###########################################

### number of genes with 1+/2+/etc mutants vs number of mutants (randomly chosen mutant subsets)

def genes_with_N_mutants(dataset, subset_sizes=100, max_N_mutants=3, N_mutants_colors=None, repeat_N_times=100, dataset_info='', 
                         total_genes=None, mappable_percent=None, plot_full_count=True, 
                         print_repeat_progress=10, rasterize_data=False):
    """ Plot % of all genes with N mutants for different-sized random subsets of dataset.

    Plot % of all genes that have between 1 and max_N_mutants mutants (separate line for each N); 
     N_mutants_colors should be a list of length max_N_mutants, giving the colors to plot the genes with each N mutants.
    Total_genes is the total number of genes in the genome; if None, it'll be extracted from dataset if possible, 
     or else the default value from the Phytozome chlamy v4.3 genome will be used (17114).

    Plot points for mutant subset sizes between 0 and almost the full dataset size, defined by the subset_sizes argument
      (can be either a list of subset sizes, or a number which will be used as the increment to a range); the last point 
      won't be the full dataset size unless it's included in subset_sizes (if list), or divisible by subset_sizes (if number).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).

    Redo repeat_N_times, with different random subsets each time, to see the full representation.
    If plot_full_count is True, also plot dots for the full dataset (in case it's not there due to subset_sizes value), 
     in a different style (black borders around the dots).

    Mappable_percent is the approximate % of mutants that are mappable, just to put in the xlabel; if None, don't print it. 
    """
    # TODO update to show simulated data in parallel to real data, once we have simulated data!  (Simulated data can go to higher total mutant numbers than we have in dataset)  Simulated lines should be paler colors or different markers or something.

    # TODO implement a log x-scale!  And make sure subset_sizes is log-scaled too, in that case, if it was a range to start with... Can just do it by hand and use custom subset_sizes to get around that for now.

    # nice default color schemes only defined for some max_N_mutants values
    if max_N_mutants<=3 and N_mutants_colors is None:
        N_mutants_colors = 'magenta darkcyan darkorange'.split()

    # default total_genes: extract from dataset if possible; 
    #  if not (dataset missing that attribute, dataset is a mutant list, or dataset value is None as well), use hardcoded default. 
    if total_genes is None:
        try:                    
            total_genes = dataset.total_genes_in_genome
        except AttributeError:  
            raise Exception("Can't determine the total number of genes from the dataset - provide as an argument!")

    # if subset_sizes is a number, make it a list by using it as a range increment
    if not hasattr(subset_sizes, '__len__'):
        subset_sizes = range(0, len(dataset), subset_sizes)
    else:
        subset_sizes.sort()
        subset_sizes = [x for x in subset_sizes if x <= len(dataset)]
    # Plot it all several times with new random mutant subsets, to make sure we have a good coverage of the random space
    for repeat in range(repeat_N_times):
        if print_repeat_progress:
            if not repeat % print_repeat_progress:  
                print "repeat %s/%s..."%(repeat, repeat_N_times)
        gene_counts_by_Nmutants = mutant_simulations.gene_counts_for_mutant_subsets(dataset, subset_sizes, max_N_mutants)
        for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
            # use custom colors if given, otherwise let mplt choose colors
            plot_kwargs = {} 
            if N_mutants_colors:    plot_kwargs['color'] = N_mutants_colors[N_mutants-1]
            if rasterize_data:      plot_kwargs['rasterized'] = True
            # only label the lines in the first repeat, to avoid 100 repeats in the legend!
            if repeat==0 and dataset_info != '__nolegend__':
                plot_kwargs['label'] = "genes with %s+ mutants%s"%(N_mutants, ', %s'%dataset_info if dataset_info else '')
            mplt.plot(subset_sizes, gene_counts, '.', linewidth=0, **plot_kwargs)
    # optionally add a separate plot set for the actual full mutant count, with different style (black border)
    if plot_full_count:
        gene_counts_by_Nmutants = mutant_simulations.gene_counts_for_mutant_subsets(dataset, [len(dataset)], max_N_mutants)
        for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
            plot_kwargs = {} 
            if N_mutants_colors:    plot_kwargs['color'] = N_mutants_colors[N_mutants-1]
            if rasterize_data:      plot_kwargs['rasterized'] = True
            mplt.plot(len(dataset), gene_counts, '.', linewidth=0, markeredgecolor='None',**plot_kwargs)
    # add a line at the total gene number - TODO this doesn't work right... figure out sensible xmin/xmax values, reset xlim
    mplt.legend(loc='upper left', prop=FontProperties(size='medium'))

    mplt.title('Percent of genes hit vs number of mutants sequenced\n(plotted for %s random mutant subsets)'%repeat_N_times)
    mplt.ylabel("Percent of genes hit (out of %s total chlamy nuclear genes)"%total_genes)
    # Y max should never be more than 100%, it makes the plot look weird
    if mplt.ylim()[1] > 0.8*total_genes:
        mplt.ylim(-0.01*total_genes, 1.01*total_genes)
        percent_ticks = range(0,101,10)
        mplt.yticks([x/100*total_genes for x in percent_ticks])
    mplt.yticks(mplt.yticks()[0], ["%.0f%%"%(x*100/total_genes) for x in mplt.yticks()[0]])
    mplt.xlabel("Number of mutants (randomly chosen subsets out of %s total)\n"%len(dataset) 
                +"(counting only mappable unique genomic mutants%s)"%(' - about %s%%'%mappable_percent if mappable_percent else ''))


### histogram of the number of insertions distributed over normalized gene length, for all genes combined

ORIENTATION_COLORS = {'both': 'blue', 'sense': 'green', 'antisense': 'magenta'}

def get_binned_mutant_counts(dataset, gene_positions, N_bins, orientation):
    """ Get the number of mutants for each gene split into N bins along gene length.

    Returns a gene:array_of_bin_mutant_counts dict giving the number of mutants in dataset for each gene in gene_positions, 
    split into N_bins along the length of the gene.  Dict values are 1D numpy arrays, not lists.

    If orientation is 'both', use all mutants; if 'sense'/'antisense', use only mutants with that orientation w.r.t. the gene.

    Treats sense/antisense genes correctly: if a gene is antisense position 100-200, a mutant at 101 goes in the last bin.
    """
    gene_binned_mutants = {}
    for gene,(strand,start,end) in gene_positions.iteritems():
        if orientation=='both': mutant_positions = [m.position.min_position-start for m in dataset.mutants_by_gene[gene]]
        else:                   mutant_positions = [m.position.min_position-start for m in dataset.mutants_by_gene[gene]
                                                   if m.orientation==orientation]
        if strand=='-':
            mutant_positions = [(end-start)-position for position in mutant_positions]
        gene_binned_mutants[gene] = numpy.histogram(mutant_positions, bins=N_bins, range=(0,end-start))[0] 
    return gene_binned_mutants


def insertions_over_gene_length(mutant_data, gene_positions=None, gene_effective_lengths_binned=None, N_bins=None, 
                                separate_sense_antisense=True, length_type=None, seq_info=''):
    """ Histogram of number of insertions binned over gene length (optionally normalized to effective length) for all genes combined.

    Mutant_data can be a orientation:(gene:bin_mutant_counts) dict (with values generated by get_binned_mutant_counts, 
       and orientations a subset of (both, sense, antisense) matching the separate_sense_antisense value), 
     or a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, in which case get_binned_mutant_counts will be used
      to get the binned count dict from it (which takes a while).
     If it's a dataset instance, N_bins must be given; otherwise it'll be inferred from the bin mutant count list length.

    If gene_effective_lengths_binned is given, the mutant density per bin will be calculated based on the effective bin lengths
     instead of raw bin lengths (the mutant density per bin will be calculated for each gene, then they'll be averaged together);
     gene_effective_lengths_binned should be a gene:list_of_bin_effective_lengths, generated by mutant_simulations.gene_mappability
      (2nd output) or another similar function - make sure it was run with with N_bins as the value for the split_into_bins arg, 
       and the same exclude_UTRs_in_bins value as used to generate gene_positions if given and mutant_data if it's a dictionary!)

    Gene_positions should be a gene:(start_pos,end_pos) dictionary (giving the mRNA start/end positions, or coding ones, or such)
      (like generated by gff_examine_file.gene_positions with include_chromosomes=False and include_strand=True);
     if mutant_data is already binned and gene_effective_lengths_binned is given, gene_positions isn't needed.

    If separate_sense_antisense is True, separate mutant counts by sense/antisense orientation compared to the gene,
     and plot them as two separate bar sets on the same graph.

    Length_type should be a string describing the type of length in gene_effective_lengths_binned (mappable, or total, or something);
     seq_info should be a string describing the gene sequence type (full, or excluding UTRs, or such) 
      used when generating mutant_data and gene_effective_lengths_binned.  Both will be used for text/labels only.

    The mutant densities per bin for each gene are averaged by giving equal weight to each base, not each gene:
     if a 2kb gene has counts [1,1], and an 18kb gene has counts [0,0], the average will be 0.1 mutant per kb, not 0.5.
    """
    ### process the arguments
    if separate_sense_antisense:    orientations = ['sense', 'antisense']
    else:                           orientations = ['both']
    # make sure gene_positions is given when needed, and NOT given if not needed
    if isinstance(mutant_data, mutant_analysis_classes.Insertional_mutant_pool_dataset) or gene_effective_lengths_binned is None:
        if gene_positions is None:    raise Exception("gene_positions must be provided if mutant_data is a dataset object "
                                                      +"or gene_effective_lengths_binned isn't given!")
    elif gene_positions is not None:  raise Exception("gene_positions SHOULDN'T be provided if mutant_data is an already binned "
                                                      +"dictionary and gene_effective_lengths_binned is given!")
    # convert mutant_data to an orientation:gene:bin_mutant_counts dict if it's a full dataset; 
    #  make sure N_bins is sane, set it if needed
    if isinstance(mutant_data, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        if N_bins is None:
            raise Exception("N_bins must be provided if mutant_data is a dataset object!")
        mutant_data = {orient: get_binned_mutant_counts(mutant_data, gene_positions, N_bins, orient) for orient in orientations}
    else:
        if not set(orientations).issubset(mutant_data.keys()):
            raise Exception("mutant_data provided to insertions_over_gene_length must have %s keys!"%orientations)
        if N_bins is not None:
            raise Exception("N_bins cannot be provided if mutant_data is a dictionary, since it already has bins!")
        # infer N_bins from the number of bins in mutant_data values, make sure it's sane
        N_bins = set.union(*[set([len(l) for l in mutant_data_val.values()]) for mutant_data_val in mutant_data.values()])
        if len(N_bins) > 1:
            raise Exception("mutant_data is malformed - the number of bins varies between genes or orientations! (%s)"%N_bins)
        else:
            N_bins = N_bins.pop()
    # if effective bin lengths are given, use that for bin lengths (and make sure N_bins matches the one from mutant_data or input), 
    #  otherwise just divide gene length by N_bins
    if gene_effective_lengths_binned is not None:
        gene_bin_lengths = gene_effective_lengths_binned
        N_bins_2 = set([len(l) for l in gene_bin_lengths.values()])
        if len(N_bins_2) > 1:
            raise Exception("gene_effective_lengths_binned is malformed - the number of bins varies between genes! (%s)"%N_bins_2)
        else:
            N_bins_2 = N_bins_2.pop()
        if N_bins != N_bins_2:
            raise Exception("Bin numbers for mutant_data and gene_effective_lengths_binned are different! %s, %s"%(N_bins,N_bins_2))
    else:
        gene_bin_lengths = {gene:[(end-start+1)/N_bins for _ in range(N_bins)] for (gene,(_,start,end)) in gene_positions.items()}
    ### get the actual data used for the plot
    bin_mutant_densities_dict = {}
    for orient in orientations:
        # convert mutant counts to mutant densities by dividing them by bin lengths (binned_mutants is a numpy array, easy division)
        gene_mutant_densities = {}
        for gene,binned_mutants in mutant_data[orient].items():
            # if some bin lengths are 0, I can put any value for the 0-length bin, since it'll have weight 0 in the average;
            #   using an absurdly high value to make it more obvious if the average doesn't work right
            gene_mutant_densities[gene] = numpy.array([(1000000 if bin_len==0 else m_count/bin_len) 
                                                       for (m_count,bin_len) in zip(binned_mutants,gene_bin_lengths[gene])])
        # average all the mutant bin densities together, using numpy.average over dimension 0 of the list of all bin density arrays;
        #  give equal weight to each base by weighing genes by bin length (effective or total)
        gene_order = list(gene_mutant_densities.keys())
        bin_mutant_densities_dict[orient] = numpy.average([gene_mutant_densities[g] for g in gene_order], 0, 
                                                          [gene_bin_lengths[g] for g in gene_order])
    ### make the plot
    bar_width = 0.8/len(bin_mutant_densities_dict)
    for N,orientation in enumerate(orientations):
        bin_mutant_densities = bin_mutant_densities_dict[orientation]
        mplt.bar([x + N*bar_width for x in range(N_bins)], bin_mutant_densities, bar_width, label=orientation, 
                 align='center', color=ORIENTATION_COLORS[orientation], edgecolor='none')
    # MAYBE-TODO integrate the p-values from insertions_over_gene_length_pvalues into the plot?  Easier done by hand, for now.
    if len(bin_mutant_densities_dict)>1:
        mplt.legend(loc='lower center', prop=FontProperties(size='medium'))
    mplt.xlabel('gene sequence length%s divided into %s equal sections'%(' (%s)'%seq_info if seq_info else '', N_bins))
    mplt.xticks([], [])
    mplt.xlim(-0.6, N_bins-0.4)
    # specify length type; if not given, set to TOTAL if gene_effective_lengths_binned wasn't given; otherwise we can't know
    if length_type is None:
        if gene_effective_lengths_binned is None:  length_type = ' TOTAL' 
        else:                            length_type = '' 
    else:                                length_type = ' '+length_type
    mplt.ylabel('average mutants per kb\nof%s bin length'%length_type)
    # convert yticks to kb
    mplt.yticks(mplt.yticks()[0], ["%.2g"%(x*1000) for x in mplt.yticks()[0]])


def insertions_over_gene_length_pvalues(mutant_data_by_orientation, gene_effective_lengths_binned, length_type=None, seq_info=''):
    """ Get p-values for number of insertions binned over gene length for all genes combined.

    mutant_data_by_orientation should be a orientation:(gene:bin_mutant_counts) dict 
     with orientations [both,sense,antisense], and values for each orientation generated by get_binned_mutant_counts.

    Gene_effective_lengths_binned should be a gene:list_of_bin_lengths, generated by mutant_simulations.gene_mappability (2nd output)
     or something else similar - make sure the bin number and gene length type (full, or excluding UTRs) 
      were the same when generating that as when generating mutant_data_by_orientation!

    Length_type should be a string describing the type of length in gene_effective_lengths_binned (mappable, or total, or something);
     seq_info should be a string describing the gene sequence type (full, or excluding UTRs, or such) 
      used when generating mutant_data_by_orientation and gene_effective_lengths_binned. Both will be used for text/labels only.

    Do two types of p-value calculations for each bin:
     - for total mutant counts (both orientations), use the binomial distribution to check whether they're significantly different 
       from the overall mutant density (the overall mutant density is the weighed average of intron/CDS mutant densities, 
        NOT the total mutant density - that seems like the most meaningful measure)
     - for sense/antisense, what we really want to know is whether there's a significant difference between sense and antisense 
       (compared to the overall gene-without-UTR sense/antisense proportions), not whether each one is different from 
       the overall mutant density, so use Fisher's exact test to check if the proportions are different.
    For both types of p-values, do false discovery rate (FDR) adjustment with the Benjamini-Hochberg method, 
     using the p_adjust function from the R language via rpy2.
    """
    # calculate probability of a mutant falling into each bin vs in all bins, based on bin effective lengths
    bin_lengths = sum([numpy.array(l) for l in gene_effective_lengths_binned.values()])
    N_bins = len(bin_lengths)
    total_gene_length = sum(bin_lengths)
    bin_probabilities = bin_lengths/total_gene_length
    # get the total number of mutants per bin, summed over all genes
    bin_mutant_counts = {orient: sum(mutant_data_by_orientation[orient].values()) for orient in mutant_data_by_orientation}
    # get the total number of mutants over all bins
    total_N_mutants = {orient: sum(bin_mutant_counts[orient]) for orient in bin_mutant_counts}
    ### get the p-values for whether total bin mutant counts are different from the overall count
    total_count_pvalues = [scipy.stats.binom_test(x=bin_count, n=total_N_mutants['both'], p=bin_probability) 
                           for (bin_count, bin_probability) in zip(bin_mutant_counts['both'], bin_probabilities)]
    total_count_pvalues_adjusted = statistics_utilities.FDR_adjust_pvalues(total_count_pvalues, method='BH')
    # print the data
    print "Whether total bin mutant counts are significantly different from overall, based on %s length:"%length_type
    print " (overall mutant density is %s (over full gene length %s)"%(total_N_mutants['both']/total_gene_length, seq_info) 
    for N, (adj_pval, pval, mutant_count, bin_len) in enumerate(zip(total_count_pvalues_adjusted, total_count_pvalues, 
                                                                   bin_mutant_counts['both'], bin_lengths)):
        print " - %.3g FDR-adjusted p-value (%.3g raw p-value) - bin %s, mutant density %.3g (raw mutant count %s, length %s)"%(
                                                                    adj_pval, pval, N+1, mutant_count/bin_len, mutant_count, bin_len)
    ### get the p-values for whether sense/antisense bin mutant count proportions are different from the overall proportions, 
    orientation_diff_pvalues = []
    for bin_N in range(N_bins):
        fisher_result = scipy.stats.fisher_exact([[bin_mutant_counts['sense'][bin_N], bin_mutant_counts['antisense'][bin_N]], 
                                                 [total_N_mutants['sense'], total_N_mutants['antisense']]])
        orientation_diff_pvalues.append(fisher_result[1])   # fisher output is (oddsratio,pvalue) - we just want pvalue
    orientation_diff_pvalues_adjusted = statistics_utilities.FDR_adjust_pvalues(orientation_diff_pvalues, method='BH')
    # print the data
    print "Whether bin mutant sense/antisense proportions are significantly different from overall:"
    print " (overall mutant sense:antisense counts are %s:%s (over full gene length %s)"%(total_N_mutants['sense'], 
                                                                                          total_N_mutants['antisense'], seq_info) 
    for N, (adj_pval, pval, sense_count, antisense_count, bin_len) in enumerate(zip(orientation_diff_pvalues_adjusted, 
                                orientation_diff_pvalues, bin_mutant_counts['sense'], bin_mutant_counts['antisense'], bin_lengths)):
        print " - %.3g FDR-adjusted p-value (%.3g raw p-value) - bin %s, %s:%s sense:antisense mutants (bin length %s)"%(
                                                                        adj_pval, pval, N+1, sense_count, antisense_count, bin_len)


### bar-chart of the density of insertions in different gene features, for all genes combined

def adjust_feature_mutant_counts(feature_mutant_counts, ignore_UTR_introns=2, feature_boundaries=1, overlapping_genes=1, 
                                 add_all=False, collapse_splice_variants=True, long_UTR_names=False):
    """ Make adjustments to a raw feature mutant count to make it fit desired format for insertion_density_by_feature.

    Changes made:
    - change '-' to 'intergenic'
    - if add_all, add the 'all' category (sum of all other categories - assumes categories are non-overlapping!)
    - add the 'intergenic' category ('all' minus 'gene')
    - if long_UTR_names, change 5'UTR to five_prime_UTR and same for 3'

    For feature-boundary categories: if feature_boundaries is 0, don't count them them except in all/gene categories; 
     if 1, count it as a randomly chosen of the features it's a boundary of (excepting gene_edge/mRNA_edge); 

    For overlapping-gene cases: if overlapping_genes is 0, don't count them except in all/gene categories; 
     if 1, count it as the first feature (useful because when looking at feature mappability we take the first of overlapping gnes);
     if 2, count both features; 
     if 3, count as a randomly chosen one,
     if 4, count as 'multiple'
     if 5, look at the features - if all the same, count as that feature, otherwise 'multiple'

    For multiple splice variants: if collapse_splice_variants is True, count as 'multiple'.
        (if they're the same feature, they would already be collapsed to one, so no need to worry about that)

    If ignore_UTR_introns is 0, just keep "5'UTR_intron" and "3'UTR_intron" counts separate;
     if 1, merge them with "5'UTR" and "3'UTR" respectively (which is probably sensible); 
     if 2, merge them both into "intron" (which is how they're annotated in all gff files, so this is useful for 
        comparisons to other data in which I did the parsing in a simple way.)
    """
    if feature_boundaries not in [0,1]:                 raise Exception("Unknown method feature_boundaries=%s"%feature_boundaries)
    if overlapping_genes not in [0,1,2,3,4,5]:          raise Exception("Unknown method overlapping_genes=%s"%overlapping_genes)
    if ignore_UTR_introns not in [0,1,2]:               raise Exception("Unknown method ignore_UTR_introns=%s"%ignore_UTR_introns)
    if feature_boundaries==1 and not collapse_splice_variants:
        raise Exception("Can't deal with feature boundaries without dealing with splice variants, since one can contain the other!")
    # make a new copy of the input before making any changes
    feature_mutant_counts = dict(feature_mutant_counts)
    # for intergenic, set to 0 if there weren't any (e.g. if the input only counted sense-orientation mutants
    if 'intergenic' not in feature_mutant_counts:
        try:                feature_mutant_counts['intergenic'] = feature_mutant_counts.pop('-')
        except KeyError:    feature_mutant_counts['intergenic'] = 0
    # calculate all/gene categories, add 'multiple' if needed
    all_features = sum(feature_mutant_counts.values())
    feature_mutant_counts['gene'] = all_features - feature_mutant_counts['intergenic']
    if add_all:
        feature_mutant_counts['all'] = all_features
    if overlapping_genes in [4,5] or collapse_splice_variants:
        feature_mutant_counts['multiple'] = 0
    # deal with overlapping genes
    for feature in list(feature_mutant_counts):
        if '&' in feature:  
            features = feature.split(' & ')
            if overlapping_genes==1:
                feature_mutant_counts[features[0]] += feature_mutant_counts.pop(feature)
            elif overlapping_genes==2:
                for f in features:
                    feature_mutant_counts[f] += feature_mutant_counts.pop(feature)
            elif overlapping_genes==3:
                random.shuffle(features)
                feature_mutant_counts[features[0]] += feature_mutant_counts.pop(feature)
            elif overlapping_genes==4:
                feature_mutant_counts['multiple'] += feature_mutant_counts.pop(feature)
            elif overlapping_genes==5:
                if len(set(features)) == 1:
                    feature_mutant_counts[features[0]] += feature_mutant_counts.pop(feature)
                else:
                    feature_mutant_counts['multiple'] += feature_mutant_counts.pop(feature)
    # splice variants should be counted as 'multiple' if desired - AFTER overlapping genes, since you can have them inside
    for feature in list(feature_mutant_counts):
        if collapse_splice_variants and '|' in feature:
            feature_mutant_counts['multiple'] += feature_mutant_counts.pop(feature)
    # deal with feature boundary cases (AFTER dealing with overlapping genes and splice variants, since they can be inside both)
    for feature in list(feature_mutant_counts):
        if '/' in feature:  
            if feature_boundaries==1:
                features = [x for x in feature.split('/') if x not in 'gene_edge mRNA_edge'.split()]
                random.shuffle(features)
                feature_mutant_counts[features[0]] += feature_mutant_counts.pop(feature)
            else: 
                del feature_mutant_counts[feature]
    # deal with UTR introns
    if ignore_UTR_introns==1:
        feature_mutant_counts["5'UTR"] += feature_mutant_counts.pop("5'UTR_intron")
        feature_mutant_counts["3'UTR"] += feature_mutant_counts.pop("3'UTR_intron")
    elif ignore_UTR_introns==2:
        feature_mutant_counts['intron'] += feature_mutant_counts.pop("5'UTR_intron") + feature_mutant_counts.pop("3'UTR_intron")
    if long_UTR_names:
        feature_mutant_counts["five_prime_UTR"] = feature_mutant_counts.pop("5'UTR")
        feature_mutant_counts["three_prime_UTR"] = feature_mutant_counts.pop("3'UTR")
    return feature_mutant_counts
    # TODO unit-tests!


def insertion_density_by_feature(dataset, feature_lengths, relative_to=None, add_line=False, exclude_multiple=False, 
                                 separate_sense_antisense=True, show_both_for_all=True, feature_length_type=''):
    """ Plot a bar-chart of the average density of insertions in different gene features, for all genes combined. 

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance.

    Feature_lengths should be a feature_name:length dictionary (giving total lengths, mappable lengths, or whatever else), 
     like generated by gff_examine_file.feature_total_lengths.
    Ignores insertions on feature boundaries (usually <1%) - they're not counted as a separate feature, 
     but still included in gene/all insertion counts.
    If exclude_multiple is True, case where feature=MULTIPLE_SPLICE_VARIANTS are treated that way as well; 
     otherwise they're counted as a separate feature.

    If relative_to is 'all', the y labels will be changed so the "all" bar is at 1, and the "all" category itself 
     won't be plotted (since it'll always be at 1); if add_line is True, a horizontal line will be added at 1.
    Similarly if relative_to is 'intergenic'; the y labels will be changed to relative measures putting the intergenic density
     at 1, and neither "intergenic" nor "all" will be plotted.

    If separate_sense_antisense is True, separate mutant counts by sense/antisense orientation compared to the gene,
      and plot them as two separate bar sets on the same graph (except for intergenic categories, which have no orientation);
     if show_both_for_all is True, also plot the both-orientation counts (behind the sense/antisense bars)
      to make it easier to compare the categories that have an orientation to the ones that don't.

    Feature_length_type will be used in the label to specify what kind of length was used (raw, mappable, etc), if given.
    """
    # include or exclude "MULTIPLE_SPLICE_VARIANTS" and "all" features depending on options.
    allowed_relative_to_vals = ('all', 'intergenic', None)
    if relative_to not in allowed_relative_to_vals:  
        raise Exception("relative_to must be one of %s, not %s!"%(allowed_relative_to_vals, relative_to))
    if exclude_multiple:    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR"
    else:                   features_printable = "all intergenic gene exon intron 5'UTR 3'UTR MULTIPLE_SPLICE_VARIANTS"
    if relative_to=='all':          features_printable = ' '.join(features_printable.split()[1:])
    elif relative_to=='intergenic': features_printable = ' '.join(features_printable.split()[2:])
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS').split()
    ### if not separating sense and antisense, just get the feature counts and densities and plot them as a single set of bars
    if not separate_sense_antisense:
        # get feature mutant counts - most of the categories match feature_lengths, fix some details to get the right format: 
        feature_mutant_counts = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset), 
                                                             ignore_UTR_introns=2, long_UTR_names=True)
        ### plot the mutants/length ratios
        feature_mutant_densities = [feature_mutant_counts[f]/feature_lengths[f] for f in features_order]
        mplt.bar(range(len(features_order)), feature_mutant_densities, align='center', edgecolor='none')
    ### if separating sense and antisense, things get more complicated, because most categories should be separated, 
    ###  but some can't (all, intergenic), so there will be 3 sets of bars to plot: sense, antisense, and both
    else:
        # get full feature counts by orientation and for both orientations
        feature_mutant_counts_both =      adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset), 
                                                                       ignore_UTR_introns=2, long_UTR_names=True)
        feature_mutant_counts_sense =     adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset 
                                                                                           if m.orientation=='sense'), 
                                                                       ignore_UTR_introns=2, long_UTR_names=True)
        feature_mutant_counts_antisense = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset 
                                                                                           if m.orientation=='antisense'), 
                                                                       ignore_UTR_introns=2, long_UTR_names=True)
        # which features should be split by orientation and which shouldn't
        no_orientation_features = set("all intergenic".split())
        orientation_features = set(features_order) - no_orientation_features
        # the both-orientation bars should be the same width/position as the two sense/antisense bars side-by-side
        widths_and_offsets = {'both': (0.8, 0), 'sense': (0.4, -0.2), 'antisense': (0.4, 0.2)}
        # plot in reasonable order of orientations, so the legend is in order
        for orientation in "both sense antisense".split():
            # get the mutant counts from the right orientation and the relevant set of features
            if orientation=='both':
                if show_both_for_all:
                    feature_mutant_counts = feature_mutant_counts_both
                else:
                    feature_mutant_counts = general_utilities.filter_dict_by_keys(feature_mutant_counts_both, 
                                                                                  no_orientation_features)
            elif orientation=='sense':
                feature_mutant_counts = general_utilities.filter_dict_by_keys(feature_mutant_counts_sense, orientation_features)
            elif orientation=='antisense':
                feature_mutant_counts = general_utilities.filter_dict_by_keys(feature_mutant_counts_antisense, orientation_features)
            # make a density list based on features_order, with 0 for features where the current orientation doesn't apply
            density_list = []
            for feature in features_order:
                if feature in feature_mutant_counts:    density = feature_mutant_counts[feature]/feature_lengths[feature]
                else:                                   density = 0
                density_list.append(density)
            width, offset = widths_and_offsets[orientation]
            bar_positions = [x+offset for x in range(len(density_list))]
            mplt.bar(bar_positions, density_list, width, color=ORIENTATION_COLORS[orientation], label=orientation, 
                     align='center', edgecolor='none')
        # need legend to show sense/antisense/both
        mplt.legend(loc='upper center', prop=FontProperties(size='medium'))
    ### plot labels etc are the same in both cases
    mplt.xticks(range(len(features_order)), features_printable.split())
    # LATER-TODO the labels are aligned funny, look into fixing that, but it seems tricky...
    mplt.xlim(-0.6, len(features_order)-0.4)
    # save ylim, set it later because some of the stuff below can mess it up
    ymax = mplt.ylim()[1]
    # if result should be relative to all/intergenic, just change the yticks and labels (and make ymax be min 1.5)
    if relative_to:
        if relative_to=='all':
            reference_mapp = feature_mutant_counts['all']/feature_lengths['all']
        elif relative_to=='intergenic':
            reference_mapp = feature_mutant_counts['intergenic']/feature_lengths['intergenic']
        yticks = [0,0.5,1,1.5, 2]
        mplt.yticks([x*reference_mapp for x in yticks], yticks)
        if add_line:
            mplt.hlines(reference_mapp, *mplt.xlim(), linestyles='dotted')
        ymax = max(ymax, 1.5*reference_mapp)
    # otherwise the y scale is the number of mutants per base, so multiply by 1000 so it's per kb
    else:
        mplt.yticks(mplt.yticks()[0], [x*1000 for x in mplt.yticks()[0]])
    # now set ylim
    mplt.ylim(0, ymax)
    # set ylabel
    feature_length_info = ' '+feature_length_type if feature_length_type else ''
    if relative_to:
        mplt.ylabel('mutant density as fraction of %s mutant density\n(over%s feature length)'%(
                        'overall' if relative_to=='all' else 'intergenic', feature_length_info), horizontalalignment='center')
    else:
        mplt.ylabel('average mutants per kb\nof%s feature length'%feature_length_info)

def _get_nice_feature_counts_data(dataset):
    """ Given a dataset, return nicely adjusted feature counts.
    """
    # TODO this currently isn't used anywhere!  Could use it in the three things below to avoid code duplication.
    # TODO add more options!  For what to do with border/multi-gene/etc.
    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR"
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS').split()
    # calculate probability of a mutant falling into each feature vs full genome, based on bin effective lengths
    feature_probabilities = {feature: feature_lengths[feature]/feature_lengths['all'] for feature in feature_lengths}
    # get full feature counts by orientation and for both orientations
    feature_mutant_counts = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset), 
                                                         ignore_UTR_introns=2, long_UTR_names=True)
    return features_order, feature_mutant_counts

def insertion_density_by_feature_pvalues(dataset, feature_lengths, feature_length_type=''):
    """ Check if density of insertions in each feature is different from overall, and between sense and antisense.

    Uses the same statistical methods as insertions_over_gene_length_pvalues - see docstring for that function.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance.

    Feature_lengths should be a feature_name:length dictionary (giving total lengths, mappable lengths, or whatever else), 
     like generated by gff_examine_file.feature_total_lengths.
    Feature_length_type will be used in the label to specify what kind of length was used, if given.

    Ignores insertions on feature boundaries (usually <1%) - they're not counted as a separate feature, 
     but still included in gene/all insertion counts.

    """
    # MAYBE-TODO refactor this with insertions_over_gene_length_pvalues to avoid code duplication?
    # LATER-TODO features_order and features_printable should probably be top-level constants somewhere... in gff_examine_file?
    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR"
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS').split()
    # calculate probability of a mutant falling into each feature vs full genome, based on bin effective lengths
    feature_probabilities = {feature: feature_lengths[feature]/feature_lengths['all'] for feature in feature_lengths}
    # get full feature counts by orientation and for both orientations
    feature_mutant_counts_both =      adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset), 
                                                                   ignore_UTR_introns=2, long_UTR_names=True)
    feature_mutant_counts_sense =     adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset 
                                                                                       if m.orientation=='sense'),
                                                                   ignore_UTR_introns=2, long_UTR_names=True)
    feature_mutant_counts_antisense = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset 
                                                                                       if m.orientation=='antisense'),
                                                                   ignore_UTR_introns=2, long_UTR_names=True)
    # for sense/antisense counts, get rid of features orientation doesn't apply to
    orientation_features = set(features_order) - set("all intergenic".split())
    feature_mutant_counts_sense = general_utilities.filter_dict_by_keys(feature_mutant_counts_sense, orientation_features)
    feature_mutant_counts_antisense = general_utilities.filter_dict_by_keys(feature_mutant_counts_antisense, orientation_features)
    # get the total number of mutants over all bins
    total_N_mutants_both = feature_mutant_counts_both['all']
    total_N_mutants_sense = feature_mutant_counts_sense['gene']
    total_N_mutants_antisense = feature_mutant_counts_antisense['gene']
    ### get the p-values for whether total feature mutant counts are different from the overall count
    total_count_pvalues = [scipy.stats.binom_test(x=feature_mutant_counts_both[feature], n=total_N_mutants_both, 
                                                  p=feature_probabilities[feature]) for feature in features_order if feature!='all']
    total_count_pvalues_adjusted = statistics_utilities.FDR_adjust_pvalues(total_count_pvalues, method='BH')
    # print the data
    print "Whether total feature mutant counts are significantly different from overall, based on %s length:"%feature_length_type
    print "  (overall mutant density %.3g - mutant count %s, genome length %s)"%(total_N_mutants_both/feature_lengths['all'], 
                                                                                 total_N_mutants_both, feature_lengths['all'])
    for feature, pval, adj_pval in zip(features_order[1:], total_count_pvalues, total_count_pvalues_adjusted):
        print " - %s:  %.3g FDR-adjusted p-value (%.3g raw p-value) - mutant density %.3g (mutant count %s, length %s)"%(
                                    feature, adj_pval, pval, feature_mutant_counts_both[feature]/feature_lengths[feature], 
                                    feature_mutant_counts_both[feature], feature_lengths[feature])
    ### get the p-values for whether sense/antisense bin mutant count proportions are different from the overall proportions, 
    features_order_orientation = [feature for feature in features_order if feature in orientation_features]
    orientation_diff_pvalues = [scipy.stats.fisher_exact([[feature_mutant_counts_sense[feature], 
                                                           feature_mutant_counts_antisense[feature]], 
                                                          [total_N_mutants_sense, total_N_mutants_antisense]])[1] 
                                for feature in features_order_orientation]
    orientation_diff_pvalues_adjusted = statistics_utilities.FDR_adjust_pvalues(orientation_diff_pvalues, method='BH')
    # print the data
    print "Whether feature mutant sense/antisense proportions are significantly different from overall:"
    print "  (overall mutant sense:antisense counts are %s:%s (over total gene length %s)"%(total_N_mutants_sense, 
                                                                              total_N_mutants_antisense, feature_lengths['gene']) 
    for feature, pval, adj_pval in zip(features_order_orientation, orientation_diff_pvalues, orientation_diff_pvalues_adjusted):
        print " - %s:  %.3g FDR-adjusted p-value (%.3g raw p-value) - %s:%s sense:antisense mutants (feature length %s)"%(
                                        feature, adj_pval, pval, feature_mutant_counts_sense[feature], 
                                        feature_mutant_counts_antisense[feature], feature_lengths[feature])

def insertion_density_by_feature_pvalues_all_pairs(dataset, feature_lengths, feature_length_type=''):
    """ Check if density of insertions is different between each pair of features.

    Uses the chi-squared test of independence, and FDR-adjustment.

    See insertion_density_by_feature_pvalues docstring for argument descriptions.
    """
    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR"
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS').split()
    # get full feature counts by orientation and for both orientations
    feature_mutant_counts = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset), 
                                                         ignore_UTR_introns=2, long_UTR_names=True)
    features = sorted(feature_mutant_counts.keys())
    all_feature_pairs = list(itertools.combinations(sorted(set(feature_mutant_counts.keys()) - set(['all'])), 2))
    all_pvalues_1 = []
    # for each feature, look at the number of mappable positions that did and did not have an insertion,
    #  and compare those two numbers between features to see if the proportions are the same.
    for feature1,feature2 in all_feature_pairs:
        feature1_ins_count = feature_mutant_counts[feature1]
        feature2_ins_count = feature_mutant_counts[feature2]
        feature1_no_ins_count = feature_lengths[feature1] - feature1_ins_count
        feature2_no_ins_count = feature_lengths[feature2] - feature2_ins_count
        pvalue_1 = statistics_utilities.chisquare_independence([feature1_ins_count,feature1_no_ins_count], 
                                                               [feature2_ins_count,feature2_no_ins_count])
        all_pvalues_1.append(pvalue_1)
    adj_pvalues_1 = statistics_utilities.FDR_adjust_pvalues(all_pvalues_1)
    print "P-values for difference between each pair of features in this dataset, based on %s length:"%feature_length_type
    for (feature1,feature2),p1,ap1 in zip(all_feature_pairs, all_pvalues_1, adj_pvalues_1):
        mark = '**' if p1<0.0001 else ('*' if p1<0.05 else '-')
        print " %s pvalue %.3g (FDR-adj %.3g) for difference between %s and %s"%(mark, p1, ap1, feature1, feature2)

def insertion_density_by_feature_pvalues_compare_datasets(dataset1, dataset2, subtract_2_from_1=False):
    """ Check if density of insertions is different between the two datasets, for each feature.

    Uses the chi-squared test of independence, and FDR-adjustment.
    If subtract_2_from_1 is True, subtract dataset2 from dataset1 
      (use in case dataset2 is a subset of 1, to make the data independent).

    See insertion_density_by_feature_pvalues docstring for argument descriptions.
    """
    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR"
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS').split()
    # get full feature counts
    feature_mutant_counts_1 = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset1), 
                                                           ignore_UTR_introns=2, long_UTR_names=True)
    feature_mutant_counts_2 = adjust_feature_mutant_counts(collections.Counter(m.gene_feature for m in dataset2), 
                                                           ignore_UTR_introns=2, long_UTR_names=True)
    assert set(feature_mutant_counts_1.keys()) == set(feature_mutant_counts_2.keys())
    features = sorted(feature_mutant_counts_1.keys())
    # subtract dataset2 from dataset1 if needed
    if subtract_2_from_1:
        feature_mutant_counts_1 = {f: count1 - feature_mutant_counts_2[f] for (f,count2) in feature_mutant_counts_1.items()}
    # To compare between two datasets, compare insertions in featureX vs in other features for both datasets, 
    #   rather than comparing insertions vs non-insertions in featureX, to stay independent of different total insertion numbers.
    all_pvalues_1 = []
    for feature in features:
        feature_ins_count_1 = feature_mutant_counts_1[feature]
        feature_ins_count_2 = feature_mutant_counts_2[feature]
        other_ins_count_1 = feature_mutant_counts_1['all'] - feature_ins_count_1
        other_ins_count_2 = feature_mutant_counts_2['all'] - feature_ins_count_2
        pvalue_1 = statistics_utilities.chisquare_independence([feature_ins_count_1,other_ins_count_1], 
                                                               [feature_ins_count_2,other_ins_count_2])
        all_pvalues_1.append(pvalue_1)
    adj_pvalues_1 = statistics_utilities.FDR_adjust_pvalues(all_pvalues_1)
    print "P-values for difference between mutant densities in each feature between the two datasets, based on %s length:"%\
            feature_length_type
    for feature,p1,ap1 in zip(features, all_pvalues_1, adj_pvalues_1):
        mark = '**' if p1<0.0001 else ('*' if p1<0.05 else '-')
        print " %s %s pvalue %.3g (FDR-adj %.3g)"%(mark, feature, p1, ap1)

### mappability by gene/feature

def mappability_feature_barchart(feature_total_lengths, feature_mappable_lengths):
    """ Plot a bar-chart of the average mappability of insertions in different gene features, for all genes combined. 

    Both inputs should be feature_type:length dictionaries, like generated by gff_examine_file.feature_total_lengths 
     and mutant_simulations.gene_mappability (3rd output) respectively.
    """
    features_printable = "all intergenic gene exon intron 5'UTR 3'UTR multi"
    features_order = features_printable.replace("5'", "five_prime_").replace("3'", "three_prime_").replace('exon','CDS')\
                                        .replace('multi', 'MULTIPLE_SPLICE_VARIANTS').split()
    assert set(feature_total_lengths.keys()) == set(feature_mappable_lengths.keys()) == set(features_order)
    feature_mappability_fractions = [feature_mappable_lengths[f]/feature_total_lengths[f] for f in features_order]
    mplt.bar(range(len(features_order)), feature_mappability_fractions, align='center', edgecolor='none')
    mplt.title('fraction of mappable insertion positions in total feature length')
    mplt.xticks(range(len(features_order)), features_printable.split())
    mplt.xlim(-0.6, len(features_order)-0.4)
    mplt.ylim(0,1)
    mplt.ylabel('mappable fraction')


def mappability_gene_scatterplot(gene_total_lengths, gene_mappable_lengths, logscale=True):
    """ Make a scatterplot of gene total vs mappable lengths.

    Both inputs should be feature_type:length dictionaries, like generated by gff_examine_file.gene_lengths 
     and mutant_simulations.gene_mappability (1st output) respectively.
    """
    genes_in_order = list(gene_total_lengths.keys())
    mplt.plot([gene_total_lengths[g] for g in genes_in_order], [gene_mappable_lengths[g] for g in genes_in_order], 
              'o', linewidth=0, markeredgecolor='blue', markerfacecolor='none')
    # make a line for 100% mappable (1:1)
    max_total, min_total = max(gene_total_lengths.values()), min(gene_total_lengths.values())
    max_mappable, min_mappable = max(gene_mappable_lengths.values()), min(gene_mappable_lengths.values())
    mplt.plot([0.001,max_total*10], [0.001,max_total*10], '--', color='black')
    mplt.title('fraction of mappable insertion positions in gene lengths')
    mplt.xlabel('total gene length' + (' (log)' if logscale else ''))
    mplt.ylabel('mappable gene length' + (' (log)' if logscale else ''))
    if logscale:
        mplt.xscale('symlog')
        mplt.yscale('symlog')
        mplt.xlim(min_total/2, 2*max_total)
        mplt.ylim(min_mappable/2 if min_mappable>0 else -0.7, 2*max_mappable)
    else:
        range_total, range_mappable = max_total-min_total, max_mappable-min_mappable
        mplt.xlim(min_total - range_total/10, max_total + range_total/10)
        mplt.ylim(min_mappable - range_mappable/10, max_mappable + range_mappable/10)


def mappability_gene_histogram(gene_mappable_lengths, gene_total_lengths=None, cumulative=False, log_xscale=True, N_bins=100):
    """ Make a histogram of gene mappable lengths or mappable fractions, cumulative or not..

    If gene_total_lengths is given, plot gene mappable/total length ratios; otherwise plot gene mappable lengths.
    If logscale, take the log of the values (lengths or ratios) before making the histogram. 

    If cumulative, the histogram will be cumulative, normalized to 100%, 
     and it won't actually be binned, there's not much point in doing that.

    First two args should be feature_type:length dictionaries, like generated by mutant_simulations.gene_mappability (1st output) 
     and gff_examine_file.gene_lengths respectively.
    """
    if gene_total_lengths is None:  values = gene_mappable_lengths.values()
    else:                           values = [gene_mappable_lengths[g]/gene_total_lengths[g] for g in gene_total_lengths]
    if log_xscale:              values = [-1 if x==0 else numpy.log10(x) for x in values]
    # when doing a cumulative histogram, don't actually use bins, there's really no point.
    if cumulative:  mplt.hist(values, bins=len(gene_mappable_lengths), histtype='step', cumulative=True, normed=True)
    else:           mplt.hist(values, bins=N_bins, histtype='bar', cumulative=False, normed=False)
    val_type = 'lengths' if gene_total_lengths is None else 'fractions'
    mplt.title('%shistogram of gene mappable %s'%('cumulative ' if cumulative else '', val_type))
    if cumulative:
        mplt.xlabel('gene mappable %s'%(val_type + (' (log)' if log_xscale else '')))
        mplt.ylabel('fraction of genes with that value or lower')
        mplt.ylim(0,1.01)
        xmax = mplt.xlim()[1]
        if log_xscale:  mplt.xlim(-1.1, xmax)
        else:           mplt.xlim(-0.01*xmax, xmax)
    else:
        mplt.xlabel('gene mappable %s, binned into %s bins'%(val_type + (' (log)' if log_xscale else ''), N_bins))
        mplt.ylabel('number of genes in bin')
    # if log_xscale, change xtick labels using the reverse of how the values were transform;
    #  if using fractions, we know we want the 0-1 range, so just set the ticks to that
    if log_xscale:              
        if gene_total_lengths is None:
            xticks = mplt.xticks()[0]
            mplt.xticks(xticks, [0 if x==-1 else general_utilities.int_or_float(10**x) for x in xticks])
        else:
            print "Warning: xticks will be messy if using log_xscale and providing gene_total_lengths - you probably shouldn't."


############################################## Plotting readcounts ##############################################

######### Help functions

def _extract_readcounts_from_dataset(dataset, perfect_only=False, total_and_perfect=False):
    """ Extract a readcount list from a mutant dataset. 

    Normally the returned list will contain total readcounts for each mutant;
     if perfect_only, the perfect instead of total readcounts are used; 
     if total_and_perfect, readcount_list will contain (total, perfect) tuples insteaad of single values.
    """
    if total_and_perfect:   return [(m.total_read_count, m.perfect_read_count) for m in dataset]
    elif perfect_only:      return [m.perfect_read_count for m in dataset]
    else:                   return [m.total_read_count for m in dataset]


def _get_readcount_data_if_needed(dataset, perfect_only=False, total_and_perfect=False):
    """ Return total/perfect readcount list from dataset object; if dataset is a number list already, just return it.

    For extraction details see _extract_readcounts_from_dataset.
    Of course, if the dataset is just a list of numbers, there's no way of checking whether those numbers are
     total or perfect-only, so the caller has to deal with that. 
    """
    # if dataset is a dataset object, extract the desired readcount list from it
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        return _extract_readcounts_from_dataset(dataset, perfect_only, total_and_perfect)
    # otherwise, check if dataset is a list of numbers (or number pairs if total_and_perfect): 
    #  if yes, return it, otherwise give an error.
    else:
        try:
            if total_and_perfect:
                assert all([len(x)==2 for x in dataset])
                sum(dataset[1]+dataset[2]+dataset[-1])
            else:
                sum(dataset[:10]+dataset[-10:])
            return dataset
        except TypeError, AssertionError:
            raise Exception("dataset format not recognized, in _get_readcount_data_if_needed!")


def get_all_datasets_glob(glob_pattern, split_filenames_on=None, readcounts_only=False, perfect_only=False, total_and_perfect=False):
    """ Get a name:dataset or name:readcount_list dict for all the files matching glob_pattern. 

    Files read using mutant_analysis_classes.read_mutant_file.

    By default returen a name:dataset dict, with full dataset objects. 
    If readcounts_only, return a name:readcount_list dict instead, to use less memory - 
        For extraction details see _extract_readcounts_from_dataset.

    By default dataset_name is file basename with no extension; if split_filenames_on is not None, 
     split the file basename on the value and use the first component as dataset_name.
    """
    all_datasets = {}
    for infile in glob.glob(glob_pattern):
        filename = splitext(os.path.basename(infile))[0]
        if split_filenames_on is not None:
            filename = filename.split(split_filenames_on)[0]
        dataset = mutant_analysis_classes.read_mutant_file(infile)
        if readcounts_only: all_datasets[filename] = _extract_readcounts_from_dataset(dataset, perfect_only, total_and_perfect)
        else:               all_datasets[filename] = dataset
    return all_datasets
    # MAYBE-TODO this function may be applicable to more than readcounts - if so, move it up?  Should probably remove the readcounts_only/etc options in that case.  So really that would be a different function...


######### Single-plot functions

def readcounts_sorted_plot(dataset_name_list, color_dict=None, perfect_only=False, x_max=None, y_max=None, y_min=0.7, 
                           log_y=True, legend=True, legend_kwargs=None):
    """ Basic sorted readcount plots - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        readcounts = _get_readcount_data_if_needed(dataset, perfect_only=perfect_only)
        kwargs = {} if color_dict is None else {'c': color_dict[name]}
        _ = mplt.plot(sorted(readcounts), '.', linewidth=0, label=name, markersize=2, **kwargs)
    mplt.xlabel('mutants sorted by readcount\n(independently for each dataset)')
    mplt.ylabel('readcount (%s)'%('log' if log_y else 'linear'))
    if log_y:   mplt.yscale('log')
    else:       mplt.yscale('linear')
    plotting_utilities.set_axes_limits(None, x_max, y_min, y_max)
    plotting_utilities.remove_half_frame()
    if legend:
        if legend_kwargs is None:   mplt.legend(loc=4, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                                                prop=FontProperties(size='medium'))
        else:                       mplt.legend(**legend_kwargs)
      

def readcounts_cumulative_plot(dataset_name_list, color_dict=None, linewidth_dict=None, linestyle_dict=None, perfect_only=False, 
                               x_min=None, x_max=None, y_max=None, y_min=None, log_x=True, legend=True, legend_kwargs=None):
    """ ROC-like "cumulative histogram" - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        N_mutants = len(dataset)
        # adding a "0th mutant at readcount 1" data point to make sure all lines start at readcount 1, 
        #  even if the first real mutant has 100 reads - looks better that way, and more consistent with the nature of the plot.
        readcounts = [1] + sorted(_get_readcount_data_if_needed(dataset, perfect_only=perfect_only))
        mutant_percentages = [0] + [n/N_mutants*100 for n in range(N_mutants)]
        kwargs = {} 
        if color_dict is not None:      kwargs.update({'c': color_dict[name]})
        if linewidth_dict is not None:  kwargs.update({'linewidth': linewidth_dict[name]})
        if linestyle_dict is not None:  kwargs.update({'linestyle': linestyle_dict[name]})
        _ = mplt.plot(readcounts, mutant_percentages, label=name, **kwargs)
    mplt.xlabel('readcount (%s)'%('log' if log_x else 'linear'))
    if log_x:   mplt.xscale('log')
    else:       mplt.xscale('linear')
    mplt.ylabel('% of mutants with that readcount or less')
    plotting_utilities.set_axes_limits(x_min, x_max, y_min, y_max)
    plotting_utilities.remove_half_frame()
    if legend:
        if legend_kwargs is None:   mplt.legend(loc=4, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                                                prop=FontProperties(size='medium'))
        else:                       mplt.legend(**legend_kwargs)
  

def readcounts_histogram(dataset_name_list, color_dict=None, Nbins=100, histtype='bar', perfect_only=False, log_x=True, log_y=False, 
                         readcount_max=None, y_max=None, y_min=None, no_edges=False, legend=True, legend_kwargs=None):
    """ Normal histogram (any type), linear or logscale - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        readcounts = _get_readcount_data_if_needed(dataset, perfect_only=perfect_only)
        if readcount_max:   readcounts = [c for c in readcounts if c<=readcount_max]
        if log_x:           readcounts = [math.log10(c) for c in readcounts]
        if log_x:           bin_range = (0, max(readcounts))
        else:               bin_range = (1, max(readcounts))
        kwargs = {} 
        if color_dict is not None:  kwargs.update({'facecolor': color_dict[name]})
        if no_edges:                kwargs.update({'edgecolor': 'none'})
        # using range to make sure that if we have one dataset with some 1-read mutants and one with only 10+-read mutants, 
        #  they both get the same bins, rather than N bins stretched over the (1,max) and (10,max) ranges respectively.
        _ = mplt.hist(readcounts, histtype=histtype, linewidth=0.3, bins=Nbins, range=bin_range, label=name, log=log_y, 
                      align='mid', **kwargs)
    mplt.xlabel('readcount (%s) (range binned into %s bins)'%(('log' if log_x else 'linear'), Nbins))
    mplt.ylabel('number of mutants with that readcount (%s)'%('log' if log_y else 'linear'))
    # TODO if x_log, change the xticklabels (and maybe add ticks?) to look like log!
    plotting_utilities.set_axes_limits(None, readcount_max, y_min, y_max)
    plotting_utilities.remove_half_frame()
    if legend:
        mplt.legend(loc=1, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                    prop=FontProperties(size='medium'))
  

######### Multi-plot functions

# this is OLD, probably doesn't work - MAYBE-TODO make it work if I actually want to?
def sample_row_multi_plots(sample_name, sample_N, total_samples, first_cumulative=True, 
                                xmax_ymax_Nbins_logx_logy_list=[(None,None,100,True,False),(None,None,100,False,False)], 
                                histtype='bar', if_xlabels=True):
    """ Plot one sample readcount info with multiple methods: cumulative histogram, and any number of normal histograms with different settings. """
    N_plots = len(xmax_ymax_Nbins_logx_logy_list) + int(first_cumulative)
    curr_plot = sample_N * N_plots + 1
    if first_cumulative:
        mplt.subplot(total_samples, N_plots, curr_plot)
        plot_data_cumulative([sample_name], legend=False)
        if not if_xlabels:
            mplt.xlabel('')
            plotting_utilities.remove_xticklabels()
        mplt.ylabel(sample_name)
        curr_plot += 1
    for xmax,ymax,Nbins,log_x,log_y in xmax_ymax_Nbins_logx_logy_list:
        mplt.subplot(total_samples, N_plots, curr_plot)
        plot_data_hist([sample_name], Nbins=Nbins, logx=log_x, logy=log_y, maxcount=xmax, histtype=histtype, legend=False)
        if ymax is not None:
            mplt.ylim((0,ymax))
        if not if_xlabels:
            mplt.xlabel('')
            plotting_utilities.remove_xticklabels()
        if curr_plot != sample_N * N_plots + 1:
            mplt.ylabel('')
            plotting_utilities.remove_yticklabels()
        else:
            mplt.ylabel(sample_name)
        curr_plot += 1



######################################## Plotting adjacent mutant data ###########################################

def _get_singles_from_counter(val_count_dict, max_val):
    """ Given a val:count dict, return a val list with each val present count times. (for creating histograms). """
    val_singles_list = []
    for val,count in val_count_dict.items():
        if max_val is None or val <= max_val:
            val_singles_list.extend([val]*count)
    return val_singles_list


colors_by_adjacent_category = {'adjacent-same-strand': 'red', 'adjacent-opposite-both': 'cyan', 'same-pos-opposite': 'orange', 
                               'adjacent-opposite-toward': 'green', 'adjacent-opposite-away': 'blue' }


def adjacent_distance_histogram(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
                                     incl_same_pos_opposite=False, max_distance=None, N_bins=100, symlog_y=False, symthresh=100):
    """ Step-histogram of the number of adjacent mutant pairs by distance, separated by type. 

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
     (Be careful with incl_same_pos_opposite - those are always distance 0, so it's hard to show them meaningfully in comparison
       to the ones with a wider ditance range unless each distance has its own bin, like with max_distance 100 and N_bins 101). 
    Only include mutant pairs with distance at most max_distance, if given.
    If symlog_y, the y axis will be symlog-scale (linear up to symthresh, exponential afterward). 
    """
    datasets, labels, colors = [], [], []
    if incl_same_strand:    
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_same_strand_dict, max_distance))
        labels.append('adjacent-same-strand')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_separate:
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_opposite_toward_dict, max_distance))
        labels.append('adjacent-opposite-toward')
        colors.append(colors_by_adjacent_category[labels[-1]])
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_opposite_away_dict, max_distance))
        labels.append('adjacent-opposite-away')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_both:
        adjacent_opposite_both = general_utilities.add_dicts_of_ints(dataset.summary.adjacent_opposite_toward_dict, 
                                                                     dataset.summary.adjacent_opposite_away_dict)
        datasets.append(_get_singles_from_counter(adjacent_opposite_both, max_distance))
        labels.append('adjacent-opposite-both')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_same_pos_opposite:
        datasets.append([0]*dataset.summary.same_position_opposite)
        labels.append('same-pos-opposite')
        colors.append(colors_by_adjacent_category[labels[-1]])
        if max_distance>=N_bins:
            print("Warning: if each distance doesn't get its own bin (N_bins = max_distance+1), "
                  +"the same-pos-opposite count won't be fairly represented!")
    xmin = min(sum(datasets, []))
    xmax = max(sum(datasets, []))
    xspread = xmax - xmin
    # have to do bins by hand - TODO or do I?
    #bins = [(low + i/N_bins*spread) for i in range(N_bins+1)]
    hist_data, bin_edges, hist_patches = mplt.hist(datasets, label=labels, color=colors, bins=N_bins, histtype='step')
    plotting_utilities.remove_half_frame()
    mplt.xlabel('distance between the two mutants in a pair (bp)' + (' (binned into %s bins)'%N_bins if (xspread+1)>N_bins else ''))
    mplt.ylabel('number of mutant pairs with given distance')
    # sometimes mplt.hist gets the axis limits wrong, so fix them by hand
    ymax = max(sum([list(d) for d in hist_data], []))
    mplt.xlim(xmin-xspread/N_bins, xmax)
    mplt.ylim(0-ymax/100, ymax*1.05)
    mplt.legend(title='pair categories by relative orientation:', prop=FontProperties(size='medium'))


def adjacent_dist1_barchart(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
                              incl_same_pos_opposite=True, logscale_x=False):
    """ Horizontal bar-plot of the number of 0-1bp adjacent mutant pairs, separated by type. 

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
    If logscale_x, the x axis will be log-scale. 
    """
    counts, labels, colors = [], [], []
    if incl_same_strand:    
        counts.append(dataset.summary.adjacent_same_strand_dict[1])
        labels.append('adjacent-same-strand')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_both:
        counts.append(dataset.summary.adjacent_opposite_toward_dict[1] + dataset.summary.adjacent_opposite_away_dict[1])
        labels.append('adjacent-opposite-both')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_separate:
        counts.append(dataset.summary.adjacent_opposite_toward_dict[1])
        labels.append('adjacent-opposite-toward')
        colors.append(colors_by_adjacent_category[labels[-1]])
        counts.append(dataset.summary.adjacent_opposite_away_dict[1])
        labels.append('adjacent-opposite-away')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_same_pos_opposite:
        counts.append(dataset.summary.same_position_opposite)
        labels.append('same-pos-opposite')
        colors.append(colors_by_adjacent_category[labels[-1]])
    # use everything reversed, so that the first one is on top instead of bottom (convert to list, matplotlib chokes on iterators)
    mplt.barh(range(len(counts)), list(reversed(counts)), color=list(reversed(colors)), align='center', log=logscale_x)
    mplt.yticks(range(len(labels)), list(reversed(labels)))
    mplt.xlabel('number of mutant pairs in given category')
    plotting_utilities.remove_half_frame()


def _get_sorted_ratios_from_dict(ratio_distance_dict, min_distance, max_distance):
    """ Given a distance:readcount_pair_list dict, return a list of readcount_pair ratios for distances between min and max. 
    """
    ratio_list = []
    for dist in ratio_distance_dict.keys():
        if min_distance <= dist <= max_distance:
            for x,y in ratio_distance_dict[dist]:
                x,y = sorted([x,y])
                ratio_list.append(y/x)
    return sorted(ratio_list)


def adjacent_readcount_ratio_plot(dataset, distance_cutoffs, distance_linestyles=None, incl_same_strand=True, 
                           incl_opposite_both=True, incl_opposite_separate=True, incl_same_pos_opposite=True, logscale_y=True):
    """ Plot sorted readcount ratios of adjacent mutant pairs, separated by type and optionally distance categories.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    Distance_cutoffs should be a list of numbers: [1,10,1000] will result in plotting separate lines for mutant pairs with 
     distance 1, 2-10, and 11-1000 for each type (except same-pos-opposite, since those always have a distance of 0)
    Distance_linestyles should be a list of valid matplotlib.pyplot.plot linestyle values, one for each cutoff in distance_cutoffs.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
    If logscale_y, the y axis will be log-scale. 
    """
    # LATER-TODO do I actually want sorted readcount ratios, or would a histogram or something be better?  But histograms don't have linestyles... I could just make histogram-like data and use mplt.plot to plot it, though.
    if distance_linestyles is None:
        if len(distance_cutoffs)==2:    distance_linestyles = ['-', ':']
        elif len(distance_cutoffs)==3:  distance_linestyles = ['-', '--', ':']
        elif len(distance_cutoffs)==4:  distance_linestyles = ['-', '--', '-.', ':']
        else:   raise Exception("If there aren't 2-4 distance_cutoffs, need to provide distance_linestyles!")
    ratio_lists, labels, colors, linestyles = [], [], [], []
    distance_mins = [0] + [x+1 for x in distance_cutoffs[:-1]]
    if incl_same_strand:    
        category = 'adjacent-same-strand'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_same_strand_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_opposite_both:
        category = 'adjacent-opposite-both'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(sorted(sum([_get_sorted_ratios_from_dict(D, distance_min, distance_max) 
                                           for D in (dataset.summary.adjacent_opposite_toward_readcounts_dict, 
                                                     dataset.summary.adjacent_opposite_away_readcounts_dict,)], [])))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_opposite_separate:
        category = 'adjacent-opposite-toward'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_opposite_toward_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
        category = 'adjacent-opposite-away'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_opposite_away_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_same_pos_opposite:
        category = 'same-pos-opposite'
        ratio_lists.append(_get_sorted_ratios_from_dict({0: dataset.summary.same_position_opposite_readcounts}, 0, 0))
        labels.append('%s (always dist 0)'%(category))
        colors.append(colors_by_adjacent_category[category])
        linestyles.append(distance_linestyles[0])
    for ratio_list,label,color,linestyle in zip(ratio_lists, labels, colors, linestyles):
       mplt.plot([(y+1)/len(ratio_list) for y in range(len(ratio_list))], ratio_list, 
                 color=color, linestyle=linestyle, label='%s - %s pairs'%(label, len(ratio_list)))
    if logscale_y:   mplt.yscale('log')
    plotting_utilities.remove_half_frame()
    mplt.ylabel('readcount ratio between the two mutants in a pair')
    mplt.xlabel('all adjacent mutant pairs, sorted by readcount ratio, normalized to 100%')
    smallfont = FontProperties(size='smaller')
    mplt.legend(prop=smallfont, loc=2)
    mplt.title('Ratios between the readcounts of adjacent mutant pairs,\nby category and distance')


################################################# Plotting RISCC data ##################################################

# TODO check that those don't have any missing requirements - I moved them from mutant_Carette.py.

def distance_histogram(distance_datasets, labels=None, colors=None, linestyles=None, N_bins=50, histtype='step', normalized=True, 
                       outfile=None, filetypes='svg png', x_max=None, sample_info=''):
    """ Plot histogram of the distances in the given datasets; optionally save to file. 
    """
    mplt.figure()
    max_y = 0
    for N, distances in enumerate(distance_datasets):
        if not len(distances):
            print "Dataset %s contains no data! Skipping."%N
            continue
        plot_kwargs = {}
        if labels:      plot_kwargs['label'] = "%s (%s)"%(labels[N], len(distances))
        if colors:      plot_kwargs['color'] = colors[N]
        if linestyles:  plot_kwargs['linestyle'] = linestyles[N]
        y_values, _, _ = mplt.hist(distances, bins=N_bins, histtype=histtype, normed=str(normalized), **plot_kwargs)
        max_y = max(max_y, max(y_values))
    mplt.xlabel('highest distance between cassette-side and genome-side reads')
    mplt.ylim(-max_y/20, max_y*1.1)
    if x_max is not None:
        mplt.xlim(mplt.xlim()[0], x_max)
    if normalized:  mplt.ylabel('fraction of mutants')
    else:           mplt.ylabel('number of mutants')
    if labels:
        mplt.legend()
    mplt.title("Carette confirmed distance distribution" + "\n(%s)"%sample_info if sample_info else "")
    if outfile:
        plotting_utilities.savefig(outfile, filetypes)
    mplt.close()

def plot_confirmed_dist_vs_percent(dataset, dataset_name, min_genomic_reads=1, min_conf_reads=0, mutant_filter=None, 
                                   xmax=None, max_allowed_distance=3000, markersize=4, alpha=0.3, color='black'):
    """ Make a max-confirmed-distance vs %-confirming-reads scatterplot.
    """
    D = max_allowed_distance
    if mutant_filter:   mutants = [m for m in dataset if mutant_filter(m)]
    else:               mutants = dataset
    # grab the data from old or new version of mutant datasets
    try:
        filtered_data = [(m.RISCC_max_confirmed_distance(D), m.RISCC_N_confirming_seqs(D), m.RISCC_N_non_confirming_seqs(D)) 
                         for m in mutants]
    except AttributeError:
        filtered_data = [(m.Carette_max_confirmed_distance(D), m.Carette_N_confirming_reads(D), m.Carette_N_non_confirming_reads(D)) 
                         for m in mutants]
    filtered_data_2 = [x for x in filtered_data if x[1]+x[2] >= min_genomic_reads and x[1] >= min_conf_reads]
    print "Filtering data: %s total mutants, %s passed general filter, %s passed readcount filter."%(len(dataset), 
                                                                                     len(filtered_data), len(filtered_data_2))
    max_conf_len = [x[0] for x in filtered_data_2]
    percent_wrong_reads = [x[2]/(x[1]+x[2])*100 for x in filtered_data_2]
    mplt.plot(max_conf_len, percent_wrong_reads, 
              marker='.', markerfacecolor=color, markersize=markersize, linestyle='None', markeredgecolor='None', alpha=alpha)
    mplt.xlabel('max distance between cassette-side sequence and matching genome-side read')
    mplt.ylabel("% reads that don't match the cassette-side sequence")
    mplt.ylim(-1, 101)
    if xmax is None:    xmax = mplt.xlim()[1]
    mplt.xlim(-10, xmax)
    mplt.title('Data on Carette genome-side reads confirming the cassette-side sequence,\n dataset %s, min %s total and %s confirming reads.'%(dataset_name, min_genomic_reads, min_conf_reads))

################################################# Testing etc ##################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
