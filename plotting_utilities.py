#! /usr/bin/env python2.7
"""
Various plotting utilities I wrote (usually for matplotlib) - see docstring for each function for what it does.
 -- Weronika Patena, 2010-2012
"""

# standard library
import itertools
import unittest
import os
# other packages
import numpy
import matplotlib
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties

# For useful tricks see ~/computers_and_programming/matplotlib_notes_and_tricks.txt file.


def savefig(figname, extensions=None, padding=0.2, dpi=300, dont_check_figname=False, **kwargs):
    """ Save current figure as figname if figname has an extension, or as figname.ext for each extension, with given padding/dpi. 
    
    Extensions can be a list, or a space-separated string that will be split into a list, or None.  If None:
     - if figname looks like it has an extension (basename contains '.'), just use figname
     - otherwise use figname.png
   
    Uses bbox_inches='tight' to get rid of extra padding and make sure no text/elements are outside the figure.
    """
    default_extension = 'png'
    try:                    extensions = extensions.split()
    except AttributeError:  pass
    kwargs.update(dict(bbox_inches='tight', pad_inches=padding, dpi=dpi))
    if extensions is not None:
        for extension in extensions:    mplt.savefig('%s.%s'%(figname, extension), **kwargs)
    elif os.path.splitext(figname)[1]:  mplt.savefig(figname, **kwargs)
    else:                               mplt.savefig('%s.%s'%(figname, default_extension), **kwargs)


################################ EASY FUNCTIONS FOR SPECIFIC PLOT TYPES ###################################

def stacked_bar_plot(list_of_sample_category_lists, sample_names=[], bar_width=0.7, normalize=False, colors='bgrcmy', **kwargs):
    """ Plot list_of_sample_category_lists as a stacked bar plot (all categories per sample on top of each other). 
    Return list of plot_bar objects, to use in legend (like this: "mplt.legend(plot_bar_list, name_list)"). 
    """
    if not len(set([len(category_list) for category_list in list_of_sample_category_lists])) == 1:
        raise ValueError("All lists in list_of_sample_category_lists must be the same length!")
    if sample_names and not len(sample_names)==len(list_of_sample_category_lists):
        raise ValueError("list_of_sample_category_lists and sample_names must be the same length!")
    N_samples = len(list_of_sample_category_lists)
    N_categories = len(list_of_sample_category_lists[0])
    if 'align' not in kwargs:   kwargs['align'] = 'center'
    if not sample_names:
        sample_names = ['' for _ in range(N_samples)]
    positions = range(N_samples)
    category_bars_for_legend = []
    bar_bottoms = [0 for _ in sample_names]
    if normalize:
        for (i, data) in enumerate (list_of_sample_category_lists):
            data_sum = sum(data)
            norm_data = [x/data_sum*100 for x in data]
            list_of_sample_category_lists[i] = norm_data
    for category_N, color in zip(range(N_categories), itertools.cycle(colors)):
        category_values = [sample_category_list[category_N] for sample_category_list in list_of_sample_category_lists]
        plot_bars = mplt.bar(positions, category_values, bottom=bar_bottoms, color=color, width=bar_width, **kwargs)
        bar_bottoms = [x+y for (x,y) in zip(bar_bottoms,category_values)]
        category_bars_for_legend.append(plot_bars[0])
    mplt.xticks(positions, sample_names)
    if normalize:
        mplt.ylim(0,100)
    return category_bars_for_legend



################################ COSMETIC MODIFICATIONS TO EXISTING PLOTS ###################################

# NOTE: mplt.gca() is get_current_axis, for when I was to act on it directly with ax.* methods and I don't have an ax arg given;
#       mplt.sca(ax) is set_current_axis, for when I have an ax input and I want to act on it indirectly with mplt.* methods.
# The mplt.* methods are interactive, but the ax.* methods etc need an mplt.draw() to show the results.

# MAYBE-TODO could make the "if ax is None:  ax = mplt.gca(); else:  mplt.sca(ax)" line into a decorator, since it shows up everywhere

def legend(*args, **kwargs):
    """ make the legend with medium instead of large font size. """
    mplt.legend(*args, prop=FontProperties(size='medium'), **kwargs)

def remove_legend(ax=None):
    """ Remove legend for ax or the current axes (detected with gca()). """
    # from Scipy matplotlib cookbook - http://www.scipy.org/Cookbook/Matplotlib/Legend
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.legend_ = None
    # alternative version here - http://stackoverflow.com/questions/5735208/remove-the-legend-on-a-matplotlib-figure
    # ax.legend().set_visible(False)
    # yes, the draw is actually needed in this case!
    mplt.draw()


def remove_xticklabels(ax=None):
    """ Remove x tick labels (leaving the ticks unchanged); acts on ax, or current axes if ax is None. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    xlim = mplt.xlim()
    xticks = mplt.xticks()[0]
    mplt.xticks(xticks, [''] * len(xticks))
    mplt.xlim(xlim)

def remove_yticklabels(ax=None):
    """ Remove y tick labels (leaving the ticks unchanged); acts on ax, or current axes if ax is None. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    lim = mplt.ylim()
    yticks = mplt.yticks()[0]
    mplt.yticks(yticks, [''] * len(yticks))
    mplt.ylim(lim)


def set_axes_limits(x_min=None, x_max=None, y_min=None, y_max=None, ax=None):
    """ Set whichever limits aren't None, keeping the others the same. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    # MAYBE-TODO could add an option to only increase (take max of current and new value), or only decrease, etc...
    if ax is None:  ax = mplt.gca()
    if x_min is None:   x_min = mplt.xlim()[0]
    if x_max is None:   x_max = mplt.xlim()[1]
    mplt.xlim((x_min, x_max))
    if y_min is None:   y_min = mplt.ylim()[0]
    if y_max is None:   y_max = mplt.ylim()[1]
    mplt.ylim((y_min, y_max))


def color_plot_frame(ax=None, color='grey', color_frame=True, color_ticks=True, color_ticklabels=True): 
    """ Change the color of the frame/ticks/ticklabels of ax (a matplotlib.axes.AxesSubplot object) to color. """
    # source: http://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    if color_frame:         mplt.setp(ax.spines.values(), color=color)
    if color_ticks:         mplt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)
    if color_ticklabels:    mplt.setp([ax.get_xticklabels(), ax.get_yticklabels()], color=color)

# more on modifying frames here: http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html

def remove_frame_and_background(ax=None):
    """ Remove plot frame, including x/y ticks and ticklabels, and the background. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    # this removes the frame, but also the background
    ax.set_frame_on(False)
    mplt.xticks([], [])
    mplt.yticks([], [])
    mplt.draw()

def remove_frame(ax=None):
    """ Remove plot frame, including x/y ticks and ticklabels. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    mplt.xticks([], [])
    mplt.yticks([], [])
    mplt.draw()

def remove_half_frame(ax=None):
    """ Remove the top and right sides of the plot frame, including x/y ticks. """
    # source: http://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib#925141
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    mplt.draw()


################################ OTHER ADDITIONS ###################################

######### color lists
# Trying to make lists of reasonably distinct colors, for various color numbers

# this is the html basic 16-color set from https://en.wikipedia.org/wiki/Web_colors#HTML_color_names, with some changes: different order, with no white, with yellow changed to gold because yellow is too pale, blue changed to dodgerblue because blue is too dark.
colors_15a = "black silver gray dodgerblue red lime gold aqua fuchsia maroon navy green olive teal purple".split()

# sets I just made up, trying to get them to be easy to distinguish
colors_6a = 'black blue darkturquoise forestgreen mediumorchid red'.split()
colors_6b = 'blue darkturquoise forestgreen darkorange mediumorchid red'.split()
colors_7a = 'black blue darkturquoise forestgreen darkorange mediumorchid red'.split()
colors_10a = 'black gray blue darkturquoise forestgreen darkgoldenrod saddlebrown darkorange red mediumorchid'.split()
colors_12a = 'black gray blue darkturquoise forestgreen yellowgreen darkgoldenrod darkorange red darkred mediumorchid blueviolet'.split()
colors_20a = "green springgreen yellowgreen darkturquoise teal skyblue dodgerblue blue navy coral red darkred plum magenta darkviolet darkorange darkgoldenrod darkgray black dimgray".split()
colors_20b = "black dimgray darkgray navy blue dodgerblue skyblue darkturquoise teal green springgreen yellowgreen darkgoldenrod darkorange coral red darkred magenta plum darkviolet ".split()
colors_21a = "black dimgray darkgray navy blue dodgerblue skyblue darkturquoise teal green springgreen yellowgreen darkgoldenrod gold darkorange coral red darkred magenta plum darkviolet ".split()

# MAYBE-TODO write functions to sort color-lists by hue (rainbow) and by brightness?

######### colormaps (continuous color ranges to use for heatmaps)

def all_colormaps_image(save_as='~/computers_and_programming/colormaps_all_matplotlib.png', close=False):
    """ Make an image with all the colormaps shown; optionally save to file. """
    # code from http://matplotlib.org/examples/pylab_examples/show_colormaps.html
    a = numpy.linspace(0, 1, 256).reshape(1,-1)
    a = numpy.vstack((a,a))
    # Get a list of the colormaps in matplotlib.  Ignore the ones that end with '_r' because these are simply 
    #  reversed versions of ones that don't end with '_r'
    maps = sorted(m for m in mplt.cm.datad if not m.endswith("_r"))
    nmaps = len(maps) + 1
    fig = mplt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i,m in enumerate(maps):
        ax = mplt.subplot(nmaps, 1, i+1)
        mplt.axis("off")
        mplt.imshow(a, aspect='auto', cmap=mplt.get_cmap(m), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], m, fontsize=10, horizontalalignment='right')
    mplt.draw()
    if save_as:     
        savefig(os.path.expanduser(save_as))
    if close:
        mplt.close()

def add_colormap_to_matplotlib(cmap_data, cmap_name, add_reverse_too=True):
    """ Add a colormap to the matplotlib colormap list. """
    # Based on code in https://github.com/dmcdougall/matplotlib/blob/85a26dca55c788c97bbf312b0732db222df55530/lib/matplotlib/cm.py
    # Basically the point is to get matplotlib.cm.get_cmap(cmap_name) to return a cmap instance.
    cm = matplotlib.cm
    cm.datad[cmap_name] = cmap_data
    cm.cmap_d[cmap_name] = cm._generate_cmap(cmap_name, cm.LUTSIZE)
    if add_reverse_too:
        reverse_name = cmap_name + '_r'
        matplotlib.cm.datad[reverse_name] = matplotlib.cm._reverse_cmap_spec(cmap_data)
        cm.cmap_d[reverse_name] = cm._generate_cmap(reverse_name, cm.LUTSIZE)

# CMRmap, from https://github.com/matplotlib/matplotlib/pull/496/files
#  (I think it should be in matplotlib already, but it's not in mine!)
_CMRmap_data = {'red'   : ( (0.000, 0.00, 0.00),
                            (0.125, 0.15, 0.15),
                            (0.250, 0.30, 0.30),
                            (0.375, 0.60, 0.60),
                            (0.500, 1.00, 1.00),
                            (0.625, 0.90, 0.90),
                            (0.750, 0.90, 0.90),
                            (0.875, 0.90, 0.90),
                            (1.000, 1.00, 1.00) ),
                'green' : ( (0.000, 0.00, 0.00),
                            (0.125, 0.15, 0.15),
                            (0.250, 0.15, 0.15),
                            (0.375, 0.20, 0.20),
                            (0.500, 0.25, 0.25),
                            (0.625, 0.50, 0.50),
                            (0.750, 0.75, 0.75),
                            (0.875, 0.90, 0.90),
                            (1.000, 1.00, 1.00) ),
                'blue':   ( (0.000, 0.00, 0.00),
                            (0.125, 0.50, 0.50),
                            (0.250, 0.75, 0.75),
                            (0.375, 0.50, 0.50),
                            (0.500, 0.15, 0.15),
                            (0.625, 0.00, 0.00),
                            (0.750, 0.10, 0.10),
                            (0.875, 0.50, 0.50),
                            (1.000, 1.00, 1.00) )}
add_colormap_to_matplotlib(_CMRmap_data, 'CMRmap')

### trying to make a green_heat colormap, same as gist_heat but green instead of red. 
# that's NOT just white-green-black - there's yellow involved too!  Maybe if I switched red with green in that map, I'd get a good green one... But they're functions rather than value lists, so it's complicated... I could try it with 'hot' instead, that has numbers, but it's not as nice.

def switch_colormap_colors(source_cmap_name, color_cycle, new_cmap_name, add_reverse_too=True):
    """ Make a colormap that's source_cmap_name with colors switched according to color_cycle, and save it as new_cmap_name. 
    
    Color_cycle must be a list of 2 or 3 colors (out of blue, green, red):
     - if it's 2 colors, those two will be switched and the third left constant
     - if it's 3 colors, move color1 data to color2, color2 data to color3, and color3 data to color1.
    """
    cm = matplotlib.cm
    all_colors = set('blue green red'.split())
    for color in color_cycle:
        if color not in all_colors:
            raise ValueError("Illegal color %s in color_cycle - should be one of (%s)!"%(color, ', '.join(all_colors)))
    if len(color_cycle) < 2:
        raise ValueError("Color_cycle must contain 2 or 3 distinct colors! Passed %s: %s!"%(len(color_cycle),', '.join(color_cycle)))
    if len(set(color_cycle)) != len(color_cycle):
        raise ValueError("All colors in color_cycle must be different! (%s, %s)"%(color1, color2))
    source_cmap_data = dict(cm.datad[source_cmap_name])
    # make a new copy of the source cmap data for the new cmap (so that colors not in color_cycle stay the same).
    new_cmap_data = dict(source_cmap_data)
    for from_color, to_color in zip(color_cycle, color_cycle[1:] + [color_cycle[0]]):
        new_cmap_data[to_color] = source_cmap_data[from_color]
    add_colormap_to_matplotlib(new_cmap_data, new_cmap_name, add_reverse_too)

# Make all the different-color variants of the gist_heat colormap (or just some of them if some are commented out)
switch_colormap_colors('gist_heat', 'red green'.split(), 'wp_GnYl_heat')
switch_colormap_colors('gist_heat', 'red blue'.split(), 'wp_BuCy_heat')
switch_colormap_colors('gist_heat', 'red green blue'.split(), 'wp_GnCy_heat')
switch_colormap_colors('gist_heat', 'red blue green'.split(), 'wp_BuPu_heat')
switch_colormap_colors('gist_heat', 'green blue'.split(), 'wp_RdPu_heat')

# Making other colormap variants based on gist_heat - single-tone ones this time, but similar otherwise.
gist_heat = mplt.cm.datad['gist_heat']
gh_red = gist_heat['red']
gh_green = gist_heat['green']
gh_blue = gist_heat['blue']
add_colormap_to_matplotlib({'red': gh_blue, 'green': gh_red, 'blue': gh_blue}, 'wp_green1_heat')
add_colormap_to_matplotlib({'red': gh_blue, 'green': gh_red, 'blue': gh_red}, 'wp_teal1_heat')
add_colormap_to_matplotlib({'red': gh_green, 'green': gh_red, 'blue': gh_green}, 'wp_green2_heat')
add_colormap_to_matplotlib({'red': gh_green, 'green': gh_red, 'blue': gh_red}, 'wp_teal2_heat')

# Explanation of making colormaps (esp. LinearSegmentedColormap): http://matplotlib.org/examples/pylab_examples/custom_cmap.html
# Another link on modifying colormaps: http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations

# Getting a bit more complicated by taking fractions of the values (I have to admit I don't quite understand how it works)
# Notes: try plotting gh_red, gh_green and gh_blue in the -0.1 to 1.1 range... I THINK what's going on is that only x in the 0-1 range are considered, and that values below 0 count as 0, and values above 1 count as 1.  000 is black, 111 is white.  So:
# - to make colors lighter, keep the value at 1 the same but raise the values for lower x; 
# - to make colors darker straightforwardly, keep the value at 0 the same but lower the values for higher x. 
# - to make colors "greyer", so that they can get "darker" but still start at black and end at white, I basically want to add a greyscale gradient to all colors, I think...  Anyway, those are just based on vague ideas and trial-and-error:
add_colormap_to_matplotlib({'red': gh_blue, 'green': lambda x: x, 'blue': gh_blue}, 'wp_green3_heat')
add_colormap_to_matplotlib({'red': gh_green, 'green': lambda x: x, 'blue': gh_blue}, 'wp_GnYl2_heat')
add_colormap_to_matplotlib({'red': gh_blue, 'green': lambda x: x, 'blue': lambda x: x}, 'wp_teal3_heat')
# teal4 is basically an average between green3 and teal3:
#  (the last lambda is basically "lambda x: (min(0,gh_green(x)) + x)/2", but in a way that works if x is a numpy array)
#  (might be better to use power functions instead of all those segmented lines, by the way... gh_green is about x**3, gh_blue about x**7.)
add_colormap_to_matplotlib({'red': gh_blue, 'green': lambda x: x, 'blue': lambda x: (numpy.int8(gh_blue(x)>0)*gh_blue(x) + x)/2}, 'wp_teal4_heat')
# MAYBE-TODO all these functions are linear anyway, so I think I could pretty well emulate this with LinearSegmentedColormap and it'd be easier...

def transform_colormap_indices(source_cmap_name, start, end, new_cmap_name, add_reverse_too=True):
    """ Make a colormap that's just the slice of source_cmap_name between start and end; save it as new_cmap_name.

    NOT PERFECT: currently the nearest indices outside the chosen range will just get moved to start/end positions without scaling!

    This will only work if the values in the source cmap color dict are (position, color1, color2) tuples, not functions!
    """
    cm = matplotlib.cm
    if not 0 <= start < end <= 1:
        raise ValueError("Start/end must satisfy 0 <= start < end <= 1! %s, %s are bad."%(start, end))
    source_cmap_data = dict(cm.datad[source_cmap_name])
    for val in source_cmap_data.values():
        try:                len(val)
        except TypeError:   
            raise ValueError("Source_cmap values must be (position, color, color) lists! Found %s instead."%val)
    # make new colormap, transforming all positions 
    new_cmap_data = {}
    for key, pos_color_list in source_cmap_data.items():
        curr_data = []
        for position, color1, color2 in pos_color_list:
            new_position = min(max(0, position-start) / (end-start), 1)
            curr_data.append((new_position, color1, color2))
            # TODO this is wrong - I need to scale the indices properly (letting them go below 0 and above 1), then remove the unneeded ones, and scale the COLOR VALUES in the nearest outside ones instead of just shifting them, since that squishes the edge colormap! Like, if I end up with (-.5, 0, 0) and (.5, 1, 1), the first should change into (0, .5, .5), NOT (0, 0, 0) like currently.
        # this may have resulted in multiple positions of 0 or 1 - get rid of those
        N_zeros = len([1 for pos,_,_ in curr_data if pos==0])
        if N_zeros>1:   del curr_data[:N_zeros-1]
        N_ones = len([1 for pos,_,_ in curr_data if pos==1])
        if N_ones>1:    del curr_data[-N_ones+1:]
        # make sure the result is sane, add to new_cmap_data
        assert curr_data[0][0] == 0
        assert curr_data[-1][0] == 1
        new_cmap_data[key] = curr_data
    # make a new copy of the source cmap data for the new cmap (so that colors not in color_cycle stay the same).
    add_colormap_to_matplotlib(new_cmap_data, new_cmap_name, add_reverse_too)

### Making different versions of the jet colormap, starting with lighter blue and ending with medium red instead of dark red
# First just take the middle and remove the sides: I want something that starts around 1/7 of jet, and ends at about 10/11.
transform_colormap_indices('jet', 1/7, 9.5/11, 'wp_light_jet1')
transform_colormap_indices('jet', 4/10.5, 9.5/11, 'wp_light_jet2')

# MAYBE-TODO add other nice colormaps, from elsewhere or mine:
#  - https://www.mathworks.com/matlabcentral/fileexchange/26026

################## OLD FUNCTIONS, moved from general_utilities.py

def plot_function_by_window_size(data,window_size_list,function,figsize=(),title='',xlabel='',xlim=None,ylim=None,yticks=None,ylabel=''):
    """ Return a plot of function on data, with a subplot for each window_size_list value. """
    if figsize: fig = mplt.figure(figsize=figsize)
    else:       fig = mplt.figure()
    for i in range(len(window_size_list)):
        window_size = window_size_list[i]
        mplt.subplot(len(window_size_list),1,i+1)
        if i==0 and title:                          mplt.title(title)
        if i==len(window_size_list)-1 and xlabel:   mplt.xlabel(xlabel)
        function_values, indices = function(data, window_size, return_indices=True)
        mplt.plot(indices, function_values, 'b,',linestyle='none') 
        if ylim:    mplt.ylim(ylim[0],ylim[1])
        if xlim:    mplt.xlim(xlim[0],xlim[1])
        if yticks:  mplt.yticks(yticks[0],yticks[1] if len(yticks)>1 else yticks[0])
        if ylabel:  mplt.ylabel(ylabel%window_size)
    return fig
    # TODO how would one even test this?


### Convert data into "linlog" scale for plotting (linear for some small range around 0, log after that)
# NOTE: there used to be an old in-progress implementation of this here, called convert_data_to_linlog, but I removed it on 2012-05-30 -  the right thing to do is use the "symlog" scale in xscale/yscale, with the appropriate linthreshx/linthreshy (also allows negative values):
# - http://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
# - http://matplotlib.sourceforge.net/examples/pylab_examples/symlog_demo.html


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    print "NO TESTS FOR THIS FILE"
    #unittest.main()
