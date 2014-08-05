'''
Created on Jul 25, 2014

@author: kunal
'''

from collections import OrderedDict
from math import log10

import pandas as pd

import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import pylab
import brewer2mpl
import seaborn as sns

sns.set_style("whitegrid")
sns.set("paper", {"lines.linewidth": 2,
                  "font.size": 10,
                  "axes.labelsize": 11,
                  "axes.xtick.labelsize": 10,
                  "axes.ytick.labelsize": 10})

# Axes size
item_labels_size = 10


def beautify(original_plot):
    def new_plot(*args, **kwargs):
        fig, ax = original_plot(*args, **kwargs)
        ax = adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # --- Added this line --- #
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontsize(item_labels_size)
        return fig, ax
    return new_plot


def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            # spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine
    return ax

@beautify
def count_box_plot(count_df, colors, ax=None, fig=None):
    groups = count_df.columns
    rects = OrderedDict()
    if not fig:
        fig = plt.figure(figsize=(6, 12))
    if not ax:
        ax = fig.add_subplot(111)

    width = 0.25  # the width of the bars
    for ix, group in enumerate(groups):
        values = count_df[group].values
        ind = range(0, len(values))
        rects[group] = ax.bar([val + (ix * width) for val in ind], values,
                              width, color=colors[group], ecolor='black')

    ax.set_xticks([val + width for val in ind])
    measurements = count_df.index
    ax.set_xticklabels(measurements)
    ax.grid(False, which='major', axis='x')
    legend = ax.legend(rects.values(), rects.keys(), loc=2,
                       prop={'size': item_labels_size})
    frame = legend.get_frame()
    frame.set_facecolor('w')
    ax.set_xlim(0 - width, ax.get_xlim()[1] + width)
    return fig, ax


def get_count_df(stats):
    variants = ['ins', 'del', 'mnp', 'snp']
    count_df = pd.concat([stat['counts'][variants] for stat in stats], axis=1)
    count_df = count_df.applymap(log10)
    count_df.index = ['Insertion', 'Deletion', 'MNP', 'SNP']
    return count_df


def plot_counts(fig, ax, count_df, colors):
    fig, ax = count_box_plot(count_df, colors, ax, fig)
    ax.set_ylabel('Number of Variants (log10)')  # , fontdict={'size':10})
    return fig, ax

def my_rstyle(ax):
    """Styles an axes to appear like ggplot2
    Must be called after all plot and axis manipulation operations have been carried out (needs to know final tick spacing)
    """
    # set the style of the major and minor grid lines, filled blocks

    ax.grid(True, 'major', color='gray', linestyle='-', linewidth=0.7, alpha=0.5)
    ax.patch.set_facecolor('white')
    ax.set_axisbelow(True)

    # restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)

    # remove the minor tick lines
    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)

    # only show bottom left ticks, pointing out of axis
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # set x, y limit

    if ax.legend_ <> None:
        lg = ax.legend_
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_alpha(0.5)
