#!/usr/bin/env python
"""
    Given 2 files, creates a 2-way Venn diagram of intersections using matplotlib.  

    Numbers are placed on the diagram.  If you don't have matplotlib installed.
    try venn_gchart.py to use the Google Chart API instead.

    The values in the diagram assume:

        * unstranded intersections
        * no features that are nested inside larger features
"""

import argparse
import sys
import os
import pybedtools
import matplotlib
matplotlib.use('Agg')

def venn_mpl(a, b, colors=None, outfn=None, labels=None):
    """
    *a* and *b* are filenames to BED-like files.

    *colors* is a list of matplotlib colors for the Venn diagram circles.

    *outfn* is the resulting output file.  This is passed directly to
    fig.savefig(), so you can supply extensions of .png, .pdf, or whatever your
    matplotlib installation supports.

    *labels* is a list of labels to use for each of the files; by default the
    labels are ['a','b']
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError:
        sys.stderr.write('matplotlib is required to make a Venn diagram with %s\n' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)

    if colors is None:
        colors = ['r','b']

    radius = 6.0
    center = 0.0
    offset = radius / 2

    if labels is None:
        labels = ['a','b']

    circle_a = Circle(xy = (center-offset, center+offset), radius=radius, edgecolor=colors[0], label=labels[0])
    circle_b = Circle(xy = (center+offset, center+offset), radius=radius, edgecolor=colors[1], label=labels[1])


    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)

    for circle in (circle_a, circle_b):
        circle.set_facecolor('none')
        circle.set_linewidth(3)
        ax.add_patch(circle)

    ax.axis('tight')
    ax.axis('equal')
    ax.set_axis_off()


    kwargs = dict(horizontalalignment='center')

    # Unique to A
    ax.text( center-2*offset, center+offset, str((a - b).count()), **kwargs)

    # Unique to B
    ax.text( center+2*offset, center+offset, str((b - a).count()), **kwargs)

    # A and B 
    ax.text( center, center+2*offset-0.5*offset, str((a + b).count()), **kwargs)

    ax.legend(loc='best')

    fig.savefig(outfn)

    plt.close(fig)

def main():
    """Create a 3-way Venn diagram, using matplotlib"""
    op = argparse.ArgumentParser(description=__doc__, prog=sys.argv[0])
    op.add_argument('-a', help='File to use for the left-most circle')
    op.add_argument('-b', help='File to use for the right-most circle')
    op.add_argument('--labels',
                  help='Optional comma-separated list of '
                       'labels for a and b', default='a,b')
    op.add_argument('--colors', default='r,b',
                  help='Comma-separated list of matplotlib-valid colors '
                       'for circles a and b.  E.g., --colors=r,b')
    op.add_argument('-o', default='out.png', 
                  help='Output file to save as.  Extension is '
                       'meaningful, e.g., out.pdf, out.png, out.svg.  Default is "%(default)s"')
    op.add_argument('--test', action='store_true', help='run test, overriding all other options.')
    options = op.parse_args()

    reqd_args = ['a','b']
    if not options.test:
        for ra in reqd_args:
            if not getattr(options,ra):
                op.print_help()
                sys.stderr.write('Missing required arg "%s"\n' % ra)
                sys.exit(1)

    venn_mpl(a=options.a, b=options.b, 
             colors=options.colors.split(','),
             labels=options.labels.split(','), 
             outfn=options.o)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
