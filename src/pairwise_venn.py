"""Produce pairwise venn-diagrom from given bed files.

Takes two bed files and produces a venn-diagram.


Usage:
   pairwise_venn.py [options] <out_figure> <bed_file>...

Options:
    --name1=<name1>  Name for the first group [default: group1].
    --name2=<name2>   Name for the second group [default: group2].
    --name3=<name3>   Name for the second group [default: group3].
    --color1=<color1>     Color for the first group [default: white].
    --color2=<color2>     Color for the second group [default: white].
    --color3=<color3>     Color for the shared group [default: white].
"""

from docopt import docopt
import matplotlib
matplotlib.use('Agg')
from matplotlib_venn import venn2
from matplotlib_venn import venn2_circles
from matplotlib import pyplot as plt
import pybedtools


def two_way_venn(bedfiles, names, colors):
    Site1 = pybedtools.BedTool(bedfiles[0])
    Site2 = pybedtools.BedTool(bedfiles[1])
    Int = Site1.intersect(Site2, wa=True)

    Sets=(len(Site1)-len(Int), len(Site2)-len(Int), len(Int))

    fig = plt.figure(figsize=(5,5))
    v = venn2(subsets=Sets,  set_labels = names)
    v.get_patch_by_id('10').set_color(colors[0])
    v.get_patch_by_id('01').set_color(colors[1])
    v.get_patch_by_id('11').set_color(colors[2])

    c = venn2_circles(subsets=Sets, linestyle='solid')

    return fig

def three_way_venn(bedfiles, names, colors):
    Site1 = pybedtools.BedTool(bedfiles[0])
    Site2 = pybedtools.BedTool(bedfiles[1])
    Site3 = pybedtools.BedTool(bedfiles[2])

    Int12 = Site1.intersect(Site2, wa=True)
    Int23 = Site2.intersect(Site3, wa=True)
    Int31 = Site3.intersect(Site1, wa=True)
    Int123 = Int12.intersect(Site3, wa=True)

#    Sets=
# main

if __name__ == '__main__':
    # reading arguments
    arguments = docopt(__doc__)
    bedfiles = arguments['<bed_file>']
    names = [arguments['--name1'], arguments['--name2']]
    colors = [arguments['--color1'], arguments['--color2'], arguments['--color3']]

    fig = two_way_venn(bedfiles, names, colors)

    OutName = arguments['<out_figure>']
    print("Saving file: " + OutName)
    fig.savefig(OutName, dpi=100, bbox_inches="tight")


