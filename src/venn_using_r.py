"""Produce pairwise venn-diagrom from given bed files.

Takes two bed files and produces a venn-diagram.


Usage:
   venn_using_r.py [options] <out_figure> <bed_file>...

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
import pybedtools as pb


if __name__ == '__main__':
    # reading arguments
    arguments = docopt(__doc__)
    bedfiles = arguments['<bed_file>']
    names = [arguments['--name1'], arguments['--name2']]
    colors = [arguments['--color1'], arguments['--color2'], arguments['--color3']]

    OutName = arguments['<out_figure>']

    pb.contrib.venn_maker.venn_maker(
            beds = bedfiles,
            names = names,
            additional_args=['euler.d=TRUE',
                'scaled=TRUE'],
            figure_filename=OutName,
            run=True
            )

