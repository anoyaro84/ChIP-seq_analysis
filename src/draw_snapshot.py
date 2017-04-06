"""Plot readcount profile of defined genomic locations.

Takes in a ChIP-seq data in bam format and genomic locations in bed format,
and produces a figure shows read count profiles in the top genomic locations.


Usage:
   draw_snapshot.py [options] <bed_file> <figure_out> <bam_file>...

Options:
    --red=<red>  red of RGB color code [default: 1].
    --green=<green>  green of RGB color code [default: 0].
    --blue=<blue>  blue of RGB color code [default: 0].
"""

from docopt import docopt
import pybedtools
import matplotlib

def draw_snapshot(sites, bamfiles, color="black", min_y=30, Nsite=5):
    import metaseq
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import Grid


    Nsites_use = min(Nsite, len(sites))

    # take read counts from samples
    ip_signals = [None]*len(bamfiles)
    ip_arrays = [None]*len(bamfiles)
    for i, bamfile in zip(range(len(bamfiles)), bamfiles):
        ip_signals[i] = metaseq.genomic_signal(bamfile, 'bam')
        ip_arrays[i] = [None]*Nsites_use
        for k in range(Nsites_use):
            ip_arrays[i][k] = ip_signals[i].local_coverage(sites[k])

    # draw figure
    fig = plt.figure(figsize=(100, 20))

    grid = Grid(fig, 142, nrows_ncols=(len(bamfiles), Nsites_use),
            axes_pad=0.05, direction="row",
            add_all=True, share_all=False, label_mode="all")

    ymaxs = [None]*len(bamfiles)

    for i in range(len(bamfiles)):
        for k in range(Nsites_use):
            grid[i*Nsites_use+k].bar(ip_arrays[i][k][0],
                    ip_arrays[i][k][1], color=color, edgecolor=color)
            xmin, xmax, ymin, ymax = grid[i*Nsites_use+k].axis()
            ymaxs[i] = max(ymaxs[i], ymax)


    for i in range(len(bamfiles)):
        for k in range(Nsites_use):
            xmin, xmax, ymin, ymax = grid[i*Nsites_use+k].axis()
            grid[i*Nsites_use+k].axis([xmin,xmax,ymin, ymaxs[i]])
            grid[i*Nsites_use+k].get_xaxis().set_visible(False)
            grid[i*Nsites_use+k].get_yaxis().set_visible(False)
            grid[i*Nsites_use+k].annotate(
                    bamfiles[i].split('/')[-1].split('.')[0]
                    +" [0-"+str(max(ymax, min_y))+"]",
                    xy=(0,1), xytext=(10, -10),
                    va='top', xycoords='axes fraction',
                    textcoords='offset points',
                    fontsize=25)
            if i==0:
                chrom = str(sites[k]).split('\t')[0]
                start = str(sites[k]).split('\t')[1]
                end = str(sites[k]).split('\t')[2].split('\n')[0]
                grid[i*Nsites_use+k].set_title(
                        "Location: " + chrom + " " + start + "-" + end
                        )

    return fig



if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    bedfile = arguments['<bed_file>']
    bamfiles = arguments['<bam_file>']
    outfile = arguments['<figure_out>']
    red = float(arguments['--red'])
    green = float(arguments['--green'])
    blue = float(arguments['--blue'])

    #printing information
    print("Reading coverage of coordinates specified by: " + bedfile)


    matplotlib.use('Agg')

    sites=pybedtools.BedTool(bedfile)
    fig = draw_snapshot(sites, bamfiles, color=(red, green, blue),
            min_y=30, Nsite=5)
    fig.savefig(outfile, dpi=100, bbox_inches="tight")

