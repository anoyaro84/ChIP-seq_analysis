"""Generate a coverage scatter plot that compares two ChIP-seq data.

Take two bam files along with a bed file and generate a scatter plot. If multiple bedfiles are given (separated by ','), union of the sites are used.

Usage:
    scatter_coverage.py [options] <outfile> <bedfiles> <bamfiles>...

Options:
    --hlsites=<hlsites>    Obtain sties from a bed file and highlight the points in scatter plots that overlap with the sites. [default: None].
    --index_name=<index_name>    Index of feature name in the bed file for highlighting (specified by hlsites). [default: 4].
    --measure=<measure>   Coverage measures. FPKM or CPM. [default: FPKM].
    --title=<title>   Title of the plot. [default: ScatterPlot].
    --kind=<kind>   Type of plot, all options in jointplot of seaborn supported (e.g. reg, scatter) [default: scatter].
"""

from docopt import docopt
import matplotlib
matplotlib.use('Agg')
import ast
import numpy as np
import pybedtools as pb
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import concurrent.futures as cf
import matplotlib.patches as mpatches
import pysam

def readcount(site, pathtobam, order):
    read = pb.BedTool.coverage(site, pathtobam, counts=True)
    return order, read

def flagstats(bamfile, order):
    from subprocess import Popen, PIPE
    p = Popen(['samtools', 'flagstat', bamfile], stdout=PIPE)
    for line in p.stdout:
        if "mapped" in line:
            nmapped = line.rstrip().split()
            nmapped = int(nmapped[0])
            break
    return nmapped, order

def create_array(Bedfiles, Bamfiles, measure, max_workers=15):
    PyBedfiles = dict()
    colname = [None] * len(Bamfiles)

    # create bedfiles
    print("Collecting bedfiles... ")
    for i,Bedfile in zip(range(len(Bedfiles)),Bedfiles):
        print("Obtaining " + Bedfile)
        PyBedfiles[Bedfile] = pb.BedTool(Bedfile)

    for i, Bamfile in enumerate(Bamfiles):
        colname[i] = (Bamfile.split('/')[-1]).split('.')[0]

    print("Obtaining unions of binding sites")
    UnionSite = PyBedfiles[Bedfile]
    for Bedfile in PyBedfiles:
        UnionSite = UnionSite.cat(PyBedfiles[Bedfile])

    # reading bam reads
    bamreads = [None]*len(Bamfiles)
    futures = []
    print("Calculating read counts from bam files")
    with cf.ThreadPoolExecutor(max_workers=max_workers) as e:
        for order in range(len(Bamfiles)):
            futures.append(e.submit(readcount, UnionSite, Bamfiles[order], order))

        for future in cf.as_completed(futures):
            order, read = future.result()
            bamreads[order] = read

    counts = np.zeros((len(UnionSite), len(Bamfiles)))
    # measuring total number of reads

    Nreads = [pysam.AlignmentFile(bam).mapped for bam in Bamfiles]
    print(Nreads)

    if measure == "FPKM":
        print("Calculating FPKM")
        for i in range(len(Bamfiles)):
            for j, site in enumerate(bamreads[i]):
                counts[j,i] = np.log(1+float(site[-1])*(1000000000/Nreads[i])/float(site.length))

    elif measure == "CPM":
        print("Calculating CPM")
        for i in range(len(Bamfiles)):
            for j, site in enumerate(bamreads[i]):
                counts[j,i] = np.log(1+float(site[-1])*1000000.0/float(Nreads[i]))

    return counts, colname, UnionSite


def draw_scatter(x, y, xname, yname, sites, index_highlight, title, kind="scatter"):
    # generate figure
    f = plt.figure(figsize=(11,9))
    sns.set(style="white", color_codes=True)

    Table = pd.concat([pd.Series(x), pd.Series(y)], axis=1)
    Table.columns = [xname, yname]
    grid = sns.jointplot(xname, yname, data=Table, kind=kind, color='k')
   # grid.ax_joint.plot([-1,max(x)*1.1],[-1,max(y)*1.1], 'r--')
    grid.ax_marg_x.set_title(title)

    HL_names = list(set(index_highlight.values()))
    colors = sns.color_palette(n_colors=len(HL_names)).as_hex()
    col_dict = dict()
    for i, name in enumerate(HL_names):
        col_dict[name] = colors[i]

    for col,ind in zip(colors,index_highlight):
        grid.ax_joint.plot(x[ind], y[ind], marker= 'o', color = col_dict[index_highlight[ind]])


    markers = []
    for col in col_dict:
        markers.append(mpatches.Patch(color=col_dict[col], label=col))

    print(markers)

    plt.legend(handles=markers, loc=2)

    grid.ax_joint.set_xlabel(xname)
    grid.ax_joint.set_ylabel(yname)

    return grid


# site = sites (pybedtool) used for scatter plot
# hlsite = sites (pybedtool) to be highlighted
# index_name = index to feature names in hlsite
def index_hlsearch(sites, hlsites, index_name):
    index_highlight = dict()
    for hlsite in hlsites:
        hit = sites.all_hits(hlsite)
        if len(hit) > 0:
            hit = hit[0]
            for i,site in enumerate(sites):
                if site[0] == hit[0] and site[1] == hit[1] and site[2] == hit[2]:
                    index_highlight[i] = hlsite[index_name]
    return index_highlight


# main function start
if __name__ == '__main__':
    arguments = docopt(__doc__)
    Bedfiles = arguments['<bedfiles>'].split(',')
    Bamfiles = arguments['<bamfiles>']
    Outfile = arguments['<outfile>']
    measure = str(arguments['--measure'])
    kind = str(arguments['--kind'])
    title = str(arguments['--title'])

    # highlight
    hlsites = arguments['--hlsites']

    if measure not in ['FPKM', 'CPM']:
        raise   ValueError("Unknown measure: " + str(measure))

    print("Coverage measure: " + str(measure))
    print("Calculating coverages...")
    counts, colname, UnionSite = create_array(Bedfiles, Bamfiles, measure, max_workers=2)

    # identify sites to be highlighted:
    if hlsites == "None":
        print("No features are highlighted")
        index_highlight = dict()
    else:
        print("Features overlap with sites in " + str(hlsites) + " will be highlighted")
        hlsites = pb.BedTool(hlsites)
        index_name = int(arguments['--index_name'])
        index_highlight = index_hlsearch(UnionSite, hlsites, index_name)


    print("Producing scatter plot")
    fig = draw_scatter(counts[:,0], counts[:,1], colname[0], colname[1], UnionSite, index_highlight, title, kind=kind)
    print("Saving figure at :" + Outfile)
    fig.savefig(Outfile, dpi=100, bbox_inches="tight")


