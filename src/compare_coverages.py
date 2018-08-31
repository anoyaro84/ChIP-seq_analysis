"""Generate correlation heatmap given number of bam and bed files.

Take multiple bam and bed files, and generate heat map indicating similairty between coverages.If multiple bedfiles are given, union is taken. To give multiple bedfiles, they are separated with comma.

Usage:
    compare_coverages.py [options] <outfile> <bedfiles> <bamfiles>...

Options:
    --color=<col_scheme>    Color scheme for matplotlib. For more options, visit http://matplotlib.org/examples/color/colormaps_reference.html [default: Reds].
    --plottype=<plottype>  Type of visualization. Either one of heatmap, LTheatmap (lower-triangle heatmap) and clustermap (heatmap with clustering) [default: heatmap].
    --grad_max=<grad_max>   upper bound of color gradient [default: 1.0].
    --grad_min=<grad_min>   lower bound of color gradient [default: 0.0].
    --labels=<labels>   Names to be used for samples. Should be delimited with comma. If not specified, it is derived from file names. [default: None].
    --cores=<n_cores>   number of cores to be used [default: 5].
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

def readcount(site, pathtobam, order):
    read = pb.BedTool.coverage(site, pathtobam)
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


def create_array(Bedfiles, Bamfiles, max_workers=15):
    mat = np.zeros((len(Bedfiles), len(Bedfiles)))
    PyBedfiles = dict()
    colname = [None] * len(Bedfiles)

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

    # measuring total number of reads
    print("Calculating FPKM")
    Nreads = [None] * len(Bamfiles)
    futures = []
    with cf.ThreadPoolExecutor(max_workers=max_workers) as e:
        for order in range(len(Bamfiles)):
            futures.append(e.submit(flagstats, Bamfiles[order], order))

        for future in cf.as_completed(futures):
            read, order = future.result()
            Nreads[order] = float(read)


    counts = np.zeros((len(UnionSite), len(Bamfiles)))
    for i in range(len(Bamfiles)):
        for j, site in enumerate(bamreads[i]):
            counts[j,i] = int(site.name)*(1000000000/Nreads[i])/float(site.length)

    mat = np.corrcoef(np.log2(counts.T+1))


    return mat, colname


# main function start

if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    Bedfiles = str(arguments['<bedfiles>']).split(',')
    Bamfiles = arguments['<bamfiles>']
    Outfile = arguments['<outfile>']
    Cmap = arguments['--color']
    PlotType = arguments['--plottype']
    vmax = float(arguments['--grad_max'])
    vmin = float(arguments['--grad_min'])
    N_cores = int(arguments['--cores'])
    Names = arguments['--labels']

    if Names != "None":
        Names = Names.split(',')
        if len(Names) != len(Bamfiles):
            raise ValueError('Wrong number of labels are given.' +
                    ' Number of labels=' + str(len(Names)) +
                    ' and Number of bam files=' + str(len(Bamfiles)))

            if not PlotType in ['heatmap', 'LTheatmap','clustermap']:
                raise ValueError('Unknown plot type is given: ' + PlotType)


    print("Type of heatmap: " + PlotType)
    print("Color gradient range: " + str(vmin) + "~" + str(vmax))
    print("Color scheme: " + str(Cmap))
    print("Number of cores given: " + str(N_cores))

    # create co-occupancy map
    mat, name = create_array(Bedfiles, Bamfiles, max_workers=N_cores)
    if Names == "None":
        Names = name


    # produce heatmap
    print("Producing heatmap at " + Outfile)
    pdmat = pd.DataFrame(data=mat, columns=Names)
    mask = np.zeros_like(pdmat, dtype=np.bool)

    if PlotType == 'LTheatmap':
        mask[np.triu_indices_from(mask)] = True
    else:
        mask[:] = False

    #plt.rcParams['font.family'] = 'Arial'
    #plt.rcParams['font.size'] = 10

    #print(pdmat)

    fig = plt.figure(figsize=(10,10))
    if PlotType in ['heatmap', 'LTheatmap']:
        ax = sns.heatmap(pdmat, mask=mask, cmap=Cmap, square=True, linewidths=1,
                vmin=vmin, vmax=vmax, xticklabels=Names, yticklabels=Names)
        plt.xticks(rotation='vertical')
        plt.yticks(rotation='horizontal')
        fig.savefig(Outfile, dpi=100)

    else:
        cm = sns.clustermap(pdmat, linewidths=1, cmap=Cmap,
                vmin=vmin, vmax=vmax)
        plt.setp(cm.ax_heatmap.get_yticklabels(), rotation='horizontal')
        plt.setp(cm.ax_heatmap.get_xticklabels(), rotation='vertical')
        cm.savefig(Outfile, dpi=100)


