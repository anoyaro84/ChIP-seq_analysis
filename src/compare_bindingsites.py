"""Generate co-occupancy heatmap given number of bed files.

Take multiple bed files and generate heat map indicating similairty between binding sites.

Usage:
    compare_bindingsites.py [options] <outfile> <bedfiles>...

Options:
    --color=<col_scheme>    Color scheme for matplotlib. For more options, visit http://matplotlib.org/examples/color/colormaps_reference.html [default: Reds].
    --plottype=<plottype>  Type of visualization. Either one of heatmap, LTheatmap (lower-triangle heatmap), clustermap (heatmap with clustering) and matrix (save text file with similarity measured) [default: heatmap].
    --measure=<sim_measure> Similarity measure. Can either be correlation (pearson correaltion in co-occupancy profile) or overlap (overlap rate in intervals) [default: correlation].
    --grad_max=<grad_max>   Upper bound of color gradient [default: 1.0].
    --grad_min=<grad_min>   Lower bound of color gradient [default: 0.0].
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

def create_array(Bedfiles, overlap=True):
    mat = np.zeros((len(Bedfiles), len(Bedfiles)))
    PyBedfiles = dict()
    colname = [None] * len(Bedfiles)

    # create bedfiles
    print("Collecting bamfiles... ")
    for i,Bedfile in zip(range(len(Bedfiles)),Bedfiles):
        print("Obtaining " + Bedfile)
        PyBedfiles[Bedfile] = pb.BedTool(Bedfile)
        colname[i] = (Bedfile.split('/')[-1]).split('.')[0]

    # measure overlapping ratio
    if overlap:
        print("Calculate overlap ratio between bed files")
        for i, Bedfile1 in zip(range(len(Bedfiles)),Bedfiles):
            for j, Bedfile2 in zip(range(len(Bedfiles)),Bedfiles):
                if i==j:
                    mat[i,j]=1
                elif j>i:
                    mat[i,j] = len(PyBedfiles[Bedfile1].intersect(PyBedfiles[Bedfile2]))
                    mat[i,j] = mat[i,j]/min(len(PyBedfiles[Bedfile1]), len(PyBedfiles[Bedfile2]))
                    mat[j,i] = mat[i,j]

    # measure similarity in co-occupancy profiles
    else:
        print("obtaining unions of bed files")
        Union = pb.BedTool(Bedfile)
        for Bedfile in Bedfiles:
            Union = Union.cat(PyBedfiles[Bedfile])

        print("Total " + str(len(Union)) + " sites were identified")
        print("Identifing occupancy patterns of the files")

        occupancy = np.zeros((len(Union), len(Bedfiles)))

        for i, Bedfile in zip(range(len(Bedfiles)), Bedfiles):
            tmp = Union.coverage(PyBedfiles[Bedfile], wo=True)
            score = [int(site.name)!=0 for site in tmp]
            occupancy[:,i] = np.array(score)

        mat = np.corrcoef(occupancy.T)


    return mat, colname


if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    Bedfiles = arguments['<bedfiles>']
    Outfile = arguments['<outfile>']
    Cmap = arguments['--color']
    PlotType = arguments['--plottype']
    vmax = float(arguments['--grad_max'])
    vmin = float(arguments['--grad_min'])
    measure = arguments['--measure']

    if not PlotType in ['heatmap', 'LTheatmap', 'clustermap', 'matrix']:
        raise ValueError('Unknown plot type is given: ' + PlotType)

    if not measure in ['correlation', 'overlap']:
        raise ValueError('Unknown similarity measure is given: ' + measure)

    print("Type of heatmap: " + PlotType)
    print("Color gradient range: " + str(vmin) + "~" + str(vmax))
    print("Color scheme: " + str(Cmap))

    # create co-occupancy map
    mat, name = create_array(Bedfiles, overlap=measure == 'overlap')

    # produce heatmap
    pdmat = pd.DataFrame(data=mat, columns=name)
    mask = np.zeros_like(pdmat, dtype=np.bool)

    if PlotType == 'LTheatmap':
        mask[np.triu_indices_from(mask)] = True
    else:
        mask[:] = False

    #plt.rcParams['font.family'] = 'Arial'
    #plt.rcParams['font.size'] = 10

    #print(pdmat)

    if PlotType in ['heatmap', 'LTheatmap']:
        print("Producing heatmap at " + Outfile)
        fig = plt.figure(figsize=(10,10))
        ax = sns.heatmap(pdmat, mask=mask, cmap=Cmap, square=True, linewidths=1,
           vmin=vmin, vmax=vmax, xticklabels=name, yticklabels=name)
        plt.xticks(rotation='vertical')
        plt.yticks(rotation='horizontal')
        fig.savefig(Outfile, dpi=100)

    elif PlotType == 'clustermap':
        print("Producing heatmap at " + Outfile)
        fig = plt.figure(figsize=(10,10))
        cm = sns.clustermap(pdmat, linewidths=1, cmap=Cmap,
                  vmin=vmin, vmax=vmax)
        plt.setp(cm.ax_heatmap.get_yticklabels(), rotation='horizontal')
        plt.setp(cm.ax_heatmap.get_xticklabels(), rotation='vertical')
        cm.savefig(Outfile, dpi=100)
    else:
        print("Saving matrix at " + Outfile)
        np.savetxt(Outfile, mat, delimiter='\t')

