"""Retrieve coverage of bam files given a bed file and bamfile.

Takes in a ChIP-seq data in bam format and genomic locations in bed format,
and produces a matrix contains read count profiles in genomic locations.
The output file is generated with xarray.


Usage:
   coverage_sites.py [options] <bed_file> <netcdf_out> <bam_file>...

Options:
    --isbigwig=<isbigwig>   Bigwig is provided instead of bam   [default: False].
    --isTable=<istable> Indicator parameter that list of bams are given in a table. If True, it will load the table and find IDs to locate bam or bigwig files [default: False].
    --colname_table=<colname>   Column name to identify IDs from the table, if isTable=True [default: IDs].
    --prefix_bam=<prefix>   Common prefix for bam_file to locate them [default: None].
    --suffix_bam=<suffix>   Common suffix for bam_file to locate them [default: None].
    --window=<window-size>  Size of genomic locations [default: 1000].
    --threads=<thred-num>   Number of threads for calculating coverage [default: 1].
    --Ref_ver=<ref_ver>     Reference genome version [default: hg19].
    --fragment=<frag-size>  Size of fragment size by which each read is extended toward 3'end (used in metaseq function) [default: None].
    --binsize=<bin-size>    Number of bins to which coverages in each genomic coordiates are summurized. If None is given, number of bins is the same as window size [default: None].
"""

from docopt import docopt
import pybedtools
import xarray as xr
import ast
import pandas as pd

# function for obtaining center of genomic locations
def midpoint_generator(bedfile):
    from pybedtools.featurefuncs import midpoint
    for interval in bedfile:
        yield midpoint(interval)


# function for calculating coverage
def coverage(Bedfile, Bamfiles, Nproc, bins=None, fragSize=None):
    import metaseq
    import numpy as np

    ip_array = []
    for Bamfile in Bamfiles:
        print("Calculating coverages from : " + Bamfile)
        ip_signal = metaseq.genomic_signal(Bamfile, 'bam')
        ip_array.append(ip_signal.array(Bedfile, bins=bins, fragment_size=fragSize, processes=Nproc))

    return np.asarray(ip_array)

# function for calculating 
def coverage_bw(Bedfile, BigWigs, bins=None):
    import subprocess as sp
    import numpy as np
    Bedfile.saveas('tmp_bed.bed')

    ip_array = []
    for BigWig in BigWigs:
        print("Calculating coverages from : " + BigWig)
        command = ['bwtool', 'extract', 'bed', 'tmp_bed.bed', BigWig, 'tmp.bed']
        print(' '.join(command))
        p = sp.Popen(command, stdout=sp.PIPE, stdin=sp.PIPE)
        p.wait()

        out = pybedtools.BedTool('tmp.bed')

        if bins == None:
            mat = np.vstack(
                    [np.array(str(a[-1]).split(',')) for a in out]
                    )
            mat[mat=='NA'] = '0'
            mat.astype(float)
            ip_array.append(mat)

    return np.asarray(ip_array)

# main 
if __name__ == '__main__':
    # reading arguments
    arguments = docopt(__doc__)
    bedfile = arguments['<bed_file>']
    bamfiles = arguments['<bam_file>']
    WinSize = int(arguments['--window'])
    Nthread = int(arguments['--threads'])
    genome_ver = arguments['--Ref_ver']
    IsBigWig = arguments['--isbigwig'] in ['True', 'true']
    Istable = arguments['--isTable'] in ['True', 'true']
    prefix = arguments['--prefix_bam']
    suffix = arguments['--suffix_bam']

    if str(prefix) == 'None':
        prefix = ''
    if str(suffix) == 'None':
        suffix = ''

    #printing information
    print("Reading coverage of coordinates specified by: " + bedfile)
    print("With +/- " + str(WinSize) + "bp of window size")
    print("Reference genome :" + genome_ver)
    print("Number of threads: " + str(Nthread))

    if Istable:
        tableName = bamfiles[0]
        colname = str(arguments['--colname_table'])
        print("loading bam files from the table: " + tableName)
        table = pd.read_csv(tableName, sep='\t', header=0)
        if colname in table.columns:
            bamfiles = [prefix+id+suffix for id in table[colname].tolist()]
        else:
            raise ValueError("The column name for IDs (" + colname + ") is not found in the table!")
    else:
        bamfiles = [prefix+id+suffix for id in bamfiles]


    FragSize = None
    if ast.literal_eval(arguments['--fragment']) != None:
        FragSize = int(arguments['--fragment'])
        print("Fragment size specified: " + str(FragSize) + "bp stretch in the 3'end direction for each read")

    BinSize = None
    if ast.literal_eval(arguments['--binsize']) != None:
        BinSize = int(arguments['--binsize'])
        print("Bin size specified: Signals in each genomic coordinates will be summarized to " + str(BinSize) + " bins")


    #identify sites
    Sites = pybedtools.BedTool(bedfile)
    MidSites = pybedtools.BedTool(midpoint_generator(Sites))
    Sites = MidSites.slop(b=WinSize, genome=genome_ver)

    #calculate covrage
    if not IsBigWig:
        array = coverage(Sites, Bamfiles=bamfiles, Nproc=Nthread, bins=BinSize, fragSize=FragSize)
    else:
        array = coverage_bw(Sites, bamfiles, BinSize)
    print("Calculating coverage was completed")

    #save output
    Out = xr.Dataset({'Coverage': xr.DataArray(array, dims = ['Sample', 'Coordinate', 'Position'],
        coords = {'Sample':bamfiles})})

    Outfile = arguments['<netcdf_out>']
    print('Saving output to :' + Outfile)
    Out.to_netcdf(str(Outfile))

