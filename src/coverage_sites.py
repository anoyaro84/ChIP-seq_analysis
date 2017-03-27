"""Retrieve coverage of bam files given a bed file and bamfile.

Takes in a ChIP-seq data in bam format and genomic locations in bed format,
and produces a matrix contains read count profiles in genomic locations.
The output file is generated with xarray.


Usage:
   coverage_sites.py [options] <bed_file> <netcdf_out> <bam_file>...

Options:
    --window=<window-size>  Size of genomic locations [default: 100].
    --threads=<thred-num>   Number of threads for calculating coverage [default: 1].
    --Ref_ver=<ref_ver>     Reference genome version [default: hg19].
    --fragment=<frag-size>  Size of fragment size by which each read is extended toward 3'end (used in metaseq function) [default: None].
    --binsize=<bin-size>    Number of bins to which coverages in each genomic coordiates are summurized. If None is given, number of bins is the same as window size [default: None].
"""

from docopt import docopt
import pybedtools
import xarray as xr
import ast


# function for obtaining center of genomic locations
def midpoint_generator(bedfile):
    from pybedtools.featurefuncs import midpoint
    for interval in bedfile:
        yield midpoint(interval)


# function for calculating coverage
def coverage(Bedfile, Bamfiles, Nproc, bins=None, fragSize=None):
    import os
    import metaseq
    import numpy as np

    ip_array = []
    for Bamfile in Bamfiles:
        print("Calculating coverages from : " + Bamfile)
        ip_signal = metaseq.genomic_signal(Bamfile, 'bam')
        ip_array.append(ip_signal.array(Bedfile, bins=bins, fragment_size=fragSize, processes=Nproc))

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

    #printing information
    print("Reading coverage of coordinates specified by: " + bedfile)
    print("With +/- " + str(WinSize) + "bp of window size")
    print("Reference genome :" + genome_ver)
    print("Number of threads: " + str(Nthread))

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
    array = coverage(Sites, Bamfiles=bamfiles, Nproc=Nthread, bins=BinSize, fragSize=FragSize)
    print("Calculating coverage was completed")

    #save output
    Out = xr.Dataset({'Coverage': xr.DataArray(array, dims = ['Sample', 'Coordinate', 'Position'],
        coords = {'Sample':bamfiles})})

    Outfile = arguments['<netcdf_out>']
    print('Saving output to :' + Outfile)
    Out.to_netcdf(str(Outfile))

