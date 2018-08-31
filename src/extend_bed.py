"""Extend sites in bed files to given size.

Takes in a bed file and size of sites and save into a bed file.

Usage:
   extend_bed.py [options] <bed_file> <out_file>

Options:
    --window=<window-size>  Size of genomic locations [default: 1000].
    --Ref_ver=<ref_ver>     Reference genome version [default: hg19].
"""

from docopt import docopt
import pybedtools

# function for obtaining center of genomic locations
def midpoint_generator(bedfile):
    from pybedtools.featurefuncs import midpoint
    for interval in bedfile:
        yield midpoint(interval)


# main 

if __name__ == '__main__':
    # reading arguments
    arguments = docopt(__doc__)
    bedfile = arguments['<bed_file>']
    WinSize = int(arguments['--window'])
    genome_ver = arguments['--Ref_ver']
    outfile = arguments['<out_file>']

    print("Identifying +/- " + str(WinSize) + "bp from center of sites in " + str(bedfile))
    print("Reference genome :" + genome_ver)

    #identify sites
    Sites = pybedtools.BedTool(bedfile)
    MidSites = pybedtools.BedTool(midpoint_generator(Sites))
    Sites = MidSites.slop(b=WinSize, genome=genome_ver)

    print("Save at " + str(outfile))
    Sites.saveas(outfile)

