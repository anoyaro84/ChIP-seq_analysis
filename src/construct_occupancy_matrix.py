"""Generate co-occupancy matrix given number of bed files.

Taking reference bed and the other multiple bed files, and construct a binary matrix indicating occupancy of the bed files in the reference regions.

Usage:
    construct_occupancy_matrix.py [options] <outfile> <refbed> <otherbed>...

Options:
"""

from docopt import docopt
import pybedtools as pb
import pandas as pd
import numpy as np


def co_occupancy(Refbed, Bedfiles):
   Reference = pb.BedTool(Refbed)

   print("Total " + str(len(Reference)) + " sites in the reference bed: " + Refbed)

   occupancy = np.zeros((len(Reference), len(Bedfiles)))
   for bed,i in zip(Bedfiles, range(len(Bedfiles))):
       print("Checking occupancy of " + bed + " in Reference bed")
       tmp = Reference.coverage(pb.BedTool(bed), wo=True)
       score = [int(site.name)!=0 for site in tmp]
       occupancy[:, i] = np.array(score)

   df = pd.concat([Reference.to_dataframe(), pd.DataFrame(occupancy)], axis=1)
   return df


if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    Refbed = arguments['<refbed>']
    Outmat = arguments['<outfile>']
    Bedfiles = arguments['<otherbed>']

    df = co_occupancy(Refbed, Bedfiles)

    print("Saving outcome at " + Outmat)
    df.to_csv(Outmat, sep='\t', header=False, index=False)
