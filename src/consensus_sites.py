"""Identify consensus binding sites givn multiple bed files.

Take multiple bed files and generate a consensus binding sites.

Usage:
    consensus_sites.py [options] <outfile> <bedfiles>...

Options:
    --thrs=<thrs>    An integer for threshold of frequency for consensus sites. If not specified, half of the number of bed files is used [default: Half].
"""

from docopt import docopt
import ast
import numpy as np
import pybedtools as pb


def create_consensus(Bedfiles, thr):
    pyBedfiles = dict()
    for i,Bedfile in zip(range(len(Bedfiles)), Bedfiles):
        print("Obtaining " + Bedfile)
        pyBedfiles[Bedfile] = pb.BedTool(Bedfile)

    Union = pyBedfiles[Bedfile]
    for Bedfile in Bedfiles:
        Union = Union.cat(pyBedfiles[Bedfile])

    Counts=np.zeros((len(Union),))
    for Bedfile in Bedfiles:
        intersect = Union.intersect(pyBedfiles[Bedfile], c=True)
        Counts = Counts + np.array([int(interval[3]) for interval in intersect])

    Consensus = Union.at(np.where(Counts>thr)[0])
    return Consensus


if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    Bedfiles = arguments['<bedfiles>']
    Outfile = arguments['<outfile>']
    thrs = arguments['--thrs']

    if thrs == "Half":
        thrs = np.floor(len(Bedfiles)/2)
    else:
        try:
            thrs = float(thrs)
        except ValueError:
            print("Should give numeric values for threshold.")
            raise

    Consensus = create_consensus(Bedfiles, thrs)
    print("Saving output at " + Outfile)
    Consensus.saveas(Outfile)

