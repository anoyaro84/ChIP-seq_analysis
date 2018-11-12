"""Generate coverage matrix given a bed files and multiple bam files.

Taking reference bed and multiple bam files, and construct a coverage matrix. If multiple bed files are given (delimited with ","), union of them will be taken.

Usage:
    construct_coverage_matrix.py [options] <outfile> <bedfiles> <bamfiles>...

Options:
    --measure=<measure>   Either Raw, FPKM or CPM used to calculate coverage. [default: FPKM].
    --cores=<cores>  Maximum number of jobs to be excuted in parallel. [default: 10].
    --sorted=<sorted>   To run memory-efficient coverage calculation. Pre-sorted bams are required [default: False].
"""

from docopt import docopt
import ast
import numpy as np
import pybedtools as pb
import pandas as pd
import concurrent.futures as cf
import os
from tempfile import NamedTemporaryFile as temp

def readcount(site, pathtobam, order, sorted):
    read = pb.BedTool.coverage(site, pathtobam, counts=True, sorted=sorted)
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


def create_array(Bedfiles, Bamfiles, measure='FPKM', max_workers=15, sorted=False):
    mat = np.zeros((len(Bedfiles), len(Bedfiles)))
    PyBedfiles = dict()

    # create bedfiles
    print("Collecting bedfiles... ")
    for i,Bedfile in zip(range(len(Bedfiles)),Bedfiles):
        print("Obtaining " + Bedfile)
        PyBedfiles[Bedfile] = pb.BedTool(Bedfile)

    if len(Bedfiles) > 1:
        print("Obtaining unions of binding sites")
        UnionSite = PyBedfiles[Bedfile]
        for Bedfile in PyBedfiles:
            UnionSite = UnionSite.cat(PyBedfiles[Bedfile])
    else:
        UnionSite = PyBedfiles[Bedfiles[0]]


    if sorted:
        print("Sorted bam/bed files are given. Extracting chromosome names")
        #tmpfile = UnionSite._tmp()
        #os.system('sort {0} -k1,1 -k2,2n > {1}'.format(UnionSite.fn, tmpfile))
        with temp('w') as f:
            command = """samtools view -H """ + Bamfiles[0] + """ | grep SQ | cut -f 2,3 | awk '{sub(/^SN:/,""); sub(/LN:/,""); print;}' >  """ + f.name
            os.system(command)
            UnionSite = UnionSite.sort(faidx=f.name)
        #UnionSite = pb.BedTool(tmpfile)

    # reading bam reads
    bamreads = [None]*len(Bamfiles)
    futures = []
    print("Calculating read counts from bam files")
    with cf.ThreadPoolExecutor(max_workers=max_workers) as e:
        for order in range(len(Bamfiles)):
            futures.append(e.submit(readcount, UnionSite, Bamfiles[order], order, sorted))

        for future in cf.as_completed(futures):
            order, read = future.result()
            bamreads[order] = read

    # measuring total number of reads
    if measure in ['FPKM', 'CPM']:
        print('Obtaining library depth..')
        Nreads = [None] * len(Bamfiles)
        futures = []
        with cf.ThreadPoolExecutor(max_workers=max_workers) as e:
            for order in range(len(Bamfiles)):
                futures.append(e.submit(flagstats, Bamfiles[order], order))

            for future in cf.as_completed(futures):
                read, order = future.result()
                Nreads[order] = float(read)

    print("Calculating " + measure)
    counts = np.zeros((len(UnionSite), len(Bamfiles)))
    for i in range(len(Bamfiles)):
        for j, site in enumerate(bamreads[i]):
            if measure == 'FPKM':
                counts[j,i] = np.log2(float(site[-1])*(1000000000/Nreads[i])/float(site.length)+1)
            elif measure == 'CPM':
                counts[j,i] = np.log2(float(site[-1])*(1000000000/Nreads[i])+1)
            else:
                counts[j,i] = int(site[-1])

    df =  pd.concat([UnionSite.to_dataframe(), pd.DataFrame(counts, columns=Bamfiles)], axis=1)

    return df


if __name__ == '__main__':
    # reading argument
    arguments = docopt(__doc__)
    Bedfiles = str(arguments['<bedfiles>']).split(',')
    Bamfiles = arguments['<bamfiles>']
    Outfile = arguments['<outfile>']
    sorted = str(arguments['--sorted']) in ['True', 'true']

    Measure = str(arguments['--measure'])
    Cores = int(arguments['--cores'])

    if not Measure in ['FPKM', 'CPM', 'Raw']:
        raise ValueError("Unknown measure: " + Measure + " , should be either of FPKM and CPM")

    df = create_array(Bedfiles, Bamfiles, measure=Measure, max_workers=Cores, sorted=sorted==True)
    df.to_csv(Outfile, sep='\t', header=True, index=False)


