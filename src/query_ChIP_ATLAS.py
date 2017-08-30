"""Query ChIP-ATLAS data. 

Usage:
    query_ChIP_ATLAS [options] <celltype> <prefix>

Options:
    --table=<table> ChIP-atlas data table. In default, it will take latest table.   [default: http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab]
    --datatype=<dtype>  Data type to be obatiend. Either of bigwig, bed and both.    [default: both]
    --threshold=<thr> Threshold used for selecting peak lists.    [default: 05]
    --filtermeta=<filter>   Include features with specific string pattern.  [default: No]
"""

import docopt
import pandas as pd
import urllib as url
import os

if __name__=='__main__':
    arguments = docopt.docopt(__doc__)

    tableName = str(arguments['--table'])
    print("Parsing table: " + tableName)

    table = pd.read_csv(tableName, sep='\t', usecols=range(0,13),
            header=None, names = ['ID', 'assembly', 'AntigenClass',
                'Antigen', 'CellTypeClass', 'CellType', 'CellDescription',
                'ProcessLog', 'Title', 'MetaDataByAuthor',
                'MetaData2', 'MetaData3', 'MetaData4'
                ])

    celltype = str(arguments['<celltype>'])
    print("Filtering table element by cell type: " + celltype)
    ind_cellType = table.CellType == celltype

    if not str(arguments['--filtermeta']).strip() == 'No':
        filter = str(arguments['--filtermeta']).strip().split(',')
        print("Filtering table with the following filters: " + ','.join(filter))

        ind_filter = table.CellType == -1
        for string in filter:
            print(string)
            ind_filter = ind_filter | (table.MetaDataByAuthor.str.contains(string) | \
                table.MetaData2.str.contains(string) | \
                table.MetaData3.str.contains(string) | \
                table.MetaData4.str.contains(string))
        ind_final = ind_cellType & ind_filter
    else:
        ind_final = ind_cellType

    print(str(ind_final.sum()) + " entries remained")
    table = table[ind_final]
    
    prefix = str(arguments['<prefix>'])
    print("Saving filtered table at " + prefix + "/table.tab")

    table.to_csv(prefix + "/table.tab", sep='\t')


    datatype = str(arguments['--datatype'])
    thr = str(arguments['--threshold'])
    # obtaining bed files
    if datatype in  ['both', 'bed']:
        if not os.path.exists(prefix + '/bed'):
           os.makedirs(prefix + '/bed')
        print("obtaining bed files with Q value threshold 10e-" + thr + " from the table.")
        for ID, genome  in zip(table.ID, table.assembly):
            address = 'http://dbarchive.biosciencedbc.jp/kyushu-u/' +\
                    genome + "/eachData/bed" + thr + '/' + ID + "." + thr + '.bed'
            local = prefix + '/bed/' + ID + '.' + thr + '.bed'
            if not os.path.isfile(local):
 #               try:
                url.urlretrieve(address, local)
#                except url.ContentTooShortError:
#                    print(ID)
#                    raise

    # obtaining bigwig
    if datatype in  ['both', 'bigwig']:
        if not os.path.exists(prefix + '/bigwig'):
            os.makedirs(prefix + '/bigwig')
        print("obtaining bigwig files from the table.")
        for ID, genome  in zip(table.ID, table.assembly):
            address = 'http://dbarchive.biosciencedbc.jp/kyushu-u/' +\
                    genome + "/eachData/bw/" + ID + ".bw"
            local = prefix + '/bigwig/' + ID + '.bw'
            if not os.path.isfile(local):
 #               try:
                url.urlretrieve(address, local)
 #               except url.ContentTooShortError:
 #                   print(ID)

