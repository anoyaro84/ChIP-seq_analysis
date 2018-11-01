from subprocess import PIPE, Popen
from urllib.request import urlopen, Request
from bs4 import BeautifulSoup
import urllib.request
import os

def get_paths(IDs, PATH_DATA, ext="bam"):
    # open specified URL
    Files = dict()
    FullPath = dict()

    if isinstance(PATH_DATA, str):
        PATH_DATA = PATH_DATA.split()

    for PATH in list(PATH_DATA):
        print("reading path:")
        print(PATH)

        response = Request(PATH)
        html = urlopen(response)

        print("..success")
        print("")

        # find the link for files
        print("find specified files...")
        soup = BeautifulSoup(html.read(), "html.parser")

        for link in soup.find_all('a'):
            name = link.get('href')
            for ID in IDs:
                if ID in name and name.endswith(ext):
                    Files[ID] = name
                    FullPath[Files[ID].split('.')[0]] = PATH + name

    # check if all specified files can be found
    for ID in IDs:
        if not ID in Files.keys():
            raise Exception(str(ID) + "is not found at the path!")

    return Files, FullPath


def query_SRX(GSM_ACC, path_ncbitoolkit='/home/NFS/users/yo.kim/lib/softwares/edirect/'):
    p1 = Popen([path_ncbitoolkit + "esearch", "-db", "gds", "-query", '"' + GSM_ACC  +'"'],stdout=PIPE)
    p2 = Popen([path_ncbitoolkit + "efetch", "-format", "docsum"], stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    p3 = Popen([path_ncbitoolkit + "xtract", "-pattern", "ExtRelations", "-element", "TargetObject"],
                stdin = p2.stdout, stdout = PIPE)
    p2.stdout.close()
    tmp = p3.communicate()[0].decode("utf-8").split('\n')
    return tmp[len(tmp)-2]

def query_SRR(GSM_ACC, path_ncbitoolkit='/home/NFS/users/yo.kim/lib/softwares/edirect/'):
    p1 = Popen([path_ncbitoolkit + "esearch", "-db", "sra", "-query", '"' + GSM_ACC  +'"'],stdout=PIPE)
    p2 = Popen([path_ncbitoolkit + "efetch", "-format", "docsum"], stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    p3 = Popen([path_ncbitoolkit + "xtract", "-pattern", "DocumentSummary", "-element", "Run@acc"],
                stdin = p2.stdout, stdout = PIPE)
    p2.stdout.close()
    tmp = p3.communicate()[0].decode("utf-8").split('\t')
    return list(set(tmp)-set(['']))


