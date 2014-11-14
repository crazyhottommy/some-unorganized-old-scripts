# this script demonstrates the usage of IntervalTree from the bx-python library
# https://bitbucket.org/james_taylor/bx-python/src/8600c7055ed08e87c6951e65ccce6771b38d0acf/lib/bx/intervals/intersection.pyx?at=default
# see post here for Interval Trees  http://informatics.malariagen.net/2011/07/07/using-interval-trees-to-query-genome-annotations-by-position/
# http://blog.nextgenetics.net/?e=45
# http://hackmap.blogspot.com/2008/11/python-interval-tree.html
# http://psaffrey.wordpress.com/2011/04/
# http://www.biostars.org/p/2244/
# Interval Tree data structure makes the query of Genomic regions very fast once built (it takes time and memory to buid the tree)
# 2013.10.08 Tommy Tang




import csv
from bx.intervals.intersection import IntervalTree, Interval


def index_gtf(gtf_file_path):
    # dictionary mapping chromosome names to interval trees
    genome = dict()
    #parse the annotations file (Gtf) and build the interval trees
    with open(gtf_file_path, "r") as annotations_file:
        reader = csv.reader(annotations_file, delimiter = '\t')
        for row in reader:
            if len(row) == 9 and not row[0].startswith('##'):
                seqid = row[0]
                start = int(row[3])
                end  = int(row[4])
                tree = None
                # build one interval tree per chromosome 
                if seqid in genome:
                    tree = genome[seqid]
                else:
#first time we've encoutered this chromosome, creat an interval tree
                    tree = IntervalTree()
                    genome[seqid] = tree
#index the feature
                tree.add(start, end, tuple(row))
    return genome


    
    # useage annotations=  genome['chr1'].find(590000,590100)   



    # another python packages gffutils https://github.com/daler/gffutils
    # one can parse the gtf file with HTSeq GFF_Reader
    # I believe biopython also has module to do that.
    # there are just way too many packages there, know what you want, choose one fit your need
    # do not re-invent the wheels ( I am not that good...even I want to)


    ## use the pybedtools http://pythonhosted.org/pybedtools/
import sys
from bedtools import BedTool
location_file = sys.argv[1]
genes_file = sys.argv[2]

for nearest in BedTool(location_file).closest(genes_file):
    print nearest
