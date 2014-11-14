#! /usr/bin/env python
# fast fastq file parser, extract ~10,000 reads with their sequences and quality 
#from a fastq file (~70milion reads) based on names  of the reads. HengLi in Princeton https://github.com/lh3/seqtk
# has a wrapper for this kind of task.
# the following code is from http://www.biostars.org/p/10353/  07/03/13
# this demonstrates the usage of set, a data structure that is much faster than list
# when you have two files, put the information from one file into a container, loop over the other file.
# usage: cat file.fastq | python parse_fastq.py id_file.txt > selected.fastq

import sys
# get filename from parameter
idfile = sys.argv[1]

# load ids in a set with  a set comprehension 

ids = set( x.strip() for x in open(idfile) )

# read the fastq file
handle = sys.stdin

while ids:
    #parse fastq
    idline = handle.readline()
    seq   = handle.readline()
    spacer = handle.readline()
    quals = handle.readline()

    id_name = idline[:-1] # except the newline \n
    if id_name in ids:
        #print fastq
        sys.stdout.write( '%s%s%s%s%' % ( idline, seq, spacer, quals) )
        ids.remove(id_name)
