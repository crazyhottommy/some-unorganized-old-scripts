#! /usr/bin/env python

#this script coun chars, words, lines, paragraphs from file(s)   02/14/2013

import sys
import glob

Usage=''' 
this program counts the character, words, lines and paragrahs
of file or files. It supports wild card when specifying the files
Usage:
count_chars_words_lines_paragraphs.py *.txt
'''

if len(sys.argv)==1:
    print Usage
         
else: 
    FileList=sys.argv[1:]
    c, w, l, p = 0, 0, 0, 0
    for FileName in FileList:
        for File in glob.glob(FileName):      #wild card serach for files  "*.txt"
            s=open(File).read()
            wc= len(s),len(s.split()), len(s.split('\n')), \
                len(s.split('\n\n'))
            print '\t'.join(map(str, wc)), '\t' + File
            c, w, l, p = c+wc[0], w+wc[1], l+wc[2], p+wc[3]
    wc=(c, w, l, p)
    print '\t'.join(map(str, wc)), '\tTOTAL'  #map function returns a list 

   
    #s=sys.stdin.read()
    #wc= len(s), len(s.split()), len(s.split('\n')),\
    # len(s.split('\n\n'))
    #print '\t'.join(map(str, wc)), '\tSTDIN'
