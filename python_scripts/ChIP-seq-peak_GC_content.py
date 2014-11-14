#!  /usr/bin/env python 

# this program extracts sequences from ChIP-seq peaks (bed intervals after
# MACS peak call) by using the coordinates in pygr worldbase module.
# bedtools can also do the same job, but need to download the whole genome sequences in local
# plus 3kb upstream and downstream sequences. ( all the sequences are the same length, 6001 bp )
# calculate the GC content 
# x axis [-3000,-2999,.....0.......2999,3000]
# y axis is the %GC content for each position
# 05/05/13 by  Tommy Tang



def GC_content(ifile):

    # ifile 1 is a bed file, with first column is the chromosome name, second column is the start, third column is the end
    # it uses pygr worldbase module to extract genomic sequences of cooresponding sequences
    # and the first base of each sequence form a new sequence, the second base of each sequence for the second new sequence
    # GC content is calculated for those new sequences.

    from pygr import worldbase

    hg19=worldbase.Bio.Seq.Genome.HUMAN.hg19() # creat hg19 to represent the whole genome

    halfwinwidth=3000

    SeqList=[]

    f=open(ifile, "r")
    for line in f:
        LineList= line.split()
        Seq= hg19[LineList[0]][ int(LineList[1])-halfwinwidth : int(LineList[2])+halfwinwidth+1] 
        # genome sequences can be accessed by a dictionary  with LineList[0]='chr1' as a key
        # and then slice the sequence according to the start and end coordinates 
        # LineList[1]----start   LineList[2]--------end
        SeqList.append(str(Seq))
        # Seq is a pygr sequence object, str(Seq) gives the actual sequence

    #print SeqList

    NewSeqList = []  # creat a new seq list, NewSeq in index1 is a sequence comprised of
    # the first base of all the old Seq. NewSeq in index 2 is a sequence comprised of
    # the second base of all the old Seq...

    for i in range(len(SeqList[0])):
        NewSeq=""
        for Seq in SeqList:
            NewSeq += Seq[i]
        NewSeqList.append(NewSeq)
    
    #print NewSeqList[0]

    #calculate the GC% 

    GC_content=[]

    for NewSeq in NewSeqList:
        TotalG= NewSeq.count("G") + NewSeq.count("g") # seqence my contain lower case, it also may contain "N"
        TotalC= NewSeq.count("C") + NewSeq.count("c")
        GC_content_NewSeq = 100 * float ((TotalG + TotalC))/ len(NewSeq)
        GC_content.append(GC_content_NewSeq)
    return GC_content
    


halfwinwidth=3000
    
from matplotlib import pyplot

x_axis= range(-halfwinwidth,halfwinwidth+1)
GC1=GC_content("/home/tommy/CHIA_PET/HRE/HIFsites_overlap_polll_position.bed")
GC2=GC_content("/home/tommy/CHIA_PET/HRE/HIFsites_not_overlap_polll_position.bed")


line1= pyplot.plot(x_axis, GC1 , "ro", label= "involved in looping")
line2= pyplot.plot(x_axis, GC2, "b^", label= "not involved in looping")

pyplot.xlabel("distance related to HRE peak summit bp")
pyplot.ylabel("GC%")
pyplot.title("GC content at HRE")


pyplot.show()
    
