#! /usr/bin/env python


#this program reads in the Refgene file downloaded from UCSC
#extract the chr_name and TSS position information for downstream analysis

import re # import the regular expression module 


'''the ifile looks like this:
790	NR_026775	chr6	+	26924771	26991753	26991753	26991753	3	26924771,26936236,26991032,	26924962,26936301,26991753,	0	LINC00240	unk	unk	-1,-1,-1,
216	NM_007188	chr7	+	150725509	150744869	150725602	150742436	16	150725509,150730691,150731359,150731591,150731810,150732667,150732968,150733159,150733630,150737354,150737584,150737914,150738185,150739047,150741057,150742295,	150725697,150731004,150731515,150731686,150731916,150732829,150733054,150733257,150733736,150737388,150737721,150738009,150738319,150739195,150741308,150744869,	0	ABCB8	cmpl	cmpl	0,2,0,0,2,0,0,2,1,2,0,2,1,0,1,0,'''


ifile=open("/home/tommy/Datasets/refGene_hg19.txt", "r")
ofile=open("/home/tommy/Datasets/refgene_clean.txt","r+")
for line in ifile:
    line_list=line.split()
    if line_list[3]=="+":    # if plus strand, TSS is at column 5
        indices=(1,2,3,4,12)
        newline=  "\t".join([line_list[i] for i in indices]) #extract column 2,3,4,5,13  list comprehension 
        ofile.write(newline+"\n")
    else:
        indices=(1,2,3,5,12) # if minus strand, TSS is at column 6
        newline= "\t".join([line_list[i] for i in indices])
        ofile.write(newline+"\n")


print "Done"
ifile.close()
ofile.close()

''' after it was done the ofile looks like this:
NR_026775	chr6	+	26924771	LINC00240
NM_007188	chr7	+	150725509	ABCB8
NM_001242928	chr14	+	74353317	ZNF410
NM_197964	chr7	+	139025877	C7orf55
NM_197962	chr1	-	193074608	GLRX2
NM_197961	chr15	-	65810035	DPP8
NR_024233	chr11	+	129872518	LINC00167
NM_007109	chr6_cox_hap2	+	2641059	TCF19'''

ofile2=open ("/home/tommy/Datasets/refgene_clean.txt","r")
for line in ofile2:
    line_list=line.split()
    SearchString="(\w+\d+)(_.+)"   #find chromosome names like chr6_cox_hap2
    Result=re.search(SearchString,line_list[1])
    if Result is not None:
        print Result.group(1) 
        print Result.group(2)
    else:
        pass
    
 
ofile2.close() 
