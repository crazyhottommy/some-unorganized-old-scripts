# this script is used to extract the promoter regions that interact with a distal enhancer
# in this case, the distal HREs based on the RNA-polll CHIA-PET data
# we have three files: RNA polll CHIA-PET, HRE bed file and the promoter region file

#  RNApoll anchor1                              RNApoll anchor2
# ###############---------------------------------################
#          ###
#          HRE                                                       |---------------->
#                                                             _______|  
#                                                                   TSS                                      
# 10.14.13 modified to use the GenomicArrayOfsets class from the HTSeq package
# bx-python has an Interval Tree class, pybedtools can also be used.
# I am not doing simply interval intersection, the CHIA-PET data has two anchors(intervals)

          

import HTSeq
import itertools


class Interaction:
    '''a new class to represent the genomic interaction data, two \
    intervals interact with
    each other.
    the GenomicInterval object is from the HTSeq package'''
    def __init__(self, iv1, iv2, tag_number):
        self.iv1 = iv1
        self.iv2 = iv2
        self.tag_number = tag_number

    def __str__(self):
        rep = "interaction object\n"
        rep += "two interacting intervals: " + str(self.iv1) + " and "\
               + str(self.iv2) + "\n" + "with tag number" + str(self.tag_number)
        return rep
        
    def get_the_other_iv(self, interval_to_eval ):
        if self.iv1==interval_to_eval:
            return self.iv2
        elif self.iv2==interval_to_eval:
            return self.iv1
        else:
            return None

interaction_pairs = set()
with open("/home/tommy/CHIA_PET/pol2_intra_interaction_tag2.txt", "r") as f:
    for line in f:
        line_list = line.split()
        iv1 = HTSeq.GenomicInterval( line_list[0], int(line_list[1]),\
                                     int(line_list[2]), ".")
        iv2 = HTSeq.GenomicInterval( line_list[3], int(line_list[4]),\
                                     int(line_list[5]), ".")
        tag_number = line_list[6]
        interaction_pair = Interaction(iv1, iv2, tag_number)
        interaction_pairs.add(interaction_pair)

# build the Genomic Array

ga=HTSeq.GenomicArrayOfSets("auto", stranded=False)
for interaction_pair in interaction_pairs:
    ga[interaction_pair.iv1] += interaction_pair.iv2
    ga[interaction_pair.iv2] += interaction_pair.iv1
  
    
HRE_ivs = {} # initiate an empty dictionary, HRE_id is the key, iv is the value

with open("/home/tommy/CHIA_PET/HIF1_ChIP-seq_pvalue1e-4_peaks.bed","r") as f:
    for line in f:
        line_list = line.split()
        HRE_iv = HTSeq.GenomicInterval( line_list[0], int(line_list[1]),\
                                        int(line_list[2]), ".")
        HRE_id = line_list[3]
        HRE_ivs[HRE_id]=HRE_iv

# have a look at the dictionary
for key,value in itertools.islice(HRE_ivs.items(),10):
    print key,value
    
# access the genomic array by intervals(HREs)
# total_set contains all the anchors with the other coresponding anchor
# overlapping with HRE     
#############################################################      
total_set= set()
for HRE_id, HRE_iv in HRE_ivs.items():
    biggest_set = None
    for iv, step in ga[HRE_ivs[HRE_id]].steps():       
        if biggest_set is None:
            biggest_set = step.copy()
        else:
            biggest_set.update(step)
    total_set.update(biggest_set)
len(total_set)
# Now it is much faster than just brute comparing each Interval.
##############################################################    

ofile=open("/home/tommy/CHIA_PET/HBS_overlap_ChIA_PET.txt","w")
for HRE_id, HRE_iv in HRE_ivs.items():
    biggest_set = None
    for iv, step in ga[HRE_ivs[HRE_id]].steps():
        if biggest_set is None:
            biggest_set = step.copy()
        else:
            biggest_set.update(step)
    if len(biggest_set) > 0:
        #print biggest_set
        for anchor in biggest_set:
            #print anchor
            ofile.write(HRE_id + "\t" + HRE_iv.chrom  + "\t" + str(HRE_iv.start) + "\t" \
                   + str(HRE_iv.end) + "\t" + anchor.chrom +"\t" + str(anchor.start) +"\t" + str(anchor.end) +"\n")
ofile.close()

    

gtffile = HTSeq.GFF_Reader("/home/tommy/Datasets/hg19_UTR_exon.gtf")

for feature in itertools.islice(gtffile, 100):
    if feature.type == "5UTR":
        print feature.attr["gene_id"], feature.attr["transcript_id"], \
          feature.iv.start_d_as_pos

# buid a genomic array to hold the tss plus upstream and downstream
tsspos = HTSeq.GenomicArrayOfSets("auto", stranded = False)
upstream = 3000
downstream = 1000
for feature in gtffile:
    if feature.type == "5UTR":
        p = feature.iv.start_d_as_pos
        if p.strand == "+":
            window = HTSeq.GenomicInterval( p.chrom, p.pos - upstream,\
            p.pos + downstream, ".")
        else:
            window = HTSeq.GenomicInterval( p.chrom, p.pos - downstream,\
            p.pos + upstream, ".")
        tsspos[window] += feature.attr["gene_id"]
        # add gene name to this window 

# count how many genes overlap with the other anchor
gene_names = set()
for the_other_anchor in total_set:
    biggest_set = None
    for iv, step in tsspos[the_other_anchor].steps():
        if biggest_set is None:
            biggest_set = step.copy()
        else:
            biggest_set.update(step)
    gene_names.update(biggest_set)

    
ofile1= open("genes_with_distal_HRE.txt","w")

for gene in sorted(gene_names):
    ofile1.write(gene+'\n')

ofile1.close()

ofile2 = open("genes_with_promo_HRE.txt","w")
for gene in sorted(promoter_genes):
    ofile2.write(gene + "\n")

             
if __name__ == "__main__":
    main()
