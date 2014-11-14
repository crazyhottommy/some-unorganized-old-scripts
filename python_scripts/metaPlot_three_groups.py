# This code was modified from the tss plot code, that it can plot any other ChIP-seq signal
# at other genomic positions. In this case, it is the HRE. HIF1a ChIP-seq data
# is available, peaks were called by MACS in Galaxy, generated a bed file. the middle point
# of each peak is used as the center of the plot. all HREs are included
# 04/10/13

def TSS_Profile(ifile1,ifile2):
    '''read in three files, ifile1 is the sortedbamfile prepared by samtool
    ifile2 is the Genomic position  file with three columns: chr, position, strand'''
    
    import HTSeq
    import numpy
    import itertools

    sortedbamfile=HTSeq.BAM_Reader(ifile1)
   
    HRE_file=open(ifile2)
    halfwinwidth=3000
    fragmentsize=200

        
    HREpos=set() 
    for line in HRE_file:
        linelist=line.split()
        HREpos.add(HTSeq.GenomicPosition(linelist[0],int(linelist[1]),'.'))# creat Genomic position objects by HTSeq
        # if there is a blank line, linelist[1] will get an index out of range error
        # make sure no blank lines, you can write a yeild Generator no_blanklines()
        #   def nonblan_lines(f):
        #       for l in f:
        #           line=l.rstrip()
        #           if line:
        #               yield line
            
    for HRE in itertools.islice(HREpos,10):
        print HRE  # print out 10 HRE postions 

   
            
              
    profile=numpy.zeros(2*halfwinwidth, dtype='i')
    for p in HREpos:
        try:
            window=HTSeq.GenomicInterval(p.chrom, p.pos-halfwinwidth-fragmentsize,p.pos+halfwinwidth + fragmentsize,".")
            for almnt in sortedbamfile[window]:
                almnt.iv.length=fragmentsize
                if p.strand==".":
                    start_in_window=almnt.iv.start- p.pos +halfwinwidth
                    end_in_window  =almnt.iv.end  - p.pos +halfwinwidth
                else:
                    start_in_window=p.pos+halfwinwidth-almnt.iv.end
                    end_in_window =p.pos+halfwinwidth-almnt.iv.start 
                start_in_window=max(start_in_window,0)
                end_in_window=min(end_in_window, 2*halfwinwidth)
                if start_in_window >= 2*halfwinwidth or end_in_window <0:
                    continue
                profile[start_in_window : end_in_window] +=1
        except:
           continue
    return profile
    ifile1.close()
    ifile2.close()

    # it is very slow to count the reads in this way, better to use the samtools flagstat
def aligned_counts(ifile1):
    '''count how many alignments are aligned back to genome, ifile1 is a sorted bam file'''
    import HTSeq
    sortedbamfile= HTSeq.BAM_Reader(ifile1)
    aligned_counts=0
    unaligned_counts=0
    for almnt in sortedbamfile:
        if almnt.aligned:
            aligned_counts+= 1
        else:
            unaligned_counts+=1
    print "number of aligned tags of %s is %d " % (ifile1, aligned_counts)
    print "number of unaligned tags of %s is %d "% (ifile1, unaligned_counts)
    return aligned_counts

input_bamfile=("/home/tommy/Downloads/Mcf7pol2.sorted.bam")
    
    #counts=aligned_counts(input_bamfile)
counts=25850138
halfwinwidth=3000

profile1=TSS_Profile(input_bamfile,\
                     "/home/tommy/CHIA_PET/TSS_-3kb_+1kb/looping_distal_HRE.position")
profile1_normalized= (profile1 *1000000.0/counts) /903.0

profile2=TSS_Profile(input_bamfile,\
                     "/home/tommy/CHIA_PET/TSS_-3kb_+1kb/non_looping_distal_HRE.position")
profile2_normalized= (profile2*1000000.0/counts) /2989.0  

profile3=TSS_Profile(input_bamfile,\
                     "/home/tommy/CHIA_PET/TSS_-3kb_+1kb/HRE_at_promoter.position")
profile3_normalized = (profile3*1000000.0/counts)/952.0

from matplotlib import pyplot
import numpy


line1=pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), profile1_normalized, color="red",label="distal_looping_HRE")
line2=pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), profile2_normalized, color="blue",label="distal_non_looping_HRE")
line3=pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), profile3_normalized, color="green",label="promoter_HRE")
pyplot.legend()
pyplot.xlabel("distance related to HRE bp")
pyplot.ylabel("tag density")
pyplot.title("MCF7 RNApol2 around HRE ")


pyplot.show()
