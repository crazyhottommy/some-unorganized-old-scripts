# This code was modified from the tss plot code, that it can plot any other ChIP-seq signal
# at other genomic positions. In this case, it is the ERE. ERa  ChIP-seq data
# GSM594602
# is available, peaks were called by MACS, generated a bed file. the middle point
# of each peak is used as the center of the plot
# 04/18/13 modified for ERE

def TSS_Profile(ifile1,ifile2):
    '''read in three files, ifile1 is the sortedbamfile prepared by samtool
    ifile2 is the Genomic position  file with three columns: chr, position, strand'''
    
    import HTSeq
    import numpy
    import itertools

    sortedbamfile=HTSeq.BAM_Reader(ifile1)
   
    ERE_file=open(ifile2)
    halfwinwidth=3000
    fragmentsize=200

        
    EREpos=set() 
    for line in ERE_file:
        linelist=line.split()
        EREpos.add(HTSeq.GenomicPosition(linelist[0],int(linelist[1]),'.'))# creat Genomic position objects by HTSeq
        # if there is a blank line, linelist[1] will get an index out of range error
        # make sure no blank lines, you can write a yeild Generator no_blanklines()
        #   def nonblan_lines(f):
        #       for l in f:
        #           line=l.rstrip()
        #           if line:
        #               yield line
            
    for ERE in itertools.islice(EREpos,10):
        print ERE  # print out 10 ERE postions 

   
            
              
    profile=numpy.zeros(2*halfwinwidth, dtype='i')
    for p in EREpos:
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

input_bamfile=("/home/tommy/CHIA_PET/Mcf7Jund.sorted.bam")
    
counts=aligned_counts(input_bamfile)

halfwinwidth=3000

profile1=TSS_Profile(input_bamfile,\
                     "/home/tommy/CHIA_PET/ER/IHM001F/ER_hg19_involved_in_looping.position")
profile1_normalized= (profile1 *1000000.0/counts) /4884.0

profile2=TSS_Profile(input_bamfile,\
                     "/home/tommy/CHIA_PET/ER/IHM001F/ER_hg19_not_involved_in_looping.position")
profile2_normalized= (profile2*1000000.0/counts) /9509.0               


from matplotlib import pyplot
import numpy


line1=pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), profile1_normalized, color="red",label="invovled in looping")
line2=pyplot.plot(numpy.arange(-halfwinwidth, halfwinwidth), profile2_normalized, color="blue",label="not invovled in looping")

pyplot.legend()
pyplot.xlabel("distance related to ERE bp")
pyplot.ylabel("tag density")
pyplot.title("MCF7 JunD enrichment around ERE ")


pyplot.show()
