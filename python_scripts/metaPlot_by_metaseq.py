# From http://www.biostars.org/p/83800/:
 
# "What I want to do is to plot reads of my histone marks (in bam file)
# around TSS with CpG and TSS without CpG (Essentially a coverage profile)."
 
 
# to install metaseq and dependencies:
# 1. get the metaseq source:
#
#     git clone https://github.com/daler/metaseq.git
#
# 2. change to dev branch:
#
#    cd metaseq
#    git checkout v0.5dev
#
# 3. Install numpy and cython first if you don't already have them.
#
#    pip install numpy
#    pip install cython
#
# 4. Install metaseq using pip.  If you don't have a lot of the scientific 
#    Python packages or genomic Python packages, this may take a while.
#  
#    pip install .
 
 
import metaseq
import pybedtools
import numpy as np
from matplotlib import pyplot as plt
 
bam = metaseq.genomic_signal('Mcf7Max.sorted.bam', 'bam')
cpg = pybedtools.BedTool('cpg.bed')
tss = pybedtools.BedTool('HIF_sites_invovled_in_looping_not_at_promoter.bed')
 
# extend by 5 kb up/downstream
tss = tss.slop(b=5000, g=pybedtools.chromsizes('hg19'))
 
tss_with_cpg = tss.intersect(cpg, u=True)
tss_without_cpg = tss.intersect(cpg, v=True)
 
# change this to as many CPUs as you have in order to run in parallel
processes = 1
 
# each read will be extended 3' to a total size of this many bp
fragment_size = 200
 
# the region +/-5kb around each TSS will be split into a total of 100 bins,
# change as needed
bins = 100
 
x = np.linspace(-5000, 5000, bins)
 
# most of the work happens here
y1 = bam.array(tss, bins=bins, processes=processes, fragment_size=fragment_size)
y2 = bam.array(tss_without_cpg, bins=bins, processes=processes, fragment_size=fragment_size)
 
 
plt.plot(x, y1.mean(axis=0), label='with cpg')
plt.plot(x, y2.mean(axis=0), label='without cpg')
plt.legend(loc='best')
plt.xlabel('Distance from TSS (bp)')
plt.ylabel('Mean read density')
plt.show()
