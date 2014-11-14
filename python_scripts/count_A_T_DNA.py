DNASeq = raw_input("Enter a DNA sequence:")
print 'Sequence:', DNASeq
Seqlength=float(len(DNASeq))
print 'Sequence Length:', Seqlength

NumberA=DNASeq.count('A')
NumberC=DNASeq.count('C')
NumberG=DNASeq.count('G')
NumberT=DNASeq.count('T')
print 'A:', NumberA/Seqlength
print 'C:', NumberC/Seqlength
