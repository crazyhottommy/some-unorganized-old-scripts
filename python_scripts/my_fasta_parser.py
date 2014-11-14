#! /usr/bin/env python

#this program process the fasta file, read the file line by line
# based on the ">", process one record at a time. There is no 
# blank line between records.



ifile=open("/home/tommy/playground/FPexamples.fta","r")


while 1:#this while loop is for the accumulating of Records
        Records=[]
        line=ifile.readline()
        if not line.startswith(">"):
            raise TypeError ("Not a FASTA file: %r" % line)
        title=line[1:]   # get the sequence name
        sequence_lines=[]
        while 1: 
            line=ifile.readline().rstrip() #this while loop is for the accumulating of one Record
            if line.startswith(">"):
                break
            sequence_lines.append(line)
            last_position=ifile.tell()
            if line=="":
                break
        Record=">"+ title + "".join(sequence_lines)
        print Record
        ifile.seek(last_position) #repostion one line back
        Records.append(Record)
        
        if line=="": #end of the file
            break
print Records       
