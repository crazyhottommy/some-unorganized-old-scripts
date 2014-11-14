def TSV_index(file_name):
    '''this function reads in a tab-delimited file line by line, split
    each line into a list, and index each line as an element of a big
    list of the whole file'''
    
    ifile=open(file_name,"r")
    file_list=[]
    for line in ifile:
        line_list=line.split()
        file_list.append(line_list)
    return file_list
    ifile.close()  #always a good habit to close your file properly
    
def count_line(file_name):
    '''count how many lines of this file'''
    ifile=open(file_name,"r")
    for i,l in enumerate(ifile):#enumerate(sequence [,start=0]) Return
                                    # an enumerate object,sequence must be a 
                                    #sequence, an iterator or some other objec which supports iteration
        pass
    return i+1
    ifile.close()
    
