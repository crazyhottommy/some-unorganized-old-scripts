fileinput = True
while fileinput == True:
    filename = raw_input('Enter file name:')
    if len(filename) > 0:
        try:
            dnafile = open(filename, 'r')
            fileinput = False
        except:
            print 'File does not exist, please try again'
    else:
        sys.exit()
