# File:     hadoop reducer file
# Function: inverted index, from each file in input folder, output inverted 
#           index
# Implementation: 
#           Get information from mapper, count number of each and file, store 
#           them to a dict, and output the sorted inverted index to output file

import sys

dic = {}
def read_mapper_output(file):
    for line in file:
        yield line.rstrip().split('\t', 1)

data = read_mapper_output(sys.stdin)

for current_word, file_name in data:
    p = dic.get(tuple((current_word, file_name)), 0)
    p = p + 1;
    dic[(current_word, file_name)] = p
    
dictk = sorted(dic.iteritems(), key=lambda d:d[0])

for ll in dictk:
    print "%s%s%s:%d" % (ll[0][0], '\t', ll[0][1], ll[1])
