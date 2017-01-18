# File:     hadoop mapper file
# Function: inverted index, from each file in input folder, output inverted 
#           index
# Implementation: 
#           Get information from input folder, replace punctuations with 
#           space, send word and file name information to reducer

import sys
import os

def getFileName():
    if 'map_input_file' in os.environ:
        return os.environ['map_input_file']
    else:
        return 'none'

for line in sys.stdin:
        words = line.strip().lstrip().rstrip().replace("\"", "").\
            replace("*", "").replace("#", "").replace("$", "").\
            replace("=", "").replace("(", "").replace(")", "").\
            replace("-", "").replace("!", "").replace(":", "").\
            replace("\'", "").replace("?", "").replace("/", "").\
            replace("&", "").replace(".", "").replace(",", "").\
            replace("`", "").replace("\x1a", "").split()
        
        for word in words:
            print "%s%s%s" % (word, '\t', getFileName())