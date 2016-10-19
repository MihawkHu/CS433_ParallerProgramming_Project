# compare two input file, judge whether they are same
# in order to judge whether the results got by dijkstra.c 
#   and dijkstra_paraller.c are the same
# this file will be called in ./run_test.sh

import sys
import filecmp

if __name__ != '__main__':
    print("Error1: Not itself.")
    
f1_path = str(sys.argv[1])
f2_path = str(sys.argv[2])

if filecmp.cmp(f1_path, f2_path):
    print("Correct!")
else:
    print("Wrong!")
