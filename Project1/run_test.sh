#!/bin/bash
rm dijkstra dijkstra_parallel
mpicc -g -Wall -o dijkstra_parallel ./dijkstra_parallel.c
gcc -g -Wall -o dijkstra ./dijkstra.c

size=(10 100 300 400 500 800 1000 2000 3000 4000 5000)
core_num=(2 5 10)

# test different process number
for((p=0;p<3;p++))
do
    echo -e "\n=================================="
    echo -e "Process number: "$[core_num[p]]
    
    for((i=0;i<11;i++))
    do
        echo -e "--------------------------------"
        echo -e "Test data " ${i} "  Size: " ${size[i]}
        mpiexec -n 2 ./dijkstra_parallel "./test_data/"${size[i]}".txt" "./test_result/dijkstra_parallel/"${size[i]}".txt"
        ./dijkstra "./test_data/"${size[i]}".txt" "./test_result/dijkstra/"${size[i]}".txt"
    done
    
done

# check
echo -e "\n\n====================================="
echo -e "Check begin"
for((i=0;i<11;i++))
do
    python ./compare.py "./test_result/dijkstra/"${size[i]}".txt" "./test_result/dijkstra_parallel/"${size[i]}".txt" 
done
echo -e "Check done"
