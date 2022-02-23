#!/bin/bash
# #2D graphs
echo "---------------beginning 2-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 5 2
./randmst_noheap 1 256 5 2
./randmst_noheap 2 512 5 2
./randmst_noheap 3 1024 5 2
./randmst_noheap 4 2048 5 2
./randmst_noheap 5 4096 5 2
./randmst_noheap 6 8192 5 2
./randmst_noheap 7 16384 5 2
./randmst_noheap 8 32768 5 2
./randmst_noheap 9 65536 5 2
./randmst_noheap 10 131072 5 2
./randmst_noheap 11 262144 5 2
echo "-----------------------------------------------------------------------------"