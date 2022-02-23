#!/bin/bash
# 1D graphs
echo "---------------Reporting results from single trial for all experiments--------------"
echo "---------------beginning 0-dim experiment, n in increasing order-------------"
./randmst_noheap 0 131072 1 0
./randmst_noheap 0 262144 1 0
echo "-----------------------------------------------------------------------------"
# #2D graphs
echo "---------------beginning 2-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 1 2
./randmst_noheap 0 256 1 2
./randmst_noheap 0 512 1 2
./randmst_noheap 0 1024 1 2
./randmst_noheap 0 2048 1 2
./randmst_noheap 0 4096 1 2
./randmst_noheap 0 8192 1 2
./randmst_noheap 0 16384 1 2
./randmst_noheap 0 32768 1 2
./randmst_noheap 0 65536 1 2
./randmst_noheap 0 131072 1 2
./randmst_noheap 0 262144 1 2
echo "-----------------------------------------------------------------------------"
# #3D graphs
echo "---------------beginning 3-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 1 3
./randmst_noheap 0 256 1 3
./randmst_noheap 0 512 1 3
./randmst_noheap 0 1024 1 3
./randmst_noheap 0 2048 1 3
./randmst_noheap 0 4096 1 3
./randmst_noheap 0 8192 1 3
./randmst_noheap 0 16384 1 3
./randmst_noheap 0 32768 1 3
./randmst_noheap 0 65536 1 3
./randmst_noheap 0 131072 1 3
./randmst_noheap 0 262144 1 3
echo "-----------------------------------------------------------------------------"
# #4D graphs
echo "---------------beginning 4-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 1 4
./randmst_noheap 0 256 1 4
./randmst_noheap 0 512 1 4
./randmst_noheap 0 1024 1 4
./randmst_noheap 0 2048 1 4
./randmst_noheap 0 4096 1 4
./randmst_noheap 0 8192 1 4
./randmst_noheap 0 16384 1 4
./randmst_noheap 0 32768 1 4
./randmst_noheap 0 65536 1 4
./randmst_noheap 0 131072 1 4
./randmst_noheap 0 262144 1 4
echo "-----------------------------------------------------------------------------"