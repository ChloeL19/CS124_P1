#!/bin/bash
# 1D graphs
echo "---------------beginning 0-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 5 0
./randmst_noheap 0 256 5 0
./randmst_noheap 0 512 5 0
./randmst_noheap 0 1024 5 0
./randmst_noheap 0 2048 5 0
./randmst_noheap 0 4096 5 0
./randmst_noheap 0 8192 5 0
./randmst_noheap 0 16384 5 0
./randmst_noheap 0 32768 5 0
./randmst_noheap 0 65536 5 0
./randmst_noheap 0 131072 5 0
./randmst_noheap 0 262144 5 0
echo "-----------------------------------------------------------------------------"
# #2D graphs
echo "---------------beginning 2-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 5 2
./randmst_noheap 0 256 5 2
./randmst_noheap 0 512 5 2
./randmst_noheap 0 1024 5 2
./randmst_noheap 0 2048 5 2
./randmst_noheap 0 4096 5 2
./randmst_noheap 0 8192 5 2
./randmst_noheap 0 16384 5 2
./randmst_noheap 0 32768 5 2
./randmst_noheap 0 65536 5 2
./randmst_noheap 0 131072 5 2
./randmst_noheap 0 262144 5 2
echo "-----------------------------------------------------------------------------"
# #3D graphs
echo "---------------beginning 3-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 5 3
./randmst_noheap 0 256 5 3
./randmst_noheap 0 512 5 3
./randmst_noheap 0 1024 5 3
./randmst_noheap 0 2048 5 3
./randmst_noheap 0 4096 5 3
./randmst_noheap 0 8192 5 3
./randmst_noheap 0 16384 5 3
./randmst_noheap 0 32768 5 3
./randmst_noheap 0 65536 5 3
./randmst_noheap 0 131072 5 3
./randmst_noheap 0 262144 5 3
echo "-----------------------------------------------------------------------------"
# #4D graphs
echo "---------------beginning 4-dim experiment, n in increasing order-------------"
./randmst_noheap 0 128 5 4
./randmst_noheap 0 256 5 4
./randmst_noheap 0 512 5 4
./randmst_noheap 0 1024 5 4
./randmst_noheap 0 2048 5 4
./randmst_noheap 0 4096 5 4
./randmst_noheap 0 8192 5 4
./randmst_noheap 0 16384 5 4
./randmst_noheap 0 32768 5 4
./randmst_noheap 0 65536 5 4
./randmst_noheap 0 131072 5 4
./randmst_noheap 0 262144 5 4
echo "-----------------------------------------------------------------------------"