#!/bin/bash

# 1D graphs
echo "---------------beginning 0-dim experiment, n in increasing order-------------"
./randmst_noheap 0 131072 5 0
./randmst_noheap 0 262144 5 0
echo "-----------------------------------------------------------------------------"
# #2D graphs
echo "---------------beginning 2-dim experiment, n in increasing order-------------"
./randmst_noheap 0 65536 5 2
./randmst_noheap 0 131072 5 2
./randmst_noheap 0 262144 5 2
echo "-----------------------------------------------------------------------------"