# Makefile for randmst assignment
# inspiration from here: v\https://www.softwaretestinghelp.com/cpp-makefile-tutorial/
# Also from : https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

# ****************************************
# Variables to control Makefile operation
CC = g++
CFLAGS = -g -Wall -O0
CFLAGSFAST = -O2 -Wall

# ****************************************
all: randmst randmstfast

randmst: randmst.cpp # IDK WHAT I'M DOING
	$(CC) $(CFLAGS) -c randmst.cpp -o randmst-debug

randmstfast: randmst.cpp
	$(CC) $(CFLAGSFAST) -c randmst.cpp -o randmst
