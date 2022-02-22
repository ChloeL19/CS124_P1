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

test: unit_tests randmst_chloe

randmst: randmst.cpp # IDK WHAT I'M DOING
	$(CC) $(CFLAGS) -c randmst.cpp -o randmst-debug

randmstfast: randmst.cpp
	$(CC) $(CFLAGSFAST) -c randmst.cpp -o randmst

unit_tests: unit_tests.cpp
	$(CC) $(CFLAGSFAST) -c unit_tests.cpp -o unit_tests

unit_tests: randmst_chloe.cpp
	$(CC) $(CFLAGSFAST) -c randmst_chloe.cpp -o randmst_chloe