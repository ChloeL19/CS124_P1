#include <cstring>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <set>
using namespace std;

struct vertex {
    int dim; // the dimension of the vertex
    vector<int> vert; // the vertex values
};

// a graph is a set of vertices
struct graph {
    set<vertex> vertices;
};

/*
Generate a complete graph of n vertices, where the weight
of each edge is a real number chosen uniformly at random on
[0,1].
*/
graph* gen1D (int n) {

}

/* 
Helper function for computing euclidean distance.
*/
float euclideanDist (vertex* v1, vertex* v2) {

}

/*
Generate graph on n vertices, where vertices are points chosen
uniformly at random from inside the unit square.
*/
graph* gen2D (){

}

// TODO: fill in other functions here based on pset description

int main () {

}