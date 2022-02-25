#include <cstdint>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>
#include <map>
#include <cmath>
#include <algorithm>
#include <ctime>

//g++ -o randmst_noheap randmst_noheap.cpp -lm -I.

// define useful typedefs
using id_type = uint64_t;
using unif_distr = std::uniform_real_distribution<double>;

/*
A struct representing a vertex in the graph.
*/
struct vertex_data {
    id_type vid; // vertex id
    int v_dim; // dimension of the vertex
    double priority = -1; // weight of lightest edge connecting this node to T'
    int rank = 0; // for fibonacci heap, number of child nodes
    vertex_data* parent; // fib heap requires doubly linked list
    std::vector<vertex_data*> children; // list (vector) of children nodes
    std::vector<double> v; // holds coordinates of vertex itself
};

/*
A struct representing a complete graph, using some pseudorandom number generators to 
allow us to obtain the random weights on edges without storing them all explicitly
*/
struct CompleteGraph {
    int n_vertices;
    int seed;
    int v_dim;
    std::default_random_engine gen;
    std::vector<vertex_data> vertices;
    /* Graph constructor. Initializes all vertices in the graph.*/
    CompleteGraph(int num_vertices, int vertex_dim, int seed) {
        // initialize random number generator
        CompleteGraph::n_vertices = num_vertices;
        CompleteGraph::seed = seed;
        CompleteGraph::v_dim = vertex_dim;
        gen.seed(seed);
        unif_distr U(0.0, 1.0);
        // initialize the vector of vertices, iterate through until num_vertices
        for (int v_i = 0; v_i != num_vertices; v_i ++){
            vertex_data vertex;
            vertex.vid = v_i;
            vertex.v_dim = vertex_dim; // Q: repetitive info storage?
            // initialize vertex value depending on dimension
            if (vertex_dim > 0){
                for (int coord = 0; coord != vertex_dim; coord++){
                    vertex.v.push_back(U(gen));
                }
            }
            vertices.push_back(vertex);
        }
    }
    ~CompleteGraph() = default;

    /*
        debugging visualizations.
    */
    void printVertices(){
        printf("--------Printing the Vertices-------\n");
        for (int i = 0; i != n_vertices; i++) {
            printf("%ld\n", vertices[i].vid);
        }
    }

    /* 
        Calculates distance between two vertices in the graph. If dimension > 0, 
        calculates the Euclidean distance. If dimension = 0, generates a random number.
    */
    double dist(int vid1, int vid2){
        vertex_data v1 = vertices[vid1];
        vertex_data v2 = vertices[vid2];

        if (v1.v_dim != v2.v_dim){
            printf("Dimensions of vertices should be equal!\n");
            return -1;
        } 

        // enforce range
        if (v1.v_dim < 0 || v1.v_dim > 4){
            printf("Check dimension range. 1D graph has dimension 0.\n");
            printf("Current dimension: %d\n", v1.v_dim);
            return -1;
        }

        if (v1.v_dim == 0) {
            unif_distr U(0.0, 1.0);
            auto n = U(gen);
            return n;
        } else {
            // credit for this distance code in c++: https://www.oreilly.com/library/view/c-cookbook/0596007612/ch11s13.html
            //printf("about to calc dist\n");
            auto first = v1.v.begin();
            auto last = v1.v.end();
            auto first2 = v2.v.begin();

            double ret = 0.0;
            while (first != last) {
                double dist = (*first++) - (*first2++);
                ret += dist * dist;
            }
            return sqrt(ret);
        }
    }

    int get_num_vertices(){
        return n_vertices;
    }
};

// list of MST algorithms
double prims_mst_algorithm(CompleteGraph& G);

/* 
    Parse command line arguments and run experiment.
    Commandline expected to be of form: ./randmst 0 numpoints numtrials dimension
    - optional flag: default value 0, not using it yet
    - numpoints: number of vertices in graph
    - numtrials: number of times to calculate the weight of min spanning tree
    - dimension: the dimension of each vertex (acceptible values are: 0,2,3,4)
*/
int main(int argc, char* argv[]){
    // parse command line arguments
    if (argc < 4) {
        printf("Usage: please check params.\n");
        return -1;
    }
    char* s = argv[1];
    int seed = strtol(s, &s, 10);
    char* n_vertices = argv[2];
    int num_vertices = strtol(n_vertices, &n_vertices, 10);
    char* n_trials = argv[3];
    int num_trials = strtol(n_trials, &n_trials, 10);
    char* dim = argv[4];
    int dimension = strtol(dim, &dim, 10);

    double sum_weight = 0;

    for (int t = 0; t != num_trials; t++){
        CompleteGraph g(num_vertices, dimension, seed);
        double trial = prims_mst_algorithm(g);
        sum_weight += trial;
        seed +=1; // change seed for each trial so we get a different graph
        //printf("Trial %d complete!\n", t);
    }
    printf("%f\n", sum_weight/num_trials);

    return 0;

};

 /*
    v,w: vertices
    dist: array[n] of integer
    prev: array[n] of vertices
    S: set of vertices, initially empty
    H: priority heap of V
    H := {s : 0}
    for v ∈ V do
        dist[v] := ∞, prev[v] :=nil
    rof
    dist[s] := 0
    while H \ne 0/
        v := deletemin(h)
        S := S ∪ {v}
        for (v,w) ∈ E and w ∈ V(G) \ S do
            if dist[w] > length(v,w)
                dist[w] := length(v,w), prev[w] := v, insert(w,dist[w],H)
            fi
        rof
    end while end Prim
*/

// We use an array instead of a heap. We take advantage of our knowledge of
// the complete graph structure to make this more efficient.
double prims_mst_algorithm(CompleteGraph& G){
    double INFTY = std::numeric_limits<double>::max();
    double NIL = -1;
    int n = G.n_vertices;


    double tree_weight = 0; // for tracking MST weight
    double H[n]; // distance from one node to another, our "heapless heap"
    for (int i = 0; i != n; i++){
        H[i] = {INFTY};
    }
    H[0] = NIL;
    int v_i = 0;

    int min_h_ind = 0; // index of vertex w smallest weight
    int curr_mst_vertex = 0; // current index we are investigating
    // debugging
    //G.printVertices();

    while (v_i < n - 1) {
        double min_h = INFTY; // smallest weight that exists to outside MST
        //printf("vi: %d\n", v_i);
        for(int v_j = 0; v_j !=n; v_j ++){
            //printf("investigating vertes: %d\n", v_j);
            if (curr_mst_vertex != v_j && H[v_j] != NIL) {
                //printf("about to calc dist\n");
                //printf("the index ids: %d, %d\n", min_h_ind, v_j);
                double edge_dist = G.dist(curr_mst_vertex, v_j); // this needs to stay constant!!
                //printf("the distance %f\n", edge_dist);
                //printf("vids: %d, %d\n", min_h_ind, v_j);
                //printf("distance: %f\n", edge_dist);
                if (edge_dist < H[v_j]){ //Q: are we sure this works?????? why???
                                        // by not changing aren't we storing wrong dist val??
                    H[v_j] = edge_dist;
                }
                if (H[v_j] < min_h) { //update global minimum vertex
                    min_h = H[v_j];
                    min_h_ind = v_j;
                }
                // printf("-----Printing the Array-----\n");
                // for (int i = 0; i != n; i++){
                //     printf("%f\t", H[i]);
                // }
                // printf("\n");
            }
        }
        // printf("the min_h value being added: %f\n", min_h);
        // printf("the index of this vertex: %d\n", min_h_ind);
        tree_weight += min_h; // needs to be lowest possible
                            // connection back into the tree
        curr_mst_vertex = min_h_ind;
        ++v_i;
        H[min_h_ind] = NIL;
    }
    return tree_weight;
}; 
