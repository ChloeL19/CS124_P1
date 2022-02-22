#include <cstdint>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>
#include <map>
#include <cmath>
#include <algorithm>

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
    std::map<int, vertex_data> vertices;
    /* Graph constructor. Initializes all vertices in the graph.*/
    CompleteGraph(int num_vertices, int vertex_dim, int seed) {
        // initialize random number generator
        CompleteGraph::n_vertices = num_vertices;
        CompleteGraph::seed = seed;
        std::default_random_engine gen(seed);
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
            // if vertex_dim is zero leave the vertex coordinates as empty vector
            vertices.insert(std::pair<int, vertex_data> (v_i, vertex));
        }
    }
    ~CompleteGraph() = default;

    /* 
        Calculates distance between two vertices in the graph. If dimension > 0, 
        calculates the Euclidean distance. If dimension = 0, generates a random number.
    */
    // FIXME: only square root when adding to MST?? --> eh maybe not
    double dist(int vid1, int vid2){
        vertex_data v1 = vertices.find(vid1)->second;
        vertex_data v2 = vertices.find(vid2)->second;

        if (v1.v_dim != v2.v_dim){
            printf("Dimensions of vertices should be equal!\n");
            return -1;
        } 

        // enforce range
        if (v1.v_dim < 0 ||v1.v_dim == 1||v1.v_dim > 4){
            printf("Check dimension range. 1D graph has dimension 0.\n");
            printf("Current dimension: %d\n", v1.v_dim);
            return -1;
        }

        if (v1.v_dim == 0) {
            std::default_random_engine gen(seed);
            unif_distr U(0.0, 1.0);
            return U(gen);
        } else {
            auto it1 = v1.v.begin();
            auto it2 = v2.v.begin();
            double sumsq = 0;
            // calculate sum of squares
            for (; it1 != v1.v.end(); it1++, it2++){
                sumsq += std::pow((float)(*it1 - *it2),2.0);
            }
            // return the square root of the sum of squares
            return std::pow(sumsq, 0.5);
        }
    }

    int get_num_vertices(){
        return n_vertices;
    }
};

/*
    A fibonacci heap. This will serve as our priority queue in Prim's 
    algorithm.
    UPDATE: ultimately unused, tear.
*/
struct fibHeap {
    std::vector<vertex_data*> tree_list; // Q: or should this be a vector of pointers??
    vertex_data* min_node = NULL;
    int max_rank = 0; // the max rank of our heap

    /* 
        Insert a node into the fibonacci tree list. 
        Return a pointer to the newly inserted node.
    */
    vertex_data* insert(vertex_data* v){
        if (v->priority < 0){
            printf("Insertion failed!\n");
            return NULL;
        }
        tree_list.push_back(v);
        // update the minimum node pointer if necessary
        if (tree_list.size() == 1) {
            min_node = tree_list[0];
        }
        if (v->priority < min_node->priority){
            min_node = v;
        }
        if (v->rank > max_rank){ // update max_rank if necessary
            max_rank = v->rank;
        }
        return tree_list.back(); // last element should always be the most 
                                // recently added element
    }

    /*
        Merge the tree in the second argument to the tree in the first argument.
        Tree in fist argument should be modified in-place.
        Assume tree in first argument is later in tree_list than the second
        argument.
        Pointer to t1 should now point to the larger modified tree.
        Returns an iterator to tree_list where new tree exists.
    */
    std::vector<vertex_data*>::iterator merge(vertex_data* t1, vertex_data* t2){
        tree_list.erase(std::find(tree_list.begin(), tree_list.end(), t2));
        if (t1->priority <= t2->priority){
            t1->children.push_back(t2);
            t1->rank++;
            if (t1->rank > max_rank){
                max_rank = t1->rank;
            }
            return std::find(tree_list.begin(), tree_list.end(), t1);
        } else {
            t2->children.push_back(t1);
            t2->rank++;
            if (t2->rank > max_rank){
                max_rank = t2->rank;
            }
            //*t1 = *t2; // cannot do this if children is vector of pointers
            // we can assume t1 will always be the current pointer
            // so t2 will always be earlier in tree_list
            auto location = std::find(tree_list.begin(), tree_list.end(), t1);
            tree_list.erase(location);
            return tree_list.insert(location, t2);
        }
    }

    /*
        Printing the rank_tracker for debugging purposes
    */
    void printRankTracker(std::map<int, vertex_data*>& rt){
        printf("-------Printing the Rank Tracker--------\n");
        for (auto m : rt){
            printf("---Rank: %d, Vertex: %ld\n", m.first, m.second->vid);
        }
    };

    /* 
        Delete and return the minimum element of 
        the fibonacci heap.
    */
    vertex_data* deleteMin(){
       // erase the minimum node from tree_list
       auto min_location = std::find(tree_list.begin(), tree_list.end(), min_node);
       tree_list.erase(min_location);
       // merge children trees into the rest of the tree list
       for (auto child_ptr : fibHeap::min_node->children) {
           fibHeap::insert(child_ptr);
           // update pointer to minimum node
           if (child_ptr->priority < min_node->priority){
               min_node = child_ptr;
           }
       }
       std::map<int, vertex_data*> rank_tracker;
       for (auto rit = tree_list.begin(); rit != tree_list.end(); rit++){
           // while this root has same rank as another
                // merge with that other
            vertex_data* root = *rit;
            auto insert_ret = rank_tracker.
                                insert(std::pair<int,vertex_data*>
                                (root->children.size(), root));
            while(!insert_ret.second){
                // merge duplicate's tree with the current root, change in-place
                // find duplicate
                vertex_data* duplicate = (*insert_ret.first).second;
                // merge tree, want it to end up at current iterator position
                rit = merge(root, duplicate);
                // delete the duplicate from the rank_tracker
                rank_tracker.erase(insert_ret.first);
                insert_ret = rank_tracker.
                                insert(std::pair<int,vertex_data*>
                                (root->children.size(), root));
            } 
       }
       return min_node;
   }

   /*
        Print out the heap. Primarily for debugging purposes.
   */
    void display(){
        // iterate through tree_list
        printf("Printing the heap now.----------------------------------\n");
        printf("The max rank of the heap: %d\n", max_rank);
        for (auto t : tree_list){
            // print the number of children for each tree in parentheses
            if (t->children.empty()){
                printf("---Vertex %ld (no children)---\t", t->vid);
            } else {
                printf("---Vertex %ld (%ld children)---\t", t->vid, t->children.size());
            }
        }
        printf("\n");
    }

    /*
        Returns true if heap is empty and false otherwise.
    */
   bool empty(){
       return tree_list.empty();
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

    printf("Starting experiment . . . \n");
    for (int t = 0; t != num_trials; t++){
        CompleteGraph g(num_vertices, dimension, seed);
        double trial = prims_mst_algorithm(g);
        sum_weight += trial;
        seed +=1; // change seed for each trial so we get a different graph
        printf("Trial %d complete!\n", t);
    }
    printf("Avg trial weight: %f\n", sum_weight/num_trials);

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
    H[0] = 0;
    int min_h_ind = 0;
    int v_i = 0;
    while (v_i < n) {
        double min_h = INFTY; // global minimum of H
        for(int v_j = 0; v_j !=n; v_j ++){
            if (min_h_ind != v_j && H[v_j] != NIL) {
                double edge_dist = G.dist(min_h_ind, v_j);
                if (edge_dist < H[v_j]){ //Q: are we sure this works?????? why???
                                        // by not changing aren't we storing wrong dist val??
                    H[v_j] = edge_dist;
                }
                if (H[v_j] < min_h) { //update global minimum vertex
                    min_h = H[v_j];
                    min_h_ind = v_j;
                }
            }
        }
        tree_weight += min_h;
        v_i++;
        H[min_h_ind] = NIL;
    }
    return tree_weight;
}; 
