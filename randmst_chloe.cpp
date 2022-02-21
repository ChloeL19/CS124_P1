#include <cstdint>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>
#include <map>
#include <cmath>
#include <algorithm>


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

    /* Calculates distance between two vertices in the graph. If dimension > 0, 
    calculates the Euclidean distance. If dimension = 0, generates a random number.*/
    double dist(int vid1, int vid2){
        vertex_data v1 = vertices.find(vid1)->second;
        vertex_data v2 = vertices.find(vid2)->second;

        if (v1.v_dim != v2.v_dim){
            printf("Dimensions of vertices should be equal!\n");
            return -1;
        } 

        // enforce range
        if (v1.v_dim < 0 ||v1.v_dim == 1||v1.v_dim > 4){
            printf("Check dimension range.\n");
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
        return tree_list.back(); // last element should always be the most 
                                // recently added element
    }

    /*
        Merge the tree in the second argument to the tree in the first argument.
        Tree in fist argument should be modified in-place.
        Assume tree in first argument is later in tree_list than the second
        argument.
        Pointer to t1 should now point to the larger modified tree.
        Returns -1 on error.
    */
    void merge(vertex_data* t1, vertex_data* t2){
        tree_list.erase(std::find(tree_list.begin(), tree_list.end(), t2));
        // Q: is this truly "in-place"???
        if (t1->priority <= t2->priority){
            t1->children.push_back(t2);
            t1->rank++;
            if (t1->rank > max_rank){
                max_rank = t1->rank;
            }
            // delete the second pointer from tree_list
            //tree_list.erase(std::find(tree_list.begin(), tree_list.end(), t2));
        } else {
            t2->children.push_back(t1);
            t2->rank++;
            if (t2->rank > max_rank){
                max_rank = t2->rank;
            }
            //*t1 = *t2; // cannot do this if children is vector of pointers
            // we can assume t1 will always be the current pointer
            // so t2 will always be earlier in tree_list
            //tree_list.erase(std::find(tree_list.begin(), tree_list.end(), t2));
            auto location = std::find(tree_list.begin(), tree_list.end(), t1);
            tree_list.erase(location);
            tree_list.insert(location, t2);
        }
    }

    /* 
        Delete and return the minimum element of 
        the fibonacci heap.
    */
    vertex_data deleteMin(){
       // Q: how do we break ties again??????
       min_node = tree_list[0];
       printf("just identified min node\n");
       // merge children trees into the rest of the tree list
       for (auto child_ptr : fibHeap::min_node->children) {
           fibHeap::insert(child_ptr);
           // update pointer to minimum node if necessary
           if (child_ptr->priority < min_node->priority){
               min_node = child_ptr;
           }
       }
       printf("inserted all children into heap\n");
       // merge the nodes such that no node has the same rank
       // update max_rank as necessary
       // should we use a map for the ranks??? idk how to do
       std::map<int, vertex_data*> rank_tracker;
       printf("created rank tracker\n");
       // NOTE: i may need to make a copy of the tree_list for iterating
       for (auto root : tree_list){
           // while this root has same rank as another
                // merge with that other
            auto insert_ret = rank_tracker.
                                insert(std::pair<int,vertex_data*>
                                (root->children.size(), root));
            printf("inserted first time into rank tracker\n");
            while(!insert_ret.second){
                // merge duplicate's tree with the current root, change in-place
                vertex_data* duplicate = (*insert_ret.first).second;
                printf("found duplicate rank tree pointer\n");
                merge(root, duplicate);
                printf("merged my thing\n");
                // delete the duplicate from the rank_tracker
                rank_tracker.erase(insert_ret.first);
                printf("erased from rank_tracker\n");
                // delete the duplicate from the tree_list
                // tree_list.erase(std::find(tree_list.begin(), tree_list.end(), duplicate));
                // try inserting into the rank tracker again
                auto insert_ret = rank_tracker.
                                insert(std::pair<int,vertex_data*>
                                (root->children.size(), root));
                // FIXME: claiming there is a duplicate rank when there should
                // not be one
            } 
       }

       // Q: free the rank tracker memory??

       return *min_node;
   }
};

/*
// list of MST algorithms
double prims_mst_algorithm(const CompleteGraph& G);
*/

/* 
    Parse command line arguments and run experiment.
    Commandline expected to be of form: ./randmst 0 numpoints numtrials dimension
    - optional flag: default value 0, not using it yet
    - numpoints: number of vertices in graph
    - numtrials: number of times to calculate the weight of min spanning tree
    - dimension: the dimension of each vertex (acceptible values are: 0,2,3,4)
*/

/*
int main(int argc, char** argv){
    // TODO
    // parse command line arguments
    if (argc < 4) {
        printf("Usage: please check params.\n");
        return -1;
    }

    return 0;

}

using minheap_t = std::priority_queue<vertex_data, std::vector<vertex_data>, std::greater<vertex_data>>;
// FIXME: we need to implement priority_queue from scratch


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
/*
double prims_mst_algorithm(CompleteGraph& G){

    // get number of vertices in graph
    int n = G.get_num_vertices();

    // define NIL as the largest value id_type can take on
    constexpr id_type NIL = -1;

    // define infinity as the largest distance type we can have
    constexpr double INFTY = std::numeric_limits<double>::max();

    // define array of integers and vertices
    std::vector<id_type> prev(n, NIL);
    std::vector<double> dist(n, INFTY);

    // define the initially empty set S
    std::unordered_set<id_type> S;
    std::unordered_set<id_type> unvisited; 
    for(id_type w = 0; w < n; ++w){ // insert all vertices as "unvisited"
        unvisited.insert(w);
    }

    // define the initial heap H
    minheap_t H;

    // choose starting vertex // Q: how to make this decision? does it matter?
    id_type s = 0;
    dist[s] = 0;
    H.push(G.vertices.find(0)->second); 
    
    // main loop
    while( ! H.empty() ){

        // pop off the top vertex from the heap along with its data
        auto vdata = H.top(); H.pop();

        // add the vertex v to our set S
        S.insert(vdata.vid);

        // remove vertex from the unvisited set
        unvisited.erase(vdata.vid);

        // loop over edges incident to our vertex v, which is every vertex that is not in S
        for(auto w: unvisited){
            auto weight_vw = G.dist(vdata.vid, w);
            if( dist[w] > weight_vw ){
                dist[w] = weight_vw; prev[w] = vdata.vid;
                H.push(G.vertices.find(w)->second);
            }
        }// end for loop
    }// end while loop

    // compute the tree weight
    double tree_weight = 0;

    // loop over vertices other than s and look at the `prev` list to obtain
    // the edges we need to compute the weight of and add to our sum
    for(id_type w = 1; w < n; ++w){
        auto weight = G.dist(w, prev[w]); //Q: is there a way to store data so we don't 
                                        // recompute distance every time?
        tree_weight += weight;
    }

    // return the tree weight
    return tree_weight;

}; */
