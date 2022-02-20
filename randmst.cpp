#include <cstdint>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>


// define useful typedefs
using id_type = uint64_t;
using unif_distr = std::uniform_int_distribution<uint64_t>;

/*
A class representing a complete graph, using some pseudorandom number generators to 
allow us to obtain the random weights on edges without storing them all explicitly
*/
class CompleteGraph {
public:

    // constructor/destructor
    CompleteGraph(id_type num_vertices, int vertex_dim, id_type seed = 17):n(num_vertices), dim(vertex_dim), p(std::numeric_limits<id_type>::max() - 58)
    {

        // define the random generator instance
        std::default_random_engine gen(seed);
        unif_distr U(1, p);

        // create the random multiplier and shift
        a = U(gen);
        b = U(gen);

    }
    ~CompleteGraph() = default;

    // get edge weight
    double get_edge_weight(id_type i, id_type j) const {

    }

    // get number of vertices
    id_type get_num_vertices() const { return n; }

    // get dimension dim
    int get_dimension() const { return dim; }

private:

    // compute a random integer between 0 and p-1
    id_type randint(id_type idx){
        return (a*idx + b) % p;
    }

    // our pseudo random number generator
    id_type a, b, p;

    // set some variables to define the graph
    id_type n;
    int dim;

};

// list of MST algorithms
double prims_mst_algorithm(const CompleteGraph& G);

int main(int argc, char** argv){



    return 0;

}

struct vertex_data {
    id_type vid; // vertex id
    double cost;

    vertex_data(id_type id, double c){
        vid = id;
        cost = c;
    }
    bool operator>(const vertex_data& vd) const {
        return cost > vd.cost;
    }

};
using minheap_t = std::priority_queue<vertex_data, std:vector<vertex_data>, std::greater<vertex_data>>;

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
double prims_mst_algorithm(const CompleteGraph& G){

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
    std::set<id_type> unvisited;
    for(id_type w = 0; w < n; ++w){
        unvisited.insert(w);
    }

    // define the initial heap H
    minheap_t H;

    // choose starting vertex
    id_type s = 0;
    dist[s] = 0;
    H.push(vertex_data(s, 0));
    
    // main loop
    while( not H.empty() ){

        // pop off the top vertex from the heap along with its data
        auto vdata = H.top(); H.pop();

        // add the vertex v to our set S
        S.insert(vdata.vid);

        // remove vertex from the unvisited set
        unvisited.erase(vdata.vid);

        // loop over edges incident to our vertex v, which is every vertex that is not in S
        for(auto w: unvisited){
            auto weight_vw = G.get_edge_weight(v, w);
            if( dist[w] > weight_vw ){
                dist[w] = weight_vw; prev[w] = v;
                H.push(vertex_data(w, dist[w]));
            }
        }// end for loop
        /*for(id_type w = 0; w < n; ++w){
            if( S.count(w) == 0 ){ // if S does not contain w, then
                auto weight_vw = G.get_edge_weight(v, w);
                if( dist[w] > weight_vw ){
                    dist[w] = weight_vw; prev[w] = v;
                    H.push(vertex_data(w, dist[w]));
                }
            }
        }// end for loop*/
    }// end while loop

    // compute the tree weight
    double tree_weight = 0;

    // loop over vertices other than s and look at the `prev` list to obtain
    // the edges we need to compute the weight of and add to our sum
    for(id_type w = 1; w < n; ++w){
        auto weight = G.get_edge_weight(w, prev[w]);
        tree_weight += weight;
    }

    // return the tree weight
    return tree_weight;

}
