#include <cstdint>
#include <cmath>
#include <array>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>
#include <set>
#include <iostream>
#include <sstream>
#include <cstdio>

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
        inv_p = 1.0 / static_cast<double>(p);

        // create the random multiplier and shift
        a = U(gen);
        b = U(gen);

    }
    ~CompleteGraph() = default;

    // get edge weight
    double get_edge_weight(id_type i, id_type j) const {
        if( dim == 0 ){
            if( i > j ){
                auto tmp = i;
                i = j;
                j = tmp;
            }
            auto k = j + i*n;
            auto ri = randint(k);
            return (static_cast<double>(ri) * inv_p);
        }else{
            std::array<double, 4> xi = compute_point_location(i);
            std::array<double, 4> xj = compute_point_location(j);
            auto delta = subtract_vectors(xi, xj);
            return compute_norm(delta);
        }
    }

    // get number of vertices
    id_type get_num_vertices() const { return n; }

    // get dimension dim
    int get_dimension() const { return dim; }

private:

    // compute a random integer between 0 and p-1
    id_type randint(id_type idx) const {
        return (a*idx + b) % p;
    }

    id_type compute_dth_component_idx(id_type i, id_type d) const {
        auto dim_cast = static_cast<id_type>(dim);
        return d + dim_cast*i;
    }

    std::array<double, 4> compute_point_location(id_type i) const {

        // init empty vector
        std::array<double, 4> x{};

        // compute individual component values via the random number generator
        for(int d = 0; d < dim; ++d){
            auto k = compute_dth_component_idx(i, d);
            x[d] = static_cast<double>(randint(k)) * inv_p;
        }

        // return the resulting vector
        return x;
    }

    std::array<double, 4> subtract_vectors(const std::array<double, 4>& xi, const std::array<double, 4>& xj) const {
        std::array<double, 4> difference;
        for(int d = 0; d < dim; ++d){
            difference[d] = xi[d] - xj[d];
        }
        return difference;
    }

    double compute_norm(const std::array<double, 4>& x) const {
        double inner_prod = 0.0;
        for(int d = 0; d < dim; ++d){
            inner_prod += x[d] * x[d];
        }
        return std::sqrt(inner_prod);
    }

    // our pseudo random number generator
    id_type a, b, p;
    double inv_p;

    // set some variables to define the graph
    id_type n;
    int dim;

};

// list of MST algorithms
double prims_mst_algorithm(const CompleteGraph& G);
void run_experiments(int num_trials, id_type num_points, int dimension);

int main(int argc, char** argv){

    try {
        if( argc < 5 || argc > 5){
            std::cout << "Wrong number of inputs. Expect 5 inputs." << std::endl;
        }else{

            // extract the command line inputs
            std::stringstream ss;
            std::array<int, 4> inputs{};
            for(int i = 1; i < argc; ++i){
                ss.str(std::string(argv[i]));
                ss >> inputs[i-1];
                ss.clear();
            }

            // make use of the integer inputs
            int flag            = inputs[0];
            id_type num_points  = static_cast<id_type>(inputs[1]);
            int num_trials      = inputs[2];
            int dimension       = inputs[3];


            // if flag is 0, then run the experiments
            if( flag == 0 ){
                run_experiments(num_trials, num_points, dimension);
            }
            
            // otherwise, do some custom test work
            else{

                CompleteGraph G(num_points, dimension);

                for(int i = 0; i < num_points; ++i){
                    for(int j = i+1; j < num_points; ++j){
                        auto w = G.get_edge_weight(i, j);
                        std::cout << "weight(" << i << ", " << j << ") = " << w << std::endl;
                    }
                }// end for i

            }// end if-else

        }// end of valid inputs
    }catch(std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
    }

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
using minheap_t = std::priority_queue<vertex_data, std::vector<vertex_data>, std::greater<vertex_data>>;

double prims_mst_algorithm(const CompleteGraph& G){

    // get number of vertices in graph
    auto n = G.get_num_vertices();

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
        auto v = vdata.vid;

        // add the vertex v to our set S
        S.insert(v);

        // remove vertex from the unvisited set
        unvisited.erase(v);

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

double run_trial(id_type num_points, int dimension, int seed){
    CompleteGraph G(num_points, dimension, seed);
    return prims_mst_algorithm(G);
}

void run_experiments(int num_trials, id_type num_points, int dimension) {

    // init some important parameters
    double avg_weight = 0.0;
    double inv_n = 1.0 / static_cast<double>(num_trials);

    // init random number generator
    std::default_random_engine gen(17);
    std::uniform_int_distribution<int> U(1, 1000000);
    std::set<int> seed_history;

    // compute the number of trials
    for(int i = 0; i < num_trials; ++i){

        // loop until we get a seed we've never seen
        int seed;
        while(true){
            seed = U(gen);
            if( seed_history.contains(seed) == 0 ){
                seed_history.insert(seed);
                break;
            }
        }

        // compute the minimum spanning tree value
        auto weight = run_trial(num_points, dimension, seed);
        avg_weight += weight * inv_n;
    }

    // print the final answer
    int npoints = static_cast<int>(num_points);
    std::printf("%0.10e %i %i %i", avg_weight, npoints, num_trials, dimension);

}
