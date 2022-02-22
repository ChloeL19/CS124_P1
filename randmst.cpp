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
#include <map>

// define useful typedefs
using id_type = uint64_t;
using float_type = double;
using unif_distr = std::uniform_int_distribution<uint64_t>;
using unif_real_distr = std::uniform_real_distribution<double>;

namespace my1 {

    /*
    A class representing a complete graph, using some pseudorandom number generators to
    allow us to obtain the random weights on edges
    */
    class CompleteGraph {
    public:

        // constructor/destructor
        CompleteGraph(id_type num_vertices, int vertex_dim, id_type seed = 17):n(num_vertices), dim(vertex_dim), p(std::numeric_limits<id_type>::max() - 58)
        {

            // define the random generator instance
            std::default_random_engine gen(seed);
            std::uniform_real_distribution<float_type> U(0.0, 1.0);
            unif_distr Uint(1, p);
            inv_p = 1.0 / static_cast<float_type>(p);

            // create the random multiplier and shift
            a = Uint(gen);
            b = Uint(gen);

            // compute points
            points.resize(n);
            for(id_type i = 0; i < n; ++i){
                for(int d = 0; d < dim; ++d){
                    points[i][d] = U(gen);
                }
            }

        }
        ~CompleteGraph() = default;

        // get edge weight
        float_type get_edge_weight(id_type i, id_type j) const {
            if(dim == 0){

                /*
                // test weights to sanity check it gets the right answer. Should
                // return an MST with weight equal to (n - 1), where n is the number of vertices
                if( i > j ){
                    auto tmp = i;
                    i = j;
                    j = tmp;
                }
                if( (j - i) == 1){
                    return 1;
                }else{
                    return 10;
                }
                 */

                if( i > j ){
                    auto tmp = i;
                    i = j;
                    j = tmp;
                }
                auto k = j + i*n;
                auto ri = randint(k);
                return (static_cast<float_type>(ri) * inv_p);

            }else{

                //std::array<float_type, 4> xi = compute_point_location(i);
                //std::array<float_type, 4> xj = compute_point_location(j);
                auto& xi = points[i];
                auto& xj = points[j];
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

        std::array<float_type, 4> compute_point_location(id_type i) const {

            // init empty vector
            std::array<float_type, 4> x{};

            // compute individual component values via the random number generator
            for(int d = 0; d < dim; ++d){
                auto k = compute_dth_component_idx(i, d);
                x[d] = static_cast<float_type>(randint(k)) * inv_p;
            }

            // return the resulting vector
            return x;
        }

        std::array<float_type, 4> subtract_vectors(const std::array<float_type, 4>& xi, const std::array<float_type, 4>& xj) const {
            std::array<float_type, 4> difference;
            for(int d = 0; d < dim; ++d){
                difference[d] = xi[d] - xj[d];
            }
            return difference;
        }

        float_type compute_norm(const std::array<float_type, 4>& x) const {
            float_type inner_prod = 0.0;
            for(int d = 0; d < dim; ++d){
                inner_prod += x[d] * x[d];
            }
            return std::sqrt(inner_prod);
        }

        // our pseudo random number generator
        id_type a, b, p;
        float_type inv_p;

        // set some variables to define the graph
        id_type n;
        int dim;

        // add the points for the vertices
        std::vector<std::array<float_type, 4>> points;

    };

    struct vertex_data {
        id_type vid; // vertex id
        float_type cost;

        vertex_data():vid(0),cost(0){}

        vertex_data(id_type id, float_type c){
            vid = id;
            cost = c;
        }
        bool operator>(const vertex_data& vd) const {
            return cost > vd.cost;
        }
    };

    /*
     * This class is meant to create a min heap for our vertices, where we are given in advance
     * that there is n vertices. The idea here is we will track these n vertices and allow for updating their
     * priorities and dynamically changing the heap structure. The purpose for this is in the context of
     * the Minimum Spanning Tree (MST) implementation where we want to pop off a vertex with minimum distance to
     * the arbitrary start vertex. The fact is, if a vertex gets added to the heap with a new, smaller distance value,
     * then we should not be adding that as a new entry to the heap but just update the priority of this vertex in the
     * heap. This keeps the heap size O(n) so we do not incur a large memory cost associated with the worst case situation
     * of having O(n^2) elements in the heap. Note that if the heap is poly(n) in size, the dynamic operations in the
     * heap are still O(\log(n))
     *
     */
    class vertex_min_heap {
    public:

        // ctor/dtor
        vertex_min_heap(id_type _num_vertices):vertex_heap(_num_vertices),_size(0), pos_map(_num_vertices), num_vertices(_num_vertices) {
            constexpr float_type INFTY = std::numeric_limits<float_type>::max();
            for(id_type i = 0; i < num_vertices; ++i){
                vertex_heap[i].vid = i;
                vertex_heap[i].cost = INFTY;
                pos_map[i] = i;
            }
        }
        ~vertex_min_heap() = default;

        // main methods
        void push(const vertex_data& value) {
            constexpr float_type INFTY = std::numeric_limits<float_type>::max();
            if( vertex_heap[pos_map[value.vid]].cost == INFTY ){
                _size++;
            }
            update(value);
        }
        vertex_data pop() {
            constexpr float_type INFTY = std::numeric_limits<float_type>::max();
            auto value = vertex_heap[0];
            vertex_heap[0].cost = INFTY;
            heapify_down(0);
            _size--;
            return value;
        }
        void update(const vertex_data& value) {
            auto idx = pos_map[value.vid];
            vertex_heap[idx].cost = value.cost;
            heapify_up(idx);
            heapify_down(idx);
        }

        // methods for accessing info
        id_type size() const {
            return _size;
        }

        bool empty() const {
            return (_size == 0);
        }

        void print() const {
            for(auto& v: vertex_heap){
                std::cout << "(" << v.vid << ", " << v.cost << ") ";
            }
            std::cout << std::endl;

            std::cout << "[ ";
            for(auto& idx: pos_map){
                std::cout << idx << " ";
            }
            std::cout << " ]" << std::endl;
        }


    private:
        id_type _size, num_vertices;
        std::vector<vertex_data> vertex_heap;
        std::vector<id_type> pos_map;

        // connectivity methods
        id_type left_child(id_type i) const {
            return 2*i + 1;
        }
        id_type right_child(id_type i) const {
            return 2*i + 2;
        }
        id_type left_parent(id_type i) const {
            return (i/2);
        }
        id_type right_parent(id_type i) const {
            return (i/2) - 1;
        }
        id_type parent(id_type i) const{
            if( i % 2 == 0 ){ // right child
                return right_parent(i);
            }else{ // left child
                return left_parent(i);
            }
        }

        // heapify methods
        void swap(id_type i1, id_type i2) {
            auto v1 = vertex_heap[i1].vid;
            auto v2 = vertex_heap[i2].vid;

            // swap vertex heap values
            auto tmp = vertex_heap[i1];
            vertex_heap[i1] = vertex_heap[i2];
            vertex_heap[i2] = tmp;

            // swap ID map values
            auto tmp2 = pos_map[v1];
            pos_map[v1] = pos_map[v2];
            pos_map[v2] = tmp2;
        }
        void heapify_up(id_type idx) {
            auto pidx = parent(idx);
            while( idx > 0 && (vertex_heap[pidx] > vertex_heap[idx]) ){
                swap(idx, pidx);
                idx = pidx;
                pidx = parent(idx);
            }
        }
        void heapify_down(id_type idx) {
            auto max_index = num_vertices-1;
            auto lc = left_child(idx), rc = right_child(idx);
            while( true ){
                id_type idx_tmp = idx;

                // figure out if we need to swap the index
                if( (lc <= max_index) && (vertex_heap[idx_tmp] > vertex_heap[lc]) ){
                    idx_tmp = lc;
                }
                if( (rc <= max_index) && (vertex_heap[idx_tmp] > vertex_heap[rc]) ){
                    idx_tmp = rc;
                }
                if(idx_tmp != idx){
                    swap(idx, idx_tmp);
                }else{
                    break;
                }

                // update the indices
                idx = idx_tmp;
                lc = left_child(idx);
                rc = right_child(idx);
            }
        }


    };
}// end namespace my1

namespace my2 {
    /*
    A struct representing a vertex in the graph.
    */
    struct vertex_data {
        id_type vid; // vertex id
        int v_dim; // dimension of the vertex
        double priority = -1; // weight of lightest edge connecting this node to T'
        std::vector<double> v; // holds coordinates of vertex itself
        int rank = 0; // for fibonacci heap, number of child nodes
        vertex_data* parent; // fib heap requires doubly linked list
        std::vector<vertex_data*> children; // list (vector) of children nodes
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
        CompleteGraph(int num_vertices, int vertex_dim, int _seed):n_vertices(num_vertices), seed(_seed) {
            // initialize random number generator
            std::default_random_engine gen(seed);
            unif_real_distr U(0.0, 1.0);
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
        double dist(int vid1, int vid2) const {
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

                int i = vid1, j = vid2;
                if( i > j ){
                    auto tmp = i;
                    i = j;
                    j = tmp;
                }
                if( (j - i) == 1){
                    return 1;
                }else{
                    return 10;
                }
                

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

        int get_num_vertices() const {
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
            return tree_list.back();// last element should always be the most 
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
}// end namespace my2


// list of MST algorithms
namespace my1 {
    float_type prims_mst_algorithm(const CompleteGraph& G);
    float_type prims_mst_algorithm2(const CompleteGraph& G);
}
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

                my1::CompleteGraph G(num_points, dimension);

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

namespace my1 {

    using minheap_t = std::priority_queue<vertex_data, std::vector<vertex_data>, std::greater<vertex_data>>;

    // first implementation of Prim's algorithm using STL's priority queue DS
    float_type prims_mst_algorithm(const CompleteGraph& G){

        // get number of vertices in graph
        auto n = G.get_num_vertices();

        // define NIL as the largest value id_type can take on
        constexpr id_type NIL = -1;

        // define infinity as the largest distance type we can have
        constexpr float_type INFTY = std::numeric_limits<float_type>::max();

        // define array of integers and vertices
        std::vector<id_type> prev(n, NIL);
        std::vector<float_type> dist(n, INFTY);

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
        int counter = 0;
        while( not unvisited.empty() ){
            if( counter++ % 1000 == 0 ){
                std::cout << "Num unvisited is " << unvisited.size() << " and H.size = " << H.size() << std::endl;
            }

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
        float_type tree_weight = 0;

        // loop over vertices other than s and look at the `prev` list to obtain
        // the edges we need to compute the weight of and add to our sum
        for(id_type w = 1; w < n; ++w){
            auto weight = G.get_edge_weight(w, prev[w]);
            tree_weight += weight;
        }

        // return the tree weight
        return tree_weight;

    }

    // second implementation of Prim's algorithm using custom memory saving min heap DS
    float_type prims_mst_algorithm2(const CompleteGraph& G){

        // get number of vertices in graph
        auto n = G.get_num_vertices();

        // define NIL as the largest value id_type can take on
        constexpr id_type NIL = -1;

        // define infinity as the largest distance type we can have
        constexpr float_type INFTY = std::numeric_limits<float_type>::max();

        // define array of integers and vertices
        std::vector<id_type> prev(n, NIL);
        std::vector<float_type> dist(n, INFTY);

        // define the initially empty set S
        std::unordered_set<id_type> S;
        std::set<id_type> unvisited;
        for(id_type w = 0; w < n; ++w){
            unvisited.insert(w);
        }

        // define the initial heap H
        vertex_min_heap H(n);

        // choose starting vertex
        id_type s = 0;
        dist[s] = 0;
        H.push(vertex_data(s, 0));

        // main loop
        int counter = 0;
        while( not unvisited.empty() ){
            /*if( counter++ % 1000 == 0 ){
                std::cout << "Num unvisited is " << unvisited.size() << " and H.size = " << H.size() << std::endl;
            }*/

            // pop off the top vertex from the heap along with its data
            auto vdata = H.pop();
            auto v = vdata.vid;

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
        }// end while loop

        // compute the tree weight
        float_type tree_weight = 0;

        // loop over vertices other than s and look at the `prev` list to obtain
        // the edges we need to compute the weight of and add to our sum
        for(id_type w = 1; w < n; ++w){
            auto weight = G.get_edge_weight(w, prev[w]);
            tree_weight += weight;
        }

        // return the tree weight
        return tree_weight;

    }

    float_type run_trial(id_type num_points, int dimension, int seed){
        CompleteGraph G(num_points, dimension, seed);
        return prims_mst_algorithm2(G);
    }
}// end namespace my1

namespace my2 {

    float_type prims_mst_algorithm(CompleteGraph& G){

        // get number of vertices in graph
        auto n = G.get_num_vertices();

        // define NIL as the largest value id_type can take on
        constexpr id_type NIL = -1;

        // define infinity as the largest distance type we can have
        constexpr float_type INFTY = std::numeric_limits<float_type>::max();

        // define array of integers and vertices
        std::vector<id_type> prev(n, NIL);
        std::vector<float_type> dist(n, INFTY);

        // define the initially empty set S
        std::unordered_set<id_type> S;
        std::set<id_type> unvisited;
        for(id_type w = 0; w < n; ++w){
            unvisited.insert(w);
        }

        // define the initial heap H
        fibHeap H;

        // choose starting vertex
        id_type s = 0;
        dist[s] = 0;
        vertex_data& v = G.vertices[s];
        v.priority = dist[s];
        H.insert(&v);

        // main loop
        int counter = 0;
        while( not unvisited.empty() ){

            // pop off the top vertex from the heap along with its data
            vertex_data vdata = H.deleteMin();
            auto v = vdata.vid;

            // remove vertex from the unvisited set
            unvisited.erase(v);

            // loop over edges incident to our vertex v, which is every vertex that is not in S
            for(auto w: unvisited){
                auto weight_vw = G.dist(v, w);
                if( dist[w] > weight_vw ){
                    dist[w] = weight_vw; prev[w] = v;
                    vertex_data& wdata = G.vertices[w];
                    wdata.priority = dist[w];
                    H.insert(&wdata);
                }
            }// end for loop
        }// end while loop

        // compute the tree weight
        float_type tree_weight = 0;

        // loop over vertices other than s and look at the `prev` list to obtain
        // the edges we need to compute the weight of and add to our sum
        for(id_type w = 1; w < n; ++w){
            auto weight = G.dist(w, prev[w]);
            tree_weight += weight;
        }

        // return the tree weight
        return tree_weight;

    }

    float_type run_trial(id_type num_points, int dimension, int seed){
        CompleteGraph G(num_points, dimension, seed);
        return prims_mst_algorithm(G);
    }

}// end namespace my2

void run_experiments(int num_trials, id_type num_points, int dimension) {

    // init some important parameters
    float_type avg_weight = 0.0;
    float_type inv_n = 1.0 / static_cast<float_type>(num_trials);

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
        auto weight = my1::run_trial(num_points, dimension, seed);
        avg_weight += weight * inv_n;
    }

    // print the final answer
    int npoints = static_cast<int>(num_points);
    std::printf("%0.10e %i %i %i", avg_weight, npoints, num_trials, dimension);

}
