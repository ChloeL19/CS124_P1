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
LONG LIVE THE FIBHEAP.
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