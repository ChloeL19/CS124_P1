#include <cstdint>
#include <vector>
#include <random>
#include <limits>
#include <queue>
#include <unordered_set>
#include <map>
#include <cmath>
#include "randmst_chloe.cpp"

// Run this command to build this file:
// g++ -o unit_tests unit_tests.cpp randmst_chloe.cpp -lm -I.

int printGraph(CompleteGraph g) {
    for (int i = 0; i != g.vertices.size(); i++){
        printf("Vertex %d coordinates: \n", i);
        auto vexes = (g.vertices.find(i)->second).v;
        if (vexes.empty()){
            printf("No coordinate \n");
        } else {
            for (int j = 0; j != vexes.size(); j++){
                printf("%.2f ", vexes.at(j));
            }
            printf("\n");
        }
    }
    return 0;
}

int printHeap(fibHeap H){
    // iterate through tree_list
    printf("Printing the heap now.----------------------------------\n");
    printf("The max rank of the heap: %d\n", H.max_rank);
    for (auto t : H.tree_list){
        // print the number of children for each tree in parentheses
        if (t->children.empty()){
            printf("---Vertex %d (no children)---\t", t->vid);
        } else {
            printf("---Vertex %d (%d children)---\t", t->vid, t->children.size());
        }
    }
    printf("\n");
    return 0;
}

int main() {
    /*
        Testing Graph Construction!
    */
    // construct 1D graph
    printf("Graph 1D-------\n");
    CompleteGraph g1(5,0,0);
    printGraph(g1);
    // construct 2D graph
    printf("Graph 2D-------\n");
    CompleteGraph g2(5,2,2);
    printGraph(g2);
    // construct 3D graph
    printf("Graph 3D-------\n");
    CompleteGraph g3(5,3,3);
    printGraph(g3);
    // construct 4D graph
    printf("Graph 4D-------\n");
    CompleteGraph g4(5,4,4);
    printGraph(g4);
    // test distance for dim=0
    printf("1D distance between v1 and v2 is %f\n", g1.dist(0,1));
    // test distance for dim=2
    printf("2D distance between v1 and v2 is %f\n", g2.dist(0,1));
    // test distance for dim=3
    printf("3D distance between v1 and v2 is %f\n", g3.dist(0,1));
    // test distance for dim=4
    printf("4D distance between v1 and v2 is %f\n", g4.dist(0,1));
    // large-ass graph
    CompleteGraph g5(300000, 4, 0);
    printf("created a freaking huge graph.\n");
    printf("4D distance between v1 and v2 is %f\n", g5.dist(0,1));

    /* Testing Heap Implementation*/
    // the vertices used for testing:
    vertex_data v1; v1.priority = 2; v1.vid=1;// no children --> priority: 2
    printf("created vertex1\n");
    vertex_data v2; v2.priority = 2; v2.vid=2;// two children --> priority: 2
    vertex_data v21; vertex_data v22; v21.priority=7; v22.priority=8;// priorities: 7 & 8
    v2.children.push_back(&v21); v2.children.push_back(&v22);
    printf("created vertex2\n");
    vertex_data v3; v3.priority=1; v3.vid=3;// no children --> priority: 1
    printf("created vertex3\n");
    vertex_data v4; v4.priority=3; v4.vid=4;// three children --> priority: 3
    vertex_data v41; vertex_data v42; vertex_data v43; // --> priorities: 10,11,12
    v41.priority=10; v42.priority=11; v43.priority=12;
    v4.children.push_back(&v41); v4.children.push_back(&v42); v4.children.push_back(&v43);
    printf("created vertex4\n");
    // test insert
    fibHeap H;
    printf("created fibHeap\n");
    H.insert(&v1);
    //printf("Inserted first vertex\n");
    H.insert(&v2);
    H.insert(&v3);
    H.insert(&v4);
    printHeap(H);
    fibHeap H2; // test different pattern of vertices, w diff priorities
    // test merge function
    printf("------Tesing the Merge function-----\n");
    H.merge(&v4, &v1);
    printf("The vertex 1 child id: %d\n", v1.children[0]->vid);
    printf("The vertex 4 child id: %d\n", v4.children[0]->vid);
    printHeap(H);
    // test deletemin
    printf("-----Now testing deleteMin-----\n");
    auto vmin = H.deleteMin(); // should pop v3 and leave rest of heap structure unchanged
    printf("Popped vertex id: %d\n", vmin->vid);
    printf("Remaining heap: \n");
    printHeap(H);
    // another test of deletMmin(), should pop v2 and then merge v1 and v3
    // recreate the set of vertices lol
    printf("Creating a new heap for testing: \n");
    v1.children.clear(); v1.priority = 2; v1.vid=1;v1.rank=0;// no children --> priority: 2
    v2.children.clear(); v2.priority = 1; v2.vid=2;v2.rank=2;// two children --> priority: 2
    v21.children.clear(); v22.children.clear(); v21.priority=7; v22.priority=8;// priorities: 7 & 8
    v21.vid=21; v22.vid=22;
    v2.children.push_back(&v21); v2.children.push_back(&v22);
    v3.children.clear(); v3.priority=2; v3.vid=3; v3.rank=0;// no children --> priority: 1
    v4.children.clear(); v4.priority=3; v4.vid=4; v4.rank=3;// three children --> priority: 3
    v41.children.clear(); v42.children.clear(); v43.children.clear(); // --> priorities: 10,11,12
    v41.priority=10; v42.priority=11; v43.priority=12;
    v41.vid=41; v42.vid=42; v43.vid=43;
    v4.children.push_back(&v41); v4.children.push_back(&v42); v4.children.push_back(&v43);
    H2.insert(&v1);
    //printf("Inserted first vertex\n");
    H2.insert(&v2);
    H2.insert(&v3);
    H2.insert(&v4);
    printHeap(H2);
    auto vmin2 = H2.deleteMin(); // should pop v2 and then merge v1 and v3
    printf("Popped vertex id: %d\n", vmin2->vid);
    printf("Remaining heap: \n");
    printHeap(H2);
    /* Testing Prim's Algorithm */
}
