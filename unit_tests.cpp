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

    /* Testing Heap Implementation*/
    // test insert
    // test deletemin

    /* Testing Prim's Algorithm */
}
