#ifndef JW_FE_TESTBUILDER
#define JW_FE_TESTBUILDER

#include <Eigen/Core>

using namespace Eigen;

namespace{

  triangulateio mesh_constructor(){
    triangulateio my_tio;
    my_tio.pointlist = NULL;
    my_tio.trianglelist = NULL;
    my_tio.edgelist = NULL;

    //Fill out node list
    my_tio.numberofpoints = 30;
    my_tio.pointlist = new double[my_tio.numberofpoints * 2];
    for(int node=0; node<my_tio.numberofpoints; node++){
      my_tio.pointlist[2*node + 0] = (double)(node%5);
      my_tio.pointlist[2*node + 1] = (double)(node/5);
    }
    printf("Points filled \n");

    //Fill out triangle list
    my_tio.numberoftriangles = 40;
    my_tio.trianglelist = new int[my_tio.numberoftriangles * 3];
    int start;
    int i;

    for(start=0; start<25; start += 5){
        for(i=0; i<4; i++){
          my_tio.trianglelist[(start/5) * 24 + 6 * i + 0] = start + i + 0;
          my_tio.trianglelist[(start/5) * 24 + 6 * i + 3] = start + i + 0;

          my_tio.trianglelist[(start/5) * 24 + 6 * i + 1] = start + i + 1;
          my_tio.trianglelist[(start/5) * 24 + 6 * i + 4] = start + i + 5;

          my_tio.trianglelist[(start/5) * 24 + 6 * i + 2] = start + i + 6;
          my_tio.trianglelist[(start/5) * 24 + 6 * i + 5] = start + i + 6;
        }
    }
    printf("Triangles filled \n");

    //Fill out edge list
    my_tio.numberofedges = 18;
    my_tio.edgelist = new int[my_tio.numberofedges * 2];
    my_tio.edgelist[0] = 0;
    my_tio.edgelist[1] = 1;
    my_tio.edgelist[2] = 1;
    my_tio.edgelist[3] = 2;
    my_tio.edgelist[4] = 2;
    my_tio.edgelist[5] = 3;
    my_tio.edgelist[6] = 3;
    my_tio.edgelist[7] = 4;
    my_tio.edgelist[8] = 4;
    my_tio.edgelist[9] = 9;
    my_tio.edgelist[10] = 9;
    my_tio.edgelist[11] = 14;
    my_tio.edgelist[12] = 14;
    my_tio.edgelist[13] = 19;
    my_tio.edgelist[14] = 19;
    my_tio.edgelist[15] = 24;
    my_tio.edgelist[16] = 24;
    my_tio.edgelist[17] = 29;
    my_tio.edgelist[18] = 29;
    my_tio.edgelist[19] = 28;
    my_tio.edgelist[20] = 28;
    my_tio.edgelist[21] = 27;
    my_tio.edgelist[22] = 27;
    my_tio.edgelist[23] = 26;
    my_tio.edgelist[24] = 26;
    my_tio.edgelist[25] = 25;
    my_tio.edgelist[26] = 25;
    my_tio.edgelist[27] = 20;
    my_tio.edgelist[28] = 20;
    my_tio.edgelist[29] = 15;
    my_tio.edgelist[30] = 15;
    my_tio.edgelist[31] = 10;
    my_tio.edgelist[32] = 10;
    my_tio.edgelist[33] = 5;
    my_tio.edgelist[34] = 5;
    my_tio.edgelist[35] = 0;
    printf("Edges filled \n");

    return my_tio;
  }

}

#endif
