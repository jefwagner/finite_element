#ifndef JW_FE_TESTBUILDER
#define JW_FE_TESTBUILDER

#include <Eigen/Core>

using namespace Eigen;

namespace{

  // The goal of this is to fake a triangulateio output.
  //  -The mesh that is being created is in described on pg 54 in the notebook.
  triangulateio mesh_constructor(){
    triangulateio my_tio;
    my_tio.pointlist = NULL;
    my_tio.trianglelist = NULL;
    my_tio.segmentlist = NULL;

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
    my_tio.numberofsegments = 18;
    my_tio.segmentlist = new int[my_tio.numberofsegments * 2];
    my_tio.segmentlist[0] = 0;
    my_tio.segmentlist[1] = 1;
    my_tio.segmentlist[2] = 1;
    my_tio.segmentlist[3] = 2;
    my_tio.segmentlist[4] = 2;
    my_tio.segmentlist[5] = 3;
    my_tio.segmentlist[6] = 3;
    my_tio.segmentlist[7] = 4;
    my_tio.segmentlist[8] = 4;
    my_tio.segmentlist[9] = 9;
    my_tio.segmentlist[10] = 9;
    my_tio.segmentlist[11] = 14;
    my_tio.segmentlist[12] = 14;
    my_tio.segmentlist[13] = 19;
    my_tio.segmentlist[14] = 19;
    my_tio.segmentlist[15] = 24;
    my_tio.segmentlist[16] = 24;
    my_tio.segmentlist[17] = 29;
    my_tio.segmentlist[18] = 29;
    my_tio.segmentlist[19] = 28;
    my_tio.segmentlist[20] = 28;
    my_tio.segmentlist[21] = 27;
    my_tio.segmentlist[22] = 27;
    my_tio.segmentlist[23] = 26;
    my_tio.segmentlist[24] = 26;
    my_tio.segmentlist[25] = 25;
    my_tio.segmentlist[26] = 25;
    my_tio.segmentlist[27] = 20;
    my_tio.segmentlist[28] = 20;
    my_tio.segmentlist[29] = 15;
    my_tio.segmentlist[30] = 15;
    my_tio.segmentlist[31] = 10;
    my_tio.segmentlist[32] = 10;
    my_tio.segmentlist[33] = 5;
    my_tio.segmentlist[34] = 5;
    my_tio.segmentlist[35] = 0;
    printf("Edges filled \n");

    return my_tio;
  }

}

#endif
