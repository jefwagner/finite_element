#include "triangle.hpp"
#include "test.hpp"
#include "mesh.hpp"

void test_constructor(){
  triangulateio my_tio;

  //Fill out node list
  my_tio.numberofpoints = 30;
  my_tio.pointlist = (double *) NULL;
  for(node=0; node<my_tio.numberofpoints; node++){
    my_tio.numberofpoints[node] = ((double)(node%5),(double)(node/5));
  }

  //Fill out triangle list
  my_tio.numberoftriangles = 40;
  my_tio.trianglelist = (int *) NULL;
  for(start=0; start<5; start ++){
    for(i=start; i<25; i+=5){
      tris[6*start + (2*(i/5))].i = i;
      tris[6*start + (2*(i/5))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   + 1].i = i;

      tris[6*start + (2*(i/5))].j = i+1;
      tris[6*start + (2*(i/5)) + 1].j = i+5;

      tris[6*start + (2*(i/5))].k = i+5;
      tris[6*start + (2*(i/5)) + 1].k = i+5;
    }
  }

  //Fill out edge list
  my_tio.numberofedges = 18;
  my_tio.edgelist = (int *) NULL;
  edges[0].i = 0;
  edges[0].j= 1;
  edges[1].i = 1;
  edges[1].j = 2;
  edges[2].i = 2;
  edges[2].j = 3;
  edges[3].i = 3;
  edges[3].j = 4;
  edges[4].i = 4;
  edges[4].j = 9;
  edges[5].i = 9;
  edges[5].j = 14;
  edges[6].i = 14;
  edges[6].j = 19;
  edges[7].i = 19;
  edges[7].j = 24;
  edges[8].i = 24;
  edges[8].j = 29;
  edges[9].i = 29;
  edges[9].j = 28;
  edges[10].i = 28;
  edges[10].j = 27;
  edges[11],i = 27;
  edges[11].j = 26;
  edges[12].i = 26;
  edges[12].j = 25;
  edges[13].i = 25;
  edges[13].j = 20;
  edges[14].i = 20;
  edges[14].j = 15;
  edges[15].i = 15;
  edges[15].j = 10;
  edges[16].i = 10;
  edges[16].j = 5;
  edges[17].i = 5;
  edges[17].j = 0;

  print_status(1==0   ,"Mesh Constructor");
}
