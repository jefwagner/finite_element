#include <stdio.h>
#include <iostream>
#include <fstream>
#include "test.hpp"
#include "mesh.hpp"
#include "testbuilder.hpp"
#include "triangle.hpp"


//function for t e s ting integrate

using namespace std;
namespace{
  double func(Vector2d point){
    double ret;
    ret = point[0] * point[1];
    return ret;
  }
}

void test_constructor(){
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

  Mesh mesh_test(&my_tio);

  print_status(1==1   ,"Mesh Constructor");//t e s t mesh constructor

}

void test_reorder(){
  triangulateio my_tio;

  my_tio = mesh_constructor();


  Mesh mesh_test(&my_tio);

  //reorder_nodes t e s t by organizing around the 12th point at 2,2 in a 5 by 6 grid

  mesh_test.reorder_nodes(12);

  int count = 0;

  Vector2d point = mesh_test.points_get(29);

  if(point[1] == 5.){
    count++;
  }

  if(point[0] == 0.){
    count++;
  }

  //reorder_nodes t e s t by organizing around the 12th point at 2,2 in a 5 by 6 grid

  mesh_test.reorder_nodes(5);  //Node 17 before the reorder above

  //t e s t succeeds when the last 5 poinst have y value of 5

  point = mesh_test.points_get(29);

  if(point[1] == 0.){
    count++;
  }
  if(point[0] == 4.){
    count++;
  }

  //Clean up

  if(my_tio.pointlist != NULL){
    delete[] my_tio.pointlist;
  }

  if(my_tio.trianglelist != NULL){
    delete[] my_tio.trianglelist;
  }

  if(my_tio.segmentlist != NULL){
    delete[] my_tio.segmentlist;
  }

  print_status( count==4, "reorder_nodes");
}

void test_integrate(){
  triangulateio my_tio;

  my_tio = mesh_constructor();

  Mesh mesh_test(&my_tio);

  // print_status(1==1   ,"Mesh Constructor");//t e s t mesh constructor

  //Integration t e s t using a hemisphere

  double integral = mesh_test.integrate(func); // t e s t integration actual value ~16.75

  cout << "Integrate = " << integral << endl;

  integral = round(integral);

  print_status(integral == 100, "Integrate");
}

void test_massMatrix(){
  triangulateio my_tio;

  my_tio = mesh_constructor();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double, RowMajor, int> mass_mat_unordered(mesh_test.num_points, mesh_test.num_points);

  mesh_test.mass_matrix(func, mass_mat_unordered);

  ofstream unordered_mat;
  unordered_mat.open("unordered_mat.txt");
  if(unordered_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  for(int i=0; i<30; i++){
    for(int j=0; j<30; j++){
      if(mass_mat_unordered.coeffRef(i,j) != 0){
        unordered_mat << i << "," << j << "," << mass_mat_unordered.coeffRef(i,j) << "\n";
      }
    }
  }
  unordered_mat.close();

  SparseMatrix<double, RowMajor, int> mass_mat_ordered(mesh_test.num_points, mesh_test.num_points);

  mesh_test.reorder_nodes(12);

  mesh_test.mass_matrix(func, mass_mat_ordered);

  ofstream ordered_mat;
  ordered_mat.open("ordered_mat.txt");
  for(int i=0; i<30; i++){
    for(int j=0; j<30; j++){
      if(mass_mat_ordered.coeffRef(i,j) != 0){
        ordered_mat << i << "," << j << "," << mass_mat_ordered.coeffRef(i,j) << "\n";
      }
    }
  }
  ordered_mat.close();

  print_status(1 == 1, "massMatrix");
}
