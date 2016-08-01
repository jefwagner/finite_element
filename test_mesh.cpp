#include <stdio.h>
#include <iostream>
#include <fstream>
#include "test.hpp"
#include "mesh.hpp"
#include "region.hpp"
#include "testbuilder.hpp"
#include "triangle.hpp"


//function for t e s ting integrate

using namespace std;
namespace{
  double func(Vector2d point){
    double ret;
    ret = (point[0] * point[1]);
    return ret;
  }

  double func_k(Vector2d point){
    return 1;
  }
}

void test_constructor(){

  // The Goal is to fake a triangulateio output
  //  -The mesh that this describes is on pg 54 in the notebook.
  //  -A copy of this is in testbuilder.hpp.
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

  Mesh mesh_test(&my_tio);

  print_status(1==1   ,"Mesh Constructor");//t e s t mesh constructor

}

void test_reorder(){

  // Call the constructor to construct the t e s t mesh.
  triangulateio my_tio;

  my_tio = mesh_constructor();


  Mesh mesh_test(&my_tio);

  //reorder_nodes t e s t by organizing around the 12th point at 2,2 in a 5 by 6 grid

  mesh_test.reorder_nodes(12);

  int count = 0;

  Vector2d point = mesh_test.points_get(29);

  // T e s t if the last point listed is the farthest away.
  if(point[1] == 5.){
    count++;
  }

  if(point[0] == 0.){
    count++;
  }

  //reorder_nodes t e s t by organizing around the 5th point at 2,3 in a 5 by 6 grid

  mesh_test.reorder_nodes(5);  //Node 17 before the reorder above

  point = mesh_test.points_get(29);

  // T e s t if the last point listed is the farthes away.
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

  // Print success only if both reorders pass above.
  print_status( count==4, "reorder_nodes");
}

void test_integrate(){
  triangulateio my_tio;

  my_tio = mesh_constructor();

  Mesh mesh_test(&my_tio);

  //Integration t e s t using a hemisphere

  double integral = mesh_test.integrate(func); // t e s t integration actual value ~16.75

  integral = round(integral);

  print_status(integral == 100, "Integrate");
}

void test_massMatrix(){
  triangulateio my_tio;

  my_tio = mesh_constructor();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double> not_bound_mat_unordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_unordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.mass_matrix(func, bound_mat_unordered, not_bound_mat_unordered);

  // Nodes on the boundary

  ofstream unordered_bound_mat;
  unordered_bound_mat.open("unordered_bound_mat.txt");
  if(unordered_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  unordered_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<bound_mat_unordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(bound_mat_unordered,k); it; ++it){
      unordered_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  unordered_bound_mat.close();

  // Nodes not on the boundary

  ofstream unordered_not_bound_mat;
  unordered_not_bound_mat.open("unordered_not_bound_mat.txt");
  if(unordered_not_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  unordered_not_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<not_bound_mat_unordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(not_bound_mat_unordered,k); it; ++it){
      unordered_not_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  unordered_not_bound_mat.close();

  // Ordered Nodes

  mesh_test.bound = NULL;
  mesh_test.not_bound == NULL;

  SparseMatrix<double> not_bound_mat_ordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_ordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.reorder_nodes(12);

  mesh_test.mass_matrix(func, bound_mat_ordered, not_bound_mat_ordered);

  // Nodes on the boundary

  ofstream ordered_bound_mat;
  ordered_bound_mat.open("ordered_bound_mat.txt");
  if(ordered_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  ordered_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<bound_mat_ordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(bound_mat_ordered,k); it; ++it){
      ordered_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  ordered_bound_mat.close();

  // Nodes not on the boundary

  ofstream ordered_not_bound_mat;
  ordered_not_bound_mat.open("ordered_not_bound_mat.txt");
  if(ordered_not_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  ordered_not_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<not_bound_mat_ordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(not_bound_mat_ordered,k); it; ++it){
      ordered_not_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  ordered_not_bound_mat.close();

  print_status(1 == 1, "massMatrix");
}

void test_stiffnessMatrix(){
  triangulateio my_tio;

  my_tio = mesh_constructor();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double> not_bound_mat_unordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_unordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.stiffness_matrix(func, bound_mat_unordered, not_bound_mat_unordered);

  // Nodes on the boundary

  ofstream unordered_bound_mat;
  unordered_bound_mat.open("unordered_bound_mat.txt");
  if(unordered_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  unordered_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<bound_mat_unordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(bound_mat_unordered,k); it; ++it){
      unordered_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  unordered_bound_mat.close();

  // Nodes not on the boundary

  ofstream unordered_not_bound_mat;
  unordered_not_bound_mat.open("unordered_not_bound_mat.txt");
  if(unordered_not_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  unordered_not_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<not_bound_mat_unordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(not_bound_mat_unordered,k); it; ++it){
      unordered_not_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  unordered_not_bound_mat.close();

  mesh_test.bound = NULL;
  mesh_test.not_bound = NULL;

  SparseMatrix<double> not_bound_mat_ordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_ordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.reorder_nodes(12);

  mesh_test.stiffness_matrix(func, bound_mat_ordered, not_bound_mat_ordered);

  // Nodes on the boundary

  ofstream ordered_bound_mat;
  ordered_bound_mat.open("ordered_bound_mat.txt");
  if(ordered_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  ordered_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<bound_mat_ordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(bound_mat_ordered,k); it; ++it){
      ordered_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  ordered_bound_mat.close();

  // Nodes not on the boundary

  ofstream ordered_not_bound_mat;
  ordered_not_bound_mat.open("ordered_not_bound_mat.txt");
  if(ordered_not_bound_mat.is_open() == false){
    cout << "Unable to open file" << endl << std::flush;
  }
  ordered_not_bound_mat << mesh_test.num_points<< "\n";
  for(int k=0; k<not_bound_mat_ordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(not_bound_mat_ordered,k); it; ++it){
      ordered_not_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
    }
  }
  ordered_not_bound_mat.close();

  print_status(1 ==1, "stiffnessMatrix");
}

void test_k_matrix(){
  double testk;
  int testcount = 0;
  triangulateio my_tio;

  my_tio = mesh_constructor_simple();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double> not_bound_mat_unordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_unordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.stiffness_matrix(func_k, bound_mat_unordered, not_bound_mat_unordered);

  for(int k=0; k<bound_mat_unordered.outerSize(); ++k){
    for(SparseMatrix<double>::InnerIterator it(bound_mat_unordered,k); it; ++it){
      testk = it.value()/not_bound_mat_unordered.coeffRef(0,0);
      if(testk == -0.25){
        testcount++;
      }
    }
  }
  print_status(testcount == 4, "k_matrix");
}

void test_b_vector(){
  triangulateio my_tio;

  my_tio = mesh_constructor_simple();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double> not_bound_mat_unordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_unordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.stiffness_matrix(func_k, bound_mat_unordered, not_bound_mat_unordered);

  VectorXd w_k, w_ij, w;
  w_k = w_k_builder(func_k, mesh_test);
  VectorXd b = b_vector_builder(bound_mat_unordered, w_k);
}

void test_w_solver(){
  triangulateio my_tio;

  my_tio = mesh_constructor_simple();

  Mesh mesh_test(&my_tio);

  SparseMatrix<double> not_bound_mat_unordered(mesh_test.num_points-mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);
  SparseMatrix<double> bound_mat_unordered(mesh_test.num_edges, mesh_test.num_points-mesh_test.num_edges);

  mesh_test.stiffness_matrix(func_k, bound_mat_unordered, not_bound_mat_unordered);

  VectorXd w_k, w_ij, w;
  w_k = w_k_builder(func_k, mesh_test);
  VectorXd b = b_vector_builder(bound_mat_unordered, w_k);
  w_ij = matrix_solver(not_bound_mat_unordered, b);
  w = w_stitcher(w_k, w_ij, mesh_test);

  print_status(w[4] == 1, "w_solver");
}
