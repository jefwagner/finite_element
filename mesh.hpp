#ifndef JW_FE_MESH
#define JW_FE_MESH

#include <vector> // Vectors from the standard library

#include <Eigen/Core> // The Vector2d data type

#include <Eigen/SparseCore> // The SparseMatrix data type

#include "triangle.hpp"

using namespace Eigen;
using namespace std;

// ------------
// Convenience structs
// ------------
// The following two structs allow for a more natural structure
// of the mesh class below
/* A 'tri' struct that holds three integer indices */
/* the indices mark the corners of the triangle in a larger array of points */
struct tri{
	int i,j,k;
	int& operator[](int l){
		if( l == 1 ){
			return j;
		}
		if( l == 2 ){
			return k;
		}
		return i;
	};
};
/* An 'edge' struct that holds two integer indices */
/* the indices mark the end-points of an edge in a larger array of points */
struct edge{
	int i,j;
	int& operator[](int l){
		if( l == 1){
			return j;
		}
		return i;
	};
};

/*
 * Mesh class:
 *
 * Members:
 * - num_points : the number of points in the mesh
 * - points : an array of the points in the mesh
 * - num_tris : the number of tris in the mesh
 * - tris : an array of tris in the mesh
 * - num_edges : the number of edges in the mesh (these are all boundary points)
 * - *edges : an array of edges
 *
 * Methods:
 * - Constructor: create a Mesh from a triangulateio data structure
 * - Destructor: release all the dynamically allocated memory
 * - to_triangulateio: output a triangulateio data structure
 * - reorder_nodes: Reorder the nodes in the mesh data
 * - integrate: Integrade a function over the mesh
 * - mass_matrix: Populate a mass matrix from the mesh
 * - stiffness_matrix: Populate a stiffness matrix from the mesh
 */
class Mesh{
	int num_points;
	Vector2d *points;
	int num_tris;
	tri *tris;
	int num_edges;
	edge *edges;

public:
	Mesh(triangulateio *in);
	~Mesh();

	void to_triangulateio( triangulateio *out);
	void reorder_nodes( int);
	double integrate( double(*func)(Vector2d));
	void mass_matrix( double(*rho)(Vector2d), SparseMatrix<double, RowMajor, int> &);
	void stiffness_matrix( double(*a)(Vector2d), SparseMatrix<double, RowMajor, int> &);

	//methods created for reorder_nodes
	int find_tri(int, int);
	tri find_edge(int, int);
	void update_arr_builder(int, int *);
	void add_without_duplicating(int *, int);
	int find_negative(int *);
	int first_non_negative(int *);
	void update_arrays(int *);
	void swap_tri(int, int);
	void swap_edge(int, int);
	void tri_final_form();
	void edge_final_form();

	//test functions
	void test_fill();
	Vector2d points_get(int);
};

#endif /* JW_FE_MESH */
