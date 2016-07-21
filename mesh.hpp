#ifndef JW_FE_MESH
#define JW_FE_MESH

#include <vector> // Vectors from the standard library

#include <Eigen/Core> // The Vector2d data type

#include <Eigen/SparseCore> // The SparseMatrix data type

#include "triangle.hpp" // The triangulateio data type

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
public:

	// Start of member data
	//#############################################
	int num_points;
	Vector2d *points;
	int num_tris;
	tri *tris;
	int num_edges;
	edge *edges;

	//----------------------------------------------
	// Added for mass_mat, integrate, and stiff_mat

	// Used for Guassian Integration of a function
	//	within the above member functions
	static Vector2d xieta[4];
	static double weight[4];
	//----------------------------------------------
	// Use either depending on whether v or the
	//	gradient of v is used in the member function
	static double v[3][4];
	static Vector2d grad_v[3];
	//----------------------------------------------

	// Start of member functions
	//##############################################

	//Construction and Destructor
	Mesh(triangulateio *in);
	~Mesh();

	//----------------------------------------------
	// Member functions that were the end goal. The
	//	member functions after this section where
	//	created as a means to acheive these member
	//	functions.
	void to_triangulateio( triangulateio *out);
	void reorder_nodes( int);
	double integrate( double(*func)(Vector2d));
	void mass_matrix( double(*rho)(Vector2d), SparseMatrix<double, RowMajor, int> &, SparseMatrix<double, RowMajor, int> &);
	void stiffness_matrix( double(*a)(Vector2d), SparseMatrix<double, RowMajor, int> &, SparseMatrix<double, RowMajor, int> &);

	//----------------------------------------------
	// Methods created for reorder_nodes

	// In order, functions and their purposes:
	//	-find_tri: Sees if a node is an element of
	//		some given triangle.
	//	-find_edge: Sees if a node is an element of
	//		some given pair of edge nodes.
	//	-update_arr_builder: Builds an array of the
	//		the new order of nodes
	//	-add_without_duplicating: Adds elements to
	//		an array if and only if that element is
	//		not in that array already
	//	-find_negative: Finds the index of the first
	//		negative 1 in an array
	//	-first_non_negative: Finds the index of the
	//		first non negative 1 value in an array
	//	-update_arrays: Updates points[], tris[],
	//		and edges[] based on the update array
	//		for one cycle, See notebook pg 48 and 49
	//		for more info on ordering by cycling.
	//	-swap_tri: Goes through every triangle and
	//		if an element is in that triangle it
	//		replaces it with the replacement node minus
	//		num_points to keep track of which ones have
	//		been changed.
	//	-swap_edge: Similar to swap_tri but instead
	//		of triangles, it now uses edges
	//	-tri_final_form: Takes all the elememnts of
	//		tris and adds num_points to undo the
	//		distinguishing of changed elements.
	//	-edge_final_form;	Similar to tri_final_form
	//		but uses edges instead of tris.

	int find_tri(int, int);
	int find_edge(int, int);
	void update_arr_builder(int, int *);
	void add_without_duplicating(int *, int);
	int find_negative(int *);
	int first_non_negative(int *);
	void update_arrays(int *);
	void swap_tri(int, int);
	void swap_edge(int, int);
	void tri_final_form();
	void edge_final_form();

	//----------------------------------------------
	// Methods created for stiffness_matrix and
	// 	mass_matrix.
	int find_not_bound(int , int *);
	int find_bound(int, int *);

	//----------------------------------------------
	// Methods created for stiffness_matrix:
	//	-jacobian: Creates a 2 by 2 matrix that is
	//		is equal to the Jacobian.
	Matrix2d jacobian(Vector2d, Vector2d, Vector2d);

	//----------------------------------------------
	// Test functions:
	//	-test_fill: Fills a simple 5 by 6 grid
	//		described on pg 54 in the notebook.
	//	-points: A getter function for points.
	//		note: This was created before points where
	//			made public, This function is now pointless
	void test_fill();
	Vector2d points_get(int);
};

#endif /* JW_FE_MESH */
