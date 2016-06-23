#ifndef JW_FE_MESH
#define JW_FE_MESH

#include <vector> // Vectors from the standard library

#include <Eigen/Core> // The Vector2d data type
#include <Eigen/SparseCore> // The SparseMatrix data type

extern "C" {
	#include <triangle.h> // The triangulateio data type
}

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
	int operator[](int l){
		if( l == 0){
			return i;
		}
		if( l == 1 ){
			return j;
		}
		if( l == 2 ){
			return k;
		}
		 return -1;
	};
}
/* An 'edge' struct that holds two integer indices */
/* the indices mark the end-points of an edge in a larger array of points */
struct edge{
	int i,j;
	int operator[](int l){
		if( l == 0){
			return i;
		}
		if( l == 1){
			return j;
		}
		return -1;
	};
}

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
	double integrate( double(*)(Vector2d));
	void mass_matrix( double(*)(Vector2d), &SparseMatrix);
	void stiffness_matrix( double(*)(Vector2d), &SparseMatrix);

	//methods created for reorder_nodes
	tri Find_tri(int, int);
	tri find_edge(int, int);
	void updateArrBuilder(int, int *);
	void add_without_duplicating(int *, int);
	int Find_negative(int *);
	int FirstNonNegative(int *);
	void updateArrays(int *);
	void swap_tri(int, int);
	void swap_edge(int, int);
	void triFinalForm();
	void edgeFinalForm();

	//test functions
	void testfill();
}

#endif /* JW_FE_MESH */
