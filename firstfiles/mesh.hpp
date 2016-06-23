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

//Constructor => fills data in Mesh, assuming triangle has already run

Mesh::Mesh(triangulatio *in){
	int i;

	//setting the values of the data to values from triangle

	num_points = in->numberofpoints;
	num_tris = in->numberoftriangles;
	num_edges = in->numberofedges

	//creating arrays

	points = new Vector2d[num_points];
	tris = new tri[num_tris];
	edges = new edge[num_edges];

	//Now to fill them
	//-----------------------------
	//First is the point list, assuming 2 doubles labeled x and y

	for(i=0; i<(2*num_points); i++){
		if(i%2 == 0){
			x = in->pointlist[i];
		} else{
			y = in->pointlist[i];
			points[(i-1)/2] = (x,y);
		}
	}

	//-------------------------------------
	//Now the triangle index, assuming 3 ints labeled, i,j,k

	for(i=0; i<(3*num_tris; i++)){
		if(i%3 == 0){
			tris[i/3].i = in->trianglelist[i];
		} else if(i%3 == 1){
			tris[(i-1)/3].j = in->trianglelist[i];
		} else {
			tris[(i-2)/3].k = in->trianglelist[i];
		}
	}

	//--------------------------------------
	//Now the edge index, assuming 3 ints labeled, i,j

	for(i=0; i<(2*num_edges); i++){
		if(i%2 == 0){
			edges[i/2] = in->edgelist[i];
		} else{
			edges[(i-1)/2] = in->edgelist[i]
		}
	}
}

//Deconstructor => frees all arrays created in Constructor

Mesh::~Mesh(){
	if(points != NULL){
		free(points);
	}

	if(tris != NULL){
		free(tris);
	}

	if(edges != NULL){
		free(edges);
	}
}

//to_triangulateio => creates a triangulateio to refine the mesh

Mesh::to_triangulateio(triangulateio *out){
	int i;

	//First to initialize the arays and variables

	out->numberofpoints = 0;
	out->numberofpointattributes = 0;
	out->pointlist = (double *) NULL;
	out->pointmarkerlist = (int *) NULL;
	out->pointattributelist = (double *) NULL;

	out->numberofsegments = 0;
	out->segmentlist = (int *) NULL;
	out->segmentmarkerlist = (int *) NULL;

	out->numberoftriangles = 0;
	out->numberofcorners = 0;
	out->numberoftriangleattributes = 0;

	out->trianglelist = (int *) NULL;
	out->triangleattributelist = (double *) NULL;
	out->trianglearealist = (double *) NULL;

	out->numberofholes = 0;
	out->holelist = (double *) NULL;

	out->numberofregions = 0;
	out->regionlist = (double *) NULL;

	out->numberofedges = 0;
	out->edgelist = (int *) NULL;
	out->edgemarkerlist = (int *) NULL;
	out->normlist = (double *) NULL;

	//Now to fill them back up again

	out->numberofpoints = num_points;
	out->numberoftriangles = num_tris;
	out->numberofedges = num_edges;

	out->pointlist = new double[num_points*2];
	out->trianglelist = new int[num_tris*2];
	out->edgelist = new int[num_edges*2];

	for(i=0; i<num_points; i++){
		points[i][0] = pointlist[2*i];
		points[i][1] = pointlist[2*i + 1];
	}

	for(i=0; i<num_tris; i++){
		tris[i].i = trianglelist[3*i];
		tris[i+1].j = trianglelist[3*i + 1];
		tris[i+3].k = trianglelist[3*i + 2];
	}

	for(i+0; i<num_edges; i++){
		edges[i].i = edgelist[2*i];
		edges[i].j = edgelist[2*i + 1];
	}
}

//reorder_nodes => will take nodes and reorganize around a single node

//The plan for this is on the board/notebook

//find the -1 in an array
Mesh::Find_negative(int *list){
	int k;

	for(k=0; k<size(list); k++){
		if(list[k] == -1){
			return k;
		}
	}
}

//Add without duplicating
Mesh::add_without_duplicating(int *list, int n){
	int p;
	bool count = true;

	for(p=0; p<size(list) && list[p] != -1; p++){
		if(list[p] == n){                            //Check for double count
			count false;
			break;
		}
		if(p+1 == size(list)){                       //Check for hitting the end
			count false;
		}
	}

	if(count){                                     //If none of the checks above
		list[Find_negative(list)] = n;               //come up positive add element
	}
}

//triangle search function
Mesh::Find_tri(int n, int i) {
	int j;
	int k;
	tri ret;
	for(j=i+1; j<num_tris; j++){
		for(k=0; k<3; k++){
			if(tris[j][k] == n){
				ret.i = j;
				ret.j = j;
				ret.k = k;
				return ret;
			}
		}
		ret.i = -1;
		ret.j = -1;
		ret.k = -1;
		return ret;
	}
}

//edge search function
Mesh::Find_edge(int n, int i){
	int j;
	int k;
	tri ret;
	for(j=i+1; j<num_edges; j++){
		for(k=0; k<2; k++){
			if(edges[j][k] == n){
				ret.i = j;
				ret.j = j;
				ret.k = k;
				return ret;
			}
		}
		ret.i = -1;
		ret.j =-1;
		ret.k = -1;
		return ret;
	}
}

//Fills the updateArr array of the new order of elements
Mesh::updateArrBuilder(int n, int *updateArr){
	int i = 0;
	while(i != -1){
		tri location = Find_tri(n, i);
		int k;

		for(k=0; k<3; k++){
			if(tris[location.i][k] != n){
				add_without_duplicating(tris[location.i][k]);
			}
		}
		i = location.i
	}
}

//Find first non negative 1 element in an array
Mesh::FirstNonNegative(int *Arr){
	int i;

	for(i=0; i<num_points; i++){
		if(Arr[i] != -1){
			return i;
		}
	}

	return -1;
}

//update the arrays for ONE cycle
Mesh::updateArrays(int *updateArr){

	//Pulling an element out to do the swapping
	int finali = FirstNonNegative(updateArr);
	Vector2d finalval = points[finali];
	int element;
	int replacement;

	//The swapping and filling the array with -1 as it does so
	while(element != -1){
		replacement = updateArr[element];

		points[element] = points[replacement];
		swap_tri(element, replacement);
		swap_edge(element, replacement);

		updateArr[element] = -1;
		element = replacement
	}

	//Sliding the first element back in and filling the updateArr with -1 in its place
	points[element] = finalval;
	swap_tri(element, finali);
	swap_edge(element, finali);

	updateArr[element] = -1;
}

//Swaps all nodes that have the value element and replaces it wilh replace
Mesh::swap_tri(int element, int replace){
	int i = 0;

	while(i != -1){
		tri location = Find_tri(element, i);
		i = location.i;
		tris[location.j][location.k] = replace - num_points;
	}
}

Mesh::swap_edge(int element, int replace){
	int i = 0;

	while(i != -1){
		tri location = Find_edge(element, i);
		i = location.i;
		tris[location.j][location.k] = replace - num_points;
	}
}

//Final Form of tris
Mesh::triFinalForm(){
	int i;
	int j;

	for(i=0; i<num_tris; i++){
		for(j=0; j<3; j++){
			tris[i][j] += num_points;
		}
	}
}

Mesh::edgeFinalForm(){
	int i;
	int j;

	for(i=0; i<num_edges; i++){
		for(j=0; j<2; j++){
			edges[i][j] += num_points;
		}
	}
}

Mesh::reorder_nodes(int n){

	//Initialize variables
	int i;
	int start = 0;
	int end = 1;

	//Create and fill the update list
	int * updateArr = new int[num_points];
	updateArr[0] = n;
	for(i=1; i<num_points; i++){
		updateArr[i] = -1;
	}


	//Actual search for nodes as described in notebook pg 49
	while(end != num_points){
		for(i=start; i<end; i++){
			updateArrBuilder(i,updateArr);
		}
		start = end;
		end = Find_negative(updateArr);
	}

	//Wlll go through each cycle and change order of arrays
	while(FirstNonNegative(updateArr) != -1){
		updateArrays(updateArr);
	}

	free(updateArr);
}

//Integrate a function over the entire Mesh
Mesh::integrate(double(*func)(Vector2d)){

	//Declaring some important variables
	int l;
	int i;
	int j;
	int k;

	Vector2d p0;
	Vector2d p1;
	Vector2d p2;

	int kk;

	double kelement;
	double allElements;

	//Filling the xieta array and weight array that will be used in this guassian
	//	integration where the first point in xieta is xi and the second is eta
	Vector2d * xieta = new Vector2d[4];
	xieta[0] = (1./3., 1./3.);
	xieta[1] = (1./5., 1./5.);
	xieta[2] = (1./5., 3./5.);
	xieta[3] = (3./5., 1./5.);

	double * weight = new double[4];
	weight[0] = -27./48.;
	weight[1] = 25./48.;
	weight[2] = 25./48.;
	weight[3] = 25./48.;

	//The actual sum, solving for the sum of each element and then summing over
	//	all elements
	for(l=0; l<num_tris; l++){
		i = tris[l].i;
		j = tris[l].j;
		k = tris[l].k;

		p0 = points[i];
		p1 = points[j];
		p2 = points[k];

		FromRefTri xy(p0, p1, p2);
		kelement = 0;

		for(kk=0; kk<4; kk++){
			kelement += func(xy(xieta[kk]))*weight[k];
		}

		allElements += kelement * xy.jac();
	}

	free(xieta);
	free(weight);
	return allElements;
}

Mesh::testfill(){
  int node;
  int start;
  int i;

  num_points = 30;
  num_tris = 40;
  num_edges = 18;

  points = new Vector2d[num_points];
  tris = new tri[num_tris];
  edges = new edge[num_edges];

  for(node=0; node<num_points; node++){
    points[node] = ((double)(node%5),(double)(node/5));
  }

  for(start=0; start<5; start ++){
    for(i=start; i<25; i+=5){
      tris[6*start + (2*(i/5))].i = i;
      tris[6*start + (2*(i/5)) + 1].i = i;

      tris[6*start + (2*(i/5))].j = i+1;
      tris[6*start + (2*(i/5)) + 1].j = i+5;

      tris[6*start + (2*(i/5))].k = i+5;
      tris[6*start + (2*(i/5)) + 1].k = i+5;
    }
  }

  edges[0]=(0,1);
  edges[1]=(1,2);
  edges[2]=(2,3);
  edges[3]=(3,4);
  edges[4]=(4,9);
  edges[5]=(9,14);
  edges[6]=(14,19);
  edges[7]=(19,24);
  edges[8]=(24,29);
  edges[9]=(29,28);
  edges[10]=(28,27);
  edges[11]=(27,26);
  edges[12]=(26,25);
  edges[13]=(25,20);
  edges[14]=(20,15);
  edges[15]=(15,10);
  edges[16]=(10,5);
  edges[17]=(5,0);
}

Mesh::mass_matrix(double(*Rho)(Vector2d), &massMat){

	//Filling the xieta array and weight array that will be used in this guassian
	//	integration where the first point in xieta is xi and the second is eta
	Vector2d * xieta = new Vector2d[4];
	xieta[0] = (1./3., 1./3.);
	xieta[1] = (1./5., 1./5.);
	xieta[2] = (1./5., 3./5.);
	xieta[3] = (3./5., 1./5.);

	double * weight = new double[4];
	weight[0] = -27./48.;
	weight[1] = 25./48.;
	weight[2] = 25./48.;
	weight[3] = 25./48.;

	double * v = {(1./3., 3./5., 1./5., 1./5.),
								(1./3., 1./5., 1./5., 3./5.),
								(1./3., 1./5., 3./5., 1./5.)}

	int ii, jj, kk;
	double nodeElement;

	for(int t=0; t<num_tris; t++){
		ii = tris[t].i;
		jj = tris[t].j;
		kk = tris[t].k;

		Vector2d p0 = points[ii];
		Vector2d p1 = points[jj];
		Vector2d p2 = points[kk];

		FromRefTri RefTri(p0, p1, p2);
		nodeElement = 0;
		for(int i=0; i<3; i++){
			for(int k=0; k<3; k++){
				for(int k=0; k<4; k++){
					nodeElement += weight[k] * Rho(RefTri(xieta[k])) * v[i][k] * v[j][k];
				}
			}
		}
		massMat.coeffRef(i,j) += nodeElement * RefTri.jac()
	}

	free(xieta);
	free(weight);
}

#endif /* JW_FE_MESH */
