#include "mesh.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>

//-----------------------------------------------
// Public arrays
//	-xieta: Xi ->  Vector2d[0]
//					Eta -> Vector2d[1]
//					From Guassian Lagrangian quadrature.
//	-weight: Coresponding weight to Xi-Eta pair
//	-v[i][j]: i refers to which v, 0,1,2, from
//							pg 18 in notebook.
//						j refers to corresponding Xi-Eta pair
//	-grad_v: gradient of v_i with respet to Xi and
//		Eta, respectively.

Vector2d Mesh::xieta[4] = {
	Vector2d(1./3., 1./3.),
	Vector2d(1./5., 1./5.),
	Vector2d(1./5., 3./5.),
	Vector2d(3./5., 1./5.)
};

double Mesh::weight[4] = {-27./48., 25./48., 25./48.,	25./48.};

double Mesh::v[3][4] = {{1./3., 3./5., 1./5., 1./5.},
												{1./3., 1./5., 1./5., 3./5.},
												{1./3., 1./5., 3./5., 1./5.}};
Vector2d Mesh::grad_v[3] = {
	Vector2d(-1, -1),
	Vector2d(0, 1),
	Vector2d(1, 0)};

//###############################################
// Constructor -> fills data in Mesh, assuming
//	triangle has already run

Mesh::Mesh(triangulateio *in){

	//setting the values of the data to values from triangle

	num_points = in->numberofpoints;
	num_tris = in->numberoftriangles;
	num_edges = in->numberofsegments;

	//creating arrays

	points = new Vector2d[num_points];
	tris = new tri[num_tris];
	edges = new edge[num_edges];

	bound = NULL;
	not_bound = NULL;

	//Now to fill them
	//----------------------------------------------
	//First is the point list, assuming 2 doubles
	//	labeled x and y

	for(int i=0; i<num_points; i++){
		points[i](0) = in->pointlist[2*i+0];
		points[i](1) = in->pointlist[2*i+1];
	}

	//----------------------------------------------
	//Now the triangle index, assuming 3 ints
	//	labeled i,j,k

	for(int t=0; t<num_tris; t++){
		tris[t].i = in->trianglelist[3*t+0];
		tris[t].j = in->trianglelist[3*t+1];
		tris[t].k = in->trianglelist[3*t+2];
	}

	//----------------------------------------------
	//Now the edge index, assuming 3 ints labeled
	//	i,j

	for(int i=0; i<num_edges; i++){
		edges[i].i = in->segmentlist[2*i+0];
		edges[i].j = in->segmentlist[2*i+1];
	}
}

//###############################################
// Deconstructor -> frees all arrays created
//	in Constructor

Mesh::~Mesh(){
	// Free memory associated with points if not
	//	already freed.
	if(points != NULL){
		delete points;
	}

	// Free memory associated with tris if not
	//	already freed.
	if(tris != NULL){
		delete tris;
	}

	// Free memory associated with edges if not
	//	already freed.
	if(edges != NULL){
		delete edges;
	}

	if(not_bound != NULL){
		delete[] not_bound;
	}

	if(bound != NULL){
		delete[] bound;
	}
}

//###############################################
// to_triangulateio -> creates a triangulateio
//	to refine the mesh by passing it back to
//	triangle.c

void Mesh::to_triangulateio(triangulateio *out){
	int i;

	// First to initialize the triangulateio
	//	 arays and variables

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

	//Now to fill them

	out->numberofpoints = num_points;
	out->numberoftriangles = num_tris;
	out->numberofsegments = num_edges;

	out->pointlist = new double[num_points*2];
	out->trianglelist = new int[num_tris*2];
	out->edgelist = new int[num_edges*2];

	for(i=0; i<num_points; i++){
		points[i][0] = out->pointlist[2*i];
		points[i][1] = out->pointlist[2*i + 1];
	}

	for(i=0; i<num_tris; i++){
		tris[i].i = out->trianglelist[3*i];
		tris[i+1].j = out->trianglelist[3*i + 1];
		tris[i+3].k = out->trianglelist[3*i + 2];
	}

	for(i=0; i<num_edges; i++){
		edges[i].i = out->edgelist[2*i];
		edges[i].j = out->edgelist[2*i + 1];
	}
}

//###############################################
// reorder_nodes -> will take nodes and reorganize
//	around a single node. The General plan for
//	this begins on pg 48 i nthe notebook.
// The actual reorder_nodes is at the bottom.
//	These are the supportion fucntions

//-----------------------------------------------
// find_negative -> find the index of the first
//	negative 1 in an array
int Mesh::find_negative(int *arr){
	int k;

	// Go throught the array of length, num_points,
	//	and find the first negative 1 and return
	//	the index of that negative 1.
	for(k=0; k<num_points; k++){
		if(arr[k] == -1){
			return k;
		}
	}

	// If no negative 1's are found, return
	//	negative 1.
	return -1;
}

//-----------------------------------------------
// add_without_duplicating -> add an element to
//	to an array if and only if that element is
//	not already in that array.
void Mesh::add_without_duplicating(int arr[], int n){
	int p;
	bool count = true;

	// First goes through the update array until it
	//	either hits the end of the array or
	//	encounters a negative 1.

	for(p=0; (p<num_points) && (arr[p] != -1); p++){

		// Checks to see if that element is in the array.
		if(arr[p] == n){
			count = false;
			break;
		}

		// Checks if it hits the end of the array.
		if(p == num_points-1){
			count = false;
		}
	}

	// IF the element is not in the array and if the
	//	checker never hits the end of the array, then
	//	replace the first negative 1 with the element.
	if(count){
		int index = find_negative(arr);
		arr[index] = n;
	}
}

//-----------------------------------------------
// find_tri -> find the index of a node in a
//	triangle if it is an element of a triangle,
//	else return -1.
int Mesh::find_tri(int n,int t) {
	int k;

	// Go through each element in triangle, t, and
	//	check if node, n, is an element. If it is
	//	return the index of the element in the
	//	triangle. If it isn't, return -1.
	for(k=0; k<3; k++){
		if(tris[t][k] == n){
			return t;
		}
	}
	return -1;
}

//-----------------------------------------------
// find_edge -> find the index of a node in an
//	edge if it is an element of an edge, else
//	return -1.
int Mesh::find_edge(int n, int e){
	int k;

	// Go through each element in edge, e, and
	//	check if node, n, is an element. If it is
	//	return the index of the element in the
	//	edge. If it isn't, return -1.
	for(k=0; k<2; k++){
		if(edges[e][k] == n){
			return e;
		}
	}
	return -1;
}

//-----------------------------------------------
// update_arr_builder -> builds an array of the
//	new order that the nodes should be in.
void Mesh::update_arr_builder(int n, int * update_arr){
	int t;
	int location;
	int k;

	// Go through each triangle and find the triangles
	//	with node, n.
	for(t=0; t<num_tris; t++){
		location = find_tri(n, t);
		if(location != -1){
			for(k=0; k<3; k++){

				// For triangles containing node, n, add the
				//	elements to the new order list if and
				//	only if those elements aren't already
				//	contained in the new order list.
				if((tris[location][k] != n) && (tris[location][k] != -1)){
					add_without_duplicating(update_arr, tris[location][k]);
				}
			}
		}
	}
}

//-----------------------------------------------
// first_non_negative -> finds the index of first
//	non-negative 1 value in an array.
int Mesh::first_non_negative(int *arr){
	int i;

	// Go through an array of length num_points and
	//	find the first non-negative 1 index, else
	//	return -1.
	for(i=0; i<num_points; i++){
		if(arr[i] != -1){
			return i;
		}
	}

	return -1;
}

//-----------------------------------------------
// update_arrays -> updates the points, tris,
//	and edges according to one cycle built from
//	the update array. For info on the cycles,
//	look in the notebook pg 49-52.
void Mesh::update_arrays(int *update_arr){

	// This is done with a cycle built from the
	//	update array. It is found by first taking
	//	the first non-negative 1 element of the
	//	update array and saving the index and
	//	the points value.
	int final_i = first_non_negative(update_arr);
	Vector2d final_val = points[final_i];

	// Next, it finds the next element in the
	//	cycle by looking at the saved index above
	//	and going to update_arr[saved index]. This
	//	wll yield another value that then becomes
	//	our replacement index.
	int element = final_i;
	int replacement = update_arr[element];

	// This is donw until the replacement value
	//	circles back around to initial index that
	//	was saved. This marks the completion of a
	//	cycle. A single update_arr can have several
	//	of these cycles. This only completes one
	//	cycle, so this funtion must be repeated
	//	until the update_arr is filled with -1.
	//	Values that have already been cycled
	//	through are replaced with -1's in the
	//	update_arr. It then makes the proper
	//	replacements in points, tris, and edges.
	while(replacement != final_i){

		points[element] = points[replacement];
		swap_tri(element, replacement);
		swap_edge(element, replacement);

		update_arr[element] = -1;
		element = replacement;
		replacement = update_arr[element];
	}

	// The final replacement must be done outside
	//	of the while loop due to the nature of
	//	cycles. Without this step, one would find
	//	themselves deleting the information of the
	//	first index.
	points[element] = final_val;
	swap_tri(element, final_i);
	swap_edge(element, final_i);

	update_arr[element] = -1;
}

//-----------------------------------------------
// swap_tri -> swaps a given node in tris with
//	a replacement node minus num_points to
//	distinguish changed nodes apart from unchanged
//	nodes.
void Mesh::swap_tri(int replace, int element){
	int i;
	int t;

	// Goes through each triangle and searches for
	//	some element and replaces it with the
	//	replacement minus num_points. This is done
	//	to distinguish changed elements apart.
	for(t=0; t<num_tris; t++){
		for( i=0; i<3; i++){
			if(tris[t][i] == element){
				tris[t][i] = replace - num_points;
			}
		}
	}
}

//-----------------------------------------------
// swap_edge -> swaps a given node in edges with
//	a replacement node minus num_points to
//	distinguish changed nodes apart from unchanged
//	nodes.
void Mesh::swap_edge(int replace, int element){
	int i;
	int e;

	// Goes through each edge and searches for
	//	some element and replaces it with the
	//	replacement minus num_points. This is done
	//	to distinguish changed elements apart.
	for(e=0; e<num_edges; e++){
		for (i=0; i<2; i++){
			if(edges[e][i] == element){
				edges[e][i] = replace - num_points;
			}
		}
	}
}

//-----------------------------------------------
// tri_final_form -> adds num_points to all nodes
//	in tris. This will be applied once all nodes
//	have been changed, so all nodes are equal to
//	their replacements.
void Mesh::tri_final_form(){
	int i;
	int j;

	// Before this function is called, all the
	//	triangle elements have been replaced. They
	//	need to have num_points added to them to
	//	undo the changes done in swap_tri.
	for(i=0; i<num_tris; i++){
		for(j=0; j<3; j++){
			tris[i][j] += num_points;
		}
	}
}

//-----------------------------------------------
// edge_final_form -> adds num_points to all nodes
//	in edges. This will be applied once all nodes
//	have been changed, so all nodes are equal to
//	their replacements.
void Mesh::edge_final_form(){
	int i;
	int j;

	// Before this function is called, all the
	//	edge elements have been replaced. They
	//	need to have num_points added to them to
	//	undo the changes done in swap_edge.
	for(i=0; i<num_edges; i++){
		for(j=0; j<2; j++){
			edges[i][j] += num_points;
		}
	}
}

//-----------------------------------------------
// reorder_nodes -> bringing all the above
//	functions together.
void Mesh::reorder_nodes(int n){
	int i;
	int start = 0;
	int end = 1;

	// Create and fill the update list with negative
	//	1's except for the first point which is the
	//	node that the other nodes are being reordered
	//	around.
	int * update_arr = new int[num_points];
	update_arr[0] = n;
	for(i=1; i<num_points; i++){
		update_arr[i] = -1;
	}

	// Building the array of the new order of nodes.
	while(end != -1){
		for(i=start; i<end; i++){
			update_arr_builder(update_arr[i],update_arr);
		}
		start = end;
		end = find_negative(update_arr);
	}

	// This will first to check to see if there are
	//	any cycles left by chekcing for any non-negative
	//	1's in the update_arr. When it finds a cycle,
	//	i.e. a non-negative element, it calls update_arrays.
	//	It will continue this until all cycles are
	//	completed.
	while(first_non_negative(update_arr) != -1){
		update_arrays(update_arr);
	}

	delete[] update_arr;

	// Finally, it puts tris and edges in their final
	//	form.
	tri_final_form();

	edge_final_form();
}

//###############################################
// integrate -> use Guassian Lagrangian
//	numberical integration technique to integrate
//	over a surface for some given function, func.
double Mesh::integrate(double(*func)(Vector2d)){
	int l;
	int i;
	int j;
	int k;

	Vector2d p0;
	Vector2d p1;
	Vector2d p2;

	int kk;

	double k_element;
	double all_elements;

	// This numerically integrates using the Gaussian
	//	Legrangian method described in detail on
	//	pg 53 in the notebook.
	all_elements = 0.;
	for(l=0; l<num_tris; l++){
		i = tris[l].i;
		j = tris[l].j;
		k = tris[l].k;

		p0 = points[i];
		p1 = points[j];
		p2 = points[k];

		// This method relies that the points are in
		//	Xi-Eta space. This is what FromRefTri
		//	defined in utils.hpp will do. This uses
		//	a number of class defined variables
		//	including the weights and the matching
		//	Xi-Eta pairs. More details on pg 53 in
		//	the notebook.
		FromRefTri xy(p0, p1, p2);
		k_element = 0.;

		for(kk=0; kk<4; kk++){
			k_element += func(xy(xieta[kk]))*weight[kk];
		}

		all_elements += 0.5 * k_element * std::abs(xy.jac());
	}

	return all_elements;
}

//###############################################
// Functions for mass_matrix and stiffness_matrix

//-----------------------------------------------
// find_not_bound -> find the index of a node in
//	not_bound if it is an element of not_bound,
//	else return -1.

int Mesh::find_not_bound(int n){
	int b;

	// Go through each element in not bound, b,
	//	and check if node, n, is an element. If it
	//	is, return the index of the element in
	//	not_bound. If it isn't, return -1.

	for(b=0; b<num_points - num_edges; b++){
		if(not_bound[b] == n){
			return b;
		}
	}
	return -1;
}

//-----------------------------------------------
// find_bound -> find the index of a node in
//	bound if it is an element of bound,
//	else return -1.

int Mesh::find_bound(int n){
	int b;

	// Go through each element in bound, b,
	//	and check if node, n, is an element. If it
	//	is, return the index of the element in
	//	bound. If it isn't, return -1.

	for(b=0; b<num_edges; b++){
		if(bound[b] == n){
			return b;
		}
	}
	return -1;
}

//###############################################
// mass_matrix -> build a matrix that will use
//	the density function to solve for a
//	Sparse Matrix that will be used to find the
//	work done.
void Mesh::mass_matrix(double(*rho)(Vector2d), SparseMatrix<double> &bound_mat, SparseMatrix<double> &not_bound_mat){
	int ii;
	int jj;
	int kk;
	double node_element;
	bool on_edge, add_to;
	int row, col;

	// Building array of points not on the boundary
	// 	 array and the points on the boundary array.
	// It will also check to see if bound and
	// 	not_bound are initialized and if they are
	// 	not, it will initialize them.

	if(not_bound == NULL){
		not_bound = new int[num_points-num_edges];
	}
	if(bound == NULL){
		bound = new int[num_edges];
	}

	for(int i=0; i<num_points-num_edges; i++){
		not_bound[i] = -1;
	}
	for(int i=0; i<num_edges; i++){
		bound[i] = -1;
	}

	for(int p=0; p<num_points; p++){
		on_edge = false;
		for(int e=0; e<num_edges; e++){
			if(find_edge(p,e) != -1){
				on_edge = true;
				break;
			}
		}
		if(!on_edge){
			for(int i=0; i<num_points - num_edges; i++){
				if(not_bound[i] == -1){
					not_bound[i] = p;
					break;
				}
			}
		} else {
			for(int i=0; i<num_edges; i++){
				if(bound[i] == -1){
					bound[i] = p;
					break;
				}
			}
		}
	}


	// This has a very similar set up to the Guassian
	//	Legrangian integration described above. This
	//	is described in detail on pg 55 in the
	//	notebook.
	for(int t=0; t<num_tris; t++){
		ii = tris[t].i;
		jj = tris[t].j;
		kk = tris[t].k;

		Vector2d p0 = points[ii];
		Vector2d p1 = points[jj];
		Vector2d p2 = points[kk];

		// This is done in Xi-Eta space. FromRefTri
		//	as defined in utils.hpp will do the
		//	conversion.
		FromRefTri ref_tri(p0, p1, p2);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				node_element = 0;
				on_edge = false;
				add_to = false;
				for(int k=0; k<4; k++){
					node_element += .5 * weight[k] * rho(ref_tri(xieta[k])) * v[i][k] * v[j][k];
				}

				for(int e=0; e<num_edges; e++){
					if((find_edge(tris[t][j], e) != -1) && (find_not_bound(tris[t][i]) != -1)){
							on_edge = true;
							col = find_bound(tris[t][j]);
							row = find_not_bound(tris[t][i]);
							break;
					} else if((find_not_bound(tris[t][i]) != -1) && (find_not_bound(tris[t][j]) != -1)){
						add_to = true;
						row = find_not_bound(tris[t][i]);
						col = find_not_bound(tris[t][j]);
						break;
					}
				}

				if(on_edge){
					bound_mat.coeffRef(col, row) += node_element * std::abs(ref_tri.jac());
				} else if(add_to){
					not_bound_mat.coeffRef(col, row) += node_element * std::abs(ref_tri.jac());
				}
			}
		}
		// cout << "Triangles Processed(Mass): " << t << "/" << num_tris << endl << std::flush;
	}
}

//###############################################
// The next function was created as a supporting
//	funciton for stiffness_matrix.

//-----------------------------------------------
// jacobian -> creates a 2 by 2 jacobian matrix
//	based on 3 points in x-y space, not Xi-Eta
//	space.
Matrix2d Mesh::jacobian(Vector2d p0, Vector2d p1, Vector2d p2){
	double grad_a;
	double grad_b;
	double grad_c;
	double grad_d;
	double grad_const;

	// Builds each element of the jacobian as described
	//	on pg 57-58 in the notebook.

	grad_a = p1[0] - p0[0];
	grad_b = p2[0] - p0[0];
	grad_c = p1[1] - p0[1];
	grad_d = p2[1] - p0[1];
	grad_const = (grad_a * grad_d) - (grad_b * grad_c);

	Matrix2d jacobian;

	jacobian(0,0) = grad_d/grad_const;
	jacobian(1,0) = -grad_b/grad_const;
	jacobian(0,1) = -grad_c/grad_const;
	jacobian(1,1) = grad_a/grad_const;

	return jacobian;
}

//-----------------------------------------------
// stiffness_matrix -> build a matrix that will use
//	the stiffness function to solve for a
//	Sparse Matrix that will be used to find the
//	work done.
void Mesh::stiffness_matrix(double(*stiff)(Vector2d), SparseMatrix<double> &bound_mat, SparseMatrix<double> &not_bound_mat){
	int ii;
	int jj;
	int kk;
	double node_element;
	bool on_edge, add_to;
	int row, col;

	// Building array of points not on the boundary
	// 	 array and the points on the boundary array.
	// It will also check to see if bound and
	// 	not_bound are initialized and if they are
	// 	not, it will initialize them.

	if(not_bound == NULL){
		not_bound = new int[num_points-num_edges];
	}
	if(bound == NULL){
		bound = new int[num_edges];
	}

	for(int i=0; i<num_points-num_edges; i++){
		not_bound[i] = -1;
	}
	for(int i=0; i<num_edges; i++){
		bound[i] = -1;
	}

	for(int p=0; p<num_points; p++){
		on_edge = false;
		for(int e=0; e<num_edges; e++){
			if(find_edge(p,e) != -1){
				on_edge = true;
				break;
			}
		}
		if(!on_edge){
			for(int i=0; i<num_points - num_edges; i++){
				if(not_bound[i] == -1){
					not_bound[i] = p;
					break;
				}
			}
		} else {
			for(int i=0; i<num_edges; i++){
				if(bound[i] == -1){
					bound[i] = p;
					break;
				}
			}
		}
	}

	// Built similar to integrate or mass_matrix.
	for(int t=0; t<num_tris; t++){
		ii = tris[t].i;
		jj = tris[t].j;
		kk = tris[t].k;

		Vector2d p0 = points[ii];
		Vector2d p1 = points[jj];
		Vector2d p2 = points[kk];

		// This must be done in Xi-Eta space.
		FromRefTri ref_tri(p0, p1, p2);
		Matrix2d jacob = jacobian(p0, p1, p2);

		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				node_element = 0;
				on_edge = false;
				add_to = false;

				// grad_v is transposed in the first part
				//	because Eigen is a column major lib,
				//	so Vector2d is a 2 by 1 and not a
				//	1 by 2.
				double gradient = (grad_v[j].transpose() * jacob) * (jacob.transpose() * grad_v[i]);

				for(int k=0; k<4; k++){

					node_element += .5 * weight[k] * stiff(ref_tri(xieta[k])) * gradient;
				}

				if((find_bound(tris[t][j]) != -1) && (find_not_bound(tris[t][i]) != -1)){
					// cout << "i: " << tris[t][i] << endl;
					// cout << "j: " << tris[t][j] << endl;
					on_edge = true;
					col = find_bound(tris[t][j]);
					row = find_not_bound(tris[t][i]);
				} else if((find_not_bound(tris[t][j]) != -1) && (find_not_bound(tris[t][i]) != -1)){
					add_to = true;
					row = find_not_bound(tris[t][i]);
					col = find_not_bound(tris[t][j]);
				}

				if(on_edge){
					bound_mat.coeffRef(col, row) += node_element * std::abs(ref_tri.jac());
				} else if(add_to){
					not_bound_mat.coeffRef(col, row) += node_element * std::abs(ref_tri.jac());
				}
			}
		}
	}
}

//###############################################
// points_get -> This function is now pointless
//	because points is a public variable.
Vector2d Mesh::points_get(int i){
	// Typical getter.
	return points[i];
}
