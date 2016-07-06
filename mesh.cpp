#include "mesh.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>

//Added for integrate, mass_matrix stiffness_matrix
//Filling the xieta array and weight array that will be used in this guassian
//	integration where the first point in xieta is xi and the second is eta

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
double Mesh::grad_v[3][2] = {{-1, -1},
															{1, 0},
															{0, 1}};

//Constructor => fills data in Mesh, assuming triangle has already run

Mesh::Mesh(triangulateio *in){

	//setting the values of the data to values from triangle

	num_points = in->numberofpoints;
	num_tris = in->numberoftriangles;
	num_edges = in->numberofsegments;

	//creating arrays

	points = new Vector2d[num_points];
	tris = new tri[num_tris];
	edges = new edge[num_edges];

	//Now to fill them
	//-----------------------------
	//First is the point list, assuming 2 doubles labeled x and y

	for(int i=0; i<num_points; i++){
		points[i](0) = in->pointlist[2*i+0];
		points[i](1) = in->pointlist[2*i+1];
	}

	//-------------------------------------
	//Now the triangle index, assuming 3 ints labeled, i,j,k

	for(int t=0; t<num_tris; t++){
		tris[t].i = in->trianglelist[3*t+0];
		tris[t].j = in->trianglelist[3*t+1];
		tris[t].k = in->trianglelist[3*t+2];
	}

	//--------------------------------------
	//Now the edge index, assuming 3 ints labeled, i,j

	for(int i=0; i<num_edges; i++){
		edges[i].i = in->segmentlist[2*i+0];
		edges[i].j = in->segmentlist[2*i+1];
	}
}

//Deconstructor => frees all arrays created in Constructor

Mesh::~Mesh(){
	if(points != NULL){
		delete points;
	}

	if(tris != NULL){
		delete tris;
	}

	if(edges != NULL){
		delete edges;
	}
}

//to_triangulateio => creates a triangulateio to refine the mesh

void Mesh::to_triangulateio(triangulateio *out){
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

//reorder_nodes => will take nodes and reorganize around a single node

//The plan for this is on the board/notebook

//find the -1 in an array
int Mesh::find_negative(int *arr){
	int k;

	for(k=0; k<num_points; k++){
		if(arr[k] == -1){
			return k;
		}
	}

	return -1;
}

//Add without duplicating
void Mesh::add_without_duplicating(int arr[], int n){
	int p;
	bool count = true;



	for(p=0; (p<num_points) && (arr[p] != -1); p++){
		if(arr[p] == n){                            //Check for double count
			count = false;
			break;
		}
		if(p == num_points-1){                      //Check for hitting the end
			count = false;
		}
	}
	if(count){                                     //If none of the checks above
		int index = find_negative(arr);
		arr[index] = n;     						            //come up positive add element
	}
}

//triangle search function
int Mesh::find_tri(int n,int t) {
	int k;
	for(k=0; k<3; k++){
		if(tris[t][k] == n){
			return t;
		}
	}
	return -1;
}

//edge search function
int Mesh::find_edge(int n, int e){
	int k;

	for(k=0; k<2; k++){
		if(edges[e][k] == n){
			return e;
		}
	}
	return -1;
}

//Fills the update_arr array of the new order of elements
void Mesh::update_arr_builder(int n, int * update_arr){
	int t;
	int location;
	int k;

	for(t=0; t<num_tris; t++){
		location = find_tri(n, t);
		if(location != -1){
			for(k=0; k<3; k++){
				if((tris[location][k] != n) && (tris[location][k] != -1)){
					add_without_duplicating(update_arr, tris[location][k]);
				}
			}
		}
	}
}

//Find first non negative 1 element in an array
int Mesh::first_non_negative(int *arr){
	int i;

	for(i=0; i<num_points; i++){
		if(arr[i] != -1){
			return i;
		}
	}

	return -1;
}

//update the arrays for ONE cycle
void Mesh::update_arrays(int *update_arr){

	//Pulling an element out to do the swapping
	int final_i = first_non_negative(update_arr);
	Vector2d final_val = points[final_i];

	int element = final_i;
	int replacement = update_arr[element];

	//The swapping and filling the array with -1 as it does so
	while(replacement != final_i){

		points[element] = points[replacement];
		swap_tri(element, replacement);
		swap_edge(element, replacement);

		update_arr[element] = -1;
		element = replacement;
		replacement = update_arr[element];
	}

	//Sliding the first element back in and filling the update_arr with -1 in its place
	points[element] = final_val;
	swap_tri(element, final_i);
	swap_edge(element, final_i);

	update_arr[element] = -1;
}

//Swaps all nodes that have the value element and replaces it wilh replace
void Mesh::swap_tri(int replace, int element){
	int i;
	int t;

	for(t=0; t<num_tris; t++){
		for( i=0; i<3; i++){
			if(tris[t][i] == element){
				tris[t][i] = replace - num_points;
			}
		}
	}
}

void Mesh::swap_edge(int replace, int element){
	int i;
	int e;

	for(e=0; e<num_edges; e++){
		for (i=0; i<2; i++){
			if(edges[e][i] == element){
				edges[e][i] = replace - num_points;
			}
		}
	}
}


//Final Form of tris
void Mesh::tri_final_form(){
	int i;
	int j;

	for(i=0; i<num_tris; i++){
		for(j=0; j<3; j++){
			tris[i][j] += num_points;
		}
	}
}

void Mesh::edge_final_form(){
	int i;
	int j;

	for(i=0; i<num_edges; i++){
		for(j=0; j<2; j++){
			edges[i][j] += num_points;
		}
	}
}

void Mesh::reorder_nodes(int n){

	//Initialize variables
	int i;
	int start = 0;
	int end = 1;

	//Create and fill the update list
	int * update_arr = new int[num_points];
	update_arr[0] = n;
	for(i=1; i<num_points; i++){
		update_arr[i] = -1;
	}

	//Actual search for nodes as described in notebook pg 49
	while(end != -1){
		for(i=start; i<end; i++){
			update_arr_builder(update_arr[i],update_arr);
		}
		start = end;
		end = find_negative(update_arr);
	}

	//Wlll go through each cycle and change order of arrays
	while(first_non_negative(update_arr) != -1){
		update_arrays(update_arr);
	}

	delete[] update_arr;

	tri_final_form();

	edge_final_form();
}

//Integrate a function over the entire Mesh
double Mesh::integrate(double(*func)(Vector2d)){

	//Declaring some important variables
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

	//The actual sum, solving for the sum of each element and then summing over
	//	all elements
	all_elements = 0.;
	for(l=0; l<num_tris; l++){
		i = tris[l].i;
		j = tris[l].j;
		k = tris[l].k;

		p0 = points[i];
		p1 = points[j];
		p2 = points[k];

		FromRefTri xy(p0, p1, p2);
		k_element = 0.;

		for(kk=0; kk<4; kk++){
			k_element += func(xy(xieta[kk]))*weight[kk];
		}

		all_elements += 0.5 * k_element * std::abs(xy.jac());
	}

	return all_elements;
}



void Mesh::mass_matrix(double(*rho)(Vector2d), SparseMatrix<double, RowMajor, int> &mass_mat){

	int ii;
	int jj;
	int kk;
	double node_element;

	for(int t=0; t<num_tris; t++){
		ii = tris[t].i;
		jj = tris[t].j;
		kk = tris[t].k;

		Vector2d p0 = points[ii];
		Vector2d p1 = points[jj];
		Vector2d p2 = points[kk];
		FromRefTri ref_tri(p0, p1, p2);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				node_element = 0;
				for(int k=0; k<4; k++){
					node_element += weight[k] * rho(ref_tri(xieta[k])) * v[i][k] * v[j][k];
				}
				mass_mat.coeffRef(tris[t][i], tris[t][j]) += node_element * ref_tri.jac();
			}
		}
	}
}

Vector2d Mesh::gradient_v(int i,  Vector2d p0, Vector2d p1, Vector2d p2){
	double grad_a;
	double grad_b;
	double grad_c;
	double grad_d;
	double grad_const;

	grad_a = p1[0] - p0[0];
	grad_b = p2[0] - p0[0];
	grad_c = p1[1] - p0[1];
	grad_d = p2[1] - p0[1];
	grad_const = 1./((grad_a * grad_d) - (grad_b * grad_c));

	double xix = grad_d/grad_const;
	double etax = -grad_c/grad_const;
	double xiy = -grad_b/grad_const;
	double etay = grad_a/grad_const;

	Vector2d ret;

	ret[0] = (grad_v[i][0] * xix) + (grad_v[i][1] * etax);
	ret[1] = (grad_v[i][0] * xiy) + (grad_v[i][1] * etay);

	return ret;
}

double Mesh::dot(Vector2d a, Vector2d b){
	return (a[0]*b[0]) + (a[1]*b[1]);
}

void Mesh::stiffness_matrix(double(*stiff)(Vector2d), SparseMatrix<double, RowMajor, int> &stiff_mat){
	int ii;
	int jj;
	int kk;
	double node_element;

	for(int t=0; t<num_tris; t++){
		ii = tris[t].i;
		jj = tris[t].j;
		kk = tris[t].k;

		Vector2d p0 = points[ii];
		Vector2d p1 = points[jj];
		Vector2d p2 = points[kk];

		FromRefTri ref_tri(p0, p1, p2);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				node_element = 0;
				for(int k=0; k<4; k++){
					node_element += weight[k] * stiff(ref_tri(xieta[k])) * dot(gradient_v(i, p0, p1,p2), gradient_v(j, p0, p1, p2));
				}
				stiff_mat.coeffRef(tris[t][i], tris[t][j]) += node_element * ref_tri.jac();
			}
		}
	}
}

Vector2d Mesh::points_get(int i){
	return points[i];
}
