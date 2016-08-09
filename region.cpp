#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
//#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
//#include <Eigen/CholmodSupport>

#include "mesh.hpp"
#include "triangle.hpp"
#include "utils.hpp"

#define PI 3.14159265358979323846E+0

using namespace std;

typedef struct triangulateio tio;

namespace{

	void add_point(double *pl, int *n, double x, double y){
		pl[2*(*n)+0] = x;
		pl[2*(*n)+1] = y;
		(*n)++;
	}

	void free_if_not_null( void *ptr){
		if( ptr != NULL){
			free(ptr);
		}
	}

	double func(Vector2d point){
		double ret;
		ret = point[0] * point[1];
		return ret;
	}
}


int num_circ(double r){
	double dx = 0.1;
	return static_cast<int>(2*PI*r/dx);
}

tio * malloc_tio(){
	tio* in = (tio *) malloc(sizeof(tio));

	in->numberofpoints = 0;
	in->numberofpointattributes = 0;
	in->pointlist = (double *) NULL;
	in->pointattributelist = (double *) NULL;
	in->pointmarkerlist = (int *) NULL;
	in->numberoftriangles = 0;
	in->numberofcorners = 0;
	in->numberoftriangleattributes = 0;
	in->trianglelist = (int *) NULL;
	in->triangleattributelist = (double *) NULL;
	in->trianglearealist = (double *) NULL;
	in->numberofsegments = 0;
	in->segmentlist = (int *) NULL;
	in->segmentmarkerlist = (int *) NULL;
	in->numberofholes = 0;
	in->holelist = (double *) NULL;
	in->numberofregions = 0;
	in->regionlist = (double *) NULL;
	in->numberofedges = 0;
	in->edgelist = (int *) NULL;
	in->edgemarkerlist = (int *) NULL;
	in->normlist = (double *) NULL;

	return in;
}

tio * malloc_pslg(int num_points){
	// Allocate the triangulateio object to hold the
	// planar straight line graph
	tio* in = malloc_tio();

	// Allocate the memory to hold the PSLG
	in->numberofpoints = num_points;
	in->pointlist = (double *) malloc(2*in->numberofpoints*sizeof(double));
	in->pointmarkerlist = (int *) malloc(in->numberofpoints*sizeof(double));
	in->numberofsegments = num_points;
	in->segmentlist = (int *) malloc(2*in->numberofsegments*sizeof(int));
	in->segmentmarkerlist = (int *) malloc(in->numberofsegments*sizeof(int));
	in->numberofholes = 2;
	in->holelist = (double *) malloc(2*in->numberofholes*sizeof(double));

	return in;
}

void free_tio(tio *in){
	free_if_not_null( (void *) in->pointlist);
	free_if_not_null( (void *) in->pointattributelist);
	free_if_not_null( (void *) in->pointmarkerlist);
	free_if_not_null( (void *) in->trianglelist);
	free_if_not_null( (void *) in->triangleattributelist);
	free_if_not_null( (void *) in->trianglearealist);
	free_if_not_null( (void *) in->neighborlist);
	free_if_not_null( (void *) in->segmentlist);
	free_if_not_null( (void *) in->segmentmarkerlist);
	free_if_not_null( (void *) in->holelist);
	free_if_not_null( (void *) in->regionlist);
	free_if_not_null( (void *) in->edgelist);
	free_if_not_null( (void *) in->edgemarkerlist);
	free_if_not_null( (void *) in->normlist);
}

void fill_pslg(tio *in, double d, double r1, double r2){
	int n = 0;
	int first, i;
	int nc1, nc2;

	double *pl;
	int *pml, *sl, *sml;

	double dth;

	pl = in->pointlist;
	pml = in->pointmarkerlist;
	sl = in->segmentlist;
	sml = in->segmentmarkerlist;

	first = n;

	// Outer Shape of the mesh: an elipse
	// 	-The shape comes from the last
	// 		2 parameters in add_point()

	for(i=0; i<56; i++){
		double theta = 2*PI*i/56;
		add_point(pl, &n, 1.5*d*cos(theta), d*sin(theta));
	}
	// Four corners of surrounding box

	// add_point(pl, &n, 1.5*d, d);
	// add_point(pl, &n, 1.5*d, -d);
	// add_point(pl, &n, -1.5*d, -d);
	// add_point(pl, &n, -1.5*d, d);
	for( i=first; i<n; i++){
		pml[i] = 1;
		sl[2*i+0] = i;
		sl[2*i+1] = i+1;
		sml[i] = 1;
	}
	sl[2*(n-1)+1] = first;
	// The first hole
	first=n;
	nc1 = num_circ(r1);
	dth = 2.*PI/nc1;
	for( i=0; i<nc1; i++){
		add_point(pl, &n, -d/2.+cos(dth*i), sin(dth*i));
	}
	for( i=first; i<n; i++){
		pml[i] = 1;
		sl[2*i+0] = i;
		sl[2*i+1] = i+1;
		sml[i] = 1;
	}
	sl[2*(n-1)+1] = first;
	in->holelist[2*0+0] = -d/2.;
	in->holelist[2*0+1] = 0;
	// the second hole
	first=n;
	nc2 = num_circ(r2);
	dth = 2.*PI/nc2;
	for( i=0; i<nc2; i++){
		add_point(pl, &n, +d/2.+cos(dth*i), sin(dth*i));
	}
	for( i=first; i<n; i++){
		pml[i] = 1;
		sl[2*i+0] = i;
		sl[2*i+1] = i+1;
		sml[i] = 1;
	}
	sl[2*(n-1)+1] = first;
	in->holelist[2*1+0] = d/2.;
	in->holelist[2*1+1] = 0.;
}

void print_mesh( fstream &f, Mesh &m){
	int i;

	f << "Points " << m.num_points << endl;
	for( i=0; i<m.num_points; i++){
		f << m.points[i][0] << " " << m.points[i][1] << endl;
	}
	f << endl;

	f << "Triangles " << m.num_tris << endl;
	for( i=0; i<m.num_tris; i++){
		f << m.tris[i][0] << " " << m.tris[i][1] << " " << m.tris[i][2] << endl;
	}
	f << endl;
}

void print_mat(Mesh &m, SparseMatrix<double> &not_bound_mat_ordered, SparseMatrix<double> &bound_mat_ordered){

	// Ordered Mass Matrix

	// Nodes on the boundary

  ofstream ordered_bound_mat;
  ordered_bound_mat.open("ordered_bound_mat.txt");
  if(ordered_bound_mat.is_open() == false){
    cout << "Unable to open ordered_bound_mat" << endl << std::flush;
  }
    ordered_bound_mat << m.num_points<< "\n";
		for(int i=0; i<bound_mat_ordered.outerSize(); i++){
			for(SparseMatrix<double>::InnerIterator it(bound_mat_ordered, i); it; ++it){
				ordered_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
			}
		}
  ordered_bound_mat.close();

	// Nodes not on the boundary

	ofstream ordered_not_bound_mat;
	ordered_not_bound_mat.open("ordered_not_bound_mat.txt");
	if(ordered_not_bound_mat.is_open() == false){
		cout << "Unable to open ordered_not_bound_mat" << endl << std::flush;
	}
		ordered_not_bound_mat << m.num_points<< "\n";
		for(int i=0; i<not_bound_mat_ordered.outerSize(); i++){
			for(SparseMatrix<double>::InnerIterator it(not_bound_mat_ordered, i); it; ++it){
				ordered_not_bound_mat << it.col() << " " << it.row() << " " << it.value() << endl;
			}
		}
	ordered_not_bound_mat.close();
}

VectorXd w_k_builder(double(*func)(Vector2d, double), Mesh &m, double d){
	VectorXd ret(m.num_edges);
	for(int i=0; i<m.num_edges; i++){
		int index = m.bound[i];
		Vector2d point = m.points[index];
		ret[i] = func(point, d);
	}
	return ret;
}

VectorXd b_vector_builder(SparseMatrix<double> &a_ik, VectorXd w_k){
	return -1. * a_ik.transpose() * w_k;
}

VectorXd matrix_solver(SparseMatrix<double> &a_ij, VectorXd b){
	VectorXd ret;

	SimplicialLLT<SparseMatrix<double> > solver;
	solver.compute(a_ij);
	ret = solver.solve(b);
	return ret;
}

VectorXd w_stitcher(VectorXd w_ik, VectorXd w_ij, Mesh &m){
	int i=0;
	int j=0;
	int b=0;
	VectorXd ret(m.num_points);

	for(b=0; b<m.num_points; b++){
		if(b == m.bound[i]){
			ret[b] = w_ik[i];
			i++;
		} else if(b == m.not_bound[j]){
			ret[b] = w_ij[j];
			j++;
		} else {
			cout << "Error in w_switcher at point p: " << b << endl;
		}
	}
	return ret;
}

void print_w(VectorXd w, Mesh &m, fstream &w_printed){
	if(w_printed.is_open() == false){
		cout << "Unable to open " << w_printed << endl;
	}

	w_printed << m.num_points << endl;

	for(int i=0; i<m.num_points; i++){
		w_printed << m.points[i][0] << " " << m.points[i][1] << " " << w[i] << endl;
	}

	w_printed << m.num_tris << endl;

	for(int t=0; t<m.num_tris; t++){
		w_printed << m.tris[t].i << " " << m.tris[t].j << " " << m.tris[t].k << endl;
	}

	w_printed.close();
}

double energy(Mesh &m, VectorXd w, double gamma, double rho){
	int ii,jj,kk;
	double final_soln = 0.0;

	for(int t=0; t<m.num_tris; t++){
		ii = m.tris[t].i;
		jj = m.tris[t].j;
		kk = m.tris[t].k;

		Vector2d p0 = m.points[ii];
		Vector2d p1 = m.points[jj];
		Vector2d p2 = m.points[kk];

		FromRefTri ref_tri(p0, p1, p2);
		Matrix2d jacob = m.jacobian(p0, p1, p2);

		// gradient of u term

		double element = 0.0;

		double gradient = 0.0;
		double u_sq = 0.0;

		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				int i_index = m.tris[t][i];
				int j_index = m.tris[t][j];

				gradient += (m.grad_v[i].transpose() * jacob).dot(jacob.transpose() * m.grad_v[j]) * (w[i_index] * w[j_index]);

				for(int k=0; k<4; k++){
					u_sq += m.v[i][k] * m.v[j][k] * w[i_index] * w[j_index] * m.weight[k];
				}

			}
		}
		element += gamma * (sqrt(1. + gradient) - 1.0) + (rho * .5 * u_sq);
		final_soln += .5 * element * std::abs(ref_tri.jac());

		// // u sqared term
		// element = 0.0;
		// for(int k=0; k<4; k++){
		// 	double u_sq = 0;
		//
		// 	for(int i=0; i<3; i++){
		// 		for(int j=0; j<3; j++){
		// 			int i_index = m.tris[t][i];
		// 			int j_index = m.tris[t][j];
		//
		// 			u_sq += m.v[i][k] * m.v[j][k] * w[i_index] * w[j_index];
		// 		}
		// 	}
		// 	element += rho * .5 * u_sq * m.weight[k];
		// }
		//
		// final_soln += .5 * element * std::abs(ref_tri.jac());
	}

	return final_soln;
}
