#include <fstream>
#include <cmath>
#include <cstdlib>

#include "mesh.hpp"
#include "triangle.hpp"

#define PI 3.14159265358979323846E+0

using namespace std;

typedef struct triangulateio tio;

namespace{

	void add_point(double *pl, int *n, double x, double y){
		pl[2*(*n)+0] = x;
		pl[2*(*n)+1] = y;
		(*n)++;
	}

	void add_point_elipse(double *pl, int *n, double a, double b, int size){
		double theta = 2*PI * ((*n)/size);
		double x = a*cos(theta);
		double y = b*sin(theta);
		pl[2*(*n)+0] = x;
		pl[2*(*n)+1] = y;
		(*n)++;
	}

	void free_if_not_null( void *ptr){
		if( ptr != NULL){
			free(ptr);
		}
	}
}


int num_circ(double r){
	double dx = 0.01;
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
	// Four corners of surrounding box

	// // Making the box into an elipse, the code for the box is still there
	// // size indicates resolution of elipse, higher size, higher resolution
	// int size = 4;

	first = n;

	// for(i=first; i<size; i++){
	// 	add_point_elipse(pl, &n, d, 1.5 * d, size);
	// }

	add_point(pl, &n, 1.5*d, d);
	add_point(pl, &n, 1.5*d, -d);
	add_point(pl, &n, -1.5*d, -d);
	add_point(pl, &n, -1.5*d, d);
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
