#include "region.hpp"
#include "test.hpp"
#include "triangle.hpp"
#include <sstream>
#include <string>
#include <iomanip>

double unitary(Vector2d x){
	return 1.0;
}

namespace{
	double func(Vector2d point, double d){
		if(nearbyint(sqrt(((point[0] + d/2.) * (point[0] + d/2.)) + (point[1] * point[1]))) == 1){
			return 0.1;
		} else if(nearbyint(sqrt(((point[0] - d/2.) * (point[0] - d/2.)) + (point[1] * point[1]))) == 1){
			return 0.1;
		}
		return 0.0;
	}

}

//-----------------------------------------------
// T e s t mesh class on a larger scale. 30 nodes
// 	vs couple thousand nodes.
void test_region(){
	tio *in;
	tio *out;
	tio *vor;

	// Initial Conditions
	double r = 5.;
	double gamma = 1.;
	double d = 4.;
	double dd = .01;
	double rho = 0.5;

	// Storing I.C.
	Mesh_Struct force(r, gamma, rho, d, dd);

	// Mesh Construction
	int num_points = 56 + num_circ(1.0) + num_circ(1.0);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, 1.0, 1.0);
	out = malloc_tio();
	vor = (tio *) NULL;
	char *triswitches = "pq30a0.1z";
	triangulate(triswitches, in, out, vor);
	Mesh m(out);

	// Print for plotting purposes
	fstream f("region_mesh.txt", fstream::out);
	print_mesh( f, m);

	// Reorder and Analyze
	int center = 12;

	SparseMatrix<double> not_bound_mat_ordered(m.num_points-m.num_edges, m.num_points-m.num_edges);
	SparseMatrix<double> bound_mat_ordered(m.num_edges, m.num_points-m.num_edges);

	m.reorder_nodes(center);

	m.mass_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);
	m.stiffness_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);

	// Build bound_mat_ordered nad not_bound_mat_ordered.

	print_mat(m, not_bound_mat_ordered, bound_mat_ordered);

	// Solving for w.

	VectorXd w_k, w_ij, w;
	w_k = w_k_builder(func, m, d);
	VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
	w_ij = matrix_solver(not_bound_mat_ordered, b);
	w = w_stitcher(w_k, w_ij, m);
	fstream f_w("w_vector.txt", fstream::out);

	print_w(w, m, f_w);

	print_status( 1==1, "Region");
}

//-----------------------------------------------
// T e s t energy calculation on a larger scale
void test_energy(){
	tio *in;
	tio *out;
	tio *vor;

	// Initial Conditions
	double r = 5.;
	double gamma = 1.;
	double d = 4.;
	double dd = .01;
	double rho = 0.5;

	// Storing I.C.
	Mesh_Struct force(r, gamma, rho, d, dd);

	// Build mesh
	int num_points = 56 + num_circ(1.0) + num_circ(1.0);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, 1.0, 1.0);
	out = malloc_tio();
	vor = (tio *) NULL;
	char *triswitches = "pq30a0.1z";
	triangulate(triswitches, in, out, vor);
	Mesh m(out);

	// Reorder and Analyze
	int center = 12;

	SparseMatrix<double> not_bound_mat_ordered(m.num_points-m.num_edges, m.num_points-m.num_edges);
	SparseMatrix<double> bound_mat_ordered(m.num_edges, m.num_points-m.num_edges);

	m.reorder_nodes(center);

	m.mass_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);
	m.stiffness_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);

	VectorXd w_k, w_ij, w;
	w_k = w_k_builder(func, m, 1.0);
	VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
	w_ij = matrix_solver(not_bound_mat_ordered, b);
	w = w_stitcher(w_k, w_ij, m);

	// Solve energy for the mesh
	energy(m, w, gamma, rho);

	print_status( 1==1, "energy");
}

//-----------------------------------------------
// T e s t finding the force from the slope of
// 	the energy vs distance plot.
void test_force(){
	tio *in;
	tio *out;
	tio *vor;

	// Initial Conditions
	double r = 5.;
	double gamma = 1.;
	double rho = 0.;
	double d = 2.1;
	double dd = .01;

	// Storing I.C.
	fstream energy_f("Energy_plot.txt", fstream::out);

	// Creating multiple meshes at several distances
	// 	-Storing the mesh
	// 	-Storing the Energy
	while(d<=10.0){
		Mesh_Struct force(r, gamma, rho, d, dd);

		int num_points = 56 + num_circ(1.0) + num_circ(1.0);

		in = malloc_pslg(num_points);
		fill_pslg(in, d, 1.0, 1.0);
		out = malloc_tio();
		vor = (tio *) NULL;
		char *triswitches = "pq30a0.1z";
		triangulate(triswitches, in, out, vor);
		Mesh m(out);

		int center = 12;

		SparseMatrix<double> not_bound_mat_ordered(m.num_points-m.num_edges, m.num_points-m.num_edges);
		SparseMatrix<double> bound_mat_ordered(m.num_edges, m.num_points-m.num_edges);

		m.reorder_nodes(center);

		m.mass_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);
		m.stiffness_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);

		VectorXd w_k, w_ij, w;
		w_k = w_k_builder(func, m, d);
		VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
		w_ij = matrix_solver(not_bound_mat_ordered, b);
		w = w_stitcher(w_k, w_ij, m);

		stringstream strs;
		strs << "data/w_vector_d_" << d << "_rho_" << rho << "_gamma_" << gamma << ".txt";
		string str = strs.str();
		cout << str << endl;
		fstream f_w(str.c_str(), fstream::out);
		print_w(w, m, f_w);

		double energy_soln = energy(m, w, gamma, rho);

		energy_f << d << " " << energy_soln << endl;

		d += dd;
	}

	energy_f.close();
}
