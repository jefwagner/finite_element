#include "region.hpp"
#include "test.hpp"
#include "triangle.hpp"

double unitary(Vector2d x){
	return 1.0;
}

namespace{
	double func(Vector2d point, double d){
		if(round(sqrt(((point[0] + d/2.) * (point[0] + d/2.)) + (point[1] * point[1]))) == 1){
			return 1.0;
		} else if(round(sqrt(((point[0] - d/2.) * (point[0] - d/2.)) + (point[1] * point[1]))) == 1){
			return 1.0;
		}
		return 0.0;
	}

}

void test_region(){
	tio *in;
	tio *out;
	tio *vor;

	double r = 5.;
	double gamma = 1.;
	double d = 4.;
	double dd = .01;
	double rho = 0.5;

	Mesh_Struct force(r, gamma, rho, d, dd);

	int num_points = 56 + num_circ(r) + num_circ(r);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, r, r);
	out = malloc_tio();
	vor = (tio *) NULL;
	char *triswitches = "pq30a0.1z";
	triangulate(triswitches, in, out, vor);
	Mesh m(out);

	fstream f("region_mesh.txt", fstream::out);
	print_mesh( f, m);

	int center = 12;

	SparseMatrix<double> not_bound_mat_ordered(m.num_points-m.num_edges, m.num_points-m.num_edges);
	SparseMatrix<double> bound_mat_ordered(m.num_edges, m.num_points-m.num_edges);

	m.reorder_nodes(center);

	m.mass_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);
	cout << "Mass Matrix done!!" << endl;
	m.stiffness_matrix(unitary, bound_mat_ordered, not_bound_mat_ordered);
	cout << "Stiffness Matrix done!!" << endl;

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

void test_energy(){
	tio *in;
	tio *out;
	tio *vor;

	double r = 5.;
	double gamma = 1.;
	double d = 4.;
	double dd = .01;
	double rho = 0.5;

	Mesh_Struct force(r, gamma, rho, d, dd);

	int num_points = 56 + num_circ(r) + num_circ(r);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, r, r);
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
	w_k = w_k_builder(func, m, r);
	VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
	w_ij = matrix_solver(not_bound_mat_ordered, b);
	w = w_stitcher(w_k, w_ij, m);

	energy(m, w, gamma, rho);

	print_status( 1==1, "energy");
}

void test_force(){
	tio *in;
	tio *out;
	tio *vor;

	double r = 5.;
	double gamma = 1.;
	double rho = 0.5;
	double d = 2.2;
	double dd = .1;

	fstream energy_f("Energy_plot.txt", fstream::out);

	while(d<=5.0){
		Mesh_Struct force(r, gamma, rho, d, dd);

		int num_points = 56 + num_circ(r) + num_circ(r);

		in = malloc_pslg(num_points);
		fill_pslg(in, d, r, r);
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
		w_k = w_k_builder(func, m, r);
		VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
		w_ij = matrix_solver(not_bound_mat_ordered, b);
		w = w_stitcher(w_k, w_ij, m);

		double energy_soln = energy(m, w, gamma, rho);

		energy_f << d << " " << energy_soln << endl;

		d += dd;
	}

	energy_f.close();
}
