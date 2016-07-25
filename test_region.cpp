#include "region.hpp"
#include "test.hpp"
#include "triangle.hpp"

double unitary(Vector2d x){
	return 1.0;
}

namespace{
	double func(Vector2d point){
		return (1-point[0]-point[1]);
	}
}

void test_region(){
	tio *in;
	tio *out;
	tio *vor;

	double r1 = 1;
	double r2 = 1.2;
	double d = 2.5;

	int num_points = 56 + num_circ(1.) + num_circ(1.2);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, r1, r2);
	out = malloc_tio();
	vor = (tio *) NULL;
	char *triswitches = "pq30a0.1z";
	triangulate(triswitches, in, out, vor);
	Mesh m(out);

	// SparseMatrix<double, RowMajor, int> sp1(m.num_points, m.num_points);
	// SparseMatrix<double, RowMajor, int> sp2(m.num_points, m.num_points);

	fstream f("region_mesh.txt", fstream::out);
	print_mesh( f, m);

	int center = 12;
	// std::cout << "Start of Unorderd Mat" << std::endl;
	// std::cout << std::endl;
	// m.mass_matrix(unitary, sp1);
	// std::cout << std::endl;
	// std::cout << "end of Unorderd Mat" << std::endl;
	//
	// m.reorder_nodes(0);
	//
	// std::cout << "Start of orderd Mat" << std::endl;
	// std::cout << std::endl;
	// m.mass_matrix(unitary, sp2);
	// std::cout << std::endl;
	// std::cout << "end of orderd Mat" << std::endl;


	SparseMatrix<double> not_bound_mat_ordered(m.num_points-m.num_edges, m.num_points-m.num_edges);
	SparseMatrix<double> bound_mat_ordered(m.num_edges, m.num_points-m.num_edges);

	m.reorder_nodes(center);

	m.mass_matrix(func, bound_mat_ordered, not_bound_mat_ordered);
	cout << "Mass Matrix done!!" << endl;
	m.stiffness_matrix(func, bound_mat_ordered, not_bound_mat_ordered);
	cout << "Stiffness Matrix done!!" << endl;

	print_mat(m, not_bound_mat_ordered, bound_mat_ordered);

	VectorXd w_k, w_ij, w;
	w_k = w_k_builder(func, m);
	VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
	w_ij = matrix_solver(not_bound_mat_ordered, b);
	w = w_stitcher(w_k, w_ij, m);

	// print_stiff_mat(m, center);

	print_status( 1==1, "Region");
}
