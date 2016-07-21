#include "region.hpp"
#include "test.hpp"
#include "triangle.hpp"

double unitary(Vector2d x){
	return 1.0;
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


	print_mat(m, center);

	// print_stiff_mat(m, center);

	print_status( 1==1, "Region");
}
