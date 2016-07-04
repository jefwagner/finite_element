#include "region.hpp"
#include "test.hpp"
#include "triangle.hpp"

void test_region(){
	tio *in;
	tio *out;
	tio *vor;

	double r1 = 1;
	double r2 = 1.2;
	double d = 2.5;

	int num_points = 4 + num_circ(1.) + num_circ(1.2);

	in = malloc_pslg(num_points);
	fill_pslg(in, d, r1, r2);
	out = malloc_tio();
	vor = (tio *) NULL;
	char *triswitches = "pq30a0.1z";
	triangulate(triswitches, in, out, vor);
	Mesh m(out);

	fstream f("region_mesh.txt", fstream::out);
	print_mesh( f, m);

	print_status( 1==1, "Region");
}
