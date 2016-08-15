#include "region.hpp"
#include "triangle.hpp"
#include "utils.hpp"
#include "mesh.hpp"
#include "test.hpp"

double unitary(Vector2d x){
	return 1.0;
}

namespace{
	double func(Vector2d point, double d){
		if(round(sqrt(((point[0] + d/2.) * (point[0] + d/2.)) + (point[1] * point[1]))) == 1){
			return -0.1;
		} else if(round(sqrt(((point[0] - d/2.) * (point[0] - d/2.)) + (point[1] * point[1]))) == 1){
			return -0.1;
		}
		return 0.0;
	}

}

int main(){
  tio *in;
  tio *out;
  tio *vor;

  double r = 5.;
  double gamma = 1.;
  double rho = 0.;
  double d = 2.1;
  double dd = .01;

  fstream energy_f("Energy_plot.txt", fstream::out);

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
    w_k = w_k_builder(func, m, 1.0);
    VectorXd b = b_vector_builder(bound_mat_ordered, w_k);
    w_ij = matrix_solver(not_bound_mat_ordered, b);
    w = w_stitcher(w_k, w_ij, m);

    double energy_soln = energy(m, w, gamma, rho);

    energy_f << d << " " << energy_soln << endl;

    d += dd;
  }

  energy_f.close();

  return 0;
}