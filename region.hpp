#ifndef JW_FE_REGION
#define JW_FE_REGION

#include <fstream>
#include "mesh.hpp"
#include "triangle.hpp"

typedef struct triangulateio tio;

class Mesh_Struct{
public:
  double length_unit;
  double acceleration_unit;
  double surface_tension;
  double energy_unit;
  double time_unit;
  double mass_unit;
  double force_unit;
  double d;
  double dd;
  double new_d;
  double density;

  Mesh_Struct(double r, double gamma, double g, double dist, double change_dist){
    length_unit = r;
    acceleration_unit = g;
    surface_tension = gamma;
    energy_unit = gamma * r * r;
    time_unit = sqrt(r / g);
    mass_unit = g * gamma * r;
    force_unit = gamma * r;
    d = dist/r;
    dd = change_dist;
    new_d = d + dd;
    density = mass_unit / length_unit;
  }
};

int num_circ(double r);
tio * malloc_tio();
tio * malloc_pslg(int num_points);
void free_tio(tio *in);
void fill_pslg(tio *in, double d, double r1, double r2);
void print_mesh( fstream &f, Mesh &m);

// Added for SparseMatrix

void print_mat(Mesh &, SparseMatrix<double> &, SparseMatrix<double> &);

// Added to evaluate Ax = b

VectorXd w_k_builder(double(*func)(Vector2d, double), Mesh &, double);
VectorXd b_vector_builder(SparseMatrix<double> &, VectorXd);
VectorXd matrix_solver(SparseMatrix<double> &, VectorXd);
VectorXd w_stitcher(VectorXd, VectorXd, Mesh &);
void print_w(VectorXd, Mesh &, fstream &w_printed);
double energy(double, double, Mesh &, VectorXd);

#endif /* JW_FE_REGION */
