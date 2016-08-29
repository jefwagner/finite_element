#ifndef JW_FE_REGION
#define JW_FE_REGION

#include <fstream>
#include "mesh.hpp"
#include "triangle.hpp"

typedef struct triangulateio tio;

class Mesh_Struct{
public:
  // The Base dimensionless variables
  double length_unit;
  double acceleration_unit;
  double energy_density_unit;

  // The terms that are manually adjusted by user in main file.
  double d;
  double dd;
  double new_d;

  // The dimensionless variables that are dependent on those above.
  double energy_unit;
  double time_unit;
  double mass_unit;
  double force_unit;
  double rho;

  // Dimensionless variables constructor. This stores everything specific to
  //  the surface. It stores the size of the system and the location of the holes
  Mesh_Struct(double r, double gamma, double density, double dist, double change_dist){
    length_unit = r * .000001;                // In terms of micrometers
    acceleration_unit = 9.8;
    energy_density_unit = gamma * .0728;      // In terms of percent surface tension of water (N/m)

    d = dist;
    dd = change_dist;
    new_d = d + dd;

    energy_unit = energy_density_unit * length_unit * length_unit;
    time_unit = sqrt(r / acceleration_unit);
    mass_unit = acceleration_unit * energy_density_unit * length_unit;
    force_unit = energy_density_unit * length_unit;

    rho = density * mass_unit / (length_unit * length_unit * length_unit);    // In terms of density of water
  }
};

// Functions necessary for building the surface and feeding it to Triangle
int num_circ(double r);
tio * malloc_tio();
tio * malloc_pslg(int num_points);
void free_tio(tio *in);
void fill_pslg(tio *in, double d, double r1, double r2);
void print_mesh( fstream &f, Mesh &m);

// Plot SparseMatrix
void print_mat(Mesh &, SparseMatrix<double> &, SparseMatrix<double> &);

// Added to evaluate Ax = b
VectorXd w_k_builder(double(*func)(Vector2d, double), Mesh &, double);
VectorXd b_vector_builder(SparseMatrix<double> &, VectorXd);
VectorXd matrix_solver(SparseMatrix<double> &, VectorXd);
VectorXd w_stitcher(VectorXd, VectorXd, Mesh &);
void print_w(VectorXd, Mesh &, fstream &w_printed);
double energy(Mesh &, VectorXd, double, double);

#endif /* JW_FE_REGION */
