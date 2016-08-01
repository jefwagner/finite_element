#ifndef JW_FE_REGION
#define JW_FE_REGION

#include <fstream>
#include "mesh.hpp"
#include "triangle.hpp"

typedef struct triangulateio tio;

int num_circ(double r);
tio * malloc_tio();
tio * malloc_pslg(int num_points);
void free_tio(tio *in);
void fill_pslg(tio *in, double d, double r1, double r2);
void print_mesh( fstream &f, Mesh &m);

// Added for SparseMatrix

void print_mat(Mesh &, SparseMatrix<double> &, SparseMatrix<double> &);

// Added to evaluate Ax = b

VectorXd w_k_builder(double(*func)(Vector2d), Mesh &);
VectorXd b_vector_builder(SparseMatrix<double> &, VectorXd);
VectorXd matrix_solver(SparseMatrix<double> &, VectorXd);
VectorXd w_stitcher(VectorXd, VectorXd, Mesh &);
void print_w(VectorXd, Mesh &);




#endif /* JW_FE_REGION */
