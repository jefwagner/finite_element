#ifndef JW_FE_REGION
#define JW_FE_REGION

#include <fstream>
#include "triangle.hpp"
#include "mesh.hpp"

typedef struct triangulateio tio;

int num_circ(double r);
tio * malloc_tio();
tio * malloc_pslg(int num_points);
void free_tio(tio *in);
void fill_pslg(tio *in, double d, double r1, double r2);
void print_mesh( fstream &f, Mesh &m);


#endif /* JW_FE_REGION */