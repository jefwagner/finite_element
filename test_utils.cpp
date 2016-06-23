#include <Eigen/Core>

#include "utils.hpp"
#include "test.hpp"

using namespace std;
using namespace Eigen;

void test_cross(){
	Vector2d u(1.,1.), v(1.,-1.);
	print_status( cross(u,v) == -2., "Vector2d cross" );
}

void test_comp(){
	Vector2d u(1.,1.), v(1.,-1.), w(1.,1.);
	print_status( !comp_eq(u,v) && comp_eq(u,w), "Vector2d comp_eq");
}

void test_FromRefTri(){
	Vector2d a(2.,1.5), b(1.,1.), c(2.,1.);
	Vector2d p0(0.,0.), p1(1.,0.), p2(0.,1.);
	FromRefTri f(a,b,c);
	int status = comp_eq(a, f(p0));
	status &= comp_eq(b, f(p1));
	status &= comp_eq(c, f(p2));
	print_status( status, "FromRefTri functor");
}
