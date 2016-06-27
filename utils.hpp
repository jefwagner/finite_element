#ifndef JW_FE_UTILS
#define JW_FE_UTILS

#include <Eigen/Core> // The Vector2d data type

using namespace Eigen;

/*
 * Cross product for Vector2d types
 */
double cross(Vector2d u, Vector2d v){
	return u(0)*v(1)-u(1)*v(0);
}

/*
 * Compare Vector2d objects
 */
int comp_eq( Vector2d u, Vector2d v){
	return (u(0)==v(0) && u(1)==v(1));
}

/*
 * FromRefTri functor.
 *
 * This functor defines the transformation from the
 * nu = (xi,eta) space of the reference triangle to
 * the x,y space of the triangle given by the three
 * points a,b,c.
 *
 * (0,0) should map to point a
 * (1,0) should map to point b
 * (0,1) should map to point c
 */
class FromRefTri{
	Vector2d a, b, c;

public:
	FromRefTri(Vector2d p0, Vector2d p1, Vector2d p2){
		a = p0;
		b = p1;
		c = p2;
	}
	double jac(){
		return cross(a,b)+cross(b,c)+cross(c,a);
	}
	Vector2d operator()(Vector2d nu){
		return (b-a)*nu(0) + (c-a)*nu(1) + a;
	}
};

#endif /*JW_FE_UTILS */
