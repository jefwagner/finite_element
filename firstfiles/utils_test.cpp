#include <iostream>
#include <string>

#include <Eigen/Core>

#include 'utils.hpp'

int print_status( int status, string s){
	if( status){
		cout << " + " << s << " Passes." << endl;
		return 0;
	}else{
		cout << " - " << s << " Fails! " << endl;
		return 1;
	}
}

void test_cross(){
	Vector2d u(1.,1.), v(1.,-1.);
	print_status( cross(u,v) == -2., "Vector2d cross" );
}

void test_FromRefTri(){
	Vector2d a(2.,1.5), b(1.,1.), c(2.,1.);
	Vector2d p0(0.,0.), p1(1.,0.), p2(0.,1.);
	FromRefTri f(a,b,c);
	int status = ( a.array() == f(p0).array() );
	status &= ( b.array() == f(p1).array() );
	status &= ( c.array() == f(p2).array() );
	print_status( status, "FromRefTri functor");
}

int main(){
	cout << endl;
	cout << "Testing utilility functions in utils.hpp" << endl;
	cout << "----------------------------------------" << endl;
	test_cross();
	test_FromRefTri();
	cout << "----------------------------------------" << endl;
	cout << endl;
}