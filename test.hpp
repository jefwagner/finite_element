#ifndef JW_FE_TEST
#define JW_FE_TEST

#include <iostream>
#include <string>

static int print_status( int status, const char *s){
	if( status){
		std::cout << " + " << s << " Passes." << std::endl;
		return 0;
	}else{
		std::cout << " - " << s << " Fails! " << std::endl;
		return 1;
	}
}

#endif /* JW_FE_TEST */
