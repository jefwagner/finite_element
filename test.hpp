#include <iostream>
#include <string>

int print_status( int status, std::string s){
	if( status){
		std::cout << " + " << s << " Passes." << std::endl;
		return 0;
	}else{
		std::cout << " - " << s << " Fails! " << std::endl;
		return 1;
	}
}
