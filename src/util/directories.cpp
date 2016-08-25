#include "directories.h"

namespace latentregpp {
	
	void setwd(std::string path) {
		prev_path = path;
	}

	std::string getwd(void){
		return prev_path;
	}
}
