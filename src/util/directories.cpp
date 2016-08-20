#include "directories.h"

namespace irtpp {
	
	void setwd(std::string path) {
		prev_path = path;
	}

	std::string getwd(void){
		return prev_path;
	}
}
