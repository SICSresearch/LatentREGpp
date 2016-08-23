#include "directories.h"

namespace lrpp {
	
	void setwd(std::string path) {
		prev_path = path;
	}

	std::string getwd(void){
		return prev_path;
	}
}
