/**
 * directories.h
 *
 *  Created on: 19/08/2016
 *      Author: sics
 */

#ifndef UTIL_DIRECTORIES_H_
#define UTIL_DIRECTORIES_H_

#include <string>

namespace latentregpp {
	
	static std::string prev_path = "";
	
	void setwd(std::string path);
	std::string getwd(void);

}

#endif /* UTIL_DIRECTORIES_H_ */
