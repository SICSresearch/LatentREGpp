/**
 * input.h
 *
 *  Created on: 14/04/2016
 *      Author: Milder
 */

#ifndef UTIL_INPUT_H_
#define UTIL_INPUT_H_

#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "matrix.h"

namespace latentregpp {

template<class T>
/**
 * Class for handle input data from \.csv files with a
 * matrix content, by default it is separating the data
 * by semicolon character but it can be modify through a constructor.
 * This class uses T template, so you can define the data type for input
 * the data, this should be the same data type of m matrix data type.
 * Use char data type is more efficient in memory.
 * @see matrix
 */
class input {
	private:

		std::string delimiter; /**< Delimiter to split each file line*/

	public:

		/**
		 * Default constructor for input class.
		 * Set delimiter with semicolon character.
		 */
		input();

		/**
		 * Constructor for input class.
		 * @param del a string delimiter.
		 */
		input(std::string);

		/**
		 * Destructor for input class.
		 */
		virtual ~input();

		/**
		 * Imports matrixes from a csv.
		 * @param filename a string with source file path to import.
		 * @param m a matrix instance for save the data.
		 */
		bool import_data(std::string filename, matrix<T>& m);

		/**
		 * Imports matrixes from a csv.
		 * @param filename a string with source file path to import.
		 * @param m a vector instance from std library for save the data.
		 */
		bool import_data(std::string filename, std::vector<T>& m);

		/**
		 * Gets the delimiter used for inputting.
		 */
		char get_delimiter() const;

		/**
		 * Sets the delimiter for inputting text matrices.
		 */
		void set_delimiter(char);
};

template<class T>
input<T>::input() {
	delimiter = ",; ";
}

template<class T>
input<T>::input(std::string delimiter) {
	this->delimiter = delimiter;
}

template<class T>
input<T>::~input() {
}

template<class T>
bool input<T>::import_data(std::string filename, matrix<T>& m) {
	std::string line;
	std::ifstream file(filename.c_str());
	if (file.is_open()) {
		while (std::getline(file, line)) {
			std::vector<T> splitted;
			std::string::const_iterator start = line.begin();
			std::string::const_iterator end = line.end();
			std::string::const_iterator next = std::find_first_of(start, end, delimiter.begin(), delimiter.end());
			while (next != end) {
				std::string to_add(start, next);
				if ( !to_add.empty() )
					splitted.push_back(strtold(to_add.c_str(), NULL));
				start = next + 1;
				next = std::find_first_of(start, end, delimiter.begin(), delimiter.end());
			}
			std::string to_add(start, next);
			if ( !to_add.empty() )
				splitted.push_back(strtold(to_add.c_str(), NULL));
			m.add_row(&splitted);
		}
		file.close();
		return true;
	}
	return false;
}


template<class T>
bool input<T>::import_data(std::string filename, std::vector<T>& m) {
	std::string line;
	std::ifstream file(filename.c_str());
	if (file.is_open()) {
		while (std::getline(file, line)) {
			std::vector<T> splitted;
			std::string::const_iterator start = line.begin();
			std::string::const_iterator end = line.end();
			std::string::const_iterator next = std::find_first_of(start, end, delimiter.begin(), delimiter.end());
			while (next != end) {
				std::string to_add(start, next);
				if ( !to_add.empty() )
					splitted.push_back(strtold(to_add.c_str(), NULL));
				start = next + 1;
				next = std::find_first_of(start, end, delimiter.begin(), delimiter.end());
			}
			std::string to_add(start, next);
			if ( !to_add.empty() )
				splitted.push_back(strtold(to_add.c_str(), NULL));
			m.insert(m.end(), splitted.begin(), splitted.end());
		}
		file.close();
		return true;
	}
	return false;
}

template<class T>
char input<T>::get_delimiter() const {
	return delimiter;
}

template<class T>
void input<T>::set_delimiter(char del) {
	this->delimiter = delimiter;
}

} /* namespace latentregpp */

#endif /* UTIL_INPUT_H_ */
