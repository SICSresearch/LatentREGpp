/**
 * matrix.h
 *
 *  Created on: 14/04/2016
 *      Author: Milder
 */

#ifndef UTIL_MATRIX_H_
#define UTIL_MATRIX_H_

#include <vector>
#include <iostream>
#include <fstream>

namespace latentregpp {

template<class T>
/**
 * Matrix class used to store information.
 */
class matrix {
	private:

	public:
		std::vector<std::vector<T> > data;
		matrix();

		/**
		 * Creates a new matrix with size rows x columns.
		 * @param rows the number of rows in matrix.
		 * @param cols the number of columns in matrix.
		 */
		matrix(int rows, int cols);

		/**
		 * Creates a new matrix specifying only the number of rows,
		 * rows can be added after using add_row(std::vector) method.
		 * @param rows a integer with rows number.
		 * @see add_row()
		 */
		matrix(int rows);

		/**
		 * Adds a new row to the current matrix.
		 * @param new_row a pointer from std vector library, this will be the new row in matrix.
		 */
		void add_row(std::vector<T> *new_row);

		/**
		 * Adds a new row to the current matrix.
		 * @param new_row a std vector library, this will be the new row in matrix.
		 */
		void add_row(std::vector<T> new_row);

		/**
		 * Adds a new row to the current matrix setting size.
		 * @param new_row a T template for new_row.
		 * @param size the size for row.
		 */
		void add_row(T* new_row, int size);

		/**
		 * Adds an empty row to matrix by default with std vector library and T template.
		 */
		void add_row();

		/**
		 * Adds a empty row with specified size.
		 * @param size the size for new row.
		 */
		void add_row(int);

		/**
		 * Clean the number of the specific row, by default this row is set with zeros.
		 * WARNING: be sure about the index row. No throw exception implemented for index array out of bounds.
		 * @param idx the index of row.
		 */
		void reset_row(int);

		/**
		 * Clean all matrix with zeros by default.
		 */
		void reset();

		/**
		 * Adds an element to the last row.
		 * @param e the T element.
		 */
		void add_element(T e);

		/**
		 * Adds an element to a specific row.
		 * @param row the number of row to add element.
		 * @param e the T element to add.
		 */
		void add_element(int row, T e);

		/**
		 * Destructor for matrix class.
		 */
		virtual ~matrix();

		/**
		 * Getter method, returns the data.
		 * @return the data. it is a vector of vectors from std library
		 * */
		std::vector<std::vector<T> > get_data();

		/**
		 * Returns the i-th row in the matrix.
		 * WARNING: be sure about the index row. No throw exception implemented for index array out of bounds.
		 * @param i the index of row to show.
		 * @return std vector according the number of row
		 * */
		std::vector<T> get_row(int i);

		/**
		 * Returns the i-th row pointer in the matrix.
		 * WARNING: be sure about the index row. No throw exception implemented for index array out of bounds.
		 * @param i the index of row to show.
		 * @return std vector according the number of row
		 * */
		std::vector<T>* get_pointer_row(int i);

		/**
		 * Returns number of rows of the matrix.
		 * @return an integer with number of row in matrix
		 * */
		int rows();

		/**
		 * Returns number of columns of a specific row.
		 * WARNING: be sure about the index row. No throw exception implemented for index array out of bounds.
		 * @param row the index number of row matrix.
		 * @return an integer with number of columns in ro.
		 * */
		int columns(int row);

		/**
		 * Exports the matrix to a .csv file
		 * @param filename where the matrix will be saved
		 * @param sep separator of .csv file
		 * */
		void export_to_csv (std::string, char sep = ';');

		/**
		 * Accessing operator for an element.
		 * Overload parenthesis operator.
		 */
		inline T &operator()(const int row, const int col);
};

template<class T>
matrix<T>::matrix(int rows) {
	data = std::vector<std::vector<T> >(rows);
}

template<class T>
matrix<T>::matrix() {
}

template<class T>
matrix<T>::~matrix() {
}

template<class T>
matrix<T>::matrix(int rows, int cols) {
	data = std::vector<std::vector<T> >(rows, std::vector<T>(cols));
}

template<class T>
void matrix<T>::add_row(std::vector<T> *new_row) {
	data.push_back(*new_row);
}

template<class T>
void matrix<T>::add_row(std::vector<T> new_row) {
	data.push_back(new_row);
}

template<class T>
void matrix<T>::add_row(T* new_row, int size) {
	data.push_back(std::vector<T>(new_row, new_row + size));
}

template<class T>
void matrix<T>::add_row() {
	data.push_back(std::vector<T>());
}

template<class T>
void matrix<T>::reset() {
	for ( int i = 0; i < data.size(); ++i )
		for ( int j = 0; j < data[i].size(); ++j )
			data[i][j] = 0;
}


template<class T>
void matrix<T>::add_row(int size) {
	data.push_back(std::vector<T>(size));
}

template<class T>
void matrix<T>::reset_row(int idx) {
	for ( unsigned int i = 0; i < data[idx].size(); ++i )
		data[idx][i] = 0;
}

template<class T>
void matrix<T>::add_element(T e) {
	data.back().push_back(e);
}

template<class T>
void matrix<T>::add_element(int row, T e) {
	data[row].push_back(e);
}

template<class T>
int matrix<T>::rows() {
	return data.size();
}

template<class T>
int matrix<T>::columns(int row) {
	if ( (unsigned)row >= data.size() )
		return -1;
	return data[row].size();
}

template<class T>
std::vector<std::vector<T> > matrix<T>::get_data() {
	return data;
}

template<class T>
std::vector<T> matrix<T>::get_row(int i) {
	return data[i];
}

template<class T>
std::vector<T>* matrix<T>::get_pointer_row(int i) {
	return &data[i];
}

template<class T>
inline T &matrix<T>::operator()(const int r, const int c) {
	return data.at(r).at(c);
}

template<class T>
std::ostream& operator<<(std::ostream &out, matrix<T> &m) {
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.columns(i); j++)
			out << (double)m(i, j) << '\t';
		out << '\n';
	}
	return (out);
}

template<class T>
void matrix<T>::export_to_csv(std::string filename, char sep) {
	std::ofstream out(filename);
	for ( auto row : data ) {
		bool first = true;
		for ( auto col : row ) {
			if ( !first ) out << sep;
			out << col;
			first = false;
		}
		out << '\n';
	}
	out.close();
}

} /* namespace latentregpp */

#endif /* UTIL_MATRIX_H_ */
