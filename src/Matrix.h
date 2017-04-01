/*
 * Matrix.h
 *
 *  Created on: 12/ott/2015
 *      Author: finalfire
 */

#ifndef MATRIX_H_
#define MATRIX_H_

template<typename T>
class Matrix {
private:
	T* m;
	unsigned rows, columns;

public:
	Matrix(unsigned r, unsigned c) : rows(r), columns(c) { m = new T[this->rows * this->columns]; }
	Matrix(unsigned r, unsigned c, T d) : rows(r), columns(c) {
		m = new T[this->rows * this->columns];

		for (size_t i = 0; i < this->rows; i++)
			for (size_t j = 0; j < this->columns; j++)
				(*this)(i,j) = d;
	}
	~Matrix() { delete[] (T*) m; }

	T& operator()(const unsigned& i, const unsigned& j) { return this->m[i*this->columns + j]; }				// set
	const T& operator()(const unsigned& i, const unsigned& j) const { return this->m[i*this->columns + j]; }	// access

	const unsigned r() const { return this->rows; }
	const unsigned c() const { return this->columns; }
};



#endif /* MATRIX_H_ */
