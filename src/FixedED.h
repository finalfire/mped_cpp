/*
 * FixedED.h
 *
 *  Created on: 22/nov/2016
 *      Author: finalfire
 */

#ifndef FIXEDED_H_
#define FIXEDED_H_

//#include <memory>
#include "Matrix.h"

template<typename T>
class FixedED {
private:
	//typedef Matrix<T> MatrixT;
	//std::shared_ptr<MatrixT> matrix;
	Matrix<T>* matrix;
public:
	FixedED(unsigned r, unsigned c) {
		//matrix = (new MatrixT(r, c, 0));
		matrix = new Matrix<T>(r, c, 0);

		// setup the matrix here
		for (size_t i = 0; i < r; ++i) (*matrix)(i,0) = i;
		for (size_t j = 0; j < c; ++j) (*matrix)(0,j) = j;
	}

	~FixedED() { delete matrix; }

	unsigned edit_distance_matching_schema_enhanced(const std::vector<unsigned>& a, const std::vector<unsigned>& b, const size_t& al, const size_t& bl,
			unsigned* sg1, unsigned* sg2, const size_t& sg1l, const size_t& sg2l, const matching_schema<bool>& m) {

		unsigned* sig1_index = new unsigned[sg1l]; for (size_t i = 0; i < sg1l; ++i) { sig1_index[sg1[i]] = i; }
		unsigned* sig2_index = new unsigned[sg2l]; for (size_t i = 0; i < sg2l; ++i) { sig2_index[sg2[i]] = i; }

		for (size_t i = 1; i < al+1; ++i)
			for (size_t j = 1; j < bl+1; ++j) {
				(*matrix)(i,j) = min(
						(*matrix)(i-1,j) + 1,																		// deletion
						(*matrix)(i,j-1) + 1,																		// insertion
						(*matrix)(i-1,j-1) + 1 * m.ms[sig1_index[a[i-1]]][sig2_index[b[j-1]]]		// if in the matching schema there's a false, they match
				);
			}

		return (*matrix)(al,bl);
	}
};

#endif /* FIXEDED_H_ */
