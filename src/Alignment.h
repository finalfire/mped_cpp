#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <vector>

template <typename T>
struct Alignment {
	int distance;

	/* REMEMBER
	 * The end of the alignment for the vectors a and b is
	 * defined by -2 in the i-th index of the vectors.
	 * They contain only one occurrence of -2.
	 * */
	std::vector<T> a;
	std::vector<T> b;

	Alignment(unsigned d, size_t res) : distance(d) {
		a.reserve(res);
		b.reserve(res);
	}
};

#endif /* ALIGNMENT_H_ */
