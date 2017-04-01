#ifndef MATCHINGSCHEMA_H_
#define MATCHINGSCHEMA_H_

#include "Utility.h"

template <typename T>
struct matching_schema {
	size_t n, m;
	size_t p1, p2;
	T** ms;
	bool** chi;

	// x is the default boolean value for the matrix
	matching_schema(size_t a, size_t b, size_t c, size_t d, const T x, bool def) : n(a), m(b), p1(c), p2(d) {
		ms = new T*[n];
		for (size_t i = 0; i < n; ++i) ms[i] = new T[m];

		chi = new bool*[n];
		for (size_t i = 0; i < n; ++i) chi[i] = new bool[m];

		init(x, def);
	}

	~matching_schema() {
		for (size_t i = 0; i < n; ++i)
			delete[] ms[i];
		delete[] ms;

		for (size_t i = 0; i < n; ++i)
			delete[] chi[i];
		delete[] chi;
	}

	void init(T b, bool d) {
		// setup the matching schema
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				ms[i][j] = b;
		// setup the constraints
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				chi[i][j] = d;
	}

	// e.g., if value == false, ms[i][j] == false means that i and j are in match
	void set_identity(std::string& sigma1, std::string& sigma2, T value) {
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				if (sigma1[i] == sigma2[j])
					ms[i][j] = value;
	}

	// e.g., if value == false, ms[i][j] == false means that i and j are in match
	void set_general(const std::string& sigma1, const std::string& sigma2, T value) {
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (p1 !=0 && p2 !=0 && (i/p2 == j/p1))
					ms[i][j] = value;
				else
					ms[i][j] = !value;
	}

	void set_general_fibre(const std::vector<std::string>& sigma1, const std::vector<std::string>& sigma2, T value) {
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (p1 !=0 && p2 !=0 && (i/p2 == j/p1))
					ms[i][j] = value;
				else
					ms[i][j] = !value;
	}

	// TODO: CONSTRAINTS
	/*void set_constraints(std::map<char, int>& sig1, std::map<char, int>& sig2, std::vector<p_constr>& v, bool mode) {
		for (size_t i = 0; i < v.size(); ++i)
			ms[sig1[v[i].first]][sig2[v[i].second]] = mode;
	}*/

	void set_constraint(size_t i, size_t j, bool t) { chi[i][j] = t; }

	void print_matching_schema(const std::string& sigma1, const std::string& sigma2) {
		std::cout << "\t";
		for (int j = 0; j < m; ++j)
			std::cout << sigma2[j] << "\t";
		std::cout << std::endl;
		for (int i = 0; i < n; ++i) {
			std::cout << sigma1[i] << "\t";
			for (int j = 0; j < m; ++j)
				std::cout << ms[i][j] << "\t";
			std::cout << std::endl;
		}
	}

	void shift(int startRow, int startCol, T value){


		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				ms[i][j] = true;

		for (int i = (startRow); i < n; ++i)
			for (int j = (startCol); j < m; ++j)
				if (p1 !=0 && p2 !=0 && ((i-(startRow))/p2 == (j-(startCol))/p1))
					ms[i][j] = value;

	
		/*for (int i=0; i<n; i++){
			for (int j=0; j<m; j++){
				std::cout<<ms[i][j]<<" | ";
			}
			std::cout<<"\n";
		}
		std::cout<<"\n";*/
	}
};

#endif /* MATCHINGSCHEMA_H_ */
