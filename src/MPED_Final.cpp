#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include "Alignment.h"
#include "EditDistance.h"
#include "FixedED.h"
#include "MatchingSchema.h"
#include "Utility.h"

/* Definitions */
#define endl '\n'

/* Consts */
const unsigned short _ASCII_LEN = 255 - 0;
const bool _DEBUG = false;
const std::string _HC_ARG("hc");
const std::string _BRUTEFORCE_ARG("ex");
const std::string _SPECIFIC_PMS("specific-permutations");
const std::string _SPECIFIC_MMS("specific-matrix");

/* Utility functions */
void extract_sigma(std::string&, std::string&);
void define_mapping(std::string&, std::map<char,int>&);
void extract_sigma_fibre(std::vector<std::string>&, std::vector<std::string>&);
void define_mapping_fibre(std::vector<std::string>&, std::map<std::string,int>&);
void print_alignment(const Alignment<int>&, const std::string&, const std::string&, const matching_schema<bool>&, int&, bool);

/* Solvers */
int hill_climbing(const std::vector<unsigned>&, const std::vector<unsigned>&, const size_t&, const size_t&,
		const std::vector<unsigned>&, const std::vector<unsigned>&, const size_t&, const size_t&, const size_t&, matching_schema<bool>&, edit_distance&);
int bruteforce(const std::vector<unsigned>&, const std::vector<unsigned>&, const size_t&, const size_t&,
		const std::vector<unsigned>&, const std::vector<unsigned>&, const size_t&, const size_t&, const matching_schema<bool>&, edit_distance&);

int main(int argc, char *argv[]) {
	std::ios_base::sync_with_stdio(false);

	// arguments: [hc|ex] p1 p2 [specific-permutations|specific-matrix]
	// heuristic
	std::string heuristic(argv[1]);
	// p1 and p2
	size_t p1 = argc > 3 ? fast_atoi(argv[2]) : 1;
	size_t p2 = argc > 3 ? fast_atoi(argv[3]) : 1;
	// specific matching schema
	bool specific_perm = false;
	bool specific_matrix = false;
	if (argc > 4) {
		std::string specific_temp(argv[4]);
		// if we have in input a specific ms with permutations
		if (specific_temp == _SPECIFIC_PMS)
			specific_perm = true;
		// if we have in input a specific ms with the whole matrix
		if (specific_temp == _SPECIFIC_MMS)
			specific_matrix = true;
	}


	// TODO: implement constraints managing
	bool has_constraints = false;
	bool default_constraints_mode = false;


	/* (1) read strings, extract sigmas and constraints*/
	std::string s1, s2;
	std::string sigma1(""), sigma2("");

	read_stdin(s1, s2);
	extract_sigma(s1, sigma1);
	extract_sigma(s2, sigma2);

	size_t sigma1l = sigma1.size();
	size_t sigma2l = sigma2.size();
	size_t s1l = s1.size();
	size_t s2l = s2.size();

	std::vector<p_constr> constraints;
	if (has_constraints)
		read_constraints(constraints);


	/* define the mapping from char -> int */
	std::map<char, int> map1;
	std::map<char, int> map2;

	define_mapping(sigma1, map1);
	define_mapping(sigma2, map2);

	/* integer representations of strings and sigmas */
	std::vector<unsigned> s1i(s1l);
	std::vector<unsigned> s2i(s1l);
	std::vector<unsigned> sigma1i(sigma1l);
	std::vector<unsigned> sigma2i(sigma2l);

	std::iota(sigma1i.begin(), sigma1i.end(), 0);				// encoding for sigma1
	std::iota(sigma2i.begin(), sigma2i.end(), 0);				// encoding for sigma2
	for (size_t i = 0; i < s1l; ++i) s1i[i] = map1[s1[i]];		// encoding for s1
	for (size_t i = 0; i < s2l; ++i) s2i[i] = map2[s2[i]];		// encoding for s2


	if (_DEBUG) {
		std::cout << "p1: " << p1 << ", p2: " << p2 << endl;
		std::cout << "s len: " << s1.size() << ", " << s2.size() << endl;
		std::cout << "strings: " << s1 << ", " << s2 << endl;
		std::cout << "sigma len: " << sigma1.size() << ", " << sigma2.size() << endl;
		std::cout << "sigmas: " << sigma1 << ", " << sigma2 << endl;
		std::cout << "int rep (s1): ";
		for (size_t i = 0; i < s1l; ++i) std::cout << s1i[i]; std::cout << endl;
		std::cout << "int rep (s2): ";
		for (size_t i = 0; i < s2l; ++i) std::cout << s2i[i]; std::cout << endl;
		std::cout << "int rep (sigma1): ";
		for (size_t i = 0; i < sigma1l; ++i) std::cout << sigma1i[i]; std::cout << endl;
		std::cout << "int rep (sigma2): ";
		for (size_t i = 0; i < sigma2l; ++i) std::cout << sigma2i[i]; std::cout << endl;
	}


	/* identity (classical) matching schema */
	matching_schema<bool> ms(sigma1l, sigma2l, p1, p2, true, default_constraints_mode);
	ms.set_general(sigma1, sigma2, false);

	// TODO: CONSTRAINTS
	//if (has_constraints)
	//	ms.set_constraints(map1, map2, constraints, !default_constraints_mode);

	if (_DEBUG) {
		if (has_constraints) {
			std::cout << "Constraints: ";
			for(std::vector<p_constr>::size_type i = 0; i < constraints.size(); ++i)
				std::cout << constraints[i].first << ", " << constraints[i].second << endl;
		}
		ms.print_matching_schema(sigma1, sigma2);
	}

	/* create an edit distance object */
	edit_distance e;

	/* here we call the solver needed */
	int distance = -1;

	// Common execution of HC or EX
	if (!specific_perm && !specific_matrix) {
		if (heuristic == _BRUTEFORCE_ARG) {
			distance = bruteforce(s1i, s2i, s1l, s2l, sigma1i, sigma2i, sigma1l, sigma2l, ms, e);
		}
		else {
			distance = hill_climbing(s1i, s2i, s1l, s2l, sigma1i, sigma2i, sigma1l, sigma2l, p1, ms, e);
		}
	}
	// For a specific matching schema
	else {
		// input matching schema
		unsigned* sigma1_spe = new unsigned[sigma1l];
		unsigned* sigma2_spe = new unsigned[sigma2l];

		// if we have the permutations
		if (specific_perm) {
			// read string
			std::string read;

			// skip the obtained edit distance
			getline(std::cin, read);

			// read the specified matching schema
			read_specific_matchingschema(read, sigma1_spe);
			read_specific_matchingschema(read, sigma2_spe);

			// check if there is a matching schema or not
			bool there_is_a_ms = false;
			for (size_t i = 0; i < sigma1l; ++i)
				if (sigma1_spe[i]) {
					there_is_a_ms = true;
					break;
				}

			if (there_is_a_ms) {
				distance = e.edit_distance_matching_schema_enhanced(s1i, s2i, s1l, s2l, sigma1_spe, sigma2_spe, sigma1l, sigma2l, ms);
			}
			else
				distance = e.edit_distance_matching_schema(s1i, s2i, s1l, s2l, ms);

		// if we have the matrix
		} else if (specific_matrix) {
			// populate linearly
			for (size_t i = 0; i < sigma1l; ++i) sigma1_spe[i] = i;
			for (size_t i = 0; i < sigma2l; ++i) sigma2_spe[i] = i;

			// read the matrix
			read_specific_matrix(ms);

			// and now for something completely different...
			distance = e.edit_distance_matching_schema_enhanced(s1i, s2i, s1l, s2l, sigma1_spe, sigma2_spe, sigma1l, sigma2l, ms);
		}
	}

	std::cout << distance << endl;

	return 0;
}


// ========================================================================

int hill_climbing(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,  const size_t& s1l, const size_t& s2l,
		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2, const size_t& sig1l, const size_t& sig2l, const size_t& p1, matching_schema<bool>& m, edit_distance& e) {

	//std::cout << "enter the void (1)" << endl;

	unsigned d = e.edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	unsigned minDist = d;
	unsigned minMinDist = minDist;

	//std::cout << "enter the void (2)" << endl;

	// for the permutations
	unsigned* sigma1_o = new unsigned[sig1l]; std::iota(sigma1_o, sigma1_o + sig1l, 0);
	unsigned* sigma2_o = new unsigned[sig2l]; std::iota(sigma2_o, sigma2_o + sig2l, 0);
	unsigned* sigma1_t = new unsigned[sig1l]; std::iota(sigma1_t, sigma1_t + sig1l, 0);
	unsigned* sigma2_t = new unsigned[sig2l]; std::iota(sigma2_t, sigma2_t + sig2l, 0);

	//std::cout << "enter the void (3)" << endl;

	// for fixpoints
	unsigned* sigma1_min = new unsigned[sig1l]; std::iota(sigma1_min, sigma1_min + sig1l, 0);
	unsigned* sigma2_min = new unsigned[sig2l]; std::iota(sigma2_min, sigma2_min + sig2l, 0);
	unsigned* sigma1_min_min = new unsigned[sig1l]; std::iota(sigma1_min_min, sigma1_min_min + sig1l, 0);
	unsigned* sigma2_min_min = new unsigned[sig2l]; std::iota(sigma2_min_min, sigma2_min_min + sig2l, 0);

	//std::cout << "enter the void (4)" << endl;

	size_t attempts = 1, shuffle_tries = 2;
	unsigned tries = 0, k_shuffle = 0;

	//std::cout << "start: " << d << endl;

	bool improved = true;
	while (improved) {
		improved = false;

		for (size_t ip = 0; ip < sig1l; ip++) {
			for (size_t jp = ip; jp < sig1l; jp++) {

				//std::copy(sigma1_t.begin(), sigma1_t.end(), sigma1_o.begin());		// reset state
				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o);						// reset state
				std::swap(sigma1_o[ip], sigma1_o[jp]);									// swap

				if (isValid(sigma1_o, sig1l, p1)) {

					for (size_t ipp = 0; ipp < sig2l; ipp++) {
						for (size_t jpp = ipp; jpp < sig2l; jpp++) {

							//std::copy(sigma2_t.begin(), sigma2_t.end(), sigma2_o.begin());			// reset state
							std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);							// reset state
							std::swap(sigma2_o[ipp], sigma2_o[jpp]);									// swap


							//int newDistance=e.enhancedWithMatchingSchema(s1,s2,s1l,s2l,sigma1_o,sigma2_o,sigma1_t,sigma2_t,sig1l,sig2l,m,minDist);
							int newDistance=e.edit_distance_matching_schema_enhanced_with_diagonal(s1, s2, s1l, s2l, sigma1_o, sigma2_o, sig1l, sig2l, m, minDist);

							if (newDistance!=-1){
								//std::cout << "start upgrade" << endl;
								//std::cout << minDist << ' ' << newDistance << endl;
								minDist=newDistance;

								improved = true;
								//std::copy(sigma1_o.begin(), sigma1_o.end(), sigma1_min.begin());
								std::copy(sigma1_o, sigma1_o + sig1l, sigma1_min);
								//std::copy(sigma2_o.begin(), sigma2_o.end(), sigma2_min.begin());
								std::copy(sigma2_o, sigma2_o + sig2l, sigma2_min);

								//std::cout << "end upgrade" << endl;
							}

							/*int d = e.edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, sigma1_o, sigma2_o, sig1l, sig2l, m);
							std::cout << "d: " << d << endl;
							if (d < minDist) {
								minDist = d;

								improved = true;
								std::copy(sigma1_o, sigma1_o + sig1l, sigma1_min);
								std::copy(sigma2_o, sigma2_o + sig2l, sigma2_min);
							}*/

							//for (size_t x = 0; x < sig1l; ++x) std::cout << sigma1_o[x] << ' '; std::cout << endl;
							//for (size_t x = 0; x < sig2l; ++x) std::cout << sigma2_o[x] << ' '; std::cout << endl << endl;

						}
						//std::copy(sigma2_t.begin(), sigma2_t.end(), sigma2_o.begin());
						std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);
					}

				}
			}
		}

		if (improved) {
			// copy sigmaMin to sigmaOrig
			//std::copy(sigma1_min.begin(), sigma1_min.end(), sigma1_o);
			std::copy(sigma1_min, sigma1_min + sig1l, sigma1_o);
			//std::copy(sigma2_min.begin(), sigma2_min.end(), sigma2_o);
			std::copy(sigma2_min, sigma2_min + sig2l, sigma2_o);
			// copy sigmaOrig to sigmaTmp
			//std::copy(sigma1_o.begin(), sigma1_o.end(), sigma1_t.begin());
			std::copy(sigma1_o, sigma1_o + sig1l, sigma1_t);
			//std::copy(sigma2_o.begin(), sigma2_o.end(), sigma2_t.begin());
			std::copy(sigma2_o, sigma2_o + sig2l, sigma2_t);
		} else {
			if (minDist < minMinDist) {
				minMinDist = minDist;

				// copy sigmaMin to sigmaMinMin
				//std::copy(sigma1_min.begin(), sigma1_min.end(), sigma1_min_min.begin());
				std::copy(sigma1_min, sigma1_min + sig1l, sigma1_min_min);
				//std::copy(sigma2_min.begin(), sigma2_min.end(), sigma2_min_min.begin());
				std::copy(sigma2_min, sigma2_min + sig2l, sigma2_min_min);

				improved = true;
				tries = 0;
			}

			if (tries < attempts) {
				improved = true;
				tries++;

				// random swap
				// for sigma1, we try _SHUFFLE_TRIES times, then if is still not valid, we retry with the original one
				for (k_shuffle = 0; k_shuffle < shuffle_tries && !isValid(sigma1_t, sig1l, p1); ++k_shuffle)
					shuffle(sigma1_t, sig1l);
				if (k_shuffle == shuffle_tries)
					std::copy(sigma1_o, sigma1_o + sig1l, sigma1_t);
					//std::copy(sigma1_o.begin(), sigma1_o.end(), sigma1_t.begin());
				// no constraints on the shuffle for sigma2
				shuffle(sigma2_t, sig2l);

				//std::copy(sigma1_t.begin(), sigma1_t.end(), sigma1_o.begin());
				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o);
				//std::copy(sigma2_t.begin(), sigma2_t.end(), sigma2_o.begin());
				std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);

				minDist = e.edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, sigma1_o, sigma2_o, sig1l, sig2l, m);
			}
		}
	}

	delete[] sigma1_o; delete[] sigma1_t;
	delete[] sigma2_o; delete[] sigma2_t;
	delete[] sigma1_min; delete[] sigma1_min_min;
	delete[] sigma2_min; delete[] sigma2_min_min;

	return minMinDist;
}

int bruteforce(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2,  const size_t& s1l, const size_t& s2l,
		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2, const size_t& sig1l, const size_t& sig2l, const matching_schema<bool>& m, edit_distance& e) {

	unsigned distance = e.edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	unsigned current = distance;

	unsigned* perm1 = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) perm1[i] = i;
	unsigned* perm2 = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) perm2[i] = i;

	FixedED<unsigned> fixed_ed(s1l+1, s2l+1);

	do {
		do {
			current = fixed_ed.edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, perm1, perm2, sig1l, sig2l, m);

			if (current < distance)
				distance = current;

		} while(std::next_permutation(perm2, perm2+sig2l));
	} while(std::next_permutation(perm1, perm1+sig1l));

	return distance;
}


// ========================================================================


/**
 * extract_sigmas works with a bitset of length _ASCII_LEN = 255 - _ASCII_START
 */
void extract_sigma(std::string& s, std::string& e) {
	std::bitset<_ASCII_LEN> symbols;

	for (size_t i = 0; i < s.size(); ++i)
		symbols[(int) s[i]] = (symbols[(int) s[i]] || 1);

	for (size_t i = 0; i < _ASCII_LEN; ++i)
		if (symbols[i])
			e += (char) i;
}

void extract_sigma_fibre(std::vector<std::string>& s, std::vector<std::string>& e) {
	for (int i=0; i<s.size(); i++)
		e.push_back(s[i]);

	sort(e.begin(), e.end());
	e.erase(unique(e.begin(), e.end()), e.end());
}

void define_mapping(std::string& s, std::map<char,int>& m) {
	for (size_t i = 0; i < s.size(); ++i)
		m.insert(std::pair<char,int>(s[i], i));
}

void define_mapping_fibre(std::vector<std::string>& s, std::map<std::string,int>& m) {
	for (size_t i = 0; i < s.size(); ++i)
		m.insert(std::pair<std::string,int>(s[i], i));
}

void print_alignment(const Alignment<int>& alignment, const std::string& sigma1, const std::string& sigma2, const matching_schema<bool>& m, int& distance, bool is_identity) {
	// print first string
	for (size_t i = 0; alignment.a[i] != _END_ALIGNMENT; i++)
		std::cout << (alignment.a[i] != _GAP_FLAG ? sigma1[alignment.a[i]] : '-');
	std::cout << endl;

	// print second string
	for (size_t i = 0; alignment.b[i] != _END_ALIGNMENT; i++)
		std::cout << (alignment.b[i] != _GAP_FLAG ? sigma2[alignment.b[i]] : '-');
	std::cout << endl;

	// print the alignment (a string containing * or space)
	size_t index = 0;
	size_t temp = 0;
	for (; alignment.a[index] != _END_ALIGNMENT; index++) {
		// if both chars are not gaps
		if (alignment.a[index] != _GAP_FLAG && alignment.b[index] != _GAP_FLAG) {
			if (	!m.ms[alignment.a[index]][alignment.b[index]] || 									// if they are in match or
					(is_identity && (sigma1[alignment.a[index]] == sigma2[alignment.b[index]]))		) {	// the self_identity is active and they are the same char
				std::cout << '*';
				temp++;
			} else
				std::cout << ' ';
		} else
			std::cout << ' ';
	}
	std::cout << endl;

	// here we return the new edit distance including the self_identity case
	distance = index - temp;
}
