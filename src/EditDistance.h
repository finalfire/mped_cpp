#ifndef EDITDISTANCE_H_
#define EDITDISTANCE_H_

#include "Utility.h"
#include "Alignment.h"
#include "MatchingSchema.h"
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
//#include <tr1/unordered_map>

//Ho utilizzato pair_hash perchè unordered map non riusciva a definire l'hash per una key con pair<int,int>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        return h1 ^ h2;
    }
};


const short _END_ALIGNMENT = -2;
const short _GAP_FLAG = -1;

struct edit_distance {

// matrice della prima edit distance, consultata successivamente
unsigned** matrix;
//sono le coppie di indici (non della sequenza, ma posizioni nella matrice) da controllare poichè sono state soggette a variazioni con la permutazione
std::set<std::pair<int, int> >* check;
std::set<std::pair<int, int> >::iterator it;
std::set<std::pair<int, int> >::iterator it2;
//(key,value)->(ogni coppia di caratteri delle due sequenze, un set di posizioni di tale coppia(key) nella matrice)
std::unordered_map<std::pair<int,int>, std::set<std::pair<int, int> > ,pair_hash>* position;


void verify_permutation(const std::vector<unsigned short>& perm1, const std::vector<unsigned short>& perm2,
		const std::vector<unsigned short>& permOrigin1, const std::vector<unsigned short>& permOrigin2, const size_t& sig1l, const size_t& sig2l, int p){

	(*check).clear();

	for (int i=0; i<sig1l; i++){
		if (permOrigin1[i]!=perm1[i]){
			std::set<std::pair<int, int> > vO=(*position)[std::pair<int, int>(permOrigin1[i],perm2[i])];
			std::set<std::pair<int, int> > vOO=(*position)[std::pair<int, int>(perm1[i],perm2[i])];
			for (it=(vO).begin(); it!=(vO).end(); ++it){

					(*check).insert((*it));

			}
			for (it=(vOO).begin(); it!=(vOO).end(); ++it){

					(*check).insert((*it));

			}
		}
	}

	//faccio questo controllo perchè è inutile calcolare su p 1 1 in quanto ho solo la permutazione sul primo alfabeto
	if (p!=1){
		for (int i=0; i<sig2l; i++){
			if (permOrigin2[i]!=perm2[i]){
				std::set<std::pair<int, int> > vO=(*position)[std::pair<int, int>(permOrigin2[i],perm1[i])];
				std::set<std::pair<int, int> > vOO=(*position)[std::pair<int, int>(perm2[i],perm1[i])];
				for (it=(vO).begin(); it!=(vO).end(); ++it){

						(*check).insert((*it));

				}
				for (it=(vOO).begin(); it!=(vOO).end(); ++it){

						(*check).insert((*it));

				}
			}
		}
	}

	//std::cout<<check->size()<<"\n";
}

//Questa soluzione non calcola le celle che non sono cambiate rispetto alla matrice di origine
// e salta zone della matrice calcolando invece le celle nella quale i caratteri sono influenzati dalle permutazioni.
// Inoltre questa soluzione tiene conto allo stesso tempo di calcolare solo quelle celle che si trovano in un determinato range/intorno di un valore threshold;
// cosi facendo le celle che sono influenzate dalle permutazioni ma che si trovano fuori tale intorno vengono saltate perchè non contribuiscono a trovare
// una distance minore del threshold.
int incrementalSolution(const std::vector<unsigned short>& a, const std::vector<unsigned short>& b, const int& al, const int& bl, const std::vector<unsigned short>& sg1, const std::vector<unsigned short>& sg2,
		const std::vector<unsigned short>& permOrigin1, const std::vector<unsigned short>& permOrigin2, const size_t& sg1l, const size_t& sg2l, const matching_schema<bool>& m, int p, int threshold){

	//tale funzione mi va a controllare cosa è cambiato dalle precedenti permutazioni
	//	e mi dice a quali coppie della matrice guardare
	verify_permutation(sg1, sg2, permOrigin1, permOrigin2, sg1l, sg2l, p);

	int MAX = std::numeric_limits<int>::max();
	const int PASSED=9000;	//significa che l'ho analizzato nell'algoritmo qui e non è cambiato il valore
	const int NOTHING=9100; //significa che non l'ho proprio analizzato (potrebbe capitare sia che non deve essere analizzato o che ancora deve essere analizzato)

	if (abs(al-bl)>=threshold || threshold<0)
		return -1;

	unsigned** d=new unsigned*[al+1];
	for (int i=0; i<al+1; i++)
		d[i]=new unsigned[bl+1];

	for (int i=0; i<al+1; i++)
		for (int j=0; j<bl+1; j++)
			d[i][j]=NOTHING;

	int boundary=std::min(bl, threshold)+1;
	for (int i=0; i<boundary; i++)
		d[0][i]=i;

	//fill matrix
	for (int i=boundary; i<bl+1; i++)
		d[0][i]=MAX;

	for (int i=0; i<al+1; i++)
		d[i][0]=i;


	bool firstBoolean=true; //tale boolean rimane a true fino a quando non trova un primo cambiamento tra una cella della matrice attuale e di quella originale!
	for (it=(*check).begin(); it!=(*check).end(); ++it){

			int first=(*it).first;
			int second=(*it).second;

			int min=std::max(1, first-threshold);
			int max=(first>MAX-threshold)?bl:std::min(bl, first+threshold);

			if (min>max)
				return -1;

			if (min>1)
				d[first][min-1]=MAX;

			if (second<min || second>max || d[first][second]!=NOTHING)
				continue;


			if (firstBoolean){
				firstBoolean=false;
				if (!m.ms[index_of(a[first-1], sg1, sg1l)][index_of(b[second-1], sg2, sg2l)])
					d[first][second]=matrix[first-1][second-1];
				else
					d[first][second]=1+std::min(std::min(matrix[first][second-1],matrix[first-1][second]),matrix[first-1][second-1]);

				if (d[first][second]!=matrix[first][second]){

					if (first<al){
							//down
								(*check).insert(std::pair<int, int>(first+1, second));

						if (second<max){
							//diagonal
								(*check).insert(std::pair<int, int>(first+1, second+1));
							//right
								(*check).insert(std::pair<int, int>(first, second+1));
						}
					}
				}
				else{
					firstBoolean=true;
					d[first][second]=PASSED;
					}
			}else{
				if (!m.ms[index_of(a[first-1], sg1, sg1l)][index_of(b[second-1], sg2, sg2l)])
						d[first][second]=(d[first-1][second-1]!=PASSED && d[first-1][second-1]!=NOTHING)?d[first-1][second-1]:matrix[first-1][second-1];
				else
					d[first][second]=1+std::min(std::min((d[first][second-1]!=PASSED && d[first][second-1]!=NOTHING)?d[first][second-1]:matrix[first][second-1],(d[first-1][second]!=PASSED && d[first-1][second]!=NOTHING)?d[first-1][second]:matrix[first-1][second]),(d[first-1][second-1]!=PASSED && d[first-1][second-1]!=NOTHING)?d[first-1][second-1]:matrix[first-1][second-1]);

				if (d[first][second]!=matrix[first][second]){
					if (first<al){
						//down
							(*check).insert(std::pair<int, int>(first+1, second));

						if (second<max){
							//diagonal
								(*check).insert(std::pair<int, int>(first+1, second+1));
							//right
								(*check).insert(std::pair<int, int>(first, second+1));
						}
					}
				}
				else
					d[first][second]=PASSED;
			}
	}

	//verifica se la distance generata è minore del threshold, altrimenti assume che la distance è uguale a quella di origine
	unsigned my_dist;
	if (d[al][bl]!=PASSED && d[al][bl]!=NOTHING)
		my_dist= d[al][bl];
	else
	{
		for (size_t i = 0; i < al+1; ++i)
			delete[] d[i];
		delete[] d;
		return -1;
	}


	//std::cout<<my_dist<<" | Threshold"<<threshold<<" | giusto: "<<edit_distance_matching_schema_enhanced(a, b, al, bl,sg1, sg2, sg1l, sg2l,  m)<<"\n";
	//std::cout<<"\n";

	for (size_t i = 0; i < al+1; ++i)
		delete[] d[i];
	delete[] d;

	if (my_dist<threshold)
		return my_dist;

	return -1;
}


// tale funzione permette di calcolare la distanza se risulta essere inferiore al threshold
// riducendo il tempo di esecuzione a O(k*m), calcolando solo una striscia di diagonali 2k+1.
// 

int edit_distance_matching_schema_enhanced_with_diagonal(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2, int s1l, int s2l,
		unsigned* sig1, unsigned* sig2, const size_t& sig1l, const size_t& sig2l, const matching_schema<bool>& m, int threshold){

	int MAX = std::numeric_limits<int>::max();

	if (abs(s1l-s2l)>=threshold || threshold<0)
		return -1;

	// TODO: this must be fixed
	// uso la dimensione piu piccola tra le due s1l e s2l (m,n)
	/*if (s1l>s2l){
		unsigned *temp=s1;		// temp è un array come s1
		unsigned templ=s1l;		// lunghezza di s1
		s1=s2;					// s1 contiene s2 adesso
		s2=temp;				// s2 contiene s1 adesso
		s1l=s2l;				// lunghezza di s1 è quella di s2
		s2l=templ;				// lunghezza di s2 è quella di s1 (temp
	}*/

	// Permutation index
	unsigned* sig1_index = new unsigned[sig1l]; for (size_t i = 0; i < sig1l; ++i) { sig1_index[sig1[i]] = i; }
	unsigned* sig2_index = new unsigned[sig2l]; for (size_t i = 0; i < sig2l; ++i) { sig2_index[sig2[i]] = i; }

	// uso due righe e non la matrice completa
	unsigned* p=new unsigned[s1l+1];
	unsigned* d=new unsigned[s1l+1];

	int boundary=std::min(s1l, threshold)+1;
	for (int i=0; i<boundary; i++)
		p[i]=i;

	// fill array
	for (int i=boundary; i<s1l+1; i++)
		p[i]=MAX;
	for (int i=0; i<s1l+1; i++)
		d[i]=MAX;

	// iteration
	for (int j=1; j<=s2l; j++){
		unsigned charS2=s2[j-1];
		d[0]=j;

		int min=std::max(1, j-threshold);
		int max=(j>MAX-threshold)?s1l:std::min(s1l, j+threshold);

		if (min>max){
			return -1;
		}

		if (min>1){
			d[min-1]=MAX;
		}

		for (int i=min; i<=max; i++){
			if (!m.ms[sig1_index[s1[i-1]]][sig2_index[charS2]])
				d[i]=p[i-1];
			else
				d[i]=1+std::min(std::min(d[i-1],p[i]),p[i-1]);
		}

		//p=d;
		for (int c=0; c<s1l+1; c++)
			p[c]=d[c];
	}

	unsigned valuefinal=p[s1l];

	delete[] p;
	delete[] d;

	if (valuefinal<threshold)
		return valuefinal;

	return -1;
}

/*int edit_distance_matching_schema_enhanced_with_diagonal(const std::vector<unsigned>& s1, const std::vector<unsigned>& s2, int s1l, int s2l,
		std::vector<unsigned>& sig1, std::vector<unsigned>& sig2, const size_t& sig1l, const size_t& sig2l, const matching_schema<bool>& m, int threshold){
	
	int MAX = std::numeric_limits<int>::max();
	
	if (abs(s1l-s2l)>=threshold || threshold<0)
		return -1;

	// TODO: this must be fixed
	// uso la dimensione piu piccola tra le due s1l e s2l (m,n)
	//if (s1l>s2l){
	//	unsigned *temp=s1;		// temp è un array come s1
	//	unsigned templ=s1l;		// lunghezza di s1
	//	s1=s2;					// s1 contiene s2 adesso
	//	s2=temp;				// s2 contiene s1 adesso
	//	s1l=s2l;				// lunghezza di s1 è quella di s2
	//	s2l=templ;				// lunghezza di s2 è quella di s1 (temp
	//}
	
	// uso due righe e non la matrice completa
	unsigned* p=new unsigned[s1l+1]; 
	unsigned* d=new unsigned[s1l+1];
	
	int boundary=std::min(s1l, threshold)+1;
	for (int i=0; i<boundary; i++)
		p[i]=i;
	
	// fill array
	for (int i=boundary; i<s1l+1; i++)
		p[i]=MAX;
	for (int i=0; i<s1l+1; i++)
		d[i]=MAX;

	// iteration
	for (int j=1; j<=s2l; j++){
		unsigned charS2=s2[j-1];
		d[0]=j;
		
		int min=std::max(1, j-threshold);
		int max=(j>MAX-threshold)?s1l:std::min(s1l, j+threshold);
		
		if (min>max){
			return -1;
		}
		
		if (min>1){
			d[min-1]=MAX;
		}
		
		for (int i=min; i<=max; i++){
			if (!m.ms[index_of(s1[i-1], sig1, sig1l)][index_of(charS2, sig2, sig2l)])
				d[i]=p[i-1];
			else
				d[i]=1+std::min(std::min(d[i-1],p[i]),p[i-1]);
		}
		
		//p=d;
		for (int c=0; c<s1l+1; c++)
			p[c]=d[c];
	}

	unsigned valuefinal=p[s1l];

	delete[] p;
	delete[] d;

	if (valuefinal<threshold)
		return valuefinal;
	
	return -1;
}*/


unsigned edit_distance_matching_schema_enhanced(const std::vector<unsigned>& a, const std::vector<unsigned>& b, const size_t& al, const size_t& bl,
		unsigned*& sg1, unsigned*& sg2, const size_t& sg1l, const size_t& sg2l, const matching_schema<bool>& m) {

	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];

	// Permutation index
	unsigned* sig1_index = new unsigned[sg1l]; for (size_t i = 0; i < sg1l; ++i) { sig1_index[sg1[i]] = i; }
	unsigned* sig2_index = new unsigned[sg2l]; for (size_t i = 0; i < sg2l; ++i) { sig2_index[sg2[i]] = i; }

	// (1) first row and first column
	for (size_t i = 0; i < al+1; ++i) d[i][0] = i;
	for (size_t j = 0; j < bl+1; ++j) d[0][j] = j;

	// (2) fill the matrix
	for (size_t i = 1; i < al+1; ++i){
		for (size_t j = 1; j < bl+1; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,																		// deletion
					d[i][j-1] + 1,																		// insertion
					d[i-1][j-1] + 1 * m.ms[sig1_index[a[i-1]]][sig2_index[b[j-1]]]	// if in the matching schema there's a false, they match
			);
		}
	}
	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return my_dist;
}

/*unsigned edit_distance_matching_schema_enhanced(const std::vector<unsigned>& a, const std::vector<unsigned>& b, const size_t& al, const size_t& bl,
		const std::vector<unsigned>& sg1, const std::vector<unsigned>& sg2, const size_t& sg1l, const size_t& sg2l, const matching_schema<bool>& m) {

	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i < al+1; ++i) d[i][0] = i;
	for (size_t j = 0; j < bl+1; ++j) d[0][j] = j;

	// (2) fill the matrix
	for (size_t i = 1; i < al+1; ++i){
		for (size_t j = 1; j < bl+1; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,																		// deletion
					d[i][j-1] + 1,																		// insertion
					d[i-1][j-1] + 1 * m.ms[index_of(a[i-1], sg1, sg1l)][index_of(b[j-1], sg2, sg2l)]	// if in the matching schema there's a false, they match
			);
		}
	}
	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return my_dist;
}*/

unsigned edit_distance_matching_schema(const std::vector<unsigned>& a, const std::vector<unsigned>& b, const size_t& al, const size_t& bl, const matching_schema<bool>& m) {

	check=new std::set<std::pair<int, int> >();
	position=new std::unordered_map<std::pair<int,int>, std::set<std::pair<int, int> > ,pair_hash>();


	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i < al+1; ++i) d[i][0] = i;
	for (size_t j = 0; j < bl+1; ++j) d[0][j] = j;

	// (2) fill the matrix
	for (size_t i = 1; i < al+1; ++i){
		for (size_t j = 1; j < bl+1; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,							// deletion
					d[i][j-1] + 1,							// insertion
					d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]	// if in the matching schema there's a false, they match
			);
			
			std::set<std::pair<int, int> >* sett=&((*position)[std::pair<int,int>(a[i-1], b[j-1])]);
			(sett)->insert(std::pair<int, int>(i, j));
		}
	}

	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	//TODO decommentare e commentare la distruzione della matrice sotto (incremental solution)
	matrix=d;

	for (size_t i = 0; i < al; ++i) delete[] d[i];
		delete[] d;

	return my_dist;
}

Alignment<int> compute_alignment(const unsigned* a, const unsigned* b, const size_t& al, const size_t& bl, const matching_schema<bool>& m) {
	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];
	char** path = new char*[al+1];
	for (size_t i = 0; i < al + 1; ++i) path[i] = new char[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i <= al; ++i) { d[i][0] = i; path[i][0] = 'n'; }
	for (size_t j = 0; j <= bl; ++j) { d[0][j] = j; path[0][j] = 'o'; }

	// (2) fill the matrix
	for (size_t i = 1; i <= al; ++i) {
		for (size_t j = 1; j <= bl; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,								// deletion
					d[i][j-1] + 1,								// insertion
					d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]		// if in the matching schema there's a false, they match
			);


			// (3) annotating the path for the backtrace
			if (d[i][j] == (d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]))
				path[i][j] = 'd';
			else if (d[i][j] == d[i-1][j] + 1)
				path[i][j] = 'n';
			else
				path[i][j] = 'o';
		}
	}

	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	// (4) here I start the backtrace
	unsigned maxl = (al > bl) ? al : bl;
	int x = 2 * maxl - 1;

	std::vector<int> all1(2*maxl, -1);
	std::vector<int> all2(2*maxl, -1);

	int i = al;
	int j = bl;

	while (i >= 0 || j >= 0) {
		if (i == 0 && j == 0) {
			break;
		} else {
			if (path[i][j] == 'n') {
				all1[x] = a[i-1];
				all2[x] = _GAP_FLAG;
				i = i-1;
				x--;
			} else if (path[i][j] == 'd') {
				all1[x] = a[i-1];
				all2[x] = b[j-1];
				i = i-1;
				j = j-1;
				x--;
			} else {
				all1[x] = _GAP_FLAG;
				all2[x] = b[j-1];
				j = j-1;
				x--;
			}
		}
	}

	// (5) here I define the alignment
	Alignment<int> alignment(my_dist, 2*maxl-x);

	size_t k = 0;
	for (k = x+1; k <= 2*maxl-1; k++)
		alignment.a[k-x-1] = all1[k];
	alignment.a[k-x-1] = -2;
	for (k = x+1; k <= 2*maxl-1; k++)
		alignment.b[k-x-1] = all2[k];
	alignment.b[k-x-1] = -2;


	// (6) deallocate everything
	for (size_t i = 0; i < al; ++i) delete[] path[i];
	delete[] path;
	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return alignment;
}

int distance_from_alignment(const Alignment<int>& alignment, const std::string& sigma1, const std::string& sigma2, const matching_schema<bool>& m, bool is_identity) {
	// print the alignment (a string containing * or space)
	size_t index = 0;
	size_t temp = 0;
	for (; alignment.a[index] != _END_ALIGNMENT; index++) {
		// if both chars are not gaps
		if (alignment.a[index] != _GAP_FLAG && alignment.b[index] != _GAP_FLAG) {
			if (	!m.ms[alignment.a[index]][alignment.b[index]] || 									// if they are in match or
					(is_identity && (sigma1[alignment.a[index]] == sigma2[alignment.b[index]]))		) {	// the self_identity is active and they are the same char
				temp++;
			}
		}
	}

	// here we return the new edit distance including the self_identity case
	return index - temp;
}

};

#endif /* EDITDISTANCE_H_ */
