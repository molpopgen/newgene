#ifndef __CA_HPP__
#define __CA_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <vector>
#include <utility>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace Sequence;

//void ca(gsl_rng * r, int * NSAM, vector<chromosome> * sample, arg * sample_history,
void ca(int * NSAM, vector<chromosome> * sample, Sequence::arg * sample_history,
	vector<pair<int,int> > * paired, int * nlinks,int config[],
	const int & ch1, const int & ch2,
	const double & t, const int & nsam, const int & nsites);

#endif
