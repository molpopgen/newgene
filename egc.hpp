#ifndef __EGC_HPP__
#define __EGC_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <gsl/gsl_rng.h>
#include <vector>
#include <utility>

using std::vector;
using std::pair;
using Sequence::chromosome;
using Sequence::arg;

pair<int,int> pick_uniform_spot_gc(const double & random_01,
				   const int & nlinksgc,
				   std::vector<chromosome>::const_iterator sample_begin,
				   const unsigned & current_nsam);

pair<int,int> pick_uniform_spot_gc_deme(const double & random_01,
					const int & nlinksgc_deme,
					std::vector<chromosome>::const_iterator sample_begin,
					const unsigned & current_nsam,
					const int & deme);
pair<int,int> pick_uniform_spot_gc_not_deme(const double & random_01,
					    const int & nlinksgc_deme,
					    std::vector<chromosome>::const_iterator sample_begin,
					    const unsigned & current_nsam,
					    const int & not_deme);
void egc(gsl_rng * r, const pair<int,int> & two, const int & to_deme,
	 const double & t,  const int & nsam, const int & nsites,
	 const int &nlinksgc, const double & p, int * NSAM, int * nlinks,
	 int config[], vector<chromosome> * sample, Sequence::arg * sample_history, vector< pair<int,int> > * paired);

#endif 
