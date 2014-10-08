#ifndef __FIXATION_HPP__
#define __FIXATION_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <utility>
#include <gsl/gsl_rng.h>

using std::vector;
using std::pair;
using Sequence::chromosome;
using Sequence::arg;

void fixation(gsl_rng * r, int * NSAM, vector<chromosome> * sample, Sequence::arg * sample_history,
	      vector<pair<int,int> > * paired, int * nlinks, int *NLINKEDPAIRS,
	      int * nlinksbw, double * t,
	      int config[], const vector<double> & trajectory,
	      const double & pmin, const double & dt,
	      const int & nsam, const int & nsites, const int & nsites_bw,
	      const double & littler, const double & littleec, const double & p);

#endif
