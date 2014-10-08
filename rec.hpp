#ifndef __REC_HPP__
#define __REC_HPP__

#include <gsl/gsl_rng.h>
#include <vector>
#include <utility>

void rec( gsl_rng * r,
	  std::vector<std::pair<int,int> >* paired,
	  const int & current_nsam,
	  const int & NLINKEDPAIRS);

#endif
