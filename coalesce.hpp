#ifndef __COALESCE_HPP__
#define __COALESCE_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <vector>

using namespace std;
using namespace Sequence;

int coalesce(const double & time,
	     const int & ttl_nsam,
	     //const int & current_nsam,
	     const int & c1,
	     const int & c2,
	     const int & nsites,
	     int * nlinks,
	     std::vector<chromosome> * sample,
	     arg * sample_history);

#endif
