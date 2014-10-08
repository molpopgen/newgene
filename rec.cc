#include <rec.hpp>
#include <util.hpp>
#include <cassert>
#include <gsl/gsl_randist.h>

using namespace std;

void rec( gsl_rng * r,
	  vector<pair<int,int> >* paired,
	  const int & current_nsam,
	  const int & NLINKEDPAIRS)
{
  //in this model, recombination just updates paired,
  //and does not change sample size
  int chindex = int(gsl_ran_flat(r,0.,double(NLINKEDPAIRS)));
  vector< pair<int,int> >::iterator pbegin = paired->begin();
  int dummy=-1,ch = -1;
#ifndef NDEBUG
  int np=0;
#endif
  for( vector<pair<int,int> >::iterator itr = pbegin ;
       itr != (pbegin+current_nsam) ; ++itr )
    {
      if(itr->first != -1 && itr->second != -1)
	{
	  ch = int(itr-pbegin);
	  ++dummy;
#ifndef NDEBUG
	  ++np;
#endif
	  if(dummy==chindex) break;
	}
    }
  assert(np>0);
  assert(ch>=0);
  assert( (pbegin+ch)->first != -1 );
  assert( (pbegin+ch)->second != -1 );
  assert( (pbegin+(pbegin+ch)->second)->first != -1 );
  assert( (pbegin+(pbegin+ch)->second)->second != -1 );

  (pbegin+(pbegin+ch)->second)->second = -1;
  (pbegin+ch)->second = -1;
  assert( ! unbalanced_pairs(paired,paired->size()) );
  assert( ! redundant_pairs(paired,current_nsam) );
}
