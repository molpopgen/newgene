#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <utility>
#include <functional>
#include <vector>
#include <limits>
#include <Sequence/Coalescent/SimTypes.hpp>
#include <gsl/gsl_rng.h>

using std::vector;
using std::pair;
using std::binary_function;
using Sequence::chromosome;


static const double DMAX = std::numeric_limits<double>::max();

struct fpair : public binary_function< pair<int,int>, int, bool >
{
  inline bool operator()(const pair<int,int> & pii, const int & i) const
  {
    return pii.first == i;
  }
};

struct fpairbothitr : public binary_function< vector<pair<int,int> >::const_iterator,
					      vector<pair<int,int> >::const_iterator, bool >
{
  inline bool operator()(vector<pair<int,int> >::const_iterator & left,
			 vector<pair<int,int> >::const_iterator & right) const
  {
    return( (left->first==right->first && left->second==right->second) ||
	    (left->first==right->second && left->second==right->first) );
  }
};

void update_paired( const int & i, vector< pair<int,int> > * paired );
//void popswitch( chromosome * c, int config[] );
void popswitch( vector<chromosome>::iterator  c, int config[] );
void popswitch( vector<chromosome>::iterator  c, int config[], const int & todeme );
int count_linked_pairs( const vector<pair<int,int> > & paired );
pair<int,int> pick2( gsl_rng * r, const int & current_nsam, 
		     const vector<pair<int,int> > & paired);
pair<int,int> pick2_not_in_deme( gsl_rng * r, vector<chromosome>::const_iterator sbegin,
				 vector<pair<int,int> >::const_iterator pbegin,
				 const int & current_nsam,const int & deme_nsam,
				 const int & not_deme );
pair<int,int> pick2_different_demes( gsl_rng * r, const int & NSAM,
				     const int config[],
				     const vector<chromosome> & sample,
				     const vector<pair<int,int> > & paired);

#ifndef NDEBUG
bool invalid_pairs(const vector<pair<int,int> > * paired,
		   const vector<chromosome> * sample);
bool demes_ok( const vector<chromosome> * sample, const int & NSAM, const int config[2] );
bool unbalanced_pairs( const vector<pair<int,int> > * paired,
		       const int & NSAM );
bool redundant_pairs( const vector<pair<int,int> > * paired ,
		      const int & NSAM );
#endif
#endif
