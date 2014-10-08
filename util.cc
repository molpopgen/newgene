#include <boost/bind.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <util.hpp>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <boost/bind.hpp>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

void update_paired( const int & i, vector< pair<int,int> > * paired )
{
  paired->erase( paired->begin() + i);
  for(vector< pair<int,int> >::iterator itr = paired->begin();
      itr != paired->end() ;++itr)
    {
      if( itr->second > i ) itr->second--;
      if( itr->first > i ) itr->first--;
    }
}

//void popswitch( chromosome * c, int config[] )
void popswitch( vector<chromosome>::iterator  c, int config[] )
{
  config[c->pop]--;
  config[!(c->pop)]++;
  c->pop = !(c->pop);
}

void popswitch( vector<chromosome>::iterator  c, int config[], const int & todeme )
{
  config[c->pop]--;
  //config[!(c->pop)]++;
  config[todeme]++;
  //c->pop = !(c->pop);
  c->pop = todeme;
}
pair<int,int> pick2( gsl_rng * r, const int & current_nsam, 
		     const vector<pair<int,int> > & paired)
{
  if(current_nsam==2) return std::make_pair(0,1);
  int i = int(gsl_ran_flat(r,0.,current_nsam));
  int j;
  do
    {
      j = int(gsl_ran_flat(r,0.,current_nsam));
    }
  while( j==i || ( paired[i].second==j && paired[j].second==i) );
  if(paired[i].second != -1)
    assert(paired[i].second != j);
  if(paired[j].second != -1)
    assert(paired[j].second != i);

  if(i>j) std::swap(i,j);
  return std::make_pair(i,j);
}

pair<int,int> pick2_different_demes( gsl_rng * r, const int & NSAM,
				     const int config[],
				     const vector<chromosome> & sample,
				     const vector<pair<int,int> > & paired)
{
  pair<int,int> two;
  do
    {
      two = pick2(r,NSAM,paired);
    }
  while( sample[two.first].pop == sample[two.second].pop );
  assert(two.first != two.second);
  assert(sample[two.first].pop != sample[two.second].pop );
  assert(two.first < two.second );
  assert(two.first<NSAM);
  assert(two.second<NSAM);
  return two;
}

//FIX -- does not pick 2 that are not linked to each other!!
pair<int,int> pick2_not_in_deme( gsl_rng * r, vector<chromosome>::const_iterator sbegin,
				 vector<pair<int,int> >::const_iterator pbegin,
				 const int & current_nsam,const int & deme_nsam,
				 const int & not_deme )
{
    assert( deme_nsam > 0 );
    assert( not_deme >= 0 );
    assert( deme_nsam <= current_nsam );
    int c1=0,c2=0;
    do
      {
	std::pair<int,int> two = Sequence::pick2(boost::bind(gsl_ran_flat,r,_1,_2),deme_nsam);
	assert(two.first != two.second);
	c1=c2=0;
	for(int i=0,j=0;i<current_nsam;++i)
	  {
	    if( (sbegin+i)->pop != not_deme && j == two.first )
	      c1 = i;
	    else if ( (sbegin+i)->pop != not_deme && j == two.second )
	      c2 = i;
	    if( (sbegin+i)->pop != not_deme ) ++j;
	  }
      }
    while( (pbegin+c1)->second == c2 && (pbegin+c2)->second == c1 );
    assert(c1!=c2);
    assert( (sbegin+c1)->pop != not_deme );
    assert( (sbegin+c2)->pop != not_deme );
    return std::make_pair(std::min(c1,c2),std::max(c1,c2));
}
int count_linked_pairs( const vector<pair<int,int> > & paired )
{
  vector< vector<pair<int,int> >::const_iterator > seen;
  for( vector<pair<int,int> >::const_iterator itr = paired.begin() ;
       itr != paired.end() ; ++itr )
    {
      if(itr->first != -1 && itr->second != -1)
	{
	  if ( find_if(seen.begin(),seen.end(),boost::bind(fpairbothitr(),_1,itr)) == seen.end() )
	    {
	      seen.push_back(itr);
	    }
	}
    }
  return ( int(seen.size()) );
}

#ifndef NDEBUG
bool invalid_pairs(const vector<pair<int,int> > * paired,
		   const vector<chromosome> * sample)
{
  for(unsigned i=0;i<paired->size();++i)
    {
      vector<pair<int,int> >::const_iterator itr = find_if(paired->begin(),paired->end(),boost::bind(fpair(),_1,i));
      assert( itr != paired->end() );
      // IF there is linkage b/w genes, make sure that the two chromose are in different populations
      if( itr->first != -1 && itr->second != -1 )
	if( (sample->begin()+itr->first)->pop == (sample->begin()+itr->second)->pop ) return true;
    }
  return false;
}

bool demes_ok( const vector<chromosome> * sample, const int & NSAM, const int config[2] )
{
  int c1=0,c2=0;
  vector<chromosome>::const_iterator sbegin = sample->begin();
  for(int i=0;i<NSAM;++i)
    {
      if( (sbegin+i)->pop==0 ) c1++;
      else if ( (sbegin+i)->pop == 1 ) c2++;
      //else abort();
    }
  if ( ! (c1==config[0] && c2==config[1] ) )
    cerr << c1 << ' ' << config[0] << ' ' << c2 << ' ' << config[1] << '\n';
  return( c1==config[0] && c2==config[1] );
}



bool unbalanced_pairs( const vector<pair<int,int> > * paired,
		       const int & NSAM )
{
  vector< pair<int,int> >::const_iterator pbeg = paired->begin();
  for( int i=0 ; i < NSAM ; ++i )
    {
      if ( (pbeg+i)->second != -1 )
	{
	  if ( (pbeg+(pbeg+i)->second)->second == -1 ||
	       (pbeg+(pbeg+i)->second)->second != (pbeg+i)->first )
	    {
	      cerr << "pair " << (i+1) << " of " << NSAM << " unbalanced " 
		   << (pbeg+i)->first << ' ' << (pbeg+i)->second << ' ' 
		   << (pbeg+(pbeg+i)->second)->first << ' ' << (pbeg+(pbeg+i)->second)->second << '\n';
	      return true;
	    }
	}
    }
  return false;
}

bool redundant_pairs( const vector<pair<int,int> > * paired ,
		      const int & NSAM )
{
  vector< pair<int,int> >::const_iterator pbeg = paired->begin();
  for( int i=0 ; i < NSAM ; ++i )
    {
      if ( (pbeg+i)->second != -1 )
	{
	  for(int j=i+1;j<NSAM;++j)
	    {
	      if ( (pbeg+j)->second == (pbeg+i)->second )
		return true;
	    }
	}
    }
  return false;;
}

#endif
