#include <Sequence/Coalescent/Recombination.hpp>
#include <egc.hpp>
#include <coalesce.hpp>
#include <util.hpp>
#include <cassert>
#include <iostream>
#include <utility>
#include <cmath>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

int conditional_tract_length(gsl_rng * r, double p, int i, int L)
/*
  Returns a random deviate from a truncated geometric distribution.
  The tract length returned is such that it begins at position i,
  and includes, at most, up to position L
*/
{
  int lim1 = L-i+1;
  double b = 1.-pow(1.-p,lim1);
  double lnq = log(1.-p);
  return 1+int(floor(log(1-gsl_rng_uniform(r)*b)/lnq));
}

pair<int,int> pick_uniform_spot_gc_deme(const double & random_01,
					const int & nlinksgc_deme,
					std::vector<chromosome>::const_iterator sample_begin,
					const unsigned & current_nsam,
					const int & deme)
{
  
#ifndef NDEBUG
  int check=0,i=0;
  while( i < int(current_nsam))
    {
      if((sample_begin+i)->pop == deme )
	check+=(sample_begin+i)->links()+1;
      ++i;
    }
  if(check!=nlinksgc_deme) cerr << check << ' ' << nlinksgc_deme << '\n';
  assert(check==nlinksgc_deme);
#endif
  
  int pos = int(random_01*double(nlinksgc_deme))+1;
  int recombinant = 0,len=0;
  bool flag = false;
  while( (sample_begin+recombinant) <
	 (sample_begin+current_nsam) )
    {
      if( (sample_begin+recombinant)->pop == deme )
	{
	  len = (sample_begin+recombinant)->links()+1;
	  if(pos <= len)
	    {
	      flag=true;
	      break;
	    }
	  pos -= len;
	}
      if(!flag)
	recombinant++;
    }
  //FIX -- DOUBLE CHECK THIS
  int rpos = (sample_begin+recombinant)->begin()->beg + pos - 1;
  assert(rpos >= (sample_begin+recombinant)->first() );
  assert(rpos <= (sample_begin+recombinant)->last() );
  //rpos has the interpretation of being the first position on
  //chromosome recombinant which is involved in the conversion event
  return std::make_pair(recombinant,rpos);
}

pair<int,int> pick_uniform_spot_gc_not_deme(const double & random_01,
					const int & nlinksgc_deme,
					std::vector<chromosome>::const_iterator sample_begin,
					const unsigned & current_nsam,
					const int & notdeme)
{
  /*
#ifndef NDEBUG
  int check=0,i=0;
  while( i < int(current_nsam))
    {
      check+=(sample_begin+i)->links()+1;
      ++i;
    }
  if(check!=nlinksgc) cerr << check << ' ' << nlinksgc << '\n';
  assert(check==nlinksgc);
#endif
  */
  int pos = int(random_01*double(nlinksgc_deme))+1;
  int recombinant = 0,len=0;
  bool flag = false;
  while( (sample_begin+recombinant) <
	 (sample_begin+current_nsam) )
    {
      if( (sample_begin+recombinant)->pop != notdeme )
	{
	  len = (sample_begin+recombinant)->links()+1;
	  if(pos <= len)
	    {
	      flag=true;
	      break;
	    }
	  pos -= len;
	}
      if(!flag)
	recombinant++;
    }
  //FIX -- DOUBLE CHECK THIS
  int rpos = (sample_begin+recombinant)->begin()->beg + pos - 1;
  assert(rpos >= (sample_begin+recombinant)->first() );
  assert(rpos <= (sample_begin+recombinant)->last() );
  //rpos has the interpretation of being the first position on
  //chromosome recombinant which is involved in the conversion event
  return std::make_pair(recombinant,rpos);
}
pair<int,int> pick_uniform_spot_gc(const double & random_01,
				   const int & nlinksgc,
				   std::vector<chromosome>::const_iterator sample_begin,
				   const unsigned & current_nsam)
{
#ifndef NDEBUG
  int check=0,i=0;
  while( i < int(current_nsam))
    {
      check+=(sample_begin+i)->links()+1;
      ++i;
    }
  if(check!=nlinksgc) cerr << check << ' ' << nlinksgc << '\n';
  assert(check==nlinksgc);
#endif
  int pos = int(random_01*double(nlinksgc))+1;
  int recombinant = 0,len=0;
  while( (sample_begin+recombinant) <
	 (sample_begin+current_nsam) )
    {
      len = (sample_begin+recombinant)->links()+1;
      if(pos <= len)break;
      pos -= len;
      recombinant++;
    }
  //FIX -- DOUBLE CHECK THIS
  int rpos = (sample_begin+recombinant)->begin()->beg + pos - 1;
  assert(rpos >= (sample_begin+recombinant)->first() );
  assert(rpos <= (sample_begin+recombinant)->last() );
  //rpos has the interpretation of being the first position on
  //chromosome recombinant which is involved in the conversion event
  return std::make_pair(recombinant,rpos);
}

//FIX -- must be changed so that the function is told what pop it is going to
void egc(gsl_rng * r, const pair<int,int> & two, const int & to_deme, const double & t,  const int & nsam, const int & nsites, 
	 const int & nlinksgc, const double & p,
	 int * NSAM, int * nlinks, int config[], 
	 vector<chromosome> * sample, arg * sample_history, vector< pair<int,int> > * paired)
{
  vector<chromosome>::iterator sbegin = sample->begin();
  vector< pair<int,int> >::iterator pbegin = paired->begin();
  //cerr << *nlinks << ' ' << *NSAM << ' ' << nlinksgc  << '\n';
  
  int tl = conditional_tract_length(r,p,two.second,nsites);
  int deme = (sbegin+two.first)->pop;
  if( two.second == (sbegin+two.first)->first() )
    {
      if( two.second + tl >= (sbegin+two.first)->last() )
	{
	  //cerr << "EGC1 ";

	  if( (pbegin+two.first)->second != -1 )
	    {
	      ( pbegin + (pbegin+two.first)->second )->second = -1;
	      ( pbegin + two.first )->second = -1;
	      assert( ! unbalanced_pairs(paired,*NSAM) );
	    }
	  //popswitch( &*(sbegin+two.first), config );
	  popswitch( (sbegin+two.first), config,to_deme );
	  assert( demes_ok(sample,*NSAM,config) );
	}
      else
	{
	  //cerr << "EGC2 ";
	  //FIX -- DOUBLE CHECK THIS
	  (*nlinks) -= crossover(*NSAM,two.first,two.second+tl-1,sample,sample_history);
	  (*NSAM)++;
	  config[deme]++;
	  assert( demes_ok(sample,*NSAM,config) );

	  sbegin = sample->begin();
	  assert( (sbegin+two.first)->nsegs > 0 );
	  assert( (sbegin+*NSAM-1)->nsegs > 0 );


	  paired->push_back( make_pair(*NSAM-1,-1) );
	  pbegin = paired->begin();
	  if( (pbegin+two.first)->second != -1)
	    {
	      (pbegin+(pbegin+two.first)->second)->second = *NSAM-1;
	      swap( (pbegin+*NSAM-1)->second, (pbegin+two.first)->second );
	      assert( ! unbalanced_pairs(paired,*NSAM) );
	    }
	  //popswitch( &*(sbegin+two.first),config);
	  popswitch( (sbegin+two.first),config,to_deme);
	  assert( demes_ok(sample,*NSAM,config) );
	}
    }
  else
    {
      //cerr << "EGC3 ";
      //need to check for "do nothing events", where 
      //the tract length is contained entirely within
      //ancestral material
      Sequence::chromosome::const_iterator seg = (sbegin+two.first)->begin();//sample[two.first].begin();
      //for( ; two.second >= seg->end ; ++seg );
      for( ; two.second > seg->end ; ++seg );
      /*
      if(seg == (sbegin+two.first)->end())
	{
	  cerr << two.second << ' ' << *(sbegin+two.first) << '\n';
	}
      */
      assert(seg != (sbegin+two.first)->end());
      bool in1 = (two.second >=seg->beg)?true:false;
      seg=(sbegin+two.first)->begin();
      //for( ; two.second+tl >= seg->end && seg < (sbegin+two.first)->end() ; ++seg );
      for( ; two.second+tl > seg->end && seg < (sbegin+two.first)->end() ; ++seg );

      bool in2 = (two.second+tl >=seg->beg && seg < (sbegin+two.first)->end() )?true:false;
      /*
      cerr << "in: "<< in1 << ' ' << in2 << '\n';
      cerr << two.second << ' ' << two.second+tl << ' '
	   << *(sbegin+two.first) << '\n';
      */
      //assert(seg != (sbegin+two.first)->end());			     
      if(in1||in2)
	{
	  //FIX -- DOUBLE CHECK THIS
	  (*nlinks) -= crossover(*NSAM,two.first,two.second-1,sample,sample_history);
	  sbegin = sample->begin();
	  (*NSAM)++;
	  config[deme]++;
	  assert( demes_ok(sample,*NSAM,config) );
	  assert( (sbegin+two.first)->nsegs > 0 );
	  assert( (sbegin+*NSAM-1)->nsegs > 0 );
	  paired->push_back( make_pair(*NSAM-1,-1) );
	  pbegin = paired->begin();
	  if( two.second+tl >= (sbegin+*NSAM-1)->last())
	    {
	      //cerr << "a ";
	      //popswitch( &*(sbegin+*NSAM-1),config );
	      popswitch( (sbegin+*NSAM-1),config,to_deme );
	      assert( demes_ok(sample,*NSAM,config) );
	    }
	  else
	    {
	      //cerr << "b ";
	      //(*nlinks) -= crossover(*NSAM,*NSAM-1,two.second+tl-1,sample,sample_history);
	      //FIX -- DOUBLE CHECK THIS
	      (*nlinks) -= crossover(*NSAM,*NSAM-1,two.second+tl,sample,sample_history);
	      sbegin = sample->begin();
	      (*NSAM)++;
	      config[deme]++;
	      assert( demes_ok(sample,*NSAM,config) );
	      assert( (sbegin+*NSAM-2)->nsegs > 0 );
	      assert( (sbegin+*NSAM-1)->nsegs > 0 );
	      paired->push_back( make_pair(*NSAM-1,-1) );
	      //cerr << *NSAM << ' ' << (paired->end()-1)->first 
	      //<< ' '  << (paired->end()-1)->second << '\n';
	      pbegin = paired->begin();
	      //popswitch( &*(sbegin+(*NSAM)-2),config );
	      popswitch( (sbegin+(*NSAM)-2),config,to_deme );
	      assert( demes_ok(sample,*NSAM,config) );

	      int rv = ::coalesce(t,2*nsam,two.first,*NSAM-1,nsites,nlinks,sample,sample_history);
	      sbegin = sample->begin();
	      (*NSAM) -= rv;
	      config[deme] -= rv;
	      assert( demes_ok(sample,*NSAM,config) );
	      update_paired( *(NSAM)-1,paired );
	      if(rv==2)
		update_paired( two.first,paired );
	      assert(! invalid_pairs(paired,sample) );
	    }
	}
    }
  assert( ! unbalanced_pairs(paired,*NSAM) );
}
