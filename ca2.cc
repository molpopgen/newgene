#include <Sequence/Coalescent/Coalesce.hpp>
#include <ca.hpp>
#include <coalesce.hpp>
#include <util.hpp>


void ca(int * NSAM, vector<chromosome> * sample, arg * sample_history, 
	vector<pair<int,int> > * paired, int * nlinks, int config[],
	const int & ch1, const int & ch2,
	const double & t, const int & nsam, const int & nsites)
{
  pair<int,int> two(ch1,ch2),othertwo;
  int rv,rv2,ch;
  vector<chromosome>::iterator sbegin = sample->begin();
  vector<pair<int,int> >::iterator pbegin = paired->begin();
  int firstflag,secondflag;
#ifndef NDEBUG
  if(*NSAM==2 && (pbegin+1)->second == 0)
    assert( pbegin->second != 1);
#endif

  int firstdeme = (sbegin+two.first)->pop, seconddeme = (sbegin+two.second)->pop;
  bool firstlinked = ( (pbegin+two.first)->second != -1 ) ? true : false;
  bool secondlinked = ( (pbegin+two.second)->second != -1 ) ? true : false;
 
  if ( firstdeme == seconddeme ) //both chromosomes are frome the same locus
    {
      if( firstlinked && secondlinked )
	{
	  othertwo.first = min( (pbegin+two.first)->second,(pbegin+two.second)->second );
	  othertwo.second = max( (pbegin+two.first)->second,(pbegin+two.second)->second );
	  seconddeme = (sbegin+othertwo.first)->pop;
	  
	  assert( (sbegin+two.first)->pop == (sbegin+two.second)->pop );
	  assert( (sbegin+othertwo.first)->pop == (sbegin+othertwo.second)->pop );
	  assert( (sbegin+two.first)->pop != (sbegin+othertwo.second)->pop );

	  //update linkage
	  (pbegin + two.first)->second = (pbegin + two.second)->second = othertwo.first;
	  (pbegin + othertwo.first)->second = (pbegin + othertwo.second)->second = two.first;

	  //do coalescent events
	  if( two.second > othertwo.second )
	    {
	      rv = ::coalesce(t,2*nsam,two.first,two.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[firstdeme] -= rv;
	      firstflag=secondflag=0;
	      if( rv ==2 )
		{
		  if(two.first < othertwo.first )
		    {
		      firstflag=secondflag=1;
		    }
		  else if (two.first < othertwo.second)
		    {
		      secondflag=1;
		    }
		}
	      assert ( (sample->begin()+othertwo.first-firstflag)->pop ==
		       (sample->begin()+othertwo.second-secondflag)->pop );
	      rv2 = ::coalesce(t,2*nsam,othertwo.first-firstflag,
			       othertwo.second-secondflag,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv2;
	      config[seconddeme] -= rv2;
	      if(rv==2)
		{
		  (pbegin+othertwo.first)->second = -1;
		  (pbegin+othertwo.second)->second = -1;
		}
	      if (rv2==2)
		{
		  (pbegin+two.first)->second = -1;
		  (pbegin+two.second)->second = -1;
		}
	      update_paired( two.second,paired );
	      if(rv==2)
		{
		  if(two.first < othertwo.second)
		    {
		      update_paired(othertwo.second,paired);
		      update_paired(two.first,paired);
		    }
		  else
		    {
		      update_paired(two.first,paired);
		      update_paired(othertwo.second,paired);
		    }
		}
	      else
		{
		  update_paired( othertwo.second,paired );
		}
	      if(rv2==2)
		{
		  update_paired(othertwo.first,paired);
		}
	    }
	  else
	    {
	      rv = ::coalesce(t,2*nsam,othertwo.first,othertwo.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[seconddeme] -= rv;
	      firstflag=secondflag=0;
	      if(rv==2)
		{
		  if(othertwo.first < two.first)
		    firstflag=secondflag=1;
		  else if (othertwo.first < two.second)
		    secondflag=1;
		}
	      assert ( (sample->begin()+two.first-firstflag)->pop ==
		       (sample->begin()+two.second-secondflag)->pop );
	      rv2 = ::coalesce(t,2*nsam,two.first-firstflag,
			       two.second-secondflag,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv2;
	      config[firstdeme] -= rv2;

	      if(rv==2)
		{
		  (pbegin+two.first)->second = -1;
		  (pbegin+two.second)->second = -1;
		}
	      if (rv2==2)
		{
		  (pbegin+othertwo.first)->second = -1;
		  (pbegin+othertwo.second)->second = -1;
		}
	      update_paired( othertwo.second,paired );
	      if(rv==2)
		{
		  if(othertwo.first < two.second)
		    {
		      update_paired( two.second,paired );
		      update_paired(othertwo.first,paired);
		    }
		  else
		    {
		      update_paired(othertwo.first,paired);
		      update_paired( two.second,paired );
		    }
		}
	      else
		{
		  update_paired( two.second,paired );
		}
	      if(rv2==2)
		update_paired(two.first,paired);
	    }
	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
      else if ( firstlinked || secondlinked )
	{
	  if( firstlinked )
	    {
	      (pbegin+two.second)->second = (pbegin+two.first)->second;
	    }
	  else
	    {
	      (pbegin+two.first)->second = (pbegin+two.second)->second;
	      (pbegin + (pbegin+two.second)->second)->second = two.first;
	      (pbegin +two.second)->second = (pbegin+two.first)->second;
	    }
	  rv = ::coalesce(t,2*nsam,two.first,two.second,nsites,nlinks,sample,sample_history);
	  *NSAM -= rv;
	  config[firstdeme] -= rv;
	  //cerr << "rv = " << rv << '\n';


	  update_paired(two.second,paired);
	  if(rv==2)
	    {
	      (paired->begin()+(paired->begin()+two.first)->second)->second = -1;
	      update_paired(two.first,paired);
	    }

	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
      else if ( !firstlinked && !secondlinked )
	{	  
	  rv = ::coalesce(t,2*nsam,two.first,two.second,nsites,nlinks,sample,sample_history);
	  *NSAM -= rv;
	  config[firstdeme] -= rv;

	  update_paired(two.second,paired);
	  if(rv==2)
	    update_paired(two.first,paired);
	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
    }
  else //chromosomes from different loci
    {
      if( firstlinked && secondlinked )
	{
	  othertwo.first = two.second;
	  othertwo.second = (pbegin+two.first)->second;
	  two.second = (pbegin+othertwo.first)->second;
	  if(two.first>two.second) swap(two.first,two.second);
	  if(othertwo.first>othertwo.second) swap(othertwo.first,othertwo.second);

	  assert( (sbegin+two.first)->pop == (sbegin+two.second)->pop );
	  assert( (sbegin+two.first)->pop == firstdeme );
	  assert( (sbegin+othertwo.first)->pop == (sbegin+othertwo.second)->pop );
	  assert( (sbegin+othertwo.first)->pop == seconddeme );

	  //update linkage
	  (pbegin+two.first)->second = (pbegin+two.second)->second = othertwo.first;
	  (pbegin+othertwo.first)->second = (pbegin+othertwo.second)->second = two.first;

	  if(two.second > othertwo.second)
	    {	 
	      rv = ::coalesce(t,2*nsam,two.first,two.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[firstdeme] -= rv;
	      firstflag=secondflag=0;
	      if(rv == 2 )
		{
		  if(two.first < othertwo.first)
		    {
		      firstflag=secondflag=1;
		    }
		  else if(two.first < othertwo.second)
		    {
		      secondflag=1;
		    }
		}
	      assert ( demes_ok(sample,*NSAM,config) );
	      assert ( (sample->begin()+othertwo.first-firstflag)->pop ==
		       (sample->begin()+othertwo.second-secondflag)->pop );
	      rv2 = ::coalesce(t,2*nsam,othertwo.first-firstflag,
			       othertwo.second-secondflag,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv2;
	      config[seconddeme] -= rv2;
	      assert ( demes_ok(sample,*NSAM,config) );
	      if(rv==2)
		{
		  (paired->begin()+othertwo.first)->second = -1;
		  (paired->begin()+othertwo.second)->second = -1;
		}
	      if(rv2==2)
		{
		  (paired->begin()+two.first)->second = -1;
		  (paired->begin()+two.second)->second = -1;
		}

	      update_paired( two.second,paired );
	      if(rv == 2)
		{
		  if( two.first < othertwo.second )
		    {
		      update_paired( othertwo.second,paired );
		      update_paired(two.first,paired);
		    }
		  else
		    {
		      update_paired(two.first,paired);
		      update_paired( othertwo.second,paired );
		    }
		}
	      else
		{
		  update_paired( othertwo.second,paired );
		}
	      if(rv2==2)
		{
		  update_paired( othertwo.first,paired );
		}
	    }
	  else
	    {
	      rv = ::coalesce(t,2*nsam,othertwo.first,othertwo.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[seconddeme] -= rv;
	      firstflag=secondflag=0;
	      if(rv == 2)
		{
		  if( othertwo.first < two.first )
		    {
		      firstflag=secondflag=1;
		    }
		  else if ( othertwo.first < two.second )
		    {
		      secondflag=1;
		    }
		}
	      assert ( demes_ok(sample,*NSAM,config) );
	      assert ( (sample->begin()+two.first-firstflag)->pop ==
		       (sample->begin()+two.second-secondflag)->pop );
	      rv2 = ::coalesce(t,2*nsam,two.first-firstflag,two.second-secondflag,
			       nsites,nlinks,sample,sample_history);
	      *NSAM -= rv2;
	      config[firstdeme] -= rv2;
	      assert ( demes_ok(sample,*NSAM,config) );
	      if( rv == 2 )
		{
		  (paired->begin()+two.first)->second = -1;
		  (paired->begin()+two.second)->second = -1;
		}
	      if( rv2 == 2 )
		{
		  (paired->begin()+othertwo.first)->second = -1;
		  (paired->begin()+othertwo.second)->second = -1;
		}

	      update_paired( othertwo.second,paired );
	      if( rv == 2 )
		{
		  if(othertwo.first < two.second)
		    {
		      update_paired(two.second,paired);
		      update_paired(othertwo.first,paired);
		    }
		  else
		    {
		      update_paired(othertwo.first,paired);
		      update_paired(two.second,paired);
		    }
		}
	      else
		{
		  update_paired( two.second,paired );
		}
	      if (rv2 == 2)
		{
		  update_paired(two.first,paired);
		}
	    }

	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
      else if ( firstlinked || secondlinked )
	{
	  if ( firstlinked )
	    {
	      ch = two.first;
	      two.first = (pbegin+two.first)->second;
	      if(two.first>two.second) swap(two.first,two.second);
	      
	      assert( (sbegin+two.first)->pop == (sbegin+two.second)-> pop );
	      assert( (sbegin+two.first)->pop == seconddeme );
	    }
	  else
	    {
	      ch = two.second;
	      two.second = (pbegin+two.second)->second;
	      if(two.first>two.second) swap(two.first,two.second);

	      assert( (sbegin+two.first)->pop == (sbegin+two.second)-> pop );
	      assert( (sbegin+two.first)->pop == firstdeme );
	    }

	  //update linkage
	  (pbegin+ch)->second = two.first;
	  (pbegin+two.first)->second = (pbegin+two.second)->second = ch;
	  
	  rv = ::coalesce(t,2*nsam,two.first,two.second,nsites,nlinks,sample,sample_history);
	  *NSAM -= rv;
	  config[(firstlinked) ? seconddeme : firstdeme] -= rv;
	  if(rv == 2)
	    ( paired->begin() + ch )->second = -1;
	  update_paired( two.second,paired );
	  if(rv==2)
	    {
	      update_paired( two.first,paired );
	    }
	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
      else if ( !firstlinked && !secondlinked )
	{
	  (pbegin + two.first)->second = two.second;
	  (pbegin + two.second)->second = two.first;

	  assert ( ! invalid_pairs(paired,sample) );
	  assert ( ! unbalanced_pairs(paired,*NSAM) );
	  assert ( ! redundant_pairs(paired,*NSAM) );
	  assert ( demes_ok(sample,*NSAM,config) );
	}
    }
}
