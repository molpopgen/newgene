#include <Sequence/Coalescent/Coalesce.hpp>
#include <ca.hpp>
#include <coalesce.hpp>
#include <util.hpp>
#include <gsl/gsl_randist.h>
#include <iostream>

void ca(gsl_rng * r, int * NSAM, pair<int,int> * two,vector<chromosome> * sample, arg * sample_history, vector<pair<int,int> > * paired,
	int * nlinks, int config[],
	const double & t, const int & nsam, const int & nsites)
{
  //pair<int,int> two = pick2(r,*NSAM,*paired), othertwo;
  pair<int,int> othertwo;
  int rv;
  vector<chromosome>::iterator sbegin = sample->begin();
  vector<pair<int,int> >::iterator pbegin = paired->begin();
#ifndef NDEBUG
  if( (pbegin+two->first)->second != -1 )
    {
      assert( (sbegin+two->first)->pop != (sbegin+(pbegin+two->first)->second)->pop);
    }
  if( (pbegin+two->second)->second != -1 )
    {
      assert( (sbegin+two->second)->pop != (sbegin+(pbegin+two->second)->second)->pop);
    }
#endif

  int deme1 = (sbegin+two->first)->pop;
  int deme2 = (sbegin+two->second)->pop;
  if( deme1 == deme2 ) //picked 2 of same gene
    {
      //possibilities:
      //a. both genes still linked to something
      //b. only two->first linked to something
      //c. only two->second linked to something
      //d. neither linked to anything
      //if( vpii1->second != -1 && vpii2->second != -1 ) //case a
      if( (pbegin+two->first)->second != -1 && (pbegin+two->second)->second != -1 ) //case a
	{
	  cerr << "A1 ";
	  assert(! unbalanced_pairs(paired,*NSAM));

	  assert( (sbegin+(pbegin+two->first)->second)->pop == !deme1 );
	  assert( (sbegin+(pbegin+two->second)->second)->pop == !deme1 );

	  othertwo = make_pair( min( (pbegin+two->first)->second,(pbegin+two->second)->second ),
				max( (pbegin+two->first)->second,(pbegin+two->second)->second ) );
	  cerr << two->first << '(' << (sbegin+two->first)->pop << ") "<< two->second <<'(' << (sbegin+two->second)->pop << ") " << othertwo.first << '(' << (sbegin+othertwo.first)->pop << ") " << othertwo.second << '(' << (sbegin+othertwo.second)->pop << ") ";	 
	  //just in case
	  assert( othertwo.first != othertwo.second );
	  assert( two->first != othertwo.first );
	  assert( two->first != othertwo.second );
	  assert( two->second != othertwo.first );
	  assert( two->second != othertwo.second );

	  (pbegin + two->first)->second = othertwo.first;
	  (pbegin + two->second)->second = othertwo.first;
	  (pbegin + othertwo.first)->second = two->first;   
	  (pbegin + othertwo.second)->second = two->first;   

	  if ( deme1 == 0  )
	    {
	      rv = ::coalesce(t,2*nsam,*NSAM,othertwo.first,othertwo.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[!deme1] -= rv;
	      cerr << rv << ' ';
	      rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;			    
	      config[deme1] -= rv;
	      cerr << rv << '\n';
	      update_paired( max(othertwo.second,two->second),paired);//update_paired( othertwo.second,paired );
	      update_paired( min(othertwo.second,two->second),paired);//update_paired( two->second,paired );
	      sbegin = sample->begin();
	      pbegin = paired->begin();
	      assert ( ! unbalanced_pairs(paired,*NSAM) );
	      assert ( ! invalid_pairs(paired,sample) );
	      assert ( demes_ok(sample,*NSAM,config) );
	    }
	  else
	    {
	      rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;			    
	      config[deme1] -= rv;		    
	      cerr << rv << ' ';
	      rv = ::coalesce(t,2*nsam,*NSAM,othertwo.first,othertwo.second,nsites,nlinks,sample,sample_history);
	      *NSAM -= rv;
	      config[!deme1] -= rv;
	      cerr << rv << '\n';
	      update_paired( max(othertwo.second,two->second),paired);//update_paired( othertwo.second,paired );
	      update_paired( min(othertwo.second,two->second),paired);//update_paired( two->second,paired );
	      sbegin = sample->begin();
	      pbegin = paired->begin();
	      assert ( ! unbalanced_pairs(paired,*NSAM) );
	      assert ( ! invalid_pairs(paired,sample) );
	      assert ( demes_ok(sample,*NSAM,config) );
	    }
	}
      else if( (pbegin+two->first)->second == -1 && (pbegin+two->second)->second == -1 ) //case d
	{
	  cerr << "A2\n";
	  assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      
	  assert( deme1 == (sbegin+two->second)->pop);	      
	  assert(two->first < two->second);

	  rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	  sbegin=sample->begin();
	  *NSAM -= rv;
	  config[deme1] -= rv;

	  update_paired(two->second,paired);

	  if(rv==2)
	    update_paired(two->first,paired);

	  assert( ! unbalanced_pairs(paired,*NSAM) );
	  assert( ! invalid_pairs(paired,sample) );
	  if(rv==2) cerr << "rv == 2 and it passed?\n";
	  pbegin=paired->begin();
	}
      else 
	{
	  if ( (pbegin+two->first)->second != -1 && (pbegin+two->second)->second == -1) // case b
	    {
	      cerr << "A3\n";
	    }
	  else if ( (pbegin+two->first)->second == -1 && (pbegin+two->second)->second != -1 ) // case c
	    {
	      cerr << "A4\n";
	      //need to update linkage info for pointer 1
	      (pbegin+two->first)->second = (pbegin+two->second)->second;
	      ( pbegin + (pbegin+two->second)->second )->second = two->first;
	    }
	  //need to do the coalsecence of the two chromosomes
	  assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      
	  rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	  sbegin=sample->begin();
	  *NSAM -= rv;
	  config[deme1] -= rv;
	  update_paired(two->second,paired);
	  cerr << "rv = " << rv << '\n';
	  assert( ! unbalanced_pairs(paired,*NSAM) );
	  assert( ! invalid_pairs(paired,sample) );
	  if(rv==2) cerr << "rv == 2 and it passed?\n";
	  pbegin=paired->begin();
	}
    }
  else //picked 1 of gene 1 and 1 of gene 2
    {
      //the possibilities are:
      //a. both genes are still linked to something in the other gene
      //   subcase 1. two genes are linked to each other -- this is an error, and a new pick2 must be written!
      //   subcase 2. two genes linked to different different chromosomes
      //b. only gene 1 linked
      //c. only gene 2 linked
      //d. neither copy is linked to anything
      //if( vpii1->second != -1 && vpii2->second != -1 ) //case a
      if( (pbegin+two->first)->second != -1 && (pbegin+two->second)->second != -1)
	{
	  cerr << "B1\n";

	  othertwo.first = two->second;
	  two->second = (pbegin+othertwo.first)->second;
	  othertwo.second = (pbegin+two->first)->second;
	  deme1 = (sbegin+two->first)->pop;
	  deme2 = (sbegin+othertwo.first)->pop;
	  assert(deme1!=deme2);
	  if(two->first>two->second) swap(two->first,two->second);
	  if(othertwo.first>othertwo.second) swap(othertwo.first,othertwo.second);

	  assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop );
	  assert( (sbegin+othertwo.first)->pop == (sbegin+othertwo.second)->pop );
	  assert( (sbegin+two->first)->pop != (sbegin+othertwo.first)->pop );
	  assert( othertwo.first != othertwo.second );
	  assert( two->first != othertwo.first );
	  assert( two->first != othertwo.second );
	  assert( two->second != othertwo.first );
	  assert( two->second != othertwo.second );

	  //step 1. update the linkage relationships
	  (pbegin+two->first)->second = othertwo.first;
	  (pbegin+othertwo.first)->second = two->first;

	  //step 2. do the coalescences themselves
	  if( othertwo.first > two->second ) 
	    {
	      assert((sbegin+othertwo.first)->pop == (sbegin+othertwo.second)->pop);	      
	      rv = ::coalesce(t,2*nsam,*NSAM,othertwo.first,othertwo.second,nsites,
			      nlinks,sample,sample_history);
	      sbegin=sample->begin();
	      *NSAM -= rv;
	      config[!deme1] -= rv;

	      assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      
	      rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,
			      nlinks,sample,sample_history);
	      sbegin=sample->begin();
	      *NSAM -= rv;
	      config[deme1] -= rv;

	      //step 3. update "paired"
	      cerr << "erase a\n";
	      update_paired(othertwo.second,paired);
	      update_paired(two->second,paired);
	      pbegin=paired->begin();
	      cerr << "rv = " << rv << '\n';
	      assert( ! invalid_pairs(paired,sample) );
	      if(rv==2) cerr << "rv == 2 and it passed?\n";
	    }
	  else
	    {		
	      assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      	      
	      rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,
			      nlinks,sample,sample_history);
	      sbegin=sample->begin();
	      *NSAM -= rv;
	      config[deme1] -= rv;

	      assert( (sbegin+othertwo.first)->pop == (sbegin+othertwo.second)->pop);	      
	      rv = ::coalesce(t,2*nsam,*NSAM,othertwo.first,othertwo.second,nsites,
			      nlinks,sample,sample_history);
	      sbegin=sample->begin();
	      *NSAM -= rv;
	      config[!deme1] -= rv;

	      //step 3. update "paired"	
	      update_paired(two->second,paired);
	      //assert( ! invalid_pairs(paired,sample) );
	      update_paired(othertwo.second,paired);
	      pbegin=paired->begin();
	      assert( ! invalid_pairs(paired,sample) );
	      if(rv==2) cerr << "rv == 2 and it passed?\n";
	    }
	  assert( ! unbalanced_pairs(paired,*NSAM) );  //THIS FAILS (RARE)
	  assert( ! invalid_pairs(paired,sample) );
	}
      else if ( (pbegin+two->first)->second != -1 && (pbegin+two->second)->second == -1 ) //case b
	{
	  cerr << "B2\n";
	  assert(two->first<two->second);
	  bool first = ( two->first < two->second ) ? true : false;
	  int ch = two->first;
	  int deme = (first==true) ? (sbegin+two->first)->pop : (sbegin+two->second)->pop;
	  two->first = (pbegin+ch)->second;
	  if(two->first>two->second) swap(two->first,two->second);
	  (pbegin+ch)->second = two->first;
	  (pbegin+two->first)->second = ch;
	  assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      
	  rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	  sbegin=sample->begin();
	  *NSAM -= rv;
	  config[!deme] -= rv;
	  update_paired(two->second,paired);
	  cerr << "rv = " << rv << '\n';
	  assert( ! invalid_pairs(paired,sample) );
	  pbegin=paired->begin();
	}
      //else if( vpii1->second == -1 && vpii2->second != -1 ) //case c
      else if( (pbegin+two->first)->second == -1 && (pbegin+two->second)->second != -1 ) //case c
	{
	  cerr << "B3\n";
	  int deme = (sbegin+two->first)->pop;//assert( (sbegin+two->first)->pop == 0);
	  //assert( (sbegin+two->second)->pop == 1);
	  int ch = two->second;
	  two->second = (pbegin+ch)->second;
	  if(two->first>two->second) swap(two->first,two->second);
	  (pbegin+ch)->second = two->first;
	  (pbegin+two->first)->second = ch;
	  assert(two->first >=0 && two->first < *NSAM);
	  assert(two->second >= 0 && two->second < *NSAM);
	  assert(ch>=0&&ch < *NSAM);
	  assert( (sbegin+two->first)->pop == (sbegin+two->second)->pop);	      
	  rv = ::coalesce(t,2*nsam,*NSAM,two->first,two->second,nsites,nlinks,sample,sample_history);
	  sbegin=sample->begin();
	  *NSAM -= rv;
	  config[deme] -=rv;//config[0] -= rv;
	  update_paired(two->second,paired);

	  pbegin=paired->begin();
	  cerr << "rv = " << rv << '\n';
	  assert( ! unbalanced_pairs(paired,*NSAM) );
	  assert( ! invalid_pairs(paired,sample) );
	  if(rv==2) cerr << "rv == 2 and it passed?\n";
	}
      else if( (pbegin+two->first)->second != -1 && (pbegin+two->second)->second != -1 ) //case d
	//"heals" a recombination event, no change in sample sizes
	{
	  cerr << "B4\n";
	  (pbegin+two->first)->second = (pbegin+two->second)->first;
	  (pbegin+two->second)->second = (pbegin+two->first)->first;
	  assert( ! unbalanced_pairs(paired,*NSAM) );
	}
    }
}
