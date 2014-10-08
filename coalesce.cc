#include <Sequence/Coalescent/Coalesce.hpp>
//#include <boost/bind.hpp>
#include <coalesce.hpp>
//#include <util.hpp>
//#include <algorithm>
//#include <gsl/gsl_randist.h>
//#include <iostream>
 
int coalesce(const double & time,
	     const int & ttl_nsam,
	     //const int & current_nsam,
	     const int & c1,
	     const int & c2,
	     const int & nsites,
	     int * nlinks,
	     std::vector<chromosome> * sample,
	     arg * sample_history)
/*!
  This is a modification of the routine in libsequence.  In that library,
  chromosomes are swapped to the end of the array, as is done in Hudson's "ms" program.
  That imposes some mind-bending book-keeping on data structures for this model,
  so we rely on simple erasing here instead.
*/
{
  bool yes1,yes2;
  int ch1=(c1<c2)?c1:c2, ch2=(c2>c1)?c2:c1;

  assert( (sample->begin()+ch1)->nsegs>0 );
  assert( (sample->begin()+ch2)->nsegs>0 );

  std::vector<chromosome>::iterator sbegin=sample->begin();

  chromosome::iterator ch1beg = (sbegin+ch1)->begin(),
    ch2beg=(sbegin+ch2)->begin();
  unsigned seg1=0,seg2=0;

  segment * tsp = (segment *)malloc(sample_history->size()*sizeof(segment));
  int tseg = -1;

  //iterate over marginal histories
  int k=0,nsegs=int(sample_history->size());
  arg::iterator imarg = sample_history->begin(),
    jmarg=imarg;
  jmarg++;
  for ( ; k<nsegs ; ++imarg,++jmarg,++k )
    {
      //ask if chromosomes ch1 and ch2 have segments
      //that are part of the i-th marginal history
      yes1 = isseg(ch1beg,(sbegin+ch1)->nsegs,imarg->beg,&seg1);
      yes2 = isseg(ch2beg,(sbegin+ch2)->nsegs,imarg->beg,&seg2);
      if( yes1 || yes2 )
	{
	  tseg++;
	  (tsp+tseg)->beg = imarg->beg;
	  (tsp+tseg)->end = (k<(nsegs-1)) ? jmarg->beg-1 : nsites-1;

	  if(yes1 && yes2)
	    {
	      imarg->nnodes++;
	      if( imarg->nnodes >= (2*ttl_nsam-2))
		{
		  tseg--;
		}
	      else
		{
		  (tsp+tseg)->desc = imarg->nnodes;
		}
	      marginal::iterator mi = imarg->begin();
	      (mi+imarg->nnodes)->time = time;
	      (mi+(ch1beg+seg1)->desc)->abv = imarg->nnodes;
	      (mi+(ch2beg+seg2)->desc)->abv = imarg->nnodes;
	      assert( (mi+(ch1beg+seg1)->desc)->abv <= int(2*ttl_nsam-2) );
	      assert( (mi+(ch2beg+seg2)->desc)->abv <= int(2*ttl_nsam-2) );
	    }
	  else
	    {
	      (tsp+tseg)->desc = (yes1==true) ? (ch1beg+seg1)->desc : (ch2beg+seg2)->desc;
	      assert( (tsp+tseg)->desc < 2*ttl_nsam-1 );
	    }
	}
    }
  *nlinks -= (sbegin+ch1)->links();
  int flag=0;
  int should_erase = -1;
  if(tseg < 0)
    {
      free(tsp);
      should_erase = ch1;
      flag=1;
      assert( (sbegin+ch1)->nsegs>0 );
      assert( (sbegin+ch2)->nsegs>0 );
    }
  else
    {
      assert( (sbegin+ch1) < sample->end() );
      (sbegin+ch1)->assign_allocated_segs(tsp,tseg+1);
      *nlinks += (sbegin+ch1)->links();
    }

  *nlinks -= (sbegin+ch2)->links();
  if(should_erase != -1)
    {
      sample->erase(sbegin+should_erase);
      sbegin=sample->begin();
    }
  sample->erase(sbegin+ch2-flag);
  return ((tseg<0)?2:1);  
}
