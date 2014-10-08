//validation that the "erasing" version of the coalescent function is not bogus

#include <Sequence/Coalescent/Coalescent.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <ctime>
#include <set>
#include <boost/bind.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

  int coalesce(const double & time,
	       const int & ttl_nsam,
	       const int & current_nsam,
	       const int & c1,
	       const int & c2,
	       const int & nsites,
	       int * nlinks,
	       std::vector<chromosome> * sample,
	       Sequence::arg * sample_history)
  /*!
    @brief Common ancestor routine for coalescent simulation.  Merges
    chromosome segments and updates marginal trees.

    Common ancestor routine for coalescent simulation. This routine performs
    the merging of two lineages by a coalescent event.  Such merges usually 
    require two sorts of operations.  The first is an update to the segments
    contained in a chromosome, and the second is an update of the nodes
    on a marginal tree.
    \param time the time at which the coalecent event is occuring
    \param ttl_nsam the total sample size being simulated
    \param current_nsam the current sample size in the simulation
    \param c1 the array index of the first chromosome involved in the coalescent event
    \param c2 the array index of the second chromosome involved in the coalescent event
    \param nsites the total mutational length of the region begin simulated.  In the
    language of Hudson (1983), this is the number of infinitely-many-alleles loci in the
    simulation.
    \param nlinks a pointer to the number of "links" currently in the simulation.
    A link is the region between two sites, such that a chromosome currently with k sites has 
    k-1 links
    \param sample a pointer to the vector of chromosomes which makes up the sample
    \param sample_history a pointer to the ancestral recombination graph
    \return the decrease in current_nsam due to the coalescent event.  Usually, 
    the return value is 1.  Sometimes, however, it is two, when the two chromosomes
    being merged have no ancestral material on the same marginal tree.
    \ingroup coalescent
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
    Sequence::arg::iterator imarg = sample_history->begin(),
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
    int should_erase = -1;//,should_swap2=-1;
    if(tseg < 0)
      {
	free(tsp);
	//(sbegin+ch1)->swap_with(*(sbegin+current_nsam-1));
	should_erase = ch1;
	/*
	if(ch2 == current_nsam-1)
	  {
	    cerr << "setting ch2=ch1\n";
	    ch2=ch1;
	  }
	*/
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
	//cerr << "needing to erase\n";
	sample->erase(sbegin+should_erase);
	sbegin=sample->begin();
      }
    //cerr << "erasing ch2\n";
    sample->erase(sbegin+ch2-flag);//current_nsam-1-flag);
    //(sbegin+ch2)->swap_with(*(sbegin+current_nsam-1-flag));
    return ((tseg<0)?2:1);  
  }

int main(int argc, char **argv)
{
  unsigned seed = atoi(argv[1]);
  int nsam = 30;
  double rho = 25;
  int nsites = 1000;
  int howmany = 1000000;
  double theta = 25;
  double littler = rho/double(nsites-1);
  vector<chromosome> isample = init_sample(vector<int>(1,30),nsites);
  marginal imarg = init_marginal(nsam);
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  pair<int,int> two;
  while(howmany--)
    {
      vector<chromosome> sample(isample);
      Sequence::arg sample_history(1,imarg);
      int NSAM=nsam,nlinks=NSAM*(nsites-1),dummy;
      double t=0;
      while(NSAM>1)
	{
	  dummy=0;
	  for(unsigned i=0;i<NSAM;++i)
	    {
	      dummy+=sample[i].links();
	    }
	  assert(dummy==nlinks);
	  double tcoal = (NSAM>1) ? gsl_ran_exponential(r,1./(double(NSAM)*double(NSAM-1))) : 1000.;
	  double trec = (rho>0.) ? gsl_ran_exponential(r,1./(littler*double(nlinks))):1000.;

	  if(tcoal<trec)
	    {
	      t+=tcoal;
	      two = pick2(boost::bind(gsl_ran_flat,r,_1,_2),NSAM);
	      NSAM -= Sequence::coalesce(t,nsam,NSAM,two.first,two.second,nsites,&nlinks,&sample,&sample_history);
	      //NSAM -= ::coalesce(t,nsam,NSAM,two.first,two.second,nsites,&nlinks,&sample,&sample_history);
	    }
	  else
	    {
	      t+=trec;
	      two = pick_uniform_spot(gsl_rng_uniform(r),nlinks,sample.begin(),NSAM);
	      nlinks -= crossover(NSAM,two.first,two.second,&sample,&sample_history);
	      NSAM++;
	    }
	}
      SimData d = infinite_sites_sim_data(boost::bind(gsl_ran_poisson,r,_1),
					  boost::bind(gsl_ran_flat,r,_1,_2),
					  nsites,sample_history,theta);
      cout << d << endl;
    }
}
