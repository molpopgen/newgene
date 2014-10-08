#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/PolyTableFunctions.hpp>
//#include <Sequence/Coalescent/SimTypes.hpp>
//#include <Sequence/Coalescent/Initialize.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <limits>
//#include <cmath>
#include <ctime>
#include <set>
#include <boost/bind.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <ca.hpp>
#include <util.hpp>
#include <rec.hpp>
#include <egc.hpp>
#include <fixation.hpp>
#include <constants.hpp>

using namespace std;
using namespace Sequence;

enum EVENT {COAL,REC,CONV,FIXATION};

ostream & operator<<( ostream & o, const EVENT & e)
{
  switch (e)
    {
    case COAL:
      o << "COAL";
      break;
    case REC:
      o << "REC";
      break;
    case CONV:
      o << "CONV";
      break;
    case FIXATION:
      o << "FIXATION";
      break;
    }
  return o;
}

int main(int argc, char **argv)
{
  if ( argc != 11 )
    {
      cerr << "usage:\n"
	   << "dcoal nsam nreps 4Nu 4Nr 4Nc nsites nsites_bw tlen tau seed\n"
	   << "where:\n"
	   << "\tnsam = sample size\n"
	   << "\tnreps = # of simulated samples to generate\n"
	   << "\t4Nu = theta, the population mutation rate\n"
	   << "\t4Nr = rho, the population recombination rate\n"
	   << "\t4Nc = the population ectopic gene conversion rate\n"
	   << "\tnsites = length of locus in base pairs\n"
	   << "\tnsites_bw = # of base pairs between loci\n"
	   << "\ttlen = mean length of an ectopic conversion tract\n"
	   << "\ttau = fixation time of duplicate, in units of 4N generations\n"
	   << "\tseed = random number seed (unsigned integer)\n"
	   << "The output is in \"ms\" format, with the first nsam\n"
	   << "lines corresponding to the ancestral gene, and the next nsam\n"
	   << "lines corresponding to the new gene\n";
      exit(1);
    }
  unsigned argn=1;
  const int nsam = atoi(argv[argn++]);
  unsigned howmany = atoi(argv[argn++]);
  const double theta = atof(argv[argn++]);
  const double rho = atof(argv[argn++]);
  const double econv = atof(argv[argn++]);    //ectopic conversion rate, i.e. 4Nc, where c is the rate per generation
  const int nsites = atoi(argv[argn++]);
  const int nsites_bw = atoi(argv[argn++]);   //distance b/w duplicates, in base pairs, must be >= 0
  const unsigned tlen = atoi(argv[argn++]);   //mean tract length of a conversion event
  const double tau = atof(argv[argn++]);      //time at which duplicate has fixed in the population, units of 4N generations
  const unsigned seed = atoi(argv[argn++]);

  copy(argv,argv+argc,ostream_iterator<char *>(cout, " "));
  cout << endl;
  const double p = 1./double(tlen);
  const double littler = rho/double(nsites-1);
  const double littleec = econv/double(nsites-1);

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  vector<chromosome> isample = init_sample(vector<int>(2,nsam),nsites);
  marginal imarg = init_marginal(2*nsam);
  double rcoal,rrecbw,rgc,tmin,ttemp,t;
  //pair<int,int> two,othertwo;
  int NSAM,ENSAM,NLINKEDPAIRS,nA0,n0B,nAB,nlinks,nlinksbw,nlinksgc,config[2];
  EVENT event;
  vector< pair<int,int> > __paired(2*nsam,make_pair(0,0) );
  typedef vector<pair<int,int> >::iterator vpii_iter;
  vpii_iter vpii1,vpii2;
  for(int i=0;i<nsam;++i)
    {
      __paired[i]=make_pair(i,i+nsam);
    }
  for(int i=nsam;i<2*nsam;++i)
    {
      __paired[i]=make_pair(i,i-nsam);
    }
  assert( ! invalid_pairs(&__paired,&isample) );

  vector<double> trajectory;

  while(howmany--)
    {
      vector<chromosome> sample(isample);
      vector<pair<int,int> > paired(__paired);
      Sequence::arg sample_history(1,imarg);
      NSAM = 2*nsam;
      NLINKEDPAIRS=nsam;
      config[0]=config[1]=nsam;
      ENSAM =  nsam;
      nAB = nsam;
      nA0 = 0;
      n0B = 0;
      nlinksbw = NLINKEDPAIRS*(nsites_bw-1);
      nlinks = NSAM*(nsites-1);
      nlinksgc = nlinks+NSAM; //because every chromo has at least 1 site
      t = 0.;
      bool fixed = false;
      //while(NSAM > 1)
      while(ENSAM > 1 || !fixed)
	{	
	  /*  
#ifndef NDEBUG
	  int nA0 = 0,n0B = 0,nAB = 0;
	  for(int i=0;i<NSAM;++i)
	    {
	      if(sample[i].pop==0)
		{
		  if(paired[i].second == -1)
		    {
		      nA0++;
		    }
		}	      
	      else
		{
		  if(paired[i].second == -1)
		    {
		      n0B++;
		    }
		}
	      nAB = (NSAM-nA0-n0B)/2;
	    }
	  NLINKEDPAIRS=count_linked_pairs(paired);
	  //cerr << "dcoal update: " << nA0 << ' ' << n0B << ' ' << nAB << ' ' << NSAM << ' ' << NLINKEDPAIRS << '\n';
	  if(nAB!=NLINKEDPAIRS)
	    {
	      for(int i=0;i<NSAM;++i)
		{
		  cerr << paired[i].first << " -> " << paired[i].second << '\n';
		}
	    }
	  assert(nAB==NLINKEDPAIRS);
#endif
	  */
	  rcoal = double(ENSAM*(ENSAM-1));
	  rrecbw = littler*nlinksbw;
	  nlinksgc=nlinks+NSAM;
	  int sum=0;
#ifndef NDEBUG
	  for( int i=0;i<NSAM;++i )
	    {
	      sum += sample[i].links();
	      assert( sample[i].last() >= sample[i].first() );
	      for( chromosome::iterator itr = sample[i].begin() ; itr < sample[i].end()-1 ; ++itr)
		{
		  assert(itr->end>=itr->beg);
		  assert((itr+1)->beg > itr->end);
		}
	    }
	  assert(sum==nlinks);
#endif
	  rgc = littleec*double(nlinksgc);

	  tmin = gsl_ran_exponential(r,1/rcoal);
	  event = COAL;
	  ttemp = (rrecbw>0. && t < tau) ? gsl_ran_exponential(r,1/(rrecbw)) : DMAX;
	  assert( ttemp >= 0. );
	  if(ttemp<tmin)
	    {
	      tmin=ttemp;
	      event=REC;
	    }
	  ttemp = (rgc > 0. && t < tau) ? gsl_ran_exponential(r,1/rgc) : DMAX;
	  assert( ttemp >= 0. );
	  if(ttemp<tmin)
	    {
	      tmin=ttemp;
	      event=CONV;
	    }
	  if(ENSAM==1 && !fixed)
	    {
	      //assert(t<tau);
	      tmin = tau-t;
	      event = FIXATION;
	    }
#ifndef NDEBUG
	  else
	    {
	      assert( tmin < DMAX );
	    }
#endif
	  if( (t<tau && t + tmin >=tau) || (t==tau&&tau==0.) )
	    event = FIXATION;
	  //cerr << "next event is " << event << " after duration " << tmin << '\n';
	  if(event == FIXATION)
	    {
	      //enter fixation phase
	      t = tau;

	      ConditionalTrajNeutral( boost::bind(gsl_rng_uniform,r),&trajectory,
				      dtp,1.,0. );

	      fixation(r, &NSAM, &sample, &sample_history,
		       &paired, &nlinks, &NLINKEDPAIRS,
		       &nlinksbw, &t,config, trajectory,
		       pmin,  dt, nsam, nsites,  nsites_bw,
		       littler, littleec,p);
	      fixed=true;
	      NLINKEDPAIRS = count_linked_pairs(paired);
	      ENSAM = config[0]+config[1]-NLINKEDPAIRS;//NLINKEDPAIRS + (config[0]-NLINKEDPAIRS) + (config[1]-NLINKEDPAIRS);
	      nAB = NLINKEDPAIRS;
	      nA0 = config[0] - NLINKEDPAIRS;
	      n0B = config[1] - NLINKEDPAIRS;
	      assert(NLINKEDPAIRS >= 0);
	      nlinksbw = NLINKEDPAIRS*(nsites_bw-1);
	      assert( demes_ok(&sample,NSAM,config) );
	      assert( int(paired.size()) == NSAM );
	    }
	  else
	    {
	      t += tmin;
	      if( event == COAL ) 
		{
		  pair<int,int> two;
		  two = pick2(r,NSAM,paired);
		  ca(&NSAM,&sample,&sample_history,&paired,&nlinks,config,
		     two.first,two.second,t,nsam,nsites);
		  assert(!unbalanced_pairs(&paired,NSAM));
		}
	      else if (event == REC) //recombination between loci
		{
		  rec(r,&paired,NSAM,NLINKEDPAIRS);
		}
	      else //ectopic gene conversion
		{
		  pair<int,int> two = pick_uniform_spot_gc(gsl_rng_uniform(r),nlinksgc,sample.begin(),NSAM);
		  egc(r,two,!(sample[two.first].pop),t,nsam,nsites,nlinksgc,p,&NSAM,&nlinks,config,&sample,&sample_history,&paired);
		  assert(!unbalanced_pairs(&paired,NSAM));
		}
	      NLINKEDPAIRS = count_linked_pairs(paired);
	      ENSAM = config[0]+config[1]-NLINKEDPAIRS;//NLINKEDPAIRS + (config[0]-NLINKEDPAIRS) + (config[1]-NLINKEDPAIRS);
	      nAB = NLINKEDPAIRS;
	      nA0 = config[0] - NLINKEDPAIRS;
	      n0B = config[1] - NLINKEDPAIRS;
	      assert(NLINKEDPAIRS >= 0);
	      nlinksbw = NLINKEDPAIRS*(nsites_bw-1);
	      assert( demes_ok(&sample,NSAM,config) );
	      assert( int(paired.size()) == NSAM );
	    }
	  if ( NSAM < int(sample.size()) / 2 )
	    {
	      sample.erase(sample.begin()+NSAM+1,sample.end());
	    }
	}

      //cerr << t << '\n';
      //cerr  << (sample_history.begin()->begin()+4*nsam-3)->time << '\n';;
      /*
#ifndef NDEBUG

      //make sure first marginal tree is sane
      vector<unsigned> nodecounts(4*nsam-2,0);
      arg::iterator itr = sample_history.begin();
      //cout << t << ' '  << (itr->begin()+4*nsam-3)->time << '\n';;
      //cout << *itr << '\n';
      for( unsigned i=0; i< 4*nsam-2 ; ++i )
	{
	  if( (itr->begin()+i)->abv != -1 )
	    {
	      nodecounts[ (itr->begin()+i)->abv ]++;
	    }
	}
      for(unsigned i=0;i<nodecounts.size();++i)
	{
	  assert(nodecounts[i]<=2);
	}

      //count total time leading to gene 1 and gene 2 in the first marginal tree
      vector<int> desc;
      for( unsigned i=0; i< 4*nsam-2 ; ++i )
        {
	  desc = get_all_descendants( itr->begin(),2*nsam,i );
	  if(desc.size()==nsam)break;
        }
      for(unsigned i=0;i<desc.size();++i)
	cout << desc[i] << ' ';
      cout << '\n';
      exit(1);
#endif
*/
      SimData d = infinite_sites_sim_data(boost::bind(gsl_ran_poisson,r,_1),
					  boost::bind(gsl_ran_flat,r,_1,_2),
					  nsites,sample_history,theta);
      cout << d << endl;
      
    }    
}
