#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <limits>
#include <ctime>
#include <set>
#include <boost/bind.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;


int main(int argc, char **argv)
{
  unsigned argn=1;
  const int nsam = atoi(argv[argn++]);
  unsigned howmany = atoi(argv[argn++]);
  const double theta = atof(argv[argn++]);
  const double tau = atof(argv[argn++]);      //time at which duplicate has fixed in the population, units of 4N generations
  const unsigned seed = atoi(argv[argn++]);

  copy(argv,argv+argc,ostream_iterator<char *>(cout, " "));
  cout << endl;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  vector<chromosome> isample = init_sample(vector<int>(1,nsam),1);
  marginal imarg = init_marginal(nsam);
  double rcoal;
  //pair<int,int> two,othertwo;
  int NSAM;
  pair<int,int> two;
  while(howmany--)
    {
      int nlinks = 0;
      vector<chromosome> sample(isample);
      Sequence::arg sample_history(1,imarg);
      vector<double> trajectory;
      NSAM = nsam;
      double t = 0.;
      bool fixed = false;
      while(NSAM > 1)
	{	
	  rcoal = double(NSAM*(NSAM-1));
	  double tcoal = gsl_ran_exponential(r,1./rcoal);
	  
	  if( t <= tau && t+tcoal >= tau )
	    {
	      t=tau;
	      const int N = 10000;
	      const int k = 25;
	      const double dtp = 1./double(2*k*N);
	      ConditionalTrajNeutral( boost::bind(gsl_rng_uniform,r),&trajectory,
				      dtp,1.,0. );
	      
	      const double dt = 1./double(4*k*N);
	      const double pmin = 1./double(2*N);

	      int gen = 0,jump,nevents;
	      double pcoal;
	      while(trajectory[gen] >= pmin && NSAM > 1)
		{
		  //cerr << gen << ' ' << trajectory.size() << '\n';
		  double pr=1.;
		  jump=0;
		  double rdm = gsl_rng_uniform(r);
		  while( pr>rdm && trajectory[gen]>=pmin )
		    {
		      /*
		      cerr << gen << ' ' << trajectory.size() << ' ' 
			   << pr << ' ' << rdm << ' ' << NSAM << '\n';
		      */
		      nevents=0;
		      pcoal = dt*double(NSAM*(NSAM-1))/trajectory[gen];
		      pr *= (1.-pcoal);
		      if(pcoal>=1.)nevents++;
		      if(!nevents)
			{
			  ++gen;
			  ++jump;
			}
		    }
		  if(pcoal>0.)
		    {
		      t += double(jump)*dt;
		      two = pick2(boost::bind(gsl_ran_flat,r,_1,_2),NSAM);
		      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,
				       1,&nlinks,&sample,&sample_history);
		    }
		}
	    }
	  else
	    {
	      t += tcoal;
	      two = pick2(boost::bind(gsl_ran_flat,r,_1,_2),NSAM);
	      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,
			       1,&nlinks,&sample,&sample_history);
	    }
	  //cerr << "tmrca = " << t << '\n';
	  if ( NSAM < int(sample.size()) / 2 )
	    {
	      sample.erase(sample.begin()+NSAM+1,sample.end());
	    }
	}
      SimData d = infinite_sites_sim_data(boost::bind(gsl_ran_poisson,r,_1),
					  boost::bind(gsl_ran_flat,r,_1,_2),
					  1,sample_history,theta);
      cout << d << endl;      
    }    
}
