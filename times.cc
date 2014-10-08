//compare fixation times to coalescence times
#include <Sequence/Coalescent/Trajectories.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <boost/bind.hpp>
using namespace std;
using namespace Sequence;


int main(int argc, char **argv)
{
  unsigned seed = atoi(argv[1]);
  vector<double> trajectory;
  int N = 10000;
  int k = 25;
  double dtp = 1./double(2*k*N);
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  for(unsigned i=0;i<1000;++i)
    {
      ConditionalTrajNeutral( boost::bind(gsl_rng_uniform,r),&trajectory,
				      dtp,1.,0. );
      double tmrca=0.;
      for(int j=2;j<=50;++j)
	{
	  tmrca += gsl_ran_exponential(r,1./(double(j*(j-1))));
	}
      cout << tmrca << ' ' << trajectory.size()/double(4*k*N) << '\n';
    }
}
