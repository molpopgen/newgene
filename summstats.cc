#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/FST.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <iostream>

using namespace std;
using namespace Sequence;

int main(int argc, char **argv )
{
  if(argc != 2)
    {
      cerr << "usage:\n"
	   << "dcoal or cnvcoal | summstats nsam\n"
	   << "where:\n"
	   << "\tnsam = sample size of the ancestral gene\n";
      exit(1);
    }
  int nsam=atoi(argv[1]);

  SimData d;
  int rv;
  vector<unsigned> config(2,nsam);
  cout << "hbk pib fixed shared priv1 priv2 s1 pi1 d1 h1 s2 pi2 d2 h2\n";
  while( (rv = d.fromfile(stdin)) != EOF )
    {
      SimData d0,d1;
      if(config[1] != d.size()-nsam) config[1]=d.size()-nsam;
      assert(config[0]+config[1] == d.size());
      if(!d.empty())
	{
	  d0.assign(&*d.pbegin(),d.numsites(),&*d.begin(),nsam);
	  d1.assign(&*d.pbegin(),d.numsites(),&*d.begin()+nsam,d.size()-nsam);
	  RemoveInvariantColumns(&d0);
	  RemoveInvariantColumns(&d1);
	}
      FST ad(&d,2,&config[0]);
      PolySIM ad0(&d0),ad1(&d1);
      assert(d0.numsites()==ad0.NumPoly());      
      set<double> shared = ad.shared(0,1),
	fixed = ad.fixed(0,1);
      pair<set<double>, set<double> > priv = ad.Private(0,1);
      
      cout << ad.HBK() << ' ' << ad.piB() << ' '
	   << fixed.size() << ' ' <<  shared.size() << ' '
	   << priv.first.size() << ' ' << priv.second.size() << ' '
	   << ad0.NumPoly() << ' ' << ad0.ThetaPi() << ' ' <<  ad0.TajimasD() << ' ' << ad0.ThetaH() << ' '
	   << ad1.NumPoly() << ' ' << ad1.ThetaPi() << ' ' <<  ad1.TajimasD() << ' ' << ad1.ThetaH() << endl;
    }
}
