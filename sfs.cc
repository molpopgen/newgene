#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/FST.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>
#include <utility>

using namespace std;
using namespace Sequence;

void output_normalized( const vector<unsigned> & v, const unsigned & n )
{
  vector<double> __v(v.begin(),v.end());
  transform(__v.begin(),__v.end(),__v.begin(),
	    bind2nd(divides<double>(),double(n)));
  transform(__v.begin(),__v.end(),__v.begin(),
	    bind2nd(divides<double>(),accumulate(__v.begin(),__v.end(),0.)));
  copy(__v.begin(),__v.end(),
       ostream_iterator<double>(cout," "));
  cout << '\n';
}

typedef set<double>::const_iterator sdci;
typedef SimData::const_site_iterator SDci;

int main(int argc, char **argv )
{
  int nsam=atoi(argv[1]);
  int nsam2 = nsam;
  if(argc==3)
    nsam2=atoi(argv[2]);
  SimData d,d1,d2;
  int rv,c1,c2;
  vector<unsigned> config(1,nsam);
  config.push_back(nsam2);
  vector<double> sharedp(nsam-1,0.),privp(nsam-1,0.);
  double nfixed=0;
  unsigned nreps=0;
  unsigned tmuts = 0;
  while( (rv = d.fromfile(stdin)) != EOF )
    {
      tmuts += d.numsites();
      ++nreps;
      FST fst(&d,2,&config[0]);
      set<double> shared = fst.shared(0,1);
      set<double> fixed = fst.fixed(0,1);
      pair< set<double>,set<double> > priv = fst.Private(0,1);
      /*
      d1.assign(&*d.pbegin(),d.numsites(),&d[0],nsam);
      d2.assign(&*d.pbegin(),d.numsites(),&d[nsam],nsam);
      RemoveInvariantColumns(&d1);
      RemoveInvariantColumns(&d2);
      cerr << "d1.S = " << d1.numsites() << ',' << "d2.S = " << d2.numsites() << '\n';
      cerr << d1 << '\n';
      cerr << fixed.size() << ' ' << shared.size() << ' '
	   << priv.first.size() << ' ' << priv.second.size() << ' ' << d.numsites() << ' '
	   << (fixed.size()+shared.size()+priv.first.size()+priv.second.size()) << '\n';
      */
      nfixed += double(fixed.size())/double(d.numsites());

      for( sdci i = shared.begin() ; i != shared.end() ; ++i )
	{
	  for( SDci j = d.sbegin() ; j != d.send() ; ++j )
	    {
	      if( *i == j->first ) //same position
		{
		  c1 = count(j->second.begin(),j->second.begin()+nsam,'1');
		  c2 = count(j->second.begin()+nsam,j->second.end(),'1');
		  if(c1>0)
		    {
		      sharedp[c1-1] += 1./double(d.numsites());
		    }
		  if(c2>0)
		    {
		      sharedp[c2-1] += 1./double(d.numsites());
		    }
		}
	    }
	}

      for( sdci i = priv.first.begin() ; i != priv.first.end() ; ++i )
	{
	  for( SDci j = d.sbegin() ; j != d.send() ; ++j )
	    {
	      if( *i == j->first ) //same position
		{
		  c1 = count(j->second.begin(),j->second.begin()+nsam,'1');
		  if(c1>0)
		    {
		      privp[c1-1] += 1./double(d.numsites());
		    }
		}
	    }
	}

      for( sdci i = priv.second.begin() ; i != priv.second.end() ; ++i )
	{
	  for( SDci j = d.sbegin() ; j != d.send() ; ++j )
	    {
	      if( *i == j->first ) //same position
		{
		  c1 = count(j->second.begin()+nsam,j->second.end(),'1');
		  if(c1>0)
		    {
		      privp[c1-1] += 1./double(d.numsites());
		    }
		}
	    }
	}
    }
  transform(sharedp.begin(),sharedp.end(),sharedp.begin(),
	    bind2nd(divides<double>(),double(nreps)));
  transform(privp.begin(),privp.end(),privp.begin(),
	    bind2nd(divides<double>(),double(nreps)));

  cout << nfixed/double(nreps) << ' ';
  copy(sharedp.begin(),sharedp.end(),
       ostream_iterator<double>(cout," "));
  copy(privp.begin(),privp.end(),
       ostream_iterator<double>(cout," "));
  cout << '\n';
}
