#include <Sequence/SimData.hpp>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>

using namespace Sequence;
using namespace std;

int main( int argc, char **argv )
{
  unsigned nsam = atoi(argv[1]);
  

  vector< vector<double> > jsfs(nsam+1, vector<double>(nsam+1,0.));
  unsigned nreps = 0;
  int rv;
  SimData d;
  unsigned sum=0;
  while( (rv=d.fromfile(stdin)) != EOF )
    {
      sum += d.numsites();
      for( SimData::const_site_iterator i = d.sbegin() ; 
	   i != d.send() ; ++i )
	{
	  int c1 = count(i->second.begin(),i->second.begin()+nsam,'1');
	  int c2 = count(i->second.begin()+nsam,i->second.end(),'1');
	  //cerr << c1 << ' ' << c2 << '\n';
	  jsfs[c1][c2] += 1.;
	}
    }
  cout.precision(10);
  for(unsigned i=0;i<nsam+1;++i)
    {
      /*
      copy(jsfs[i].begin(), jsfs[i].end(),
	   ostream_iterator<double>(cerr," "));
      cerr << '\n';
      */
      transform( jsfs[i].begin(), jsfs[i].end(), jsfs[i].begin(),
		 bind2nd(divides<double>(),double(sum)) );
      copy(jsfs[i].begin(), jsfs[i].end(),
	   ostream_iterator<double>(cout," "));
      cout << '\n';
    }
}
