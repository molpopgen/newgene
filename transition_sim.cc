#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

vector<vector<double> > fill_transition_matrix(const unsigned & twoN, 
					       const double & c,
					       const double & r);
unsigned n(const unsigned & i);
bool absorbing( const unsigned & i );
unsigned iterate(gsl_rng * r,unsigned * state,const vector<vector<double> > & tm);

int main( int argc, char ** argv )
{
  int arg = 1;
  const unsigned twoN = atoi(argv[arg++]);
  const double c = atof(argv[arg++]);
  const double r = atof(argv[arg++]);
  unsigned howmany = atoi(argv[arg++]);
  const unsigned seed = atoi(argv[arg++]);

  vector< vector<double> > tm = fill_transition_matrix(twoN,c,r);

  gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,seed);

  unsigned steps,nsam;
  while( howmany-- )
    {
      unsigned state = 3; //initial state for the sample is (2,0,0)
      nsam = n(state);
      unsigned step = 0, tstep = 0;

      while (! absorbing(state) )
	{
	  //	  cerr << state << ' ';
	  steps = iterate(rng,&state,tm);
	  //cerr << state << ' ' << steps << '\n';
	  step += steps;
	  tstep += nsam*steps;
	  nsam = n(state);
	}
      cout << step << ' ' << tstep << endl;
    }
}

vector<vector<double> > fill_transition_matrix(const unsigned & twoN, 
					       const double & c,
					       const double & r)
{
  vector<vector<double> > tm(21,vector<double>(21,0.) );
  const double tc = 2.*c;
  const double thc = 3.*c;
  const double fc = 4.*c;
  const double tr = 2.*r;

  tm[0][0] = 1. - tc - r;
  tm[0][4] = tm[0][5] = c;
  tm[0][8] = r;
  tm[1][1] = 1.;
  tm[2][2] = 1.;
  tm[3][0] = 1./double(twoN);
  tm[3][3] = 1. - tm[3][1] - tr - fc;
  tm[3][9] = tr;
  tm[3][10] = tm[3][11] = tc;
  tm[4][1] = tm[3][0];
  tm[4][4] = 1. - tm[4][1] - tc;
  tm[4][8] = tc;
  tm[5][2] = tm[3][0];
  tm[5][5] = 1. - tm[5][2] - tc;
  tm[5][8] = tc;
  tm[6][0] = tm[3][0];
  tm[6][6] = 1. - tm[6][0] - r - thc;
  tm[6][7] = tm[6][12] = tm[6][14] = c;
  tm[6][13] = r;
  tm[7][0] = tm[3][0];
  tm[7][6] = tm[7][13] = tm[7][15] = c;
  tm[7][7] = 1. - tm[7][0] - thc - r;
  tm[7][12] = r;
  tm[8][0] = tm[3][0];
  tm[8][4] = tm[8][5] = c;
  tm[8][8] = 1. - tm[8][0] - tc;
  tm[9][3] = tm[9][6] = tm[9][7] = tm[3][0];
  tm[9][9] = 1. - 3./double(twoN) - fc - r;
  tm[9][10] = tm[9][11] = tm[9][17] = tm[9][18] = c;
  tm[9][16] = r;
  tm[10][6] = 3./double(twoN);
  tm[10][9] = tc;
  tm[10][10] = 1. - tm[10][6] - fc - r;
  tm[10][16] = tm[10][19] = c;
  tm[10][17] = r;
  tm[11][7] = 3./double(twoN);
  tm[11][8] = tc;
  tm[11][11] = 1.-tm[11][7] - fc - r;
  tm[11][16] = tm[11][20] = c;
  tm[11][18] = r;
  tm[12][7] = 2./double(twoN);
  tm[12][8] = 1./double(twoN);
  tm[12][12] = 1. - tm[12][7] - tm[12][8] - thc;
  tm[12][13] = tc;
  tm[12][14] = c;
  tm[13][6] = 2./double(twoN);
  tm[13][8] = tm[3][0];
  tm[13][12] = tc;
  tm[13][13] = 1.-3./double(twoN)-thc;
  tm[13][14] = c;
  tm[14][4] = 3./double(twoN);
  tm[14][13] = thc;
  tm[14][14] = 1. - tm[14][4] - tm[14][13];
  tm[15][3] = tm[14][4];
  tm[15][12] = tm[14][13];
  tm[15][15] = tm[14][14];
  tm[16][9] = 8./double(twoN);
  tm[16][12] = tm[16][13] = 2./double(twoN);
  tm[16][16] = 1. - 12./double(twoN) - fc;
  tm[16][17] = tm[16][18] = tc;
  tm[17][10] = 4./double(twoN);
  tm[17][13] = 8./double(twoN);
  tm[17][16] = thc;
  tm[17][17] = 1. - 12./double(twoN) - fc;
  tm[17][19] = c;
  tm[18][11] = tm[17][10];
  tm[18][12] = tm[17][13];
  tm[18][16] = thc;
  tm[18][18] = tm[17][17];
  tm[18][20] = c;
  tm[19][14] = 12./double(twoN);
  tm[19][17] = fc;
  tm[19][19] = 1. - tm[19][14] - tm[19][17];
  tm[20][15]=tm[19][14];
  tm[20][18] = tm[19][17];
  tm[20][20] = tm[19][19];

  return tm;
}

unsigned n(const unsigned & i)
{
  if ( i < 3 ) return 1;
  else if ( i < 9 ) return 2;
  else if ( i < 16 ) return 3;
  return 4;
}

bool absorbing( const unsigned & i )
{
  return (i==1 || i == 2);
}

unsigned iterate(gsl_rng * r,unsigned * state,const vector<vector<double> > & tm)
{
  unsigned oldstate = *state;
  //tm[i][i] = P(stay at state i)
  //therefore, P(leave state i) = 1-tm[i][i]
  //the time to leave is therefore geometric with probability 1-tm[i][i]
  double p = 1.-tm[*state][*state];
  unsigned steps = gsl_ran_geometric(r, p );

  //now, figure out which state we jump to
  vector<double>::const_iterator itr,itr2=tm[*state].begin();
  double rdm = gsl_rng_uniform(r),sum=0.;
  //cerr << "current state is " << *state << ' ';
  while( (itr = find_if(itr2,tm[*state].end(),bind2nd(greater<double>(),0.))) < tm[*state].end() )
    {
      unsigned dist = unsigned(std::distance(tm[*state].begin(),itr));
      //cerr << "dist = " << dist << ' ';
      if ( dist == *state ) itr2=itr+1;
      else
	{
	  sum += *itr;
	  //cerr << "( " << sum << "," << p <<") ";
	  if( rdm <= sum/p ) 
	    {
	      assert(tm[*state][dist] > 0.);
	      *state = dist;
	      //cerr << ", changed to : " << dist << " in ";
	      break;
	    }
	  else itr2=itr+1;
	}
    }
  //cerr << steps << ' ' << " steps\n";
#ifndef NDEBUG
  if(*state == oldstate)
    {
      copy(tm[*state].begin(),tm[*state].end(),ostream_iterator<double>(cerr, " "));
      cerr << '\n';
    }
#endif
  assert(*state != oldstate);
  return steps;
}
