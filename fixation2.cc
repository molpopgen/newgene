#include<iostream>
#include <Sequence/Coalescent/Coalescent.hpp>
#include <boost/bind.hpp>
#include <fixation.hpp>
#include <ca.hpp>
#include <rec.hpp>
#include <egc.hpp>
#include <util.hpp>
#include <kahan.hpp>
#include <numeric>
#include <cmath>

#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

enum EVENT {COALAM,COALNOTAM,RECABTOAM,RECABTOA,RECATOAM,RECAMTOA,CONVAB,CONVBA,CONVBAM,NEVENTS};

void fixation(gsl_rng * r, int * NSAM, vector<chromosome> * sample, Sequence::arg * sample_history,
	      vector<pair<int,int> > * paired, int * nlinks, int *NLINKEDPAIRS,
	      int * nlinksbw, double * t,
	      int config[], const vector<double> & trajectory,
	      const double & pmin, const double & dt,
	      const int & nsam, const int & nsites, const int & nsites_bw,
	      const double & littler, const double & littleec, const double & p)
{ 
  int gen = 0,nevents,jump,nlinksgc,nlinksgc0=0,nlinksgc1=0;
  double tfixation = 0.;
  pair<int,int> two;
  EVENT event=NEVENTS; 
  double rdm,pr,rcoal0,rcoal1,rcoal,rrecbw,rgc=-1.,rgc01=-1.,rgc10=-1.;
  //  sum=
  rrecbw=rcoal=rcoal1=rcoal0=-1.;
  kahan_sum<double> sum;
  vector<chromosome>::iterator sbegin = sample->begin();
  vector<pair<int,int> >::iterator pbegin = paired->begin();
  for( int i=0;i<*NSAM;++i )
    {
      if( (sbegin+i)->pop == 0 )
	{
	  nlinksgc0 += (sbegin+i)->links()+1;
	}
      else
	{
	  nlinksgc1 += (sbegin+i)->links()+1;
	}
    }

  int config_s[3];
  config_s[0]=config_s[1]=config_s[2]=0;

  if(trajectory[0]<1.)
    /*
      CNV
    */
    {
      for(int i=0;i<*NSAM;++i)
	{
	  if( (sbegin+i)->pop == 0 )
	    {
	      if( (pbegin+i)->second == -1 )
		{
		  (sbegin+i)->pop = 0;
		  config_s[0]++;	  
		}
	      else
		{
		  (sbegin+i)->pop = 1;
		  config_s[1]++;
		}
	    }
	  else if ( (sbegin+i)->pop == 1)
	    {
	      (sbegin+i)->pop=2;
	      config_s[2]++;
	    }
	}
    }
  else
    /*all chromosomes are linked to the duplicated background
      at the start:
      deme 0 is empty
      deme 1 is A- + the A component of AB pairs
      deme 2 is -B + the B component of AB pairs
    */
    {
      for(int i=0;i<*NSAM;++i)
	{ 
	  if( (sbegin + i)->pop == 0 )
	    {
	      (sbegin+i)->pop=1;
	      config_s[1]++;
	    }
	  else if ( (sbegin+i)->pop == 1)
	    {
	      (sbegin+i)->pop=2;
	      config_s[2]++;
	    }
	}
    }

  int ENSAM =config_s[0]+config_s[1]+config_s[2]-*NLINKEDPAIRS;
  int nAplus = config_s[1]-*NLINKEDPAIRS;

  double rates[NEVENTS];
  while( trajectory[gen] >= pmin )
    {
      pr=1.;
      jump=0;
      rdm=gsl_ran_flat(r,0.,1.);
      ENSAM = config_s[0]+config_s[1]+config_s[2]-*NLINKEDPAIRS;
      while( pr > rdm && trajectory[gen]>=pmin )
      {
	  nevents=0;
  //coalescence in the ancestral gene, ancestry unlinked to fixation event
	  rates[COALAM] = (config_s[0]>1) ? dt*double(config_s[0]*(config_s[0]-1))/(1.-trajectory[gen])  : 0.;
	  //coalescence amongst chromosomes linked to fixation event
	  rates[COALNOTAM] = (ENSAM-config_s[0] > 1) ? dt*double( (ENSAM-config_s[0])*(ENSAM-config_s[0]-1) )/trajectory[gen] : 0.;
	  //cerr << trajectory[gen] << ' ' << rates[COALAM] << ' ' << rates[COALNOTAM] << '\n';
	  nAplus = config_s[1]-*NLINKEDPAIRS;
	  assert(nAplus >= 0);

	  //Recombination in an AB pair, A has ancestry in portion of population unlinked to fixation
	  rates[RECABTOAM] = dt*littler*double(nsites_bw-1)*double(*NLINKEDPAIRS)*(1.-trajectory[gen]);
	  //Recombination in an AB pair, A has ancestry in portion of population linked to fixation
	  rates[RECABTOA] = dt*littler*double(nsites_bw-1)*double(*NLINKEDPAIRS)*trajectory[gen];
	  //recombination from A- background onto A+
	  rates[RECAMTOA] = dt*littler*double(nsites_bw-1)*double(config_s[0])*trajectory[gen];
	  //recombination from A+ background onto A-
	  rates[RECATOAM] = dt*littler*double(nsites_bw-1)*double(config_s[1]-*NLINKEDPAIRS)*(1.-trajectory[gen]);
	  nlinksgc=*nlinks+(*NSAM);
	  //conversion from gene A to gene B
	  rates[CONVAB] = dt*littleec*double(nlinksgc0)*trajectory[gen];
	  //conversion from duplicate to ancestral gene, ancestry of ancestral gene not linked to fixation
	  rates[CONVBA] = dt*littleec*double(nlinksgc1)*trajectory[gen];
	  //conversion from duplicate to ancestral gene, ancestry of ancestral gene is linked to fixation
	  rates[CONVBAM] = dt*littleec*double(nlinksgc1)*(1.-trajectory[gen]);
	  assert(nlinksgc == nlinksgc0+nlinksgc1);

	  sum = kahan_sum<double>();
	  for(unsigned i=0;i<int(NEVENTS);++i)
	    sum+=rates[i];
	  if(sum>=1.) nevents++;
	  pr *= (1.-sum);
	  if(nevents==0)
	    {
	      ++gen;
	      ++jump;
	    }
      }
      if(sum>0.)
	{
          double ttemp = double(jump)*dt;
          *t += ttemp;
          tfixation += ttemp;
	  event=NEVENTS;
	  kahan_sum<double> tsum;
	  rdm=gsl_rng_uniform(r);
	  for(int i = 0 ; i < NEVENTS ; ++i)
	    {
	      tsum += rates[i]/sum;
	      if( rdm <= tsum && rates[i]>0.)
		{
		  event = EVENT(i);
		  break;
		}
	    }
	  assert(event != NEVENTS);

	  if( event == COALAM )
	    {
	      two = pick2_in_deme( boost::bind(gsl_ran_flat,r,_1,_2), *sample,
				   *NSAM,config_s[0],0 );
	      assert ( (sbegin+two.first)->pop == 0 );
	      assert( (sbegin+two.first)->pop == (sbegin+two.second)->pop );
	      assert( two.first != two.second );
	      assert( !( (pbegin+two.first)->second == two.second 
			 && (pbegin+two.second)->second == two.first) );
	      ca(NSAM,sample,sample_history,paired,nlinks,config_s,
		 two.first,two.second,*t,nsam,nsites);
	      //The ca function can result in containers being resized,
	      //which will invalidate iterators, so we need to re-assign
	      sbegin=sample->begin();
	      pbegin=paired->begin();
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if (event == COALNOTAM)
	    {
	      two = pick2_not_in_deme( r, sbegin,pbegin,
				       *NSAM,config_s[1]+config_s[2],0 );
	      assert( (sbegin+two.first)->pop != 0 );
	      assert( (sbegin+two.second)->pop != 0 );

	      assert( two.first != two.second );
	      assert( !( (pbegin+two.first)->second == two.second 
			 && (pbegin+two.second)->second == two.first) );
	      ca(NSAM,sample,sample_history,paired,nlinks,config_s,
		 two.first,two.second,*t,nsam,nsites);
	      //The ca function can result in containers being resized,
	      //which will invalidate iterators, so we need to re-assign
	      sbegin=sample->begin();
	      pbegin=paired->begin();
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if ( event == RECABTOAM || event == RECABTOA )
	    {
	      //recombination in an AB just updates linkage relationships,
	      //and does not change sample size
	      int chindex = int(gsl_ran_flat(r,0.,double(*NLINKEDPAIRS)));
	      int dummy=-1,ch=-1;
	      for( vector<pair<int,int> >::iterator itr = pbegin ;
		   itr != (pbegin+*NSAM) ; ++itr )
		{
		  if(itr->first != -1 && itr->second != -1)
		    {
		      ch = int(itr-pbegin);
		      ++dummy;
		      if(dummy==chindex) break;
		    }
		}
	      assert( (sbegin+ch)->pop == 1 || (sbegin+ch)->pop == 2);

	      if (event == RECABTOAM)
		{
		  if( (sbegin+ch)->pop== 1 )
		    {
		      (sbegin+ch)->pop=0;
		    }
		  else if ((sbegin+ch)->pop == 2 )
		    {
		      (sbegin + (pbegin+ch)->second)->pop=0;
		    }
		  else
		    {
		      assert(false);
		    }
		  config_s[0]++;
		  config_s[1]--;
		}
	      (pbegin+(pbegin+ch)->second)->second = -1;
	      (pbegin+ch)->second = -1;
	    }
	  else if( event == RECATOAM )
	    {
	      int chindex;
	      int dummy=-1,ch = -1;
	      chindex = int(gsl_ran_flat(r,0.,double(nAplus)));
	      vector<chromosome>::iterator citr = sbegin;
	      for( vector<pair<int,int> >::iterator itr = pbegin ;
		   itr != (pbegin+*NSAM) ; ++itr,++citr )
		{
		  if(citr->pop==1 && itr->second == -1)
		    {
		      ch = int(itr-pbegin);
		      ++dummy;
		      if(dummy==chindex) break;
		    }
		}
	      assert(ch>=0);
	      assert( (sbegin+ch)->pop == 1 );
	      assert( (pbegin+ch)->second == -1 );
	      (sbegin+ch)->pop=0;
	      config_s[1]--;
	      config_s[0]++;
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if (event == RECAMTOA)
	    {
	      assert(config_s[0]>0);
	      int chindex = int(gsl_ran_flat(r,0.,double(config_s[0])));
	      int dummy=-1,ch=-1;
	      vector<chromosome>::iterator citr = sbegin;
	      for( vector<pair<int,int> >::iterator itr = pbegin ;
		   itr != (pbegin+*NSAM) ; ++itr,++citr )
		{
		  if(citr->pop==0 && itr->second == -1)
		    {
		      ch = int(itr-pbegin);
		      ++dummy;
		      if(dummy==chindex) break;
		    }
		}
	      assert(ch>=0);
	      assert( (sbegin+ch)->pop == 0 );
	      assert( (pbegin+ch)->second == -1 );
	      (sbegin+ch)->pop=1;
	      config_s[1]++;
	      config_s[0]--;
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if (event == CONVAB)
	    {	
	      two = pick_uniform_spot_gc_not_deme( gsl_rng_uniform(r), 
						   nlinksgc0, sbegin, *NSAM, 2 );
	      assert( (sbegin+two.first)->pop != 2 );
	      egc(r,two,2,*t,nsam,nsites,nlinksgc,p,NSAM,nlinks,config_s,sample,sample_history,paired);
	      sbegin=sample->begin();
	      pbegin=paired->begin();
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if (event == CONVBA || event == CONVBAM )//duplicate to ancestor
	    {
	      two = pick_uniform_spot_gc_deme( gsl_rng_uniform(r),
					       nlinksgc1, sbegin, *NSAM, 2 );
	      assert( (sbegin+two.first)->pop == 2 );
	      int todeme = (event==CONVBA) ? 1 : 0;
	      egc(r,two,todeme,*t,nsam,nsites,nlinksgc,p,NSAM,nlinks,config_s,sample,sample_history,paired);
	      sbegin=sample->begin();
	      pbegin=paired->begin();
	      assert(!unbalanced_pairs(paired,*NSAM));
	    } 
	  else 
	    {
	      assert(false);
	    }
#ifndef NDEBUG
	  int c1=0,c2=0,c3=0,pop;
	  for( int i = 0 ; i <  *NSAM ; ++i )
	    {
	      pop=(sbegin+i)->pop;
	      switch(pop)
		{
		case 0:
		  c1++;
		  break;
		case 1:
		  c2++;
		  break;
		case 2:
		  c3++;
		  break;
		}
	    }
	  assert(c1 == config_s[0]);
	  assert(c2 == config_s[1]);
	  assert(c3 == config_s[2]);
	  assert(c1+c2+c3 == *NSAM);
#endif
	  nlinksgc0=nlinksgc1=0;
	  for( int i=0;i<*NSAM;++i )
	    {
	      if( (sbegin+i)->pop != 2 )
		{
		  nlinksgc0 += (sbegin+i)->links()+1;
		}
	      else
		{
		  nlinksgc1 += (sbegin+i)->links()+1;
		}
	    }
	  *NLINKEDPAIRS = count_linked_pairs(*paired);
	  assert(*NLINKEDPAIRS >= 0);
	  *nlinksbw = *NLINKEDPAIRS*(nsites_bw-1);
	  assert( demes_ok(sample,*NSAM,config_s) );
	  assert( int(paired->size()) == *NSAM );
	  assert( config_s[0]+config_s[1]+config_s[2] == *NSAM );
	  }
    }

  //one final update of the sample sizes, etc.
  *NLINKEDPAIRS = count_linked_pairs(*paired);
  *nlinksbw = *NLINKEDPAIRS*(nsites_bw-1);
  ENSAM =config_s[0]+config_s[1]+config_s[2]-*NLINKEDPAIRS;

  //A guard against quick fixation, similar to what happens if a sweep is strong.
  //Testing suggests this is rare (<= 1/10000 replicates)
  while( ENSAM - config_s[0] > 1 )
    {
      two = pick2_not_in_deme( r, sbegin,pbegin,
			       *NSAM,config_s[1]+config_s[2],0 );
      ca(NSAM,sample,sample_history,paired,nlinks,config_s,
	 two.first,two.second,*t,nsam,nsites);
      *NLINKEDPAIRS = count_linked_pairs(*paired);
      ENSAM = config_s[0]+config_s[1]+config_s[2]-*NLINKEDPAIRS;
    }

  for( int i=0;i<*NSAM;++i )
    {
      (pbegin+i)->second = -1;
      (sbegin+i)->pop=0;
    }
  config[0]=*NSAM;
  config[1]=0;
}
