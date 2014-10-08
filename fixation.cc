#include<iostream>
#include <Sequence/Coalescent/Coalescent.hpp>
#include <boost/bind.hpp>
#include <fixation.hpp>
#include <ca.hpp>
#include <rec.hpp>
#include <egc.hpp>
#include <util.hpp>
#include <iostream>
#include <numeric>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

//FIX--rates of conversion and recombination are not right!!  Must be adjusted by population size!!
//Quadruple-check rates of coalescence
void fixation(gsl_rng * r, int * NSAM, vector<chromosome> * sample, arg * sample_history,
	      vector<pair<int,int> > * paired, int * nlinks, int *NLINKEDPAIRS,
	      int * nlinksbw, double * t,
	      int config[], const vector<double> & trajectory,
	      const double & pmin, const double & dt, const int & k, const int & N,
	      const int & nsam, const int & nsites, const int & nsites_bw,
	      const double & littler, const double & littleec, const double & p,const bool & unlinked)
{
 
  int gen = 0,nevents,jump,nlinksgc,nlinksgc0=0,nlinksgc1=0;
  double tfixation = 0.;
  pair<int,int> two;
  //cerr << "gen = " << gen << ' ' << trajectory.size() << '\n';
  EVENT event; 
  //while( (*NSAM>1 && !(*NSAM==2 && *NLINKEDPAIRS==1) )
  //&& gen < trajectory.size() && trajectory[gen] >= pmin )
  //vector<double> rates(3);
  double sum,rdm,pr,rcoal0,rcoal1,rcoal01,rcoal,rrecbw,rgc,rgc01,rgc10;
  sum=rrecbw=rcoal=rcoal1=rcoal0=rcoal01=-1.;
  //while( gen < trajectory.size() && trajtory[gen] >= pmin )
  //bool flag=false;
  //double dt2 = 1/double(2*k*N);
  for( int i=0;i<*NSAM;++i )
    {
      if( (sample->begin()+i)->pop == 0 )
	{
	  nlinksgc0 += (sample->begin()+i)->links()+1;
	}
      else
	{
	  nlinksgc1 += (sample->begin()+i)->links()+1;
	}
    }
  //  cerr << *nlinks << ' ' << nlinksgc0 << ' ' << nlinksgc1 << ' ' <<  *NSAM << '\n';
  while( trajectory[gen] >= pmin )
    {
      //cerr << "times: " << gen << ' ' << *t << ' ' << tfixation << '\n';
      //cerr <<"sanity : "<< *NSAM << ' ' << *NLINKEDPAIRS << ' ' << config[0] << ' ' << config[1] << endl;
      assert(config[0]>=0);
      assert(config[1]>=0);
      pr=1.;
      jump=0;
      rdm=gsl_ran_flat(r,0.,1.);
      //cerr << tfixation << ' '<< config[0] << ' ' << config[1] << '\n';
      while( pr > rdm && trajectory[gen]>=pmin )
	{
	  nevents=0;
	  /*
	   rcoal0 = ( config[0] > 1 ) ? 
	    dt*double(config[0]*(config[0]-1)) : 0.;
	   rcoal1 = ( config[1] > 1 ) ? 
	    dt*double(config[1]*(config[1]-1))/trajectory[gen] : 0.;
	  */
	   rcoal0 = ( config[0] > 1 ) ? 
	    dt*double(config[0]*(config[0]-1))/2. : 0.;
	   rcoal1 = ( config[1] > 1 ) ? 
	     dt*double(config[1]*(config[1]-1))/(2.*trajectory[gen]) : 0.;
	   //the factor of 2 below is to put the rate in units of 4N generations
	   //       double rcoal01 = (config[0] && config[1]) ?
	   // 	2.*dt*double(config[0]*config[1]-*NLINKEDPAIRS)/trajectory[gen] : 0.;
	   /*
	     rcoal01 = (unlinked) ? 0 : ((config[0] && config[1]) ?
	     2.*dt*double(config[0]*config[1]-*NLINKEDPAIRS)/trajectory[gen] : 0.);
	   */
	   rcoal01 = (unlinked) ? 0 : ((config[0] && config[1]) ?
				       dt*double(config[0]*config[1]-*NLINKEDPAIRS)/trajectory[gen] : 0.);
	   rcoal = rcoal0+rcoal1+rcoal01;
	   
	   //FIXED--adjust for change in frequency
	   rrecbw = dt*2.*littler*(*nlinksbw);
	   nlinksgc=*nlinks+(*NSAM);
	   assert(nlinksgc == nlinksgc0+nlinksgc1);
	   //rgc = dt*2.*littleec*double(nlinksgc);
	   rgc10 = dt*2.*littleec*double(nlinksgc0);
	   rgc01 = dt*trajectory[gen]*2.*littleec*double(nlinksgc1);
	   rgc=rgc10+rgc01;
	   
	   sum = rcoal + rrecbw + rgc;
	   if(sum>=1) nevents++;
	   pr *= (1.-sum);
	   if(nevents==0)
	     {
	       ++gen;
	       ++jump;
	     }
	}
      if(sum > 0.)
	{
	  double ttemp = double(jump)*dt/2.;
	  *t += ttemp;
	  tfixation += ttemp;//double(jump)*dt2;//ttemp;
	  rdm = gsl_ran_flat(r,0.,1.);
	  if( rdm <= rcoal/sum )
	    {
	      event = COAL;
	    }
	  else if (rdm <= (rcoal+rrecbw)/sum)
	    {
	      event = REC;
	    }
	  else
	    {
	      event = CONV;
	    }
	
	  //cerr << "rates: " << rcoal << ' ' <<  rrecbw << ' ' << rgc << '\n';
	  /*
	    double tmin = gsl_ran_exponential(r,1./rcoal);
	    event = COAL;
	    double ttemp = (rrecbw > 0.) ? gsl_ran_exponential(r,1./(rrecbw)) : DMAX;
	    if(ttemp<tmin)
	    {
	    tmin=ttemp;
	    event = REC;
	    }
	    ttemp =  (rgc > 0.) ? gsl_ran_exponential(r,1./rgc) : DMAX;
	    assert( ttemp >= 0. );
	    if(ttemp<tmin)
	    {
	    tmin=ttemp;
	    event=CONV;
	    }
	    //cerr << "tmin = : " << tmin << ' ' << dt*tmin << '\n';
	    *t += dt*tmin;
	    tfixation += dt*tmin;
	  */
	  if( event == COAL )
	    {
	      two = make_pair(-1,-1);
	      pair<int,int> oldtwo(two);
	      //two.first=two.second=-1;
	      double rdm = gsl_rng_uniform(r);
	      //cerr << "sumcheck " << rdm << ' ' << rcoal0/rcoal << ' ';
	      //Sequence::pick2 in deme is not the right function to use b/c of linkage?
	      //cerr << tfixation << ' ' << config[0] << ' '<<  config[1] << ' ';
	      if( rdm <= rcoal0/rcoal )
		{
		  //coalescence in ancestral gene
		  //cerr << " deme0 ";
		  two = pick2_in_deme( boost::bind(gsl_ran_flat,r,_1,_2), *sample,
				       *NSAM,config[0],0 );
		  //cerr << "0\n";
		  assert( (sample->begin()+two.first)->pop == (sample->begin()+two.second)->pop );
		}
	      else
		{
		  //cerr << (rcoal0+rcoal1)/rcoal << ' ';
		  if (rdm <= (rcoal0+rcoal1)/rcoal)
		    {
		      //cerr << " deme1 ";
		      //coalescence in new copy
		      two = pick2_in_deme( boost::bind(gsl_ran_flat,r,_1,_2), *sample,
					   *NSAM,config[1],1 );
		      assert( (sample->begin()+two.first)->pop == 1 );
		      assert( (sample->begin()+two.first)->pop == (sample->begin()+two.second)->pop );
		      //cerr << "1\n";
		    }
		  else
		    {
		      //cerr << " deme01 ";
		      //coalescence "between copies"
		      //cerr << sum/rcoal << ' ';
		      two = pick2_different_demes(r,*NSAM,config,*sample,*paired);

		      assert( (sample->begin()+two.first)->pop != (sample->begin()+two.second)->pop );
		      //cerr << "2\n";
		      /*
			cerr << "//\n"<<*NSAM << ' ' << config[0] << ' ' << config[1] <<'\n';
			for( unsigned i=0;i<*NSAM;++i )
			{
			cerr << i << ' '
			<< (paired->begin()+i)->first << ' ' << (paired->begin()+i)->second << '\n';
			}
		      */
		    }
		}
	      //cerr << '\n';
	      //cerr << "donesumcheck\n";
	      assert(two!=oldtwo);
	      assert(two.first != -1 );
	      assert( two.first != two.second );
	      /*
		cerr <<"who's who: " << two.first << ' ' << two.second << ' '
		<< (paired->begin()+two.first)->second << ' ' << (paired->begin()+two.second)->second << '\n';
	      */
	      assert( !( (paired->begin()+two.first)->second == two.second 
			 && (paired->begin()+two.second)->second == two.first) );
	      if(unlinked)
		assert( (sample->begin()+two.first)->pop == (sample->begin()+two.second)->pop);
	      ca(NSAM,sample,sample_history,paired,nlinks,config,
		 two.first,two.second,*t,nsam,nsites);
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if ( event == REC )
	    {
	      rec(r,paired,*NSAM,*NLINKEDPAIRS);
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  else if (event == CONV )
	    {	
	      if( gsl_rng_uniform(r) <= rgc01/rgc )
		{
		  //pick an A-bearing chromosome
		  two = pick_uniform_spot_gc_deme( gsl_rng_uniform(r), 
						       nlinksgc0, sample->begin(), *NSAM, 0 );
		  assert( (sample->begin()+two.first)->pop == 0 );
		}
	      else
		{
		  two = pick_uniform_spot_gc_deme( gsl_rng_uniform(r),
						   nlinksgc1, sample->begin(), *NSAM, 1 );
		  assert( (sample->begin()+two.first)->pop == 1 );
		} 
	      egc(r,two,todeme,*t,nsam,nsites,nlinksgc,p,NSAM,nlinks,config,sample,sample_history,paired);
	      assert(!unbalanced_pairs(paired,*NSAM));
	    }
	  //gen = int( tfixation*double(4*k*N) );
	  /*
	    if(gen >= trajectory.size() )
	    {
	    cerr << gen << ' ' << trajectory.size() << '\n';
	    }
	  */
	  //THE NEXT LINES FORCE THE DUPS TO BE UNLINKED
	  /*	  
	  for(unsigned i=0;i<paired->size();++i)
	    (*paired)[i].second=-1;
	  */
	  nlinksgc0=nlinksgc1=0;
	  for( int i=0;i<*NSAM;++i )
	    {
	      if( (sample->begin()+i)->pop == 0 || (sample->begin()+i)->pop == 1)
		{
		  nlinksgc0 += (sample->begin()+i)->links()+1;
		}
	      else
		{
		  nlinksgc1 += (sample->begin()+i)->links()+1;
		}
	    }
	  *NLINKEDPAIRS = count_linked_pairs(*paired);
	  if(unlinked)assert(*NLINKEDPAIRS==0);
	  assert(*NLINKEDPAIRS >= 0);
	  *nlinksbw = *NLINKEDPAIRS*(nsites_bw-1);
	  //cerr << paired->size() << ' ' << *NSAM << ' ' << config[0] << ' ' << config[1] << '\n';
	  assert( demes_ok(sample,*NSAM,config) );
	  assert( int(paired->size()) == *NSAM );
	  //cerr << "gen = " << gen << ' ' << trajectory[gen] << ' ' << trajectory.size() << '\n';
	}
    }
  if (config[1]>1) //guard against v. quick fixations (similar to strong selection fix/hack for sweeps)
    {
      while(config[1]>1)
	{
	  two = pick2_in_deme( boost::bind(gsl_ran_flat,r,_1,_2), *sample,
			       *NSAM,config[1],1 );
	  ca(NSAM,sample,sample_history,paired,nlinks,config,
	     two.first,two.second,*t,nsam,nsites);
	}
    }
  
  //cerr << "done!\n";
  //cerr << "premig: " << config[0]<<' '<<config[1]<<'\n';
  for( int i=0;i<*NSAM;++i )
    {
      (paired->begin()+i)->second = -1;
      if( (sample->begin()+i)->pop==1) 
	{
	  (sample->begin()+i)->pop=0;
	  config[0]++;
	  config[1]--;
	}
    }
  assert(config[1]==0);
  /*
    cerr << "exiting at time " << *t << ' ' 
    << *NSAM << ' '  << config[0] << ' ' << config[1] << '\n';
  */
  //exit(10);
}
