#include <cmath>
#include <iostream>

//using std::cerr;
using namespace std;

double ai(const int &i, const int &n)
/*
  Equation 6 of Tajima 1990
*/
{
  return double(2*(n+1))/double( (i+1)*(i+2)*(n-1) );
}

double ETtot(const int & n)
// Total time on tree, standard coalescent
{
  double t=0.;
  for(int i=1;i<n;++i)
    {
      t+= 1./double(i);
    }
  return t;
}

double ETtotFix(const int &  n)
// Total time on tree, conditional on fixation
{
  double t=0.;
  for(int i=2;i<=n;++i)
    {
      t += 1./double(i+1);
    }
  return t;
}

double ETiFix(const int &i, const int & n)

{
  return( 1./double(i+1) - 1./double(n) );
}

double ETFix(const int & n)
/*
  Unlabeled equation, top or p 499, Tajima 1990.
  Only considering time, not 2Nu
*/
{
  double t=0.;
  for(int i=1;i<n-1;++i)
    {
      t += ai(i,n)*ETiFix(i,n);
    }
  return t;
}

double Dnominator(const int & n, const int & S)
{
  double a1=ETtot(n);
  double a2=0.;
  for(int i=1;i<n;++i)
    {
      a2+=1./double(i*i);
    }
  double b1 = (n+1)/(3*(n-1));
  double b2 = (2.*double(pow(double(n),2.)+n+3))/double(9*n*(n-1));
  double c1 = b1-1./a1;
  double c2 = b2 - double(n+1)/(a1*n) + a2/pow(a1,2);
  double e1 = c1/a1;
  double e2 = c2/(pow(a1,2) + a2);
  return( pow(e1*double(S)+e2*double(S*(S-1)),0.5) );
}

double EDFix(const int &n,const double &theta)
{
  double ES=ETtotFix(n)*theta;
  double EThetaW=ES/ETtot(n);
  double EPi = ETFix(n)*theta;
  return( (EPi-EThetaW)/Dnominator(n,ES) );
}

double lambdai(const int & i)
{
  return double(i*(i+1))/2.;
}

double ziT(const int &i, const int &n, const double & T)
//Equation 20 of Tajima 1990
{
  double rv=0.;//,rv2=0.;
  double li=lambdai(i);
  //double a,b,c,d;
  for(int j=1;j<n;++j)
    {
      if(i!=j)
	{
	  double lj = lambdai(j);
	  /*
	  a = double(2*j+1)*pow(-1.,j+1);
	  b = (lj/(li-lj));
	  c = (exp(-1.*lj*T)-exp(-1.*li*T))*exp(-2.*lj);
	  d= double(2.*i+1)*pow(-1.,i+1)*li*T*exp(-1.*li*(T+2.));
	  */
	  rv += double(2*j+1)*pow(-1.,j+1)*(lj/(li-lj))*
	    (exp(-1.*lj*T)-exp(-1.*li*T))*exp(-2.*lj) +
	    double(2.*i+1)*pow(-1.,i+1)*li*T*exp(-1.*li*(T+2.));
	  
	  //rv2 += exp(log(a)+log(b)+log(c))+d;
	}
    }
  //cerr << rv << ' ' << rv2 << '\n';
  return rv;
}

double yT(const int &n,const double & t)
//Equation 4 of Tajima 1990
{
  double rv=0.,li;
  for(int i=1;i<n;++i)
    {
      li=lambdai(i);
      rv+=double(2*i+1)*pow(-1.,i+1)*li*exp(-1.*li*(t+2.));
    }
  return rv;
}

double ETGivenNumerator(const int &i,const int & n, const double & T)
//Equation 21 of Tajima 1990, numerator only
{
  double t=0.;
  for( int j=i+1;j<n;++j)
    {
      t+=ziT(j,n,T);
    }
  return t;
}

double ETkGivenT(const int &n,const double & T)
//Equation 22 of Tajima 1990, time component only
{
  double rv=0.;
  for(int i=1;i<n-1;++i)
    {
      rv+=ai(i,n)*ETGivenNumerator(i,n,T);
    }
  return rv/yT(n,T);
}
