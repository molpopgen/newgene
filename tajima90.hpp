#ifndef __TAJIMA90__
#define __TAJIMA90__

double ai(const int &i, const int &n);
double ETtot(const int & n);
double ETtotFix(const int &  n);
double ETiFix(const int &i, const int & n);
double ETFix(const int & n);
double Dnominator(const int & n, const int & S);
double EDFix(const int &n,const double &theta);
double ziT(const int &i, const int &n, const double & T);
double yT(const int &n,const double & t);
double ETGivenNumerator(const int &i,const int & n, const double & T);
double ETkGivenT(const int &n,const double & T);
#endif
