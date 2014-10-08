#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

static const int N = 10000;
//static const int k = 25;
static const int k = 2;
static const double dtp = 1./double(2*k*N);
static const double dt = 1./double(4*k*N);
static const double pmin = 1./double(2*N);

#endif
