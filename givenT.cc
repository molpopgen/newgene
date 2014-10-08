#include <tajima90.hpp>
#include <iostream>

using namespace std;

int main( int argc, char **argv)
{
 //  const int n = atoi(argv[1]);
//   const double T = atof(argv[2]);
//   const double theta = atof(argv[3]);
//   double EPi = ETkGivenT(n,T)*theta;
  int twoN=20;
  for(double i=0.01;i<150;i+=0.1)
    {
      cout << i << ' ' << ETkGivenT(twoN,i) << endl;
    }
}
