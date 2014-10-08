#include<vector>
#include <iostream>

using namespace std;

int main(int argc,char **argv)
{
  vector<int> v;
  for(unsigned i=0;i<1000000;++i)
    {
      v.insert(v.begin(),i);
      if(i%10000==0.) cerr << v.size()<<'\n';
    }
  cerr << v.size() << '\n';
}
