#include <iostream>

#include "ZZFactoring.h"

using namespace NTL;


// simple test program for factoring integers using ECM method
int main(int argc, char* argv[]) {

  if (argc<=1) {
    std::cerr<<"Usage:"<<std::endl
	     <<"  factor number_to_factor"<<std::endl;
    return 1;
  }

  ZZ n;
  conv(n,argv[1]);

  // set to true to get more output
  bool verbose=false;

  ZZ nc;
  vec_pair_ZZ_long factorization;

  double start = GetTime();
  factor(factorization,n,ZZ::zero(),0,verbose);
  double end = GetTime();

  // write out answer
  std::cout<<n<<':';
  for (long i=0; i<factorization.length(); ++i)
    if (factorization[i].b==1)
      std::cout<<' '<<factorization[i].a;
    else
      std::cout<<' '<<factorization[i].a<<'^'<<factorization[i].b;
  std::cout<<std::endl;

  // display statistics
  if (end>start)
    std::cout<<(end-start)<<" seconds"<<std::endl;

  return 0;
}
