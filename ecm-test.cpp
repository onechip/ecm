#include <iostream>
#include <time.h>
#include "ZZFactoring.h"

using namespace NTL;

// multiply out factorization to get n
void mul(ZZ& n, const vec_pair_ZZ_long& factorization) {
  n=1;
  for (long i=0; i<factorization.length(); ++i)
    n*=power(factorization[i].a,factorization[i].b);
}


// simple test program for factoring integers using ECM method
int main(int argc, char* argv[]) {

  // initialize random number generator
  SetSeed(to_ZZ(time(0)));

  // number of random integers to factor
  long count=16;

  // size of integers (in bits)
  long bits=128;

  // set to true to get more output
  bool verbose=false;

  ZZ n,nc;
  vec_pair_ZZ_long factorization;
  double min = HUGE_VAL;
  double max = 0;
  double start = GetTime();
  double prev = start;
  double end = start;

  // factor count integers
  for (long i=0; i<count; ++i) {
    RandomBits(n,bits);
    std::cout<<"n="<<n<<std::endl;
    factor(factorization,n,ZZ::zero(),0,verbose);
    std::cout<<factorization<<std::endl;
    mul(nc,factorization);
    if (n!=nc) {
      std::cout<<"check="<<nc<<std::endl;
      Error("FAILED FACTORIZATION");
    }
    std::cout<<std::endl;
    end = GetTime();
    if (end-prev<min) min=end-prev;
    if (end-prev>max) max=end-prev;
    prev=end;
  }

  // display statistics
  if (min<0) min=0;
  std::cout<<"min="<<min<<"  avg="<<((end-start)/count)<<"  max="<<max
	   <<"  (seconds)"<<std::endl;

  return 0;
}

