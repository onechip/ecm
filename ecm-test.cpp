#include <stdio.h>
#include <time.h>

#include "ZZFactoring.h"


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
    cout<<"n="<<n<<endl;
    factor(factorization,n,ZZ::zero(),0,verbose);
    cout<<factorization<<endl;
    mul(nc,factorization);
    if (n!=nc) {
      cout<<"check="<<nc<<endl;
      Error("FAILED FACTORIZATION");
    }
    cout<<endl;
    end = GetTime();
    if (end-prev<min) min=end-prev;
    if (end-prev>max) max=end-prev;
    prev=end;
  }

  // display statistics
  if (min<0) min=0;
  cout<<"min="<<min<<"  avg="<<((end-start)/count)<<"  max="<<max
      <<"  (seconds)"<<endl;

  return 0;
}
