#include <cmath>

#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_long.h>
#include "ZZFactoring.h"
#include "EC_p.h"

NTL_START_IMPL;

/* Implementation of factoring methods.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */


// factors is kept sorted by p
void addFactor(vec_pair_ZZ_long& factors, const ZZ& p, long exponent=1) {
  // fast path: factors.length()==0
  if (factors.length()==0) {
    factors.SetLength(1);
    factors[0].a = p;
    factors[0].b = exponent;
    return;
  }

  // fast path: p>=factors[factors.length()-1].a
  if (p>=factors[factors.length()-1].a) {
    if (p==factors[factors.length()-1].a)
      factors[factors.length()-1].b += exponent;
    else {
      factors.SetLength(factors.length()+1);
      factors[factors.length()-1].a = p;
      factors[factors.length()-1].b = exponent;
    }
    return;
  }

  // binary search to find location to insert
  long low=0;
  long high=factors.length();
  while (high>low) {
    long mid = (low+high)/2;
    if (factors[mid].a<p)
      low=mid+1;
    else // (p<=factors[mid].a)
      high=mid;
  }
  if ((low<factors.length())&&(factors[low].a==p))
    factors[low].b += exponent;
  else {
    // insert factor
    factors.SetLength(factors.length()+1);
    for (long i=factors.length()-1; i>low; --i)
      factors[i] = factors[i-1];
    factors[low].a = p;
    factors[low].b = exponent;
  }
}



/****************  PollardRho  ****************/

// factor n into a*b using Pollard Rho method
// pre-condition: n>1
inline void PollardRho(ZZ& a, ZZ& b, const ZZ& n, 
		       const ZZ& bnd=ZZ::zero()) {
  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(n);

  ZZ d;
  ZZ_p x1;
  random(x1);
  ZZ_p x2(x1);

  ZZ end(IsZero(bnd)?5*SqrRoot(SqrRoot(n)):2*SqrRoot(bnd));
  for (; !IsZero(end); --end) {
    x1 = x1*x1 + 1;
    x2 = x2*x2 + 1;
    x2 = x2*x2 + 1;
    GCD(d,n,rep(x2-x1));
    if ((d>1)&&(d<n)) {
      a=d; b=n/d;
      return;
    }
  }
  // failure
  a=1; b=n;
}

void PollardRho(ZZ& q, const ZZ& n,
		const ZZ& bnd, double failure_prob,
		bool verbose) {
  Error("PollardRho: not implemented");
}

// recursively factor n and add results to factors
// pre-condition: factors is of zero-length or 
//                a properly sorted list of factors
// pre-condition: n>1
/*
void PollardRho(vec_pair_ZZ_long& factors, const ZZ& n,
		const ZZ& _bnd=ZZ::zero(), 
		long deterministic=0, long verbose=0) {
  bool pn = deterministic ? ProvePrime(n) : ProbPrime(n);
  if (pn) {
    addFactor(factors,n);
    return;
  }

  ZZ bnd(n<=_bnd ? ZZ::zero() : _bnd);

  ZZ a,b;
  do {
    PollardRho(a,b,n,bnd);
    if (!IsOne(a)&&!IsOne(b))
      break;
    if (!deterministic||!IsZero(bnd)) {
      addFactor(factors,n);
      return;
    }
  } while (true);

  PollardRho(factors,a,bnd,deterministic,verbose);
  PollardRho(factors,b,bnd,deterministic,verbose);
}
*/



/****************  ECM  ****************/

// get optimal parameters
void ECM_parameters(long& B1, long& B2, double& prob, long& D,
		    long pbits, long nbits) {
  // smoothness bounds
  double logp = pbits*M_LN2;
  double b1 = 1.25*exp(0.77*sqrt(logp)*sqrt(log(logp)));
  B1 = b1<NTL_SP_BOUND ? (long)b1 : NTL_SP_BOUND;
  double b2 = 70*exp(0.78*sqrt(logp)*sqrt(log(logp)));
  B2 = b2<NTL_SP_BOUND ? (long)b2 : NTL_SP_BOUND;

  // memory use in stage 2 is 3*D*nbits/8
  D = (long)sqrt((B2-B1)/2.0);
  if (2*D>=B1) D = (B1-1)/2;
  if (D<2) D=2;

  // all curve orders are divisible by 12
  logp-=log(12.0); 

  // probability that curve order is B1-smooth
  double logB1 = log((double)B1);
  double u = logp/logB1;
  prob = std::pow(u,-u);
  // add in probability that order has one larger factor (up to B2)
  // numerical integration: time is O(log(B2/B1))
  long min_x = B1;
  double min_y = std::pow(u-1,1-u)/logB1/B1;
  while (min_x<B2) {
    long max_x = 2*min_x;
    if (max_x>B2 || min_x>max_x) max_x=B2;
    double logx = log((double)max_x);
    double v = logx/logB1;
    double max_y = std::pow(u-v,v-u)/logx/max_x;
    prob += (max_x-min_x)*(min_y+max_y)/2;
    min_x = max_x;
    min_y = max_y;
  }
  if (prob>1) prob=1;
}

// returns a non-trivial factor q of ZZ_p::modulus(), or 1 otherwise
void ECM_stage_one(ZZ& q, EC_p& Q, PrimeSeq& seq, long bound) {
  long sbound = (long)sqrt((double)bound);
  seq.reset(0);
  long p = seq.next();
  for (; p<=sbound; p=seq.next()) {
    long pp,t=p;
    do { pp=t; t*=p; } while (t>pp && t<=bound); // we might overflow t here
    mul(Q,Q,pp);
  }
  for (; p<=bound; p=seq.next())
    mul(Q,Q,p);
  if (!IsZero(Q))
    GCD(q,ZZ_p::modulus(),rep(Q.Z));
  else
    set(q);
}

// returns a non-trivial factor q of ZZ_p::modulus(), or 1 otherwise
void ECM_stage_two(ZZ& q, EC_p& Q, PrimeSeq& seq, long bound, long D) {
  long B1 = seq.next();
  if (B1<=2*D) {
    // no primes to work with
    set(q);
    return;
  }

  EC_p R;
  mul(R,B1,Q);
  // check R for divisor
  if (IsZero(R)) {
    set(q);
    return;
  }
  GCD(q,ZZ_p::modulus(),rep(R.Z));
  if (!IsOne(q))
    return;
  EC_p T;
  mul(T,B1-2*D,Q);

  // compute point multiples S[d]=2*d*Q
  EC_p S[D+1];
  S[0]=Q;
  doub(S[1],Q);
  doub(S[2],S[1]);
  for (long d=3; d<=D; ++d)
    addh(S[d],S[d-1],S[1],S[d-2]);
  ZZ_p beta[D+1];
  for (long d=0; d<=D; ++d)
    mul(beta[d],S[d].X,S[d].Z);

  ZZ_p g,t,t2;
  set(g);
  ZZ_p alpha;
  long r=B1;
  long p=seq.next();
  do {
    mul(alpha,R.X,R.Z);
    do {
      long delta = (p-r)/2;
      if (delta>D) break;
      //g *= (R.X-S[delta].X) * (R.Z+S[delta].Z) - alpha + beta[delta];
      sub(t,R.X,S[delta].X);
      add(t2,R.Z,S[delta].Z);
      t*=t2;
      t-=alpha;
      t+=beta[delta];
      g*=t;
      // next prime
      p = seq.next();
      if (p==0) {
	// ran out of primes (should never happen)
	p=NTL_MAX_LONG;
	break;
      }
    } while (p<=bound);
    if (p>bound)
      break;
    addh(T,R,S[D],T);
    swap(R,T);
    r+=2*D;
  } while (true);
  if (!IsZero(g))
    GCD(q,ZZ_p::modulus(),rep(g));
  else
    set(q);
}


// choose random curve
// if IsZero(Z), then X is non-trivial factor of ZZ_p::modulus()
void ECM_random_curve(EC_pCurve& curve, ZZ_p& X, ZZ_p& Z) {
  ZZ sigma;
  RandomBnd(sigma,ZZ_p::modulus()-6);
  sigma+=6;
  ZZ_p u,v;
  u = sqr(to_ZZ_p(sigma))-5;
  v = 4*to_ZZ_p(sigma);
  ZZ_p C,Cd;
  C = (v-u)*sqr(v-u)*(3*u+v);
  Cd = 4*u*sqr(u)*v;
  // make sure Cd is invertible
  ZZ Cinv;
  if (InvModStatus(Cinv,rep(Cd),ZZ_p::modulus())!=0) {
    conv(X,Cinv);
    clear(Z);
    return;
  }
  C*=to_ZZ_p(Cinv);
  C-=2;

  // random curve
  ZZ_pX f;
  SetCoeff(f,3);
  SetCoeff(f,2,C);
  SetCoeff(f,1);
  conv(curve,f);
  curve.SetRepresentation(curve.MONTGOMERY);

  // initial point
  mul(X,u,sqr(u));
  mul(Z,v,sqr(v));
}


// test one curve (q is factor or 1 in case of failure)
void ECM_one_curve(ZZ& q, PrimeSeq& seq, long B1, long B2, long D) {
  EC_pBak bak;
  bak.save();

  // random curve and point
  EC_pCurve curve;
  ZZ_p X,Z;
  ECM_random_curve(curve,X,Z);
  if (IsZero(Z)) {
    q=rep(X);
    return;
  }
  EC_p::init(curve);

  // initial point
  EC_p Q;
  Q.X = X;
  Q.Z = Z;

  // attempt to factor
  ECM_stage_one(q,Q,seq,B1);
  if (IsOne(q)&&!IsZero(Q))
    ECM_stage_two(q,Q,seq,B2,D);
}


void ECM(ZZ& q, const ZZ& N, long ncurves, long B1, long B2, long D, 
	 bool verbose) {
  // initialize ZZ_p
  ZZ_pBak bak;
  bak.save();
  ZZ_p::init(N);

  PrimeSeq seq;

  for (long i=0; i<ncurves; ++i) {
    if (verbose) {
      if (ncurves<NTL_MAX_LONG)
	cout<<"ECM: "<<NumBits(N)<<" bits; "
	    <<i<<"/"<<ncurves<<" curves\r"<<flush;
      else
	cout<<"ECM: "<<NumBits(N)<<" bits; "
	    <<i<<" curves\r"<<flush;
    }

    // try to find a factor
    ECM_one_curve(q,seq,B1,B2,D);
    if (!IsOne(q)) {
      if (verbose) {
	cout<<"ECM: "<<NumBits(N)<<" bits; "
	    <<(i+1)<<" curves; found "<<q<<"  "<<endl;
      }
      return;
    }
  }
  if (verbose)
    cout<<"ECM: "<<NumBits(N)<<" bits; FAILED!        "<<endl;
  set(q);
}


// ECM primitive
void ECM(ZZ& q, const ZZ& n, const ZZ& _bnd, double failure_prob,
	 bool verbose) {
  ZZ bnd;
  SqrRoot(bnd,n);
  if (_bnd>0 && _bnd<bnd)
    bnd=_bnd;

  long B1,B2,D;
  double prob;
  ECM_parameters(B1,B2,prob,D,NumBits(bnd),NumBits(n));
  if (verbose)
    cout<<"ECM: bnd="<<NumBits(bnd)<<"bits  B1="<<B1<<"  B2="<<B2
	<<"  D="<<D<<"  prob="<<prob<<endl;
  
  if (failure_prob<=0) {
    // run "forever"
    if (verbose)
      cout<<"ECM: expected ncurves="<<(1/prob)<<endl;
    ECM(q,n,NTL_MAX_LONG,B1,B2,D,verbose); 
  }
  else if (failure_prob>=1) {
    // run on one curve
    ECM(q,n,1,B1,B2,D,verbose); 
  }
  else {
    // run just the right number of times
    // failure_prob = (1-prob)^ncurves;
    double ncurves_d;
    if (prob<0.1)
      ncurves_d = to_double(log(failure_prob)/log1p(to_RR(-prob)));
    else
      ncurves_d = log(failure_prob)/log(1-prob);
    long ncurves = ncurves_d<=1 ? 1 :
      (ncurves_d>=NTL_MAX_LONG ? NTL_MAX_LONG : (long)ncurves_d);
    if (verbose)
      cout<<"ECM: good_prob="<<prob<<"  failure_prob="<<failure_prob
	  <<"  ncurves="<<ncurves<<endl;
    ECM(q,n,ncurves,B1,B2,D,verbose); 
  }
}


/****************  ProbPrime and power detection  ****************/

// same as ProbPrime(), but with no trial division stage
bool ProbPrime_notd(const ZZ& n, long NumTrials=10) {
  // first try W == 2....the exponentiation
  // algorithm runs slightly faster in this case
  ZZ W;
  W = 2;
  if (MillerWitness(n, W))
    return false;

  while (--NumTrials>=0) {
    do {
      RandomBnd(W,n);
    } while (IsZero(W));
    // W == 0 is not a useful candidate for a witness!
    if (MillerWitness(n, W))
      return false;
  }
  return true;
}

/* Quickly find factor of n in the case where n is a prime power (and not
 * a power of 2).  n>2, odd.  Returns q, a factor of n, or 1 in case of
 * failure.
 */
void PrimePowerTest(ZZ& q, const ZZ& n) {
  // n-1 == 2^k * m, m odd
  ZZ n1,m;
  sub(n1,n,1);
  m=n1;
  long k = MakeOdd(m);
  set(q);
  
  ZZ z;
  conv(z,2);
  PowerMod(z,z,m,n);
  do {
    if (IsOne(z) || z==n1) return;
    SqrMod(z,z,n);
  } while (--k>0);
  if (IsOne(z)) return;
  z*=2;
  z-=2;
  GCD(q,z,n);
  if (q==0 || q==n)
    set(q);
}



/****************  TrialDivision  ****************/

// trial division primitive
void TrialDivision(vec_pair_ZZ_long& factors, ZZ& q, const ZZ& n, long bnd) {
  factors.SetLength(0);

  if (&q!=&n) q=n;
  if (bnd==0) {
    bnd=10000;  // should probably be higher
  }

  PrimeSeq s;
  ZZ d;
  for (long p=s.next(); (p>0 && p<=bnd); p=s.next()) {
    if (DivRem(d,q,p)==0) {
      long e=1;
      q=d;
      while (DivRem(d,q,p)==0) {
	++e;
	q=d;
      }
      addFactor(factors,to_ZZ(p),e);
      if (IsOne(q))
	return;
    }
    if (d<=p) {
      // q must be prime
      addFactor(factors,q);
      set(q);
      return;
    }
  }
}


/****************  general factoring  ****************/

// recursive factoring (precondition: n composite)
// currently only uses ECM
void factor_r(vec_pair_ZZ_long& factors, const ZZ& _n,
	      const ZZ& bnd, double failure_prob, bool verbose) {
  ZZ q;
  ZZ n(_n);
  do {
    // attempt to factor n
    if (bnd>0)
      ECM(q,n,bnd,failure_prob,verbose);
    else
      ECM(q,n,ZZ::zero(),0,verbose);

    if (IsOne(q)) {
      // give up
      addFactor(factors,n);
      return;
    }

    // compute other factor
    div(n,n,q);
    if (n<q) swap(n,q);

    // q is small factor, n is large factor
    if (ProbPrime_notd(q))
      addFactor(factors,q);
    else 
      factor_r(factors,q,bnd,failure_prob,verbose);

    // check if n is still composite
    if (ProbPrime_notd(n)) {
      addFactor(factors,n);
      return;
    }
  } while (true);
}


// general purpose factoring method
void factor(vec_pair_ZZ_long& factors, const ZZ& _n,
            const ZZ& bnd, double failure_prob, 
	    bool verbose) {
  ZZ n(_n);
  if (n<=1) {
    abs(n,n);
    if (n<=1) {
      factors.SetLength(0);
      return;
    }
  }

  // upper bound on size of smallest prime factor
  ZZ upper_bound;
  SqrRoot(upper_bound,n);
  if (bnd>0 && bnd<upper_bound)
    upper_bound=bnd;

  // figure out appropriate lower_bound for trial division
  long B1,B2,D;
  double prob;
  ECM_parameters(B1,B2,prob,D,NumBits(upper_bound),NumBits(n));
  ZZ lower_bound;
  conv(lower_bound,max(B2,1<<14));
  if (lower_bound>upper_bound)
    lower_bound=upper_bound;
  
  // start factoring with trial division
  TrialDivision(factors,n,n,to_long(lower_bound));
  if (IsOne(n))
    return;
  if (upper_bound<=lower_bound || ProbPrime_notd(n)) {
    addFactor(factors,n);
    return;
  }

  /* n is composite and smallest prime factor is assumed to be such that
   *     lower_bound < factor <= upper_bound
   *
   * Ramp-up to searching for factors of size upper_bound.  This is a good
   * idea in cases where we have no idea what size factors N might have,
   * but we don't want to spend too much time doing this.
   */
  
  for(lower_bound<<=4; lower_bound<upper_bound; lower_bound<<=4) {
    ZZ q;
    ECM(q,n,lower_bound,1,verbose); // one curve only
    if (!IsOne(q)) {
      div(n,n,q);
      if (n<q) swap(n,q);
      // q is small factor, n is large factor
      if (ProbPrime_notd(q))
	addFactor(factors,q);
      else
	factor_r(factors,q,bnd,failure_prob,verbose);
      if (ProbPrime_notd(n)) {
	addFactor(factors,n);
	return;
      }
      // new upper_bound
      SqrRoot(upper_bound,n);
      if (bnd>0 && bnd<upper_bound)
	upper_bound=bnd;
    }
  }

  // search for factors of size bnd
  factor_r(factors,n,bnd,failure_prob,verbose);
}



/****************  ProvePrime  ****************/

bool ProvePrime(const ZZ& _n) {
  ZZ n(_n);
  if (n<0)
    abs(n,n);
  if (n<=1)
    return 0;

  if (n<=1000000) {
    // n is small so use trial division to check primality
    long ln = to_long(n);
    long end = to_long(SqrRoot(n));
    PrimeSeq s;
    for (long p=s.next(); p<=end; p=s.next())
      if ((ln%p)==0)
	return 0;
    return 1;
  }

  // check small primes
  PrimeSeq s;
  for (long p=s.next(); p<1000; p=s.next())
    if (divide(n,p))
      return 0;

  // obviously, something is missing here!

  return ProbPrime(n);
}

NTL_END_IMPL;
