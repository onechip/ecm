#ifndef NTL_ZZFactoring__H
#define NTL_ZZFactoring__H

#include <NTL/ZZ.h>
#include "pair_ZZ_long.h"

NTL_OPEN_NNS;

/*
Copyright (C) 2003 Chris Studholme

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/



/**************************************************************************\

                                 Factoring

\**************************************************************************/


/* General purpose integer factoring using best known methods (namely,
 * Lenstra's elliptic curve method).  If N<0, the absolute value of N will 
 * be factored.  If N==0 or N==1, factors will be set to zero length.
 *
 * This method is intended to be used for factoring or testing the smoothness
 * of composites with no particular structure.  Although this method is capable
 * of factoring a structured composite (such as an RSA modulus), other methods
 * may be more efficient.
 *
 * When finished, the factors vector will contain the factors and exponents 
 * sorted in increasing order by factor.  
 * 
 * If bnd>0, the method only runs long enough to find prime factors less than
 * bnd.  In this case, the method may inadvertantly find larger prime factors;
 * however, it may also miss some prime factors less than bnd.  The
 * failure_prob parameter indicates the desired probability of failing to find
 * factors less than bnd.  Smaller values of failure_prob will result in longer
 * running times.  
 *
 * If bnd==0, the method always runs until the complete prime factorization
 * of N has been found (ie. failure_prob has no effect).  Large primes in the 
 * factorization are tested using ProbPrime().
 */
void factor(vec_pair_ZZ_long& factors, const ZZ& N,
            const ZZ& bnd=ZZ::zero(), double failure_prob=0.1, 
	    bool verbose=false);


/* Elliptic curve method primitive.  The parameters are:
 *
 *   N>1          - the composite to factor (relatively prime to 6)
 *   NumBits(bnd) - the size of the factor sought (0 for SqrRoot(N))
 *   failure_prob - probability of failing to find desired factor
 *                  (1 indicates that only one curve is tested,
 *                   0 indicates run until factor is found)
 *
 * If N is prime, an attempt to factor N is still made and this method may
 * run forever (especially if failure_prob=0).  If bnd is set less than 
 * SqrRoot(N) and failure_prob=0, then the method will run until some 
 * non-trivial factor of N is found; furthermore, if N has no factors < bnd,
 * the method may take a very long time as the relavent algorithm parameters
 * have been optimized for finding factors < bnd.
 *
 * The result is q which will either be a non-trivial factor of N,
 * or 1 if no such factor was found.  q may be larger than bnd.  failure_prob
 * indicates that if N has a factor < bnd, the probability of finding a non-
 * trivial factor of N is at least 1-failure_prob.
 */
void ECM(ZZ& q, const ZZ& N, const ZZ& bnd=ZZ::zero(), 
	 double failure_prob=0, bool verbose=false);


/* Pollard Rho factoring primitive.  The parameters are:
 *
 *   N>1          - the number to factor
 *   NumBits(bnd) - the size of the factor sought (0 for SqrRoot(N))
 *   failure_prob - probability of failing to find desired factor
 *                  (0 indicates run until a factor is found)
 *
 * If N is prime, an attempt to factor N is still made and this method may
 * run forever (especially if failure_prob=0).  
 *
 * The result is q which will either be a non-trivial factor of N,
 * or 1 if no such factor was found.  q may be larger than bnd.  failure_prob
 * indicates that if N has a factor < bnd, the probability of finding a non-
 * trivial factor of N is at least 1-failure_prob.
 */
void PollardRho(ZZ& q, const ZZ& n, const ZZ& bnd=ZZ::zero(), 
		double failure_prob=0, bool verbose=false);


/* Trial division factoring primitive.  The parameters are:
 *
 *   N>1 - the number to factor
 *   bnd - maximum prime to check (need not be prime)
 *
 * The result is prime factors in factors and any remaining large factor (may
 * be prime) in q (or 1 if there is no remaining factor).
 *
 * If bnd==0, an appropriate bound is chosen as a function of N.
 *
 * All factors returned in factors are guarenteed to be prime.
 *
 * If n<=(bnd^2), factors will contain the complete prime factorization of N
 * and q=1.
 */
void TrialDivision(vec_pair_ZZ_long& factors, ZZ& q,
		   const ZZ& N, long bnd=0);



/**************************************************************************\

                             Primality Proving

\**************************************************************************/

/* Prove that n is prime.  If n<0, then the absolute value of n is tested.
 * If n is 0 or 1, then 0 is returned.
 *
 * Not implemented; just calls ProbPrime() for now.
 */
bool ProvePrime(const ZZ& n);




NTL_CLOSE_NNS;

#endif
