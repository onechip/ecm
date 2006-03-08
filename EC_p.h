#ifndef NTL_EC_p__H
#define NTL_EC_p__H

#include <NTL/ZZ_pX.h>
#include <ostream>

/*
Copyright (C) 2006 Chris Studholme

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

*****************************************************

MODULE: EC_p

SUMMARY:

  Elliptic curves over ZZ_p.

*/


#define NTLX_OPT_RETURN(T,x) return x

NTL_OPEN_NNS;

// classes defined here
class EC_pCurve;
class EC_p;
class EC_pBak;


/* Class representing an elliptic curve defined by the equation:
 *
 *   y^2 = f(x)
 *
 * where f(x) is a polynomial in x.  Currently, only f(x) monic and degree 3
 * is supported.
 *
 * Points can have a variety of representations (a_i are coefficients of f):
 *
 *   Affine coordinates:
 *     Y^2 = X^3 + a_2*X^2 + a_1*X + a_0   
 *       Z = 1
 *
 *   Projective coordinates:
 *     Y^2*Z = X^3 + a_2*X^2*Z + a_1*X*Z^2 + a_0*Z^3
 *
 *   Jacobian (modified) projective coordinates:
 *     Y^2 = X^3 + a_2*X^2*Z^2 + a_1*X*Z^4 + a_0*Z^6
 *
 *   Montgomery coordinates:
 *     same as projective but without Y coordinate
 *
 * The default is projective coordinates.
 */
class EC_pCurve { 
public:
  // options for point representation
  static const long AFFINE;
  static const long PROJECTIVE;
  static const long JACOBIAN;
  static const long MONTGOMERY;

  friend class EC_p;
  friend void conv(EC_pCurve&, const ZZ_pX&);

public:  // should be protected
  // point representation
  long rep;
  // curve parameters
  ZZ_pX f;

public:
  // creates trivial (singular) curve f(x)=0
  // with projective representation of points
  EC_pCurve() {
    rep=PROJECTIVE;
    SetCoeff(f,3);
  }

  // from polynomial
  EC_pCurve(const ZZ_pX& f) {
    conv(*this,f);
  }

  // copy constructor
  EC_pCurve(const EC_pCurve& other) {
    rep=other.rep;
    f=other.f;
  }

  // destructor
  ~EC_pCurve() {}

  // assignment
  inline EC_pCurve& operator=(const EC_pCurve& other) {
    rep=other.rep;
    f=other.f;
    return *this;
  }

  // set modified projective coordinates
  inline void SetRepresentation(long rep) {
    this->rep = rep;
  }

  // get representation
  inline long representation() const {
    return rep;
  }

  // polynomial for y^2
  inline const ZZ_pX& y2() const {
    return f;
  }
};

// resets point representation to projective coordinates
void conv(EC_pCurve& curve, const ZZ_pX& f);

// test if curve is singular
bool IsSingular(const EC_pCurve& curve);

// change curve to equivalent curve, but with a_2=0
void MakeStandard(EC_pCurve& curve);

// access to coefficients of y^2 (ie. f)
inline const ZZ_p& coeff(const EC_pCurve& curve, long i) {
  return coeff(curve.y2(),i);
}

// set coefficient of y^2 (ie. f)
inline void SetCoeff(EC_pCurve& curve, long i, const ZZ_p& a) {
  if (i<0 || i>2)
    Error("SetCoeff: index out of range");
  SetCoeff(curve.f,i,a);
}

/* Class to represent a point on an elliptic curve over ZZ_p.
 */
class EC_p {
protected:
  static EC_pCurve* c;
  static EC_p* inf;

  friend class EC_pBak;

public:
  // coordinates of point (see c.rep for point representation)
  ZZ_p X,Y,Z;

public:
  // initialize to point at infinity
  EC_p();

  // copy constructor
  EC_p(const EC_p& other) : 
    X(other.X),Y(other.Y),Z(other.Z) {}

  /* Initialize to 0 and swap with other. 
   */ 
  EC_p(EC_p& other, NTL_NNS INIT_TRANS_TYPE) :
    X(other.X,INIT_TRANS),Y(other.Y,INIT_TRANS),Z(other.Z,INIT_TRANS) {}

  // destructor
  ~EC_p() {}

  // assignment
  EC_p& operator=(const EC_p& a) {
    X=a.X; Y=a.Y; Z=a.Z;
    return *this;
  }

  // compute affine coordinates
  void affine(ZZ_p& x, ZZ_p& y) const;

  // compute affine coordinate (x only)
  void affine(ZZ_p& x) const;

  // adjust so Z=1
  void makeAffine();

  // test if this is a point on the curve
  bool isValid() const;

  // set curve to work with
  static void init(const EC_pCurve& curve) {
    if (inf) delete inf;
    if (c) delete c;
    c = new EC_pCurve(curve);
    inf = new EC_p();
  }

  // read-only reference to current curve
  inline static const EC_pCurve& curve() {
    return *c;
  }

  // read-only reference to point at infinity
  inline static const EC_p& zero() {
    return *inf;
  }

  // read-only reference to point at infinity
  inline static const EC_p& infinity() {
    return *inf;
  }
};


// tools
void swap(EC_p& a, EC_p& b);

void clear(EC_p& x);

inline bool IsZero(const EC_p& x) {
  return IsZero(x.Z);
}

inline bool IsInfinity(const EC_p& x) {
  return IsZero(x.Z);
}

inline bool IsValid(const EC_p& x) {
  return x.isValid();
}

inline void MakeAffine(EC_p& x) {
  x.makeAffine();
}

bool operator==(const EC_p& a, const EC_p& b);
inline bool operator!=(const EC_p& a, const EC_p& b) {
  return !(a==b);
}


// stream operators
std::ostream& operator<<(std::ostream&, const EC_pCurve&);
std::ostream& operator<<(std::ostream&, const EC_p&);


// random point
void random(EC_p& p);


// arithmetic
void negate(EC_p& x, const EC_p& a);
void add(EC_p& x, const EC_p& a, const EC_p& b);
void sub(EC_p& x, const EC_p& a, const EC_p& b);
void doub(EC_p& x, const EC_p& a);
inline EC_p doub(const EC_p& a) {
  EC_p x;  doub(x,a);  NTLX_OPT_RETURN(EC_p,x);
}

// for Montgomery representation: computes x=a+b where c=a-b
void addh(EC_p& x, const EC_p& a, const EC_p& b, const EC_p& c);

inline EC_p operator+(const EC_p& a, const EC_p& b) {
  EC_p x;  add(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p operator-(const EC_p& a, const EC_p& b) {
  EC_p x;  sub(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p operator-(const EC_p& a) {
  EC_p x;  negate(x,a);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p& operator+=(EC_p& x, const EC_p& a) { 
  add(x,x,a);  return x;
}
inline EC_p& operator-=(EC_p& x, const EC_p& a) { 
  sub(x,x,a);  return x;
}


// point multiplication
void mul(EC_p& x, const EC_p& a, const ZZ& b);
void mul(EC_p& x, const EC_p& a, long b);
inline void mul(EC_p& x, const ZZ& a, const EC_p& b) {
  mul(x,b,a);
}
inline void mul(EC_p& x, long a, const EC_p& b) {
  mul(x,b,a);
}

inline EC_p operator*(const EC_p& a, const ZZ& b) {
  EC_p x;  mul(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p operator*(const EC_p& a, long b) {
  EC_p x;  mul(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p operator*(const ZZ& a, const EC_p& b) {
  EC_p x;  mul(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p operator*(long a, const EC_p& b) {
  EC_p x;  mul(x,a,b);  NTLX_OPT_RETURN(EC_p,x);
}
inline EC_p& operator*=(EC_p& x, const ZZ& a) {
  mul(x,x,a);  return x;
}
inline EC_p& operator*=(EC_p& x, long a) {
  mul(x,x,a);  return x;
}


/* Class to backup curve.
 */
class EC_pBak {
protected:
  EC_pCurve* backup;
  bool auto_restore;

public:
  EC_pBak();

  // calls restore() if a curve has been backed-up
  ~EC_pBak();

  // save copy of current curve
  void save();

  // restore backup curve to use
  void restore();

private:
  EC_pBak(const EC_pBak&);  // copy disabled
  void operator=(const EC_pBak&);  // assignment disabled
};


NTL_CLOSE_NNS;

#endif
