#include "EC_p.h"


/* Implementation of points on elliptic curve.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */


/**************** EC_pCurve methods ****************/

const long EC_pCurve::AFFINE=0;
const long EC_pCurve::PROJECTIVE=1;
const long EC_pCurve::JACOBIAN=2;
const long EC_pCurve::MONTGOMERY=3;

// assignment to polynomial
void conv(EC_pCurve& curve, const ZZ_pX& f) {
  if (deg(f)!=3 || !IsOne(coeff(f,3)))
    Error("conv: f must be monic and degree 3");
  curve.rep=EC_pCurve::PROJECTIVE;
  curve.f=f;
}

// test if curve is singular
bool IsSingular(const EC_pCurve& c) {
  // test if 4A^3 + 27B^2 - 18ABC - (AC)^2 + 4BC^3 == 0
  const ZZ_p& A = coeff(c.y2(),1);
  const ZZ_p& B = coeff(c.y2(),0);
  const ZZ_p& C = coeff(c.y2(),2);
  ZZ_p A2,C2,BC;
  sqr(A2,A);
  sqr(C2,C);
  mul(BC,B,C);
  return IsZero(27*sqr(B) + A*(4*A2-18*BC) + C2*(4*BC-A2));
}

// set to equivalent curve, but with a_2=0
void MakeStandard(EC_pCurve& curve) {
  const ZZ_p& C = coeff(curve.y2(),2);
  if (!IsZero(C)) {
    const ZZ_p& A = coeff(curve.y2(),1);
    const ZZ_p& B = coeff(curve.y2(),0);
    ZZ_p D(-C/3);
    SetCoeff(curve,0,B+D*(D*(D+C)+A));
    SetCoeff(curve,1,A+D*(3*D+2*C));
    SetCoeff(curve,2,ZZ_p::zero());
  }
}

ostream& operator<<(ostream& stream, const EC_pCurve& curve) {
  return stream<<"["<<curve.y2()<<" "<<curve.representation()<<"]";
}


/**************** EC_p methods ****************/

EC_pCurve* EC_p::c = NULL;
EC_p* EC_p::inf = NULL;

EC_p::EC_p() : Y(ZZ_p_NoAlloc),Z(ZZ_p_NoAlloc) {
  switch (c->rep) {
  case EC_pCurve::AFFINE:
    Y._ZZ_p__rep.SetSize(ZZ_p::ModulusSize());
    set(Y);
    break;
  case EC_pCurve::PROJECTIVE:
  case EC_pCurve::JACOBIAN:
    Y._ZZ_p__rep.SetSize(ZZ_p::ModulusSize());
    set(Y);
    Z._ZZ_p__rep.SetSize(ZZ_p::ModulusSize());
    break;
  case EC_pCurve::MONTGOMERY:
    Z._ZZ_p__rep.SetSize(ZZ_p::ModulusSize());
    break;
  default:
    Error("EC_p: unknown point representation");
  }
}

// test if this is a point on the curve
bool EC_p::isValid() const {
  if (IsZero(Z))
    return true;
  const ZZ_p& A = coeff(c->f,1);
  const ZZ_p& B = coeff(c->f,0);
  const ZZ_p& C = coeff(c->f,2);
  switch (c->rep) {
  case EC_pCurve::AFFINE:
    // Y^2 = X^3 + CX^2 + AX + B
    return IsZero(eval(c->f,X)-sqr(Y));
  case EC_pCurve::PROJECTIVE: {
    // (Y^2)Z = X^3 + C(X^2)Z + AX(Z^2) + B(Z^3)
    ZZ_p Z2;
    sqr(Z2,Z);
    return IsZero(X*(X*(X+C*Z)+A*Z2)+(B*Z2-sqr(Y))*Z);
  }
  case EC_pCurve::MONTGOMERY:
    // Y^2 = X^3/Z + C(X^2) + AXZ + B(Z^2) and test for quadratic residue
    return Jacobi(rep(X*(X*(X*inv(Z)+C)+A*Z)+B*sqr(Z)),ZZ_p::modulus())!=-1;
  case EC_pCurve::JACOBIAN: {
    // Y^2 = X^3 + C(X^2)(Z^2) + AX(Z^4) + B(Z^6)
    ZZ_p Z2,Z4;
    sqr(Z2,Z);
    sqr(Z4,Z2);
    return IsZero(X*(X*(X+C*Z2)+A*Z4)+B*Z2*Z4-sqr(Y));
  }
  }
  Error("EC_p: unknown point representation");
  return false;
}

// compute affine coordinates
void EC_p::affine(ZZ_p& x, ZZ_p& y) const {
  if (IsZero(Z))
    Error("EC_p: affine representation of point at infinity");
  switch (c->rep) {
  case EC_pCurve::AFFINE:
    x=X; y=Y;
    break;
  case EC_pCurve::PROJECTIVE:
    if (IsOne(Z)) {
      x=X; y=Y;
    }
    else {
      ZZ_p Zinv;
      inv(Zinv,Z);
      mul(x,X,Zinv);
      mul(y,Y,Zinv);
    }
    break;
  case EC_pCurve::MONTGOMERY:
    // should compute some y here
    Error("EC_p: affine() noY not implemented");
    break;
  case EC_pCurve::JACOBIAN:
    if (IsOne(Z)) {
      x=X; y=Y;
    }
    else {
      ZZ_p Zinv,Zinv2;
      inv(Zinv,Z);
      sqr(Zinv2,Zinv);
      mul(x,X,Zinv2);
      Zinv*=Zinv2;
      mul(y,Y,Zinv);
    }
    break;
  default:
    Error("EC_p: unknown point representation");
  }
}

// compute affine coordinate (x only)
void EC_p::affine(ZZ_p &x) const {
  if (IsZero(Z))
    Error("EC_p: affine representation of point at infinity");
  switch (c->rep) {
  case EC_pCurve::AFFINE:
    x=X;
    break;
  case EC_pCurve::PROJECTIVE:
  case EC_pCurve::MONTGOMERY:
    if (IsOne(Z))
      x=X;
    else {
      ZZ_p Zinv;
      inv(Zinv,Z);
      mul(x,X,Zinv);
    }
    break;
  case EC_pCurve::JACOBIAN:
    if (IsOne(Z))
      x=X;
    else {
      ZZ_p Zinv,Zinv2;
      inv(Zinv,Z);
      sqr(Zinv2,Zinv);
      mul(x,X,Zinv2);
    }
    break;
  default:
    Error("EC_p: unknown point representation");
  }
}

// adjust so Z=1
void EC_p::makeAffine() {
  if (IsOne(Z)||IsZero(Z))
    return;
  // need inverse of Z here
  ZZ_p Zinv;
  inv(Zinv,Z);
  switch (c->rep) {
  case EC_pCurve::AFFINE:
    Error("EC_p: affine point with Z!=1");
    break;
  case EC_pCurve::PROJECTIVE:
    Y*=Zinv;
  case EC_pCurve::MONTGOMERY:
    X*=Zinv;
    break;
  case EC_pCurve::JACOBIAN: {
    ZZ_p Zinv2;
    sqr(Zinv2,Zinv);
    X*=Zinv2;
    Y*=Zinv2;
    Y*=Zinv;
  }
    break;
  default:
    Error("EC_p: unknown point representation");
  }
  set(Z);
}

// swap points
void swap(EC_p& a, EC_p& b) {
  swap(a.X,b.X);
  swap(a.Y,b.Y);
  swap(a.Z,b.Z);
}

// set to point at infinity
void clear(EC_p& x) {
  clear(x.X); clear(x.Z);
  set(x.Y);
}

// test for equality
bool operator==(const EC_p& a, const EC_p& b) {
  if (IsZero(a.Z))
    return IsZero(b.Z);
  if (IsZero(b.Z))
    return false;
  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
    return a.X==b.X && a.Y==b.Y;
  case EC_pCurve::PROJECTIVE:
    if (a.Z==b.Z) return a.X==b.X && a.Y==b.Y;
    return a.X*b.Z==b.X*a.Z && a.Y*b.Z==b.Y*a.Z;
  case EC_pCurve::MONTGOMERY:
    if (a.Z==b.Z) return a.X==b.X;
    return a.X*b.Z==b.X*a.Z;
  case EC_pCurve::JACOBIAN: {
    if (a.Z==b.Z) return a.X==b.X && a.Y==b.Y;
    ZZ_p aZ2,bZ2;
    sqr(aZ2,a.Z);
    sqr(bZ2,b.Z);
    if (a.X*bZ2!=b.X*aZ2)
      return false;
    aZ2*=a.Z;  // aZ2 = a.Z^3
    bZ2*=b.Z;  // bZ2 = b.Z^3
    return a.Y*bZ2==b.Y*aZ2;
  }
  }
  Error("EC_p: unknown point representation");
  return false;
}

// stream out
ostream& operator<<(ostream& stream, const EC_p& x) {
  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
    return stream<<"("<<x.X<<" "<<x.Y<<")";
  case EC_pCurve::PROJECTIVE:
    return stream<<"["<<x.X<<" "<<x.Y<<" "<<x.Z<<"]";
  case EC_pCurve::MONTGOMERY:
    return stream<<"["<<x.X<<" : "<<x.Z<<"]";
  case EC_pCurve::JACOBIAN:
    return stream<<"<"<<x.X<<" "<<x.Y<<" "<<x.Z<<">";
  }
  Error("EC_p: unknown point representation");
  return stream;
}


/**************** EC_pBak methods ****************/

EC_pBak::EC_pBak() {
  backup=NULL;
  auto_restore=false;
}

EC_pBak::~EC_pBak() {
  if (auto_restore)
    restore();
  if (backup) delete backup;
}

void EC_pBak::save() {
  if (backup) delete backup;
  backup = EC_p::c ? new EC_pCurve(*EC_p::c) : NULL;
  auto_restore = true;
}

void EC_pBak::restore() {
  if (backup)
    EC_p::init(*backup);
  else {
    if (EC_p::c) {
      delete EC_p::c;
      EC_p::c = NULL;
    }
  }
  auto_restore=false;
}

EC_pBak::EC_pBak(const EC_pBak&) {
  Error("EC_pBak: copy not allowed");
} 
void EC_pBak::operator=(const EC_pBak&) {
  Error("EC_pBak: assignment not allowed");
}


/**************** random() ****************/


// random point
void random(EC_p& p) {
  set(p.Z);
  ZZ_p Y2;
  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
  case EC_pCurve::PROJECTIVE:
  case EC_pCurve::JACOBIAN:
    do {
      random(p.X);
      // y^2 = x*(x*(x+C)+A) + B
      eval(Y2,EC_p::curve().y2(),p.X);
      // check if Y2 is a square
      long j = Jacobi(rep(Y2),ZZ_p::modulus());
      if (j!=-1) {
	if (j==0)
	  clear(p.Y);
	else {
	  ZZ y;
	  SqrRootMod(y,rep(Y2),ZZ_p::modulus());
	  conv(p.Y,y);
	}
	return;
      }
    } while (true);

  case EC_pCurve::MONTGOMERY:
    do {
      random(p.X);
      // y^2 = x*(x*(x+C)+A)+B
      eval(Y2,EC_p::curve().y2(),p.X);
      // check if Y2 is a square
      if (Jacobi(rep(Y2),ZZ_p::modulus())!=-1)
	return;
    } while (true);

  default:
    Error("EC_p: unknown point representation");
  }
}


/**************** point addition ****************/


// negate point
void negate(EC_p& x, const EC_p& a) {
  if (&x!=&a) {
    x.X = a.X;
    x.Z = a.Z;
  }
  if (EC_p::curve().rep!=EC_pCurve::MONTGOMERY)
    negate(x.Y,a.Y);
}

// double point
void doub(EC_p& x, const EC_p& a) {
  if (IsZero(a.Z)) {
    clear(x);
    return;
  }

  const ZZ_p& A = coeff(EC_p::curve().y2(),1);
  const ZZ_p& B = coeff(EC_p::curve().y2(),0);
  const ZZ_p& C = coeff(EC_p::curve().y2(),2);

  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
    if (IsZero(a.Y)) {
      clear(x);
      return;
    }
    if (!IsZero(C))
      Error("EC_p: double() C!=0 not supported");
    Error("EC_p: double() affine coordinates not supported");
    break;

  case EC_pCurve::PROJECTIVE: {
    if (IsZero(a.Y)) {
      clear(x);
      return;
    }
    if (!IsZero(C))
      Error("EC_p: double() C!=0 not supported");
    ZZ_p l,m,v,c2;
    conv(c2,2);
    // mu = 3X^2 + AZ^2
    conv(m,3);
    sqr(l,a.X);
    mul(m,m,l);
    sqr(l,a.Z);
    mul(l,l,A);
    add(m,m,l);
    // lambda = 2XY
    // nu = 2YZ
    mul(l,c2,a.Y);
    mul(v,l,a.Z);
    mul(l,l,a.X);
    // Y = mu*(3*lambda*nu - mu^2) - 2*Y^2*nu^2
    sqr(x.Y,a.Y);
    mul(l,l,v);
    ZZ_p mm,vv,tmp;
    sqr(mm,m);
    sqr(vv,v);
    mul(x.Y,c2,x.Y);
    mul(x.Y,x.Y,vv);
    conv(tmp,3);
    mul(tmp,tmp,l);
    sub(tmp,tmp,mm);
    mul(tmp,tmp,m);
    sub(x.Y,tmp,x.Y);
    // X = nu*(mu^2 - 2*lambda*nu)
    mul(x.X,c2,l);
    sub(x.X,mm,x.X);
    mul(x.X,x.X,v);
    // Z = nu^3
    mul(x.Z,v,vv);
    break;
  }

  case EC_pCurve::JACOBIAN: {
    if (IsZero(a.Y)) {
      clear(x);
      return;
    }
    if (!IsZero(C))
      Error("EC_p: double() C!=0 not supported");
    ZZ_p M,S,tmp;
    // M = 3X^2 + AZ^4
    sqr(M,a.X);
    conv(tmp,3);
    mul(M,M,tmp);
    sqr(tmp,a.Z);
    sqr(tmp,tmp);
    mul(tmp,tmp,A);
    add(M,M,tmp);
    // S = 4XY^2
    sqr(S,a.Y);
    mul(S,S,a.X);
    conv(tmp,4);
    mul(S,S,tmp);
    // new Z = 2YZ
    conv(tmp,2);
    mul(x.Z,tmp,a.Z);
    mul(x.Z,x.Z,a.Y);
    // new X = M^2 - 2S
    ZZ_p tmp2;
    sqr(tmp2,M);
    conv(tmp,2);
    mul(tmp,tmp,S);
    sub(x.X,tmp2,tmp);
    // new Y = M(S - new X) - 8Y^4
    sqr(tmp2,a.Y);
    sqr(tmp2,tmp2);
    conv(tmp,8);
    mul(tmp2,tmp,tmp2);
    sub(tmp,S,x.X);
    mul(tmp,tmp,M);
    sub(x.Y,tmp,tmp2);
    break;
  }

  case EC_pCurve::MONTGOMERY: {
    if (IsZero(B)) {
      ZZ_p X2;
      sqr(X2,a.X);
      ZZ_p newX,newZ;
      //newZ = 4*a.Z*a.X*(X2+(C*a.X+A*a.Z)*a.Z);
      mul(newZ,C,a.X);
      mul(newX,A,a.Z);
      newZ+=newX;
      newZ*=a.Z;
      newZ+=X2;
      newZ*=a.X;
      newZ*=a.Z;
      newZ*=4;
      //newX = sqr(X2-A*sqr(a.Z));
      sqr(newX,a.Z);
      newX*=A;
      sub(newX,X2,newX);
      sqr(newX,newX);
      // swap is faster than assignment
      swap(x.X,newX);
      swap(x.Z,newZ);
    }
    else {
      // general solution
      ZZ_p X2;
      sqr(X2,a.X);
      ZZ_p newX,newZ;
      newX = sqr(X2-A*sqr(a.Z)) - 4*B*(2*a.X+C*a.Z)*a.Z*sqr(a.Z);
      newZ = 4*a.Z*(a.X*X2 + C*X2*a.Z + A*a.X*sqr(a.Z) + B*a.Z*sqr(a.Z));
      swap(x.X,newX);
      swap(x.Z,newZ);
    }
    break;
  }

  default:
    Error("EC_p: unknown point representation");
  }
}

// add two points
void add(EC_p& x, const EC_p& a, const EC_p& b) {
  if (IsZero(a.Z)) {
    if (&x!=&b)
      x=b;
    return;
  }
  if (IsZero(b.Z)) {
    if (&x!=&a)
      x=a;
    return;
  }
  if (&a==&b) {
    doub(x,a);
    return;
  }

  const ZZ_p& C = coeff(EC_p::curve().y2(),2);
  if (!IsZero(C))
    Error("EC_p: add() C!=0 not supported");

  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
    Error("EC_p: add() affine coordinates not supported");
    break;

  case EC_pCurve::PROJECTIVE: {
    ZZ_p aU,bU;
    mul(aU,b.X,a.Z);
    mul(bU,a.X,b.Z);
    ZZ_p aS,bS;
    mul(aS,b.Y,a.Z);
    mul(bS,a.Y,b.Z);
    ZZ_p W,R;
    sub(W,aU,bU);
    sub(R,aS,bS);
    if (IsZero(W)) {
      if (IsZero(R))
	doub(x,a);
      else
	clear(x);
      return;
    }
    ZZ_p T,M;
    add(T,aU,bU);
    add(M,aS,bS);
    // a.Z*b.Z
    ZZ_p ZZ;
    mul(ZZ,a.Z,b.Z);
    // R^2*a.Z*b.Z
    ZZ_p RRZZ;
    sqr(RRZZ,R);
    RRZZ*=ZZ;
    // W^2*T and W^3
    ZZ_p WWW,WWT;
    sqr(WWT,W);
    mul(WWW,WWT,W);
    WWT*=T;
    // new X = W * (R^2*a.Z*b.Z - W^2*T)
    sub(x.X,RRZZ,WWT);
    x.X*=W;
    // new Z = W^3 * a.Z * b.Z
    mul(x.Z,WWW,ZZ);
    // new Y = (R*(3*W^2*T - R^2*a.Z*b.Z) - W^3*M) / 2
    ZZ_p tmp;
    conv(tmp,3);
    mul(x.Y,tmp,WWT);
    conv(tmp,2);
    RRZZ*=tmp;
    x.Y-=RRZZ;
    x.Y*=R;
    WWW*=M;
    x.Y-=WWW;
    x.Y/=tmp;
    break;
  }

  case EC_pCurve::JACOBIAN: {
    ZZ_p aZ,bZ;
    sqr(aZ,a.Z);
    sqr(bZ,b.Z);
    ZZ_p aU,bU;
    mul(aU,b.X,aZ);
    mul(bU,a.X,bZ);
    ZZ_p aS,bS;
    mul(aZ,aZ,a.Z);  // aZ = a.Z^3
    mul(bZ,bZ,b.Z);  // bZ = b.Z^3
    mul(aS,b.Y,aZ);
    mul(bS,a.Y,bZ);
    ZZ_p W,R;
    sub(W,aU,bU);
    sub(R,aS,bS);
    if (IsZero(W)) {
      if (IsZero(R))
	doub(x,a);
      else
	clear(x);
      return;
    }
    ZZ_p T,M;
    add(T,aU,bU);
    add(M,aS,bS);
    // set new Z (x.Z = a.Z * b.Z * W)
    mul(x.Z,a.Z,b.Z);
    mul(x.Z,x.Z,W);
    // some temps
    mul(M,M,W);
    sqr(W,W);
    mul(M,M,W);  // M = MW^3
    mul(T,T,W);  // T = TW^2
    // set new X (x.X = R^2 - TW^2)
    ZZ_p tmp;
    sqr(tmp,R);
    sub(x.X,tmp,T);
    // set new Y (x.Y = ((TW^2 - 2x.X)R - MW^3) / 2)
    conv(tmp,2);
    mul(tmp,tmp,x.X);
    sub(T,T,tmp);
    mul(T,T,R);
    sub(T,T,M);
    conv(tmp,2);
    div(x.Y,T,tmp);
    break;
  }

  case EC_pCurve::MONTGOMERY:
    Error("EC_p: use addh() for montgomery coordinates");
    break;

  default:
    Error("EC_p: unknown point representation");
  }
}

// add for MONTGOMERY representation (c=a-b)
void addh(EC_p& x, const EC_p& a, const EC_p& b, const EC_p& c) {
  switch (EC_p::curve().rep) {
  case EC_pCurve::AFFINE:
  case EC_pCurve::PROJECTIVE:
  case EC_pCurve::JACOBIAN:
    add(x,a,b);
    break;

  case EC_pCurve::MONTGOMERY: {
    const ZZ_p& A = coeff(EC_p::curve().y2(),1);
    const ZZ_p& B = coeff(EC_p::curve().y2(),0);
    const ZZ_p& C = coeff(EC_p::curve().y2(),2);
    if (IsZero(B)) {
      ZZ_p newX,newZ,t;
      //newX = c.Z * sqr(a.X*b.X-A*a.Z*b.Z))
      mul(newX,a.X,b.X);
      mul(t,a.Z,b.Z);
      t*=A;
      newX-=t;
      sqr(newX,newX);
      newX*=c.Z;
      //newZ = c.X * sqr(a.X*b.Z-b.X*a.Z))
      mul(newZ,a.X,b.Z);
      mul(t,b.X,a.Z);
      newZ-=t;
      sqr(newZ,newZ);
      newZ*=c.X;
      // swap is faster than assignment
      swap(x.X,newX);
      swap(x.Z,newZ);
    }
    else {
      // general solution
      ZZ_p ZZ;
      mul(ZZ,a.Z,b.Z);
      ZZ_p newX,newZ;
      newX = c.Z*(sqr(a.X*b.X-A*ZZ)-4*B*(a.X*b.Z+b.X*a.Z+C*ZZ)*ZZ);
      newZ = c.X*sqr(a.X*b.Z-b.X*a.Z);
      swap(x.X,newX);
      swap(x.Z,newZ);
    }
    break;
  }

  default:
    Error("EC_p: unknown point representation");
  }
}

// subtract two points
void sub(EC_p& x, const EC_p& a, const EC_p& b) {
  EC_p nb;
  negate(nb,b);
  add(x,a,nb);
}


/**************** multiplication ****************/


// add a point to itself several times
void mul(EC_p& x, const EC_p& a, const ZZ& n) {
  if (n<=1) {
    if (n<0)
      mul(x,-a,-n);
    else if (IsZero(n))
      clear(x);
    else if (&x!=&a) // n==1
      x=a;
    return;
  }

  if (EC_p::curve().rep==EC_pCurve::MONTGOMERY) {
    if (n==2) {
      doub(x,a);
      return;
    }
    EC_p u(a);
    EC_p t;
    doub(t,a);
    for (long j=NumBits(n)-2; j>0; --j) {
      if (bit(n,j)!=0) {
	addh(u,t,u,a);
	doub(t,t);
      }
      else {
	addh(t,u,t,a);
	doub(u,u);
      }
    }
    if (IsOdd(n))
      addh(x,u,t,a);
    else
      doub(x,u);
    return;
  }

  // default implementation
  ZZ m(3*n);
  EC_p result(a);
  for (long e=NumBits(m)-2; e>0; --e) {
    doub(result,result);
    if (bit(m,e)==0) {
      if (bit(n,e)!=0) 
	sub(result,result,a);
    }
    else if (bit(n,e)==0) 
      add(result,result,a);
  }
  swap(x,result);
}

void mul(EC_p& x, const EC_p& a, long n) {
  if (n<=1) {
    if (n<0)
      mul(x,-a,-n);
    else if (n==0)
      clear(x);
    else if (&x!=&a) // n==1
      x=a;
    return;
  }

  switch (EC_p::curve().rep) {
  case EC_pCurve::MONTGOMERY:
    if (n==2) {
      doub(x,a);
      return;
    }
    EC_p u(a);
    EC_p t;
    doub(t,a);
    for (long j=1<<(NumBits(n)-2); j>1; j>>=1) {
      if (n&j) {
	addh(u,t,u,a);
	doub(t,t);
      }
      else {
	addh(t,u,t,a);
	doub(u,u);
      }
    }
    if (n&1)
      addh(x,u,t,a);
    else
      doub(x,u);
    return;
  }

  // make sure n is not too large (3*n must fit in a long)
  if (n>=(1L<<(NTL_BITS_PER_LONG-3))) {
    mul(x,a,to_ZZ(n));
    return;
  }

  // default implementation
  long m=3*n;
  long mask=1<<(NumBits(m)-2);
  EC_p result(a);
  do {
    doub(result,result);
    if ((m&mask)==0) {
      if ((n&mask)!=0)
	sub(result,result,a);
    }
    else if ((n&mask)==0)
      add(result,result,a);
    mask>>=1;
  } while (mask>1);
  swap(x,result);
}

