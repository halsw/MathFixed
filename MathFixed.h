/*
 * This file is part of the MathFixed library
 * Usage: A template library for the implementation
 *        of math functions for use by fixed point types
 * Version 1.0.1
 * Developed by Evan https://github.com/halsw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 * Functions:
 *   fxibits() gets or sets the integer bits for a type
 *   fxfbits() gets or sets the fractional bits for a type
 *   fxnan() gets the representation of NaN
 *   fxisnan() tests if argument is not a number
 *   fxisinf() tests if argument is infinity, but here just a copy of fxisnan()
 *   fxabs() the absolute value
 *   fxmax() the greater of two numbers
 *   fxmin() the lesser of two numbers
 *   fxmaxs() sort two numbers descending and return difference
 *   fxmins() sort two numbers ascending and return difference
 *   fxsq() the square of a number
 *   fxsqrt() the square root
 *   fxhypot() the hypotenuse of two numbers
 *   fxcbrt() the cubic root
 *   fxfloor() the immediately smaller integer of a number
 *   fxceil() the immediately larger integer of a number
 *   fxround() the closest integer of a number
 *   fxtrunc() the integer part of a number
 *   fxdiv() the integer division of two numbers
 *   fxmod() the remainder of the integer division 
 *   fxmodf() get the fractional and integer part of a number
 *   fxfrexp() get the mantissa and exponent(base 2) of a number
 *   fxrandom() get a random number between two limits
 *   fxsin() the sine
 *   fxcos() the cosine
 *   fxtan() the tangent
 *   fxcot() the cotangent
 *   fxatan2() the inverse tangent of the ratio of two numbers
 *   fxatan() the inverse tangent
 *   fxasin() the inverse sine
 *   fxacos() the inverse cosine
 *   fxexp() the natural exponential
 *   fxlog2() the base 2 logarithm
 *   fxlog() the natural logarithm
 *   fxlog10() the base 10 logarithm
 *   fxpow() raise a number to given exponent
 *   fxsinh() the hyperbolic sine
 *   fxsinh() the hyperbolic cosine
 *   fxtanh() the hyperbolic tangent
 */
 
#ifndef FIXEDMATH_H
#define FIXEDMATH_H
#include <limits.h> 

template <class T>
  unsigned char fxibits(unsigned char integer_bits=0) {
    static unsigned char ibits = 0;
    if (integer_bits) ibits=integer_bits;
    return ibits;
  }

template <class T>
  unsigned char fxfbits(unsigned char fractional_bits=0) {
    static unsigned char fbits = 8;
    if (fractional_bits) fbits=fractional_bits;
    return fbits;
  }

template <class T>
  T fxnan() {
     static T r[8]={0.0};
     if (r[0]!=0) return r[0];
     long long* p= (long long*)r;
     unsigned char n = fxibits<T>()+fxfbits<T>();
     *p = -1;
     if ( n & 1 ) 
       *p<<=n;
     return r[0];
  }

template <>
  float fxnan() {
     float r;
     long* p = (long*)&r;
     *p = 0xFFC00000;
     return r;
  }

template <>
  double fxnan() {
     float r;
     long long* p = (long long*)&r;
     *p = 0xFFC0000080000000;
     return r;
  }

template <class T>
  inline bool fxisnan(T x) {
    return x == fxnan<T>();
  }
template <> inline bool fxisnan(double x) {return isnan(x);}    
template <> inline bool fxisnan(float x) {return isnan(x);}    

template <class T>
  inline T fxisinf(T x) {
    return x == fxnan<T>();
  }
template <> inline double fxisinf(double x) {return isinf(x);}    
template <> inline float fxisinf(float x) {return isinf(x);}    
    
template <class T>
  inline T fxabs(T x) {
    return x<0 && x!=fxnan<T>()?-x:x;
  }
template <> inline double fxabs(double x) {return abs(x);}    
template <> inline float fxabs(float x) {return abs(x);}    

template <class T>
  T fxmax(T x, T y) {
    if (x == fxnan<T>()) return x;
    if (y == fxnan<T>()) return y;
    return x<y?y:x;
  }
template <> inline double fxmax(double x, double y) {return max(x,y);}    
template <> inline float fxmax(float x, float y) {return max(x,y);}    

template <class T>
  T fxmin(T x, T y) {
    if (x == fxnan<T>()) return x;
    if (y == fxnan<T>()) return y;
    return x>y?y:x;
  }
template <> inline double fxmin(double x, double y) {return min(x,y);}    
template <> inline float fxmin(float x, float y) {return min(x,y);}    

template <class T>
  T fxmaxs(T *x, T *y) {
    T t=*x;
    if (t == fxnan<T>()) return t;
    if (*y == fxnan<T>()) return *y;
    if (t < *y) {
      *x = *y;
      *y = t;
    }
    return *x - *y;
  }

template <class T>
  T fxmins(T *x, T *y) {
    T t=*x;
    if (t == fxnan<T>()) return t;
    if (*y == fxnan<T>()) return *y;
    if (t > *y) {
      *x = *y;
      *y = t;
    }
    return *y - *x;
  }

template <class T>
  inline T fxsq(T x) {
    return x == fxnan<T>()?x:x*x;
  }
template <> inline double fxsq(double x) {return sq(x);}    
template <> inline float fxsq(float x) {return sq(x);}    

template <class T>
  T fxsqrt(T x) {
  long n = fxfbits<T>();
  if (x == fxnan<T>()) return x;
  T r = 0.5 * x;
  T p = 0.0;
  if ( x < 0) return fxnan<T>();
  do {
    p = r;
    r = 0.5 * ( r + x/r );
    p -= r;
    if ( p<1.0 ) n--;
  } while (p != 0.0 && n);
  return r;
}    
template <> inline double fxsqrt(double x) {return sqrt(x);}    
template <> inline float fxsqrt(float x) {return sqrt(x);}    

template <class T>
  inline T fxhypot(T x, T y) {
  if (x == fxnan<T>()) return x;
  if (y == fxnan<T>()) return y;
  return fxsqrt(x*x + y*y);
}    
template <> inline double fxhypot(double x, double y) {return hypot(x,y);}    
template <> inline float fxhypot(float x, float y) {return hypot(x,y);}    

template <class T>
  T fxcbrt (T x) {
    T r=x / 3.0, e = 1.0;  
    long n = fxfbits<T>();
    if (x == fxnan<T>()) return x;
    if ( x < 0) return fxnan<T>();
    while (e != 0.0 && n) {
      e = ((x / r) / r - r) / 3.0;
      r += e;
      if (fxabs(e)<1.0) n--;
    };
    return(r);
};
template <> inline double fxcbrt(double x) {return cbrt(x);}    
template <> inline float fxcbrt(float x) {return cbrt(x);}    


template <class T>
  T fxfloor(T x) {
     T r[8];
     long long* p = (long long*)r;
     long long mask = -1;
     if (x == fxnan<T>()) return x;
     mask = -1;
     mask<<=fxfbits<T>();
     r[0] = x;
     *p&=mask;     
     return r[0];
  }  
template <> inline double fxfloor(double x) {return floor(x);}    
template <> inline float fxfloor(float x) {return floor(x);}    

template <class T>
  T fxceil(T x) {
    T r = fxfloor(x);
    if (x == r) return(r);
    return r + static_cast<T>(1);
  }
template <> inline double fxceil(double x) {return ceil(x);}    
template <> inline float fxceil(float x) {return ceil(x);}    

template <class T>
  T fxround(T x) {
    T f = fxfloor(x);
    T r;
    if (x == f) return(f);
    r = f+static_cast<T>(1);
    if ( (x-f) < (r-x) ) return f;
    return r;  
  }
template <> inline double fxround(double x) {return round(x);}    
template <> inline float fxround(float x) {return round(x);}    

template <class T>
  T fxtrunc(T x) {
    return x<0?fxceil(x):fxfloor(x);
  }
template <> inline double fxtrunc(double x) {return trunc(x);}    
template <> inline float fxtrunc(float x) {return trunc(x);}    

template <class T>
  T fxdiv(T x, T y) {
    if (x == fxnan<T>()) return x;
    if (y == fxnan<T>() || y == 0) return fxnan<T>();
    return fxtrunc(x/y);
  }

template <class T>
  inline T fxmod(T x, T y) {
    T r = fxdiv(x,y);
    if (r == fxnan<T>()) return r;
    return x - r*y;  
  }
template <> inline double fxmod(double x, double y) {return fmod(x,y);}    
template <> inline float fxmod(float x, float y) {return fmod(x, y);}    

template <class T>
  T fxmodf(T x, T* iptr) {
    if (x == fxnan<T>()) return x;
    *iptr = fxtrunc(x);
    return x - *iptr;  
  }
template <> inline double fxmodf(double x, double* iptr) {return modf(x,iptr);}    

template <class T>
  T fxfrexp(T x, int* expn) {
    int n = 0;
    if (x==0.0 || x == fxnan<T>()) {
      *expn = 0;
      return x;
    }
    if (x<0) return -fxfrexp(-x, expn);
    if (x>=0.5 && x<1.0) {
      *expn = 0;
      return x;
    }
    if (x<0.5)
      while (x<0.5) {
        x*=2;
        n--;
      }
    else
      while (x>1.0) {
        x/=2;
        n++;
      }
    *expn = n;
    return x;
  }
template <> inline double fxfrexp(double x, int* expn) {return frexp(x,expn);}    
template <> inline float fxfrexp(float x, int* expn) {return frexp(x,expn);}    

template <class T>
  T fxrandom(T x, T y = 0.0) {
    T r[8];
    int *p;
    r[0] = fxmins(&x, &y);
    if (r[0] == fxnan<T>()) return r[0];
    if (r[0] == 0.0) return x;
    p = (int*)r;
    for (int i = 0; i<4; i++)
      if (p[i]) {
        p[i] = random(p[i]);
        break;
      } else
        p[i] = random(65536);
    return x+r[0];
  }
template <>
  double fxrandom(double x, double y = 0.0) {
    double r = fxmins(&x, &y);
    r *= random(0x7FFFFFFF)/2147483647.0;
    return x+r;
  }
template <>
  float fxrandom(float x, float y = 0.0) {
    float r = fxmins(&x, &y);
    r *= random(0x7FFFFF)/8388607.0;
    return x+r;
  }

template <class T>
  T fxsin (T x) {
    const T A  = 1.27323954;
    const T B  = 0.405284735;
    const T C  = 0.255;
    T r;
    if (x == fxnan<T>()) return x;
    if (x>0)
      while (x>PI) x-=TWO_PI;
    else  
      while (x<-PI) x+=TWO_PI;
    if ( x<0 )
      r=x;
    else
      r=-x;
    r = x*(A - B*r);
    if ( r<0 )
      r *=-C*(r + 1.0) + 1.0;
    else
      r *= C*(r - 1.0) + 1.0;
    return r;
  }
template <> inline double fxsin(double x) {return sin(x);}    
template <> inline float fxsin(float x) {return sin(x);}    

template <class T>
  inline T fxcos (T x) {
    return x == fxnan<T>()?x:fxsin(x+HALF_PI);
  }
template <> inline double fxcos(double x) {return cos(x);}    
template <> inline float fxcos(float x) {return cos(x);}    

template <class T>
  inline T fxtan (T x) {
    return x == fxnan<T>() || fxmod(x+HALF_PI,static_cast<T>(PI))==0.0?fxnan<T>():fxsin(x)/fxsin(x+HALF_PI);
  }
template <> inline double fxtan(double x) {return tan(x);}    
template <> inline float fxtan(float x) {return tan(x);}    

template <class T>
  inline T fxcot (T x) {
    return x == fxnan<T>() || fxmod(x,static_cast<T>(PI))==0.0?fxnan<T>():fxsin(x+HALF_PI)/fxsin(x);
  }
template <> inline double fxcot(double x) {return 1.0/tan(x);}    
template <> inline float fxcot(float x) {return 1.0/tan(x);}    

template <class T>
  T fxatan2 (T y, T x) {
    T t0, t1, t2, t3, t4, t5;
    if (x == fxnan<T>()) return x;
    if (y == fxnan<T>()) return y;
    t5 = fxabs(x);
    t2 = fxabs(y);
    t0 = fxmax(t5, t2);
    t1 = fxmin(t5, t2);
    t3 = 1.0 / t0;
    t3 *= t1;
    t4 = t3 * t3;
    t0 = -0.013480470;
    t0 *= t4 ;
    t0 += 0.057477314;
    t0 *= t4;
    t0 -= 0.121239071;
    t0 *= t4;
    t0 += 0.195635925;
    t0 *= t4;
    t0 -= 0.332994597;
    t0 *= t4;
    t0 += 0.999995630;
    t3 *= t0;
    if (t2 > t5) t3 = HALF_PI - t3;
    if (x < 0) t3 = PI - t3;
    if (y < 0) t3 = - t3;
    return t3;
  }  
template <> inline double fxatan2(double y, double x) {return atan2(y,x);}    
template <> inline float fxatan2(float y, float x) {return atan2(y,x);}    

template <class T>
  inline T fxatan(T x) {
    return fxatan2(x,static_cast<T>(1.0));
  }
template <> inline double fxatan(double x) {return atan(x);}    
template <> inline float fxatan(float x) {return atan(x);}    

template <class T>
  T fxasin(T x) {
    T r = -0.0187293;
    bool n=false;
    if (x == fxnan<T>()) return x;
    if (x<0) {
      x = -x;
      n = true;
    }  
    r *= x;
    r += 0.0742610;
    r *= x;
    r -= 0.2121144;
    r *= x;
    r += 1.5707288;
    r = HALF_PI - fxsqrt(1.0 - x)*r;
    if (n) return(-r);
    return r;
  }
template <> inline double fxasin(double x) {return asin(x);}    
template <> inline float fxasin(float x) {return asin(x);}    

template <class T>
  T fxacos(T x) {
    T r = -0.0187293;
    bool n=false;
    if (x == fxnan<T>()) return x;
    if (x<0) {
      x = -x;
      n = true;
    }  
    r *= x;
    r += 0.0742610;
    r *= x;
    r -= 0.2121144;
    r *= x;
    r += 1.5707288;
    r *= fxsqrt(1.0-x);
    if (n) return PI - r;
    return r;
  }
template <> inline double fxacos(double x) {return acos(x);}    
template <> inline float fxacos(float x) {return acos(x);}    
  
template <class T>
  T fxexp (T x) {
    static T f[6] = { 0.0,\
      7.389056098930650227230427460575,\
      54.598150033144239078110261202861,\
      2980.9579870417282747435920994529,\
      8886110.5205078726367630237407815,\
      78962960182680.695160978022635108};
    int n;    
    T r = 1.0, t = x;
    if (x == fxnan<T>()) return x;
    if (x == 0.0) return 1.0;
    if (x<0.0) return 1.0/fxexp(-x);
    if (x>=1.0) {
      if (f[0] == 0.0) {
        f[0] = 2.7182818284590452353602874713527;
        for (n=1; n<6; n++) if (f[n]<=f[n-1]) break; 
        for (n=n; n<6; n++) f[n] = fxnan<T>();       
      }
      n = 0;
      while (t>=1.0) {
        if (fxmod(t,(T)2.0)>=1.0)
          if (f[n] == fxnan<T>())
            return fxnan<T>();
          else  
            r*=f[n];
        t /= 2;
        n++;
      }
      t = fxmod(x,(T)1.0);
    }  
    n = 0;  
    while ((x=t/2) != 0.0) {
      t=x;
      n++;
    }
    t += 1.0;
    while (n-- > 0) t*=t; 
    return r*t;
  }
template <> inline double fxexp(double x) {return exp(x);}    
template <> inline float fxexp(float x) {return exp(x);}

template <class T>
  T fxlog2(T x) {
    int expn;
    if (x == fxnan<T>() || x<=0.0) return fxnan<T>();
    x = fxfrexp(x, &expn);
    return expn + x*(2.8854 - 0.7213*x) - 2.1640;
  }

template <class T>
  inline T fxlog(T x) {
    T r=fxlog2(x);
    return r == fxnan<T>() ? r:0.69314718056*r;
  }
template <> inline double fxlog(double x) {return log(x);}    
template <> inline float fxlog(float x) {return log(x);}

template <class T>
  inline T fxlog10(T x) {
    T r=fxlog2(x);
    return r == fxnan<T>() ? r:0.30103*r;
  }
template <> inline double fxlog10(double x) {return log10(x);}    
template <> inline float fxlog10(float x) {return log10(x);}

template <class T>
  inline T fxpow (T x, T y) {
    T t = 1.0;
    T b = x;
    T p = 0.0;
    bool ovf = fxabs(x)>1.0;
    if ( x == fxnan<T>() ) return x;
    if ( y == fxnan<T>() ) return y;
    if ( x == 0.0 ) return 0.0;
    if ( y == 0.0 || x == 1.0 ) return 1.0;
    if (y<0) { 
      t = fxpow(x,-y);
      return t==0.0 ? fxnan<T>() : 1.0/t;
    }  
    while ( y >= 1.0 && t != 0.0 ) {
      if (ovf && fxabs(b)<p) return fxnan<T>(); 
      if ( fxmod(y,static_cast<T>(2)) >= static_cast<T>(1) ) t*=b;
      if (ovf) p = fxabs(b);
      b *= b;
      y /= 2.0;
    }
    if ( y == 0.0 || t == 0.0 ) return t;
    if ( x < 0.0 ) return fxnan<T>();
    return t*fxexp(y*fxlog(x));
  }
template <> inline double fxpow(double x, double y) {return pow(x,y);}    
template <> inline float fxpow(float x, float y) {return pow(x,y);}

template <class T>
  inline T fxsinh(T x) {
    return ( x == fxnan<T>() ) ? fxnan<T>() : (fxexp(x) - fxexp(-x)) / 2.0;
  }
template <> inline double fxsinh(double x) {return sinh(x);}    
template <> inline float fxsinh(float x) {return sinh(x);}

template <class T>
  inline T fxcosh(T x) {
    return ( x == fxnan<T>() ) ? fxnan<T>() : (fxexp(x) + fxexp(-x)) / 2.0;
  }
template <> inline double fxcosh(double x) {return cosh(x);}    
template <> inline float fxcosh(float x) {return cosh(x);}

template <class T>
  inline T fxtanh(T x) {
    if ( x == fxnan<T>() ) return x;
    x=fxexp(2*x);
    return (x - 1.0)/(x + 1.0);
  }
template <> inline double fxtanh(double x) {return tanh(x);}    
template <> inline float fxtanh(float x) {return cosh(x);}

#endif
