/*
 * This file is part of the MathFixed library
 * Usage: A template library for the implementation
 *        of math functions for use by fixed point types
 * Version 1.1.0
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
 *   fxpgmread() read PROGMEM stored values
 *   fxsize() gets number of real numbers needed for representation
 *   fxnan() gets the representation of NaN
 *   fxisnan() tests if argument is not a number
 *   fxisinf() tests if argument is infinity, but here just a copy of fxisnan()
 *   fxabs() the absolute value
 *   fxalmost() set the almost equal comparison tolerance if argument >1 then tolerance = 1/argument
 *   fxequ() almost equal comparison (set by fxalmost tollerance)
 *   fxsign() the sign (-1, 1)
 *   fxcopysign() copies the sign of second argument to the first one
 *   fxmax() the greater of two numbers
 *   fxmin() the lesser of two numbers
 *   fxswap() swap two numbers
 *   fxmaxs() sort two numbers descending and return difference
 *   fxmins() sort two numbers ascending and return difference
 *   fxsq() the square of a number
 *   fxsqrt() the square root
 *   fxhypot() the hypotenuse of two numbers
 *   fxnorm() the hypotenuse of an array of numbers
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
 *   fxcosh() the hyperbolic cosine
 *   fxtanh() the hyperbolic tangent
 *   fxconj() the conjugate of a number (used for compatibility with Complex numbers & Quaternions)
 *   
 * fxm namespace:  
 *   bool fxm::readpgm, variable that controls fxpgmread() access to program memory
 *   constexpr fxm::fxinfo_s info, structure that holds basic info about a used fixed point type
 *   fxm::pgm_read() function to access program memory
 *   fxm::fracbits() constexpr function returning the fractional bits of fixed point types
 *   fxm::intbits() constexpr function returning the integer bits of fixed point types
 *   fxm::zero() constexpr zero reprepresenation of fixed point types
 *   fxm::one() constexpr unit(1) reprepresenation of fixed point types
 *   fxm::accuracy() constexpr smallest in magnitude number of fixed point types 
 *   fxm::convert() constexpr conversion of standard types to fixed point types 
 *   
 * Defines:
 *   fx() macro to store fixed point types in program memory
 *   TFixed the default type used by the library (can also be floating type)
 *   FIXEDMATH_MAX the maximum in bit count integer that contains all used fixed point types
 *   
 */
 
#ifndef FIXEDMATH_H
#define FIXEDMATH_H

#ifndef TFixed
#define TFixed float
#endif

#ifndef FIXEDMATH_MAX
#define FIXEDMATH_MAX int32_t
#endif

#include <limits.h> 


#define fx( x ) fxm::convert<TFixed>( x )

#ifdef PROGMEM

namespace fxm { 
  bool readpgm = true;
  template <class T>  inline T pgm_read(const void* p); 
}

template <class T>
inline T fxpgmread(const T* x) { 
  if (!x)  return (T)(fxm::readpgm = !fxm::readpgm);
  return fxm::readpgm ? fxm::pgm_read<T>(x) : ((T*)x)[0]; 
}

#else

template <class T>
inline T fxpgmread(const void* x) { 
  if (!x) return 0;
  return ((T*)x)[0]; 
}

#define PROGMEM
#define pgm_read_byte( p ) (char*)( p )[0]
#define pgm_read_word( p ) (uint16_t*)( p )[0]
#define pgm_read_dword( p ) (uint32_t*)( p )[0]
#define memcpy_P( d, s, l ) memcpy( d, s, l )
#endif

namespace fxm {
  template <class T> union f {
    FIXEDMATH_MAX bin;
    T val;
  };  

  template <class T> struct fxinfo_s {
    uint8_t intp;
    uint8_t frac;
    f<T> nan;
    T min;
  };  

  template <class T>
    inline T pgm_read(const void* p) {
      T val;
      memcpy_P(&val, p, sizeof(T));
      return val;
    }
  
  template <class T> constexpr T shift(int8_t n) { return n==1 ? 0.5 : (n>0 ? shift<T>(n-1)*(T)0.5 : shift<T>(n+1)/(T)0.5); }

  template <class T, int8_t N=1> constexpr uint8_t fracbits() {return shift<T>(N) == (T)0 ? N-1 : fracbits<T,N+1>();  }
  template <> constexpr uint8_t fracbits<float,1>() {return 0; }
  template <> constexpr uint8_t fracbits<double,1>() {return 0; }

  template <class T, int8_t N=0> constexpr uint8_t intbits() {return shift<T>(N) <= shift<T>(N+1) ? -N : intbits<T,N-1>();  }
  template <> constexpr uint8_t intbits<float,0>() {return 0; }
  template <> constexpr uint8_t intbits<double,0>() {return 0; }

  template <class T> constexpr f<T> nan() {
   return (f<T>){.bin=(intbits<T>()+fracbits<T>()) ? ((int32_t)-1) << ((intbits<T>()+fracbits<T>()) & 1)*(intbits<T>()+fracbits<T>()) : 0xFFC00000};
  }

  template <class T> constexpr f<T> zero() {
   return f<T>{.bin=0};
  }

  template <class T> constexpr f<T> one() {
   return f<T>{.bin=(FIXEDMATH_MAX)1 << fracbits<T>()};
  }

  template <class T> constexpr f<T> accuracy() {
   return f<T>{.bin=1};
  }

  template <class T> constexpr T convert( float x ) {
   return fracbits<T>() ? x : (FIXEDMATH_MAX)round(x *((FIXEDMATH_MAX)1 << fracbits<T>()));
  }

  template <class T> constexpr fxinfo_s<T> info  = {
    .intp = intbits<T>(),
    .frac = fracbits<T>(),
    .nan = nan<T>()
  };

}

template <class T=TFixed> constexpr size_t fxsize(T x) { return 1;}

template <class T=TFixed> inline T fxnan() { return fxm::info<T>.nan.val;}

template <class T=TFixed> inline T& fxnan( T& x ) { return x = fxm::info<T>.nan.val;}

template <class T=TFixed> inline bool fxisnan(T x) {return x == fxm::info<T>.nan.val;}

template <class T=TFixed> inline bool fxisnan(void* x) {return ((T*)x)[0] == fxm::info<T>.nan.val;}

template <> inline bool fxisnan(double x) {return isnan(x);}    
template <> inline bool fxisnan(float x) {return isnan(x);}    

template <class T=TFixed> inline bool fxisinf(T x) {return x == fxm::info<T>.nan.val;}

template <> inline bool fxisinf(double x) {return isinf(x);}    
template <> inline bool fxisinf(float x) {return isinf(x);}    

template <class T=TFixed>
  inline T fxabs(T x) {
    return x<0 && x!=fxm::info<T>.nan.val?-x:x;
  }
template <> inline double fxabs(double x) {return abs(x);}    
template <> inline float fxabs(float x) {return abs(x);}    

template <class T=TFixed>
  inline T fxalmost(T e = fxm::info<T>.nan.val ) {
    static T epsilon = 0.0;
    if ( fxisnan(e)  ) return epsilon;
    epsilon = fxabs(e);
    if (epsilon>static_cast<T>(1.0)) epsilon = static_cast<T>(1.0) / epsilon;
  }

template <class T=TFixed>
  inline bool fxequ(T x, T y = static_cast<T>(0.0)) {
    if ( fxisnan(x) || fxisnan(y) ) return false;
    return fxabs(x-y) <= fxalmost<T>();
  }

template <class T=TFixed>
  inline T fxsign(T x) {
    return x<0.0 ? static_cast<T>(-1.0) : static_cast<T>(1.0);
  }

template <class T=TFixed>
  inline T fxcopysign(T x, T y) {
    return x<0.0 ^ y<0.0 ? -x : x;
  }

template <class T=TFixed>
  T fxmax(T x, T y) {
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val) return y;
    return x<y?y:x;
  }
template <> inline double fxmax(double x, double y) {return max(x,y);}    
template <> inline float fxmax(float x, float y) {return max(x,y);}    

template <class T=TFixed>
  T fxmin(T x, T y) {
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val) return y;
    return x>y?y:x;
  }
template <> inline double fxmin(double x, double y) {return min(x,y);}    
template <> inline float fxmin(float x, float y) {return min(x,y);}    

template <class T=TFixed>
  inline void fxswap(T& x, T& y) {
    T t=x;
    x = y;
    y = t;
  }

template <class T=TFixed>
  T fxmaxs(T& x, T& y) {
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val) return y;
    if (x < y) fxswap(x, y);
    return x - y;
  }

template <class T=TFixed>
  T fxmins(T& x, T& y) {
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val) return y;
    if (x > y) fxswap(x, y);
    return y - x;
  }

template <class T=TFixed>
  inline T fxsq(T x) {
    return x == fxm::info<T>.nan.val?x:x*x;
  }
template <> inline double fxsq(double x) {return sq(x);}    
template <> inline float fxsq(float x) {return sq(x);}    

template <class T=TFixed>
  T fxsqrt(T x) {
  long n = fxm::info<T>.frac;
  if (x == fxm::info<T>.nan.val) return x;
  T r = 0.5 * x;
  T p = 0.0;
  if ( x < 0) return fxm::info<T>.nan.val;
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

template <class T=TFixed>
  T fxhypot(T x, T y=1.0) {
  T ax,ay;
  if (x == fxm::info<T>.nan.val) return x;
  if (y == fxm::info<T>.nan.val) return y;
  ax=fxabs(x);
  ay=fxabs(y);
  fxmaxs(ax, ay);
  if (ax == 0.0) return 0.0;
  ay /= ax;
  return ax*fxsqrt(1.0+ay*ay);
}    
template <> inline double fxhypot(double x, double y=1.0) {return hypot(x,y);}    
template <> inline float fxhypot(float x, float y=1.0) {return hypot(x,y);}    

template <class T=TFixed>
  T fxnorm(const T x[], uint8_t n=3) {
    uint8_t i,m;
    T nrm = (T)1.0, val = (T)0.0;
    for( i=0; i<n; i++ ) {
      if ( fxisnan(x[i]) ) return x[i];
      if ( fxabs(x[i]) <= val ) continue;
      val = fxabs(x[i]);
      m=i;
    }
    if ( val==0.0 ) return 0.0;
    for( i=0; i<n; i++ ) if ( i!=m )
      nrm += fxsq(x[i]/val);
    return val*fxsqrt(nrm);
  }    


template <class T=TFixed>
  T fxcbrt (T x) {
    T r=x / 3.0, e = 1.0;  
    long n = fxm::info<T>.frac;
    if (x == fxm::info<T>.nan.val) return x;
    if ( x < 0) return fxm::info<T>.nan.val;
    while (e != 0.0 && n) {
      e = ((x / r) / r - r) / 3.0;
      r += e;
      if (fxabs(e)<1.0) n--;
    };
    return(r);
};
template <> inline double fxcbrt(double x) {return cbrt(x);}    
template <> inline float fxcbrt(float x) {return cbrt(x);}    


template <class T=TFixed>
  T fxfloor(T x) {
     T r[8];
     FIXEDMATH_MAX* p = (FIXEDMATH_MAX*)r;
     FIXEDMATH_MAX mask = -1;
     if (x == fxm::info<T>.nan.val) return x;
     mask = -1;
     mask <<= fxm::info<T>.frac;
     r[0] = x;
     *p&=mask;     
     return r[0];
  }  
template <> inline double fxfloor(double x) {return floor(x);}    
template <> inline float fxfloor(float x) {return floor(x);}    

template <class T=TFixed>
  T fxceil(T x) {
    T r = fxfloor(x);
    if (x == r) return(r);
    return r + static_cast<T>(1);
  }
template <> inline double fxceil(double x) {return ceil(x);}    
template <> inline float fxceil(float x) {return ceil(x);}    

template <class T=TFixed>
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

template <class T=TFixed>
  T fxtrunc(T x) {
    return x<0?fxceil(x):fxfloor(x);
  }
template <> inline double fxtrunc(double x) {return trunc(x);}    
template <> inline float fxtrunc(float x) {return trunc(x);}    

template <class T=TFixed>
  T fxdiv(T x, T y) {
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val || y == 0) return fxm::info<T>.nan.val;
    return fxtrunc(x/y);
  }

template <class T=TFixed>
  inline T fxmod(T x, T y) {
    T r = fxdiv(x,y);
    if (r == fxm::info<T>.nan.val) return r;
    return x - r*y;  
  }
template <> inline double fxmod(double x, double y) {return fmod(x,y);}    
template <> inline float fxmod(float x, float y) {return fmod(x, y);}    

template <class T=TFixed>
  T fxmodf(T x, T* iptr) {
    if (x == fxm::info<T>.nan.val) return x;
    *iptr = fxtrunc(x);
    return x - *iptr;  
  }
template <> inline double fxmodf(double x, double* iptr) {return modf(x,iptr);}    

template <class T=TFixed>
  T fxfrexp(T x, int* expn) {
    int n = 0;
    if (x==0.0 || x == fxm::info<T>.nan.val) {
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

template <class T=TFixed>
  T fxrandom(T x, T y = 0.0) {
    T r=fxmins(x, y);
    r *= static_cast<T>(random(0x7FFFFFFF)/2147483647.0);
    return x+r;
  }

template <class T=TFixed>
  T fxsin (T x) {
    const T A  = 1.27323954;
    const T B  = 0.405284735;
    const T C  = 0.255;
    T r;
    if (x == fxm::info<T>.nan.val) return x;
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

template <class T=TFixed>
  inline T fxcos (T x) {
    return x == fxm::info<T>.nan.val?x:fxsin(x+HALF_PI);
  }
template <> inline double fxcos(double x) {return cos(x);}    
template <> inline float fxcos(float x) {return cos(x);}    

template <class T=TFixed>
  inline T fxtan (T x) {
    return x == fxm::info<T>.nan.val || fxmod(x+HALF_PI,static_cast<T>(PI))==0.0?fxm::info<T>.nan.val:fxsin(x)/fxsin(x+HALF_PI);
  }
template <> inline double fxtan(double x) {return tan(x);}    
template <> inline float fxtan(float x) {return tan(x);}    

template <class T=TFixed>
  inline T fxcot (T x) {
    return x == fxm::info<T>.nan.val || fxmod(x,static_cast<T>(PI))==0.0?fxm::info<T>.nan.val:fxsin(x+HALF_PI)/fxsin(x);
  }
template <> inline double fxcot(double x) {return 1.0/tan(x);}    
template <> inline float fxcot(float x) {return 1.0/tan(x);}    

template <class T=TFixed>
  T fxatan2 (T y, T x) {
    T t0, t1, t2, t3, t4, t5;
    if (x == fxm::info<T>.nan.val) return x;
    if (y == fxm::info<T>.nan.val) return y;
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

template <class T=TFixed>
  inline T fxatan(T x) {
    return fxatan2(x,static_cast<T>(1.0));
  }
template <> inline double fxatan(double x) {return atan(x);}    
template <> inline float fxatan(float x) {return atan(x);}    

template <class T=TFixed>
  T fxasin(T x) {
    T r = -0.0187293;
    bool n=false;
    if (x == fxm::info<T>.nan.val) return x;
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

template <class T=TFixed>
  T fxacos(T x) {
    T r = -0.0187293;
    bool n=false;
    if (x == fxm::info<T>.nan.val) return x;
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
  
template <class T=TFixed>
  T fxexp (T x) {
    static T f[6] = { 0.0,\
      7.389056098930650227230427460575,\
      54.598150033144239078110261202861,\
      2980.9579870417282747435920994529,\
      8886110.5205078726367630237407815,\
      78962960182680.695160978022635108};
    int n;    
    T r = 1.0, t = x;
    if (x == fxm::info<T>.nan.val) return x;
    if (x == 0.0) return 1.0;
    if (x<0.0) return 1.0/fxexp(-x);
    if (x>=1.0) {
      if (f[0] == 0.0) {
        f[0] = 2.7182818284590452353602874713527;
        for (n=1; n<6; n++) if (f[n]<=f[n-1]) break; 
        for (n=n; n<6; n++) f[n] = fxm::info<T>.nan.val;       
      }
      n = 0;
      while (t>=1.0) {
        if (fxmod(t,(T)2.0)>=1.0)
          if (f[n] == fxm::info<T>.nan.val)
            return fxm::info<T>.nan.val;
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

template <class T=TFixed>
  T fxlog2(T x) {
    int expn;
    if (x == fxm::info<T>.nan.val || x<=0.0) return fxm::info<T>.nan.val;
    x = fxfrexp(x, &expn);
    return expn + x*(2.8854 - 0.7213*x) - 2.1640;
  }

template <class T=TFixed>
  inline T fxlog(T x) {
    T r=fxlog2(x);
    return r == fxm::info<T>.nan.val ? r:0.69314718056*r;
  }
template <> inline double fxlog(double x) {return log(x);}    
template <> inline float fxlog(float x) {return log(x);}

template <class T=TFixed>
  inline T fxlog10(T x) {
    T r=fxlog2(x);
    return r == fxm::info<T>.nan.val ? r:0.30103*r;
  }
template <> inline double fxlog10(double x) {return log10(x);}    
template <> inline float fxlog10(float x) {return log10(x);}

template <class T=TFixed>
  inline T fxpow (T x, T y) {
    T t = 1.0;
    T b = x;
    T p = 0.0;
    bool ovf = fxabs(x)>1.0;
    if ( x == fxm::info<T>.nan.val ) return x;
    if ( y == fxm::info<T>.nan.val ) return y;
    if ( x == 0.0 ) return 0.0;
    if ( y == 0.0 || x == 1.0 ) return 1.0;
    if (y<0) { 
      t = fxpow(x,-y);
      return t==0.0 ? fxm::info<T>.nan.val : 1.0/t;
    }  
    while ( y >= 1.0 && t != 0.0 ) {
      if (ovf && fxabs(b)<p) return fxm::info<T>.nan.val; 
      if ( fxmod(y,static_cast<T>(2)) >= static_cast<T>(1) ) t*=b;
      if (ovf) p = fxabs(b);
      b *= b;
      y /= 2.0;
    }
    if ( y == 0.0 || t == 0.0 ) return t;
    if ( x < 0.0 ) return fxm::info<T>.nan.val;
    return t*fxexp(y*fxlog(x));
  }
template <> inline double fxpow(double x, double y) {return pow(x,y);}    
template <> inline float fxpow(float x, float y) {return pow(x,y);}

template <class T=TFixed>
  inline T fxsinh(T x) {
    return ( x == fxm::info<T>.nan.val ) ? fxm::info<T>.nan.val : (fxexp(x) - fxexp(-x)) / 2.0;
  }
template <> inline double fxsinh(double x) {return sinh(x);}    
template <> inline float fxsinh(float x) {return sinh(x);}

template <class T=TFixed>
  inline T fxcosh(T x) {
    return ( x == fxm::info<T>.nan.val ) ? fxm::info<T>.nan.val : (fxexp(x) + fxexp(-x)) / 2.0;
  }
template <> inline double fxcosh(double x) {return cosh(x);}    
template <> inline float fxcosh(float x) {return cosh(x);}

template <class T=TFixed>
  inline T fxtanh(T x) {
    if ( x == fxm::info<T>.nan.val ) return x;
    x=fxexp(2*x);
    return (x - 1.0)/(x + 1.0);
  }
template <> inline double fxtanh(double x) {return tanh(x);}    
template <> inline float fxtanh(float x) {return cosh(x);}

template <class T=TFixed>
  inline T fxconj(T x) { return x;  }
  
#ifndef FX_PRECISION
#define FX_PRECISION 2
#endif
#ifndef FX_WIDTH
#define FX_WIDTH 32
#endif

#define endl eol()
uint8_t fxPrecision = FX_PRECISION;
uint8_t fxWidth = 0;
uint8_t fxIndent = 0;
uint8_t fxPos = 0;
#define endl eol()
class eol {
  public:
  uint8_t n;
  inline eol():n(1) {}
  inline eol(uint8_t lf):n(lf) {}
};
class spc {
  public:
  uint8_t n;
  inline spc():n(0) {}
  inline spc(uint8_t s):n(s) {}
};
class setprecision {
  public:
  inline setprecision(uint8_t p) {if (p!=255) fxPrecision=p;}
};
class setw : public setprecision {
  public:
  inline setw(uint8_t w, uint8_t p=255) : setprecision(p) {if (w!=255) fxWidth=w;}
};
class indent : public setw {
  public:
  indent(uint8_t i=255, uint8_t w=255, uint8_t p=255):setw(255,255) { if (p==255 && i<255) { fxWidth=i; if (w!=255) fxPrecision=w; fxIndent = fxPos;} else if (i==255) fxIndent = fxPos; else fxIndent = i; }
};

Stream& operator << (Stream& s, const char x);

template <class T> Stream& operator << (Stream& s, const T &x) {
  if (fxPos < fxIndent) s << spc(fxIndent-fxPos);
  if (fxWidth) {
    char buf[FX_WIDTH];
    dtostrf((double)x, fxWidth, fxPrecision, buf);
    fxPos += s.print(buf);
  } else 
    fxPos += s.print((float)x,fxPrecision);
  return s; 
};


template <> inline Stream& operator << (Stream& s, const setprecision &x) { return s; };
template <> inline Stream& operator << (Stream& s, const setw &x) { return s; };
template <> inline Stream& operator << (Stream& s, const indent &x) { return s; };
template <>  Stream& operator << (Stream& s, const eol &x) { fxPrecision = FX_PRECISION; fxWidth = 0; fxPos = 0; fxIndent = 0; for (uint8_t i=0; i<x.n; i++) s.write('\n'); return s; };
template <>  Stream& operator << (Stream& s, const spc &x) { if (!x.n) {s.write(' '); fxPos++;} else do { s.write(' '); fxPos++;} while (fxPos % x.n); return s; };

Stream& operator << (Stream& s, char x) {
  if (fxPos < fxIndent) s << spc(fxIndent-fxPos);
  if (fxWidth) {
    char buf[FX_WIDTH];
    memset(buf,32,fxWidth-1);
    buf[fxWidth] = x;
    buf[fxWidth+1] = x;
    fxPos += s.print(buf);
  } else { 
    s.write(x);
    fxPos++;
  }
  if (x == '\n' || x == '\r') fxPos = 0;
  return s;   
}

Stream& operator << (Stream& s, const char f[]) {
  while( char c=*(f++) ) { 
    s.write(c) ;
    if (c == '\n' || c == '\r') fxPos = 0; else fxPos++;
  }
  return s; 
};

Stream& operator << (Stream& s, const __FlashStringHelper* f) {
  char* buf=reinterpret_cast<const char*>(f);
  while(char c=pgm_read_byte(buf++)) {
    s.write(c) ;
    if (c == '\n' || c == '\r') fxPos = 0; else fxPos++;
  }  
  return s; 
};

#endif
