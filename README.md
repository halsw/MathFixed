# Fixed Point Math Library for Arduino and Teensy
A library with template functions that provides Math.h alternatives for fixed point types

Note that this library does not implement fixed point types but these are passed as parameters  
The library also redirects calls to the original Math.h functions if it is used with the **double** and **float** types

The library should be initialized for every type used using the **fxibits()** and **fxfbits()** functions to define respectively the integer and fractional number of bits of the type used.

It implements also a NaN value for every type that coincides with the +INF and -INF values that are not implemented.

The library isn't optimized yet for fast execution and has only been tested the **FixedPoints** library, which needs to be installed before running the example .ino file

## Functions
The implemented functions are:
**fxibits()** gets or sets the integer bits for a type
**fxfbits()** gets or sets the fractional bits for a type
**fxnan()** gets the representation of NaN
**fxisnan()** tests if argument is not a number
**fxisinf()** tests if argument is infinity, but here just a copy of **fxisnan()**
**fxabs()** the absolute value
**fxmax()** the greater of two numbers
**fxmin()** the lesser of two numbers
**fxmaxs()** sort two numbers descending and return difference
**fxmins()** sort two numbers ascending and return difference
**fxsq()** the square of a number
**fxsqrt()** the square root
**fxhypot()** the hypotenuse of two numbers
**fxcbrt()** the cubic root
**fxfloor()** the immediately smaller integer of a number
**fxceil()** the immediately larger integer of a number
**fxround()** the closest integer of a number
**fxtrunc()** the integer part of a number
**fxdiv()** the integer division of two numbers
**fxmod()** the remainder of the integer division 
**fxmodf()** get the fractional and integer part of a number
**fxfrexp()** get the mantissa and exponent(base 2) of a number
**fxrandom()** get a random number between two limits
**fxsin()** the sine
**fxcos()** the cosine
**fxtan()** the tangent
**fxcot()** the cotangent
**fxatan2()** the inverse tangent of the ratio of two numbers
**fxatan()** the inverse tangent
**fxasin()** the inverse sine
**fxacos()** the inverse cosine
**fxexp()** the natural exponential
**fxlog2()** the base 2 logarithm
**fxlog()** the natural logarithm
**fxlog10()** the base 10 logarithm
**fxpow()** raise a number to given exponent
**fxsinh()** the hyperbolic sine
**fxsinh()** the hyperbolic cosine
**fxtan()** the hyperbolic tangent

## Example File
There example .ino file only tests the above functions for speed and not accuracy  