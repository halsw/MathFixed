/*
 * This file is part of the MathFixed library
 * Usage: Tests speed performance of its functions
 * Dependecies: FixedPoints library
 * 
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
 */
#include "MathFixed.h"
#include <FixedPoints.h>

//The number of integer bits
#define TFIXED_INT 7

//The number of fractional bits
#define TFIXED_FRAC 8

//Define the type used for calculations
#define TFIXED SFixed<TFIXED_INT, TFIXED_FRAC>

//The period of each cycle
#define PERIOD_MS 10

//The number of cycles for speed testing
#define CYCLES 1000

#define PROFILE(X,V) static long x##X = 0;\ 
 long t##X = micros();\
  X V;\
  t##X = micros() - t##X;\
  x##X += (t##X - x##X)/sample;\

#define PPRINT(X) Serial.print(#X ": ");\
  Serial.print(x##X);\
  Serial.println("us");\
  x##X = 0.0;

unsigned int wait() {
  static unsigned int load = 0;
  static unsigned long loopMS=PERIOD_MS;
  unsigned long now=millis();
  if (now>loopMS) {
    load = 100<<8;
    loopMS = now + PERIOD_MS;
  } else {
    load += (25600 - (loopMS-now)*(25600/PERIOD_MS) - load)>>3;
    while (loopMS>millis());
    loopMS += PERIOD_MS;
  }
  return(load);
}

void setup() {
  Serial.begin(115200);
  randomSeed(analogRead(0)); //assuming A0 is not connected
  fxibits<TFIXED>(TFIXED_INT); //Define integer bits in the library
  fxfbits<TFIXED>(TFIXED_FRAC);//Define fractional bits in the library
  Serial.print("MathFixed library test, please wait ");
  Serial.print((unsigned long)CYCLES*PERIOD_MS/1000);
  Serial.println("s between updates...");  
}

void loop() {
  static long sample = 1;  
  TFIXED x = fxrandom((TFIXED)-127,(TFIXED)127);//Generate first argument for function testing
  TFIXED y = fxrandom((TFIXED)127); //Generate second argument for function testing
  TFIXED r;
  int load, t;
  PROFILE(fxisnan,(x) );
  PROFILE(fxabs,  (x) );
  PROFILE(fxmin,  (x,y) );
  PROFILE(fxmax,  (x,y) );
  PROFILE(fxmins, (&x,&y) );
  PROFILE(fxmaxs, (&x,&y) );
  PROFILE(fxsq,   (x) );
  PROFILE(fxsqrt, (x) );
  PROFILE(fxhypot,(x,y) );
  PROFILE(fxcbrt, (x) );
  PROFILE(fxfloor,(x) );
  PROFILE(fxceil, (x) );
  PROFILE(fxround,(x) );
  PROFILE(fxtrunc,(x) );
  PROFILE(fxdiv,  (x,y) );
  PROFILE(fxmod,  (x,y) );
  PROFILE(fxmodf, (x,&r) );
  PROFILE(fxfrexp,(x,&t) );
  PROFILE(fxrandom,(x,y) );
  PROFILE(fxsin,  (x) );
  PROFILE(fxcos,  (x) );
  PROFILE(fxtan,  (x) );
  PROFILE(fxcot,  (x) );
  PROFILE(fxatan2,(x,y) );
  PROFILE(fxatan, (x) );
  PROFILE(fxasin, (x) );
  PROFILE(fxacos, (x) );
  PROFILE(fxexp, (x) );
  PROFILE(fxlog2, (y) );
  PROFILE(fxlog,  (y) );
  PROFILE(fxlog10,(y) );
  PROFILE(fxpow,  (y,x) );
  PROFILE(fxsinh, (x) );
  PROFILE(fxcosh, (x) );
  PROFILE(fxtanh, (x) );
  load=wait()>>8;
  if (sample++ > CYCLES) {
    sample = 1;
    Serial.print("Load(%):");
    Serial.println(load);
    PPRINT(fxisnan);
    PPRINT(fxabs);
    PPRINT(fxmax);
    PPRINT(fxmin);
    PPRINT(fxmaxs);
    PPRINT(fxmins);
    PPRINT(fxsq);
    PPRINT(fxsqrt);
    PPRINT(fxhypot);
    PPRINT(fxcbrt);
    PPRINT(fxfloor);
    PPRINT(fxceil);
    PPRINT(fxround);
    PPRINT(fxtrunc);
    PPRINT(fxdiv);
    PPRINT(fxmod);
    PPRINT(fxmodf);
    PPRINT(fxfrexp);
    PPRINT(fxrandom);
    PPRINT(fxsin);
    PPRINT(fxcos);
    PPRINT(fxtan);
    PPRINT(fxcot);
    PPRINT(fxatan2);
    PPRINT(fxatan);
    PPRINT(fxasin);
    PPRINT(fxacos);
    PPRINT(fxexp);
    PPRINT(fxlog2);
    PPRINT(fxlog);
    PPRINT(fxlog10);
    PPRINT(fxpow);
    PPRINT(fxsinh);
    PPRINT(fxcosh);
    PPRINT(fxtanh);
  }
}
