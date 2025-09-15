#ifndef _XSAT_H_
#define _XSAT_H_
//MUST BE NEAGATIVE!!
#define EPSILON 0
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

double DLE(double x,double y);
double DLT(double x,double y);
double DGT(double x,double y);
double DGE(double x,double y);
double DEQ(double x,double y);
double DNE(double x,double y);

// Float32 functions
double DLE_f32(float x, float y);
double DLT_f32(float x, float y);
double DGT_f32(float x, float y);
double DGE_f32(float x, float y);
double DEQ_f32(float x, float y);
double DNE_f32(float x, float y);
double ulp_f32(float x, float y);

double  BAND(double x,double y);
double  BOR(double x,double y);
float TR32(double x);
double MAX(double a, double b);

static inline uint64_t ordered_bits(double d) {
  uint64_t u; memcpy(&u, &d, sizeof u);
  return (u & (1ULL<<63)) ? ~u : (u ^ (1ULL<<63));
}

double ulp(double x, double y) {
  if (isnan(x) || isnan(y)) return NAN;
  if (isinf(x) && isinf(y)) {
    return (signbit(x) != signbit(y)) ? INFINITY : NAN;
  }
  if ((isinf(x) && !isinf(y)) || (isinf(y) && !isinf(x))) {
    return INFINITY;
  }
  if (x == y) return 0.0;
  uint64_t a = ordered_bits(x), b = ordered_bits(y);
  uint64_t d = (a >= b) ? (a - b) : (b - a);
  return (double)d;
}

// ============= FLOAT32 SUPPORT =============
static inline uint32_t ordered_bits_f32(float f) {
  uint32_t u;
  memcpy(&u, &f, sizeof u);
  return (u & (1U<<31)) ? ~u : (u ^ (1U<<31));
}

// Float32 ULP distance
double ulp_f32(float x, float y) {
  if (isnan(x) || isnan(y)) return NAN;
  if (isinf(x) && isinf(y)) {
    return (signbit(x) != signbit(y)) ? INFINITY : NAN;
  }
  if ((isinf(x) && !isinf(y)) || (isinf(y) && !isinf(x))) {
    return INFINITY;
  }
  if (x == y) return 0.0;
  uint32_t a = ordered_bits_f32(x), b = ordered_bits_f32(y);
  uint32_t d = (a >= b) ? (a - b) : (b - a);
  return (double)d;  // Return as double for consistency
}

// Float32-specific comparison functions
double DLE_f32(float x, float y) {
  return x <= y ? 0.0 : ulp_f32(x, y);
}

double DLT_f32(float x, float y) {
  return x < y ? 0.0 : ulp_f32(x, y) + 1;
}

double DGE_f32(float x, float y) {
  return DLE_f32(y, x);
}

double DGT_f32(float x, float y) {
  return DLT_f32(y, x);
}

double DEQ_f32(float x, float y) {
  return ulp_f32(x, y);
}

double DNE_f32(float x, float y) {
  return (x == y) ? 1.0 : 0.0;
}

double DLE(double x,double y){
    return x<=y?0.0:ulp(x,y);
}
double DLT(double x,double y){
  return x<y?0.0:ulp(x,y)+1;
}

double DGE(double x,double y)  {
   return DLE(y,x);
}
double  DGT(double x,double y)  {
  return DLT(y,x);
}
double DEQ(double x, double y){
  return ulp(x,y);
}
double  DNE(double x,double y) {
  return  (x==y)?1.0:0.0;
}

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
double MAX(double a, double b) {
  return (((a) > (b)) ? (a) : (b));
}
double  BAND(double x,double y){return x+y;}
double  BOR(double x,double y){return x*y;}

//#endif
float TR32(double x){
  return (float)x;
}
#endif
