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

// Float32-specific comparison functions
double DLE_f32(float x, float y) {
  return x <= y ? 0.0 : (x - y) * (x - y);
}

double DLT_f32(float x, float y) {
  return x < y ? 0.0 : (x - y) * (x - y) + 1;
}

double DGE_f32(float x, float y) {
  return DLE_f32(y, x);
}

double DGT_f32(float x, float y) {
  return DLT_f32(y, x);
}

double DEQ_f32(float x, float y) {
  return (x - y) * (x - y);
}

double DNE_f32(float x, float y) {
  return (x == y) ? 1.0 : 0.0;
}

double DLE(double x,double y){
  return x <= y ? 0.0 : (x - y) * (x - y);
}
double DLT(double x,double y){
  return x < y ? 0.0 : (x - y) * (x - y) + 1;
}

double DGE(double x,double y)  {
   return DLE(y,x);
}
double  DGT(double x,double y)  {
  return DLT(y,x);
}
double DEQ(double x, double y){
  return (x - y) * (x - y);
}
double  DNE(double x,double y) {
  return (x == y) ? 1.0 : 0.0;
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
