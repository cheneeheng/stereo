/*********************************************************************
 * klt_util.h
 *********************************************************************/

#ifndef _KLT_UTIL_H_
#define _KLT_UTIL_H_

#include <ippi.h>
//#include "pyramid.h"


#define MAX_KERNEL_WIDTH 	71


typedef struct  {
  int ncols;
  int nrows;
  float *data;
}  _KLT_FloatImageRec, *_KLT_FloatImage;

typedef struct  {
  float xDisp;
  float yDisp;
}  _KLT_FloatPixelDispRec, *_KLT_FloatPixelDisp;

typedef struct  {
  int width;
  float data[MAX_KERNEL_WIDTH];
}  ConvolutionKernel;


float _KLTInterpolate(  
  float x, 
  float y, 
  _KLT_FloatImage img);

int _KLTTrackFeature(
  float x1,  
  float y1,
  float *x2, 
  float *y2,
  _KLT_FloatImage img1, 
  _KLT_FloatImage gradx1,
  _KLT_FloatImage grady1,
  _KLT_FloatImage img2, 
  _KLT_FloatImage gradx2,
  _KLT_FloatImage grady2,
  int width,           
  int height,
  int max_iterations,
  float small,         
  float th,            
  float max_residue,   
  int lighting_insensitive,
  float *residue);

_KLT_FloatImage _KLTCreateFloatImage(
  int ncols, 
  int nrows);

void _KLTFreeFloatImage(
  _KLT_FloatImage);
	
void _KLTPrintSubFloatImage(
  _KLT_FloatImage floatimg,
  int x0, int y0,
  int width, int height);

void _KLTWriteFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename);

/* for affine mapping */
void _KLTWriteAbsFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename,float scale);

void _KLTWritePatchToPPM(
  _KLT_FloatImage img,
  char *filename, 
  KLT_locType x,
  KLT_locType y,
  int patch_width,
  int patch_height,
  short colorChannel);
  
float _minEigenvalue(float gxx, float gxy, float gyy);


#endif


