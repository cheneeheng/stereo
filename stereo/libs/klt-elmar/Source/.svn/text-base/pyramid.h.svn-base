/*********************************************************************
 * pyramid.h
 *********************************************************************/

#ifndef _PYRAMID_H_
#define _PYRAMID_H_

#include <ippi.h>
#include "klt_util.h"

typedef struct  {
  int subsampling;
  int nLevels;
  _KLT_FloatImage *img;
  int *ncols, *nrows;
}  _KLT_PyramidRec, *_KLT_Pyramid;


typedef struct  {
  _KLT_Pyramid pyramid;
  _KLT_Pyramid pyramid_gradx;
  _KLT_Pyramid pyramid_grady;
  IppiRect *roi;
}  _KLT_PatchRec, *_KLT_Patch;

typedef struct  {
  int nFeatures;
  _KLT_PatchRec *patches;
}  _KLT_PatchListRec, *_KLT_PatchList;

_KLT_Pyramid _KLTCreatePyramid(
  int ncols,
  int nrows,
  int subsampling,
  int nlevels);

void _KLTComputePyramid(
  _KLT_FloatImage floatimg, 
  _KLT_Pyramid pyramid,
  float sigma_fact,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel, 
  float *sigma_last);


void _KLTFreePyramid(
  _KLT_Pyramid pyramid);

_KLT_Patch _KLTCreatePatch();

void _KLTDeletePatchImages(
  _KLT_Patch patch);

void _KLTFreePatch(
  _KLT_Patch patch);


#endif
