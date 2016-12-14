/*********************************************************************
 * convolve.h
 *********************************************************************/

#ifndef _CONVOLVE_H_
#define _CONVOLVE_H_

#include <ippi.h>
#include "klt.h"
#include "klt_util.h"

void _KLTToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  _KLT_FloatImage floatimg);

void _KLTRoiToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  _KLT_FloatImage floatimg,
  int nFeatures,
  IppiRect rois[]);

void _KLTPatchToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  IppiRect *roi,
  _KLT_FloatImage patch);

void _KLTComputeGradients(
  _KLT_FloatImage img,
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel, 
  float *sigma_last,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady);


void _KLTGetKernelWidths(
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel,
  //int *gauss_width,
  //int *gaussderiv_width,
  float *sigma_last);
  
void _KLTComputeSmoothedImage(
  _KLT_FloatImage img,
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel, 
  float *sigma_last,
  _KLT_FloatImage smooth);

void _KLTRoiComputeSmoothedImage(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth,
  int nFeatures,
  IppiRect rois[]);

void computeKernelWidths(
  float sigma, 
  int *gaussWidth, 
  int *gaussderivWidth);

IppiRect *_KLTComputePatch(
  int roi_width, 
  int roi_height, 
  KLT_Feature feature, 
  int ncols, 
  int nrows );



#endif
