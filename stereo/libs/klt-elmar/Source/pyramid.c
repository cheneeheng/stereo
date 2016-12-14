/*********************************************************************
 * pyramid.c
 *
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <stdlib.h>		/* malloc() ? */
#include <string.h>		/* memset() ? */
#include <math.h>		/* */
#ifdef HAVE_IPP
	#include <ippi.h>  
	#include <ippcore.h>
#endif

/* Our includes */
#include "base.h"
#include "error.h"
#include "convolve.h"	/* for computing pyramid */
#include "pyramid.h"
#include "klt_util.h"   /* printing */
#include "klt.h" 


/*********************************************************************
 *
 */

_KLT_Pyramid _KLTCreatePyramid(
  int ncols,
  int nrows,
  int subsampling,
  int nlevels)
{
  _KLT_Pyramid pyramid;
  int nbytes = sizeof(_KLT_PyramidRec) +	
    nlevels * sizeof(_KLT_FloatImage *) +
    nlevels * sizeof(int) +
    nlevels * sizeof(int);
  int i;

  if (subsampling != 2 && subsampling != 4 && 
      subsampling != 8 && subsampling != 16 && subsampling != 32)
    KLTError("(_KLTCreatePyramid)  Pyramid's subsampling must "
             "be either 2, 4, 8, 16, or 32");

     
  /* Allocate memory for structure and set parameters */
  pyramid = (_KLT_Pyramid)  malloc(nbytes);
  if (pyramid == NULL)
    KLTError("(_KLTCreatePyramid)  Out of memory");
     
  /* Set parameters */
  pyramid->subsampling = subsampling;
  pyramid->nLevels = nlevels;
  pyramid->img = (_KLT_FloatImage *) (pyramid + 1);
  pyramid->ncols = (int *) (pyramid->img + nlevels);
  pyramid->nrows = (int *) (pyramid->ncols + nlevels);

  /* Allocate memory for each level of pyramid and assign pointers */
  for (i = 0 ; i < nlevels ; i++)  {
    pyramid->img[i] =  _KLTCreateFloatImage(ncols, nrows);
    pyramid->ncols[i] = ncols;  pyramid->nrows[i] = nrows;
    ncols /= subsampling;  nrows /= subsampling;
  }

  return pyramid;
}


/*********************************************************************
 *
 */

void _KLTFreePyramid(
  _KLT_Pyramid pyramid)
{
  int i;

  /* Free images */
  for (i = 0 ; i < pyramid->nLevels ; i++)
    _KLTFreeFloatImage(pyramid->img[i]);

  /* Free structure */
  free(pyramid);
}


_KLT_Patch _KLTCreatePatch()
{
	return (_KLT_Patch)malloc(sizeof(_KLT_PatchRec));
}

void _KLTDeletePatchImages(
  _KLT_Patch patch)
{
	if(patch->pyramid)
	{
		_KLTFreePyramid(patch->pyramid);
		_KLTFreePyramid(patch->pyramid_gradx);
		_KLTFreePyramid(patch->pyramid_grady);
		free(patch->roi);
	}
}

void _KLTFreePatch(
  _KLT_Patch patch)
{
	if(patch)
	{
		_KLTDeletePatchImages(patch);
		free(patch);
	}
}


/*********************************************************************
 * TODO: sometimes works, but not always??? IppiResize doesn't work properly!?!
 */
#ifdef HAVE_IPP
	void _KLTComputePyramid(
	  _KLT_FloatImage img, 
	  _KLT_Pyramid pyramid,
	  float sigma_fact,
	  ConvolutionKernel *gauss_kernel,
	  ConvolutionKernel *gaussderiv_kernel, 
	  float *sigma_last)
	{
	  _KLT_FloatImage currimg, tmpimg;
	  int ncols = img->ncols, nrows = img->nrows;
	  int subsampling = pyramid->subsampling;
	  float sigma = subsampling * sigma_fact;  // empirically determined
	  int oldncols;
	  int i;
	  float subfactor = 1.0/(float)subsampling;
		
	  if (subsampling != 2 && subsampling != 4 && 
	      subsampling != 8 && subsampling != 16 && subsampling != 32)
	    KLTError("(_KLTComputePyramid)  Pyramid's subsampling must "
	             "be either 2, 4, 8, 16, or 32");
	
	  assert(pyramid->ncols[0] == img->ncols);
	  assert(pyramid->nrows[0] == img->nrows);
	
	  // Copy original image to level 0 of pyramid 
	  memcpy(pyramid->img[0]->data, img->data, ncols*nrows*sizeof(float));

/*	    char fname[80];
	    sprintf(fname, "helpPyramida.pgm");
		_KLTWriteFloatImageToPGM(img, fname);
*/
	  currimg = img;
	  for (i = 1 ; i < pyramid->nLevels ; i++)  {
	    tmpimg = _KLTCreateFloatImage(ncols, nrows);
	    _KLTComputeSmoothedImage(currimg, sigma, gauss_kernel, gaussderiv_kernel, sigma_last, tmpimg);

/*	    sprintf(fname, "helpPyramidb%d.pgm", i);
		_KLTWriteFloatImageToPGM(tmpimg, fname);
*/	
	    // Subsample 
	    oldncols = ncols;
	    ncols /= subsampling;  nrows /= subsampling;
	    IppiSize size_src={tmpimg->ncols,tmpimg->nrows};
	    IppiSize size_dest={pyramid->img[i]->ncols,pyramid->img[i]->nrows};
	    //IppiRect roi={0,0,tmpimg->ncols-1,tmpimg->nrows-1};
	    IppiRect roi={0,0,tmpimg->ncols,tmpimg->nrows};
	//    fprintf(stderr,"\nRETURN : %s\n",ippGetStatusString(ippiResize_32f_C1R((const Ipp32f*)tmpimg->data,
	//        size_src,tmpimg->ncols*sizeof(float),roi,
	//	(Ipp32f*)pyramid->img[i]->data,pyramid->img[i]->ncols*sizeof(float),
	//	size_dest,subfactor,subfactor,IPPI_INTER_NN))); //nearest neighbor
	//    fflush(stderr);
	    ippiResize_32f_C1R((const Ipp32f*)tmpimg->data,
	            size_src,tmpimg->ncols*sizeof(float),roi,
	    	(Ipp32f*)pyramid->img[i]->data,pyramid->img[i]->ncols*sizeof(float),
	    	size_dest,subfactor,subfactor,IPPI_INTER_NN); //nearest neighbor
	
	    // Reassign current image 
	    currimg = pyramid->img[i];

/*	    sprintf(fname, "helpPyramidc%d.pgm", i);
	    _KLTWriteFloatImageToPGM(pyramid->img[i], fname);
*/					
	    _KLTFreeFloatImage(tmpimg);
	  }
	}
#else
	void _KLTComputePyramid(
	  _KLT_FloatImage img, 
	  _KLT_Pyramid pyramid,
	  float sigma_fact,
	  ConvolutionKernel *gauss_kernel,
	  ConvolutionKernel *gaussderiv_kernel, 
	  float *sigma_last)
	{
	  _KLT_FloatImage currimg, tmpimg;
	  int ncols = img->ncols, nrows = img->nrows;
	  int subsampling = pyramid->subsampling;
	  int subhalf = subsampling / 2;
	  float sigma = subsampling * sigma_fact;  // empirically determined 
	  int oldncols;
	  int i, x, y;
		
	  if (subsampling != 2 && subsampling != 4 && 
	      subsampling != 8 && subsampling != 16 && subsampling != 32)
	    KLTError("(_KLTComputePyramid)  Pyramid's subsampling must "
	             "be either 2, 4, 8, 16, or 32");
	
	  assert(pyramid->ncols[0] == img->ncols);
	  assert(pyramid->nrows[0] == img->nrows);
	
	  // Copy original image to level 0 of pyramid 
	  memcpy(pyramid->img[0]->data, img->data, ncols*nrows*sizeof(float));
	
	  currimg = img;
	  for (i = 1 ; i < pyramid->nLevels ; i++)  {
	    tmpimg = _KLTCreateFloatImage(ncols, nrows);
	    _KLTComputeSmoothedImage(currimg, sigma, gauss_kernel, gaussderiv_kernel, sigma_last, tmpimg);
	
	
	    // Subsample 
	    oldncols = ncols;
	    ncols /= subsampling;  nrows /= subsampling;
	    for (y = 0 ; y < nrows ; y++)
	      for (x = 0 ; x < ncols ; x++)
	        pyramid->img[i]->data[y*ncols+x] = 
	          tmpimg->data[(subsampling*y+subhalf)*oldncols +
	                      (subsampling*x+subhalf)];
	
	    // Reassign current image 
	    currimg = pyramid->img[i];
					
	    _KLTFreeFloatImage(tmpimg);
	  }
	}
#endif






