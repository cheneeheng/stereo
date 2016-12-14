/*********************************************************************
 * convolve.c
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <math.h>
#include <stdlib.h>   /* malloc(), realloc() */
#include <string.h>
#ifdef HAVE_IPP
	#include <ippi.h>  
	#include <ippcore.h>
#endif


/* Our includes */
#include "base.h"
#include "error.h"
#include "convolve.h"
#include "klt_util.h"   /* printing */





/* Kernels */
//static ConvolutionKernel gauss_kernel;
//static ConvolutionKernel gaussderiv_kernel;
//static float sigma_last = -10.0;


/*********************************************************************
 * _KLTToFloatImage
 *
 * Given a pointer to image data (probably unsigned chars), copy
 * data to a float image.
 */

void _KLTToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  _KLT_FloatImage floatimg)
{

  /* Output image must be large enough to hold result */
  assert(floatimg->ncols >= ncols);
  assert(floatimg->nrows >= nrows);

  floatimg->ncols = ncols;
  floatimg->nrows = nrows;


#ifdef HAVE_IPP
  IppiSize roiSize={ncols,nrows};
  ippiConvert_8u32f_C1R((const Ipp8u*)img,ncols*sizeof(unsigned char),(Ipp32f *)(floatimg->data),
          ncols*sizeof(float), roiSize);
#else
  KLT_PixelType *ptrend = img + ncols*nrows;
  float *ptrout = floatimg->data;

  while (img < ptrend)  *ptrout++ = (float) *img++;
#endif
    
}



/*********************************************************************
 * _KLTPatchToFloatImage
 *
 * Given a pointer to image data (probably unsigned chars), copy
 * data to a float image, but only in a specific ROI
 */

void _KLTPatchToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  IppiRect *roi,
  _KLT_FloatImage patch)
{
   /* Output image must be large enough to hold result */
  //assert(floatimg->ncols >= ncols);
  //assert(floatimg->nrows >= nrows);
  //IppiSize roiSize={(int)(roi_width+1),(int)(roi_height+1)}; //+1 because of rounding errors - x,y are float
  
  
  //memset((void *)floatimg->data, 0, ncols*nrows*sizeof(float));

  int windowOffset=((int)roi->y)*ncols+((int)roi->x);
#ifdef HAVE_IPP
  IppiSize roiSize={roi->width,roi->height}; 
  ippiConvert_8u32f_C1R((const Ipp8u*)(img+windowOffset),
		  ncols*sizeof(unsigned char),
		  (Ipp32f *)(patch->data),
		  roiSize.width*sizeof(float), 
		  roiSize);
#else
  KLT_PixelType *ptrroi = img+windowOffset;
  KLT_PixelType *ptrrowend = ptrroi+roi->width;
  int rowDiff = ncols - roi->width;
  float *ptrout = patch->data;
  float *ptrfloatend = ptrout+roi->width*roi->height;
  while(ptrout<ptrfloatend)
  {
	  while (ptrroi < ptrrowend)  *ptrout++ = (float) *ptrroi++;
	  ptrroi+=rowDiff;
	  ptrrowend += ncols;
  }
#endif
  
}



/*********************************************************************
 * _computeKernels
 */

static void _computeKernels(
  float sigma,
  ConvolutionKernel *gauss,
  ConvolutionKernel *gaussderiv,
  float *sigma_last)
{
  const float factor = 0.01f;   /* for truncating tail */
  int i;

  assert(MAX_KERNEL_WIDTH % 2 == 1);
  assert(sigma >= 0.0);

  /* Compute kernels, and automatically determine widths */
  {
    const int hw = MAX_KERNEL_WIDTH / 2;
    float max_gauss = 1.0f, max_gaussderiv = (float) (sigma*exp(-0.5f));
	
    /* Compute gauss and deriv */
    for (i = -hw ; i <= hw ; i++)  {
      gauss->data[i+hw]      = (float) exp(-i*i / (2*sigma*sigma));
      gaussderiv->data[i+hw] = -i * gauss->data[i+hw];
    }

    /* Compute widths */
    gauss->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gauss->data[i+hw] / max_gauss) < factor ; 
         i++, gauss->width -= 2);
    gaussderiv->width = MAX_KERNEL_WIDTH;
    for (i = -hw ; fabs(gaussderiv->data[i+hw] / max_gaussderiv) < factor ; 
         i++, gaussderiv->width -= 2);
    if (gauss->width == MAX_KERNEL_WIDTH || 
        gaussderiv->width == MAX_KERNEL_WIDTH)
      KLTError("(_computeKernels) MAX_KERNEL_WIDTH %d is too small for "
               "a sigma of %f", MAX_KERNEL_WIDTH, sigma);
  }

  /* Shift if width less than MAX_KERNEL_WIDTH */
  for (i = 0 ; i < gauss->width ; i++)
    gauss->data[i] = gauss->data[i+(MAX_KERNEL_WIDTH-gauss->width)/2];
  for (i = 0 ; i < gaussderiv->width ; i++)
    gaussderiv->data[i] = gaussderiv->data[i+(MAX_KERNEL_WIDTH-gaussderiv->width)/2];
  /* Normalize gauss and deriv */
  {
    const int hw = gaussderiv->width / 2;
    float den;
			
    den = 0.0;
    for (i = 0 ; i < gauss->width ; i++)  den += gauss->data[i];
    for (i = 0 ; i < gauss->width ; i++)  gauss->data[i] /= den;
    den = 0.0;
    for (i = -hw ; i <= hw ; i++)  den -= i*gaussderiv->data[i+hw];
    for (i = -hw ; i <= hw ; i++)  gaussderiv->data[i+hw] /= den;
  }

  *sigma_last = sigma;
}
	

/*********************************************************************
 * _KLTGetKernelWidths
 *
 */

void _KLTGetKernelWidths(
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel,
  //int *gauss_width,
  //int *gaussderiv_width,
  float *sigma_last)
{
  _computeKernels(sigma, gauss_kernel, gaussderiv_kernel, sigma_last);
  //*gauss_width = gauss_kernel.width;
  //*gaussderiv_width = gaussderiv_kernel.width;
}


/*********************************************************************
 * _convolveImageHoriz
 */

static void _convolveImageHoriz(
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout)
{
  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

#ifdef HAVE_IPP
  IppiSize roisize_src={imgin->ncols,imgin->nrows};
  IppiSize roisize_kernel={kernel.width,1};
  memset(imgout->data,0,imgout->ncols*imgout->nrows*sizeof(float));
  ippiConvValid_32f_C1R((const Ipp32f*)imgin->data,
                       imgin->ncols*sizeof(float),
                       roisize_src,(const Ipp32f*)(kernel.data),
	               kernel.width*sizeof(float), roisize_kernel,
		       (Ipp32f*)imgout->data+((int)(kernel.width/2)),
		       imgout->ncols*sizeof(float));
#else
  float *ptrrow = imgin->data;           /* Points to row's first pixel */
  register float *ptrout = imgout->data, /* Points to next output pixel */
    *ppp;
  register float sum;
  register int radius = kernel.width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;
  /* For each row, do ... */
  for (j = 0 ; j < nrows ; j++)  {

    /* Zero leftmost columns */
    for (i = 0 ; i < radius ; i++)
      *ptrout++ = 0.0;

    /* Convolve middle columns with kernel */
    for ( ; i < ncols - radius ; i++)  {
      ppp = ptrrow + i - radius;
      sum = 0.0;
      for (k = kernel.width-1 ; k >= 0 ; k--)
        sum += *ppp++ * kernel.data[k];
      *ptrout++ = sum;
    }

    /* Zero rightmost columns */
    for ( ; i < ncols ; i++)
      *ptrout++ = 0.0;

    ptrrow += ncols;
  }
#endif
}




/*********************************************************************
 * _convolveImageVert
 */

static void _convolveImageVert(
  _KLT_FloatImage imgin,
  ConvolutionKernel kernel,
  _KLT_FloatImage imgout)
{
  /* Kernel width must be odd */
  assert(kernel.width % 2 == 1);

  /* Must read from and write to different images */
  assert(imgin != imgout);

  /* Output image must be large enough to hold result */
  assert(imgout->ncols >= imgin->ncols);
  assert(imgout->nrows >= imgin->nrows);

#ifdef HAVE_IPP
  IppiSize roisize_src={imgin->ncols,imgin->nrows};
  IppiSize roisize_kernel={1,kernel.width};
  memset(imgout->data,0,imgout->ncols*imgout->nrows*sizeof(float));
  ippiConvValid_32f_C1R((const Ipp32f*)imgin->data,
                       imgin->ncols*sizeof(float),
                       roisize_src,(const Ipp32f*)kernel.data,
	               sizeof(float), roisize_kernel,
		       (Ipp32f*)imgout->data+((int)(kernel.width/2))*
		           imgout->ncols, imgout->ncols*sizeof(float));
#else
  float *ptrcol = imgin->data;            /* Points to row's first pixel */
  register float *ptrout = imgout->data,  /* Points to next output pixel */
    *ppp;
  register float sum;
  register int radius = kernel.width / 2;
  register int ncols = imgin->ncols, nrows = imgin->nrows;
  register int i, j, k;

  /* For each column, do ... */
  for (i = 0 ; i < ncols ; i++)  {

    /* Zero topmost rows */
    for (j = 0 ; j < radius ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    /* Convolve middle rows with kernel */
    for ( ; j < nrows - radius ; j++)  {
      ppp = ptrcol + ncols * (j - radius);
      sum = 0.0;
      for (k = kernel.width-1 ; k >= 0 ; k--)  {
        sum += *ppp * kernel.data[k];
        ppp += ncols;
      }
      *ptrout = sum;
      ptrout += ncols;
    }

    /* Zero bottommost rows */
    for ( ; j < nrows ; j++)  {
      *ptrout = 0.0;
      ptrout += ncols;
    }

    ptrcol++;
    ptrout -= nrows * ncols - 1;
  }  
#endif
}




/*********************************************************************
 * _convolveSeparate
 */

static void _convolveSeparate(
  _KLT_FloatImage imgin,
  ConvolutionKernel horiz_kernel,
  ConvolutionKernel vert_kernel,
  _KLT_FloatImage imgout)
{
  /* Create temporary image */
  _KLT_FloatImage tmpimg;
  tmpimg = _KLTCreateFloatImage(imgin->ncols, imgin->nrows); 
  
  /* Do convolution */
  _convolveImageHoriz(imgin, horiz_kernel, tmpimg);

  _convolveImageVert(tmpimg, vert_kernel, imgout);

  /* Free memory */
  _KLTFreeFloatImage(tmpimg);
}


/*********************************************************************
 * _convolveRoiSeparate
 */


	
/*********************************************************************
 * _KLTComputeGradients
 */

void _KLTComputeGradients(
  _KLT_FloatImage img,
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel,  
  float *sigma_last,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady)
{
				
  /* Output images must be large enough to hold result */
  assert(gradx->ncols >= img->ncols);
  assert(gradx->nrows >= img->nrows);
  assert(grady->ncols >= img->ncols);
  assert(grady->nrows >= img->nrows);

  /* Compute kernels, if necessary */
  if (fabs(sigma - *sigma_last) > 0.05)
    _computeKernels(sigma, gauss_kernel, gaussderiv_kernel, sigma_last);
	
  _convolveSeparate(img, *gaussderiv_kernel, *gauss_kernel, gradx); 		
  _convolveSeparate(img, *gauss_kernel, *gaussderiv_kernel, grady);

}



	

/*********************************************************************
 * _KLTComputeSmoothedImage
 */

void _KLTComputeSmoothedImage(
  _KLT_FloatImage img,
  float sigma,
  ConvolutionKernel *gauss_kernel,
  ConvolutionKernel *gaussderiv_kernel, 
  float *sigma_last,
  _KLT_FloatImage smooth)
{
  /* Output image must be large enough to hold result */
  assert(smooth->ncols >= img->ncols);
  assert(smooth->nrows >= img->nrows);

  /* Compute kernel, if necessary; gauss_deriv is not used */
  if (fabs(sigma - *sigma_last) > 0.05)
    _computeKernels(sigma, gauss_kernel, gaussderiv_kernel, sigma_last);

  _convolveSeparate(img, *gauss_kernel, *gauss_kernel, smooth);

}





/* Computes the ROI for the feature patch
 * 
 * 
 */
IppiRect *_KLTComputePatch(
					int roi_width, 
					int roi_height, 
					KLT_Feature feature, 
					int ncols, 
					int nrows )
{
	int half_roi_width  = roi_width/2;
	int half_roi_height  = roi_height/2;
	IppiRect *roi;
	
	roi = (IppiRect *)malloc(sizeof(IppiRect));
	if( (int)feature->x < 0 || (int)feature->x >=ncols ||
	    (int)feature->y < 0 || (int)feature->y >=nrows )
	{
		roi->x=0;
		roi->y=0;
		roi->width=0;
		roi->height=0;
	}
	else
	{
		roi->x=(int)feature->x-half_roi_width;
		roi->y=(int)feature->y-half_roi_height;
		roi->width=roi_width;
		roi->height=roi_height;
		if (roi->x < 0)
			roi->width-=roi->x, roi->x = 0;
		if (roi->y < 0)
			roi->height-=roi->y, roi->y = 0;
		if (roi->x+roi->width > ncols)
			roi->width = ncols - roi->x;
		if (roi->y+roi->height > nrows)
			roi->height = nrows - roi->y;
	}
	
	return roi;
}
