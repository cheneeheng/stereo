/*********************************************************************
 * klt_util.c
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>  /* malloc() */
#include <math.h>		/* fabs() */
#include <ipp.h>

/* Our includes */
#include "base.h"
#include "error.h"
#include "pnmio.h"
#include "klt.h"
#include "klt_util.h"



/*********************************************************************
 * _KLTInterpolate
 * 
 * Given a point (x,y) in an image, computes the bilinear interpolated 
 * gray-level value of the point in the image.  
 */

float _KLTInterpolate(
  float x, 
  float y, 
  _KLT_FloatImage img)
{
  int xt = (int) x;  /* coordinates of top-left corner */
  int yt = (int) y;
  float ax = x - xt;
  float ay = y - yt;
  float *ptr = img->data + (img->ncols*yt) + xt;

#ifndef _DNDEBUG
  if (xt<0 || yt<0 || xt>=img->ncols-1 || yt>=img->nrows-1) {
    fprintf(stderr, "(xt,yt)=(%d,%d)  imgsize=(%d,%d)\n"
            "(x,y)=(%f,%f)  (ax,ay)=(%f,%f)\n",
            xt, yt, img->ncols, img->nrows, x, y, ax, ay);
    fflush(stderr);
  }
#endif

  assert (xt >= 0 && yt >= 0 && xt <= img->ncols - 2 && yt <= img->nrows - 2);

  return ( (1-ax) * (1-ay) * *ptr +
           ax   * (1-ay) * *(ptr+1) +
           (1-ax) *   ay   * *(ptr+(img->ncols)) +
           ax   *   ay   * *(ptr+(img->ncols)+1) );                      
}


/*********************************************************************/

float _KLTComputeSmoothSigma(
  KLT_TrackingContext tc)
{
  return (tc->smooth_sigma_fact * max(tc->window_width, tc->window_height));
}


/*********************************************************************
 * _KLTCreateFloatImage
 */

_KLT_FloatImage _KLTCreateFloatImage(
  int ncols,
  int nrows)
{
  _KLT_FloatImage floatimg;
  int nbytes = sizeof(_KLT_FloatImageRec) +  ncols * nrows * sizeof(float);// +  ncols * nrows * sizeof(float);
  floatimg = (_KLT_FloatImage) malloc(nbytes);//ippMalloc(nbytes);//malloc(nbytes);
  //floatimg.data=ippiMalloc_32f_C1(ncols,nrows,(ncols*sizeof(float)));  //you must use ippifree()!!
  if (floatimg == NULL)
    KLTError("(_KLTCreateFloatImage)  Out of memory");
  floatimg->ncols = ncols;
  floatimg->nrows = nrows;
  floatimg->data = (float *)  (floatimg + 1);

  return(floatimg);
}


/*********************************************************************
 * _KLTFreeFloatImage
 */

void _KLTFreeFloatImage(
  _KLT_FloatImage floatimg)
{
  free(floatimg);
}


/*********************************************************************
 * _KLTPrintSubFloatImage
 */

void _KLTPrintSubFloatImage(
  _KLT_FloatImage floatimg,
  int x0, int y0,
  int width, int height)
{
  int ncols = floatimg->ncols;
  int offset;
  int i, j;

  assert(x0 >= 0);
  assert(y0 >= 0);
  assert(x0 + width <= ncols);
  assert(y0 + height <= floatimg->nrows);

  fprintf(stderr, "\n");
  for (j = 0 ; j < height ; j++)  {
    for (i = 0 ; i < width ; i++)  {
      offset = (j+y0)*ncols + (i+x0);
      fprintf(stderr, "%6.2f ", *(floatimg->data + offset));
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}
	

/*********************************************************************
 * _KLTWriteFloatImageToPGM
 */

void _KLTWriteFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename)
{
  int npixs = img->ncols * img->nrows;
  float mmax = -999999.9f, mmin = 999999.9f;
  float fact;
  float *ptr;
  uchar *byteimg, *ptrout;
  int i;

  /* Calculate minimum and maximum values of float image */
  ptr = img->data;
  for (i = 0 ; i < npixs ; i++)  {
    mmax = max(mmax, *ptr);
    mmin = min(mmin, *ptr);
    ptr++;
  }
	
  /* Allocate memory to hold converted image */
  byteimg = (uchar *) malloc(npixs * sizeof(uchar));

  /* Convert image from float to uchar */
  fact = 255.0f / (mmax-mmin);
  ptr = img->data;
  ptrout = byteimg;
  for (i = 0 ; i < npixs ; i++)  {
    *ptrout++ = (uchar) ((*ptr++ - mmin) * fact);
  }

  /* Write uchar image to PGM */
  pgmWriteFile(filename, byteimg, img->ncols, img->nrows);

  /* Free memory */
  free(byteimg);
}


/*********************************************************************
 * _KLTWritePatchToPGM
 */

void _KLTWritePatchToPPM(
  _KLT_FloatImage img,
  char *filename, 
  KLT_locType x,
  KLT_locType y,
  int patch_width,
  int patch_height,
  short colorChannel) //0 red, 1 green, 2 blue
{
  int npixs = img->ncols * img->nrows;
  float mmax = -999999.9f, mmin = 999999.9f;
  float fact;
  float *ptr;
  uchar *byteimg, *ptrout;
  uchar *redimg, *grnimg, *bluimg;
  uchar *colorChannelP, *blackChannelP1, *blackChannelP2;
  int i;
  int xx, yy;
  int ncols = img->ncols;
  int nbytes = npixs * sizeof(uchar);
  int offset;
  
  /* Calculate minimum and maximum values of float image */
  ptr = img->data;
  for (i = 0 ; i < npixs ; i++)  {
    mmax = max(mmax, *ptr);
    mmin = min(mmin, *ptr);
    ptr++;
  }
	
  /* Allocate memory to hold converted image */
  byteimg = (uchar *) malloc(nbytes);

  /* Convert image from float to uchar */
  fact = 255.0f / (mmax-mmin);
  ptr = img->data;
  ptrout = byteimg;
  for (i = 0 ; i < npixs ; i++)  {
    *ptrout++ = (uchar) ((*ptr++ - mmin) * fact);
  }
  
  /* Allocate memory for component images */
  redimg = (uchar *)  malloc(nbytes);
  grnimg = (uchar *)  malloc(nbytes);
  bluimg = (uchar *)  malloc(nbytes);
  if (redimg == NULL || grnimg == NULL || bluimg == NULL)
    KLTError("(KLTWriteFeaturesToPPM)  Out of memory\n");

  /* Copy grey image to component images */
  if (sizeof(KLT_PixelType) != 1)
    KLTWarning("(KLTWriteFeaturesToPPM)  KLT_PixelType is not uchar");
  memcpy(redimg, byteimg, nbytes);
  memcpy(grnimg, byteimg, nbytes);
  memcpy(bluimg, byteimg, nbytes);

  switch(colorChannel)
  {
	case 1:
    	colorChannelP = grnimg;
    	blackChannelP1 = redimg;
    	blackChannelP2 = bluimg;
    	break;
  	case 2:
    	colorChannelP = bluimg;
    	blackChannelP1 = grnimg;
    	blackChannelP2 = redimg;
    	break;
  	case 0:
    default:
    	colorChannelP = redimg;
    	blackChannelP1 = grnimg;
    	blackChannelP2 = bluimg;
    	break;
  }
  
  
  //draw roi-rectangle
  //upper line
  xx=x-patch_width/2-1;
  yy=y-patch_height/2-1;
  offset= yy * ncols + xx;
  for(i=0; i<patch_width+2; i++, offset++){
	  if(offset>=0 && offset<nbytes){
		  *(redimg + offset)=255;
		  *(blackChannelP1 + offset)=0;
		  *(blackChannelP2 + offset)=0;
	  }
  }
  //side line left
  offset= yy * ncols + xx;
  for(i=0; i<patch_height; i++){
	  offset+= ncols;
	  if(offset>=0 && offset<nbytes){
		  *(colorChannelP + offset)=255;
		  *(blackChannelP1 + offset)=0;
		  *(blackChannelP2 + offset)=0;
  	  }
  }	  
  //side line right
  xx=x+patch_width/2+1;
  offset= yy * ncols + xx;
  for(i=0; i<patch_height; i++){
	  offset+= ncols;
	  if(offset>=0 && offset<nbytes){
		  *(colorChannelP + offset)=255;
		  *(blackChannelP1 + offset)=0;
		  *(blackChannelP2 + offset)=0;
      }
  }	  
  //lower line
  xx=x-patch_width/2-1;
  yy=y+patch_height/2+1;
  offset= yy * ncols + xx;
  for(i=0; i<patch_width+2; i++, offset++){
	  if(offset>=0 && offset<nbytes){
		  *(colorChannelP + offset)=255;
		  *(blackChannelP1 + offset)=0;
		  *(blackChannelP2 + offset)=0;
      }
  }

  
  //draw feature position
  offset= (int)y * ncols + (int)x;
  if(offset>=0 && offset<nbytes){
	  *(colorChannelP + offset)=255;
	  *(blackChannelP1 + offset)=0;
	  *(blackChannelP2 + offset)=0;
  }

  

  /* Write uchar image to PGM */
  //ppmWriteFileRGB(filename, byteimg, img->ncols, img->nrows);
  ppmWriteFileRGB(filename, redimg, grnimg, bluimg, img->ncols, img->nrows);

  /* Free memory */
  free(byteimg);
  free(redimg);
  free(grnimg);
  free(bluimg);
}
/*********************************************************************
 * _KLTWriteFloatImageToPGM
 */

void _KLTWriteAbsFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename,float scale)
{
  int npixs = img->ncols * img->nrows;
  float fact;
  float *ptr;
  uchar *byteimg, *ptrout;
  int i;
  float tmp;
	
  /* Allocate memory to hold converted image */
  byteimg = (uchar *) malloc(npixs * sizeof(uchar));

  /* Convert image from float to uchar */
  fact = 255.0f / scale;
  ptr = img->data;
  ptrout = byteimg;
  for (i = 0 ; i < npixs ; i++)  {
    tmp = (float) (fabs(*ptr++) * fact);
    if(tmp > 255.0) tmp = 255.0;
    *ptrout++ =  (uchar) tmp;
  }

  /* Write uchar image to PGM */
  pgmWriteFile(filename, byteimg, img->ncols, img->nrows);

  /* Free memory */
  free(byteimg);
}
