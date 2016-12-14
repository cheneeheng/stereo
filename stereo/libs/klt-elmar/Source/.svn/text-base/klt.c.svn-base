/*********************************************************************
 * klt.c
 *
 * Kanade-Lucas-Tomasi tracker
 *********************************************************************/

/* Standard includes */
#include <assert.h>
#include <math.h>    /* logf() */
#include <stdlib.h>  /* malloc() */

/* Our includes */
#include "base.h"
#include "convolve.h"
#include "error.h"
#include "klt.h"
#include "pyramid.h"


static const int mindist = 10;
static const int window_size = 7;
static const int min_eigenvalue = 1;
static const float min_determinant = 0.01f;
static const float min_displacement = 0.1f;
static const int max_patchDisplacement = 10;
static const int max_iterations = 10;
static const float max_residue = 10.0f;
static const float grad_sigma = 1.0f;
static const float smooth_sigma_fact = 0.1f;
static const float pyramid_sigma_fact = 0.9f;
static const int sequentialMode = 0;
static const KLT_BOOL lighting_insensitive = FALSE;
/* for affine mapping*/
static const int affineConsistencyCheck = -1;
static const int affine_window_size = 15;
static const int affine_max_iterations = 10;
static const float affine_max_residue = 10.0;
static const float affine_min_displacement = 0.02f;
static const float affine_max_displacement_differ = 1.5f;

static const KLT_BOOL smoothBeforeSelecting = TRUE;
static const KLT_BOOL bindToCorner = FALSE;
static const KLT_BOOL writeInternalImages = FALSE;
static const int search_range = 15;
static const int nSkippedPixels = 0;

//static const KLT_BOOL linearPrediction = FALSE;
static const KLT_BOOL patchMode = FALSE;


//for stereo tracking
static const int stereo_window_size = 7;
static const float stereo_min_determinant = 0.01f;
static const float stereo_min_displacement = 0.1f;
static const int stereo_max_iterations = 10;
static const float stereo_max_residue = 10.0f;
static const KLT_BOOL stereo_lighting_insensitive = FALSE;
static const int stereo_nPyramidLevels = 1;
static const int stereo_yOffset = 0;
static const unsigned int stereo_yOffsetTolerance = 5;
static const unsigned int stereo_xOffset = 50;
static const int stereo_min_pixelDisplacement = 5;
static const int stereo_max_pixelDisplacement = 200;
static const unsigned int stereo_max_patchDisplacement = 5;
static const int stereo_writeInternalImages=0;
static const KLT_BOOL stereo_rejectMultipleCorrespondences = TRUE;
static const float stereo_max_distanceFromEpipolarLine = 0.1;
//static const int stereoMode = 1; //tracking features


extern int KLT_verbose;


/*********************************************************************
 * _createArray2D
 *
 * Creates a two-dimensional array.
 *
 * INPUTS
 * ncols:      no. of columns
 * nrows:      no. of rows
 * nbytes:     no. of bytes per entry
 *
 * RETURNS
 * Pointer to an array.  Must be coerced.
 *
 * EXAMPLE
 * char **ar;
 * ar = (char **) createArray2D(8, 5, sizeof(char));
 */

static void** _createArray2D(int ncols, int nrows, int nbytes)
{
  char **tt;
  int i;

  tt = (char **) malloc(nrows * sizeof(void *) +
                        ncols * nrows * nbytes);
  if (tt == NULL)
    KLTError("(createArray2D) Out of memory");

  for (i = 0 ; i < nrows ; i++)
    tt[i] = ((char *) tt) + (nrows * sizeof(void *) +
                             i * ncols * nbytes);

  return((void **) tt);
}


/*********************************************************************
 * KLTCreateTrackingContext
 *
 */

KLT_TrackingContext KLTCreateTrackingContext()
{
  KLT_TrackingContext tc;

  /* Allocate memory */
  tc = (KLT_TrackingContext)  malloc(sizeof(KLT_TrackingContextRec));

  /* Set values to default values */
  tc->mindist = mindist;
  tc->window_width = window_size;
  tc->window_height = window_size;
  tc->sequentialMode = sequentialMode;
  tc->smoothBeforeSelecting = smoothBeforeSelecting;
  tc->bindToCorner = bindToCorner;
  tc->writeInternalImages = writeInternalImages;
  tc->lighting_insensitive = lighting_insensitive;
  tc->min_eigenvalue = min_eigenvalue;
  tc->min_determinant = min_determinant;
  tc->max_iterations = max_iterations;
  tc->min_displacement = min_displacement;
  tc->max_patchDisplacement = max_patchDisplacement;
  tc->max_residue = max_residue;
  tc->grad_sigma = grad_sigma;
  tc->smooth_sigma_fact = smooth_sigma_fact;
  tc->pyramid_sigma_fact = pyramid_sigma_fact;
  tc->nSkippedPixels = nSkippedPixels;
  //tc->patchlist_last = NULL;
  tc->pyramid_last = NULL;
  tc->pyramid_last_gradx = NULL;
  tc->pyramid_last_grady = NULL;
  /* for affine mapping */
  tc->affineConsistencyCheck = affineConsistencyCheck;
  tc->affine_window_width = affine_window_size;
  tc->affine_window_height = affine_window_size;
  tc->affine_max_iterations = affine_max_iterations;
  tc->affine_max_residue = affine_max_residue;
  tc->affine_min_displacement = affine_min_displacement;
  tc->affine_max_displacement_differ = affine_max_displacement_differ;

  tc->pixelDisplacements_last = NULL;
  tc->patchMode = patchMode;

  /* Change nPyramidLevels and subsampling */
  KLTChangeTCPyramid(tc, search_range);

  /* Update border, which is dependent upon  */
  /* smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling */
  KLTUpdateTCBorder(tc);

  return(tc);
}



/*********************************************************************
 * KLTCreateTrackingContext
 *
 */

KLT_StereoTrackingContext KLTCreateStereoTrackingContext()
{
  KLT_StereoTrackingContext stc;

  /* Allocate memory */
  stc = (KLT_StereoTrackingContext)  malloc(sizeof(KLT_StereoTrackingContextRec));

  /* Set values to default values */
  stc->window_width = stereo_window_size;
  stc->window_height = stereo_window_size;
  stc->lighting_insensitive = stereo_lighting_insensitive;
  stc->min_determinant = stereo_min_determinant;
  stc->max_iterations = stereo_max_iterations;
  stc->min_displacement = stereo_min_displacement;
  stc->max_residue = stereo_max_residue;
  stc->max_patchDisplacement = stereo_max_patchDisplacement;

  stc->nPyramidLevels = stereo_nPyramidLevels;

  stc->yOffset = stereo_yOffset;
  stc->yOffsetTolerance = stereo_yOffsetTolerance;
  stc->xOffset = stereo_xOffset;
  stc->min_pixelDisplacement = stereo_min_pixelDisplacement;
  stc->max_pixelDisplacement = stereo_max_pixelDisplacement;

  stc->writeInternalImages = stereo_writeInternalImages;

  stc->rejectMultipleCorrespondences = stereo_rejectMultipleCorrespondences;

  stc->get_distanceFromEpipolarLine=NULL;
  stc->max_distanceFromEpipolarLine= stereo_max_distanceFromEpipolarLine;
//  stc->stereoMode = stereoMode;

  return(stc);
}

/*********************************************************************
 * KLTCreateFeatureList
 *
 */

KLT_FeatureList KLTCreateFeatureList(
  int nFeatures)
{
  KLT_FeatureList fl;
  KLT_Feature first;
  int nbytes = sizeof(KLT_FeatureListRec) +
    nFeatures * sizeof(KLT_Feature) +
    nFeatures * sizeof(KLT_FeatureRec);
  int i;

  /* Allocate memory for feature list */
  fl = (KLT_FeatureList)  malloc(nbytes);

  /* Set parameters */
  fl->nFeatures = nFeatures;

  /* Set pointers */
  fl->feature = (KLT_Feature *) (fl + 1);
  first = (KLT_Feature) (fl->feature + nFeatures);
  for (i = 0 ; i < nFeatures ; i++) {
    fl->feature[i] = first + i;
    fl->feature[i]->patch = NULL;
    fl->feature[i]->aff_img = NULL;           /* initialization fixed by Sinisa Segvic */
    fl->feature[i]->aff_img_gradx = NULL;
    fl->feature[i]->aff_img_grady = NULL;
  }
  /* Return feature list */
  return(fl);
}


/*********************************************************************
 * KLTCreateFeatureHistory
 *
 */

KLT_FeatureHistory KLTCreateFeatureHistory(
  int nFrames)
{
  KLT_FeatureHistory fh;
  KLT_Feature first;
  int nbytes = sizeof(KLT_FeatureHistoryRec) +
    nFrames * sizeof(KLT_Feature) +
    nFrames * sizeof(KLT_FeatureRec);
  int i;

  /* Allocate memory for feature history */
  fh = (KLT_FeatureHistory)  malloc(nbytes);

  /* Set parameters */
  fh->nFrames = nFrames;

  /* Set pointers */
  fh->feature = (KLT_Feature *) (fh + 1);
  first = (KLT_Feature) (fh->feature + nFrames);
  for (i = 0 ; i < nFrames ; i++)
    fh->feature[i] = first + i;

  /* Return feature history */
  return(fh);
}


/*********************************************************************
 * KLTCreateFeatureTable
 *
 */

KLT_FeatureTable KLTCreateFeatureTable(
  int nFrames,
  int nFeatures)
{
  KLT_FeatureTable ft;
  KLT_Feature first;
  int nbytes = sizeof(KLT_FeatureTableRec);
  int i, j;

  /* Allocate memory for feature history */
  ft = (KLT_FeatureTable)  malloc(nbytes);

  /* Set parameters */
  ft->nFrames = nFrames;
  ft->nFeatures = nFeatures;

  /* Set pointers */
  ft->feature = (KLT_Feature **)
    _createArray2D(nFrames, nFeatures, sizeof(KLT_Feature));
  first = (KLT_Feature) malloc(nFrames * nFeatures * sizeof(KLT_FeatureRec));
  for (j = 0 ; j < nFeatures ; j++)
    for (i = 0 ; i < nFrames ; i++)
      ft->feature[j][i] = first + j*nFrames + i;

  /* Return feature table */
  return(ft);
}


/*********************************************************************
 * KLTPrintTrackingContext
 */

void KLTPrintTrackingContext(
  KLT_TrackingContext tc)
{
  fprintf(stderr, "\n\nTracking context:\n\n");
  fprintf(stderr, "\tmindist = %d\n", tc->mindist);
  fprintf(stderr, "\twindow_width = %d\n", tc->window_width);
  fprintf(stderr, "\twindow_height = %d\n", tc->window_height);
  fprintf(stderr, "\tsequentialMode = %d\n", tc->sequentialMode);
  fprintf(stderr, "\tsmoothBeforeSelecting = %s\n",
          tc->smoothBeforeSelecting ? "TRUE" : "FALSE");
  fprintf(stderr, "\twriteInternalImages = %s\n",
          tc->writeInternalImages ? "TRUE" : "FALSE");

  fprintf(stderr, "\tmin_eigenvalue = %d\n", tc->min_eigenvalue);
  fprintf(stderr, "\tmin_determinant = %f\n", tc->min_determinant);
  fprintf(stderr, "\tmin_displacement = %f\n", tc->min_displacement);
  fprintf(stderr, "\tmax_iterations = %d\n", tc->max_iterations);
  fprintf(stderr, "\tmax_residue = %f\n", tc->max_residue);
  fprintf(stderr, "\tgrad_sigma = %f\n", tc->grad_sigma);
  fprintf(stderr, "\tsmooth_sigma_fact = %f\n", tc->smooth_sigma_fact);
  fprintf(stderr, "\tpyramid_sigma_fact = %f\n", tc->pyramid_sigma_fact);
  fprintf(stderr, "\tnSkippedPixels = %d\n", tc->nSkippedPixels);
  fprintf(stderr, "\tborderx = %d\n", tc->borderx);
  fprintf(stderr, "\tbordery = %d\n", tc->bordery);
  fprintf(stderr, "\tnPyramidLevels = %d\n", tc->nPyramidLevels);
  fprintf(stderr, "\tsubsampling = %d\n", tc->subsampling);
  fprintf(stderr, "\tpatchMode = %d\n", tc->patchMode);

  //fprintf(stderr, "\n\tpatchlist_last = %s\n", (tc->patchlist_last!=NULL) ?
  //        "points to old patchlist" : "NULL");

  fprintf(stderr, "\n\tpixelDisplacements_last = %s\n", (tc->pixelDisplacements_last!=NULL) ?
          "points to old pixel displacement" : "NULL");

  fprintf(stderr, "\n\tpyramid_last = %s\n", (tc->pyramid_last!=NULL) ?
          "points to old image" : "NULL");
  fprintf(stderr, "\tpyramid_last_gradx = %s\n",
          (tc->pyramid_last_gradx!=NULL) ?
          "points to old image" : "NULL");
  fprintf(stderr, "\tpyramid_last_grady = %s\n",
          (tc->pyramid_last_grady!=NULL) ?
          "points to old image" : "NULL");

  fprintf(stderr, "\n\n");
}


/*********************************************************************
 * KLTChangeTCPyramid
 *
 */

void KLTChangeTCPyramid(
  KLT_TrackingContext tc,
  int search_range)
{
  float window_halfwidth;
  float subsampling;

  /* Check window size (and correct if necessary) */
  if (tc->window_width % 2 != 1) {
    tc->window_width = tc->window_width+1;
    KLTWarning("(KLTChangeTCPyramid) Window width must be odd.  "
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height % 2 != 1) {
    tc->window_height = tc->window_height+1;
    KLTWarning("(KLTChangeTCPyramid) Window height must be odd.  "
               "Changing to %d.\n", tc->window_height);
  }
  if (tc->window_width < 3) {
    tc->window_width = 3;
    KLTWarning("(KLTChangeTCPyramid) Window width must be at least three.  \n"
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height < 3) {
    tc->window_height = 3;
    KLTWarning("(KLTChangeTCPyramid) Window height must be at least three.  \n"
               "Changing to %d.\n", tc->window_height);
  }
  window_halfwidth = min(tc->window_width,tc->window_height)/2.0f;

  subsampling = ((float) search_range) / window_halfwidth;

  if (subsampling < 1.0)  {		/* 1.0 = 0+1 */
    tc->nPyramidLevels = 1;
  } else if (subsampling <= 3.0)  {	/* 3.0 = 2+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 2;
  } else if (subsampling <= 5.0)  {	/* 5.0 = 4+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 4;
  } else if (subsampling <= 9.0)  {	/* 9.0 = 8+1 */
    tc->nPyramidLevels = 2;
    tc->subsampling = 8;
  } else {
    /* The following lines are derived from the formula:
       search_range =
       window_halfwidth * \sum_{i=0}^{nPyramidLevels-1} 8^i,
       which is the same as:
       search_range =
       window_halfwidth * (8^nPyramidLevels - 1)/(8 - 1).
       Then, the value is rounded up to the nearest integer. */
    float val = (float) (log(7.0*subsampling+1.0)/log(8.0));
    tc->nPyramidLevels = (int) (val + 0.99);
    tc->subsampling = 8;
  }

  tc->max_patchDisplacement=search_range;
}


/*********************************************************************
 * NOTE:  Manually must ensure consistency with _KLTComputePyramid()
 */

static float _pyramidSigma(
  KLT_TrackingContext tc)
{
  return (tc->pyramid_sigma_fact * tc->subsampling);
}


/*********************************************************************
 * Updates border, which is dependent upon
 * smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling
 */

void KLTUpdateTCBorder(
  KLT_TrackingContext tc)
{
  float val;
  int pyramid_gauss_hw;
  int smooth_gauss_hw;
  //int gauss_width, gaussderiv_width;
  int num_levels = tc->nPyramidLevels;
  int n_invalid_pixels;
  int window_hw;
  int ss = tc->subsampling;
  int ss_power;
  int border;
  int i;
  ConvolutionKernel gauss_kernel;
  ConvolutionKernel gaussderiv_kernel;
  float sigma_last = -10.0;

  /* Check window size (and correct if necessary) */
  if (tc->window_width % 2 != 1) {
    tc->window_width = tc->window_width+1;
    KLTWarning("(KLTUpdateTCBorder) Window width must be odd.  "
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height % 2 != 1) {
    tc->window_height = tc->window_height+1;
    KLTWarning("(KLTUpdateTCBorder) Window height must be odd.  "
               "Changing to %d.\n", tc->window_height);
  }
  if (tc->window_width < 3) {
    tc->window_width = 3;
    KLTWarning("(KLTUpdateTCBorder) Window width must be at least three.  \n"
               "Changing to %d.\n", tc->window_width);
  }
  if (tc->window_height < 3) {
    tc->window_height = 3;
    KLTWarning("(KLTUpdateTCBorder) Window height must be at least three.  \n"
               "Changing to %d.\n", tc->window_height);
  }
  window_hw = max(tc->window_width, tc->window_height)/2;

  /* Find widths of convolution windows */
  _KLTGetKernelWidths(_KLTComputeSmoothSigma(tc),&gauss_kernel, &gaussderiv_kernel, &sigma_last);
  //gauss_width = gauss_kernel.width;
  //gaussderiv_width = gaussderiv_kernel.width;
  //smooth_gauss_hw = gauss_width/2;
  smooth_gauss_hw = gauss_kernel.width/2;

  _KLTGetKernelWidths(_pyramidSigma(tc),&gauss_kernel, &gaussderiv_kernel, &sigma_last);
  //gauss_width = gauss_kernel.width;
  //gaussderiv_width = gaussderiv_kernel.width;
  //pyramid_gauss_hw = gauss_width/2;
  pyramid_gauss_hw = gauss_kernel.width/2;

  /* Compute the # of invalid pixels at each level of the pyramid.
     n_invalid_pixels is computed with respect to the ith level
     of the pyramid.  So, e.g., if n_invalid_pixels = 5 after
     the first iteration, then there are 5 invalid pixels in
     level 1, which translated means 5*subsampling invalid pixels
     in the original level 0. */
  n_invalid_pixels = smooth_gauss_hw;
  for (i = 1 ; i < num_levels ; i++)  {
    val = ((float) n_invalid_pixels + pyramid_gauss_hw) / ss;
    n_invalid_pixels = (int) (val + 0.99999);  /* Round up */
  }

  /* ss_power = ss^(num_levels-1) */
  ss_power = 1;
  for (i = 1 ; i < num_levels ; i++)
    ss_power *= ss;

  /* Compute border by translating invalid pixels back into */
  /* original image */
  border = (n_invalid_pixels + window_hw) * ss_power;

  tc->borderx = border;
  tc->bordery = border;
}


/*********************************************************************
 * KLTFreeTrackingContext
 * KLTFreeStereoTrackingContext
 * KLTFreeFeatureList
 * KLTFreeFeatureHistory
 * KLTFreeFeatureTable
 */

void KLTFreeTrackingContext(
  KLT_TrackingContext tc)
{
  /*if (tc->patchlist_last)
  {
	int j=0;
	for(j=0;j<tc->patchlist_last->nFeatures;j++)
	{
		if(tc->patchlist_last->patches[j].pyramid)
		{
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid);
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid_gradx);
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid_grady);
			free(tc->patchlist_last->patches[j].roi);
		}
	}
	free(tc->patchlist_last->patches);
	free(tc->patchlist_last);
  }*/
  if(tc->pixelDisplacements_last)
	  free(tc->pixelDisplacements_last);
  if (tc->pyramid_last)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last);
  if (tc->pyramid_last_gradx)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_gradx);
  if (tc->pyramid_last_grady)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_grady);

  free(tc);
  tc=NULL;
}

void KLTFreeStereoTrackingContext(
  KLT_StereoTrackingContext stc)
{
  free(stc);
  stc=NULL;
}


void KLTFreeFeatureList(
  KLT_FeatureList fl)
{
  /* for affine mapping */
  int indx;
  for (indx = 0 ; indx < fl->nFeatures ; indx++)  {
    /* free image and gradient  */
    _KLTFreeFloatImage(fl->feature[indx]->aff_img);
    _KLTFreeFloatImage(fl->feature[indx]->aff_img_gradx);
    _KLTFreeFloatImage(fl->feature[indx]->aff_img_grady);
    fl->feature[indx]->aff_img = NULL;
    fl->feature[indx]->aff_img_gradx = NULL;
    fl->feature[indx]->aff_img_grady = NULL;
    _KLTFreePatch(fl->feature[indx]->patch);
  }

  free(fl);
}

void KLTFreeFeatureHistory(
  KLT_FeatureHistory fh)
{
  free(fh);
}

void KLTFreeFeatureTable(
  KLT_FeatureTable ft)
{
  free(ft->feature[0][0]);  /* this plugs a memory leak found by Stefan Wachter */
  free(ft->feature);
  free(ft);
}




/*********************************************************************
 * KLTStopSequentialMode
 */

void KLTFreeStoredImages(
  KLT_TrackingContext tc)
{
  /*if(tc->patchlist_last)
  {
	  int j=0;
	  for(j=0;j<tc->patchlist_last->nFeatures;j++)
	  {
		if(tc->patchlist_last->patches[j].pyramid)
		{
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid);
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid_gradx);
			_KLTFreePyramid(tc->patchlist_last->patches[j].pyramid_grady);
			free(tc->patchlist_last->patches[j].roi);
		}
	  }
	  free(tc->patchlist_last->patches);
	  free(tc->patchlist_last);
	  tc->patchlist_last=NULL;
  }*/
  if(tc->pixelDisplacements_last)
  {
  	  free(tc->pixelDisplacements_last);
  	  tc->pixelDisplacements_last=NULL;
  }

  if (tc->pyramid_last)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last),tc->pyramid_last = NULL;
  if (tc->pyramid_last_gradx)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_gradx),tc->pyramid_last_gradx = NULL;
  if (tc->pyramid_last_grady)
    _KLTFreePyramid((_KLT_Pyramid) tc->pyramid_last_grady),tc->pyramid_last_grady = NULL;
}


/*********************************************************************
 * KLTStopSequentialMode
 */

void KLTStopSequentialMode(
  KLT_TrackingContext tc)
{
  tc->sequentialMode = FALSE;

  KLTFreeStoredImages(tc);
}


/*********************************************************************
 * KLTCountRemainingFeatures
 */

int KLTCountRemainingFeatures(
  KLT_FeatureList fl)
{
  int count = 0;
  int i;

  for (i = 0 ; i < fl->nFeatures ; i++)
    if (fl->feature[i]->val >= 0)
      count++;

  return count;
}


/*********************************************************************
 * KLTCountRemainingFeatures
 * counts also features where the number of iterations was not sufficient
 * -> softer restriction
 */

int KLTExtCountRemainingFeatures(
  KLT_FeatureList fl)
{
  int count = 0;
  int i;

  for (i = 0 ; i < fl->nFeatures ; i++)
    if (fl->feature[i]->val >= 0 || fl->feature[i]->val==-3)
      count++;

  return count;
}


static void _copyFeature(KLT_Feature f1, KLT_Feature f2)
{
	f2->x=f1->x;
	f2->y=f1->y;
	f2->val=f1->val;
	/* for affine mapping */
	f2->aff_img=f1->aff_img;
	f2->aff_img_gradx=f1->aff_img_gradx;
	f2->aff_img_grady=f1->aff_img_grady;
	f2->aff_x=f1->aff_x;
	f2->aff_y=f1->aff_y;
	f2->aff_Axx=f1->aff_Axx;
	f2->aff_Ayx=f1->aff_Ayx;
	f2->aff_Axy=f1->aff_Axy;
	f2->aff_Ayy=f1->aff_Ayy;
}

int KLTGetRemainingFeatures(
  KLT_FeatureList fl, KLT_FeatureList fl_valid)
{
  int count = 0;
  int i;

  for (i = 0 ; i < fl->nFeatures ; i++)
  {
    if (fl->feature[i]->val >= 0)
    {
    	_copyFeature(fl->feature[i], fl_valid->feature[count]);
        count++;
    }
  }

  return count;
}

/*********************************************************************
 * KLTSetVerbosity
 */

void KLTSetVerbosity(
  int verbosity)
{
  KLT_verbose = verbosity;
}



/*********************************************************************
 * KLTMoveFeature
 */
void KLTMoveFeature(KLT_Feature featureSrc, KLT_Feature featureDest)
{
	featureDest->x = featureSrc->x;
	featureDest->y = featureSrc->y;
	featureDest->val = featureSrc->val;
	if(featureDest->patch)
		_KLTFreePatch(featureDest->patch);
	featureDest->patch = featureSrc->patch;
	featureSrc->patch = NULL;

	featureDest->aff_Axx = featureSrc->aff_Axx;
	featureDest->aff_Axy = featureSrc->aff_Axy;
	featureDest->aff_Ayy = featureSrc->aff_Ayy;
	featureDest->aff_Ayx = featureSrc->aff_Ayx;
	if(featureDest->aff_img)
		_KLTFreeFloatImage(featureDest->aff_img);
	featureDest->aff_img = featureSrc->aff_img;
	featureSrc->aff_img = NULL;
	if(featureDest->aff_img_gradx)
		_KLTFreeFloatImage(featureDest->aff_img_gradx);
	featureDest->aff_img_gradx = featureSrc->aff_img_gradx;
	featureSrc->aff_img_gradx = NULL;
	if(featureDest->aff_img_grady)
		_KLTFreeFloatImage(featureDest->aff_img_grady);
	featureDest->aff_img_grady = featureSrc->aff_img_grady;
	featureSrc->aff_img_grady = NULL;
	featureDest->aff_x = featureSrc->aff_x;
	featureDest->aff_y = featureSrc->aff_y;
}


/*********************************************************************
 * KLTMoveFeatureList
 */
void KLTMoveFeatureList(KLT_FeatureList featurelistSrc, KLT_FeatureList featurelistDest)
{
	int j;
	if(!featurelistDest || featurelistDest->nFeatures != featurelistSrc->nFeatures)
	{
		KLTWarning("(KLTCopyFeatureList) Feature cannot be copied cause destination "
				"list doesn't exit or is not of the same size");
		return;
	}

	for(j=0;j<featurelistSrc->nFeatures; j++)
	{
		KLTMoveFeature(featurelistSrc->feature[j], featurelistDest->feature[j]);
	}
}



/*********************************************************************
 * KLTCopyFeatureList
 */
int KLTCopyFeatureList(KLT_FeatureList featurelistSrc, KLT_FeatureList featurelistDest)
{
	int j;
	if(!featurelistDest || featurelistDest->nFeatures != featurelistSrc->nFeatures)
	{
		KLTWarning("(KLTCopyFeatureList) Feature cannot be copied cause destination "
				"list doesn't exit or is not of the same size");
		return -1;
	}

	for(j=0;j<featurelistSrc->nFeatures; j++)
	{
		KLT_Feature featureSrc = featurelistSrc->feature[j];
		KLT_Feature featureDest = featurelistDest->feature[j];

		featureDest->x = featureSrc->x;
		featureDest->y = featureSrc->y;
		featureDest->val = featureSrc->val;

		featureDest->aff_Axx = featureSrc->aff_Axx;
		featureDest->aff_Axy = featureSrc->aff_Axy;
		featureDest->aff_Ayy = featureSrc->aff_Ayy;
		featureDest->aff_Ayx = featureSrc->aff_Ayx;
		featureDest->aff_x = featureSrc->aff_x;
		featureDest->aff_y = featureSrc->aff_y;
	}

	return 0;
}


/*********************************************************************
 * KLTExtCopyFeatureList
 */
/*void KLTExtCopyFeatureList(KLT_FeatureList featurelistSrc, KLT_FeatureList featurelistDest, KLT_TrackingContext tc, int ncols, int nrows)
{
	if(KLTCopyFeatureList(featurelistSrc,featurelistDest))
		return;
}*/

void KLTExtPropagateFeatureList(KLT_FeatureList featurelistSrc, KLT_FeatureList featurelistDest, KLT_TrackingContext tc, int ncols, int nrows)
{
	if(KLTCopyFeatureList(featurelistSrc,featurelistDest))
		return;

	if(tc->sequentialMode==2 && tc->pixelDisplacements_last != NULL) //linear feature propagation
	{
		int nFeatures=featurelistDest->nFeatures;

		if (KLT_verbose >= 1)  {
			fprintf(stderr,  "(KLT_EXT) Copying with linear feature propagation.\n");
		}
		int j;

		for( j=0;j<nFeatures; j++)
		{
			KLT_Feature feature=featurelistDest->feature[j];
			if(feature->val==0)
			{
				float help = feature->x + tc->pixelDisplacements_last[j].xDisp;
				//float hWin = (float)tc->window_width/2.0;
				if(help>=tc->borderx && help<(ncols-tc->borderx))
					feature->x = help;

				help = feature->y + tc->pixelDisplacements_last[j].yDisp;
				//hWin = (float)tc->window_height/2.0;
				if(help>=tc->bordery && help<(nrows-tc->bordery))
					feature->y = help;
			}
		}
        //KLTWriteFeatureListToPPM(fl_help, img2, ncols, nrows, "feat2_guessed.ppm");

	}
}


void KLTCopyTrackingContext(
  KLT_TrackingContext tc_from,
  KLT_TrackingContext tc_to)
{

  tc_to->mindist = tc_from->mindist;
  tc_to->window_width = tc_from->window_width;
  tc_to->window_height = tc_from->window_height;
  tc_to->sequentialMode = tc_from->sequentialMode;
  tc_to->smoothBeforeSelecting = tc_from->smoothBeforeSelecting;
  tc_to->writeInternalImages = tc_from->writeInternalImages;
  tc_to->lighting_insensitive = tc_from->lighting_insensitive;
  tc_to->min_eigenvalue = tc_from->min_eigenvalue;
  tc_to->min_determinant = tc_from->min_determinant;
  tc_to->max_iterations = tc_from->max_iterations;
  tc_to->min_displacement = tc_from->min_displacement;
  tc_to->max_patchDisplacement = tc_from->max_patchDisplacement;
  tc_to->max_residue = tc_from->max_residue;
  tc_to->grad_sigma = tc_from->grad_sigma;
  tc_to->smooth_sigma_fact = tc_from->smooth_sigma_fact;
  tc_to->pyramid_sigma_fact = tc_from->pyramid_sigma_fact;
  tc_to->nSkippedPixels = tc_from->nSkippedPixels;
  //tc->patchlist_last = NULL;
  tc_to->pyramid_last = NULL;
  tc_to->pyramid_last_gradx = NULL;
  tc_to->pyramid_last_grady = NULL;
  /* for affine mapping */
  tc_to->affineConsistencyCheck = tc_from->affineConsistencyCheck;
  tc_to->affine_window_width = tc_from->affine_window_width;
  tc_to->affine_window_height = tc_from->affine_window_height;
  tc_to->affine_max_iterations = tc_from->affine_max_iterations;
  tc_to->affine_max_residue = tc_from->affine_max_residue;
  tc_to->affine_min_displacement = tc_from->affine_min_displacement;
  tc_to->affine_max_displacement_differ = tc_from->affine_max_displacement_differ;

  tc_to->pixelDisplacements_last = NULL;
  tc_to->patchMode = tc_from->patchMode;




}

