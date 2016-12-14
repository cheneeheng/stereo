/*********************************************************************
 * stereo_util.c
 *
 * utilities for stereo camera setups
 *********************************************************************/



#include <math.h>
#include <stdlib.h>  /* malloc(), qsort() */
#include <ippi.h> /*IppiRect*/
#include <assert.h>
#include <stdio.h>

#include "klt.h"
#include "klt_util.h"
#include "convolve.h"
#include "stereo_util.h"



#define MIN(a,b)(a<b?a:b)
#define MAX(a,b)(a>b?a:b)


extern int KLT_verbose;

/*********************************************************************
 * _compareFeatureLoc
 *
 * compare function for the qsort algorithm
 * compares two KLT_SortElement structs respectively to there y value
 */

static int _compareFeatureLoc(const void *a, const void *b)
{
  KLT_locType v1 =  ((KLT_SortElement) a)->y;
  KLT_locType v2 =  ((KLT_SortElement) b)->y;

  if (v1 > v2)  return(1);
  else if (v1 < v2)  return(-1);
  else return(0);
}


/*********************************************************************
 * _createSortList
 *
 * creates a KLT_SortList which is used to get a sorted feature list
 *
 */

static KLT_SortList _createSortList(
		int nFeatures)
{
	int nbytes = sizeof(KLT_SortListRec)+
				nFeatures*sizeof(KLT_SortElement)+
				nFeatures*sizeof(KLT_SortElementRec);
	int i;
	KLT_SortElement first;

	KLT_SortList dl=(KLT_SortList)malloc(nbytes);
	dl->nFeatures = nFeatures;

	/* Set pointers */
	dl->element = (KLT_SortElement*)(dl+1);
	first = (KLT_SortElement) (dl->element + nFeatures);
	for (i = 0 ; i < nFeatures ; i++) {
	    dl->element[i] = first + i;
    	dl->element[i]->index=-1;
    	dl->element[i]->y=-1;
    }

	return dl;
}



/*********************************************************************
 * _freeSortList
 *
 * frees the allocated memory of a KLT_SortList
 *
 */

static void _freeSortList(KLT_SortList dl)
{
	free(dl);
}



/*********************************************************************
 * _featureList2sortList
 *
 * copies the necessary values of a feature list into a sortlist
 *
 */

static void _featureList2sortList(KLT_FeatureList fl, KLT_SortList sl)
{
	int j;
	int nFeatures=sl->nFeatures;

	assert(nFeatures == fl->nFeatures);

	for(j=0;j<sl->nFeatures;j++)
	{
		sl->element[j]->index=j;
		sl->element[j]->y=fl->feature[j]->y;
	}


}


/*********************************************************************
 * _computePyramid
 *
 * computes the image pyramid necessary for stereo tracking for the
 * normal image mode
 *
 */

static void _computePyramid(
					KLT_PixelType *img,
					int ncols,
					int nrows,
					KLT_TrackingContext tc)
{
	float subsampling = (float) tc->subsampling;
	int nPyramidLevels=tc->nPyramidLevels;
	int i;
	ConvolutionKernel gauss_kernel;
	ConvolutionKernel gaussderiv_kernel;
	float sigma_last = -10.0;

	_KLT_FloatImage floatimg = _KLTCreateFloatImage(ncols, nrows);
	_KLT_FloatImage tmpimg = _KLTCreateFloatImage(ncols, nrows);
	_KLTToFloatImage(img, ncols, nrows, tmpimg);
	_KLTComputeSmoothedImage(tmpimg, _KLTComputeSmoothSigma(tc), &gauss_kernel, &gaussderiv_kernel, &sigma_last, floatimg);

	_KLT_Pyramid pyramid = _KLTCreatePyramid(ncols, nrows, (int) subsampling, nPyramidLevels);
	_KLTComputePyramid(floatimg, pyramid, tc->pyramid_sigma_fact, &gauss_kernel, &gaussderiv_kernel, &sigma_last);
	_KLT_Pyramid pyramid_gradx = _KLTCreatePyramid(ncols, nrows, (int) subsampling, nPyramidLevels);
	_KLT_Pyramid pyramid_grady = _KLTCreatePyramid(ncols, nrows, (int) subsampling, nPyramidLevels);
	for (i = 0 ; i < nPyramidLevels ; i++)
		_KLTComputeGradients(pyramid->img[i], tc->grad_sigma,
		&gauss_kernel, &gaussderiv_kernel, &sigma_last,
		pyramid_gradx->img[i],
		pyramid_grady->img[i]);


	if(tc->pyramid_last!=0)
		_KLTFreePyramid(tc->pyramid_last);
	if(tc->pyramid_last_gradx!=0)
		_KLTFreePyramid(tc->pyramid_last_gradx);
	if(tc->pyramid_last_grady!=0)
		_KLTFreePyramid(tc->pyramid_last_grady);

	tc->pyramid_last = pyramid;
	tc->pyramid_last_gradx = pyramid_gradx;
	tc->pyramid_last_grady = pyramid_grady;

	_KLTFreeFloatImage(tmpimg);
	_KLTFreeFloatImage(floatimg);
}



/*********************************************************************
 * _computePatchPyramid
 *
 * computes the image pyramid necessary for stereo tracking for the
 * patch mode
 *
 */

static void _computePatchPyramid(
					KLT_PixelType *img,
					KLT_FeatureList featurelist,
					int ncols,
					int nrows,
					int roi_width,
					int roi_height,
					KLT_TrackingContext tc)
{
	int indx;
	float subsampling = (float) tc->subsampling;
	int nPyramidLevels=tc->nPyramidLevels;
	int nFeatures=featurelist->nFeatures;
	ConvolutionKernel gauss_kernel;
	ConvolutionKernel gaussderiv_kernel;
	float sigma_last = -10.0;


	for(indx=0;indx<nFeatures; indx++)
	{
		_KLT_Patch featurePatch=featurelist->feature[indx]->patch;
		if(featurePatch)
			_KLTDeletePatchImages(featurePatch);
		else
		{
			featurePatch=_KLTCreatePatch();
			featurelist->feature[indx]->patch = featurePatch;
		}

		int i;
		IppiRect *roi;
		roi =_KLTComputePatch(roi_width, roi_height, featurelist->feature[indx], ncols, nrows);
		int width=roi->width;
		int height=roi->height;
		_KLT_FloatImage tmpimg = _KLTCreateFloatImage(width, height);
		_KLT_FloatImage floatimg = _KLTCreateFloatImage(width, height);
		_KLTPatchToFloatImage(img, ncols, nrows, roi, tmpimg);
		_KLTComputeSmoothedImage(tmpimg, _KLTComputeSmoothSigma(tc), &gauss_kernel, &gaussderiv_kernel, &sigma_last, floatimg);
		_KLT_Pyramid pyramid = _KLTCreatePyramid(width, height, (int) subsampling, nPyramidLevels);
		_KLTComputePyramid(floatimg, pyramid, tc->pyramid_sigma_fact, &gauss_kernel, &gaussderiv_kernel, &sigma_last);
		_KLT_Pyramid pyramid_gradx = _KLTCreatePyramid(width, height, (int) subsampling, nPyramidLevels);
		_KLT_Pyramid pyramid_grady = _KLTCreatePyramid(width, height, (int) subsampling, nPyramidLevels);
		for (i = 0 ; i < nPyramidLevels ; i++)
			_KLTComputeGradients(pyramid->img[i], tc->grad_sigma,
					&gauss_kernel, &gaussderiv_kernel, &sigma_last,
					pyramid_gradx->img[i],
					pyramid_grady->img[i]);



		featurePatch->pyramid = pyramid;
		featurePatch->pyramid_gradx = pyramid_gradx;
		featurePatch->pyramid_grady = pyramid_grady;
		featurePatch->roi = roi;

		_KLTFreeFloatImage(tmpimg);
		_KLTFreeFloatImage(floatimg);
	}

}



/*********************************************************************
 * _computeIntensityDifference
 *
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images.
 */

static void _computeIntensityDifferenceValue(
  _KLT_FloatImage img1,   /* images */
  _KLT_FloatImage img2,
  float x1, float y1,     /* center of window in 1st img */
  float x2, float y2,     /* center of window in 2nd img */
  int width, int height,  /* size of window */
  double *imgdiff)   /* output */
{
  register int hw = width/2, hh = height/2;
  float g1, g2;
  register int i, j;

  /* Compute values */
  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)  {
      g1 = _KLTInterpolate(x1+i, y1+j, img1);
      g2 = _KLTInterpolate(x2+i, y2+j, img2);
      *imgdiff += fabs(g1 - g2);
    }
}



/*********************************************************************
 * _computeIntensityDifferenceLightingInsensitive
 *
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images; normalizes for overall gain and bias.
 */

static void _computeIntensityDifferenceLightingInsensitiveValue(
  _KLT_FloatImage img1,   /* images */
  _KLT_FloatImage img2,
  float x1, float y1,     /* center of window in 1st img */
  float x2, float y2,     /* center of window in 2nd img */
  int width, int height,  /* size of window */
  double *imgdiff)   /* output */
{
  register int hw = width/2, hh = height/2;
  float g1, g2, sum1_squared = 0, sum2_squared = 0;
  register int i, j;

  float sum1 = 0, sum2 = 0;
  float mean1, mean2,alpha,belta;
  /* Compute values */
  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)  {
      g1 = _KLTInterpolate(x1+i, y1+j, img1);
      g2 = _KLTInterpolate(x2+i, y2+j, img2);
      sum1 += g1;    sum2 += g2;
      sum1_squared += g1*g1;
      sum2_squared += g2*g2;
   }
  mean1=sum1_squared/(width*height);
  mean2=sum2_squared/(width*height);
  alpha = (float) sqrt(mean1/mean2);
  mean1=sum1/(width*height);
  mean2=sum2/(width*height);
  belta = mean1-alpha*mean2;

  for (j = -hh ; j <= hh ; j++)
    for (i = -hw ; i <= hw ; i++)  {
      g1 = _KLTInterpolate(x1+i, y1+j, img1);
      g2 = _KLTInterpolate(x2+i, y2+j, img2);
      *imgdiff += fabs(g1- g2*alpha-belta);
    }
}


static void _writeInternalImages(	KLT_TrackingContext tc_left,
									KLT_TrackingContext tc_right,
									KLT_StereoTrackingContext stc,
									KLT_FeatureList fl_left,
									KLT_FeatureList fl_right)
{
	char fname[80];
	char dir[]="InternalImages/";
	int i, indx;
	for (i = 0 ; i < stc->nPyramidLevels ; i++)
	{
		if(tc_left->patchMode)
		{
			for (indx = 0 ; indx < fl_left->nFeatures ; indx++)
			{
				sprintf(fname, "%sl_p%d_kltimg_tf_i%d.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_left->feature[indx]->patch->pyramid->img[i], fname);
				sprintf(fname, "%sl_p%d_kltimg_tf_i%d_gx.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_left->feature[indx]->patch->pyramid_gradx->img[i], fname);
				sprintf(fname, "%sl_p%d_kltimg_tf_i%d_gy.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_left->feature[indx]->patch->pyramid_grady->img[i], fname);
			}

		}
		else
		{
			sprintf(fname, "%sl_kltimg_tf_i%d.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_left->pyramid_last)->img[i],fname);
			sprintf(fname, "%sl_kltimg_tf_i%d_gx.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_left->pyramid_last_gradx)->img[i],fname);
			sprintf(fname, "%sl_kltimg_tf_i%d_gy.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_left->pyramid_last_grady)->img[i],fname);
		}

		if(tc_right->patchMode)
		{
			for (indx = 0 ; indx < fl_right->nFeatures ; indx++)
			{
				sprintf(fname, "%sr_p%d_kltimg_tf_i%d.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_right->feature[indx]->patch->pyramid->img[i], fname);
				sprintf(fname, "%sr_p%d_kltimg_tf_i%d_gx.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_right->feature[indx]->patch->pyramid_gradx->img[i], fname);
				sprintf(fname, "%sr_p%d_kltimg_tf_i%d_gy.pgm", dir, indx, i);
				_KLTWriteFloatImageToPGM(fl_right->feature[indx]->patch->pyramid_grady->img[i], fname);
			}

		}
		else
		{
			sprintf(fname, "%sr_kltimg_tf_i%d.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_right->pyramid_last)->img[i],fname);
			sprintf(fname, "%sr_kltimg_tf_i%d_gx.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_right->pyramid_last_gradx)->img[i],fname);
			sprintf(fname, "%sr_kltimg_tf_i%d_gy.pgm", dir, i);
			_KLTWriteFloatImageToPGM(((_KLT_Pyramid)tc_right->pyramid_last_grady)->img[i],fname);
		}
	}



}





static void _horizontalSearch(	KLT_TrackingContext tc_left,
								KLT_TrackingContext tc_right,
								KLT_StereoTrackingContext stc,
								KLT_FeatureList fl_left,
								KLT_FeatureList fl_right,
								KLT_FeatureList fl_result)
{
	int j;
	int nFeaturesToTrack = fl_left->nFeatures;
	int nFeaturesRight = fl_right->nFeatures;
	unsigned int yOffsetTolerance = stc->yOffsetTolerance;
	int matches[nFeaturesRight];

	if(stc->rejectMultipleCorrespondences)
	{
		for(j=0;j<nFeaturesRight;j++)
			matches[j]=0;
	}

	KLT_SortList sl = _createSortList(nFeaturesRight);

	_featureList2sortList(fl_right, sl);

	  qsort((sl->element[0]), nFeaturesRight, sizeof(KLT_SortElementRec), _compareFeatureLoc);

	  /*for(j=0;j<nFeaturesRight;j++)
	  {
		  fprintf(stderr, "j: %d, index: %d, y: %f\n", j, sl->element[j]->index, sl->element[j]->y  );
	  }*/

	  for(j=0; j<nFeaturesToTrack; j++)
	  {
		  int lowerIndex, upperIndex;
		  int i=0;
		  _KLT_Patch patch_right_i;
		  float best_x=-1, best_y=-1;
		  int val, best_val, best_index_i;
		  float best_residue=stc->max_residue+1;
		  double imgdiff;
		  double best_imgdiff=-1;
		  //int best_relatedFeature=0;
		  float x_left=fl_left->feature[j]->x;
		  float x_left_fl=x_left; //used to recognize features which are not matchable (not reachable by trinagulation)
		  float y_left=fl_left->feature[j]->y;
		  float y_right_est=y_left+stc->yOffset;
		  _KLT_FloatImage img_left, gradx_left, grady_left, img_right, gradx_right, grady_right;

		  if(tc_left->patchMode)
		  {
			  _KLT_Patch patch_left_i=fl_left->feature[j]->patch;
			  x_left-=patch_left_i->roi->x;
			  y_left-=patch_left_i->roi->y;
			  img_left=patch_left_i->pyramid->img[0];
			  gradx_left=patch_left_i->pyramid_gradx->img[0];
			  grady_left=patch_left_i->pyramid_grady->img[0];
	 	  }
		  else
		  {
			  img_left=((_KLT_Pyramid)(tc_left->pyramid_last))->img[0];
			  gradx_left=((_KLT_Pyramid)(tc_left->pyramid_last_gradx))->img[0];
			  grady_left=((_KLT_Pyramid)(tc_left->pyramid_last_grady))->img[0];
		  }

		  while( (i<nFeaturesRight) && (sl->element[i]->y < (y_right_est-yOffsetTolerance)) )
			i++;
		  lowerIndex= (i==0) ? i : i--;

		  //search upper index
		  while( (i<nFeaturesRight) && (sl->element[i]->y < (y_right_est+yOffsetTolerance)) )
			i++;
		  upperIndex = --i;//(i==nFeatures) ? --i : i;


		  best_val = KLT_NO_STEREO_MATCH;
		  for(; i>=lowerIndex; i--)
		  {
			  int index_i=sl->element[i]->index;
			  float residue;

			  float x_right=fl_right->feature[index_i]->x;
			  float y_right=fl_right->feature[index_i]->y;

			  //check if feature is out of range for a possible match
			  if( (x_right > (x_left_fl-stc->min_pixelDisplacement)) || (x_right < (x_left_fl-stc->max_pixelDisplacement)) )
			  //if( x_right > (x_left_fl-yOffset) )
			  {
				  val = KLT_OUT_OF_STEREO_RANGE;
				  if (KLT_verbose >= 1)
					  fprintf(stderr,"j: %d, i: %d, out of stereo\n" ,j,index_i);
				  continue;
			  }

			  if(tc_right->patchMode)
			  {
				  patch_right_i=fl_right->feature[index_i]->patch;
				  x_right-=patch_right_i->roi->x;
				  y_right-=patch_right_i->roi->y;
				  img_right=patch_right_i->pyramid->img[0];
				  gradx_right=patch_right_i->pyramid_gradx->img[0];
				  grady_right=patch_right_i->pyramid_grady->img[0];
			  }
			  else
			  {
				  img_right=((_KLT_Pyramid)(tc_right->pyramid_last))->img[0];
				  gradx_right=((_KLT_Pyramid)(tc_right->pyramid_last_gradx))->img[0];
				  grady_right=((_KLT_Pyramid)(tc_right->pyramid_last_grady))->img[0];
			  }


			  /*  //variant where only the intensity differences are compared
			  imgdiff=0;
			  if (stc->lighting_insensitive) {
			    _computeIntensityDifferenceLightingInsensitiveValue(img_left, img_right, x_left, y_left,
			    													x_right, y_right,
			    												tc_left->window_width, tc_left->window_height,
			    												&imgdiff);
			  } else {
			    _computeIntensityDifferenceValue(img_left, img_right, x_left, y_left,
			    									x_right, y_right,
			    									tc_left->window_width, tc_left->window_height,
			    									&imgdiff);
			  }
			  */

			  val=_KLTTrackFeature( x_left, y_left, &x_right, &y_right,
					  img_left,gradx_left,grady_left,
					  img_right,gradx_right,gradx_right,
					  stc->window_width,stc->window_height,
					  stc->max_iterations,
					  stc->min_determinant,
					  stc->min_displacement,
					  stc->max_residue,
					  stc->lighting_insensitive,
					  &residue);

			  if (KLT_verbose >= 1)
				  fprintf(stderr,"j: %d, i: %d, val: %d\n" ,j,index_i,val);

			  //if( val==KLT_TRACKED && residue<best_residue)
			  //often the feature position is jumping between two locations, so that the iterations
			  //are never enough although the match would be quite good
			  if( (val==KLT_TRACKED || val==KLT_MAX_ITERATIONS) && residue<best_residue)
			  {
					/*if(imgdiff < 1200 && (best_imgdiff < 0 || imgdiff<best_imgdiff) )
					 {
					 best_imgdiff=imgdiff;
					 best_val=(int)imgdiff;*/
				  best_index_i=index_i;
				  best_val=val;
				  best_residue=residue;
				  if(tc_right->patchMode)
		  		  {
		  			  best_x=x_right + patch_right_i->roi->x;
		  			  best_y=y_right + patch_right_i->roi->y;
		  		  }
		  		  else
		  		  {
		  			  best_x=x_right;
		  			  best_y=y_right;
		  		  }
			  }
		  }

		  if(stc->rejectMultipleCorrespondences && best_val != KLT_NO_STEREO_MATCH)
		  {

			  if(matches[best_index_i] == 0) // correspondence not yet found
			  	  matches[best_index_i] = j;
			  else if(matches[best_index_i] == -1) // correspondence many times found
				  best_val = KLT_MULTIPLE_STEREO_MATCH;
			  else // correspondence second time found
			  {
				  fl_result->feature[matches[best_index_i]]->val=KLT_MULTIPLE_STEREO_MATCH;
				  matches[best_index_i] = -1;

			  }
		  }

		  fl_result->feature[j]->val=best_val;
		  fl_result->feature[j]->x=best_x;
		  fl_result->feature[j]->y=best_y;
	  }

	  _freeSortList(sl);

}


static void _epipolarSearch(
		KLT_TrackingContext tc_left,
		KLT_TrackingContext tc_right,
		KLT_StereoTrackingContext stc,
		KLT_FeatureList fl_left,
		KLT_FeatureList fl_right,
		KLT_FeatureList fl_result) {
	int j, i;
	int nFeaturesToTrack = fl_left->nFeatures;
	int nFeaturesRight = fl_right->nFeatures;
	unsigned int yOffsetTolerance = stc->yOffsetTolerance;
	int matches[nFeaturesRight];

	if(stc->rejectMultipleCorrespondences)
	{
		for(j=0;j<nFeaturesRight;j++)
			matches[j]=0;
	}


	for (j=0; j<nFeaturesToTrack; j++) {
		_KLT_Patch patch_right_i;
		float x_right=0, y_right=0, best_x=-1, best_y=-1;
		int val, best_val, best_i;
		float best_residue=stc->max_residue+1;
		double imgdiff;
		double best_imgdiff=-1;
		//int best_relatedFeature=0;
		float x_left=fl_left->feature[j]->x;
		float x_left_fl=x_left;
		float y_left=fl_left->feature[j]->y;
		float y_left_fl=y_left;
		_KLT_FloatImage img_left, gradx_left, grady_left, img_right,
				gradx_right, grady_right;
		float y_right_est=y_left_fl+stc->yOffset;

		if (tc_left->patchMode) {
			_KLT_Patch patch_left_i=fl_left->feature[j]->patch;
			x_left-=patch_left_i->roi->x;
			y_left-=patch_left_i->roi->y;
			img_left=patch_left_i->pyramid->img[0];
			gradx_left=patch_left_i->pyramid_gradx->img[0];
			grady_left=patch_left_i->pyramid_grady->img[0];
		} else {
			img_left=((_KLT_Pyramid)(tc_left->pyramid_last))->img[0];
			gradx_left=((_KLT_Pyramid)(tc_left->pyramid_last_gradx))->img[0];
			grady_left=((_KLT_Pyramid)(tc_left->pyramid_last_grady))->img[0];
		}

		best_val = KLT_NO_STEREO_MATCH;

		for (i=0; i<nFeaturesRight; i++)
		{
			float residue;

			x_right=fl_right->feature[i]->x;
			y_right=fl_right->feature[i]->y;

			//fprintf(stderr,"j: %d, i: %d --> %f\n" ,j,i, stc->get_distanceFromEpipolarLine(x_left_fl, y_left_fl, x_right, y_right));

			//check if feature is out of range for a possible match
			if (y_right < (y_right_est-yOffsetTolerance) ||
				y_right > (y_right_est+yOffsetTolerance) ||
				x_right > (x_left_fl-stc->min_pixelDisplacement) ||
				x_right < (x_left_fl-stc->max_pixelDisplacement) ||
				stc->get_distanceFromEpipolarLine(x_left_fl, y_left_fl,x_right, y_right) > stc->max_distanceFromEpipolarLine)
			{
				val = KLT_OUT_OF_STEREO_RANGE;
				if (KLT_verbose >= 1)
					fprintf(stderr,"j: %d, i: %d, out of stereo\n" ,j,i);
				continue;
			}

			//fprintf(stderr,"yOffsetTolerance: %d, y_right_est: %f, min: %d, max: %d\n" ,yOffsetTolerance,y_right_est, stc->min_pixelDisplacement, stc->min_pixelDisplacement);
			//fprintf(stderr,"left_x: %f, left_y: %f, right_x: %f, right_y: %f\n" ,x_left_fl,y_left_fl, x_right, y_right);

			if (tc_right->patchMode) {
				patch_right_i=fl_right->feature[i]->patch;
				x_right-=patch_right_i->roi->x;
				y_right-=patch_right_i->roi->y;
				img_right=patch_right_i->pyramid->img[0];
				gradx_right=patch_right_i->pyramid_gradx->img[0];
				grady_right=patch_right_i->pyramid_grady->img[0];
			} else {
				img_right=((_KLT_Pyramid)(tc_right->pyramid_last))->img[0];
				gradx_right=((_KLT_Pyramid)(tc_right->pyramid_last_gradx))->img[0];
				grady_right=((_KLT_Pyramid)(tc_right->pyramid_last_grady))->img[0];
			}

			/*  //variant where only the intensity differences are compared
			 imgdiff=0;
			 if (stc->lighting_insensitive) {
			 _computeIntensityDifferenceLightingInsensitiveValue(img_left, img_right, x_left, y_left,
			 x_right, y_right,
			 tc_left->window_width, tc_left->window_height,
			 &imgdiff);
			 } else {
			 _computeIntensityDifferenceValue(img_left, img_right, x_left, y_left,
			 x_right, y_right,
			 tc_left->window_width, tc_left->window_height,
			 &imgdiff);
			 }
			 */

			val=_KLTTrackFeature(x_left, y_left, &x_right, &y_right, img_left,
					gradx_left, grady_left, img_right, gradx_right,
					gradx_right, stc->window_width, stc->window_height,
					stc->max_iterations, stc->min_determinant,
					stc->min_displacement, stc->max_residue,
					stc->lighting_insensitive, &residue);

			if (KLT_verbose >= 1)
				fprintf(stderr,"j: %d, i: %d, val: %d\n" ,j,i,val);

			//if( val==KLT_TRACKED && residue<best_residue)
			//often the feature position is jumping between two locations, so that the iterations
			//are never enough although the match would be quite good
			if ( (val==KLT_TRACKED || val==KLT_MAX_ITERATIONS) && residue<best_residue) {
				/*if(imgdiff < 1200 && (best_imgdiff < 0 || imgdiff<best_imgdiff) )
				 {
				 best_imgdiff=imgdiff;
				 best_val=(int)imgdiff;*/
				best_i=i;
				best_val=val;
				best_residue=residue;
				if (tc_right->patchMode) {
					best_x=x_right + patch_right_i->roi->x;
					best_y=y_right + patch_right_i->roi->y;
				} else {
					best_x=x_right;
					best_y=y_right;
				}
			}
		}

		if(stc->rejectMultipleCorrespondences && best_val != KLT_NO_STEREO_MATCH)
		{
			if (matches[best_i] == 0) // correspondence not yet found
				matches[best_i] = j;
			else if (matches[best_i] == -1) // correspondence many times found
				best_val = KLT_MULTIPLE_STEREO_MATCH;
			else // correspondence seconde time found
			{
				fl_result->feature[matches[best_i]]->val=KLT_MULTIPLE_STEREO_MATCH;
				matches[best_i] = -1;
			}
		}

		fl_result->feature[j]->val=best_val;
		fl_result->feature[j]->x=best_x;
		fl_result->feature[j]->y=best_y;

	}

}






/*********************************************************************
 * _matchFeatures
 *
 * matches the features in a stereo image pair
 *
 */

static void _matchFeatures(
		KLT_PixelType *img_left,
		KLT_PixelType *img_right,
        int ncols,
		int nrows,
		KLT_TrackingContext tc_left,
		KLT_TrackingContext tc_right,
		KLT_StereoTrackingContext stc,
		KLT_FeatureList fl_left,
		KLT_FeatureList fl_right,
		KLT_FeatureList fl_result)
{
  int roiSrc_width  = tc_left->window_width+2*tc_left->borderx;
  int roiSrc_height = tc_left->window_height+2*tc_left->bordery;
  int roiDest_width  = roiSrc_width  + 2*stc->max_patchDisplacement;
  int roiDest_height = roiSrc_height + 2*stc->max_patchDisplacement;


  //_KLT_FloatImage helpImg_left = _KLTCreateFloatImage(ncols, nrows);
  //_KLT_FloatImage helpImg_right = _KLTCreateFloatImage(ncols, nrows);
  //_KLTToFloatImage(img_left, ncols, nrows, helpImg_left);
  //_KLTToFloatImage(img_right, ncols, nrows, helpImg_right);
  //_KLTFreeFloatImage(helpImg_left);
  //_KLTFreeFloatImage(helpImg_right);



  if(tc_left->patchMode)
	  _computePatchPyramid(img_left, fl_left, ncols, nrows, roiSrc_width, roiSrc_height ,tc_left);
  else
	  _computePyramid(img_left, ncols, nrows, tc_left);

  if(tc_right->patchMode)
	  _computePatchPyramid(img_right, fl_right, ncols, nrows, roiDest_width, roiDest_height ,tc_right);
  else
	  _computePyramid(img_right, ncols, nrows, tc_right);


  if (stc->writeInternalImages)
	  _writeInternalImages(tc_left, tc_right, stc, fl_left, fl_right);




  if(stc->get_distanceFromEpipolarLine == NULL)
	  _horizontalSearch(tc_left, tc_right, stc, fl_left, fl_right, fl_result);
  else
	  _epipolarSearch(tc_left, tc_right, stc, fl_left, fl_right, fl_result);














}




/*********************************************************************
 * KLTExtFindStereoMatches
 *
 * matches the features in a stereo image pair
 * therefor following constrainst have to be taken into account:
 * - this function works only with a coplanar stereo camera pair with
 *   cameras mounted on the same high
 * - the left image has to be the left image
 * - it is not possible to use for the left and the right image the
 *   same tracking context, even if the parameters are the same
 *
 * The resulting feature list contains the coordinates of the features
 * in the right image with the indices corresponding to the reature
 * list of the left image.
 * Experiments have shown, that even if the maximum iteration number
 * exceeds, the coordinates may be sense- and therfor useful.
 *
 */

void KLTExtFindStereoMatches(
	  KLT_TrackingContext tc_left,
	  KLT_PixelType *img_left,
	  KLT_FeatureList fl_left,
	  IppiRect *roi_left,
	  KLT_TrackingContext tc_right,
	  KLT_PixelType *img_right,
	  KLT_FeatureList fl_right,
	  IppiRect *roi_right,
	  int ncols,
	  int nrows,
	  KLT_StereoTrackingContext tc_stereo,
	  KLT_FeatureList fl_result,
	  int width_dispersionFactor,
  	  int height_dispersionFactor)
{
	KLT_BOOL newRoi_left=FALSE;
	KLT_BOOL newRoi_right=FALSE;
	unsigned int offset=tc_stereo->xOffset;
	int oldnPL_left = tc_left->nPyramidLevels;
	int oldnPL_right = tc_right->nPyramidLevels;
	tc_left->nPyramidLevels = tc_stereo->nPyramidLevels;
	tc_right->nPyramidLevels = tc_stereo->nPyramidLevels;

	//tc_left->writeInternalImages=TRUE;
	//tc_right->writeInternalImages=TRUE;

	KLTUpdateTCBorder(tc_left);
	KLTUpdateTCBorder(tc_right);





	//KLTPrintTrackingContext(tc_left);

	//muss nicht mehr sein
	//assert(fl_left->nFeatures == fl_right->nFeatures);

	if(roi_left==NULL)
	{
		newRoi_left = TRUE;
		roi_left = (IppiRect*)malloc(sizeof(IppiRect));
		//set up left ROI
		roi_left->x=offset;
		roi_left->y=0;
		roi_left->width=ncols-offset;
		roi_left->height=nrows;
	}

	if(roi_right==NULL)
	{
		newRoi_right = TRUE;
		roi_right = (IppiRect*)malloc(sizeof(IppiRect));
		//set up right ROI
		roi_right->x=0;
		roi_right->y=0;
		roi_right->width=ncols-offset;
		roi_right->height=nrows;
	}

	KLTExtSelectGoodFeatures(tc_left, img_left, ncols, nrows, fl_left, roi_left, width_dispersionFactor, height_dispersionFactor);
	KLTExtSelectGoodFeatures(tc_right, img_right, ncols, nrows, fl_right, roi_right, width_dispersionFactor, height_dispersionFactor);

	if(newRoi_left)
		free(roi_left), roi_left=NULL;
	if(newRoi_right)
		free(roi_right), roi_right=NULL;

	_matchFeatures(img_left, img_right, ncols, nrows, tc_left, tc_right, tc_stereo, fl_left, fl_right, fl_result);

	//set back the old nPyramidLevel values
	if(tc_left->nPyramidLevels != oldnPL_left || !tc_left->sequentialMode)
	{
		KLTFreeStoredImages(tc_left);
		tc_left->nPyramidLevels = oldnPL_left;
		KLTUpdateTCBorder(tc_left);
		if(tc_left->patchMode)
		{
			int roiSrc_width  = tc_left->window_width+2*tc_left->borderx;
  			int roiSrc_height = tc_left->window_height+2*tc_left->bordery;
  			//TODO: implement it more efficiently that it doesn't have to delete the first level images, which have to be redone here afterwards again
  			//or at least do it so, that it is only necessary for the valid correspondences
	  		_computePatchPyramid(img_left, fl_left, ncols, nrows, roiSrc_width, roiSrc_height ,tc_left);
		}
	}

	if(tc_right->nPyramidLevels != oldnPL_right || !tc_right->sequentialMode)
	{
		KLTFreeStoredImages(tc_right);
		tc_right->nPyramidLevels = oldnPL_right;
		KLTUpdateTCBorder(tc_right);
		if(tc_right->patchMode)
		{
			int roiSrc_width  = tc_right->window_width+2*tc_left->borderx;
  			int roiSrc_height = tc_right->window_height+2*tc_left->bordery;
	  		_computePatchPyramid(img_right, fl_right, ncols, nrows, roiSrc_width, roiSrc_height ,tc_right);
		}
	}


		/*KLTWriteFeatureListToColoredPPM(fl_left, img_left, ncols, nrows, "feat_left.ppm");
	KLTWriteFeatureListToColoredPPM(fl_right, img_right, ncols, nrows, "feat_right.ppm");
	KLTWriteFeatureListToColoredPPM(fl_result, img_right, ncols, nrows, "feat_result.ppm");*/


}


