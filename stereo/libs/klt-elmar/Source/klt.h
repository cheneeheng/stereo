/*********************************************************************
 * klt.h
 *
 * Kanade-Lucas-Tomasi tracker
 *********************************************************************/

#ifndef _KLT_H_
#define _KLT_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef float KLT_locType;
typedef unsigned char KLT_PixelType;

#define KLT_BOOL int

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#ifndef NULL
#define NULL  0
#endif

#define KLT_REFIND          		 1
#define KLT_TRACKED          		 0
#define KLT_NOT_FOUND        		-1
#define KLT_SMALL_DET        		-2
#define KLT_MAX_ITERATIONS   		-3
#define KLT_OOB              		-4
#define KLT_LARGE_RESIDUE    		-5
#define KLT_NO_STEREO_MATCH	 		-6
#define KLT_OUT_OF_STEREO_RANGE	 	-7
#define KLT_MULTIPLE_STEREO_MATCH	-8

#define KLTEXT_STORE_RESULT_IN_SOURCE_LIST		1
#define KLTEXT_STORE_RESULT_IN_DESTINATION_LIST	0

#include <ippi.h>
#include "klt_util.h" /* for affine mapping */
#include "stereo_util.h" /* for stereo tracking support */
#include "pyramid.h"

/*******************
 * Structures
 */

typedef struct  {
  /* Available to user */
  int mindist;			/* min distance b/w features */
  int window_width, window_height;
  int sequentialMode;	/* whether to save most recent image to save time - if set to 2 it means that a linearPrediction of the feature positions is applied */
  /* can set to TRUE manually, but don't set to */
  /* FALSE manually */
  KLT_BOOL smoothBeforeSelecting;	/* whether to smooth image before */
  KLT_BOOL bindToCorner;	/* whether to bind to the next corner in each tracking step */    
  /* selecting features */
  int writeInternalImages;	/* whether to write internal images - if set to 1 overwrite old image, if set to 2-9 writes every image as sequence, if set > 10 writes only the image of feature (writeInternalImages-10) as sequence */
  /* tracking features */
  KLT_BOOL lighting_insensitive;  /* whether to normalize for gain and bias (not in original algorithm) */
  
  /* Available, but hopefully can ignore */
  int min_eigenvalue;		/* smallest eigenvalue allowed for selecting */
  float min_determinant;	/* th for determining lost */
  float min_displacement;	/* th for stopping tracking when pixel changes little */
  int max_iterations;		/* th for stopping tracking when too many iterations */
  float max_residue;		/* th for stopping tracking when residue is large */
  float grad_sigma;
  float smooth_sigma_fact;
  float pyramid_sigma_fact;
  int nSkippedPixels;		/* # of pixels skipped when finding features */
  int borderx;			/* border in which features will not be found */
  int bordery;
  int nPyramidLevels;		/* computed from search_ranges */
  int subsampling;		/* 		" */

  /*for extended KLT*/
  int max_patchDisplacement;	/* amount of pixels to enlarge searching patch on each side */

  
  /* for affine mapping */ 
  int affine_window_width, affine_window_height;
  int affineConsistencyCheck; /* whether to evaluates the consistency of features with affine mapping 
				 -1 = don't evaluates the consistency
				 0 = evaluates the consistency of features with translation mapping
				 1 = evaluates the consistency of features with similarity mapping
				 2 = evaluates the consistency of features with affine mapping
			      */
  int affine_max_iterations;  
  float affine_max_residue;
  float affine_min_displacement;        
  float affine_max_displacement_differ; /* th for the difference between the displacement calculated 
					   by the affine tracker and the frame to frame tracker in pel*/

  /* User must not touch these */
  KLT_BOOL patchMode;
  _KLT_FloatPixelDispRec *pixelDisplacements_last;
  void *pyramid_last;
  void *pyramid_last_gradx;
  void *pyramid_last_grady;
}  KLT_TrackingContextRec, *KLT_TrackingContext;


typedef struct  {
  KLT_locType x;
  KLT_locType y;
  int val;	
  _KLT_Patch patch;

  /* for affine mapping */
  _KLT_FloatImage aff_img; 
  _KLT_FloatImage aff_img_gradx;
  _KLT_FloatImage aff_img_grady;
  KLT_locType aff_x;
  KLT_locType aff_y;
  KLT_locType aff_Axx;
  KLT_locType aff_Ayx;
  KLT_locType aff_Axy;
  KLT_locType aff_Ayy;
}  KLT_FeatureRec, *KLT_Feature;

typedef struct  {
  int nFeatures;
  KLT_Feature *feature;
}  KLT_FeatureListRec, *KLT_FeatureList;

typedef struct  {
  int nFrames;
  KLT_Feature *feature;
}  KLT_FeatureHistoryRec, *KLT_FeatureHistory;

typedef struct  {
  int nFrames;
  int nFeatures;
  KLT_Feature **feature;
}  KLT_FeatureTableRec, *KLT_FeatureTable;



typedef struct  {
  /* Available to user */
  int window_width, window_height;
  /* can set to TRUE manually, but don't set to */
  /* tracking features */
  KLT_BOOL lighting_insensitive;  /* whether to normalize for gain and bias (not in original algorithm) */
  
  /* Available, but hopefully can ignore */
  float min_determinant;	/* th for determining lost */
  float min_displacement;	/* th for stopping tracking when pixel changes little */
  int max_iterations;		/* th for stopping tracking when too many iterations */
  float max_residue;		/* th for stopping tracking when residue is large */

  unsigned int max_patchDisplacement; /* amount of pixels to enlarge searching patch on each side */
  
  int nPyramidLevels;		
  
  int yOffset;			/*offset between rows, which features can belong to*/
  unsigned int yOffsetTolerance;			/*offset between rows, which features can belong to*/
  unsigned int xOffset;			/*offset at the left resp. right side of the image's borders, only used, if tc->left and/or tc->right == 0*/
 
  int min_pixelDisplacement; //minimum displacement in pixels
  int max_pixelDisplacement; //maximum displacement in pixels
    
  int writeInternalImages;
  
  KLT_BOOL rejectMultipleCorrespondences;	// whether to reject multiple correspondences or not     

  float (*get_distanceFromEpipolarLine)(float x_left, float y_left, float x_right, float y_right); //function pointer to the epipolar distance calculation routine
  float max_distanceFromEpipolarLine; //maximum distance from epipolar line 
}  KLT_StereoTrackingContextRec, *KLT_StereoTrackingContext;




/*******************
 * Functions
 */

/* Create */
KLT_TrackingContext KLTCreateTrackingContext(void);
KLT_StereoTrackingContext KLTCreateStereoTrackingContext(void);
KLT_FeatureList KLTCreateFeatureList(
  int nFeatures);
KLT_FeatureHistory KLTCreateFeatureHistory(
  int nFrames);
KLT_FeatureTable KLTCreateFeatureTable(
  int nFrames,
  int nFeatures);

/* Free */
void KLTFreeTrackingContext(
  KLT_TrackingContext tc);
void KLTFreeStereoTrackingContext(
  KLT_StereoTrackingContext stc);
void KLTFreeFeatureList(
  KLT_FeatureList fl);
void KLTFreeFeatureHistory(
  KLT_FeatureHistory fh);
void KLTFreeFeatureTable(
  KLT_FeatureTable ft);

/* Processing */
void KLTSelectGoodFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img,
  int ncols,
  int nrows,
  KLT_FeatureList fl);
void KLTExtSelectGoodFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img, 
  int ncols, 
  int nrows,
  KLT_FeatureList fl,
  IppiRect *roi,
  int width_dispersionFactor,
  int height_dispersionFactor);
void KLTExtSelectGoodFeaturesInSameImage(
  KLT_TrackingContext tc,
  KLT_PixelType *img, 
  int ncols, 
  int nrows,
  KLT_FeatureList fl,
  IppiRect *roi);
void KLTTrackFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img1,
  KLT_PixelType *img2,
  int ncols,
  int nrows,
  KLT_FeatureList fl);
void KLTExtTrackFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img1,
  KLT_PixelType *img2,
  int ncols,
  int nrows,
  KLT_FeatureList flSrc,
  KLT_FeatureList flDest);
void KLTReplaceLostFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img,
  int ncols,
  int nrows,
  KLT_FeatureList fl);
void KLTExtReplaceLostFeatures(
  KLT_TrackingContext tc,
  KLT_PixelType *img, 
  int ncols, 
  int nrows,
  KLT_FeatureList fl,
  IppiRect *roi);

/* Utilities */
int KLTCountRemainingFeatures(
  KLT_FeatureList fl);
int KLTExtCountRemainingFeatures(
  KLT_FeatureList fl);
void KLTPrintTrackingContext(
  KLT_TrackingContext tc);
void KLTChangeTCPyramid(
  KLT_TrackingContext tc,
  int search_range);
void KLTUpdateTCBorder(
  KLT_TrackingContext tc);
void KLTStopSequentialMode(
  KLT_TrackingContext tc);
void KLTSetVerbosity(
  int verbosity);
float _KLTComputeSmoothSigma(
  KLT_TrackingContext tc);
int KLTGetRemainingFeatures(
  KLT_FeatureList fl, 
  KLT_FeatureList fl_valid);
void KLTFreeStoredImages(
  KLT_TrackingContext tc);

/* Storing/Extracting Features */
void KLTStoreFeatureList(
  KLT_FeatureList fl,
  KLT_FeatureTable ft,
  int frame);
void KLTExtractFeatureList(
  KLT_FeatureList fl,
  KLT_FeatureTable ft,
  int frame);
void KLTStoreFeatureHistory(
  KLT_FeatureHistory fh,
  KLT_FeatureTable ft,
  int feat);
void KLTExtractFeatureHistory(
  KLT_FeatureHistory fh,
  KLT_FeatureTable ft,
  int feat);

/* Writing/Reading */
void KLTWriteFeatureListToPPM(
  KLT_FeatureList fl,
  KLT_PixelType *greyimg,
  int ncols,
  int nrows,
  char *filename);
void KLTWriteStereoResultToPPM(
  KLT_FeatureList fl_left, 
  KLT_PixelType *img_left, 
  KLT_FeatureList fl_right, 
  KLT_PixelType *img_right, 
  KLT_FeatureList fl_result, 
  int ncols, 
  int nrows, 
  char *filename);
void KLTWriteFeatureListToColoredPPM(
  KLT_FeatureList featurelist,
  KLT_PixelType *greyimg,
  int ncols,
  int nrows,
  char *filename);
void KLTWriteFeatureList(
  KLT_FeatureList fl,
  char *filename,
  char *fmt);
void KLTWriteFeatureHistory(
  KLT_FeatureHistory fh,
  char *filename,
  char *fmt);
void KLTWriteFeatureTable(
  KLT_FeatureTable ft,
  char *filename,
  char *fmt);
KLT_FeatureList KLTReadFeatureList(
  KLT_FeatureList fl,
  char *filename);
KLT_FeatureHistory KLTReadFeatureHistory(
  KLT_FeatureHistory fh,
  char *filename);
KLT_FeatureTable KLTReadFeatureTable(
  KLT_FeatureTable ft,
  char *filename);

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
  int height_dispersionFactor);


void KLTMoveFeature(
  KLT_Feature featureSrc, 
  KLT_Feature featureDest);

void KLTMoveFeatureList(
  KLT_FeatureList featurelistSrc, 
  KLT_FeatureList featurelistDest);

int KLTCopyFeatureList(
  KLT_FeatureList featurelistSrc, 
  KLT_FeatureList featurelistDest);

/*inline void KLTExtCopyFeatureList(
  KLT_FeatureList featurelistSrc, 
  KLT_FeatureList featurelistDest,
  KLT_TrackingContext tc,
  int ncols, 
  int nrows);*/

void KLTExtPropagateFeatureList(
  KLT_FeatureList featurelistSrc, 
  KLT_FeatureList featurelistDest,
  KLT_TrackingContext tc,
  int ncols, 
  int nrows);

void KLTCopyTrackingContext(
  KLT_TrackingContext tc_from, 
  KLT_TrackingContext tc_to);



#ifdef __cplusplus
}
#endif

#endif






