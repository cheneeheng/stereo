#ifndef STEREO_UTIL_H_
#define STEREO_UTIL_H_

//#include "klt.h"



typedef struct  {
  double y;
  int index;
}  KLT_SortElementRec, *KLT_SortElement;


typedef struct  {
  KLT_SortElement *element;
  int nFeatures;
}  KLT_SortListRec, *KLT_SortList;




#endif /*STEREO_UTIL_H_*/
