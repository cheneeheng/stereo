#include <exception>
#include "config.h"
#include <XVImageIO.h>
#include <XVImageSeq.h>
#ifndef DLR
#include <XVDig1394.h>
#endif
//#include "XVAVI.h"
#include <XVWindowX.h>
#include <XVMatrix.h>
#include "klt.h"
#include "nav_config.h"

#ifdef RPC_EXTENSION
#include "ScannerApplication.h"
#else
  class ScannerApplication{ //dummy class
   public:
    ScannerApplication(bool flag,float interv,u_int Port){ };
    ~ScannerApplication(){};
    XVImageRGB<XV_RGB> get_image(int which){
               static XVImageRGB<XV_RGB> im(0,0);
             return im;}
  };
#endif
   
typedef enum {TR_UNKNOWN,TR_LOCKED,TR_EXPLORING} TrackingState;

typedef struct
{
  TrackingState st; 
  XVColVector   T;  
  XVMatrix      R;  
  XVColVector	D;
  KLT_FeatureRec state;
}NewRecord;


class ImageTracker:public exception
{
   private:
      XVVideo<XVImageRGB<XV_RGB24> >
      				*dev;
      ScannerApplication        *rpc_handle;
      int			active_index; //for template
      int			num_points;
      int			which; //which camera?
      XVImageScalar<unsigned char> gray_image;
      XVImageScalar<unsigned char> tmp_image;
      KLT_TrackingContext 	tc;
      KLT_FeatureList           fl;
      KLT_FeatureTable 		ft;
      CameraParams		cam_param;
      Ipp8u			*DistortBuffer;
      short int			left_index;
      bool			initialized;
      bool			pose_valid;
      XVMatrix			current_robot_pose;
   public:
      				ImageTracker(const Config &config,int which);
      void			step(XVDrawWindowX<XV_RGB>*left_win,bool track);
      void			replace();				      		
#ifdef DLR
      bool			get_robot_pose(XVMatrix &pose)
      					{pose=current_robot_pose;
					 return pose_valid;};
#endif
      XVImageRGB<XV_RGB24> &      frame();
      void			next_frame();
      KLT_Feature		state(int which) { return fl->feature[which];};
#ifdef DLR
      void			set_pose(XVMatrix &pose)
      			        { if(rpc_handle)
			          rpc_handle->report_modeller_pose(pose);};
#endif
      int			num_active() {return
                                  KLTCountRemainingFeatures(fl);};
      void			step_track(XVImageRGB<XV_RGB24> &image);
      virtual			~ImageTracker() throw()
				   {if(rpc_handle) delete rpc_handle;
				    if(dev) delete dev;};
};
