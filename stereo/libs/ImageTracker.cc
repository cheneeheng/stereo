#include <iostream>
#include "ippi.h"
#include "ippcv.h"
#include "XVWindowX.h"
#include "ImageTracker.h"



XVImageRGB<XV_RGB24>  &
ImageTracker::frame(void)
{
    static XVImageRGB<XV_RGB> rpc_image;

    return dev->frame(left_index);
}

void 
ImageTracker::step_track(XVImageRGB<XV_RGB24> &image)
{

    tmp_image.resize(image);
    gray_image.resize(image);
    RGBtoScalar(image,tmp_image);
    IppiSize roi={tmp_image.Width(), tmp_image.Height()};
    ippiUndistortRadial_8u_C1R(tmp_image.data(),tmp_image.Width(),
                               gray_image.lock(),gray_image.Width(),
			       roi,cam_param.f[0], cam_param.f[1],
			       cam_param.C[0],cam_param.C[1],
			       cam_param.kappa[0],cam_param.kappa[1],
			       DistortBuffer);

    gray_image.unlock();
    KLTTrackFeatures(tc, (KLT_PixelType*)gray_image.data(), 
                         (KLT_PixelType*)gray_image.data(),
			 gray_image.Width(), gray_image.Height(), fl);
}

void 
ImageTracker::next_frame()
{
  if(dev)
  {
    // switch to next buffer
     // schedule next frame
    left_index^=1;
    dev->initiate_acquire(left_index^1);
    // wait for previous one
    dev->wait_for_completion(left_index);
  }
}

void 
ImageTracker::step(XVDrawWindowX<XV_RGB>*left_win, bool track)
{
  XVImageRGB<XV_RGB> rpc_image; 

if(track)
{
  // put images from framegrabber on the screen
  if(left_win)
   left_win->CopySubImage(dev->frame(left_index));
  step_track(dev->frame(left_index));
 if(left_win)
  for(int k=0;k<num_points;k++)
     if(fl->feature[k]->val>=0)
         left_win->drawEllipse(fl->feature[k]->x-5,
                           fl->feature[k]->y-5, 11,11,"red");
 }

next_frame();
}

void ImageTracker::replace(void)
{

  KLTReplaceLostFeatures(tc, (KLT_PixelType*)gray_image.data(), 
                             gray_image.Width(), gray_image.Height(), fl);

}

ImageTracker::ImageTracker(const Config  &config,int twhich)
{
 dev=NULL,rpc_handle=NULL;pose_valid=false;which=twhich;
 cam_param=config.camera_params[which];
 XVImageRGB<XV_RGB> rpc_image; 
 // open input devices
  switch(config.type)
  {
     case FILE_INTERFACE:
       dev=new XVAVI<XVImageRGB<XV_RGB24> >(config.file[which].mask);
       break;
     default:
       cerr << "unknown streaming type" << endl;
       throw 10;
  }
   num_points=config.num_points;
   if(!dev && !rpc_handle)
   {
     cerr << "Couldn't open input source" << endl;
     throw 11;
   }
  if(dev)
  {
     left_index=0; initialized=false;
     dev->initiate_acquire(left_index);
     dev->initiate_acquire(left_index^1);
     dev->wait_for_completion(left_index);
  }
  
   int tr_count=0;
   int width=SCALE_FACTOR*dev->frame(left_index).Width();
   int height= SCALE_FACTOR*dev->frame(left_index).Height();
   if(which) return;
   tc = KLTCreateTrackingContext();
   fl = KLTCreateFeatureList(num_points);
   tc->sequentialMode = TRUE;
   tc->writeInternalImages = FALSE;
   tc->affineConsistencyCheck = -1;
   tc->sequentialMode = true;
   tc->bindToCorner=true;
   tc->patchMode=true;
   tc->mindist=25;
   tc->writeInternalImages = FALSE;
   tc->window_width=7;
   tc->window_height=7;
   tc->affineConsistencyCheck = -1; 
   KLTChangeTCPyramid(tc, 10);
   KLTUpdateTCBorder(tc);

#ifdef CHECK
   XVDrawWindowX<XV_RGB> window(320,240);
   window.map();
   window.CopySubImage(dev->frame(left_index));
#endif
   RGBtoScalar(dev->frame(left_index),gray_image);
   XVImageScalar<u_char> undist_image(dev->frame(left_index));
   {
     IppiSize roi={dev->frame(left_index).Width(),
                   dev->frame(left_index).Height()};
     int dist_buf_size;
     ippiUndistortGetSize(roi,&dist_buf_size);
     DistortBuffer=new Ipp8u[dist_buf_size];
   }
   {
     IppiSize roi={gray_image.Width(),gray_image.Height()};
     ippiUndistortRadial_8u_C1R((Ipp8u*)gray_image.data(), gray_image.Width(),
                               undist_image.lock(),undist_image.Width(),
			       roi,cam_param.f[0], cam_param.f[1],
			       cam_param.C[0],cam_param.C[1],
			       cam_param.kappa[0],cam_param.kappa[1],
			       DistortBuffer);

     undist_image.unlock();
   }
   {
     IppiRect roi={296,160,300,300};
     KLTExtSelectGoodFeatures(tc, (KLT_PixelType*)undist_image.data(), 
                              width, height, fl,&roi,2,2);
   }
#ifdef CHECK
   for(int j=0;j<num_points;j++)
   {   
     window.fillEllipse(fl->feature[j]->x-5,fl->feature[j]->y-5,
                      10,10,"red");
     tr_count++;
   }   
   window.swap_buffers();
   window.flush();
   sleep(20);
#endif
   dev->initiate_acquire(left_index^1);

   KLTSetVerbosity(0);

}
