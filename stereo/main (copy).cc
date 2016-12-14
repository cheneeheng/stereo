#define HIGH_SPEED_ID 302662258
#define LOW_SPEED_ID  235376646

#ifndef Sqr
#define Sqr(a) ((a)*(a))
#endif

//#define FREQ_CAM
//#define FREQ_STEREO
//#define FREQ_TRACKER

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/calib3d.hpp>

//#include <opencv2/core/utility.hpp>
//#include <opencv2/core/cuda.hpp>
//#include <opencv2/cudastereo.hpp>

/*
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
*/

#include <dc1394/control.h>

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <semaphore.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include <cuda.h>
#include <Stereo.h>

#include "klt.h"
#include "Vgps.h"

#include <XVImageSeq.h>
#include <XVMpeg.h>
#include <XVImageIO.h>
#include <XVColorSeg.h>
#include <XVBlobFeature.h>
#include <XVTracker.h>
#include <XVWindowX.h>
#include "ippi.h"
#include "ippcv.h"
#include "ippcc.h"

// threads
sem_t mutex_t1,mutex_t2,mutex_t3;
sem_t lock_t1,lock_t2,lock_t3;

// images
cv::Mat img_left;
cv::Mat img_right;
std::vector<cv::Mat> img_right_container(4);
cv::Mat point3D;
cv::Mat disp;

const int  	      num_disparities=128; //def = 128
static int16_t        mGamma[num_disparities*16];

// pose estimation
typedef struct{
     XVMatrix           R;
     XVColVector	T;
}TransfMatrix;


//==============================================================================================================================================================
/*
pcl::PointCloud<pcl::PointXYZRGB>::Ptr MatToPoinXYZ(cv::Mat OpenCVPointCloud, cv::Mat OpenCVImage)
         {

             //char pr=100, pg=100, pb=100;
             pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);//(new pcl::pointcloud<pcl::pointXYZ>);

             cv::Vec3f single_point_obj; cv::Vec3b single_point_obj_rgb; float tmp = 5.0;

             for(int i=0;i<OpenCVPointCloud.size().height;i++)
             {
             for(int ii=0;ii<OpenCVPointCloud.size().width;ii++)
             {
               single_point_obj = OpenCVPointCloud.at<cv::Vec3f>(i,ii);
               single_point_obj_rgb = OpenCVImage.at<cv::Vec3b>(i,ii); 

                pcl::PointXYZRGB point;
                point.x = single_point_obj[0];
                point.y = -single_point_obj[1];
                point.z = single_point_obj[2] < 100 ? single_point_obj[2] : 0 ;
                point.z = point.z > 0 ? point.z : 0 ;   
                point.b = point.z != 0 ? single_point_obj_rgb[0]: 0 ;
                point.g = point.z != 0 ? single_point_obj_rgb[1]: 0 ;
                point.r = point.z != 0 ? single_point_obj_rgb[2]: 0 ;

                // when color needs to be added:
                //uint8_t pr = 0, pg = 255, pb = 0; // Example: Green color
                //uint32_t rgb = (static_cast<uint32_t>(pr) << 16 | static_cast<uint32_t>(pg) << 8 | static_cast<uint32_t>(pb));
                //point.rgb = *reinterpret_cast<float*>(&rgb);

                point_cloud_ptr -> points.push_back(point);
             }
             }

             point_cloud_ptr->width = (int)point_cloud_ptr->points.size();
             point_cloud_ptr->height = 1;

             return point_cloud_ptr;

         }
*/
//==============================================================================================================================================================
static void calc_pseudo(cv::Mat depth_image,cv::Mat  &disp_depth)
{
  u_char      *ptr=disp_depth.data;
  uint16_t  *depth=(uint16_t*)depth_image.data;
  for(int i=0;i<depth_image.rows*depth_image.cols; ++i )
  {
    if(depth[i]>=num_disparities*16) 
      {ptr[3*i]=ptr[3*i+1]=ptr[3*i+2]=0;continue;} // TOO FAR
    int pval = mGamma[depth[i]];
    int lb = pval & 0xff;
			switch ( pval >> 8 )
			{
				case 0:
					ptr[3*i+0] = 255;
					ptr[3*i+1] = 255-lb;
					ptr[3*i+2] = 255-lb;
					break;
				case 1:
					ptr[3*i+0] = 255;
					ptr[3*i+1] = lb;
					ptr[3*i+2] = 0;
					break;
				case 2:
					ptr[3*i+0] = 255-lb;
					ptr[3*i+1] = 255;
					ptr[3*i+2] = 0;
					break;
				case 3:
					ptr[3*i+0] = 0;
					ptr[3*i+1] = 255;
					ptr[3*i+2] = lb;
					break;
				case 4:
					ptr[3*i+0] = 0;
					ptr[3*i+1] = 255-lb;
					ptr[3*i+2] = 255;
					break;
				case 5:
					ptr[3*i+0] = 0;
					ptr[3*i+1] = 0;
					ptr[3*i+2] = 255-lb;
					break;
                                // TOO NEAR
				default:
					ptr[3*i+0] = 100; 
					ptr[3*i+1] = 100;
					ptr[3*i+2] = 100;
					break;
    }
  }
}

//====================================================================================================================================
void* grab(void* arg)
{

//  C API
//  CvCapture *cap1, *cap2, *cam_left, *cam_right;
//  cap1=cvCreateCameraCapture(cv::CAP_FIREWIRE+0);
//  cap2=cvCreateCameraCapture(cv::CAP_FIREWIRE+1);
//  if((u_long)cvGetCaptureProperty(cap1,CV_CAP_PROP_GUID)==LOW_SPEED_ID)
//     cam_left=cap1,cam_right=cap2;
//  else
//     cam_left=cap2,cam_right=cap1;
//  cvSetCaptureProperty(cam_left,CV_CAP_PROP_FPS,30);
//  cvSetCaptureProperty(cam_right,CV_CAP_PROP_ISO_SPEED,800);
//  cvSetCaptureProperty(cam_right,CV_CAP_PROP_EXPOSURE,200);
//  cvSetCaptureProperty(cam_right,CV_CAP_PROP_GAIN,500);
//  cvSetCaptureProperty(cam_right,CV_CAP_PROP_FPS,120);
//  cvSetCaptureProperty(cam_right,CV_CAP_PROP_MODE,4);
//  while(1)
//  {
//    for(int i=0;i<10;i++)
//    {
//      cvGrabFrame(cam_left);
//      img_left = cv::cvarrToMat(cvRetrieveFrame(cam_left),true);
//
//      for(int ii=0;ii<4;ii++)
//      {
//        cvGrabFrame(cam_right);
//        img_right_container[ii] = cv::cvarrToMat(cvRetrieveFrame(cam_right),true);
//      }
//    }
//  }

  // C++ API
  cv::VideoCapture cap1(cv::CAP_FIREWIRE+0), cap2(cv::CAP_FIREWIRE+1);
  cv::VideoCapture left_cam, right_cam;
  if((u_long)cap1.get(CV_CAP_PROP_GUID)==LOW_SPEED_ID)
     left_cam=cap1,right_cam=cap2;
  else
     left_cam=cap2,right_cam=cap1;
  left_cam.set(CV_CAP_PROP_FPS,30);
  right_cam.set(CV_CAP_PROP_ISO_SPEED,800);
  right_cam.set(CV_CAP_PROP_EXPOSURE,300);
  right_cam.set(CV_CAP_PROP_GAIN,350);
  right_cam.set(CV_CAP_PROP_FPS,120);
  right_cam.set(CV_CAP_PROP_MODE,4);

  struct timeval start_time, end_time;

  cv::Mat img_left_,img_right_;

  while(true)
  {

#ifdef FREQ_CAM
    gettimeofday(&start_time, NULL);
    for(int i = 0;i<50;i++)
    {
#endif

    sem_wait(&mutex_t1);
    right_cam >> img_right_; img_right_container[0] = img_right_;
    right_cam >> img_right_; img_right_container[1] = img_right_;
    right_cam >> img_right_; img_right_container[2] = img_right_;
    right_cam >> img_right_; img_right_container[3] = img_right_;
    left_cam  >> img_left_;
    img_left_.copyTo(img_left);
    img_right_.copyTo(img_right);
    cv::imshow("left",img_left_),
    cv::imshow("right",img_right_);
    cvWaitKey(1);
    sem_post(&mutex_t1);
    // lock3 : to prevent tracking thread from taking the same resource
    sem_post(&lock_t3);
       
#ifdef FREQ_CAM
    }
    gettimeofday(&end_time, NULL);
    std::cout << 50/ ((end_time.tv_sec - start_time.tv_sec) + 
                      (end_time.tv_usec- start_time.tv_usec) * 1e-6) 
	      << " [Hz]"<< std::endl;
#endif

  }

  return 0;
}

//====================================================================================================================================
void* stereo(void* arg)
{

  for ( int i = 0; i < num_disparities*16 ; ++i )
  {
    float v = i/(num_disparities*16.0);
    v = powf(v, 3)* 6;
    mGamma[i] = v*6*256;
  }

  cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(0,0,0);
  cv::Rect roi1, roi2;

  cv::FileStorage fs("./intrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open in.yml" << std::endl;
    return 0;
  }
  cv::Mat M1, D1, M2, D2;
  fs["M1"] >> M1;
  fs["D1"] >> D1;
  fs["M2"] >> M2;
  fs["D2"] >> D2;

  fs.open("./extrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open ex.yml" << std::endl;
    return 0;
  }    
  cv::Mat R, T, R1, P1, R2, P2, Q;
  fs["R1"] >> R1;
  fs["P1"] >> P1;
  fs["R2"] >> R2;
  fs["P2"] >> P2;
  fs["Q"] >> Q;    

  cv::Size size_(640,480);
//  cv::stereoRectify( M1, D1, M2, D2, size_, R, T, R1, R2, P1, P2, Q,
//                     cv::CALIB_ZERO_DISPARITY, -1, size_, &roi1, &roi2 );
    
  cv::Mat map11, map12, map21, map22;
  cv::Mat img_left_r, img_right_r, img_left_r_g;
  cv::initUndistortRectifyMap(M1, D1, R1, P1, size_, CV_16SC2, map11, map12);
  cv::initUndistortRectifyMap(M2, D2, R2, P2, size_, CV_16SC2, map21, map22);
    
  int cn = 1;
  //int numberOfDisparities = 128; //multiple of 16
  int sgbmWinSize = 5; //filter window size
  sgbm->setBlockSize(sgbmWinSize);
  // left disparity
  //sgbm->setMinDisparity(-64);
  //sgbm->setNumDisparities(64+num_disparities);
  // right disparity
  sgbm->setMinDisparity(-num_disparities);
  sgbm->setNumDisparities(num_disparities+64);
  sgbm->setUniquenessRatio(10); // ratio to determine which points match
  sgbm->setSpeckleWindowSize(100); // dots in image filter
  sgbm->setSpeckleRange(1); // implicitly multiplied by 16
  sgbm->setDisp12MaxDiff(1); // left-right disparity check allowed pixel
  sgbm->setPreFilterCap(63);
  sgbm->setP1(8*cn*sgbmWinSize*sgbmWinSize);
  sgbm->setP2(32*cn*sgbmWinSize*sgbmWinSize);
  //sgbm->setMode(cv::StereoSGBM::MODE_SGBM);

  cv::Mat img_left_local, img_right_local;
  cv::Mat disp_image(size_,CV_8UC3);

  sleep(1); //to allow image reading first
  struct timeval start_time, end_time;

  int w = 640, h = 480;
  stereoInit(w,h);
// extern "C" void stereoUpload( const unsigned char *left, const unsigned char *right );
// extern "C" void stereoProcess();
// extern "C" void stereoDownload( float *disparityLeft, float *disparityRight );

  float disparityLeft[640*480], disparityRight[640*480]; 

//  pcl::visualization::CloudViewer viewer("Cloud Viewer");

  cv::Mat im(size_,CV_8UC1),im2(size_,CV_8UC1),im3(size_,CV_16SC1),im4(size_,CV_8UC1);

//  cv::Mat iml,imr;
//  iml = cv::imread("Left_01.png");
//  imr = cv::imread("Right_01.png");

  while(true)
  {

#ifdef FREQ_STEREO
    gettimeofday(&start_time, NULL);
    for(int i = 0;i<50;i++)
    {
#endif
    
    sem_wait(&mutex_t2);
    img_left.copyTo(img_left_local);
    img_right.copyTo(img_right_local);
    sem_post(&mutex_t2);
    //remap(iml , img_left_r, map11, map12, cv::INTER_LINEAR);
    //remap(imr, img_right_r, map21, map22, cv::INTER_LINEAR);
    remap(img_left_local , img_left_r, map11, map12, cv::INTER_LINEAR);
    remap(img_right_local, img_right_r, map21, map22, cv::INTER_LINEAR);
    cvtColor(img_left_r, img_left_r_g, CV_RGB2GRAY, 1);
    //cvtColor(img_right_r, img_right_r, CV_RGB2GRAY, 1);
    //sgbm->compute(img_left_r_g,img_right_r,disp);
    sgbm->compute(img_right_r,img_left_r_g,disp);
    calc_pseudo(-disp,disp_image);
    //disp.convertTo(disp_image, CV_8U, 255/(num_disparities*16.));
    cv::imshow("disparity1",disp_image); cvWaitKey(1);

/*
    stereoUpload(img_left_r_g.data, img_right_r.data);
    stereoProcess();
    stereoDownload(disparityLeft, disparityRight);
    uint16_t  *depth = (uint16_t*)im3.data;
    for(int i=0;i<640*480;i++)
    {
      depth[i] = disparityLeft[i]*255;
      if (disparityLeft[i] != 0)
        {std::cout << depth[i]<< "  " << disparityLeft[i] << "  " << depth2[i] << std::endl;}
    }
    im3.convertTo(im, CV_8U, 255/(num_disparities*16.));
    //calc_pseudo(im3,im);
    cv::imshow("disparity_cuda2",im);  cvWaitKey(1);
*/

    //for(int i=0;i<640*480;i++)
    //  im3.data[i] = disparityRight[i];
    //im3.convertTo(im2, CV_8U, 255/(num_disparities*16.));
    //calc_pseudo(im3,im2);
    //imshow("test2",im2);

    // 3D PROJECTION
    //cv::reprojectImageTo3D(disp/16.0,point3D,Q,false);
    
//    std::vector<cv::Point3f> t1(1);
//    std::vector<cv::Point3f> t2(1);
//    cv::Point3f t3(320.0,240.0,disp.at<int16_t> (240,320)/16.0);
//    t1[0] = t3;
//    cv::perspectiveTransform(t1, t2, Q);
//    std::cout << point3D.at<cv::Vec3f> (240,320) << "   ";
//    fflush(stdout);
//    std::cout << t2[0] << std::endl;
//    fflush(stdout);

//    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud;
//    cloud = MatToPoinXYZ(point3D,img_left_r);
//    viewer.showCloud(cloud);


#ifdef FREQ_STEREO
    }
    gettimeofday(&end_time, NULL);
    std::cout << 50/ ((end_time.tv_sec - start_time.tv_sec) + 
                      (end_time.tv_usec- start_time.tv_usec) * 1e-6) 
	      << " [Hz]"<< std::endl;
#endif

//    sgbm->compute(img_left_r_g,img_right_r,disp);
//    calc_pseudo(disp,disp_image);
//    cv::imshow("disparity2",disp_image);
//    cvWaitKey(1);
  }

  return 0;
}

//====================================================================================================================================
void* tracker(void* arg)
{

  sleep(2); // make thread slower than the first thread

  cv::FileStorage fs("./intrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open in.yml" << std::endl;
    return 0;
  }
  cv::Mat M2;
  fs["M2"] >> M2; // right cam matrix

  fs.open("./extrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open ex.yml" << std::endl;
    return 0;
  }    
  cv::Mat Q;
  fs["Q"] >> Q;  

  int width,height,posx,posy;
  int num_points = 4;

  //XVMatrix vec_1(3,num_points);
  //XVMatrix vec_2(3,num_points);
  XVMatrix vec_1_tmp(3,num_points);
  XVMatrix vec_2_tmp(3,num_points);
  XVColVector D1(num_points);
  XVColVector D2(num_points);

  KLT_TrackingContext   tc;
  KLT_FeatureList       fl;
  KLT_FeatureTable      ft;
  tc = KLTCreateTrackingContext();
  fl = KLTCreateFeatureList(num_points);
  tc->sequentialMode = TRUE;
  tc->writeInternalImages = FALSE;
  tc->affineConsistencyCheck = -1;
  tc->sequentialMode = true;
  tc->bindToCorner=true;
  tc->patchMode=true;
  tc->mindist=30;
  tc->writeInternalImages = FALSE;
  tc->window_width=9;
  tc->window_height=9;
  tc->affineConsistencyCheck = -1; 
  KLTChangeTCPyramid(tc, 10);
  KLTUpdateTCBorder(tc);

  width  = 100;
  height = 100;
  posx   = 260;
  posy   = 150;
  IppiRect roi={posx,posy,width,height};

  // showing tracking region
  cv::Mat rgb_tmp = cv::Mat::zeros(480,640, CV_8UC1);
  cv::Mat rgb_tmp2 = cv::Mat::zeros(480,640, CV_8UC1);
  cv::Mat tmp_img = img_right_container[3];
  tmp_img.rowRange(posy,posy+height).copyTo(rgb_tmp.rowRange(posy,posy+height));
  rgb_tmp.colRange(posx,posx+width).copyTo(rgb_tmp2.colRange(posx,posx+width));
  cv::imshow("tracking_region",rgb_tmp2); cvWaitKey(1); 

  XVImageScalar<unsigned char> gray_image(640,480);
  XVImageScalar<unsigned char> gray_image_tmp(640,480);
  XVDrawWindowX<XV_RGB> window(649,480,0,0);
  window.map();
  window.swap_buffers();
  window.flush();
  gray_image.remap((unsigned char*)tmp_img.data,false);
  window.CopySubImage(gray_image);
  window.swap_buffers();
  window.flush();
  
  bool flag = true;
  while(flag)
  {
    flag = false;
    KLTExtSelectGoodFeatures(tc, (KLT_PixelType*)gray_image.data(), 640, 480, fl, &roi, 2, 2);
    for(int j=0;j<num_points;j++) 
      for(int jj=0;jj<num_points;jj++) 
        if(jj!=j)
          if(fl->feature[j]->y == fl->feature[jj]->y && fl->feature[j]->x == fl->feature[jj]->x)
            flag = true;
    if(KLTCountRemainingFeatures(fl)!=num_points) flag = true;
    tmp_img = img_right_container[3];
    tmp_img.rowRange(posy,posy+height).copyTo(rgb_tmp.rowRange(posy,posy+height));
    rgb_tmp.colRange(posx,posx+width).copyTo(rgb_tmp2.colRange(posx,posx+width));
    cv::imshow("tracking_region",rgb_tmp2); cvWaitKey(1); 
    gray_image.remap((unsigned char*)tmp_img.data,false);
  }
    //KLTExtReplaceLostFeatures(tc, (KLT_PixelType*)gray_image.data(), 640, 480, fl, &roi);

  for(int j=0;j<num_points;j++)
    if(fl->feature[j]->val>=0)
      window.fillEllipse(fl->feature[j]->x-5, fl->feature[j]->y-5, 10, 10, "red");  
  window.swap_buffers();
  window.flush();

  cv::Mat disp_; disp.copyTo(disp_);

  cv::Vec3f point,point1,point2;
  std::vector<cv::Vec3f> point_(num_points);
  //std::vector<cv::Point3f> point_(num_points);
  std::vector<cv::Point3f> p1_(1),p2_(1);
  int counter = 0;
  //for(int j=0;j<num_points;j++) if(fl->feature[j]->val>=0) counter += 1;

  cv::reprojectImageTo3D(disp_/16.0,point3D,Q,true);

  for(int j=0;j<num_points;j++)
//  if(fl->feature[j]->val>=0)
  {
    //cv::Point3f tmp(fl->feature[j]->x,fl->feature[j]->y,-disp_.at<int16_t> (fl->feature[j]->y,fl->feature[j]->x)/16.0);
    //p1_[0] = tmp;
    //cv::perspectiveTransform(p1_, p2_, Q);
    //point_[j] = p2_[0];
    //vec_1[0][j]=p2_[0].x;
    //vec_1[1][j]=p2_[0].y;
    //vec_1[2][j]=p2_[0].z;
    point = point3D.at<cv::Vec3f> (fl->feature[j]->y,fl->feature[j]->x);
//    vec_1[0][j]=point[0];
//    vec_1[1][j]=point[1];
//    vec_1[2][j]=point[2];         
//    D1[j]=sqrt(Sqr(vec_1[0][j])+Sqr(vec_1[1][j])+Sqr(vec_1[2][j]));
//    vec_1[0][j]/=D1[j];
//    vec_1[1][j]/=D1[j];
//    vec_1[2][j]/=D1[j];
    point_[j] = point;
  }

  KLTSetVerbosity(0);
  
  TransfMatrix Trans;
  TransfMatrix Trans2;
  XVMatrix  R(3,3);
  XVColVector T(3);
  Trans.R=R;
  Trans.T=T;
  Trans2.R=R;
  Trans2.T=T;
  Vgps vgps(Trans.R,Trans.T);

  struct timeval start_time, end_time;

  std::vector<cv::Mat> img_container(4);

  while(true)
  {

#ifdef FREQ_TRACKER
    gettimeofday(&start_time, NULL);
    for(int i = 0;i<50;i++)
    {
#endif

      sem_wait(&lock_t3);
      sem_wait(&mutex_t3);
      img_container = img_right_container;
      sem_post(&mutex_t3);

      for(int ii=0;ii<4;ii++)
      {
        gray_image.remap((unsigned char*)img_container[ii].data,false);
        KLTTrackFeatures(tc, (KLT_PixelType*)gray_image.data(), (KLT_PixelType*)gray_image.data(), 640, 480, fl);
        window.CopySubImage(gray_image);
        for(int k=0;k<num_points;k++)
          if(fl->feature[k]->val>=0)
            window.fillEllipse(fl->feature[k]->x-5, fl->feature[k]->y-5, 10, 10, "green");
        window.swap_buffers();
        window.flush();
        //while(KLTCountRemainingFeatures(fl)!=num_points)
        //  KLTExtReplaceLostFeatures(tc, (KLT_PixelType*)gray_image.data(), 640, 480, fl, &roi);
      }

      counter = 0;
      
      //printf("x: %f y: %f z: %f\n",vec_1[0][0],vec_1[1][0],vec_1[2][0]);

      XVMatrix vec_1(3,num_points);
      XVMatrix vec_2(3,num_points);
     
      //std::cout << M2.at<double>(0,2) << M2.at<double>(1,2) << M2.at<double>(0,0) << M2.at<double>(1,1) << std::endl;

      for(int j=0;j<num_points;j++)
      //if(fl->feature[j]->val>=0)
      {
        //if(fl->feature[j]->val>=0) counter++;
        //cv::Point3f tmp(fl->feature[j]->x,fl->feature[j]->y,disp_.at<int16_t> (fl->feature[j]->y,fl->feature[j]->x)/16.0);
        //printf("xpix: %f ypix: %f M2x: %f M2y: %f\n",fl->feature[j]->x,fl->feature[j]->y,M2.at<double>(0,0),M2.at<double>(1,1));
        //p1_[0] = tmp;
        //cv::perspectiveTransform(p1_, p2_, Q);
        //vec_2[0][j]=p2_[0].x;
        //vec_2[1][j]=p2_[0].y;
        //vec_2[2][j]=p2_[0].z;
        //point = point3D.at<cv::Vec3f>(fl->feature[j]->y,fl->feature[j]->x);
        //printf("x: %f y: %f z: %f x: %f y: %f z: %f\n",p2_[0].x,p2_[0].y,p2_[0].z,point[0],point[1],point[2]);
        //vec_2[0][j]=point[0];
        //vec_2[1][j]=point[1];
        //vec_2[2][j]=point[2];
        vec_2[0][j]=fl->feature[j]->x - M2.at<double>(0,2);
        vec_2[1][j]=fl->feature[j]->y - M2.at<double>(1,2);
        vec_2[2][j]=(M2.at<double>(0,0) + M2.at<double>(1,1))/2;
        D2[j]=sqrt(Sqr(vec_2[0][j])+Sqr(vec_2[1][j])+Sqr(vec_2[2][j]));
        vec_2[0][j]/=D2[j];
        vec_2[1][j]/=D2[j];
        vec_2[2][j]/=D2[j];
        //printf("x: %f y: %f z: %f\n",vec_2[0][j],vec_2[1][j],vec_2[2][j]);

          //vec_1[0][j]=point_[j].x;
          //vec_1[1][j]=point_[j].y;
          //vec_1[2][j]=point_[j].z;
          vec_1[0][j]=point_[j][0];
          vec_1[1][j]=point_[j][1];
          vec_1[2][j]=point_[j][2];   
          D1[j]=sqrt(Sqr(vec_1[0][j])+Sqr(vec_1[1][j])+Sqr(vec_1[2][j]));      
          vec_1[0][j]/=D1[j];
          vec_1[1][j]/=D1[j];
          vec_1[2][j]/=D1[j];
          //printf("x: %f y: %f z: %f\n",vec_1[0][j],vec_1[1][j],vec_1[2][j]);

      }
      //std::cout << counter << std::endl;

      //printf("x: %f y: %f z: %f\n",vec_1[0][0],vec_1[1][0],vec_1[2][0]);

/*      printf("Points1 :\n");
      printf("%f %f\n%f %f\n%f %f\n",
           vec_1[0][0],vec_2[0][0],
           vec_1[1][0],vec_2[1][0],
           vec_1[2][0],vec_2[2][0]);
*/

      printf("%f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f\n",
           vec_1[0][0],vec_2[0][0],vec_1[0][1],vec_2[0][1],vec_1[0][2],vec_2[0][2],vec_1[0][3],vec_2[0][3],
           vec_1[1][0],vec_2[1][0],vec_1[1][1],vec_2[1][1],vec_1[1][2],vec_2[1][2],vec_1[1][3],vec_2[1][3],
           vec_1[2][0],vec_2[2][0],vec_1[2][1],vec_2[2][1],vec_1[2][2],vec_2[2][2],vec_1[2][3],vec_2[2][3]);

//      vec_1_tmp = vec_1;
      vgps.find_pose(vec_2,vec_1,D1,Trans.R,Trans.T,1200);
      Trans2.T=Trans.R.t()*Trans.T;

      //if(!vgps.validate_set(vec_2,vec_1,Trans.R,Trans.T))
        disp.copyTo(disp_); 

//      if(vgps.validate_set(vec_2,vec_1,Trans.R,Trans.T))
//{
      //printf("Feature Count: %d   D1: %f   D2: %f   ",KLTCountRemainingFeatures(fl),D1[0],D2[0]);
      //std::cout << "T limit(1e-3) : " << Trans.T.ip();
      std::cout << "Validate : " << vgps.validate_set(vec_2,vec_1,Trans.R,Trans.T) << std::endl;
//      printf("Points :\n");
//      printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",
//           vec_1[0][0],vec_1[0][1],vec_1[0][2],vec_1[0][3],
//           vec_1[1][0],vec_1[1][1],vec_1[1][2],vec_1[1][3],
//           vec_1[2][0],vec_1[2][1],vec_1[2][2],vec_1[2][3]);

//      printf("Points2 :\n");
//      printf("%f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f\n",
//           vec_1[0][0],vec_2[0][0],vec_1[0][1],vec_2[0][1],vec_1[0][2],vec_2[0][2],vec_1[0][3],vec_2[0][3],
//           vec_1[1][0],vec_2[1][0],vec_1[1][1],vec_2[1][1],vec_1[1][2],vec_2[1][2],vec_1[1][3],vec_2[1][3],
//           vec_1[2][0],vec_2[2][0],vec_1[2][1],vec_2[2][1],vec_1[2][2],vec_2[2][2],vec_1[2][3],vec_2[2][3]);

      printf("[R T] :\n");
      printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",
           Trans.R[0][0],Trans.R[0][1],Trans.R[0][2],Trans2.T[0],
           Trans.R[1][0],Trans.R[1][1],Trans.R[1][2],Trans2.T[1],
           Trans.R[2][0],Trans.R[2][1],Trans.R[2][2],Trans2.T[2]);

//}

#ifdef FREQ_TRACKER
    }

 
      printf("Feature Count: %d   D1: %f   ",KLTCountRemainingFeatures(fl),D1[0]);
      //std::cout << "T limit(1e-3) : " << Trans.T.ip();
      std::cout << "Validate : " << vgps.validate_set(vec_2,vec_1,Trans.R,Trans.T) << std::endl;
/*
      printf("[R T] :\n");
      printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",
           Trans.R[0][0],Trans.R[0][1],Trans.R[0][2],Trans.T[0],
           Trans.R[1][0],Trans.R[1][1],Trans.R[1][2],Trans.T[1],
           Trans.R[2][0],Trans.R[2][1],Trans.R[2][2],Trans.T[2]);
*/

    //printf("Feature Count: %d\n",KLTCountRemainingFeatures(fl));
    gettimeofday(&end_time, NULL);
    std::cout << (50 * 4) / ((end_time.tv_sec - start_time.tv_sec) + 
                      (end_time.tv_usec- start_time.tv_usec) * 1e-6) 
	      << " [Hz]"<< std::endl;
#endif

  }

  return 0;
}


//==============================================================================================================================================================
int main(int argc, char **argv)
{

  cv::namedWindow("tracking_region"); cvMoveWindow("tracking_region",0,600);
  cv::namedWindow("left");
  cv::namedWindow("right");
  cv::namedWindow("disparity1"); cvMoveWindow("disparity1",0,600);
  //cv::namedWindow("disparity2"); cvMoveWindow("disparity2",1000,0);

  // Start multithread
  pthread_t thread_grab, thread_stereo, thread_tracker;
//  pthread_t thread_grab, thread_stereo;

  sem_init(&lock_t1, 0, 1);
  sem_init(&lock_t2, 0, 0);
  sem_init(&lock_t3, 0, 0);
  sem_init(&mutex_t1, 0, 1);
  sem_init(&mutex_t2, 0, 1);
  sem_init(&mutex_t3, 0, 1);

  pthread_attr_t attr;
  cpu_set_t cpus;
  pthread_attr_init(&attr);

  CPU_ZERO(&cpus);
  CPU_SET(1, &cpus);
  pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
  pthread_create(&thread_grab, &attr, grab, NULL);

  CPU_ZERO(&cpus);
  CPU_SET(2, &cpus);
  pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
  pthread_create(&thread_stereo, &attr, stereo, NULL);

  CPU_ZERO(&cpus);
  CPU_SET(3, &cpus);
  pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
  pthread_create(&thread_tracker, &attr, tracker, NULL);

  pthread_join(thread_grab, NULL);
  pthread_join(thread_stereo, NULL);
  pthread_join(thread_tracker, NULL);

  return 0;
}
