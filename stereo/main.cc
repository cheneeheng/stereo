#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#include <sys/time.h>
#include <dc1394/control.h>
#include <iostream>
#include <opencv2/videoio.hpp>
#include "opencv2/calib3d.hpp"

#define HIGH_SPEED_ID 302662258
#define LOW_SPEED_ID  235376646

const int  	      num_disparities=128;
static int16_t        mGamma[num_disparities*16];

//==============================================================================================================================================================
static void calc_pseudo(cv::Mat depth_image,cv::Mat  disp_depth)
{
		  u_char      *ptr=disp_depth.data;
		  uint16_t  *depth=(uint16_t*)depth_image.data;
		  for(int i=0;i<depth_image.rows*depth_image.cols; ++i )
		  {
			if(depth[i]>=num_disparities*16) 
			               {ptr[3*i]=ptr[3*i+1]=ptr[3*i+2]=0;continue;}
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

				default:
					ptr[3*i+0] = 0;
					ptr[3*i+1] = 0;
					ptr[3*i+2] = 0;
					break;
			}
		  }
}


//==============================================================================================================================================================
int main(int argc, char **argv)
{
  CvCapture *cap1,*cap2, *left_cam,*right_cam;
  struct timeval start_time, end_time;

  bool        active = true;
  static char buffer_left[80],buffer_right[80];
  int         index = 1;

  cap1=cvCreateCameraCapture(cv::CAP_FIREWIRE+0);
  cap2=cvCreateCameraCapture(cv::CAP_FIREWIRE+1);
  if((u_long)cvGetCaptureProperty(cap1,CV_CAP_PROP_GUID)==LOW_SPEED_ID)
     left_cam=cap1,right_cam=cap2;
  else
     left_cam=cap2,right_cam=cap1;
  cvSetCaptureProperty(left_cam,CV_CAP_PROP_FPS,30);

  cvSetCaptureProperty(right_cam,CV_CAP_PROP_ISO_SPEED,800);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_EXPOSURE,300);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_GAIN,500);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_FPS,120);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_MODE,4);

  cv::Mat im1,im2;
  cv::namedWindow("left");
  cv::namedWindow("right");
  cvGrabFrame(left_cam);
  cvGrabFrame(right_cam);

  cv::namedWindow("left");
  cv::namedWindow("right");
  cv::namedWindow("disparity");
  cvMoveWindow("disparity",650,0);


  cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(0,0,0);
  //cv::StereoSGBM sgbm;
  cv::Rect roi1, roi2;

  cv::FileStorage fs("../calib/images/intrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open in.yml" << std::endl;
    return -1;
  }
  cv::Mat M1, D1, M2, D2;
  fs["M1"] >> M1;
  fs["D1"] >> D1;
  fs["M2"] >> M2;
  fs["D2"] >> D2;

  fs.open("../calib/images/extrinsics.yml", CV_STORAGE_READ);
  if(!fs.isOpened())
  {
    std::cerr << "Couldn't open ex.yml" << std::endl;
    return -1;
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
  cv::Mat img1r, img2r;
  cv::initUndistortRectifyMap(M1, D1, R1, P1, size_, CV_16SC2, map11, map12);
  cv::initUndistortRectifyMap(M2, D2, R2, P2, size_, CV_16SC2, map21, map22);

  for ( int i = 0; i < num_disparities*16 ; ++i )
  {
    float v = i/(num_disparities*16.0);
    v = powf(v, 3)* 6;
    mGamma[i] = v*6*256;
  }

  int numberOfDisparities = 80;
  int SADWindowSize = 3;
  int sgbmWinSize;
  sgbmWinSize = SADWindowSize;
  sgbm->setPreFilterCap(63);
  sgbm->setBlockSize(sgbmWinSize);
    
  int cn = 3;
    
//  sgbm.P1 = 8*cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
//  sgbm.P2 = 32*cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
//  sgbm.minDisparity = 0;
//  sgbm.numberOfDisparities = 80;
//  sgbm.uniquenessRatio = 10;
//  sgbm.speckleWindowSize = 100;
//  sgbm.speckleRange = 32;
//  sgbm.disp12MaxDiff = 1;
//  sgbm.fullDP =1;
    
  sgbm->setP1(8*cn*sgbmWinSize*sgbmWinSize);
  sgbm->setP2(32*cn*sgbmWinSize*sgbmWinSize);
  sgbm->setMinDisparity(0);
  sgbm->setNumDisparities(numberOfDisparities);
  sgbm->setUniquenessRatio(10);
  sgbm->setSpeckleWindowSize(100);
  sgbm->setSpeckleRange(32);
  sgbm->setDisp12MaxDiff(1);

  cv::Mat disp;
  cv::Mat disp_image(size_,CV_8UC3);
  cv::Mat point3D;


   do 
   {
     im1=cv::cvarrToMat(cvRetrieveFrame(left_cam),true);
     im2=cv::cvarrToMat(cvRetrieveFrame(right_cam),true);
     cvGrabFrame(left_cam), cvGrabFrame(right_cam);
     remap(im1, img1r, map11, map12, cv::INTER_LINEAR);
     remap(im2, img2r, map21, map22, cv::INTER_LINEAR);
//     sgbm->compute(img1r,img2r,disp);

//     calc_pseudo(disp,disp_image);
//     cv::reprojectImageTo3D(disp,point3D,Q,true);

//(point3D.at<cv::Vec3f>(t_fl->feature[i]->y,
//                                         t_fl->feature[i]->x)[2]>3 ||
//                  fabs(point3D.at<cv::Vec3f>(t_fl->feature[i]->y,
//                                       t_fl->feature[i]->x)[2]-
//                     point3D.at<cv::Vec3f>(t_fl->feature[i]->y-1,
//                                   t_fl->feature[i]->x)[2])>0.4 ||
//              fabs(point3D.at<cv::Vec3f>(t_fl->feature[i]->y,
//                                   t_fl->feature[i]->x)[2]-


     cv::imshow("left",img1r);
     cv::imshow("right",img2r);
//     cv::imshow("disparity",disp_image);
     switch(cvWaitKey(1))
     {
	case 0x1b:
	   active=false;
	   break;
        default:
	  break;
     }
   }while(active);


   return 0;
}
