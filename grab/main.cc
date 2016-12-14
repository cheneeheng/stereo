#include <sys/time.h>
#include <dc1394/control.h>
#include <iostream>

#include "opencv2/opencv.hpp"
#include "opencv2/videoio.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#define HIGH_SPEED_ID 302662258
#define LOW_SPEED_ID  235376646

int main(int argc, char **argv)
{

  bool calib = false;
  bool write = false;
  cv::CommandLineParser parser(argc, argv, "{calib|false|},{write|false|}");
  calib = parser.get<bool>("calib"); 
  write = parser.get<bool>("write"); 

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
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_GAIN,450);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_FPS,120);
  cvSetCaptureProperty(right_cam,CV_CAP_PROP_MODE,4);

  cv::Mat im1,im2;
  cv::namedWindow("left"); cv::moveWindow("left",0,0);
  cv::namedWindow("right"); cv::moveWindow("right",0,550);
  cv::namedWindow("left2"); cv::moveWindow("left2",650,0);
  cv::namedWindow("right2"); cv::moveWindow("right2",650,550);
  cvGrabFrame(left_cam);
  cvGrabFrame(right_cam);

if(calib)
{
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
  cv::Mat map11, map12, map21, map22;
  cv::Mat img1r, img2r;
  cv::initUndistortRectifyMap(M1, D1, R1, P1, size_, CV_16SC2, map11, map12);
  cv::initUndistortRectifyMap(M2, D2, R2, P2, size_, CV_16SC2, map21, map22);

  do 
  {
    im1=cv::cvarrToMat(cvRetrieveFrame(left_cam),true);
    im2=cv::cvarrToMat(cvRetrieveFrame(right_cam),true);
     cvGrabFrame(left_cam);
     cvGrabFrame(right_cam);
     cv::imshow("left",im1);
     cv::imshow("right",im2);

     remap(im1, img1r, map11, map12, cv::INTER_LINEAR);
     remap(im2, img2r, map21, map22, cv::INTER_LINEAR);
     cv::imshow("left2",img1r);
     cv::imshow("right2",img2r);

    cvWaitKey(1);

/*
     switch(cvWaitKey(1))
     {
	case 0x1b:
	   active=false;
	   break;
        default:
          sprintf(buffer_left, "../calib/images/Left_%02d.png", index);
          cv::imwrite(buffer_left,im1);
          sprintf(buffer_right,"../calib/images/Right_%02d.png",index);
          cv::imwrite(buffer_right,im2);
	  std::cout <<"grabbed " << index << std::endl;
	  index++;
        case -1:
	  break;
     }
*/  }while(active);

}
else
{

  do 
  {
    im1=cv::cvarrToMat(cvRetrieveFrame(left_cam),true);
    im2=cv::cvarrToMat(cvRetrieveFrame(right_cam),true);
     cvGrabFrame(left_cam);
     cvGrabFrame(right_cam);
     cv::imshow("left",im1);
     cv::imshow("right",im2);

if(write){
     switch(cvWaitKey(1))
     {
	case 0x1b:
	   active=false;
	   break;
        default:
          sprintf(buffer_left, "../calib/images/Left_%02d.png", index);
          cv::imwrite(buffer_left,im1);
          sprintf(buffer_right,"../calib/images/Right_%02d.png",index);
          cv::imwrite(buffer_right,im2);
	  std::cout <<"grabbed " << index << std::endl;
	  index++;
        case -1:
	  break;
     }
}
else
  cvWaitKey(1);




  }while(active);

}

  return 0;
}









