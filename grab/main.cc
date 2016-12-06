#include <sys/time.h>
#include <dc1394/control.h>
#include <iostream>

#include "opencv2/opencv.hpp"
#include "opencv2/videoio.hpp"

#define HIGH_SPEED_ID 302662258
#define LOW_SPEED_ID  235376646

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

  do 
  {
    im1=cv::cvarrToMat(cvRetrieveFrame(left_cam),true);
    im2=cv::cvarrToMat(cvRetrieveFrame(right_cam),true);
     cvGrabFrame(left_cam);
     cvGrabFrame(right_cam);
     cv::imshow("left",im1);
     cv::imshow("right",im2);

     switch(cvWaitKey(1))
     {
	case 0x1b:
	   active=false;
	   break;
        default:
          sprintf(buffer_left,"Left_%02d.png",index);
          cv::imwrite(buffer_left,im1);
          sprintf(buffer_right,"Right_%02d.png",index);
          cv::imwrite(buffer_right,im2);
	  std::cout <<"grabbed " << index << std::endl;
	  index++;
        case -1:
	  break;
     }
  }while(active);

  return 0;
}









