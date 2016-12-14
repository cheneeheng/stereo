cmake_minimum_required(VERSION 2.8)
find_package( OpenCV REQUIRED )

include_directories( ~ /usr/local/src/XVision2/include /usr/local/intel/ipp/5.2/em64t/include /usr/include/ffmpeg ~/KINECT/TRACKER/klt-elmar/Source ~/KINECT/TRACKER/vgps)

link_directories( /usr/local/src/XVision2/lib /usr/local/intel/ipp/5.2/em64t/sharedlib ~/KINECT/TRACKER/klt-elmar/Release ~/KINECT/TRACKER/vgps)

add_executable( main main.cc )

target_link_libraries( main ${OpenCV_LIBS} XVSeg XVCons XVTrack XVDevs XVTools XVImages png Xext X11 ippccem64t ippcvem64t ippiem64t pthread klt-elmar vgps)

