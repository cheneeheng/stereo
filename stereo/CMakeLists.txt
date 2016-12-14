cmake_minimum_required(VERSION 2.8)
find_package( OpenCV REQUIRED )
find_package(VTK REQUIRED)
find_package(PCL REQUIRED)
##set constant
##set flag
SET(GCC_COVERAGE_COMPILE_FLAGS "-DHAVE_IPP -g")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
include_directories( ~/stereo/stereo/libs /usr/local/include /usr/local/cuda/include /usr/local/src/XVision2/include /usr/local/intel/ipp/5.2/em64t/include /usr/include/ffmpeg ~/stereo/stereo/libs/klt-elmar/Source ${PCL_INCLUDE_DIRS})
include(${VTK_USE_FILE})
add_definitions(${PCL_DEFINITIONS})
link_directories( /usr/local/lib /usr/local/cuda/lib64 /usr/local/src/XVision2/lib /usr/local/intel/ipp/5.2/em64t/sharedlib ~/stereo/stereo/libs/klt-elmar/Release ~/stereo/stereo/libs/vgps ${PCL_LIBRARY_DIRS})
add_executable( stereo main.cc )
target_link_libraries( stereo ${OpenCV_LIBS} ${VTK_LIBRARIES} ${PCL_LIBRARIES} cuda cudart Stereo XVSeg XVCons XVTrack XVDevs XVTools XVImages png Xext X11 ippccem64t ippcvem64t ippiem64t pthread klt-elmar vgps)



#cmake_minimum_required(VERSION 2.8)
#find_package( OpenCV REQUIRED )
#set constant
#set flag
#SET(GCC_COVERAGE_COMPILE_FLAGS "-DHAVE_IPP -g")
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
#include_directories( ~/stereo/stereo/libs /usr/local/include /usr/local/cuda/include /usr/local/src/XVision2/include /usr/local/intel/ipp/5.2/em64t/include /usr/include/ffmpeg ~/stereo/stereo/libs/klt-elmar/Source)
#link_directories( /usr/local/lib /usr/local/cuda/lib64 /usr/local/src/XVision2/lib /usr/local/intel/ipp/5.2/em64t/sharedlib ~/stereo/stereo/libs/klt-elmar/Release ~/stereo/stereo/libs/vgps)
#add_executable( stereo main.cc )
#target_link_libraries( stereo ${OpenCV_LIBS} cuda cudart Stereo XVSeg XVCons XVTrack XVDevs XVTools XVImages png Xext X11 ippccem64t ippcvem64t ippiem64t pthread klt-elmar vgps)


