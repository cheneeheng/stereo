PROG = stereo_calib
CXXFLAGS = -O3 -I/usr/local/include 
LIBS= -L/usr/local/lib -lopencv_core -lopencv_highgui -lopencv_features2d -lopencv_calib3d -lopencv_imgproc -lopencv_imgcodecs
OBJS= stereo_calib.o

all: $(OBJS)
	$(CXX) -o $(PROG) $(OBJS) $(LIBS)
clean: ;	rm -f *.o $(PROG)
