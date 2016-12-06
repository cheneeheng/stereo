PROG =   grab
OBJS =   main.o
CXXFLAGS = -I/usr/local/include
LFLAGS   = -L/usr/local/lib -lopencv_core -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video

all: $(OBJS)
	$(CXX) -o $(PROG) $(OBJS) $(LFLAGS)

clean: ;	rm -f $(OBJS) $(PROG)
