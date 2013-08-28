###############################################################################
# Example makefile for portaudiov19 + fftw3 + openGL projects.
# ihsan@kehribar.me wrote this file.
###############################################################################
ifeq ($(shell uname), Darwin)
	LIBS = `/opt/local/bin/pkg-config --libs portaudio-2.0 fftw3`
	LIBS += -framework GLUT -framework OpenGL
	CFLAGS = `/opt/local/bin/pkg-config --cflags gl portaudio-2.0 fftw3`
else
	LIBS = `pkg-config --libs portaudio-2.0 fftw3`
	CFLAGS = `pkg-config --cflags portaudio-2.0 fftw3`
	LIBS += -lGL -lglut -lGLU
endif

TARGET = test

.PHONY: clean

$(TARGET):
	g++ $(TARGET).c -o $(TARGET) $(CFLAGS) $(LIBS)

clean:
	rm -f *.o *~ ./$(TARGET) 
