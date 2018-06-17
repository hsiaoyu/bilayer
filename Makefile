CC = g++

#flags#
CFLAG = -c
COMPILE_FLAGS = -std=c++11 -O2
INCLUDES = -I/var/local/libigl/include/ -I/var/local/eigen 
LIBS = -L /h1/hsiaoyu/program/usr/local/lib  -L /h1/hsiaoyu/program/usr/local/lib64/ -lglfw -lGLEW -lGLU -lGL

SOURCES = simSetup.cpp main.cpp energy1.cpp energy2.cpp energy3.cpp energy4.cpp energy5.cpp \
Dd01.cpp Dd02.cpp Dd03.cpp Dd11.cpp Dd12.cpp Dd13.cpp Dd21.cpp Dd22.cpp Dd23.cpp \
Ddtop01.cpp Ddtop02.cpp Ddtop03.cpp Ddtop11.cpp Ddtop12.cpp Ddtop13.cpp Ddtop21.cpp Ddtop22.cpp Ddtop23.cpp \
Dv01.cpp Dv02.cpp Dv03.cpp Dv11.cpp Dv12.cpp Dv13.cpp Dv21.cpp Dv22.cpp Dv23.cpp \
DDd01.cpp DDd02.cpp DDd03.cpp DDd11.cpp DDd12.cpp DDd13.cpp DDd21.cpp DDd22.cpp DDd23.cpp \
DDdtop01.cpp DDdtop02.cpp DDdtop03.cpp DDdtop11.cpp DDdtop12.cpp DDdtop13.cpp DDdtop21.cpp DDdtop22.cpp DDdtop23.cpp \
DDv01.cpp DDv02.cpp DDv03.cpp DDv11.cpp DDv12.cpp DDv13.cpp DDv21.cpp DDv22.cpp DDv23.cpp \

OBJECTS = $(SOURCES:.cpp=.o) 
EXECUTABLE = bilayer

.PHONY:clean

all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

%.o : %.cpp
	$(CC) $(INCLUDES) $(COMPILE_FLAGS) $(CFLAG) $< -o $@

clean:
	rm *o $(EXECUTABLE)
