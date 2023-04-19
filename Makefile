CC=g++
NOWARNTAG=-Wno-varargs -Wno-deprecated -Wno-unused-function -Wno-unused-variable -Wno-format -Wno-sign-compare -Wno-unused-parameter -Wno-attributes
CPPFLAGS=-std=c++11 -g -O3 ${NOWARNTAG}
LINKLIB=-lz -lpthread -lhts
BWASRC=$(wildcard ./src/bwalib/*.cpp)
GETSRC=$(wildcard ./src/*.cpp)
BWAOBJ=$(patsubst %.cpp,%.o,$(BWASRC))
GETOBJ=$(patsubst %.cpp,%.o,$(GETSRC))
SRCDIR=./src/
BWADIR=./src/bwalib/
PROG=getools


all:$(PROG)

$(PROG):$(BWAOBJ) $(GETOBJ)
	$(CC) $(CPPFLAGS) $^ -o $@ $(LINKLIB)

%.o:$(BWADIR)%.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

%.o:$(SRCDIR)%.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

.PHONY:all clean cleanlocal

cleanlocal:
	rm -f $(PROG) $(BWAOBJ) $(GETOBJ)

clean:cleanlocal
