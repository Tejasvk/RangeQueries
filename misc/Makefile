appname := a.out

CXX := g++
#CXXFLAGS := -Wall -g -std=c++17 -O3 
CXXFLAGS :=  -pipe -std=c++17 -O3 

#srcfiles := $(shell find . -maxdepth 1 -name "main.cpp")
srcfiles := \
	$(wildcard sources/*.cpp) \
	$(wildcard main.cpp) 
	
	
	 

objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend
