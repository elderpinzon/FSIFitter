CXX=$(shell root-config --cxx)
ROO=rootcint

ROOTFLAGS=$(shell root-config --cflags)
CFLAGS=-c -Wall -v -g -std=gnu++0x -O2 $(ROOTFLAGS)

LDFLAGS=-g -O2 $(shell root-config --ldflags) -DVERBOSE
LDLIBS=$(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -std=gnu++0x

INCLUDE=-I$(shell root-config --incdir) -I./
SOFLAGS=-fPIC -shared $(INCLUDE)

OCT_INCS=$(shell mkoctfile -p ALL_CXXFLAGS)
OCT_LIB_DIR=$(shell mkoctfile -p LFLAGS)
OCT_LIBS=$(shell mkoctfile -p OCTAVE_LIBS)

HEADERS=ExternalDataSet.hxx FSIFitFCN.hxx FSIParameterScan.hxx FSIChi2Grid.hxx
SOURCES=ExternalDataSet.cxx FSIFitFCN.cxx FSIParameterScan.cxx FSIChi2Grid.cxx
#DICT=dict.cc
LIBPISCAT=libPiScat.so

MAIN=PiScatFit2016.cxx
OBJ=$(MAIN:.cc=.o)
EXECUTABLE=PiScatFit2016

all: $(EXECUTABLE) $(LIBPISCAT)

# $(EXECUTABLE): $(MAIN) $(LIBPISCAT)
# 	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) -DVERBOSE

PiScatFit2016: PiScatFit2016.cxx $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

$(LIBPISCAT): $(SOURCES) $(DICT)
	$(CXX) $(SOFLAGS) $(OCT_INCS) -o $@ $(LDFLAGS) $(LDLIBS) $(OCT_LIB_DIR) $(OCT_LIBS) $^

#$(DICT): $(HEADERS)
#	$(ROO) -f $@ $(CFLAGS) -p $^

clean:
	rm -f PiScatFit2016.exe *.o *.so dict.*
