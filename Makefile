CXX=$(shell root-config --cxx)
ROO=rootcint

ROOTFLAGS=$(shell root-config --cflags)
CFLAGS=-c -Wall -v -g -std=gnu++0x -O2 $(ROOTFLAGS)

LDFLAGS=-g -O2 $(shell root-config --ldflags) -DVERBOSE
LDLIBS=$(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -std=gnu++0x

INCLUDE=-I$(shell root-config --incdir) -I./
SOFLAGS=-fPIC -shared $(INCLUDE)

HEADERS=ExternalDataSet.hxx FSIFitFCN.hxx FSIParameterScan.hxx FSIChi2Grid.hxx NuclExternalDataSet.hxx NuclFSIFitFCN.hxx NuclFSIParameterScan.hxx NuclFSIChi2Grid.hxx #DUETDataSet.hxx
SOURCES=ExternalDataSet.cxx FSIFitFCN.cxx FSIParameterScan.cxx FSIChi2Grid.cxx NuclExternalDataSet.cxx NuclFSIFitFCN.cxx NuclFSIParameterScan.cxx NuclFSIChi2Grid.cxx #DUETDataSet.cxx
#DICT=dict.cc
LIBPISCAT=libPiScat.so
LIBNUCLSCAT=libNuclScat.so

MAIN=PiScatFit2016.cxx NuclScatFit2016.cxx
OBJ=$(MAIN:.cc=.o)
EXECUTABLE=PiScatFit2016 NuclScatFit2016

all: $(EXECUTABLE) $(LIBPISCAT) $(LIBNUCLSCAT)

# $(EXECUTABLE): $(MAIN) $(LIBPISCAT)
# 	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) -DVERBOSE

PiScatFit2016: PiScatFit2016.cxx $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) -DVERBOSE

NuclScatFit2016: NuclScatFit2016.cxx $(LIBNUCLSCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) -DVERBOSE


$(LIBPISCAT): $(SOURCES) $(DICT)
	$(CXX) $(SOFLAGS) -o $@ $(LDFLAGS) $(LDLIBS) $^

$(LIBNUCLSCAT): $(SOURCES) $(DICT)
	$(CXX) $(SOFLAGS) -o $@ $(LDFLAGS) $(LDLIBS) $^

#$(DICT): $(HEADERS)
#	$(ROO) -f $@ $(CFLAGS) -p $^

clean:
	rm -f PiScatFit2016.exe NuclScatFit2016.exe *.o *.so dict.*
