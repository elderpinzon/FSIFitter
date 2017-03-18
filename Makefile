CXX=$(shell root-config --cxx)
ROO=rootcint

ROOTFLAGS=$(shell root-config --cflags)
CFLAGS=-c -Wall -v -g -std=gnu++0x -O2 $(ROOTFLAGS)

LDFLAGS=-g -O2 $(shell root-config --ldflags) -DVERBOSE
#LDLIBS=$(shell root-config --glibs) -lMinuit -lMathMore -lMinuit2 -std=gnu++0x
LDLIBS=$(shell root-config --glibs) -lMinuit -lMinuit2 -std=gnu++0x

INCLUDE=-I$(shell root-config --incdir) -I./
SOFLAGS=-fPIC -shared $(INCLUDE)

OCT_INCS=$(shell mkoctfile -p ALL_CXXFLAGS)
OCT_LIB_DIR=$(shell mkoctfile -p LFLAGS)
OCT_LIBS=$(shell mkoctfile -p OCTAVE_LIBS)

HEADERS=ExternalDataSet.hxx FSIFitFCN.hxx FSIParameterScan.hxx FSIChi2Grid.hxx FSIFitterUtils.hxx InterpolatedCrossSections.hxx ModelPrediction.hxx TN032Envelopes.hxx ThrowParms.hxx InterpolatedCrossSectionsOctave.hxx FSIChi2NotGridNorm.hxx FSIFitFCNNorm.hxx
SOURCES=ExternalDataSet.cxx FSIFitFCN.cxx FSIParameterScan.cxx FSIChi2Grid.cxx FSIFitterUtils.cxx InterpolatedCrossSections.cxx ModelPrediction.cxx TN032Envelopes.cxx ThrowParms.cxx InterpolatedCrossSectionsOctave.cxx FSIChi2GridNorm.cxx FSIFitFCNNorm.cxx

LIBPISCAT=libPiScat.so

MAIN=ScatFit.cxx
EXECUTABLE=ScatFit

MAIN7=ScatFitNorm.cxx
EXECUTABLE7=ScatFitNorm

MAIN2=ConstructFiniteGrid.cxx
EXECUTABLE2=ConstructFiniteGrid

MAIN3=RunCrossSectionInterpolation.cxx
EXECUTABLE3=RunCrossSectionInterpolation

MAIN6=RunCrossSectionInterpolationOctave.cxx
EXECUTABLE6=RunCrossSectionInterpolationOctave

MAIN4=PlotResult.cxx
EXECUTABLE4=PlotResult

MAIN5=MergeCovariances.cxx
EXECUTABLE5=MergeCovariances

all: $(EXECUTABLE) $(LIBPISCAT)

ScatFit: $(MAIN) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

ScatFitNorm: $(MAIN7) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

ConstructFiniteGrid: $(MAIN2) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

RunCrossSectionInterpolation: $(MAIN3) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

RunCrossSectionInterpolationOctave: $(MAIN6) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

PlotResult: $(MAIN4) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

MergeCovariances: $(MAIN5) $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOSE

PreComputeGridDUET: PreComputeGridDUET.cxx $(LIBPISCAT)
	$(CXX) -o $@.exe $^ $(LDFLAGS) $(INCLUDE) $(LDLIBS) $(OCT_INCS) $(OCT_LIB_DIR) $(OCT_LIBS) -DVERBOS

$(LIBPISCAT): $(SOURCES) $(DICT)
	$(CXX) $(SOFLAGS) $(OCT_INCS) -o $@ $(LDFLAGS) $(LDLIBS) $(OCT_LIB_DIR) $(OCT_LIBS) $^

clean:
	rm -f $(EXECUTABLE).exe $(EXECUTABLE2).exe $(EXECUTABLE3).exe $(EXECUTABLE4).exe *.o *.so dict.*
