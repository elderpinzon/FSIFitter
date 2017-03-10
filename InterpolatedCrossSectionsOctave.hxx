#ifndef INTERPOLATEDCROSSSECTIONSOCTAVE_HXX
#define INTERPOLATEDCROSSSECTIONSOCTAVE_HXX

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

// Octave libraries
#include <octave/config.h>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>

#include "TTree.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TMultiDimFit.h"
#include "TPaveText.h"

#include "FSIParameterScan.hxx"
#include "ExternalDataSet.hxx"
#include "FSIFitterUtils.hxx"
#include "FSIFitterConfig.hxx"

class InterpolatedCrossSectionsOctave{
  
private:

  TString filename;
  Int_t thisTarget;
  Int_t thisPID;
  std::vector<Double_t> allMoms;

  int nPars;
  int nStepsOctave[6];

  // Octave objects
  Matrix *x_vec[nFSIpars];  // Axes points to build mesh
  Array<double> *v_reac; // Mesh of cross-sections
  Array<double> *v_qe; // Mesh of cross-sections
  Array<double> *v_abs; // Mesh of cross-sections
  Array<double> *v_cx; // Mesh of cross-sections
  Array<double> *v_abscx; // Mesh of cross-sections


public:

  InterpolatedCrossSectionsOctave();
  InterpolatedCrossSectionsOctave(TString fFileName,Int_t target,Int_t pid);
  ~InterpolatedCrossSectionsOctave();
  void BuildGrid();
  void BuildMomentumIndices();
  Double_t GetSplinedGridPoint(Int_t thisXSec,const std::vector<double> &x);
  Int_t CarbonMomToIndex(Double_t thisMom);
  Int_t MomToIndex(Double_t thisMom);
  std::vector<Double_t> GetAllMoms(){return allMoms;};
  Int_t GetNumberMomenta(){return allMoms.size();}
};

#endif  
