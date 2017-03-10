#ifndef FSICHI2GRIDNORM_HXX
#define FSICHI2GRIDNORM_HXX

#include <iostream>
#include <string>
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
#include "ModelPrediction.hxx"
#include "InterpolatedCrossSectionsOctave.hxx"

class FSIChi2GridNorm{
  
private:

  // Vector of datasets to be used in the chisquare calculation
  std::vector< ExternalDataSet::ExternalDataSet > fAllDataSets;
  std::vector< ExternalDataSet::ExternalDataSet > fNonDUETDataSets;
  int nPointsUsed;

  TMatrixD fullCovariance;
  TMatrixD fullCovarianceInv;
  TMatrixD fullCorrelation;
  TVectorD allErrors;
  TVectorD allMomenta;
  TVectorD allSigmaData;
  TVectorD allTargets;
  std::vector<TString> allTargetsString;
  TVectorD allPids;
  TVectorD allTypes;
  
  // Name of the file with the result of NEUT grid scan
  TString FileName;

  // A tree for storage/reading and its necessary variables
  TFile *foutTree;
  TTree *ChiSquareFiniteGrid;
  Double_t fqe;
  Double_t fabs;
  Double_t fcx;
  Double_t fFSIPars[nFSIpars];
  Double_t fchisquare;
  Double_t freac;
  Double_t felas;
  Double_t ftot;

  // The TMultiDimFit object
  TMultiDimFit *multifit;

  // Octave objects
  Matrix *x_vec[nFSIpars];  // Axes points to build mesh
  Array<double> *v_vec; // Mesh of chi^2

  // Interpolated cross setions
  InterpolatedCrossSectionsOctave *aInterpolatedCrossSections;
  InterpolatedCrossSectionsOctave *aInterpolatedCrossSection[6][2];
  

public:

  FSIChi2GridNorm(){};
  ~FSIChi2GridNorm(){};
  FSIChi2GridNorm(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets);

  // Getters and Setters
  void SetGridScanFileName(TString fileName){FileName = fileName;};
  Double_t GrabInterpolatedCrossSectionOctave(Int_t thisNuclei, Int_t thisPid, Int_t thisType, std::vector<double> fsipars);
  Double_t GetChi2(std::vector<Double_t> fsipars,std::vector<Double_t> normpars);

  void LoadInterpolatedCrossSections();
  
  bool nuclFit;
  bool kOctave;
  bool kMultiDimFit;

};

#endif  
