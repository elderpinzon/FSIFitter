#ifndef FSICHI2GRID_HXX
#define FSICHI2GRID_HXX

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

class FSIChi2Grid{
  
private:

  // Vector of datasets to be used in the chisquare calculation
  std::vector< ExternalDataSet::ExternalDataSet > fAllDataSets;
  Int_t nPointsUsed;

  // Auxiliary vectors and matrices
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
  Double_t fFSIPars[nFSIparsFitted];
  Double_t fchisquare;
  Double_t freac;
  Double_t felas;
  Double_t ftot;

  // The TMultiDimFit object
  TMultiDimFit *multifit;

  // Octave objects
  Matrix *x_vec[nFSIpars];  // Axes points to build mesh
  Array<double> *v_vec; // Mesh of chi^2


public:

  FSIChi2Grid(){};
  ~FSIChi2Grid(){};
  FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets);

  // Getters and Setters
  void SetOctaveInterpolation(){kOctave=true;};
  void SetMultiDimFitInterpolation(){kMultiDimFit=true;};
  void SetGridScanFileName(TString fileName){FileName = fileName;};
  TTree* GetFiniteGrid(){return ChiSquareFiniteGrid;};
  Double_t GetInterpolatedGridPoint(const std::vector<Double_t> &x);
  Double_t GetSplinedGridPoint(const std::vector<Double_t> &x);


  // Functions to buld or run things
  void BuildCovariance();
  void FastBuildFiniteGrid();
  void BuildOctaveGrid();
  void BuildTMultiDimFitGrid();
  void RunTMultiDimFitInterpolation();
  void InterpolateFiniteGrid();
  void PlotChiSquare(const std::vector<Double_t> &x_center, int index);
  Double_t CalculateModelChiSquare(ModelPrediction model);  
  bool CheckInsideGridBoundary(const std::vector<Double_t> &x);
  
  bool nuclFit;
  bool kOctave;
  bool kMultiDimFit;

};

#endif  
