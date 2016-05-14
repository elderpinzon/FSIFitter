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
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>

#include "TTree.h"
#include "TVector.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include <iostream>
//#include <TMultiDimFit.h>
#include <TMultiDimFit.h>

#include "FSIParameterScan.hxx"
#include "ExternalDataSet.hxx"
#include "FSIFitterUtils.hxx"

// Temporarily setting grid size here - need better way of doing it!
const int nFSIpars = 3;

class FSIChi2Grid{
  
private:

  // Vector of datasets to be used in the chisquare calculation
  std::vector< ExternalDataSet::ExternalDataSet > fAllDataSets;

  // NEUT MC parameter scan to be used
  //FSIParameterScan fFSIParameterScan;

  // A tree for storage/reading and its necessary variables
  TTree *ChiSquareFiniteGrid;
  Double_t fqe;
  Double_t fabs;
  Double_t fcx;
  Double_t fchisquare;

  // The TMultiDimFit object
  TMultiDimFit *multifit;

  // Octave objects
  Matrix *x_vec[nFSIpars];  // Axes points to build mesh
  Array<double> *v_vec; // Mesh of chi^2


public:

  FSIChi2Grid(){};
  ~FSIChi2Grid(){};
  FSIChi2Grid(ExternalDataSet SingleDataSet);
  FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets);  
  //FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets, FSIParameterScan &aFSIParameterScan);
  FSIChi2Grid(TTree &BuiltGridTree);
  void ImportDataSets(ExternalDataSet AllDataSets);
  void BuildFiniteGrid();
  TTree* GetFiniteGrid();
  void InterpolateFiniteGrid(bool kUseInputGrid);
  void UseExistingGrid();
  double GetInterpolatedGridPoint(const std::vector<double> &x);
  double GetSplinedGridPoint(const std::vector<double> &x);
  
};

#endif  
