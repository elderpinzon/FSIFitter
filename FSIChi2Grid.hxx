#ifndef PIONCHI2GRID_HXX
#define PIONCHI2GRID_HXX

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

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
};

#endif  
