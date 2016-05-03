#ifndef NUCLCHI2GRID_HXX
#define NUCLCHI2GRID_HXX

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

#include "NuclFSIParameterScan.h"
#include "NuclExternalDataSet.hxx"

class NuclFSIChi2Grid{
  
private:

  // Vector of datasets to be used in the chisquare calculation
  std::vector< NuclExternalDataSet::NuclExternalDataSet > fAllDataSets;

  // NEUT MC parameter scan to be used
  //FSIParameterScan fFSIParameterScan;

  // A tree for storage/reading and its necessary variables
  TTree *ChiSquareFiniteGrid;
  Double_t freac;
  Double_t felas;
  Double_t ftot;
  Double_t fchisquare;

  // The TMultiDimFit object
  TMultiDimFit *multifit;

public:

  NuclFSIChi2Grid(){};
  ~NuclFSIChi2Grid(){};
  NuclFSIChi2Grid(NuclExternalDataSet SingleDataSet);
  NuclFSIChi2Grid(std::vector< NuclExternalDataSet::NuclExternalDataSet > AllDataSets);  
  //NuclFSIChi2Grid(std::vector< NuclExternalDataSet::NuclExternalDataSet > AllDataSets, FSIParameterScan &aFSIParameterScan);
  NuclFSIChi2Grid(TTree &BuiltGridTree);
  void ImportDataSets(NuclExternalDataSet AllDataSets);
  void BuildFiniteGrid();
  TTree* GetFiniteGrid();
  void InterpolateFiniteGrid(bool kUseInputGrid);
  void UseExistingGrid();
  double GetInterpolatedGridPoint(const std::vector<double> &x);
};

#endif  
