#include <iostream>
#include <cmath>
#include <initializer_list>

#include "TVirtualFitter.h"
#include "TStyle.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/Util.h"


#include "ExternalDataSet.hxx"
#include "FSIParameterScan.h"
#include "FSIFitFCN.hxx"
#include "FSIChi2Grid.hxx"

// Running options
bool kUseInputGrid = false;
bool kBuildGridFromScratch = false;

// Input file with MC scan
std::string inputfileScan="";

// Input file with existing grid
std::string inputfileGrid="";

// Output file with fit results
std::string outputfileFit ="output/fit_result.root";

// Vector of data sets for easier manipulation
std::vector< ExternalDataSet::ExternalDataSet > AllDataSets;

// FSIChi2Grid object
FSIChi2Grid::FSIChi2Grid aFSIChi2Grid;
