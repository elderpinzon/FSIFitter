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


#include "NuclExternalDataSet.hxx"
#include "NuclFSIParameterScan.h"
#include "NuclFSIFitFCN.hxx"
#include "NuclFSIChi2Grid.hxx"

// Running options
bool kUseInputGrid = false;
bool kBuildGridFromScratch = false;

// Input file with MC scan
std::string inputfileScan="";

// Input file with existing grid
std::string inputfileGrid="";

// Output file with fit results
std::string outputfileFit ="output/fit_result_nucl.root";

// Vector of data sets for easier manipulation
std::vector< NuclExternalDataSet::NuclExternalDataSet > AllDataSets;

// NuclFSIChi2Grid object
NuclFSIChi2Grid::NuclFSIChi2Grid aFSIChi2Grid;
