#include <iostream>
#include <cmath>
#include <initializer_list>

#include "TVirtualFitter.h"
#include "TStyle.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TMinuit.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TStopwatch.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
// #include "Math/Util.h"


#include "ExternalDataSet.hxx"
#include "FSIParameterScan.hxx"
#include "FSIFitFCN.hxx"
#include "FSIFitFCNNorm.hxx"
#include "FSIChi2Grid.hxx"
#include "FSIChi2GridNorm.hxx"
#include "FSIFitterUtils.hxx"
#include "ThrowParms.hxx"
#include "TN032Envelopes.hxx"
#include "ModelPrediction.hxx"
#include "InterpolatedCrossSectionsOctave.hxx"

// Running options
bool kBuildGridFromScratch = true; //default option
bool kUseInputGrid = false;
bool kUseInputMultiDimFit = false;
bool nuclFit=0;
bool kOctave = false;
int kPiece = -1;
int kTotalPieces = -1;
int dataFit = 0; // 0: {C}, 1: {C,Al,O}, 2:{C,O,Al,Fe,Cu,Pb}
int nFitPars = 6;
bool fitNorm = false;
Double_t scalingFactor = 1.0;

// Input file with MC scan
std::string inputfileScan="";

// Input file with existing grid
std::string inputfileGrid="";

// Output file with fit results
std::string outputfileFit ="output/fit_result.root";

// Vector of data sets for easier manipulation
std::vector< ExternalDataSet::ExternalDataSet > AllDataSets;
std::vector< ExternalDataSet::ExternalDataSet > NonDUETDataSets;
int nDOF = 0;

// FSIChi2Grid object
FSIChi2Grid::FSIChi2Grid *aFSIChi2Grid;
FSIChi2GridNorm *aFSIChi2GridNorm;

// Clock for time-keeping
TStopwatch stopwatch;
TStopwatch stopwatch_2;

//***********************************************************************************
// Define list of data sets to be added and compile them as ExternalDataSets
// objects in a vector to be used later
//***********************************************************************************
void AddDataSets(){
  std::vector< std::string > DataSetsToBeAdded;
  if (nuclFit){
    std::cout<<"Nucleon Fit data---"<<std::endl;
    DataSetsToBeAdded = {
              // "c_reac_p_ch55",
              "c_reac_p_mi54",
              "c_reac_p_re72",
              "c_tot_p_ca54",
              "c_tot_p_ja78",
              "c_tot_p_ma53",
              "c_tot_p_sc79",
              "c_elas_p_147",
              "c_elas_p_131",
              "c_elas_p_58",
              "c_elas_p_33",
              "c_elas_p_126",
              "c_elas_p_24",
              "c_elas_p_102",
              "c_elas_p_124"
            };
  }
}

//***********************************************************************************
// Alternative (more automatic) way of adding datasets
//***********************************************************************************
void AddDataSets(TString Nuclei){

  TSystemDirectory dir("pion_data", "pion_data");
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".csv") && fname.BeginsWith(Nuclei)) {
	fname = fname.ReplaceAll(".csv","");
	std::cout << "\nDataset to be added: " << fname << std::endl;
	ExternalDataSet::ExternalDataSet fDataSet(fname);
	AllDataSets.push_back(fDataSet);
	if(!(fname.Contains("DUET"))) NonDUETDataSets.push_back(fDataSet);
	else nDOF += 10;
      }
    }
  }
}
