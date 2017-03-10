#ifndef FSIFITTERCONFIG_HXX_
#define FSIFITTERCONFIG_HXX_

#include <iostream>
#include "TString.h"

const TString parNames[7] = {"FEFQE","FEFABS","FEFCX","FEFINEL","FEFQEH","FEFCXH","FEFALL"};
const Int_t nFSIpars = 7;
const Double_t nomPars[7] = {0.9,1.25,0.8,1.0,1.8,1.8,1.0};

// Larger grid corresponding to scan_all_c_o_al_fe_cu_pb_170318.root
const TString ScanFileName = "input/scan_all_c_o_al_fe_cu_pb_170318.root";
const Int_t nVariations = 17*17*16*13*11;
const Int_t nSteps[7] = {17,17,16,13,11,7,7};
const Double_t FSIParsMin[7]  = {0.10,0.35,0.10,0.20,0.80,1.20,0.70};
const Double_t FSIParsMax[7]  = {1.70,1.951,1.60,2.60,2.80,2.40,1.30}; //1.551 is a trick
const Double_t FSIParsStep[7] = {0.10,0.10,0.10,0.20,0.20,0.20,0.10};

// // Smaller grid corresponding to scan_all_c_o_al_fe_cu_pb_extra_c.root
// const TString ScanFileName = "input/scan_all_c_o_al_fe_cu_pb_extra_c.root";
// const Int_t nVariations =7*7*7*7*7;
// const int nSteps[7] = {7,7,7,7,7,7,7};
// const Double_t FSIParsMin[7]  = {0.60,0.95,0.50,0.40,1.20,1.20,0.70};
// const Double_t FSIParsMax[7]  = {1.20,1.551,1.10,1.60,2.40,2.40,1.30}; //1.551 is a trick
// const Double_t FSIParsStep[7] = {0.10,0.10,0.10,0.20,0.20,0.20,0.10};

// Grid actually being fitted (removed FEFCXH and FEFALL)
const int nFSIparsFitted = 5;
const Int_t nMaxPointsUsed = 190; // current max size with all nuclei

const std::vector<TString> IndivXsecs = {"xqe","xabs","xcx","xdcx","xhadr"};
const std::vector<TString> xsName = {"REAC","QE","ABS","CX","ABS+CX"};
const std::vector<TString> iTypeString = {"reac","inel","abs","cx","abscx"};//"dcx"
const std::vector<TString> iTypeString2 = {"rea","ine","abs","scx","abscx"};
const std::vector<TString> iTypeString3 = {"Reactive","Quasi-Elastic","Absorption (ABS)","Single Charge Exchange (CX)","ABS+CX"};//"dcx"

const std::vector<Int_t> NucleiA = {6,8,13,26,29,82};
const std::vector<TString> Nuclei = {"c","o","al","fe","cu","pb"};
const std::vector<TString> Nuclei2 = {"c12","o","al","fe","cu","pb"};
const std::vector<TString> pionLatex = {"#pi^{+}","#pi^{-}"};
const std::vector<Int_t> pids = {211,-211};  
const std::vector<TString> pionType = {"piP","piM"};
const std::vector<TString> pionType2 = {"pip","pim"};
const std::vector<int> pionCode = {211,-211};

//////////////////////////////////////////////////////////////////////////
// This is a placeholder of a class that would read in the config above //
// To be implemented in the future (maybe)
//////////////////////////////////////////////////////////////////////////
struct ConfigParams {
  
  const int param1;
  const int param2;
  // ...
  
  // Define a constructor that will load stuff from a configuration file.
  // Tune this to meet your needs, such as parsing the data from a startup command line.
  ConfigParams(const std::string & configFileName);
  int getConfigFromFile(const std::string & configFileName, const std::string & key);
  
};

  // The global that you will use to fetch the configurations during runtime.
extern const ConfigParams configParams;

#endif
