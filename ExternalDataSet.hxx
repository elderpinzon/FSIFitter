#ifndef PIONEXTERNALDATASET_HXX
#define PIONEXTERNALDATASET_HXX

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

#include "TTree.h"
#include "TVector.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"

#include "FSIParameterScan.hxx"

class ExternalDataSet{
  
private:
  
  std::string FileName;
  TString SetName;
  TString Nuclei;
  TString pionType;
  TString intrTypeString;
  Int_t intrType;
  Int_t nDataPoints;
  TTree *MiniMCScan; //A TTree with only the MC scans relevant for this data set
  TVectorD vec_data_Mom;
  TVectorD vec_data_Xsec;
  TVectorD vec_data_Xsec_error;
  TGraphErrors *fDataPointsGraph;
  TMatrixDSym m_data_cov;
  TMatrixDSym m_data_cor;
  TMatrixDSym m_data_cov_inv;
  double fDataSetChiSquare;

  
public:

  ExternalDataSet();
  ExternalDataSet(TString fFileName);
  //ExternalDataSet& operator=(const ExternalDataSet&);
  void Initialize(); 
  void ParseDataSetName();
  void LoadData(); //Load data into TGraph and TVectors
  void BuildDataDiagonalMatrices(); //Only for non-DUET data
  int IntrNameToInt(TString fIntrType);  
  TGraphErrors* GetTGraph(){ return fDataPointsGraph;};
  TVectorD GetMomVector(){ return vec_data_Mom;};
  Int_t GetIntrType(){ return intrType;};
  TString GetIntrTypeString(){ return intrTypeString;};
  TString GetPionType(){ return pionType;};
  TString GetNucleus(){ return Nuclei;};
  TString GetSetName(){ return SetName;};
  Int_t GetNumberOfDataPoints(){ return nDataPoints;};
  double CalculateChiSquare(Double_t this_abs, Double_t this_cx, Double_t this_qe, FSIParameterScan &aFSIParameterScan);
  double CalculateChiSquare(Double_t this_abs, Double_t this_cx, Double_t this_qe);
  double GetChiSquare(){ return fDataSetChiSquare;};
  void PrintParameterSet(double f_qe, double f_abs, double f_cx);
  void FillMiniTreeFromMCScan();
};

#endif  
