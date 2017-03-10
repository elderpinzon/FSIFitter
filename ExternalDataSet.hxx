#ifndef EXTERNALDATASET_HXX
#define EXTERNALDATASET_HXX

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
#include "FSIFitterUtils.hxx"

const int maxDataPoints = 13;

class ExternalDataSet{
  
private:
  
  TString FileName;
  TString SetName;
  TString NucleiString;
  Int_t Nuclei;
  TString pionTypeString;
  Int_t pionType;
  TString nucleonType;
  TString intrTypeString;
  Int_t intrType;
  Int_t nDataPoints;

  TVectorD vec_data_Mom;
  TVectorD vec_data_Xsec;
  TVectorD vec_data_Xsec_error;
  TGraphErrors *fDataPointsGraph;
  TMatrixD m_data_cov;
  TMatrixD m_data_cor;
  TMatrixD m_data_cov_inv;

  
public:

  ExternalDataSet();
  ExternalDataSet(TString fFileName);
  void Initialize();
  void ParseDataSetName();
  void LoadData(); //Load data into TGraph and TVectors
  void BuildDataDiagonalMatrices(); //Only for non-DUET data
  void LoadDUETCovariance(); // Only for DUET separate measurement

  // Getters
  TGraphErrors* GetTGraph(){ return fDataPointsGraph;};
  TVectorD GetMomVector(){ return vec_data_Mom;};
  TVectorD GetXsecVector(){ return vec_data_Xsec;};
  TVectorD GetXsecErrorVector(){ return vec_data_Xsec_error;};
  TMatrixD GetCovariance(){return m_data_cov;};
  TMatrixD GetCovarianceInv(){return m_data_cov_inv;};
  Int_t GetIntrType(){ return intrType;};
  TString GetIntrTypeString(){ return intrTypeString;};
  Int_t GetPionType(){ return pionType;};
  TString GetPionTypeString(){ return pionTypeString;};
  Int_t GetNucleus(){ return Nuclei;};
  TString GetNucleusString(){ return NucleiString;};
  TString GetSetName(){ return SetName;};
  TString GetFileName(){ return FileName;};
  Int_t GetNumberOfDataPoints(){ return nDataPoints;};

};

#endif  
