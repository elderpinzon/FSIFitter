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
  std::string SetName;
  std::string Nuclei;
  std::string pionType;
  int intrType;
  int nDataPoints;
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
  ExternalDataSet(std::string fFileName);
  //ExternalDataSet& operator=(const ExternalDataSet&);
  void ParseDataSetName();
  int IntrNameToInt(std::string fIntrType);
  std::string GetSetName();
  void LoadData(); //Load data into TGraph and TVectors
  void BuildDataDiagonalMatrices(); //Only for non-DUET data
  void Initialize(); 
  int GetNumberOfDataPoints();
  double CalculateChiSquare(Double_t this_abs, Double_t this_cx, Double_t this_qe, FSIParameterScan &aFSIParameterScan);
  double CalculateChiSquare(Double_t this_abs, Double_t this_cx, Double_t this_qe);
  double GetChiSquare();
  void PrintParameterSet(double f_qe, double f_abs, double f_cx);
  void FillMiniTreeFromMCScan();
  
};

#endif  
