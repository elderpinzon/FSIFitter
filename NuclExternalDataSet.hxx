#ifndef NUCLEXTERNALDATASET_HXX
#define NUCLEXTERNALDATASET_HXX

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

#include "TTree.h"
#include "TVector.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"

#include "NuclFSIParameterScan.h"

class NuclExternalDataSet{
  
private:
  
  std::string FileName;
  std::string SetName;
  std::string Nuclei;
  std::string nucleonType;
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

  NuclExternalDataSet();
  NuclExternalDataSet(std::string fFileName);
  NuclExternalDataSet& operator=(const NuclExternalDataSet&);
  void ParseDataSetName();
  int IntrNameToInt(std::string fIntrType);
  std::string GetSetName();
  void LoadData(); //Load data into TGraph and TVectors
  void BuildDataDiagonalMatrices(); //Only for non-DUET data
  void Initialize(); 
  int GetNumberOfDataPoints();
  double CalculateChiSquare(Double_t this_reac, Double_t this_elas, Double_t this_tot, NuclFSIParameterScan &aNuclFSIParameterScan);
  double CalculateChiSquare(Double_t this_reac, Double_t this_elas, Double_t this_tot);
  double GetChiSquare();
  void PrintParameterSet(double f_reac, double f_elas, double f_tot);
  void FillMiniTreeFromMCScan();
  
};

#endif  
