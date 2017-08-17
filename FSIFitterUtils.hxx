#ifndef FSIFITTERUTILS_HXX_
#define FSIFITTERUTILS_HXX_

#include <iostream>
#include <vector>

#include <TMultiDimFit.h>
#include <TGraphErrors.h>
#include "TMath.h"

#include "FSIFitterConfig.hxx"

namespace FSIFitterUtils{

  double BuildInterpolatedFunctionAndEvaluate(TMultiDimFit *aTMultiDimFit, const std::vector<double> &x0);

  TGraphErrors* MergeGraphs(TCollection* li);
  TGraphErrors* MergeGraphsIntoEnvelope(TGraphErrors *gr_min, TGraphErrors *gr_max);
  void PrintParameterSet(double FSIPars[nFSIpars]);

  TMatrixD CovarianceToCorrelation(TMatrixD* h); 
  void MergeMatrices(TMatrixD &matrix1, TMatrixD matrix2);
  Int_t IntrNameToInt(TString fIntrType);
  Int_t NucleiNameToInt(TString fNuclei);
  Int_t PionTypeToInt(TString fPionType);
  Int_t TurnParamsToIndex(Double_t *kindex);
  void TurnIndexToParams(Int_t index,Double_t *kindex);
};

//***********************************************************************************
// Class to create a TF1 out of TMultiDimFits. Useful for plotting
//***********************************************************************************
class TF1FromMultiDimFits{
  
private:
  
  std::vector <TMultiDimFit*> *fVectorMultiDimFits;
  std::vector<double> fFSIPars;

public:
  
  TF1FromMultiDimFits(std::vector<TMultiDimFit*> *vectorMultiDimFits, std::vector<double> FSIPars): fVectorMultiDimFits(vectorMultiDimFits),fFSIPars(FSIPars) {};;
  double operator()(double *x, double *) const;

};


#endif
