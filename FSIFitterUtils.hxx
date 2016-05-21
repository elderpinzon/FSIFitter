#ifndef FSIFITTERUTILS_HXX_
#define FSIFITTERUTILS_HXX_

#include <iostream>

#include <TMultiDimFit.h>
#include <TGraphErrors.h>

namespace FSIFitterUtils{

  double BuildInterpolatedFunctionAndEvaluate(TMultiDimFit *aTMultiDimFit, const std::vector<double> &x);

  TGraphErrors* MergeGraphs(TCollection* li);
  TGraphErrors* MergeGraphsIntoEnvelope(TGraphErrors *gr_min, TGraphErrors *gr_max);

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
