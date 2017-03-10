#ifndef _PISCAT_FCN_H
#define _PISCAT_FCN_H

#include <iostream>
#include <vector>
#include <fstream>

// Minuit headers
#include "TVirtualFitter.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "Math/IFunction.h"

#include "FSIChi2Grid.hxx"

class FSIFitFCN : public ROOT::Minuit2::FCNBase {
  
 public:

  FSIFitFCN(){};
  ~FSIFitFCN(){};
  FSIFitFCN(FSIChi2Grid *aFSIChi2Grid);
  void SetOctaveInterpolation(){ kOctave = true; kMultiDimFit = false;};
  void SetMultiDimFitInterpolation() {kMultiDimFit = true; kOctave = false;};
  double operator() (const std::vector<double> & x) const;
  double operator() (const double * par) const;
  double Up() const {return 1.;}
  void SetScalingFactor(Double_t factor){scalingFactor = factor;};
  void SetNPARS(Int_t NPARS){_NPARS = NPARS;};

 private:
  
  FSIChi2Grid *fFSIChi2Grid;
  Double_t scalingFactor;
  Int_t _NPARS;
  bool kOctave;
  bool kMultiDimFit;
  
};
#endif
