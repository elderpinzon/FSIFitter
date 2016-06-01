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
  FSIFitFCN(FSIChi2Grid *aFSIChi2Grid, bool fnuclFit, bool fOctave, bool fMultiDimFit);
  void UseOctaveInterpolation(){ kOctave = true; kMultiDimFit = false;};
  void UseTMultiDimFitInterpolation() {kMultiDimFit = true; kOctave = false;};
  double operator() (const std::vector<double> & x) const;
  double Up() const {return 1.;}

 private:
  
  FSIChi2Grid *fFSIChi2Grid;
  bool nuclFit;
  bool kOctave;
  bool kMultiDimFit;
  
};
#endif
