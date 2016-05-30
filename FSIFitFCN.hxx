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

  FSIFitFCN(FSIChi2Grid *aFSIChi2Grid): fFSIChi2Grid(aFSIChi2Grid), kOctave(false), kMultiDimFit(false){};
  void UseOctaveInterpolation(){ kOctave = true; kMultiDimFit = false;};
  void UseTMultiDimFitInterpolation() {kMultiDimFit = true; kOctave = false;};
  double operator() (const std::vector<double> & x) const;
  double Up() const {return 1.;}

 private:
  
  FSIChi2Grid *fFSIChi2Grid;
  bool kOctave;
  bool kMultiDimFit;
  
};
#endif
