#ifndef _FSIFITFCNNORM_H
#define _FSIFITFCNNORM_H

#include <iostream>
#include <vector>
#include <fstream>

// Minuit headers
#include "TVirtualFitter.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "Math/IFunction.h"

#include "FSIChi2GridNorm.hxx"

class FSIFitFCNNorm : public ROOT::Minuit2::FCNBase {
  
 public:

  FSIFitFCNNorm(){};
  ~FSIFitFCNNorm(){};
  FSIFitFCNNorm(FSIChi2GridNorm *aFSIChi2Grid);
  double operator() (const std::vector<double> & x) const;
  double operator() (const double * par) const;
  double Up() const {return 1.;}
  void SetNPARS(Int_t NPARS){_NPARS = NPARS;};

 private:
  
  FSIChi2GridNorm *fFSIChi2GridNorm;
  Int_t _NPARS;
  
};
#endif
