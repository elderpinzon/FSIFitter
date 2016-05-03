#ifndef _NUCLSCAT_FCN_H
#define _NUCLSCAT_FCN_H

#include <iostream>
#include <vector>
#include <fstream>

// Minuit headers
#include "TVirtualFitter.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "Math/IFunction.h"

#include "NuclFSIChi2Grid.hxx"

class NuclFSIFitFCN : public ROOT::Minuit2::FCNBase {
  
 public:

  NuclFSIFitFCN(NuclFSIChi2Grid *aFSIChi2Grid): fNuclFSIChi2Grid(aFSIChi2Grid){};
  
  double operator() (const std::vector<double> & x) const;
  double Up() const {return 1.;}

 private:
  
  NuclFSIChi2Grid *fNuclFSIChi2Grid;
  
};
#endif
