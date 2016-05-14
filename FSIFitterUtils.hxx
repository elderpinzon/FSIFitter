#ifndef FSIFITTERUTILS_HXX_
#define FSIFITTERUTILS_HXX_

#include <TMultiDimFit.h>

namespace FSIFitterUtils{

  double BuildInterpolatedFunctionAndEvaluate(TMultiDimFit *aTMultiDimFit, const std::vector<double> &x);

};

#endif
