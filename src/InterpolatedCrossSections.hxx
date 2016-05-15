#ifndef INTERPOLATEDCROSSSECTIONS_HXX
#define INTERPOLATEDCROSSSECTIONS_HXX

#include <string>

#include "TMultiDimFit.h"
#include "FSIParameterScan.hxx"
#include "FSIFitterUtils.hxx"

std::vector< std::string > xsecs = {"xqe","xabs","xcx","xdcx","xhadr"};

class InterpolatedCrossSections{

private:

  std::string filename;
  // One for each cross section (qe, abs, cx, dcx, hadr)
  TMultiDimFit *xsecmultidimfit[5];

public:

  InterpolatedCrossSections();
  InterpolatedCrossSections(std::string fFileName);
  void RunMultiDimFitInterpolations();
  double GetInterpolatedCrossSection(int xs, const std::vector<double> &x);

};

#endif
