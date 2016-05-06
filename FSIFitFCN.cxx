//***********************************************************************************
// piScatFCN.cxx
// Minimization function for TMinuit2
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "FSIFitFCN.hxx"

double FSIFitFCN::operator() (const std::vector<double> &x) const {

  // Simply get the result of the interpolation evaluated at this grid point
  double inter_chi2 = fFSIChi2Grid->GetInterpolatedGridPoint(x);
  double splined_chi2 = fFSIChi2Grid->GetSplinedGridPoint(x);

  for(int i = 0; i < (int)x.size(); i++)
    std::cout << "\tx[" << i <<"]: "<<x[i];

  std::cout << "\tinter_chi2: " << inter_chi2 << " \t splined_chi2: " << splined_chi2 << std::endl;
  
  //return inter_chi2;
  return splined_chi2;
  
}

