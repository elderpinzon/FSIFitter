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
  double chi2 = fFSIChi2Grid->GetInterpolatedGridPoint(x);

  for(int i = 0; i < (int)x.size(); i++)
    std::cout << "\tx[" << i <<"]: "<<x[i];

  std::cout << "\tchi2: " << chi2 << std::endl;
  
  return chi2;
  
}

