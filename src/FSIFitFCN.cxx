//***********************************************************************************
// piScatFCN.cxx
// Minimization function for TMinuit2
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "FSIFitFCN.hxx"

double FSIFitFCN::operator() (const std::vector<double> &x) const {

  for(int i = 0; i < (int)x.size(); i++)
    std::cout << "\tx[" << i <<"]: "<<x[i];
  std::cout << std::endl;

  // Simply get the result of the desired interpolation evaluated at this grid point
  if(kMultiDimFit)  return fFSIChi2Grid->GetInterpolatedGridPoint(x);
  else if(kOctave)  return fFSIChi2Grid->GetSplinedGridPoint(x);
  else{
    std::cout << "FSIFitFCN: Must select either Octave or TMultiDimFit!" << std::endl;
    std::exit(-1);
  }

}

