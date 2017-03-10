//***********************************************************************************
// piScatFCN.cxx
// Minimization function for TMinuit2
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "FSIFitFCN.hxx"

FSIFitFCN::FSIFitFCN(FSIChi2Grid *aFSIChi2Grid): fFSIChi2Grid(aFSIChi2Grid){
  
  scalingFactor = 1.0;
  
}

double FSIFitFCN::operator() (const std::vector<double> &x) const {
  
  // Print the parameters to follow progression of the minimization
  for(int i = 0; i < (int)x.size(); i++)
    printf(" %4.4f",x[i]);
  
  double value = -1;
  
  // Simply get the result of the desired interpolation evaluated at this grid point
  if(!kOctave) {
    value = fFSIChi2Grid->GetInterpolatedGridPoint(x); 
    std::cout << " MultiDimFit fcn: " << value << std::endl;
  }else if(kOctave) {
    value = fFSIChi2Grid->GetSplinedGridPoint(x); 
    std::cout << " octave fcn: " << value << std::endl;
  }else{
    std::cout << "FSIFitFCN: Must select either Octave or TMultiDimFit!" << std::endl;
    std::exit(-1);
  }
  
  // Apply scaling factor if needed
  value /= scalingFactor;

  return value;

}

// Overload function to make this work with new Minuit
double FSIFitFCN::operator()( const double *  par ) const { 
  std::vector<double> p(par, par+_NPARS); 
  return (*this)(p); 
}
