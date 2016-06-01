//***********************************************************************************
// piScatFCN.cxx
// Minimization function for TMinuit2
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "FSIFitFCN.hxx"

FSIFitFCN::FSIFitFCN(FSIChi2Grid *aFSIChi2Grid, bool fnuclFit, bool fOctave, bool fMultiDimFit): fFSIChi2Grid(aFSIChi2Grid), kOctave(false), kMultiDimFit(false){

  nuclFit=fnuclFit;
  kOctave=fOctave;
  kMultiDimFit=fMultiDimFit;
}

double FSIFitFCN::operator() (const std::vector<double> &x) const {

  for(int i = 0; i < (int)x.size(); i++)
    std::cout << "\tx[" << i <<"]: "<<x[i];
  std::cout << std::endl;

  // Simply get the result of the desired interpolation evaluated at this grid point
  if(kMultiDimFit) {
  	 // std::cout << "MultiDimFit fcn" << std::endl;
  	return fFSIChi2Grid->GetInterpolatedGridPoint(x);
  }else if(kOctave) {
  	// std::cout << "octave fcn" << std::endl;
  	return fFSIChi2Grid->GetSplinedGridPoint(x);
  }else{
    std::cout << "FSIFitFCN: Must select either Octave or TMultiDimFit!" << std::endl;
    std::exit(-1);
  }

}

