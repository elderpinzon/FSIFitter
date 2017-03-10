//***********************************************************************************
// piScatFCN.cxx
// Minimization function for TMinuit2
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "FSIFitFCNNorm.hxx"

FSIFitFCNNorm::FSIFitFCNNorm(FSIChi2GridNorm *aFSIChi2GridNorm): fFSIChi2GridNorm(aFSIChi2GridNorm){}

double FSIFitFCNNorm::operator() (const std::vector<double> &x) const {

  std::cout << "FSIFitFCNNorm:: called with " 
	    << x.size() << " parameters" << std::endl;
  
  // Separate FSI and normalization
  std::vector<Double_t> fsi;
  std::vector<Double_t> norm;
  for(int i = 0; i < (int)x.size(); i++){
    
    if(i<nFSIparsFitted) fsi.push_back(x[i]);
    else norm.push_back(x[i]);
    
    printf(" %4.4f\n",x[i]);
    
  }
  
  return fFSIChi2GridNorm->GetChi2(fsi,norm);
  
}

// Overload function to make this work with new Minuit
double FSIFitFCNNorm::operator()( const double *  par ) const { 
  std::vector<double> p(par, par+_NPARS); 
  return (*this)(p); 
}
