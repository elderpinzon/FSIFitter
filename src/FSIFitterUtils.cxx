#include "FSIFitterUtils.hxx"

//***********************************************************************************
// Build resulting function from the results of given TMultiDim object and evaluate at x
// Most of this code was copied from TMultiDimFit::MakeRealCode()
//***********************************************************************************
double FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(TMultiDimFit *aTMultiDimFit, const std::vector<double> &x){
  
  // Grab all necessary pieces
  Int_t gNCoefficients = aTMultiDimFit->GetNCoefficients();
  Int_t gNVariables = aTMultiDimFit->GetNVariables();
  const Int_t *gPower = aTMultiDimFit->GetPowers();
  const Int_t *gPowerIndex = aTMultiDimFit->GetPowerIndex();
  const TVectorD *gCoefficient = aTMultiDimFit->GetCoefficients();
  const TVectorD *gXMax = aTMultiDimFit->GetMaxVariables();
  const TVectorD *gXMin = aTMultiDimFit->GetMinVariables();
  const double gDMean =  aTMultiDimFit->GetMeanQuantity();
  double returnValue = gDMean;
  
  // Build the polynomial and evaluate it
  int i = 0, j = 0, k = 0;
  for (i = 0; i < gNCoefficients ; i++) {
    // Evaluate the ith term in the expansion
    double term = (*gCoefficient)[i];
    for (j = 0; j < gNVariables; j++) {
      // Evaluate the polynomial in the jth variable.
      int power = gPower[gPowerIndex[i] * gNVariables + j];
      double p1 = 1, p2 = 0, p3 = 0, r = 0;
      double v =  1 + 2. / ((*gXMax)[j] - (*gXMin)[j]) * (x[j] - (*gXMax)[j]);
      // what is the power to use!
      switch(power) {
      case 1: r = 1; break;
      case 2: r = v; break;
      default:
	p2 = v;
	for (k = 3; k <= power; k++) {
	  p3 = p2 * v;
	  p1 = p2; p2 = p3;
	}
	r = p3;
      }
      // multiply this term by the poly in the jth var
      term *= r;
    }
    // Add this term to the final result
    returnValue += term;
  }
  
  return returnValue;

}
