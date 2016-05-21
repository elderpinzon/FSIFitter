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

//***********************************************************************************
// Merges a collection of TGraphErrors into a single ordered one. Useful for plotting
// TGraphErrors::Merge() seems to lose the error info when merging
//***********************************************************************************
TGraphErrors* FSIFitterUtils::MergeGraphs(TCollection* li){

  TGraphErrors *merged = new TGraphErrors();
  TIter next(li);

  while (TObject* o = next()) {

    TGraph *g = dynamic_cast<TGraph*>(o);
    if (!g) {
      std::cout << "Merge Cannot merge - an object which doesn't inherit from TGraph found in the list" << std::endl;
      exit(-1);
    }

    int n0 = merged->GetN();
    int n1 = n0+g->GetN();

    merged->Set(n1);
    Double_t * x = g->GetX();
    Double_t * y = g->GetY();
    Double_t * ex = g->GetEX();
    Double_t * ey = g->GetEY();
    for (Int_t i = 0 ; i < g->GetN(); i++) {
      //std::cout << "x[i] " << x[i] << "ex[i]: " << ex[i] << std::endl;
      merged->SetPoint(n0+i, x[i], y[i]);
      merged->SetPointError(n0+i, ex[i], ey[i]);
    }
  }
  //std::cout << merged->GetN() << std::endl;
  return merged;
}

//***********************************************************************************
// Joins two TGraphErrors back to back. Useful for plotting error envelopes
//***********************************************************************************
TGraphErrors* FSIFitterUtils::MergeGraphsIntoEnvelope(TGraphErrors *gr_min, TGraphErrors *gr_max){

  // Fill the envelope
  if(gr_min->GetN() != gr_max->GetN()){
    std::cout << "gr_min and gr_max have different binning!" << std::endl;
    std::exit(-1);
  }

  TGraphErrors *gr_shade = new TGraphErrors(2*gr_min->GetN());

  for(int ipoint=0; ipoint<gr_min->GetN(); ipoint++){
    
    double x_min,y_min,x_max,y_max;
    //Get min and max envelopes
    gr_min->GetPoint(ipoint,x_min,y_min);
    gr_max->GetPoint(gr_min->GetN()-ipoint-1,x_max,y_max);
    
    gr_shade->SetPoint(ipoint,x_min,y_min);
    gr_shade->SetPoint(gr_min->GetN()+ipoint,x_max,y_max);
    
  }
  
  gr_shade->SetFillStyle(3013);
  gr_shade->SetFillColor(16);
  return gr_shade;

}

//***********************************************************************************
// Turns a set of TMultiDimFits into a single TF1
//***********************************************************************************
double TF1FromMultiDimFits::operator()(double *x, double *) const{

  std::vector<double> tmp;
  tmp.push_back(x[0]);
  tmp.insert(tmp.end(), fFSIPars.begin(), fFSIPars.end());
  
  double calc = 0;
  for(int i = 0; i < (int)fVectorMultiDimFits->size(); i++)
    calc += FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(fVectorMultiDimFits->at(i),tmp); 
  
  return calc;

}
