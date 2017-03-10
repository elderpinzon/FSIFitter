#include "FSIFitterUtils.hxx"


//***********************************************************************************
// Convert interaction type to integer for easiness later on
//***********************************************************************************
Int_t FSIFitterUtils::IntrNameToInt(TString fIntrType){
  
  int fIntrTypeNumber = -1;
  
  // if (nuclFit){
  //   if(fIntrType == "reac")   fIntrTypeNumber = 1;
  //   if(fIntrType == "elas")   fIntrTypeNumber = 2;
  //   if(fIntrType == "tot")    fIntrTypeNumber = 3;
  // }else{
  if(fIntrType == "reac")   fIntrTypeNumber = 1;
  if(fIntrType == "inel")   fIntrTypeNumber = 2;
  if(fIntrType == "abs")    fIntrTypeNumber = 3;
  if(fIntrType == "cx")     fIntrTypeNumber = 4;
  if(fIntrType == "abscx")  fIntrTypeNumber = 5;
  if(fIntrType == "dcx")    fIntrTypeNumber = 6;
  //}

  return fIntrTypeNumber;

}

//***********************************************************************************
// Convert target nuclei name to Z
//***********************************************************************************
Int_t FSIFitterUtils::NucleiNameToInt(TString fNuclei){
  
  int fNucleiNumber = -1;
  
  if(fNuclei == "c")   fNucleiNumber = 6;
  if(fNuclei == "o")   fNucleiNumber = 8;
  if(fNuclei == "al")   fNucleiNumber = 13;
  if(fNuclei == "fe")   fNucleiNumber = 26;
  if(fNuclei == "cu")   fNucleiNumber = 29;
  if(fNuclei == "pb")   fNucleiNumber = 82;

  return fNucleiNumber;

}

//***********************************************************************************
// Convert pion type string to int
//***********************************************************************************
Int_t FSIFitterUtils::PionTypeToInt(TString fPionType){
  
  int fPionTypeNumber = -1;
  
  if(fPionType == "piP")   fPionTypeNumber = 211;
  if(fPionType == "piM")   fPionTypeNumber = -211;

  return fPionTypeNumber;

}

//***********************************************************************************
// Convert parameter indeces into single index
//***********************************************************************************
Int_t FSIFitterUtils::TurnParamsToIndex(Double_t *kindex){

  Double_t index = 0;
  for(Int_t i = 0; i < nFSIparsFitted; i++){
    Int_t factor = 1.0;
    for(Int_t j = i+1; j < nFSIparsFitted; j++){
      factor *= nSteps[j];
    }
    
    index += kindex[i] * factor;
    //printf("%s: %.2f %d\n",parNames[i].Data(),kindex[i],factor);
  }

  // Add 0.5 for proper rounding
  Int_t roundedIndex = index+0.5;
  return roundedIndex;
}

//***********************************************************************************
// Convert single index into parameter indeces
//***********************************************************************************
void FSIFitterUtils::TurnIndexToParams(Int_t index, Double_t *kindex){

  for(Int_t i = 0; i < nFSIparsFitted; i++){
    
    // Start from full index
    kindex[i] = index;

    // Do modulus a few times
    for(Int_t j = 0; j < i; j++){
    
      Int_t factor_modulus = 1;      

      // Determine the factor for this step
      for(Int_t k = j+1; k < nFSIparsFitted; k++){
	factor_modulus *= nSteps[k];
      }
      
      // Take modulus
      kindex[i] = (Int_t)kindex[i] % factor_modulus;
      printf("factor_modulus: %d kindex[%d]: %.2f\n",factor_modulus,i,kindex[i]);
    }
    
    // Finally do division
    Int_t factor = 1;
    for(Int_t j = i+1; j < nFSIparsFitted; j++){
      factor *= nSteps[j];
    }
    
    kindex[i] = (Int_t)kindex[i]/factor;
    printf("factor: %d kindex[%d]: %.2f\n",factor,i,kindex[i]);

    // Note: kindex is returned by reference
  }
}

//***********************************************************************************
// Build resulting function from the results of given TMultiDim object and evaluate at x
// Most of this code was copied from TMultiDimFit::MakeRealCode()
//***********************************************************************************
double FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(TMultiDimFit *aTMultiDimFit, const std::vector<double> &x0){

  // HACK!!!!
  // Remove FEFCXH manually here
  std::vector<double> x;
  x = x0;
  // x.erase(x.begin() + 5);
  
  // std::cout << "\t... FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate ";
  // for(int i = 0; i <x.size(); i++){
  //   std::cout << " " << x[i];
  // }
  // std::cout << std::endl;

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
// Print the parameter set.
//***********************************************************************************
void FSIFitterUtils::PrintParameterSet(double FSIPars[nFSIpars]){
  
  for(int par = 0; par<nFSIpars; par++)
    std::cout << FSIPars[par] << " ";
  std::cout << std::endl;
  
}
//***********************************************************************************
// Calculate correlation matrix from covariance matrix
//***********************************************************************************
TMatrixD FSIFitterUtils::CovarianceToCorrelation(TMatrixD* h){

  TMatrixD h_diag;
  h_diag.ResizeTo((*h));
  TMatrixD correlation;
  correlation.ResizeTo((*h));
    
  for(Int_t i=0; i<h->GetNrows(); i++){
    for(Int_t j=0; j<h->GetNcols(); j++){
      if(i==j) h_diag[i][j] = 1/TMath::Sqrt((*h)[i][j]);
      else  h_diag[i][j] = 0;
    }
  }
  
  correlation = h_diag * (*h) * h_diag;
  
  return correlation;

}

//***********************************************************************************
// Merge covariances by appending one to one corner of the other
//***********************************************************************************
void FSIFitterUtils::MergeMatrices(TMatrixD &matrix1, TMatrixD matrix2){

  TMatrixD newMatrix;
  newMatrix.ResizeTo(matrix1.GetNcols()+matrix2.GetNcols(),matrix1.GetNrows()+matrix2.GetNrows());

  for(int ix = 0; ix < newMatrix.GetNcols(); ix++){
    for(int iy = 0; iy < newMatrix.GetNrows(); iy++){

      if(ix < matrix1.GetNcols() && iy < matrix1.GetNrows())
        newMatrix(ix,iy) = matrix1(ix,iy);
      else if(ix >= matrix1.GetNcols() && iy >= matrix1.GetNrows())
        newMatrix(ix,iy) = matrix2(ix - matrix1.GetNcols(), iy - matrix1.GetNcols());
      else
        newMatrix(ix,iy) = 0;
    }
  }

  // Copy newMatrix back to matrix1 since that one gets passed by reference with the result
  matrix1.ResizeTo(matrix1.GetNcols()+matrix2.GetNcols(),matrix1.GetNrows()+matrix2.GetNrows());
  matrix1 = newMatrix;

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

  
