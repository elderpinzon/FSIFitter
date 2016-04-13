#include "FSIChi2Grid.hxx"
#include "TH1.h"
#include <iomanip>

//***********************************************************************************
// Definition of a grid object. A vector of the datasets to be used for the chisquare
// calculation is an argument
//***********************************************************************************
//FSIChi2Grid::FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets, FSIParameterScan &aFSIParameterScan) : fFSIParameterScan(aFSIParameterScan){
FSIChi2Grid::FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets){
  
  std::cout << "\nStarting Grid Object. Will use these datasets:" << std::endl;
  for (int dataset = 0; dataset < (int)AllDataSets.size(); dataset++){
    std::cout << "\t"<< dataset+1 << ". " << AllDataSets[dataset].GetSetName() << std::endl;
  }

  //fFSIParameterScan = aFSIParameterScan;
  //std::cout << "\nWill use MC parameter scan with # entries: " << fFSIParameterScan.fChain->GetEntries() << std::endl;
  //fFSIParameterScan.Test();

  // Copy vector of data set to internal data set vector
  fAllDataSets = AllDataSets;

  // Set-up TTree for storage of the finite grid
  ChiSquareFiniteGrid = new TTree("Chi2Grid","Finite Chi Square Grid Tree");
  ChiSquareFiniteGrid->Branch("qe",&fqe,"qe/D");
  ChiSquareFiniteGrid->Branch("abs",&fabs,"abs/D");
  ChiSquareFiniteGrid->Branch("cx",&fcx,"cx/D");
  ChiSquareFiniteGrid->Branch("chisquare",&fchisquare,"chisquare/D");

}

//***********************************************************************************
// Alternative definition of a grid object. 
// Uses previously calculated grid TTree to save time
//***********************************************************************************
FSIChi2Grid::FSIChi2Grid(TTree &BuiltGridTree) : ChiSquareFiniteGrid(&BuiltGridTree){

  std::cout << "\nStarting Grid Object. Will use TTree with # entries: " 
	    << ChiSquareFiniteGrid->GetEntries() << std::endl;

}

//***********************************************************************************
// Load entries of TMultiDimFit object using exisitng TTree
//***********************************************************************************
void FSIChi2Grid::UseExistingGrid(){
  
  // Set up the branches to be read from existing TTree
  ChiSquareFiniteGrid->SetBranchAddress("qe",&fqe);
  ChiSquareFiniteGrid->SetBranchAddress("abs",&fabs);
  ChiSquareFiniteGrid->SetBranchAddress("cx",&fcx);
  ChiSquareFiniteGrid->SetBranchAddress("chisquare",&fchisquare);

  // Loop over entries and fill TMultiDimFit object
  for(int ientry = 0; ientry < ChiSquareFiniteGrid->GetEntries() ; ientry++){
    
    ChiSquareFiniteGrid->GetEntry(ientry);
    
    std::cout << "Adding qe:" << fqe 
	      << " abs: " << fabs 
	      << " cx: " << fcx 
	      << " chisquare : " << fchisquare
	      <<" to TMultiDimFit object"
	      << std::endl;
    
    // Add this set of parameters and its chisquare to the TMultiDim object
    Double_t fsi_pars[3];
    fsi_pars[0] = fqe;
    fsi_pars[1] = fabs;
    fsi_pars[2] = fcx;
    multifit->AddRow(fsi_pars,fchisquare);

  }

  std::cout<< std::endl << "FSIChi2Grid::UseExisitingGrid()... done" << std::endl;

}

//***********************************************************************************
// The finite chisquare grid will be built here. Each entry is stored in a TTree
// Most importantly and entry for each parameter set is added to the TMultiDim object
//*********************************************************************************** 
void FSIChi2Grid::BuildFiniteGrid(){

  std::cout << "Building finite chisquare grid" << std::endl;
  
  // Right now the boundaries and step sizes are hard-coded
  // This should be fixed somehow
  for(Double_t qe = 0.75; qe < 1.2; qe = qe+0.05){
    for(Double_t abs = 1.0; abs < 1.55; abs = abs+0.05){
      for(Double_t cx = 0.50; cx < 1.05; cx = cx+0.05){

	// To collect individual contributions to chi-square from datasets
	Double_t totalChiSquare = 0;
	
	// Loop through datasets and calculate chisquare for a parameter set
	for (int dataset = 0; dataset < (int)fAllDataSets.size(); dataset++){
	  //double this_chisquare = fAllDataSets[dataset].CalculateChiSquare(qe,abs,cx,fFSIParameterScan);
	  double this_chisquare = fAllDataSets[dataset].CalculateChiSquare(qe,abs,cx);
	  std::cout << "Individual dataset: " << fAllDataSets[dataset].GetSetName() 
		    << " has contribution to chisquare: " << this_chisquare << std::endl;
	  totalChiSquare += this_chisquare;
	}
	
	// Weird conversion to store elements in TTree with grid info
	fqe = qe;
	fabs = abs;
	fcx = cx;
	fchisquare = totalChiSquare;
	ChiSquareFiniteGrid->Fill();

	std::cout << "\nAdding qe:" << fqe
		  << " abs: " << fabs
		  << " cx: " << fcx
		  << " chisquare : " << fchisquare
		  <<" to TMultiDimFit object\n"
		  << std::endl;
	
	// Add this set of parameters and its chisquare to the TMultiDim object
	Double_t fsi_pars[3];
	fsi_pars[0] = qe;
	fsi_pars[1] = abs;
	fsi_pars[2] = cx;
	multifit->AddRow(fsi_pars,fchisquare);
	
      }
    }
  }
  
  // Store chisquare finite grid in output file
  TFile *fout = new TFile("output/finite_chi2_grid.root","RECREATE");
  ChiSquareFiniteGrid->Write();
  fout->Close();

  std::cout<< std::endl << "FSIChi2Grid::BuildGrid()... done" << std::endl;

}
//***********************************************************************************
// Just return the TTree with the finite grid
//***********************************************************************************
TTree* FSIChi2Grid::GetFiniteGrid(){
  
  return ChiSquareFiniteGrid;

}

//***********************************************************************************
// This is the method that is called in the main executable
// The TMultiDim object is initialized and either BuildFiniteGrid() or UseExisitingGrid()
// are used to fill it with the multi-dimensional points to be interpolated.
// The resulting function is built in FSIChi2Grid::InterpolatedGrid()
//***********************************************************************************
void FSIChi2Grid::InterpolateFiniteGrid(bool kUseInputGrid){

  std::cout << "\nBuilding interpolated chisquare grid" << std::endl;

  // Create MultiDim fitter
  multifit = new TMultiDimFit(3,TMultiDimFit::kMonomials,"kv");
  //Int_t mPowers[] = {10,10,10,10};
  Int_t mPowers[] = {5,5,5,5};
  multifit->SetMaxPowers(mPowers);
  multifit->SetMaxFunctions(10000);
  multifit->SetMaxStudy(10000);
  multifit->SetMaxTerms(100);
  multifit->SetPowerLimit(10);
  multifit->SetMinRelativeError(0.001);

  // Build Finite Grid here and fill the rows of the MultiDimFit
  if(kUseInputGrid) UseExistingGrid();
  else BuildFiniteGrid();
  
  // Print out the statistics of input sample
  multifit->Print("s");
  
  // Run interpolation
  multifit->MakeHistograms();
  multifit->FindParameterization("");
  //multifit->Fit("M");
  
  // Print and store interpolation result
  multifit->Print("prcfkm");
  multifit->MakeMethod("output/chi2");

  // Store interpolation results in output file
  TFile *fout = new TFile("output/interpolation_results.root","RECREATE");
  multifit->Write();
  const TMatrixD *inter_corr = multifit->GetCorrelationMatrix();
  inter_corr->Write("inter_corr");
  inter_corr->Print();
  fout->Close();

  std::cout<< std::endl << "FSIChi2Grid::InterpolateFiniteGrid()... done" << std::endl;
  
}

//***********************************************************************************
// Build resulting function from the results of the TMultiDim object
// Most of this code was copied from TMultiDimFit::MakeRealCode()
//***********************************************************************************
double FSIChi2Grid::GetInterpolatedGridPoint(const std::vector<double> &x){
  
  // Grab all necessary pieces
  Int_t gNCoefficients = multifit->GetNCoefficients();
  Int_t gNVariables = multifit->GetNVariables();
  const Int_t *gPower = multifit->GetPowers();
  const Int_t *gPowerIndex = multifit->GetPowerIndex();
  const TVectorD *gCoefficient = multifit->GetCoefficients();
  const TVectorD *gXMax = multifit->GetMaxVariables();
  const TVectorD *gXMin = multifit->GetMinVariables();
  const double gDMean =  multifit->GetMeanQuantity();
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
