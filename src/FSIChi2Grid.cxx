#include "FSIChi2Grid.hxx"
#include "TStopwatch.h" // remove later
#include "TROOT.h"
#include "TCanvas.h" //remove later

//***********************************************************************************
// Definition of a grid object. A vector of the datasets to be used for the chisquare
// calculation is an argument
//***********************************************************************************
FSIChi2Grid::FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets){

  kOctave = false;
  kMultiDimFit = false;
  nPointsUsed = 0;
  
  std::cout << "\nStarting Grid Object. Will use these datasets:" << std::endl;
  for (int dataset = 0; dataset < (int)AllDataSets.size(); dataset++){
    std::cout << "\t"<< dataset+1 << ". " << AllDataSets[dataset].GetFileName() << std::endl;
    nPointsUsed += AllDataSets[dataset].GetNumberOfDataPoints();
  }
  
  std::cout << "... " << AllDataSets.size() << " data sets with "
	    << nPointsUsed << " points used in total." << std::endl;
  
  // Copy vector of data set to internal data set vector
  fAllDataSets = AllDataSets;
  
  // Store chisquare finite grid in output file
  foutTree = new TFile("output/finite_chi2_grid.root","RECREATE");

  // Set-up TTree for storage of the finite grid
  ChiSquareFiniteGrid = new TTree("Chi2Grid","Finite Chi Square Grid Tree");

}

//***********************************************************************************
// Builds a giant covariance by stacking all the individual "made-up" covariances for
// each data set. There are only correlation among poins on each dataset
//***********************************************************************************
void FSIChi2Grid::BuildCovariance(){

  std::cout << "\n--- Building full covariance matrix" 
	    << std::endl;
  
  // Intermediate tools
  int pointIndex = 0;
  allErrors.ResizeTo(nPointsUsed);
  allMomenta.ResizeTo(nPointsUsed);
  allSigmaData.ResizeTo(nPointsUsed);
  allTargets.ResizeTo(nPointsUsed);
  allPids.ResizeTo(nPointsUsed);
  allTypes.ResizeTo(nPointsUsed);
  
  for (int dataset = 0; dataset < (int)fAllDataSets.size(); dataset++){

    // Add covariance for this dataset to the full covariance;
    FSIFitterUtils::MergeMatrices(fullCovariance,fAllDataSets[dataset].GetCovariance());

    // Need to keep track of the information on each entry of the matrix
    for(int point = 0; point < fAllDataSets[dataset].GetXsecErrorVector().GetNoElements(); point++){
      allMomenta[pointIndex] = fAllDataSets[dataset].GetMomVector()[point];
      allSigmaData[pointIndex] = fAllDataSets[dataset].GetXsecVector()[point];
      allErrors[pointIndex] = fAllDataSets[dataset].GetXsecErrorVector()[point];
      allTargets[pointIndex] = fAllDataSets[dataset].GetNucleus();
      allPids[pointIndex] = fAllDataSets[dataset].GetPionType();
      allTypes[pointIndex] = fAllDataSets[dataset].GetIntrType();
      allTargetsString.push_back(fAllDataSets[dataset].GetNucleusString());
      
      // UGLY HACK TO DEAL WITH DUET SEPARATE
      // First 5 points are CX
      if(fAllDataSets[dataset].GetSetName() == "DUETSeparate" && point < 5) 
	allTypes[pointIndex] = 4;
      // Last 5 points are ABS
      if(fAllDataSets[dataset].GetSetName() == "DUETSeparate" && point >= 5)
	allTypes[pointIndex] = 3;
      
      pointIndex++;
    }
  }

  // Caclulate correlation from covariance
  fullCorrelation.ResizeTo(nPointsUsed,nPointsUsed);
  fullCorrelation = FSIFitterUtils::CovarianceToCorrelation(&fullCovariance);
  std::cout << "fullCovariance nrows: " << fullCovariance.GetNrows()
	    << " nPointsUsed: " << nPointsUsed
	    << std::endl;
  fullCovarianceInv.ResizeTo(nPointsUsed,nPointsUsed);
  fullCovarianceInv = fullCovariance;
  fullCovarianceInv.Invert();

  // // Store in a file for checks
  // TFile *fcov = new TFile("fcov.root","RECREATE");
  // fullCovariance.Write("fullCovariance");
  // fullCovarianceInv.Write("fullCovarianceInverse");
  // fullCorrelation.Write("fullCorrelation");
  // fcov->Close();

  // Check that it can be Cholesky decomposed
  TDecompChol decomp(fullCovariance);
  if (!decomp.Decompose()){
    std::cout << "Decomposition failed for fullCovariance" << std::endl;
    std::exit(1);
  }
  
  std::cout << "--- done building full covariance matrix" << std::endl;

}

//***********************************************************************************
// The finite chisquare grid will be built here. Each entry is stored in a TTree
// Most importantly and entry for each parameter set is added to the TMultiDim object
//***********************************************************************************
// Fills TVectorD's for each parameter index and then multiples them with the
// covariance matrix to calculate the chi2's.
//*********************************************************************************** 
// Builds only with a 5D parameter grid
//*********************************************************************************** 
void FSIChi2Grid::FastBuildFiniteGrid(){

  std::cout << "FSIChi2Grid::FastBuildFiniteGrid: Building finite chisquare grid" << std::endl;

  TStopwatch s;
  s.Start();

  // Define output TTree branches
  ChiSquareFiniteGrid->Branch("FEFQE",&fFSIPars[0],"FEFQE/D");
  ChiSquareFiniteGrid->Branch("FEFABS",&fFSIPars[1],"FEFABS/D");
  ChiSquareFiniteGrid->Branch("FEFCX",&fFSIPars[2],"FEFCX/D");
  ChiSquareFiniteGrid->Branch("FEFINEL",&fFSIPars[3],"FEFINEL/D");
  ChiSquareFiniteGrid->Branch("FEFQEH",&fFSIPars[4],"FEFQEH/D");
  ChiSquareFiniteGrid->Branch("chisquare",&fchisquare,"chisquare/D");
    
  FSIParameterScan aFSIParameterScan(FileName,nuclFit);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "FSIChi2Grid::FastBuildFiniteGrid from grid scan: " << FileName 
	    << " with " << nentries << " entries."
	    << std::endl;

  // Initialize array in the heap
  Float_t** allSigmaMC = new Float_t*[nVariations];
  for(int i = 0; i < nVariations; ++i)
    allSigmaMC[i] = new Float_t[nMaxPointsUsed];

  // Make sure that the number of data point is under the limit
  if(nPointsUsed > nMaxPointsUsed){
    std::cout << "FATAL ERROR: nPointsUsed>nMaxPointsUsed "
	      << nPointsUsed << std::endl;
    std::exit(1);
  }

  // Loop through the tree
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;
    
    if(jentry % 100000 == 0)
      printf("... aFSIParameterScan entry #%d. %.2f%% done...\n",(int)jentry,(Double_t)jentry/nentries*100.);
    
    aFSIParameterScan.GetEntry(jentry);

    // Doing 5D fit so ignore some entries
    if(TMath::Abs(aFSIParameterScan.FEFCXH-1.8) > 0.01) continue;
    if(TMath::Abs(aFSIParameterScan.FEFALL-1.0) > 0.01) continue;

    // Get individual indices
    Double_t kindex[7] = {0};
    kindex[0] = (aFSIParameterScan.FEFQE - FSIParsMin[0])/FSIParsStep[0];
    kindex[1] = (aFSIParameterScan.FEFABS - FSIParsMin[1])/FSIParsStep[1];
    kindex[2] = (aFSIParameterScan.FEFCX - FSIParsMin[2])/FSIParsStep[2];
    kindex[3] = (aFSIParameterScan.FEFINEL - FSIParsMin[3])/FSIParsStep[3];
    kindex[4] = (aFSIParameterScan.FEFQEH - FSIParsMin[4])/FSIParsStep[4];
    kindex[5] = (aFSIParameterScan.FEFCXH - FSIParsMin[5])/FSIParsStep[5];
    kindex[6] = (aFSIParameterScan.FEFALL - FSIParsMin[6])/FSIParsStep[6];

    // Turn it into a single global index
    // Note: this function will turn the 7-D grid into a nFSIParsFitted-D one
    Int_t int_index = FSIFitterUtils::TurnParamsToIndex(kindex);
    
    // Loop over allMomenta, allTargets, allPids to decide if should fill TVectorD
    for(int i = 0; i<nPointsUsed; i++){
      
      if(TMath::Abs(aFSIParameterScan.target - allTargets[i]) > 0.01) continue;
      if(TMath::Abs(aFSIParameterScan.pid - allPids[i]) > 0.01) continue;
      if(TMath::Abs(aFSIParameterScan.mom - allMomenta[i]) > 0.01) continue;
      
      // If we are still here we got it! Now pick the proper xsec
      Double_t properSigma = -1;
      // If reactive data point
      if(allTypes[i] == 1) properSigma = aFSIParameterScan.xreac;
      // If quasi-elastic data point
      if(allTypes[i] == 2) properSigma = aFSIParameterScan.xqe;
      // If ABS data point
      if(allTypes[i] == 3) properSigma = aFSIParameterScan.xabs;
      // If CX data point
      if(allTypes[i] == 4) properSigma = aFSIParameterScan.xcx;
      // If ABS+CX data point
      if(allTypes[i] == 5) properSigma = aFSIParameterScan.xcx + aFSIParameterScan.xabs;
      
      // If LE point manually vary the HE parameters
      if(aFSIParameterScan.mom < 400){

	// Hardcoded: 5D grid
	//for(int vindex5 = 0; vindex5 < 7; vindex5++){
	  for(int vindex4 = 0; vindex4 < nSteps[4]; vindex4++){
	    for(int vindex3 = 0; vindex3 < nSteps[3]; vindex3++){

	      Double_t vindex[7] = {0};
	      memcpy(vindex, kindex, sizeof(kindex));
	      vindex[4] = vindex4;
	      vindex[3] = vindex3;
	    	    	      
	      Int_t _int_index = FSIFitterUtils::TurnParamsToIndex(vindex);
	      allSigmaMC[_int_index][i] = properSigma;
	      //}
	    }
	  }
      } // end of LE loop
      
	// If HE point manually vary the LE parameters
      else if(aFSIParameterScan.mom > 500){
	
	for(int vindex2 = 0; vindex2 < nSteps[2]; vindex2++){
	  for(int vindex1 = 0; vindex1 < nSteps[1]; vindex1++){
	    for(int vindex0 = 0; vindex0 < nSteps[0]; vindex0++){
	      
	      Double_t vindex[7] = {0};
	      memcpy(vindex, kindex, sizeof(kindex));
	      vindex[2] = vindex2;
	      vindex[1] = vindex1;
	      vindex[0] = vindex0;
	      
	      Int_t _int_index = FSIFitterUtils::TurnParamsToIndex(vindex);
              allSigmaMC[_int_index][i] = properSigma;
	    }
	  }
	}
      } // end of HE loop
      
	// In LE-HE transition region, no need to vary manually
      else{
	
	allSigmaMC[int_index][i] = properSigma;
	
      } // end of LE-HE loop
      
    }
    
  }
  
  // Actual chi2 calculation for each parameter combo
  Array<octave_idx_type> index_array(dim_vector((octave_idx_type)nFSIparsFitted));
  for(Int_t set = 0; set < nVariations; set++){

    // From global index to individual FSI indices
    Double_t *kIndices;
    kIndices = new Double_t[nFSIparsFitted];

    FSIFitterUtils::TurnIndexToParams(set,kIndices);

    // From individual indeces set to FSI parameter set
    for(Int_t i = 0; i<nFSIparsFitted; i++){
      fFSIPars[i] = kIndices[i] * FSIParsStep[i] + FSIParsMin[i];
      index_array(i) = kIndices[i];
    }
    
    // Some auxiliary vectors (in TMatrixD form for convenience)
    TMatrixD vec_diff_Xsec(nPointsUsed,1);
    TMatrixD vec_diff_Xsec_transpose(1,nPointsUsed);
    
    for(int point = 0; point < nPointsUsed; point++){
      vec_diff_Xsec[point][0] = allSigmaData[point] - allSigmaMC[set][point];

      // Make sure that the values are sensible
      if(allSigmaData[point] <= 0 || allSigmaMC[set][point] <= 0){
	std::cout << "Missing digma. point: " << point << " set: " << set 
		  << " allSigmaData: " << allSigmaData[point]
		  << " allSigmaMC: " << allSigmaMC[set][point] << std::endl;
	std::exit(1);
      }
    }
    vec_diff_Xsec_transpose.Transpose(vec_diff_Xsec);

    // Calculate chi2 from filled vectors and covariance matrix
    TMatrixD res(1,1);
    res = vec_diff_Xsec_transpose * fullCovarianceInv * vec_diff_Xsec;
    Double_t _kChiSquare = res[0][0];
    if(set%1000==0) std::cout << "\n_kChiSquare[" << set << "]: " << _kChiSquare << std::endl;
    
    // Fill TMultiDimFit
    if(kMultiDimFit)
      multifit->AddRow(fFSIPars,_kChiSquare);
    
    // Or the octave matrix
    if(kOctave)
      (*v_vec)(index_array) = _kChiSquare;
    
    fchisquare = _kChiSquare;
    ChiSquareFiniteGrid->Fill();
    
  }

  foutTree->cd();
  ChiSquareFiniteGrid->Write();

  std::cout<< std::endl << "FSIChi2Grid::FastBuildGrid()... done in " 
	   << s.RealTime() << " seconds. Stored grid with entries: " 
	   << ChiSquareFiniteGrid->GetEntries() << std::endl;
  
}


//***********************************************************************************
// This is the method that is called in the main executable
// The TMultiDim object is initialized and either BuildFiniteGrid() or UseExisitingGrid()
// are used to fill it with the multi-dimensional points to be interpolated.
// The resulting function is built in FSIChi2Grid::InterpolatedGrid()
//***********************************************************************************
void FSIChi2Grid::InterpolateFiniteGrid(){

  std::cout << "\nFSIChi2Grid::Building interpolated chisquare grid" << std::endl;

  // For octave need to first set up the grid elements
  if (kOctave)
    BuildOctaveGrid();
  // Or for TMultiDimFit set up the TMultiDimFit object
  else if(kMultiDimFit)
    BuildTMultiDimFitGrid();

  // Actually fill the grid by reading from input file
  FastBuildFiniteGrid();
  
  // For octave nothing else needs to happen beacuse it is a spline
  // interpolation. GetSplinedGridPoint() can simply be called
  
  // For TMultiDimFit we need to actually run the interpolation algorithm
  // that will get the polynomial functions
  
  if(kMultiDimFit){
    
    RunTMultiDimFitInterpolation();
  }
  
  std::cout<< std::endl << "FSIChi2Grid::InterpolateFiniteGrid()... done" << std::endl;
}

//***********************************************************************************
// Set up the octave elements that hold the grid
//***********************************************************************************
void FSIChi2Grid::BuildOctaveGrid(){
  
  // Octave initializations
  string_vector argv (2);
  argv(0) = "FSIOctaveInterpolation";
  argv(1) = "-q";
  octave_main (2, argv.c_str_vec (), 1);
  
  for (int idim=0; idim<nFSIparsFitted; idim++){
    
    std::cout << "FSIChi2Grid::Building octave array of parameters: " << idim << " dim: " << nSteps[idim] << std::endl;
    
    // Vector with FSI pars
    x_vec[idim] = new Matrix(nSteps[idim],1);
    for (octave_idx_type i = 0; i < nSteps[idim]; i++){
      (*x_vec[idim])(i,0) = FSIParsMin[idim] + i*FSIParsStep[idim];
    }
    //std::cout << (*x_vec[idim]);
  }
  
  // Hardcoded 5 parameters here!!
  dim_vector dim_vec(nSteps[0],nSteps[1],nSteps[2],nSteps[3],nSteps[4]);
  
  // Initialize with large chi^2 in case some failed in filling tree
  v_vec = new Array<double>(dim_vec,999);
  //std::cout << "v_vec at the start: " << (*v_vec) << std::endl; 
  
}

//***********************************************************************************
// Set up the TMultiDimFit object
//***********************************************************************************
void FSIChi2Grid::BuildTMultiDimFitGrid(){
  
  // Create TMultiDimFit fitter
  multifit = new TMultiDimFit(nFSIparsFitted,TMultiDimFit::kMonomials,"kv");

  // These parameters have been "manually" tuned to give the best
  // interpolation result, but feel free to play with them
  Int_t mPowers[7];
  mPowers = {5,5,5,5,5,5,5};
  multifit->SetMaxPowers(mPowers);
  multifit->SetMaxFunctions(5000);
  multifit->SetMaxStudy(300000);
  multifit->SetMaxTerms(1000);
  multifit->SetMinRelativeError(0.001);

}

//***********************************************************************************
// Run interpolation algorithm (FindParameterization) for the TMultiDimFit object
//***********************************************************************************
void FSIChi2Grid::RunTMultiDimFitInterpolation(){
  
  std::cout << "Starting TMultiDimFit interpolation of Chi2 grid..." << std::endl;
  TStopwatch s;
  s.Start();
  
  // Print out the statistics of input sample
  multifit->Print("s");
  
  // Run interpolation
  std::cout << "multifit->MakeHistograms()" << std::endl;
  multifit->MakeHistograms();
  std::cout << "multifit->FindParameterization()" << std::endl;
  multifit->FindParameterization("");
  
  // Print and store interpolation result
  multifit->Print("prcfkm");
  multifit->MakeMethod("InterMultiDimChi2");

  s.Stop();
  
  std::cout << "Finished interpolation of Chi2 grid in "
	    << s.RealTime() << " seconds." << std::endl;
  
}

//***********************************************************************************
// Check if point is within the grid
//***********************************************************************************
bool FSIChi2Grid::CheckInsideGridBoundary(const std::vector<double> &x){
  
  bool IsOutside = false;
  for(int i = 0; i < (int)x.size(); i++){
    if(x[i] > FSIParsMax[i] || x[i] < FSIParsMin[i])
      IsOutside = true;
  }

  return IsOutside;

}

//***********************************************************************************
// Evaluate interpolated function at required point in n-dimensional space
//***********************************************************************************
Double_t FSIChi2Grid::GetInterpolatedGridPoint(const std::vector<Double_t> &x){

  if(CheckInsideGridBoundary(x))
    return 999.9;

  return FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(multifit,x);

}

//***********************************************************************************
// Get point using interpn spline function of GNU-Octave
//***********************************************************************************
Double_t FSIChi2Grid::GetSplinedGridPoint(const std::vector<Double_t> &x){

  Matrix *x_point[nFSIparsFitted];
  const int nPointsToInterp = 1;
  for (int idim=0; idim<nFSIparsFitted; idim++) {
    x_point[idim] = new Matrix(nPointsToInterp,1);
    (*x_point[idim])(0,0) = x[idim];
  }

  octave_value_list in;
  for (int idim=0; idim<nFSIparsFitted; idim++) {
    in.append( octave_value( *x_vec[idim] ) );
  }

  in.append( octave_value (*v_vec) );
  
  for (int idim=0; idim<nFSIparsFitted; idim++)
    in.append( octave_value( *x_point[idim] ) );
  
  in.append( octave_value("spline"));  // Default is linear
  in.append( octave_value(999));       // 999 outside the grid
  octave_value_list out_v_point = feval ("interpn",in,nPointsToInterp);
  
  Matrix v_point = Matrix(1,1);
  if (!error_state && out_v_point.length () > 0) {
    v_point = out_v_point(0).matrix_value ();
  }
  double splined_chi2 = v_point(0);
  // std::cout << "splined_chi2  : " << splined_chi2 << std::endl;
  
  return splined_chi2;

}

//***********************************************************************************
// Used to plot slices of the interpolated chi2 surface.
// Note: probably no longer works
//***********************************************************************************
void FSIChi2Grid::PlotChiSquare(const std::vector<Double_t> &x_center, Int_t set){

  gROOT->ProcessLine(".x ~/.rootlogon.C");

  Double_t true_chi2;
  ChiSquareFiniteGrid->SetBranchAddress("chisquare",&true_chi2);
  ChiSquareFiniteGrid->Print();
  Double_t kindex[nFSIparsFitted];
  
  TCanvas *c_inter = new TCanvas("","",0,0,2000,800);
  c_inter->Divide(4,2);
  TPaveText *pt = new TPaveText(0.15,0.15,0.95,.87);
  pt->AddText("Around point");
  pt->AddLine(.0,.875,1.,.875);

  for(Int_t k = 0; k<nFSIparsFitted; k++){

    std::cout << "Plotting for set: " << set
	      << " par " << parNames[k].Data()
	      <<std::endl;

    std::vector<Double_t> x = x_center;

    // Fill "true graph" using the ChiSquareFiniteGrid tree
    // This is the finite grid before interpolation
    TGraph *gTrue = new TGraph(nSteps[k]);
    for(Int_t i = 0; i<nSteps[k]; i++){
      
      // For x-axis grid
      x[k] = FSIParsMin[k] + FSIParsStep[k]*i;

      // Get index
      for(Int_t j = 0; j<nFSIpars; j++){
	kindex[j] = (x[j] - FSIParsMin[j])/FSIParsStep[j];
      }
      
      Int_t int_index = FSIFitterUtils::TurnParamsToIndex(kindex);
      printf("int_index: %d",int_index);
      ChiSquareFiniteGrid->GetEntry(int_index);
      ChiSquareFiniteGrid->Show(int_index);
      gTrue->SetPoint(i,x[k],true_chi2);
      //gTrue->SetPoint(i,x[k],fchisquare);
      
    }
    gTrue->Print();
    
    // Number of steps for TGraphs with interpolation
    Int_t pSteps = 100;

    // TGraph with multidimfit interpolation (and many more points)
    TGraph *gMultiDim = new TGraph(pSteps);

    std::cout << "Got multidimfit with variables: " << multifit->GetNVariables() << std::endl;
    for(int i = 0; i<=pSteps; i++){
      
      x[k] = FSIParsMin[k] + (FSIParsMax[k] - FSIParsMin[k])/(double)pSteps * i;
      //std::vector<double> y2 = {x[0],x[1],x[2],x[3],x[4]};;
      //double value2 = GetInterpolatedGridPoint(y);
      double value2 = GetInterpolatedGridPoint(x);
      std::cout << "y[" << k <<"]: " << x[k] << " value2: " << value2 << std::endl;
      gMultiDim->SetPoint(i,x[k],value2);
      
    }

    // TGraph with octave interpolation (and many more points)
    TGraph *gOctave = new TGraph(pSteps);
    
    for(Int_t i = 0; i<=pSteps; i++){
      
      x[k] = FSIParsMin[k] + (FSIParsMax[k] - FSIParsMin[k])/(double)pSteps * i;
      if(x[k] == 1.551) x[k] = 1.55; // ugly hack to deal with earlier ugly hack
      Double_t value = GetSplinedGridPoint(x);
      //std::cout << "x[" << k <<"] : " << x[k] << " value: " << value << std::endl;
      gOctave->SetPoint(i,x[k],value);
      
      // double value2 = GetInterpolatedGridPoint(x);
      std::cout << "x[0] : " << x[0] << " value2: " << value << std::endl;
      // gMultiDim->SetPoint(i,x[0],value2);
      
    }
    gOctave->SetLineWidth(2);
    gOctave->SetTitle(Form("#chi2 vs. %s",parNames[k].Data()));
    gMultiDim->SetLineWidth(2);
    gMultiDim->SetLineColor(2);
    
    c_inter->cd(k+1);
    gOctave->Draw("AL");
    gMultiDim->Draw("Lsame");
    gTrue->Draw("Psame");

    pt->AddText(Form("%s: %.2f",parNames[k].Data(),x_center[k]));
		
  }

  c_inter->cd(8);
  pt->Draw();

  
  c_inter->SaveAs(Form("/home/elder/share/tmp/frame_inter_%d.png",set));
  c_inter->SaveAs(Form("/home/elder/share/tmp/frame_inter_%d.eps",set));
  
}

//***********************************************************************************
// Takes a ModelPrediction and calculates the same chi2 that is used for the fit.
// This is meant to give some comparison between the NEUT best fit and other models
// Note: Only works for "pip-carbon" now. Need to fix that. Main problem is that only have
// model predicitons for pip carbon
//***********************************************************************************
Double_t FSIChi2Grid::CalculateModelChiSquare(ModelPrediction model){
  
  std::cout << "Calculating chi-square for model " << model.GetModelName()
	    << std::endl;


  TMatrixD vec_diff_Xsec(nPointsUsed,1);
  TMatrixD vec_diff_Xsec_transpose(1,nPointsUsed);
  
  for(Int_t point = 0; point < nPointsUsed; point++){

    Double_t sigmaModel = model.GetCrossSectionValue(allTargetsString[point], "pip", iTypeString2[allTypes[point]-1], allMomenta[point]);
      
    vec_diff_Xsec[point][0] = allSigmaData[point] - sigmaModel;

    std::cout << allMomenta[point] << " " 
	      << allTypes[point] << " "
	      << sigmaModel << " "
	      << allSigmaData[point] << " "
	      << vec_diff_Xsec[point][0] << " "
	      << std::endl;
    
  }

  vec_diff_Xsec_transpose.Transpose(vec_diff_Xsec);
  
  TMatrixD res(1,1);
  res = vec_diff_Xsec_transpose * fullCovarianceInv * vec_diff_Xsec;
  Double_t modelChiSquare = res[0][0];
  std::cout << "\nmodelChiSquare: " << modelChiSquare << std::endl;
  
  return modelChiSquare;

}
