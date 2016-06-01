#include "FSIChi2Grid.hxx"

//***********************************************************************************
// Definition of a grid object. A vector of the datasets to be used for the chisquare
// calculation is an argument
//***********************************************************************************
//FSIChi2Grid::FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets, FSIParameterScan &aFSIParameterScan) : fFSIParameterScan(aFSIParameterScan){
FSIChi2Grid::FSIChi2Grid(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets, bool fnuclFit, bool fOctave, bool fMultiDimFit){

  // Setting grid size here temporarily
  //nFSIPars = 3;
  nuclFit=fnuclFit;
  kOctave=fOctave;
  kMultiDimFit=fMultiDimFit;

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
  if (nuclFit){
    ChiSquareFiniteGrid->Branch("reac",&freac,"reac/D");
    ChiSquareFiniteGrid->Branch("elas",&felas,"elas/D");
    ChiSquareFiniteGrid->Branch("tot",&ftot,"tot/D");
  }else{
    ChiSquareFiniteGrid->Branch("qe",&fqe,"qe/D");
    ChiSquareFiniteGrid->Branch("abs",&fabs,"abs/D");
    ChiSquareFiniteGrid->Branch("cx",&fcx,"cx/D");   
  }
  ChiSquareFiniteGrid->Branch("chisquare",&fchisquare,"chisquare/D");

}

//***********************************************************************************
// Alternative definition of a grid object. 
// Uses previously calculated grid TTree to save time
//***********************************************************************************
FSIChi2Grid::FSIChi2Grid(TTree &BuiltGridTree, bool fnuclFit, bool fOctave, bool fMultiDimFit) : ChiSquareFiniteGrid(&BuiltGridTree){
  
  nuclFit=fnuclFit;
  kOctave=fOctave;
  kMultiDimFit=fMultiDimFit;

  std::cout << "\nStarting Grid Object. Will use TTree with # entries: " 
	    << ChiSquareFiniteGrid->GetEntries() << std::endl;

}

//***********************************************************************************
// Load entries of TMultiDimFit object using exisitng TTree
//***********************************************************************************
void FSIChi2Grid::UseExistingGrid(){
  
  // Set up the branches to be read from existing TTree
  if (nuclFit){
    ChiSquareFiniteGrid->SetBranchAddress("reac",&freac);
    ChiSquareFiniteGrid->SetBranchAddress("elas",&felas);
    ChiSquareFiniteGrid->SetBranchAddress("tot",&ftot);
  }else{
    ChiSquareFiniteGrid->SetBranchAddress("qe",&fqe);
    ChiSquareFiniteGrid->SetBranchAddress("abs",&fabs);
    ChiSquareFiniteGrid->SetBranchAddress("cx",&fcx);   
  }
  ChiSquareFiniteGrid->SetBranchAddress("chisquare",&fchisquare);

  // Loop over entries and fill TMultiDimFit object
  for(int ientry = 0; ientry < ChiSquareFiniteGrid->GetEntries() ; ientry++){
    
    ChiSquareFiniteGrid->GetEntry(ientry);
  if (nuclFit){
    std::cout << "Adding reac:" << freac 
        << " elas: " << felas 
        << " tot: " << ftot 
        << " chisquare : " << fchisquare
        <<" to TMultiDimFit object"
        << std::endl;
  }else{
    std::cout << "Adding qe:" << fqe 
        << " abs: " << fabs 
        << " cx: " << fcx 
        << " chisquare : " << fchisquare
        <<" to TMultiDimFit object"
        << std::endl;
  }  
    
    // Convert back to indices (this is dumb)
    int i_qe = (fqe-0.75)/0.05;
    int i_abs = (fabs-1.0)/0.05;
    int i_cx = (fcx-0.5)/0.05;

    // Add this set of parameters and its chisquare to the TMultiDim object
    Double_t fsi_pars[3];
  if (nuclFit){
    fsi_pars[0] = freac;
    fsi_pars[1] = felas;
    fsi_pars[2] = ftot;
  }else{
    fsi_pars[0] = fqe;
    fsi_pars[1] = fabs;
    fsi_pars[2] = fcx;
  }  
    multifit->AddRow(fsi_pars,fchisquare);
    (*v_vec)(i_qe,i_abs,i_cx) = fchisquare;
    
  }

  std::cout << "v_vec at the end: " << (*v_vec) << std::endl; 

  std::cout<< std::endl << "FSIChi2Grid::UseExisitingGrid()... done" << std::endl;

}

//***********************************************************************************
// The finite chisquare grid will be built here. Each entry is stored in a TTree
// Most importantly and entry for each parameter set is added to the TMultiDim object
//*********************************************************************************** 
void FSIChi2Grid::BuildFiniteGrid(){

  std::cout << "Building finite chisquare grid" << std::endl;
  
  if (nuclFit){
    for(Double_t reac = 0.9; reac < 1.5; reac = reac+0.1){
      for(Double_t elas = 1.0; elas < 1.6; elas = elas+0.1){
        for(Double_t tot = 0.7; tot < 1.3; tot = tot+0.1){

    // To collect individual contributions to chi-square from datasets
    Double_t totalChiSquare = 0;
    
    // Loop through datasets and calculate chisquare for a parameter set
    for (int dataset = 0; dataset < (int)fAllDataSets.size(); dataset++){
      //double this_chisquare = fAllDataSets[dataset].CalculateChiSquare(qe,abs,cx,fNucleonNEUTParameterScan);
      double this_chisquare = fAllDataSets[dataset].CalculateChiSquare(reac,elas,tot);
      std::cout << "Individual dataset: " << fAllDataSets[dataset].GetSetName() 
          << " has contribution to chisquare: " << this_chisquare << std::endl;
      totalChiSquare += this_chisquare;
    }
    
    // Weird conversion to store elements in TTree with grid info
    freac = reac;
    felas = elas;
    ftot = tot;
    fchisquare = totalChiSquare;
    ChiSquareFiniteGrid->Fill();
    
    std::cout << "\nAdding reac:" << freac
        << " elas: " << felas
        << " tot: " << ftot
        << " chisquare : " << fchisquare
        <<" to TMultiDimFit object\n"
        << std::endl;
    // Add this set of parameters and its chisquare to the TMultiDim object
    Double_t fsi_pars[3];
    fsi_pars[0] = reac;
    fsi_pars[1] = elas;
    fsi_pars[2] = tot;
    multifit->AddRow(fsi_pars,fchisquare);
    
        }
      }
    }
  }else{

   /*
  // Was trying to simplify loops but gave up :(
  int fsi_iter[nFSIpars];
  double fsi_par[nFSIpars];

  for (int idim=0; idim<nFSIpars; idim++){

    for(fsi_iter[idim] = 0; fsi_iter[idim] < nSteps[idim]; fsi_iter[idim]++){

      fsi_par[idim] = grid_start[idim] + grid_step[idim] * fsi_iter[idim];
      printf("%d %d %.2f\n",idim,fsi_iter[idim],fsi_par[idim]);
    }
  }
  */

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
	Double_t fsi_pars[nFSIpars];
	fsi_pars[0] = qe;
	fsi_pars[1] = abs;
	fsi_pars[2] = cx;
	multifit->AddRow(fsi_pars,fchisquare);

  if (kOctave){
	
	// Convert back to indices (this is dumb)
	int i_qe = (qe-0.75)/0.05;
	int i_abs = (abs-1.0)/0.05;
	int i_cx = (cx-0.5)/0.05;

	// Set chi2 in mesh
	std::cout << "\nSetting grid point i_qe:" << i_qe
                  << " i_abs: " << i_abs
                  << " i_cx: " << i_cx
                  << " chisquare : " << fchisquare
                  <<" to Octave mesh\n"
                  << std::endl;
	(*v_vec)(i_qe,i_abs,i_cx) = fchisquare;
	
      }
    }
  }
}
  if (kOctave) std::cout << "v_vec at the end: " << (*v_vec) << std::endl; 

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
if (kOctave){
  // Grid initializations
  double grid_step[nFSIpars] = {0.05,0.05,0.05};  
  double grid_start[nFSIpars] = {0.75,1.0,0.5};
  //double grid_stops[nFSIpars] = {1.15,1.5,1.0};
  int nSteps[nFSIpars] = {9,11,11};
  
  // Octave initializations
  string_vector argv (2);
  argv(0) = "FSIOctaveInterpolation";
  argv(1) = "-q";
  octave_main (2, argv.c_str_vec (), 1);

  for (int idim=0; idim<nFSIpars; idim++){
    
    // Number of step grid points
    //nSteps[idim] = (grid_stops[idim]-grid_start[idim])/grid_step[idim] + 1;

    std::cout << "Building octave array of parameters: " << idim << " dim: " << nSteps[idim] << std::endl;

    // Vector with FSI pars
    x_vec[idim] = new Matrix(nSteps[idim],1);
    for (octave_idx_type i = 0; i < nSteps[idim]; i++){
      std::cout << grid_start[idim] + i*grid_step[idim] << std::endl;
      (*x_vec[idim])(i,0) = grid_start[idim] + i*grid_step[idim];
    }
    std::cout << (*x_vec[idim]);
  }
  
  //dim_vector dim_vec(nSteps,nSteps);       // Start off 2D
  //dim_vec.resize(nFSIpars,nSteps);         // Then increase
  dim_vector dim_vec(nSteps[0],nSteps[1],nSteps[2]);
  v_vec = new Array<double>(dim_vec,999); // Initialize with large chi^2 in case some failed in filling tree
  std::cout << "v_vec at the start: " << (*v_vec) << std::endl; 

} 
// else if (kMultiDimFit){

  // Create MultiDim fitter
  multifit = new TMultiDimFit(nFSIpars,TMultiDimFit::kMonomials,"kv");
  Int_t mPowers[4];
  // Create MultiDim fitter
  if (nuclFit) mPowers = {10,10,10,10};
  else mPowers = {5,5,5,5};
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
// } else {
//     std::cout << "FSIFitFCN: Must select either Octave or TMultiDimFit!" << std::endl;
//     std::exit(-1);
// }
  std::cout<< std::endl << "FSIChi2Grid::InterpolateFiniteGrid()... done" << std::endl;
  
}

//***********************************************************************************
// Evaluate interpolated function at required point in n-dimensional space
//***********************************************************************************
double FSIChi2Grid::GetInterpolatedGridPoint(const std::vector<double> &x){

  return FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(multifit,x);

}

//***********************************************************************************
// Get point using interpn spline function of GNU-Octave
//***********************************************************************************
double FSIChi2Grid::GetSplinedGridPoint(const std::vector<double> &x){

  // Interpolate to 0.75 1.0 0.5
  Matrix *x_point[nFSIpars];
  const int nPointsToInterp = 1;
  for (int idim=0; idim<nFSIpars; idim++) {
    x_point[idim] = new Matrix(nPointsToInterp,1);
    (*x_point[idim])(0,0) = x[idim];
    // std::cout << "x_point:[ " << idim << "]: " <<  (*x_point[idim]) << std::endl;
  }
  // std::cout << "octave_value_list----------------" << std::endl;
  octave_value_list in;
  for (int idim=0; idim<nFSIpars; idim++) {
    in.append( octave_value( *x_vec[idim] ) );
  }

  in.append( octave_value (*v_vec) );
  
  for (int idim=0; idim<nFSIpars; idim++)
    in.append( octave_value( *x_point[idim] ) );

  in.append( octave_value("spline"));  // Default is linear
  octave_value_list out_v_point = feval ("interpn",in,nPointsToInterp);
  
  Matrix v_point = Matrix(1,1);
  if (!error_state && out_v_point.length () > 0) {
    v_point = out_v_point(0).matrix_value ();
  }
  double splined_chi2 = v_point(0);
  // std::cout << "splined_chi2  : " << splined_chi2 << std::endl;

  return splined_chi2;

}
