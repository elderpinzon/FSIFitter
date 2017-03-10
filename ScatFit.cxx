//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "ScatFit.hxx"

// Print the cmd line syntax
void Usage(){
  std::cout << "Cmd line syntax should be:" << std::endl;
  std::cout << "./ScatFit.exe -d index_datasets -o outputfile  (OPTIONAL for octave:  -t 0)" << std::endl;
}


//***********************************************************************************
// Parse the arguments for running the fit
//***********************************************************************************
void ParseArgs(int argc, char **argv){

  std::cout << "Parsing arguments ";
  
  int nargs = 1; 
  if(argc<(nargs*2+1)){ Usage(); exit(1); }
  for(int i = 1; i < argc; i+=2){
    if(std::string(argv[i]) == "-t") kOctave = true;
    else if(std::string(argv[i]) == "-o") outputfileFit = (argv[i+1]);
    else if(std::string(argv[i]) == "-d") dataFit = atoi(argv[i+1]);
    else if(std::string(argv[i]) == "-f") nFitPars = atoi(argv[i+1]);
    else {  
      std::cout << "Invalid argument:" << argv[i] << " "<< argv[i+1] << std::endl;
      Usage();
      exit(1);
    }
  } 
}

int main(int argc, char* argv[]){

  ParseArgs(argc,argv);
  
  std::cout << "-----------------------------------------------------------" << std::endl;
  if (nuclFit) std::cout << "Welcome to the nucleon scattering external data fitter - v2016" << std::endl;
  else std::cout << "Welcome to the pion scattering external data fitter - v2016" << std::endl;
  std::cout << "  (Now Including DUET data with correlation information!!)"    << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  std::cout << "kOctave: " << kOctave
	    << " output: " << outputfileFit
	    << " datafit: " << dataFit
	    << " nFitPars: " << nFitPars
	    << std::endl;
  
  if(kOctave) std::cout << "Octave" << std::endl;
  else std::cout << "TMultiDimFit" << std::endl;
  
  // Set data sets to be used. Will be added to AllDataSets vector
  if(dataFit == 0){
    AddDataSets("c_");
  }
  
  if(dataFit == 3){
    AddDataSets("o_");
  }
  
  if(dataFit == 4){
    AddDataSets("al_");
  }
  
  if(dataFit == 5){
    AddDataSets("fe_");
  }
  
  if(dataFit == 6){
    AddDataSets("cu_");
  }
  
  if(dataFit == 7){
    AddDataSets("pb_");
  }
  
  if(dataFit == 8){
    AddDataSets("fe_");
    AddDataSets("cu_");
    AddDataSets("pb_");
  }
  
  if(dataFit == 1){
    AddDataSets("c_");
    AddDataSets("o_");
    AddDataSets("al_");
  }
  
  if(dataFit == 2){
    AddDataSets("c_");
    AddDataSets("o_");
    AddDataSets("al_");
    AddDataSets("fe_");
    AddDataSets("cu_");
    AddDataSets("pb_");
  }
      
  // Define Chi2 Grid class from datasets
  aFSIChi2Grid = new FSIChi2Grid(AllDataSets);
  
  // Set interpolation method
  if(kOctave)
      aFSIChi2Grid->SetOctaveInterpolation();
  else
    aFSIChi2Grid->SetMultiDimFitInterpolation();
  
  // Build "made-up" covariance
  aFSIChi2Grid->BuildCovariance();

  // Grid scan to be used
  aFSIChi2Grid->SetGridScanFileName(ScanFileName);

  // Calculating model chi-squares  
  Bool_t calcModels = false;
  if(calcModels){
    
    ModelPrediction::ModelPrediction geant4("geant4_xs.root");
    aFSIChi2Grid->CalculateModelChiSquare(geant4);
    
    ModelPrediction::ModelPrediction genie_hA("genie_hA_xs.root");
    genie_hA.SetGeVtoMeV();
    aFSIChi2Grid->CalculateModelChiSquare(genie_hA);
    
    ModelPrediction::ModelPrediction genie_hA2014("genie_hA2014_xs.root");
    genie_hA2014.SetGeVtoMeV();
    aFSIChi2Grid->CalculateModelChiSquare(genie_hA2014);
    
    ModelPrediction::ModelPrediction genie_hN2015("genie_hN2015.root");
    genie_hN2015.SetGeVtoMeV();
    aFSIChi2Grid->CalculateModelChiSquare(genie_hN2015);
    
    ModelPrediction::ModelPrediction nuwro("nuwro_scattering_xs.root");
    nuwro.SetGeVtoMeV();
    aFSIChi2Grid->CalculateModelChiSquare(nuwro);
    
    ModelPrediction::ModelPrediction fluka("FLUKA_pipC.root");
    fluka.SetGeVtoMeV();
    aFSIChi2Grid->CalculateModelChiSquare(fluka);

    return 0;
  }
  
  // Interpolate finite grid
  aFSIChi2Grid->InterpolateFiniteGrid();

  // Plot slices of the Chi2 grid for visualization (and debugging)
  Bool_t plotGridSlices = false;
  if(plotGridSlices){

    // Since the goal is to compare the interpolations. Enable both
    aFSIChi2Grid->SetOctaveInterpolation();
    aFSIChi2Grid->SetMultiDimFitInterpolation();
    
    // Slices around these points
    std::vector<Double_t>x[8];
    x[0] = {0.9,1.25,0.9,0.8,1.8,1.8,1.1};
    x[1] = {0.8,1.05,0.8,0.6,1.4,1.8,1.3};
    x[2] = {1.1,1.05,1.1,1.2,1.4,1.8,0.9};

    for(Int_t i = 0; i < 1; i++)
      aFSIChi2Grid->PlotChiSquare(x[i],i);
    return 0;
  }
  
  // Alright now run the fit
  std::cout << "\nWill now run TMinuit2 using the interpolated grid" << std::endl;
  
  // Define FCN function. 
  // The grid object is passed beacuse the FCN function will use the result 
  // of the interpolation, which is stored there
  FSIFitFCN fcn(aFSIChi2Grid);
  fcn.SetNPARS(nFSIparsFitted);

  if(kOctave)
    fcn.SetOctaveInterpolation();
  else
    fcn.SetMultiDimFitInterpolation();

  // Define Minuit object
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);
  ROOT::Math::Functor f(fcn,nFSIparsFitted);
  min.SetFunction(f);
  
  // Starting point
  const double nomCenter[7] = {0.9,1.25,0.8,0.9,1.8,1.8,1.0};
  const double error[7] = {0.1,0.1,0.1,0.2,0.2,0.2,0.1};

  // FSI parameters
  for(int par = 0; par<nFSIparsFitted; par++){
    min.SetVariable(par,parNames[par].Data(),nomCenter[par],error[par]);
  }
  
  // Actually create and run minimizer
  min.Minimize();
  min.PrintResults();

  // Chi2 and NDOF
  Int_t nNonDUETDataSets = NonDUETDataSets.size();
  nDOF = nDOF + 2*nNonDUETDataSets - nFSIparsFitted;
  TVectorT<double> Chi2NDOF(3);
  Chi2NDOF[0] = min.MinValue();
  Chi2NDOF[1] = nDOF;
  Chi2NDOF[2] = (min.MinValue()/nDOF);

  // Save results as histos
  TH1D *fsiparsHist = new TH1D("fsiparsHist","FSI Parameters",nFSIparsFitted,0,nFSIparsFitted);

  // Grab best fit parameters and covariance
  const double *xs = min.X();
  TMatrixD *cov = new TMatrixD(nFSIparsFitted,nFSIparsFitted);
  TVectorT<double> BestFitPars(nFSIparsFitted);

  for (int ipar=0; ipar<nFSIparsFitted; ipar++) {
    printf("%.3f±%.3f\n", xs[ipar], sqrt(min.CovMatrix(ipar,ipar)));
    BestFitPars[ipar] = xs[ipar];

    for (int ipar2=0; ipar2<nFSIparsFitted; ipar2++)
      (*cov)(ipar,ipar2) = min.CovMatrix(ipar,ipar2);
    
    fsiparsHist->GetXaxis()->SetBinLabel(ipar+1,parNames[ipar].Data());
    fsiparsHist->SetBinContent(ipar+1,xs[ipar]);
    fsiparsHist->SetBinError(ipar+1,sqrt(min.CovMatrix(ipar,ipar)));
  }

  // Calculate correlation
  TMatrixD corr = FSIFitterUtils::CovarianceToCorrelation(cov);

  // Store results in output file
  TFile *fout = new TFile(outputfileFit.c_str(),"RECREATE");
  Chi2NDOF.Write("Chi2NDOF");
  fsiparsHist->Write();
  BestFitPars.Write("BestFitPars");
  cov->Write("cov");
  corr.Write("corr");
  
  // Also store as TH2D to set bin labels
  TH2D* cov_histo = new TH2D(*cov);
  cov_histo->SetName("cov_histo");
  cov_histo->SetMinimum(-cov_histo->GetMaximum());
  TH2D* corr_histo = new TH2D(corr);
  corr_histo->SetMinimum(-1.0);
  
  // Set bin labels for TH2D matrices
  for(int par = 0; par<nFSIparsFitted; par++){
    
    cov_histo->GetXaxis()->SetBinLabel(par+1,parNames[par].Data());
    cov_histo->GetYaxis()->SetBinLabel(par+1,parNames[par].Data());
    corr_histo->GetXaxis()->SetBinLabel(par+1,parNames[par].Data());
    corr_histo->GetYaxis()->SetBinLabel(par+1,parNames[par].Data());
  }
  
  cov_histo->Write("cov_histo");
  corr_histo->Write("corr_histo");
  
  fout->Close();
  
  // Run fit many times with random starting locations
  Bool_t runManyTimes = false;
  if(runManyTimes){
    
    // Histograms to store bestfit points
    TH1D *hParsInitial[nFSIparsFitted];
    TH1D *hParsFinal[nFSIparsFitted];
    for(int i = 0; i<nFSIparsFitted; i++){
      hParsInitial[i] = new TH1D(Form("initial_%s",parNames[i].Data()),Form("%s;Parameter Initial Value",parNames[i].Data()),100,FSIParsMin[i],FSIParsMax[i]);
      hParsFinal[i] = new TH1D(Form("final_%s",parNames[i].Data()),Form("%s;Parameter Best Fit",parNames[i].Data()),100,BestFitPars[i]*0.9,BestFitPars[i]*1.1);
    }
  
    // Loop on fitting tries
    int nTotTries = 100;
    for(int nTries = 0; nTries < nTotTries; nTries++){
      
      std::cout << "Attempt #" << nTries << std::endl;
      
      // New Minuit instance
      ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
      min.SetMaxFunctionCalls(1000000);
      min.SetMaxIterations(100000);
      min.SetTolerance(0.001);
      ROOT::Math::Functor f(fcn,nFSIparsFitted);
      min.SetFunction(f);
      
      // Random number generator
      TRandom3 *gRan = new TRandom3(0);
      gRan->SetSeed(0);
      
      // Array to store starting position for this try
      Double_t tryStart[7];
      const double error[7] = {0.1,0.1,0.1,0.2,0.2,0.2,0.1};
      
      for(int k = 0; k<7; k++){
		
	// For uniform throw in half parameter range
	Double_t random =  gRan->Uniform(-1,1);
	if(k<5) tryStart[k] = nomPars[k] + (FSIParsMax[k]- FSIParsMin[k]) * random/2.;
	else tryStart[k] = nomPars[k];
	
	std::cout << "k: " << k 
		  << " gRan: " << random
		  << " " << tryStart[k]
		  << std::endl;
      }
      
      for(int par = 0; par<nFSIparsFitted; par++){
	min.SetVariable(par,parNames[par].Data(),tryStart[par],error[par]);
	hParsInitial[par]->Fill(tryStart[par]);
      }
      
      
      // Actually create and run minimizer
      min.Minimize();
      min.PrintResults();

      // Grab best fit parameters
      const double *xs = min.X();
      TVectorT<double> BestFitPars(nFSIparsFitted);
      
      for (int ipar=0; ipar<nFSIparsFitted; ipar++){
	hParsFinal[ipar]->Fill(xs[ipar]);
      }
      
      delete gRan;
      
    }
    
    // Store output histograms to output file
    fout->cd();
    for(int k = 0; k<nFSIparsFitted; k++){
      std::cout << "Writing par: " << k << std::endl;
      hParsFinal[k]->Print();
      hParsInitial[k]->Write(Form("initial_%d",k));
      hParsFinal[k]->Write(Form("final_%d",k));
    }
    
  }
  
  fout->Close();
  
  std::cout << "\nMinimization done for ... "
	    << AllDataSets.size() << " datasets."
	    << std::endl;
  
  std::cout<<" Good Bye!" << std::endl;
  
  return 0;
}
