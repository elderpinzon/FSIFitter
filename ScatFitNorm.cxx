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
  std::cout << "genWeightsFromNRooTracker_BANFFv2_2014.exe -p psyche_inputfile -r run_period -o outputfile" << std::endl;
}

// Messy way to process cmd line arguments.

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
    else if(std::string(argv[i]) == "-n") fitNorm = true;
    else {
      std::cout << "Invalid argument:" << argv[i] << " "<< argv[i+1] << std::endl;
      Usage();
      exit(1);
    }
  } 
}

int main(int argc, char* argv[]){

  ParseArgs(argc,argv);
  
  std::cout << "kOctave: " << kOctave
	    << " output: " << outputfileFit
	    << " datafit: " << dataFit
	    << " nFitPars: " << nFitPars
	    << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  if (nuclFit) std::cout << "Welcome to the nucleon scattering external data fitter - v2016" << std::endl;
  else std::cout << "Welcome to the pion scattering external data fitter - v2016" << std::endl;
  std::cout << "  (Now Including DUET data with correlation information!!)"    << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  if(kOctave) std::cout << "Octave" << std::endl;
  else std::cout << "TMultiDimFit" << std::endl;
  
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
  
  // Define new way of doing the chi2 "grid"
  aFSIChi2GridNorm = new FSIChi2GridNorm(AllDataSets);
  aFSIChi2GridNorm->LoadInterpolatedCrossSections();

  // Loop over all datasets, but skip DUET. Will do it later
  // Set parameter names for DUET datasets
  Int_t nNonDUETDataSets = NonDUETDataSets.size();
  
  // Number of paramateres in fit. 
  // FSI + #Datasets (other than DUET)
  Int_t npar = nFSIparsFitted + nNonDUETDataSets;
  
  // Define FCN function. 
  // The grid object is passed beacuse the FCN function will use the result 
  // of the interpolation, which is stored there
  FSIFitFCNNorm fcn(aFSIChi2GridNorm);
  fcn.SetNPARS(npar);

  // Define Minuit object
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);
  ROOT::Math::Functor f(fcn,npar);
  min.SetFunction(f);

  printf("Initialized TMinuit2 with %d FSI parameters and %d normalization parameters\n",
	 5,nNonDUETDataSets);

  // Starting point
  const double nomCenter[7] = {0.9,1.25,0.8,0.9,1.8,1.8,1.0};
  const double error[7] = {0.1,0.1,0.1,0.2,0.2,0.2,0.1};
    
  // FSI parameters
  for(int par = 0; par<5; par++){
    min.SetVariable(par,parNames[par].Data(),nomCenter[par],error[par]);
  }
  
  // Normalization parameters
  Double_t fitNormError = 0;
  if(fitNorm) fitNormError = 0.1;
  
  for (int par = nFSIparsFitted; par < npar; par++){
    min.SetVariable(par,Form("%s_%s_%s_%s",
			     NonDUETDataSets[par-5].GetNucleusString().Data(),
			     NonDUETDataSets[par-5].GetIntrTypeString().Data(),
			     NonDUETDataSets[par-5].GetPionTypeString().Data(),
			     NonDUETDataSets[par-5].GetSetName().Data()),1.0,fitNormError);

    
    printf("Norm par %d for data set: %s\n",par,NonDUETDataSets[par-5].GetSetName().Data());
  }
  
  
  min.Minimize();
  min.PrintResults();
  
  // Chi2 and NDOF
  nDOF = nDOF + 2*nNonDUETDataSets - npar;
  TVectorT<double> Chi2NDOF(3);
  Chi2NDOF[0] = min.MinValue();
  Chi2NDOF[1] = nDOF;
  Chi2NDOF[2] = (min.MinValue()/nDOF);

  // Save results as histos
  TH1D *fsiparsHist = new TH1D("fsiparsHist","FSI Parameters",nFSIparsFitted,0,nFSIparsFitted);
  TH1D *normparsHist = new TH1D("normparsHist","Normalization Parameters",nNonDUETDataSets,0,nNonDUETDataSets);
  

  // Grab best fit parameters and covariance
  const double *xs = min.X();
  TMatrixD *cov = new TMatrixD(npar,npar);
  TVectorT<double> BestFitPars(npar);

  for (int ipar=0; ipar<npar; ipar++) {
    printf("%.3f±%.3f\n", xs[ipar], sqrt(min.CovMatrix(ipar,ipar)));
    BestFitPars[ipar] = xs[ipar];

    for (int ipar2=0; ipar2<npar; ipar2++)
      (*cov)(ipar,ipar2) = min.CovMatrix(ipar,ipar2);
    
    // Store FSI parameters
    if(ipar < 5){
      fsiparsHist->GetXaxis()->SetBinLabel(ipar+1,parNames[ipar].Data());
      fsiparsHist->SetBinContent(ipar+1,xs[ipar]);
      fsiparsHist->SetBinError(ipar+1,sqrt(min.CovMatrix(ipar,ipar)));
    }
    // Normalization parameters
    else{
      normparsHist->GetXaxis()
	->SetBinLabel(ipar-5+1,Form("%s_%s_%s_%s",
				    NonDUETDataSets[ipar-5].GetNucleusString().Data(),
				    NonDUETDataSets[ipar-5].GetIntrTypeString().Data(),
				    NonDUETDataSets[ipar-5].GetPionTypeString().Data(),
				    NonDUETDataSets[ipar-5].GetSetName().Data()));
      normparsHist->SetBinContent(ipar-5+1,xs[ipar]);
      normparsHist->SetBinError(ipar-5+1,sqrt(min.CovMatrix(ipar,ipar)));
    }
  }

  // // Get the covariance matrix
  // TMatrixD *cov = new TMatrixD(npar,npar,minuit->GetCovarianceMatrix());

  // Calculate correlation
  TMatrixD corr = FSIFitterUtils::CovarianceToCorrelation(cov);

  // Store results in output file
  TFile *fout = new TFile(outputfileFit.c_str(),"RECREATE");
  Chi2NDOF.Write("Chi2NDOF");
  fsiparsHist->Write();
  normparsHist->Write();
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
  for(int par = 0; par<npar; par++){

    if(par < 5){
      cov_histo->GetXaxis()->SetBinLabel(par+1,parNames[par].Data());
      cov_histo->GetYaxis()->SetBinLabel(par+1,parNames[par].Data());
      corr_histo->GetXaxis()->SetBinLabel(par+1,parNames[par].Data());
      corr_histo->GetYaxis()->SetBinLabel(par+1,parNames[par].Data());
    }
    else{

      cov_histo->GetXaxis()
	->SetBinLabel(par+1,Form("%s_%s_%s_%s",
				 NonDUETDataSets[par-5].GetNucleusString().Data(),
				 NonDUETDataSets[par-5].GetIntrTypeString().Data(),
				 NonDUETDataSets[par-5].GetPionTypeString().Data(),
				 NonDUETDataSets[par-5].GetSetName().Data()));
      cov_histo->GetYaxis()
	->SetBinLabel(par+1,Form("%s_%s_%s_%s",
				 NonDUETDataSets[par-5].GetNucleusString().Data(),
				 NonDUETDataSets[par-5].GetIntrTypeString().Data(),
				 NonDUETDataSets[par-5].GetPionTypeString().Data(),
				 NonDUETDataSets[par-5].GetSetName().Data()));
      corr_histo->GetXaxis()
	->SetBinLabel(par+1,Form("%s_%s_%s_%s",
				 NonDUETDataSets[par-5].GetNucleusString().Data(),
				 NonDUETDataSets[par-5].GetIntrTypeString().Data(),
				 NonDUETDataSets[par-5].GetPionTypeString().Data(),
				 NonDUETDataSets[par-5].GetSetName().Data()));
      corr_histo->GetYaxis()
	->SetBinLabel(par+1,Form("%s_%s_%s_%s",
				 NonDUETDataSets[par-5].GetNucleusString().Data(),
				 NonDUETDataSets[par-5].GetIntrTypeString().Data(),
				 NonDUETDataSets[par-5].GetPionTypeString().Data(),
				 NonDUETDataSets[par-5].GetSetName().Data()));
    }
  }
  
  cov_histo->Write("cov_histo");
  corr_histo->Write("corr_histo");
  fout->Close();

  return 0;

  // // Run fit many times with random starting locations
  // npar = nFSIpars;
  // // Histograms to store bestfit points
  // TH1D *hParsInitial[npar];
  // TH1D *hParsFinal[npar];
  // for(int i = 0; i<npar; i++){
  //   hParsInitial[i] = new TH1D(Form("initial_%s",parNames[i].Data()),Form("%s;Parameter Initial Value",parNames[i].Data()),100,FSIParsMin[i],FSIParsMax[i]);
  //   hParsFinal[i] = new TH1D(Form("final_%s",parNames[i].Data()),Form("%s;Parameter Best Fit",parNames[i].Data()),100,BestFitPars_out[i]*0.9,BestFitPars_out[i]*1.1);
  // }
  
  // // Loop on fittintg tries
  // int nTotTries = 0;//100;
  // for(int nTries = 0; nTries < nTotTries; nTries++){
    
  //   std::cout << "nTries: " << nTries << std::endl;
  //   npar = nFSIpars;
  //   // New TMinuit instance
  //   TFitterMinuit *aminuit = new TFitterMinuit(nFSIpars);
  //   aminuit->SetMinuitFCN(&fcn);

  //   // Random number generator
  //   TRandom3 *gRan = new TRandom3(0);
  //   gRan->SetSeed(0);
    
  //   // Array to store starting position for this try
  //   Double_t tryStart[7];
  //   const double error[7] = {0.1,0.1,0.1,0.2,0.2,0.2,0.1};
    
  //   for(int k = 0; k<7; k++){
      
  //     // // For gaussian throw
  //     // Double_t sigma = 0.2;
  //     // Double_t random = gRan->Gaus(1,sigma);
  //     // if(k<5) tryStart[k] = nomPars[k]* random;
  //     // else tryStart[k] = nomPars[k];

  //     // // For uniform throw in all parameter range
  //     // Double_t random =  gRan->Uniform();
  //     // if(k<5) tryStart[k] = FSIParsMin[k] + (FSIParsMax[k]- FSIParsMin[k]) * random;
  //     // else tryStart[k] = nomPars[k];

  //     // For uniform throw in half parameter range
  //     Double_t random =  gRan->Uniform(-1,1);
  //     if(k<5) tryStart[k] = nomPars[k] + (FSIParsMax[k]- FSIParsMin[k]) * random/2.;
  //     else tryStart[k] = nomPars[k];
      
  //     std::cout << "k: " << k 
  // 		<< " gRan: " << random
  // 		<< " " << tryStart[k]
  // 		<< std::endl;
  //   }
    
  //   for(int par = 0; par<nFSIpars; par++){
  //     if(!kOctave && par == 5) continue; // skip FEFCXH for TMultiDimFit
  //     aminuit->SetParameter(par,parNames[par].Data(),tryStart[par],error[par],0,0);
  //     hParsInitial[par]->Fill(tryStart[par]);
  //   }

  //   if(!kOctave && nFitPars == 5){aminuit->FixParameter(5); npar--;} //fix FEFALL

  //   if(kOctave){
  //     aminuit->FixParameter(5); npar--; //fix FEFCXH always for octave
  //     if(nFitPars == 5) aminuit->FixParameter(6); npar--; //fix FEFALL
  //   }
    
  //   // Actually create and run minimizer
  //   aminuit->CreateMinimizer();
  //   int iret = aminuit->Minimize();//0,100.0);
  //   if (iret != 0) {
  //     return iret;
  //   }
    
  //   // // Grab best fit parameters
  //   // Double_t *BestFitPars_out = new Double_t[npar];
  //   for (int ipar=0; ipar<npar; ipar++) {
  //   //   BestFitPars_out[ipar] = minuit->GetParameter(ipar);
  //     hParsFinal[ipar]->Fill(aminuit->GetParameter(ipar));
  //   }

  //   delete gRan;
  //   //delete aminuit;

  // }

  // TFile *fout2 = new TFile("manyfits_histos.root","RECREATE");
  // //fout->cd();
  // for(int k = 0; k<npar; k++){
  //   std::cout << "Writing par: " << k << std::endl;
  //   hParsFinal[k]->Print();
  //   hParsInitial[k]->Write(Form("initial_%d",k));
  //   hParsFinal[k]->Write(Form("final_%d",k));
  // }
  
  //fout->Close();
  
  std::cout << "\nMinimization done for ... "
	    << AllDataSets.size() << " datasets."
	    << std::endl;
    
  std::cout<<" Good Bye!" << std::endl;
    
  return 0;
}
