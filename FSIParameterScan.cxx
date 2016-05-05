//***********************************************************************************
// FSIParameterScan.cxx
// Class to read and manipulate the result of the MC scans
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#define FSIParameterScan_cxx
#include "FSIParameterScan.hxx"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// To be deletedx
#include <iostream>
#include <TMultiDimFit.h>

//***********************************************************************************
// Parse the arguments for running the fit
//***********************************************************************************
void FSIParameterScan::DetermineGridSize(){
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Double_t fsi_pars[4];
  for(int par=0; par<3; par++){
    fsi_par_low[par]   =  9999;
    fsi_par_high[par]  = -9999;
    fsi_par_step[par]  =  0;
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    fsi_pars[0] = mom;
    fsi_pars[1] = qe;
    fsi_pars[2] = abs;
    fsi_pars[3] = cx;
    
    for(int par=0; par<3; par++){
      if(fsi_pars[par]<fsi_par_low[par]) fsi_par_low[par] = fsi_pars[par];
      if(fsi_pars[par]>fsi_par_high[par]) fsi_par_high[par] = fsi_pars[par];
    }
  }
}

//***********************************************************************************
// A simple test loop that prints the 10 first and last entries for debugging purposes
//***********************************************************************************
void FSIParameterScan::Test(){

  if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry<10 || (nentries-jentry)<10)
	std::cout<<"jentry: " <<jentry << " mom: "<<mom<<" qe: "<<qe<<" abs: "<<abs<<" cx: "<<cx<<std::endl;
    }
}

//***********************************************************************************
// Try to interpolate the xsecs in 4D (fsipars + mom)
// This is not used
//***********************************************************************************

void FSIParameterScan::InterpolateXsecs()
{

  // One for each xsec: xqe, xabs, xcx, xdcx, xhadr
  TMultiDimFit *multifit[5];
  std::vector< std::string > xsecs = {"xqe","xabs","xcx","xdcx","xhadr"};
  
  for(int i = 0; i<1; i++){

    std::cout << "Interpolating xsec = " << xsecs[i] << std::endl;
    
    // Create MultiDim fitter
    multifit[i] = new TMultiDimFit(4,TMultiDimFit::kMonomials,"kv");
    Int_t mPowers[] = {10,10,10,10};
    multifit[i]->SetMaxPowers(mPowers);
    multifit[i]->SetMaxFunctions(1000);
    multifit[i]->SetMaxStudy(10000);
    ////multifit[i]->SetMaxStudy(1000);
    multifit[i]->SetMaxTerms(100);
    multifit[i]->SetPowerLimit(10);
    Double_t fsi_pars[4];
  
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      fsi_pars[0] = mom;
      fsi_pars[1] = qe;
      fsi_pars[2] = abs;
      fsi_pars[3] = cx;
      if(i == 0) multifit[i]->AddRow(fsi_pars,xqe);
      if(i == 1) multifit[i]->AddRow(fsi_pars,xabs);
      if(i == 2) multifit[i]->AddRow(fsi_pars,xcx);
      if(i == 3) multifit[i]->AddRow(fsi_pars,xdcx);
      if(i == 4) multifit[i]->AddRow(fsi_pars,xhadr);
    }
    
    // Print out the statistics of input sample
    multifit[i]->Print("s");
    
    // Run interpolation
    multifit[i]->MakeHistograms();
    multifit[i]->FindParameterization("");
    //multifit[i]->Fit("M");
    
    // Print and store interpolation result
    //multifit[i]->Print("rc");
    multifit[i]->Print("prcfkm");
    //multifit[i]->MakeCode();
    multifit[i]->MakeMethod(Form("%s",xsecs[i].c_str()));
    
    /*
    TFile *fout = new TFile("interpolation_covariance.root","RECREATE");
    const TMatrixD *inter_corr = multifit[i]->GetCorrelationMatrix();
    inter_corr->Write("inter_corr");
    inter_corr->Print();*/
    std::cout << "Chi2 " << xsecs[i] << ": " << multifit[i]->GetChi2() << std::endl;
    multifit[i]->Clear();
    //fout->Close();
  }
}
