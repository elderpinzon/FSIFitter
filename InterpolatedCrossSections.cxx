//***********************************************************************************
// InterpolatedCrossSections.cxx
// Interpolates cross sections from an input grid scan using TMultiDimFit
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "InterpolatedCrossSections.hxx"

InterpolatedCrossSections::InterpolatedCrossSections(std::string fFileName) : filename(fFileName)
{

  std::cout << "\nWill interpolate cross sections from file: "<< filename << std::endl;
  
}

InterpolatedCrossSections::~InterpolatedCrossSections(){

  for(int i = 0; i<5; i++)
    delete xsecmultidimfit[i];
 
}

void InterpolatedCrossSections::RunMultiDimFitInterpolations(){

  // Load results from FSI parameter scan
  FSIParameterScan aFSIParameterScan(filename);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;

  // Initialize MultiDimFit objects
  for(int i = 0; i<5; i++){
    xsecmultidimfit[i] = new TMultiDimFit(6,TMultiDimFit::kMonomials,"kv");
    //Int_t mPowers[] = {3,3,3,3,3,3,3};
    Int_t mPowers[] = {5,5,5,5,5,5,5};
    //Int_t mPowers[] = {7,7,7,7,7,7,7};
    //Int_t mPowers[] = {9,9,9,9,9,9,9};
    xsecmultidimfit[i]->SetMaxPowers(mPowers);
    // xsecmultidimfit[i]->SetMaxFunctions(10000);
    // xsecmultidimfit[i]->SetMaxStudy(500000);
    // xsecmultidimfit[i]->SetMaxTerms(10000);
    // xsecmultidimfit[i]->SetMaxFunctions(10);
    // xsecmultidimfit[i]->SetMaxStudy(10);
    // xsecmultidimfit[i]->SetMaxTerms(20);
    xsecmultidimfit[i]->SetMaxFunctions(100);
    xsecmultidimfit[i]->SetMaxStudy(10000);
    xsecmultidimfit[i]->SetMaxTerms(50);
    xsecmultidimfit[i]->SetPowerLimit(1);
    xsecmultidimfit[i]->SetMinRelativeError(0.001);

  }

  if (aFSIParameterScan.fChain == 0) return;

  // Loop over intput tree and add entries to MultiDimFit objects
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;

    aFSIParameterScan.GetEntry(jentry);

    // Only do pip-carbon, for FEFCXH = 1.8
    if(TMath::Abs(aFSIParameterScan.FEFCXH - 1.8) > 0.01) continue;
    if(TMath::Abs(aFSIParameterScan.FEFALL - 1.0) > 0.01) continue;
    if(aFSIParameterScan.target != 6) continue;
    if(aFSIParameterScan.pid != 211) continue;
    if(aFSIParameterScan.mom > 500) continue;
    //if(aFSIParameterScan.mom < 400) continue;
    
    // Interval variables to store input variables below
    //Double_t fsi_pars[7];
    //Double_t kFSIPars[6];

    Double_t fsi_pars[6];
    Double_t kFSIPars[5];

    fsi_pars[0] = aFSIParameterScan.mom;
    fsi_pars[1] = aFSIParameterScan.FEFQE;
    fsi_pars[2] = aFSIParameterScan.FEFABS;
    fsi_pars[3] = aFSIParameterScan.FEFCX;
    fsi_pars[4] = aFSIParameterScan.FEFINEL;
    fsi_pars[5] = aFSIParameterScan.FEFQEH;
    //fsi_pars[6] = aFSIParameterScan.FEFALL;
    //fsi_pars[6] = aFSIParameterScan.FEFCXH;
    //fsi_pars[7] = aFSIParameterScan.FEFALL;

    // IF LE point set FSIPars to varied values of HE pars
    if(aFSIParameterScan.mom < 400){

      //FEFINEL
      for(kFSIPars[3] = FSIParsMin[3]; kFSIPars[3] <= FSIParsMax[3]; kFSIPars[3] = kFSIPars[3] + 0.2){
	//FEFQEH
	for(kFSIPars[4] = FSIParsMin[4]; kFSIPars[4] <= FSIParsMax[4]; kFSIPars[4] = kFSIPars[4] + 0.2){
	  //FEFCXH
	  //for(kFSIPars[5] = FSIParsMin[5]; kFSIPars[5] <= FSIParsMax[5]; kFSIPars[5] = kFSIPars[5] + 0.2){
	    
	  fsi_pars[4] = kFSIPars[3];
	  fsi_pars[5] = kFSIPars[4];
	  //fsi_pars[6] = kFSIPars[5];
	    
	  xsecmultidimfit[0]->AddRow(fsi_pars,aFSIParameterScan.xqe);
	  xsecmultidimfit[1]->AddRow(fsi_pars,aFSIParameterScan.xabs);
	  xsecmultidimfit[2]->AddRow(fsi_pars,aFSIParameterScan.xcx);
	  xsecmultidimfit[3]->AddRow(fsi_pars,aFSIParameterScan.xdcx);
	  xsecmultidimfit[4]->AddRow(fsi_pars,aFSIParameterScan.xhadr);
	    
	  //}
	}
      }
    }
    
    
    // IF HE point set Facts to varied values of LE pars
    else if(aFSIParameterScan.mom > 500){
      //FEFQE
      for(kFSIPars[0] = FSIParsMin[0]; kFSIPars[0] <= FSIParsMax[0]; kFSIPars[0] = kFSIPars[0] + 0.1){
	//FEABS
	for(kFSIPars[1] = FSIParsMin[1]; kFSIPars[1] <= FSIParsMax[1]; kFSIPars[1] = kFSIPars[1] + 0.1){
	  //FEFCX
	  for(kFSIPars[2] = FSIParsMin[2]; kFSIPars[2] <= FSIParsMax[2]; kFSIPars[2] = kFSIPars[2] + 0.1){
	    
	    fsi_pars[1] = kFSIPars[0];
            fsi_pars[2] = kFSIPars[1];
            fsi_pars[3] = kFSIPars[2];

	    xsecmultidimfit[0]->AddRow(fsi_pars,aFSIParameterScan.xqe);
	    xsecmultidimfit[1]->AddRow(fsi_pars,aFSIParameterScan.xabs);
	    xsecmultidimfit[2]->AddRow(fsi_pars,aFSIParameterScan.xcx);
	    xsecmultidimfit[3]->AddRow(fsi_pars,aFSIParameterScan.xdcx);
	    xsecmultidimfit[4]->AddRow(fsi_pars,aFSIParameterScan.xhadr);

	    
	  }
	}
      }
    }
    
    // If LE-HE transition no need to vary manually
    else{

      xsecmultidimfit[0]->AddRow(fsi_pars,aFSIParameterScan.xqe);
      xsecmultidimfit[1]->AddRow(fsi_pars,aFSIParameterScan.xabs);
      xsecmultidimfit[2]->AddRow(fsi_pars,aFSIParameterScan.xcx);
      xsecmultidimfit[3]->AddRow(fsi_pars,aFSIParameterScan.xdcx);
      xsecmultidimfit[4]->AddRow(fsi_pars,aFSIParameterScan.xhadr);
    }
    
  }

  //TFile *fout = new TFile("output/ResultCrossSectionInterpolation.root","RECREATE");
  
  for(int i = 0; i<5; i++){
    
    if(i!=1) continue;
    
    TStopwatch s;
    s.Start();

    std::cout << "Interpolating xsec = " << IndivXsecs[i] << std::endl;    
    // Print out the statistics of input sample
    xsecmultidimfit[i]->Print("s");

    // Run interpolation
    //xsecmultidimfit[i]->MakeHistograms();
    xsecmultidimfit[i]->FindParameterization("");
    //xsecmultidimfit[i]->Fit("M");

    // Print and store interpolation result
    xsecmultidimfit[i]->Print("prcfkm");
    // xsecmultidimfit[i]->MakeCode();
    xsecmultidimfit[i]->MakeMethod(Form("%s",IndivXsecs[i].Data()));
    //xsecmultidimfit[i]->Write(Form("interp_%s",IndivXsecs[i].Data()));

    std::cout<< std::endl << "... interpolation done in " 
	   << s.RealTime() << " seconds." << std::endl;

  }
  
  //fout->Close();
  
}

double InterpolatedCrossSections::GetInterpolatedCrossSection(int xs_index, const std::vector<double> &x){
  
  return FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(xsecmultidimfit[xs_index],x);
  
}
