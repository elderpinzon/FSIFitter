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

void InterpolatedCrossSections::RunMultiDimFitInterpolations(){

  // Load results from FSI parameter scan
  FSIParameterScan aFSIParameterScan(filename);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;

  Double_t fsi_pars[4];
  
  // Initialize MultiDimFit objects
  for(int i = 0; i<5; i++){
    xsecmultidimfit[i] = new TMultiDimFit(4,TMultiDimFit::kMonomials,"k");
    //Int_t mPowers[] = {10,10,10,10};
    Int_t mPowers[] = {5,5,5,5};
    xsecmultidimfit[i]->SetMaxPowers(mPowers);
    xsecmultidimfit[i]->SetMaxFunctions(1000);
    xsecmultidimfit[i]->SetMaxStudy(10000);
    xsecmultidimfit[i]->SetMaxTerms(100);
    xsecmultidimfit[i]->SetPowerLimit(10);
  }

  if (aFSIParameterScan.fChain == 0) return;

  // Loop over intput tree and add entries to MultiDimFit objects
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;

    fsi_pars[0] = aFSIParameterScan.mom;
    fsi_pars[1] = aFSIParameterScan.qe;
    fsi_pars[2] = aFSIParameterScan.abs;
    fsi_pars[3] = aFSIParameterScan.cx;

    xsecmultidimfit[0]->AddRow(fsi_pars,aFSIParameterScan.xqe);
    xsecmultidimfit[1]->AddRow(fsi_pars,aFSIParameterScan.xabs);
    xsecmultidimfit[2]->AddRow(fsi_pars,aFSIParameterScan.xcx);
    xsecmultidimfit[3]->AddRow(fsi_pars,aFSIParameterScan.xdcx);
    xsecmultidimfit[4]->AddRow(fsi_pars,aFSIParameterScan.xhadr);

  }

  for(int i = 0; i<5; i++){

    std::cout << "Interpolating xsec = " << xsecs[i] << std::endl;    
    // Print out the statistics of input sample
    //xsecmultidimfit[i]->Print("s");

    // Run interpolation
    //xsecmultidimfit[i]->MakeHistograms();
    xsecmultidimfit[i]->FindParameterization("");
    //xsecmultidimfit[i]->Fit("M");

    // Print and store interpolation result
    //xsecmultidimfit[i]->Print("rc");
    //xsecmultidimfit[i]->Print("prcfkm");
    //xsecmultidimfit[i]->MakeCode();
    //xsecmultidimfit[i]->MakeMethod(Form("%s",xsecs[i].c_str()));
    //std::cout << "Chi2 " << xsecs[i] << ": " << xsecmultidimfit[i]->GetChi2() << std::endl;

  }
}

double InterpolatedCrossSections::GetInterpolatedCrossSection(int xs_index, const std::vector<double> &x){

  return FSIFitterUtils::BuildInterpolatedFunctionAndEvaluate(xsecmultidimfit[xs_index],x);

}
