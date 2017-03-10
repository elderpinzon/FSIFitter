#include "ModelPrediction.hxx"

ModelPrediction::ModelPrediction(TString FileName){

  std::cout<<"ModelPrediction with filename: " << FileName << std::endl;
  TString FileNameWDir = "input/" + FileName;
  inFile = new TFile(FileNameWDir,"OPEN");

  ModelName = FileName.ReplaceAll("_xs.root","");
  ModelName = ModelName.ReplaceAll("_"," ");

  kColor = 1;
  kGeVtoMeV = false;

}

ModelPrediction::~ModelPrediction(){

  inFile->Close();

}

TH1F* ModelPrediction::GetCrossSection(TString nuclei, TString type, TString xsec){

  histo = (TH1F*)inFile->Get(Form("%s_%s_mom_%s_xs",type.Data(),nuclei.Data(),xsec.Data()));
  if(!histo){
    std::cout << "...didn't find: "
	      << Form("%s_%s_mom_%s_xs",type.Data(),nuclei.Data(),xsec.Data())
	      <<" exiting before crash" << std::endl;
    std::exit(1);
  }

  histo->SetLineColor(kColor);
  histo->Smooth();

  if(kGeVtoMeV){
    // std::cout << "GeVtoMeV for "
    // 	      << Form("%s_%s_mom_%s_xs",type.Data(),nuclei.Data(),xsec.Data())
    // 	      << ModelName << std::endl;
    histo = GeVtoMeV();
  }
  
  return histo;

}

Double_t ModelPrediction::GetCrossSectionValue(TString nuclei, TString type, TString xsec,Double_t momentum){

  // Ugly hack
  if(nuclei == "c")
    nuclei = "c12";

  // First get histogram
  GetCrossSection(nuclei, type, xsec);
  
  //Now "interpolate" to get the cross section
  Int_t nBins = histo->GetNbinsX();
  Double_t maxMom = histo->GetBinCenter(nBins);
  
  if(momentum > maxMom){
    std::cout << "ModelPrediction::GetCrossSectionValue(): requesting momentum out of range: "
	      << momentum << " " << maxMom
	      << ". Returning penultimate bin value."
	      << std::endl;
    return histo->GetBinContent(nBins-1);
  }
  
  //Find two bins surrounding required momentum
  Int_t bin1 = TMath::Floor(momentum/(maxMom/nBins));
  Int_t bin2 = TMath::Ceil(momentum/(maxMom/nBins));

  // Linear interpolation
  Double_t slope = (histo->GetBinContent(bin2) - histo->GetBinContent(bin1)) / (histo->GetBinCenter(bin2) - histo->GetBinCenter(bin1));

  Double_t value = histo->GetBinContent(bin1) + slope * (momentum - histo->GetBinCenter(bin1));

  return value;

}

TH1F* ModelPrediction::GeVtoMeV(){

  Int_t nbins = histo->GetNbinsX();
  TH1F *htmp = new TH1F(histo->GetTitle(),histo->GetName(),nbins,histo->GetBinLowEdge(1)*1000.,histo->GetBinLowEdge(nbins+1)*1000.);
  for(Int_t i = 0; i<nbins; i++){
    htmp->SetBinContent(i,histo->GetBinContent(i));
  }
  htmp->SetLineColor(kColor);
  return htmp;
}
