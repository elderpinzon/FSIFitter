#include "ModelPrediction.hxx"

ModelPrediction::ModelPrediction(TString FileName){

  std::cout<<"ModelPrediction with filename: " << FileName << std::endl;
  TString FileNameWDir = "input/" + FileName;
  inFile = new TFile(FileNameWDir,"OPEN");
  ModelName = FileName.ReplaceAll("_xs.root","");

  kColor = 1;
  kGeVtoMeV = false;

}

ModelPrediction::~ModelPrediction(){

  inFile->Close();

}

TH1D* ModelPrediction::GetCrossSection(TString nuclei, TString type, TString xsec){

  char *title = Form("%s_%s_mom_%s_xs",type.Data(),nuclei.Data(),xsec.Data());

  histo = (TH1D*)inFile->Get(title);

  if(!histo){
    title = Form("%s_%s_mom_%s",type.Data(),nuclei.Data(),xsec.Data());
    histo = (TH1D*)inFile->Get(title);
  }
  
  if(!histo){
    std::cout << ModelName << "...didn't find: "
	      << title
	      << std::endl;

    // If not found, return empty histogram
    histo = new TH1D(Form("%s_%s",ModelName.Data(),title),Form("%s_%s",ModelName.Data(),title),10,-10.0,0.0);
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

TGraph* ModelPrediction::GetCrossSectionGraph(TString nuclei, TString type, TString xsec){

  char *title = Form("%s_%s_mom_%s_xs",type.Data(),nuclei.Data(),xsec.Data());
  graph = (TGraph*)inFile->Get(title);
  if(!graph){
    std::cout << ModelName << "...didn't find: "
	      << title
	      << std::endl;

    // If not found, return empty graph
    graph = new TGraph(1);
  }
  else{
    std::cout << ModelName << " found TGraph! "
	      << title << std::endl;
  }
  
  graph->SetLineColor(kColor);

  return graph;

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

TH1D* ModelPrediction::GeVtoMeV(){

  Int_t nbins = histo->GetNbinsX();
  TH1D *htmp = new TH1D(histo->GetTitle(),histo->GetName(),nbins,histo->GetBinLowEdge(1)*1000.,histo->GetBinLowEdge(nbins+1)*1000.);
  for(Int_t i = 0; i<nbins; i++){
    htmp->SetBinContent(i,histo->GetBinContent(i));
  }
  htmp->SetLineColor(kColor);
  return htmp;
}
