#ifndef MODELPREDICTION_HXX_
#define MODELPREDICTION_HXX_

#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TMath.h>

class ModelPrediction{

private:

  TString ModelName;
  TFile *inFile;
  TH1F* histo;
  Int_t kColor;
  bool kGeVtoMeV;

public:

  ModelPrediction(){};
  ModelPrediction(TString FileName);
  ~ModelPrediction();
  TString GetModelName(){return ModelName;};
  TH1F* GetCrossSection(TString nuclei, TString type, TString xsec);
  Double_t GetCrossSectionValue(TString nuclei, TString type, TString xsec,Double_t momentum);
  void SetLineColor(Int_t color){kColor = color;};
  void SetGeVtoMeV(){kGeVtoMeV = true;};
  TH1F* GeVtoMeV();
};

#endif
