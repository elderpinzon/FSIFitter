#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TGraphErrors.h>

#include <FSIFitterUtils.hxx>

class TN032Envelopes{

private:

  // Histos for TN-032 plots
  TH1D *nominal;
  TH1D *varied;
  TGraphErrors *gr_nom;
  TGraphErrors *gr_varied;
  TGraphErrors *gr_min;
  TGraphErrors *gr_max;
  TGraphErrors *gr_shade;
  Double_t kMaxScale;

public:

  TN032Envelopes(){};
  TN032Envelopes(TString nuclei, Int_t pionType, TString iTypeString);
  TGraph* GetCrossSectionNominal(){return gr_nom;};
  TGraph* GetCrossSectionShade(){ return gr_shade;};
  TGraph* GetCrossSectionMin(){return gr_min;};
  TGraph* GetCrossSectionMax(){return gr_max;};
  Double_t GetMaxScale(){return kMaxScale;};
  
};
