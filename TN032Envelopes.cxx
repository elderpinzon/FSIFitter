#include "TN032Envelopes.hxx"

//***********************************************************************************
// Builds TN032 1sigma envelopes using files in input folder
//***********************************************************************************
TN032Envelopes::TN032Envelopes(TString Nuclei, Int_t pionType, TString iTypeString){

  bool kVERBOSE = false;

  std::cout << "TN032Envelopes for nuclei: " << Nuclei 
	    << " type: " << pionType << " xsec: " << iTypeString << std::endl;

  // Loop through sets from TN-032, but first initialize nominal (set==0)
  TFile *fset = TFile::Open(Form("input/tn032_plotfiles/all_%s_%d_%d_plots.root",Nuclei.Data(),0,pionType),"OPEN");
  nominal = (TH1D*)fset->Get(Form("h%s",iTypeString.Data()));
  nominal->Print();
  gr_nom = new TGraphErrors(nominal);
  gr_min = new TGraphErrors(nominal);
  gr_max = new TGraphErrors(nominal);
  gr_shade = new TGraphErrors(gr_nom->GetN()*2);
  kMaxScale = 0;
  
  for(int set = 16; set>0; set--){

    if(kVERBOSE) std::cout<<"set: " << set << std::endl;	  
    TFile *fset = new TFile(Form("input/tn032_plotfiles/all_%s_%d_%d_plots.root",Nuclei.Data(),set,pionType),"OPEN");
    if(kVERBOSE) std::cout << "Opened file: " << fset->GetName() << std::endl;

    // Get varied cross section
    varied = (TH1D*)fset->Get(Form("h%s",iTypeString.Data()));
    varied->SetLineColor(set);
    //varied->SetLineWidth(4);
    varied->Smooth();
    if(varied->GetMaximum()*1.15 > kMaxScale) kMaxScale = varied->GetMaximum()*1.15;
    gr_varied = new TGraphErrors(varied);
    
    // Loop over momentum 
    for (int ipoint=0; ipoint<gr_varied->GetN(); ipoint++) { // Loop over momentum
      
      double x_new,y_new,x_min,y_min,x_max,y_max;
      
      gr_varied->GetPoint(ipoint,x_new,y_new);
      
      //Get min and max envelopes
      gr_min->GetPoint(ipoint,x_min,y_min);
      gr_max->GetPoint(ipoint,x_max,y_max);
      
      // If this throw is higher, use for max
      if(y_new > y_max){
	gr_max->SetPoint(ipoint,x_new,y_new);
      }
      
      // If this throw is lower, use for min
      if(y_new < y_min){
	gr_min->SetPoint(ipoint,x_new,y_new);
      }
      
    }

    fset->Close();
    
  }
  
  gr_nom->SetLineColor(kAzure+2);
  gr_min->SetLineColor(kAzure+2);
  gr_max->SetLineColor(kAzure+2);
  // Merge graphs to form uncertainty band
  gr_shade = FSIFitterUtils::MergeGraphsIntoEnvelope(gr_min,gr_max);
  gr_shade->SetFillColor(kAzure+2);
  //gr_shade->SetFillStyle(3344);
  
}
