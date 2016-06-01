#ifndef FSIParameterScan_hxx
#define FSIParameterScan_hxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <cstdlib>

class FSIParameterScan {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   bool nuclFit;

   // Declaration of leaf types
   Double_t        mom;
   Double_t        qe;
   Double_t        abs;
   Double_t        cx;
   Double_t        nall;
   Double_t        nqe;
   Double_t        nabs;
   Double_t        ncx;
   Double_t        ndcx;
   Double_t        nhadr;
   Double_t        xqe;
   Double_t        xabs;
   Double_t        xcx;
   Double_t        xdcx;
   Double_t        xhadr;
   //nucleon
   Double_t        totfact;
   Double_t        elafact;
   Double_t        spifact;
   Double_t        dpifact;
   Double_t        xsecelastic;
   Double_t        xsecspi;
   Double_t        xsecdpi;
   Double_t        xsecrxn;
   Double_t        xsectotal;   

   // List of branches
   TBranch        *b_mom;   //!
   TBranch        *b_qe;   //!
   TBranch        *b_abs;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_nall;   //!
   TBranch        *b_nqe;   //!
   TBranch        *b_nabs;   //!
   TBranch        *b_ncx;   //!
   TBranch        *b_ndcx;   //!
   TBranch        *b_nhadr;   //!
   TBranch        *b_xqe;   //!
   TBranch        *b_xabs;   //!
   TBranch        *b_xcx;   //!
   TBranch        *b_xdcx;   //!
   TBranch        *b_xhadr;   //!
   //nucleon
   TBranch        *b_totfact;   //!
   TBranch        *b_elafact;   //!
   TBranch        *b_spifact;   //!
   TBranch        *b_dpifact;   //!
   TBranch        *b_xsecelastic;   //!
   TBranch        *b_xsecspi;   //!
   TBranch        *b_xsecdpi;   //!
   TBranch        *b_xsecrxn;   //!
   TBranch        *b_xsectotal;   //!

   Double_t qe_low;
   Double_t qe_high;
   Double_t qe_step;
   Double_t abs_low;
   Double_t abs_high;
   Double_t abs_step;
   Double_t cx_low;
   Double_t cx_high;
   Double_t cx_step;

   Double_t fsi_par_low[3];
   Double_t fsi_par_high[3];
   Double_t fsi_par_step[3];

   FSIParameterScan(){};
   FSIParameterScan(std::string ScanFileName, bool fnuclFit=0);
   virtual ~FSIParameterScan();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual void     Test();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void DetermineGridSize();
   void InterpolateXsecs();
};

#endif

#ifdef FSIParameterScan_cxx
FSIParameterScan::FSIParameterScan(std::string ScanFileName, bool fnuclFit) : fChain(0) 
{
    nuclFit=fnuclFit;

  TTree *tree = 0;
  Long_t id,size,flags,modtime;
  if(!gSystem->GetPathInfo(ScanFileName.c_str(),&id,&size,&flags,&modtime)){
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(ScanFileName.c_str());
    if (!f || !f->IsOpen()) {
      f = new TFile(ScanFileName.c_str());
    }
    if (nuclFit) f->GetObject("T",tree);
    else f->GetObject("t",tree);
  }
  else{
    std::cout << "\nFSIParameterScan could not find specified parameter scan file "
	      << ScanFileName << " ... exiting now" << std::endl;
    exit(-1);
  }
    
  //std::cout << "\nReading specified parameter scan file: " << ScanFileName 
  //	    << ". File has " << tree->GetEntries() << " entries." << std::endl;
  
  
  Init(tree);
}

FSIParameterScan::~FSIParameterScan()
{
  // Not sure why this is called right before FSIChi2Grid::InterpolateFiniteGrid() :(
  //std::cout<<"~FSIParameterScan"<<std::endl;
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FSIParameterScan::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) {std::cout<<"can't get entry"<<std::endl;return 0;}
   return fChain->GetEntry(entry);
}
Long64_t FSIParameterScan::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FSIParameterScan::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   if (nuclFit){
     fChain->SetBranchAddress("mom", &mom, &b_mom);
     fChain->SetBranchAddress("totfact", &totfact, &b_totfact);
     fChain->SetBranchAddress("elafact", &elafact, &b_elafact);
     fChain->SetBranchAddress("spifact", &spifact, &b_spifact);
     fChain->SetBranchAddress("dpifact", &dpifact, &b_dpifact);
     fChain->SetBranchAddress("xsecelastic", &xsecelastic, &b_xsecelastic);
     fChain->SetBranchAddress("xsecspi", &xsecspi, &b_xsecspi);
     fChain->SetBranchAddress("xsecdpi", &xsecdpi, &b_xsecdpi);
     fChain->SetBranchAddress("xsecrxn", &xsecrxn, &b_xsecrxn);
     fChain->SetBranchAddress("xsectotal", &xsectotal, &b_xsectotal);
   }else{
     fChain->SetBranchAddress("mom", &mom, &b_mom);
     fChain->SetBranchAddress("qe", &qe, &b_qe);
     fChain->SetBranchAddress("abs", &abs, &b_abs);
     fChain->SetBranchAddress("cx", &cx, &b_cx);
     fChain->SetBranchAddress("nall", &nall, &b_nall);
     fChain->SetBranchAddress("nqe", &nqe, &b_nqe);
     fChain->SetBranchAddress("nabs", &nabs, &b_nabs);
     fChain->SetBranchAddress("ncx", &ncx, &b_ncx);
     fChain->SetBranchAddress("ndcx", &ndcx, &b_ndcx);
     fChain->SetBranchAddress("nhadr", &nhadr, &b_nhadr);
     fChain->SetBranchAddress("xqe", &xqe, &b_xqe);
     fChain->SetBranchAddress("xabs", &xabs, &b_xabs);
     fChain->SetBranchAddress("xcx", &xcx, &b_xcx);
     fChain->SetBranchAddress("xdcx", &xdcx, &b_xdcx);
     fChain->SetBranchAddress("xhadr", &xhadr, &b_xhadr);
   }
   Notify();
}

Bool_t FSIParameterScan::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FSIParameterScan::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FSIParameterScan::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FSIParameterScan_cxx
