#ifndef NuclFSIParameterScan_h
#define NuclFSIParameterScan_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <cstdlib>

class NuclFSIParameterScan {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        mom;
   Double_t        totfact;
   Double_t        elafact;
   Double_t        spifact;
   Double_t        dpifact;

   // Double_t        nall;
   // Double_t        nqe;
   // Double_t        nabs;
   // Double_t        ncx;
   // Double_t        ndcx;
   // Double_t        nhadr;
   Double_t        xsecelastic;
   Double_t        xsecspi;
   Double_t        xsecdpi;
   Double_t        xsecrxn;
   Double_t        xsectotal;

   // List of branches
   TBranch        *b_mom;   //!
   TBranch        *b_totfact;   //!
   TBranch        *b_elafact;   //!
   TBranch        *b_spifact;   //!
   TBranch        *b_dpifact;   //!
   // TBranch        *b_nall;   //!
   // TBranch        *b_nqe;   //!
   // TBranch        *b_nabs;   //!
   // TBranch        *b_ncx;   //!
   // TBranch        *b_ndcx;   //!
   // TBranch        *b_nhadr;   //!
   TBranch        *b_xsecelastic;   //!
   TBranch        *b_xsecspi;   //!
   TBranch        *b_xsecdpi;   //!
   TBranch        *b_xsecrxn;   //!
   TBranch        *b_xsectotal;   //!

   // Double_t qe_low;
   // Double_t qe_high;
   // Double_t qe_step;
   // Double_t abs_low;
   // Double_t abs_high;
   // Double_t abs_step;
   // Double_t cx_low;
   // Double_t cx_high;
   // Double_t cx_step;

   Double_t fsi_par_low[3];
   Double_t fsi_par_high[3];
   Double_t fsi_par_step[3];

   NuclFSIParameterScan(){};
   NuclFSIParameterScan(std::string ScanFileName);
   virtual ~NuclFSIParameterScan();
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

#ifdef NuclFSIParameterScan_cxx
NuclFSIParameterScan::NuclFSIParameterScan(std::string ScanFileName) : fChain(0) 
{
  TTree *tree = 0;
  Long_t id,size,flags,modtime;
  if(!gSystem->GetPathInfo(ScanFileName.c_str(),&id,&size,&flags,&modtime)){
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(ScanFileName.c_str());
    if (!f || !f->IsOpen()) {
      f = new TFile(ScanFileName.c_str());
    }
    f->GetObject("T",tree);
  }
  else{
    std::cout << "\nNuclFSIParameterScan could not find specified parameter scan file "
	      << ScanFileName << " ... exiting now" << std::endl;
    exit(-1);
  }
    
  //std::cout << "\nReading specified parameter scan file: " << ScanFileName 
  //	    << ". File has " << tree->GetEntries() << " entries." << std::endl;
  
  
  Init(tree);
}

NuclFSIParameterScan::~NuclFSIParameterScan()
{
  // Not sure why this is called right before FSIChi2Grid::InterpolateFiniteGrid() :(
  //std::cout<<"~NuclFSIParameterScan"<<std::endl;
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NuclFSIParameterScan::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) {std::cout<<"can't get entry"<<std::endl;return 0;}
   return fChain->GetEntry(entry);
}
Long64_t NuclFSIParameterScan::LoadTree(Long64_t entry)
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

void NuclFSIParameterScan::Init(TTree *tree)
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

   fChain->SetBranchAddress("mom", &mom, &b_mom);
   fChain->SetBranchAddress("totfact", &totfact, &b_totfact);
   fChain->SetBranchAddress("elafact", &elafact, &b_elafact);
   fChain->SetBranchAddress("spifact", &spifact, &b_spifact);
   fChain->SetBranchAddress("dpifact", &dpifact, &b_dpifact);
   // fChain->SetBranchAddress("nall", &nall, &b_nall);
   // fChain->SetBranchAddress("nqe", &nqe, &b_nqe);
   // fChain->SetBranchAddress("nabs", &nabs, &b_nabs);
   // fChain->SetBranchAddress("ncx", &ncx, &b_ncx);
   // fChain->SetBranchAddress("ndcx", &ndcx, &b_ndcx);
   // fChain->SetBranchAddress("nhadr", &nhadr, &b_nhadr);
   fChain->SetBranchAddress("xsecelastic", &xsecelastic, &b_xsecelastic);
   fChain->SetBranchAddress("xsecspi", &xsecspi, &b_xsecspi);
   fChain->SetBranchAddress("xsecdpi", &xsecdpi, &b_xsecdpi);
   fChain->SetBranchAddress("xsecrxn", &xsecrxn, &b_xsecrxn);
   fChain->SetBranchAddress("xsectotal", &xsectotal, &b_xsectotal);
   Notify();
}

Bool_t NuclFSIParameterScan::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NuclFSIParameterScan::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NuclFSIParameterScan::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NuclFSIParameterScan_cxx
