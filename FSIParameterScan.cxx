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
      if(nuclFit) std::cout<<"jentry: " <<jentry << " mom: "<<mom<<" tot: "<<totfact<<" elas: "<<elafact<<" spi: "<<spifact<<" dpi: "<<dpifact<<std::endl;
      //else std::cout<<"jentry: " <<jentry << " mom: "<<mom<<" qe: "<<qe<<" abs: "<<abs<<" cx: "<<cx<<std::endl;
    }
}
