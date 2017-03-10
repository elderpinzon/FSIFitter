#include "FSIChi2GridNorm.hxx"
#include "TStopwatch.h" // remove later
#include "TROOT.h"
#include "TCanvas.h" //remove later

//***********************************************************************************
// Definition of a grid object. A vector of the datasets to be used for the chisquare
// calculation is an argument
//***********************************************************************************
FSIChi2GridNorm::FSIChi2GridNorm(std::vector< ExternalDataSet::ExternalDataSet > AllDataSets){

  nPointsUsed = 0;

  std::cout << "\nStarting Grid Object. Will use these datasets:" << std::endl;
  for (Int_t dataset = 0; dataset < (int)AllDataSets.size(); dataset++){
    std::cout << "\t"<< dataset+1 << ". " << AllDataSets[dataset].GetFileName() << std::endl;
    nPointsUsed += AllDataSets[dataset].GetNumberOfDataPoints();
    if(!(AllDataSets[dataset].GetSetName().Contains("DUET")))
      fNonDUETDataSets.push_back(AllDataSets[dataset]);
  }
  
  std::cout << "... " << AllDataSets.size() << " data sets with "
	    << nPointsUsed << " points used in total." << std::endl;

  // Copy vector of data set to internal data set vector
  fAllDataSets = AllDataSets;

  // Store chisquare finite grid in output file
  foutTree = new TFile("output/finite_chi2_grid.root","RECREATE");

  // Set-up TTree for storage of the finite grid
  ChiSquareFiniteGrid = new TTree("Chi2Grid","Finite Chi Square Grid Tree");

}

void FSIChi2GridNorm::LoadInterpolatedCrossSections(){

  // Octave engine initialization
  string_vector argv (2);
  argv(0) = "FSIOctaveInterpolation";
  argv(1) = "-q";
  octave_main (2, argv.c_str_vec (), 1);

  for(Int_t iNuclei = 0; iNuclei < (int)NucleiA.size(); iNuclei++){
    for(Int_t ipid = 0; ipid < (Int_t)pids.size(); ipid++){
      aInterpolatedCrossSection[iNuclei][ipid] = new InterpolatedCrossSectionsOctave(ScanFileName,
										     NucleiA[iNuclei],
										     pids[ipid]);
      aInterpolatedCrossSection[iNuclei][ipid]->BuildGrid();
    }
  }
  
  std::cout << "Done loading cross sections to octave objects" << std::endl;
  
}

Double_t FSIChi2GridNorm::GrabInterpolatedCrossSectionOctave(Int_t thisNuclei, Int_t thisPid, Int_t thisType, std::vector<double> fsipars){
  
  for(Int_t iNuclei = 0; iNuclei < (Int_t)NucleiA.size(); iNuclei++){
    for(Int_t ipid = 0; ipid < (Int_t)pids.size(); ipid++){
      
      if(thisNuclei == NucleiA[iNuclei] && thisPid == pids[ipid])
	return aInterpolatedCrossSection[iNuclei][ipid]->GetSplinedGridPoint(thisType,fsipars);
      
    }
  }
  
  std::cout << "No interpolated cross section available for thisNuclei: " << thisNuclei
	    << " thisPid: " << thisPid 
	    << std::endl;
  
  return 0;
}

Double_t FSIChi2GridNorm::GetChi2(std::vector<Double_t> fsipars, std::vector<Double_t> normpars){

  // Global chi2
  Double_t chi2 = 0;

  // Loop over all NonDUET datasets, but skip DUET. Will do it later
  for (Int_t dataset = 0; dataset < (int)fNonDUETDataSets.size(); dataset++){

    if(fNonDUETDataSets[dataset].GetSetName() == "DUETSeparate")
      continue;
    
    for(Int_t point = 0; point < fNonDUETDataSets[dataset].GetXsecErrorVector().GetNoElements(); point++){

      // Add momentum to vector of parameters to be passed
      fsipars.insert(fsipars.begin(),fNonDUETDataSets[dataset].GetMomVector()[point]);
      
      Int_t thisType = fNonDUETDataSets[dataset].GetIntrType();
      Double_t mc_xsec = GrabInterpolatedCrossSectionOctave(fNonDUETDataSets[dataset].GetNucleus(),
							    fNonDUETDataSets[dataset].GetPionType(),
							    thisType,fsipars);

      // Scale MC prediction by normalization factor
      mc_xsec *= normpars[dataset];
      Double_t thisChi2 = (fNonDUETDataSets[dataset].GetXsecVector()[point] - mc_xsec) / fNonDUETDataSets[dataset].GetXsecErrorVector()[point];
      
      chi2 += thisChi2*thisChi2 / fNonDUETDataSets[dataset].GetNumberOfDataPoints();

      // Add pull term
      Double_t pullTerm = (normpars[dataset] - 1)/0.4;

      chi2 += pullTerm * pullTerm;

      std::cout << " dataset: " << fNonDUETDataSets[dataset].GetSetName()
		<< " npoints: " << fNonDUETDataSets[dataset].GetNumberOfDataPoints()
		<< " type: " << thisType
		<< " mom: " << fNonDUETDataSets[dataset].GetMomVector()[point]
		<< " data: " << fNonDUETDataSets[dataset].GetXsecVector()[point]
		<< "±" << fNonDUETDataSets[dataset].GetXsecErrorVector()[point]
		<< " norm: " << normpars[dataset]
		<< " mc: " << mc_xsec
		<< " thisChi2: " << thisChi2
		<< " chi2: " << chi2
		<< std::endl;

      // Remove the momentum, just in case
      fsipars.erase(fsipars.begin());
    }
  }

  // Loop over all datasets, but only work on DUET
  for (Int_t dataset = 0; dataset < (int)fAllDataSets.size(); dataset++){

    if(fAllDataSets[dataset].GetSetName() != "DUETSeparate")
      continue;

    // Will only use this for DUET
    TMatrixD vec_diff_Xsec(10,1);
    TMatrixD vec_diff_Xsec_transpose(1,10);
    
    
    for(Int_t point = 0; point < fAllDataSets[dataset].GetXsecErrorVector().GetNoElements(); point++){

      // Add momentum to vector of parameters to be passed
      fsipars.insert(fsipars.begin(),fAllDataSets[dataset].GetMomVector()[point]);
      
      Int_t thisType = fAllDataSets[dataset].GetIntrType();

      // Ugly hack. First five DUET are CX
      if(fAllDataSets[dataset].GetSetName() == "DUETSeparate" && point <5)
	thisType = 4;
      // Ugly hack. Last five DUET are ABS
      else if(fAllDataSets[dataset].GetSetName() == "DUETSeparate" && point >= 5)
	thisType = 3;
      
      Double_t mc_xsec = GrabInterpolatedCrossSectionOctave(fAllDataSets[dataset].GetNucleus(),
							    fAllDataSets[dataset].GetPionType(),
							    thisType,fsipars);

      vec_diff_Xsec[point][0] = fAllDataSets[dataset].GetXsecVector()[point] - mc_xsec;

      std::cout << " dataset: " << fAllDataSets[dataset].GetSetName()
		<< " npoints: " << fAllDataSets[dataset].GetNumberOfDataPoints()
		<< " type: " << thisType
		<< " mom: " << fAllDataSets[dataset].GetMomVector()[point]
		<< " data: " << fAllDataSets[dataset].GetXsecVector()[point]
		<< "±" << fAllDataSets[dataset].GetXsecErrorVector()[point]
		<< " mc: " << mc_xsec
		<< std::endl;

      // Remove the momentum, just in case
      fsipars.erase(fsipars.begin());
      
    }

    // Now adding DUET
    vec_diff_Xsec.Print();
    fAllDataSets[dataset].GetCovarianceInv().Print();
    vec_diff_Xsec_transpose.Transpose(vec_diff_Xsec);
    TMatrixD res(1,1);
    res = vec_diff_Xsec_transpose * fAllDataSets[dataset].GetCovarianceInv() * vec_diff_Xsec;
    Double_t duetChiSquare = res[0][0];
    std::cout << "\n_duetChiSquare: " << duetChiSquare << std::endl;
    chi2 += duetChiSquare;
  }

  std::cout << "Computed chi2: " << chi2 
	    << " for " << (int)fAllDataSets.size() << " datasets"
	    << std::endl;
  return chi2;

}
