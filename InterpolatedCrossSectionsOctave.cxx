//***********************************************************************************
// InterpolatedCrossSectionsOctave.cxx
// Interpolates cross sections from an input grid scan using Octave
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "InterpolatedCrossSectionsOctave.hxx"

InterpolatedCrossSectionsOctave::InterpolatedCrossSectionsOctave(TString fFileName,Int_t target, Int_t pid) : filename(fFileName), thisTarget(target), thisPID(pid)
{
  
  std::cout << "\nWill interpolate cross sections from file: "<< filename 
	    << "\ntarget: " << thisTarget
	    << "\npid: " << thisPID << std::endl;

  // Steps for FSI. Will set step for Momenta in BuildMomentumIndices();
  for(int i = 0; i<nFSIparsFitted; i++)
    nStepsOctave[i+1] = nSteps[i];
  
  //5 FSI + Momentum
  nPars = nFSIparsFitted + 1;
  
  BuildMomentumIndices();
  
  for (int idim=0; idim<nPars; idim++){
    
    std::cout << "Building octave array of parameters: " 
	      << idim << " dim: " << nStepsOctave[idim] 
	      << std::endl;
    
    // Vector with FSI pars
    x_vec[idim] = new Matrix(nStepsOctave[idim],1);
    for (octave_idx_type i = 0; i < nStepsOctave[idim]; i++){
      // First entry is momentum
      if(idim == 0){
	(*x_vec[idim])(i,0) = allMoms[i];
      }
      // Other entries are FSI parameters
      // Note: FSIParsMin and FSIParsStep are off by one index
      else{
	(*x_vec[idim])(i,0) = FSIParsMin[idim-1] + i*FSIParsStep[idim-1];
      }
    }
  }
  
  dim_vector dim_vec(nStepsOctave[0],nStepsOctave[1],nStepsOctave[2],nStepsOctave[3],nStepsOctave[4],nStepsOctave[5]);

  // Initialize with large value in case some failed in filling tree
  v_reac = new Array<double>(dim_vec,9999);
  v_qe = new Array<double>(dim_vec,9999);
  v_abs = new Array<double>(dim_vec,9999);
  v_cx = new Array<double>(dim_vec,9999);
  v_abscx = new Array<double>(dim_vec,9999);
    
}

InterpolatedCrossSectionsOctave::~InterpolatedCrossSectionsOctave(){
  
}

//***********************************************************************************
// Loop once over input file to build momentum indices to be used for the grid
//***********************************************************************************
void InterpolatedCrossSectionsOctave::BuildMomentumIndices(){
  
  allMoms.clear();
  
  // Load results from FSI parameter scan
  FSIParameterScan aFSIParameterScan(filename);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;
  
  if (aFSIParameterScan.fChain == 0) return;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;
    
    aFSIParameterScan.GetEntry(jentry);
    
    // Select target
    if(TMath::Abs(aFSIParameterScan.target - thisTarget) > 0.1) continue;
    
    // Select pid
    if(TMath::Abs(aFSIParameterScan.pid - thisPID) > 0.1) continue;
    
    // If haven't been included, add to vector
    if(!(std::find(allMoms.begin(), allMoms.end(), aFSIParameterScan.mom) != allMoms.end())) {
      allMoms.push_back(aFSIParameterScan.mom);
    }
    
  }
  
  // Sort vector by ascending order (for no particular reason)
  std::sort(allMoms.begin(), allMoms.end());
  std::cout << "allMoms size: " << allMoms.size() << ". Listed below: " << std::endl;
  for(Int_t i = 0; i < (Int_t)allMoms.size(); i++)
    printf("%.2f ",allMoms[i]);
  std::cout<<std::endl;

  // Set # steps for momentum to be used to initialize octave arrays
  // It is the first entry of the nstepsOctave array
  nStepsOctave[0] = (int)allMoms.size();
  
}

//***********************************************************************************
// Actually build the grid by looping over input grid scan file
//***********************************************************************************
void InterpolatedCrossSectionsOctave::BuildGrid(){

  // Load results from FSI parameter scan
  FSIParameterScan aFSIParameterScan(filename);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;

  if (aFSIParameterScan.fChain == 0) return;

  Int_t countVVEC = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;

    aFSIParameterScan.GetEntry(jentry);

    // Select target
    if(TMath::Abs(aFSIParameterScan.target - thisTarget) > 0.1) continue;

    // Select pid
    if(TMath::Abs(aFSIParameterScan.pid - thisPID) > 0.1) continue;

    //Hardcoded: Skip variations of FEFALL and FEFCXH
    if(TMath::Abs(aFSIParameterScan.FEFALL - 1.0) > 0.01) continue;
    if(TMath::Abs(aFSIParameterScan.FEFCXH - 1.8) > 0.01) continue;


    // Get index
    Int_t kIndices[8] = {0};
    kIndices[0] = MomToIndex(aFSIParameterScan.mom);
    // 0.1 added for round-off
    kIndices[1] = (aFSIParameterScan.FEFQE - FSIParsMin[0])/FSIParsStep[0] + 0.1;
    kIndices[2] = (aFSIParameterScan.FEFABS - FSIParsMin[1])/FSIParsStep[1] + 0.1;
    kIndices[3] = (aFSIParameterScan.FEFCX - FSIParsMin[2])/FSIParsStep[2] + 0.1;
    kIndices[4] = (aFSIParameterScan.FEFINEL - FSIParsMin[3])/FSIParsStep[3] + 0.1;
    kIndices[5] = (aFSIParameterScan.FEFQEH - FSIParsMin[4])/FSIParsStep[4] + 0.1;
    kIndices[6] = (aFSIParameterScan.FEFQEH - FSIParsMin[5])/FSIParsStep[5] + 0.1;
    kIndices[7] = (aFSIParameterScan.FEFQEH - FSIParsMin[6])/FSIParsStep[6] + 0.1;

    Array<octave_idx_type> index_array(dim_vector((octave_idx_type)nPars,1));
    
    // If LE point manually vary the HE parameters
    if(aFSIParameterScan.mom < 400){
      
      // Hardcoded: Skipping FEFCXH and FEFALL
      for(int vindex5 = 0; vindex5 < nStepsOctave[5]; vindex5++){
	for(int vindex4 = 0; vindex4 < nStepsOctave[4]; vindex4++){
	  
	  // It's one index off because mom is the first index
	  kIndices[5] = vindex5;
	  kIndices[4] = vindex4;
    	  
	  // Set index_array
	  for(Int_t i = 0; i<nPars; i++){
	    index_array(i) = kIndices[i];
	  }
	  countVVEC++;
	  (*v_reac)(index_array) = aFSIParameterScan.xreac;
	  (*v_qe)(index_array) = aFSIParameterScan.xqe;
	  (*v_abs)(index_array) = aFSIParameterScan.xabs;
	  (*v_cx)(index_array) = aFSIParameterScan.xcx;
	  (*v_abscx)(index_array) = aFSIParameterScan.xabs+aFSIParameterScan.xcx;
	  
	}
      }
    }
    
    // If HE point manually vary the LE parameters
    else if(aFSIParameterScan.mom > 500){
      
      //std::cout << "Varying LE parameters manually" << std::endl;
      
      for(int vindex3 = 0; vindex3 < nStepsOctave[3]; vindex3++){
    	for(int vindex2 = 0; vindex2 < nStepsOctave[2]; vindex2++){
    	  for(int vindex1 = 0; vindex1 < nStepsOctave[1]; vindex1++){
	    
    	    kIndices[1] = vindex1;
    	    kIndices[2] = vindex2;
    	    kIndices[3] = vindex3;
      
    	    // Set index_array
    	    for(Int_t i = 0; i<nPars; i++){
    	      index_array(i) = kIndices[i];
    	    }
    	    countVVEC++;
	    (*v_reac)(index_array) = aFSIParameterScan.xreac;
	    (*v_qe)(index_array) = aFSIParameterScan.xqe;
	    (*v_abs)(index_array) = aFSIParameterScan.xabs;
	    (*v_cx)(index_array) = aFSIParameterScan.xcx;
	    (*v_abscx)(index_array) = aFSIParameterScan.xabs+aFSIParameterScan.xcx;
	    
    	  }
    	}
      }
    }

    else{
      // Set index array
      for(Int_t i = 0; i<nPars; i++){
    	index_array(i) = kIndices[i];
      }
      countVVEC++;
      (*v_reac)(index_array) = aFSIParameterScan.xreac;
      (*v_qe)(index_array) = aFSIParameterScan.xqe;
      (*v_abs)(index_array) = aFSIParameterScan.xabs;
      (*v_cx)(index_array) = aFSIParameterScan.xcx;
      (*v_abscx)(index_array) = aFSIParameterScan.xabs+aFSIParameterScan.xcx;
      
    }

  }
    
  std::cout << "countVVEC: " << countVVEC << std::endl;
  //std::cout << "v_reac at the end: " << (*v_reac) << std::endl; 

}

//***********************************************************************************
// Get point using interpn spline function of GNU-Octave
//***********************************************************************************
Double_t InterpolatedCrossSectionsOctave::GetSplinedGridPoint(Int_t thisXSec, const std::vector<double> &x){

  // Build octave command. Start by passing the full matrix
  octave_value_list in;
  for (int idim=0; idim<nPars; idim++) {
    in.append( octave_value( *x_vec[idim] ) );
  }

  // Pass the v_vec filled with requested xsec
  if(thisXSec == 1)
    in.append( octave_value (*v_reac) );
  else if(thisXSec == 2)
    in.append( octave_value (*v_qe) );
  else if(thisXSec == 3)
    in.append( octave_value (*v_abs) );
  else if(thisXSec == 4)
    in.append( octave_value (*v_cx) );
  else if(thisXSec == 5)
    in.append( octave_value (*v_abscx) );
  else
    std::cout << "ERROR: InterpolatedCrossSectionsOctave being asked for wrong xsec" << std::endl;

  // Pass the point we want to interpolate to
  Matrix *x_point[nPars];
  const int nPointsToInterp = 1;
  for (int idim=0; idim<nPars; idim++) {
    x_point[idim] = new Matrix(nPointsToInterp,1);
    (*x_point[idim])(0,0) = x[idim];
    //std::cout << "x_point:[" << idim << "]: " <<  (*x_point[idim]) << std::endl;
  }

  for (int idim=0; idim<nPars; idim++)
    in.append( octave_value( *x_point[idim] ) );

  // Spline interpolation
  //in.append( octave_value("spline"));
  in.append( octave_value("linear"));

  // Default value outside of the grid
  in.append( octave_value(9999));

  // Evaluate command
  octave_value_list out_v_point = feval ("interpn",in,nPointsToInterp);
  
  // Get result
  Matrix v_point = Matrix(1,1);
  if (!error_state && out_v_point.length () > 0) {
    v_point = out_v_point(0).matrix_value ();
  }

  Double_t splined_chi2 = v_point(0);

  return splined_chi2;

}

//***********************************************************************************
// Get index position of momentum in allMoms array
//***********************************************************************************
Int_t InterpolatedCrossSectionsOctave::MomToIndex(Double_t thisMom){

  for(int i = 0; i<(int)allMoms.size(); i++){

    if(TMath::Abs(thisMom-allMoms[i] < 0.1))
      return i;
  }

  // If we are here, we are in trouble
  std::cout << "ERROR: InterpolatedCrossSectionsOctave::MomToIndex failed for thisMom: "
	    << thisMom << std::endl;

  return -999;

}
