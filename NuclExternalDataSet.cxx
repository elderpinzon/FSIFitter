//***********************************************************************************
// NuclExternalDataSet.cxx
// Class to store and manipulate external data
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "NuclExternalDataSet.hxx"

//***********************************************************************************
// Definition of a NuclExternalDataSet object
// Defined from the name of the csv data file
//***********************************************************************************
NuclExternalDataSet::NuclExternalDataSet(std::string fFileName) : FileName(fFileName)
{
  // This will load the data points into a TGraphErrors, TVectors and a TMatrix
  Initialize();
}

NuclExternalDataSet& NuclExternalDataSet::operator=(const NuclExternalDataSet &aDataSet ){}

//***********************************************************************************
// Parse info about the data set from the filename
//***********************************************************************************
void NuclExternalDataSet::ParseDataSetName(){
  
  std::stringstream ss(FileName);
  std::vector<std::string> parser;
  
  while( ss.good() ){
    std::string substr;
    getline( ss, substr, '_' );
    parser.push_back(substr);
  }
  
  Nuclei = parser[0];
  intrType = IntrNameToInt(parser[1]);
  nucleonType = parser[2];
  SetName = parser[3];

  std::cout << "\n%%%%% Building NuclExternalDataSet: "
      << "\nSetName: "   << SetName
      << "\nNuclei: "    << Nuclei
      << "\nIntr Type: " << intrType
      << "\nNucleon Type: " << nucleonType
      <<std::endl;
  
}

//***********************************************************************************
// Convert interaction type to integer for easiness later on
//***********************************************************************************
int NuclExternalDataSet::IntrNameToInt(std::string fIntrType){
  
  int fIntrTypeNumber = -1;
  
  if(fIntrType == "reac")   fIntrTypeNumber = 1;
  if(fIntrType == "elas")   fIntrTypeNumber = 2;
  if(fIntrType == "tot")    fIntrTypeNumber = 3;
  // if(fIntrType == "cx")     fIntrTypeNumber = 4;
  // if(fIntrType == "dcx")    fIntrTypeNumber = 5;
  // if(fIntrType == "abscx")  fIntrTypeNumber = 6;
  return fIntrTypeNumber;

}

//***********************************************************************************
// Just return the data set name
//***********************************************************************************
std::string NuclExternalDataSet::GetSetName(){
  
  return SetName;

}

//***********************************************************************************
// Load the data from the csv files into a TGraphErrors and TVectors
//***********************************************************************************
void NuclExternalDataSet::LoadData(){

  // Load data from csv file
  std::string FileNameWithDirectory = "nucleon_data/" + FileName + ".csv";
  fDataPointsGraph = new TGraphErrors(FileNameWithDirectory.c_str(),"%lf %lg %lg");
  
  // Set number of data points in this experiment
  nDataPoints = fDataPointsGraph->GetN();
  
  // Resize vector to required length
  vec_data_Mom.ResizeTo(TVector(nDataPoints));
  vec_data_Xsec.ResizeTo(TVector(nDataPoints));
  vec_data_Xsec_error.ResizeTo(TVector(nDataPoints));

  //Copy points to global variables
  for(Int_t i=0; i<nDataPoints; i++) {

    Double_t fX, fY, fError;
    fDataPointsGraph->GetPoint(i,fX,fY);
    fError = fDataPointsGraph->GetErrorY(i);
    vec_data_Mom[i] = fX;
    vec_data_Xsec[i] = fY;
    vec_data_Xsec_error[i] = fError;
    std::cout << i+1 << " " 
	      << vec_data_Mom[i] << " " 
	      << vec_data_Xsec[i] << " " 
	      << vec_data_Xsec_error[i] 
	      << std::endl;
    
  }

  std::cout << "Filled data vectors with # entries: " << vec_data_Mom.GetNoElements() << std::endl;

}

//***********************************************************************************
// Build diagonal covariance matrix for this data set
//***********************************************************************************
void NuclExternalDataSet::BuildDataDiagonalMatrices(){
  
  // Resize matrices
  m_data_cov.ResizeTo(nDataPoints,nDataPoints);
  m_data_cov_inv.ResizeTo(nDataPoints,nDataPoints);
  m_data_cor.ResizeTo(nDataPoints,nDataPoints);

  // Build matrices
  for (int ix=0; ix<nDataPoints; ix++) {
    for (int iy=0; iy<nDataPoints; iy++) {
      
      // Make diagonals
      if (ix==iy) {
	// 1 for correlation matrix
	       m_data_cor(ix,iy) = 1;

	// Error^2 for covariance matrix
	       m_data_cov(ix,iy) = vec_data_Xsec_error(ix)*vec_data_Xsec_error(iy);
      }

      else {
	       m_data_cor(ix,iy) = 0;
	       m_data_cov(ix,iy) = 0;
      }
      
      m_data_cov_inv(ix,iy) = m_data_cov(ix,iy);
      
    }
  }

  // Check that it's decomposable
  TDecompChol decomp(m_data_cov);
  if (!decomp.Decompose())
    std::cout << "Decomposition failed";

  // Check that eigenvalues > 0 
  TVectorD eigenvals(m_data_cov.GetNrows());
  m_data_cov.EigenVectors(eigenvals);
  for (int eval=0; eval<eigenvals.GetNoElements(); eval++) {
    if (eigenvals(eval)<=0)
      std::cout << "Error: Eigenvalue #" << eval << " = " << eigenvals(eval) << " < 0" << std::endl;
  }

  m_data_cov_inv.Invert();

  std::cout << "Built covariance matrices for this dataset" << std::endl;

}

//***********************************************************************************
// Just return the number of data points
//***********************************************************************************
int NuclExternalDataSet::GetNumberOfDataPoints(){
  
  return nDataPoints;

}   

//***********************************************************************************
// Calculate the chisquare for this data set for a particular set of parameters
//***********************************************************************************
//double NuclExternalDataSet::CalculateChiSquare(Double_t this_qe, Double_t this_abs, Double_t this_cx, NuclFSIParameterScan &aNuclFSIParameterScan){
double NuclExternalDataSet::CalculateChiSquare(Double_t this_reac, Double_t this_elas, Double_t this_tot){

  bool VERBOSE_LEVEL = false;//true;

  if(VERBOSE_LEVEL){
    std::cout<< "Calculating chi-square for dataset: " << FileName << std::endl;
    PrintParameterSet(this_reac,this_elas,this_tot);
  }
  
  // TVectorD for the MC
  TVector vec_MC_Xsec(nDataPoints);

  // Loop through memory resident TTree
  const Int_t nentries = MiniMCScan->GetEntries();

  double spifact;
  // double dpifact;
  double elafact;
  double totfact;
  double mom;
  double xsec_mc;
  
  MiniMCScan->SetBranchAddress("spifact",&spifact);
  MiniMCScan->SetBranchAddress("elafact",&elafact);
  MiniMCScan->SetBranchAddress("totfact",&totfact);
  MiniMCScan->SetBranchAddress("mom",&mom);
  MiniMCScan->SetBranchAddress("specific_xsec",&xsec_mc);

  // Boolean to flag if parameter set was found
  bool foundParSet = false;

  // Loop through the tree
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    MiniMCScan->GetEntry(jentry);
    
    // If a different set of parameters, move on
    if(!(TMath::Abs(spifact - this_reac) <0.01    && 
	 TMath::Abs(elafact - this_elas) < 0.01 && 
	 TMath::Abs(totfact - this_tot) < 0.01 )) continue;
    
    // Loop through data set points
    for (int point = 0; point < (int)vec_data_Mom.GetNoElements(); point++){
      
      // Compare momentum value from MC scan with momentum from this data set
      if(TMath::Abs(mom - vec_data_Mom[point])<0.05){
	
	vec_MC_Xsec[point] = xsec_mc;
	
	if(VERBOSE_LEVEL)
	  std::cout << "mom: " << mom << " found vec_MC_Xsec: " << vec_MC_Xsec[point] << std::endl;
	
	// Flag parameter set as found
	if(!foundParSet) foundParSet = true;
      }
    }
  }
  
  // Exit if a parameter set was not found. This should not happen
  if(!foundParSet){
    std::cout << "NuclExternalDataSet::CalculateChiSquare(): Warning did not find this parameter set!" 
	      << "\n ... Leaving now. " << std::endl;
    PrintParameterSet(this_reac,this_elas,this_tot);
    exit(-1);
  }

  // Check that the Data and MC Vectors match in size
  const int nbins = vec_data_Xsec.GetNoElements();
  if (nbins != vec_MC_Xsec.GetNoElements()) {
    std::cout << "Error NuclExternalDataSet::CalculateChiSquare: Number of data bins (" <<
      nbins << ") != mc bins (" << vec_MC_Xsec.GetNoElements() << ")" << std::endl;
  std:exit (1);
  }

  // Need a new nDataPoints in case some MC point is missing
  int newDataPoints = nDataPoints;

  // Actually do the chi-square math
  double sum_chi2 = 0;

  for (int ix=0; ix<nbins; ix++) {
    for (int iy=0; iy<nbins; iy++) {
      
      // Check for missing MC information. For now, skip this data point
      if(vec_MC_Xsec(ix) != 0 || vec_MC_Xsec(iy) != 0)
	sum_chi2 += (vec_data_Xsec(ix)-vec_MC_Xsec(ix)) * m_data_cov_inv(ix,iy) * (vec_data_Xsec(iy)-vec_MC_Xsec(iy));
      else{
	if(ix==iy){
	  // newDataPoints--;
	  std::cout << "... there wasn't xsec info for mom: " << vec_data_Mom[ix] 
		    << ".Not touching the chisquare and reducing the number of data points" << std::endl;
	}
      }
      
      if(VERBOSE_LEVEL && ix==iy) 
      // if(ix==iy) 
	std::cout << "mom: " << vec_data_Mom[ix]
		  << " vec_data_Xsec: "<<vec_data_Xsec(ix) 
		  << " vec_MC_Xsec: "<<vec_MC_Xsec(ix) 
		  << " vec_diff_Xsec: "<< vec_data_Xsec(ix) - vec_MC_Xsec(ix) 
		  << " m_data_cov_inv(ix): " << m_data_cov_inv(ix,ix) 
		  << " sum_chi2: " << sum_chi2
		  << std::endl;
    }
  }
  
  // Divide by the number of data points
  sum_chi2 /= newDataPoints;
  
  if(VERBOSE_LEVEL) std::cout << "chi_square: " << sum_chi2 << std::endl;
  
  fDataSetChiSquare = sum_chi2;
  
  return sum_chi2;
  
}

//***********************************************************************************
// Just return the calculate chi square
//***********************************************************************************
double NuclExternalDataSet::GetChiSquare(){
  
  if(fDataSetChiSquare == -9999){
    std::cout << "NuclExternalDataSet::GetChiSquare ERROR: Haven't computed chi-square yet!" << std::endl;
    exit(-1);
  }
  return fDataSetChiSquare;

}

//***********************************************************************************
// Print the parameter set. Should make something more general
//***********************************************************************************
void NuclExternalDataSet::PrintParameterSet(double f_reac, double f_elas, double f_tot){

  std::cout << "f_reac: " << f_reac
      << "\tf_elas: " << f_elas
      << "\tf_tot: " << f_tot
      << std::endl;
}
//***********************************************************************************
void NuclExternalDataSet::FillMiniTreeFromMCScan(){

  std::string filename = "input/fullscan0418_selec.root";
  NuclFSIParameterScan aNuclFSIParameterScan(filename);
  const Int_t nentries = aNuclFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;
  
  MiniMCScan = aNuclFSIParameterScan.fChain->CloneTree(0);
  MiniMCScan->SetDirectory(0);

  double specific_xsec;
  MiniMCScan->Branch("specific_xsec",&specific_xsec,"specific_xsec/D");
  
  // Loop through the tree
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aNuclFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;
    
    aNuclFSIParameterScan.GetEntry(jentry);
    
    // Loop through data set points
    for (int point = 0; point < (int)vec_data_Mom.GetNoElements(); point++){
      
      // Compare momentum value from MC scan with momentum from this data set
      // to determine if this entry should be saved
      // NOTE: A check for the target should be added later on!!
      if(TMath::Abs(aNuclFSIParameterScan.mom - vec_data_Mom[point])<0.1){
	
  // If reactive data point
  if(intrType == 1) specific_xsec = aNuclFSIParameterScan.xsecrxn;
  // If quasi-elastic data point
  if(intrType == 2) specific_xsec = aNuclFSIParameterScan.xsecelastic;
  // If total data point
  if(intrType == 3) specific_xsec = aNuclFSIParameterScan.xsectotal;
  // If CX data point
  // if(intrType == 4) vec_MC_Xsec[point] = aNucleonNEUTParameterScan.xcx;
  // If ABS+CX data point
  // if(intrType == 6) vec_MC_Xsec[point] = (aNucleonNEUTParameterScan.xabs + aNucleonNEUTParameterScan.xcx);
  
	MiniMCScan->Fill();
	
      }
    }
  }
  
  std::cout << "MiniMCScan has # entries: " << MiniMCScan->GetEntries()
	    << ". Entries per data point: " << MiniMCScan->GetEntries()/nDataPoints
	    << std::endl;

}
      


  
  

//***********************************************************************************
// This is called by the constructor 
//***********************************************************************************
void NuclExternalDataSet::Initialize(){

  fDataSetChiSquare = -9999;
  ParseDataSetName();
  LoadData();
  BuildDataDiagonalMatrices();
  FillMiniTreeFromMCScan();
}
