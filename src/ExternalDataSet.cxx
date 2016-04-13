//***********************************************************************************
// ExternalDataSet.cxx
// Class to store and manipulate external data
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "ExternalDataSet.hxx"

//***********************************************************************************
// Definition of a ExternalDataSet object
// Defined from the name of the csv data file
//***********************************************************************************
ExternalDataSet::ExternalDataSet(std::string fFileName) : FileName(fFileName)
{
  // This will load the data points into a TGraphErrors, TVectors and a TMatrix
  Initialize();
}

ExternalDataSet& ExternalDataSet::operator=(const ExternalDataSet &aDataSet ){}

//***********************************************************************************
// Parse info about the data set from the filename
//***********************************************************************************
void ExternalDataSet::ParseDataSetName(){
  
  std::stringstream ss(FileName);
  std::vector<std::string> parser;
  
  while( ss.good() ){
    std::string substr;
    getline( ss, substr, '_' );
    parser.push_back(substr);
  }
  
  Nuclei = parser[0];
  intrType = IntrNameToInt(parser[1]);
  pionType = parser[2];
  SetName = parser[3];

  std::cout << "\n%%%%% Building ExternalDataSet: "
	    << "\nSetName: "   << SetName
	    << "\nNuclei: "    << Nuclei
	    << "\nIntr Type: " << intrType
	    << "\nPion Type: " << pionType
	    <<std::endl;
  
}

//***********************************************************************************
// Convert interaction type to integer for easiness later on
//***********************************************************************************
int ExternalDataSet::IntrNameToInt(std::string fIntrType){
  
  int fIntrTypeNumber = -1;
  
  if(fIntrType == "reac")   fIntrTypeNumber = 1;
  if(fIntrType == "inel")   fIntrTypeNumber = 2;
  if(fIntrType == "abs")    fIntrTypeNumber = 3;
  if(fIntrType == "cx")     fIntrTypeNumber = 4;
  if(fIntrType == "dcx")    fIntrTypeNumber = 5;
  if(fIntrType == "abscx")  fIntrTypeNumber = 6;
  return fIntrTypeNumber;

}

//***********************************************************************************
// Just return the data set name
//***********************************************************************************
std::string ExternalDataSet::GetSetName(){
  
  return SetName;

}

//***********************************************************************************
// Load the data from the csv files into a TGraphErrors and TVectors
//***********************************************************************************
void ExternalDataSet::LoadData(){

  // Load data from csv file
  std::string FileNameWithDirectory = "data/" + FileName + ".csv";
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
void ExternalDataSet::BuildDataDiagonalMatrices(){
  
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
int ExternalDataSet::GetNumberOfDataPoints(){
  
  return nDataPoints;

}   

//***********************************************************************************
// Calculate the chisquare for this data set for a particular set of parameters
//***********************************************************************************
//double ExternalDataSet::CalculateChiSquare(Double_t this_qe, Double_t this_abs, Double_t this_cx, FSIParameterScan &aFSIParameterScan){
double ExternalDataSet::CalculateChiSquare(Double_t this_qe, Double_t this_abs, Double_t this_cx){

  bool VERBOSE_LEVEL = false;//true;

  if(VERBOSE_LEVEL){
    std::cout<< "Calculating chi-square for dataset: " << FileName << std::endl;
    PrintParameterSet(this_qe,this_abs,this_cx);
  }
  
  // TVectorD for the MC
  TVector vec_MC_Xsec(nDataPoints);

  // Loop through memory resident TTree
  const Int_t nentries = MiniMCScan->GetEntries();

  double qe;
  double abs;
  double cx;
  double mom;
  double xsec_mc;
  
  MiniMCScan->SetBranchAddress("qe",&qe);
  MiniMCScan->SetBranchAddress("abs",&abs);
  MiniMCScan->SetBranchAddress("cx",&cx);
  MiniMCScan->SetBranchAddress("mom",&mom);
  MiniMCScan->SetBranchAddress("specific_xsec",&xsec_mc);

  // Boolean to flag if parameter set was found
  bool foundParSet = false;

  // Loop through the tree
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    MiniMCScan->GetEntry(jentry);
    
    // If a different set of parameters, move on
    if(!(TMath::Abs(qe - this_qe) <0.01    && 
	 TMath::Abs(abs - this_abs) < 0.01 && 
	 TMath::Abs(cx - this_cx) < 0.01 )) continue;
    
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
    std::cout << "ExternalDataSet::CalculateChiSquare(): Warning did not find this parameter set!" 
	      << "\n ... Leaving now. " << std::endl;
    PrintParameterSet(this_qe,this_abs,this_cx);
    exit(-1);
  }

  // Check that the Data and MC Vectors match in size
  const int nbins = vec_data_Xsec.GetNoElements();
  if (nbins != vec_MC_Xsec.GetNoElements()) {
    std::cout << "Error ExternalDataSet::CalculateChiSquare: Number of data bins (" <<
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
	  newDataPoints--;
	  std::cout << "... there wasn't xsec info for mom: " << vec_data_Mom[ix] 
		    << ".Not touching the chisquare and reducing the number of data points" << std::endl;
	}
      }
      
      if(VERBOSE_LEVEL && ix==iy) 
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
double ExternalDataSet::GetChiSquare(){
  
  if(fDataSetChiSquare == -9999){
    std::cout << "ExternalDataSet::GetChiSquare ERROR: Haven't computed chi-square yet!" << std::endl;
    exit(-1);
  }
  return fDataSetChiSquare;

}

//***********************************************************************************
// Print the parameter set. Should make something more general
//***********************************************************************************
void ExternalDataSet::PrintParameterSet(double f_qe, double f_abs, double f_cx){

  std::cout << "f_qe: " << f_qe
	    << "\tf_abs: " << f_abs
	    << "\tf_cx: " << f_cx
	    << std::endl;
}

void ExternalDataSet::FillMiniTreeFromMCScan(){

  std::string filename = "input/scan_result_ieki_LE.root";
  FSIParameterScan aFSIParameterScan(filename);
  const Int_t nentries = aFSIParameterScan.fChain->GetEntries();
  std::cout << "MC Scan has # entries = " << nentries << std::endl;
  
  MiniMCScan = aFSIParameterScan.fChain->CloneTree(0);
  MiniMCScan->SetDirectory(0);

  double specific_xsec;
  MiniMCScan->Branch("specific_xsec",&specific_xsec,"specific_xsec/D");
  
  // Loop through the tree
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = aFSIParameterScan.LoadTree(jentry);
    if (ientry < 0) break;
    
    aFSIParameterScan.GetEntry(jentry);
    
    // Loop through data set points
    for (int point = 0; point < (int)vec_data_Mom.GetNoElements(); point++){
      
      // Compare momentum value from MC scan with momentum from this data set
      // to determine if this entry should be saved
      // NOTE: A check for the target should be added later on!!
      if(TMath::Abs(aFSIParameterScan.mom - vec_data_Mom[point])<0.05){
	
	// If reactive data point
	if(intrType == 1) specific_xsec = (aFSIParameterScan.xqe + aFSIParameterScan.xabs + aFSIParameterScan.xcx +aFSIParameterScan.xdcx +aFSIParameterScan.xhadr);
	// If quasi-elastic data point
	if(intrType == 2) specific_xsec = aFSIParameterScan.xqe;
	// If ABS data point
	if(intrType == 3) specific_xsec = aFSIParameterScan.xabs;
	// If CX data point
	if(intrType == 4) specific_xsec = aFSIParameterScan.xcx;
	// If ABS+CX data point
	if(intrType == 6) specific_xsec = (aFSIParameterScan.xabs + aFSIParameterScan.xcx);
	
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
void ExternalDataSet::Initialize(){

  fDataSetChiSquare = -9999;
  ParseDataSetName();
  LoadData();
  BuildDataDiagonalMatrices();
  FillMiniTreeFromMCScan();
}
