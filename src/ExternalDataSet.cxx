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
ExternalDataSet::ExternalDataSet(TString fFileName): FileName(fFileName)
{
  // This will load the data points into a TGraphErrors, TVectors and a TMatrix
  Initialize();
}

//***********************************************************************************
// Parse info about the data set from the filename
//***********************************************************************************
void ExternalDataSet::ParseDataSetName(){
  
  TObjArray *parser = FileName.Tokenize("_");
  
  NucleiString = ((TObjString *)(parser->At(0)))->String();
  Nuclei = FSIFitterUtils::NucleiNameToInt(NucleiString);
  intrTypeString = ((TObjString *)(parser->At(1)))->String();
  intrType = FSIFitterUtils::IntrNameToInt(intrTypeString);
  pionTypeString = ((TObjString *)(parser->At(2)))->String();
  pionType = FSIFitterUtils::PionTypeToInt(pionTypeString);
  SetName = ((TObjString *)(parser->At(3)))->String().ReplaceAll(".csv","");
  
  std::cout << "\n%%%%% Building Pion ExternalDataSet: "
	    << "\nSetName: "   << SetName
	    << "\nNuclei: "    << Nuclei
	    << "\nIntr Type: " << intrType
	    << "\nPion Type: " << pionType
	    <<std::endl;
}

//***********************************************************************************
// Load the data from the csv files into a TGraphErrors and TVectors
//***********************************************************************************
void ExternalDataSet::LoadData(){
  
  // Load data from csv file
  std::string FileNameWithDirectory = "";
  FileNameWithDirectory = "pion_data/" + FileName + ".csv";
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
  
  // Correlation factor
  Double_t correlationFactor = 0;
  //double correlationFactor = 0.4*0.4;

  std::cout << "Assuming correlationFactor: " << correlationFactor << std::endl;

  // Resize matrices
  m_data_cov.ResizeTo(nDataPoints,nDataPoints);
  m_data_cov_inv.ResizeTo(nDataPoints,nDataPoints);
  m_data_cor.ResizeTo(nDataPoints,nDataPoints);

  // Build matrices
  for (Int_t ix=0; ix<nDataPoints; ix++) {
    for (Int_t iy=0; iy<nDataPoints; iy++) {
      
      // Make diagonals
      if (ix==iy) {

	// Error^2 for covariance matrix
	m_data_cov(ix,iy) = vec_data_Xsec_error(ix)*vec_data_Xsec_error(iy) * nDataPoints;
	//m_data_cov(ix,iy) = vec_data_Xsec_error(ix)*vec_data_Xsec_error(iy);
      }
      
      else {
	//m_data_cov(ix,iy) = vec_data_Xsec_error(ix)*vec_data_Xsec_error(iy) * correlationFactor;
	m_data_cov(ix,iy) = vec_data_Xsec_error(ix)*vec_data_Xsec_error(iy) * nDataPoints * correlationFactor;
      }
      
      m_data_cov_inv(ix,iy) = m_data_cov(ix,iy);
      
    }
  }

  // Check that it's decomposable
  TDecompChol decomp(m_data_cov);
  if (!decomp.Decompose()){
    std::cout << "Decomposition failed for matrix:" << std::endl;
    m_data_cov.Print();
    std::cout << " try to invert it: " << std::endl;
    m_data_cov_inv.Invert();
    m_data_cov_inv.Print();
    std::cout << " check if we recover an unitary matrix: " << std::endl;
    TMatrixD check_unit(nDataPoints,nDataPoints);
    check_unit = m_data_cov * m_data_cov_inv;
    check_unit.Print();
    //std::exit(1);
  }
  
  
  // Check that eigenvalues > 0 
  TVectorD eigenvals(m_data_cov.GetNrows());
  m_data_cov.EigenVectors(eigenvals);
  for (Int_t eval=0; eval<eigenvals.GetNoElements(); eval++) {
    if (eigenvals(eval)<=0)
      std::cout << "Error: Eigenvalue #" << eval << " = " << eigenvals(eval) << " < 0" << std::endl;
  }
  
  m_data_cor = FSIFitterUtils::CovarianceToCorrelation(&m_data_cov);
  m_data_cov_inv.Invert();
  
  std::cout << "Built covariance matrices for this dataset" << std::endl;
  
}

//***********************************************************************************
// Load covariance martix for DUET separate measurement
//***********************************************************************************
void ExternalDataSet::LoadDUETCovariance(){
  
  // Still need to resize this matrices
  m_data_cov.ResizeTo(nDataPoints,nDataPoints);
  m_data_cov_inv.ResizeTo(nDataPoints,nDataPoints);
  m_data_cor.ResizeTo(nDataPoints,nDataPoints);
  std::cout << "Loading covariance for DUET separate measurement" << std::endl;

  // File where the covariances are
  TFile *fin = new TFile("pion_data/duet_cov_corr.root","OPEN");
  TMatrixDSym *cov = (TMatrixDSym*)fin->Get("duet_covariance");
  TMatrixDSym *cor = (TMatrixDSym*)fin->Get("duet_correlation");

  m_data_cov = (*cov);
  m_data_cov_inv = *cov;
  m_data_cor = *cor;
  m_data_cov_inv.Invert();
  
  std::cout << "Loaded covariance for DUET separate measurement" << std::endl;
  m_data_cov.Print();
  m_data_cov_inv.Print();

}

//***********************************************************************************
// This is called by the constructor 
//***********************************************************************************
void ExternalDataSet::Initialize(){

  ParseDataSetName();
  LoadData();
  BuildDataDiagonalMatrices();
  if(SetName == "DUETSeparate") LoadDUETCovariance();

}
