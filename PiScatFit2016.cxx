//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "PiScatFit2016.hxx"

//***********************************************************************************
// Parse the arguments for running the fit
//***********************************************************************************
int getArgs(int argc, char* argv[]){

  std::cout << "Parsing arguments ";

  if(argc == 1){
    
    std::cout<<" ... no options selected: running default option ";
    kBuildGridFromScratch = true;

  }

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){
      
    case 'i':
      kBuildGridFromScratch = true;
      inputfileScan = argv[2];
      ++argv; --argc;
      break;
    case 'g':
      kUseInputGrid = true;
      inputfileGrid = argv[2];
      ++argv; --argc;
      break;
    case 'o':
      outputfileFit = argv[2];
      ++argv; --argc;
      break;
    default:
      std:: cout<< "Error: Option unknown" << std::endl;
      return -1;

      ++argv; --argc;
    }
  }
  
  std::cout << "... done parsing arguments " << std::endl;
  
  return 0;
   
}

//***********************************************************************************
// Define list of data sets to be added and compile them as ExternalDataSets
// objects in a vector to be used later
//***********************************************************************************

void AddDataSets(){
  
  // Name list of data sets to be added
  // csv files with these exact name are expected in the /data folder
  std::vector< std::string > DataSetsToBeAdded = {"c_abscx_piP_Ashery",
						  "c_abscx_piP_DUET",
						  "c_abscx_piP_Navon1",
						  "c_abscx_piP_Navon2",
						  "c_abs_piP_Belloti",
						  "c_cx_piP_Ashery",
						  "c_cx_piP_Belloti2",
						  "c_cx_piP_Jones",
						  "c_cx_piP_Navon",
						  "c_inel_piP_Jones",
						  "c_inel_piP_Levenson",
						  "c_reac_piP_Ashery",
						  "c_reac_piP_Saunders",
						  "c_abscx_piM_Ashery",
						  "c_abscx_piM_Miller",
						  "c_abscx_piM_Navon",
						  "c_cx_piM_Ashery",
						  "c_cx_piM_Hilscher",
						  "c_cx_piM_Navon",
						  "c_reac_piM_Ashery",
						  "c_reac_piM_Binon"};
    
  // Loop through list of data sets. Create them and add them to the vector
  for (int dataset = 0; dataset < (int)DataSetsToBeAdded.size(); dataset++){
    
    std::cout << "\nDataset to be added: " << DataSetsToBeAdded[dataset] << std::endl;
    ExternalDataSet::ExternalDataSet fDataSet(DataSetsToBeAdded[dataset]);
    AllDataSets.push_back(fDataSet);
    
  }
  
}


int main(int argc, char* argv[]){

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Welcome to the pion scattering external data fitter - v2016" << std::endl;
  std::cout << "  (Now Including DUET data with correlation information!!)"    << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  
  int args = getArgs(argc, argv);
  if(args != 0){
    std::cout << "Usage " << std::endl;
    return 1;
  }

  // Default option: Build chisquare grid using dataset and MC Scan
  if(kBuildGridFromScratch) {

      AddDataSets();
      
      /*
      // Was trying to pass the MC scan information from here but failed
      std::string filename = "scan_result_ieki_LE.root";
      FSIParameterScan aFSIParameterScan(filename);
      FSIChi2Grid::FSIChi2Grid TestGrid(AllDataSets,aFSIParameterScan);
      */
      
      // Define Chi2 Grid class from datasets
      FSIChi2Grid::FSIChi2Grid TestGrid(AllDataSets);
      aFSIChi2Grid = TestGrid;
      
  }
  
  // Use specified grid. It will be interpolated and the ouput will go to TMinuit below
  if(kUseInputGrid){
    
    // Define grid object from previously calculated grid TTree
    std::cout << "\nWill use existing grid from specified file: " << inputfileGrid << std::endl;
    TFile *fin = TFile::Open(inputfileGrid.c_str());
    TTree *BuiltGridTree = (TTree*)fin->Get("Chi2Grid");
    FSIChi2Grid::FSIChi2Grid TestGrid(*BuiltGridTree);
    aFSIChi2Grid = TestGrid;
    
  }
  
  // Interpolate the grid
  // Specify if using existing grid or if calculating from scratch
  aFSIChi2Grid.InterpolateFiniteGrid(kUseInputGrid);

  // Load TMinuit2 library
  //gSystem->Load("libMinuit2");
  
  std::cout << "\nWill not run TMinuit2 using the interpolated grid" << std::endl;
  
  // Number of paramateres in fit. 
  // Need to make this more general
  int npar = 3;

  // Define Minuit object. 
  // Not using most up to date constructors because this way is much easier
  TFitterMinuit *minuit = new TFitterMinuit(npar);

  // Define FCN function. 
  // The grid object is passed beacuse the FCN function will use the result 
  // of the interpolation, which is stored there
  FSIFitFCN fcn(&aFSIChi2Grid);
  minuit->SetMinuitFCN(&fcn);

  // starting values
  double startX = 0.9;
  double startY = 0.9;
  
  // Set Parameters for TMinuit
  minuit->SetParameter(0,"qe",startX,0.1,0,0);
  minuit->SetParameter(1,"abs",startY,0.1,0,0);
  minuit->SetParameter(2,"cx",startY,0.1,0,0);
  minuit->SetPrintLevel(3);

  // Actually create and run minimizer
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();
  if (iret != 0) {
    return iret;
  }

  // Get and store the covariance matrix
  TMatrixD *matrix = new TMatrixD(npar,npar,minuit->GetCovarianceMatrix());
  TFile *fout = new TFile(outputfileFit.c_str(),"RECREATE");
  matrix->Write("cov_minuit");
  fout->Close();

  std::cout << std::endl << "Minimization done... Good Bye!" << std::endl;
  
  return 0;
}


//***********************************************************************************
// Interpolate the xsecs instead of the chisquare grid
// This is not used
//***********************************************************************************
void InterpolateXsecs(){
  
  // Use to run the interpolation on the reac cross section
  std::string filename = "input/scan_result_ieki_LE.root";
  FSIParameterScan aFSIParameterScan(filename);
  aFSIParameterScan.InterpolateXsecs();

}
