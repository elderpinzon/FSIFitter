//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "NuclScatFit2016.hxx"

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
  std::vector< std::string > DataSetsToBeAdded = {
              "c_reac_p_mi54",
              "c_reac_p_1re72",
              "c_reac_p_2re72",
              "c_reac_p_3re72",
              "c_reac_p_4re72",
              "c_reac_p_5re72",
              "c_reac_p_6re72",
              "c_reac_p_7re72",
              "c_reac_p_8re72",
              // "c_reac_p_ch55",
              "c_tot_p_ca54",
              // "c_tot_p_ja78",
              "c_tot_p_ma53",
              "c_tot_p_1sc79",
              "c_tot_p_2sc79",
              "c_tot_p_3sc79",
              "c_tot_p_4sc79",
              "c_tot_p_5sc79",
              "c_tot_p_6sc79",
              "c_tot_p_7sc79",
              // "c_tot_p_8sc79",
              "c_tot_p_9sc79",
              "c_tot_p_10sc79",
              "c_tot_p_11sc79",
              // "c_tot_p_12sc79",
              // "c_tot_p_13sc79",
              "c_tot_p_14sc79",
              "c_tot_p_15sc79",
              "c_tot_p_16sc79",
              "c_tot_p_17sc79",
              "c_tot_p_18sc79",
              "c_tot_p_19sc79",
              "c_tot_p_20sc79",
              "c_elas_p_147",
              "c_elas_p_131",
              "c_elas_p_58",
              "c_elas_p_33",
              "c_elas_p_126",
              "c_elas_p_24",
              "c_elas_p_102",
              "c_elas_p_124"
            };
    
  // Loop through list of data sets. Create them and add them to the vector
  for (int dataset = 0; dataset < (int)DataSetsToBeAdded.size(); dataset++){
    
    std::cout << "\nDataset to be added: " << DataSetsToBeAdded[dataset] << std::endl;
    NuclExternalDataSet::NuclExternalDataSet fDataSet(DataSetsToBeAdded[dataset]);
    AllDataSets.push_back(fDataSet);
    
  }
  
}


int main(int argc, char* argv[]){

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Welcome to the nucleon scattering external data fitter - v2016" << std::endl;
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
      NuclFSIChi2Grid::NuclFSIChi2Grid TestGrid(AllDataSets);
      aFSIChi2Grid = TestGrid;
      
  }
  
  // Use specified grid. It will be interpolated and the ouput will go to TMinuit below
  if(kUseInputGrid){
    
    // Define grid object from previously calculated grid TTree
    std::cout << "\nWill use existing grid from specified file: " << inputfileGrid << std::endl;
    TFile *fin = TFile::Open(inputfileGrid.c_str());
    TTree *BuiltGridTree = (TTree*)fin->Get("Chi2Grid");
    NuclFSIChi2Grid::NuclFSIChi2Grid TestGrid(*BuiltGridTree);
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
  NuclFSIFitFCN fcn(&aFSIChi2Grid);
  minuit->SetMinuitFCN(&fcn);

  // starting values
  double startX = 0.9;
  double startY = 1.0;
  double startZ = 0.7;
  
  // Set Parameters for TMinuit
  minuit->SetParameter(0,"reac",startX,0.1,0,0);
  minuit->SetParameter(1,"elas",startY,0.1,0,0);
  minuit->SetParameter(2,"tot",startZ,0.1,0,0);
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
// void InterpolateXsecs(){
  
//   // Use to run the interpolation on the reac cross section
//   std::string filename = "input/scan_result_ieki_LE.root";
//   FSIParameterScan aFSIParameterScan(filename);
//   aFSIParameterScan.InterpolateXsecs();

// }
