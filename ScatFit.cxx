//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "ScatFit.hxx"

//***********************************************************************************
// Parse the arguments for running the fit
//***********************************************************************************
int getArgs(int argc, char* argv[]){
  // nuclFit=0;
  std::cout << "Parsing arguments ";

  if(argc == 1){
    
    std::cout<<" ... no options selected: running default option ";
    kBuildGridFromScratch = true;
    kOctave = false;

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
    case 't':
      kOctave = true;
      ++argv; --argc;
      break;
    case 'f':
      nuclFit = true;
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
  std::vector< std::string > DataSetsToBeAdded;
  if (nuclFit){
    std::cout<<"Nucleon Fit data---"<<std::endl;
    DataSetsToBeAdded = {
              // "c_reac_p_ch55",
              "c_reac_p_mi54",
              "c_reac_p_re72",
              "c_tot_p_ca54",
              "c_tot_p_ja78",
              "c_tot_p_ma53",
              "c_tot_p_sc79",
              "c_elas_p_147",
              "c_elas_p_131",
              "c_elas_p_58",
              "c_elas_p_33",
              "c_elas_p_126",
              "c_elas_p_24",
              "c_elas_p_102",
              "c_elas_p_124"
            };
            
  } else {
    std::cout<<"Pion Fit data---"<<std::endl;
  // Name list of data sets to be added
  // csv files with these exact name are expected in the /data folder
    DataSetsToBeAdded = {
              "c_abscx_piP_Ashery",
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
						  "c_reac_piM_Binon"
            };
  }
    
  // Loop through list of data sets. Create them and add them to the vector
  for (int dataset = 0; dataset < (int)DataSetsToBeAdded.size(); dataset++){
    
    std::cout << "\nDataset to be added: " << DataSetsToBeAdded[dataset] << std::endl;
    ExternalDataSet::ExternalDataSet fDataSet(DataSetsToBeAdded[dataset], nuclFit);
    AllDataSets.push_back(fDataSet);
    
  }
  
}


int main(int argc, char* argv[]){
  int args = getArgs(argc, argv);
  if(args != 0){
    std::cout << "Usage " << std::endl;
    return 1;
  }

  if(!kBuildGridFromScratch && !kUseInputGrid){
    std::cout << "didn't set input type, -i/-g" << std::endl;
    return 1;
  }

  std::cout << "-----------------------------------------------------------" << std::endl;
  if (nuclFit) std::cout << "Welcome to the nucleon scattering external data fitter - v2016" << std::endl;
  else std::cout << "Welcome to the pion scattering external data fitter - v2016" << std::endl;
  std::cout << "  (Now Including DUET data with correlation information!!)"    << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  

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
      if (kOctave){
        FSIChi2Grid::FSIChi2Grid TestGrid(AllDataSets, nuclFit, 1, 0);
        aFSIChi2Grid = TestGrid;

      } else{
        FSIChi2Grid::FSIChi2Grid TestGrid(AllDataSets, nuclFit, 0, 1); 
        aFSIChi2Grid = TestGrid;

      }      
  }
  
  // Use specified grid. It will be interpolated and the ouput will go to TMinuit below
  if(kUseInputGrid){
    
    // Define grid object from previously calculated grid TTree
    std::cout << "\nWill use existing grid from specified file: " << inputfileGrid << std::endl;
    TFile *fin = TFile::Open(inputfileGrid.c_str());
    TTree *BuiltGridTree = (TTree*)fin->Get("Chi2Grid");
      if (kOctave){
        FSIChi2Grid::FSIChi2Grid TestGrid(*BuiltGridTree, nuclFit, 1, 0);
        aFSIChi2Grid = TestGrid;
      } else{
        FSIChi2Grid::FSIChi2Grid TestGrid(*BuiltGridTree, nuclFit, 0, 1);
        aFSIChi2Grid = TestGrid;       
      }
    // FSIChi2Grid::FSIChi2Grid TestGrid(*BuiltGridTree, nuclFit);
    
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
  std::cout << "after minuit" << std::endl;

  // Define FCN function. 
  // The grid object is passed beacuse the FCN function will use the result 
  // of the interpolation, which is stored there
  FSIFitFCN fcn(&aFSIChi2Grid, nuclFit, kOctave, (!kOctave));
  minuit->SetMinuitFCN(&fcn);
  double startX, startY, startZ;
  // starting values
  if (nuclFit){
    startX = 0.9;
    startY = 1.0;
    startZ = 0.7;
    minuit->SetParameter(0,"reac",startX,0.1,0,0);
    minuit->SetParameter(1,"elas",startY,0.1,0,0);
    minuit->SetParameter(2,"tot",startZ,0.1,0,0); 
  }else{
    startX = 0.9;
    startY = 1.2;
    startZ = 0.7;
    minuit->SetParameter(0,"qe",startX,0.1,0,0);
    minuit->SetParameter(1,"abs",startY,0.1,0,0);
    minuit->SetParameter(2,"cx",startZ,0.1,0,0);
  }

  
  // Set Parameters for TMinuit

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
  std::string filename;
  // Use to run the interpolation on the reac cross section
  if (!nuclFit) filename = "input/scan_result_ieki_LE.root";
  FSIParameterScan aFSIParameterScan(filename, nuclFit);
  aFSIParameterScan.InterpolateXsecs();

}
