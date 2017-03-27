//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "ScatFit.hxx"

int main(int argc, char* argv[]){

  gROOT->ProcessLine(".x ~/rootmaclogon.C");

  TString outdir = "/home/elder/share/tmp/octave_xsec";

  bool drawModels = false;
  bool drawTN032Envelopes = true;//false;

  // Read input best fit and covariance
  if(argc<2){
    std::cout << "Usage: ./PlotResult.exe input-file.root" << std::endl;
    std::exit(1);
  }

  TFile *fin = new TFile(argv[1],"OPEN");
  TMatrixDSym *cov = (TMatrixDSym*)fin->Get("cov");
  TVectorD *bestfit = (TVectorD*)fin->Get("BestFitPars");
  std::cout << "Plotting envelopes from following FSI covariance in file: "
	    << argv[1]
	    << std::endl;
  cov->Print();
  
  // Make output folder
  gROOT->ProcessLine(Form(".! mkdir -p %s/%s",outdir.Data(),argv[1]));

  // Output file
  TFile *fOut = new TFile(Form("%s/%s/plotting_summary.root",outdir.Data(),argv[1]),"RECREATE");

  // By default do all nuclei, polarities, and interaction channels
  Int_t nXsecs = xsName.size();
  Int_t nNuclei = Nuclei.size();
  Int_t nPID = pionType.size();
  Int_t nMaxMomenta = 50;
  printf("nNuclei: %d nXsecs: %d nPID: %d\n", nNuclei, nXsecs, nPID);
  
  // Load true (regenerated) xsecs for best pi-C fit
  // this is (was) useful to check the interpolation accuracy
  TFile *ftrue = TFile::Open("input/true_best_fit_c.root");
  TH1D *hTrue[nXsecs];
  for(Int_t i = 0; i < nXsecs; i++){
    hTrue[i] = (TH1D*)ftrue->Get(Form("h%s",xsName[i].Data()));
    hTrue[i]->SetLineWidth(1);
    hTrue[i]->SetLineColor(4);
  }

  // Link correct folder with pion_data selection
  gROOT->ProcessLine(".! ln -s pion_data_plots pion_data");

  // Add all data sets
  AddDataSets("c_");
  AddDataSets("o_");
  AddDataSets("al_");
  AddDataSets("fe_");
  AddDataSets("cu_");
  AddDataSets("pb_");
  
  // Number of FSI parameters
  Int_t npar = 5;//7;

  // Histograms to store thrown parameter values and initialization
  TH1F *histoFSIThrows[npar];
  TH1F *histoFSIThrowsUsed[npar];
  for(Int_t i = 0; i<npar; i++){
    histoFSIThrows[i] = new TH1F(parNames[i].Data(),
				 Form("%s;Thrown Parameter",parNames[i].Data()),
				 100,(*bestfit)[i]*0.3,(*bestfit)[i]*1.7);
    
    histoFSIThrowsUsed[i] = new TH1F(Form("%s-used",parNames[i].Data()),
					 Form("%s;Thrown Parameter",parNames[i].Data()),
					 100,(*bestfit)[i]*0.3,(*bestfit)[i]*1.7);
  }

  // Max scales for plotting each histogram (horrible)
  std::vector<int> histoMax = {550,340,250,100,300};

  // TGraphs to store the interpolated cross sections used
  // to draw the +- 1 sigma lines and the envelopes
  TGraphErrors *gr_bestfitOct[nNuclei][nPID][nXsecs];
  TGraphErrors *gr_minOct[nNuclei][nPID][nXsecs];
  TGraphErrors *gr_maxOct[nNuclei][nPID][nXsecs];
  TGraphErrors *gr_shadeOct[nNuclei][nPID][nXsecs];
  
  // Histos to build xs PDFS on each bin which are then fitted to get 1-sigma envelope
  TH1D *xsPDFs[nNuclei][nPID][nXsecs][nMaxMomenta];
  TCanvas *xsPDFCanvas[nNuclei][nPID][nXsecs][nMaxMomenta];

  // ThrowParms does the Cholesky decompostion of the covariance matrix
  // and gets ready to throw a set of correlated random parameters
  ThrowParms *throwObject;
  throwObject = new ThrowParms((*bestfit),(*cov));
  throwObject->SetSeed(1000);

  // Octave engine initialization
  string_vector argv_oct (2);
  argv_oct(0) = "FSIOctaveInterpolation";
  argv_oct(1) = "-q";
  octave_main (2, argv_oct.c_str_vec (), 1);

  // Initialize splined cross sections
  InterpolatedCrossSectionsOctave *aInterpolatedCrossSection[6][2];
  for(Int_t iNuclei = 0; iNuclei < nNuclei; iNuclei++){
    for(Int_t ipid = 0; ipid < (Int_t)pids.size(); ipid++){
      
      aInterpolatedCrossSection[iNuclei][ipid] = new InterpolatedCrossSectionsOctave(ScanFileName,NucleiA[iNuclei],pids[ipid]);
      aInterpolatedCrossSection[iNuclei][ipid]->BuildGrid();
      
      // Initialize TGraphs here
      for(Int_t iXsec = 0; iXsec<nXsecs; iXsec++){
	
	Int_t nMom = aInterpolatedCrossSection[iNuclei][ipid]->GetNumberMomenta();
	
	gr_bestfitOct[iNuclei][ipid][iXsec] = new TGraphErrors(nMom);
	gr_minOct[iNuclei][ipid][iXsec] = new TGraphErrors(nMom);
	gr_maxOct[iNuclei][ipid][iXsec] = new TGraphErrors(nMom);
	gr_shadeOct[iNuclei][ipid][iXsec] = new TGraphErrors(2*nMom);
	
	gr_minOct[iNuclei][ipid][iXsec]->SetLineColor(kGreen+2);
	gr_maxOct[iNuclei][ipid][iXsec]->SetLineColor(kGreen+2);
	
      }
    }
  }
  
  std::cout << "Done loading cross sections to octave objects" << std::endl;
  
  // How many valid throws to do?
  // Note: Throws outside the grid are not counted
  Int_t nthrows = 1000;
  Int_t nOutside = 0;
  Int_t tt = 0;
  
  while(tt < nthrows){
    
    std::cout<< "Throw # " << tt << std::endl;
    stopwatch.Start();
    
    bool outsideGrid = false;
    std::vector<Double_t> thrownPars;
    throwObject->ThrowSet(thrownPars);
    std::vector<Double_t> params;
    
    // Count drawing best fit as first throw
    if(tt == 0){
      std::cout << "Starting with bestfit pars: " << std::endl;
      for(Int_t i = 0; i<(Int_t)thrownPars.size(); i++)
	params.push_back((*bestfit)[i]);
    }
    // Then do throws
    else{
      
      for(Int_t i = 0; i<(Int_t)thrownPars.size(); i++){
	
	params.push_back(thrownPars[i]);

	histoFSIThrows[i]->Fill(thrownPars[i]);
	printf("%s: min:%.4f throw:%.4f max: %.4f\n",parNames[i].Data(),FSIParsMin[i],thrownPars[i],FSIParsMax[i]);
	// Check if thrown parameter is inside the Grid
	if(thrownPars[i] < FSIParsMin[i] || thrownPars[i] > FSIParsMax[i]){
	  
	  // If outside grid, must exit otherwise 999s will screw up envelopes
	  std::cout << "Outside grid! Skipping Throw!" << std::endl;
	  nOutside++;
	  outsideGrid = true;
	}
      }
    }
    
    // If outside grid, must exit otherwise 999s will screw up envelopes
    if(outsideGrid){
      std::cout << "Outside grid! Skipping Throw!" << std::endl;
      continue;
    }


    // Need to add extra entry to vector 'params' for mom value
    // as that is the format expected by InterpolatedCrossSectionOctave
    params.insert(params.begin(),0);
    

    for(Int_t i = 0; i<(Int_t)thrownPars.size(); i++){
      histoFSIThrowsUsed[i]->Fill(thrownPars[i]);
    }

    for(Int_t iNuclei = 0; iNuclei < nNuclei; iNuclei++){
      for(Int_t ipid = 0; ipid < (Int_t)pids.size(); ipid++){
	for(Int_t iType = 0; iType<nXsecs; iType++){
	  
	  std::vector<Double_t> allMoms = aInterpolatedCrossSection[iNuclei][ipid]->GetAllMoms();
	  TGraph *gxsecsOct = new TGraph((Int_t)allMoms.size());
	  
	  for(Int_t i = 0; i < (Int_t)allMoms.size(); i++){
	    
	    // The first entry should be the momentum
	    params[0] = allMoms[i];
	    stopwatch_2.Start();
	    Double_t out = aInterpolatedCrossSection[iNuclei][ipid]->GetSplinedGridPoint(iType+1,params);
	    stopwatch_2.Stop();
	    gxsecsOct->SetPoint(i,allMoms[i],out);
	    // std::cout << " Nuclei: " << NucleiA[iNuclei]
	    // 	      << " pid: " << pids[ipid]
	    // 	      << " type: " << iType
	    // 	      << " mom: " << allMoms[i]
	    // 	      << " mc: " << out
	    // 	      << ". Took " << stopwatch_2.RealTime() << " seconds."
	    // 	      << std::endl;
	    if(out == 9999)
	      std::exit(1);
	    
	    // Draw best fit
	    if(tt == 0){
	    
	      // Initialize PDF histos
	      char *title = Form("PDF_%s_%s_%s_%.2f_MeV",
				 Nuclei[iNuclei].Data(),
				 pionType[ipid].Data(),
				 xsName[iType].Data(),
				 allMoms[i]);

	      xsPDFs[iNuclei][ipid][iType][i] = new TH1D(title,title,50,0.5*out,1.5*out);
	      gr_bestfitOct[iNuclei][ipid][iType]->SetPoint(i,allMoms[i],out);

	    }

	    xsPDFs[iNuclei][ipid][iType][i]->Fill(out);
	    
	  }
	}
      }
    }
    
    // We made all the way here, so let's count this throw
    tt++;
    stopwatch.Stop();
    std::cout << "Took " << stopwatch.RealTime() << " seconds." << std::endl;
  
  }
  
  std::cout << nOutside << " throws outside the grid." << std::endl;
  
  // Fit PDFs and fill TGraphs with  \pm 1-\sigma envelopes
  // Also calculate  (nom-data)/sqrt(th*th + exp*exp) for error inflation
  TH1F *mcDataCompatibility = new TH1F("mcDataCompatibility",
				       "Data-MC compatibility;(Best fit - Data)/(th*th+exp*exp)^{1/2};Numer of Data points"
				       ,26,-5.0,5.0);

  for(Int_t iNuclei = 0; iNuclei < nNuclei; iNuclei++){
      
    for(Int_t ipid = 0; ipid<2; ipid++){
    
      std::vector<Double_t> allMoms = aInterpolatedCrossSection[iNuclei][ipid]->GetAllMoms();
      Int_t nMom = aInterpolatedCrossSection[iNuclei][ipid]->GetNumberMomenta();
      
      for(Int_t iType = 0; iType<nXsecs; iType++){
	
	for(Int_t i = 0; i < nMom; i++){
	  
	  // Fit PDF to Gaussian
	  xsPDFs[iNuclei][ipid][iType][i]->Fit("gaus","Q");
	  
	  // Store them and print them right away
	  char *title = Form("PDF_%s_%s_%s_%.2f_MeV",
			     Nuclei[iNuclei].Data(),
			     pionType[ipid].Data(),
			     xsName[iType].Data(),
			     allMoms[i]);
	  xsPDFCanvas[iNuclei][ipid][iType][i] = new TCanvas(title,title);
	  xsPDFCanvas[iNuclei][ipid][iType][i]->cd();
	  xsPDFs[iNuclei][ipid][iType][i]->Draw();
	  fOut->cd();
	  xsPDFs[iNuclei][ipid][iType][i]->Write();
	  //xsPDFCanvas[iNuclei][ipid][iType][i]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/pdf_%s.png",argv[1],title));
	  //xsPDFCanvas[iNuclei][ipid][iType][i]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/pdf_%s.eps",argv[1],title));

	  // Get mean and sigma
	  Double_t mean = xsPDFs[iNuclei][ipid][iType][i]->GetFunction("gaus")->GetParameter(1);
	  Double_t sigma = xsPDFs[iNuclei][ipid][iType][i]->GetFunction("gaus")->GetParameter(2);
	  Double_t bestfit = gr_bestfitOct[iNuclei][ipid][iType]->Eval(allMoms[i]);
	  std::cout << "Bestfit: " << bestfit << " Mean: " << mean << " sigma: " << sigma << std::endl;

	  // Just check that the mean of the PDF is indeed close to the best fit
	  if(TMath::Abs(bestfit-mean)/bestfit > 0.05)
	    std::cout << "Warning more than 5% difference" << std::endl;

	  // Set envelopes
	  gr_minOct[iNuclei][ipid][iType]->SetPoint(i,allMoms[i],mean-sigma);
	  gr_maxOct[iNuclei][ipid][iType]->SetPoint(i,allMoms[i],mean+sigma);
      
	  // Loop over datasets to calculate (nom-data)/sqrt(th*th + exp*exp)
	  // for each data point that matches this PDF
	  for (Int_t dataset = 0; dataset < (Int_t)AllDataSets.size(); dataset++){
	    
	    // Loop over points
	    for(Int_t point =  0; point < AllDataSets[dataset].GetNumberOfDataPoints(); point++){
	      
	      // Check if matches
	      if(AllDataSets[dataset].GetNucleus() == NucleiA[iNuclei] &&
		 AllDataSets[dataset].GetPionType() == pids[ipid] &&
		 AllDataSets[dataset].GetIntrType() == (iType+1) &&
		 TMath::Abs(AllDataSets[dataset].GetMomVector()[point] - allMoms[i]) < 0.1){
		
		
		Double_t comparison = (mean - AllDataSets[dataset].GetXsecVector()[point]) / TMath::Sqrt(sigma*sigma + AllDataSets[dataset].GetXsecErrorVector()[point] * AllDataSets[dataset].GetXsecErrorVector()[point]);
		mcDataCompatibility->Fill(comparison);
		printf("A: %d p: %d int: %d mom: %.2f xsec: %.2f err: %.2f mean: %.2f sigma: %.2f comp :%.4f\n",
		       AllDataSets[dataset].GetNucleus(),
		       AllDataSets[dataset].GetPionType(),
		       AllDataSets[dataset].GetIntrType(),
		       AllDataSets[dataset].GetMomVector()[point],
		       AllDataSets[dataset].GetXsecVector()[point],
		       AllDataSets[dataset].GetXsecErrorVector()[point],
		       mean, sigma, comparison);
		
	      }
	    }
	  }
	}
      }
    }
  }

  TCanvas *c_mcDataCompatibility = new TCanvas();
  c_mcDataCompatibility->cd();
  mcDataCompatibility->Draw();
  c_mcDataCompatibility->SaveAs(Form("%s/%s/mcDataCompatibility.png",outdir.Data(),argv[1]));
  c_mcDataCompatibility->SaveAs(Form("%s/%s/mcDataCompatibility.eps",outdir.Data(),argv[1]));

  fOut->cd();
  mcDataCompatibility->Write();

  ////////////////////////////////////////////////////
  // This is the main plotting part of this executable
  ////////////////////////////////////////////////////

  TCanvas *c_reac[6][2][nXsecs]; //[Nuclei][pionType][intrType]
  TCanvas *c_all[6][2]; //[Nuclei][pionType
  TCanvas *c_all_ratio[6][2]; //[Nuclei][pionType]
  Int_t counter = 0;

  // Load model predictions from other generators
  ModelPrediction::ModelPrediction geant4("geant4_xs.root");
  geant4.SetLineColor(kBlack);
  
  ModelPrediction::ModelPrediction genie_hA("genie_hA_xs.root");
  genie_hA.SetGeVtoMeV();
  genie_hA.SetLineColor(kRed);
  
  ModelPrediction::ModelPrediction genie_hA2014("genie_hA2014_xs.root");
  genie_hA2014.SetGeVtoMeV();
  genie_hA2014.SetLineColor(kBlue);
  
  ModelPrediction::ModelPrediction genie_hN2015("genie_hN2015.root");
  genie_hN2015.SetGeVtoMeV();
  genie_hN2015.SetLineColor(kYellow+1);
  
  ModelPrediction::ModelPrediction nuwro("nuwro_scattering_xs.root");
  nuwro.SetGeVtoMeV();
  nuwro.SetLineColor(kGreen+2);
  
  ModelPrediction::ModelPrediction fluka("FLUKA_pipC.root");
  fluka.SetGeVtoMeV();
  fluka.SetLineColor(kCyan);//Orange+10);
    
  for(Int_t iNuclei = 0; iNuclei < nNuclei; iNuclei++){
  
    for(Int_t ipid = 0; ipid<2; ipid++){
      
      c_all[iNuclei][ipid] = new TCanvas(Form("%s_%s",Nuclei[iNuclei].Data(),pionType[ipid].Data()),Form("%s_%s",Nuclei[iNuclei].Data(),pionType[ipid].Data()),2400,1200);
      c_all[iNuclei][ipid]->Divide(3,2);

      c_all_ratio[iNuclei][ipid] = new TCanvas(Form("%s_%s_ratio",Nuclei[iNuclei].Data(),pionType[ipid].Data()),Form("%s_%s_ratio",Nuclei[iNuclei].Data(),pionType[ipid].Data()),1000,1600);
      
      for(Int_t iType = 0; iType <nXsecs; iType++){

	std::cout << "\nWorking on nuclei: " << Nuclei[iNuclei]
		  << " pionType: " << pionType[ipid]
		  <<" intrType: " << iTypeString[iType]
		  <<std::endl;

	TList *collection = new TList;
	
	// Loop through list of data sets and add them to the TCollection
	
	for (Int_t dataset = 0; dataset < (Int_t)AllDataSets.size(); dataset++){
	  
	  if(AllDataSets[dataset].GetNucleusString() == Nuclei[iNuclei]
	     && AllDataSets[dataset].GetPionTypeString() == pionType[ipid]
	     && AllDataSets[dataset].GetIntrType() == iType+1){
	    collection->Add(AllDataSets[dataset].GetTGraph());
	    counter += AllDataSets[dataset].GetNumberOfDataPoints();
	  }
	}
	
	std::cout << "Collection has "<< collection->GetSize() << " elements." << std::endl;
	
	TGraphErrors *merged = new TGraphErrors();
	//merged->Merge(collection); // This method loses the error bars in the merging process.
	merged = FSIFitterUtils::MergeGraphs(collection);
	//merged->SetMarkerStyle(kFullCircle);
	merged->SetMarkerStyle(kCircle);
	merged->SetMarkerSize(1.6);
	//merged->SetLineWidth(2);
	merged->Sort();	

	c_reac[iNuclei][ipid][iType]= new TCanvas(Form("%s_%s_%s",Nuclei[iNuclei].Data(),pionType[ipid].Data(),iTypeString[iType].Data()),Form("%s-%s %s",pionLatex[ipid].Data(),Nuclei[iNuclei].Data(),iTypeString[iType].Data()));
	const char* hTemplateTitle = Form("%s %s %s",pionLatex[ipid].Data(),NucleiFull[iNuclei].Data(),iTypeString3[iType].Data());
	TH1F *hTemplate = new TH1F(hTemplateTitle,Form("%s;Momentum [MeV/c];#sigma [mb]",hTemplateTitle),10,0,2000);
	
	TN032Envelopes *tn032Envel = new TN032Envelopes (Nuclei[iNuclei],pionCode[ipid],iTypeString[iType]);

	// Get max scale for template histogram
	Double_t max_scale = tn032Envel->GetMaxScale();
	if(merged->GetN() > 0 && TMath::MaxElement(merged->GetN(),merged->GetY())*1.15 > max_scale)
	  max_scale = TMath::MaxElement(merged->GetN(),merged->GetY());
	
	hTemplate->SetLineColor(0);
	hTemplate->Draw();
	hTemplate->GetYaxis()->SetRangeUser(0,max_scale);
	
	if(drawTN032Envelopes){
	  tn032Envel->GetCrossSectionShade()->Draw("fsame");
	  tn032Envel->GetCrossSectionMin()->Draw("Csame");
	  tn032Envel->GetCrossSectionMax()->Draw("Csame");
	}
	
	if(drawModels){
	  tn032Envel->GetCrossSectionNominal()->Draw("Csame");
	  geant4.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hA.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hA2014.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hN2015.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  nuwro.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  fluka.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	}

	//hTrue[iType]->Draw("samehist");
	
	//gr_bestfit[iNuclei][ipid][iType]->Draw("Csame");
	gr_bestfitOct[iNuclei][ipid][iType]->Draw("Csame");
	gr_minOct[iNuclei][ipid][iType]->Draw("Csame");
	gr_maxOct[iNuclei][ipid][iType]->Draw("Csame");

	merged->Draw("Psame");

	c_reac[iNuclei][ipid][iType]->SaveAs(Form("%s/%s/%s_%s_%s.png",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data(),iTypeString[iType].Data()));
	c_reac[iNuclei][ipid][iType]->SaveAs(Form("%s/%s/%s_%s_%s.eps",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data(),iTypeString[iType].Data()));

	// Draw composite 'c_all' histogram
	c_all[iNuclei][ipid]->cd(iType+1);
	hTemplate->Draw();

	if(drawTN032Envelopes){
	  tn032Envel->GetCrossSectionShade()->Draw("fsame");
	  tn032Envel->GetCrossSectionMin()->Draw("Csame");
	  tn032Envel->GetCrossSectionMax()->Draw("Csame");
	}
	
	if(drawModels){
	  tn032Envel->GetCrossSectionNominal()->Draw("Csame");
	  geant4.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hA.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hA2014.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  genie_hN2015.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  nuwro.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	  fluka.GetCrossSection(Nuclei2[iNuclei],pionType2[ipid],iTypeString2[iType])->Draw("Chistsame");
	}

	//hTrue[iType]->Draw("samehist");

	// Create shade from the min and max histos
	gr_shadeOct[iNuclei][ipid][iType] = FSIFitterUtils::MergeGraphsIntoEnvelope(gr_minOct[iNuclei][ipid][iType],gr_maxOct[iNuclei][ipid][iType]);
	gr_shadeOct[iNuclei][ipid][iType]->SetFillColor(kRed+1);
	//gr_shadeOct[iNuclei][ipid][iType]->SetFillStyle(3344);
	gr_shadeOct[iNuclei][ipid][iType]->Draw("fsame");

	gr_bestfitOct[iNuclei][ipid][iType]->Draw("Csame");
	// gr_minOct[iNuclei][ipid][iType]->Draw("Csame");
	// gr_maxOct[iNuclei][ipid][iType]->Draw("Csame");

	merged->Draw("Psame");

	// Draw composite c_all_ratio
	c_all_ratio[iNuclei][ipid]->cd();
	const char* hTemplateTitleRatio = Form("%s_ratio",hTemplateTitle);
	TH1F *hTemplateRatio = new TH1F(hTemplateTitleRatio,Form(";Momentum [MeV/c];#sigma_{%s}/#sigma_{%s, Best fit}",xsName[iType].Data(),xsName[iType].Data()),10,0,2000);
	hTemplateRatio->SetLineColor(0);
	hTemplateRatio->GetYaxis()->SetTitleOffset(0.5);
	hTemplateRatio->GetYaxis()->SetTitleSize(0.1);
	hTemplateRatio->GetYaxis()->SetLabelSize(0.06);
	hTemplateRatio->GetXaxis()->SetTitleSize(0.07);
	hTemplateRatio->GetXaxis()->SetLabelSize(0.06);
	hTemplateRatio->GetYaxis()->SetRangeUser(0,2.0);

	TPad *pad = new TPad("pad", "pad", 0, 1.0 - (iType+1)*0.2, 1, 1.0 - iType*0.2);
	pad->SetLeftMargin(0.1);
	if(iType != 4) pad->SetBottomMargin(0); // Upper and lower plot are joined
	if(iType != 0) pad->SetTopMargin(0);
	pad->Draw();
	pad->cd();
	
	hTemplateRatio->Draw();

	// TGraphError with ratio band
	TGraphErrors *gr_shadeOctRatio = new TGraphErrors(gr_shadeOct[iNuclei][ipid][iType]->GetN());
	Double_t *nxShade = gr_shadeOct[iNuclei][ipid][iType]->GetX();
        Double_t *nyShade = gr_shadeOct[iNuclei][ipid][iType]->GetY();
	for(Int_t imom = 0; imom < gr_shadeOct[iNuclei][ipid][iType]->GetN(); imom++){
	  
	  Double_t bestfit = gr_bestfitOct[iNuclei][ipid][iType]->Eval(nxShade[imom]);
	  Double_t ratio = nyShade[imom]/bestfit;
	  gr_shadeOctRatio->SetPoint(imom,nxShade[imom],ratio);
	}
	gr_shadeOctRatio->SetFillColor(kRed+1);
	gr_shadeOctRatio->Draw("fsame");
	
	// TGraphErrors with ratio of data points to best fit
	TGraphErrors *merged_ratio = new TGraphErrors(merged->GetN());
	merged_ratio->SetMarkerStyle(kCircle);
	merged_ratio->SetMarkerSize(1.6);
	merged->Sort();	

	Double_t *nx = merged->GetX();
	Double_t *ny = merged->GetY();

	for(Int_t point = 0; point < merged->GetN(); point++){
	  
	  Double_t bestfit = gr_bestfitOct[iNuclei][ipid][iType]->Eval(nx[point]);
	  Double_t ratio = ny[point]/bestfit;
	  merged_ratio->SetPoint(point,nx[point],ratio);
	  Double_t ratio_error = merged->GetErrorY(point)/bestfit;
	  merged_ratio->SetPointError(point,0,ratio_error);
	  	  
	}
	merged_ratio->Draw("Psame");

	const char* latexTitle = Form("%s %s %s",pionLatex[ipid].Data(),NucleiFull[iNuclei].Data(),xsName[iType].Data());
	TLatex *l = new TLatex(0.7,0.8,latexTitle);
	l->SetNDC();
	l->Draw("same");

	//Store envelopes to files here
	fOut->cd();
	char *title = Form("%s_%s_%s",
			   Nuclei[iNuclei].Data(),
			   pionType[ipid].Data(),
			   xsName[iType].Data());

	gr_bestfitOct[iNuclei][ipid][iType]->Write(Form("gr_bestfitOct_%s",title));
        gr_minOct[iNuclei][ipid][iType]->Write(Form("gr_minOct_%s",title));
        gr_maxOct[iNuclei][ipid][iType]->Write(Form("gr_maxOct_%s",title));
        gr_shadeOct[iNuclei][ipid][iType]->Write(Form("gr_shadeOct_%s",title));
	gr_shadeOctRatio->Write(Form("gr_shadeOctRatio_%s",title));
	merged->Write(Form("data_%s",title));
	merged_ratio->Write(Form("dataRatio_%s",title));
	delete collection;	  

      }


      TPaveText *pt = new TPaveText(0.15,0.15,0.95,.87);

      pt->AddText("FSIFitter best fit");
      pt->AddText("FSIFitter #pm1#sigma band");
      pt->GetLineWith("FSIFitter #pm")->SetTextColor(kRed+1);

      if(drawTN032Envelopes){
	pt->AddText("TN-032 #pm1#sigma band");
	pt->GetLineWith("TN-032")->SetTextColor(kAzure+2);
      }

      if(drawModels){
	// Draw box with legend
	if(drawModels){
	  pt->AddText("Geant4 Bertini");
	  pt->GetLineWith("Geant4")->SetTextColor(kBlack);
	  pt->AddText("GENIE hA");
	  pt->GetLineWith("GENIE hA")->SetTextColor(kRed);
	  pt->AddText("GENIE hA2014");
	  pt->GetLineWith("2014")->SetTextColor(kBlue);
	  pt->AddText("GENIE hN2015");
	  pt->GetLineWith("2015")->SetTextColor(kYellow+1);
	  pt->AddText("NuWro");
	  pt->GetLineWith("NuWro")->SetTextColor(kGreen+2);
	  pt->AddText("FLUKA");
	  pt->GetLineWith("FLUKA")->SetTextColor(kCyan);//+10);
	}

	c_all[iNuclei][ipid]->cd(6);
	pt->Draw();//"same");
      }
      
      fOut->cd();
      c_all[iNuclei][ipid]->Write();
      c_all[iNuclei][ipid]->SaveAs(Form("%s/%s/%s_%s_all.png",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));
      c_all[iNuclei][ipid]->SaveAs(Form("%s/%s/%s_%s_all.eps",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));

      c_all_ratio[iNuclei][ipid]->SaveAs(Form("%s/%s/%s_%s_all_ratio.png",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));
      c_all_ratio[iNuclei][ipid]->SaveAs(Form("%s/%s/%s_%s_all_ratio.eps",outdir.Data(),argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));

    }
  }
  
  std::cout << "Total numer of data points plotted: " << counter << std::endl;
  
  
  // Draw and store throw histograms
  TCanvas *c_throws = new TCanvas("c_throws","c_throws",2400,1200);
  c_throws->Divide(3,2);
  for(Int_t j = 0; j<npar; j++){
    c_throws->cd(j+1);
    //gStyle->SetOptStat(0111);
    gStyle->SetOptFit();
    histoFSIThrows[j]->SetLineWidth(1);
    histoFSIThrows[j]->Fit("gaus");
    histoFSIThrowsUsed[j]->SetLineWidth(1);
    histoFSIThrowsUsed[j]->SetLineColor(2);
    histoFSIThrowsUsed[j]->Fit("gaus");
    histoFSIThrows[j]->Draw();
    histoFSIThrowsUsed[j]->Draw("same");
    fOut->cd();
    histoFSIThrows[j]->Write();
    histoFSIThrowsUsed[j]->Write();
  }
  
  c_throws->SaveAs(Form("%s/%s/throws.png",outdir.Data(),argv[1]));
  c_throws->SaveAs(Form("%s/%s/throws.eps",outdir.Data(),argv[1]));

  c_throws->Write();

  fOut->Close();

  return 0;
}
