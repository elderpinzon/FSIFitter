//***********************************************************************************
// PiScatFit v2016
// See README file for description and usage
// Authors:
//   Elder Pinzon <elder@yorku.ca>
//   York University
//***********************************************************************************

#include "ScatFit.hxx"
#include "InterpolatedCrossSectionsOctave.hxx"

int main(int argc, char* argv[]){

  gROOT->ProcessLine(".x ~/rootmaclogon.C");

  bool drawModels = false;
  bool drawTN032Envelopes = false;

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
  gROOT->ProcessLine(Form(".! mkdir -p /home/elder/share/tmp/octave_xsec/%s",argv[1]));

  // Output file
  TFile *fOut = new TFile(Form("/home/elder/share/tmp/octave_xsec/%s/plotting_summary.root",argv[1]),"RECREATE");

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
  TH1F *histoFSIThrowsOneSigma[npar];
  for(Int_t i = 0; i<npar; i++){
    histoFSIThrows[i] = new TH1F(parNames[i].Data(),
				 Form("%s;Thrown Parameter",parNames[i].Data()),
				 100,(*bestfit)[i]*0.3,(*bestfit)[i]*1.7);
    
    histoFSIThrowsOneSigma[i] = new TH1F(Form("%s-1sigma",parNames[i].Data()),
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
  for(Int_t iNuclei = 0; iNuclei < (Int_t)NucleiA.size(); iNuclei++){
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
  Int_t nthrows = 250;
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
      histoFSIThrowsOneSigma[i]->Fill(thrownPars[i]);
    }

    for(Int_t iNuclei = 0; iNuclei <(Int_t)Nuclei.size(); iNuclei++){
      for(Int_t ipid = 0; ipid < (Int_t)pids.size(); ipid++){
	for(Int_t iType = 0; iType<nXsecs; iType++){
	  
	  std::vector<Double_t> allMoms = aInterpolatedCrossSection[iNuclei][ipid]->GetAllMoms();
	  TGraph *gxsecsOct = new TGraph((Int_t)allMoms.size());
	  
	  for(Int_t i = 0; i < (Int_t)allMoms.size(); i++){
	    
	    // The first entry should be the momentum
	    params[0] = allMoms[i];
	    Double_t out = aInterpolatedCrossSection[iNuclei][ipid]->GetSplinedGridPoint(iType+1,params);
	    gxsecsOct->SetPoint(i,allMoms[i],out);
	    std::cout << " Nuclei: " << NucleiA[iNuclei]
		      << " pid: " << pids[ipid]
		      << " type: " << iType
		      << " mom: " << allMoms[i]
		      << " mc: " << out
		      << std::endl;

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
				       "Data-MC compatibility;(nom-data)/sqrt(th*th+exp*exp)"
				       ,60,-3.0,3.0);

  for(Int_t iNuclei = 0; iNuclei <(Int_t)Nuclei.size(); iNuclei++){
      
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
	  std::cout << "Mean: " << mean << " sigma: " << sigma << std::endl;      

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
  c_mcDataCompatibility->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/mcDataCompatibility.png",argv[1]));
  c_mcDataCompatibility->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/mcDataCompatibility.eps",argv[1]));

  fOut->cd();
  mcDataCompatibility->Write();

  ////////////////////////////////////////////////////
  // This is the main plotting part of this executable
  ////////////////////////////////////////////////////

  TCanvas *c_reac[6][2][nXsecs]; //[Nuclei][pionType][intrType]
  TCanvas *c_all[6][2]; //[Nuclei][pionType]
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
    
  for(Int_t iNuclei = 0; iNuclei <(Int_t)Nuclei.size(); iNuclei++){
  
    for(Int_t ipid = 0; ipid<2; ipid++){
      
      c_all[iNuclei][ipid] = new TCanvas(Form("%s_%s",Nuclei[iNuclei].Data(),pionType[ipid].Data()),Form("%s_%s",Nuclei[iNuclei].Data(),pionType[ipid].Data()),2400,1200);
      c_all[iNuclei][ipid]->Divide(3,2);
      
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
	const char* hTemplateTitle = Form("%s-%s %s",pionLatex[ipid].Data(),Nuclei[iNuclei].Data(),iTypeString3[iType].Data());
	TH1F *hTemplate = new TH1F(hTemplateTitle,Form("%s;Momentum[MeV/c];#sigma[mb]",hTemplateTitle),10,0,2000);
	
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

	c_reac[iNuclei][ipid][iType]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/%s_%s_%s.png",argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data(),iTypeString[iType].Data()));
	c_reac[iNuclei][ipid][iType]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/%s_%s_%s.eps",argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data(),iTypeString[iType].Data()));

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

	gr_bestfitOct[iNuclei][ipid][iType]->Draw("Csame");
	gr_minOct[iNuclei][ipid][iType]->Draw("Csame");
	gr_maxOct[iNuclei][ipid][iType]->Draw("Csame");

	merged->Draw("Psame");

	// Create shade from the min and max histos
	gr_shadeOct[iNuclei][ipid][iType] = FSIFitterUtils::MergeGraphsIntoEnvelope(gr_minOct[iNuclei][ipid][iType],gr_maxOct[iNuclei][ipid][iType]);
	gr_shadeOct[iNuclei][ipid][iType]->SetFillColor(2);
	//gr_shadeOct[iNuclei][ipid][iType]->SetFillStyle(3344);


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
	delete collection;	  

      }

      if(drawModels){
	// Draw box with legend
	c_all[iNuclei][ipid]->cd(6);
	TPaveText *pt = new TPaveText(0.15,0.15,0.95,.87);
	pt->AddText("Cascade Models");
	pt->AddLine(.0,0.87,0,0.87);
	// pt->AddText("NEUT with current #pm1#sigma band");
	// pt->GetLineWith("NEUT with")->SetTextColor(kMagenta+1);
	pt->AddText("NEUT");
	pt->GetLineWith("NEUT")->SetTextColor(kMagenta+1);
	if(drawTN032Envelopes){
	  pt->AddText("NEUT new #pm1#sigma band");
	  pt->GetLineWith("NEUT new")->SetTextColor(kCyan+1);
	}
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
	pt->Draw();//"same");
      }
	  
      c_all[iNuclei][ipid]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/%s_%s_all.png",argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));
      c_all[iNuclei][ipid]->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/%s_%s_all.eps",argv[1],Nuclei[iNuclei].Data(),pionType[ipid].Data()));

    }
  }
  
  std::cout << "Total numer of data points plotted: " << counter << std::endl;
  
  
  // Draw and store throw histograms
  TCanvas *c_throws = new TCanvas("c_throws","c_throws",1800,600);
  c_throws->Divide(3,2);
  for(Int_t j = 0; j<npar; j++){
    c_throws->cd(j+1);
    //gStyle->SetOptStat(0111);
    gStyle->SetOptFit();
    histoFSIThrows[j]->Fit("gaus");
    histoFSIThrows[j]->Draw();
    histoFSIThrowsOneSigma[j]->SetLineColor(2);
    histoFSIThrowsOneSigma[j]->Draw("same");
    fOut->cd();
    histoFSIThrows[j]->Write();
    histoFSIThrowsOneSigma[j]->Write();
  }
  
  c_throws->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/throw.png",argv[1]));
  c_throws->SaveAs(Form("/home/elder/share/tmp/octave_xsec/%s/throw.eps",argv[1]));

  fOut->Write();
  c_throws->Write();

  fOut->Close();

  return 0;
}
