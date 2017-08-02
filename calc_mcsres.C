//DEFINITIONS
//Momentum if for whole track so don't need to separate planes
//Use best plane for the pid
//Produce 2D histograms of calculated momentum vs true momentum
//For range-based and MCS momentum estimation

TH2F *h_McsTrue[3];
TH2F *h_RangeTrue[3];
TH2F *h_McsRange[3];

TH2F *h_McsTrueDiff[3];
TH2F *h_RangeTrueDiff[3];

TH1F *h_McsResolution[3];
TH1F *h_McsBias[3];
TH1F *h_RangeResolution[3];
TH1F *h_RangeBias[3];

//Diagnostic histograms
TH1F *h_SameMerged[3];

//Number of momentum bins
int nbins = 10;

// Vectors to store all histograms and graphs
std::vector<TH2*> v_AllTH2;
std::vector<TH1*> v_AllTH1;
std::vector<TGraphAsymmErrors*> v_AllTGraph;

void calc_mcsres_Begin(TTree*){

  for (int part_i = 0; part_i < 3; part_i++){
    std::string particle = "null"; double max_mom = 0;
    if(part_i == 0){ particle = "muon"; max_mom = 1.5; }
    if(part_i == 1){ particle = "pion"; max_mom = 2.5; }
    if(part_i == 2){ particle = "proton"; max_mom = 2.5; }

    //Initialise histograms and graphs
    TString h_McsTrue_name = Form("h_McsTrue_%s",particle.c_str());
    h_McsTrue[part_i] = new TH2F(h_McsTrue_name,"",100,0,max_mom,100,0,max_mom);
    h_McsTrue[part_i]->GetXaxis()->SetTitle("True Momentum (GeV)");
    h_McsTrue[part_i]->GetYaxis()->SetTitle("MCS Momentum (GeV)");
    v_AllTH2.push_back(h_McsTrue[part_i]);

    TString h_RangeTrue_name = Form("h_RangeTrue_%s",particle.c_str());
    h_RangeTrue[part_i] = new TH2F(h_RangeTrue_name,"",100,0,max_mom,100,0,max_mom);
    h_RangeTrue[part_i]->GetXaxis()->SetTitle("True Momentum (GeV)");
    h_RangeTrue[part_i]->GetYaxis()->SetTitle("Range Momentum (GeV)");
    v_AllTH2.push_back(h_RangeTrue[part_i]);

    TString h_McsRange_name = Form("h_McsRange_%s",particle.c_str());
    h_McsRange[part_i] = new TH2F(h_McsRange_name,"",100,0,max_mom,100,0,max_mom);
    h_McsRange[part_i]->GetXaxis()->SetTitle("Range Momentum (GeV)");
    h_McsRange[part_i]->GetYaxis()->SetTitle("MCS Momentum (GeV)");
    v_AllTH2.push_back(h_McsRange[part_i]);

    TString h_McsTrueDiff_name = Form("h_McsTrueDiff_%s",particle.c_str());
    h_McsTrueDiff[part_i] = new TH2F(h_McsTrueDiff_name,"",100,0,max_mom,150,-1,1);
    h_McsTrueDiff[part_i]->GetXaxis()->SetTitle("True Momentum (GeV)");
    h_McsTrueDiff[part_i]->GetYaxis()->SetTitle("(P_{MCS}-P_{true})/P_{true}");
    v_AllTH2.push_back(h_McsTrueDiff[part_i]);

    TString h_RangeTrueDiff_name = Form("h_RangeTrueDiff_%s",particle.c_str());
    h_RangeTrueDiff[part_i] = new TH2F(h_RangeTrueDiff_name,"",100,0,max_mom,150,-1,1);
    h_RangeTrueDiff[part_i]->GetXaxis()->SetTitle("True Momentum (GeV)");
    h_RangeTrueDiff[part_i]->GetYaxis()->SetTitle("(P_{Range}-P_{true})/P_{true}");
    v_AllTH2.push_back(h_RangeTrueDiff[part_i]);

    TString h_McsResolution_name = Form("h_McsResolution_%s",particle.c_str());
    h_McsResolution[part_i] = new TH1F(h_McsResolution_name,"",10,0,1.5);
    h_McsResolution[part_i]->GetXaxis()->SetTitle("Momentum (GeV)");
    h_McsResolution[part_i]->GetYaxis()->SetTitle("MCS Fractional Resolution");
    v_AllTH1.push_back(h_McsResolution[part_i]);

    TString h_McsBias_name = Form("h_McsBias_%s",particle.c_str());
    h_McsBias[part_i] = new TH1F(h_McsBias_name,"",10,0,1.5);
    h_McsBias[part_i]->GetXaxis()->SetTitle("Momentum (GeV)");
    h_McsBias[part_i]->GetYaxis()->SetTitle("MCS Fractional Bias");
    v_AllTH1.push_back(h_McsBias[part_i]);

    TString h_RangeResolution_name = Form("h_RangeResolution_%s",particle.c_str());
    h_RangeResolution[part_i] = new TH1F(h_RangeResolution_name,"",10,0,1.5);
    h_RangeResolution[part_i]->GetXaxis()->SetTitle("Momentum (GeV)");
    h_RangeResolution[part_i]->GetYaxis()->SetTitle("Range Fractional Resolution");
    v_AllTH1.push_back(h_RangeResolution[part_i]);

    TString h_RangeBias_name = Form("h_RangeBias_%s",particle.c_str());
    h_RangeBias[part_i] = new TH1F(h_RangeBias_name,"",10,0,1.5);
    h_RangeBias[part_i]->GetXaxis()->SetTitle("Momentum (GeV)");
    h_RangeBias[part_i]->GetYaxis()->SetTitle("Range Fractional Bias");
    v_AllTH1.push_back(h_RangeBias[part_i]);

    TString h_SameMerged_name = Form("h_SameMerged_%s",particle.c_str());
    h_SameMerged[part_i] = new TH1F(h_SameMerged_name,"",10,0,10);
    h_SameMerged[part_i]->GetXaxis()->SetTitle("Geant Particles");
    h_SameMerged[part_i]->GetYaxis()->SetTitle("Events");
    v_AllTH1.push_back(h_SameMerged[part_i]);
  }

}

int calc_mcsres() {

  int nMerged = 0;

  for (int i = 0; i < geant_list_size; i++){
    if (MergedId[i]==MergedId[0]) nMerged++;
  }

  //Loop over tracks in event
  for (int track_i = 0; track_i < ntracks_pmalgtrackmaker; track_i++){

    int plane_i = trkpidbestplane_pmalgtrackmaker[track_i];

    //For each track, check that it is primary and that it stops (truth first, then PID module)
    int g4i = -1;
    for (int j = 0; j < geant_list_size; j++){
      if (trkidtruth_pmalgtrackmaker[track_i][plane_i]==TrackId[j]){ g4i =j; }
    }
    //Skip if track not in geant list
    if (g4i == -1 ) continue;
    //Skip if the particle is not the primary particle
    if (MergedId[g4i]!=MergedId[0]) continue;
  
/*
    //Skip if the particle is not the last geant particle corresponding to the primary
    bool isLast = 1;
    for (int j = 0; j < geant_list_size; j++){
      if (MergedId[g4i]==MergedId[j] && j>g4i) isLast = 0;
    }
    if (!isLast) continue;
*/
    //Get the particle index
    int part_i = -1;
    double momrange = 0;
    if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==13){ 
      part_i=0;
      momrange = trkmomrange_pmalgtrackmaker[track_i];
    }
    else if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==211){ 
      part_i=1;
      momrange = trkmomrange_pmalgtrackmaker[track_i];
    }
    else if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==2212){ 
      part_i=2;
      momrange = trkmomrangep_pmalgtrackmaker[track_i];
    }
    if (part_i==-1) continue;

    double chi2 = 0;
    int npt = 0;
    for (int hit_i = 0; hit_i < ntrkhits_pmalgtrackmaker[track_i][plane_i]; hit_i++){
      //Calculate the chi2 with a horizontal line at 3 MeV/cm
      if (trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]>1000) continue;
      if (trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]>25 || trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]<0.1) continue;
        // Needs to be done over the the same range as dedx_range_pro
      double err = 0.04231+0.0001783*trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]*trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i];
      err *= trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i];
      chi2 += pow((trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]-3)/std::sqrt(0.7+pow(err,2)),2); // Chi2 with horizontal line at 3+/-0.7 MeV/cm
      npt++;
    }
    if(chi2/npt < trkpidchi_pmalgtrackmaker[track_i][plane_i]) continue;

    //Get the true momentum, the MCS momentum and the range momentum
    //Fill 2D histograms
    h_McsTrue[part_i]->Fill(P[g4i], trkmcsmom_pmalgtrackmaker[track_i]);
    h_RangeTrue[part_i]->Fill(P[g4i], momrange);
    h_McsRange[part_i]->Fill(momrange, trkmcsmom_pmalgtrackmaker[track_i]);
    h_McsTrueDiff[part_i]->Fill(P[g4i], (trkmcsmom_pmalgtrackmaker[track_i]-P[g4i])/P[g4i]);
    h_RangeTrueDiff[part_i]->Fill(P[g4i], (momrange-P[g4i])/P[g4i]);

  }

  int ipart = -1;
  if (pdg[0]==13) ipart = 0;
  else if (pdg[0]==211) ipart = 1;
  else if (pdg[0]==2212) ipart = 2;
  if (ipart!=-1) h_SameMerged[ipart]->Fill(nMerged);

  return 0;
}

void calc_mcsres_Terminate(){

  for (int part_i = 0; part_i < 3; part_i++){
    for (int j = 1; j < nbins; j++){
      TString h_McsBin_name = Form("h_McsBin_%i_%i",part_i,j);
      TString h_RangeBin_name = Form("h_RangeBin_%i_%i",part_i,j);
      int firstbin = 100*((double)(j-1)/(double)nbins) + 1;
      int lastbin = 100*((double)j/(double)nbins);
      TH1D* h_McsBin = h_McsTrueDiff[part_i]->ProjectionY(h_McsBin_name,firstbin,lastbin);
      TH1D* h_RangeBin = h_RangeTrueDiff[part_i]->ProjectionY(h_RangeBin_name,firstbin,lastbin);

      if(h_McsBin->GetEntries()<100) continue;

      //h_MomBin->GetXaxis()->SetTitle("Momentum (GeV)");
      TF1 *mcsgaus = new TF1("mg","gaus",-0.6,0.6);
      h_McsBin->Fit(mcsgaus, "RQE");
      TF1 *rangegaus = new TF1("rg","gaus",-0.6,0.6);
      h_RangeBin->Fit(rangegaus,"RQE");
      v_AllTH1.push_back(h_McsBin);
      v_AllTH1.push_back(h_RangeBin);

      h_McsResolution[part_i]->SetBinContent(j,mcsgaus->GetParameter(2));
      h_McsResolution[part_i]->SetBinError(j,mcsgaus->GetParError(2));
      h_McsBias[part_i]->SetBinContent(j,mcsgaus->GetParameter(1));
      h_McsBias[part_i]->SetBinError(j,mcsgaus->GetParError(1));

      h_RangeResolution[part_i]->SetBinContent(j,rangegaus->GetParameter(2));
      h_RangeResolution[part_i]->SetBinError(j,rangegaus->GetParError(2));
      h_RangeBias[part_i]->SetBinContent(j,rangegaus->GetParameter(1));
      h_RangeBias[part_i]->SetBinError(j,rangegaus->GetParError(1));
    }
  }

  //Write all histograms and graphs to file
  for (unsigned int i = 0; i < v_AllTH2.size(); i++){
    v_AllTH2[i]->Write();
  }
  for (unsigned int i = 0; i < v_AllTH1.size(); i++){
    v_AllTH1[i]->Write();
  }
  for (unsigned int i = 0; i < v_AllTGraph.size(); i++){
    v_AllTGraph[i]->Write();
  }

}
