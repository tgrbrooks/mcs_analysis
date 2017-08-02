//All dom brailsford's work, just playing with a few parameters
{
  TChain *chain = new TChain("analysistree/anatree");
  chain->Add("/sbnd/data/users/tbrooks/mcs_analysis/muons/muon_anatree_100-1257MeV.root");
  chain->Add("/sbnd/data/users/tbrooks/mcs_analysis/pions/pion_anatree_100-2250MeV.root");
  chain->Add("/sbnd/data/users/tbrooks/mcs_analysis/protons/proton_anatree_200-2250MeV.root");
  std::cout<<chain->GetEntries()<<std::endl;
  TFile *file = new TFile("output_mcs.root","RECREATE");
  chain->Draw("calc_mcsres.C+");
  file->Close();
}
