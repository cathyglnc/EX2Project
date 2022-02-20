#include "./libPerso.h"


void energy_drift (){

  std::cout << "\n Display the energy drift of the background (R19_*.gz) from 2022-02-03.\n";
  gStyle->SetOptStat (1111110);


  // pull data
  Long64_t entry_g, stamp_g;
  Int_t plug_g, qdc_g;
  TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;

  TChain *f_in = new TChain ("t_entry");
  f_in->Add ("./Data/fixData.root");
  f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
  f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g);
  f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g);
  f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);
  
  // data containers
  const Double_t entry_g_low = 0.;
  const Double_t entry_g_high = 1250;
  const Int_t entry_g_bin = TMath::Abs (entry_g_high - entry_g_low)/10;
  const Double_t qdc_g_low = 0.;
  const Double_t qdc_g_high = 22500;
  const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/20;
  TH2F *h_qdcVsEntry_g = new TH2F ("h_qdcVsEntry_g", "h_qdcVsEntry_g", entry_g_bin, entry_g_low, entry_g_high, qdc_g_bin, qdc_g_low, qdc_g_high);
  

  // data loop
  for (Long64_t i = 0; i < f_in->GetEntries ()/1.; i++)
    {      
      f_in->GetEntry (i);
      if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;

      h_qdcVsEntry_g->Fill (stamp_g/60.e9, qdc_g);       
    }


  // data display
  TCanvas *c_qdcVsEntry_g = new TCanvas ("c_qdcVsEntry_g", "c_qdcVsEntry_g", 700, 500);
  c_qdcVsEntry_g->Draw ();
  c_qdcVsEntry_g->cd ();

  h_qdcVsEntry_g->GetYaxis ()->SetRangeUser (19.e3, 21.e3);
  h_qdcVsEntry_g->Draw ("COLZ");

  c_qdcVsEntry_g->Print ("fix_bkgrRun_c_qdcVsEntry_g.pdf");

  std::cout << "\n";
}
