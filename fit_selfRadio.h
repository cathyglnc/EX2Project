
#include "./libPerso.h"

// experimental data
TChain *f_exp;
Long64_t entry_g, stamp_g;
Int_t plug_g, qdc_g;
TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;

// simulation
Long64_t event;
std::vector <Int_t> *id = 0;
std::vector <Int_t> *det_g = 0;
std::vector <Double_t> *ene_g = 0;
std::vector <Double_t> *pos_x = 0;
std::vector <Double_t> *pos_y = 0;
std::vector <Double_t> *pos_z = 0;
std::vector <Double_t> *mom_x = 0;
std::vector <Double_t> *mom_y = 0;
std::vector <Double_t> *mom_z = 0;
TBranch *b_event, *b_id, *b_det_g, *b_ene_g, *b_pos_x, *b_pos_y, *b_pos_z, *b_mom_x, *b_mom_y, *b_mom_z;

const Double_t qdc_g_low = 0.; // QDC chn.
const Double_t qdc_g_high = 70000.;
const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/100;
const Double_t ene_g_low = 0.; // MeV
const Double_t ene_g_high = 2.;
const Int_t ene_g_bin = TMath::Abs (ene_g_high - ene_g_low)*100;

TH1F *h_qdc_g[2];
TH1F *h_sim_g[2]; // [raw, resolution]
TH1F *h_ene_g[2]; // [crude, fitted]

Int_t peak1, peak2;
Double_t slope, offset;
