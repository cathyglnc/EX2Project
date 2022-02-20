

#include "./libPerso.h"



Double_t my_lin (Double_t *var, Double_t *par){
  // par[0]: slope
  // par[1]: offset
  
  return par[0]*var[0] + par[1];
}


Double_t my_gaus (Double_t *var, Double_t *par){
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2));
}


Double_t my_resp (Double_t *var, Double_t *par){

  Double_t *parLin = &par[0];
  Double_t *parGaus1 = &par[2];
  Double_t *parGaus2 = &par[5];
  
  return my_lin (var, parLin) + my_gaus (var, parGaus1) + my_gaus (var, parGaus2);
}


void selectTime (){

  std::cout << "\n Select a background spectrum by specifying the time window.\n" << std::endl;
  gStyle->SetOptStat (11111110);


  // set parameters
  const Int_t my_plug = 24;
  const Long64_t my_start = 12.*3600.e9; // h
  const Long64_t my_span = 1.*3600.e9; // 1 h
  
  
  // pull data
  Long64_t entry_g, stamp_g;
  Int_t plug_g, qdc_g;
  TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;

  TChain *f_in = new TChain ("t_entry");
  f_in->Add ("fixData.root"); // bkgr
  f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
  f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g); // 1..40
  f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g); // ns
  f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);
  

  // data containers
  const Double_t qdc_g_low = 0.;
  const Double_t qdc_g_high = 40.e3;
  const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/20;
    std::cout<<qdc_g_bin<<std::endl;
  TH1F *h_qdc_g = new TH1F ("h_qdc_g", "h_qdc_g", qdc_g_bin, qdc_g_low, qdc_g_high);


  // get time span
  f_in->GetEntry (0);
  std::cout << " start at " << std::setw (5) << (Int_t) (stamp_g/60.e9) << " min" << std::endl;
  f_in->GetEntry (f_in->GetEntries () - 1);
  std::cout << " stop at  " << std::setw (5) << (Int_t) (stamp_g/60.e9) << " min\n" << std::endl;  

  // data loop
  for (Long64_t i = 0; i < f_in->GetEntries ()/1.; i++)
    {     
      f_in->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0){
            //shows how many percents of the file have been read
            std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
        }
      
      if (plug_g == my_plug)
	{
	  if (stamp_g > my_start && stamp_g < my_start + my_span)
	    {
	      h_qdc_g->Fill (qdc_g);
	    }
	}
    }


  // data display
  TCanvas *c_qdc_g = new TCanvas ("c_qdc_g", "c_qdc_g", 700, 500);
  c_qdc_g->Draw ();
  c_qdc_g->cd ();

  h_qdc_g->GetXaxis ()->SetTitleOffset (0.95);
  h_qdc_g->GetYaxis ()->SetTitleOffset (0.95);
  h_qdc_g->GetXaxis ()->SetTitle ("QDC chn.");
  h_qdc_g->GetXaxis ()->SetRangeUser (18.e3, 22.e3);
  h_qdc_g->Draw ("");


  // fit
  TF1 *f_resp = new TF1 ("f_resp", my_resp, 18.e3, 21.5e3, 8);
  f_resp->SetParameter (0, -0.005); // slope
  f_resp->SetParameter (1, 110.); // offset
  f_resp->SetParameter (2, 260.); // ampl
  f_resp->SetParameter (3, 19600.); // mean
  f_resp->SetParameter (4, 190.); // sigma
  f_resp->SetParameter (5, 480.); // ampl
  f_resp->SetParameter (6, 20060.); // mean
  f_resp->SetParameter (7, 240.); // sigma

  for (Int_t i = 1; i < 8; i++) f_resp->SetParLimits (i, 0., 1.e6);
  f_resp->SetNpx (1e3);
  f_resp->SetLineColor (2);
  f_resp->SetLineWidth (2);
  
 
  h_qdc_g->Fit (f_resp, "RI");

  f_resp->Draw ("SAME");
  
  c_qdc_g->Print ("selectTime_c_qdc_g.pdf");
}
