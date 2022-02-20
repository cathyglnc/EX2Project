
#include "fit_selfRadio.h"

#define duration 21


Double_t my_reso (Double_t *var, Double_t *par){
    // relative energy resolution of the scitillators: de/e = a + b/sqrt (e)
    // de = sigma = a*e + b*sqrt (e)
    
    return par[0]*var[0] + par[1]*TMath::Sqrt (var[0]);
}


Double_t my_exp (Double_t *var, Double_t *par){
    
    // par[0]: amplitude
    // par[1]: bin
    
    return par[0]*TMath::Exp (-par[1]*var[0]);
}


Int_t fill_sim (const Int_t det_sim, const Int_t scale_sim){
    // fill the histogram from the simulation
    std::cout << " Filling spectrum from simulation." << std::endl;
    
    
    // pull data
    TChain *f_sim = new TChain ("v3");
    f_sim->Add ("./Data/la138_array_fullShell.root");
    f_sim->SetBranchAddress ("event", &event, &b_event);
    f_sim->SetBranchAddress ("ene_g", &ene_g, &b_ene_g);
    f_sim->SetBranchAddress ("det_g", &det_g, &b_det_g);
    f_sim->SetBranchAddress ("id", &id, &b_id);
    f_sim->SetBranchAddress ("pos_x", &pos_x, &b_pos_x);
    f_sim->SetBranchAddress ("pos_y", &pos_y, &b_pos_y);
    f_sim->SetBranchAddress ("pos_z", &pos_z, &b_pos_z);
    f_sim->SetBranchAddress ("mom_x", &mom_x, &b_mom_x);
    f_sim->SetBranchAddress ("mom_y", &mom_y, &b_mom_y);
    f_sim->SetBranchAddress ("mom_z", &mom_z, &b_mom_z);
    
    for (Int_t i = 0; i < 2; i++) h_sim_g[i] = new TH1F (Form ("h_sim_g[%i]", i), Form ("h_sim_g[%i]", i), ene_g_bin, ene_g_low, ene_g_high);
    h_sim_g[1]->SetLineColor (2);
    
    TF1 *f_reso = new TF1 ("f_reso", my_reso, 0., 700., 2); // LaBr3 resolution: you can extract from your 152Eu run, or take what I have here
    f_reso->SetParameters (-0.0104, 0.0233);
    Double_t fac_res = 1.; // make the resolution a little worse for higher energies; i.e. 0.79MeV-100%, 1.47MeV-130%

    // data loop
    for (Long64_t i = 0; i < f_sim->GetEntries ()/scale_sim; i++)
    {
        f_sim->GetEntry (i);
        if (i%((Long64_t) f_sim->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_sim->GetEntries () + .5) << " % done" << std::endl;
        
        
        for (Int_t j = 0; j < (Int_t) ene_g->size (); j++)
        {
            if (det_g->at (j) != det_sim) continue;
            if (ene_g->at (j) < 0.1) continue; // Suppress first peek
            
//             if (ene_g->at (j) < 0.79) fac_res = 1.; // worsen the energy resolution just for higher energies; need to see if you want to use
//             else fac_res = (1.35 - 1.0)/(1.470 - 0.790)*(ene_g->at (j) - 0.790) + 1.0;
            
            h_sim_g[0]->Fill (ene_g->at (j));
            h_sim_g[1]->Fill (gRandom->Gaus (ene_g->at (j), fac_res*f_reso->Eval (ene_g->at (j))));
        }
    }
    
    return 0;
}


Int_t fill_exp (const Int_t qdc_exp, const Long64_t start_exp, const Long64_t span_exp){
    // fill the histogram from the experiment
    std::cout << " Filling spectrum from experiment." << std::endl;
    
    // pull data
    f_exp = new TChain ("t_entry");
    f_exp->Add ("./Data/fixData.root"); // bkgr run
    f_exp->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
    f_exp->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g); // 1..40
    f_exp->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g); // ns
    f_exp->SetBranchAddress ("qdc_g"  , &qdc_g  , &b_qdc_g);
    
    for (Int_t i = 0; i < 2; i++)
    {
        h_qdc_g[i] = new TH1F (Form ("h_qdc_g[%i]", i), Form ("h_qdc_g[%i]", i), qdc_g_bin, qdc_g_low, qdc_g_high);
        h_ene_g[i] = new TH1F (Form ("h_ene_g[%i]", i), Form ("h_ene_g[%i]", i), ene_g_bin, ene_g_low, ene_g_high); // I define both exp. energy histograms already here, but fill only the first one (crude)
    }
    
    
    // data loop
    for (Long64_t i = 0; i < f_exp->GetEntries ()/1.; i++)
    {
        f_exp->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_exp->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_exp->GetEntries () + .5) << " % done" << std::endl;
        
        if (plug_g == qdc_exp)
        {
            if (stamp_g > start_exp && stamp_g < start_exp + span_exp)
            {
                for (Int_t j = 0; j < 2; j++) h_qdc_g[j]->Fill (qdc_g);
                h_ene_g[0]->Fill (slope*qdc_g + offset);
            }
        }
    }
    
    return 0;
}


Double_t my_hist_fit (Double_t *var, Double_t *par){
    // fit the bins of the experiment (linear) as well as the amplitude to the simulation
    // takes into account as well a general exponential
    // par[0]: exp. bin offset
    // par[1]: exp. bin scaling
    // par[2]: hist scaling
    // par[3]: expo. amplitude
    // par[4]: expo. bin scaling
    
    Double_t bin_exp = h_qdc_g[1]->GetXaxis ()->FindBin ((Int_t) ((var[0] - par[0])/par[1]));
    Double_t *setExpo = &par[3];
    
    //  return par[2]*h_qdc_g[1]->GetBinContent ((Int_t) bin_exp); // no exponential
    return par[2]*h_qdc_g[1]->GetBinContent ((Int_t) bin_exp) - my_exp (var, setExpo);
}


void fit_selfRadio (const Long64_t my_start, Double_t *my_ratio, Double_t *my_offset){

    std::cout << "\n";
    gStyle->SetOptStat (1111110);
    const Int_t my_det = 19; // detector position in LaBr3 array
    const Int_t my_qdc = 24; // QDC chn. from digitizer modules
    const Int_t my_scale = 1; // might want to scale statistics of simulation to experiment
    const Long64_t my_span = 1.*3600.e9; // 1 h
    const Double_t fit_low = 0.6; // MeV
    const Double_t fit_high = 1.7;

    // crude calibration (checking the h_qdc_g in c_hist)
    peak1 = 10600;
    peak2 = 20100;
    slope = (1.46 - 0.79)/(peak2 - peak1); //0.0705263e-3
    offset = 1.46 - slope*peak2; //42.42105e-3

    
    //instedad filling each time we put it as a parameter
    //fill_sim (my_det, my_scale);
    fill_exp (my_qdc, my_start, my_span);


//    // display histograms
//    TCanvas *c_hist = new TCanvas ("c_hist", "all histograms involved", 700, 500);
//    c_hist->Draw ();
//    c_hist->Divide (2, 2);
//    for (Int_t i = 0; i < 4; i++)
//      {
//        c_hist->cd (i + 1)->SetTickx ();
//        c_hist->cd (i + 1)->SetTicky ();
//
//        if (i == 0) h_sim_g[0]->Draw ("HIST"); // bare simulation
//        else if (i == 1) h_sim_g[1]->Draw ("HIST"); // simulation with detector resolution
//        else if (i == 2) h_qdc_g[0]->Draw ("HIST");
//        else h_ene_g[0]->Draw ("HIST");
//      }
//
//    c_hist->Print ("fit_selfRadio_c_hist.pdf");


    // fit simulation with fit function (= experimental data - exponential)
    // TF1 *f_fit[i] = new TF1 (Form ("f_fit[%i]", i), my_hist_fit, fit_low, fit_high, 3); // no exponential
    TF1 *f_fit = new TF1 ("f_fit", my_hist_fit, fit_low, fit_high, 5); // see parameters in fit function
    f_fit->SetLineWidth (2);
    f_fit->SetLineColor (3);
    f_fit->SetNpx ((fit_high - fit_low)*ene_g_bin/(ene_g_high - ene_g_low)); // binning from the simulations
    f_fit->SetParameters (offset, slope, 1., 5.e3, 5.);
    

    TCanvas *c_fit = new TCanvas ("c_fit", "the actual fit", 700, 500);
    c_fit->Draw ();
    c_fit->cd ()->SetTickx ();
    c_fit->cd ()->SetTicky ();
    
    h_sim_g[1]->Fit (f_fit, "RI0");
    h_sim_g[1]->Draw ("HIST");
    f_fit->Draw ("SAME");
    
    *my_ratio = f_fit->GetParameter (1);
    *my_offset = f_fit->GetParameter (0);

    TF1 *f_expo = new TF1 ("f_expo", my_exp, fit_low, fit_high, 2); // see parameters in fit function
    f_expo->SetLineWidth (2);
    f_expo->SetLineColor (4);
    f_expo->SetLineStyle (2);
    f_expo->SetNpx (1e3);
    f_expo->SetParameters (f_fit->GetParameter (3), f_fit->GetParameter (4));
     f_expo->Draw ("SAME");

    c_fit->Print ("fit_selfRadio_c_fit.pdf");


    // calibrate experimental data with fit parameters
    // data loop
    for (Long64_t i = 0; i < f_exp->GetEntries ()/1.; i++)
      {
        f_exp->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_exp->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_exp->GetEntries () + .5) << " % done" << std::endl;
        
        if (plug_g == my_qdc)
      {
        if (stamp_g > my_start && stamp_g < my_start + my_span) h_ene_g[1]->Fill (*my_ratio * qdc_g + *my_offset);
      }
      }

    TCanvas *c_calib = new TCanvas ("c_calib", "calibrated histogram", 700, 500);
    c_calib->Draw ();
    c_calib->cd ()->SetTickx ();
    c_calib->cd ()->SetTicky ();

    h_ene_g[1]->Draw ("HIST");
    f_fit->Draw ("SAME");
    f_expo->Draw ("SAME");

    c_calib->Print ("fit_selfRadio_c_calib.pdf");
    
    *my_ratio = *my_ratio*1000;
    *my_offset = *my_offset*1000;
    
    std::cout << "\n";
}
void GetAllMeans_21(){
    
    Double_t ratio[duration];
    Double_t offset[duration];
    
    TFile f("./Data/Simulation_factors.root","recreate");
    TTree *tree = new TTree("tree","SimuFact");
    tree->Branch("ratio",ratio,Form("ratio[%i]/D",duration));
    tree->Branch("offset",offset,Form("offset[%i]/D",duration));
    
    Long64_t my_start;
    
    fill_sim (19, 1.);
    
    for (Int_t i=0; i<duration; i++)
    {
        my_start = i*3600.e9; // h
        fit_selfRadio (my_start, &ratio[i], &offset[i]);
        cout<<"step : "<< i <<" ratio : "<<ratio[i]<<" offset : "<<offset[i]<<endl;
    }
    
    tree->Fill();
    tree->Write();
}

