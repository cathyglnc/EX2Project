#include "./libPerso.h"
#include "./fit_selfRadio.h"

#define duration 21


Double_t InterpolationRatio(Double_t temps){
    TFile f1("./Data/Eu152_calibrated_20220203.root");
    TTree *tree1 = f1.Get<TTree>("tree");
    Double_t ratio1;
    tree1->SetBranchAddress("ratio", &ratio1);
    
    TFile f2("./Data/Eu152_calibrated_20220204.root");
    TTree *tree2 = f2.Get<TTree>("tree");
    Double_t ratio2;
    tree2->SetBranchAddress("ratio", &ratio2);
    
    tree1->GetEntry(0);
    f1.Close();
    tree2->GetEntry(0);
    f2.Close();
    
    Double_t tempsTot = 1277.;
    return (ratio1+temps*(ratio2-ratio1)/tempsTot);
}

Double_t InterpolationOffset(Double_t temps){
    TFile f1("./Data/Eu152_calibrated_20220203.root");
    TTree *tree1 = f1.Get<TTree>("tree");
    Double_t offset1;
    tree1->SetBranchAddress("offset", &offset1);
    
    TFile f2("./Data/Eu152_calibrated_20220204.root");
    TTree *tree2 = f2.Get<TTree>("tree");
    Double_t offset2;
    tree2->SetBranchAddress("offset", &offset2);
    
    tree1->GetEntry(0);
    f1.Close();
    tree2->GetEntry(0);
    f2.Close();
    
    Double_t tempsTot = 1277.;
    return ((offset2+offset1)/2+temps*(offset2-offset1)/tempsTot);
}

void hist_superposition(){
    
    TColor *c3 = new TColor(9001,0,0,0); //black
    TColor *c4 = new TColor(9000,1,0,0); // red
    TColor *c7 = new TColor(9002,0,0,1); // blue
    
    TFile *f1 = new TFile("./Data/Eu152_calibrated_20220203.root");
    TH1 *h1 ; f1->GetObject("h_calibrated[24];1",h1);
    TFile *f2 = new TFile("./Data/Eu152_calibrated_20220204.root");
    TH1 *h2; f2->GetObject("h_calibrated[24];1",h2);
    
    // pull data
    Long64_t entry_g, stamp_g;
    Int_t plug_g, qdc_g;
    TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;
    Int_t det= 24;
    
    const Double_t qdc_g_low = -1.e3;
    const Double_t qdc_g_high = 70.e3;
    const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/50;
    
    TChain *f_in = new TChain ("t_entry");
    f_in->Add ("./Data/eu152_20220203.root");
    f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
    f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g);
    f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g);
    f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);
    
    
    // h_interpolated ------------------------------------------------------------
    Double_t ratio_inter;
    Double_t offset_inter;
    TH1F *h_interpolated = new TH1F ("h_interpolated", Form ("Number of count with respect to the Energy (detector %i)", det), qdc_g_bin, 0, 2000);
    
    //for(Int_t j=0; j<duration; j++){
        ratio_inter=InterpolationRatio(12);
        offset_inter=InterpolationOffset(12);
        for (Long64_t i = 0; i < f_in->GetEntries (); i++)
        {
            
            f_in->GetEntry (i);
            
            if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
            if (plug_g-1 == det){
                h_interpolated->Fill (qdc_g*ratio_inter + offset_inter);
            }
        }
    //}
    
    // h_calibrated ------------------------------------------------------------
    
//    Double_t offset_calib = 9.319;//04
//    Double_t ratio_calib = 0.07227;
    
    Double_t offset_calib = 9.714;//03
    Double_t ratio_calib = 0.07204;
    TH1F *h_calibrated = new TH1F ("h_calibrated", Form ("Number of count with respect to the Energy (detector %i)", det), qdc_g_bin, 0, 2000);
    
    //for(Int_t j=0; j<duration; j++){
        for (Long64_t i = 0; i < f_in->GetEntries (); i++)
        {
            
            f_in->GetEntry (i);
            
            if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
            if (plug_g-1 == det){
                h_calibrated->Fill (qdc_g*ratio_calib + offset_calib);
            }
        }
    //}
    
    
    // Simulation ---------------------------------------------------------------
    Double_t ratio[duration];
    Double_t offset[duration];
    
    TH1F *h_simu = new TH1F ("h_simu", Form ("Number of count with respect to the Energy (detector %i)", det), qdc_g_bin, 0, 2000);
    
    TFile f3("./Data/Simulation_factors.root");
    TTree *tree = f3.Get<TTree>("tree");
    tree->SetBranchAddress("ratio", ratio);
    tree->SetBranchAddress("offset", offset);
    
    tree->GetEntry(0);

    //for(Int_t j=0; j<duration; j++){
        for (Long64_t i = 0; i < f_in->GetEntries (); i++)
        {
        f_in->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
        if (plug_g-1 == det){
                h_simu->Fill (qdc_g*ratio[12] + offset[12]);
            }
        }
    //}
    
    
        

    f3.Close();
    
    TCanvas *c6 = new TCanvas("c6","Histogram superposition (calibrated on 20220203)",200,10,700,500);
    c6->Draw ();
    c6->SetTickx ();
    c6->cd ()->SetTicky ();
    h_interpolated->SetLineColor(9002);
    h_interpolated->Draw("HIST");
    h_simu->SetLineColor(9000);
    h_simu->Draw("SAME");
    h_calibrated->SetLineColor(9001);
    h_calibrated->Draw ("SAME");
    
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
       //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
       //legend->AddEntry(h1,"calibration 20220203","l");
       legend->AddEntry(h_interpolated,"interpolation","l");
       legend->AddEntry(h_simu,"simulation","l");
       legend->AddEntry(h_calibrated,"calibration 20220203","l");
       legend->Draw();

    TFile f("./Data/Superposition_calibrated.root","recreate");
    c6->Write();
}
