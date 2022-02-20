
#include "./libPerso.h"

Double_t my_gaus (Double_t *var, Double_t *par){
  // par[0]: amplitude
  // par[1]: mean
  // par[2]: sigma
  
  return par[0]*TMath::Exp (-.5*TMath::Power ((var[0] - par[1])/par[2], 2));
}

void AskUser(int Order, Double_t amp_params[7], Double_t mean_params[7], Double_t std_params[7],Double_t InfLim_params[7], Double_t SupLim_params[7]){
    printf( "-- Peak nb %d selection -- \n", Order+1);
    
    printf( "Specify Fitting Inferior Limit :");
    scanf("%lf", &InfLim_params[Order]);
    printf( "Specify Fitting Superior Limit :");
    scanf("%lf", &SupLim_params[Order]);
    
    printf( "Enter Amplitude :");
    scanf("%lf", &amp_params[Order]);
    printf("Enter Mean :");
    scanf("%lf", &mean_params[Order]);
    printf("Enter Standard Deviation :");
    scanf("%lf", &std_params[Order]);
}

// Temps entre calibration Avant et après 1277 min
// Ratio avant : 29.45 +- 0.08009
// Ratio après : 28.71 +- 0.07834

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
    tree2->GetEntry(0);
    cout<<"ratio1 : "<<ratio1<<" ratio2 : "<<ratio2<<endl;
    
    Double_t tempsTot = 1277;
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
    tree2->GetEntry(0);
    cout<<"offset1 : "<<offset1<<" offset2 : "<<offset2<<endl;
    
    Double_t tempsTot = 1277.;
    return (offset1+temps*(offset2-offset1)/tempsTot);
}

void calib_each_Time (){
    
    Double_t temps;
//    printf( "Enter time for interpolated calibration :");
//    scanf("%lf", &temps);
    
    Double_t ratio = InterpolationRatio(temps);
    Double_t offset = InterpolationOffset(temps);
    
    std::cout << "\n Simple script to loop through a standard format root file filling a histogram.\n";
    const Int_t det = 24;
    Double_t para_mean;
    
    //Définition des tableaux de paramètres
    //Detector 23
    //Double_t InfLim_params[7] = {4100,5100,5550,10100,12600,17200,18750};
    //Double_t SupLim_params[7] = {4700,5500,5950,10775,13400,18000,19750};
    
    //Detector 24
    Double_t InfLim_params[7] = {4300,5300,5800,10400,13000,17550,18750};
    Double_t SupLim_params[7] = {4900,5800,6300,11200,13800,18050,19750};
    Double_t amp_params[7] = {4400,750,780,815,620,63,460};
    Double_t mean_params[7] = {4375,5275,5750,10400,13000,17600,19189};
    Double_t std_params[7] = {200,200,200,200,200,200,220};
    
    Double_t error_mean[7];
  
    // pull data
    Long64_t entry_g, stamp_g;
    Int_t plug_g, qdc_g;
    TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;

    TChain *f_in = new TChain ("t_entry");
    f_in->Add ("./Data/eu152_20220203.root");
    f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
    f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g);
    f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g);
    f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);


    // data container
    const Double_t qdc_g_low = -1.e3;
    const Double_t qdc_g_high = 70.e3;
    const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/20;
    TH1F *h_qdc_g[40];
    for (Int_t i = 0; i < 40; i++) h_qdc_g[i] = new TH1F (Form ("h_qdc_g[%i]", i), Form ("qdc_chn %i", i + 1), qdc_g_bin, qdc_g_low, qdc_g_high);

    // data loop
    for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
        f_in->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;

            //if (i < 10) std::cout << i << " " << plug_g << " " << stamp_g << " " << qdc_g << std::endl;
      
        h_qdc_g[plug_g - 1]->Fill (qdc_g);
    }
    
//    // fit
//    TF1 *f_fit = new TF1 ("f_spec", my_gaus, 18700, 19600, 3);
//    f_fit->SetParameter (0, 450.);
//    f_fit->SetParameter (1, 19100.);
//    f_fit->SetParameter (2, 220.);
//    // f_fit->SetParameters (0.1, 400., 4500., 5., 0.5);
//    f_fit->SetNpx (1e3);
//    h_qdc_g[det]->Fit (f_fit, "RI");

//  // data display
//  TCanvas *c_qdc = new TCanvas ("c_qdc", "", 700, 500);
//  c_qdc->Draw ();
//  c_qdc->Divide (8, 5);
//  for (Int_t i = 0; i < 40; i++)
//    {
//      c_qdc->cd (i + 1)->SetTickx ();
//      c_qdc->cd (i + 1)->SetTicky ();
//      c_qdc->cd (i + 1)->SetLogy ();
//
//      h_qdc_g[i]->GetXaxis ()->SetRangeUser (0., 30.e3);
//      h_qdc_g[i]->Draw ("HIST");
//    }
//  c_qdc->Print ("loop_rootFile_c_qdc.pdf");

    // single hist display
    TCanvas *c_det = new TCanvas ("c_det", "", 700, 500);
    c_det->Draw ();
    c_det->cd ()->SetTickx ();
    c_det->cd ()->SetTicky ();

    h_qdc_g[det]->GetXaxis ()->SetRangeUser (0., 30.e3);
  
    h_qdc_g[det]->Draw ("HIST");
    Double_t mean[7];
    for(int i=0; i<8; i++){
            //AskUser(i, amp_params, mean_params, std_params, InfLim_params, SupLim_params);
        
            // fit
            TF1 *f_fit = new TF1 ("f_spec", my_gaus, InfLim_params[i], SupLim_params[i], 3);
            f_fit->SetParameter (0, amp_params[i]);
            f_fit->SetParameter (1, mean_params[i]);
            f_fit->SetParameter (2, std_params[i]);
            // f_fit->SetParameters (0.1, 400., 4500., 5., 0.5);
            f_fit->SetNpx (1e3);
            h_qdc_g[det]->Fit (f_fit, "RI");
        
            mean[i] = f_fit->GetParameter(1);
        
            //Erreur
            error_mean[i]= f_fit->GetParError(1);
        
            //Affichage
            f_fit->Draw ("SAME");
            TPaveText *pt = new TPaveText(mean_params[i]-2,mean_params[i]+2,amp_params[i],amp_params[i]+2);
            char* mean = (char*)(&mean_params[i]);
            pt->AddText(mean);
    }
    
    auto c1 = new TCanvas("c1","Graph of mean with error bars",200,10,700,500);
    c1->SetFillColor(42);
    c1->SetGrid();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    
    const Int_t n = 8;
    
    Double_t energy[7]={344.2785,411.1163,443.965,778.9040,964.064,1299.140,1408.006};
    Double_t error_energy[7]={0.0012,0.0011,0.003,0.0018,0.018,0.010,0.003};
    
    for(int j=0; j<7; j++){
            printf("Moyenne : ");
        std::cout << mean[j];
        printf("- Energy : ");
        std::cout << energy[j];
        printf("\n");
    }
    TGraphErrors *gre3 = new TGraphErrors(n, mean, energy, error_mean, error_energy);
    gre3->SetTitle("QDC evolution with respect to mean");
    gre3->Draw("a*");
    //Fit the graph with the predefined "pol3" function
    gre3->Fit("pol1"); // polX : polynome de degré X
    //Access the fit resuts
    TF1 *f3 = gre3->GetFunction("pol1");
    f3->SetLineWidth(1);
    //Double_t offset = f3->GetParameter(0);
    //Double_t ratio = f3->GetParameter("p0");
    cout << "ratio " << ratio << endl;
    cout << "offset " << offset << endl;
    
    
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Legend","C");
    legend->AddEntry("f3","Function fit","l");
    legend->AddEntry("gr","Graph with error bars","lep");
    legend->Draw();
    gStyle->SetOptStat (0);
    gStyle->SetOptFit (1);
    
    //Définition de la fonction de conversion
    TH1F *h_calibrated[40];
    for (Int_t i = 0; i < 40; i++) h_calibrated[i] = new TH1F (Form ("h_calibrated[%i]", i), Form ("energy_chn %i", i + 1), qdc_g_bin, 0, 2000);
    
    auto c2 = new TCanvas("c2","Nb of counts with respect to the energy",200,10,700,500);
    c2->Draw ();
    c2->SetTickx ();
    c2->cd ()->SetTicky ();
    //h_calibrated[det]->GetXaxis ()->SetRangeUser (0., 30.e3);
    
    cout << "Loading Scaled Histogram" << endl;
    for (Long64_t i = 0; i < f_in->GetEntries (); i++)
      {
        
        f_in->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
        
        h_calibrated[plug_g - 1]->Fill (qdc_g*ratio + offset);
      }
    
    h_calibrated[det]->Draw ("HIST");
    
//    TFile f("Eu152_calibrated_20220203.root","new");
//    h_calibrated[det]->Write();
    // Tracé de spectre en QDC
    c_det->Print ("loop_rootFile_c_det.pdf");
    
    // Tracé du fit pour conversion QDC <=> Energies
    c1->Print ("TGraph_Error.pdf");
    
  std::cout << "\n";
}
