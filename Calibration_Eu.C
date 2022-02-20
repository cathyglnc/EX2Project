
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

void Calibration_Eu (int option){

    // choosing detector
    const Int_t det = 24;
    
    // defining Parameters Arrays (with default variables obtained with AskUser)
    
    //Detector 23
    //Double_t InfLim_params[7] = {4100,5100,5550,10100,12600,17200,18750};
    //Double_t SupLim_params[7] = {4700,5500,5950,10775,13400,18000,19750};
    
    //Detector 24
    Double_t InfLim_params[7] = {4300,5300,5800,10400,13000,17550,18750};
    Double_t SupLim_params[7] = {4900,5800,6300,11200,13800,18050,19750};
    Double_t amp_params[7]= {4400,750,780,815,620,63,460};
    Double_t mean_params[7]= {4375,5275,5750,10400,13000,17600,19189};
    Double_t std_params[7]= {200,200,200,200,200,200,220};
    
    // defining storage for fit results
    Double_t mean[7];
    Double_t error_mean[7];
  
  // pull data
    Long64_t entry_g, stamp_g;
    Int_t plug_g, qdc_g;
    TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;

    // data for the first calibration
    TChain *f_in = new TChain ("t_entry");
    f_in->Add (Form("./Data/eu152_2022020%i.root",option));
    cout<<Form("./Data/eu152_2022020%i.root",option)<<endl;
    f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
    f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g);
    f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g);
    f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);

    // data container
    const Double_t qdc_g_low = -1.e3;
    const Double_t qdc_g_high = 70.e3;
    const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/50;
    TH1F *h_qdc_g[40];
    for (Int_t i = 0; i < 40; i++) h_qdc_g[i] = new TH1F (Form ("h_qdc_g[%i]", i), Form ("qdc_chn %i", i + 1), qdc_g_bin, qdc_g_low, qdc_g_high);

    // taking data from first file
    for (Long64_t i = 0; i < f_in->GetEntries (); i++)
    {
        f_in->GetEntry (i);
        if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
      
        h_qdc_g[plug_g - 1]->Fill (qdc_g);
    }

    // drawing QDC spectrum "before"
    TCanvas *c_det = new TCanvas ("c_det", "", 700, 500);
    c_det->SetTitle("Nb of count with respect to QDC");
    c_det->Draw ();
    c_det->cd ()->SetTickx ();
    c_det->cd ()->SetTicky ();
    h_qdc_g[det]->GetXaxis ()->SetRangeUser (0., 30.e3);
    h_qdc_g[det]->Draw ("HIST");
    
    // loop on each fittable peaks that will be use in the calibration
    for(int i=0; i<7; i++){
        // possibility to ask user for peaks locations if first run
        // AskUser(i, amp_params, mean_params, std_params, InfLim_params, SupLim_params);
        
        // fitting each peaks with my_gaus and given parameters
        TF1 *f_fit = new TF1 ("f_spec", my_gaus, InfLim_params[i], SupLim_params[i], 3);
        f_fit->SetParameter (0, amp_params[i]);
        f_fit->SetParameter (1, mean_params[i]);
        f_fit->SetParameter (2, std_params[i]);
        f_fit->SetNpx (1e3);
        h_qdc_g[det]->Fit (f_fit, "RI");
        
        //storing results from the fit
        mean[i] = f_fit->GetParameter(1);
        error_mean[i]= f_fit->GetParError(1);
        
        //displaying fit results
        f_fit->Draw ("SAME");
        TPaveText *pt = new TPaveText(mean_params[i]-2,mean_params[i]+2,amp_params[i],amp_params[i]+2);
        char* mean = (char*)(&mean_params[i]);
        pt->AddText(mean);
        //pt->Draw();
    }
    
    //drawing Energy function of QDC
    auto c1 = new TCanvas("c1","Graph of mean with error bars",200,10,700,500);
    c1->SetGrid();
    c1->cd ()->SetTickx ();
    c1->cd ()->SetTicky ();

    const Int_t n = 7;
    
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
    gre3->GetXaxis()->SetTitle("mean");
    gre3->GetYaxis()->SetTitle("QDC");
    gre3->GetXaxis()->CenterTitle(true);
    gre3->GetYaxis()->CenterTitle(true);
    gre3->Draw("a*");
    
    //Fit the graph with the predefined "pol3" function
    gre3->Fit("pol1"); // polX : polynome de degré X
    
    //Access the fit resuts
    TF1 *f3 = gre3->GetFunction("pol1");
    f3->SetLineWidth(1);
    Double_t offset = f3->GetParameter("p0");
    Double_t ratio = f3->GetParameter("p1");
    
    //making legend on the E/QDC graph
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("Legend","C");
    legend->AddEntry("f3","Function fit","l");
    legend->AddEntry("gr","Graph with error bars","lep");
    legend->Draw();
    gStyle->SetOptStat (0);
    gStyle->SetOptFit (1);
    gPad->SetRightMargin (.03);
    gPad->SetTopMargin (.005);
    gPad->SetLeftMargin (.14);
    gPad->SetBottomMargin (.14);
    
    //defining calibrated histogram
    TH1F *h_calibrated[40];
    for (Int_t i = 0; i < 40; i++) h_calibrated[i] = new TH1F (Form ("h_calibrated[%i]", i), Form ("energy_chn %i", i + 1), qdc_g_bin, 0, 2000);
    
    //displaying first calibrated histogram
    auto c2 = new TCanvas("c2","Nb of counts with respect to the energy",200,10,700,500);
    c2->Draw ();
    c2->cd ()->SetTickx ();
    c2->cd ()->SetTicky ();

    cout << "Loading Scaled Histogram" << endl;
    for (Long64_t i = 0; i < f_in->GetEntries (); i++)
      {
          f_in->GetEntry (i);
          if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0) std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
        
          h_calibrated[plug_g - 1]->Fill (qdc_g*ratio + offset);
      }
    h_calibrated[det]->Draw ("HIST");
    
    //storing calibrated hist in root file
    TFile f(Form("./Data/Eu152_calibrated_2022020%i.root",option),"recreate");
    h_calibrated[det]->Write();
    TTree *tree = new TTree("tree","Ratio and Offsets Storage");
    tree->Branch("ratio",&ratio,Form("ratio_2022020%i/D",option));
    tree->Branch("offset",&offset,Form("offset_2022020%i/D",option));
    tree->Fill();
    tree->Write();
    
  std::cout << "\n";
}
