

#include "./libPerso.h"

#define duration 21
#define start 0


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

void FillHistogram(TH1F *hist_hour[duration], Int_t my_plug, Int_t ratio_ref){
    //if ratio_ref is 0 then we take default ratio and offset
    //if ratio_ref is 1 then we take ratio and offset from calibration of 20020203
    //if ratio_ref is 2 then we take ratio and offset from Interpolation
    //if ratio_ref is 3 then we take ratio and offset from Simu
    
    Double_t ratio[duration];
    Double_t offset[duration];
    if(ratio_ref==1){
        TFile f("./Data/Eu152_calibrated_20220203.root");
        TTree *tree = f.Get<TTree>("tree");
        tree->SetBranchAddress("ratio", ratio);
        tree->SetBranchAddress("offset", offset);
        tree->GetEntry(0);
        f.Close();
    }
    else if(ratio_ref == 3){
        TFile f("./Data/Simulation_factors.root");
        TTree *tree = f.Get<TTree>("tree");
        tree->SetBranchAddress("ratio", ratio);
        tree->SetBranchAddress("offset", offset);
        for(Int_t i=0; i<duration; i++){
            tree->GetEntry(i);
            cout<<"step : "<<i<<" ratio : "<<ratio[i]<<" offset : "<<offset[i]<<endl;
        }
        f.Close();
        
    }
    else{
        ratio[0] = 1;
        offset[0] = 0;
    }

    const Long64_t my_span = 1.*3600.e9; // 1 h
    
    // pull data
    Long64_t entry_g, stamp_g;
    Int_t plug_g, qdc_g;
    TBranch *b_entry_g, *b_plug_g, *b_stamp_g, *b_qdc_g;
    
    TChain *f_in = new TChain ("t_entry");
    f_in->Add ("./Data/fixData.root"); // bkgr
    f_in->SetBranchAddress ("entry_g", &entry_g, &b_entry_g);
    f_in->SetBranchAddress ("plug_g"  , &plug_g  , &b_plug_g); // 1..40
    f_in->SetBranchAddress ("stamp_g", &stamp_g, &b_stamp_g); // ns
    f_in->SetBranchAddress ("qdc_g" , &qdc_g , &b_qdc_g);
    
    // data containers
    const Double_t qdc_g_low = 0.;
    const Double_t qdc_g_high = 40.e3;
    
    const Int_t qdc_g_bin = TMath::Abs (qdc_g_high - qdc_g_low)/20;
    //TH1F *hist_hour = new TH1F ("hist_hour", "hist_hour", qdc_g_bin, qdc_g_low, qdc_g_high);
    
    
    for (Int_t i = 0; i < duration; i++) hist_hour[i] = new TH1F (Form ("hist_hour[%i]", i), Form ("qdc_chn %i", i + 1), 2000, 0, 3000);
    
    // get time span
    f_in->GetEntry (0);
    std::cout << " start at " << std::setw (5) << (Int_t) (stamp_g/60.e9) << " min" << std::endl;
    f_in->GetEntry (f_in->GetEntries () - 1);
    std::cout << " stop at  " << std::setw (5) << (Int_t) (stamp_g/60.e9) << " min\n" << std::endl;
    
    Long64_t time;
    
    //hour loop
    for(Long64_t j=0; j<duration; j++){
        
        //defining the time stamp for this step
        time = j*3600.e9;
        
        if(ratio_ref ==2){
            ratio[0] = InterpolationRatio(j*60);
            offset[0] = InterpolationOffset(j*60);
            cout<<"ratio : "<<ratio[0]<<" offset : "<<offset[0]<<endl;
        }
        // data loop
        for (Long64_t i = 0; i < f_in->GetEntries ()/1.; i++)
        {
            f_in->GetEntry (i);
            
            //shows progress in file reading
            if (i > 11 && i%((Long64_t) f_in->GetEntries ()/10) == 0){
                std::cout << " " << std::setw (3) << (Int_t) ((i*100.)/f_in->GetEntries () + .5) << " % done" << std::endl;
            }
            
            //for the selected detector
            if (plug_g == my_plug)
            {
                //selection per hour
                if (stamp_g > time && stamp_g < time + my_span)
                {
                    if(ratio_ref == 3){
                        hist_hour[j]->Fill(qdc_g*ratio[j]+offset[j]);
                    }
                    else{
                        hist_hour[j]->Fill(qdc_g*ratio[0]+offset[0]);
                    }
                    
                }
            }
        }
    }
    //delete hist_hour[duration];
}

void SetParametersQDC(TF1 *fct_to_parametrize){
    fct_to_parametrize->SetParameter (0, -0.005); // slope
    fct_to_parametrize->SetParameter (1, 110.); // offset
    fct_to_parametrize->SetParameter (2, 260.); // ampl
    fct_to_parametrize->SetParameter (3, 19600.); // mean
    fct_to_parametrize->SetParameter (4, 190.); // sigma
    fct_to_parametrize->SetParameter (5, 480.); // ampl
    fct_to_parametrize->SetParameter (6, 20060.); // mean
    fct_to_parametrize->SetParameter (7, 240.); // sigma
}

void SetParametersENERGY(TF1 *fct_to_parametrize){
    fct_to_parametrize->SetParameter (0, -0.005); // slope
    fct_to_parametrize->SetParameter (1, 110.); // offset
    fct_to_parametrize->SetParameter (2, 400.); // ampl
    fct_to_parametrize->SetParameter (3, 1430.); // mean
    fct_to_parametrize->SetParameter (4, 30.); // sigma
    fct_to_parametrize->SetParameter (5, 500.); // ampl
    fct_to_parametrize->SetParameter (6, 1470.); // mean
    fct_to_parametrize->SetParameter (7, 30.); // sigma
}

void FittingStep(TH1F *hist_hour[duration], TF1 *f_resp, Double_t *mean1, Double_t *error_mean1, Double_t *mean2, Double_t *error_mean2, Double_t option){
    
    if(option==0){
        f_resp = new TF1 ("f_resp", my_resp, 18.e3, 21.5e3, 8);
        //set fit param in function of the option seleted
        SetParametersQDC(f_resp);
    }
    else{
        f_resp= new TF1 ("f_resp", my_resp, 1350, 1550, 8);
        //set fit param in function of the option seleted
        SetParametersENERGY(f_resp);
    }
    
    for (Int_t i = 1; i < 8; i++) f_resp->SetParLimits (i, 0., 1.e6);
    f_resp->SetNpx (1e3);
    f_resp->SetLineColor (2);
    f_resp->SetLineWidth (2);
    
    for (Int_t j=0;j<duration;j++)
    {
        //fitting the main peek at each hour with f_resp
        hist_hour[j]->Fit (f_resp, "RI0");
        
        mean1[j] = f_resp->GetParameter(3);
        error_mean1[j]=f_resp->GetParError(3);
        mean2[j] = f_resp->GetParameter(6);
        error_mean2[j]=f_resp->GetParError(6);
    }
}

void Cumulated_Histogram(TH1F *hist_hour[duration],TH1F *hist_cumulated){
    
    for (Int_t i=0; i<duration;i++){
        hist_cumulated->Add(hist_hour[i]);
    }
}

void FittingCumulated(TH1F *hist_cumulated, TF1 *f_resp, Double_t *mean1, Double_t *error_mean1, Double_t *mean2, Double_t *error_mean2){
    
    f_resp= new TF1 ("f_resp", my_resp, 1350, 1550, 8);
    //set fit param in function of the option seleted
    SetParametersENERGY(f_resp);
    
    for (Int_t i = 1; i < 8; i++) f_resp->SetParLimits (i, 0., 1.e6);
    f_resp->SetNpx (1e3);
    f_resp->SetLineColor (2);
    f_resp->SetLineWidth (2);
    
    hist_cumulated->Fit (f_resp, "RI0");
    
    *mean1 = f_resp->GetParameter(3);
    *error_mean1 = f_resp->GetParError(3);
    *mean2 = f_resp->GetParameter(6);
    *error_mean2 = f_resp->GetParError(6);
    
}

void meanEvolTime_21 (Int_t option){
    //different behaviours for different values of option:
    //0 each peek is fitted but nothing is calibrated
    //1 peeks are fitted and only one calibration (20220203) is used
    //2 peeks are fitted and calibration is done via Intrepolation of 2-times calibration
    //3 peeks are fitted on simulation calibration
    //4 : 1 and 2 together
    //5 : 1 and 3 together
    //6 : all together
    
    std::cout << "\n Select a background spectrum by specifying the time window.\n" << std::endl;
    gStyle->SetOptStat (11111110);
    
    // set parameters
    const Int_t my_plug = 24;
    Double_t mean1[duration];
    Double_t error_mean1[duration];
    Double_t mean2[duration];
    Double_t error_mean2[duration];
    
    Double_t mean_cumulated1;
    Double_t error_mean_cumulated1;
    Double_t mean_cumulated2;
    Double_t error_mean_cumulated2;
    
    int n = duration;
    Double_t hours1[21]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
    Double_t hours2[21]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
    for (Int_t i=0; i<21; i++){
        hours2[i]=hours1[i]-0.05;
    }
    Double_t error_hours[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    //Filling the histogram by reading the datafile
    TH1F *hist_hour1[duration];
    TH1F *hist_hour2[duration];
    TH1F *hist_cumulated1;
    TH1F *hist_cumulated2;
    hist_cumulated1 = new TH1F ("hist_cumulated1", "energy", 2000, 0, 3000);
    hist_cumulated2 = new TH1F ("hist_cumulated2", "energy", 2000, 0, 3000);
    
    TF1 *f_resp1;
    TF2 *f_resp2;
    
    
    if(option==4){
        auto c = new TCanvas("c","Nb of counts with respect to the energy",200,10,700,500);
        
        FillHistogram(hist_hour1, my_plug, 1);
        FillHistogram(hist_hour2, my_plug, 2);
        
        FittingStep(hist_hour1, f_resp1, mean1, error_mean1, mean2, error_mean2, option);
        TGraphErrors *low_points1 = new TGraphErrors(n, hours1, mean1, error_hours, error_mean1);
        low_points1->SetTitle("mean evolution with respect to time");
        TColor *c3 = new TColor(9001,1,0,0); //red
        TColor *c31 = new TColor(9002,0,0,0); //black
        TColor *c32 = new TColor(9000,0,0,1); //blue
        low_points1->SetLineColor(9002);
        low_points1->Draw("a*");
        TGraphErrors *high_points1 = new TGraphErrors(n, hours1, mean2, error_hours, error_mean2);
        high_points1->SetLineColor(9002);
        high_points1->Draw("* SAME");
        
        
        FittingStep(hist_hour2, f_resp1, mean1, error_mean1, mean2, error_mean2, option);
        TGraphErrors *low_points2 = new TGraphErrors(n, hours2, mean1, error_hours, error_mean1);
        low_points2->SetTitle("mean evolution with respect to time");
        low_points2->SetLineColor(9000);
        low_points2->Draw("* SAME");
        TGraphErrors *high_points2 = new TGraphErrors(n, hours2, mean2, error_hours, error_mean2);
        high_points2->SetLineColor(9000);
        high_points2->Draw("* SAME");
        
        
        Cumulated_Histogram(hist_hour1,hist_cumulated1);
        Cumulated_Histogram(hist_hour2,hist_cumulated2);
        
        Double_t hour_cumulated[2] = {0,duration};
        
        FittingCumulated(hist_cumulated1, f_resp1, &mean_cumulated1, &error_mean_cumulated1, &mean_cumulated2, &error_mean_cumulated2);
        Double_t mean_cumulated_line1[2] = {mean_cumulated1,mean_cumulated1};
        Double_t mean_cumulated_line2[2]= {mean_cumulated2,mean_cumulated2};
        TGraph *low_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line1);
        TGraph *high_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line2);
        
        
        FittingCumulated(hist_cumulated2, f_resp2, &mean_cumulated1, &error_mean_cumulated1, &mean_cumulated2, &error_mean_cumulated2);
        mean_cumulated_line1[0] = mean_cumulated1;
        mean_cumulated_line1[1] = mean_cumulated1;
        mean_cumulated_line2[0] = mean_cumulated2;
        mean_cumulated_line2[1] = mean_cumulated2;
        TGraph *low_line2 = new TGraph(2,hour_cumulated,mean_cumulated_line1);
        TGraph *high_line2 = new TGraph(2,hour_cumulated,mean_cumulated_line2);
        
        
        low_line1->SetLineColor(9002);
        high_line1->SetLineColor(9002);
        low_line2->SetLineColor(9000);
        high_line2->SetLineColor(9000);
        low_line1->Draw("SAME");
        low_line2->Draw("SAME");
        high_line1->Draw("SAME");
        high_line2->Draw("SAME");
        auto legend = new TLegend(0.1,0.7,0.48,0.9);
           //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend->AddEntry(low_points1,"calibration 20220203","l"); //black
        legend->AddEntry(low_points2,"interpolation","l"); // blue
           legend->Draw();
        
        TFile f("./Data/Mean_fctof_Time.root","recreate");
        TTree *tree = new TTree("tree","Mean vs Hours");
        tree->Branch("means1",mean1,Form("mean1[%i]/D",duration));
        tree->Branch("error_means1",error_mean1,Form("error_mean1[%i]/D",duration));
        tree->Branch("means2",mean2,Form("mean2[%i]/D",duration));
        tree->Branch("error_means2",error_mean2,Form("error_mean2[%i]/D",duration));
        tree->Branch("hours",hours1,Form("hour[%i]/D",duration));
        tree->Branch("error_hours",error_hours,Form("error_hour[%i]/D",duration));
        tree->Branch("mean_accumulated_1Calib",&mean_cumulated1,"mean_accumulated_1Calib/D");
        tree->Branch("mean_accumulated_Interpol",&mean_cumulated2,"mean_accumulated_Interpol/D");
        
        tree->Fill();
        tree->Write();
        
        TCanvas *c2 = new TCanvas("c2","Histogram cumulation",200,10,700,500);
        c2->Draw ();
        c2->SetTickx ();
        c2->cd ()->SetTicky ();
        hist_cumulated1->SetLineColor(9001);
        hist_cumulated2->Draw("HIST"); // interpolated
        hist_cumulated1->Draw("SAME");
       
        //f_resp1->Draw("SAME");
        //f_resp2->Draw("SAME");
        
        
    }
    else if (option==5){
        auto c1 = new TCanvas("c1","Nb of counts with respect to the energy",200,10,700,500);
        
        FillHistogram(hist_hour1, my_plug, 3);
        FillHistogram(hist_hour2, my_plug, 2);
        
        hist_hour1[0]->Draw("HIST");
        for (Int_t i = 1; i<duration; i++){
            hist_hour1[i]->Draw("SAME");
        }
            
        auto c = new TCanvas("c","Nb of counts with respect to the energy",200,10,700,500);
        
        FittingStep(hist_hour1, f_resp1, mean1, error_mean1, mean2, error_mean2, option);
        TGraphErrors *low_points1 = new TGraphErrors(n, hours1, mean1, error_hours, error_mean1);
        low_points1->SetTitle("mean evolution with respect to time");
        TColor *c3 = new TColor(9001,1,0,0); //red
        TColor *c31 = new TColor(9000,0,0,1); //blue
        low_points1->SetLineColor(9001);
        low_points1->Draw("a*");
        TGraphErrors *high_points1 = new TGraphErrors(n, hours1, mean2, error_hours, error_mean2);
        high_points1->SetLineColor(9001);
        high_points1->Draw("* SAME");
        
        
        FittingStep(hist_hour2, f_resp1, mean1, error_mean1, mean2, error_mean2, option);
        TGraphErrors *low_points2 = new TGraphErrors(n, hours2, mean1, error_hours, error_mean1);
        low_points2->SetTitle("mean evolution with respect to time");
        low_points2->SetLineColor(9000);
        low_points2->Draw("* SAME");
        TGraphErrors *high_points2 = new TGraphErrors(n, hours2, mean2, error_hours, error_mean2);
        high_points2->SetLineColor(9000);
        high_points2->Draw("* SAME");
        
        
        Cumulated_Histogram(hist_hour1,hist_cumulated1);
        Cumulated_Histogram(hist_hour2,hist_cumulated2);
        
        Double_t hour_cumulated[2] = {start,duration};
        
        FittingCumulated(hist_cumulated1, f_resp1, &mean_cumulated1, &error_mean_cumulated1, &mean_cumulated2, &error_mean_cumulated2);
        Double_t mean_cumulated_line1[2] = {mean_cumulated1,mean_cumulated1};
        Double_t mean_cumulated_line2[2]= {mean_cumulated2,mean_cumulated2};
        TGraph *low_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line1);
        TGraph *high_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line2);
        
        
        FittingCumulated(hist_cumulated2, f_resp2, &mean_cumulated1, &error_mean_cumulated1, &mean_cumulated2, &error_mean_cumulated2);
        mean_cumulated_line1[0] = mean_cumulated1;
        mean_cumulated_line1[1] = mean_cumulated1;
        mean_cumulated_line2[0] = mean_cumulated2;
        mean_cumulated_line2[1] = mean_cumulated2;
        TGraph *low_line2 = new TGraph(2,hour_cumulated,mean_cumulated_line1);
        TGraph *high_line2 = new TGraph(2,hour_cumulated,mean_cumulated_line2);
        
        
        low_line1->SetLineColor(9001);
        high_line1->SetLineColor(9001);
        low_line2->SetLineColor(9000);
        high_line2->SetLineColor(9000);
        low_line1->Draw("SAME");
        low_line2->Draw("SAME");
        high_line1->Draw("SAME");
        high_line2->Draw("SAME");
        auto legend = new TLegend(0.1,0.7,0.48,0.9);
           //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
           legend->AddEntry(low_points1,"simulation","l"); //red
           legend->AddEntry(low_points2,"interpolation","l"); // blue
           legend->Draw();
        
        TFile f("./Data/Mean_fctof_Time.root","recreate");
        TTree *tree = new TTree("tree","Mean vs Hours");
        tree->Branch("means1",mean1,Form("mean1[%i]/D",duration-start));
        tree->Branch("error_means1",error_mean1,Form("error_mean1[%i]/D",duration-start));
        tree->Branch("means2",mean2,Form("mean2[%i]/D",duration-start));
        tree->Branch("error_means2",error_mean2,Form("error_mean2[%i]/D",duration-start));
        tree->Branch("hours",hours1,Form("hour[%i]/D",duration-start));
        tree->Branch("error_hours",error_hours,Form("error_hour[%i]/D",duration-start));
        tree->Branch("mean_accumulated_1Calib",&mean_cumulated1,"mean_accumulated_1Calib/D");
        tree->Branch("mean_accumulated_Interpol",&mean_cumulated2,"mean_accumulated_Interpol/D");
        
        tree->Fill();
        tree->Write();
        
        TCanvas *c2 = new TCanvas("c2","Histogram cumulation",200,10,700,500);
        c2->Draw ();
        c2->SetTickx ();
        c2->cd ()->SetTicky ();
        hist_cumulated1->SetLineColor(9001);
        hist_cumulated2->Draw("HIST"); // interpolated
        hist_cumulated1->Draw("SAME");
        
        //f_resp1->Draw("SAME");
        //f_resp2->Draw("SAME");
        
    }
    
    else{
        FillHistogram(hist_hour1, my_plug, option);
        
        FittingStep(hist_hour1, f_resp1, mean1, error_mean1, mean2, error_mean2, option);
        TGraphErrors *low_points1 = new TGraphErrors(n, hours1, mean1, error_hours, error_mean1);
        low_points1->SetTitle("mean evolution with respect to time");
        TColor *c3 = new TColor(9001,1,0,0);
        low_points1->SetLineColor(9001);
        low_points1->Draw("a*");
        TGraphErrors *high_points1 = new TGraphErrors(n, hours1, mean2, error_hours, error_mean2);
        high_points1->SetLineColor(9001);
        high_points1->Draw("* SAME");
        
        
        Cumulated_Histogram(hist_hour1,hist_cumulated1);
        
        Double_t hour_cumulated[2] = {0,duration};
        
        FittingCumulated(hist_cumulated1, f_resp1, &mean_cumulated1, &error_mean_cumulated1, &mean_cumulated2, &error_mean_cumulated2);
        Double_t mean_cumulated_line1[2] = {mean_cumulated1,mean_cumulated1};
        Double_t mean_cumulated_line2[2]= {mean_cumulated2,mean_cumulated2};
        TGraph *low_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line1);
        TGraph *high_line1 = new TGraph(2,hour_cumulated,mean_cumulated_line2);
        
        low_line1->SetLineColor(9001);
        high_line1->SetLineColor(9001);
        low_line1->Draw("SAME");
        high_line1->Draw("SAME");
        
        
        TFile f("./Data/Simulation_results.root","recreate");
        TTree *tree = new TTree("tree","Mean vs Hours");
        tree->Branch("means1",mean1,Form("mean1[%i]/D",duration));
        tree->Branch("error_means1",error_mean1,Form("error_mean1[%i]/D",duration));
        tree->Branch("means2",mean2,Form("mean2[%i]/D",duration));
        tree->Branch("error_means2",error_mean2,Form("error_mean2[%i]/D",duration));
        tree->Branch("hours",hours1,Form("hour[%i]/D",duration));
        tree->Branch("error_hours",error_hours,Form("error_hour[%i]/D",duration));
        tree->Branch("mean_accumulated_highpeek",&mean_cumulated1,"mean_accumulated_highpeek/D");
        tree->Branch("mean_accumulated_lowpeek",&mean_cumulated2,"mean_accumulated_lowpeek/D");
        
        tree->Fill();
        tree->Write();
        
        TCanvas *c2 = new TCanvas("c2","Histogram cumulation",200,10,700,500);
        c2->Draw ();
        c2->SetTickx ();
        c2->cd ()->SetTicky ();
        hist_cumulated1->SetLineColor(9001);
        hist_cumulated1->Draw("HIST");
    }
    
    
}
