//Macro that Analyse the root files
// 
// Before calling this macro the files 
//    DigOut.root   
//
// files should be created by running
//    raw_to_root.C

#include <stdio.h>
#include <fcntl.h>
#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TNtuple.h>
#include "Riostream.h"

#include <math.h>
#include "TMath.h"
#include <TRandom.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TLine.h"

#include "ROOT/RDataFrame.hxx" //to add tree RDataFrame.hxx
#include <iostream>

#include "WaveFormFunctions.C"

class Pulse : public TObject {
public:
  float adqTime;
  float fBase;  
  float fMaxBin;
  float fMax;
  float fInt;
  float fRCharge;
  float fFrontTime_50;
  float fTailTime_50; 
  float fT0;
  float fT0_lf;
 
  Pulse() { };
  ClassDef(Pulse,1);
};

void accessTreeInfo2(){
  TFile* file = new TFile("DigOut11.root");
  TChain *tree = new TChain("tree");
  // Open the ROOT file
  tree->Add("DigOut11.root");    
  //tree variables
  int run, event; 
  Pulse *B0 = new Pulse();   
  Pulse *B1 = new Pulse();
  Pulse *B2 = new Pulse();
  Pulse *B3 = new Pulse();  
  Pulse *B4 = new Pulse();  
  Pulse *B5 = new Pulse();  
  Pulse *B6 = new Pulse();  
  Pulse *B7 = new Pulse();  
  Pulse *B8 = new Pulse();  
  Pulse *B9 = new Pulse();  
  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);  


  tree->SetBranchAddress("B0",&B0);
  tree->SetBranchAddress("B1",&B1);
  tree->SetBranchAddress("B2",&B2);
  tree->SetBranchAddress("B3",&B3);  
  tree->SetBranchAddress("B4",&B4);  
  tree->SetBranchAddress("B5",&B5);  
  tree->SetBranchAddress("B6",&B6);  
  tree->SetBranchAddress("B7",&B7);  
  tree->SetBranchAddress("B8",&B8);  
  tree->SetBranchAddress("B9",&B9); 

 

  Long64_t nEvents = tree->GetEntries();
  int iEvent;   
  
  bool adqComplete[noCh];
  for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
  }


   Int_t i, iChannels = bins_per_record+1;

   //TH2F *h_Per[noCh];

/* TRYING TO DO IT FOR ALL CHANNELS AT ONCE
  for(i=0; i<noCh; i++) {
    hName  = Form("h_Per%01d", i);
    hTitle = Form("WaveForm Persistence %01d;Time,ns; Voltage, ADC", translator[i]);
    h_Per[i] = new TH2F(hName,hTitle,iChannels,0,bin_width*(bins_per_record-1),
                                      200,6500,8500);
    h_Per[i]= (TH2F*)file->Get("h_Pers[i]"); 

  } */

//histogramas que luego se llenan con h_pers (persistencias)
  TH2F *h_Per0 = new TH2F("h_Per0", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per1 = new TH2F("h_Per1", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per2 = new TH2F("h_Per2", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per3 = new TH2F("h_Per3", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per4 = new TH2F("h_Per4", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per5 = new TH2F("h_Per5", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per6 = new TH2F("h_Per6", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per7 = new TH2F("h_Per7", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per8 = new TH2F("h_Per8", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500); 
  TH2F *h_Per9 = new TH2F("h_Per9", "title",iChannels,0,bin_width*(bins_per_record-1),200,6500,8500);  

//histogramas para los valores màximos (amplitudes)
  auto hB0Max = new TH1F("hB0Max","Amplitude Distribution 4; Amplitude, mV; Entries",200,0.,8500.);
  auto hB1Max = new TH1F("hB1Max","Amplitude Distribution 5; Amplitude, mV; Entries",200,0.,8500.);
  auto hB2Max = new TH1F("hB2Max","Amplitude Distribution 6; Amplitude, mV; Entries",200,0.,8500.);
  auto hB3Max = new TH1F("hB3Max","Amplitude Distribution 7; Amplitude, mV; Entries",200,0.,8500.);
  auto hB4Max = new TH1F("hB4Max","Amplitude Distribution 9; Amplitude, mV; Entries",200,0.,8500.);
  auto hB5Max = new TH1F("hB5Max","Amplitude Distribution 10; Amplitude, mV; Entries",200,0.,8500.);
  auto hB6Max = new TH1F("hB6Max","Amplitude Distribution 11; Amplitude, mV; Entries",200,0.,8500.);
  auto hB7Max = new TH1F("hB7Max","Amplitude Distribution 13; Amplitude, mV; Entries",200,0.,8500.);
  auto hB8Max = new TH1F("hB8Max","Amplitude Distribution 14; Amplitude, mV; Entries",200,0.,8500.);
  auto hB9Max = new TH1F("hB9Max","Amplitude Distribution 15; Amplitude, mV; Entries",200,0.,8500.);


//histogramas para las distribuciones de diferencia de tiempos de arribo (fT0)
  auto hTA0 = new TH1F("hTA0","arrTime ch 13-9; Time, ns; Entries",200,-10.,10.); //2a-1a
  auto hTA1 = new TH1F("hTA1","arrTime ch 14-10; Time, ns; Entries",200,-10.,10.); //2b-1b
  auto hTA2 = new TH1F("hTA2","arrTime ch 15-11; Time, ns; Entries",200,-10.,10.); //2c-1c
  auto hTA3 = new TH1F("hTA3","arrTime ch 6-7; Time, ns; Entries",200,-10.,10.); //2b-2a
  auto hTA4 = new TH1F("hTA4","arrTime ch 4-5; Time, ns; Entries",200,-10.,10.); //0b-0a


//llenado de los histogramas de persistencia con histogramas ya creados
  h_Per0= (TH2F*)file->Get("h_Pers0"); 
  h_Per1= (TH2F*)file->Get("h_Pers1"); 
  h_Per2= (TH2F*)file->Get("h_Pers2"); 
  h_Per3= (TH2F*)file->Get("h_Pers3"); 
  h_Per4= (TH2F*)file->Get("h_Pers4"); 
  h_Per5= (TH2F*)file->Get("h_Pers5"); 
  h_Per6= (TH2F*)file->Get("h_Pers6"); 
  h_Per7= (TH2F*)file->Get("h_Pers7"); 
  h_Per8= (TH2F*)file->Get("h_Pers8"); 
  h_Per9= (TH2F*)file->Get("h_Pers9"); 

//llenado de los histogramas de distribucion de diff of arrTIme 

  for ( Int_t i = 0; i < nEvents ; i++)
  {
    tree->GetEntry(i);
    Double_t fT0diff0 = B7->fT0 - B4->fT0;  //diferencia de tiempo ch 13-9(2a-1a)
    Double_t fT0diff1 = B8->fT0 - B5->fT0;  //diferencia de tiempo ch 14-10(2b-1b)
    Double_t fT0diff2 = B9->fT0 - B6->fT0;  //diferencia de tiempo ch 15-11(2c-1c)
    Double_t fT0diff3 = B2->fT0 - B3->fT0;  //diferencia de tiempo ch 6-7(2b-2a)
    Double_t fT0diff4 = B0->fT0 - B1->fT0;  //diferencia de tiempo ch 4-5(0b-0a)
    hTA0->Fill(fT0diff0);
    hTA1->Fill(fT0diff1);
    hTA2->Fill(fT0diff2);
    hTA3->Fill(fT0diff3);
    hTA4->Fill(fT0diff4);
  }



//llenado de los histogramas de distribuciòn de amplitudes bajo la condiciòn
  for (iEvent=0; iEvent<nEvents; iEvent++) {
    tree->GetEntry(iEvent);
    for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
    }      
    if(B0->fMax > -10000.) adqComplete[0]=true;
    if(B1->fMax > -10000.) adqComplete[1]=true;
    if(B2->fMax > -10000.) adqComplete[2]=true;
    if(B3->fMax > -10000.) adqComplete[3]=true;
    if(B4->fMax > -10000.) adqComplete[4]=true;
    if(B5->fMax > -10000.) adqComplete[5]=true;
    if(B6->fMax > -10000.) adqComplete[6]=true;
    if(B7->fMax > -10000.) adqComplete[7]=true;
    if(B8->fMax > -10000.) adqComplete[8]=true;
    if(B9->fMax > -10000.) adqComplete[9]=true;
 
    if(adqComplete[0]) hB0Max->Fill(B0->fMax);
    if(adqComplete[1]) hB1Max->Fill(B1->fMax);
    if(adqComplete[2]) hB2Max->Fill(B2->fMax);
    if(adqComplete[3]) hB3Max->Fill(B3->fMax);
    if(adqComplete[4]) hB4Max->Fill(B4->fMax);
    if(adqComplete[5]) hB5Max->Fill(B5->fMax);
    if(adqComplete[6]) hB6Max->Fill(B6->fMax);
    if(adqComplete[7]) hB7Max->Fill(B7->fMax);
    if(adqComplete[8]) hB8Max->Fill(B8->fMax);
    if(adqComplete[9]) hB9Max->Fill(B9->fMax);


  }
  
  //display arrival time difference distribution
  auto c1 = new TCanvas();
  c1->Divide(2,3);

  c1->cd(1);
  hTA0->Draw();
  c1->cd(2);
  hTA1->Draw();
  c1->cd(3);
  hTA2->Draw();
  c1->cd(4);
  hTA3->Draw();
  c1->cd(5);
  hTA4->Draw();
  
  c1->Draw();  

  //display amplitude histograms

  auto c3 = new TCanvas("c3", "Amplitude Distribution", 10,10,1400,700);
  c3->Divide(2,noCh/2);

  c3->cd(1);
  hB0Max->Draw();
  c3->cd(2);
  hB1Max->Draw();
  c3->cd(3);
  hB2Max->Draw();
  c3->cd(4);
  hB3Max->Draw();
  c3->cd(5);
  hB4Max->Draw();
  c3->cd(6);
  hB5Max->Draw();
  c3->cd(7);
  hB6Max->Draw();
  c3->cd(8);
  hB7Max->Draw();
  c3->cd(9);
  hB8Max->Draw();
  c3->cd(10);
  hB9Max->Draw();
  
  c3->Draw();


//display persistence histograms
/*   for (int i = 0; i < 10; i++){
    h_Per[i]->GetXaxis()->SetRangeUser(4000, 7500);
  } */
  
  h_Per0->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per1->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per2->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per3->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per4->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per5->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per6->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per7->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per8->GetXaxis()->SetRangeUser(3800, 6500);
  h_Per9->GetXaxis()->SetRangeUser(3800, 6500);


  auto c2 = new TCanvas("c2", "Persistence", 1200, 800);
  c2->Divide(2, noCh / 2);
  c2->cd(1);
  h_Per0->Draw("colz");
  c2->cd(2);
  h_Per1->Draw("colz");
  c2->cd(3);
  h_Per2->Draw("colz");
  c2->cd(4);
  h_Per3->Draw("colz");
  c2->cd(5);
  h_Per4->Draw("colz");
  c2->cd(6);
  h_Per5->Draw("colz");
  c2->cd(7);
  h_Per6->Draw("colz");
  c2->cd(8);
  h_Per7->Draw("colz");
  c2->cd(9);
  h_Per8->Draw("colz");
  c2->cd(10);
  h_Per9->Draw("colz");
  
  c2->Draw();


/*   for(int i=1; i<noCh+1;i++){
    h_Per[i-1]->Write();
    h_Per[i-1]->GetXaxis()->SetRangeUser(4000, 7600);
    //h_Pers[i-1]->GetYaxis()->SetRangeUser(-200, 1600);

  } */

}
