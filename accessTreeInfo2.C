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
  tree->Add("DigOut.root");
  tree->Add("DigOut1.root");
  tree->Add("DigOut2.root");
  tree->Add("DigOut3.root");
  tree->Add("DigOut4.root");
  tree->Add("DigOut5.root");
  tree->Add("DigOut6.root");
  tree->Add("DigOut7.root");
  tree->Add("DigOut8.root");
  tree->Add("DigOut9.root");
  tree->Add("DigOut10.root");
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

  tree->SetBranchAddress("B0",&B0);  //ch4
  tree->SetBranchAddress("B1",&B1);  //ch5
  tree->SetBranchAddress("B2",&B2);  //ch6 
  tree->SetBranchAddress("B3",&B3);  //ch7 
  tree->SetBranchAddress("B4",&B4);  //ch9
  tree->SetBranchAddress("B5",&B5);  //ch10
  tree->SetBranchAddress("B6",&B6);  //ch11
  tree->SetBranchAddress("B7",&B7);  //ch13
  tree->SetBranchAddress("B8",&B8);  //ch14
  tree->SetBranchAddress("B9",&B9);  //ch15

 

  Long64_t nEvents = tree->GetEntries();
  int iEvent;   
  
  bool adqComplete[noCh];
  for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
  }


  Int_t i, iChannels = bins_per_record+1;

  //histogramas para los valores màximos (amplitudes)
  auto hB0Max = new TH1F("hB0Max","Amplitude Distribution 4; Amplitude, mV; Entries",200,0.,5500.);
  auto hB1Max = new TH1F("hB1Max","Amplitude Distribution 5; Amplitude, mV; Entries",200,0.,5500.);
  auto hB2Max = new TH1F("hB2Max","Amplitude Distribution 6; Amplitude, mV; Entries",200,0.,5500.);
  auto hB3Max = new TH1F("hB3Max","Amplitude Distribution 7; Amplitude, mV; Entries",200,0.,5500.);
  auto hB4Max = new TH1F("hB4Max","Amplitude Distribution 9; Amplitude, mV; Entries",200,0.,5500.);
  auto hB5Max = new TH1F("hB5Max","Amplitude Distribution 10; Amplitude, mV; Entries",200,0.,5500.);
  auto hB6Max = new TH1F("hB6Max","Amplitude Distribution 11; Amplitude, mV; Entries",200,0.,5500.);
  auto hB7Max = new TH1F("hB7Max","Amplitude Distribution 13; Amplitude, mV; Entries",200,0.,5500.);
  auto hB8Max = new TH1F("hB8Max","Amplitude Distribution 14; Amplitude, mV; Entries",200,0.,5500.);
  auto hB9Max = new TH1F("hB9Max","Amplitude Distribution 15; Amplitude, mV; Entries",200,0.,5500.);

  //histogramas para las distribuciones de diferencia de tiempos de arribo (fT0)
  
  //auto hTA0 = new TH1F("hTA0","arrTime ch 13-9(2a-1a); Time, ns; Entries",100,-10,10); //2a-1a
  //auto hTA1 = new TH1F("hTA1","arrTime ch 14-10(2b-1b); Time, ns; Entries",100,-10,10); //2b-1b

  auto hTA0 = new TH1F("hTA0","  (13+14+15)-(9+10+11); Time, ns; Entries",100,-10,10); //2a-1a
  auto hTA1 = new TH1F("hTA1","(6+7)-(4+5); Time, ns; Entries",100,-1,1); //2b-1b
  auto hTA2 = new TH1F("hTA2","arrTime ch 15-11(2c-1c); Time, ns; Entries",100,-10,10); //2c-1c
  auto hTA3 = new TH1F("hTA3","arrTime ch 6-7(2b-2a); Time, ns; Entries",100,-10,10); //2b-2a
  auto hTA4 = new TH1F("hTA4","arrTime ch 4-5(0b-0a); Time, ns; Entries",100,-10,10); //0b-0a


  //llenado de los histogramas de persistencia con histogramas recuperados del archivo root
  TH2F *h_Per0= (TH2F*)file->Get("h_Pers0"); 
  TH2F *h_Per1= (TH2F*)file->Get("h_Pers1"); 
  TH2F *h_Per2= (TH2F*)file->Get("h_Pers2"); 
  TH2F *h_Per3= (TH2F*)file->Get("h_Pers3"); 
  TH2F *h_Per4= (TH2F*)file->Get("h_Pers4"); 
  TH2F *h_Per5= (TH2F*)file->Get("h_Pers5"); 
  TH2F *h_Per6= (TH2F*)file->Get("h_Pers6"); 
  TH2F *h_Per7= (TH2F*)file->Get("h_Pers7"); 
  TH2F *h_Per8= (TH2F*)file->Get("h_Pers8"); 
  TH2F *h_Per9= (TH2F*)file->Get("h_Pers9"); 

  //llenado de los histogramas de distribuciòn de amplitudes bajo la condiciòn
  for (iEvent=0; iEvent<nEvents; iEvent++) {
    tree->GetEntry(iEvent);

    for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
    }      

    if(B0->fMax>-10000) {adqComplete[0]=true;}
    if(B1->fMax>-10000) {adqComplete[1]=true;}
    if(B2->fMax>-10000) {adqComplete[2]=true;}
    if(B3->fMax>-10000) {adqComplete[3]=true;}
    if(B4->fMax>-10000) {adqComplete[4]=true;}
    if(B5->fMax>-10000) {adqComplete[5]=true;}
    if(B6->fMax>-10000) {adqComplete[6]=true;}
    if(B7->fMax>-10000) {adqComplete[7]=true;}
    if(B8->fMax>-10000) {adqComplete[8]=true;}
    if(B9->fMax>-10000) {adqComplete[9]=true;}

    if(adqComplete[0]) {hB0Max->Fill(B0->fMax);}
    if(adqComplete[1]) {hB1Max->Fill(B1->fMax);}
    if(adqComplete[2]) {hB2Max->Fill(B2->fMax);}
    if(adqComplete[3]) {hB3Max->Fill(B3->fMax);}
    if(adqComplete[4]) {hB4Max->Fill(B4->fMax);}
    if(adqComplete[5]) {hB5Max->Fill(B5->fMax);}
    if(adqComplete[6]) {hB6Max->Fill(B6->fMax);}
    if(adqComplete[7]) {hB7Max->Fill(B7->fMax);}
    if(adqComplete[8]) {hB8Max->Fill(B8->fMax);}
    if(adqComplete[9]) {hB9Max->Fill(B9->fMax);}

  

   double_t fT0diff0 = 0., fT0diff1 = 0., fT0diff2=0., fT0diff3=0., fT0diff4=0.;


    for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
    }

    if(B0->fT0>-10000.) adqComplete[0] = true;
    if(B1->fT0>-10000.) adqComplete[1] = true;
    if(B2->fT0>-10000.) adqComplete[2] = true;
    if(B3->fT0>-10000.) adqComplete[3] = true;
    if(B4->fT0>-10000.) adqComplete[4] = true;
    if(B5->fT0>-10000.) adqComplete[5] = true;
    if(B6->fT0>-10000.) adqComplete[6] = true;
    if(B7->fT0>-10000.) adqComplete[7] = true;
    if(B8->fT0>-10000.) adqComplete[8] = true;
    if(B9->fT0>-10000.) adqComplete[9] = true;

    double_t arrTCh0=0., arrTCh1=0., arrTCh2=0., arrTCh3=0., arrTCh4=0., arrTCh5=0., arrTCh6=0., arrTCh7=0., arrTCh9=0., arrTCh10=0., arrTCh11=0., arrTCh13=0., arrTCh14=0., arrTCh15=0;
    if(adqComplete[0]) {arrTCh4= B0->fT0;}
    if(adqComplete[1]) {arrTCh5= B1->fT0;}
    if(adqComplete[2]) {arrTCh6= B2->fT0;}
    if(adqComplete[3]) {arrTCh7= B3->fT0;}
    if(adqComplete[4]) {arrTCh9= B4->fT0;}
    if(adqComplete[5]) {arrTCh10= B5->fT0;}
    if(adqComplete[6]) {arrTCh11= B6->fT0;}
    if(adqComplete[7]) {arrTCh13= B7->fT0;}
    if(adqComplete[8]) {arrTCh14= B8->fT0;}
    if(adqComplete[9]) {arrTCh15= B9->fT0;}

    double_t C2 = 0., C1 = 0., C3 = 0., C4 = 0.;
    if(adqComplete[7] ) {
      C2 = (arrTCh13 + arrTCh14 + arrTCh15);
      C1 = (arrTCh9 + arrTCh10 + arrTCh11);
      hTA0->Fill(C2-C1);
    }
    if(adqComplete[6] ) {
      C3 = (arrTCh6);
      C4 = (arrTCh4);
      hTA1->Fill(C3-C4);
    }
    //if(adqComplete[7] && adqComplete[4]) { 
    //  fT0diff0 = (B7->fT0 - B4->fT0);  //diferencia de tiempo ch 13-9(2a-1a)
    //  hTA0->Fill(fT0diff0);
    //}
    //if(adqComplete[8] && adqComplete[5]) {
    //  fT0diff1 = (B8->fT0 - B5->fT0);   //diferencia de tiempo ch 14-10(2b-1b)
    //   hTA1->Fill(fT0diff1);
    //}
    //if(adqComplete[9] && adqComplete[6]) {
    //  fT0diff2 = (B9->fT0 - B6->fT0);  //diferencia de tiempo ch 15-11(2c-1c)
    //  hTA2->Fill(fT0diff2);
    //}
    //if(adqComplete[2] && adqComplete[3]) {
    //  fT0diff3 = (B2->fT0 - B3->fT0);  //diferencia de tiempo ch 6-7(2b-2a)
    //  hTA3->Fill(fT0diff3);
    //}
    //if(adqComplete[0] && adqComplete[1]) {
    //  fT0diff4 = (B0->fT0 - B1->fT0);  //diferencia de tiempo ch 4-5(0b-0a)
    //  hTA4->Fill(fT0diff4);
    //}
  }
  //B4);  ch9
  //B5);  ch10
  //B6);  ch11
  //B7);  ch13
  //B8);  ch14
  //B9);  ch15

  //display arrival time difference distribution (resolution time = sigma/sqrt(2))
  auto c1 = new TCanvas();
  c1->Divide(2,1);

  c1->cd(1);
  hTA0->Fit("gaus");
  hTA0->Draw();
  c1->cd(2);
  hTA1->Fit("gaus");
  hTA1->Draw();
  //c1->cd(3);
  //hTA2->Fit("gaus");
  //hTA2->Draw();
  //c1->cd(4);
  //hTA3->Fit("gaus");
  //hTA3->Draw();
  //c1->cd(5);
  //hTA4->Fit("gaus");
  //hTA4->Draw();
  
  c1->Draw();  

  //display amplitude histograms


  auto c3 = new TCanvas("c3", "Amplitude Distribution", 10,10,1400,700);
  c3->Divide(2,noCh/2);

  c3->cd(1);
  gPad->SetLogy();  
  hB0Max->SetFillStyle(3001);
  hB0Max->SetFillColor(kRed);
  hB0Max->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hB1Max->SetFillStyle(3001);
  hB1Max->SetFillColor(kRed);
  hB1Max->Draw();
  c3->cd(3);
  gPad->SetLogy();
  hB2Max->SetFillStyle(3001);
  hB2Max->SetFillColor(kRed);
  hB2Max->Draw();
  c3->cd(4);
  gPad->SetLogy();
  hB3Max->SetFillStyle(3001);
  hB3Max->SetFillColor(kRed);
  hB3Max->Draw();
  c3->cd(5);
  gPad->SetLogy();
  hB4Max->SetFillStyle(3001);
  hB4Max->SetFillColor(kRed);
  hB4Max->Draw();
  c3->cd(6);
  gPad->SetLogy();
  hB5Max->SetFillStyle(3001);
  hB5Max->SetFillColor(kRed);
  hB5Max->Draw();
  c3->cd(7);
  gPad->SetLogy();
  hB6Max->SetFillStyle(3001);
  hB6Max->SetFillColor(kRed);
  hB6Max->Draw();
  c3->cd(8);
  gPad->SetLogy();
  hB7Max->SetFillStyle(3001);
  hB7Max->SetFillColor(kRed);
  hB7Max->Draw();
  c3->cd(9);
  gPad->SetLogy();
  hB8Max->SetFillStyle(3001);
  hB8Max->SetFillColor(kRed);
  hB8Max->Draw();
  c3->cd(10);
  gPad->SetLogy();
  hB9Max->SetFillStyle(3001);
  hB9Max->SetFillColor(kRed);
  hB9Max->Draw();
  
  //c3->Draw();


  //display persistence histograms

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
  gPad->SetLogz();
  h_Per0->Draw("colz");
  c2->cd(2);
  gPad->SetLogz();
  h_Per1->Draw("colz");
  c2->cd(3);
  gPad->SetLogz();
  h_Per2->Draw("colz");
  c2->cd(4);
  gPad->SetLogz();
  h_Per3->Draw("colz");
  c2->cd(5);
  gPad->SetLogz();
  h_Per4->Draw("colz");
  c2->cd(6);
  gPad->SetLogz();
  h_Per5->Draw("colz");
  c2->cd(7);
  gPad->SetLogz();
  h_Per6->Draw("colz");
  c2->cd(8);
  gPad->SetLogz();
  h_Per7->Draw("colz");
  c2->cd(9);
  gPad->SetLogz();
  h_Per8->Draw("colz");
  c2->cd(10);
  gPad->SetLogz();
  h_Per9->Draw("colz");
  
  //c2->Draw();

}
