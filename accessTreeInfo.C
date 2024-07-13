//Macro that Analyse the root files
// 
// Before calling this macro the files 
//    DigOut.root   
//
// files should be created by running
//    column_to_root.C

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
  float fBase;  
  float fMaxBin;
  float fMax;
  float fInt;
  float fRCharge;
  float fT0_30;
 
  Pulse() { };
  ClassDef(Pulse,1);
};

const Int_t noCh = 4;

void accessTreeInfo(){
  TFile* file = new TFile("Tek000_ALL.root");
  TChain *tree = new TChain("tree");
  // Open the ROOT file
  tree->Add("Tek000_ALL.root");   //0, 1, 2, 5, 6, 8, 9
  tree->Add("Tek001_ALL.root"); 
  tree->Add("Tek002_ALL.root");
  tree->Add("Tek005_ALL.root"); 
  tree->Add("Tek006_ALL.root"); 
  tree->Add("Tek008_ALL.root"); 
  tree->Add("Tek009_ALL.root"); 


  //tree variables
  int run, event; 
  Pulse *B0 = new Pulse();   
  Pulse *B1 = new Pulse();
  Pulse *B2 = new Pulse();
  Pulse *B3 = new Pulse();  
  
  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);  
  tree->SetBranchAddress("B0",&B0);
  tree->SetBranchAddress("B1",&B1);
  tree->SetBranchAddress("B2",&B2);
  tree->SetBranchAddress("B3",&B3); 

 

  Long64_t nEvents = tree->GetEntries();
  int iEvent;   
  
  bool adqComplete[noCh];
  for(int j=0; j<noCh; j++){
      adqComplete[j]=false;
  }


   Int_t i, iChannels = bins_per_record+1;

   //TH2F *h_Per[noCh];


//histogramas que luego se llenan con h_pers (persistencias)
//  TH2F *h_Per0 = new TH2F("h_Per0", "title",iChannels,0,bin_width*(bins_per_record-1),200,-2,1); 
//  TH2F *h_Per1 = new TH2F("h_Per1", "title",iChannels,0,bin_width*(bins_per_record-1),200,-2,1); 
//  TH2F *h_Per2 = new TH2F("h_Per2", "title",iChannels,0,bin_width*(bins_per_record-1),200,-2,1); 
//  TH2F *h_Per3 = new TH2F("h_Per3", "title",iChannels,0,bin_width*(bins_per_record-1),200,-2,1); 

//histogramas para los valores màximos (amplitudes)
  auto hB0Max = new TH1F("hB0Max","Amplitude Distribution 1; Amplitude, mV; Entries",200,-8000.,8000.);
  auto hB1Max = new TH1F("hB1Max","Amplitude Distribution 2; Amplitude, mV; Entries",200,-20.,20.);
  auto hB2Max = new TH1F("hB2Max","Amplitude Distribution 3; Amplitude, mV; Entries",200,-20.,20.);
  auto hB3Max = new TH1F("hB3Max","Amplitude Distribution 4; Amplitude, mV; Entries",200,-20.,20.);


//histogramas para las distribuciones de diferencia de tiempos de arribo (fT0)
  auto hTA0 = new TH1F("hTA0","fT0_30 ch 1-2; Time, ns; Entries",200,0,2000); 
  auto hTA1 = new TH1F("hTA1","fT0_30 ch 3-4; Time, ns; Entries",200,0,2000); 



//llenado de los histogramas de persistencia con histogramas ya creados
//  h_Per0= (TH2F*)file->Get("h_Pers0"); 
//  h_Per1= (TH2F*)file->Get("h_Pers1"); 
//  h_Per2= (TH2F*)file->Get("h_Pers2"); 
//  h_Per3= (TH2F*)file->Get("h_Pers3");
  
//llenado de los histogramas de distribucion de diff of arrTIme 

  for ( Int_t i = 0; i < nEvents ; i++)
  {
    tree->GetEntry(i);
    float fT0diff0 = B0->fT0_30 - B1->fT0_30;  //diferencia de tiempo ch 1-2
    float fT0diff1 = B2->fT0_30 - B3->fT0_30;  //diferencia de tiempo ch 3-4

    hTA0->Fill(fT0diff0);
    hTA1->Fill(fT0diff1);

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

 
    if(adqComplete[0]) hB0Max->Fill(B0->fMax);
    if(adqComplete[1]) hB1Max->Fill(B1->fMax);
    if(adqComplete[2]) hB2Max->Fill(B2->fMax);
    if(adqComplete[3]) hB3Max->Fill(B3->fMax);



  }
  
  //display arrival time difference distribution
  auto c1 = new TCanvas();
  c1->Divide(2);

  c1->cd(1);
  hTA0->Draw();
  c1->cd(2);
  hTA1->Draw();

  
  c1->Draw();  

  //display amplitude histograms

  auto c3 = new TCanvas("c3", "Amplitude Distribution");
  c3->Divide(2,noCh/2);

  c3->cd(1);
  hB0Max->Draw();
  c3->cd(2);
  hB1Max->Draw();
  c3->cd(3);
  hB2Max->Draw();
  c3->cd(4);
  hB3Max->Draw();

  
  c3->Draw();


//display persistence histograms
/*   for (int i = 0; i < 10; i++){
    h_Per[i]->GetXaxis()->SetRangeUser(4000, 7500);
  } */
  
//  h_Per0->GetXaxis()->SetRangeUser(3800, 6500);
//  h_Per1->GetXaxis()->SetRangeUser(3800, 6500);
//  h_Per2->GetXaxis()->SetRangeUser(3800, 6500);
//  h_Per3->GetXaxis()->SetRangeUser(3800, 6500);



//  auto c2 = new TCanvas("c2", "Persistence", 1200, 800);
//  c2->Divide(2, noCh / 2);
//  c2->cd(1);
//  h_Per0->Draw("colz");
//  c2->cd(2);
//  h_Per1->Draw("colz");
//  c2->cd(3);
//  h_Per2->Draw("colz");
//  c2->cd(4);
//  h_Per3->Draw("colz");

  
//  c2->Draw();


/*   for(int i=1; i<noCh+1;i++){
    h_Per[i-1]->Write();
    h_Per[i-1]->GetXaxis()->SetRangeUser(4000, 7600);
    //h_Pers[i-1]->GetYaxis()->SetRangeUser(-200, 1600);

  } */

}
