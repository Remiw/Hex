//Macro to display waveforms from the digitizer CAEN DT5560SE
//Remember to change the name of the input file <<in>> "*.csv"
//The variable <<noCh>> is the number of adquired channels
//The <<translator>> array is the channels that were adquired
//Exploring two functions to determine the array size and methods

#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TLine.h"

#include "Riostream.h"

#include <iostream>
#include <fstream>
#include <string>

#include <stdlib.h> //Atof funtion

#include "WaveFormFunctions.C" //details of the record and values to obtain

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

//WfDisplay
void raw_to_root(){

   gROOT->Reset();

  DigRecord fWaveForm, fWaveFormBL, fWaveFormInv; 
  DigRecordg *pfWaveForm, *pfWaveFormBL, *pfWaveFormInv;

  pfWaveForm    = &fWaveForm;
  pfWaveFormBL  = &fWaveFormBL;
  pfWaveFormInv = &fWaveFormInv;
  
  
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

  Pulse *pPulse[noCh];

  pPulse[0] = B0;
  pPulse[1] = B1;
  pPulse[2] = B2;
  pPulse[3] = B3;  
  pPulse[4] = B4;  
  pPulse[5] = B5;  
  pPulse[6] = B6;  
  pPulse[7] = B7;  
  pPulse[8] = B8;  
  pPulse[9] = B9;  
  
  int run = 1;
  int iEvent;  

  TTree tTree("tree","Parameters of the pulses in CAEN input channels");

  tTree.Branch("run",&run,"i1/I");
  tTree.Branch("event",&iEvent,"i1/I");
  tTree.Branch("B0",&B0,10000,1);
  tTree.Branch("B1",&B1,10000,1);
  tTree.Branch("B2",&B2,10000,1);
  tTree.Branch("B3",&B3,10000,1);
  tTree.Branch("B4",&B4,10000,1);
  tTree.Branch("B5",&B5,10000,1);
  tTree.Branch("B6",&B6,10000,1);
  tTree.Branch("B7",&B7,10000,1);
  tTree.Branch("B8",&B8,10000,1);
  tTree.Branch("B9",&B9,10000,1);
  
 
  Bool_t eof = false;

  int iBLfrom[noCh];     // time window for baseline subrtaction
  int iBLto[noCh];     
  int iPULSEfrom[noCh];  // time window for pulse
  int iPULSEto[noCh];

  float fBase[noCh];     // average baseline determined beforehand       
  
  iBLfrom[0] = 50/8;        iBLto[0] =  3000/8;   iPULSEfrom[0] = 3800/8;         iPULSEto[0] = 4400/8;
  iBLfrom[1] = iBLfrom[0];  iBLto[1] = iBLto[0];  iPULSEfrom[1] = iPULSEfrom[0];  iPULSEto[1] = iPULSEto[0];
  iBLfrom[2] = iBLfrom[0];  iBLto[2] = iBLto[0];  iPULSEfrom[2] = iPULSEfrom[0];  iPULSEto[2] = iPULSEto[0];
  iBLfrom[3] = iBLfrom[0];  iBLto[3] = iBLto[0];  iPULSEfrom[3] = iPULSEfrom[0];  iPULSEto[3] = iPULSEto[0];
  iBLfrom[4] = iBLfrom[0];  iBLto[4] = iBLto[0];  iPULSEfrom[4] = iPULSEfrom[0];  iPULSEto[4] = iPULSEto[0];
  iBLfrom[5] = iBLfrom[0];  iBLto[5] = iBLto[0];  iPULSEfrom[5] = iPULSEfrom[0];  iPULSEto[5] = iPULSEto[0];
  iBLfrom[6] = iBLfrom[0];  iBLto[6] = iBLto[0];  iPULSEfrom[6] = iPULSEfrom[0];  iPULSEto[6] = iPULSEto[0];
  iBLfrom[7] = iBLfrom[0];  iBLto[7] = iBLto[0];  iPULSEfrom[7] = iPULSEfrom[0];  iPULSEto[7] = iPULSEto[0];
  iBLfrom[8] = iBLfrom[0];  iBLto[8] = iBLto[0];  iPULSEfrom[8] = iPULSEfrom[0];  iPULSEto[8] = iPULSEto[0];
  iBLfrom[9] = iBLfrom[0];  iBLto[9] = iBLto[0];  iPULSEfrom[9] = iPULSEfrom[0];  iPULSEto[9] = iPULSEto[0];

   Int_t i, iChannels = bins_per_record+1;

    //abrir archivo csv
    ifstream in;
    in.open("data10.csv");
    
    if (!in.is_open()) {
        cerr << "Error, no se puede leer el archivo" << endl;
        return;
    }

  TCanvas *c1 = new TCanvas("c1","WaveForm",10,10,1400,700);
  c1->Divide(2,noCh/2);
//   gPad->SetGrid();

//Construction the histograms for the waveforms

  TH1F *h_WF[noCh];
  TH2F *h_Pers[noCh];
  Float_t min1,min2,min3,min4,th=2700./0.488;
  Float_t time1,time2,time3,time4;
  TString hName, hTitle,hNmin,hNmin_C,hTmin;
  Float_t waveform[bins_per_record + 6];
  for(int i=0; i<bins_per_record+adqDetails+1; i++){
      waveform[i]=-10000.;
  }
  Float_t auxWF=0, auxTime=0.;
  bool boolTime=false;

  for(i=0; i<noCh; i++) {
    hName  = Form("h_%01d", i);
    hTitle = Form("WaveForm %01d;Time,ns; Voltage, ADC", translator[i]);
    h_WF[i] = new TH1F(hName,hTitle,iChannels,0,bin_width*(bins_per_record-1));
    h_WF[i]->SetStats(0);    
    hName  = Form("h_Pers%01d", i);
    hTitle = Form("WaveForm Persistence %01d;Time,ns; Voltage, ADC", translator[i]);
    h_Pers[i] = new TH2F(hName,hTitle,iChannels,0,bin_width*(bins_per_record-1),
                                      200,6500,8500);
  }
  

  TObject *obj;
  int peak=-10;
  
  int auxJ=0;
  //To fill the tree we need an specific procedure for the data was not recorded
  bool signalRecord[noCh], newRecord=false;
  double timeOfTheRecord=0., timePreviousRecord=0.;
  double adqTimeofTheRecord[noCh];
  //To double check the zero values
  for(int j=0;j<noCh;j++){
      adqTimeofTheRecord[j]=0.;
  }
    
    //Start the tree inputs
//          while (!in.eof()) {
  for (iEvent=0; iEvent < 1e4; iEvent++) {
    //fill the waveforms
        for (int i = 0; i < adqDetails+bins_per_record+1; i++) {
            waveform[i] = 0.;
            in >> auxWF;
            in.ignore();
            waveform[i] = auxWF;
        }
        if(iEvent==0) timePreviousRecord=waveform[0]; //1st event
        
        timeOfTheRecord=waveform[0];
        
        if(timePreviousRecord<timeOfTheRecord){ //new record 
//             cout << timePreviousRecord<< "/// , ///"<<timeOfTheRecord << endl;
            newRecord=true;
            timePreviousRecord=waveform[0];
        }
        
        if(newRecord){
            for (int j = 0; j < noCh; j++){
                if(!signalRecord[j]){
                    pPulse[j]->adqTime = -10000.;
                    fBase[j] = -10000.;
                    pPulse[j]->fBase = -10000.;
                    pPulse[j]->fMaxBin = -10000.;
                    pPulse[j]->fMax = -10000.;
                    pPulse[j]->fInt = -10000.;
                    pPulse[j]->fRCharge = -10000.;
                    pPulse[j]->fFrontTime_50 = -10000.;
                    pPulse[j]->fTailTime_50 = -10000.;
                    pPulse[j]->fT0 = -10000.;
                    pPulse[j]->fT0_lf = -10000.;
                }
            }
            //Fill the tree with the pulses at the same adq Time
            tTree.Fill(); 
            newRecord=false;
            for (int j = 0; j < noCh; j++){
            signalRecord[j]=false;
            }
        }
        
//         if (!in.good()) break;
//waveform entries 
//[0] - Time
//[1] - channel 
//[2] - bins_per_record
//[3] - bool
//[4] - offset
//[5] - entry[0] of the waveform
//[5+1024-1] = [adqDetails+iChannels-2] = entry[1023] =0
//[5+1024] = [adqDetails+iChannels-1] = entry[1024] =0
        
        
        if (in.good() && waveform[3]==1 && waveform[adqDetails+iChannels-1]==0){ //The adquisition is correct
                if(boolTime && waveform[0] > auxTime){
//                     cout << "time condition:"<< waveform[0] <<"," << auxTime<<endl;
                    if(waveform[0] > auxTime) boolTime=false;
                }
                if(!boolTime){
                    for (int j = 0; j < noCh; j++){
                        if (translator[j] == waveform[1]){
                            adqTimeofTheRecord[j]=waveform[0];
                            for(int k=adqDetails; k<adqDetails+iChannels-1; k++){
                                h_WF[j]->SetBinContent(k,waveform[k]);
                                h_Pers[j]->Fill(k*bin_width,waveform[k]);
                                auxJ = j;
                                fWaveForm.data[k]=waveform[k];
                                signalRecord[j]=true;
                           }
//                             peak=GetPeakPosition(fWaveForm, 400,550, -1);
//                         cout <<"peak="<< peak << endl;
//                         cout<< "min="<<GetPeak(fWaveForm, 400,550, -1)<<endl;
                           //Define the waveform values               
                           
      //Pulse data
      pPulse[j]->adqTime = adqTimeofTheRecord[0];
      fBase[j] = GetBaseLine(pfWaveForm, iBLfrom[j], iBLto[j]);
      pPulse[j]->fBase = fBase[j];
      SubtractBaseLine(pfWaveForm, pfWaveFormBL, fBase[j]);
	  InvertWaveForm(pfWaveFormBL, pfWaveFormBL);
      pPulse[j]->fMaxBin = GetPeakPosition(fWaveFormBL, iPULSEfrom[j], iPULSEto[j]);
      pPulse[j]->fMax = GetPeak(fWaveFormBL, iPULSEfrom[j], iPULSEto[j]);      
      pPulse[j]->fInt = GetIntegral(fWaveFormBL, iPULSEfrom[j], iPULSEto[j], kTRUE);
      pPulse[j]->fRCharge = GetRCharge(fWaveFormBL, iPULSEfrom[j], iPULSEto[j]);      
      pPulse[j]->fFrontTime_50 = GetFrontThresholdPosition(fWaveFormBL, iPULSEfrom[j], iPULSEto[j], 0.5);
      pPulse[j]->fTailTime_50 = GetTailThresholdPosition(fWaveFormBL, iPULSEfrom[j], iPULSEto[j], 0.5);
      pPulse[j]->fT0 = GetTime0(fWaveFormBL, iPULSEfrom[j], iPULSEto[j], 0.2, 0.5); 
      pPulse[j]->fT0_lf = GetTime0_lf(fWaveFormBL, iPULSEfrom[j], iPULSEto[j], 0.2, 0.5); 
                        }
                }
                //cout << "event=" << iEvent << endl;
                c1->cd(auxJ+1);
                h_WF[auxJ]->Draw("same");
                //gPad->Modified();
                //gPad->Update();
                }
        }
        else{
            if(!boolTime) auxTime=waveform[0];
            boolTime=true;
//             in.ignore(1000, ';0;0');
            jumpToLine(in, iEvent+1);
           // cout << "**********salto=" << iEvent << endl;
           // cout <<"***********details="<<waveform[0]<<","<<waveform[1]<<","<<waveform[2]<<","<<waveform[3]<<endl;
//             obj = gPad->WaitPrimitive(); //cout<<iEntries<<endl;
        }
//             cout <<"todos="<<waveform[0]<<","<<waveform[1]<<","<<waveform[2]<<","<<waveform[3]<<endl;
   //Para esperar por cada waveform

            
  }
  
  auto c2 = new TCanvas();
  c2->Divide(2,noCh/2);
  for(int i=1; i<noCh+1;i++){
    c2->cd(i);
    h_Pers[i-1]->Draw("colz");
  }
  
  TString OutFileName = Form("DigOut10.root");
  TFile*  file = new TFile(OutFileName,"recreate");
  tTree.Write();
  for(int i=1; i<noCh+1;i++){
    h_Pers[i-1]->Write();
    //h_Pers[i-1]->GetXaxis()->SetRangeUser(4000, 7600);
    //h_Pers[i-1]->GetYaxis()->SetRangeUser(-200, 1600);

  }
  file->Close(); 
  


}
