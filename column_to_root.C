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

#include "WaveFormFunctions.C"

const Int_t nChannels = 4;

void jumpToLine(std::istream& os, int n){
	// Clear error flags, just in case.
	os.clear();	
	// Start reading from the beginning of the file.
	os.seekg(0, std::ios::beg);
	// Skip to line n.
	for (int i = 1; i < n; ++i)
	{
		os.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
}

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

void column_to_root() {

   Int_t iChannels = bins_per_record+1;
   //Here is the name of your file
   TString inputFile = "Tek009_ALL";  //0, 1, 2, 5, 6, 8, 9
   TString readFile = inputFile + ".csv";
   ifstream in;
   in.open(readFile);
//    in.open("t_waves0.csv");

// Values for baseline and pulse region

  int iBLfrom[nChannels];     // time window for baseline subrtaction
  int iBLto[nChannels];
  int iPULSEfrom[nChannels];  // time window for pulse
  int iPULSEto[nChannels];

  float fBase[nChannels];     // average baseline determined beforehand


  iBLfrom[0] = 0;  iBLto[0] =  280./bin_width;  iPULSEfrom[0] =300./bin_width;  iPULSEto[0] = 800./bin_width;
  iBLfrom[1] = 0;  iBLto[1] =  280./bin_width;  iPULSEfrom[1] =300./bin_width;  iPULSEto[1] = 800./bin_width;
  iBLfrom[2] = 0;  iBLto[2] =  300./bin_width;  iPULSEfrom[2] =450./bin_width;  iPULSEto[2] = 700./bin_width;
  iBLfrom[3] = 0;  iBLto[3] =  iBLto[0];  iPULSEfrom[3] = iPULSEfrom[0];  iPULSEto[3] = iPULSEto[0];

//Construction the histograms for the waveforms
  TH1F *h_WF[4];
  TString hName, hTitle;

  for(int i=0; i<4; i++) {
    hName  = Form("h_%01d", i);
    hTitle = Form("WaveForm %01d", i);
    h_WF[i] = new TH1F(hName,hTitle,iChannels,0,bin_width*(bins_per_record-1));
    h_WF[i]->SetStats(0);
    h_WF[i]->SetTitle(";Time, ns; Voltage, V ");
  }
   h_WF[0]->SetLineColor(kRed); h_WF[1]->SetLineColor(kCyan); h_WF[2]->SetLineColor(kBlue); h_WF[3]->SetLineColor(kGreen);

  //The variables that we need for the rootification
  Float_t time, wf0,wf1,wf2,wf3;
  char comma;
  float peak;
  int iEvent;
    // ** The variables that we need for the correct use of the waveformfuntions
  OSC_Record fWaveForm0,fWaveForm1,fWaveForm2,fWaveForm3;
  OSC_Record fWaveFormBL0,fWaveFormBL1,fWaveFormBL2,fWaveFormBL3;
  OSC_Record *pfWaveForm[4];
  OSC_Record *pfWaveFormBL[4];

  pfWaveForm[0]    = &fWaveForm0;
  pfWaveForm[1]    = &fWaveForm1;
  pfWaveForm[2]    = &fWaveForm2;
  pfWaveForm[3]    = &fWaveForm3;
  pfWaveFormBL[0]  = &fWaveFormBL0;
  pfWaveFormBL[1]  = &fWaveFormBL1;
  pfWaveFormBL[2]  = &fWaveFormBL2;
  pfWaveFormBL[3]  = &fWaveFormBL3;
   // **

// ---- Start tree contruction ----
  Pulse *B0 = new Pulse();
  Pulse *B1 = new Pulse();
  Pulse *B2 = new Pulse();
  Pulse *B3 = new Pulse();

  Pulse *pPulse[nChannels];

  pPulse[0] = B0;
  pPulse[1] = B1;
  pPulse[2] = B2;
  pPulse[3] = B3;

  TTree tTree("tree","Parameters of pulses in Oscilloscope input channels");
  int run = 1;
  tTree.Branch("run",&run,"i1/I");//''
  tTree.Branch("event",&iEvent,"i1/I");
  tTree.Branch("B0",&B0,10000,1);
  tTree.Branch("B1",&B1,10000,1);
  tTree.Branch("B2",&B2,10000,1);
  tTree.Branch("B3",&B3,10000,1);

// ---- End tree contruction ----

  	std::string lineContent;
	
  
  //rootification start here
  int iCh, nEvents;
  nEvents=40000;
  jumpToLine(in,13);
//   std::getline(in, lineContent);
//   cout << lineContent << endl;
  bool eof = false;
  
  auto cWF = new TCanvas();
//   std::getline(in, lineContent);
//   cout << lineContent << endl;
  
  for (int iEvent=0; iEvent < nEvents; iEvent++) {
//       cout << iEvent << endl;
//   while(!in.eof()){
    //fill the waveforms
      for(int i=0; i<iChannels; i++){
        in >> time >> comma >> wf0 >> comma >> wf1 >> comma >> wf2;
        if (!in.good()) break;
//         std::getline(in, lineContent);
//         if(i<20 && iEvent ==0)cout << lineContent << endl;        
//         if(i<<20 && iEvent==0) cout << time<<"  " << wf0<<"  " << wf1<<"  " << wf2 << endl;
        h_WF[0]->SetBinContent(i,wf0);
        h_WF[1]->SetBinContent(i,wf1);
        h_WF[2]->SetBinContent(i,wf2);
        h_WF[3]->SetBinContent(i,wf3);
        fWaveForm0.data[i]=wf0;
        fWaveForm1.data[i]=wf1;
        fWaveForm2.data[i]=wf2;
        fWaveForm3.data[i]=wf3;
      }

      for (iCh=0; iCh<nChannels; iCh++) {
        fBase[iCh] = GetBaseLine(pfWaveForm[iCh], iBLfrom[iCh], iBLto[iCh]);
        pPulse[iCh]->fBase = fBase[iCh];
        SubtractBaseLine(pfWaveForm[iCh], pfWaveFormBL[iCh], fBase[iCh]);
        InvertWaveForm(pfWaveFormBL[iCh], pfWaveFormBL[iCh]);
        pPulse[iCh]->fMaxBin = GetPeakPosition(pfWaveFormBL[iCh], iPULSEfrom[iCh], iPULSEto[iCh]);
        pPulse[iCh]->fMax = GetPeak(pfWaveFormBL[iCh], iPULSEfrom[iCh], iPULSEto[iCh]);
        pPulse[iCh]->fInt = GetIntegral(pfWaveFormBL[iCh], iPULSEfrom[iCh], iPULSEto[iCh], kFALSE);
        pPulse[iCh]->fRCharge = GetRCharge(pfWaveFormBL[iCh], iPULSEfrom[iCh], iPULSEto[iCh]);
        pPulse[iCh]->fT0_30 = GetFrontThresholdPosition(pfWaveFormBL[iCh], iPULSEfrom[iCh], iPULSEto[iCh], 0.3);
      }
      
      if(iEvent==0){
        for(int i=0; i<iChannels; i++){
            h_WF[1]->SetBinContent(i,fWaveForm0.data[i]);
            h_WF[2]->SetBinContent(i,fWaveFormBL0.data[i]);
            h_WF[3]->SetBinContent(i,fWaveFormBL1.data[i]);
            h_WF[3]->SetBinContent(i,fWaveFormBL2.data[i]);
        }
      }
//       cWF->cd(1);
//       h_WF[2]->Draw("same");   
//       h_WF[3]->Draw("same");   
//       h_WF[2]->Draw("same");   
//       gPad->Modified();
//       gPad->Update();      
//     if(in.eof()) break;

//       cout<< "fBase=" <<pPulse[1]->fBase<<endl;
//       cout<< "Maxbin=" <<pPulse[1]->fMaxBin<<endl;
//       cout<< "fMax="<<pPulse[1]->fMax<<endl;
//       cout<< "fInt="<<pPulse[1]->fInt<<endl;
//       cout<< "Time arrival="<<(pPulse[1]->fT0_30)*bin_width<<endl;

      tTree.Fill();
  }

  //El archivo de salida .root
  TString outputFile = inputFile + ".root";
  TFile*  file = new TFile(outputFile,"recreate");
  tTree.Write();
  h_WF[0]->Write();
  h_WF[1]->Write();
  h_WF[2]->Write();
  h_WF[3]->Write();
  file->Close();

}

