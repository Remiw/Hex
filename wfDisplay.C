#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TLine.h"
#include "Riostream.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> // Atof function

int const bins_per_record = 1250-1;
float const bin_width = 0.16; // Sample Interval, ns

// noCh is the number of adq channels and has to be even
int const noCh = 4;
Int_t translator[noCh] = {1, 2, 3, 4};

typedef struct {
  std::array<float, bins_per_record> data;
} DigRecord;

// wf Display

void wfDisplay(int nEvents) {

  gROOT->Reset();

  // fWaveForm tiene un miembro data que es un array de punto flotante de tamaño bins_per_record.
  DigRecord fWaveForm; 

  // abrir archivo csv
  ifstream in;
  in.open("Tek000_ALL.csv");

  Float_t time, ch1, ch2, ch3, ch4;
  Int_t iChannels = bins_per_record+1;

  // Ignorar las primeras 11 líneas que no contienen datos de señal
  for (int i = 0; i < 12; i++) {
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  TCanvas *c1 = new TCanvas("c1", "WaveForm", 10, 10, 1400, 700);
  c1->Divide(2, noCh / 2);

  // Construcción de los histogramas
  TH1F *h_WF[noCh];  // señal para cada canal
  TH2F *h_Pers[noCh]; // persistencias para cada canal
  TString hName, hTitle;

  for (int j = 0; j < noCh; j++) {
    hName = Form("h_%01d", j);
    hTitle = Form("Waveform Ch %01d;time,ns;Amp,mV", translator[j]);
    h_WF[j] = new TH1F(hName, hTitle, iChannels, 0, bin_width * (bins_per_record - 1));
    h_WF[j]->SetStats(0);
    h_Pers[j] = new TH2F(Form("hPers_%01d", j), hTitle, iChannels, 0, bin_width * (bins_per_record - 1), 200, -1, 1);
  }

  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    // Llenar los histogramas de las formas de onda

    for (int i = 0; i < bins_per_record; i++) {
      in >> time;
      in.ignore(); // Ignora la coma
      in >> ch1;
      in.ignore(); // Ignora la coma
      in >> ch2;
      in.ignore(); // Ignora la coma
      in >> ch3;
      in.ignore(); // Ignora la coma
      in >> ch4;
      if (!in.good()) break;
      // time = waveform[0], ch = waveform[1]

      h_WF[0]->SetBinContent(i, ch1); // Canal 1
      h_WF[1]->SetBinContent(i, ch2); // Canal 2
      h_WF[2]->SetBinContent(i, ch3); // Canal 3
      h_WF[3]->SetBinContent(i, ch4); // Canal 4

      h_Pers[0]->Fill(i*bin_width, ch1);
      h_Pers[1]->Fill(i*bin_width, ch2);
      h_Pers[2]->Fill(i*bin_width, ch3);
      h_Pers[3]->Fill(i*bin_width, ch4);
    }

    // Actualizar todos los histogramas juntos
    for (int j = 0; j < noCh; j++) {
      c1->cd(j + 1);
      h_WF[j]->Draw();
    }
    c1->Update();
  }

  auto c2 = new TCanvas();
  c2->Divide(2, noCh / 2);

  for (int i = 0; i < noCh; i++) {
    c2->cd(i + 1);
    h_Pers[i]->Draw("colz");
  }

  c2->Update();
  in.close();
}
