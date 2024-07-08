#include "TH1.h"
#include "TF1.h"

// time, ch, bins, bool, offset; + wfdata [1024]
int const adqDetails= 5; 
int const bins_per_record = 1024-1; 
float const bin_width = 8.; // Sample Interval, ns

// noCh is the number of adq iChannels has to be odd
int const noCh = 10;
// const Int_t nChannels = 4;
//translator helps with the identification of the adq channels
Int_t translator[noCh] = {4, 5, 6, 7, 9, 10, 11, 13, 14, 15};

//baseline values
int BLMin = 8150;
int BLMax = 8220;

typedef struct {
  std::array<float,bins_per_record> data;
} DigRecord;

// we skip 2 lines when the adq was wrong with jumpToLine
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






int CopyWaveForm(DigRecord* pWFin, DigRecord* pWFout) {
  for (int i=0; i<1024; i++) pWFout->data[i] = pWFin->data[i];
  return 0;
}

void InvertWaveForm(DigRecord* pWFin, DigRecord* pWFout) {
  for (int i=0; i<1024; i++) pWFout->data[i] = -pWFin->data[i];
}

int SumUpWaveForms(DigRecord* pWFin_1, DigRecord* pWFin_2, DigRecord* pWFout) {
  for (int i=0; i<1024; i++) pWFout->data[i] = pWFin_1->data[i] + pWFin_2->data[i];
  return 0;
}


int CheckChannelRange(int base_from, int base_to) {
  if(base_from < 0)    { cout << "unexpected first bin for base " << base_from << endl; return -1; }
  if(base_from > 1024) { cout << "unexpected first bin for base " << base_from << endl; return -1; }
  if(base_to > 1024)   { cout << "unexpected last bin for base " << base_to << endl; return -1; }
  if(base_to < 0)      { cout << "unexpected last bin for base " << base_to << endl; return -1; }
  if(base_from > base_to) { 
    cout << "unexpected interval for base: " << base_from << " .. " << base_to << endl; 
    return -1;
  } 
  return 0;
}

float GetBaseLine(DigRecord* pWFin, int base_from, int base_to) {

  if(CheckChannelRange(base_from, base_to) != 0) return -100000.;

  int iBaseFrom = base_from;
  int iBaseTo = base_to;

  TH1F  *hTmp = new TH1F("hTmp","Pedestal",800,BLMin,BLMax);
  int i;
  float base = 0.;
  for (i=iBaseFrom; i<=iBaseTo; i++) { hTmp->Fill(pWFin->data[i]); }
//   Int_t iBin = hTmp->GetMaximumBin();
//   base = hTmp->GetXaxis()->GetBinCenter(iBin);
//   hTmp->Delete();
//   return base;
  Double_t iBin = hTmp->GetMean();
  hTmp->Delete();
  return iBin;

}



void SubtractBaseLine(DigRecord* pWFin, DigRecord* pWFout, float base) {
  for (int i=0; i<1024; i++) pWFout->data[i] = pWFin->data[i] - base;
  return;
}


float SubtractBaseLine(DigRecord* pWFin, DigRecord* pWFout, int base_from, int base_to) {

  float base = GetBaseLine(pWFin, base_from, base_to);
  for (int i=0; i<1024; i++) pWFout->data[i] = pWFin->data[i] - base;
  return base;

}


float GetIntegral(DigRecord rec, int from, int to, bool only_positive = kFALSE) {

  // if only_positive is true, negative part of the pulse is not added to the sum
  // in this case one has to expect slight additional baseline shift

  if(CheckChannelRange(from, to) != 0) return -1000000.;

  float sum = 0.;
  for (int i=from; i<=to; i++) { 
    sum = sum + rec.data[i];
    if(only_positive && (rec.data[i] < 0)) sum = sum - rec.data[i]; 
  }
  return sum;

}

float GetRCharge(DigRecord rec, int from, int to) {


  if(CheckChannelRange(from, to) != 0) return -1000000.;

  float sum = 0.;
  for (int i=from; i<=to; i++) { 
    sum = sum + rec.data[i];
    sum = sum - rec.data[i]; 
  }
  return sum;

}


int GetPeakPosition(DigRecord rec, int from, int to, int sign=1) {

  int iFrom = from;
  int iTo = to;

  int iMinPosition = -20;
  int iMaxPosition = -20;

  float max = -50000.;
  float min =  50000.;

  for (int i=iFrom; i<=iTo; i++) {
    if (max<rec.data[i]) { max = rec.data[i]; iMaxPosition = i; }
    if (min>rec.data[i]) { min = rec.data[i]; iMinPosition = i; }
  }

  if(sign < 0) return iMinPosition;
  return iMaxPosition;

}

float GetPeak(DigRecord rec, int from, int to, int sign=1) {

  int iPeak = GetPeakPosition(rec, from, to, sign);

  return rec.data[iPeak];

}


float GetFrontThresholdPosition(DigRecord rec, int from, int to, float threshold) {

  // WaveForm is expected to have a positive pulse and a substracted baseline
  // before calling use SubstractBaseline() and InvertWaveForm() if needed 
 
  // threshold is a fraction from the peak value, should be from 0.05 to 0.95

  if((threshold < 0.05) || (threshold > 0.95)) {
    cout << "GetFrontThresholdPosition() error. Unexpected threshold value " << threshold << endl;
    return -1.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThreshold = rec.data[iPeak]*threshold;

  float front; 
  int   i, iAbove, iBelow;

  iBelow = from;
  for (i=iPeak; i>=from; i--) { 
    if ((rec.data[i] < fThreshold) && (rec.data[i+1] >= fThreshold)) {
      iBelow = i;
      break;
    }
  }  
  iAbove = iPeak;
  for (i=iPeak; i>=from; i--) { 
    if ((rec.data[i-1] <= fThreshold) && (rec.data[i] > fThreshold)) {
      iAbove = i;
      break;
    }
  }  
  front = iBelow + (iAbove-iBelow)*(fThreshold - rec.data[iBelow])/(rec.data[iAbove]-rec.data[iBelow]);
  
  return front;
}

float GetTailThresholdPosition(DigRecord rec, int from, int to, float threshold) {

  // WaveForm is expected to have a positive pulse and a substracted baseline
  // before calling use SubstractBaseline() and InvertWaveform() if needed 
 
  // threshold is a fraction from the peak value, should be from 0.05 to 0.95


  if((threshold < 0.05) || (threshold > 0.95)) {
    cout << "GetTailThresholdPosition() error. Unexpected threshold value " << threshold << endl;
    return -1000000.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThreshold = rec.data[iPeak]*threshold;

  float tail; 
  int   i, iAbove, iBelow;

  iBelow = to;
  for (i=iPeak; i<=to; i++) { 
    if ((rec.data[i] < fThreshold) && (rec.data[i-1] >= fThreshold)) {
      iBelow = i;
      break;
    }
  }  
  iAbove = iPeak;
  for (i=iPeak; i<=to; i++) { 
    if ((rec.data[i+1] <= fThreshold) && (rec.data[i] > fThreshold)) {
      iAbove = i;
      break;
    }
  }  
  tail = iAbove + (iBelow-iAbove)*(rec.data[iAbove]-fThreshold)/(rec.data[iAbove]-rec.data[iBelow]);
  
  return tail;
}


float GetTime0 (DigRecord rec, int from, int to, float threshold_1, float threshold_2) {

  // points of two thresholds crossings are found (fTime_1,fThr_1) and (fTime_2,fThr_2)
  // extrapolation of the straight line between them to thr == 0 is returned as fTime_0. 

  // WaveForm is expected to have a positive pulse and a substracted baseline
  // before calling use SubstractBaseline() and InvertWaveform() if needed 
 
  // each threshold is a fraction from the peak value, should be from 0.05 to 0.95
      
  if((threshold_1 < 0.05) || (threshold_1 > 0.95)) {
    cout << "GetTime0() error. Unexpected threshold value " << threshold_1 << endl;
    return -1000000.;
  }

  if((threshold_2 < 0.05) || (threshold_2 > 0.95)) {
    cout << "GetTime0() error. Unexpected threshold value " << threshold_2 << endl;
    return -1000000.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThr_1 = rec.data[iPeak]*threshold_1;
  float fThr_2 = rec.data[iPeak]*threshold_2;

  float fTime_1 = GetFrontThresholdPosition(rec, from, to, threshold_1);
  float fTime_2 = GetFrontThresholdPosition(rec, from, to, threshold_2);

  return (fThr_2*fTime_1 - fThr_1*fTime_2)/(fThr_2-fThr_1); 
   
}

float GetTime0_lf (DigRecord rec, int from, int to, float threshold_1, float threshold_2) {

  // The same as GetTime0(), but makes a linear fit in the range between thresholds
      
  if((threshold_1 < 0.05) || (threshold_1 > 0.95)) {
    cout << "GetTime0() error. Unexpected threshold value " << threshold_1 << endl;
    return -1000000.;
  }

  if((threshold_2 < 0.05) || (threshold_2 > 0.95)) {
    cout << "GetTime0() error. Unexpected threshold value " << threshold_2 << endl;
    return -1000000.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThr_1 = rec.data[iPeak]*threshold_1;
  float fThr_2 = rec.data[iPeak]*threshold_2;

  float fTime_1 = GetFrontThresholdPosition(rec, from, to, threshold_1);
  float fTime_2 = GetFrontThresholdPosition(rec, from, to, threshold_2);

  TH1F *hWF = new TH1F("hWF","Temporary histogram",to - from,from,to);
  for (int i = 0; i < (to-from); i++) hWF->SetBinContent(i,rec.data[from+i]);

  TF1 *FitFunc = new TF1("FitFunc","(x-[0])*TMath::Tan([1])");
  FitFunc->SetRange(fTime_1,fTime_2);
  hWF->Fit("FitFunc","RQN");

  float t0 = FitFunc->GetParameter(0);

  hWF->Delete();
  FitFunc->Delete();

  return t0;
}

int TakePeakPosition(DigRecord* pWfFin, int from, int to, int sign = -1) {

  int iFrom = from;
  int iTo = to;

  int iMinPosition = 0;
  int iMaxPosition = 0;

  float max = -1024.;
  float min =  1024.;

  for (int i=iFrom; i<=iTo; i++) { 
    if (max<pWfFin->data[i]) { max = pWfFin->data[i]; iMaxPosition = i; }
    if (min>pWfFin->data[i]) { min = pWfFin->data[i]; iMinPosition = i; }
  }

  if(sign < 0) return iMinPosition;
  return iMaxPosition; 

}

float TakePeak(DigRecord* pWfFin, int from, int to, int sign ) {

  int iPeak = TakePeakPosition(pWfFin, from, to, sign); 
 
  return pWfFin->data[iPeak]; 

}


float GetFrontBaseLine(DigRecord* pWfFin,int from, int to, int sign) {
    
    int iBaseFrom = TakePeakPosition(pWfFin, from, to, sign)-20;
    int iBaseTo = iBaseFrom+80; 
    
    
  int i,counter1,counter2;
  counter1=0; counter2=0;
  TH1F *hTmp1 = new TH1F("hTmp1","Pedestal1",500,0.,3500.);
  for (i=0; i<=iBaseTo; i++){
      hTmp1->Fill(pWfFin->data[i]);
      counter1++;
      }
      
  Double_t iBin1 = hTmp1->GetMean();
//   cout<<"iBin1="<<iBin1<<endl;
//   cout<<"counter1="<<counter1<<endl;
    
  hTmp1->Delete();
 
  return iBin1;
  
}

float GetBiggestBaseLine(DigRecord* pWfFin,int from, int to, int sign) {

    
    int iBaseFrom = TakePeakPosition(pWfFin, from, to, sign)-20;
    int iBaseTo = iBaseFrom+80; 
    
    
  int i,counter1,counter2;
  counter1=0; counter2=0;
  TH1F *hTmp1 = new TH1F("hTmp1","Pedestal1",500,0.,3500.);
  for (i=0; i<=iBaseTo; i++){
      hTmp1->Fill(pWfFin->data[i]);
      counter1++;
      }
      
  TH1F *hTmp2 = new TH1F("hTmp2","Pedestal2",500,0.,3500.);
  for (i=iBaseFrom; i<=1024; i++){
      hTmp2->Fill(pWfFin->data[i]); 
      counter2++;
      }

  Double_t iBin1 = hTmp1->GetMean();
  Double_t iBin2 = hTmp2->GetMean();
  Double_t base;
  cout<<"iBin1="<<iBin1<<",     "<<"iBin2="<<iBin2<<endl;
  cout<<"counter1="<<counter1<<",     "<<"counter2="<<counter2<<endl;
  if(counter1>=counter2){
      base=iBin1;
  }
  else{
      base=iBin2;
  }
  cout<<"base="<<base<<endl;  
    
  hTmp1->Delete();
  hTmp2->Delete();
 
  return base;
  
}

void SubtractFrontBaseLine(DigRecord* pWfFin, DigRecord* pWFout, float base) {
  for (int i=0; i<1024; i++) pWFout->data[i] = pWfFin->data[i] - base;
  return;
}

void Normalize(DigRecord* pWFin, DigRecord* pWFout, float peak) {
  for (int i=0; i<1024; i++) pWFout->data[i] = pWFin->data[i]/peak;
  return;
}


// ---- FOR FIX THRESHOLD


float GetTdcStart(DigRecord rec, int from, int to, float threshold) {

  float fMax = GetPeak(rec, from, to);
  float fFraction = threshold/fMax;

  if(fFraction > 0.95) return -1000.;
  if(fFraction < 0.05) fFraction = 0.05;
//Note the fraction is necesarly to make usulfull the previus funtion
//but the threshold is fixed to the original value
  float fStart = GetFrontThresholdPosition(rec, from, to, fFraction);

  return fStart;

}


float GetTdcWidth(DigRecord rec, int from, int to, float threshold) {

  float fMax = GetPeak(rec, from, to);
  float fFraction = threshold/fMax;

  if(fFraction > 0.95) return -1000.;
  if(fFraction < 0.05) fFraction = 0.05;

  float fStart = GetFrontThresholdPosition(rec, from, to, fFraction);
  float fStop = GetTailThresholdPosition(rec, from, to, fFraction);

  return fStop - fStart;

}
