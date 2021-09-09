#include "MRemovePUTail.hh"
#include "QEvent.hh"
#include "QEventList.hh"
#include "QHeader.hh"
#include "QFitter.hh"
#include "QPulse.hh"
#include "QPulseInfo.hh"
#include "QEventList.hh"
#include "QPulseParameters.hh"
#include "QBaselineData.hh"
#include "QBaseType.hh"
#include "QAveragePulseHandle.hh"
#include "QRunDataHandle.hh"
#include "QAverageVector.hh"
#include "QVector.hh"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include <cmath>

REGISTER_MODULE(MRemovePUTail)

using namespace Cuore;

//double *AvgPulseFunction(double x, double *pulse, QVector pulseVec){
double AvgPulseFunction(double x, std::vector<double> pulseVec){

  double val;
  val = pulseVec[x];
  return val;

}

double AvgPulseFit(double *x, double *par){

  // parameters are: par[0] = amplitude of normal pulse, par[1] = array from QVector.GetArray(), par[2] = QVector of average pulse
  
  double xx = x[0];
  double ampOfEventPulse = par[0];
  std::vector<double> pulse;
  for (int i = 2; i<par[1]; i++){
    pulse.push_back(par[i]);
  }
  //  pulse[-1] = pulse[0]; // set fit value back to original baseline to avoid large jump in fit for last element
  double ampOfAvgPulseInt = AvgPulseFunction(xx, pulse);
  //  std::cout << "amplitude pulse is " << ampOfEventPulse << ", and amplitude of average pulse is " << ampOfAvgPulseInt << ". Returned fit value is " << ampOfEventPulse * ampOfAvgPulseInt/ampOfAvgPulseInt << std::endl;
  //  return (ampOfEventPulse * ampOfAvgPulseInt/ampOfAvgPulseInt);
  return (ampOfEventPulse * ampOfAvgPulseInt);
  
}

void MRemovePUTail::Init(QEvent& ev)
{
  ev.Require<QHeader>("DAQ", "Header");
  ev.Require<QPulseParameters>("PulseBasicParameters", "Parameters");
  ev.Require<QPulseInfo>("DAQ", "PulseInfo");
  ev.Require<QPulse>("DAQ", "Pulse");
  ev.Require<QBaselineData>("BaselineModule", "BaselineData");
  ev.Require<QDouble>("PulseBasicParameters", "MaxMinInWindow");
  
  fAvgPulseInput = GetString("AvgPulseInput", "default");
  fDataset = GetInt("Dataset", -1, false);
  fThisRun = GetInt("Run", -1, false);
  fMaxTau = GetDouble("MaxTau", -1, false);
  
  ev.Add<QVector>("testPUPulse").SetWrite(GetBool("SaveSamples", false));
  
}

void MRemovePUTail::Do(QEvent& ev, const QEventList& neighbours)
{
  
  Debug("event number %i", ev.GetReadNumber());
  std::cout << Form("event number %i", ev.GetReadNumber()) << std::endl;
  
  double hardMax = fMaxTau; 
  double maxTail;
  if (!neighbours.Empty()) {
    
    //Get trigger time of the channel
    const QHeader& evHeader = ev.Get<QHeader>("DAQ", "Header");
    const unsigned long long trigTime = evHeader.GetTime().GetFromStartRunNs();
    const QPulseInfo& pInf = ev.Get<QPulseInfo>("DAQ", "PulseInfo");
    const QPulse& pulseData = ev.Get<QPulse>("DAQ", "Pulse");
    
    const int run = evHeader.GetRun();
    
    QRunDataHandle rHan(run);
    GlobalData().Get<QRunData>("", &rHan, "");
    
    QVector pulse = pulseData.GetSamples();
    Debug("other event contains %d elements ", neighbours.Size());
    
    for(size_t iNeighbour = 0; iNeighbour < neighbours.Size(); iNeighbour++){
      QEvent otherEv = neighbours[iNeighbour];
      Debug("we're on iNeighbour %d", iNeighbour);
      const QHeader& otherHeader = otherEv.Get<QHeader>("DAQ", "Header");
      const QPulseInfo& otherpInf = otherEv.Get<QPulseInfo>("DAQ", "PulseInfo");
      const unsigned long long otherTrigTime = otherHeader.GetTime().GetFromStartRunNs();
      
      
      if(otherTrigTime > trigTime && ev.GetReadNumber()!= 3){ //only events with this selection can cause pile-up
	
	
	Debug("run over pulse %d", ev.GetReadNumber());
	
	int ch = pInf.GetChannelId();
	
	if(fDataset < 0){
	  GlobalHandle<QInt> dHan("Dataset");
	  GlobalData().Get("", &dHan, "");
	  fDataset = dHan.Get();
	}
	
	Debug("Dataset is %d", fDataset);
	
	
	QPulseParameters otherEvParams = otherEv.Get<QPulseParameters>("PulseBasicParameters", "Parameters");
	QPulse otherEvPulseData = otherEv.Get<QPulse>("DAQ", "Pulse");
	double decayTimeNs = otherEvParams.fSlowDecayTime*1E9;
	double riseTimeNs = otherEvParams.fRiseTime*1E9;
	QDouble amplitude = otherEv.Get<QDouble>("PulseBasicParameters", "MaxMinInWindow");
	
	Debug("trigTime is %llu, decayTime is %d, riseTime is %d, amplitude is %d, and other trigTime is %llu", (double)trigTime, decayTimeNs, riseTimeNs,(double)amplitude, (double)otherTrigTime);
	
	double tauSpacing = (trigTime - (otherTrigTime+riseTimeNs))/decayTimeNs;
	Debug("Tau for event is %d", tauSpacing);
	Debug("Test debug diana");
	
	QAveragePulseHandle apHan(ch);
	apHan.SetDataset(fDataset);
	GlobalData().Get(GetString("APOwner", "AveragePulses",false),&apHan, fAvgPulseInput);
	Debug("got global handle Average Pulse");
	
	QAverageVector avgPulse = apHan.Get();
	QVector otherPulse = otherEvPulseData.GetSamples();
	
	double otherPulseArray[otherPulse.Size()];
	double xvaluesForOtherPulse[otherPulse.Size()];
	for(int i = 0; i<otherPulse.Size(); i++){
	  otherPulseArray[i] = otherPulse[i];
	  xvaluesForOtherPulse[i] = i;
	}
	
	Debug("got other pulse sample info");
	
	avgPulse.ShiftReal(otherPulse.GetMaxIndex());
	
	
	double fitparms[2+avgPulse.Size()];
	double ampdoub = (double)(amplitude);
	double ampdoubsize = (double)(avgPulse.Size());
	
	fitparms[0] = ampdoub;
	fitparms[1] = ampdoubsize;
	
	for (int i = 0; i < avgPulse.Size(); i++){
	  //	  for (int y = 0; i < otherPulse.Size(); i++){
	  //	  if(otherPulse.Size()>avgPulse.Size()){
	  //	    while(i<otherPulse.Size() && i>avgPulse.Size()){
	  //	      fitparms[i+2]=avgPulse[0];
	  //	    }
	  //	  if (avgPulse[i] = avgPulse[otherPulse.Size()-1])
	  //	    fitparms[2+i] = avgPulse[0];
	  //	  else
	    fitparms[2+i] = avgPulse[i];
	    //	  }

	}

	std::cout << "fit parameter size is " << sizeof(fitparms)/sizeof(fitparms[0]) << ", and size of other event pulse is " << otherPulse.Size() << std::endl;
		 	
	double samplingFreq = rHan.Get().GetChannelRunData(ch).fSamplingFrequency;
	Debug("got sampling frequency");
	
	TF1 *apfit = new TF1("avg", AvgPulseFit, 0., double(pulse.Size()-1), (pulse.Size()+2));
	apfit->SetParameters(fitparms);
	Debug("define fit");
	
	
	TGraph *otherPgraph = new TGraph(otherPulse.Size(), xvaluesForOtherPulse, otherPulseArray);
	Debug("got graph of pulse");
	
	TCanvas *pulseCan = new TCanvas("pulseCan", "", 500, 500);
	pulseCan->cd();
	otherPgraph->Draw();
	apfit->Draw("histsame");
	pulseCan->SaveAs("testpulse.root");
	pulseCan->Close();

	TCanvas *avgCan = new TCanvas("avgCan", "", 500, 500);
	avgCan->cd();
	apfit->Draw();
	avgCan->SaveAs("testavg.root");
	avgCan->Close();
	
	otherPgraph->Fit(apfit);
	Debug("fit pulse");
	double calculatedAmp = apfit->GetParameter(0);
	
	QVector avgPulsePrime = calculatedAmp * avgPulse.ShiftReal(tauSpacing/samplingFreq);
	
	
	for(int i=0; i<pulse.Size(); i++){
	  pulse[i] -= avgPulsePrime[i];
	  
	}
	
	ev.Get<QVector>("testPUPulse") = pulse;
	
	//	delete fitparms;
	delete apfit;
	
	
      }
    }
    /*   TGraph *testPulseGraph = pulse.GetGraph();
	 TCanvas *can = new TCanvas(Form("pulse%d", ev), "", 500, 500);
	 can->cd();
	 testPulseGraph->Draw();
	 delete can;
	 delete testPulseGraph; */
  }
}


void MRemovePUTail::Done()
{
  // anything saved in do loop
}
