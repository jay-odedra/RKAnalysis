#define afterWeights_cxx
#include "afterWeights.h"
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void afterWeights::Loop()
{
  if (fChain == 0) return;

  // Analysis BDT distributions before weight
  TH1F *bdtOPFPF = new TH1F("bdtOPFPF","bdtOPFPF", 30, -15., 15.);
  TH1F *bdtGPFPF = new TH1F("bdtGPFPF","bdtGPFPF", 30, -15., 15.);
  bdtOPFPF->Sumw2();
  bdtGPFPF->Sumw2();
  TH1F *bdtOPFLPT = new TH1F("bdtOPFLPT","bdtOPFLPT", 30, -15., 15.);
  TH1F *bdtGPFLPT = new TH1F("bdtGPFLPT","bdtGPFLPT", 30, -15., 15.);
  bdtOPFLPT->Sumw2();
  bdtGPFLPT->Sumw2();

  // Analysis BDT distributions after weight
  TH1F *wbdtOPFPF = new TH1F("wbdtOPFPF","wbdtOPFPF", 30, -15., 15.);
  TH1F *wbdtGPFPF = new TH1F("wbdtGPFPF","wbdtGPFPF", 30, -15., 15.);
  wbdtOPFPF->Sumw2();
  wbdtGPFPF->Sumw2();
  TH1F *wbdtOPFLPT = new TH1F("wbdtOPFLPT","wbdtOPFLPT", 30, -15., 15.);
  TH1F *wbdtGPFLPT = new TH1F("wbdtGPFLPT","wbdtGPFLPT", 30, -15., 15.);
  wbdtOPFLPT->Sumw2();
  wbdtGPFLPT->Sumw2();

  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Loop over Bs 
    for (unsigned int ii=0; ii<fit_mass->size(); ii++) {

      // signal events only
      if (bmatchMC->at(ii)==0) continue;

      // Fix outliers
      float theWeightG = analysisBdtGWithSyst->at(ii);
      float theWeightO = analysisBdtOWithSyst->at(ii);
      if (theWeightG>5) { cout << "Out G: " << theWeightG << endl; theWeightG=5; }
      if (theWeightO>5) { cout << "Out O: " << theWeightO << endl; theWeightO=5; }
      
      // PFPF
      if (probe_pfmvaId->at(ii)<20 && tag_pfmvaId->at(ii)<20) {
	bdtOPFPF->Fill(analysisBdtO->at(ii));
	wbdtOPFPF->Fill(analysisBdtO->at(ii),theWeightO);
	bdtGPFPF->Fill(analysisBdtG->at(ii));
	wbdtGPFPF->Fill(analysisBdtG->at(ii),theWeightG);

      } else { // PFLPT
	bdtOPFLPT->Fill(analysisBdtO->at(ii));
	wbdtOPFLPT->Fill(analysisBdtO->at(ii),theWeightO);
	bdtGPFLPT->Fill(analysisBdtG->at(ii));
	wbdtGPFLPT->Fill(analysisBdtG->at(ii),theWeightG);
      }

    } // Loop over Bs

  } // Loop over events


  // Compute efficiencies for cut @5, check [chiara]
  cout << "Low edge for bin 21 = " << wbdtOPFPF->GetBinLowEdge(21) << endl;
  
  cout << "---------------------------------" << endl;
  Double_t error = 1.;
  cout << "PFPF: Eff Otto, no weight = "      << bdtOPFPF->IntegralAndError(21, 30, error, "")/bdtOPFPF->IntegralAndError(1, 30, error, "")     << endl;
  cout << "PFPF: Eff Otto, with weight = "    << wbdtOPFPF->IntegralAndError(21, 30, error, "")/wbdtOPFPF->IntegralAndError(1, 30, error, "")   << endl;
  cout << "PFPF: Eff George, no weight = "    << bdtGPFPF->IntegralAndError(21, 30, error, "")/bdtGPFPF->IntegralAndError(1, 30, error, "")     << endl;
  cout << "PFPF: Eff George, with weight = "  << wbdtGPFPF->IntegralAndError(21, 30, error, "")/wbdtGPFPF->IntegralAndError(1, 30, error, "")   << endl;
  cout << "PFLPT: Eff Otto, no weight = "     << bdtOPFLPT->IntegralAndError(21, 30, error, "")/bdtOPFLPT->IntegralAndError(1, 30, error, "")   << endl;
  cout << "PFLPT: Eff Otto, with weight = "   << wbdtOPFLPT->IntegralAndError(21, 30, error, "")/wbdtOPFLPT->IntegralAndError(1, 30, error, "") << endl;
  cout << "PFLPT: Eff George, no weight = "   << bdtGPFLPT->IntegralAndError(21, 30, error, "")/bdtGPFLPT->IntegralAndError(1, 30, error, "")   << endl;
  cout << "PFLPT: Eff George, with weight = " << wbdtGPFLPT->IntegralAndError(21, 30, error, "")/wbdtGPFLPT->IntegralAndError(1, 30, error, "") << endl;
  cout << "---------------------------------" << endl;


  gStyle->SetOptStat(0);

  TCanvas c1("c1","",1);
  bdtOPFPF->SetTitle("");
  wbdtOPFPF->SetTitle("");
  bdtOPFPF->GetXaxis()->SetTitle("BDT [Otto]");
  wbdtOPFPF->GetXaxis()->SetTitle("BDT [Otto]");
  bdtOPFPF->SetLineColor(2);
  wbdtOPFPF->SetLineColor(1);
  bdtOPFPF->Scale(1./bdtOPFPF->Integral());
  wbdtOPFPF->Scale(1./wbdtOPFPF->Integral());   
  wbdtOPFPF->Draw();
  bdtOPFPF->Draw("same");
  c1.SaveAs("otto_PFPF.png");

  TCanvas c2("c2","",1);
  bdtGPFPF->SetTitle("");
  wbdtGPFPF->SetTitle("");
  bdtGPFPF->GetXaxis()->SetTitle("BDT [George]");
  wbdtGPFPF->GetXaxis()->SetTitle("BDT [George]");
  bdtGPFPF->SetLineColor(2);
  wbdtGPFPF->SetLineColor(1);
  bdtGPFPF->Scale(1./bdtGPFPF->Integral());
  wbdtGPFPF->Scale(1./wbdtGPFPF->Integral());   
  wbdtGPFPF->Draw();
  bdtGPFPF->Draw("same");
  c2.SaveAs("george_PFPF.png");

  TCanvas c1b("c1b","",1);
  bdtOPFLPT->SetTitle("");
  wbdtOPFLPT->SetTitle("");
  bdtOPFLPT->GetXaxis()->SetTitle("BDT [Otto]");
  wbdtOPFLPT->GetXaxis()->SetTitle("BDT [Otto]");
  bdtOPFLPT->SetLineColor(2);
  wbdtOPFLPT->SetLineColor(1);
  bdtOPFLPT->Scale(1./bdtOPFLPT->Integral());
  wbdtOPFLPT->Scale(1./wbdtOPFLPT->Integral());   
  bdtOPFLPT->Draw();
  wbdtOPFLPT->Draw("same");
  c1b.SaveAs("otto_PFLPT.png");

  TCanvas c2b("c2b","",1);
  bdtGPFLPT->SetTitle("");
  wbdtGPFLPT->SetTitle("");
  bdtGPFLPT->GetXaxis()->SetTitle("BDT [George]");
  wbdtGPFLPT->GetXaxis()->SetTitle("BDT [George]");
  bdtGPFLPT->SetLineColor(2);
  wbdtGPFLPT->SetLineColor(1);
  bdtGPFLPT->Scale(1./bdtGPFLPT->Integral());
  wbdtGPFLPT->Scale(1./wbdtGPFLPT->Integral());   
  bdtGPFLPT->Draw();
  wbdtGPFLPT->Draw("same");
  c2b.SaveAs("george_PFLPT.png");
}
