#define mcVsSplots_cxx

// Data vs MC signal distributions

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void drawTH1pair(TH1* h1, TH1* h2, 
		 const string& xAxisNameTmp = "", const string& yAxisName = "Events", 
		 float lumi=-1, const string& canvasName = "default", 
		 const string& outputDIR = "./", int mycolor=2,
		 const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC", const string& outputFILE = "outFile.root") 
{
  string xAxisName = "";
  string separator = "::";
  Bool_t setXAxisRangeFromUser = false;
  Double_t xmin = 0;
  Double_t xmax = 0;

  size_t pos = xAxisNameTmp.find(separator);
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos); 
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    separator = ",";
    pos = xrange.find(separator);
    string numString = ""; 
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + separator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  if (yAxisName == "a.u.") {
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
  }
  else if (lumi>-1) {
    h1->Scale(lumi/h1->Integral());
    h2->Scale(lumi/h2->Integral());
  }

  // To deal with splots
  for (int ii=0; ii<h1->GetNbinsX(); ii++){
    int iip1 = ii+1;
    if (h1->GetBinContent(iip1)<=0) { 
      h1->SetBinContent(iip1,0);
      h1->SetBinError(iip1,0);
    } 
    if (h2->GetBinContent(iip1)<=0) {
      h2->SetBinContent(iip1,0);
      h2->SetBinError(iip1,0);
    } 
  }
  
  h1->SetStats(0);
  h2->SetStats(0);

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->SetTitle("");
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),h2->GetMaximum()) * 1.2);
  if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->Draw("EP");

  h2->SetTitle("");
  h2->SetLineColor(mycolor);
  h2->SetLineWidth(2);
  h2->Draw("hist E same");

  TLegend leg2 (0.15,0.65,0.45,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg2.AddEntry(h2,legEntry2.c_str(),"L");

  leg2.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.,4.);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  frame->GetXaxis()->SetTitleSize(0.05);
  
  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  TH1D* den = (TH1D*) h2->Clone("den");
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++) den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(mycolor);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());

  if (yAxisName == "a.u.") h1->GetYaxis()->SetRangeUser(max(0.0001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  else h1->GetYaxis()->SetRangeUser(max(0.001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  canvas->SetLogy(0);

  delete canvas;
  frame->Reset("ICES");

  TFile fileOut(outputFILE.c_str(), "UPDATE");
  fileOut.cd();
  ratio->Write( ("ratio_" + canvasName).c_str());
  h2->Write( ("MC_" + canvasName).c_str());
  fileOut.Close();
}

void mcVsSplots(int wantPFPF, int doOtto)
{
  // Input files (after TnP selection)
  TFile *fileMC = new TFile("../files_first/mcDistributionsWithTnP_forComparisonWithSPlots.root");
  TFile *fileSPlots;
  if (wantPFPF==1) fileSPlots = new TFile("../files_first/myFileSPlots_PFPF.root");
  if (wantPFPF==0 && doOtto==1) fileSPlots = new TFile("../files_first/myFileSPlots_PFLPT_Otto.root");
  if (wantPFPF==0 && doOtto==0) fileSPlots = new TFile("../files_first/myFileSPlots_PFLPT_George.root");
  
  // MC histos 
  TH1F *mc_bdtO, *mc_bdtG;
  if (wantPFPF==1) mc_bdtO = (TH1F*)fileMC->Get("bdtOHsPFPF");
  if (wantPFPF==1) mc_bdtG = (TH1F*)fileMC->Get("bdtGHsPFPF");
  if (wantPFPF==0) mc_bdtO = (TH1F*)fileMC->Get("bdtOHsPFLPT");
  if (wantPFPF==0) mc_bdtG = (TH1F*)fileMC->Get("bdtGHsPFLPT");

  // s-Plots output 
  TH1F *data_bdtO = (TH1F*)fileSPlots->Get("h1_theAnalysisBdtO_B__theAnalysisBdtO");
  TH1F *data_bdtG = (TH1F*)fileSPlots->Get("h1_theAnalysisBdtG_B__theAnalysisBdtG");

  // Rebin
  mc_bdtO->Rebin();
  mc_bdtG->Rebin();
  data_bdtO->Rebin();
  data_bdtG->Rebin();

  // In sPlots, put bins with <0 weight to zero
  int nBinsData = data_bdtO->GetNbinsX();
  for (int ii=1; ii<=nBinsData; ii++) {
    if (data_bdtO->GetBinContent(ii)<0) data_bdtO->SetBinContent(ii,0);
    if (data_bdtG->GetBinContent(ii)<0) data_bdtG->SetBinContent(ii,0);
  }

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 

  if (wantPFPF==0 && doOtto==1) drawTH1pair(data_bdtO, mc_bdtO,  "BDT [Otto]",  "a.u.",1.,"ottoBdt","./",2,"Data","MC");
  if (wantPFPF==0 && doOtto==0) drawTH1pair(data_bdtG, mc_bdtG,  "BDT [George]","a.u.",1.,"georgeBdt","./",2,"Data","MC");
  if (wantPFPF==1) {
    drawTH1pair(data_bdtO, mc_bdtO,  "BDT [Otto]",  "a.u.",1.,"ottoBdt","./",2,"Data","MC");
    drawTH1pair(data_bdtG, mc_bdtG,  "BDT [George]","a.u.",1.,"georgeBdt","./",2,"Data","MC");
  }
}
