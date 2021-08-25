#define nanoVsTnp_cxx

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>

using namespace std;

// To check the effect of tnp selection on signal distributions in MC
// To be run on the outcome of mcPlots

void nanoVsTnp()
{
  // Files
  TFile fileTnP("../files_first/mcDistributionsWithTnP.root");
  TFile fileNoTnP("../files_first/mcDistributionsNoTnP.root");

  // Histos
  TH1F *bdtOFromTnP    = (TH1F*)fileTnP.Get("bdtOHs");
  TH1F *bdtGFromTnP    = (TH1F*)fileTnP.Get("bdtGHs");
  TH1F *tagPfIdFromTnP = (TH1F*)fileTnP.Get("tagPfIdHs");
  TH1F *probePfIdFromTnP = (TH1F*)fileTnP.Get("probePfIdHs");
  TH1F *tagIdFromTnP   = (TH1F*)fileTnP.Get("tagIdHs");
  TH1F *probeIdFromTnP = (TH1F*)fileTnP.Get("probeIdHs");
  TH1F *tagPtFromTnP   = (TH1F*)fileTnP.Get("tagPtHs");
  TH1F *probePtFromTnP = (TH1F*)fileTnP.Get("probePtHs");
  TH1F *kPtFromTnP     = (TH1F*)fileTnP.Get("kPtHs");
  bdtOFromTnP-> Scale(1./bdtOFromTnP->Integral());   
  bdtGFromTnP-> Scale(1./bdtGFromTnP->Integral());   
  tagPfIdFromTnP-> Scale(1./tagPfIdFromTnP->Integral());   
  probePfIdFromTnP-> Scale(1./probePfIdFromTnP->Integral());   
  tagIdFromTnP-> Scale(1./tagIdFromTnP->Integral());   
  probeIdFromTnP-> Scale(1./probeIdFromTnP->Integral());   
  tagPtFromTnP-> Scale(1./tagPtFromTnP->Integral());   
  probePtFromTnP-> Scale(1./probePtFromTnP->Integral());   
  kPtFromTnP-> Scale(1./kPtFromTnP->Integral());   

  TH1F *bdtONoTnP    = (TH1F*)fileNoTnP.Get("bdtOHs");
  TH1F *bdtGNoTnP    = (TH1F*)fileNoTnP.Get("bdtGHs");
  TH1F *tagPfIdNoTnP = (TH1F*)fileNoTnP.Get("tagPfIdHs");
  TH1F *probePfIdNoTnP = (TH1F*)fileNoTnP.Get("probePfIdHs");
  TH1F *tagIdNoTnP   = (TH1F*)fileNoTnP.Get("tagIdHs");
  TH1F *probeIdNoTnP = (TH1F*)fileNoTnP.Get("probeIdHs");
  TH1F *tagPtNoTnP   = (TH1F*)fileNoTnP.Get("tagPtHs");
  TH1F *probePtNoTnP = (TH1F*)fileNoTnP.Get("probePtHs");
  TH1F *kPtNoTnP     = (TH1F*)fileNoTnP.Get("kPtHs");
  bdtONoTnP-> Scale(1./bdtONoTnP->Integral());   
  bdtGNoTnP-> Scale(1./bdtGNoTnP->Integral());   
  tagPfIdNoTnP-> Scale(1./tagPfIdNoTnP->Integral());   
  probePfIdNoTnP-> Scale(1./probePfIdNoTnP->Integral());   
  tagIdNoTnP-> Scale(1./tagIdNoTnP->Integral());   
  probeIdNoTnP-> Scale(1./probeIdNoTnP->Integral());   
  tagPtNoTnP-> Scale(1./tagPtNoTnP->Integral());   
  probePtNoTnP-> Scale(1./probePtNoTnP->Integral());   
  kPtNoTnP-> Scale(1./kPtNoTnP->Integral());   


  // --------------------------------------------------------
  // Cosmetics
  bdtOFromTnP -> SetLineWidth(2);  
  bdtONoTnP   -> SetLineWidth(2);  
  bdtOFromTnP -> SetLineColor(2);  
  bdtONoTnP   -> SetLineColor(4);  
  bdtOFromTnP -> SetTitle("");
  bdtONoTnP   -> SetTitle("");
  bdtOFromTnP -> GetXaxis()->SetTitle("BDT [Otto]");
  bdtONoTnP   -> GetXaxis()->SetTitle("BDT [Otto]");
  //
  bdtGFromTnP -> SetLineWidth(2);  
  bdtGNoTnP   -> SetLineWidth(2);  
  bdtGFromTnP -> SetLineColor(2);  
  bdtGNoTnP   -> SetLineColor(4);  
  bdtGFromTnP -> SetTitle("");
  bdtGNoTnP   -> SetTitle("");
  bdtGFromTnP -> GetXaxis()->SetTitle("BDT [George]");
  bdtGNoTnP   -> GetXaxis()->SetTitle("BDT [George]");
  //
  tagPfIdFromTnP -> SetLineWidth(2);  
  tagPfIdNoTnP   -> SetLineWidth(2);  
  tagPfIdFromTnP -> SetLineColor(2);  
  tagPfIdNoTnP   -> SetLineColor(4);  
  tagPfIdFromTnP -> SetTitle("");
  tagPfIdNoTnP   -> SetTitle("");
  tagPfIdFromTnP -> GetXaxis()->SetTitle("tag eleID [PF]");
  tagPfIdNoTnP   -> GetXaxis()->SetTitle("tag eleID [PF]");
  //
  probePfIdFromTnP -> SetLineWidth(2);  
  probePfIdNoTnP   -> SetLineWidth(2);  
  probePfIdFromTnP -> SetLineColor(2);  
  probePfIdNoTnP   -> SetLineColor(4);  
  probePfIdFromTnP -> SetTitle("");
  probePfIdNoTnP   -> SetTitle("");
  probePfIdFromTnP -> GetXaxis()->SetTitle("probe eleID [PF]");
  probePfIdNoTnP   -> GetXaxis()->SetTitle("probe eleID [PF]");
  //
  tagIdFromTnP -> SetLineWidth(2);  
  tagIdNoTnP   -> SetLineWidth(2);  
  tagIdFromTnP -> SetLineColor(2);  
  tagIdNoTnP   -> SetLineColor(4);  
  tagIdFromTnP -> SetTitle("");
  tagIdNoTnP   -> SetTitle("");
  tagIdFromTnP -> GetXaxis()->SetTitle("tag eleID [LPT]");
  tagIdNoTnP   -> GetXaxis()->SetTitle("tag eleID [LPT]");
  //
  probeIdFromTnP -> SetLineWidth(2);  
  probeIdNoTnP   -> SetLineWidth(2);  
  probeIdFromTnP -> SetLineColor(2);  
  probeIdNoTnP   -> SetLineColor(4);  
  probeIdFromTnP -> SetTitle("");
  probeIdNoTnP   -> SetTitle("");
  probeIdFromTnP -> GetXaxis()->SetTitle("probe eleID [LPT]");
  probeIdNoTnP   -> GetXaxis()->SetTitle("probe eleID [LPT]");
  //
  tagPtFromTnP -> SetLineWidth(2);  
  tagPtNoTnP   -> SetLineWidth(2);  
  tagPtFromTnP -> SetLineColor(2);  
  tagPtNoTnP   -> SetLineColor(4);  
  tagPtFromTnP -> SetTitle("");
  tagPtNoTnP   -> SetTitle("");
  tagPtFromTnP -> GetXaxis()->SetTitle("tag pT [GeV]");
  tagPtNoTnP   -> GetXaxis()->SetTitle("tag pT [GeV]");
  //
  probePtFromTnP -> SetLineWidth(2);  
  probePtNoTnP   -> SetLineWidth(2);  
  probePtFromTnP -> SetLineColor(2);  
  probePtNoTnP   -> SetLineColor(4);  
  probePtFromTnP -> SetTitle("");
  probePtNoTnP   -> SetTitle("");
  probePtFromTnP -> GetXaxis()->SetTitle("probe pT [GeV]");
  probePtNoTnP   -> GetXaxis()->SetTitle("probe pT [GeV]");
  //
  kPtFromTnP -> SetLineWidth(2);  
  kPtNoTnP   -> SetLineWidth(2);  
  kPtFromTnP -> SetLineColor(2);  
  kPtNoTnP   -> SetLineColor(4);  
  kPtFromTnP -> SetTitle("");
  kPtNoTnP   -> SetTitle("");
  kPtFromTnP -> GetXaxis()->SetTitle("k pT [GeV]");
  kPtNoTnP   -> GetXaxis()->SetTitle("k pT [GeV]");


  // -------------------------------------------------------------
  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.15,0.65,0.60,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(kPtFromTnP, "With tnp", "lp");
  leg->AddEntry(kPtNoTnP,   "No tnp", "lp");

  TCanvas c1("c1","c1",1);
  bdtOFromTnP->DrawNormalized("hist");
  bdtONoTnP->DrawNormalized("samehist");
  leg->Draw();
  c1.SaveAs("ottoBdt_naniVsTnp.png");

  TCanvas c2("c2","c2",1);
  bdtGFromTnP->DrawNormalized("hist");
  bdtGNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c2.SaveAs("georgeBdt_naniVsTnp.png");

  TCanvas c3("c3","c3",1);
  tagPfIdFromTnP->DrawNormalized("hist");
  tagPfIdNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c3.SaveAs("tagPfId_naniVsTnp.png");

  TCanvas c4("c4","c4",1);
  probePfIdFromTnP->DrawNormalized("hist");
  probePfIdNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c4.SaveAs("probePfId_naniVsTnp.png");

  TCanvas c3a("c3a","c3a",1);
  tagIdNoTnP->DrawNormalized("hist");
  tagIdFromTnP->DrawNormalized("samehist");
  leg->Draw();
  c3a.SaveAs("tagId_naniVsTnp.png");

  TCanvas c4a("c4a","c4a",1);
  probeIdFromTnP->DrawNormalized("hist");
  probeIdNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c4a.SaveAs("probeId_naniVsTnp.png");

  TCanvas c5("c5","c5",1);
  tagPtNoTnP->DrawNormalized("hist");
  tagPtFromTnP->DrawNormalized("samehist");
  leg->Draw();
  c5.SaveAs("tagPt_naniVsTnp.png");

  TCanvas c6("c6","c6",1);
  probePtFromTnP->DrawNormalized("hist");
  probePtNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c6.SaveAs("probePt_naniVsTnp.png");

  TCanvas c7("c7","c7",1);
  kPtFromTnP->DrawNormalized("hist");
  kPtNoTnP->DrawNormalized("samehist");
  leg->Draw();
  c7.SaveAs("kPt_naniVsTnp.png");
}
