#define mcPlots_cxx
#include "mcPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// To be run on MC tnp formatted ntuples to get distributions for 
// - comparison with sPlots (checkWrtData=1)
// - comparison between distribution pre/post tnp selection
// Select signal and background based on match with MC-truth

void mcPlots::Loop(int checkWrtData)
{
  if (fChain == 0) return;

  if (checkWrtData==1) cout << "Check with data" << endl;
  else cout << "Check with MC" << endl;

  // Analysis BDT distributions to final comparisons
  TH1F *bdtOHs = new TH1F("bdtOHs","bdtOHs", 30, -15., 15.);
  TH1F *bdtGHs = new TH1F("bdtGHs","bdtGHs", 30, -15., 15.);
  bdtOHs->Sumw2();
  bdtGHs->Sumw2();
  TH1F *bdtOHb = new TH1F("bdtOHb","bdtOHb", 30, -15., 15.);
  TH1F *bdtGHb = new TH1F("bdtGHb","bdtGHb", 30, -15., 15.);
  bdtOHb->Sumw2();
  bdtGHb->Sumw2();

  // Analysis BDT distributions to final comparisons with sPlots only
  TH1F *bdtOHsPFPF = new TH1F("bdtOHsPFPF","bdtOHsPFPF", 30, -15., 15.);
  TH1F *bdtGHsPFPF = new TH1F("bdtGHsPFPF","bdtGHsPFPF", 30, -15., 15.);
  bdtOHsPFPF->Sumw2();
  bdtGHsPFPF->Sumw2();
  TH1F *bdtOHsPFLPT = new TH1F("bdtOHsPFLPT","bdtOHsPFLPT", 30, -15., 15.);
  TH1F *bdtGHsPFLPT = new TH1F("bdtGHsPFLPT","bdtGHsPFLPT", 30, -15., 15.);
  bdtOHsPFLPT->Sumw2();
  bdtGHsPFLPT->Sumw2();

  // Other distributions
  TH1F *tagPfIdHs   = new TH1F("tagPfIdHs",  "tagPfIdHs",  40,-10,10);
  TH1F *probePfIdHs = new TH1F("probePfIdHs","probePfIdHs",40,-10,10);
  TH1F *tagIdHs     = new TH1F("tagIdHs",    "tagIdHs",    40,-20,20);
  TH1F *probeIdHs   = new TH1F("probeIdHs",  "probeIdHs",  40,-20,20);
  tagPfIdHs->Sumw2();      
  probePfIdHs->Sumw2();      
  tagIdHs->Sumw2();      
  probeIdHs->Sumw2();      
  TH1F *tagPfIdHb   = new TH1F("tagPfIdHb",  "tagPfIdHb",  40,-10,10);
  TH1F *probePfIdHb = new TH1F("probePfIdHb","probePfIdHb",40,-10,10);
  TH1F *tagIdHb     = new TH1F("tagIdHb",    "tagIdHb",    40,-20,20);
  TH1F *probeIdHb   = new TH1F("probeIdHb",  "probeIdHb",  40,-20,20);
  tagPfIdHb->Sumw2();      
  probePfIdHb->Sumw2();      
  tagIdHb->Sumw2();      
  probeIdHb->Sumw2();      
  // 
  TH1F *tagPtHs   = new TH1F("tagPtHs",  "tagPtHs",  20,0.,10.);
  TH1F *probePtHs = new TH1F("probePtHs","probePtHs",20,0.,10.);
  TH1F *kPtHs     = new TH1F("kPtHs",    "kPtHs",    20,0.,10.);
  tagPtHs->Sumw2();
  probePtHs->Sumw2();
  kPtHs->Sumw2();
  TH1F *tagPtHb   = new TH1F("tagPtHb",  "tagPtHb",  20,0.,10.);
  TH1F *probePtHb = new TH1F("probePtHb","probePtHb",20,0.,10.);
  TH1F *kPtHb     = new TH1F("kPtHb",    "kPtHb",    20,0.,10.);
  tagPtHb->Sumw2();
  probePtHb->Sumw2();
  kPtHb->Sumw2();

  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // HLT applied in all cases at dumper level

    // JPsi range 
    if (pair_mass<2.9 || pair_mass>3.2) continue;

    // Bmass range
    if (B_mass<4.5 || B_mass>6) continue;

    // Safety cuts
    if (tagMvaId>20.5)     continue;
    if (probeMvaId>20.5)   continue;
    if (tagPfmvaId>20.5)   continue;
    if (probePfmvaId>20.5) continue;
    // 
    if (theAnalysisBdtO>15)  continue;
    if (theAnalysisBdtO<-15) continue;
    if (theAnalysisBdtG>15)  continue;
    if (theAnalysisBdtG<-15) continue;

    // Fill distributions for real Bs
    if (B_matchMC==1) {

      bdtOHs -> Fill(theAnalysisBdtO);
      bdtGHs -> Fill(theAnalysisBdtG);
      
      tagPtHs->Fill(tagPt);
      probePtHs->Fill(probePt);
      kPtHs->Fill(kPt);
      
      if (tagPfmvaId<20)   tagPfIdHs->Fill(tagPfmvaId);
      if (probePfmvaId<20) probePfIdHs->Fill(probePfmvaId);
      if (tagMvaId<20)     tagIdHs->Fill(tagMvaId);
      if (probeMvaId<20)   probeIdHs->Fill(probeMvaId);

      // Only for comparison with sPlots
      if (checkWrtData==1) { 

	if (pair_mass<3.0 || pair_mass>3.2) continue;

	if (probePfmvaId<20 && tagPfmvaId<20) { // PFPF
	  bdtOHsPFPF -> Fill(theAnalysisBdtO);
	  bdtGHsPFPF -> Fill(theAnalysisBdtG);

	} else {   // PFLPT
	  if (pair_mass<3.05 || pair_mass>3.15) continue;
	  if (theAnalysisBdtO>0) bdtOHsPFLPT -> Fill(theAnalysisBdtO);
	  if (theAnalysisBdtG>-4.5) bdtGHsPFLPT -> Fill(theAnalysisBdtG);
	}

      }


    } else {

      bdtOHb -> Fill(theAnalysisBdtO);
      bdtGHb -> Fill(theAnalysisBdtG);
      
      tagPtHb->Fill(tagPt);
      probePtHb->Fill(probePt);
      kPtHb->Fill(kPt);
      
      if (tagPfmvaId<20)   tagPfIdHb->Fill(tagPfmvaId);
      if (probePfmvaId<20) probePfIdHb->Fill(probePfmvaId);
      if (tagMvaId<20)     tagIdHb->Fill(tagMvaId);
      if (probeMvaId<20)   probeIdHb->Fill(probeMvaId);
    }     

  } // Loop over entries
  
  
  // Cosmetics
  bdtOHs->SetLineWidth(2);
  bdtOHs->SetLineColor(6);
  bdtGHs->SetLineWidth(2);
  bdtGHs->SetLineColor(6);
  bdtOHb->SetLineWidth(2);
  bdtOHb->SetLineColor(4);
  bdtGHb->SetLineWidth(2);
  bdtGHb->SetLineColor(4);

  tagPfIdHs->SetLineWidth(2);   
  tagPfIdHs->SetLineColor(6);  
  probePfIdHs->SetLineWidth(2);   
  probePfIdHs->SetLineColor(6);  
  tagIdHs->SetLineWidth(2);   
  tagIdHs->SetLineColor(6); 
  probeIdHs->SetLineWidth(2);   
  probeIdHs->SetLineColor(6);  
  tagPfIdHb->SetLineWidth(2);   
  tagPfIdHb->SetLineColor(4);  
  probePfIdHb->SetLineWidth(2);   
  probePfIdHb->SetLineColor(4);  
  tagIdHb->SetLineWidth(2);   
  tagIdHb->SetLineColor(4); 
  probeIdHb->SetLineWidth(2);   
  probeIdHb->SetLineColor(4);  

  tagPtHs->SetLineWidth(2);  
  tagPtHs->SetLineColor(6);    
  probePtHs->SetLineWidth(2);  
  probePtHs->SetLineColor(6);    
  kPtHs->SetLineWidth(2);  
  kPtHs->SetLineColor(6);    
  tagPtHb->SetLineWidth(2);  
  tagPtHb->SetLineColor(4);    
  probePtHb->SetLineWidth(2);  
  probePtHb->SetLineColor(4);    
  kPtHb->SetLineWidth(2);  
  kPtHb->SetLineColor(4);    

  // Save output
  TFile myfile("myFileMC.root","RECREATE");
  bdtOHs->Write();
  bdtGHs->Write();
  if (checkWrtData==1) {
    bdtOHsPFLPT->Write();
    bdtGHsPFLPT->Write();
    bdtOHsPFPF->Write();
    bdtGHsPFPF->Write();
  }
  tagIdHs->Write();  
  probeIdHs->Write();   
  probePfIdHs->Write();   
  tagPfIdHs->Write();    
  tagPtHs->Write(); 
  probePtHs->Write();  
  kPtHs->Write(); 
  myfile.Close();

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.10,0.65,0.45,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(bdtOHs, "Matching", "lp");
  leg->AddEntry(bdtOHb, "Not matching", "lp");
  //
  TLegend *leg2;
  leg2 = new TLegend(0.50,0.65,0.95,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(bdtOHs, "Matching", "lp");
  leg2->AddEntry(bdtOHb, "Not matching", "lp");


  TCanvas c1("c1","c",1);
  bdtOHs->SetTitle("");
  bdtOHb->SetTitle("");
  bdtOHs->GetXaxis()->SetTitle("BDT [Otto]");
  bdtOHb->GetXaxis()->SetTitle("BDT [Otto]");
  bdtOHs->DrawNormalized("hist");
  bdtOHb->DrawNormalized("samehist");
  leg->Draw();
  c1.SaveAs("ottoBdt.png");

  TCanvas c2("c2","c",1);
  bdtGHs->SetTitle("");
  bdtGHb->SetTitle("");
  bdtGHs->GetXaxis()->SetTitle("BDT [George]");
  bdtGHb->GetXaxis()->SetTitle("BDT [George]");
  bdtGHs->DrawNormalized("hist");
  bdtGHb->DrawNormalized("samehist");
  leg->Draw();
  c2.SaveAs("georgeBdt.png");

  TCanvas c3("c3","c",1);
  tagIdHs->SetTitle("");
  tagIdHb->SetTitle("");
  tagIdHs->GetXaxis()->SetTitle("tag eleID [LPT]");
  tagIdHb->GetXaxis()->SetTitle("tag eleID [LPT]");
  tagIdHb->DrawNormalized("hist");
  tagIdHs->DrawNormalized("samehist");
  leg->Draw();
  c3.SaveAs("tagLptId.png");

  TCanvas c4("c4","c",1);
  probeIdHs->SetTitle("");
  probeIdHb->SetTitle("");
  probeIdHs->GetXaxis()->SetTitle("probe eleID [LPT]");
  probeIdHb->GetXaxis()->SetTitle("probe eleID [LPT]");
  probeIdHb->DrawNormalized("hist");
  probeIdHs->DrawNormalized("samehist");
  leg->Draw();
  c4.SaveAs("probeLptId.png");

  TCanvas c3a("c3a","c",1);
  tagPfIdHs->SetTitle("");
  tagPfIdHb->SetTitle("");
  tagPfIdHs->GetXaxis()->SetTitle("tag eleID [PF]");
  tagPfIdHb->GetXaxis()->SetTitle("tag eleID [PF]");
  tagPfIdHb->DrawNormalized("hist");
  tagPfIdHs->DrawNormalized("samehist");
  leg->Draw();
  c3a.SaveAs("tagPFId.png");

  TCanvas c4a("c4a","c",1);
  probePfIdHs->SetTitle("");
  probePfIdHb->SetTitle("");
  probePfIdHs->GetXaxis()->SetTitle("probe eleID [PF]");
  probePfIdHb->GetXaxis()->SetTitle("probe eleID [PF]");
  probePfIdHb->DrawNormalized("hist");
  probePfIdHs->DrawNormalized("samehist");
  leg->Draw();
  c4a.SaveAs("probePfId.png");

  TCanvas c5("c5","c",1);
  tagPtHs->SetTitle("");
  tagPtHb->SetTitle("");
  tagPtHs->GetXaxis()->SetTitle("tag pt");
  tagPtHb->GetXaxis()->SetTitle("tag pt");
  tagPtHb->DrawNormalized("hist");
  tagPtHs->DrawNormalized("samehist");
  leg2->Draw();
  c5.SaveAs("tagPt.png");

  TCanvas c6("c6","c",1);
  probePtHs->SetTitle("");
  probePtHb->SetTitle("");
  probePtHs->GetXaxis()->SetTitle("probe pt");
  probePtHb->GetXaxis()->SetTitle("probe pt");
  probePtHb->DrawNormalized("hist");
  probePtHs->DrawNormalized("samehist");
  leg2->Draw();
  c6.SaveAs("probePt.png");

  TCanvas c7("c7","c",1);
  kPtHs->SetTitle("");
  kPtHb->SetTitle("");
  kPtHs->GetXaxis()->SetTitle("k pt");
  kPtHb->GetXaxis()->SetTitle("k pt");
  kPtHb->DrawNormalized("hist");
  kPtHs->DrawNormalized("samehist");
  leg2->Draw();
  c7.SaveAs("kPt.png");
}
