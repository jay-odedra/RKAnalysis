#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>
#include "TLorentzVector.h"

using namespace std;

void tnpForBTreeFormat(const char* filename) {

  cout << "Formatting " << filename << endl;  

  // Original ntuple   
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(TString("/tmp/crovelli/")+TString(filename));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("TaPtree");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }
  
  fileOrig->cd();
  if (!treeOrig) {
    cout << "TaPTree not existing !" << endl; 
    return;    
  }

  treeOrig->SetMakeClass(0);
  cout << "TreeOrig->Size = "<< treeOrig->GetEntries() << endl;

  // number of entries saved in the first tree
  int nentriesOrig = treeOrig->GetEntries();   

  // Tree for the final format
  TFile *fileNew = TFile::Open(TString("/tmp/crovelli/FormattedTnPForB_")+TString(filename),"recreate");
  fileNew->ls();
  fileNew->cd();
  TDirectory *myDir = (TDirectory*)fileNew->mkdir("tnpAna");
  myDir->cd();
  TTree *treeNew = new TTree("fitter_tree","reduced tree for T&P with B");
  treeNew->SetAutoSave(-99999999999);
  treeNew->SetAutoFlush(-99999999999);

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  // original tree leaves
  Int_t           theRun = 0;
  Int_t           theEvent = 0;
  Int_t           nvtx = 0;
  Int_t           sampleID = 0;
  Float_t         rho = 0;
  Float_t         pu_weight = 0;
  Float_t         perEveW = 0;
  Int_t           hlt9 = 0;
  Float_t         trg_muon_pt = 0;
  Int_t           selectedBSize = 0;
  vector<float>   *tag_pt = nullptr;
  vector<float>   *tag_eta = nullptr;
  vector<float>   *tag_phi = nullptr;
  vector<bool>    *tag_isPF = nullptr;
  vector<bool>    *tag_isPFOverlap = nullptr;
  vector<bool>    *tag_isLowPt = nullptr;
  vector<float>   *tag_mvaId = nullptr;
  vector<float>   *tag_pfmvaId = nullptr;
  vector<int>     *tag_convveto = nullptr;
  vector<float>   *mll_fullfit = nullptr;
  vector<float>   *mll_raw = nullptr; 
  vector<float>   *fit_mass = nullptr; 
  vector<float>   *fit_Bpt = nullptr; 
  vector<float>   *fit_Bcos2D = nullptr; 
  vector<float>   *fit_Bsvprob = nullptr; 
  vector<float>   *fit_Bxysig = nullptr; 
  vector<bool>    *bmatchMC = nullptr; 
  vector<float>   *K_pt = nullptr; 
  vector<float>   *K_eta = nullptr; 
  vector<float>   *K_phi = nullptr; 
  vector<float>   *probe_pt = nullptr;
  vector<float>   *probe_eta = nullptr;
  vector<float>   *probe_phi = nullptr;
  vector<bool>    *probe_isPF = nullptr;
  vector<bool>    *probe_isPFOverlap = nullptr;
  vector<bool>    *probe_isLowPt = nullptr;
  vector<float>   *probe_mvaId = nullptr;
  vector<float>   *probe_pfmvaId = nullptr;
  vector<int>     *probe_convveto = nullptr;
  vector<float>   *probe_drm = nullptr;
  vector<float>   *tag_drm = nullptr;
  vector<float>   *bPtOrder = nullptr;  
  vector<float>   *b2DOrder = nullptr;  
  vector<float>   *bXYOrder = nullptr;  
  vector<float>   *bMassOrder = nullptr;  
  vector<float>   *p4Trk = nullptr;  
  vector<float>   *analysisBdtG = nullptr;  
  vector<float>   *analysisBdtO = nullptr;  
  vector<float>   *LKdz = nullptr;  
  vector<float>   *L1L2dr = nullptr;  
  vector<float>   *LKdr = nullptr;  
  vector<float>   *Kiso = nullptr;  
  vector<float>   *BBDPhi = nullptr;  
  vector<float>   *BTrkdxy2 = nullptr;  
  vector<float>   *Dmass = nullptr;  
  vector<float>   *Dmass_flip = nullptr;  
  vector<float>   *Dmass_ll = nullptr;  
  vector<float>   *Dmass_ll_flip = nullptr;  
  vector<float>   *B_iso04_rel = nullptr;  
  vector<float>   *eleDR = nullptr;  
  vector<float>   *k_DCASig = nullptr;  
  vector<float>   *k_dzTrg = nullptr;  
  vector<float>   *k_dxy_sig = nullptr;  
  vector<float>   *k_iso04_rel = nullptr;  
  vector<float>   *k_svip2d = nullptr;  
  vector<float>   *k_svip3d = nullptr;  
  vector<float>   *tag_dxy_sig = nullptr;  
  vector<float>   *tag_dzTrg = nullptr;  
  vector<float>   *tag_iso04_rel = nullptr;  
  vector<float>   *probe_dxy_sig = nullptr;  
  vector<float>   *probe_dzTrg = nullptr;  
  vector<float>   *probe_iso04_rel = nullptr;  
  vector<float>   *llkDR = nullptr;  
  vector<float>   *ptAsym = nullptr;  
  Int_t           selectedPairsSize = 0;  
  
  // List of branches - original tree
  TBranch        *b_theRun;   //!
  TBranch        *b_theEvent;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_sampleID;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_hlt9;
  TBranch        *b_trg_muon_pt;   //! 
  TBranch        *b_selectedBSize;   //!
  TBranch        *b_tag_pt;   //!
  TBranch        *b_tag_eta;   //!
  TBranch        *b_tag_phi;   //!
  TBranch        *b_tag_isPF;   //!
  TBranch        *b_tag_isPFOverlap;   //!
  TBranch        *b_tag_isLowPt;   //!
  TBranch        *b_tag_mvaId;   //!
  TBranch        *b_tag_pfmvaId;   //!
  TBranch        *b_tag_convveto;   //!
  TBranch        *b_mll_fullfit;   //!  
  TBranch        *b_mll_raw;   //!
  TBranch        *b_fit_mass;   //!
  TBranch        *b_fit_Bpt;   //!
  TBranch        *b_fit_Bcos2D;   //!
  TBranch        *b_fit_Bsvprob;   //!
  TBranch        *b_fit_Bxysig;   //!
  TBranch        *b_bmatchMC;   //!
  TBranch        *b_K_pt;   //!
  TBranch        *b_K_eta;   //!
  TBranch        *b_K_phi;   //!
  TBranch        *b_probe_pt;   //!
  TBranch        *b_probe_eta;   //!
  TBranch        *b_probe_phi;   //!
  TBranch        *b_probe_isPF;   //!
  TBranch        *b_probe_isPFOverlap;   //!
  TBranch        *b_probe_isLowPt;   //!
  TBranch        *b_probe_mvaId;   //!
  TBranch        *b_probe_pfmvaId;   //!
  TBranch        *b_probe_convveto;   //!
  TBranch        *b_probe_drm;   //!
  TBranch        *b_tag_drm;   //!
  TBranch        *b_bPtOrder;   //!
  TBranch        *b_b2DOrder;   //!
  TBranch        *b_bXYOrder;   //!
  TBranch        *b_bMassOrder;   //!
  TBranch        *b_p4Trk;   //!
  TBranch        *b_analysisBdtG;   //!
  TBranch        *b_analysisBdtO;   //!
  TBranch        *b_LKdz;   //!
  TBranch        *b_L1L2dr;   //!
  TBranch        *b_LKdr;   //!
  TBranch        *b_Kiso;   //!
  TBranch        *b_BBDPhi;   //!
  TBranch        *b_BTrkdxy2;   //!
  TBranch        *b_Dmass;   //!
  TBranch        *b_Dmass_flip;   //!
  TBranch        *b_Dmass_ll;   //!
  TBranch        *b_Dmass_ll_flip;   //!
  TBranch        *b_B_iso04_rel;   //!
  TBranch        *b_eleDR;   //!
  TBranch        *b_k_DCASig;   //!
  TBranch        *b_k_dzTrg;   //!
  TBranch        *b_k_dxy_sig;   //!
  TBranch        *b_k_iso04_rel;   //!
  TBranch        *b_k_svip2d;   //!
  TBranch        *b_k_svip3d;   //!
  TBranch        *b_tag_dxy_sig;   //!
  TBranch        *b_tag_dzTrg;   //!
  TBranch        *b_tag_iso04_rel;   //!
  TBranch        *b_probe_dxy_sig;   //!
  TBranch        *b_probe_dzTrg;   //!
  TBranch        *b_probe_iso04_rel;   //!
  TBranch        *b_llkDR;   //!
  TBranch        *b_ptAsym;   //!
  TBranch        *b_selectedPairsSize;   //!
  
  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("theRun", &theRun, &b_theRun);
  treeOrig->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
  treeOrig->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  treeOrig->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
  treeOrig->SetBranchAddress("rho", &rho, &b_rho);
  treeOrig->SetBranchAddress("hlt9", &hlt9, &b_hlt9);
  treeOrig->SetBranchAddress("selectedBSize", &selectedBSize, &b_selectedBSize);
  treeOrig->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
  treeOrig->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
  treeOrig->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
  treeOrig->SetBranchAddress("tag_isPF", &tag_isPF, &b_tag_isPF);
  treeOrig->SetBranchAddress("tag_isPFOverlap", &tag_isPFOverlap, &b_tag_isPFOverlap);
  treeOrig->SetBranchAddress("tag_isLowPt", &tag_isLowPt, &b_tag_isLowPt);
  treeOrig->SetBranchAddress("tag_mvaId", &tag_mvaId, &b_tag_mvaId);
  treeOrig->SetBranchAddress("tag_pfmvaId", &tag_pfmvaId, &b_tag_pfmvaId);
  treeOrig->SetBranchAddress("tag_convveto", &tag_convveto, &b_tag_convveto);
  treeOrig->SetBranchAddress("mll_fullfit", &mll_fullfit, &b_mll_fullfit);
  treeOrig->SetBranchAddress("mll_raw", &mll_raw, &b_mll_raw);
  treeOrig->SetBranchAddress("fit_mass", &fit_mass, &b_fit_mass);
  treeOrig->SetBranchAddress("fit_Bpt", &fit_Bpt, &b_fit_Bpt);
  treeOrig->SetBranchAddress("fit_Bcos2D", &fit_Bcos2D, &b_fit_Bcos2D);
  treeOrig->SetBranchAddress("fit_Bsvprob", &fit_Bsvprob, &b_fit_Bsvprob);
  treeOrig->SetBranchAddress("fit_Bxysig", &fit_Bxysig, &b_fit_Bxysig);
  treeOrig->SetBranchAddress("bmatchMC", &bmatchMC, &b_bmatchMC);
  treeOrig->SetBranchAddress("K_pt", &K_pt, &b_K_pt);
  treeOrig->SetBranchAddress("K_eta", &K_eta, &b_K_eta);
  treeOrig->SetBranchAddress("K_phi", &K_phi, &b_K_phi);
  treeOrig->SetBranchAddress("probe_pt", &probe_pt, &b_probe_pt);
  treeOrig->SetBranchAddress("probe_eta", &probe_eta, &b_probe_eta);
  treeOrig->SetBranchAddress("probe_phi", &probe_phi, &b_probe_phi);
  treeOrig->SetBranchAddress("probe_isPF", &probe_isPF, &b_probe_isPF);
  treeOrig->SetBranchAddress("probe_isPFOverlap", &probe_isPFOverlap, &b_probe_isPFOverlap);
  treeOrig->SetBranchAddress("probe_isLowPt", &probe_isLowPt, &b_probe_isLowPt);
  treeOrig->SetBranchAddress("probe_mvaId", &probe_mvaId, &b_probe_mvaId);
  treeOrig->SetBranchAddress("probe_pfmvaId", &probe_pfmvaId, &b_probe_pfmvaId);
  treeOrig->SetBranchAddress("probe_convveto", &probe_convveto, &b_probe_convveto);
  treeOrig->SetBranchAddress("probe_drm", &probe_drm, &b_probe_drm);
  treeOrig->SetBranchAddress("tag_drm", &tag_drm, &b_tag_drm);
  treeOrig->SetBranchAddress("bPtOrder", &bPtOrder, &b_bPtOrder);
  treeOrig->SetBranchAddress("b2DOrder", &b2DOrder, &b_b2DOrder);
  treeOrig->SetBranchAddress("bXYOrder", &bXYOrder, &b_bXYOrder);
  treeOrig->SetBranchAddress("bMassOrder", &bMassOrder, &b_bMassOrder);
  treeOrig->SetBranchAddress("p4Trk", &p4Trk, &b_p4Trk);
  treeOrig->SetBranchAddress("analysisBdtG", &analysisBdtG, &b_analysisBdtG);
  treeOrig->SetBranchAddress("analysisBdtO", &analysisBdtO, &b_analysisBdtO);
  treeOrig->SetBranchAddress("LKdz", &LKdz, &b_LKdz);
  treeOrig->SetBranchAddress("L1L2dr", &L1L2dr, &b_L1L2dr);
  treeOrig->SetBranchAddress("LKdr", &LKdr, &b_LKdr);
  treeOrig->SetBranchAddress("Kiso", &Kiso, &b_Kiso);
  treeOrig->SetBranchAddress("BBDPhi", &BBDPhi, &b_BBDPhi);
  treeOrig->SetBranchAddress("BTrkdxy2", &BTrkdxy2, &b_BTrkdxy2);
  treeOrig->SetBranchAddress("Dmass", &Dmass, &b_Dmass);
  treeOrig->SetBranchAddress("Dmass_flip", &Dmass_flip, &b_Dmass_flip);
  treeOrig->SetBranchAddress("Dmass_ll", &Dmass_ll, &b_Dmass_ll);
  treeOrig->SetBranchAddress("Dmass_ll_flip", &Dmass_ll_flip, &b_Dmass_ll_flip);
  treeOrig->SetBranchAddress("B_iso04_rel", &B_iso04_rel, &b_B_iso04_rel);
  treeOrig->SetBranchAddress("eleDR", &eleDR, &b_eleDR);
  treeOrig->SetBranchAddress("k_DCASig", &k_DCASig, &b_k_DCASig);
  treeOrig->SetBranchAddress("k_dzTrg", &k_dzTrg, &b_k_dzTrg);
  treeOrig->SetBranchAddress("k_dxy_sig", &k_dxy_sig, &b_k_dxy_sig);
  treeOrig->SetBranchAddress("k_iso04_rel", &k_iso04_rel, &b_k_iso04_rel);
  treeOrig->SetBranchAddress("k_svip2d", &k_svip2d, &b_k_svip2d);
  treeOrig->SetBranchAddress("k_svip3d", &k_svip3d, &b_k_svip3d);
  treeOrig->SetBranchAddress("tag_dxy_sig", &tag_dxy_sig, &b_tag_dxy_sig);
  treeOrig->SetBranchAddress("tag_dzTrg", &tag_dzTrg, &b_tag_dzTrg);
  treeOrig->SetBranchAddress("tag_iso04_rel", &tag_iso04_rel, &b_tag_iso04_rel);
  treeOrig->SetBranchAddress("probe_dxy_sig", &probe_dxy_sig, &b_probe_dxy_sig);
  treeOrig->SetBranchAddress("probe_dzTrg", &probe_dzTrg, &b_probe_dzTrg);
  treeOrig->SetBranchAddress("probe_iso04_rel", &probe_iso04_rel, &b_probe_iso04_rel);
  treeOrig->SetBranchAddress("llkDR", &llkDR, &b_llkDR);
  treeOrig->SetBranchAddress("ptAsym", &ptAsym, &b_ptAsym);
  treeOrig->SetBranchAddress("selectedPairsSize", &selectedPairsSize, &b_selectedPairsSize);
  
  // New variables
  Int_t     theNvtx;
  Float_t   tagPt;
  Float_t   tagEta;
  Float_t   tagPhi;
  Float_t   tagPfmvaId;
  Float_t   tagMvaId;
  Float_t   tagDxySig;
  Float_t   tagDzTrg;
  Float_t   tagIso04Rel;
  Float_t   probePt;
  Float_t   probeEta;
  Float_t   probePhi;
  Float_t   probePfmvaId;
  Float_t   probeMvaId;
  Float_t   probeDxySig;
  Float_t   probeDzTrg;
  Float_t   probeIso04Rel;
  Float_t   elesDr;
  Float_t   kPt;
  Float_t   kEta;
  Float_t   kPhi;
  Float_t   B_isoRel;
  Float_t   B_cos2d;
  Float_t   B_pt;
  Float_t   B_mass;
  Float_t   B_xysig;
  Float_t   B_svprob;
  Int_t     B_matchMC;
  Float_t   kDCASig;
  Float_t   kDzTrg;
  Float_t   kIso04Rel;
  Float_t   kSvip2d;
  Float_t   kSvip3d;
  Float_t   theLlkDR;
  Float_t   thePtAsym;
  Float_t   theLKdz;
  Float_t   theL1L2dr;
  Float_t   theLKdr;
  Float_t   theKiso;
  Float_t   theBBDPhi;
  Float_t   theBTrkdxy2;
  Float_t   theAnalysisBdtG;
  Float_t   theAnalysisBdtO;
  Float_t   pair_mass;
  Float_t   weight;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    
    theTreeNew->Branch("theNvtx", &theNvtx, "theNvtx/I");
    
    theTreeNew->Branch("tagPt",&tagPt,"tagPt/F");
    theTreeNew->Branch("tagEta",&tagEta,"tagEta/F");
    theTreeNew->Branch("tagPhi",&tagPhi,"tagPhi/F");
    theTreeNew->Branch("tagPfmvaId",&tagPfmvaId,"tagPfmvaId/F");
    theTreeNew->Branch("tagMvaId",&tagMvaId,"tagMvaId/F");
    theTreeNew->Branch("tagDxySig",&tagDxySig,"tagDxySig/F");
    theTreeNew->Branch("tagDzTrg",&tagDzTrg,"tagDzTrg/F");
    theTreeNew->Branch("tagIso04Rel",&tagIso04Rel,"tagIso04Rel/F");

    theTreeNew->Branch("probePt",&probePt,"probePt/F");
    theTreeNew->Branch("probeEta",&probeEta,"probeEta/F");
    theTreeNew->Branch("probePhi",&probePhi,"probePhi/F");
    theTreeNew->Branch("probePfmvaId",&probePfmvaId,"probePfmvaId/F");
    theTreeNew->Branch("probeMvaId",&probeMvaId,"probeMvaId/F");
    theTreeNew->Branch("probeDxySig",&probeDxySig,"probeDxySig/F");
    theTreeNew->Branch("probeDzTrg",&probeDzTrg,"probeDzTrg/F");
    theTreeNew->Branch("probeIso04Rel",&probeIso04Rel,"probeIso04Rel/F");

    theTreeNew->Branch("elesDr",&elesDr,"elesDr/F");

    theTreeNew->Branch("kPt",&kPt,"kPt/F");
    theTreeNew->Branch("kEta",&kEta,"kEta/F");
    theTreeNew->Branch("kPhi",&kPhi,"kPhi/F");
    theTreeNew->Branch("kDCASig",&kDCASig,"kDCASig/F");
    theTreeNew->Branch("kDzTrg",&kDzTrg,"kDzTrg/F");
    theTreeNew->Branch("kIso04Rel",&kIso04Rel,"kIso04Rel/F");
    theTreeNew->Branch("kSvip2d",&kSvip2d,"kSvip2d/F");
    theTreeNew->Branch("kSvip3d",&kSvip3d,"kSvip3d/F");

    theTreeNew->Branch("B_mass", &B_mass, "B_mass/F");
    theTreeNew->Branch("B_cos2d", &B_cos2d, "B_cos2d/F");
    theTreeNew->Branch("B_pt", &B_pt, "B_pt/F");
    theTreeNew->Branch("B_svprob", &B_svprob, "B_svprob/F");
    theTreeNew->Branch("B_xysig", &B_xysig, "B_xysig/F");
    theTreeNew->Branch("B_isoRel", &B_isoRel, "B_isoRel/F");
    theTreeNew->Branch("B_matchMC", &B_matchMC, "B_matchMC/I");

    theTreeNew->Branch("theLlkDR", &theLlkDR, "theLlkDR/F");
    theTreeNew->Branch("thePtAsym", &thePtAsym, "thePtAsym/F");
    theTreeNew->Branch("theLKdz", &theLKdz, "theLKdz/F");
    theTreeNew->Branch("theL1L2dr", &theL1L2dr, "theL1L2dr/F");
    theTreeNew->Branch("theLKdr", &theLKdr, "theLKdr/F");
    theTreeNew->Branch("theKiso", &theKiso, "theKiso/F");
    theTreeNew->Branch("theBBDPhi", &theBBDPhi, "theBBDPhi/F");
    theTreeNew->Branch("theBTrkdxy2", &theBTrkdxy2, "theBTrkdxy2/F");

    theTreeNew->Branch("theAnalysisBdtG", &theAnalysisBdtG, "theAnalysisBdtG/F");
    theTreeNew->Branch("theAnalysisBdtO", &theAnalysisBdtO, "theAnalysisBdtO/F");
    
    theTreeNew->Branch("pair_mass", &pair_mass, "pair_mass/F");
    theTreeNew->Branch("weight", &weight, "weight/F");
  }

  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    // Trigger
    if (hlt9==0) continue;

    // Loop over Bs
    for (unsigned int ii=0; ii<fit_mass->size(); ii++) {

      // further selection on electrons (chiara, eventualmente da stringere cosi' nel dumper)
      if (tag_pt->at(ii)<1.5)   continue;   
      if (probe_pt->at(ii)<1.5) continue;   
      if (K_pt->at(ii)<1.5)     continue;   

      // further selection on eleID (chiara, eventualmente da stringere cosi' nel dumper)     
      if ( tag_pfmvaId->at(ii)<0 )   continue;
      if ( probe_pfmvaId->at(ii)<0 ) continue;
      if ( tag_mvaId->at(ii)<0 )     continue;
      if ( probe_mvaId->at(ii)<0 )   continue;

      // e+e- invariant mass selection (chiara, eventualmente da stringere cosi' nel dumper)     
      if (mll_raw->at(ii)<2.9 || mll_raw->at(ii)>3.2) continue;  

      // save new variables, making flat tree
      theNvtx = nvtx;     
      //
      B_mass    = fit_mass->at(ii);
      B_pt      = fit_Bpt->at(ii);
      B_cos2d   = fit_Bcos2D->at(ii);
      B_svprob  = fit_Bsvprob->at(ii);
      B_xysig   = fit_Bxysig->at(ii);
      if (sampleID>0) {
	B_matchMC = bmatchMC->at(ii); 
      } else {
	B_matchMC = -1;
      }
      B_isoRel  = B_iso04_rel->at(ii);  
      //
      tagPt       = tag_pt->at(ii);
      tagEta      = tag_eta->at(ii);
      tagPhi      = tag_phi->at(ii);
      tagPfmvaId  = tag_pfmvaId->at(ii);
      tagMvaId    = tag_mvaId->at(ii);
      tagDxySig   = tag_dxy_sig->at(ii);
      tagDzTrg    = tag_dzTrg->at(ii);
      tagIso04Rel = tag_iso04_rel->at(ii);
      //
      probePt       = probe_pt->at(ii);
      probeEta      = probe_eta->at(ii);
      probePhi      = probe_phi->at(ii);
      probePfmvaId  = probe_pfmvaId->at(ii);
      probeMvaId    = probe_mvaId->at(ii);
      probeDxySig   = probe_dxy_sig->at(ii);
      probeDzTrg    = probe_dzTrg->at(ii);
      probeIso04Rel = probe_iso04_rel->at(ii);
      //
      kPt  = K_pt->at(ii);
      kEta = K_eta->at(ii);
      kPhi = K_phi->at(ii);
      //
      theLKdz   = LKdz->at(ii); 
      theL1L2dr = L1L2dr->at(ii);   
      theLKdr   = LKdr->at(ii);   
      theKiso   = Kiso->at(ii); 
      theBBDPhi = BBDPhi->at(ii); 
      theBTrkdxy2 = BTrkdxy2->at(ii);  
      //
      kDCASig   = k_DCASig->at(ii);    
      kDzTrg    = k_dzTrg->at(ii); 
      kIso04Rel = k_iso04_rel->at(ii); 
      kSvip2d   = k_svip2d->at(ii); 
      kSvip3d   = k_svip3d->at(ii); 
      theLlkDR  = llkDR->at(ii); 
      thePtAsym = ptAsym->at(ii); 
      //
      theAnalysisBdtG = analysisBdtG->at(ii);
      theAnalysisBdtO = analysisBdtO->at(ii);
      //
      pair_mass = mll_raw->at(ii);
      elesDr = eleDR->at(ii);   
      //
      if (theRun==1) weight = perEveW*pu_weight;
      else weight = 1.;
      
      treeNew->Fill();
    }
  }
  
  // new format
  treeNew->Write();
  fileNew->Close();
  fileNew->ls();
  
  fileOrig->cd();
  fileOrig->Close();  
}

