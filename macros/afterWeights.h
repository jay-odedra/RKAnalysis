//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 25 16:28:04 2021 by ROOT version 6.12/07
// from TTree TaPtree/TaPtree
// found on file: ppp.root
//////////////////////////////////////////////////////////

#ifndef afterWeights_h
#define afterWeights_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class afterWeights {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
   // Declaration of leaf types
   Int_t           theRun;
   Int_t           theEvent;
   Int_t           nvtx;
   Int_t           sampleID;
   Float_t         rho;
   Int_t           hlt9;
   Float_t         trg_muon_pt;
   Int_t           selectedBSize;
   vector<float>   *tag_pt;
   vector<float>   *tag_eta;
   vector<float>   *tag_phi;
   vector<bool>    *tag_isPF;
   vector<bool>    *tag_isPFOverlap;
   vector<bool>    *tag_isLowPt;
   vector<float>   *tag_mvaId;
   vector<float>   *tag_pfmvaId;
   vector<bool>    *tag_convveto;
   vector<float>   *tag_ptMc;
   vector<float>   *probe_ptMc;
   vector<bool>    *tag_matchMcFromJPsi;
   vector<bool>    *tag_matchMc;
   vector<float>   *mll_fullfit;
   vector<float>   *mll_raw;
   vector<float>   *fit_mass;
   vector<float>   *fit_Bpt;
   vector<float>   *fit_Bcos2D;
   vector<float>   *fit_Bsvprob;
   vector<float>   *fit_Bxysig;
   vector<bool>    *bmatchMC;
   vector<float>   *K_pt;
   vector<float>   *K_eta;
   vector<float>   *K_phi;
   vector<float>   *probe_pt;
   vector<float>   *probe_eta;
   vector<float>   *probe_phi;
   vector<bool>    *probe_isPF;
   vector<bool>    *probe_isPFOverlap;
   vector<bool>    *probe_isLowPt;
   vector<float>   *probe_mvaId;
   vector<float>   *probe_pfmvaId;
   vector<float>   *probe_dxySig;
   vector<float>   *probe_dzSig;
   vector<bool>    *probe_convveto;
   vector<bool>    *probe_matchMcFromJPsi;
   vector<bool>    *probe_matchMc;
   vector<float>   *probe_drm;
   vector<float>   *tag_drm;
   vector<float>   *bPtOrder;
   vector<float>   *b2DOrder;
   vector<float>   *bXYOrder;
   vector<float>   *bMassOrder;
   vector<float>   *p4Trk;
   vector<float>   *analysisBdtG;
   vector<float>   *analysisBdtO;
   vector<float>   *analysisBdtGWithSyst;
   vector<float>   *analysisBdtOWithSyst;
   vector<float>   *LKdz;
   vector<float>   *L1L2dr;
   vector<float>   *LKdr;
   vector<float>   *Kiso;
   vector<float>   *BBDPhi;
   vector<float>   *BTrkdxy2;
   vector<float>   *Dmass;
   vector<float>   *Dmass_flip;
   vector<float>   *Dmass_ll;
   vector<float>   *Dmass_ll_flip;
   vector<float>   *B_iso04_rel;
   vector<float>   *eleDR;
   vector<float>   *k_DCASig;
   vector<float>   *k_dzTrg;
   vector<float>   *k_dxy_sig;
   vector<float>   *k_iso04_rel;
   vector<float>   *k_svip2d;
   vector<float>   *k_svip3d;
   vector<float>   *tag_dxy_sig;
   vector<float>   *tag_dzTrg;
   vector<float>   *tag_iso04_rel;
   vector<float>   *probe_dxy_sig;
   vector<float>   *probe_dzTrg;
   vector<float>   *probe_iso04_rel;
   vector<float>   *llkDR;
   vector<float>   *ptAsym;
   Int_t           selectedPairsSize;

   // List of branches
   TBranch        *b_theRun;   //!
   TBranch        *b_theEvent;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_sampleID;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_hlt9;   //!
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
   TBranch        *b_tag_ptMc;   //!
   TBranch        *b_probe_ptMc;   //!
   TBranch        *b_tag_matchMcFromJPsi;   //!
   TBranch        *b_tag_matchMc;   //!
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
   TBranch        *b_probe_dxySig;   //!
   TBranch        *b_probe_dzSig;   //!
   TBranch        *b_probe_convveto;   //!
   TBranch        *b_probe_matchMcFromJPsi;   //!
   TBranch        *b_probe_matchMc;   //!
   TBranch        *b_probe_drm;   //!
   TBranch        *b_tag_drm;   //!
   TBranch        *b_bPtOrder;   //!
   TBranch        *b_b2DOrder;   //!
   TBranch        *b_bXYOrder;   //!
   TBranch        *b_bMassOrder;   //!
   TBranch        *b_p4Trk;   //!
   TBranch        *b_analysisBdtG;   //!
   TBranch        *b_analysisBdtO;   //!
   TBranch        *b_analysisBdtGWithSyst;   //!
   TBranch        *b_analysisBdtOWithSyst;   //!
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

   afterWeights(TTree *tree=0);
   virtual ~afterWeights();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef afterWeights_cxx
afterWeights::afterWeights(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../files_first/mcDumperNani_noSelApartHLT__withWeights.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../files_first/mcDumperNani_noSelApartHLT__withWeights.root");
    }
    f->GetObject("TaPtree",tree);
    
  }
  Init(tree);
}

afterWeights::~afterWeights()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t afterWeights::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t afterWeights::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void afterWeights::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer
  tag_pt = 0;
  tag_eta = 0;
  tag_phi = 0;
  tag_isPF = 0;
  tag_isPFOverlap = 0;
  tag_isLowPt = 0;
  tag_mvaId = 0;
   tag_pfmvaId = 0;
   tag_convveto = 0;
   tag_ptMc = 0;
   probe_ptMc = 0;
   tag_matchMcFromJPsi = 0;
   tag_matchMc = 0;
   mll_fullfit = 0;
   mll_raw = 0;
   fit_mass = 0;
   fit_Bpt = 0;
   fit_Bcos2D = 0;
   fit_Bsvprob = 0;
   fit_Bxysig = 0;
   bmatchMC = 0;
   K_pt = 0;
   K_eta = 0;
   K_phi = 0;
   probe_pt = 0;
   probe_eta = 0;
   probe_phi = 0;
   probe_isPF = 0;
   probe_isPFOverlap = 0;
   probe_isLowPt = 0;
   probe_mvaId = 0;
   probe_pfmvaId = 0;
   probe_dxySig = 0;
   probe_dzSig = 0;
   probe_convveto = 0;
   probe_matchMcFromJPsi = 0;
   probe_matchMc = 0;
   probe_drm = 0;
   tag_drm = 0;
   bPtOrder = 0;
   b2DOrder = 0;
   bXYOrder = 0;
   bMassOrder = 0;
   p4Trk = 0;
   analysisBdtG = 0;
   analysisBdtO = 0;
   analysisBdtGWithSyst = 0;
   analysisBdtOWithSyst = 0;
   LKdz = 0;
   L1L2dr = 0;
   LKdr = 0;
   Kiso = 0;
   BBDPhi = 0;
   BTrkdxy2 = 0;
   Dmass = 0;
   Dmass_flip = 0;
   Dmass_ll = 0;
   Dmass_ll_flip = 0;
   B_iso04_rel = 0;
   eleDR = 0;
   k_DCASig = 0;
   k_dzTrg = 0;
   k_dxy_sig = 0;
   k_iso04_rel = 0;
   k_svip2d = 0;
   k_svip3d = 0;
   tag_dxy_sig = 0;
   tag_dzTrg = 0;
   tag_iso04_rel = 0;
   probe_dxy_sig = 0;
   probe_dzTrg = 0;
   probe_iso04_rel = 0;
   llkDR = 0;
   ptAsym = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("theRun", &theRun, &b_theRun);
   fChain->SetBranchAddress("theEvent", &theEvent, &b_theEvent);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("hlt9", &hlt9, &b_hlt9);
   fChain->SetBranchAddress("trg_muon_pt", &trg_muon_pt, &b_trg_muon_pt);
   fChain->SetBranchAddress("selectedBSize", &selectedBSize, &b_selectedBSize);
   fChain->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
   fChain->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
   fChain->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
   fChain->SetBranchAddress("tag_isPF", &tag_isPF, &b_tag_isPF);
   fChain->SetBranchAddress("tag_isPFOverlap", &tag_isPFOverlap, &b_tag_isPFOverlap);
   fChain->SetBranchAddress("tag_isLowPt", &tag_isLowPt, &b_tag_isLowPt);
   fChain->SetBranchAddress("tag_mvaId", &tag_mvaId, &b_tag_mvaId);
   fChain->SetBranchAddress("tag_pfmvaId", &tag_pfmvaId, &b_tag_pfmvaId);
   fChain->SetBranchAddress("tag_convveto", &tag_convveto, &b_tag_convveto);
   fChain->SetBranchAddress("tag_ptMc", &tag_ptMc, &b_tag_ptMc);
   fChain->SetBranchAddress("probe_ptMc", &probe_ptMc, &b_probe_ptMc);
   fChain->SetBranchAddress("tag_matchMcFromJPsi", &tag_matchMcFromJPsi, &b_tag_matchMcFromJPsi);
   fChain->SetBranchAddress("tag_matchMc", &tag_matchMc, &b_tag_matchMc);
   fChain->SetBranchAddress("mll_fullfit", &mll_fullfit, &b_mll_fullfit);
   fChain->SetBranchAddress("mll_raw", &mll_raw, &b_mll_raw);
   fChain->SetBranchAddress("fit_mass", &fit_mass, &b_fit_mass);
   fChain->SetBranchAddress("fit_Bpt", &fit_Bpt, &b_fit_Bpt);
   fChain->SetBranchAddress("fit_Bcos2D", &fit_Bcos2D, &b_fit_Bcos2D);
   fChain->SetBranchAddress("fit_Bsvprob", &fit_Bsvprob, &b_fit_Bsvprob);
   fChain->SetBranchAddress("fit_Bxysig", &fit_Bxysig, &b_fit_Bxysig);
   fChain->SetBranchAddress("bmatchMC", &bmatchMC, &b_bmatchMC);
   fChain->SetBranchAddress("K_pt", &K_pt, &b_K_pt);
   fChain->SetBranchAddress("K_eta", &K_eta, &b_K_eta);
   fChain->SetBranchAddress("K_phi", &K_phi, &b_K_phi);
   fChain->SetBranchAddress("probe_pt", &probe_pt, &b_probe_pt);
   fChain->SetBranchAddress("probe_eta", &probe_eta, &b_probe_eta);
   fChain->SetBranchAddress("probe_phi", &probe_phi, &b_probe_phi);
   fChain->SetBranchAddress("probe_isPF", &probe_isPF, &b_probe_isPF);
   fChain->SetBranchAddress("probe_isPFOverlap", &probe_isPFOverlap, &b_probe_isPFOverlap);
   fChain->SetBranchAddress("probe_isLowPt", &probe_isLowPt, &b_probe_isLowPt);
   fChain->SetBranchAddress("probe_mvaId", &probe_mvaId, &b_probe_mvaId);
   fChain->SetBranchAddress("probe_pfmvaId", &probe_pfmvaId, &b_probe_pfmvaId);
   fChain->SetBranchAddress("probe_dxySig", &probe_dxySig, &b_probe_dxySig);
   fChain->SetBranchAddress("probe_dzSig", &probe_dzSig, &b_probe_dzSig);
   fChain->SetBranchAddress("probe_convveto", &probe_convveto, &b_probe_convveto);
   fChain->SetBranchAddress("probe_matchMcFromJPsi", &probe_matchMcFromJPsi, &b_probe_matchMcFromJPsi);
   fChain->SetBranchAddress("probe_matchMc", &probe_matchMc, &b_probe_matchMc);
   fChain->SetBranchAddress("probe_drm", &probe_drm, &b_probe_drm);
   fChain->SetBranchAddress("tag_drm", &tag_drm, &b_tag_drm);
   fChain->SetBranchAddress("bPtOrder", &bPtOrder, &b_bPtOrder);
   fChain->SetBranchAddress("b2DOrder", &b2DOrder, &b_b2DOrder);
   fChain->SetBranchAddress("bXYOrder", &bXYOrder, &b_bXYOrder);
   fChain->SetBranchAddress("bMassOrder", &bMassOrder, &b_bMassOrder);
   fChain->SetBranchAddress("p4Trk", &p4Trk, &b_p4Trk);
   fChain->SetBranchAddress("analysisBdtG", &analysisBdtG, &b_analysisBdtG);
   fChain->SetBranchAddress("analysisBdtO", &analysisBdtO, &b_analysisBdtO);
   fChain->SetBranchAddress("analysisBdtGWithSyst", &analysisBdtGWithSyst, &b_analysisBdtGWithSyst);
   fChain->SetBranchAddress("analysisBdtOWithSyst", &analysisBdtOWithSyst, &b_analysisBdtOWithSyst);
   fChain->SetBranchAddress("LKdz", &LKdz, &b_LKdz);
   fChain->SetBranchAddress("L1L2dr", &L1L2dr, &b_L1L2dr);
   fChain->SetBranchAddress("LKdr", &LKdr, &b_LKdr);
   fChain->SetBranchAddress("Kiso", &Kiso, &b_Kiso);
   fChain->SetBranchAddress("BBDPhi", &BBDPhi, &b_BBDPhi);
   fChain->SetBranchAddress("BTrkdxy2", &BTrkdxy2, &b_BTrkdxy2);
   fChain->SetBranchAddress("Dmass", &Dmass, &b_Dmass);
   fChain->SetBranchAddress("Dmass_flip", &Dmass_flip, &b_Dmass_flip);
   fChain->SetBranchAddress("Dmass_ll", &Dmass_ll, &b_Dmass_ll);
   fChain->SetBranchAddress("Dmass_ll_flip", &Dmass_ll_flip, &b_Dmass_ll_flip);
   fChain->SetBranchAddress("B_iso04_rel", &B_iso04_rel, &b_B_iso04_rel);
   fChain->SetBranchAddress("eleDR", &eleDR, &b_eleDR);
   fChain->SetBranchAddress("k_DCASig", &k_DCASig, &b_k_DCASig);
   fChain->SetBranchAddress("k_dzTrg", &k_dzTrg, &b_k_dzTrg);
   fChain->SetBranchAddress("k_dxy_sig", &k_dxy_sig, &b_k_dxy_sig);
   fChain->SetBranchAddress("k_iso04_rel", &k_iso04_rel, &b_k_iso04_rel);
   fChain->SetBranchAddress("k_svip2d", &k_svip2d, &b_k_svip2d);
   fChain->SetBranchAddress("k_svip3d", &k_svip3d, &b_k_svip3d);
   fChain->SetBranchAddress("tag_dxy_sig", &tag_dxy_sig, &b_tag_dxy_sig);
   fChain->SetBranchAddress("tag_dzTrg", &tag_dzTrg, &b_tag_dzTrg);
   fChain->SetBranchAddress("tag_iso04_rel", &tag_iso04_rel, &b_tag_iso04_rel);
   fChain->SetBranchAddress("probe_dxy_sig", &probe_dxy_sig, &b_probe_dxy_sig);
   fChain->SetBranchAddress("probe_dzTrg", &probe_dzTrg, &b_probe_dzTrg);
   fChain->SetBranchAddress("probe_iso04_rel", &probe_iso04_rel, &b_probe_iso04_rel);
   fChain->SetBranchAddress("llkDR", &llkDR, &b_llkDR);
   fChain->SetBranchAddress("ptAsym", &ptAsym, &b_ptAsym);
   fChain->SetBranchAddress("selectedPairsSize", &selectedPairsSize, &b_selectedPairsSize);
   Notify();
}

Bool_t afterWeights::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void afterWeights::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t afterWeights::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef afterWeights_cxx
