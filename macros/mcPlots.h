//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 24 15:14:56 2021 by ROOT version 6.12/07
// from TTree fitter_tree/reduced tree for T&P with B
// found on file: /eos/cms/store/user/crovelli/LowPtEle/TnpDataB/March21noRegression/FormattedTnPForB_March21_BuToKJpsi_Toee_v2.root
//////////////////////////////////////////////////////////

#ifndef mcPlots_h
#define mcPlots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class mcPlots {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           theNvtx;
   Float_t         tagPt;
   Float_t         tagEta;
   Float_t         tagPhi;
   Float_t         tagPfmvaId;
   Float_t         tagMvaId;
   Float_t         tagDxySig;
   Float_t         tagDzTrg;
   Float_t         tagIso04Rel;
   Float_t         probePt;
   Float_t         probeEta;
   Float_t         probePhi;
   Float_t         probePfmvaId;
   Float_t         probeMvaId;
   Float_t         probeDxySig;
   Float_t         probeDzTrg;
   Float_t         probeIso04Rel;
   Float_t         elesDr;
   Float_t         kPt;
   Float_t         kEta;
   Float_t         kPhi;
   Float_t         kDCASig;
   Float_t         kDzTrg;
   Float_t         kIso04Rel;
   Float_t         kSvip2d;
   Float_t         kSvip3d;
   Float_t         B_mass;
   Float_t         B_cos2d;
   Float_t         B_pt;
   Float_t         B_svprob;
   Float_t         B_xysig;
   Float_t         B_isoRel;
   Int_t           B_matchMC;
   Float_t         theLlkDR;
   Float_t         thePtAsym;
   Float_t         theLKdz;
   Float_t         theL1L2dr;
   Float_t         theLKdr;
   Float_t         theKiso;
   Float_t         theBBDPhi;
   Float_t         theBTrkdxy2;
   Float_t         theAnalysisBdtG;
   Float_t         theAnalysisBdtO;
   Float_t         pair_mass;
   Float_t         weight;

   // List of branches
   TBranch        *b_theNvtx;   //!
   TBranch        *b_tagPt;   //!
   TBranch        *b_tagEta;   //!
   TBranch        *b_tagPhi;   //!
   TBranch        *b_tagPfmvaId;   //!
   TBranch        *b_tagMvaId;   //!
   TBranch        *b_tagDxySig;   //!
   TBranch        *b_tagDzTrg;   //!
   TBranch        *b_tagIso04Rel;   //!
   TBranch        *b_probePt;   //!
   TBranch        *b_probeEta;   //!
   TBranch        *b_probePhi;   //!
   TBranch        *b_probePfmvaId;   //!
   TBranch        *b_probeMvaId;   //!
   TBranch        *b_probeDxySig;   //!
   TBranch        *b_probeDzTrg;   //!
   TBranch        *b_probeIso04Rel;   //!
   TBranch        *b_elesDr;   //!
   TBranch        *b_kPt;   //!
   TBranch        *b_kEta;   //!
   TBranch        *b_kPhi;   //!
   TBranch        *b_kDCASig;   //!
   TBranch        *b_kDzTrg;   //!
   TBranch        *b_kIso04Rel;   //!
   TBranch        *b_kSvip2d;   //!
   TBranch        *b_kSvip3d;   //!
   TBranch        *b_B_mass;   //!
   TBranch        *b_B_cos2d;   //!
   TBranch        *b_B_pt;   //!
   TBranch        *b_B_svprob;   //!
   TBranch        *b_B_xysig;   //!
   TBranch        *b_B_isoRel;   //!
   TBranch        *b_B_matchMC;   //!
   TBranch        *b_theLlkDR;   //!
   TBranch        *b_thePtAsym;   //!
   TBranch        *b_theLKdz;   //!
   TBranch        *b_theL1L2dr;   //!
   TBranch        *b_theLKdr;   //!
   TBranch        *b_theKiso;   //!
   TBranch        *b_theBBDPhi;   //!
   TBranch        *b_theBTrkdxy2;   //!
   TBranch        *b_theAnalysisBdtG;   //!
   TBranch        *b_theAnalysisBdtO;   //!
   TBranch        *b_pair_mass;   //!
   TBranch        *b_weight;   //!

   mcPlots(TTree *tree=0);
   virtual ~mcPlots();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int checkWrtData);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mcPlots_cxx
mcPlots::mcPlots(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/user/crovelli/LowPtEle/TnpDataB/March21noRegression/FormattedTnPForB_March21_BuToKJpsi_Toee_v2.root");
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../files_first/FormattedTnPForB_mcDumperNani_noSelApartHLT.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/eos/cms/store/user/crovelli/LowPtEle/TnpDataB/March21noRegression/FormattedTnPForB_March21_BuToKJpsi_Toee_v2.root");
      //f = new TFile("../files_first/FormattedTnPForB_mcDumperNani_noSelApartHLT.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("/eos/cms/store/user/crovelli/LowPtEle/TnpDataB/March21noRegression/FormattedTnPForB_March21_BuToKJpsi_Toee_v2.root:/tnpAna");
    //TDirectory * dir = (TDirectory*)f->Get("../files_first/FormattedTnPForB_mcDumperNani_noSelApartHLT.root:/tnpAna");
    dir->GetObject("fitter_tree",tree);
    
  }
  Init(tree);
}

mcPlots::~mcPlots()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t mcPlots::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t mcPlots::LoadTree(Long64_t entry)
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

void mcPlots::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("theNvtx", &theNvtx, &b_theNvtx);
   fChain->SetBranchAddress("tagPt", &tagPt, &b_tagPt);
   fChain->SetBranchAddress("tagEta", &tagEta, &b_tagEta);
   fChain->SetBranchAddress("tagPhi", &tagPhi, &b_tagPhi);
   fChain->SetBranchAddress("tagPfmvaId", &tagPfmvaId, &b_tagPfmvaId);
   fChain->SetBranchAddress("tagMvaId", &tagMvaId, &b_tagMvaId);
   fChain->SetBranchAddress("tagDxySig", &tagDxySig, &b_tagDxySig);
   fChain->SetBranchAddress("tagDzTrg", &tagDzTrg, &b_tagDzTrg);
   fChain->SetBranchAddress("tagIso04Rel", &tagIso04Rel, &b_tagIso04Rel);
   fChain->SetBranchAddress("probePt", &probePt, &b_probePt);
   fChain->SetBranchAddress("probeEta", &probeEta, &b_probeEta);
   fChain->SetBranchAddress("probePhi", &probePhi, &b_probePhi);
   fChain->SetBranchAddress("probePfmvaId", &probePfmvaId, &b_probePfmvaId);
   fChain->SetBranchAddress("probeMvaId", &probeMvaId, &b_probeMvaId);
   fChain->SetBranchAddress("probeDxySig", &probeDxySig, &b_probeDxySig);
   fChain->SetBranchAddress("probeDzTrg", &probeDzTrg, &b_probeDzTrg);
   fChain->SetBranchAddress("probeIso04Rel", &probeIso04Rel, &b_probeIso04Rel);
   fChain->SetBranchAddress("elesDr", &elesDr, &b_elesDr);
   fChain->SetBranchAddress("kPt", &kPt, &b_kPt);
   fChain->SetBranchAddress("kEta", &kEta, &b_kEta);
   fChain->SetBranchAddress("kPhi", &kPhi, &b_kPhi);
   fChain->SetBranchAddress("kDCASig", &kDCASig, &b_kDCASig);
   fChain->SetBranchAddress("kDzTrg", &kDzTrg, &b_kDzTrg);
   fChain->SetBranchAddress("kIso04Rel", &kIso04Rel, &b_kIso04Rel);
   fChain->SetBranchAddress("kSvip2d", &kSvip2d, &b_kSvip2d);
   fChain->SetBranchAddress("kSvip3d", &kSvip3d, &b_kSvip3d);
   fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
   fChain->SetBranchAddress("B_cos2d", &B_cos2d, &b_B_cos2d);
   fChain->SetBranchAddress("B_pt", &B_pt, &b_B_pt);
   fChain->SetBranchAddress("B_svprob", &B_svprob, &b_B_svprob);
   fChain->SetBranchAddress("B_xysig", &B_xysig, &b_B_xysig);
   fChain->SetBranchAddress("B_isoRel", &B_isoRel, &b_B_isoRel);
   fChain->SetBranchAddress("B_matchMC", &B_matchMC, &b_B_matchMC);
   fChain->SetBranchAddress("theLlkDR", &theLlkDR, &b_theLlkDR);
   fChain->SetBranchAddress("thePtAsym", &thePtAsym, &b_thePtAsym);
   fChain->SetBranchAddress("theLKdz", &theLKdz, &b_theLKdz);
   fChain->SetBranchAddress("theL1L2dr", &theL1L2dr, &b_theL1L2dr);
   fChain->SetBranchAddress("theLKdr", &theLKdr, &b_theLKdr);
   fChain->SetBranchAddress("theKiso", &theKiso, &b_theKiso);
   fChain->SetBranchAddress("theBBDPhi", &theBBDPhi, &b_theBBDPhi);
   fChain->SetBranchAddress("theBTrkdxy2", &theBTrkdxy2, &b_theBTrkdxy2);
   fChain->SetBranchAddress("theAnalysisBdtG", &theAnalysisBdtG, &b_theAnalysisBdtG);
   fChain->SetBranchAddress("theAnalysisBdtO", &theAnalysisBdtO, &b_theAnalysisBdtO);
   fChain->SetBranchAddress("pair_mass", &pair_mass, &b_pair_mass);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t mcPlots::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void mcPlots::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t mcPlots::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef mcPlots_cxx
