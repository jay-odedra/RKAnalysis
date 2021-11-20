#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "../FastForest/include/fastforest.h"

#include <cmath>
#include <iostream> 

#include "../include/TestMva.hh"    

using namespace std;

TestMva::TestMva(TTree *tree)     
  : BParkBase(tree) {        

  // Chiara: to be set by hand  
  sampleID = 1;   // 0 = data, >=1 MC  
}

TestMva::~TestMva() {

  // output
  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
}     

void TestMva::Loop() {

  if (fChain == 0) return;

  // Load Analysis MVA weights
  std::string bdtfile = "../data/xgbmodel_kee_final_12B_0.txt";
  std::vector<std::string> feat = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11"};
  const auto fastForest = fastforest::load_txt(bdtfile.c_str(), feat);

  // Loop over events
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%50000==0) cout << jentry << endl;

    // Event info     
    theRun   = run;
    theEvent = event;
    theSampleID = sampleID;

    // # Vertices
    nvtx = PV_npvs;
    
    // B candidates
    if (nBToKEE<=0) continue;

    for (u_int thisB=0; thisB<nBToKEE; thisB++) {

      float thisBmass   = BToKEE_fit_mass[thisB];
      float thisBpt     = BToKEE_fit_pt[thisB];
      float thisBcos    = BToKEE_fit_cos2D[thisB];
      float thisBsvprob = BToKEE_svprob[thisB];
      float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcB(thisB);      

      // Triplet
      float k_pt     = BToKEE_fit_k_pt[thisB];   
      float ele1_pt  = BToKEE_fit_l1_pt[thisB];
      float ele2_pt  = BToKEE_fit_l2_pt[thisB];
      float k_eta    = BToKEE_fit_k_eta[thisB];   
      float ele1_eta = BToKEE_fit_l1_eta[thisB];
      float ele2_eta = BToKEE_fit_l2_eta[thisB];
      float k_phi    = BToKEE_fit_k_phi[thisB];   
      float ele1_phi = BToKEE_fit_l1_phi[thisB];
      float ele2_phi = BToKEE_fit_l2_phi[thisB];

      // Indices
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];
      int k_idx    = BToKEE_kIdx[thisB];  

      // Vertices
      float kvx = ProbeTracks_vx[k_idx];
      float kvy = ProbeTracks_vy[k_idx];
      float kvz = ProbeTracks_vz[k_idx];

      float l1vx = Electron_vx[ele1_idx];
      float l1vy = Electron_vy[ele1_idx];
      float l1vz = Electron_vz[ele1_idx];

      float l2vx = Electron_vx[ele2_idx];
      float l2vy = Electron_vy[ele2_idx];
      float l2vz = Electron_vz[ele2_idx];

      // Vertex distances
      float l1kDz = fabs(kvz-l1vz);
      float l2kDz = fabs(kvz-l2vz);
      float LKdz = l1kDz;
      if (l2kDz<LKdz) LKdz = l2kDz;

      // Vectors
      TVector3 l1v3;
      TVector3 l2v3;
      TVector3 kv3;
      l1v3.SetPtEtaPhi(ele1_pt,ele1_eta,ele1_phi);
      l2v3.SetPtEtaPhi(ele2_pt,ele2_eta,ele2_phi);
      kv3.SetPtEtaPhi(k_pt,k_eta,k_phi);
	
      // DeltaR
      float l1kDr = l1v3.DeltaR(kv3);
      float l2kDr = l2v3.DeltaR(kv3);
      float LKdr = l1kDr;
      if (l2kDr<LKdr) LKdr = l2kDr;
      float L1L2dr = l1v3.DeltaR(l2v3);

      // Isolation
      float L1iso = BToKEE_l1_iso03[thisB];
      float L2iso = BToKEE_l2_iso03[thisB];
	
      // Analysis MVA
      std::vector<float> vecBdt = {thisBsvprob, thisBxysig, ele1_pt, ele2_pt, k_pt, thisBcos, thisBpt, LKdz, LKdr, L1L2dr, L1iso, L2iso};
      float scoreBdt = fastForest(vecBdt.data());

      // Tree fill
      Bmass.push_back(thisBmass);
      Bpt.push_back(thisBpt);
      Bcos2D.push_back(thisBcos);
      Bsvprob.push_back(thisBsvprob);
      Bxysig.push_back(thisBxysig);
      BmatchMC.push_back(isThisAMcB);
      
      Ele1Pt.push_back(ele1_pt);
      Ele1Eta.push_back(ele1_eta);
      Ele1Phi.push_back(ele1_phi);
      Ele1IsPF.push_back(Electron_isPF[ele1_idx]);           
      Ele1IsLowPt.push_back(Electron_isLowPt[ele1_idx]);
      Ele1MvaId.push_back(Electron_mvaId[ele1_idx]);
      Ele1PfmvaId.push_back(Electron_pfmvaId[ele1_idx]);
      Ele1Iso03.push_back(L1iso);
      
      Ele2Pt.push_back(ele2_pt);
      Ele2Eta.push_back(ele2_eta);
      Ele2Phi.push_back(ele2_phi);
      Ele2IsPF.push_back(Electron_isPF[ele2_idx]);           
      Ele2IsLowPt.push_back(Electron_isLowPt[ele2_idx]);
      Ele2MvaId.push_back(Electron_mvaId[ele2_idx]);
      Ele2PfmvaId.push_back(Electron_pfmvaId[ele2_idx]);
      Ele2Iso03.push_back(L2iso);
      
      KPt.push_back(k_pt);
      KEta.push_back(k_eta);
      KPhi.push_back(k_phi);
      
      dzL1K.push_back(l1kDz);
      dzL2K.push_back(l2kDz);
      dzLK.push_back(LKdz);
      
      drL1K.push_back(l1kDr);
      drL2K.push_back(l2kDr);
      drLK.push_back(LKdr);

      drL1L2.push_back(L1L2dr);
      
      analysisBdt.push_back(scoreBdt);

      if (sampleID>0) {     // MC      
	bool l1mc = (Electron_genPartIdx[ele1_idx]>-0.5);
	bool l2mc = (Electron_genPartIdx[ele2_idx]>-0.5);
	Ele1MatchMc.push_back(l1mc);
	Ele2MatchMc.push_back(l2mc);
      } else {
	Ele1MatchMc.push_back(0);
	Ele2MatchMc.push_back(0);
      }

    } // Loop over good Bs
      
    // At least one tag and one probe
    if (drL1K.size()<=0) continue;

    // Filling the output tree
    outTree_->Fill();

    // Cleaning all vectors used for the output tree, ready for a new entry
    Bmass.clear();
    Bpt.clear();
    Bcos2D.clear();
    Bsvprob.clear();
    Bxysig.clear();
    BmatchMC.clear();

    Ele1Pt.clear();  
    Ele1Eta.clear();  
    Ele1Phi.clear();  
    Ele1IsPF.clear();  
    Ele1IsLowPt.clear();  
    Ele1MvaId.clear();  
    Ele1PfmvaId.clear();  
    Ele1Iso03.clear();

    Ele2Pt.clear();  
    Ele2Eta.clear();  
    Ele2Phi.clear();  
    Ele2IsPF.clear();  
    Ele2IsLowPt.clear();  
    Ele2MvaId.clear();  
    Ele2PfmvaId.clear();  
    Ele2Iso03.clear();

    KPt.clear();
    KEta.clear();
    KPhi.clear();
      
    dzL1K.clear();
    dzL2K.clear();
    dzLK.clear();
    drL1K.clear();
    drL2K.clear();
    drLK.clear();
    drL1L2.clear();

    analysisBdt.clear();

    Ele1MatchMc.clear();
    Ele2MatchMc.clear();
  }
}

bool TestMva::isMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);
  bool RK_res = okMatch && RK_res1 && RK_res2;

  return RK_res;
}

void TestMva::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
};


void TestMva::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;

  outTree_->Branch("theRun", &theRun, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  
  outTree_->Branch("Ele1Pt",  "std::vector<float>", &Ele1Pt);  
  outTree_->Branch("Ele1Eta", "std::vector<float>", &Ele1Eta);  
  outTree_->Branch("Ele1Phi", "std::vector<float>", &Ele1Phi);  
  outTree_->Branch("Ele1IsPF", "std::vector<bool>", &Ele1IsPF);  
  outTree_->Branch("Ele1IsLowPt", "std::vector<bool>", &Ele1IsLowPt);  
  outTree_->Branch("Ele1MvaId", "std::vector<float>", &Ele1MvaId);  
  outTree_->Branch("Ele1PfmvaId", "std::vector<float>", &Ele1PfmvaId);  
  outTree_->Branch("Ele1Iso03", "std::vector<float>", &Ele1Iso03);  
  outTree_->Branch("Ele1MatchMc", "std::vector<bool>", &Ele1MatchMc);

  outTree_->Branch("Ele2Pt",  "std::vector<float>", &Ele2Pt);  
  outTree_->Branch("Ele2Eta", "std::vector<float>", &Ele2Eta);  
  outTree_->Branch("Ele2Phi", "std::vector<float>", &Ele2Phi);  
  outTree_->Branch("Ele2IsPF", "std::vector<bool>", &Ele2IsPF);  
  outTree_->Branch("Ele2IsLowPt", "std::vector<bool>", &Ele2IsLowPt);  
  outTree_->Branch("Ele2MvaId", "std::vector<float>", &Ele2MvaId);  
  outTree_->Branch("Ele2PfmvaId", "std::vector<float>", &Ele2PfmvaId);  
  outTree_->Branch("Ele2Iso03", "std::vector<float>", &Ele2Iso03);  
  outTree_->Branch("Ele2MatchMc", "std::vector<bool>", &Ele2MatchMc);

  outTree_->Branch("Bmass",   "std::vector<float>", &Bmass);  
  outTree_->Branch("Bpt",     "std::vector<float>", &Bpt);  
  outTree_->Branch("Bcos2D",  "std::vector<float>", &Bcos2D);  
  outTree_->Branch("Bsvprob", "std::vector<float>", &Bsvprob);  
  outTree_->Branch("Bxysig",  "std::vector<float>", &Bxysig);  
  outTree_->Branch("BmatchMC", "std::vector<bool>", &BmatchMC);  

  outTree_->Branch("KPt",  "std::vector<float>", &KPt);  
  outTree_->Branch("KEta", "std::vector<float>", &KEta);  
  outTree_->Branch("KPhi", "std::vector<float>", &KPhi);  

  outTree_->Branch("dzL1K",  "std::vector<float>", &dzL1K);
  outTree_->Branch("dzL2K",  "std::vector<float>", &dzL2K);
  outTree_->Branch("dzLK",   "std::vector<float>", &dzLK);
  outTree_->Branch("drL1K",  "std::vector<float>", &drL1K);
  outTree_->Branch("drL2K",  "std::vector<float>", &drL2K);
  outTree_->Branch("drLK",   "std::vector<float>", &drLK);
  outTree_->Branch("drL1L2", "std::vector<float>", &drL1L2);

  outTree_->Branch("analysisBdt", "std::vector<float>", &analysisBdt);
}

