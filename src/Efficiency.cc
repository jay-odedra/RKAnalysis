#include "../include/Efficiency.hh"
#include "../FastForest/include/fastforest.h"
#include "TLorentzVector.h"
#include <chrono>
#include <ctime>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::Loop() {
  std::cout << "Loop ..." << std::endl;

  if (fChain == 0) return;

  // Otto's BDT
//  std::cout << "Loading BDT model ..." << std::endl;
//  std::string bdtfileOttoPFPF = "./data/otto_model.txt";
//  std::vector<std::string> featOttoPFPF = {"f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10",
//					   "f11","f12","f13","f14","f15","f16","f17","f18","f19","f20",
//					   "f21", "f22","f23","f24","f25"};
//  const auto fastForestOttoPFPF = fastforest::load_txt(bdtfileOttoPFPF.c_str(), featOttoPFPF);

  // Some init
  auto start = std::chrono::system_clock::now();
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Starting event loop with " << nentries << " entries ..."<< std::endl;

  // Event loop 
  //nentries=1000;
  for (Long64_t jentry=0; jentry<nentries; ++jentry) {

    // Initialise variables
    initVars();

    // Timing
    Long64_t ientry = this->LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    //if (jentry<35050) continue; if (jentry>35050) break; //@@ debug
    //if (jentry>50001) break; //@@ debug
    if (jentry>0 && jentry%10000==0) {
      auto now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = now-start;
      std::chrono::duration<double> predicted_duration = elapsed_seconds * (nentries*1.)/(jentry*1.);
      std::chrono::system_clock::time_point end = start + std::chrono::seconds((int)predicted_duration.count());
      std::time_t end_time = std::chrono::system_clock::to_time_t(end);
      std::string tmp = std::ctime(&end_time); 
      tmp.resize(tmp.length()-1); // remove trailing \n
      std::cout << "Event: " << jentry
		<< " Elapsed: " << elapsed_seconds.count() << "s "
		<< " ETA: " << tmp
		<< std::endl;
    }
    
    // Scalars
    theRun_ = run;
    theEvent_ = event;
    
    // Triggering muon pT
    bool nTriggerMuon=0;
    int idTrgMu=-1;
    float ptTrgMu=0.;
    if(nMuon>0){
      for (u_int iM=0; iM<nMuon; iM++) {
	if(Muon_isTriggering[iM]) {
	  nTriggerMuon=nTriggerMuon+1;
	  if(Muon_pt[iM]>ptTrgMu) {
	    ptTrgMu=Muon_pt[iM];
	    idTrgMu=iM;
	  }
	}
      }
    }
    trg_muon_pt_ = ptTrgMu;
    
    // Trigger HLT paths
    hlt7_ip4_  = (int)HLT_Mu7_IP4;
    hlt8_ip3_  = (int)HLT_Mu8_IP3;
    hlt9_ip6_  = (int)HLT_Mu9_IP6;
    hlt12_ip6_ = (int)HLT_Mu12_IP6;
    
    // Find B->Kee decay 
    Daughters daughters;
    bool found = isBToKEE(daughters);
    is_bkee_ = found;
    
    // Determine if B->Kee decay in acceptance
    in_acc_ = inAcceptance(daughters.first) ? 1 : 0;

    // Set (by reference) the GEN pt,eta
    int e1_gen_idx = -1;
    int e2_gen_idx = -1;
    bool ok = setElePtEta(daughters,
			  e1_gen_idx,e1_gen_pt_,e1_gen_eta_,
			  e2_gen_idx,e2_gen_pt_,e2_gen_eta_);

    // If GEN electrons not found, continue (shouldn't happen, skip event?!)
    if (!ok) { std::cout << "ERROR!" << std::endl; continue; } // outTree_->Fill();

    // Safety limit on number of BToKEE candidates
    //int maxBToKEE = nBToKEE < 1000 ? nBToKEE : 1000;
    int maxBToKEE = nBToKEE;
    if ( maxBToKEE > 1000 ) {
      std::cout << "Too many BToKEE candidates (" << nBToKEE 
		<< "). Truncating ..." << std::endl;
      maxBToKEE = 1000;
    }
    if (verbose_>3) std::cout << "nBToKEE: " << nBToKEE
			      << " maxBToKEE: " << maxBToKEE
			      << std::endl;

    for (u_int iB=0; iB<maxBToKEE; ++iB) {
      if (verbose_>3) std:: cout << "iB: " << iB << std::endl;
      int e1_reco_idx = -1;
      int e2_reco_idx = -1;
      e1_reco_pt_ = -10.;
      e2_reco_pt_ = -10.;
      e1_reco_eta_ = -10.;
      e2_reco_eta_ = -10.;
      if ( recoCand(iB,
		    e1_gen_idx,e2_gen_idx,
		    e1_reco_idx,e1_reco_pt_,e1_reco_eta_,
		    e2_reco_idx,e2_reco_pt_,e2_reco_eta_) ) {
	isMatched_ = 1;
	break;
      }
    }

//    int matched = -1;
//    int e1_reco_idx = -1;
//    int e2_reco_idx = -1;
//    float e1_reco_pt = -10.;
//    float e2_reco_pt = -10.;
//    float e1_reco_eta = -10.;
//    float e2_reco_eta = -10.;
//    for (u_int iB=0; iB<maxBToKEE; ++iB) {
//      if (verbose_>3) std:: cout << "iB: " << iB << std::endl;
//      if ( recoCand(iB,
//		    e1_gen_idx,e2_gen_idx,
//		    e1_reco_idx,e1_reco_pt,e1_reco_eta,
//		    e2_reco_idx,e2_reco_pt,e2_reco_eta) ) {
//	matched = iB;
//	break;
//      }
//    }
//    if (matched>=0){
//      isMatched_ = 1;
//      e1_reco_pt_ = e1_reco_pt;
//      e2_reco_pt_ = e2_reco_pt; 
//      e1_reco_eta_ = e1_reco_eta;
//      e2_reco_eta_ = e2_reco_eta;
//    }

    // Fill the output tree
    outTree_->Fill();
    
  } // Event loop
  
}

////////////////////////////////////////////////////////////////////////////////
//
Efficiency::Efficiency(TChain* tree, int isMC, std::string output) :
  EfficiencyBase((TTree*)tree, isMC) 
{
  prepareOutputs(output);   
  initVars();
}

////////////////////////////////////////////////////////////////////////////////
//
Efficiency::~Efficiency() {
  outFile_->cd();
  h_entries_->Write();
  h_selection_->Write();
  outTree_->Write();
  outFile_->Close();
}     

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::initVars() {

  // Scalars
  theRun_=0;
  theEvent_=0;
  nvtx_=0;

  // Trigger
  trg_muon_pt_=0.;
  hlt7_ip4_=0;
  hlt8_ip3_=0;
  hlt9_ip6_=0;
  hlt12_ip6_=0;

  // Decay and acceptance
  is_bkee_=0;
  in_acc_=0;

  // RECO-to-GEN matching
  isMatched_=0;

  // GEN pt,eta
  e1_gen_pt_=-10.;
  e2_gen_pt_=-10.;
  e1_gen_eta_=-10.;
  e2_gen_eta_=-10.;

  // RECO pt,eta
  e1_reco_pt_=-10.;
  e2_reco_pt_=-10.;
  e1_reco_eta_=-10.;
  e2_reco_eta_=-10.;

  // Pre-selection andBDT
  ip3d_=-10.;
  cos2d_=-10.;
  bdt_=-10.;

}

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::prepareOutputs(std::string filename) {
  std::cout << "Preparing outputs ..." << std::endl;

  filename_=filename;
  std::string outname = filename_+".root";
  std::cout << "output: " << outname << std::endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");
  
  bookOutputTree();
  bookOutputHistos();

  // loading weights for pileup if needed
  //if (dopureweight_) SetPuWeights(puWFileName_);
};

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::bookOutputTree() {
  std::cout << "Booking tree ..." << std::endl;

  outTree_ = new TTree("tree","tree");

  // Scalars
  outTree_->Branch("theRun", &theRun_, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent_, "theEvent/I");    

  // Trigger
  outTree_->Branch("trg_muon_pt", &trg_muon_pt_, "trg_muon_pt/F");
  outTree_->Branch("HLT_Mu7_IP4", &hlt7_ip4_, "HLT_Mu7_IP4/I");
  outTree_->Branch("HLT_Mu8_IP3", &hlt8_ip3_, "HLT_Mu8_IP3/I");
  outTree_->Branch("HLT_Mu9_IP6", &hlt9_ip6_, "HLT_Mu9_IP6/I");
  outTree_->Branch("HLT_Mu12_IP6", &hlt12_ip6_, "HLT_Mu12_IP6/I");

  // Decay and acceptance
  outTree_->Branch("isBToKEE", &is_bkee_, "isBToKEE/I");
  outTree_->Branch("inAcc", &in_acc_, "inAcc/I");

  // RECO-to-GEN matching
  outTree_->Branch("isMatched", &isMatched_, "isMatched/I");  

  // GEN pt,eta
  outTree_->Branch("e1_gen_pt", &e1_gen_pt_, "e1_gen_pt/F");  
  outTree_->Branch("e2_gen_pt", &e2_gen_pt_, "e2_gen_pt/F");  
  outTree_->Branch("e1_gen_eta", &e1_gen_eta_, "e1_gen_eta/F");  
  outTree_->Branch("e2_gen_eta", &e2_gen_eta_, "e2_gen_eta/F");  

  //RECO pt,eta
  outTree_->Branch("e1_reco_pt", &e1_reco_pt_, "e1_reco_pt/F");  
  outTree_->Branch("e2_reco_pt", &e2_reco_pt_, "e2_reco_pt/F");  
  outTree_->Branch("e1_reco_eta", &e1_reco_eta_, "e1_reco_eta/F");  
  outTree_->Branch("e2_reco_eta", &e2_reco_eta_, "e2_reco_eta/F");  

  // Pre-selection and BDT
  outTree_->Branch("ip3d", &ip3d_, "ip3d/F");  
  outTree_->Branch("cos2d", &cos2d_, "cos2d/F");  
  outTree_->Branch("bdt", &bdt_, "bdt/F");  

}

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::bookOutputHistos() {
  std::cout << "Booking histos ..." << std::endl;
  h_entries_   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection_ = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
}

////////////////////////////////////////////////////////////////////////////////
//
bool Efficiency::isBToKEE(Efficiency::Daughters& daughters) {
  daughters = Daughters();
  Cands cands;
  if (verbose_>3) std::cout << "### isBToKEE: find daughters to Bs ... " << std::endl;
  for ( int idx = 0; idx < nGenPart; ++idx ) {
    int d_idx = idx;
    int d_pdg = GenPart_pdgId[d_idx];
    int m_idx = GenPart_genPartIdxMother[d_idx];
    int m_pdg = GenPart_pdgId[m_idx];
    if ( ( abs(m_pdg) == 521 && abs(d_pdg) == 11 ) ||   // find daughter electrons
	 ( abs(m_pdg) == 521 && abs(d_pdg) == 321 ) ) { // find daughter kaon
      if ( cands.find(m_idx) == cands.end() ) {
	cands[m_idx] = std::make_pair(Ints(),Ints());//Daughters();//
      }
      cands[m_idx].first.push_back(d_idx);
      cands[m_idx].second.push_back(d_pdg);
      if (verbose_>3) std::cout << " d_idx: " << d_idx 
				<< " d_pdg: " << d_pdg
				<< " m_idx: " << m_idx
				<< " m_pdg: " << m_pdg
				<< std::endl;
    }
  }
  Cands::iterator iter;
  if (verbose_>3) std::cout << "### isBToKEE: loop through Bs and their daughters ... " << std::endl;
  for (iter = cands.begin(); iter != cands.end(); ++iter) {
    Ints::iterator begin = iter->second.second.begin();
    Ints::iterator end = iter->second.second.end();
    if (verbose_>3) std::cout << " m_idx: " << iter->first
			      << " m_pdg: " << GenPart_pdgId[iter->first]
			      << " m_ndaughters: " << iter->second.second.size()
			      << std::endl;
    for ( int i = 0; i < iter->second.second.size(); ++i) {
      int d_idx = iter->second.first[i];
      int d_pdg = iter->second.second[i];
      int m_idx = GenPart_genPartIdxMother[d_idx];
      int m_pdg = GenPart_pdgId[m_idx];
      if (verbose_>3) std::cout << "  d_idx: " << d_idx
				<< " d_pdg: " << d_pdg
				<< " m_idx: " << m_idx
				<< " m_pdg: " << m_pdg
				<< std::endl;
    }
    if ( ( iter->second.second.size() == 3 ) &&
	 ( std::find(begin,end,11) != end ) &&
	 ( std::find(begin,end,-11) != end ) &&
	 ( ( std::find(begin,end,321) != end ) ||
	   ( std::find(begin,end,-321) != end ) ) ) {
      if (verbose_>3) std::cout << " B->KEE decay!"
				<< std::endl;
      daughters = iter->second; // returns (by ref) the first B->Kee decay found
      return true; 
    } 
    if (verbose_>3) std::cout << std::endl;
  }
  return false;
}
 
////////////////////////////////////////////////////////////////////////////////
//
bool Efficiency::inAcceptance(std::vector<int>& indices) {
  bool in_acc = true;
  if (verbose_>3) std::cout << "### inAcceptance: check if daughters in acceptance ..." << std::endl;
  for ( int i = 0; i < indices.size(); ++i ) {
    int d_idx = indices[i];
    int d_pdg = GenPart_pdgId[d_idx];
    float d_pt = GenPart_pt[d_idx];
    float d_eta = GenPart_eta[d_idx];
    //if (abs(d_pdg)==11) continue; //@@ ignore kaon !!!
    if (verbose_>3) std::cout << " i: " << i
			      << " d_idx: " << d_idx
			      << " d_pt: " << d_pt
			      << " d_eta: " << d_eta
			      << std::endl;
    if ( d_pt < 0.5 || fabs(d_eta) > 2.5 ) { in_acc = false; }
  }
  if (verbose_>3) std::cout << " IN ACC? " << in_acc
			    << std::endl;
  return in_acc;
}

////////////////////////////////////////////////////////////////////////////////
//
bool Efficiency::setElePtEta(Efficiency::Daughters& daughters,
			     int& e1_gen_idx,
			     float& e1_gen_pt,
			     float& e1_gen_eta,
			     int& e2_gen_idx,
			     float& e2_gen_pt,
			     float& e2_gen_eta) {
  Ints indices;
  for ( int i = 0; i < daughters.first.size(); ++i ) {
    int d_idx = daughters.first[i];
    int d_pdg = daughters.second[i];
    if ( abs(d_pdg) == 11 ) { indices.push_back(d_idx); }
  }
  if (indices.size()!=2) return false;

  int d1_idx = indices[0];
  int d2_idx = indices[1];
  float d1_pt = GenPart_pt[d1_idx];
  float d2_pt = GenPart_pt[d2_idx];

  e1_gen_idx = d1_pt > d2_pt ? d1_idx : d2_idx; // find leading electron
  e2_gen_idx = d1_pt > d2_pt ? d2_idx : d1_idx; // find sub-leading electron
  if (e1_gen_idx==e2_gen_idx) return false;

  e1_gen_pt = GenPart_pt[e1_gen_idx];
  e2_gen_pt = GenPart_pt[e2_gen_idx];
  e1_gen_eta = GenPart_eta[e1_gen_idx];
  e2_gen_eta = GenPart_eta[e2_gen_idx];
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool Efficiency::recoCand(int theB,
			  int& e1_gen_idx,int& e2_gen_idx,
			  int& e1_reco_idx,float& e1_reco_pt,float& e1_reco_eta,
			  int& e2_reco_idx,float& e2_reco_pt,float& e2_reco_eta) {
  
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

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

  bool found = ( ( ele1_genPartIdx >= 0
		   && ele2_genPartIdx >= 0
		   && k_genPartIdx >= 0
		   ) &&
		 ( abs(ele1_genPdgId) == 11
		   && abs(ele2_genPdgId) == 11
		   && abs(k_genPdgId) == 321 
		   ) &&
		 //( ele1_genPdgId * ele1_genPdgId < 0 ) &&
		 ( abs(ele1_genMotherPdgId) == 521
		   && abs(ele2_genMotherPdgId) == 521
		   && abs(k_genMotherPdgId) == 521
		   ) );

  if (found) {
//    std::cout << "FOUND " 
//	      << e1_gen_idx << " " 
//	      << ele1_genPartIdx << " " 
//	      << e2_gen_idx << " " 
//	      << ele2_genPartIdx << " " 
//	      << std::endl;
    if (ele1_genPartIdx==e1_gen_idx){
      e1_reco_idx = ele1_idx;
      e1_reco_pt  = Electron_pt[ele1_idx];
      e1_reco_eta = Electron_eta[ele1_idx];
    } else if (ele1_genPartIdx==e2_gen_idx){
      e2_reco_idx = ele1_idx;
      e2_reco_pt  = Electron_pt[ele1_idx];
      e2_reco_eta = Electron_eta[ele1_idx];
    }
    if (ele2_genPartIdx==e1_gen_idx){
      e1_reco_idx = ele2_idx;
      e1_reco_pt  = Electron_pt[ele2_idx];
      e1_reco_eta = Electron_eta[ele2_idx];
    } else if (ele2_genPartIdx==e2_gen_idx){
      e2_reco_idx = ele2_idx;
      e2_reco_pt  = Electron_pt[ele2_idx];
      e2_reco_eta = Electron_eta[ele2_idx];
    }
  }

  if (verbose_>3) std::cout << "recoCand: " << found <<  std::endl;
  return found;

}
