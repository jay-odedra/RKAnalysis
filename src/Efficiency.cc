#include "../include/Efficiency.hh"
#include "TLorentzVector.h"
#include <chrono>
#include <ctime>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::Loop() {
  std::cout << "Loop ..." << std::endl;

  if (fChain == 0) return;

  fastforest::FastForest fastForestOttoPFPF;
  fastforest::FastForest fastForestOttoPFLP;
  loadModels(fastForestOttoPFPF,fastForestOttoPFLP);
  
  // Some init
  auto start = std::chrono::system_clock::now();
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Starting event loop with " << nentries << " entries ..."<< std::endl;

  // Event loop 
  int cntr1 = 0;
  int cntr2 = 0;
  int cntr3 = 0;
  Long64_t jentry=0;
  for (; jentry<nentries; ++jentry) {

    std::vector<float> analysisBdtO;
    std::vector<float> analysisBdtG;

    // Limit events to process
    //if (jentry>5000) break;

    // Limit to event range
    //if (jentry<19181) continue; if (jentry>19181) break;

    // Initialise variables
    initVars();

    // Load tree and count events
    Long64_t ientry = this->LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    cntr1++;

    // Timing per interval
    int interval = 1000;
    if (jentry%interval==0) { timing(nentries,jentry,start); }
    
    // Scalars
    theRun_ = run;
    theEvent_ = event;
    
    // Triggering muon pT
    int nTriggerMuon=0;
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
			  e2_gen_idx,e2_gen_pt_,e2_gen_eta_,
			  e12_gen_dr_);

    // If GEN electrons not found, continue (shouldn't happen, skip event?!)
    if (!ok) { std::cout << "ERROR!" << std::endl; continue; } // outTree_->Fill();

    // Safety limit on number of BToKEE candidates
    int maxBToKEE = EfficiencyBase::nBToKEE_max_;
    if ( nBToKEE < maxBToKEE ) { maxBToKEE = nBToKEE; }
    else {
      std::cout << "Too many BToKEE candidates (" << nBToKEE 
		<< "). Truncating to " << maxBToKEE
		<< " ..." << std::endl;
      cntr3++;
    }
    if (verbose_>3) std::cout << "nBToKEE: " << nBToKEE
			      << " maxBToKEE: " << maxBToKEE
			      << std::endl;

    // Loop through all RECO'ed B->Kee candidates and find the GEN-matched triplet
    for (u_int iB=0; iB<maxBToKEE; ++iB) {
      if (verbose_>4) std:: cout << "iB: " << iB << std::endl;
      int e1_reco_idx = -1;
      int e2_reco_idx = -1;
      e1_reco_pt_ = -10.;
      e2_reco_pt_ = -10.;
      e1_reco_eta_ = -10.;
      e2_reco_eta_ = -10.;
      e12_reco_dr_ = -10.;
      e1_reco_pf_ = 0;
      e2_reco_pf_ = 0;
      e1_reco_lowpt_ = 0;
      e2_reco_lowpt_ = 0;
      e1_reco_overlap_ = 0;
      e2_reco_overlap_ = 0;
      if ( recoCand(iB,
		    e1_gen_idx,e2_gen_idx,
		    e1_reco_idx,e1_reco_pt_,e1_reco_eta_,
		    e1_reco_pf_,e1_reco_lowpt_,e1_reco_overlap_,
		    e2_reco_idx,e2_reco_pt_,e2_reco_eta_,
		    e2_reco_pf_,e2_reco_lowpt_,e2_reco_overlap_,
		    e12_reco_dr_) ) {
	// Pre-selection and BDT
	ip3d_ = BToKEE_k_svip3d[iB];
	cos2d_ = BToKEE_fit_cos2D[iB];
	bdt_ = evaluateModels(iB,fastForestOttoPFPF,fastForestOttoPFLP);
	mll_ = BToKEE_mll_fullfit[iB];
	// Mark if matched and break loop
	isMatched_ = 1;
	h_cand_->Fill(iB<1000?iB:1000,1.); // overflows into final bin
	break;
      }
    }

    // Fill the output tree
    cntr2++;
    outTree_->Fill();

  } // Event loop

  // Final timing
  timing(nentries,jentry,start);
  
  std::cout << "Summary:" << std::endl
	    << "  Number of events processed: " << cntr1 << std::endl
	    << "  Number of events written to tree: " << cntr2 << std::endl
	    << "  Number of events with too many B->Kee candidates: " << cntr3 << std::endl;

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
  h_cand_->Write();
  outTree_->Write();
  std::cout << "Written TTree '" << outTree_->GetName()
	    << "' to file '" << outFile_->GetName() 
	    << "' ..."
	    << std::endl;
  outFile_->Close();
}     

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::timing( int nentries, int jentry, auto start ) {
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
  e12_gen_dr_=-10.;

  // RECO pt,eta
  e1_reco_pt_=-10.;
  e2_reco_pt_=-10.;
  e1_reco_eta_=-10.;
  e2_reco_eta_=-10.;
  e12_reco_dr_=-10.;

  // RECO algo
  e1_reco_pf_=0;
  e2_reco_pf_=0;
  e1_reco_lowpt_=0;
  e2_reco_lowpt_=0;
  e1_reco_overlap_=0;
  e2_reco_overlap_=0;

  // Pre-selection and BDT
  ip3d_=-10.;
  cos2d_=-10.;
  bdt_=-10.;
  mll_=0.;

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
  outTree_->Branch("e12_gen_dr", &e12_gen_dr_, "e12_gen_dr/F");

  //RECO pt,eta
  outTree_->Branch("e1_reco_pt", &e1_reco_pt_, "e1_reco_pt/F");
  outTree_->Branch("e2_reco_pt", &e2_reco_pt_, "e2_reco_pt/F");
  outTree_->Branch("e1_reco_eta", &e1_reco_eta_, "e1_reco_eta/F");
  outTree_->Branch("e2_reco_eta", &e2_reco_eta_, "e2_reco_eta/F");
  outTree_->Branch("e12_reco_dr", &e12_reco_dr_, "e12_reco_dr/F");

  //RECO algo
  outTree_->Branch("e1_reco_pf", &e1_reco_pf_, "e1_reco_pf/I");
  outTree_->Branch("e2_reco_pf", &e2_reco_pf_, "e2_reco_pf/I");
  outTree_->Branch("e1_reco_lowpt", &e1_reco_lowpt_, "e1_reco_lowpt/I");
  outTree_->Branch("e2_reco_lowpt", &e2_reco_lowpt_, "e2_reco_lowpt/I");
  outTree_->Branch("e1_reco_overlap", &e1_reco_overlap_, "e1_reco_overlap/I");
  outTree_->Branch("e2_reco_overlap", &e2_reco_overlap_, "e2_reco_overlap/I");

  // Pre-selection and BDT
  outTree_->Branch("ip3d", &ip3d_, "ip3d/F");
  outTree_->Branch("cos2d", &cos2d_, "cos2d/F");
  outTree_->Branch("bdt", &bdt_, "bdt/F");
  outTree_->Branch("mll", &mll_, "mll/F");

}

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::bookOutputHistos() {
  std::cout << "Booking histos ..." << std::endl;
  h_entries_   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection_ = new TH1F("h_selection","Selection breakdown", 8, -0.5, 7.5);
  h_cand_ = new TH1F("h_cand","Candidate index", 1001, -0.5, 1000.5);
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
			     float& e2_gen_eta,
			     float& e12_gen_dr) {
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
  float e1_gen_phi = GenPart_phi[e1_gen_idx];
  float e2_gen_phi = GenPart_phi[e2_gen_idx];
  e12_gen_dr = DeltaR(e1_gen_eta,
		      e1_gen_phi,
		      e2_gen_eta,
		      e2_gen_phi);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
bool Efficiency::recoCand(int theB,
			  int& e1_gen_idx,int& e2_gen_idx,
			  int& e1_reco_idx,float& e1_reco_pt,float& e1_reco_eta,
			  int& e1_reco_pf,int& e1_reco_lowpt,int& e1_reco_overlap,
			  int& e2_reco_idx,float& e2_reco_pt,float& e2_reco_eta,
			  int& e2_reco_pf,int& e2_reco_lowpt,int& e2_reco_overlap,
			  float& e12_reco_dr
			  ) {
  if (verbose_>4) std::cout << "Efficiency::recoCand:" 
			    << " e1_gen_idx: " << e1_gen_idx
			    << " e2_gen_idx: " << e2_gen_idx
			    << std::endl;
  
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int kaon_idx  = BToKEE_kIdx[theB];
  if (verbose_>4) std::cout << " ele1_idx: " << ele1_idx
			    << " ele2_idx: " << ele2_idx
			    << " kaon_idx: " << kaon_idx
			    << std::endl;
  if (ele1_idx == ele2_idx) {
    std::cout << "Same reco'ed electrons in the B->Ke candidate!!!"
	      << " theB: " << theB
	      << " ele1_idx: " << ele1_idx
	      << std::endl;
    return false;
  }

  int ele1_isPF = Electron_isPF[ele1_idx];
  int ele2_isPF = Electron_isPF[ele2_idx];
  int ele1_isLowPt = Electron_isLowPt[ele1_idx];
  int ele2_isLowPt = Electron_isLowPt[ele2_idx];
  int ele1_isPFoverlap = Electron_isPFoverlap[ele1_idx];
  int ele2_isPFoverlap = Electron_isPFoverlap[ele2_idx];
  if (verbose_>4) std::cout << " ele1_isPF: " << ele1_isPF
			    << " ele2_isPF: " << ele2_isPF
			    << std::endl
			    << " ele1_isLowPt: " << ele1_isLowPt
			    << " ele2_isLowPt: " << ele2_isLowPt
			    << std::endl
			    << " ele1_isPFoverlap: " << ele1_isPFoverlap
			    << " ele2_isPFoverlap: " << ele2_isPFoverlap
			    << std::endl;
  
  int ele1_genPartIdx = Electron_genPartIdx[ele1_idx];
  int ele2_genPartIdx = Electron_genPartIdx[ele2_idx];
  int kaon_genPartIdx = ProbeTracks_genPartIdx[kaon_idx];
  if (verbose_>4) std::cout << " ele1_genPartIdx: " << ele1_genPartIdx
			    << " ele2_genPartIdx: " << ele2_genPartIdx
			    << " kaon_genPartIdx: " << kaon_genPartIdx
			    << std::endl;
  if (ele1_genPartIdx >= EfficiencyBase::nGenPart_max_ || 
      ele2_genPartIdx >= EfficiencyBase::nGenPart_max_ || 
      kaon_genPartIdx >= EfficiencyBase::nGenPart_max_ ) {
    if (verbose_>4) std::cout << "Incorrect GenPart index!!!"
			      << " ele1_genPartIdx: " << ele1_genPartIdx
			      << " ele2_genPartIdx: " << ele2_genPartIdx
			      << " kaon_genPartIdx: " << kaon_genPartIdx
			      << std::endl; 
    return false;
  }

  int ele1_genMotherIdx = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele2_genMotherIdx = GenPart_genPartIdxMother[ele2_genPartIdx];
  int kaon_genMotherIdx = GenPart_genPartIdxMother[kaon_genPartIdx];
  if (verbose_>4) std::cout << " ele1_genMotherIdx: " << ele1_genMotherIdx
			    << " ele2_genMotherIdx: " << ele2_genMotherIdx
			    << " kaon_genMotherIdx: " << kaon_genMotherIdx
			    << std::endl;
  
  int ele1_genGMotherIdx = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele2_genGMotherIdx = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int kaon_genGMotherIdx = GenPart_genPartIdxMother[kaon_genMotherIdx];
  if (verbose_>4) std::cout << " ele1_genGMotherIdx: " << ele1_genGMotherIdx
			    << " ele2_genGMotherIdx: " << ele2_genGMotherIdx
			    << " kaon_genGMotherIdx: " << kaon_genGMotherIdx
			    << std::endl;
  
  int ele1_genPdgId = GenPart_pdgId[ele1_genPartIdx];
  int ele2_genPdgId = GenPart_pdgId[ele2_genPartIdx];
  int kaon_genPdgId = GenPart_pdgId[kaon_genPartIdx];
  if (verbose_>4) std::cout << " ele1_genPdgId: " << ele1_genPdgId
			    << " ele2_genPdgId: " << ele2_genPdgId
			    << " kaon_genPdgId: " << kaon_genPdgId
			    << std::endl;
  
  int ele1_genMotherPdgId = GenPart_pdgId[ele1_genMotherIdx];
  int ele2_genMotherPdgId = GenPart_pdgId[ele2_genMotherIdx];
  int kaon_genMotherPdgId = GenPart_pdgId[kaon_genMotherIdx];
  if (verbose_>4) std::cout << " ele1_genMotherPdgId: " << ele1_genMotherPdgId
			    << " ele2_genMotherPdgId: " << ele2_genMotherPdgId
			    << " kaon_genMotherPdgId: " << kaon_genMotherPdgId
			    << std::endl;
  
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];
  int kaon_genGMotherPdgId = GenPart_pdgId[kaon_genGMotherIdx];
  if (verbose_>4) std::cout << " ele1_genGMotherPdgId: " << ele1_genGMotherPdgId
			    << " ele2_genGMotherPdgId: " << ele2_genGMotherPdgId
			    << " kaon_genGMotherPdgId: " << kaon_genGMotherPdgId
			    << std::endl;
  
  bool found = ( ( ele1_genPartIdx >= 0
		   && ele2_genPartIdx >= 0
		   && kaon_genPartIdx >= 0
		   ) &&
		 ( abs(ele1_genPdgId) == 11
		   && abs(ele2_genPdgId) == 11
		   && abs(kaon_genPdgId) == 321 
		   ) &&
		 //( ele1_genPdgId * ele1_genPdgId < 0 ) &&
		 ( abs(ele1_genMotherPdgId) == 521
		   && abs(ele2_genMotherPdgId) == 521
		   && abs(kaon_genMotherPdgId) == 521
		   ) );
  
  if (found) {
    float e1_reco_phi = 0.;
    float e2_reco_phi = 0.;
    if (verbose_>4) std::cout << "FOUND " 
			      << e1_gen_idx << " " 
			      << ele1_genPartIdx << " " 
			      << e2_gen_idx << " " 
			      << ele2_genPartIdx << " " 
			      << std::endl;
    if (e1_gen_idx==ele1_genPartIdx){
      e1_reco_idx = ele1_idx;
      e1_reco_pt  = Electron_pt[ele1_idx];
      e1_reco_eta = Electron_eta[ele1_idx];
      e1_reco_phi = Electron_phi[ele1_idx];
      e1_reco_pf  = Electron_isPF[ele1_idx];
      e1_reco_lowpt = Electron_isLowPt[ele1_idx];
      e1_reco_overlap = Electron_isPFoverlap[ele1_idx];
    } 
    if (e1_gen_idx==ele2_genPartIdx){
      e1_reco_idx = ele2_idx;
      e1_reco_pt  = Electron_pt[ele2_idx];
      e1_reco_eta = Electron_eta[ele2_idx];
      e1_reco_phi = Electron_phi[ele2_idx];
      e1_reco_pf  = Electron_isPF[ele2_idx];
      e1_reco_lowpt = Electron_isLowPt[ele2_idx];
      e1_reco_overlap = Electron_isPFoverlap[ele2_idx];
    }
    if (e2_gen_idx==ele1_genPartIdx){
      e2_reco_idx = ele1_idx;
      e2_reco_pt  = Electron_pt[ele1_idx];
      e2_reco_eta = Electron_eta[ele1_idx];
      e2_reco_phi = Electron_phi[ele1_idx];
      e2_reco_pf  = Electron_isPF[ele1_idx];
      e2_reco_lowpt = Electron_isLowPt[ele1_idx]; //@@ was e1_reco_lowpt (bug!)
      e2_reco_overlap = Electron_isPFoverlap[ele1_idx];
    }
    if (e2_gen_idx==ele2_genPartIdx){
      e2_reco_idx = ele2_idx;
      e2_reco_pt  = Electron_pt[ele2_idx];
      e2_reco_eta = Electron_eta[ele2_idx];
      e2_reco_phi = Electron_phi[ele2_idx];
      e2_reco_pf  = Electron_isPF[ele2_idx];
      e2_reco_lowpt = Electron_isLowPt[ele2_idx];
      e2_reco_overlap = Electron_isPFoverlap[ele2_idx];
    }
    e12_reco_dr = DeltaR(e1_reco_eta,
			 e1_reco_phi,
			 e2_reco_eta,
			 e2_reco_phi);
  }

  if (verbose_>4) std::cout << "recoCand: " << found <<  std::endl;
  return found;

}

////////////////////////////////////////////////////////////////////////////////
//
float Efficiency::DeltaR(float eta1,
			 float phi1,
			 float eta2,
			 float phi2) {
  float PI=3.1415972;
  float deta=eta1-eta2;
  float dphi=phi1-phi2;
  if(dphi>PI){ dphi-=2.0*PI; }
  else if(dphi<=-PI){ dphi+=2.0*PI; }
  return TMath::Sqrt(deta*deta+dphi*dphi);
}


////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::loadModels(fastforest::FastForest& fastForestOttoPFPF,
			    fastforest::FastForest& fastForestOttoPFLP) {
  std::cout << "Efficiency::loadModels ..." << std::endl;
  
  // PFPF
  std::string bdtfileOttoPFPF = "./data/otto_model.txt";
  std::vector<std::string> featOttoPFPF = {
    "f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13",
    "f14","f15","f16","f17","f18","f19","f20","f21", "f22","f23","f24","f25"
  };
  fastForestOttoPFPF = fastforest::load_txt(bdtfileOttoPFPF.c_str(), featOttoPFPF);

  // PFLP
  std::string bdtfileOttoPFLP = "./data/otto_model_pflp.txt";
  std::vector<std::string> featOttoPFLP = {
    "f0","f1", "f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13",
    "f14","f15","f16","f17","f18","f19","f20","f21", "f22","f23","f24","f25","f26","f27"
  };
  fastForestOttoPFLP = fastforest::load_txt(bdtfileOttoPFLP.c_str(), featOttoPFLP);

}


////////////////////////////////////////////////////////////////////////////////
//
float Efficiency::evaluateModels(u_int thisB,
				 fastforest::FastForest& fastForestOttoPFPF,
				 fastforest::FastForest& fastForestOttoPFLP) {
  //std::cout << "Efficiency::evaluateModels ..." << std::endl;
  
  int ele1_idx = BToKEE_l1Idx[thisB];
  int ele2_idx = BToKEE_l2Idx[thisB];
  int k_idx    = BToKEE_kIdx[thisB];
  
  float k_pt     = BToKEE_fit_k_pt[thisB];   
  float ele1_pt  = BToKEE_fit_l1_pt[thisB];
  float ele2_pt  = BToKEE_fit_l2_pt[thisB];

  float k_eta    = BToKEE_fit_k_eta[thisB];   
  float ele1_eta = BToKEE_fit_l1_eta[thisB];
  float ele2_eta = BToKEE_fit_l2_eta[thisB];

  float k_phi    = BToKEE_fit_k_phi[thisB];   
  float ele1_phi = BToKEE_fit_l1_phi[thisB];
  float ele2_phi = BToKEE_fit_l2_phi[thisB];
  
  TLorentzVector ele1TLV(0,0,0,0);
  ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.000511);
  TLorentzVector ele2TLV(0,0,0,0);
  ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.000511);
  TLorentzVector kTLV(0,0,0,0);
  kTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.493677);
  
  float thisBmass   = BToKEE_fit_mass[thisB];
  float thisBpt     = BToKEE_fit_pt[thisB];
  float thisBcos    = BToKEE_fit_cos2D[thisB];
  float thisBsvprob = BToKEE_svprob[thisB];
  float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
  
  float BToKEE_fit_l1_normpt=BToKEE_fit_l1_pt[thisB]/BToKEE_fit_mass[thisB];
  float BToKEE_fit_l2_normpt=BToKEE_fit_l2_pt[thisB]/BToKEE_fit_mass[thisB];
  float BToKEE_l1_dxy_sig=(Electron_dxy[ele1_idx]) /Electron_dxyErr[ele1_idx];
  float BToKEE_l2_dxy_sig=(Electron_dxy[ele2_idx]) /Electron_dxyErr[ele2_idx];
  float BToKEE_fit_k_normpt=BToKEE_fit_k_pt[thisB] /BToKEE_fit_mass[thisB];
  float BToKEE_k_DCASig=ProbeTracks_DCASig[k_idx];
  float BToKEE_k_dxy_sig=ProbeTracks_dxyS[k_idx];
  float BToKEE_fit_normpt=BToKEE_fit_pt[thisB] /BToKEE_fit_mass[thisB];
  float BToKEE_l_xy_sig = (BToKEE_l_xy[thisB]) /BToKEE_l_xy_unc[thisB];
  float BToKEE_eleDR= DeltaR(ele1_eta,ele1_phi,ele2_eta,ele2_phi);
  TLorentzVector dll=ele1TLV+ele2TLV;
  
  float BToKEE_llkDR=dll.DeltaR(kTLV);
  float BToKEE_l1_iso04_rel=BToKEE_l1_iso04[thisB]/BToKEE_fit_l1_pt[thisB];
  float BToKEE_l2_iso04_rel=BToKEE_l2_iso04[thisB]/BToKEE_fit_l2_pt[thisB];
  float BToKEE_k_iso04_rel = BToKEE_k_iso04[thisB] / BToKEE_fit_k_pt[thisB];
  float BToKEE_b_iso04_rel =BToKEE_b_iso04[thisB]/BToKEE_fit_pt[thisB];
  TVector3 diele_p3 = dll.Vect();
  TVector3 k_p3 = kTLV.Vect();
  TVector3 pv2sv_p3(PV_x-BToKEE_vtx_x[thisB], PV_y-BToKEE_vtx_y[thisB], PV_z-BToKEE_vtx_z[thisB]);
  
  float BToKEE_ptAsym = ( (diele_p3.Cross(pv2sv_p3)).Mag() - (k_p3.Cross(pv2sv_p3)).Mag() ) 
    / ( (diele_p3.Cross(pv2sv_p3)).Mag() + (k_p3.Cross(pv2sv_p3)).Mag() );
  
  float BToKEE_l1_dzTrg=Electron_dzTrg[ele1_idx];
  float BToKEE_l2_dzTrg=Electron_dzTrg[ele2_idx];
  float BToKEE_k_dzTrg=ProbeTracks_dzTrg[k_idx];
  
  float BToKEE_l1_pfmvaId_lowPt=Electron_pfmvaId[ele1_idx];
  if(ele1_pt>5) BToKEE_l1_pfmvaId_lowPt=20;
  
  float BToKEE_l2_pfmvaId_lowPt=Electron_pfmvaId[ele2_idx];
  if(ele2_pt>5) BToKEE_l2_pfmvaId_lowPt=20;
  
  float BToKEE_l1_pfmvaId_highPt=Electron_pfmvaId[ele1_idx];
  if(ele1_pt<=5) BToKEE_l1_pfmvaId_highPt=20;
  
  float BToKEE_l2_pfmvaId_highPt=Electron_pfmvaId[ele2_idx];
  if(ele2_pt<=5) BToKEE_l2_pfmvaId_highPt=20;
  
  float scoreBdtOtto=0.;
  if(Electron_isPF[ele1_idx]&&Electron_isPF[ele2_idx]){
    std::vector<float> vecBdtOtto = {BToKEE_b_iso04_rel, BToKEE_eleDR, thisBcos, BToKEE_fit_k_normpt,
				     BToKEE_fit_l1_normpt,BToKEE_fit_l2_normpt, BToKEE_fit_normpt,
				     BToKEE_k_DCASig, BToKEE_k_dzTrg, BToKEE_k_iso04_rel, 
				     BToKEE_k_svip2d[thisB], BToKEE_k_svip3d[thisB],
				     BToKEE_l1_dxy_sig, BToKEE_l1_dzTrg, BToKEE_l1_iso04_rel, 
				     //l1_mvaId,
				     BToKEE_l1_pfmvaId_highPt, BToKEE_l1_pfmvaId_lowPt,
				     BToKEE_l2_dxy_sig,BToKEE_l2_dzTrg,BToKEE_l2_iso04_rel,
				     // l2_mvaId,
				     BToKEE_l2_pfmvaId_highPt, BToKEE_l2_pfmvaId_lowPt, 
				     thisBxysig,BToKEE_llkDR, BToKEE_ptAsym, thisBsvprob};
    return fastForestOttoPFPF(vecBdtOtto.data());
  } else {
    std::vector<float> vecBdtOtto = {BToKEE_b_iso04_rel, BToKEE_eleDR, thisBcos, BToKEE_fit_k_normpt,
				     BToKEE_fit_l1_normpt,BToKEE_fit_l2_normpt, BToKEE_fit_normpt,
				     BToKEE_k_DCASig, BToKEE_k_dzTrg, BToKEE_k_iso04_rel, 
				     BToKEE_k_svip2d[thisB], BToKEE_k_svip3d[thisB],
				     BToKEE_l1_dxy_sig, BToKEE_l1_dzTrg, BToKEE_l1_iso04_rel, 
				     Electron_mvaId[ele1_idx],
				     BToKEE_l1_pfmvaId_highPt, BToKEE_l1_pfmvaId_lowPt,
				     BToKEE_l2_dxy_sig,BToKEE_l2_dzTrg,BToKEE_l2_iso04_rel,
				     Electron_mvaId[ele2_idx],
				     BToKEE_l2_pfmvaId_highPt, BToKEE_l2_pfmvaId_lowPt, 
				     thisBxysig,BToKEE_llkDR, BToKEE_ptAsym, thisBsvprob};
    return fastForestOttoPFLP(vecBdtOtto.data());
  }

  return -10.;
  
}

