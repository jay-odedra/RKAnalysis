#include "../include/Efficiency.hh"
#include "TLorentzVector.h"
#include <chrono>
#include <ctime>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::Loop() {
  std::cout << "Loop ..." << std::endl;

  //bool resonant = true;

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
    //if (jentry>1000) break;

    // Limit to event range
    //if (jentry>19181) break; if (jentry<19181) continue;
    //if (jentry>4) break; if (jentry!=4) continue;

    // Debug
    //std::cout << std::endl << "##### Event number: " << jentry << std::endl;

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
    theLumi_ = lumi;
    theEvent_ = event;

    Rho_fixedGridRhoAll_ =   Rho_fixedGridRhoAll;
    Rho_fixedGridRhoFastjetAll_ =   Rho_fixedGridRhoFastjetAll;
    Rho_fixedGridRhoFastjetCentral_ =   Rho_fixedGridRhoFastjetCentral;
    Rho_fixedGridRhoFastjetCentralCalo_ =   Rho_fixedGridRhoFastjetCentralCalo;
    Rho_fixedGridRhoFastjetCentralChargedPileUp_ =   Rho_fixedGridRhoFastjetCentralChargedPileUp;
    Rho_fixedGridRhoFastjetCentralNeutral_ =   Rho_fixedGridRhoFastjetCentralNeutral;
    
//    // Triggering muon pT
//    int nTriggerMuon=0;
//    int idTrgMu=-1;
//    float ptTrgMu=0.;
//    float etaTrgMu=-10.;
//    if(nMuon>0){
//      for (u_int iM=0; iM<nMuon; iM++) {
//	if(Muon_isTriggering[iM]) {
//	  nTriggerMuon=nTriggerMuon+1;
//	  if(Muon_pt[iM]>ptTrgMu) {
//	    ptTrgMu=Muon_pt[iM];
//	    etaTrgMu=Muon_eta[iM];
//	    idTrgMu=iM;
//	  }
//	}
//      }
//    }
//    trg_muon_pt_ = ptTrgMu;
//    trg_muon_eta_ = etaTrgMu;
    
    // Trigger HLT paths
//    hlt7_ip4_  = (int)HLT_Mu7_IP4;
//    hlt8_ip3_  = (int)HLT_Mu8_IP3;
//    hlt9_ip6_  = (int)HLT_Mu9_IP6;
//    hlt12_ip6_ = (int)HLT_Mu12_IP6;

    hlt_10p0_ = (int)HLT_DoubleEle10_eta1p22_mMax6;
    hlt_9p5_ = (int)HLT_DoubleEle9p5_eta1p22_mMax6;
    hlt_9p0_ = (int)HLT_DoubleEle9_eta1p22_mMax6;
    hlt_8p5_ = (int)HLT_DoubleEle8p5_eta1p22_mMax6;
    hlt_8p0_ = (int)HLT_DoubleEle8_eta1p22_mMax6;
    hlt_7p5_ = (int)HLT_DoubleEle7p5_eta1p22_mMax6;
    hlt_7p0_ = (int)HLT_DoubleEle7_eta1p22_mMax6;
    hlt_6p5_ = (int)HLT_DoubleEle6p5_eta1p22_mMax6;
    hlt_6p0_ = (int)HLT_DoubleEle6_eta1p22_mMax6;
    hlt_5p5_ = (int)HLT_DoubleEle5p5_eta1p22_mMax6;
    hlt_5p0_ = (int)HLT_DoubleEle5_eta1p22_mMax6;
    hlt_4p5_ = (int)HLT_DoubleEle4p5_eta1p22_mMax6;
    hlt_4p0_ = (int)HLT_DoubleEle4_eta1p22_mMax6;
    
    l1_11p0_ = (int)L1_DoubleEG11_er1p2_dR_Max0p6;
    l1_10p5_ = (int)L1_DoubleEG10p5_er1p2_dR_Max0p6;
    l1_10p0_ = (int)L1_DoubleEG10_er1p2_dR_Max0p6;
    l1_9p5_ = (int)L1_DoubleEG9p5_er1p2_dR_Max0p6;
    l1_9p0_ = (int)L1_DoubleEG9_er1p2_dR_Max0p7;
    l1_8p5_ = (int)L1_DoubleEG8p5_er1p2_dR_Max0p7;
    l1_8p0_ = (int)L1_DoubleEG8_er1p2_dR_Max0p7;
    l1_7p5_ = (int)L1_DoubleEG7p5_er1p2_dR_Max0p7;
    l1_7p0_ = (int)L1_DoubleEG7_er1p2_dR_Max0p8;
    l1_6p5_ = (int)L1_DoubleEG6p5_er1p2_dR_Max0p8;
    l1_6p0_ = (int)L1_DoubleEG6_er1p2_dR_Max0p8;
    l1_5p5_ = (int)L1_DoubleEG5p5_er1p2_dR_Max0p8;
    l1_5p0_ = (int)L1_DoubleEG5_er1p2_dR_Max0p9;
    l1_4p5_ = (int)L1_DoubleEG4p5_er1p2_dR_Max0p9;
    l1_4p0_ = (int)L1_DoubleEG4_er1p2_dR_Max0p9;

    // JSON filters
    if (verbose_>1) {
      std::cout << "[Efficiency::Loop]"
		<< " Run: " << theRun_
		<< " LS: " << theLumi_
		<< " Event: " << theEvent_ << std::endl;
    }
    if (verbose_>2) {
      std::cout << "[Efficiency::Loop]"
		<< " The following JSONs identify a valid run/lumi: " << std::endl;
    }
//      int idx = 0;
//      for ( auto& filter : jsonFilters_ ) {
//	bool ok = filter.isGoodRunLS(theRun_,theLumi_);
//	if (ok) jsonFlags_[idx] = true;
//	idx++;
//      }
    for ( auto& filter : jsonFilters_ ) {
      bool ok = filter.isGoodRunLS(theRun_,theLumi_);
      if (ok) {
	std::string name = JsonFilter::jsonFileName( filter.jsonFilePath() );
  if      (name == "trigger_OR")                   {tmp0 = 1; }
  else if (name == "L1_10p5_HLT_5p0_Excl_Final") {tmp1 = 1;}
  else if (name == "L1_10p5_HLT_6p5_Excl_Final") {tmp2 = 1;}
  else if (name == "L1_11p0_HLT_6p5_Excl_Final") {tmp3 = 1;}
  else if (name == "L1_4p5_HLT_4p0_Excl_Final") {tmp4 = 1;}
  else if (name == "L1_5p0_HLT_4p0_Excl_Final") {tmp5 = 1;}
  else if (name == "L1_5p5_HLT_4p0_Excl_Final") {tmp6 = 1;}
  else if (name == "L1_5p5_HLT_6p0_Excl_Final") {tmp7 = 1;}
  else if (name == "L1_6p0_HLT_4p0_Excl_Final") {tmp8 = 1;}
  else if (name == "L1_6p5_HLT_4p5_Excl_Final") {tmp9 = 1;}
  else if (name == "L1_7p0_HLT_5p0_Excl_Final") {tmp10 = 1;}
  else if (name == "L1_7p5_HLT_5p0_Excl_Final") {tmp11 = 1;}
  else if (name == "L1_8p0_HLT_5p0_Excl_Final") {tmp12 = 1;}
  else if (name == "L1_8p5_HLT_5p0_Excl_Final") {tmp13 = 1;}
  else if (name == "L1_8p5_HLT_5p5_Excl_Final") {tmp14 = 1;}
  else if (name == "L1_9p0_HLT_6p0_Excl_Final") {tmp15 = 1;}
	//else { std::cout << "FILTER NOT FOUND!!!" << name << std::endl; }
      }
    }

    int e1_gen_idx = -1;
    int e2_gen_idx = -1;
    if (_isMC) {
    
      // Find B->Kee decay
      Daughters daughters;
      bool found = isBToKEE(daughters);//, resonant);
      is_bkee_ = found;
      
      // Determine if B->Kee decay in acceptance
      in_acc_ = inAcceptance(daughters.first) ? 1 : 0;
      
      // Set (by reference) the GEN pt,eta
      bool ok = setElePtEta(daughters,
			    e1_gen_idx,e1_gen_pt_,e1_gen_eta_,
			    e2_gen_idx,e2_gen_pt_,e2_gen_eta_,
			    e12_gen_dr_);

      // If GEN electrons not found, continue (shouldn't happen, skip event?!)
      if (!ok) { 
	std::cout << "ERROR!"
		  << " daughters.second.size(): " << daughters.second.size()
		  << ", ";
	for ( const auto& d : daughters.second ) { std::cout << d << " "; }
	std::cout << std::endl; 
	continue; 
      } // outTree_->Fill();

    } else { // is data ...
      is_bkee_ = true;
      in_acc_ = true;
    }

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
      e1_reco_pf_ = -10;
      e2_reco_pf_ = -10;
      e1_reco_lowpt_ = -10;
      e2_reco_lowpt_ = -10;
      e1_reco_overlap_ = -10;
      e2_reco_overlap_ = -10;
      e1_reco_loose_ = -10;
      e2_reco_loose_ = -10;
      e1_reco_medium_ = -10;
      e2_reco_medium_ = -10;
      e1_reco_tight_ = -10;
      e2_reco_tight_ = -10;

      bool found = recoCand(iB,
			    e1_gen_idx,e2_gen_idx,
			    e1_reco_idx,e1_reco_pt_,e1_reco_eta_,
			    e1_reco_pf_,e1_reco_lowpt_,e1_reco_overlap_,
			    e1_reco_loose_,e1_reco_medium_,e1_reco_tight_,
			    e2_reco_idx,e2_reco_pt_,e2_reco_eta_,
			    e2_reco_pf_,e2_reco_lowpt_,e2_reco_overlap_,
			    e2_reco_loose_,e2_reco_medium_,e2_reco_tight_,
			    e12_reco_dr_);//, resonant);
      //std::cout << "TEST " << _isMC << " " << found << std::endl;
      if ( found ) {
	// Pre-selection and BDT
	ip3d_ = BToKEE_k_svip3d[iB];
	cos2d_ = BToKEE_fit_cos2D[iB];
	bdt_ = evaluateModels(iB,fastForestOttoPFPF,fastForestOttoPFLP);
	mll_ = BToKEE_mll_fullfit[iB];
	// B candidates
	b_mass_ = BToKEE_fit_mass[iB];
	b_mass_err_ = BToKEE_fit_massErr[iB];
	b_pt_ = BToKEE_fit_pt[iB];
	b_l1_pt_ = BToKEE_fit_l1_pt[iB];
	b_l2_pt_ = BToKEE_fit_l2_pt[iB];
	b_k_pt_ = BToKEE_fit_k_pt[iB];
	b_lxy_ = BToKEE_l_xy[iB];
	b_lxyerr_ = BToKEE_l_xy_unc[iB];
	b_svprob_ = BToKEE_svprob[iB];

  // all bdt vars
  
  int ele1_idx = BToKEE_l1Idx[iB];
  int ele2_idx = BToKEE_l2Idx[iB];
  int k_idx    = BToKEE_kIdx[iB];
  
  float k_pt     = BToKEE_fit_k_pt[iB];   
  float ele1_pt  = BToKEE_fit_l1_pt[iB];
  float ele2_pt  = BToKEE_fit_l2_pt[iB];

  float k_eta    = BToKEE_fit_k_eta[iB];   
  float ele1_eta = BToKEE_fit_l1_eta[iB];
  float ele2_eta = BToKEE_fit_l2_eta[iB];

  float k_phi    = BToKEE_fit_k_phi[iB];   
  float ele1_phi = BToKEE_fit_l1_phi[iB];
  float ele2_phi = BToKEE_fit_l2_phi[iB];
  
  TLorentzVector ele1TLV(0,0,0,0);
  ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.000511);
  TLorentzVector ele2TLV(0,0,0,0);
  ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.000511);
  TLorentzVector kTLV(0,0,0,0);
  kTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.493677);
  
  float thisBmass   = BToKEE_fit_mass[iB];
  float thisBpt     = BToKEE_fit_pt[iB];
  float thisBcos    = BToKEE_fit_cos2D[iB];
  float thisBsvprob = BToKEE_svprob[iB];
  float thisBxysig  = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];
  
  BToKEE_fit_l1_normpt_=BToKEE_fit_l1_pt[iB]/BToKEE_fit_mass[iB];
  BToKEE_fit_l2_normpt_=BToKEE_fit_l2_pt[iB]/BToKEE_fit_mass[iB];
  BToKEE_l1_dxy_sig_=(Electron_dxy[ele1_idx]) /Electron_dxyErr[ele1_idx];
  BToKEE_l2_dxy_sig_=(Electron_dxy[ele2_idx]) /Electron_dxyErr[ele2_idx];
  BToKEE_fit_k_normpt_=BToKEE_fit_k_pt[iB] /BToKEE_fit_mass[iB];
  BToKEE_k_DCASig_=ProbeTracks_DCASig[k_idx];
  BToKEE_k_dxy_sig_=ProbeTracks_dxyS[k_idx];
  BToKEE_fit_normpt_=BToKEE_fit_pt[iB] /BToKEE_fit_mass[iB];
  BToKEE_l_xy_sig_ = (BToKEE_l_xy[iB]) /BToKEE_l_xy_unc[iB];
  BToKEE_eleDR_= DeltaR(ele1_eta,ele1_phi,ele2_eta,ele2_phi);
  TLorentzVector dll=ele1TLV+ele2TLV;
  
  BToKEE_llkDR_=dll.DeltaR(kTLV);
  BToKEE_l1_iso04_rel_=BToKEE_l1_iso04[iB]/BToKEE_fit_l1_pt[iB];
  BToKEE_l2_iso04_rel_=BToKEE_l2_iso04[iB]/BToKEE_fit_l2_pt[iB];
  BToKEE_k_iso04_rel_ = BToKEE_k_iso04[iB] / BToKEE_fit_k_pt[iB];
  BToKEE_b_iso04_rel_ =BToKEE_b_iso04[iB]/BToKEE_fit_pt[iB];
  TVector3 diele_p3 = dll.Vect();
  TVector3 k_p3 = kTLV.Vect();
  TVector3 pv2sv_p3(PV_x-BToKEE_vtx_x[iB], PV_y-BToKEE_vtx_y[iB], PV_z-BToKEE_vtx_z[iB]);
  
  BToKEE_ptAsym_ = ( (diele_p3.Cross(pv2sv_p3)).Mag() - (k_p3.Cross(pv2sv_p3)).Mag() ) 
    / ( (diele_p3.Cross(pv2sv_p3)).Mag() + (k_p3.Cross(pv2sv_p3)).Mag() );
  
  BToKEE_l1_dzTrg_=Electron_dzTrg[ele1_idx];
  BToKEE_l2_dzTrg_=Electron_dzTrg[ele2_idx];
  BToKEE_k_dzTrg_=ProbeTracks_dzTrg[k_idx];
  float valpfmvaidinit=20;
  if(_isMC) valpfmvaidinit=20;
  BToKEE_l1_pfmvaId_lowPt_=Electron_pfmvaId[ele1_idx];
  BToKEE_l1_pfmvaId_lowPt_=valpfmvaidinit;
  
  BToKEE_l2_pfmvaId_lowPt_=Electron_pfmvaId[ele2_idx];
  BToKEE_l2_pfmvaId_lowPt_=valpfmvaidinit;
  
  BToKEE_l1_pfmvaId_highPt_=Electron_pfmvaId[ele1_idx];
  BToKEE_l1_pfmvaId_highPt_=valpfmvaidinit;
  
  BToKEE_l2_pfmvaId_highPt_=Electron_pfmvaId[ele2_idx];
  BToKEE_l2_pfmvaId_highPt_=valpfmvaidinit;

  // new vars from bparking
  BToKEE_b_iso03_ = BToKEE_b_iso03[iB];
  BToKEE_b_iso03_dca_ = BToKEE_b_iso03_dca[iB];
  BToKEE_b_iso03_dca_tight_ = BToKEE_b_iso03_dca_tight[iB];
  BToKEE_b_iso04_dca_ = BToKEE_b_iso04_dca[iB];
  BToKEE_b_iso04_dca_tight_ = BToKEE_b_iso04_dca_tight[iB];
  BToKEE_fit_eta_ =BToKEE_fit_eta[iB];
  //BToKEE_fit_k_eta_ =BToKEE_fit_k_eta[iB];
  //BToKEE_fit_k_phi_ = BToKEE_fit_k_phi[iB];
  BToKEE_fit_phi_ = BToKEE_fit_phi[iB];
  BToKEE_k_iso03_ =BToKEE_k_iso03[iB];
  BToKEE_k_iso03_dca_ =BToKEE_k_iso03_dca[iB];
  BToKEE_k_iso03_dca_tight_ = BToKEE_k_iso03_dca_tight[iB];
  BToKEE_k_iso04_dca_ =BToKEE_k_iso04_dca[iB];
  BToKEE_k_iso04_dca_tight_ = BToKEE_k_iso04_dca_tight[iB];
  BToKEE_k_svip2d_ = BToKEE_k_svip2d[iB];
  BToKEE_k_svip2d_err_ = BToKEE_k_svip2d_err[iB];
  BToKEE_k_svip3d_err_ = BToKEE_k_svip3d_err[iB];
  BToKEE_l1_iso03_dca_ = BToKEE_l1_iso03_dca[iB];
  BToKEE_l1_iso03_dca_tight_ = BToKEE_l1_iso03_dca_tight[iB];
  BToKEE_l1_iso04_ = BToKEE_l1_iso04[iB];
  BToKEE_l1_iso04_dca_ = BToKEE_l1_iso04_dca[iB];
  BToKEE_l1_iso04_dca_tight_ = BToKEE_l1_iso04_dca_tight[iB];
  BToKEE_l2_iso03_ = BToKEE_l2_iso03[iB];
  BToKEE_l2_iso03_dca_ = BToKEE_l2_iso03_dca[iB];
  BToKEE_l2_iso03_dca_tight_ = BToKEE_l2_iso03_dca_tight[iB];
  BToKEE_l2_iso04_ = BToKEE_l2_iso04[iB];
  BToKEE_l2_iso04_dca_ = BToKEE_l2_iso04_dca[iB];
  BToKEE_l2_iso04_dca_tight_ = BToKEE_l2_iso04_dca_tight[iB];
  BToKEE_maxDR_ = BToKEE_maxDR[iB];
  BToKEE_minDR_ = BToKEE_minDR[iB];
  BToKEE_b_n_isotrk_ = BToKEE_b_n_isotrk[iB];
  BToKEE_b_n_isotrk_dca_ = BToKEE_b_n_isotrk_dca[iB];
  BToKEE_b_n_isotrk_dca_tight_ = BToKEE_b_n_isotrk_dca_tight[iB];
  BToKEE_k_n_isotrk_ = BToKEE_k_n_isotrk[iB]; 
  BToKEE_k_n_isotrk_dca_ =   BToKEE_k_n_isotrk_dca[iB]; 
  BToKEE_k_n_isotrk_dca_tight_ =   BToKEE_k_n_isotrk_dca_tight[iB]; 
  BToKEE_l1_n_isotrk_ =   BToKEE_l1_n_isotrk[iB]; 
  BToKEE_l1_n_isotrk_dca_ =   BToKEE_l1_n_isotrk_dca[iB]; 
  BToKEE_l1_n_isotrk_dca_tight_ =   BToKEE_l1_n_isotrk_dca_tight[iB]; 
  BToKEE_l2_n_isotrk_ =    BToKEE_l2_n_isotrk[iB]; 
  BToKEE_l2_n_isotrk_dca_=   BToKEE_l2_n_isotrk_dca[iB]; 
  BToKEE_l2_n_isotrk_dca_tight_ =  BToKEE_l2_n_isotrk_dca_tight[iB];
  Electron_fBrem_l1_ = Electron_fBrem[ele1_idx];
  Electron_fBrem_l2_ = Electron_fBrem[ele2_idx];
  Electron_ip3d_l1_ = Electron_ip3d[ele1_idx];
  Electron_ip3d_l2_ = Electron_ip3d[ele2_idx];
  Electron_pfRelIso_l1_ = Electron_pfRelIso[ele1_idx];
  Electron_pfRelIso_l2_ = Electron_pfRelIso[ele2_idx];
  Electron_sip3d_l1_ = Electron_sip3d[ele1_idx];
  Electron_sip3d_l2_ = Electron_sip3d[ele2_idx];
  Electron_trkRelIso_l1_ = Electron_trkRelIso[ele1_idx];
  Electron_trkRelIso_l2_ = Electron_trkRelIso[ele2_idx];
  ProbeTracks_dzS_ =ProbeTracks_dzS[k_idx];
  ProbeTracks_eta_ = ProbeTracks_eta[k_idx];
  ProbeTracks_nValidHits_ = ProbeTracks_nValidHits[k_idx];




	// Mark if matched and break loop
	isMatched_ = 1;
	h_cand_->Fill(iB<1000?iB:1000,1.); // overflows into final bin
	break;
      }
    }

//    std::cout << "TEST jsonNames_.size(): " << jsonNames_.size()
//	      << " jsonFlags_.size(): " << jsonFlags_.size()
//	      << " flags: ";
//    for ( auto const& flag : jsonFlags_ ) { std::cout << flag << " "; }
//    std::cout << std::endl;

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
Efficiency::Efficiency(TChain* tree, int isMC, int mode, std::string output) :
  EfficiencyBase((TTree*)tree, isMC),
  _mode(mode),
  jsonFilters_(),
  jsonNames_(),
  jsonFlags_(),
  jsonPtrs_(),
  tmp0(0),
  tmp1(0),
  tmp2(0),
  tmp3(0),
  tmp4(0),
  tmp5(0),
  tmp6(0),
  tmp7(0),
  tmp8(0),
  tmp9(0),
  tmp10(0),
  tmp11(0),
  tmp12(0),
  tmp13(0),
  tmp14(0),
  tmp15(0)
{

  // Obtain JSON file paths (hacked for now!)
  std::string json_path = "/eos/cms/store/group/phys_bphys/DiElectronX/test/trigger/JSON/Eras_CDEF_to361091/";
  std::vector<std::string> json_files;
  if (!_isMC) {
    json_files.push_back(json_path+"trigger_OR.json");
    json_files.push_back(json_path+"L1_10p5_HLT_5p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_10p5_HLT_6p5_Excl_Final.json");
    json_files.push_back(json_path+"L1_11p0_HLT_6p5_Excl_Final.json");
    json_files.push_back(json_path+"L1_4p5_HLT_4p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_5p0_HLT_4p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_5p5_HLT_4p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_5p5_HLT_6p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_6p0_HLT_4p0_Excl_Final.json");
    json_files.push_back(json_path+"L1_6p5_HLT_4p5_Excl_Final.json");
    json_files.push_back(json_path+"L1_7p0_HLT_5p0_Excl_Final.json") ;
    json_files.push_back(json_path+"L1_7p5_HLT_5p0_Excl_Final.json") ;
    json_files.push_back(json_path+"L1_8p0_HLT_5p0_Excl_Final.json") ;
    json_files.push_back(json_path+"L1_8p5_HLT_5p0_Excl_Final.json") ;
    json_files.push_back(json_path+"L1_8p5_HLT_5p5_Excl_Final.json") ;
    json_files.push_back(json_path+"L1_9p0_HLT_6p0_Excl_Final.json"); 
  }

  // Initialise JSON filters, etc
  if (verbose_>0) {
    std::cout << "[Efficiency::Efficiency]"
	      << " Found " << json_files.size()
	      << " JSON file paths..." << std::endl;
  }
  for ( auto const& file : json_files ) {
    if (verbose_>0) {
      std::cout << "  [Efficiency::Efficiency]"
		<< " Parsing JSON file path: " << file << std::endl;
    }
    JsonFilter filter(file,verbose_);
    filter.fillRunLSMap();
    jsonFilters_.push_back(filter);
    std::string name = JsonFilter::jsonFileName( filter.jsonFilePath() );
    jsonNames_.push_back(name);
    if (verbose_>0) {
      std::cout << "  [Efficiency::Efficiency]"
		<< "Parsed JSON file name: " << name << std::endl;
    }
  }
//  std::vector<bool>::iterator iter = jsonFlags_.begin();
//  for ( auto const& name : jsonNames_ ) {
//    std::advance(iter,1);
//    jsonPtrs_.push_back(&*iter);
//  }
  jsonFlags_.clear();
  jsonFlags_.resize(jsonNames_.size(),false);
//  for ( auto& flag : jsonFlags_ ) {
//    jsonPtrs_.push_back(&flag);
//  }

  // Output trees (after JSON file parsing!)
  prepareOutputs(output);
  
  // Initialise vars (after JSON file parsing!)
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
  std::cout << "Processed event " << jentry << " out of " << nentries << " (" << int(100.*jentry/nentries) << "%). "
	    << "Time elapsed: " << elapsed_seconds.count() << " s. "
	    << "Time remaining: " << (predicted_duration.count() - elapsed_seconds.count() ) << " s. "
	    << "Predicted finish time: " << tmp
	    << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
//
void Efficiency::initVars() {

  // JSON flags
  jsonFlags_.clear();
  jsonFlags_.resize(jsonNames_.size(),false);
  tmp0=0;
  tmp1=0;
  tmp2=0;
  tmp3=0;
  tmp4=0;
  tmp5=0;
  tmp6=0;
  tmp7=0;
  tmp8=0;
  tmp9=0;
  tmp10=0;
  tmp11=0;
  tmp12=0;
  tmp13=0;
  tmp14=0;
  tmp15=0;

  // Scalars
  theRun_=0;
  theLumi_=0;
  theEvent_=0;
  nvtx_=0;
  Rho_fixedGridRhoAll_ = -1000.;
  Rho_fixedGridRhoFastjetAll_ = -1000.;
  Rho_fixedGridRhoFastjetCentral_ = -1000.;
  Rho_fixedGridRhoFastjetCentralCalo_ = -1000.;
  Rho_fixedGridRhoFastjetCentralChargedPileUp_ = -1000.;
  Rho_fixedGridRhoFastjetCentralNeutral_ = -1000.;

  // Trigger
//  trg_muon_pt_=0.;
//  trg_muon_eta_=-10.;
//  hlt7_ip4_=0;
//  hlt8_ip3_=0;
//  hlt9_ip6_=0;
//  hlt12_ip6_=0;

  hlt_10p0_=0;
  hlt_9p5_=0;
  hlt_9p0_=0;
  hlt_8p5_=0;
  hlt_8p0_=0;
  hlt_7p5_=0;
  hlt_7p0_=0;
  hlt_6p5_=0;
  hlt_6p0_=0;
  hlt_5p5_=0;
  hlt_5p0_=0;
  hlt_4p5_=0;
  hlt_4p0_=0;
  
  l1_11p0_=0;
  l1_10p5_=0;
  l1_10p0_=0;
  l1_9p5_=0;
  l1_9p0_=0;
  l1_8p5_=0;
  l1_8p0_=0;
  l1_7p5_=0;
  l1_7p0_=0;
  l1_6p5_=0;
  l1_6p0_=0;
  l1_5p5_=0;
  l1_5p0_=0;
  l1_4p5_=0;
  l1_4p0_=0;
  
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

  // RECO ID
  e1_reco_loose_=0;
  e1_reco_medium_=0;
  e1_reco_tight_=0;
  e2_reco_loose_=0;
  e2_reco_medium_=0;
  e2_reco_tight_=0;

  // Pre-selection and BDT
  ip3d_=-10.;
  cos2d_=-10.;
  bdt_=-10.;
  mll_=0.;

  // B candidates
  b_mass_ = -10.;
  b_mass_err_ = -10.;
  b_pt_ = -10.;
  b_l1_pt_ = -10.;
  b_l2_pt_ = -10.;
  b_k_pt_ = -10.;
  b_cos2D_ = -10.;
  b_lxy_ = -10.;
  b_lxyerr_ = -10.;
  b_svprob_ = -10.;

  // bdt vars
  BToKEE_fit_l1_normpt_ = -1000.;
  BToKEE_fit_l2_normpt_ = -1000.;
  BToKEE_l1_dxy_sig_ = -1000.;
  BToKEE_l2_dxy_sig_ = -1000.;
  BToKEE_fit_k_normpt_ = -1000.;
  BToKEE_k_DCASig_ = -1000.;
  BToKEE_k_dxy_sig_ = -1000.;
  BToKEE_fit_normpt_ = -1000.;
  BToKEE_l_xy_sig_ = -1000.;
  BToKEE_eleDR_ = -1000.;
  BToKEE_llkDR_ = -1000.;
  BToKEE_l1_iso04_rel_ = -1000.;
  BToKEE_l2_iso04_rel_ = -1000.;
  BToKEE_k_iso04_rel_ = -1000.;
  BToKEE_b_iso04_rel_ = -1000.;
  BToKEE_ptAsym_ = -1000.;
  BToKEE_l1_dzTrg_ = -1000.;
  BToKEE_l2_dzTrg_ = -1000.;
  BToKEE_k_dzTrg_ = -1000.;
  BToKEE_l1_pfmvaId_lowPt_ = -1000.;
  BToKEE_l2_pfmvaId_lowPt_ = -1000.;
  BToKEE_l1_pfmvaId_highPt_ = -1000.;
  BToKEE_l2_pfmvaId_highPt_ = -1000.;

  //new vars

  BToKEE_b_iso03_= -1000.;
  BToKEE_b_iso03_dca_= -1000.;
  BToKEE_b_iso03_dca_tight_= -1000.;
  BToKEE_b_iso04_dca_ = -1000.;
  BToKEE_b_iso04_dca_tight_ = -1000.;
  BToKEE_fit_eta_ = -1000.;
  //BToKEE_fit_k_eta_ = -1000.; ///doddgydydgdygdgdyd
  //BToKEE_fit_k_phi_ = -1000.;
//  //
  BToKEE_fit_phi_ = -1000.;
  BToKEE_k_iso03_ = -1000.;
  BToKEE_k_iso03_dca_ = -1000.;
  BToKEE_k_iso03_dca_tight_ = -1000.;
  BToKEE_k_iso04_dca_ = -1000.;
  BToKEE_k_iso04_dca_tight_ = -1000.;
  BToKEE_k_svip2d_ = -1000.;
  BToKEE_k_svip2d_err_= -1000.;
  BToKEE_k_svip3d_err_= -1000.;
  BToKEE_l1_iso03_dca_ = -1000.;
  BToKEE_l1_iso03_dca_tight_ = -1000.;
  BToKEE_l1_iso04_ = -1000.;
  BToKEE_l1_iso04_dca_ = -1000.;
  BToKEE_l1_iso04_dca_tight_ = -1000.;
  BToKEE_l2_iso03_= -1000.;
  BToKEE_l2_iso03_dca_ = -1000.;
  BToKEE_l2_iso03_dca_tight_ = -1000.;
  BToKEE_l2_iso04_ = -1000.;
  BToKEE_l2_iso04_dca_ = -1000.;
  BToKEE_l2_iso04_dca_tight_= -1000.;
  BToKEE_maxDR_ = -1000.;
  BToKEE_minDR_ = -1000.;
  BToKEE_b_n_isotrk_ = -1000.;
  BToKEE_b_n_isotrk_dca_ = -1000.;
  BToKEE_b_n_isotrk_dca_tight_ = -1000.;
  BToKEE_k_n_isotrk_ = -1000.;
  BToKEE_k_n_isotrk_dca_= -1000.;
  BToKEE_k_n_isotrk_dca_tight_ = -1000.;
  BToKEE_l1_n_isotrk_ = -1000.;
  BToKEE_l1_n_isotrk_dca_ = -1000.;
  BToKEE_l1_n_isotrk_dca_tight_ = -1000.;
  BToKEE_l2_n_isotrk_ = -1000.;
  BToKEE_l2_n_isotrk_dca_= -1000.;
  BToKEE_l2_n_isotrk_dca_tight_ = -1000.;
  Electron_fBrem_l1_ = -1000.;
  Electron_fBrem_l2_ = -1000.;
  Electron_ip3d_l1_ = -1000.;
  Electron_ip3d_l2_ = -1000.;
  Electron_pfRelIso_l1_ = -1000.;
  Electron_pfRelIso_l2_ = -1000.;
  Electron_sip3d_l1_= -1000.;
  Electron_sip3d_l2_ = -1000.;
  Electron_trkRelIso_l1_ = -1000.;
  Electron_trkRelIso_l2_ = -1000.;
  ProbeTracks_dzS_ = -1000.;
  ProbeTracks_eta_ = -1000.;
  ProbeTracks_nValidHits_ = -1000.;
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
  outTree_->Branch("theLumi", &theLumi_, "theLumi/I");
  outTree_->Branch("theEvent", &theEvent_, "theEvent/I");

//  int idx = 0;
//  for ( auto const& name : jsonNames_ ) {
//    if (verbose_>1) std::cout << "[Efficiency::bookOutputTree]"
//			      << " Adding branch " << name << std::endl;
//    outTree_->Branch(name.c_str(),
//		     &(jsonFlags_[idx]),
//		     std::string(name+"/O").c_str());
//    idx++;
//  }

//  int idx = 0;
//  for ( auto const& flag : jsonFlags_ ) {
//    std::string name = jsonNames_[idx];
//    if (verbose_>1) std::cout << "[Efficiency::bookOutputTree]"
//			      << " Adding branch " << name << std::endl;
//    outTree_->Branch("JSON_"+name).c_str(),
//		     &flag,
//		     std::string("JSON_"+name+"/O").c_str());
//    idx++;
//  }

  outTree_->Branch("JSON_trigger_OR", &tmp0,  "JSON_trigger_OR/I");
  outTree_->Branch("JSON_L1_10p5_HLT_5p0_Excl_Final",&tmp1 , "JSON_L1_10p5_HLT_5p0_Excl_Final");
  outTree_->Branch("JSON_L1_10p5_HLT_6p5_Excl_Final",&tmp2 , "JSON_L1_10p5_HLT_6p5_Excl_Final");
  outTree_->Branch("JSON_L1_11p0_HLT_6p5_Excl_Final",&tmp3 , "JSON_L1_11p0_HLT_6p5_Excl_Final");
  outTree_->Branch("JSON_L1_4p5_HLT_4p0_Excl_Final", &tmp4,  "JSON_L1_4p5_HLT_4p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_5p0_HLT_4p0_Excl_Final", &tmp5,  "JSON_L1_5p0_HLT_4p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_5p5_HLT_4p0_Excl_Final", &tmp6,  "JSON_L1_5p5_HLT_4p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_5p5_HLT_6p0_Excl_Final", &tmp7,  "JSON_L1_5p5_HLT_6p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_6p0_HLT_4p0_Excl_Final", &tmp8,  "JSON_L1_6p0_HLT_4p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_6p5_HLT_4p5_Excl_Final", &tmp9,  "JSON_L1_6p5_HLT_4p5_Excl_Final/I");
  outTree_->Branch("JSON_L1_7p0_HLT_5p0_Excl_Final", &tmp10,  "JSON_L1_7p0_HLT_5p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_7p5_HLT_5p0_Excl_Final", &tmp11,  "JSON_L1_7p5_HLT_5p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_8p0_HLT_5p0_Excl_Final", &tmp12,  "JSON_L1_8p0_HLT_5p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_8p5_HLT_5p0_Excl_Final", &tmp13,  "JSON_L1_8p5_HLT_5p0_Excl_Final/I");
  outTree_->Branch("JSON_L1_8p5_HLT_5p5_Excl_Final", &tmp14,  "JSON_L1_8p5_HLT_5p5_Excl_Final/I");
  outTree_->Branch("JSON_L1_9p0_HLT_6p0_Excl_Final", &tmp15,  "JSON_L1_9p0_HLT_6p0_Excl_Final/I");
  
  // Trigger
//  outTree_->Branch("trg_muon_pt", &trg_muon_pt_, "trg_muon_pt/F");
//  outTree_->Branch("trg_muon_eta", &trg_muon_eta_, "trg_muon_eta/F");
//  outTree_->Branch("HLT_Mu7_IP4", &hlt7_ip4_, "HLT_Mu7_IP4/I");
//  outTree_->Branch("HLT_Mu8_IP3", &hlt8_ip3_, "HLT_Mu8_IP3/I");
//  outTree_->Branch("HLT_Mu9_IP6", &hlt9_ip6_, "HLT_Mu9_IP6/I");
//  outTree_->Branch("HLT_Mu12_IP6", &hlt12_ip6_, "HLT_Mu12_IP6/I");

  outTree_->Branch("HLT_DoubleEle10p0", &hlt_10p0_, "HLT_DoubleEle10p0/I");
  outTree_->Branch("HLT_DoubleEle9p5", &hlt_9p5_, "HLT_DoubleEle9p5/I");
  outTree_->Branch("HLT_DoubleEle9p0", &hlt_9p0_, "HLT_DoubleEle9p0/I");
  outTree_->Branch("HLT_DoubleEle8p5", &hlt_8p5_, "HLT_DoubleEle8p5/I");
  outTree_->Branch("HLT_DoubleEle8p0", &hlt_8p0_, "HLT_DoubleEle8p0/I");
  outTree_->Branch("HLT_DoubleEle7p5", &hlt_7p5_, "HLT_DoubleEle7p5/I");
  outTree_->Branch("HLT_DoubleEle7p0", &hlt_7p0_, "HLT_DoubleEle7p0/I");
  outTree_->Branch("HLT_DoubleEle6p5", &hlt_6p5_, "HLT_DoubleEle6p5/I");
  outTree_->Branch("HLT_DoubleEle6p0", &hlt_6p0_, "HLT_DoubleEle6p0/I");
  outTree_->Branch("HLT_DoubleEle5p5", &hlt_5p5_, "HLT_DoubleEle5p5/I");
  outTree_->Branch("HLT_DoubleEle5p0", &hlt_5p0_, "HLT_DoubleEle5p0/I");
  outTree_->Branch("HLT_DoubleEle4p5", &hlt_4p5_, "HLT_DoubleEle4p5/I");
  outTree_->Branch("HLT_DoubleEle4p0", &hlt_4p0_, "HLT_DoubleEle4p0/I");

  outTree_->Branch("L1_DoubleEG11p0", &l1_11p0_, "L1_DoubleEG11p0/I");
  outTree_->Branch("L1_DoubleEG10p5", &l1_10p5_, "L1_DoubleEG10p5/I");
  outTree_->Branch("L1_DoubleEG10p0", &l1_10p0_, "L1_DoubleEG10p0/I");
  outTree_->Branch("L1_DoubleEG9p5", &l1_9p5_, "L1_DoubleEG9p5/I");
  outTree_->Branch("L1_DoubleEG9p0", &l1_9p0_, "L1_DoubleEG9p0/I");
  outTree_->Branch("L1_DoubleEG8p5", &l1_8p5_, "L1_DoubleEG8p5/I");
  outTree_->Branch("L1_DoubleEG8p0", &l1_8p0_, "L1_DoubleEG8p0/I");
  outTree_->Branch("L1_DoubleEG7p5", &l1_7p5_, "L1_DoubleEG7p5/I");
  outTree_->Branch("L1_DoubleEG7p0", &l1_7p0_, "L1_DoubleEG7p0/I");
  outTree_->Branch("L1_DoubleEG6p5", &l1_6p5_, "L1_DoubleEG6p5/I");
  outTree_->Branch("L1_DoubleEG6p0", &l1_6p0_, "L1_DoubleEG6p0/I");
  outTree_->Branch("L1_DoubleEG5p5", &l1_5p5_, "L1_DoubleEG5p5/I");
  outTree_->Branch("L1_DoubleEG5p0", &l1_5p0_, "L1_DoubleEG5p0/I");
  outTree_->Branch("L1_DoubleEG4p5", &l1_4p5_, "L1_DoubleEG4p5/I");
  outTree_->Branch("L1_DoubleEG4p0", &l1_4p0_, "L1_DoubleEG4p0/I");

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
//  outTree_->Branch("e1_reco_pf", &e1_reco_pf_, "e1_reco_pf/I");
//  outTree_->Branch("e2_reco_pf", &e2_reco_pf_, "e2_reco_pf/I");
//  outTree_->Branch("e1_reco_lowpt", &e1_reco_lowpt_, "e1_reco_lowpt/I");
//  outTree_->Branch("e2_reco_lowpt", &e2_reco_lowpt_, "e2_reco_lowpt/I");
//  outTree_->Branch("e1_reco_overlap", &e1_reco_overlap_, "e1_reco_overlap/I");
//  outTree_->Branch("e2_reco_overlap", &e2_reco_overlap_, "e2_reco_overlap/I");

  //RECO ID
  outTree_->Branch("e1_reco_loose", &e1_reco_loose_, "e1_reco_loose/I");
  outTree_->Branch("e2_reco_loose", &e2_reco_loose_, "e2_reco_loose/I");
  outTree_->Branch("e1_reco_medium", &e1_reco_medium_, "e1_reco_medium/I");
  outTree_->Branch("e2_reco_medium", &e2_reco_medium_, "e2_reco_medium/I");
  outTree_->Branch("e1_reco_tight", &e1_reco_tight_, "e1_reco_tight/I");
  outTree_->Branch("e2_reco_tight", &e2_reco_tight_, "e2_reco_tight/I");

  // Pre-selection and BDT
  outTree_->Branch("ip3d", &ip3d_, "ip3d/F");
  outTree_->Branch("cos2d", &cos2d_, "cos2d/F");
  outTree_->Branch("bdt", &bdt_, "bdt/F");
  outTree_->Branch("mll", &mll_, "mll/F");

  // B candidates
  outTree_->Branch("b_mass", &b_mass_, "b_mass/F");
  outTree_->Branch("b_mass_err", &b_mass_err_, "b_mass_err/F");
  outTree_->Branch("b_pt", &b_pt_, "b_pt/F");
  outTree_->Branch("b_l1_pt", &b_l1_pt_, "b_l1_pt/F");
  outTree_->Branch("b_l2_pt", &b_l2_pt_, "b_l2_pt/F");
  outTree_->Branch("b_k_pt", &b_k_pt_, "b_k_pt/F");
  outTree_->Branch("b_cos2D", &b_cos2D_, "b_cos2D/F");
  outTree_->Branch("b_lxy", &b_lxy_, "b_lxy/F");
  outTree_->Branch("b_lxyerr", &b_lxyerr_, "b_lxyerr/F");
  outTree_->Branch("b_svprob", &b_svprob_, "b_svprob/F");
  //bdt vars
  outTree_->Branch("BToKEE_fit_l1_normpt", &BToKEE_fit_l1_normpt_, "BToKEE_fit_l1_normpt/F");
  outTree_->Branch("BToKEE_fit_l2_normpt", &BToKEE_fit_l2_normpt_, "BToKEE_fit_l2_normpt/F");
  outTree_->Branch("BToKEE_l1_dxy_sig", &BToKEE_l1_dxy_sig_, "BToKEE_l1_dxy_sig/F");
  outTree_->Branch("BToKEE_l2_dxy_sig", &BToKEE_l2_dxy_sig_, "BToKEE_l2_dxy_sig/F");
  outTree_->Branch("BToKEE_fit_k_normpt", &BToKEE_fit_k_normpt_, "BToKEE_fit_k_normpt/F");
  outTree_->Branch("BToKEE_k_DCASig", &BToKEE_k_DCASig_, "BToKEE_k_DCASig/F");
  outTree_->Branch("BToKEE_k_dxy_sig", &BToKEE_k_dxy_sig_, "BToKEE_k_dxy_sig/F");
  outTree_->Branch("BToKEE_fit_normpt", &BToKEE_fit_normpt_, "BToKEE_fit_normpt/F");
  outTree_->Branch("BToKEE_l_xy_sig", &BToKEE_l_xy_sig_, "BToKEE_l_xy_sig/F");
  outTree_->Branch("BToKEE_eleDR", &BToKEE_eleDR_, "BToKEE_eleDR/F");
  outTree_->Branch("BToKEE_llkDR", &BToKEE_llkDR_, "BToKEE_llkDR/F");
  outTree_->Branch("BToKEE_l1_iso04_rel", &BToKEE_l1_iso04_rel_, "BToKEE_l1_iso04_rel/F");
  outTree_->Branch("BToKEE_l2_iso04_rel", &BToKEE_l2_iso04_rel_, "BToKEE_l2_iso04_rel/F");
  outTree_->Branch("BToKEE_k_iso04_rel", &BToKEE_k_iso04_rel_, "BToKEE_k_iso04_rel/F");
  outTree_->Branch("BToKEE_b_iso04_rel", &BToKEE_b_iso04_rel_, "BToKEE_b_iso04_rel/F");
  outTree_->Branch("BToKEE_ptAsym", &BToKEE_ptAsym_, "BToKEE_ptAsym/F");
  outTree_->Branch("BToKEE_l1_dzTrg", &BToKEE_l1_dzTrg_, "BToKEE_l1_dzTrg/F");
  outTree_->Branch("BToKEE_l2_dzTrg", &BToKEE_l2_dzTrg_, "BToKEE_l2_dzTrg/F");
  outTree_->Branch("BToKEE_k_dzTrg", &BToKEE_k_dzTrg_, "BToKEE_k_dzTrg/F");
  outTree_->Branch("BToKEE_l1_pfmvaId_lowPt", &BToKEE_l1_pfmvaId_lowPt_, "BToKEE_l1_pfmvaId_lowPt/F");
  outTree_->Branch("BToKEE_l2_pfmvaId_lowPt", &BToKEE_l2_pfmvaId_lowPt_, "BToKEE_l2_pfmvaId_lowPt/F");
  outTree_->Branch("BToKEE_l1_pfmvaId_highPt", &BToKEE_l1_pfmvaId_highPt_, "BToKEE_l1_pfmvaId_highPt/F");
  outTree_->Branch("BToKEE_l2_pfmvaId_highPt", &BToKEE_l2_pfmvaId_highPt_, "BToKEE_l2_pfmvaId_highPt/F");
  outTree_->Branch("BToKEE_b_iso03",            &BToKEE_b_iso03_, "BToKEE_b_iso03/F");
  outTree_->Branch("BToKEE_b_iso03_dca",            &BToKEE_b_iso03_dca_, "BToKEE_b_iso03_dca/F");
  outTree_->Branch("BToKEE_b_iso03_dca_tight",            &BToKEE_b_iso03_dca_tight_, "BToKEE_b_iso03_dca_tight/F");
  outTree_->Branch("BToKEE_b_iso04_dca",            &BToKEE_b_iso04_dca_ , "BToKEE_b_iso04_dca_/F");
  outTree_->Branch("BToKEE_b_iso04_dca_tight",            &BToKEE_b_iso04_dca_tight_ , " BToKEE_b_iso04_dca_tight_/F");
  outTree_->Branch("BToKEE_fit_eta",            &BToKEE_fit_eta_ , "BToKEE_fit_eta_/F");
  //outTree_->Branch("BToKEE_fit_k_eta",            &BToKEE_fit_k_eta_ , "BToKEE_fit_k_eta_/F"); ///error
  //outTree_->Branch("BToKEE_fit_k_phi",            &BToKEE_fit_k_phi_ , "BToKEE_fit_k_phi_/F"); ///error
  outTree_->Branch("BToKEE_fit_phi",            &BToKEE_fit_phi_ , "BToKEE_fit_phi_/F");
  outTree_->Branch("BToKEE_k_iso03",            &BToKEE_k_iso03_ , "BToKEE_k_iso03_/F");
  outTree_->Branch("BToKEE_k_iso03_dca",            &BToKEE_k_iso03_dca_ , "BToKEE_k_iso03_dca_/F");
  outTree_->Branch("BToKEE_k_iso03_dca_tight",            &BToKEE_k_iso03_dca_tight_ , "BToKEE_k_iso03_dca_tight_/F");
  outTree_->Branch("BToKEE_k_iso04_dca",            &BToKEE_k_iso04_dca_ , "BToKEE_k_iso04_dca_/F");
  outTree_->Branch("BToKEE_k_iso04_dca_tight",            &BToKEE_k_iso04_dca_tight_ , "BToKEE_k_iso04_dca_tight_/F");
  outTree_->Branch("BToKEE_k_svip2d",            &BToKEE_k_svip2d_ , "BToKEE_k_svip2d_/F");
  outTree_->Branch("BToKEE_k_svip2d_err",            &BToKEE_k_svip2d_err_, "BToKEE_k_svip2d_err/F");
  outTree_->Branch("BToKEE_k_svip3d_err",            &BToKEE_k_svip3d_err_, "BToKEE_k_svip3d_err/F");
  outTree_->Branch("BToKEE_l1_iso03_dca",            &BToKEE_l1_iso03_dca_ , "BToKEE_l1_iso03_dca_/F");
  outTree_->Branch("BToKEE_l1_iso03_dca_tight",            &BToKEE_l1_iso03_dca_tight_ , "BToKEE_l1_iso03_dca_tight_/F");
  outTree_->Branch("BToKEE_l1_iso04",            &BToKEE_l1_iso04_ , "BToKEE_l1_iso04_/F");
  outTree_->Branch("BToKEE_l1_iso04_dca",            &BToKEE_l1_iso04_dca_ , "BToKEE_l1_iso04_dca_/F");
  outTree_->Branch("BToKEE_l1_iso04_dca_tight",            &BToKEE_l1_iso04_dca_tight_ , "BToKEE_l1_iso04_dca_tight_/F");
  outTree_->Branch("BToKEE_l2_iso03",            &BToKEE_l2_iso03_, "BToKEE_l2_iso03/F");
  outTree_->Branch("BToKEE_l2_iso03_dca",            &BToKEE_l2_iso03_dca_ , "BToKEE_l2_iso03_dca_/F");
  outTree_->Branch("BToKEE_l2_iso03_dca_tight",            &BToKEE_l2_iso03_dca_tight_ , "BToKEE_l2_iso03_dca_tight_/F");
  outTree_->Branch("BToKEE_l2_iso04",            &BToKEE_l2_iso04_ , "BToKEE_l2_iso04_/F");
  outTree_->Branch("BToKEE_l2_iso04_dca",            &BToKEE_l2_iso04_dca_ , "BToKEE_l2_iso04_dca_/F");
  outTree_->Branch("BToKEE_l2_iso04_dca_tight",            &BToKEE_l2_iso04_dca_tight_, "BToKEE_l2_iso04_dca_tight/F");
  outTree_->Branch("BToKEE_maxDR",            &BToKEE_maxDR_ , "BToKEE_maxDR_/F");
  outTree_->Branch("BToKEE_minDR",            &BToKEE_minDR_ , "BToKEE_minDR_/F");
  outTree_->Branch("BToKEE_b_n_isotrk",            &BToKEE_b_n_isotrk_ , "  BToKEE_b_n_isotrk_/F");
  outTree_->Branch("BToKEE_b_n_isotrk_dca",            &BToKEE_b_n_isotrk_dca_ , "BToKEE_b_n_isotrk_dca_/F");
  outTree_->Branch("BToKEE_b_n_isotrk_dca_tight",            &BToKEE_b_n_isotrk_dca_tight_ , "BToKEE_b_n_isotrk_dca_tight_/F");
  outTree_->Branch("BToKEE_k_n_isotrk",            &BToKEE_k_n_isotrk_ , "BToKEE_k_n_isotrk_/F");
  outTree_->Branch("BToKEE_k_n_isotrk_dca",            &BToKEE_k_n_isotrk_dca_, "BToKEE_k_n_isotrk_dca/F");
  outTree_->Branch("BToKEE_k_n_isotrk_dca_tight",            &BToKEE_k_n_isotrk_dca_tight_ , "BToKEE_k_n_isotrk_dca_tight_/F");
  outTree_->Branch("BToKEE_l1_n_isotrk",            &BToKEE_l1_n_isotrk_ , "BToKEE_l1_n_isotrk_/F");
  outTree_->Branch("BToKEE_l1_n_isotrk_dca",            &BToKEE_l1_n_isotrk_dca_ , "BToKEE_l1_n_isotrk_dca_/F");
  outTree_->Branch("BToKEE_l1_n_isotrk_dca_tight",            &BToKEE_l1_n_isotrk_dca_tight_ , "BToKEE_l1_n_isotrk_dca_tight_/F");
  outTree_->Branch("BToKEE_l2_n_isotrk",            &BToKEE_l2_n_isotrk_ , "BToKEE_l2_n_isotrk_/F");
  outTree_->Branch("BToKEE_l2_n_isotrk_dca",            &BToKEE_l2_n_isotrk_dca_, "BToKEE_l2_n_isotrk_dca/F");
  outTree_->Branch("BToKEE_l2_n_isotrk_dca_tight",            &BToKEE_l2_n_isotrk_dca_tight_ , "BToKEE_l2_n_isotrk_dca_tight_/F");
  outTree_->Branch("Electron_fBrem_l1",            &Electron_fBrem_l1_ , "Electron_fBrem_l1_/F");
  outTree_->Branch("Electron_fBrem_l2",            &Electron_fBrem_l2_ , "Electron_fBrem_l2_/F");
  outTree_->Branch("Electron_ip3d_l1",            &Electron_ip3d_l1_ , "Electron_ip3d_l1_/F");
  outTree_->Branch("Electron_ip3d_l2",            &Electron_ip3d_l2_ , "Electron_ip3d_l2_/F");
  outTree_->Branch("Electron_pfRelIso_l1",            &Electron_pfRelIso_l1_ , "Electron_pfRelIso_l1_/F");
  outTree_->Branch("Electron_pfRelIso_l2",            &Electron_pfRelIso_l2_ , "Electron_pfRelIso_l2_/F");
  outTree_->Branch("Electron_sip3d_l1",            &Electron_sip3d_l1_, "Electron_sip3d_l1/F");
  outTree_->Branch("Electron_sip3d_l2",            &Electron_sip3d_l2_ , "Electron_sip3d_l2_/F");
  outTree_->Branch("Electron_trkRelIso_l1",            &Electron_trkRelIso_l1_ , "Electron_trkRelIso_l1_/F");
  outTree_->Branch("Electron_trkRelIso_l2",            &Electron_trkRelIso_l2_ , "Electron_trkRelIso_l2_/F");
  outTree_->Branch("Rho_fixedGridRhoAll",            &Rho_fixedGridRhoAll_ , "Rho_fixedGridRhoAll_/F");
  outTree_->Branch("Rho_fixedGridRhoFastjetAll",            &Rho_fixedGridRhoFastjetAll_ , "Rho_fixedGridRhoFastjetAll_/F");
  outTree_->Branch("Rho_fixedGridRhoFastjetCentral",            &Rho_fixedGridRhoFastjetCentral_ , "Rho_fixedGridRhoFastjetCentral_/F");
  outTree_->Branch("Rho_fixedGridRhoFastjetCentralCalo",            &Rho_fixedGridRhoFastjetCentralCalo_ , "Rho_fixedGridRhoFastjetCentralCalo_/F");
  outTree_->Branch("Rho_fixedGridRhoFastjetCentralChargedPileUp",            &Rho_fixedGridRhoFastjetCentralChargedPileUp_ , "Rho_fixedGridRhoFastjetCentralChargedPileUp_/F");
  outTree_->Branch("Rho_fixedGridRhoFastjetCentralNeutral",            &Rho_fixedGridRhoFastjetCentralNeutral_ , "Rho_fixedGridRhoFastjetCentralNeutral_/F");
  outTree_->Branch("ProbeTracks_dzS",            &ProbeTracks_dzS_, "ProbeTracks_dzS/F");
  outTree_->Branch("ProbeTracks_eta",            &ProbeTracks_eta_, "ProbeTracks_eta/F");
  outTree_->Branch("ProbeTracks_nValidHits",            &ProbeTracks_nValidHits_, "ProbeTracks_nValidHits/F");
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
bool Efficiency::isBToKEE(Efficiency::Daughters& daughters) { //, bool resonant) {
  daughters = Daughters();
  Cands cands;
  if (verbose_>3) std::cout << "### isBToKEE: find daughters to Bs for mode " << _mode << "... " << std::endl;
  for ( int idx = 0; idx < nGenPart; ++idx ) {
    int d_idx = idx;
    int d_pdg = GenPart_pdgId[d_idx];
    int m_idx = GenPart_genPartIdxMother[d_idx];
    int m_pdg = GenPart_pdgId[m_idx];
    int g_idx = GenPart_genPartIdxMother[m_idx];
    int g_pdg = GenPart_pdgId[g_idx];
    int gg_idx = GenPart_genPartIdxMother[g_idx];
    int gg_pdg = GenPart_pdgId[gg_idx];
    if (verbose_>3) std::cout << " GEN particle provenance:"
			      << " idx: " << idx
			      << " d_idx: " << d_idx
			      << " (d_pdg: " << d_pdg
			      << ") m_idx: " << m_idx
			      << " (m_pdg: " << m_pdg
			      << ") g_idx: " << g_idx
			      << " (g_pdg: " << g_pdg
			      << ") gg_idx: " << gg_idx
			      << " (gg_pdg: " << gg_pdg
			      << ")"
			      << std::endl;
    if ( abs(m_pdg) == 521 && abs(d_pdg) == 321 ) { // find daughter kaon
      if ( cands.find(m_idx) == cands.end() ) {
	cands[m_idx] = std::make_pair(Ints(),Ints());
      }
      cands[m_idx].first.push_back(d_idx);
      cands[m_idx].second.push_back(d_pdg);
      if (verbose_>3) std::cout << " Kaon provenance:"
				<< " d_idx: " << d_idx
				<< " d_pdg: " << d_pdg
				<< " m_idx: " << m_idx
				<< " m_pdg: " << m_pdg
				<< std::endl;
    }
    if ( _mode == 1 ) {
      if ( abs(m_pdg) == 521 && abs(d_pdg) == 11 ) { // find daughter electrons
	if ( cands.find(m_idx) == cands.end() ) {
	  cands[m_idx] = std::make_pair(Ints(),Ints());
	}
	cands[m_idx].first.push_back(d_idx);
	cands[m_idx].second.push_back(d_pdg);
	if (verbose_>3) std::cout << " Electron provenance:"
				  << " d_idx: " << d_idx
				  << " d_pdg: " << d_pdg
				  << " m_idx: " << m_idx
				  << " m_pdg: " << m_pdg
				  << std::endl;
      }
    } else if ( _mode == 2 ) {
      if ( abs(g_pdg) == 521 && abs(m_pdg) == 443 && abs(d_pdg) == 11 ) { // find daughter eles via J/psi
	if ( cands.find(g_idx) == cands.end() ) {
	  cands[g_idx] = std::make_pair(Ints(),Ints());
	}
	cands[g_idx].first.push_back(d_idx);
	cands[g_idx].second.push_back(d_pdg);
	if (verbose_>3) std::cout << " Electron provenance:"
				  << " d_idx: " << d_idx
				  << " d_pdg: " << d_pdg
				  << " m_idx: " << m_idx
				  << " m_pdg: " << m_pdg
				  << " g_idx: " << g_idx
				  << " g_pdg: " << g_pdg
				  << std::endl;
      }
    } else if ( _mode == 3 ) {
      //std::cout << "B --> K Psi2S --> K ee: not implemented yet!" << std::endl;
      if ( abs(g_pdg) == 521 && abs(m_pdg) == 100443 && abs(d_pdg) == 11 ) { // find daughter eles via Psi(2D) and J/psi
	if ( cands.find(g_idx) == cands.end() ) {
	  cands[g_idx] = std::make_pair(Ints(),Ints());
	}
	cands[g_idx].first.push_back(d_idx);
	cands[g_idx].second.push_back(d_pdg);
	if (verbose_>3) std::cout << " Electron provenance:"
				  << " d_idx: " << d_idx
				  << " d_pdg: " << d_pdg
				  << " m_idx: " << m_idx
				  << " m_pdg: " << m_pdg
				  << " g_idx: " << g_idx
				  << " g_pdg: " << g_pdg
				  << " gg_idx: " << gg_idx
				  << " gg_pdg: " << gg_pdg
				  << std::endl;
      }
    } else {
      std::cout << "Mode not known!" << std::endl;
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
      if (verbose_>3) std::cout << " d_idx: " << d_idx
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
			  int& e1_reco_loose,int& e1_reco_medium,int& e1_reco_tight,
			  int& e2_reco_idx,float& e2_reco_pt,float& e2_reco_eta,
			  int& e2_reco_pf,int& e2_reco_lowpt,int& e2_reco_overlap,
			  int& e2_reco_loose,int& e2_reco_medium,int& e2_reco_tight,
			  float& e12_reco_dr //, bool resonant
			  ) {
  if (_isMC && verbose_>4) std::cout << "Efficiency::recoCand:" 
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

  if (_isMC) {

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

    bool found = false;
    if ( _mode == 1 ) { 
      found = ( ( ele1_genPartIdx >= 0 &&
		  ele2_genPartIdx >= 0 &&
		  kaon_genPartIdx >= 0
		  ) &&
		( abs(ele1_genPdgId) == 11 &&
		  abs(ele2_genPdgId) == 11 &&
		  abs(kaon_genPdgId) == 321 
		  ) &&
		//( ele1_genPdgId * ele1_genPdgId < 0 ) &&
		( abs(ele1_genMotherPdgId) == 521 &&
		  abs(ele2_genMotherPdgId) == 521 &&
		  abs(kaon_genMotherPdgId) == 521
		  ) );
    } else if ( _mode == 2 ) { 
      found = ( ( ele1_genPartIdx >= 0 &&
		  ele2_genPartIdx >= 0 &&
		  kaon_genPartIdx >= 0
		  ) &&
		( abs(ele1_genPdgId) == 11 &&
		  abs(ele2_genPdgId) == 11 &&
		  abs(kaon_genPdgId) == 321 
		  ) &&
		//( ele1_genPdgId * ele1_genPdgId < 0 ) &&
		( abs(ele1_genGMotherPdgId) == 521 &&
		  abs(ele2_genGMotherPdgId) == 521 &&
		  abs(ele1_genMotherPdgId) == 443 && // via J/psi
		  abs(ele2_genMotherPdgId) == 443 && // via J/psi
		  abs(kaon_genMotherPdgId) == 521
		  ) );
    } else if ( _mode == 3 ) {
      //std::cout << "B --> K Psi2S --> K ee: not implemented yet!" << std::endl;
      found = ( ( ele1_genPartIdx >= 0 &&
		  ele2_genPartIdx >= 0 &&
		  kaon_genPartIdx >= 0
		  ) &&
		( abs(ele1_genPdgId) == 11 &&
		  abs(ele2_genPdgId) == 11 &&
		  abs(kaon_genPdgId) == 321 
		  ) &&
		//( ele1_genPdgId * ele1_genPdgId < 0 ) &&
		( abs(ele1_genGMotherPdgId) == 521 &&
		  abs(ele2_genGMotherPdgId) == 521 &&
		  abs(ele1_genMotherPdgId) == 100443 && // via Psi(2S)
		  abs(ele2_genMotherPdgId) == 100443 && // via Psi(2S)
		  abs(kaon_genMotherPdgId) == 521
		  ) );
    } else {
      std::cout << "Mode not known!!" << std::endl;
    }
    
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
	e1_reco_loose = (int)Electron_LooseID[ele1_idx];
	e1_reco_medium = (int)Electron_MediumID[ele1_idx];
	e1_reco_tight = (int)Electron_TightID[ele1_idx];
      } 
      if (e1_gen_idx==ele2_genPartIdx){
	e1_reco_idx = ele2_idx;
	e1_reco_pt  = Electron_pt[ele2_idx];
	e1_reco_eta = Electron_eta[ele2_idx];
	e1_reco_phi = Electron_phi[ele2_idx];
	e1_reco_pf  = Electron_isPF[ele2_idx];
	e1_reco_lowpt = Electron_isLowPt[ele2_idx];
	e1_reco_overlap = Electron_isPFoverlap[ele2_idx];
	e1_reco_loose = (int)Electron_LooseID[ele2_idx];
	e1_reco_medium = (int)Electron_MediumID[ele2_idx];
	e1_reco_tight = (int)Electron_TightID[ele2_idx];
      }
      if (e2_gen_idx==ele1_genPartIdx){
	e2_reco_idx = ele1_idx;
	e2_reco_pt  = Electron_pt[ele1_idx];
	e2_reco_eta = Electron_eta[ele1_idx];
	e2_reco_phi = Electron_phi[ele1_idx];
	e2_reco_pf  = Electron_isPF[ele1_idx];
	e2_reco_lowpt = Electron_isLowPt[ele1_idx]; //@@ was e1_reco_lowpt (bug!)
	e2_reco_overlap = Electron_isPFoverlap[ele1_idx];
	e2_reco_loose = (int)Electron_LooseID[ele1_idx];
	e2_reco_medium = (int)Electron_MediumID[ele1_idx];
	e2_reco_tight = (int)Electron_TightID[ele1_idx];
      }
      if (e2_gen_idx==ele2_genPartIdx){
	e2_reco_idx = ele2_idx;
	e2_reco_pt  = Electron_pt[ele2_idx];
	e2_reco_eta = Electron_eta[ele2_idx];
	e2_reco_phi = Electron_phi[ele2_idx];
	e2_reco_pf  = Electron_isPF[ele2_idx];
	e2_reco_lowpt = Electron_isLowPt[ele2_idx];
	e2_reco_overlap = Electron_isPFoverlap[ele2_idx];
	e2_reco_loose = (int)Electron_LooseID[ele2_idx];
	e2_reco_medium = (int)Electron_MediumID[ele2_idx];
	e2_reco_tight = (int)Electron_TightID[ele2_idx];
      }
      e12_reco_dr = DeltaR(e1_reco_eta,
			   e1_reco_phi,
			   e2_reco_eta,
			   e2_reco_phi);
    }

    if (verbose_>4) std::cout << "recoCand found: " << found <<  std::endl;
    return found;
    
  } else { // is data

    e1_reco_idx = ele1_idx;
    e1_reco_pt  = Electron_pt[ele1_idx];
    e1_reco_eta = Electron_eta[ele1_idx];
    float e1_reco_phi = Electron_phi[ele1_idx];
    e1_reco_pf  = Electron_isPF[ele1_idx];
    e1_reco_lowpt = Electron_isLowPt[ele1_idx];
    e1_reco_overlap = Electron_isPFoverlap[ele1_idx];
    e1_reco_loose = (int)Electron_LooseID[ele1_idx];
    e1_reco_medium = (int)Electron_MediumID[ele1_idx];
    e1_reco_tight = (int)Electron_TightID[ele1_idx];
    
    e2_reco_idx = ele2_idx;
    e2_reco_pt  = Electron_pt[ele2_idx];
    e2_reco_eta = Electron_eta[ele2_idx];
    float e2_reco_phi = Electron_phi[ele2_idx];
    e2_reco_pf  = Electron_isPF[ele2_idx];
    e2_reco_lowpt = Electron_isLowPt[ele2_idx];
    e2_reco_overlap = Electron_isPFoverlap[ele2_idx];
    e2_reco_loose = (int)Electron_LooseID[ele2_idx];
    e2_reco_medium = (int)Electron_MediumID[ele2_idx];
    e2_reco_tight = (int)Electron_TightID[ele2_idx];
    
    e12_reco_dr = DeltaR(e1_reco_eta,
			 e1_reco_phi,
			 e2_reco_eta,
			 e2_reco_phi);

    if (verbose_>4) std::cout << "recoCand found: " << true <<  std::endl;
    return true;
    
  }
  
  std::cerr << "recoCand: shouldn't get here!!!" <<  std::endl;
  return true;

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

