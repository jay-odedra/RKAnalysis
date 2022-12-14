#ifndef Efficiency_h
#define Efficiency_h

#include "../FastForest/include/fastforest.h"
#include "../PhysicsTools/TriggerLuminosity/interface/JsonFilter.h"
#include "EfficiencyBase.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <string>

using namespace std;

class Efficiency : public EfficiencyBase {

public:
  
  Efficiency(TChain* tree, int isMC, int mode, std::string filename);
  virtual ~Efficiency();
  
  void Loop();
  
  void timing(int nentries, int jentry, auto start);
  void initVars();
  void prepareOutputs(std::string filename);
  void bookOutputTree();
  void bookOutputHistos();

  // Decay and acceptance
  typedef std::vector<int> Ints;
  typedef std::pair<Ints,Ints> Daughters; // pair< vector<d_idx>, vector<d_pdgid> >
  typedef std::map<int,Daughters> Cands; // map< m_idx, Daughters >
  bool isBToKEE(Daughters& daughters); //, bool resonant = false);
  bool inAcceptance(std::vector<int>& indices);
  bool setElePtEta(Daughters& daughters,
		   int& e1_gen_idx,
		   float& e1_gen_pt,
		   float& e1_gen_eta,
		   int& e2_gen_idx,
		   float& e2_gen_pt,
		   float& e2_gen_eta,
		   float& e12_gen_dr);
  bool recoCand(int theB,
		int& e1_gen_idx,int& e2_gen_idx,
		int& e1_reco_idx,float& e1_reco_pt,float& e1_reco_eta,
		int& e1_reco_pf,int& e1_reco_lowpt,int& e1_reco_overlap,
		int& e1_reco_loose,int& e1_reco_medium,int& e1_reco_tight,
		int& e2_reco_idx,float& e2_reco_pt,float& e2_reco_eta,
		int& e2_reco_pf,int& e2_reco_lowpt,int& e2_reco_overlap,
		int& e2_reco_loose,int& e2_reco_medium,int& e2_reco_tight,
		float& e12_reco_dr); //,bool resonant = false);
  float DeltaR(float eta1,
	       float phi1,
	       float eta2,
	       float phi2);
  void loadModels(fastforest::FastForest& fastForestOttoPFPF,
		  fastforest::FastForest& fastForestOttoPFLP);
  float evaluateModels(u_int thisB,
		       fastforest::FastForest& fastForestOttoPFPF,
		       fastforest::FastForest& fastForestOttoPFLP);
  
private:

  int _mode;
  
  // Misc
  int verbose_=0;
  std::string filename_="";

  // JSON files
  std::vector<JsonFilter> jsonFilters_;
  std::vector<std::string> jsonNames_;
  std::vector<bool> jsonFlags_;
  std::vector<bool*> jsonPtrs_;
  int tmp0;
  int tmp1;
  int tmp2;
  int tmp3;
  int tmp4;
  int tmp5;
  int tmp6;
  int tmp7;
  int tmp8;
  int tmp9;
  int tmp10;
  int tmp11;
  int tmp12;
  int tmp13;
  int tmp14;
  int tmp15;

  // Output
  TFile* outFile_=nullptr;
  TTree* outTree_=nullptr;
  TH1F* h_entries_=nullptr;
  TH1F* h_selection_=nullptr;
  TH1F* h_cand_=nullptr;

  // Scalars
  int theRun_;
  int theLumi_;
  int theEvent_;
  int nvtx_;
  float Rho_fixedGridRhoAll_ ;
  float Rho_fixedGridRhoFastjetAll_ ;
  float Rho_fixedGridRhoFastjetCentral_ ;
  float Rho_fixedGridRhoFastjetCentralCalo_ ;
  float Rho_fixedGridRhoFastjetCentralChargedPileUp_ ;
  float Rho_fixedGridRhoFastjetCentralNeutral_ ;

  // Trigger
//  float trg_muon_pt_;
//  float trg_muon_eta_;
//  int hlt7_ip4_;
//  int hlt8_ip3_;
//  int hlt9_ip6_;
//  int hlt12_ip6_;

  int hlt_10p0_;
  int hlt_9p5_;
  int hlt_9p0_;
  int hlt_8p5_;
  int hlt_8p0_;
  int hlt_7p5_;
  int hlt_7p0_;
  int hlt_6p5_;
  int hlt_6p0_;
  int hlt_5p5_;
  int hlt_5p0_;
  int hlt_4p5_;
  int hlt_4p0_;

  int l1_11p0_;
  int l1_10p5_;
  int l1_10p0_;
  int l1_9p5_;
  int l1_9p0_;
  int l1_8p5_;
  int l1_8p0_;
  int l1_7p5_;
  int l1_7p0_;
  int l1_6p5_;
  int l1_6p0_;
  int l1_5p5_;
  int l1_5p0_;
  int l1_4p5_;
  int l1_4p0_;

  // Decay and acceptance
  int is_bkee_;
  int in_acc_;

  // RECO-to-GEN matching
  int isMatched_;

  // GEN pt,eta
  float e1_gen_pt_;
  float e2_gen_pt_;
  float e1_gen_eta_;
  float e2_gen_eta_;
  float e12_gen_dr_;

  // RECO pt,eta
  float e1_reco_pt_;
  float e2_reco_pt_;
  float e1_reco_eta_;
  float e2_reco_eta_;
  float e12_reco_dr_;

  // RECO algo
  int e1_reco_pf_;
  int e2_reco_pf_;
  int e1_reco_lowpt_;
  int e2_reco_lowpt_;
  int e1_reco_overlap_;
  int e2_reco_overlap_;

  // RECO ID
  int e1_reco_loose_;
  int e1_reco_medium_;
  int e1_reco_tight_;

  int e2_reco_loose_;
  int e2_reco_medium_;
  int e2_reco_tight_;

  // Pre-selection and BDT
  float ip3d_;
  float cos2d_;
  float bdt_;
  float mll_;

  // B candidates
  float b_mass_;
  float b_mass_err_;
  float b_pt_;
  float b_l1_pt_;
  float b_l2_pt_;
  float b_k_pt_;
  float b_cos2D_;
  float b_lxy_;
  float b_lxyerr_;
  float b_svprob_;
  //bdt vars
  float BToKEE_fit_l1_normpt_;
  float BToKEE_fit_l2_normpt_;
  float BToKEE_l1_dxy_sig_;
  float BToKEE_l2_dxy_sig_;
  float BToKEE_fit_k_normpt_;
  float BToKEE_k_DCASig_;
  float BToKEE_k_dxy_sig_;
  float BToKEE_fit_normpt_;
  float BToKEE_l_xy_sig_;
  float BToKEE_eleDR_;
  float BToKEE_llkDR_;
  float BToKEE_l1_iso04_rel_;
  float BToKEE_l2_iso04_rel_;
  float BToKEE_k_iso04_rel_;
  float BToKEE_b_iso04_rel_;
  float BToKEE_ptAsym_;
  float BToKEE_l1_dzTrg_;
  float BToKEE_l2_dzTrg_;
  float BToKEE_k_dzTrg_;
  float BToKEE_l1_pfmvaId_lowPt_;
  float BToKEE_l2_pfmvaId_lowPt_;
  float BToKEE_l1_pfmvaId_highPt_;
  float BToKEE_l2_pfmvaId_highPt_;
  //new vars
  float BToKEE_b_iso03_;
  float BToKEE_b_iso03_dca_;
  float BToKEE_b_iso03_dca_tight_;
  float BToKEE_b_iso04_dca_ ;
  float BToKEE_b_iso04_dca_tight_ ;
  float BToKEE_fit_eta_ ;
  float BToKEE_fit_k_eta_ ;
  float BToKEE_fit_k_phi_ ;
  float BToKEE_fit_phi_ ;
  float BToKEE_k_iso03_ ;
  float BToKEE_k_iso03_dca_ ;
  float BToKEE_k_iso03_dca_tight_ ;
  float BToKEE_k_iso04_dca_ ;
  float BToKEE_k_iso04_dca_tight_ ;
  float BToKEE_k_svip2d_ ;
  float BToKEE_k_svip2d_err_;
  float BToKEE_k_svip3d_err_;
  float BToKEE_l1_iso03_dca_ ;
  float BToKEE_l1_iso03_dca_tight_ ;
  float BToKEE_l1_iso04_ ;
  float BToKEE_l1_iso04_dca_ ;
  float BToKEE_l1_iso04_dca_tight_ ;
  float BToKEE_l2_iso03_;
  float BToKEE_l2_iso03_dca_ ;
  float BToKEE_l2_iso03_dca_tight_ ;
  float BToKEE_l2_iso04_ ;
  float BToKEE_l2_iso04_dca_ ;
  float BToKEE_l2_iso04_dca_tight_;
  float BToKEE_maxDR_ ;
  float BToKEE_minDR_ ;
  float BToKEE_b_n_isotrk_ ;
  float BToKEE_b_n_isotrk_dca_ ;
  float BToKEE_b_n_isotrk_dca_tight_ ;
  float BToKEE_k_n_isotrk_ ;
  float BToKEE_k_n_isotrk_dca_;
  float BToKEE_k_n_isotrk_dca_tight_ ;
  float BToKEE_l1_n_isotrk_ ;
  float BToKEE_l1_n_isotrk_dca_ ;
  float BToKEE_l1_n_isotrk_dca_tight_ ;
  float BToKEE_l2_n_isotrk_ ;
  float BToKEE_l2_n_isotrk_dca_;
  float BToKEE_l2_n_isotrk_dca_tight_ ;
  float Electron_fBrem_l1_ ;
  float Electron_fBrem_l2_ ;
  float Electron_ip3d_l1_ ;
  float Electron_ip3d_l2_ ;
  float Electron_pfRelIso_l1_ ;
  float Electron_pfRelIso_l2_ ;
  float Electron_sip3d_l1_;
  float Electron_sip3d_l2_ ;
  float Electron_trkRelIso_l1_ ;
  float Electron_trkRelIso_l2_ ;
  float ProbeTracks_dzS_;
  float ProbeTracks_eta_; 
  float ProbeTracks_nValidHits_;
};

#endif
