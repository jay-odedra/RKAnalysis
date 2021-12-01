#ifndef Efficiency_h
#define Efficiency_h

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
  
  Efficiency(TChain* tree, int isMC, std::string filename);
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
  bool isBToKEE(Daughters& daughters);
  bool inAcceptance(std::vector<int>& indices);
  bool setElePtEta(Daughters& daughters,
		   int& e1_gen_idx,
		   float& e1_gen_pt,
		   float& e1_gen_eta,
		   int& e2_gen_idx,
		   float& e2_gen_pt,
		   float& e2_gen_eta);
  bool recoCand(int theB,
		int& e1_gen_idx,int& e2_gen_idx,
		int& e1_reco_idx,float& e1_reco_pt,float& e1_reco_eta,
		int& e1_reco_pf,int& e1_reco_lowpt,int& e1_reco_overlap,
		int& e2_reco_idx,float& e2_reco_pt,float& e2_reco_eta,
		int& e2_reco_pf,int& e2_reco_lowpt,int& e2_reco_overlap);
    
private:
  
  // Misc
  int verbose_=0;
  std::string filename_="";

  // Output
  TFile* outFile_=nullptr;
  TTree* outTree_=nullptr;
  TH1F* h_entries_=nullptr;
  TH1F* h_selection_=nullptr;
  TH1F* h_cand_=nullptr;

  // Scalars
  int theRun_;
  int theEvent_;
  int nvtx_;

  // Trigger
  float trg_muon_pt_;
  int hlt7_ip4_;
  int hlt8_ip3_;
  int hlt9_ip6_;
  int hlt12_ip6_;

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

  // RECO pt,eta
  float e1_reco_pt_;
  float e2_reco_pt_;
  float e1_reco_eta_;
  float e2_reco_eta_;

  // RECO algo
  int e1_reco_pf_;
  int e2_reco_pf_;
  int e1_reco_lowpt_;
  int e2_reco_lowpt_;
  int e1_reco_overlap_;
  int e2_reco_overlap_;

  // Pre-selection and BDT
  float ip3d_;
  float cos2d_;
  float bdt_;

};

#endif
