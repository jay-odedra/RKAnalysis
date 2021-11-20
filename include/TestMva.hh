#ifndef TestMva_h
#define TestMva_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class TestMva : public BParkBase{
public:

  //! constructor
  TestMva(TTree *tree=0);
  //! destructor
  virtual ~TestMva();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  
private:
  
  // Analysis methods
  bool isMcB( int myB );
  void bookOutputTree();

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;

  // dataset name
  std::string _datasetName;      
  
  // settings
  int sampleID;

  //---output tree branches variables
  int    theRun;
  int    theEvent;
  int    theSampleID;       
  int    nvtx;
  //
  vector <float> Ele1Pt={};   
  vector <float> Ele1Eta={};
  vector <float> Ele1Phi={};
  vector <bool>  Ele1IsPF={};
  vector <bool>  Ele1IsLowPt={};
  vector <float> Ele1MvaId={};
  vector <float> Ele1PfmvaId={};
  vector <float> Ele1Iso03={};
  vector <bool>  Ele1MatchMc={};
  //
  vector <float> Ele2Pt={};   
  vector <float> Ele2Eta={};
  vector <float> Ele2Phi={};
  vector <bool>  Ele2IsPF={};
  vector <bool>  Ele2IsLowPt={};
  vector <float> Ele2MvaId={};
  vector <float> Ele2PfmvaId={};
  vector <float> Ele2Iso03={};
  vector <bool>  Ele2MatchMc={};
  //
  vector <float> Bmass={};
  vector <float> Bpt={};
  vector <float> Bcos2D={};
  vector <float> Bsvprob={};
  vector <float> Bxysig={};
  vector <bool> BmatchMC={};
  //
  vector <float> KPt={};
  vector <float> KEta={};
  vector <float> KPhi={};
  //
  vector <float> dzL1K={};
  vector <float> dzL2K={};
  vector <float> dzLK={};
  vector <float> drL1K={};
  vector <float> drL2K={};
  vector <float> drLK={};
  vector <float> drL1L2={};

  vector <float> analysisBdt={};
};

#endif
