// Based on ROOT/tutorials/roostats/rs301_splot.C

// To be run on data tnp formatted ntuples 
// - Fit B mass as a reference distribution
// - Extract analysis BDT distributions for signal and background with the sPlots technique

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "RooStats/SPlot.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include <iostream>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*, float, float, int, int);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*, int, int);
void MakeHistos(RooWorkspace*, int, int);
void getDataSet(const char *, RooWorkspace*, float, float, int, int);

void B_splot(int wantPFPF, int checkO)
{
  // set range of observable
  Float_t lowRange  = 0.;   
  Float_t highRange = 10.;   
  if (wantPFPF) {
    lowRange  = 4.5;   
    highRange = 6.0;   
  } else {
    lowRange  = 4.8;   
    highRange = 5.8;   
  }
    
  // Create a new workspace to manage the project.
  RooWorkspace* wspace = new RooWorkspace("myWS");
  
  // add the signal and background models to the workspace.
  // Inside this function you will find a description our model.
  AddModel(wspace, lowRange, highRange, wantPFPF, checkO);
  
  // add dataset from converted root tree
  getDataSet("/eos/cms/store/user/crovelli/LowPtEle/TnpDataB/March21noRegression/FormattedTnPForB_March21_ParkingBPH1and2_Run2018D_part1and2.root", wspace, lowRange, highRange, wantPFPF, checkO);
  
  // inspect the workspace if you wish
  wspace->Print();
  
  // make a new dataset with sWeights added for every event.
  DoSPlot(wspace);

  // Make some plots showing the discriminating variable and
  // the control variable after unfolding.
  MakePlots(wspace, wantPFPF, checkO);
  
  // Save variables in histos
  MakeHistos(wspace, wantPFPF, checkO);

  // cleanup
  delete wspace;
}

// Signal and background fit models
void AddModel(RooWorkspace* ws, float lowRange, float highRange, int wantPFPF, int checkO){
  
  // make a RooRealVar for the observables
  RooRealVar B_mass("B_mass", "M_{B}", lowRange, highRange,"GeV");

  // --------------------------------------
  // signal model
  std::cout << "make B model" << std::endl;

  RooRealVar m1("m1", "B Mass", 5.278, 5.27, 5.28);                  
  RooRealVar alpha1("alpha1", "alpha1",  1.85, 0.0, 10.);          

  float sigma1UpperLimit = 0.08;
  float sigma1LowerLimit = 0.02;
  float sigma1Init = 0.05;
  if (wantPFPF==0 && checkO==0) {
    sigma1UpperLimit = 0.04;
    sigma1Init = 0.02;
  }
  RooRealVar sigma1("sigma1", "sigma1",  sigma1Init, sigma1LowerLimit, sigma1UpperLimit);
  
  float n1UpperLimit = 10.0;
  float n1LowerLimit = 0.0;
  float n1Init = 5.0;
  if (wantPFPF==0) {
    n1Init = 3.5;
    n1UpperLimit = 4.0;
    n1LowerLimit = 3.0;
  }
  RooRealVar n1("n1", "N1", n1Init, n1LowerLimit, n1UpperLimit);

  RooCBShape bModel("bModel","bModel",B_mass,m1,sigma1,alpha1,n1);

  RooRealVar bYield("bYield","fitted yield for Bs", 2000 , 1., 500000) ;              

  
  // --------------------------------------
  // background model
  std::cout << "make background model" << std::endl;

  float alphaUpperLimit = 0.;
  if (wantPFPF==0 && checkO==0) alphaUpperLimit = 0.5;
  RooRealVar alpha("alpha", "Decay const for background mass spectrum", -0.1, -1., alphaUpperLimit,"1/GeV");      
  
  RooExponential bkgModel("bkgModel", "bkg Mass Model", B_mass, alpha);

  RooRealVar bkgYield("bkgYield","fitted yield for background", 10000 , 1., 5000000) ;             
  
  // --------------------------------------
  // combined model
  std::cout << "make full model" << std::endl; 
  RooAddPdf model("model","jpsi+background models",
		  RooArgList(bModel, bkgModel),
		  RooArgList(bYield,bkgYield)); 

  std::cout << "import model" << std::endl;
  ws->import(model);
}

// Add s-weights to the dataset
void DoSPlot(RooWorkspace* ws){
  
  std::cout << "Calculate sWeights" << std::endl;
  
  // get what we need out of the workspace to do the fit
  RooAbsPdf* model = ws->pdf("model");
  RooRealVar* bYield = ws->var("bYield");
  RooRealVar* bkgYield = ws->var("bkgYield");
  RooDataSet* data = (RooDataSet*) ws->data("data");

  // fit the model to the data.
  model->fitTo(*data, Extended() );


  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  // This *could* be done with the lines below, however this is taken care of
  // by the RooStats::SPlot constructor (or more precisely the AddSWeight method).
  RooRealVar* m1  = ws->var("m1");
  RooRealVar* sigma1 = ws->var("sigma1");
  RooRealVar* alpha1 = ws->var("alpha1");
  RooRealVar* n1 = ws->var("n1");
  RooRealVar* alpha  = ws->var("alpha");
  m1->setConstant();   
  sigma1->setConstant();   
  alpha1->setConstant();   
  n1->setConstant();   
  alpha->setConstant();   

  RooMsgService::instance().setSilentMode(false);

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
					       *data, model, RooArgList(*bYield,*bkgYield) );


  // Check that our weights have the desired properties
  std::cout << "Check SWeights:" << std::endl;

  std::cout << std::endl <<  "Yield of Bs is "
	    << bYield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("bYield") << std::endl;
  
  std::cout << "Yield of background is "
	    << bkgYield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("bkgYield") << std::endl
	    << std::endl;
  
  for(Int_t i=0; i < 10; i++) {
    std::cout << "B Weight = " << sData->GetSWeight(i,"bYield")
	      << ", bkg Weight = " << sData->GetSWeight(i,"bkgYield")
	      << ", Total Weight = " << sData->GetSumOfEventSWeight(i)
	      << std::endl;
  }
  
  std::cout << std::endl;

  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  ws->import(*data, Rename("dataWithSWeights"));
}

// Control plots
void MakePlots(RooWorkspace* ws, int wantPFPF, int checkO){
  
  // Here we make plots of the discriminating variable (B_mass) after the fit
  // and of the control variable (ID, or whatever) after unfolding with sPlot.
  std::cout << std::endl;
  std::cout << "make plots" << std::endl;

  // make our canvas
  TCanvas* cdataO = new TCanvas("sPlotO","sPlot demo", 400, 600);
  cdataO->Divide(1,3);
  TCanvas* cdataG = new TCanvas("sPlotG","sPlot demo", 400, 600);
  cdataG->Divide(1,3);

  // get what we need out of the workspace
  RooAbsPdf* model    = ws->pdf("model");
  RooAbsPdf* bModel   = ws->pdf("bModel");
  RooAbsPdf* bkgModel = ws->pdf("bkgModel");
  RooRealVar* theAnalysisBdtO = ws->var("theAnalysisBdtO");         
  RooRealVar* theAnalysisBdtG = ws->var("theAnalysisBdtG");         
  RooRealVar* B_mass = ws->var("B_mass");


  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
  
  // this shouldn't be necessary, need to fix something with workspace
  // do this to set parameters back to their fitted values.
  model->fitTo(*data, Extended() );

  // Plot B_mass for data with full model and individual components overlaid
  RooPlot* frame = B_mass->frame() ;
  data->plotOn(frame ) ;
  model->plotOn(frame) ;
  model->plotOn(frame,Components(*bModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame->SetTitle("Fit of model to discriminating variable");

  cdataO->cd(1);
  frame->Draw() ;
  cdataG->cd(1);
  frame->Draw() ;

  // ------------------------------------------------------------
  // Now use the sWeights to show our variable distribution for B and background.
  //
  // Plot our variable for the B component.
  // Do this by plotting all events weighted by the sWeight for the B component.
  // The SPlot class adds a new variable that has the name of the corresponding
  // yield + "_sw".
   
  // create weighted data set
  RooDataSet * dataw_B = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bYield_sw") ;
  
  RooPlot* frame2o = theAnalysisBdtO->frame() ;      
  dataw_B->plotOn(frame2o, DataError(RooAbsData::SumW2) ) ;
  frame2o->SetTitle("Analysis BDT for signal, Otto");
  cdataO->cd(2);
  frame2o->Draw() ;

  RooPlot* frame2g = theAnalysisBdtG->frame() ;      
  dataw_B->plotOn(frame2g, DataError(RooAbsData::SumW2) ) ;
  frame2g->SetTitle("Analysis BDT for signal, George");
  cdataG->cd(2);
  frame2g->Draw() ;

  
  // Plot interesting variables for background
  RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;

  RooPlot* frame3o = theAnalysisBdtO->frame() ;          
  dataw_bkg->plotOn(frame3o,DataError(RooAbsData::SumW2) ) ;
  frame3o->SetTitle("Analysis BDT for background, Otto");
  cdataO->cd(3);
  frame3o->Draw() ;

  RooPlot* frame3g = theAnalysisBdtG->frame() ;          
  dataw_bkg->plotOn(frame3g,DataError(RooAbsData::SumW2) ) ;
  frame3g->SetTitle("Analysis BDT for background, George");
  cdataG->cd(3);
  frame3g->Draw() ;
  
  if (wantPFPF==1 || checkO==1) cdataO->SaveAs("SPlotOtto.png");
  if (wantPFPF==1 || checkO==0) cdataG->SaveAs("SPlotGeorge.png");


  // Fit variable
  TCanvas* cdata2 = new TCanvas("cdata2","data fit", 1);
  RooPlot* frame4 = B_mass->frame() ;
  data->plotOn(frame4, Binning(100) ) ;
  model->plotOn(frame4) ;
  model->plotOn(frame4,Components(*bModel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame4,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame4->SetTitle("Fit of model to discriminating variable");
  frame4->Draw() ;
  cdata2->SaveAs("Fit.png");
}

// Control plots
void MakeHistos(RooWorkspace* ws, int wantPFPF, int checkO){
  
  gStyle->SetOptStat(0);

  std::cout << std::endl;
  std::cout << "save histos" << std::endl;

  RooRealVar* theAnalysisBdtO = ws->var("theAnalysisBdtO");                
  RooRealVar* theAnalysisBdtG = ws->var("theAnalysisBdtG");                

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

  // create weighted data set
  RooDataSet * dataw_b   = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bYield_sw") ;
  RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
  
  // convert to TH1
  TH1 *h1_theAnalysisBdtO_B = dataw_b->createHistogram("h1_theAnalysisBdtO_B",*theAnalysisBdtO,Binning(30));         
  TH1 *h1_theAnalysisBdtO_bkg  = dataw_bkg->createHistogram("h1_theAnalysisBdtO_bkg",*theAnalysisBdtO,Binning(30));           
  TH1 *h1_theAnalysisBdtG_B = dataw_b->createHistogram("h1_theAnalysisBdtG_B",*theAnalysisBdtG,Binning(30));         
  TH1 *h1_theAnalysisBdtG_bkg  = dataw_bkg->createHistogram("h1_theAnalysisBdtG_bkg",*theAnalysisBdtG,Binning(30));           

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();

  h1_theAnalysisBdtO_B->Write();        
  h1_theAnalysisBdtO_bkg->Write();        
  h1_theAnalysisBdtG_B->Write();        
  h1_theAnalysisBdtG_bkg->Write();        

  if (wantPFPF==1 || checkO==1) {
    TCanvas* ch1 = new TCanvas("ch1","ch1", 1);
    h1_theAnalysisBdtO_B->SetLineWidth(2);
    h1_theAnalysisBdtO_bkg->SetLineWidth(2);
    h1_theAnalysisBdtO_B->SetLineColor(6);
    h1_theAnalysisBdtO_bkg->SetLineColor(4);
    h1_theAnalysisBdtO_B->SetTitle("");
    h1_theAnalysisBdtO_bkg->SetTitle("");
    h1_theAnalysisBdtO_B->DrawNormalized("hist");
    h1_theAnalysisBdtO_bkg->DrawNormalized("samehist");
    ch1->SaveAs("BdtOH.png");
  }

  if (wantPFPF==1 || checkO==0){
    TCanvas* ch2 = new TCanvas("ch2","ch2", 1);
    h1_theAnalysisBdtG_B->SetLineWidth(2);
    h1_theAnalysisBdtG_bkg->SetLineWidth(2);
    h1_theAnalysisBdtG_B->SetLineColor(6);
    h1_theAnalysisBdtG_bkg->SetLineColor(4);
    h1_theAnalysisBdtG_B->SetTitle("");
    h1_theAnalysisBdtG_bkg->SetTitle("");
    h1_theAnalysisBdtG_B->DrawNormalized("hist");
    h1_theAnalysisBdtG_bkg->DrawNormalized("samehist");
    ch2->SaveAs("BdtGH.png");
  }
}

// Convert ROOT tree in RooDataset
void getDataSet(const char *rootfile, RooWorkspace *ws, float lowRange, float highRange, int wantPFPF, int checkO) {    
  
  cout << "roofitting file " << rootfile << endl;
  
  // fit variables
  RooRealVar B_mass("B_mass", "M_{inv}", lowRange, highRange,"GeV");
  // 
  // discriminating variables
  RooRealVar pair_mass("pair_mass", "pair_mass", 2.9, 3.2, "");           
  RooRealVar probeMvaId("probeMvaId", "probeMvaId", -0.5, 20.5, "");           
  RooRealVar tagMvaId("tagMvaId", "tagMvaId", -0.5, 20.5, "");           
  RooRealVar probePfmvaId("probePfmvaId", "probePfmvaId", -0.5, 20.5, "");           
  RooRealVar tagPfmvaId("tagPfmvaId", "tagPfmvaId", -0.5, 20.5, "");           
  //
  // BDT output
  RooRealVar theAnalysisBdtO("theAnalysisBdtO", "theAnalysisBdtO", -15., 15., "");           
  RooRealVar theAnalysisBdtG("theAnalysisBdtG", "theAnalysisBdtG", -15., 15., "");           

  RooArgSet setall(B_mass,pair_mass,probeMvaId,tagMvaId,probePfmvaId,tagPfmvaId,theAnalysisBdtO,theAnalysisBdtG);

  TFile *file = TFile::Open(rootfile);
  TTree *tree = (TTree*)file->Get("tnpAna/fitter_tree");

  RooDataSet *data = new RooDataSet("data","data",tree,setall,0); 

  // Inclusive
  if (wantPFPF==1) {   // PF-PF
    data = (RooDataSet*)data->reduce("pair_mass>3 && pair_mass<3.2 && probePfmvaId<20 && tagPfmvaId<20");
  } else {             // PF-LPT
    if (checkO==1) 
      data = (RooDataSet*)data->reduce("pair_mass>3.05 && pair_mass<3.15 && (probeMvaId<20 || tagMvaId<20) && theAnalysisBdtO>0");
    // data = (RooDataSet*)data->reduce("pair_mass>3.05 && pair_mass<3.15 && (probeMvaId<20 || tagMvaId<20) && ( (probeMvaId<20 && probeMvaId>1) || (tagMvaId<20 && tagMvaId>1) ) && theAnalysisBdtO>0");
    else
      //data = (RooDataSet*)data->reduce("pair_mass>3.05 && pair_mass<3.15 && (probeMvaId<20 || tagMvaId<20) && ( (probeMvaId<20 && probeMvaId>1) || (tagMvaId<20 && tagMvaId>1) ) && theAnalysisBdtG>-6");
      data = (RooDataSet*)data->reduce("pair_mass>3.05 && pair_mass<3.15 && (probeMvaId<20 || tagMvaId<20) && theAnalysisBdtG>-4.5");
  }
  data->Print();

  ws->import(*data);
}
