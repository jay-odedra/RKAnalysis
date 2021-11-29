#include "../include/Efficiency.hh"
#include "TChain.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc,
	 char **argv) {

  if (argc <2 ) {
    std::cerr << "Usage: "
	      << argv[0]
	      << " [file_list.txt] [outputFile] [1=MC,0=data] [nmin] [nmax] [ptB] [svprob] [co2D]"
	      << std::endl;
    return 1;
  }
  
  char inputFileName[500];
  strcpy(inputFileName,argv[1]);

  char outputFileName[500];
  if (argc <= 2 ) strcpy(outputFileName,argv[1]);
  else strcpy(outputFileName,argv[2]);

  int isMC=1;
  isMC=atoi(argv[3]);

  int nmin=0;
  int nmax=10000;
  if(argc>4){
    nmin=atoi(argv[4]);
    nmax=atoi(argv[5]);
  }

  float ptB=0.;
  float svprob=0.;
  float cos2D=0.;
  if(argc>6){
    ptB=atof(argv[6]);
    svprob=atof(argv[7]);
    cos2D=atof(argv[8]);
  }

  std::cout << " your request: isMC=" << isMC
	    << " nmin=" << nmin
	    << " nmax=" << nmax
	    << " cuts ptb=" << ptB
	    << " svprob=" << svprob
	    << " cos2D=" << cos2D
	    << std::endl;
  
  // -------------------------
  // Loading the file
  TChain* theChain = new TChain("Events");

  char Buffer[5000];
  char MyRootFile[10000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  char tmpFileName[256];
  std::vector<std::string> filesToRemove;
  int nfiles=0;
  int i=0;
  while( !(inputFile->eof()) && i<10000 ){
    inputFile->getline(Buffer,500);
    if(nfiles>=nmin&& nfiles<nmax){
      if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
	{
	  nfiles=nfiles+1;
	  sscanf(Buffer,"%s",MyRootFile);
	  theChain->Add(TString(MyRootFile));
	  std::cout << "chaining " << MyRootFile << std::endl;
	}
    }
    i=i+1;
  }
  inputFile->close();
  delete inputFile;
  std::cout << "we will process " << nfiles << " files"
	    << std::endl;
  
  Efficiency tnp(theChain, isMC, outputFileName);
  //tnp.SetCuts(ptB,svprob,cos2D);
  tnp.Loop();
  
  return 0;

}
