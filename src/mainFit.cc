/********
 *
 *Author: Baran Bodur
 *Date: 2022-06-28
 *Description: Code for toy MC to study fitting the nue16O based on histograms obtained from
 *	drawer.cc. Uses unbinned liklelihood, where model PDFs are built from the input MC tree
 *	This version is designed to incorporate systematics via weighted KDEs
 *	Arguments: 
 *		1) Input card file name
 *		2) Output filename (without extensitons)
 *
 ********/

// C++ libs
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

// ROOT Libraries
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

// Project Libraries
#include "utils.h"
#include "vectorUtils.h"
#include "organizer.h"
#include "experiment.h"
#include "fitKdePdf.h"
#include "mockData.h"
#include "fitGlobal.h"

// namespaces and typedefs
using namespace std;
using namespace TMath;

// Start Main Function
int main(int argc, char * argv[]) {
	
	cout << "Death is but a feather, duty heavier than mountain..." << endl;
	if(argc < 3) {
		cout << "Not enough arguments. Usage: ./mainFit [cardFile] [outName]" << endl;
	}
	setRootStyle(); // set root style
	TString cardName = (TString) argv[1];
	TString outName = (TString) argv[2];
	Organizer org1("org1",cardName,outName); // create organizer via input card

	// Save 1.0 (or -1.0) sigma systematic changes for expIndex (SK4)
	//expIndex=1000 means all experiments
	//org1.printSysHists(1000,1.0);
	
	//org1.getChi2SysDiff(0, 3.0);
	//org1.getChi2SysDiff(0, 2.0);
	//org1.getChi2SysDiff(0, 1.0);
	//org1.getChi2SysDiff(0, 0.0);
	//org1.getChi2SysDiff(0, -1.0);
  //org1.getChi2SysDiff(0, -2.0);
  //org1.getChi2SysDiff(0, -3.0);
	//org1.checkSysEffect(0,2);
	
	// Generate mock data
	org1.genMockData();

	// Make fits to mock data
	org1.makeAllFits();
	
	// Make only profile fit...
	//org1.onlyProfileFit();
	
	// Make fits to mock data
	//org1.reweightFit();
	
	// Scan all pars with mock data
	//org1.printParScans(outName);


	org1.closeOutputs();	
	return 0;
} // The end
