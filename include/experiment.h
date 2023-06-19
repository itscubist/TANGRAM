/********
*Author: Baran Bodur
*Date:  2022-07-26
*Description: Holds info about experiment, its interaction PDFs and mock (or real) data
*
********/

#ifndef EXPERIMENT_INCLUDED
#define EXPERIMENT_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
// ROOT headers
#include "TMath.h"
#include "TText.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "THStack.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
// Project Headers
#include "fitKdePdf.h"
#include "organizer.h"
#include "mockData.h"

class MockData;
class Organizer;

class Experiment{
public:
	// Constructor
	Experiment(Organizer *orgIn, unsigned int expIndex);

	//Destructor
	
	// Public Functions
	void genMockData(std::vector< std::vector<double> > usedSysVecIn);
	void reweightOneSys(unsigned int sysIndex, double nSigma, unsigned int mockIndexIn);
	void updateAllPdfs();
	void saveCurPdfs(TString saveName);

	void updateProjHists(double scale);
	void saveCurProj(TString saveName);

	void printSysHists(double sysShift=1.0);
	void checkSysEffect(unsigned int sysIndex, unsigned int pdfIndex);

	THStack* makeStackedHist(std::vector<double> intFracs, std::vector<double> sysVals, 
		TString type, unsigned int nX, unsigned int eX, TString nameStart);
	TH1D* makeTotalHist(std::vector<double> intFracs, std::vector<double> sysVals, 
		TString type, unsigned int nX, unsigned int eX, TString nameStart);
	
	// Public Variables
	unsigned int expIndex;
	std::vector<FitKdePdf*> pdfs;
	std::vector<MockData*> mocks;

private:
	// Private Functions
	Organizer* org;
	
	// Private Variables


};
#endif
