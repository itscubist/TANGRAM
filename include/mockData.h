/********
*Author: Baran Bodur
*Date: 2022-07-26
*Description: Holds and generates mock data given interaction PDFs to do fitting studies
*
********/

#ifndef MOCKDATA_INCLUDED
#define MOCKDATA_INCLUDED

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
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom3.h"
// Project Headers
#include "organizer.h"
#include "experiment.h"
#include "fitKdePdf.h"
#include "utils.h"

class Organizer;
class FitKdePdf;
class Experiment;

struct DataPoint{
	double eVis;
	double gTag;
	int nTag;
	double weight;
	double cosZen;
	double azi;
};

void printDataPoint(DataPoint dP);

class MockData{
public:
	// Constructor
	MockData(Organizer *orgIn, unsigned int expIndexIn, unsigned int mockIndexIn);
	
	//Destructor

	// Public Functions
	unsigned int genMockData(std::vector<double> usedSysVecIn);
	void makeProjHists();
	void saveProjHists();
	std::vector<DataPoint> MockData::makeAsimovDataset(unsigned int intIndIn, unsigned int nEvents);
	std::vector<double> MockData::shuffleOrder(std::vector<double> inVec);
	unsigned int readData(); // read data from outside

	// For reweighting
	void reweightMockData(std::vector<double> oriSysVec, std::vector<double> newSysVec);
	void reweightOneSys(unsigned int sysIndex, double nSigma);
	
	// Public Variables
	Organizer *org;
	unsigned int expIndex;
	unsigned int mockIndex;
	TString mockName;
	
	unsigned int nPoints;
	double totalWeight;
	std::vector<DataPoint> dataPoint;
	std::vector<unsigned int> byTypeEvents;
	std::vector<double> simTrueEvents;
	std::vector<double> simPoisEvents; 
	std::vector<double> simExpEvents; 
	std::vector<double> simSystematics; 

	// Data projection hists:
	std::vector<TH1D*> evisProjs;
	std::vector< std::vector<TH1D*> > gtagProjs;
	TH1D* hNeutProj;
	TH2D* hSysVsChiProb;
	TH1D* hIdealEvis;
	

private: 
	// Private Functions
	
	// Private Variables
	
};
#endif
