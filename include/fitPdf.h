/********
*Author: Baran Bodur
*Date: 2021-07-08
*Description: Make PDFs from input tree, to be used in fitting later
*
********/
#ifndef FITPDF_INCLUDED
#define FITPDF_INCLUDED

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
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom3.h"
// Project Headers
#include "mockData.h"

struct DataPoint;
class MockData;
class FitTool;

class FitPdf{
public:
	// Constructor from input file and histograms in it (estimated via kde)
	FitPdf(TString inName, TFile* inFile, TString intName, int wallOpt, 
			int histOpt, int inNeutType, int inGammaType, int inDirType, int inCol=2);
	// Consturctor from input tree and appropriate cut
	FitPdf(TString inName, TTree* inTree, TString inCut, int inNeutType, 
			int inGammaType, int inCol=2); 
	
	//Destrcutor
	~FitPdf();

	// Public Functions
	float getProb(float eVis, float gTag, int nTag); // get value of pdf at point (interpolates)
	float getProb(float eVis, float gTag, int nTag, float cosZen, float azi); // get value of pdf at point (interpolates)
	float getProb(DataPoint inPoint); // get value of pdf at point (interpolates)
	DataPoint getRandPoint(); // get a random data point from the pdf 
	int getRandNtag(); // get a random ntag (not the whole data point)
	void getRandEneAtNtag(float &eVis, float &gTag, int nTag); // get random eVis/gTag at given nTag	
	void getRandDirAtEne(float &cosZen, float &azi, float eVis);
	void makeSmoothDists(TString smoothOpt, unsigned int nTimes); // Provide smoothing for pdfs
	void savePlots(TFile* inFile); // save histograms
	float scaleExpEvents(float scale); // scale expected events
	float scaleSimEvents(float scale); // scale expected events
	
	// To get projections
	TH1D* getEneProj(int nTag, float scale=1.0); // get eVis projection of histograms
	TH1D* getGtagProj(int nTag, float scale=1.0); // get gTag projection of histograms
	TH1D* getNtagProj(float scale=1.0); // get nTag projection of histograms
	TH1D* getCosZenProj(float scale=1.0); // get cosZen projection of histograms
	TH1D* getAziProj(float scale=1.0); // get azi projection of histograms
	void makeProjHists(float scale=1.0); // Stores projection histograms of this pdf


	// Public Variables
	float expEvents;
	float simEvents;
	TString pdfName;
	int pdfColor;
	int neutInfType; // :0 no neutron info (pdf is 2d), 1: neutron info
	int gammaInfType; // :0 no gamma tagging used, 1: gamma tagging used
	int dirInfType; // :0 no dir info, 1: dir info is used

	// Public projection hists
	TH1D *hEvisProjN0, *hEvisProjN1, *hEvisProjNall;
	TH1D *hGtagProjN0, *hGtagProjN1, *hGtagProjNall;
	TH1D *hNtagProj;
	TH1D *hCosZenProj, *hAziProj;
	std::vector<TH1D*> projHists;
	
	// Smoothed histogram range
	int hsEvisBins;
	double hsEvisLow, hsEvisUp;
	int hsGtagBins;
	double hsGtagLow, hsGtagUp;
	
	// histograms to store smoothed info
	TH2F* hsEgN0;
	TH2F* hsEgN1;
	TH2F* hsEgNall;

	// dir histograms
	TH2F* hDir2d;
	TH2F* hDir2dEvisR[4];

private:
	// Private Variables
	TRandom3* randAlThor; //need rng
	float neutProb[2]; // prob of 0 and 1 neutron respectively
	TH1F* hNeutProb;

	// raw histogram range
	int hrEvisBins;
	double hrEvisLow, hrEvisUp;
	int hrGtagBins;
	double hrGtagLow, hrGtagUp;

	// histograms to store raw info
	TH2F* hrEgN0;
	TH2F* hrEgN1;
	TH2F* hrEgNall;
	
	// temporary histogram
	TH2F* hTemp;

};
#endif
