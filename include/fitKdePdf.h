/********
*Author: Baran Bodur
*Date: 2022-06-16
*Description: Make KDE based PDFs from input tree, to be used in fitting later
*	 Systematics can be taken into account by changing event weights used in preparing the KDE
*
********/
#ifndef FITKDEPDF_INCLUDED
#define FITKDEPDF_INCLUDED

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
#include "TKDE.h"
// Project Headers
#include "organizer.h"
#include "experiment.h"
#include "mockData.h"

class Organizer;
class FitKdePdf;
struct DataPoint;

// Struct to hold many KDE setup information (vectors etc.)
struct KdeHolder {
	// Descriptors
	TString name;
	bool evisOrGtag;
	int neutronType;
	int evisSubbin;
	// Info With 0 systematics
	unsigned int evCount;
	double evTotalWeight0;
	std::vector<double> evX0;
	std::vector<double> evWeight0;
	TKDE* kdePdf0;
	TH1D* histPdf0;
	TH1D* histRaw;
	// Systematic Weights From Tree
	unsigned int nRelScaleSys;
	unsigned int nRelShiftSys;
	std::vector<bool> relScaleSys; // Which scale systematics are relevant for this KDE
	std::vector<bool> relShiftSys; // Which shift systematics are relevant for this KDE
	std::vector< std::vector<double> > evSysWeight;	// for each scale syst affecting this KDE
	std::vector<double> curSys;
	// Info With Current Systematics
	double evTotalWeight;
	std::vector<double> evX;
	std::vector<double> evWeight;
	TKDE* kdePdf;
	TH1D* histPdf;
	TH1D* histProj;
}; 


class FitKdePdf{
public:
	// Constructor from input file and histograms in it (estimated via kde)
	FitKdePdf(Organizer *orgIn, unsigned int expIndexIn, unsigned int intIndexIn);
	
	//Destructor

	// Public Functions for Setup
	std::vector<bool> getRelSys(unsigned int &nRelSys, 
			TString scaleOrShift="any", TString sysType="any"); 
	void fillKdeHolders(); // loops over experiment tree to fill releveant info to kde holders
	void initKdeHolders(); // initialize kde holder with previously filled info
	void updateAllPdfs(TString upOpt=""); // update all kde holders, ntag prob and exp events with sys
	void saveCurPdfs(TString saveName); // save current pdf hists into root output file
	int updateKdeHolder(KdeHolder &khIn, TString upOpt=""); // update given holder with new sys
	int updateNtagAndNev(TString upOpt=""); // update ntag prob and exp event with sys
	void fillWhichKdeInd();
	// For evis subbins
	unsigned int getEvisSubbinForGtag(double evis); // which gtag kde should be called for given evis
	TString getGtagSubbinRange(unsigned int subbin); // get a descriptive name for subbin

	// For projection histograms:
	void updateProjHists(double scale);
	void saveCurProj(TString saveName);
		
	// To get probabilities, random values etc...
	double getProb(double eVis, double gTag, int nTag); // get value of pdf at point (interpolates)
	double getProb(double eVis, double gTag, int nTag, double cosZen, double azi); // with dir 
	double getProb(DataPoint inPoint); // get value of pdf at point (interpolates)
	DataPoint getRandPoint(); // get a random data point from the pdf 
	int getRandNtag(); // get a random ntag (not the whole data point)
	double getRandEneAtNtag(int nTag); // get random eVis/gTag at given nTag	
	double getRandGtagAtEneAndNtag(double eVis, int nTag); // get random eVis/gTag at given nTag	
	void getRandDirAtEne(double &cosZen, double &azi, double eVis);
	void savePlots(TFile* inFile); // save histograms
	
	double scaleExpEvents(double scale); // scale expected events
	double scaleSimEvents(double scale); // scale expected events
	
	// To get projections
	TH1D* getEneProj(int nTag, float scale=1.0); // get eVis projection of histograms
	TH1D* getGtagProj(int nTag, float scale=1.0); // get gTag projection of histograms
	TH1D* getNtagProj(float scale=1.0); // get nTag projection of histograms
	TH1D* getCosZenProj(float scale=1.0); // get cosZen projection of histograms
	TH1D* getAziProj(float scale=1.0); // get azi projection of histograms
	void makeProjHists(float scale=1.0); // Stores projection histograms of this pdf


	// Public Variables for Setup
	Organizer *org;
	TString pdfName;
	unsigned int expIndex;
	unsigned int intIndex;
	
	// ** For Fundamental KDE Vectors
	unsigned int nRelScaleSys;
	unsigned int nRelShiftSys;
	std::vector<bool> relScaleSys;
	std::vector<bool> relShiftSys;
	std::vector<KdeHolder> evisKdes; // vector is over neutron count: Nall, N0, N1 
	std::vector< std::vector<KdeHolder> > gtagKdes; // vector is over neutron count and evis subbin
	std::vector<unsigned int> ntagInds;
	std::vector< std::vector<unsigned int> > evisInds;

	// For neutron tagging and expected events
	double simEvents;
	double expEvents0;
	double expEvents;
	double neutProb0[2]; // prob of 0 and 1 neutron respectively
	double neutProb[2]; // prob of 0 and 1 neutron respectively
	TH1D* hNeutProb0;
	TH1D* hNeutProb;
	TH1D* hNeutProj;
	unsigned int nRelScaleSysNeut;
	unsigned int nRelShiftSysNeut;
	std::vector<bool> relScaleSysNeut;
	std::vector<bool> relShiftSysNeut;
	std::vector<double> curSysNeut;
	unsigned int nRelScaleSysNev;
	unsigned int nRelShiftSysNev;
	std::vector<bool> relScaleSysNev;
	std::vector<bool> relShiftSysNev;
	std::vector<double> curSysNev;
	unsigned int neutEvCount;
	std::vector<double> neutEvX0;
	std::vector<double> neutEvWeight0;
	std::vector<int> ntagVec0;
	std::vector< std::vector<double> > neutSysWeights;
	std::vector< std::vector<double> > nevSysWeights;


	// Public projection hists
	TH1D *hEvisProjN0, *hEvisProjN1, *hEvisProjNall;
	TH1D *hGtagProjN0, *hGtagProjN1, *hGtagProjNall;
	TH1D *hNtagProj;
	TH1D *hCosZenProj, *hAziProj;
	std::vector<TH1D*> projHists;
	
	// dir histograms
	TH2F* hDir2d;
	TH2F* hDir2dEvisR[4];

private:
	// Private Variables

};
#endif
