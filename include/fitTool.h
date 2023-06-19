/********
*Author: Baran Bodur
*Date: 2021-07-13
*Description: Make unbinned fit given dataset and pdfs
*
********/
#ifndef FITTOOL_INCLUDED
#define FITTOOL_INCLUDED

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
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// Project Headers
#include "fitPdf.h"
#include "mockData.h"

class FitPdf;
class MockData;
class FitTool;


//Minimized function and minimizer defined globally for access
ROOT::Math::GradFunctor fitFunc;
ROOT::Math::Minimizer* minFarshaw;

class FitTool{
public:
	// Constructor
	FitTool(TString inName, const unsigned int nPdf, FitPdf *inPdf[]);

	//Destructor
	~FitTool();

	// Public Functions
	
	// set (mock-)data to be fitted and parameters
	void resetFit(); // sets fit function, init vals, steps again
	void setData(MockData *inData); // sets data to be fitted
	void setInitPars(double *inPars);
	void setFitStep(double *inSteps);
	void setSimPars(); // from input pdfs and mock data
	
	// controls fitting
	double makeFit();

	// To look at many fit results
	void printLastFit();
	void prepFitHists();
	void fillFitHists();
	void saveFitHists(TFile* outFile);

	// Likelihood Projections
	void makeLhProfilePlots();
	void saveLhProfilePlots();

	// Projections histograms of a single fit
	THStack* makeStackedHist(const double *inPars, unsigned int hNo, TString nameStart); // stack his	
	void getStackedHists(); // make stacked hists with fit/sim/pois pars
	void saveStackedHists(TFile* outFile); // save stacked hists
	void delStackedHists(); // delete prev stacked hists (to make way for newer ones)
	void printFitRes(TString outPdfName); // Prints per fit plots to pdf
	
	
	// Likelihoods and gradients for fitting
	double getLikelihood(const double *inPars);
	double getLhGradient(const double *inPars, unsigned int coord);
	
	// Public Variables
	TString name;
	MockData* dataToFit;
	unsigned int nType;
	std::vector<double> fitSteps;
	std::vector<double> curPars;
	std::vector<double> initPars; // fit start point
	std::vector<double> fitPars; // fit end result
	std::vector<double> simTruePars; // simulation truth
	std::vector<double> simPoisPars; // simulated with poisson fluc
	float fitMinLh, simTrueLh, simPoisLh;

	// Histograms to judge the average fit performance
	std::vector<TH1F*> h1FitDist;
	std::vector< std::vector<TH2F*> > h2FitDist;
	TH1F* hFitLh, *hSimPoisLh, *hSimTrueLh, *hLhFitMinPois, *hLhFitMinTrue;

	// Fit Likelihood Graphs
	std::vector< std::vector<TGraph*> > fitLhProfiles;
	
	unsigned int stackHistCtr;
	// Stacked projected histograms per fit with fit,sim and pois params
	std::vector< std::vector<THStack*> > hFitProj;
	std::vector< std::vector<THStack*> > hSimProj;
	std::vector< std::vector<THStack*> > hPoisProj;
	std::vector< std::vector<TH1D*> > hDataProj;
	std::vector< std::vector<double> > fitParsVec;
	std::vector< std::vector<double> > simTrueParsVec;
	std::vector< std::vector<double> > simPoisParsVec;
	std::vector<double> fitMinLhVec, simTrueLhVec, simPoisLhVec;
	
	// Stacked projected histograms per fit
	//std::vector<THStack*> hfEvisProjN0, hfEvisProjN1, hfEvisProjNall;
	//std::vector<THStack*> hfGtagProjN0, hfGtagProjN1, hfGtagProjNall;
	//std::vector<THStack*> hfNtagProj;
		
	// Stacked projected histograms per fit
	//std::vector<THStack*> hsEvisProjN0, hsEvisProjN1, hsEvisProjNall;
	//std::vector<THStack*> hsGtagProjN0, hsGtagProjN1, hsGtagProjNall;
	//std::vector<THStack*> hsNtagProj;

	// Stacked projected histograms per fit
	//std::vector<THStack*> hpEvisProjN0, hpEvisProjN1, hpEvisProjNall;
	//std::vector<THStack*> hpGtagProjN0, hpGtagProjN1, hpGtagProjNall;
	//std::vector<THStack*> hpNtagProj;

	// To store PDF pointers publicly
	std::vector<FitPdf*> pdfs;

private:
	// Private Variables
	TRandom3* randAlThor; //need rng
	
	// Flags for internal checkcs
	bool dataWasSet = false;
	bool histsReady = false;
	
	// For functors definitions within class. some member function pointers
	typedef  double (FitTool::*LhFn)(const double *inPars);
	typedef  double (FitTool::*GradFn)(const double *inPars, unsigned int coord);
	LhFn getLh;
	GradFn getGrad;




};
#endif
