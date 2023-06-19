/********
*Author: Baran Bodur
*Date: 2022-03-24
*Description: Combined unbinned fit for multiple SK phases (or experiments in general)
* 
*
********/
#ifndef FITTOOLMULTIEXP_INCLUDED
#define FITTOOLMULTIEXP_INCLUDED

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
#include "fitTool.h"
#include "mockData.h"

class FitPdf;
class MockData;
class FitTool;
class FitToolMultiExp;


//Minimized function and minimizer defined globally for access
ROOT::Math::GradFunctor fitFuncMultiExp;
ROOT::Math::Minimizer* minFarshawMultiExp;

class FitToolMultiExp{
public:
	// Constructor
	FitToolMultiExp(TString inName, const unsigned int nPdf, const unsigned int nFitIn, 
			FitTool *inFitter[]);

	//Destructor
	~FitToolMultiExp();

	// Public Functions
	
	// set (mock-)data to be fitted and parameters
	void resetFit(); // sets fit function, init vals, steps again
	void setData(const unsigned int nFitIn, MockData *inData[]); // sets data to be fitted
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

	// Likelihood Profiles
	double profileFit(int varIndex, std::vector<double> profVals);
	void saveProfileGraphs(TFile* outFile);


	// Projections histograms of a single fit
	void getStackedHists(); // make stacked hists with fit/sim/pois pars
	void saveStackedHists(TFile* outFile); // save stacked hists
	void delStackedHists(); // delete prev stacked hists (to make way for newer ones)
	void printFitRes(TString outPdfName); // Prints per fit plots to pdf
	
	
	// Likelihoods and gradients for fitting
	double getLikelihood(const double *inPars);
	double getLhGradient(const double *inPars, unsigned int coord);
	
	// Public Variables
	TString name;
	unsigned int nType;
	unsigned int nFitter;

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
	
	// Stacked projected histograms per fit with fit,sim and pois params
	std::vector< std::vector<double> > fitParsVec;
	std::vector< std::vector<double> > simTrueParsVec;
	std::vector< std::vector<double> > simPoisParsVec;
	std::vector<double> fitMinLhVec, simTrueLhVec, simPoisLhVec;
	
	unsigned int stackHistCtr;

	// Profiled 2DeltaNLL curves
	std::vector<TGraph*> profileNllPlots;
	unsigned int profileNllCtr;

private:
	// Private Variables
	std::vector<FitTool*> fitters;
	
	// Flags for internal checks
	bool dataWasSet = false;
	bool histsReady = false;
	
	// For functors definitions within class. some member function pointers
	typedef  double (FitToolMultiExp::*LhFn)(const double *inPars);
	typedef  double (FitToolMultiExp::*GradFn)(const double *inPars, unsigned int coord);
	LhFn getLhMultiExp;
	GradFn getGradMultiExp;




};
#endif
