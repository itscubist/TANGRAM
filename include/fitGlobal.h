/********
*Author: Baran Bodur
*Date: 2022-07-26
*Description: Performs multi experiment / multi PDF fitting in the presence of systematics
*
********/

#ifndef FITGLOBAL_INCLUDED
#define FITGLOBAL_INCLUDED

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
#include "TLegend.h"
#include "TBox.h"
#include "TLine.h"
#include "TText.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// Project Headers
#include "organizer.h"
#include "experiment.h"
#include "fitKdePdf.h"
#include "mockData.h"

struct FitResult{
	std::vector<double> fitFracs;
	std::vector<double> fitFracErrs;
	std::vector<double> fitSys;
	std::vector<double> fitSysErrs;
	std::vector<double> simTrueFracs;
	std::vector<double> simTrueSys;
	std::vector<double> simPoisFracs;
	std::vector<double> simPoisSys;
	double fitLh;
	double fitMainFixedLh;
	double simTrueLh;
	double simPoisLh;
	double fitChi2;
	double simTrueChi2;
	double simPoisChi2;
	// For testing
	double nDataEvents;
	double fitEvisIntegral;
	double simEvisIntegral;
	double fitEvisBeginDiff;
	double simEvisBeginDiff;
};


class FitGlobal{
public:
	// Constructor
	FitGlobal(Organizer *orgIn);

	//Destructor

	// Public Functions
	
	// set (mock-)data to be fitted and parameters
	void resetFit(); // sets fit function, init vals, steps again
	void setMockIndex(unsigned int m); // sets mock data index to be fitted
	void setSimPars(); // from input pdfs and mock data
	
	// controls fitting
	double makeFit(bool fixSys, bool fixMainPar=false);
	double makeProfileFit(unsigned int varIndex, std::vector<double> profVals, bool fixSys);

	// utility
	std::vector<double> combIntAndSys(std::vector<double> intVec, std::vector<double> sysVec);
	void sepIntAndSys(const double *inPars, std::vector<double> &intVec, std::vector<double> &sysVec);

	// To look at many fit results
	void printLastFit();
	void prepFitHists();
	void fillFitHists();
	void saveFitHists();
	void saveProfGraphs();

	void printFitGraphs(unsigned int fitNo, TString type);

	// To Get Binned Chi Square Value and Print Fit Result Into Root File
	void FitGlobal::getFitChiAndHists(unsigned int fitNo, TString type); 
	
	// likelihood scan
	std::vector<double> scanParameter(unsigned int index, std::vector<double> inValues, 
			bool isSim=true);
	void scanAll();
	void scanAfterFit();

	// Likelihoods and gradients for fitting (Wei at the end means, it works for weighted data)
	double getLikelihood(const double *inPars);
	double getLhGradient(const double *inPars, unsigned int coord);
	double getLikelihoodWei(const double *inPars);
	double getLhGradientWei(const double *inPars, unsigned int coord);

	
	// Public Variables
	Organizer *org;
	unsigned int nDim;
	unsigned int mockIndex;
	std::vector<double> fitCurValues;
	std::vector<double> sysCurValues;
	std::vector<double> intCurFracs;
	std::vector<FitResult> fitResults;
	std::vector< std::vector<FitResult> > fitResultsPerExp;
	
	// Histograms to judge the average fit performance
	// 1 entry from each fit
	std::vector<TH1F*> h1FitDist, h1FitDiffDist;
	std::vector< std::vector<TH2F*> > h2FitDist;
	TH1F *hFitLh, *hFitFixedLh, *hSimPoisLh, *hSimTrueLh, *hLhFitMinPois, *hLhFitMinTrue, 
			 *hLhFitMinFitFixed;
	std::vector<TH2F*> h2FitDiffDist;
	std::vector<TH2F*> h2FitVsTrueDist;
	
	TH1F *hFitChi2, *hSimTrueChi2, *hChi2FitMinTrue;
	TH1F *hNDataEvents, *hFitEvisIntegral, *hSimEvisIntegral, *hFitEvisBeginDiff, *hSimEvisBeginDiff;
	std::vector<TH2F*> h2Chi2DiffDist;
	
	// For profiling
	std::vector<TGraph*> profGraphs;
	unsigned int nProfCtr;
	
	// Minimzer and minimized function
	ROOT::Math::GradFunctor fitFunc;
	//ROOT::Math::Functor fitFunc;
	ROOT::Math::Minimizer* minFarshaw;

private:
	// Private Functions
	
	// Private Variables
	bool isMockData;
	bool dataIsSet;
	bool histsReady;
	
	// For functors definitions within class. some member function pointers
	typedef  double (FitGlobal::*LhFn)(const double *inPars);
	typedef  double (FitGlobal::*GradFn)(const double *inPars, unsigned int coord);
	LhFn getLh;
	GradFn getGrad;
	
};
#endif
