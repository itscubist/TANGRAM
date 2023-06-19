/********
*Author: Baran Bodur
*Date: 2022-06-28
*Description: Organizer for unbinned fitter, reads card file, stores all the necessary parameters.
*
********/
#ifndef ORGANIZER_INCLUDED
#define ORGANIZER_INCLUDED

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
#include "TCanvas.h"
#include "TRandom.h"
#include "TRandom3.h"
// Project Headers
#include "experiment.h"
#include "fitKdePdf.h"
#include "mockData.h"
#include "fitGlobal.h"

class MockData;
class FitKdePdf;
class Experiment;
class FitGlobal;
class Organizer;

class Organizer{
public:
	// Constructor
	Organizer(TString inName, TString cardFile, TString outFileName);
	
	//Destructor
	//~Organizer();

	//Public Functions
	void readCardFile(TString cardFile);
	void genMockData();
	void reweightOneSys(unsigned int sysIndex, double nSigma, unsigned int mockIndexIn);
	void reweightFit();
	void updateCurSys(std::vector<double> inSys);
	void updateAllPdfs();
	void saveCurPdfs(TString saveName);
	void onlyProfileFit();
	void makeAllFits();
	void updateProjHists(double scale);
	void saveCurProj(TString saveName);

	void printSysHists(unsigned int expIndex, double sysShift=1.0); 
	void printParScans(TString outFileName);
	void checkSysEffect(unsigned int sysIndex, unsigned int pdfIndex);
	void getChi2SysDiff(unsigned int sysNo, double sysVal);

	void openOutputs(TString outFileName);
	void closeOutputs();


	// Public Variables
	
	// Fixed!!
	// Read from input card to use in the code globally
	
	// Generic
	TString cardFileName;
	TString orgName;
	int genRandData;
	bool nEventsPoisson;
	bool randSystematics;
	bool fixSystematics;
	bool doMainFixedFit;
	bool doReweightFit;
	unsigned int reweightFitSysIndex;
	bool saveBinnedFitResults;
	int randSeed;
	unsigned int nSubFit;
	unsigned int nFit;
	unsigned int nGraph;
	unsigned int nProfile;
	unsigned int nProfSubFit;
	int profIndex;
	double profStart;
	double profEnd;
	double profStep;
	bool makeSysEffectPlots;
	// Experiment Related
	unsigned int nExp;
	std::vector<TString> expNames;
	std::vector<TString> expFiles;
	std::vector<TString> expTrees;
	std::vector<bool> neutronInfo;
	std::vector<bool> gammaInfo;
	std::vector<bool> dirInfo;
	std::vector<double> evisLow;
	std::vector<double> evisUp;
	std::vector<double> evisLowBuf;
	std::vector<double> evisUpBuf;
	std::vector<int> evisBinsNorm;
	std::vector<int> evisProjRebinFactor;
	std::vector<double> gtagLow;
	std::vector<double> gtagUp;
	std::vector<double> gtagLowBuf;
	std::vector<double> gtagUpBuf;
	std::vector<int> gtagBinsNorm;
	std::vector<int> gtagProjRebinFactor;
	std::vector<double> reweightFactors;
	std::vector<TString> dataFiles;
	std::vector<TString> dataTrees;
	// Interaction Related
	unsigned int nInt;
	std::vector<TString> intNames;
	std::vector<int> intColors;
	std::vector<double> intVarFraction; 
	std::vector<double> intInitValues; 
	std::vector<double> intInitSigmas; 
	std::vector<int> intIsFixed; 
	std::vector<int> evisSubbinsForGtag; 
	std::vector<int> intNameCutLow;
	std::vector<int> intNameCutUp;
	std::vector<int> nmueCutLow;
	std::vector<int> nmueCutUp;
	std::vector<double> nuEneCutLow;
	std::vector<double> nuEneCutUp;
	std::vector<double> wallCutLow;
	std::vector<double> wallCutUp;
	std::vector<double> intStepSizeRatio;
	std::vector<double> intStepSize;
	std::vector<TString> evisKdeOptStr;
	std::vector<double> evisKdeRho;
	std::vector<double> evisKdeLow;
	std::vector<TString> gtagKdeOptStr;
	std::vector<double> gtagKdeRho;
	std::vector<double> gtagKdeLow;
	std::vector<int> intScan;
	// Systematic Related
	unsigned int nSys;
	double sysWeightMinScale;
	double sysWeightMaxScale;
	double sysWeightMinShift;
	double sysWeightMaxShift;
	double sysLowLimit;
	double sysUpLimit;
	std::vector<TString> sysNames;
	std::vector<bool> shiftOrScale; 
	std::vector<TString> weightNames;
	std::vector<bool> effNeutron;
	std::vector<bool> effGamma;
	std::vector<bool> effEnergy;
	std::vector<bool> effNevents;
	std::vector< std::vector<bool> > effExp; 
	std::vector< std::vector<bool> > effInt; 
	std::vector<bool> isFixed; 
	std::vector<double> sysSigmas; 
	std::vector<double> sysPosCorr; 
	std::vector<double> sysNegCorr; 
	std::vector<double> sysInitValues; 
	std::vector<double> sysInitSigmas; 
	std::vector<double> sysGenValues; 
	std::vector<double> sysStepSizeRatio; 
	std::vector<double> sysStepSize; 
	std::vector<int> sysScan;
	
	//Derived options
	bool makeScan;

	// Reweight store:
	std::vector<double> nueO_0sigma, nueO_p1sigma, nueO_m1sigma;
	std::vector<double> sys_0sigma, sys_p1sigma, sys_m1sigma;

	// Variable!!
	TFile* outFile;
	TString outPdfFile;
	TCanvas* canny;
	TString outPdfFileEvisSys, outPdfFileGtagSys, outPdfFileNtagSys;
	TString outPdfFileScans;
	TCanvas* cannyEvisSys, *cannyGtagSys, *cannyNtagSys;
	TCanvas* cannyScan;
	std::vector<Experiment*> exps; // Experiments
	FitGlobal* fitter;

	
	std::vector<double> sysCurValues;
	std::vector<double> intCurValues;

private:
	
	//Private Variables 
	
	//Private Functions


};
#endif
