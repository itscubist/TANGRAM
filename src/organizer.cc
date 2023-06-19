/********
*Author: Baran Bodur
*Date: 2022-06-28
*Description: Make PDFs from input tree, to be used in fitting later
*
********/

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
// ROOT headers
#include "TMath.h"
#include "TString.h"

// project headers
#include "organizer.h"
#include "mockData.h"
#include "experiment.h"
#include "fitGlobal.h"
#include "vectorUtils.h"

using namespace std;

// Constructor
Organizer::Organizer(TString inName, TString cardFile, TString outFileName) {
	
	// Set name and read card file
	orgName = inName;
	cout << cardFile << endl;
	cout << outFileName << endl;
	readCardFile(cardFile);
	gRandom = new TRandom3(randSeed);
	cout << "Random Seed: " << randSeed << endl;
	cout << "Rand al'Thor says: " << gRandom->Rndm() << endl;
	printVector(effExp,"Affected Experiments Per Systematic");
	printVector(effInt,"Affected Interactions Per Systematic");
	sysCurValues = sysInitValues;
	intCurValues = intInitValues;

	openOutputs(outFileName);

	// Create and init experiments
	exps.resize(nExp);
	for(unsigned int e=0; e<nExp; e++) {
		exps[e] = new Experiment(this,e);
	}

	// Create fitter
	fitter = new FitGlobal(this);

}

void Organizer::readCardFile(TString cardFile) {
	
	ifstream orgCard(cardFile); // Open organizer card file
	string line; // dummy to read lines in the card
	string name, xx_name, value;
	unsigned int marker;
	unsigned int vecIndex;
	while(!orgCard.eof()){ // Read until file ends

		// Read necessary info from the card file 1 by 1
		getline(orgCard,line); // Get next line
		if(orgCard.eof()) break; // Break if File Ended
		if(line[0]=='#') continue; // Skip comments or empty lines
		marker = line.find(' '); // Find space, where value starts
		value = line.substr(marker+1); // Get value
		name = line.substr(0,marker); // Get full parameter name
		cout << name << " - " << value <<  endl; // Print out read line
		xx_name = name.substr(3); // Get parameter name for numbered entries
		if( isdigit(name[0]) && isdigit(name[1]) ) {
			vecIndex = stoi(name.substr(0,2))-1; // Corresponding vector index for numbered entries
			cout << xx_name << " - " << vecIndex <<  endl; // Print out read name in detail
		}
		else {
			vecIndex = 0;
		}
		
		// ** Assign entries based on card
		
		// Initial Section 
		if(name=="CARD_NAME") cardFileName = value;
		if(name=="GEN_RAND_DATA") genRandData = stoi(value);
		if(name=="N_EVENTS_POISSON") nEventsPoisson = stoi(value);
		if(name=="RAND_SYSTEMATICS") randSystematics = stoi(value);
		if(name=="FIX_SYSTEMATICS") fixSystematics = stoi(value);
		if(name=="DO_MAIN_FIXED_FIT") doMainFixedFit = stoi(value);
		if(name=="DO_REWEIGHT_FIT") doReweightFit = stoi(value);
		if(name=="REWEIGHT_FIT_SYS_INDEX") reweightFitSysIndex = stoi(value);
		if(name=="SAVE_BINNED_FIT_RESULTS") saveBinnedFitResults = stoi(value);
		if(name=="RAND_SEED") randSeed = stoi(value);
		if(name=="N_SUBFIT") nSubFit = stoi(value);
		if(name=="N_FIT") nFit = stoi(value);
		if(name=="N_GRAPH") nGraph = stoi(value);
		if(name=="N_PROFILE") nProfile = stoi(value);
		if(name=="N_PROF_SUBFIT") nProfSubFit = stoi(value);
		if(name=="PROF_INDEX") profIndex = stoi(value);
		if(name=="PROF_START") profStart = stod(value);
		if(name=="PROF_END") profEnd = stod(value);
		if(name=="PROF_STEP") profStep = stod(value);
		if(name=="MAKE_SYS_EFFECT_PLOTS") makeSysEffectPlots = stoi(value);

		// Experiment Section
		if(name=="N_EXP") { // nExp should alwasy be first in this section, resizes vectors correctly
			nExp = stoi(value);
			expNames.resize(nExp);
			expFiles.resize(nExp);
			expTrees.resize(nExp);
			neutronInfo.resize(nExp);
			gammaInfo.resize(nExp);
			dirInfo.resize(nExp);
			evisLow.resize(nExp);
			evisUp.resize(nExp);
			evisLowBuf.resize(nExp);
			evisUpBuf.resize(nExp);
			evisBinsNorm.resize(nExp);
			evisProjRebinFactor.resize(nExp);
			gtagLow.resize(nExp);
			gtagUp.resize(nExp);
			gtagLowBuf.resize(nExp);
			gtagUpBuf.resize(nExp);
			gtagBinsNorm.resize(nExp);
			gtagProjRebinFactor.resize(nExp);
			reweightFactors.resize(nExp);
			dataFiles.resize(nExp);
			dataTrees.resize(nExp);
		}
		if(xx_name=="EXP_NAME") expNames[vecIndex] = (TString) value;
		if(xx_name=="EXP_FILE") expFiles[vecIndex] = (TString) value;
		if(xx_name=="EXP_TREE") expTrees[vecIndex] = (TString) value;
		if(xx_name=="NEUTRON_INFO") neutronInfo[vecIndex] = stoi(value);
		if(xx_name=="GAMMA_INFO") gammaInfo[vecIndex] = stoi(value);
		if(xx_name=="DIR_INFO") dirInfo[vecIndex] = stoi(value);
		if(xx_name=="EVIS_LOW") evisLow[vecIndex] = stod(value);
		if(xx_name=="EVIS_UP") evisUp[vecIndex] = stod(value);
		if(xx_name=="EVIS_LOW_BUF") evisLowBuf[vecIndex] = stod(value);
		if(xx_name=="EVIS_UP_BUF") evisUpBuf[vecIndex] = stod(value);
		if(xx_name=="EVIS_BINS_NORM") evisBinsNorm[vecIndex] = stoi(value);
		if(xx_name=="EVIS_PROJ_REBIN_FACTOR") evisProjRebinFactor[vecIndex] = stoi(value);
		if(xx_name=="GTAG_LOW") gtagLow[vecIndex] = stod(value);
		if(xx_name=="GTAG_UP") gtagUp[vecIndex] = stod(value);
		if(xx_name=="GTAG_LOW_BUF") gtagLowBuf[vecIndex] = stod(value);
		if(xx_name=="GTAG_UP_BUF") gtagUpBuf[vecIndex] = stod(value);
		if(xx_name=="GTAG_BINS_NORM") gtagBinsNorm[vecIndex] = stoi(value);
		if(xx_name=="GTAG_PROJ_REBIN_FACTOR") gtagProjRebinFactor[vecIndex] = stoi(value);
		if(xx_name=="REWEIGHT_FACTOR") reweightFactors[vecIndex] = stod(value);
		if(xx_name=="DATA_FILE") dataFiles[vecIndex] = (TString) value;
		if(xx_name=="DATA_TREE") dataTrees[vecIndex] = (TString) value;

		// Interaction section
		if(name=="N_INT") { // nInt should always be first in this section, resizes vectors correctly
			nInt = stoi(value);
			intNames.resize(nInt);
			intColors.resize(nInt);
			intVarFraction.resize(nInt);
			intInitValues.resize(nInt);
			intInitSigmas.resize(nInt);
			intIsFixed.resize(nInt);
			evisSubbinsForGtag.resize(nInt);
			intNameCutLow.resize(nInt);
			intNameCutUp.resize(nInt);
			nmueCutLow.resize(nInt);
			nmueCutUp.resize(nInt);
			nuEneCutLow.resize(nInt);
			nuEneCutUp.resize(nInt);
			wallCutLow.resize(nInt);
			wallCutUp.resize(nInt);
			intStepSizeRatio.resize(nInt);
			intStepSize.resize(nInt);
			evisKdeOptStr.resize(nInt);
			evisKdeRho.resize(nInt);
			evisKdeLow.resize(nInt);
			gtagKdeOptStr.resize(nInt);
			gtagKdeRho.resize(nInt);
			gtagKdeLow.resize(nInt);
			intScan.resize(nInt);
		}
		if(xx_name=="INT_NAME") intNames[vecIndex] = (TString) value;
		if(xx_name=="INT_COLOR") intColors[vecIndex] = stoi(value);
		if(xx_name=="INT_VAR_FRAC") intVarFraction[vecIndex] = stod(value);
		if(xx_name=="INT_INIT_VAL") intInitValues[vecIndex] = stod(value);
		if(xx_name=="INT_INIT_SIGMA") intInitSigmas[vecIndex] = stod(value);
		if(xx_name=="INT_IS_FIXED") intIsFixed[vecIndex] = stoi(value);
		if(xx_name=="EVIS_SUBBINS_FOR_GTAG") evisSubbinsForGtag[vecIndex] = stoi(value);
		if(xx_name=="INT_NAME_CUT_LOW") intNameCutLow[vecIndex] = stoi(value);
		if(xx_name=="INT_NAME_CUT_UP") intNameCutUp[vecIndex] = stoi(value);
		if(xx_name=="NMUE_CUT_LOW") nmueCutLow[vecIndex] = stoi(value);
		if(xx_name=="NMUE_CUT_UP") nmueCutUp[vecIndex] = stoi(value);
		if(xx_name=="NU_ENE_CUT_LOW") nuEneCutLow[vecIndex] = stod(value);
		if(xx_name=="NU_ENE_CUT_UP") nuEneCutUp[vecIndex] = stod(value);
		if(xx_name=="WALL_CUT_LOW") wallCutLow[vecIndex] = stod(value);
		if(xx_name=="WALL_CUT_UP") wallCutUp[vecIndex] = stod(value);
		if(xx_name=="INT_STEP_SIZE_RATIO") intStepSizeRatio[vecIndex] = stod(value);
		if(xx_name=="EVIS_KDE_OPT_STR") evisKdeOptStr[vecIndex] = (TString) value;
		if(xx_name=="EVIS_KDE_RHO") evisKdeRho[vecIndex] = stod(value);
		if(xx_name=="EVIS_KDE_LOW") evisKdeLow[vecIndex] = stod(value);
		if(xx_name=="GTAG_KDE_OPT_STR") gtagKdeOptStr[vecIndex] = (TString) value;
		if(xx_name=="GTAG_KDE_RHO") gtagKdeRho[vecIndex] = stod(value);
		if(xx_name=="GTAG_KDE_LOW") gtagKdeLow[vecIndex] = stod(value);
		if(xx_name=="INT_SCAN") intScan[vecIndex] = stoi(value);
		
		// Systematics section
		if(name=="N_SYS") { // nSys should always be first in this section, resizes vectors correctly
			nSys = stoi(value);
			sysNames.resize(nSys);
			shiftOrScale.resize(nSys);
			weightNames.resize(nSys);
			effNeutron.resize(nSys);
			effGamma.resize(nSys);
			effEnergy.resize(nSys);
			effNevents.resize(nSys);
			effExp.resize(nSys);
			effInt.resize(nSys);
			isFixed.resize(nSys);
			sysSigmas.resize(nSys);
			sysPosCorr.resize(nSys);
			sysNegCorr.resize(nSys);
			sysInitValues.resize(nSys);
			sysInitSigmas.resize(nSys);
			sysGenValues.resize(nSys);
			sysStepSizeRatio.resize(nSys);
			sysStepSize.resize(nSys);
			sysScan.resize(nSys);
		}
		if(name=="SYS_WEIGHT_MIN_SCALE") sysWeightMinScale = stod(value);
		if(name=="SYS_WEIGHT_MAX_SCALE") sysWeightMaxScale = stod(value);
		if(name=="SYS_WEIGHT_MIN_SHIFT") sysWeightMinShift = stod(value);
		if(name=="SYS_WEIGHT_MAX_SHIFT") sysWeightMaxShift = stod(value);
		if(name=="SYS_LOW_LIMIT") sysLowLimit = stod(value);
		if(name=="SYS_UP_LIMIT") sysUpLimit = stod(value);
		if(xx_name=="SYS_NAME") sysNames[vecIndex] = (TString) value;
		if(xx_name=="SHIFT_OR_SCALE") shiftOrScale[vecIndex] = stoi(value);
		if(xx_name=="WEIGHT_NAME") weightNames[vecIndex] = (TString) value;
		if(xx_name=="EFF_NEUTRON") effNeutron[vecIndex] = stoi(value);
		if(xx_name=="EFF_GAMMA") effGamma[vecIndex] = stoi(value);
		if(xx_name=="EFF_ENERGY") effEnergy[vecIndex] = stoi(value);
		if(xx_name=="EFF_NEVENTS") effNevents[vecIndex] = stoi(value);
		if(xx_name=="EFF_EXP") { // This one has nExp arguments so requires further work
			vector<bool> temp; temp.resize(nExp);
			if(nExp==1) temp[0] = stoi(value);
			else {
				unsigned int start = 0, end;
				for(unsigned int i = 0; i<nExp; i++) {
					end = value.find(' ',start);	
					temp[i] = stoi(value.substr(start,start-end));
					start = end + 1;
				}
			}
			effExp[vecIndex] = temp;
		}
		if(xx_name=="EFF_INT") {
			vector<bool> temp; temp.resize(nInt);
			if(nInt==1) temp[0] = stoi(value);
			else {
				unsigned int start = 0, end;
				for(unsigned int i = 0; i<nInt; i++) {
					end = value.find(' ',start);		
					temp[i] = stoi(value.substr(start,start-end));
					start = end + 1;
				}
			}
			effInt[vecIndex] = temp;
		}
		if(xx_name=="IS_FIXED") isFixed[vecIndex] = stoi(value);
		if(xx_name=="SYS_SIGMA") sysSigmas[vecIndex] = stod(value);
		if(xx_name=="SYS_POS_CORR") sysPosCorr[vecIndex] = stod(value);
		if(xx_name=="SYS_NEG_CORR") sysNegCorr[vecIndex] = stod(value);
		if(xx_name=="SYS_INIT_VAL") sysInitValues[vecIndex] = stod(value);
		if(xx_name=="SYS_INIT_SIGMA") sysInitSigmas[vecIndex] = stod(value);
		if(xx_name=="SYS_GEN_VAL") sysGenValues[vecIndex] = stod(value);
		if(xx_name=="SYS_STEP_SIZE_RATIO") sysStepSizeRatio[vecIndex] = stod(value);
		if(xx_name=="SYS_SCAN") sysScan[vecIndex] = stoi(value);


	} // End of file read

	// Fill derived quantities
	makeScan == false;
	if( getSum(sysScan) + getSum(intScan) > 0 ) makeScan = true;
	sysStepSize = sysSigmas*sysStepSizeRatio;
	intStepSize = intInitValues*intStepSizeRatio;
	
	orgCard.close();
	return;	
}

// Generates requested number of mock datasets for all experiment/interations
// Systematics used in generation is consistent (even if randSystematic is true) within a dataset
void Organizer::genMockData() {
	
	// Prepare systematics vectors to be used:
	vector< vector<double> > usedSysVec;
	usedSysVec.resize(nFit);
	for(unsigned int f=0; f<nFit; f++) {
		usedSysVec[f].resize(nSys);
		for(unsigned int s=0; s<nSys; s++) {
			if(randSystematics == false || isFixed[s]==1 ) usedSysVec[f][s] = sysGenValues[s];
			else if(randSystematics == true) { 
				//usedSysVec[f][s] = gRandom->Gaus(sysGenValues[s], sysSigmas[s]);
				double sysVal = gRandom->Gaus(sysGenValues[s], sysSigmas[s]);
				while(sysVal<-2.0*sysSigmas[s] || sysVal>2.0*sysSigmas[s]) {
					sysVal = gRandom->Gaus(sysGenValues[s], sysSigmas[s]);
				}
				usedSysVec[f][s] = sysVal;
			}
			cout << "Fit: "<< f << " " << sysNames[s] << " GenMean: " << sysGenValues[s]
			     << " GenValue: " << usedSysVec[f][s] << endl;
		}
	}

	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->genMockData(usedSysVec);
	}
}

// Reweights 1 systematic with given index and sigma to all systematics at 0 case
void Organizer::reweightOneSys(unsigned int sysIndex, double nSigma, unsigned int mockIndexIn) {
	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->reweightOneSys(sysIndex, nSigma, mockIndexIn);
	}
}

// Does reweight fit:
void Organizer::reweightFit() {
	double nSigma = 1.0;
	double nSigmaM = -1.0;
	unsigned int sysIndex = reweightFitSysIndex;
	fitter->prepFitHists();
	for(unsigned int f=0; f<nFit; f++) {
		//Fit before reweight
		fitter->resetFit();
		fitter->setMockIndex(f);
		fitter->makeFit(fixSystematics, false);
		nueO_0sigma.push_back( fitter->fitResults[f].fitFracs[0] );
		sys_0sigma.push_back( fitter->fitResults[f].fitSys[reweightFitSysIndex] );
		fitter->fillFitHists();
		//Fit after reweight
		reweightOneSys(sysIndex, nSigma, f);
		fitter->resetFit();
		fitter->setMockIndex(f);
		fitter->makeFit(fixSystematics, false);
		nueO_p1sigma.push_back( fitter->fitResults[f].fitFracs[0] );
		sys_p1sigma.push_back( fitter->fitResults[f].fitSys[reweightFitSysIndex] );
		//Fit after reweight
		reweightOneSys(sysIndex, nSigmaM, f);
		fitter->resetFit();
		fitter->setMockIndex(f);
		fitter->makeFit(fixSystematics, false);
		nueO_m1sigma.push_back( fitter->fitResults[f].fitFracs[0] );
		sys_m1sigma.push_back( fitter->fitResults[f].fitSys[reweightFitSysIndex] );
	}
	fitter->saveFitHists();
}

// Updates current systematics
void Organizer::updateCurSys(vector<double> inSys) {
	if(inSys.size() != nSys) {
		cout << "WARNING! In vector size: " << inSys.size() << " does not match nSys: " << nSys << endl;
		return;
	}
	for(unsigned int s=0; s<nSys; s++) {
		sysCurValues[s] = inSys[s];
		if( sysCurValues[s]!=sysCurValues[s] ) {
			cout << "WARNING! Systematic: " << sysNames[s] << " Value: " << sysCurValues[s] << endl; 
			cout << "Assigning 0 to this systematic!" << endl;
			sysCurValues[s] = 0;
		}
	}
}

// Updates all pdfs with current systematics
void Organizer::updateAllPdfs() {
	cout << "Updating all PDFs using systematics: " << endl;
	for(unsigned int s=0; s<nSys; s++) {
		cout << "Systematic: " << sysNames[s] << " Value: " << sysCurValues[s] << endl; 
	}
	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->updateAllPdfs();
	}
}

// Only Profile Fit
void Organizer::onlyProfileFit() {
	vector<double> profSteps = fillLinear(profStart,profEnd+0.5*profStep,profStep);
	printVector(profSteps,"Profile Step Vector");
	for(unsigned int f=0; f<nFit; f++) {
		fitter->resetFit();
		fitter->setMockIndex(f);
		if(f<nProfile) {
			//fitter->makeProfileFit(profIndex,profSteps,true);	//fixed systematics at init val
			//fitter->resetFit();
			fitter->makeProfileFit(profIndex,profSteps,false); //free systematics
			saveCurPdfs( Form("mockFit%d",f) );
			updateProjHists(1.0);
			saveCurProj( Form("mockFit%d",f) );
		}
	}
	fitter->saveProfGraphs();
}

// Make fits to all mock data
void Organizer::makeAllFits() {
	
	
	vector<double> profSteps = fillLinear(profStart,profEnd+0.5*profStep,profStep);
	printVector(profSteps,"Profile Step Vector");
	fitter->prepFitHists();
	for(unsigned int f=0; f<nFit; f++) {
		fitter->resetFit();
		fitter->setMockIndex(f);
		if(f<nProfile) {
			fitter->makeProfileFit(profIndex,profSteps,fixSystematics);	
		}
		fitter->makeFit(fixSystematics, false);
		if(doMainFixedFit==true) fitter->makeFit(fixSystematics, true);
		if(f<nGraph) {
			fitter->printFitGraphs(f,"fit");
			fitter->printFitGraphs(f,"sim");
			fitter->printFitGraphs(f,"pois");	
		}
		saveCurPdfs( Form("mockFit%d",f) );
		updateProjHists(1.0);
		saveCurProj( Form("mockFit%d",f) );
		
		// If requested make binned chi2 test with data (both sim and fit)
		if( saveBinnedFitResults == true ) {
			fitter->getFitChiAndHists(f,"fit");
			fitter->getFitChiAndHists(f,"sim");
			fitter->getFitChiAndHists(f,"pois");
			cout << endl << "*** PRINTING CHI2 FROM BINNED CHECK! *** " << endl << endl;
			// Calculate NDF:
			double ndf = 0;
			for(unsigned int e=0; e<nExp; e++) {
				ndf += (double) evisBinsNorm[e]/evisProjRebinFactor[e] ;
				ndf += (double) gtagBinsNorm[e]/gtagProjRebinFactor[e] ;
				ndf -= 1.0;
				if(neutronInfo[e]==true) {
					ndf += (double) evisBinsNorm[e]/evisProjRebinFactor[e] ;
					ndf += (double) gtagBinsNorm[e]/gtagProjRebinFactor[e] ;	
					ndf -= 1.0;
				}
			}
			cout << " NDF: " << ndf << endl; 
			cout << " Fit Chi2: " << fitter->fitResults[f].fitChi2 << endl;
			cout << " Sim True Chi2: " << fitter->fitResults[f].simTrueChi2 << endl;
			cout << " Sim Pois Chi2: " << fitter->fitResults[f].simPoisChi2 << endl;
			cout << endl << "**************************************** " << endl << endl;
		}
		
		fitter->fillFitHists();
		// Scan requested parameters around the best fit point
		if( makeScan == true ) fitter->scanAfterFit();
	
	} // end of loop
	

	fitter->saveFitHists();
	fitter->saveProfGraphs();
}
	
// Save current PDF histograms for all experiment interactions
void Organizer::saveCurPdfs(TString saveName) {
	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->saveCurPdfs(saveName);
	}
}

// Updates all projection hists with current PDFs
void Organizer::updateProjHists(double scale) {
	cout << "Updating all projection histograms" << endl;
	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->updateProjHists(scale);
	}
}

// Save current Proj histograms for all experiment interactions
void Organizer::saveCurProj(TString saveName) {
	for(unsigned int e=0; e<nExp; e++) {
		exps[e]->saveCurProj(saveName);
	}
}

// Open root file and pdf file (canvas) to write
void Organizer::openOutputs(TString outFileName) {
	outFile = new TFile(outFileName+".root","RECREATE");
	outPdfFile = outFileName+".pdf";
	
	// Calculate columns and rows, and create TCanvas accordingly
	unsigned int nRows = 4;
	unsigned int nCols = 0;
	for(unsigned int e=0;e<nExp;e++) {
		if(neutronInfo[e]==false) nCols += 1;
		if(neutronInfo[e]==true) nCols += 3;
	}
	unsigned int widthX1 = 800, heightX1 = 600;
	canny = new TCanvas("canny","canny",widthX1*nCols,heightX1*nRows);
	canny->Divide(nCols,nRows);
	canny->Print(outPdfFile+"[");


	// For systematic plots:
	if( makeSysEffectPlots == true) {
		outPdfFileEvisSys = outFileName+"_EvisSys1sigma.pdf";
		outPdfFileGtagSys = outFileName+"_GtagSys1sigma.pdf";
		outPdfFileNtagSys = outFileName+"_NtagSys1sigma.pdf";
		unsigned int nRowsSys = 4;
		unsigned int nColsSys = 2;
		unsigned int widthX2 = 800, heightX2 = 1200;
		cannyEvisSys = new TCanvas("cannyEvisSys","cannyEvisSys",widthX2,heightX2);
		cannyEvisSys->Divide(nColsSys,nRowsSys);
		cannyEvisSys->Print(outPdfFileEvisSys+"[");
		cannyGtagSys = new TCanvas("cannyGtagSys","cannyGtagSys",widthX2,heightX2);
		cannyGtagSys->Divide(nColsSys,nRowsSys);
		cannyGtagSys->Print(outPdfFileGtagSys+"[");
		cannyNtagSys = new TCanvas("cannyNtagSys","cannyNtagSys",widthX2,heightX2);
		cannyNtagSys->Divide(nColsSys,nRowsSys);
		cannyNtagSys->Print(outPdfFileNtagSys+"[");
	}
	
	// For scanning
	outPdfFileScans = outFileName+"_ParScans.pdf";
	unsigned int nRowsScan = 10;
	unsigned int nColsScan = 5;
	unsigned int widthX3 = 800, heightX3 = 1200;
	cannyScan = new TCanvas("cannyScan","cannyScan",widthX3,heightX3);
	cannyScan->Divide(nColsScan,nRowsScan);
	if( makeScan==true ) cannyScan->Print(outPdfFileScans+"[");


}
	
// Close root file and pdf file (canvas) to write
void Organizer::closeOutputs() {
	outFile->Close();
	canny->Print(outPdfFile+"]");
	delete canny;
	if( makeSysEffectPlots == true ) {
		cannyEvisSys->Print(outPdfFileEvisSys+"]");
		cannyGtagSys->Print(outPdfFileGtagSys+"]");
		cannyNtagSys->Print(outPdfFileNtagSys+"]");
		delete cannyEvisSys;
		delete cannyGtagSys;
		delete cannyNtagSys;
	}
	if( makeScan==true ) {
		cannyScan->Print(outPdfFileScans+"]");
		delete cannyScan;
	}
	
	cout << "Closed and deleted output tools" << endl;
	// Print reweight results if doing reweight fitting
	if(doReweightFit == 1) {
		// Get mean of all reweight fit results
		double sys_0sigma_mean = getMean(sys_0sigma);
		double sys_p1sigma_mean = getMean(sys_p1sigma);
		double sys_m1sigma_mean = getMean(sys_m1sigma);
		double nueO_0sigma_mean = getMean(nueO_0sigma);
		double nueO_p1sigma_mean = getMean(nueO_p1sigma);
		double nueO_m1sigma_mean = getMean(nueO_m1sigma);

		cout << endl << "****************************" << endl << endl;
		cout << "Reweight Sys Name: " << sysNames[reweightFitSysIndex] << endl;
		cout << "sys fit at 0 sigma: " << sys_0sigma_mean << endl;
		cout << "sys fit at +1 sigma: " << sys_p1sigma_mean << endl;
		cout << "sys fit at -1 sigma: " << sys_m1sigma_mean << endl;
		cout << endl << " ***** Now its effect on nueO result *****" << endl << endl; 
		cout << "nueO fit at 0 sigma: " << nueO_0sigma_mean << endl;
		cout << "nueO fit at +1 sigma: " << nueO_p1sigma_mean << endl;
		cout << "nueO fit at -1 sigma: " << nueO_m1sigma_mean << endl;

		// Calculate the mean of nueO results and then the effect:
		double plusPerDiffOfMean = 100*(nueO_p1sigma_mean - nueO_0sigma_mean)/nueO_0sigma_mean;
		double minusPerDiffOfMean = 100*(nueO_m1sigma_mean - nueO_0sigma_mean)/nueO_0sigma_mean;
		double maxPerDiffOfMean = abs(plusPerDiffOfMean) > abs(minusPerDiffOfMean) ?
			abs(plusPerDiffOfMean) : abs(minusPerDiffOfMean);

		cout << "***** 100*(<nueO fit+-> - <nueO fit0>)/<nueO fit0>: *****" << endl << endl;
		cout << "Percent diff at +1 sigma: " << plusPerDiffOfMean << endl;
		cout << "Percent diff at -1 sigma: " << minusPerDiffOfMean << endl;
		cout << "Max absolute percent effect: " << maxPerDiffOfMean << endl;
		

		// First calculate the effect and then get the mean
		vector<double> plusPerDiff = (nueO_p1sigma - nueO_0sigma)/nueO_0sigma;
		vector<double> minusPerDiff = (nueO_m1sigma - nueO_0sigma)/nueO_0sigma;

		double plusPerDiff_mean = 100 * getMean(plusPerDiff);
		double minusPerDiff_mean = 100 * getMean(minusPerDiff);

		double maxPerDiff_mean = abs(plusPerDiff_mean) > abs(minusPerDiff_mean) ? 
			abs(plusPerDiff_mean) : abs(minusPerDiff_mean);

		cout << "***** <Percent Effect on nueO fit>: *****" << endl << endl;
		cout << "Percent diff at +1 sigma: " << plusPerDiff_mean << endl;
		cout << "Percent diff at -1 sigma: " << minusPerDiff_mean << endl;
		cout << "Max absolute percent effect: " << maxPerDiff_mean << endl;
	}
}

// Print 1 sigma systematic changes for 1 experiment
void Organizer::printSysHists(unsigned int expIndex, double sysShift) {
	if(expIndex==1000) {
		for(unsigned int e=0;e<nExp;e++) {
			exps[e]->printSysHists(sysShift);
		}
	}
	else {
		exps[expIndex]->printSysHists(sysShift);
	}
}

void Organizer::printParScans(TString outFileName) {
	//For variable scans:
	outPdfFileScans = outFileName+"_ParScans.pdf";
	unsigned int nRowsScan = 10;
	unsigned int nColsScan = 5;
	unsigned int widthX3 = 800, heightX3 = 1200;
	cannyScan = new TCanvas("cannyScan","cannyScan",widthX3,heightX3);
	cannyScan->Divide(nColsScan,nRowsScan);
	cannyScan->Print(outPdfFileScans+"[");
	
	for(unsigned int f=0; f<nFit; f++) {
		fitter->resetFit();
		fitter->setMockIndex(f);
		fitter->scanAll();
	}

	cannyScan->Print(outPdfFileScans+"]");
	delete cannyScan;
}

void Organizer::checkSysEffect(unsigned int sysIndex, unsigned int pdfIndex) {
	for(unsigned int e=0;e<nExp;e++) {
		exps[e]->checkSysEffect(sysIndex, pdfIndex);
	}
}

// To Get Binned Chi2 Between Nominal and Varied Systematic
void Organizer::getChi2SysDiff(unsigned int sysNo, double sysVal) {
	
	vector< double > intFracs; intFracs.resize(nInt);
	vector< double > sysVals; sysVals.resize(nSys);
	for(unsigned int i=0;i<nInt;i++) intFracs[i] = intInitValues[i];	
	for(unsigned int s=0;s<nSys;s++) sysVals[s] = sysGenValues[s];
	sysVals[sysNo] = sysVal*sysSigmas[sysNo];
	
	vector< double > intFracs0; intFracs0.resize(nInt);
	vector< double > sysVals0; sysVals0.resize(nSys);
	for(unsigned int i=0;i<nInt;i++) intFracs0[i] = intInitValues[i];	
	for(unsigned int s=0;s<nSys;s++) sysVals0[s] = sysGenValues[s];

	TString nameStart, nameStart0;
	nameStart = Form("hTotal_Sys%d_%fSigma",sysNo,sysVal);
	nameStart0 = Form("hTotal_Sys%d_0Sigma",sysNo);
	
	// Temporaries for looping over
	vector<TH1D*> tempHists;
	vector<TH1D*> tempHists0;
	unsigned int thInd = 0;
	double chi2Val = 0;
	double ndf = 0;
	unsigned int nBins = 0;
	
	// First loop for eVis histograms
	for(unsigned int e=0;e<nExp;e++) {
		unsigned int nMax;
		unsigned int nMin;
		if(neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(neutronInfo[e]==true) {nMin=1; nMax = 3;}
		for(unsigned int nn=nMin; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "eVis";
			tempHists.push_back( 
					exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			tempHists0.push_back( 
					exps[e]->makeTotalHist( intFracs0, sysVals0, hType, tN, tE, nameStart0 ) );
			nBins = tempHists[thInd]->GetNbinsX(); 
			for(unsigned int b=0;b<nBins;b++) {
				double theoVal = tempHists0[thInd]->GetBinContent(b+1);
				double poisVal = tempHists[thInd]->GetBinContent(b+1);
				double poisErrSqr = 1;
				if(poisVal>0 && theoVal>0) {
					chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
				}	
			}
		//cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
			thInd++;
		}
	}
	// Second loop to fill 2nd row with gTag histograms
	for(unsigned int e=0;e<nExp;e++) {
		if(gammaInfo[e]==false) continue;
		unsigned int nMax;
		unsigned int nMin;
		if(neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(neutronInfo[e]==true) {nMin=1; nMax = 3;}
		for(unsigned int nn=nMin; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "gTag";
			tempHists.push_back( 
					exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			tempHists0.push_back( 
					exps[e]->makeTotalHist( intFracs0, sysVals0, hType, tN, tE, nameStart0 ) );
			nBins = tempHists[thInd]->GetNbinsX();
			for(unsigned int b=0;b<nBins;b++) {
				double theoVal = tempHists0[thInd]->GetBinContent(b+1);
				double poisVal = tempHists[thInd]->GetBinContent(b+1);
				double poisErrSqr = 1;
				if(poisVal>0 && theoVal>0) {
					chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
				}	
			}
		//cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
			thInd++;
		}
	}
	// Third loop to fill 3rd row with nTag histograms
	for(unsigned int e=0;e<nExp;e++) {
		unsigned int nMax;
		unsigned int nMin;
		if(neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(neutronInfo[e]==true) {nMin=1; nMax = 3;}
		unsigned int tN = 0;
		unsigned int tE = 0;
		TString hType = "nTag";
		tempHists.push_back(
				exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
		tempHists0.push_back(
				exps[e]->makeTotalHist( intFracs0, sysVals0, hType, tN, tE, nameStart0 ) );
		nBins = tempHists[thInd]->GetNbinsX(); 
		for(unsigned int b=0;b<nBins;b++) {
			double theoVal = tempHists0[thInd]->GetBinContent(b+1);
			double poisVal = tempHists[thInd]->GetBinContent(b+1);
			double poisErrSqr = 1;
			if(poisVal>0 && theoVal>0) {
				chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
			}	
		}
	//cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
		thInd++;
	}

	outFile->cd();
	for(unsigned int v=0;v<tempHists.size();v++) {
		tempHists0[v]->Write();
		delete tempHists0[v];
		tempHists[v]->Write();
		delete tempHists[v];
	}

	cout << endl << "*** SYSTEMATIC: " << sysNames[sysNo] << " ***" << endl << endl;
	cout << " Value of systematic: " << sysVal << endl;
	cout << " Chi2 To systematic at 0: " << chi2Val << endl;
	cout << endl << endl;

} // end of getChi2SysDiff
