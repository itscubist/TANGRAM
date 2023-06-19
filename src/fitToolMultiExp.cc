/********
*Author: Baran Bodur
*Date: 2022-03-24
*Description: To fit data from multiple SK phases together
*
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
#include "fitPdf.h"
#include "mockData.h"
#include "fitTool.h"
#include "fitToolMultiExp.h"

using namespace std;

// Constructor
FitToolMultiExp::FitToolMultiExp(TString inName, const unsigned int nPdf, 
		const unsigned int nFitIn, FitTool *inFitter[]) {
	
	// Set namd & pdf no
	name = inName;
	nFitter = nFitIn;
	nType = nPdf;


	// Reserve spaces for vectors
	fitSteps.resize(nType);
	curPars.resize(nType);
	initPars.resize(nType);
	
	h1FitDist.reserve(nType);

	// Save pdfs in this object
	for(unsigned int f=0; f<nFitter; f++) {
		fitters.push_back(inFitter[f]);			
	}
	for(unsigned int p=0; p<nType; p++) {
		fitPars.push_back(0);
		simTruePars.push_back(0);
		simPoisPars.push_back(0);
	}

	// Create minimizer object
	minFarshawMultiExp = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
	minFarshawMultiExp->SetMaxFunctionCalls(100000); // max calls, reduce if it takes too much time 
	minFarshawMultiExp->SetMaxIterations(30000); // not used
	minFarshawMultiExp->SetTolerance(0.05); // tolerance for EDM
	minFarshawMultiExp->SetStrategy(2); // 0:coarse 1: normal 2: detailed 
	minFarshawMultiExp->SetPrintLevel(0); 

	
	// Set initial pars and fit steps (temporary can be set later externally
	double tempInit[nPdf];
	double tempStep[nPdf];
	for(unsigned int p=0; p<nType; p++) {
		tempInit[p] = 1;
		tempStep[p] = 0.01;
	}
	setInitPars(tempInit);
	setFitStep(tempStep);
	setSimPars(); // based on pdfs ( sets only simTruePars, not simPoisPars)

	// Create functor to get value & gradient of the function to be minimized at given point
	getLhMultiExp = &FitToolMultiExp::getLikelihood;
	getGradMultiExp = &FitToolMultiExp::getLhGradient;

	// functor to be used in minimization	
	fitFuncMultiExp = ROOT::Math::GradFunctor(this,getLhMultiExp,getGradMultiExp,nType); 
	//cout << "Dimensions: " <<fitFuncMultiExp.NDim() << endl;

	// Bind the fit function to the minimizer
	minFarshawMultiExp->SetFunction(fitFuncMultiExp);
	
	// Bind variables
	for(unsigned int p = 0; p < nType; p++) {
			minFarshawMultiExp->SetLowerLimitedVariable(p,(string)fitters[0]->pdfs[p]->pdfName,
					tempInit[p],tempStep[p],0.0);
	}

	//Prepare fit histograms
	prepFitHists();
	
	stackHistCtr = 0;
	profileNllCtr = 0;

} // end of constructor

//Destructor
FitToolMultiExp::~FitToolMultiExp() {
	/*
	delete minFarshawMultiExp;
	
	if(histsReady == true) {
		delete hFitLh;
		delete hSimPoisLh;
		delete hSimTrueLh;
		delete hLhFitMinPois;
		delete hLhFitMinTrue;
		for(unsigned int h = 0; h < nType; h++) {
			delete h1FitDist[h];
			for(unsigned int h2 = 0; h2 < h; h2++) {
				delete h2FitDist[h][h2];
			}
		} // end of for
	} // end of if histograms ready

	delStackedHists(); // delete stacked hists
	*/
} // end of destructor

// Tool to reset fit function, minimizer and fit step and initial guess
void FitToolMultiExp::resetFit() { // sets fit function, init vals, steps again
	
	fitFuncMultiExp = ROOT::Math::GradFunctor(this,getLhMultiExp,getGradMultiExp,nType);	
	// Bind the fit function to the minimizer
	minFarshawMultiExp->SetFunction(fitFuncMultiExp);
	// Bind variables
	for(unsigned int p = 0; p < nType; p++) { // set vars again
			minFarshawMultiExp->SetLowerLimitedVariable(p,(string)fitters[0]->pdfs[p]->pdfName,
					initPars[p],fitSteps[p],0.0);
	} // end of loop over fit vars
} // end of resetFit

// Set data to individual fitters
void FitToolMultiExp::setData(const unsigned int nFitIn, MockData *inData[]) {
	if( nFitIn != nFitter ) {
		cout << " Data/Fitter Dimensions Do No Match, Fitting will not be successful!"; 
		return;
	}
	for(unsigned int f=0; f<nFitter; f++) {
		fitters[f]->setData(inData[f]);	// set each fitter with its data	
	}
	dataWasSet = true; // set dataWasSet flag
	setSimPars(); // set simulated parameters
} // end of setData

// set initial parameters
void FitToolMultiExp::setInitPars(double *inPars) {
	for(unsigned int p=0; p<nType; p++) initPars[p] = inPars[p];	
} // end of setInitPars

// set fit step in directions
void FitToolMultiExp::setFitStep(double *inSteps) {
	for(unsigned int p=0; p<nType; p++) fitSteps[p] = inSteps[p];	
} // end of setFitStep

// setSimPars from input pdfs (simTruth) and mock data (sim Result with fluctuations)
void FitToolMultiExp::setSimPars() { 	
	for(unsigned int p=0; p<nType; p++) {
		//simTruePars[p] = pdfs[p]->expEvents; 
		simTruePars[p] = 1.0; // every type of interaction gets expected # of events
		simPoisPars[p] = 0.0; 
		double denomPoisPars = 0;
		if(dataWasSet == true) {
			for(unsigned int f=0; f<nFitter; f++) {
				simPoisPars[p] += fitters[f]->dataToFit->byTypeEvents[p];
				denomPoisPars += fitters[f]->pdfs[p]->expEvents;
			}
			simPoisPars[p] /= denomPoisPars; // something around 1, # of events/expected events per int.
		}
	}
	simTrueLh = getLikelihood(simTruePars.data()); // get likelihood with true sim Params
	if(dataWasSet == true) simPoisLh = getLikelihood(simPoisPars.data()); //with fluctuated params
} // end of setSimPars

// run the fit (save result) : makeFit
double FitToolMultiExp::makeFit() {
	
	double initVal = 1.0;
	bool refitFlag = false;
	double refitTh = 0.05;
	double refitInit[nType];
		
	for(unsigned int p=0; p<nType; p++) minFarshawMultiExp->SetVariableValue(p,initVal);
	minFarshawMultiExp->Minimize();
	const double *fitResult = minFarshawMultiExp->X();
	const double *fitErr = minFarshawMultiExp->Errors();
	for(unsigned int p=0; p<nType; p++) { // set fit results
		fitPars[p] = fitResult[p];	
		refitInit[p] = fitPars[p];
		if( fitPars[p] < refitTh ) {
			refitFlag = true; 
			refitInit[p] = initVal;
		}
		cout << fitters[0]->pdfs[p]->pdfName << " fit To: " << fitPars[p] << " +/- "
			   << fitErr[p] << endl;
		for(unsigned int f=0; f<nFitter; f++) { // set fit results for individual fitters
			fitters[f]->fitPars[p] = fitters[f]->pdfs[p]->expEvents * fitPars[p];
		}
	}
	for(unsigned int f=0; f<nFitter; f++) { // set fit results for individual fitters
		fitters[f]->fitMinLh = fitters[f]->getLikelihood(fitters[f]->fitPars.data()); 
	}
	cout << "FitParsVec size: " << fitPars.size() << endl;
	//cout << "True Pars 1: " << *(simTruePars.data()+3) << endl;
	cout << "# of function calls for min:" << minFarshawMultiExp->NCalls() << endl;
	fitMinLh = minFarshawMultiExp->MinValue();

	refitFlag = false; // prevent refitting (hardcoded)
	if( refitFlag == true ) { // refit if too close to boundary value
		cout << " ** Starting refitting, because a variable was too close to the boundary!" << endl;
		for(unsigned int p=0; p<nType; p++) minFarshawMultiExp->SetVariableValue(p,refitInit[p]);
		minFarshawMultiExp->Minimize();
		fitResult = minFarshawMultiExp->X();
		for(unsigned int p=0; p<nType; p++) { // set fit results
			fitPars[p] = fitResult[p];	
			cout << fitters[0]->pdfs[p]->pdfName << " fit To: " << fitPars[p] << " +/- " 
					 << fitErr[p] << endl;
			for(unsigned int f=0; f<nFitter; f++) { // set fit results for individual fitters
				fitters[f]->fitPars[p] = fitters[f]->pdfs[p]->expEvents * fitPars[p];
			}
		}
		for(unsigned int f=0; f<nFitter; f++) { // set fit results for individual fitters
			fitters[f]->fitMinLh = fitters[f]->getLikelihood(fitters[f]->fitPars.data()); 
		}
		cout << "FitParsVec size: " << fitPars.size() << endl;
		//cout << "True Pars 1: " << *(simTruePars.data()+3) << endl;
		cout << "# of function calls for 2nd min:" << minFarshawMultiExp->NCalls() << endl;
		fitMinLh = minFarshawMultiExp->MinValue();
	} // end of refit if

	return fitMinLh;
} // end of makeFit


// Fix the parameter of interest at various values, and fit the other parameters, 
// effectively profiling over the other parameters
double FitToolMultiExp::profileFit(int varIndex, vector<double> profVals) {
	double initVal = 1.0;
	unsigned int nVal = profVals.size();
	vector<double> nllVals;
	double minNll = 1e10;
	unsigned int minIndex = 0;
	for(unsigned int v = 0; v < nVal; v++){
		for(unsigned int p=0; p<nType; p++) {
			if( p == varIndex ) minFarshawMultiExp->SetVariableValue(p, profVals[v]);
			else minFarshawMultiExp->SetVariableValue(p,initVal);
		}
		minFarshawMultiExp->FixVariable(varIndex);	
		minFarshawMultiExp->Minimize();
		nllVals.push_back( minFarshawMultiExp->MinValue() );
		minFarshawMultiExp->ReleaseVariable(varIndex);	
		
		if( minFarshawMultiExp->MinValue() < minNll ) {
			minNll = minFarshawMultiExp->MinValue();
			minIndex = v;
		}
		
		const double *fitResult = minFarshawMultiExp->X();
		cout << "Profile NLL Res Over Variable: " << fitters[0]->pdfs[varIndex]->pdfName 
					<< " Value: " << profVals[v] << " Profile NLL At Val: " << nllVals[v] << endl;
		for(unsigned int p=0; p<nType; p++) { // set fit results
			cout << fitters[0]->pdfs[p]->pdfName << " fit To: " << fitResult[p] << endl; 
		}
	}

	vector<double> twoDelNll;
	for(unsigned int v = 0; v < nVal; v++) {
		twoDelNll.push_back( 2*( nllVals[v] - minNll ) );
	}
	TGraph* gTwoDelNll = new TGraph(nVal,profVals.data(),twoDelNll.data());
	gTwoDelNll->SetName("Profile Likelihood");
	gTwoDelNll->SetTitle("Profile Likelihood");
	gTwoDelNll->GetXaxis()->SetTitle(fitters[0]->pdfs[varIndex]->pdfName);
	gTwoDelNll->GetYaxis()->SetTitle("2*DeltaNLL (Profiled)");
	profileNllPlots.push_back(gTwoDelNll);
	profileNllCtr++;
	
	return minNll;
}

void FitToolMultiExp::saveProfileGraphs(TFile* outFile) {
	outFile->cd();
	// Loop over and save
	for(unsigned int ps = 0 ; ps < profileNllPlots.size(); ps++) {
		profileNllPlots[ps]->Write();	
	}
}


// print info about last fit result
void FitToolMultiExp::printLastFit() {
	return; // will write this later
} // end of printLastFit

// prepare fit histograms to be filled
void FitToolMultiExp::prepFitHists() {
	
	// Initialize histograms
	
	hFitLh = new TH1F("hFitLh_"+name,"hFitLh",20000,0,20000);
	hSimPoisLh = new TH1F("hSimPoisLh_"+name,"hSimPoisLh",20000,0,20000);
	hSimTrueLh = new TH1F("hSimTrueLh_"+name,"hSimTrueLh",20000,0,20000);
	hLhFitMinPois = new TH1F("hLhFitMinPois_"+name,"hLhFitMinPois",1000,-50,50);
	hLhFitMinTrue = new TH1F("hLhFitMinTrue_"+name,"hLhFitMinTrue",1000,-50,50);
	// Loop over histograms
	for(unsigned int h = 0; h < nType; h++) {
		TString hName1 = "h1FitDist_" + name + "_" + fitters[0]->pdfs[h]->pdfName;	
		h1FitDist[h] = new TH1F(hName1,hName1,500,0,5);
		TString hAxis1 = "Fitted " + fitters[0]->pdfs[h]->pdfName + " Events";
		h1FitDist[h]->GetXaxis()->SetTitle(hAxis1); 
		h1FitDist[h]->SetLineWidth(2); 
		h1FitDist[h]->SetLineColor(4); 
		vector<TH2F*> tempVec;
		for(unsigned int h2 = 0; h2 < h; h2++) {
			TString hName2 = "h1FitDist_" + name + "_X_" + fitters[0]->pdfs[h]->pdfName + "_Y_" + 
					fitters[0]->pdfs[h2]->pdfName; 
			TString hAxis2 = "Fitted " + fitters[0]->pdfs[h2]->pdfName + " Events";
			TH2F* tempH = new TH2F(hName2,hName2,100,0,5,
					100,0,5);
			tempH->GetXaxis()->SetTitle(hAxis1); 
			tempH->GetYaxis()->SetTitle(hAxis2); 
			tempH->SetDrawOption("COLZ"); 
			tempVec.push_back(tempH);
		}
		h2FitDist.push_back(tempVec);
	}
	
	histsReady = true; // set hist ready flag

} // end of prepFitHists

// fill fit histograms after a fit, for a toy-MC study where many fits are performed
void FitToolMultiExp::fillFitHists() {
	if(histsReady == false) return; // if histograms not set then nothing to fill

	hFitLh->Fill(fitMinLh);
	hSimPoisLh->Fill(simPoisLh);
	hSimTrueLh->Fill(simTrueLh);
	hLhFitMinPois->Fill(fitMinLh - simPoisLh);
	hLhFitMinTrue->Fill(fitMinLh - simTrueLh);

	for(unsigned int h = 0; h < nType; h++) {
		h1FitDist[h]->Fill(fitPars[h]);
		for(unsigned int h2 = 0; h2 < h; h2++) {
			h2FitDist[h][h2]->Fill(fitPars[h],fitPars[h2]);
		}
	}

} // end of fillFitHists

// save fit hists to a root file
void FitToolMultiExp::saveFitHists(TFile* outFile) {
	if(histsReady == false) return;
	outFile->cd();
	// Loop over and save
	hFitLh->Write();
	hSimPoisLh->Write();
	hSimTrueLh->Write();
	hLhFitMinPois->Write();
	hLhFitMinTrue->Write();	
	for(unsigned int h = 0; h < nType; h++) {
		h1FitDist[h]->Write();
		for(unsigned int h2 = 0; h2 < h; h2++) {
			h2FitDist[h][h2]->Write();
		}
	}
} // end of save fit hists

void FitToolMultiExp::getStackedHists() { // make stacked hists with fit/sim/pois pars
	
	// This is responsibility of each sub fitter
	for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
		fitters[f]->getStackedHists();
	}
		
	// Also save the fit results:
	fitParsVec.push_back(fitPars);	
	simTrueParsVec.push_back(simTruePars);	
	simPoisParsVec.push_back(simPoisPars);	
	fitMinLhVec.push_back(fitMinLh);
	simTrueLhVec.push_back(simTrueLh);
	simPoisLhVec.push_back(simPoisLh);

	stackHistCtr++;
} // end of getStackedHists

void FitToolMultiExp::saveStackedHists(TFile* outFile) { // save stacked hists
	
	// This is responsibility of each sub fitter
	for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
		fitters[f]->saveStackedHists(outFile);
	}
	
}// end of saveStackedHists

void FitToolMultiExp::delStackedHists() { // delete prev stacked hists
		
	// This is responsibility of each sub fitter
	for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
		fitters[f]->delStackedHists();
	}
	
	fitParsVec.clear();
	simTrueParsVec.clear();
	simPoisParsVec.clear();
	fitMinLhVec.clear();
	simTrueLhVec.clear();
	simPoisLhVec.clear();

	stackHistCtr = 0;

} // end of delStackedHists

void FitToolMultiExp::printFitRes(TString outPdfName) {
	
	// This is responsibility of each sub fitter
	for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
		TString fitPdfName = outPdfName + Form("_subExp%i.pdf",f); 
		fitters[f]->printFitRes(fitPdfName);
	}	
}

// ** Member function variation of likelihood/gradient functions

// calculate likelihood given parameters
double FitToolMultiExp::getLikelihood(const double *inPars){
	//cout << nType << endl;
	//cout << inPars[0] << endl;
	//cout << "In getLikelihood!" << endl;
	if(dataWasSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	else { // only if data is set
		double lh=0;
		vector<double> subFitPars; 
		subFitPars.resize(nType);
		for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
			for(unsigned int p=0; p<nType; p++) { // loop over pdfs
				subFitPars[p] = fitters[f]->pdfs[p]->expEvents * inPars[p];
			}
			lh += fitters[f]->getLikelihood( subFitPars.data() );	
		}
		
		// Temporary modification to add restriction, improve later!
		/*
		float p3Exp = 1.0;
		float p4Exp = 1.0;
		float p3Err =  0.11* p3Exp;
		float p4Err =  0.03* p4Exp;
		if( nType > 3 ) lh += 0.5*pow((inPars[3] - p3Exp),2)/pow(p3Err,2);
		if( nType > 4 ) lh += 0.5*pow((inPars[4] - p4Exp),2)/pow(p4Err,2);
		*/
		
		return lh;
	} //end of dataWasSet else
} // end of getLikelihood

// calculate likelihood gradient given parameters and direction
double FitToolMultiExp::getLhGradient(const double *inPars, unsigned int coord) {
	//cout << "In getLhGradient" << endl;
	if(dataWasSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	else { // only if data is set
		double grad=0.0; // init gradient from 1
		vector<double> subFitPars; 
		subFitPars.resize(nType);
		for(unsigned int f=0; f<nFitter; f++) { // for individual fitters
			for(unsigned int p=0; p<nType; p++) { // loop over pdfs
				subFitPars[p] = fitters[f]->pdfs[p]->expEvents * inPars[p];
			}
			grad += fitters[f]->pdfs[coord]->expEvents * 
				fitters[f]->getLhGradient( subFitPars.data(), coord );	
		}

		// Temporary modification to add restriction, improve later!
		/*
		float p3Exp = 1.0;
		float p4Exp = 1.0;
		float p3Err = 0.11 * p3Exp;
		float p4Err = 0.03 * p4Exp;
		if( coord == 3 ) grad += (inPars[3] - p3Exp)/pow(p3Err,2); 
		if( coord == 4 ) grad += (inPars[4] - p4Exp)/pow(p4Err,2); 
		*/

		return grad;	
	} //end of dataWasSet else

} // end of getLhGradient

