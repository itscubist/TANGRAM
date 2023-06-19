/********
*Author: Baran Bodur
*Date: 2022-07-26
*Description: Performs global fitting to multi exp/multi PDF in the presence of systematics
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
#include "experiment.h"
#include "organizer.h"
#include "fitKdePdf.h"
#include "mockData.h"
#include "fitGlobal.h"
#include "vectorUtils.h"

using namespace std;

// Functions go here
	
// Constructor
FitGlobal::FitGlobal(Organizer *orgIn) {
	org = orgIn;
	nDim = org->nInt + org->nSys;

	// Resize vectors
	h1FitDist.resize(nDim);
	h1FitDiffDist.resize(nDim);
	h2FitDiffDist.resize(nDim);
	h2Chi2DiffDist.resize(nDim);
	h2FitVsTrueDist.resize(nDim);
	fitCurValues.resize(nDim);
	sysCurValues.resize(org->nSys);
	intCurFracs.resize(org->nInt);
	fitResults.resize(org->nFit);
	fitResultsPerExp.resize(org->nFit);
	for(unsigned int f=0; f<org->nFit; f++ ) {
		fitResultsPerExp[f].resize(org->nExp);
		for(unsigned int e=0; e<org->nExp; e++ ) {
			fitResultsPerExp[f][e].fitFracs.resize(org->nInt);
			fitResultsPerExp[f][e].fitFracErrs.resize(org->nInt);
			fitResultsPerExp[f][e].fitSys.resize(org->nSys);
			fitResultsPerExp[f][e].fitSysErrs.resize(org->nSys);
			fitResultsPerExp[f][e].simTrueFracs.resize(org->nInt);
			fitResultsPerExp[f][e].simTrueSys.resize(org->nSys);
			fitResultsPerExp[f][e].simPoisFracs.resize(org->nInt);
			fitResultsPerExp[f][e].simPoisSys.resize(org->nSys);
			if(e==0) {
				fitResults[f].fitFracs.resize(org->nInt);
				fitResults[f].fitFracErrs.resize(org->nInt);
				fitResults[f].fitSys.resize(org->nSys);
				fitResults[f].fitSysErrs.resize(org->nSys);
				fitResults[f].simTrueFracs.resize(org->nInt);
				fitResults[f].simTrueSys.resize(org->nSys);
				fitResults[f].simPoisFracs.resize(org->nInt);
				fitResults[f].simPoisSys.resize(org->nSys);	
			}
		}
	}
	
	// Create minimizer object: 
	// Note: Possibly update the card reader to read these options from card later
	minFarshaw = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
	minFarshaw->SetMaxFunctionCalls(5000); // max calls, reduce if it takes too much time 
	// set max function calls increase 3000 -> 5000, 20220908!
	minFarshaw->SetMaxIterations(30000); // not used
	minFarshaw->SetTolerance(0.05); // tolerance for EDM
	minFarshaw->SetStrategy(2); // 0:coarse 1: normal 2: detailed 
	minFarshaw->SetPrintLevel(1); 
	
	
	// Create functor to get value & gradient of the function to be minimized at given point
	getLh = &FitGlobal::getLikelihood;
	getGrad = &FitGlobal::getLhGradient;
	//getLh = &FitGlobal::getLikelihoodWei;
	//getGrad = &FitGlobal::getLhGradientWei;
	
	// Functor to be used in minimization	
	fitFunc = ROOT::Math::GradFunctor(this,getLh,getGrad,nDim); 
	//fitFunc = ROOT::Math::Functor(this,getLh,nDim); 
	cout << "Dimensions: " << fitFunc.NDim() << endl;

	// Bind the fit function to the minimizer
	minFarshaw->SetFunction(fitFunc);
	
	// Bind variables: Interaction fractions
	for(unsigned int p = 0; p < org->nInt; p++) {
		string tempStr = (string) org->intNames[p]; 
		cout << "Adding interaction fraction variable to fitter: " << tempStr << endl;
		minFarshaw->SetLowerLimitedVariable(p, tempStr, org->intInitValues[p],org->intStepSize[p],0.0);
	}
	// Bind variables: Systematics
	for(unsigned int s = 0; s < org->nSys; s++) {
		string tempStr = (string) org->sysNames[s]; 
		cout << "Adding systematic variable to fitter: " << tempStr << endl;
		minFarshaw->SetLimitedVariable(s+org->nInt,tempStr, org->sysInitValues[s],
				org->sysStepSize[s],org->sysSigmas[s]*org->sysLowLimit,org->sysSigmas[s]*org->sysUpLimit);
	}

	dataIsSet = false;
	histsReady = false;
	nProfCtr = 0;
	cout << "Global fitter is initialized!" << endl;

} // end of constructor

void FitGlobal::resetFit() {
	fitFunc = ROOT::Math::GradFunctor(this,getLh,getGrad,nDim); 
	//fitFunc = ROOT::Math::Functor(this,getLh,nDim); 
	// Bind the fit function to the minimizer
	minFarshaw->SetFunction(fitFunc);
	// Bind variables: Interaction fractions
	for(unsigned int p = 0; p < org->nInt; p++) {
			minFarshaw->SetLowerLimitedVariable(p,(string)org->intNames[p],
					org->intInitValues[p],org->intStepSize[p],0.0);
	}
	// Bind variables: Systematics
	for(unsigned int s = 0; s < org->nSys; s++) {
			minFarshaw->SetLimitedVariable(s+org->nInt,(string)org->sysNames[s], org->sysInitValues[s],
					org->sysStepSize[s],org->sysSigmas[s]*org->sysLowLimit,org->sysSigmas[s]*org->sysUpLimit);
		}
}

double FitGlobal::makeFit(bool fixSys, bool fixMainPar) {
	
	if(dataIsSet==false) return 0;

	// Do multiple fits starting from different initial positions if requested
	TRandom3* tempRand = new TRandom3();
	tempRand->SetSeed(775028171);

	if(org->nSubFit<1) org->nSubFit = 1;
	std::vector<FitResult> tempFitResults;
	tempFitResults.resize(org->nSubFit);
	std::vector<double> edms;
	edms.resize(org->nSubFit);
	std::vector<int> nFunctionCalls;
	nFunctionCalls.resize(org->nSubFit);

	if(fixMainPar==true) {
		cout << " This is a fit where true value of 1st parameter: " << org->intNames[0]
			   << " is fixed to its initial value!" << endl;
	}

	for(unsigned int f=0; f<org->nSubFit; f++) {

		cout << " This is the " << f+1 << " subfit of " << org->nSubFit << " subfits!" << endl;
		cout << " Initial values for this subfit are: " << endl;
		if(f==0) { // if f=0 then use the initial values given
			// Make sure initial conditions are right
			for(unsigned int p=0; p<org->nInt; p++) {
				minFarshaw->SetVariableValue(p,org->intInitValues[p]);
				cout << org->intNames[p] << " init value is: " << org->intInitValues[p] << endl;
				if(org->intIsFixed[p]==1) minFarshaw->FixVariable(p);
				if(fixMainPar==true && p==0) minFarshaw->FixVariable(p);
			}
			for(unsigned int s=0; s<org->nSys; s++) {
				minFarshaw->SetVariableValue(s+org->nInt, org->sysInitValues[s]);
				cout << org->sysNames[s] << " init value is: " << org->sysInitValues[s] << endl;
				if(fixSys == true) minFarshaw->FixVariable(s+org->nInt); // fix systematics if asked
				if(org->isFixed[s]==1) minFarshaw->FixVariable(s+org->nInt); // if this systematic is fixed
			}
		}
		else {
			for(unsigned int p=0; p<org->nInt; p++) {
				//double tInitVal = gRandom->Gaus(org->intInitValues[p],org->intInitSigmas[p]);
				double tInitVal = tempRand->Gaus(org->intInitValues[p],org->intInitSigmas[p]);
				if(tInitVal<=0) tInitVal = 0.25;
				cout << org->intNames[p] << " init value is: " << tInitVal << endl;
				minFarshaw->SetVariableValue(p,tInitVal);
				if(org->intIsFixed[p]==1) {
					minFarshaw->SetVariableValue(p, org->intInitValues[p]);
					minFarshaw->FixVariable(p);
					tInitVal = org->intInitValues[p];
				}
				if(fixMainPar==true && p==0) {
					minFarshaw->SetVariableValue(p, org->intInitValues[p]);
					minFarshaw->FixVariable(p);
					tInitVal = org->intInitValues[p];
				}
			}
			for(unsigned int s=0; s<org->nSys; s++) {
				//double tInitVal = gRandom->Gaus(org->sysInitValues[s],org->sysInitSigmas[s]);
				double tInitVal = tempRand->Gaus(org->sysInitValues[s],org->sysInitSigmas[s]);
				minFarshaw->SetVariableValue(s+org->nInt, tInitVal);
				if(org->isFixed[s]==1 || fixSys==true) { // if this one or all systematics are fixed
					minFarshaw->SetVariableValue(s+org->nInt, org->sysInitValues[s]); //set to value in card
					minFarshaw->FixVariable(s+org->nInt); // fix
					tInitVal = org->sysInitValues[s];
				}
				cout << org->sysNames[s] << " init value is: " << tInitVal << endl;
			}	
		} // Initial values for the subfit are determined 
	
		// Minimize and get result
		minFarshaw->Minimize();
		const double *tempFitResult = minFarshaw->X();
		const double *tempFitErr = minFarshaw->Errors();
		nFunctionCalls[f] = minFarshaw->NCalls();
		edms[f] = minFarshaw->Edm();	
		cout << "subfit: " << f+1 << " is fitted to (mock) data: " << mockIndex << endl; 

		// Store the result in temporary holders
		tempFitResults[f].fitFracs.resize(org->nInt);
		tempFitResults[f].fitFracErrs.resize(org->nInt);
		tempFitResults[f].fitSys.resize(org->nSys);
		tempFitResults[f].fitSysErrs.resize(org->nSys);
		for(unsigned int p=0; p<org->nInt; p++) { // set fit results
			tempFitResults[f].fitFracs[p] = tempFitResult[p];
			tempFitResults[f].fitFracErrs[p] = tempFitErr[p];
			cout << org->intNames[p] << " fit To: " << tempFitResults[f].fitFracs[p] << " +/- "
					<< tempFitResults[f].fitFracErrs[p] << endl;
		}
		for(unsigned int s=0; s<org->nSys; s++) {
			tempFitResults[f].fitSys[s] = tempFitResult[s+org->nInt];
			tempFitResults[f].fitSysErrs[s] = tempFitErr[s+org->nInt];
			cout << org->sysNames[s] << " fit To: " << tempFitResults[f].fitSys[s] << " +/- "
					<< tempFitResults[f].fitSysErrs[s] << endl;
		}

		// Get likelihood and print info about this sub fit:
		vector<double> fitLhInVec = 
			combIntAndSys(tempFitResults[f].fitFracs, tempFitResults[f].fitSys);
		tempFitResults[f].fitLh = getLikelihood( fitLhInVec.data() );
		//tempFitResults[f].fitLh = getLikelihoodWei( fitLhInVec.data() );
			
		cout << "Min Likelihood: " << tempFitResults[f].fitLh << endl;
		cout << "# of function calls for min:" << minFarshaw->NCalls() << endl;
		cout << "Effective Distance To Minimum: " << minFarshaw->Edm() << endl;
	
	} // end of first loop over subfits  
	
	// Second loop over subfits to summarize results and select minimum
	unsigned int minIndex = 0;
	double minNll = 100000000;
	double minEdm = 1000;
	for(unsigned int f=0; f<org->nSubFit; f++) {
		cout << "Subfit: " << f << " NLL: " << 	tempFitResults[f].fitLh << " Function Calls: " 
			   << nFunctionCalls[f] << " EDM: " << edms[f] << endl; 
		if( tempFitResults[f].fitLh < minNll ) { // determine the fit with min NLL
			minNll = tempFitResults[f].fitLh;
			minIndex = f;
			minEdm = edms[f];
		}
		else if( tempFitResults[f].fitLh == minNll && edms[f] < minEdm) {	// if same NLL, use min EDM
			minNll = tempFitResults[f].fitLh;
			minIndex = f;
			minEdm = edms[f];
		}
	} // end of 2nd loop over subfits to determine the minimum point
	cout << "Min NLL is: " << minNll << " from subfit: " << minIndex << endl;
	
	// Now fill standard holders based on the minimum:
	if(fixMainPar==false) {
		for(unsigned int p=0; p<org->nInt; p++) { // set fit results
			fitResults[mockIndex].fitFracs[p] = tempFitResults[minIndex].fitFracs[p];
			fitResults[mockIndex].fitFracErrs[p] = tempFitResults[minIndex].fitFracErrs[p];
			cout << org->intNames[p] << " fit To: " << fitResults[mockIndex].fitFracs[p] << " +/- "
					<< fitResults[mockIndex].fitFracErrs[p] << endl;
			for(unsigned int e=0;e<org->nExp;e++) {
				fitResultsPerExp[mockIndex][e].fitFracs[p] = fitResults[mockIndex].fitFracs[p];
				fitResultsPerExp[mockIndex][e].fitFracErrs[p] = fitResults[mockIndex].fitFracErrs[p];
			}
		}
		for(unsigned int s=0; s<org->nSys; s++) {
			fitResults[mockIndex].fitSys[s] = tempFitResults[minIndex].fitSys[s];
			fitResults[mockIndex].fitSysErrs[s] = tempFitResults[minIndex].fitSysErrs[s];
			cout << org->sysNames[s] << " fit To: " << fitResults[mockIndex].fitSys[s] << " +/- "
					<< fitResults[mockIndex].fitSysErrs[s] << endl;
			for(unsigned int e=0;e<org->nExp;e++) {
				fitResultsPerExp[mockIndex][e].fitSys[s] = fitResults[mockIndex].fitSys[s];
				fitResultsPerExp[mockIndex][e].fitSysErrs[s] = fitResults[mockIndex].fitSysErrs[s];
			}
		}

		fitResults[mockIndex].fitLh = tempFitResults[minIndex].fitLh;
		fitResults[mockIndex].fitMainFixedLh = tempFitResults[minIndex].fitLh;
		for(unsigned int e=0;e<org->nExp;e++) {
			fitResultsPerExp[mockIndex][e].fitLh = fitResults[mockIndex].fitLh;
			fitResultsPerExp[mockIndex][e].fitMainFixedLh = fitResults[mockIndex].fitLh;
		}
	}
	// Only fill fitMainFixedLh, if adoing a specific fit with first interaction parameter fixed
	else if(fixMainPar==true) {
		fitResults[mockIndex].fitMainFixedLh = tempFitResults[minIndex].fitLh;	
	}
	
	delete tempRand;
	return fitResults[mockIndex].fitLh;




	// OLD VERSION COMMENTED OUT DURING TESTING
/*
	// Make sure initial conditions are right
	for(unsigned int p=0; p<org->nInt; p++) minFarshaw->SetVariableValue(p, org->intInitValues[p]);
	for(unsigned int s=0; s<org->nSys; s++) {
		minFarshaw->SetVariableValue(s+org->nInt, org->sysInitValues[s]);
		if(fixSys == true) minFarshaw->FixVariable(s+org->nInt); // fix systematics if asked
	}
	
	// Minimize and get result
	minFarshaw->Minimize();
	const double *fitResult = minFarshaw->X();
	const double *fitErr = minFarshaw->Errors();

	cout << "Fitted to (mock) data: " << mockIndex << endl; 
	
	for(unsigned int p=0; p<org->nInt; p++) { // set fit results
		fitResults[mockIndex].fitFracs[p] = fitResult[p];
		fitResults[mockIndex].fitFracErrs[p] = fitErr[p];
		cout << org->intNames[p] << " fit To: " << fitResults[mockIndex].fitFracs[p] << " +/- "
			   << fitResults[mockIndex].fitFracErrs[p] << endl;
		for(unsigned int e=0;e<org->nExp;e++) {
			fitResultsPerExp[mockIndex][e].fitFracs[p] = fitResults[mockIndex].fitFracs[p];
			fitResultsPerExp[mockIndex][e].fitFracErrs[p] = fitResults[mockIndex].fitFracErrs[p];
		}
	}
	for(unsigned int s=0; s<org->nSys; s++) {
		fitResults[mockIndex].fitSys[s] = fitResult[s+org->nInt];
		fitResults[mockIndex].fitSysErrs[s] = fitErr[s+org->nInt];
		cout << org->sysNames[s] << " fit To: " << fitResults[mockIndex].fitSys[s] << " +/- "
			   << fitResults[mockIndex].fitSysErrs[s] << endl;
		for(unsigned int e=0;e<org->nExp;e++) {
			fitResultsPerExp[mockIndex][e].fitSys[s] = fitResults[mockIndex].fitSys[s];
			fitResultsPerExp[mockIndex][e].fitSysErrs[s] = fitResults[mockIndex].fitSysErrs[s];
		}
	}

	// Get likelihood:
	vector<double> fitLhInVec = 
		combIntAndSys(fitResults[mockIndex].fitFracs, fitResults[mockIndex].fitSys);
	fitResults[mockIndex].fitLh = getLikelihood( fitLhInVec.data() );
	for(unsigned int e=0;e<org->nExp;e++) {
		fitResultsPerExp[mockIndex][e].fitLh = fitResults[mockIndex].fitLh;
	}
		
	cout << "Min Likelihood: " << fitResults[mockIndex].fitLh << endl;
	//cout << "# of function calls for min:" << minFarshaw->NCalls() << endl;
	//cout << "Effective Distance To Minimum: " << minFarshaw->Edm() << endl;

	return fitResults[mockIndex].fitLh;
	*/
}
	
double FitGlobal::makeProfileFit(unsigned int varIndex, vector<double> profVals, bool fixSys) {
	
	// Do multiple fits starting from different initial positions if requested
	TRandom3* tempRand = new TRandom3();
	tempRand->SetSeed(775028172);
	if(org->nSubFit<1) org->nSubFit = 1;

	unsigned int nVal = profVals.size();
	cout << "Starting profiling over var with index: " << varIndex << " Steps: " << nVal << endl;
	vector<double> nllVals;
	double minNll = 1e10;
	unsigned int minIndex = 0;
	for(unsigned int v = 0; v < nVal; v++){ // loop over profile values
		
		cout << "Profiling, current value: " << profVals[v] << endl;

		// For each value make multiple subfits (if subfit >1)
		std::vector<FitResult> tempFitResults;
		tempFitResults.resize(org->nSubFit);
		std::vector<double> edms;
		edms.resize(org->nSubFit);
		std::vector<int> nFunctionCalls;
		nFunctionCalls.resize(org->nSubFit);
		
		for(unsigned int f=0; f<org->nSubFit; f++) { // loop over subfit for each profile value

			cout << " This is the " << f+1 << " subfit of " << org->nSubFit << " subfits!" << endl;
			cout << " Initial values for this subfit are: " << endl;
			
			if(f==0) { // if f=0 then use the initial values given
				// Make sure initial conditions are right
				for(unsigned int p=0; p<org->nInt; p++) {
					if( p == varIndex ) {
						minFarshaw->SetVariableValue(p, profVals[v]);
						minFarshaw->FixVariable(p);
						cout << org->intNames[p] << " prof value is: " << profVals[v] << endl;
					}
					else{
						minFarshaw->SetVariableValue(p,org->intInitValues[p]);
						cout << org->intNames[p] << " init value is: " << org->intInitValues[p] << endl;
						if(org->intIsFixed[p]==1) minFarshaw->FixVariable(p);
					}
				}
				for(unsigned int s=0; s<org->nSys; s++) {
					if( s+org->nInt == varIndex ) {
						minFarshaw->SetVariableValue(s+org->nInt, profVals[v]);
						minFarshaw->FixVariable(s+org->nInt);
						cout << org->sysNames[s] << " prof value is: " << profVals[v] << endl;
					}
					else {
						minFarshaw->SetVariableValue(s+org->nInt, org->sysInitValues[s]);
						cout << org->sysNames[s] << " init value is: " << org->sysInitValues[s] << endl;
						if(fixSys == true) minFarshaw->FixVariable(s+org->nInt); // fix systematics if asked
						if(org->isFixed[s]==1) minFarshaw->FixVariable(s+org->nInt); // if this sys is fixed
					}
				}
			} //(end of if)
			else {
				for(unsigned int p=0; p<org->nInt; p++) {
					if( p == varIndex ) {
						minFarshaw->SetVariableValue(p, profVals[v]);
						minFarshaw->FixVariable(p);
						cout << org->intNames[p] << " prof value is: " << profVals[v] << endl;
					}
					else {
						//double tInitVal = gRandom->Gaus(org->intInitValues[p],org->intInitSigmas[p]);
						double tInitVal = tempRand->Gaus(org->intInitValues[p],org->intInitSigmas[p]);
						if(tInitVal<=0) tInitVal = 0.25;
						cout << org->intNames[p] << " init value is: " << tInitVal << endl;
						minFarshaw->SetVariableValue(p,tInitVal);
						if(org->intIsFixed[p]==1) {
							minFarshaw->SetVariableValue(p, org->intInitValues[p]);
							minFarshaw->FixVariable(p);
							tInitVal = org->intInitValues[p];
						}
					}
				}
				for(unsigned int s=0; s<org->nSys; s++) {
					if( s+org->nInt == varIndex ) {
						minFarshaw->SetVariableValue(s+org->nInt, profVals[v]);
						minFarshaw->FixVariable(s+org->nInt);
						cout << org->sysNames[s] << " prof value is: " << profVals[v] << endl;
					}
					else {
						//double tInitVal = gRandom->Gaus(org->sysInitValues[s],org->sysInitSigmas[s]);
						double tInitVal = tempRand->Gaus(org->sysInitValues[s],org->sysInitSigmas[s]);
						minFarshaw->SetVariableValue(s+org->nInt, tInitVal);
						if(org->isFixed[s]==1 || fixSys==true) { // if this one or all systematics are fixed
							minFarshaw->SetVariableValue(s+org->nInt, org->sysInitValues[s]); // to val. in card
							minFarshaw->FixVariable(s+org->nInt); // fix
							tInitVal = org->sysInitValues[s];
						}
						cout << org->sysNames[s] << " init value is: " << tInitVal << endl;
					}
				}	
			} // Initial values for the subfit are determined (end of else) 
		
			// Minimize all subfits and choose the minimum:
			minFarshaw->Minimize();
			const double *tempFitResult = minFarshaw->X();
			const double *tempFitErr = minFarshaw->Errors();
			nFunctionCalls[f] = minFarshaw->NCalls();
			edms[f] = minFarshaw->Edm();	
			cout << "subfit: " << f+1 << " is fitted to (mock) data: " << mockIndex << endl; 

			// Store the result in temporary holders
			tempFitResults[f].fitFracs.resize(org->nInt);
			tempFitResults[f].fitFracErrs.resize(org->nInt);
			tempFitResults[f].fitSys.resize(org->nSys);
			tempFitResults[f].fitSysErrs.resize(org->nSys);
			for(unsigned int p=0; p<org->nInt; p++) { // set fit results
				tempFitResults[f].fitFracs[p] = tempFitResult[p];
				tempFitResults[f].fitFracErrs[p] = tempFitErr[p];
				cout << org->intNames[p] << " fit To: " << tempFitResults[f].fitFracs[p] << " +/- "
						<< tempFitResults[f].fitFracErrs[p] << endl;
			}
			for(unsigned int s=0; s<org->nSys; s++) {
				tempFitResults[f].fitSys[s] = tempFitResult[s+org->nInt];
				tempFitResults[f].fitSysErrs[s] = tempFitErr[s+org->nInt];
				cout << org->sysNames[s] << " fit To: " << tempFitResults[f].fitSys[s] << " +/- "
						<< tempFitResults[f].fitSysErrs[s] << endl;
			}

			// Get likelihood and print info about this sub fit:
			vector<double> fitLhInVec = 
				combIntAndSys(tempFitResults[f].fitFracs, tempFitResults[f].fitSys);
			tempFitResults[f].fitLh = getLikelihood( fitLhInVec.data() );
			//tempFitResults[f].fitLh = getLikelihoodWei( fitLhInVec.data() );
				
			cout.precision(9);
			cout << "Min Likelihood: " << tempFitResults[f].fitLh << endl;
			cout << "# of function calls for min:" << minFarshaw->NCalls() << endl;
			cout << "Effective Distance To Minimum: " << minFarshaw->Edm() << endl;
			cout.precision(6);
		
		} // end of first loop over subfits  
		
		// Second loop over subfits to summarize results and select minimum
		unsigned int minIndexS = 0;
		double minNllS = 1e10;
		double minEdmS = 1e10;
		for(unsigned int f=0; f<org->nSubFit; f++) {
			cout.precision(9);
			cout << "Subfit: " << f << " NLL: " << 	tempFitResults[f].fitLh << " Function Calls: " 
					<< nFunctionCalls[f] << " EDM: " << edms[f] << endl; 
			cout.precision(6);
			if( tempFitResults[f].fitLh < minNllS ) { // determine the fit with min NLL
				minNllS = tempFitResults[f].fitLh;
				minIndexS = f;
				minEdmS = edms[f];
			}
			else if( tempFitResults[f].fitLh == minNllS && edms[f] < minEdmS) {	//same NLL, use min EDM
				minNllS = tempFitResults[f].fitLh;
				minIndexS = f;
				minEdmS = edms[f];
			}
		} // end of 2nd loop over subfits to determine the minimum point
		cout.precision(9);
		cout << "Min NLL is: " << minNllS << " from subfit: " << minIndexS << endl;
		cout.precision(6);
		
		// Now min NLL to profiled NLLs vector
		nllVals.push_back( minNllS );
		minFarshaw->ReleaseVariable(varIndex);	
		
		if( minFarshaw->MinValue() < minNll ) {
			minNll = minFarshaw->MinValue();
			minIndex = v;
		}
		
		cout.precision(9);
		cout << "Profile NLL Res Over Variable: " << org->intNames[varIndex] 
					<< " Value: " << profVals[v] << " Profile NLL At Val: " << nllVals[v] << endl;
		cout.precision(6);
		for(unsigned int p=0; p<org->nInt; p++) { // set fit results
			cout << org->intNames[p] << " fit To: " << tempFitResults[minIndexS].fitFracs[p] << endl; 
		}
		for(unsigned int s=0; s<org->nSys; s++) { // set fit results
			cout << "Systematic: " << org->sysNames[s] << " fit To: " 
				   << tempFitResults[minIndexS].fitSys[s] << endl; 
		}

		/* // Old profile init	
		for(unsigned int p=0; p<nDim; p++) {
			if( p == varIndex ) minFarshaw->SetVariableValue(p, profVals[v]);
			else if( p < org->nInt ) minFarshaw->SetVariableValue(p,org->intInitValues[p]);
			else if( p >= org->nInt ) {
			  minFarshaw->SetVariableValue(p,org->sysInitValues[p-org->nInt]);
				if(fixSys==true) {
			  	minFarshaw->SetVariableValue(p,org->sysInitValues[p-org->nInt]);
					minFarshaw->FixVariable(p);	
				}
			}
		} // end of loop over dimension
	
		minFarshaw->FixVariable(varIndex);	
		minFarshaw->Minimize();
		nllVals.push_back( minFarshaw->MinValue() );
		minFarshaw->ReleaseVariable(varIndex);	
		

		if( minFarshaw->MinValue() < minNll ) {
			minNll = minFarshaw->MinValue();
			minIndex = v;
		}
		
		const double *fitResult = minFarshaw->X();
		cout << "Profile NLL Res Over Variable: " << org->intNames[varIndex] 
					<< " Value: " << profVals[v] << " Profile NLL At Val: " << nllVals[v] << endl;
		for(unsigned int p=0; p<org->nInt; p++) { // set fit results
			cout << org->intNames[p] << " fit To: " << fitResult[p] << endl; 
		}
		for(unsigned int s=0; s<org->nSys; s++) { // set fit results
			cout << "Systematic: " << org->sysNames[s] << " fit To: " << fitResult[s+org->nInt] << endl; 
		}
		*/

	} // end of loop over prof vals
	
	vector<double> twoDelNll;
	for(unsigned int v = 0; v < nVal; v++) {
		twoDelNll.push_back( 2*( nllVals[v] - minNll ) );
	}

	TGraph* gTwoDelNll = new TGraph(nVal,profVals.data(),twoDelNll.data());
	TString gName = Form("ProfileLikelihood_Fit%d",nProfCtr);
	if( fixSys==true )  gName = gName + "_fixedSys";
	if( fixSys==false )  gName = gName + "_withSys";
	gTwoDelNll->SetName(gName);
	gTwoDelNll->SetTitle("Profile Likelihood");
	gTwoDelNll->GetXaxis()->SetTitle(org->intNames[varIndex] + " Fraction");
	gTwoDelNll->GetYaxis()->SetTitle("2*DeltaNLL (Profiled)");
	profGraphs.push_back(gTwoDelNll);

  // Increase prof counter and make sure to release all vars
	nProfCtr++;
	for(unsigned int p=0; p<nDim; p++) {
		minFarshaw->ReleaseVariable(p);
	}
	return minNll;

} // end of profile likehood function

// Loop over and save graphs
void FitGlobal::saveProfGraphs() {
	org->outFile->cd();
	for(unsigned int p=0;p<nProfCtr;p++) {
		profGraphs[p]->Write();
	}

}

// Set data index to fit:
void FitGlobal::setMockIndex(unsigned int m) {
	dataIsSet = true;
	for(unsigned int e=0;e<org->nExp;e++) {
		if( !(m < org->exps[e]->mocks.size()) ) dataIsSet = false;
	}
	if(dataIsSet==true) {
		mockIndex = m;
		setSimPars();
		isMockData = true;
	}
	else {
		cout << "!!! *** ERROR: cannot set data, data with the requested index " << m 
			   << "does not exist for some experiments!" << endl;
	}
}

// Set simulated pars from input pdfs and mock data
void FitGlobal::setSimPars() { 
	if(dataIsSet==false) {
		cout << "No mock data was set, cannot determine simulated parameters!" << endl;
		return;
	}

	vector<double> expSum;
	expSum.resize(org->nInt);
	for(unsigned int e=0;e<org->nExp;e++) {
		for(unsigned int p=0;p<org->nInt;p++) {
			if(e==0) { // To only do once in the beginning
				fitResults[mockIndex].simTrueFracs[p] = 0;	
				fitResults[mockIndex].simPoisFracs[p] = 0;	
				expSum[p] = 0;
			}
			fitResults[mockIndex].simTrueFracs[p] += org->exps[e]->mocks[mockIndex]->simTrueEvents[p];
			fitResults[mockIndex].simPoisFracs[p] += org->exps[e]->mocks[mockIndex]->simPoisEvents[p];
			expSum[p] += org->exps[e]->mocks[mockIndex]->simExpEvents[p];


			fitResultsPerExp[mockIndex][e].simTrueFracs[p] = 
				org->exps[e]->mocks[mockIndex]->simTrueEvents[p] / 
				org->exps[e]->mocks[mockIndex]->simExpEvents[p];
			fitResultsPerExp[mockIndex][e].simPoisFracs[p] = 
				org->exps[e]->mocks[mockIndex]->simPoisEvents[p] / 
				org->exps[e]->mocks[mockIndex]->simExpEvents[p];
		}
		for(unsigned int s=0;s<org->nSys;s++) {
			if(e==0) { // To only do once in the beginning: systematics should be same for all exp
				fitResults[mockIndex].simTrueSys[s] = org->exps[e]->mocks[mockIndex]->simSystematics[s];
				fitResults[mockIndex].simPoisSys[s] = org->exps[e]->mocks[mockIndex]->simSystematics[s];
			}
			fitResultsPerExp[mockIndex][e].simTrueSys[s] = 
				org->exps[e]->mocks[mockIndex]->simSystematics[s];
			fitResultsPerExp[mockIndex][e].simPoisSys[s] = 
				org->exps[e]->mocks[mockIndex]->simSystematics[s];
		}
		// Assign 0 for these as there is no likelihood function per signle experiment as of yet
		fitResultsPerExp[mockIndex][e].simTrueLh = 0;
		fitResultsPerExp[mockIndex][e].simPoisLh = 0;	
	}
	
	// Divide to get fraction:
	for(unsigned int p=0;p<org->nInt;p++) {
		fitResults[mockIndex].simTrueFracs[p] /= expSum[p];
		fitResults[mockIndex].simPoisFracs[p] /= expSum[p];
	}
	
	// Get likelihood:
	vector<double> simTrueLhInVec = 
		combIntAndSys(fitResults[mockIndex].simTrueFracs, fitResults[mockIndex].simTrueSys);
	fitResults[mockIndex].simTrueLh = getLikelihood( simTrueLhInVec.data() );
	//fitResults[mockIndex].simTrueLh = getLikelihoodWei( simTrueLhInVec.data() );
	
	vector<double> simPoisLhInVec = 
		combIntAndSys(fitResults[mockIndex].simPoisFracs, fitResults[mockIndex].simPoisSys);
	fitResults[mockIndex].simPoisLh = getLikelihood( simPoisLhInVec.data() );
	//fitResults[mockIndex].simPoisLh = getLikelihoodWei( simPoisLhInVec.data() );

}

// Prepare fit histograms
void FitGlobal::prepFitHists() {

	// Initialize histograms	
	hFitLh = new TH1F("hFitLh_all","hFitLh",30000,0,30000);
	hSimPoisLh = new TH1F("hSimPoisLh_all","hSimPoisLh",30000,0,30000);
	hSimTrueLh = new TH1F("hSimTrueLh_all","hSimTrueLh",30000,0,30000);
	hFitFixedLh = new TH1F("hFixFixedLh_all","hFitFixedLh",30000,0,30000);
	hLhFitMinPois = new TH1F("hLhFitMinPois_all","hLhFitMinPois",1000,-50,50);
	hLhFitMinTrue = new TH1F("hLhFitMinTrue_all","hLhFitMinTrue",1000,-50,50);
	hLhFitMinFitFixed = new TH1F("hLhFitMinFitFixed_all","hLhFitMinFitFixed",1000,-50,50);
	
	hFitChi2 = new TH1F("hFitChi2","hFitChi2",1000,0,1000);
	hSimTrueChi2 = new TH1F("hSimTrueChi2","hSimTrueChi2",1000,0,1000);
	hChi2FitMinTrue = new TH1F("hChi2FitMinTrue","hChi2FitMinTrue",1000,-50,50);
	
	hNDataEvents = new TH1F("hNDataEvents","hNDataEvents",10000,0,10000);
	hFitEvisIntegral = new TH1F("hFitEvisIntegral","hFitEvisIntegral",10000,0,10000);
	hSimEvisIntegral = new TH1F("hSimEvisIntegral","hSimEvisIntegral",10000,0,10000);
	hFitEvisBeginDiff = new TH1F("hFitEvisBeginDiff","hFitEvisBeginDiff",2000,-1000,1000);
	hSimEvisBeginDiff = new TH1F("hSimEvisBeginDiff","hSimEvisBeginDiff",2000,-1000,1000);

	// Loop over histograms
	for(unsigned int h = 0; h < nDim; h++) {
		TString hName1, hName2, hAxis1, hAxis2;
		unsigned int xBins, xBins2, yBins2;
		float xMin, xMax, yMin, yMax;
		TString hNameDiff1, hAxisDiff1;
		unsigned int xDiffBins;
		float xDiffMin, xDiffMax;
		TString hNameFitMinPois2, hNameChi2FitMinTrue2;
		TString hNameFitVsTrue, hAxisTrue, hAxisFit;
		if( h < org->nInt ) {
			// For standard plots
			hName1 = "h1FitDist_" + org->intNames[h];	
			xBins = 500;
			xBins2 = 50;
			xMin = 0;
			xMax = 5;
			hAxis1 = "Fitted " + org->intNames[h] + " Fraction";
			// For fitted - true plots
			hNameDiff1 = "h1FitDiffDist_" + org->intNames[h];	
			xDiffBins = 1000;
			xDiffMin = -5; 
			xDiffMax = 5; 
			hAxisDiff1 = "Fitted - True " + org->intNames[h] + " Fraction";
			// For 2D plots vs fitMinPois
			hNameFitMinPois2 = "hLhFitMinPois_Vs_Fitted-True_" + org->intNames[h] + "_Fraction";
			hNameChi2FitMinTrue2 = "hChi2FitMinTrue_Vs_Fitted-True_" + org->intNames[h] + "_Fraction";
			hNameFitVsTrue = "hFitVsTrueValue_" + org->intNames[h];
			hAxisFit = "Fitted " + org->intNames[h] + " Fraction";
			hAxisTrue = "True " + org->intNames[h] + " Fraction";

		}
		else if( h >= org->nInt ) {
			// For standard plots
			hName1 = "h1FitDist_" + org->sysNames[h-org->nInt];	
			xBins = 1000;
			xBins2 = 100;
			xMin = -5;
			xMax = 5;
			hAxis1 = "Fitted " + org->sysNames[h-org->nInt] + " Systematic (Sigma)";	
			// For fitted - true plots
			hNameDiff1 = "h1FitDiffDist_" + org->sysNames[h-org->nInt];	
			xDiffBins = 1000;
			xDiffMin = -5; 
			xDiffMax = 5; 
			hAxisDiff1 = "(Fitted - True)/Sigma " + org->sysNames[h-org->nInt] + " Systematic (Sigma)";	
			// For 2D plots vs fitMinPois
			hNameFitMinPois2 = "hLhFitMinPois_Vs_Fitted-True_" + org->sysNames[h-org->nInt] + 
				" Systematic (Sigma)";
			hNameChi2FitMinTrue2 = "hChi2FitMinTrue_Vs_Fitted-True_" + org->sysNames[h-org->nInt] + 
				" Systematic (Sigma)";
			hNameFitVsTrue = "hFitVsTrueValue_" + org->sysNames[h-org->nInt];
			hAxisFit = "Fitted " + org->sysNames[h-org->nInt] + " Systematic (Sigma)"; 
			hAxisTrue = "True " + org->sysNames[h-org->nInt] + " Systematic (Sigma)"; 
		}
		h1FitDist[h] = new TH1F(hName1,hName1,xBins,xMin,xMax);
		h1FitDist[h]->GetXaxis()->SetTitle(hAxis1); 
		h1FitDist[h]->SetLineWidth(2); 
		h1FitDist[h]->SetLineColor(4); 
		
		h1FitDiffDist[h] = new TH1F(hNameDiff1,hNameDiff1,xDiffBins,xDiffMin,xDiffMax);
		h1FitDiffDist[h]->GetXaxis()->SetTitle(hAxisDiff1); 
		h1FitDiffDist[h]->SetLineWidth(2); 
		h1FitDiffDist[h]->SetLineColor(4); 

		h2FitDiffDist[h] = new TH2F(hNameFitMinPois2,hNameFitMinPois2,
				xDiffBins,xDiffMin,xDiffMax,1000,-50,50);
		h2FitDiffDist[h]->GetXaxis()->SetTitle(hAxisDiff1); 
		h2FitDiffDist[h]->GetYaxis()->SetTitle("Fitted-True NLL"); 
		
		h2FitVsTrueDist[h] = new TH2F(hNameFitVsTrue,hNameFitVsTrue,
				xBins2,xMin,xMax,xBins2,xMin,xMax);
		h2FitVsTrueDist[h]->GetXaxis()->SetTitle(hAxisFit); 
		h2FitVsTrueDist[h]->GetYaxis()->SetTitle(hAxisTrue); 
	
		h2Chi2DiffDist[h] = new TH2F(hNameChi2FitMinTrue2,hNameChi2FitMinTrue2,
				xDiffBins,xDiffMin,xDiffMax,1000,-50,50);
		h2Chi2DiffDist[h]->GetXaxis()->SetTitle(hAxisDiff1); 
		h2Chi2DiffDist[h]->GetYaxis()->SetTitle("Fitted-True Binned Chi2"); 
		
		// 2D histograms
		vector<TH2F*> tempVec;
		vector<TH2F*> tempVec2;
		for(unsigned int h2 = 0; h2 < h; h2++) {
			if( h2 < org->nInt ) {
				hName2 = hName1 + "_vs_" + org->intNames[h2]; 
				yBins2 = 50;
				yMin = 0;
				yMax = 5;
				hAxis2 = "Fitted " + org->intNames[h2] + " Fraction";
			}
			else if( h2 >= org->nInt ) {
				hName2 = hName1 + "_vs_" + org->sysNames[h2-org->nInt]; 
				yBins2 = 100;
				yMin = -5;
				yMax = 5;
				hAxis2 = "Fitted " + org->sysNames[h2-org->nInt] + " Systematic (Sigma)";
			}
			TH2F* tempH = new TH2F(hName2,hName2,xBins2,xMin,xMax,yBins2,yMin,yMax);
			tempH->GetXaxis()->SetTitle(hAxis1); 
			tempH->GetYaxis()->SetTitle(hAxis2); 
			tempH->SetDrawOption("COLZ"); 
			tempVec.push_back(tempH);
		} // end of 2nd loop over fit var for 2d
		h2FitDist.push_back(tempVec);
	} // end of 1st loop over fit var
	histsReady = true;
} // end of prepare fit histograms function

// Fill fit histograms
void FitGlobal::fillFitHists() {
	
	if(histsReady == false) return;

	hFitLh->Fill(fitResults[mockIndex].fitLh);
	hFitFixedLh->Fill(fitResults[mockIndex].fitMainFixedLh);
	hSimPoisLh->Fill(fitResults[mockIndex].simPoisLh);
	hSimTrueLh->Fill(fitResults[mockIndex].simTrueLh);
	hLhFitMinPois->Fill(fitResults[mockIndex].fitLh - fitResults[mockIndex].simPoisLh);
	hLhFitMinTrue->Fill(fitResults[mockIndex].fitLh - fitResults[mockIndex].simTrueLh);
	hLhFitMinFitFixed->Fill(fitResults[mockIndex].fitLh - fitResults[mockIndex].fitMainFixedLh);

	hFitChi2->Fill(fitResults[mockIndex].fitChi2);
	hSimTrueChi2->Fill(fitResults[mockIndex].simTrueChi2);
	hChi2FitMinTrue->Fill(fitResults[mockIndex].fitChi2 - fitResults[mockIndex].simPoisChi2);
	hNDataEvents->Fill(fitResults[mockIndex].nDataEvents);
	hFitEvisIntegral->Fill(fitResults[mockIndex].fitEvisIntegral);
	hSimEvisIntegral->Fill(fitResults[mockIndex].simEvisIntegral);
	hFitEvisBeginDiff->Fill(fitResults[mockIndex].fitEvisBeginDiff);
	hSimEvisBeginDiff->Fill(fitResults[mockIndex].simEvisBeginDiff);

	for(unsigned int h = 0; h < nDim; h++) {
		double fillX, fillY;
		double fillDiffX;
		if( h < org->nInt ) {
			fillX = fitResults[mockIndex].fitFracs[h];
			fillDiffX = fillX - fitResults[mockIndex].simTrueFracs[h];
			fillY = fitResults[mockIndex].simTrueFracs[h];
		}
		else {
			fillX = fitResults[mockIndex].fitSys[h-org->nInt]/org->sysSigmas[h-org->nInt];
			fillDiffX = fillX-fitResults[mockIndex].simTrueSys[h-org->nInt]/org->sysSigmas[h-org->nInt]; 
			fillY = fitResults[mockIndex].simTrueSys[h-org->nInt]/org->sysSigmas[h-org->nInt];
		}
		h1FitDist[h]->Fill(fillX);
		h1FitDiffDist[h]->Fill(fillDiffX);
		h2FitDiffDist[h]->Fill(fillDiffX, 
				(fitResults[mockIndex].fitLh - fitResults[mockIndex].simPoisLh) ); 
		h2Chi2DiffDist[h]->Fill(fillDiffX, 
				(fitResults[mockIndex].fitChi2 - fitResults[mockIndex].simTrueChi2) );
		h2FitVsTrueDist[h]->Fill(fillX,fillY);
		for(unsigned int h2 = 0; h2 < h; h2++) {
			if( h2 < org->nInt ) fillY = fitResults[mockIndex].fitFracs[h2];
			else fillY = fitResults[mockIndex].fitSys[h2-org->nInt]/org->sysSigmas[h2-org->nInt];
			h2FitDist[h][h2]->Fill(fillX,fillY);
		}
	}

} // end of fill fit histograms function

// Save fit histograms
void FitGlobal::saveFitHists() {
	if(histsReady == false) return;
	org->outFile->cd();
	// Loop over and save
	hFitLh->Write();
	hFitFixedLh->Write();
	hSimPoisLh->Write();
	hSimTrueLh->Write();
	hLhFitMinPois->Write();
	hLhFitMinFitFixed->Write();	
	// Mostly for testing
	hFitChi2->Write();
	hSimTrueChi2->Write();
	hChi2FitMinTrue->Write();
	hNDataEvents->Write();
	hFitEvisIntegral->Write();
	hSimEvisIntegral->Write();
	hFitEvisBeginDiff->Write();
	hSimEvisBeginDiff->Write();
	for(unsigned int h = 0; h < nDim; h++) {
		h1FitDist[h]->Write();
		h1FitDiffDist[h]->Write();
		h2FitDiffDist[h]->Write();
		h2Chi2DiffDist[h]->Write();
		h2FitVsTrueDist[h]->Write();
		for(unsigned int h2 = 0; h2 < h; h2++) {
			h2FitDist[h][h2]->Write();
		}
	}
	
} // end of save fit histograms

// To Print Fit Result Into Pdf File
void FitGlobal::printFitGraphs(unsigned int fitNo, TString type) {
	
	// Acquire required info based on fitNo and plotting type
	vector< double > intFracs; intFracs.resize(org->nInt);
	vector< double > sysVals; sysVals.resize(org->nSys);
	double lhVal;
	TString nameStart;
	TString textStart;
	if(type == "fit") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].fitFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].fitSys[s];
		lhVal = fitResults[fitNo].fitLh; 
		nameStart = Form("hStack_Fit%d",fitNo);
		textStart = "Fit Min Lh: "; 
	}
	else if(type == "pois") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].simPoisFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].simPoisSys[s];
		lhVal = fitResults[fitNo].simPoisLh; 
		nameStart = Form("hStack_SimPois%d",fitNo);
		textStart = "Sim Pois Lh: "; 	
	}
	else if(type == "sim") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].simTrueFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].simTrueSys[s];
		lhVal = fitResults[fitNo].simTrueLh; 
		nameStart = Form("hStack_SimTrue%d",fitNo);
		textStart = "Sim True Lh: ";	
	}
	else {
		cout << "WARNING! In printFitGraphs, graph type: " << type << " is not an option!" << endl;
		return;
	}
	
	// Calculate columns and rows, and create TCanvas accordingly
	unsigned int nRows = 4;
	unsigned int nCols = 0;
	for(unsigned int e=0;e<org->nExp;e++) {
		if(org->neutronInfo[e]==false) nCols += 1;
		if(org->neutronInfo[e]==true) nCols += 3;
	}

	// Temporaries for looping over
	vector<THStack*> tempStacks;
	unsigned tsInd = 0;
	unsigned int cRow;
	unsigned int cCol;
	double lrMargin = 0.05;

	// First loop to fill first row with eVis histograms
	cRow = 0;
	cCol = 1;
	for(unsigned int e=0;e<org->nExp;e++) {
		unsigned int nMax;
		if(org->neutronInfo[e]==false) nMax = 1;
		if(org->neutronInfo[e]==true) nMax = 3;
		for(unsigned int nn=0; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "eVis";
			tempStacks.push_back( 
					org->exps[e]->makeStackedHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			org->canny->cd( cRow*nCols + cCol );
			org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
			org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
			double maxY = 1.25*org->exps[e]->mocks[fitNo]->evisProjs[tN]->GetMaximum();
			org->exps[e]->mocks[fitNo]->evisProjs[tN]->GetYaxis()->SetRangeUser(0,maxY);
			org->exps[e]->mocks[fitNo]->evisProjs[tN]->Draw();
			tempStacks[tsInd]->Draw("same hist");
			tsInd++;
			org->exps[e]->mocks[fitNo]->evisProjs[tN]->Draw("same");	
			cCol += 1;
		}
	}
	// Second loop to fill 2nd row with gTag histograms
	cRow = 1;
	cCol = 1;
	for(unsigned int e=0;e<org->nExp;e++) {
		if(org->gammaInfo[e]==false) continue;
		unsigned int nMax;
		if(org->neutronInfo[e]==false) nMax = 1;
		if(org->neutronInfo[e]==true) nMax = 3;
		for(unsigned int nn=0; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "gTag";
			tempStacks.push_back( 
					org->exps[e]->makeStackedHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			org->canny->cd( cRow*nCols + cCol );
			org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
			org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
			double maxY = 1.25*org->exps[e]->mocks[fitNo]->gtagProjs[tN][tE]->GetMaximum();
			org->exps[e]->mocks[fitNo]->gtagProjs[tN][tE]->GetYaxis()->SetRangeUser(0,maxY);
			org->exps[e]->mocks[fitNo]->gtagProjs[tN][tE]->Draw();
			tempStacks[tsInd]->Draw("same hist");
			tsInd++;
			org->exps[e]->mocks[fitNo]->gtagProjs[tN][tE]->Draw("same");	
			cCol += 1;
		}
	}
	// Third loop to fill 3rd row with nTag histograms
	cRow = 2;
	cCol = 1;
	for(unsigned int e=0;e<org->nExp;e++) {
		unsigned int nMax;
		if(org->neutronInfo[e]==false) nMax = 1;
		if(org->neutronInfo[e]==true) nMax = 3;
		unsigned int tN = 0;
		unsigned int tE = 0;
		TString hType = "nTag";
		tempStacks.push_back(
				org->exps[e]->makeStackedHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
		org->canny->cd( cRow*nCols + cCol );
		org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
		org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
		double maxY = 1.25*org->exps[e]->mocks[fitNo]->hNeutProj->GetMaximum();
		org->exps[e]->mocks[fitNo]->hNeutProj->GetYaxis()->SetRangeUser(0,maxY);
		org->exps[e]->mocks[fitNo]->hNeutProj->Draw();
		tempStacks[tsInd]->Draw("same hist");
		tsInd++;
		org->exps[e]->mocks[fitNo]->hNeutProj->Draw("same");	
		cCol += nMax;
	}

	// Legend
	cRow = 3;
	cCol = 1;
	TLegend* leggy = new TLegend(0,0.25,1,0.75,"","nbNDC");
	leggy->SetBorderSize(0);
	leggy->SetTextFont(43);
	leggy->SetEntrySeparation(0.2);
	vector<TBox*> intBoxes; intBoxes.resize(org->nInt);
	for(unsigned int i=0;i<org->nInt;i++) {
		intBoxes[i] = new TBox(0.1,0.1,0.2,0.2);
		intBoxes[i]->SetLineColor(org->intColors[i]);
		intBoxes[i]->SetFillColor(org->intColors[i]);
		leggy->AddEntry(intBoxes[i],org->intNames[i]);
		
	}
	TLine* line1 = new TLine(0.1,0.1,0.2,0.2);
	line1->SetLineColor(6);
	leggy->AddEntry(line1,Form("MockData_%d",fitNo));
	org->canny->cd( cRow*nCols + cCol );
	org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
	org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
	leggy->Draw();

	if(nCols>1) {
		cRow = 3;
		cCol = 2;
		double tX = 0.4;
		double tY = 0.9;
		double dY = 0.05;
		double dyCtr = 0;
		TString textStr;
		TText* texty = new TText();
		texty->SetTextFont(8);
		org->canny->cd( cRow*nCols + cCol )->Clear();
		org->canny->cd( cRow*nCols + cCol );
		org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
		org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);

		textStr = textStart + Form("%.2f",lhVal);
		texty->DrawTextNDC(tX,tY-dY*dyCtr,textStr);
		dyCtr +=1;
		for(unsigned int i=0; i < org->nInt; i++) {
			textStr = org->intNames[i] + " : " + Form("%.2f",intFracs[i]);
			texty->DrawTextNDC(tX,tY-dY*dyCtr,textStr);
			dyCtr +=1;
			if(dyCtr>15) { // start writin to new column if current one is full
				cCol++;
				if(cCol<nCols) { // nothing to do if not enough columns
					dyCtr = 0;
					org->canny->cd( cRow*nCols + cCol )->Clear();
					org->canny->cd( cRow*nCols + cCol );
					org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
					org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
				}
			}
		}			
		
		for(unsigned int s=0; s < org->nSys; s++) {
			textStr = org->sysNames[s] + " : " + Form("%.2f",sysVals[s]);
			texty->DrawTextNDC(tX,tY-dY*dyCtr,textStr);
			dyCtr +=1;
			if(dyCtr>15) { // start writing to new column if current one is full
				cCol++;
				if(cCol<nCols) { // nothing to do if not enough columns
					dyCtr = 0;
					org->canny->cd( cRow*nCols + cCol )->Clear();
					org->canny->cd( cRow*nCols + cCol );
					org->canny->cd( cRow*nCols + cCol )->SetLeftMargin(lrMargin);
					org->canny->cd( cRow*nCols + cCol )->SetRightMargin(lrMargin);
				}
			}
		}			
	
	} // end of if more than 1 column

	
	org->canny->Print(org->outPdfFile);
	for(unsigned int v=0;v<tempStacks.size();v++) delete tempStacks[v];
	cout << "Printed stacked plots, fitNo: " << fitNo << " of type: " << type << endl;

} // end of printFitGraphs

// To Get Binned Chi Square Value and Print Fit Result Into Root File
void FitGlobal::getFitChiAndHists(unsigned int fitNo, TString type) {
	
	// Acquire required info based on fitNo and plotting type
	vector< double > intFracs; intFracs.resize(org->nInt);
	vector< double > sysVals; sysVals.resize(org->nSys);
	double lhVal;
	TString nameStart;
	TString textStart;
	if(type == "fit") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].fitFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].fitSys[s];
		nameStart = Form("hTotal_Fit%d",fitNo);
	}
	else if(type == "pois") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].simPoisFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].simPoisSys[s];
		nameStart = Form("hTotal_SimPois%d",fitNo);
	}
	else if(type == "sim") {
		for(unsigned int i=0;i<org->nInt;i++) intFracs[i] = fitResults[fitNo].simTrueFracs[i];	
		for(unsigned int s=0;s<org->nSys;s++) sysVals[s] = fitResults[fitNo].simTrueSys[s];
		nameStart = Form("hTotal_SimTrue%d",fitNo);
	}
	else {
		cout << "WARNING! In printFitGraphs, graph type: " << type << " is not an option!" << endl;
		return;
	}
	
	// Temporaries for looping over
	vector<TH1D*> tempHists;
	unsigned int thInd = 0;
	double chi2Val = 0;
	double ndf = 0;
	unsigned int nBins = 0;

	double nDataEvents = 0;
	double evisIntegral = 0;
	double evisBeginDiff = 0;
	
	// First loop for eVis histograms
	for(unsigned int e=0;e<org->nExp;e++) {
		unsigned int nMax;
		unsigned int nMin;
		if(org->neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(org->neutronInfo[e]==true) {nMin=1; nMax = 3;}
		for(unsigned int nn=nMin; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "eVis";
			tempHists.push_back( 
					org->exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			nBins = tempHists[thInd]->GetNbinsX(); 
			for(unsigned int b=0;b<nBins;b++) {
				double theoVal = tempHists[thInd]->GetBinContent(b+1);
				double poisVal = org->exps[e]->mocks[fitNo]->evisProjs[tN]->GetBinContent(b+1);
				double poisErrSqr = poisVal;
				if(poisVal>0 && theoVal>0) {
					chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
				}	
				if(b<5) {
					evisBeginDiff += (poisVal - theoVal); 
				}
			}
			nDataEvents += org->exps[e]->mocks[fitNo]->evisProjs[tN]->Integral();
			evisIntegral += tempHists[thInd]->Integral(); 	
			cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
			thInd++;
		}
	}
	// Second loop to fill 2nd row with gTag histograms
	for(unsigned int e=0;e<org->nExp;e++) {
		if(org->gammaInfo[e]==false) continue;
		unsigned int nMax;
		unsigned int nMin;
		if(org->neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(org->neutronInfo[e]==true) {nMin=1; nMax = 3;}
		for(unsigned int nn=nMin; nn<nMax; nn++) {
			unsigned int tN = nn;
			unsigned int tE = 0;
			TString hType = "gTag";
			tempHists.push_back( 
					org->exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
			nBins = tempHists[thInd]->GetNbinsX(); 
			for(unsigned int b=0;b<nBins;b++) {
				double theoVal = tempHists[thInd]->GetBinContent(b+1);
				double poisVal = org->exps[e]->mocks[fitNo]->gtagProjs[tN][tE]->GetBinContent(b+1);
				double poisErrSqr = poisVal;
				if(poisVal>0 && theoVal>0) {
					chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
				}	
			}
			cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
			thInd++;
		}
	}
	// Third loop to fill 3rd row with nTag histograms
	for(unsigned int e=0;e<org->nExp;e++) {
		unsigned int nMax;
		unsigned int nMin;
		if(org->neutronInfo[e]==false) { nMin=0; nMax = 1;}
		if(org->neutronInfo[e]==true) {nMin=1; nMax = 3;}
		unsigned int tN = 0;
		unsigned int tE = 0;
		TString hType = "nTag";
		tempHists.push_back(
				org->exps[e]->makeTotalHist( intFracs, sysVals, hType, tN, tE, nameStart ) );
		nBins = tempHists[thInd]->GetNbinsX(); 
		for(unsigned int b=0;b<nBins;b++) {
			double theoVal = tempHists[thInd]->GetBinContent(b+1);
			double poisVal = org->exps[e]->mocks[fitNo]->hNeutProj->GetBinContent(b+1);
			double poisErrSqr = poisVal;
			if(poisVal>0 && theoVal>0) {
				chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
			}	
		}
		cout << nameStart << " " << hType << " nBins: " << nBins << " chi2Val: " << chi2Val << endl; 
		thInd++;
	}
	
	fitResults[fitNo].nDataEvents = nDataEvents;
	if(type == "fit") {
		fitResults[fitNo].fitChi2 = chi2Val;
		fitResults[fitNo].fitEvisIntegral = evisIntegral;
		fitResults[fitNo].fitEvisBeginDiff = evisBeginDiff;
	}
	else if(type == "pois") {
		fitResults[fitNo].simPoisChi2 = chi2Val;
		fitResults[fitNo].simEvisIntegral = evisIntegral;
		fitResults[fitNo].simEvisBeginDiff = evisBeginDiff;
	}
	else if(type == "sim") {
		fitResults[fitNo].simTrueChi2 = chi2Val;
	}

	org->outFile->cd();
	for(unsigned int v=0;v<tempHists.size();v++) {
		tempHists[v]->Write();
		delete tempHists[v];
	}
	cout << "Saved total plots, fitNo: " << fitNo << " of type: " << type << endl;

} // end of printFitGraphs


// Combine interaction and systematic variables into a single vector
vector<double> FitGlobal::combIntAndSys(vector<double> intVec, vector<double> sysVec) {
	vector<double> combVec;
	combVec.resize(nDim);
	for(unsigned int i=0; i<org->nInt; i++) combVec[i] = intVec[i];
	for(unsigned int s=org->nInt; s<nDim; s++) combVec[s] = sysVec[s-org->nInt];
	return combVec;
}
	
// separate interaction and systematic variables from array pointer to two vectors
void FitGlobal::sepIntAndSys(const double *inPars, vector<double> &intVec, vector<double> &sysVec) {
	intVec.resize(org->nInt);
	sysVec.resize(org->nSys);
	for(unsigned int d=0; d<nDim; d++) {
		if( d<org->nInt ) intVec[d] = inPars[d];
		else sysVec[d-org->nInt] = inPars[d];
	}	
}

// gives the likelihood while changing one parameter and keeping others at true values
vector<double> FitGlobal::scanParameter(unsigned int index, vector<double> inValues, bool isSim) {
	// Get true values (with pois fluctuation)
	vector<double> inVec;
	if( isSim==true ) {
		inVec = combIntAndSys(fitResults[mockIndex].simPoisFracs, fitResults[mockIndex].simPoisSys);
	}
	else { 
		inVec = combIntAndSys(fitResults[mockIndex].fitFracs, fitResults[mockIndex].fitSys);
	}
	vector<double> scanResult; scanResult.resize(inValues.size());
	for(unsigned int i=0; i<inValues.size(); i++) {
		inVec[index] = inValues[i];
		scanResult[i] = getLikelihood( inVec.data() );
		//scanResult[i] = getLikelihoodWei( inVec.data() );
	}
	return scanResult;
} // end of scan parameter function

// scan requested parameters (entered from card) around best fit point
// call only after performing a fit
void FitGlobal::scanAfterFit() {
	double stepInt = 0.1;
	double stepSys = 0.1;
	vector<TGraph*> scanGraphs;
	double totalNllDiff = 0;
	unsigned int scannedCtr = 0;
	for(unsigned int d=0; d<nDim; d++) {
		vector<double> outValues;
		vector<double> inValues;
		bool isScanned = false;
		TString varName;
		if( d<org->nInt ) {
			if( org->intScan[d]==true ) {
				isScanned = true;
				scannedCtr++;
			}
			inValues = fillLinear(0.0,3.0,stepInt);
			varName = Form("scanGraph_Mock%d",mockIndex) + org->intNames[d];
		}
		else {
			if( org->sysScan[d-org->nInt]==true ) {
				isScanned = true;	
				scannedCtr++;
			}
			stepSys = stepInt*org->sysSigmas[d-org->nInt];
			inValues = fillLinear( -3.0*org->sysSigmas[d-org->nInt],
				3.0*org->sysSigmas[d-org->nInt], stepSys);
			varName = Form("scanGraph_Mock%d",mockIndex) + org->sysNames[d-org->nInt];
		}
		if( isScanned == true ) {
			outValues = scanParameter(d, inValues, false);
		
			// Print to output to see
			printVector(inValues, (string)varName);
			printVector(outValues, "Likelihood Value");
		
			unsigned int minIndex = getMinIndex(outValues);
			double minLh = outValues[minIndex];
			double minPos = inValues[minIndex];
			double defLh;
			double lhDiff;
		
			TGraph* tempG = new TGraph(inValues.size(), inValues.data(), outValues.data() );
			if( d<org->nInt ) {
				defLh = tempG->Eval(fitResults[mockIndex].simPoisFracs[d]);
			}
			else {
				defLh = tempG->Eval(fitResults[mockIndex].simPoisFracs[d-org->nInt]);
			}
			lhDiff = defLh - minLh;
			totalNllDiff += lhDiff;
			cout << "Min Pos: " << minPos << " minLh: " << minLh << " trueLh: " << defLh << endl;
			cout << "NLL Diff for this var: " << lhDiff << endl;
			
			tempG->SetTitle(varName);
			tempG->SetName(varName);
			tempG->SetLineColor(4);
			tempG->SetLineWidth(2);
			tempG->GetXaxis()->SetTitle(varName);
			tempG->GetYaxis()->SetTitle("NLL Value");
			scanGraphs.push_back(tempG);
		} // end of if scanned
	} // end of for loop over vars
	
	// Print graphs
	for(unsigned int s=0; s<scannedCtr; s++) {
		if(s>=50) break; // cannot print more than 50 for now
		org->cannyScan->cd(s+1);
		scanGraphs[s]->Draw("APC");
	}
	org->cannyScan->Print(org->outPdfFileScans);
}


// calls scanParameter for all parameters in a reasonable range, and prints results into graphs
void FitGlobal::scanAll() {
	// init scan size and other pars
	double stepInt = 0.2;
	double stepSys;
	vector<TGraph*> scanGraphs;
	double totalNllDiff = 0;
	for(unsigned int d=0; d<nDim; d++) {
		// Create input vector and do the scan
		vector<double> outValues;
		vector<double> inValues;
		TString varName;
		if( d<org->nInt ) {
			inValues = fillLinear(0.0,2.0,stepInt);
			varName = "scanGraph_" + org->intNames[d];
		}
		else {
			stepSys = 0.25*org->sysSigmas[d-org->nInt];
			inValues = fillLinear( -3.0*org->sysSigmas[d-org->nInt],
				3.0*org->sysSigmas[d-org->nInt], stepSys);
			varName = "scanGraph_" + org->sysNames[d-org->nInt];
		}
		outValues = scanParameter(d, inValues);
		
		// Print to output to see
		printVector(inValues, (string)varName);
		printVector(outValues, "Likelihood Value");
		
		unsigned int minIndex = getMinIndex(outValues);
		double minLh = outValues[minIndex];
		double minPos = inValues[minIndex];
		double defLh;
		double lhDiff;
		
		TGraph* tempG = new TGraph(inValues.size(), inValues.data(), outValues.data() );
		if( d<org->nInt ) {
			defLh = tempG->Eval(1.0);
		}
		else {
			defLh = tempG->Eval(0.0);
		}
		lhDiff = defLh - minLh;
		totalNllDiff += lhDiff;
		cout << "Min Pos: " << minPos << " minLh: " << minLh << " trueLh: " << defLh << endl;
		cout << "NLL Diff for this var: " << lhDiff << endl;
		cout << "Total NLL Diff so far: " << totalNllDiff << endl << endl;
		
		
		tempG->SetTitle(varName);
		tempG->SetName(varName);
		tempG->SetLineColor(4);
		tempG->SetLineWidth(2);
		tempG->GetXaxis()->SetTitle(varName);
		tempG->GetYaxis()->SetTitle("NLL Value");
		scanGraphs.push_back(tempG);
		//Print the scan
	}
	
	// Print graphs
	for(unsigned int d=0; d<nDim; d++) {
		if(d>=50) break; // cannot print more than 50 for now
		org->cannyScan->cd(d+1);
		scanGraphs[d]->Draw("APC");
	}
	org->cannyScan->Print(org->outPdfFileScans);
	

} // end of scan all function


// Definition of likelihood function used for fitting
double FitGlobal::getLikelihood(const double *inPars) {
	// Checks and init:
	if(dataIsSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	
	// Separate int fractions and systematics, update all PDFs with systematics:
	vector<double> intFracs;
	vector<double> sysVals;
	sepIntAndSys(inPars, intFracs, sysVals);
	org->updateCurSys(sysVals);
	org->updateAllPdfs();
	
	// Init necessary tools to get likelihood
	double lh = 0;
	double extendedTerm = 0;
	double probTerm = 0;
	double logProbTerm = 0;
	double pullTerm = 0;
	double constTerm = 0; // technicall can be ignored, but good to check if likelihood is positive

	// Sum over experiments and interactions
	for(unsigned int e=0;e<org->nExp;e++) {
		double dPoints = (double) org->exps[e]->mocks[mockIndex]->nPoints;
		constTerm += dPoints*log(dPoints) - dPoints;
		for(unsigned int d=0;d<org->exps[e]->mocks[mockIndex]->nPoints;d++) {
			probTerm = 0;
			for(unsigned int p=0;p<org->nInt;p++) {
				if(d==0) extendedTerm += intFracs[p] * org->exps[e]->pdfs[p]->expEvents;
				probTerm += intFracs[p] * org->exps[e]->pdfs[p]->expEvents * 
					org->exps[e]->pdfs[p]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);

			}
			logProbTerm -= log(probTerm);
		}
	}	

	// Now loop over systematics to add pull terms:
	for(unsigned int s=0; s<org->nSys; s++) {
		// CENTERED...
		pullTerm += 0.5*pow( org->sysCurValues[s]/org->sysSigmas[s], 2.0);
		//pullTerm += 0.5*pow( (org->sysCurValues[s]-fitResults[mockIndex].simTrueSys[s]) 
		//	/ org->sysSigmas[s], 2.0);
		constTerm += log( sqrt(2*TMath::Pi())*org->sysSigmas[s] );
	}

	lh = extendedTerm + logProbTerm + pullTerm + constTerm;

	cout.precision(9);
	//cout << "First sys value is: " << sysVals[0] << endl;
	cout << "Likelihood is: " << lh << endl;
	//cout << "SK4 non-nueCC Event count: " << org->exps[3]->pdfs[2]->expEvents << endl;
	//cout << "Extended Term is: " << extendedTerm << endl;
	cout.precision(6);

	return lh;
}

// Definition of likelihood gradient used in fitting
double FitGlobal::getLhGradient(const double *inPars, unsigned int coord) {
	if(dataIsSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	vector<double> intFracs;
	vector<double> sysVals;
	sepIntAndSys(inPars, intFracs, sysVals);

	double grad = 0;
	// Interaction fraction gradients can be calculated analytically
	if(coord < org->nInt) { 
		
		// Set systematic values and update pdfs
		org->updateCurSys(sysVals);
		org->updateAllPdfs();
		// Init
		double extendedTerm = 0;
		double probTermUp = 0;
		double probTermBot = 0;
		double probTerm = 0;
		
		// Sum over experiments and interactions
		for(unsigned int e=0;e<org->nExp;e++) {
			extendedTerm += org->exps[e]->pdfs[coord]->expEvents;
			for(unsigned int d=0;d<org->exps[e]->mocks[mockIndex]->nPoints;d++) {
				// Top part of prob term is over interaction we are differentiating over only
				probTermUp = org->exps[e]->pdfs[coord]->expEvents * 
					org->exps[e]->pdfs[coord]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);
				// Bottom part of prob term is summed over all ineractions:
				probTermBot = 0;
				for(unsigned int p=0;p<org->nInt;p++) {
					probTermBot += intFracs[p] * org->exps[e]->pdfs[p]->expEvents *
						org->exps[e]->pdfs[p]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);
				}
				probTerm -= probTermUp / probTermBot;
			}
		}
		// Pull terms drop of at differentiation
		grad = extendedTerm + probTerm;
		
		cout << "Lh gradient for dimension: " << org->intNames[coord] << " is: " << grad << endl;
	}
	// Systematic gradients change PDFs numerically, so need to differentiate numerically
	else {
		// Determine systematic index
		unsigned int sInd = coord - org->nInt;
		// Get likelihood half a step left (negative in the systematic)
		vector<double> sysVecLeft = sysVals; 
		sysVecLeft[sInd] -= org->sysStepSize[sInd]/2.0;
		vector<double> inVecLeft = combIntAndSys(intFracs, sysVecLeft);
		double lhLeft = getLikelihood( inVecLeft.data() );
		// Get likelihood half a step right (positive in the systematic)
		vector<double> sysVecRight = sysVals;
		sysVecRight[sInd] += org->sysStepSize[sInd]/2.0;
		vector<double> inVecRight = combIntAndSys(intFracs, sysVecRight);
		double lhRight = getLikelihood( inVecRight.data() );
		// Calculate gradient numerically
		grad = (lhRight-lhLeft) / org->sysStepSize[sInd];
		
		cout << "Lh gradient for dimension: " << org->sysNames[sInd] << " is: " << grad << endl;
	}
	
	return grad;
}


// Definition of likelihood function used for fitting in case of weighted data
double FitGlobal::getLikelihoodWei(const double *inPars) {
	// Checks and init:
	if(dataIsSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	
	// Separate int fractions and systematics, update all PDFs with systematics:
	vector<double> intFracs;
	vector<double> sysVals;
	sepIntAndSys(inPars, intFracs, sysVals);
	org->updateCurSys(sysVals);
	org->updateAllPdfs();
	
	// Init necessary tools to get likelihood
	double lh = 0;
	double extendedTerm = 0;
	double probTerm = 0;
	double logProbTerm = 0;
	double pullTerm = 0;
	double constTerm = 0; // technicall can be ignored, but good to check if likelihood is positive

	// Sum over experiments and interactions
	for(unsigned int e=0;e<org->nExp;e++) {
		double dPoints = (double) org->exps[e]->mocks[mockIndex]->nPoints;
		double dWeights = (double) org->exps[e]->mocks[mockIndex]->totalWeight;
		constTerm += dWeights*log(dWeights) - dWeights;
		for(unsigned int d=0;d<org->exps[e]->mocks[mockIndex]->nPoints;d++) {
			probTerm = 0;
			for(unsigned int p=0;p<org->nInt;p++) {
				if(d==0) extendedTerm += intFracs[p] * org->exps[e]->pdfs[p]->expEvents;
				probTerm += intFracs[p] * org->exps[e]->pdfs[p]->expEvents * 
					org->exps[e]->pdfs[p]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);

			}
			logProbTerm -= org->exps[e]->mocks[mockIndex]->dataPoint[d].weight * log(probTerm);
		}
	}	

	// Now loop over systematics to add pull terms:
	for(unsigned int s=0; s<org->nSys; s++) {
		// *********** TEMPORARY CHANGE FOR TEST, TO TRY RANDOM SYS BUT MAKING PULL ALSO
		// CENTERED...
		pullTerm += 0.5*pow( org->sysCurValues[s]/org->sysSigmas[s], 2.0);
		//pullTerm += 0.5*pow( (org->sysCurValues[s]-fitResults[mockIndex].simTrueSys[s]) 
		//	/ org->sysSigmas[s], 2.0);
		constTerm += log( sqrt(2*TMath::Pi())*org->sysSigmas[s] );
	}

	lh = extendedTerm + logProbTerm + pullTerm + constTerm;

	cout << "Likelihood is: " << lh << endl;

	return lh;
}

// Definition of likelihood gradient used in fitting
double FitGlobal::getLhGradientWei(const double *inPars, unsigned int coord) {
	if(dataIsSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	vector<double> intFracs;
	vector<double> sysVals;
	sepIntAndSys(inPars, intFracs, sysVals);

	double grad = 0;
	// Interaction fraction gradients can be calculated analytically
	if(coord < org->nInt) { 
		
		// Set systematic values and update pdfs
		org->updateCurSys(sysVals);
		org->updateAllPdfs();
		// Init
		double extendedTerm = 0;
		double probTermUp = 0;
		double probTermBot = 0;
		double probTerm = 0;
		
		// Sum over experiments and interactions
		for(unsigned int e=0;e<org->nExp;e++) {
			extendedTerm += org->exps[e]->pdfs[coord]->expEvents;
			for(unsigned int d=0;d<org->exps[e]->mocks[mockIndex]->nPoints;d++) {
				// Top part of prob term is over interaction we are differentiating over only
				probTermUp = org->exps[e]->pdfs[coord]->expEvents * 
					org->exps[e]->pdfs[coord]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);
				// Bottom part of prob term is summed over all ineractions:
				probTermBot = 0;
				for(unsigned int p=0;p<org->nInt;p++) {
					probTermBot += intFracs[p] * org->exps[e]->pdfs[p]->expEvents *
						org->exps[e]->pdfs[p]->getProb(org->exps[e]->mocks[mockIndex]->dataPoint[d]);
				}
				probTerm -= org->exps[e]->mocks[mockIndex]->dataPoint[d].weight * 
										(probTermUp / probTermBot);
			}
		}
		// Pull terms drop of at differentiation
		grad = extendedTerm + probTerm;
		
		cout << "Lh gradient for dimension: " << org->intNames[coord] << " is: " << grad << endl;
	}
	// Systematic gradients change PDFs numerically, so need to differentiate numerically
	else {
		// Determine systematic index
		unsigned int sInd = coord - org->nInt;
		// Get likelihood half a step left (negative in the systematic)
		vector<double> sysVecLeft = sysVals; 
		sysVecLeft[sInd] -= org->sysStepSize[sInd]/2.0;
		vector<double> inVecLeft = combIntAndSys(intFracs, sysVecLeft);
		double lhLeft = getLikelihoodWei( inVecLeft.data() );
		// Get likelihood half a step right (positive in the systematic)
		vector<double> sysVecRight = sysVals;
		sysVecRight[sInd] += org->sysStepSize[sInd]/2.0;
		vector<double> inVecRight = combIntAndSys(intFracs, sysVecRight);
		double lhRight = getLikelihoodWei( inVecRight.data() );
		// Calculate gradient numerically
		grad = (lhRight-lhLeft) / org->sysStepSize[sInd];
		
		cout << "Lh gradient for dimension: " << org->sysNames[sInd] << " is: " << grad << endl;
	}
	
	return grad;
}
