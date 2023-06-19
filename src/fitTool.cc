/********
*Author: Baran Bodur
*Date: 2021-07-13
*Description: Fit Pdfs to data by using the likelihood functions defined here
*	as well as the minimizer. It also creates some fit hists/graphs in pdf and root format
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

using namespace std;

// Constructor
FitTool::FitTool(TString inName, const unsigned int nPdf, FitPdf *inPdf[]) {
	// Set namd & pdf no
	name = inName;
	nType = nPdf;

	// Reserve spaces for vectors
	fitSteps.resize(nType);
	curPars.resize(nType);
	initPars.resize(nType);
	
	//fitPars.reserve(nType);
	//simTruePars.reserve(nType);
	//simPoisPars.reserve(nType);
	
	h1FitDist.reserve(nType);

	// Save pdfs in this object
	for(unsigned int p=0; p<nType; p++) {
			pdfs.push_back(inPdf[p]);			
			fitPars.push_back(0);
			simTruePars.push_back(0);
			simPoisPars.push_back(0);
	}

	// Create minimizer object

	
	minFarshaw = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
	minFarshaw->SetMaxFunctionCalls(100000); // max calls, reduce if it takes too much time 
	minFarshaw->SetMaxIterations(30000); // not used
	minFarshaw->SetTolerance(0.05); // tolerance for EDM
	minFarshaw->SetStrategy(2); // 0:coarse 1: normal 2: detailed 
	minFarshaw->SetPrintLevel(0); 
	
	// Set initial pars and fit steps (temporary can be set later externally
	double tempInit[nPdf];
	double tempStep[nPdf];
	for(unsigned int p=0; p<nType; p++) {
		tempInit[p] = 100;
		tempStep[p] = 1;
	}
	setInitPars(tempInit);
	setFitStep(tempStep);
	setSimPars(); // based on pdfs ( sets only simTruePars, not simPoisPars)

	// Create functor to get value & gradient of the function to be minimized at given point
	getLh = &FitTool::getLikelihood;
	getGrad = &FitTool::getLhGradient;

	fitFunc = ROOT::Math::GradFunctor(this,getLh,getGrad,nType); // functor to be used	
	//fitFunc = ROOT::Math::GradFunctor(&getLikelihood,&getLhGradient,nType); // functor to be used	
	//cout << "Dimensions: " <<fitFunc.NDim() << endl;
	
	// Bind the fit function to the minimizer
	minFarshaw->SetFunction(fitFunc);
	
	// Bind variables
	for(unsigned int p = 0; p < nType; p++) {
			minFarshaw->SetLowerLimitedVariable(p,(string)pdfs[p]->pdfName,
					tempInit[p],tempStep[p],0.0);
	}

	//Prepare fit histograms
	prepFitHists();
	
	stackHistCtr = 0;

} // end of constructor

//Destructor
FitTool::~FitTool() {
/*
	delete minFarshaw;
	
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

void FitTool::resetFit() { // sets fit function, init vals, steps again
	
	fitFunc = ROOT::Math::GradFunctor(this,getLh,getGrad,nType); // functor to be used	
	// Bind the fit function to the minimizer
	minFarshaw->SetFunction(fitFunc);
	// Bind variables
	for(unsigned int p = 0; p < nType; p++) { // set vars again
			minFarshaw->SetLowerLimitedVariable(p,(string)pdfs[p]->pdfName,
					initPars[p],fitSteps[p],0.0);
	} // end of loop over fit vars
} // end of resetFit



// setData to be used in the fit
void FitTool::setData(MockData *inData) {
	dataToFit = inData; // set data pointer
	dataWasSet = true; // set dataWasSet flag
	setSimPars();
} // end of setData

// set initial parameters
void FitTool::setInitPars(double *inPars) {
	for(unsigned int p=0; p<nType; p++) initPars[p] = inPars[p];	
} // end of setInitPars

// set fit step in directions
void FitTool::setFitStep(double *inSteps) {
	for(unsigned int p=0; p<nType; p++) fitSteps[p] = inSteps[p];	
} // end of setFitStep

// setSimPars from input pdfs (simTruth) and mock data (sim Result with fluctuations)
void FitTool::setSimPars() { 
	for(unsigned int p=0; p<nType; p++) {
		simTruePars[p] = pdfs[p]->expEvents; 
		if(dataWasSet == true) simPoisPars[p] = dataToFit->byTypeEvents[p];
	}
	simTrueLh = getLikelihood(simTruePars.data()); // get likelihood with true sim Params
	if(dataWasSet == true) simPoisLh = getLikelihood(simPoisPars.data()); //with fluctuated params
} // end of setSimPars

// run the fit (save result) : makeFit
double FitTool::makeFit() {
	minFarshaw->Minimize();
	const double *fitResult = minFarshaw->X();
	for(unsigned int p=0; p<nType; p++) fitPars[p] = fitResult[p];	
	cout << "FitParsVec size: " << fitPars.size() << endl;
	//cout << "True Pars 1: " << *(simTruePars.data()+3) << endl;
	cout << "# of function calls for min:" << minFarshaw->NCalls() << endl; ;
	fitMinLh = minFarshaw->MinValue();
	return fitMinLh;
} // end of makeFit


// print info about last fit result
void FitTool::printLastFit() {
	return; // will write this later
} // end of printLastFit

// prepare fit histograms to be filled
void FitTool::prepFitHists() {


	// Initialize histograms
	hFitLh = new TH1F("hFitLh_"+name,"hFitLh",10000,0,10000);
	hSimPoisLh = new TH1F("hSimPoisLh_"+name,"hSimPoisLh",10000,0,10000);
	hSimTrueLh = new TH1F("hSimTrueLh_"+name,"hSimTrueLh",10000,0,10000);
	hLhFitMinPois = new TH1F("hLhFitMinPois_"+name,"hLhFitMinPois",1000,-50,50);
	hLhFitMinTrue = new TH1F("hLhFitMinTrue_"+name,"hLhFitMinTrue",1000,-50,50);
	// Loop over histograms
	for(unsigned int h = 0; h < nType; h++) {
		TString hName1 = "h1FitDist_" +name + "_" + pdfs[h]->pdfName;	
		h1FitDist[h] = new TH1F(hName1,hName1,500,0,10*pdfs[h]->expEvents);
		TString hAxis1 = "Fitted " + pdfs[h]->pdfName + " Events";
		h1FitDist[h]->GetXaxis()->SetTitle(hAxis1); 
		h1FitDist[h]->SetLineWidth(2); 
		h1FitDist[h]->SetLineColor(4); 
		vector<TH2F*> tempVec;
		for(unsigned int h2 = 0; h2 < h; h2++) {
			TString hName2 = "h1FitDist_" + name + "_X_" + pdfs[h]->pdfName + "_Y_" + pdfs[h2]->pdfName; 
			TString hAxis2 = "Fitted " + pdfs[h2]->pdfName + " Events";
			TH2F* tempH = new TH2F(hName2,hName2,100,0,10*pdfs[h]->expEvents,
					100,0,pdfs[h2]->expEvents);
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
void FitTool::fillFitHists() {
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
void FitTool::saveFitHists(TFile* outFile) {
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

void FitTool::getStackedHists() { // make stacked hists with fit/sim/pois pars
	unsigned int nProjHists = pdfs[0]->projHists.size();
	vector<THStack*> htFit, htSim, htPois;
	vector<TH1D*> htData;
	for(unsigned int h=0; h < nProjHists; h++) {
		htFit.push_back( makeStackedHist(fitPars.data(),h,Form("hStackFit%d_",stackHistCtr) 
					+name ) );
		htSim.push_back( makeStackedHist(simTruePars.data(),h,Form("hStackSim%d_",stackHistCtr) 
					+name ) );
		htPois.push_back( makeStackedHist(simPoisPars.data(),h,Form("hStackPois%d_",stackHistCtr)
					+name) );
		htData.push_back( dataToFit->projHists[h] );
	} // end of histogram type loop
	hFitProj.push_back(htFit);
	hSimProj.push_back(htSim);
	hPoisProj.push_back(htPois);
	hDataProj.push_back(htData);

	// Also save the fit results:
	fitParsVec.push_back(fitPars);	
	simTrueParsVec.push_back(simTruePars);	
	simPoisParsVec.push_back(simPoisPars);	
	fitMinLhVec.push_back(fitMinLh);
	simTrueLhVec.push_back(simTrueLh);
	simPoisLhVec.push_back(simPoisLh);

	stackHistCtr++;
} // end of getStackedHists

void FitTool::saveStackedHists(TFile* outFile) { // save stacked hists
	unsigned int nProjHists = pdfs[0]->projHists.size();
	unsigned int nFit = hFitProj.size();
	outFile->cd();
	for(unsigned int f=0; f < nFit; f++) {
		for(unsigned int h=0; h < nProjHists; h++) {
			hFitProj[f][h]->Write();
			hSimProj[f][h]->Write();
			hPoisProj[f][h]->Write();
		}
	}

}// end of saveStackedHists

void FitTool::delStackedHists() { // delete prev stacked hists
	unsigned int nProjHists = pdfs[0]->projHists.size();
	unsigned int nFit = hFitProj.size();
	for(unsigned int f=0; f < nFit; f++) {
		for(unsigned int h=0; h < nProjHists; h++) {
			delete hFitProj[f][h];
			delete hSimProj[f][h];
			delete hPoisProj[f][h];
		}
	}
	hFitProj.clear();
	hSimProj.clear();
	hPoisProj.clear();
	hDataProj.clear(); // do not delete these before, because they are associated with data
	fitParsVec.clear();
	simTrueParsVec.clear();
	simPoisParsVec.clear();
	fitMinLhVec.clear();
	simTrueLhVec.clear();
	simPoisLhVec.clear();

	stackHistCtr = 0;

} // end of delStackedHists


// Make Stacked Plot With Given Parameters, Projection and Name
THStack* FitTool::makeStackedHist(const double *inPars, unsigned int hNo, TString nameStart) {	
	// Sleect names for hNo
	TString nameEnds[9] = {"projEvisN0","projEvisN1","projEvisNall",
		"projGtagN0","projGtagN1","projGtagNall","projNtag","projCosZen","projAzi"};
	TString xTitles[9] = {"Visible Energy N0 (MeV)", "Visible Energy N1(MeV)", 
		"Visible Energy Nall (MeV)","Gamma Tag Par. N0","Gamma Tag Par. N1","Gamma Tag Par. Nall",
		"Tagged Neutrons","Cos(Zenith)","Azimuth Angle (Degrees)"};
	TString sName = nameStart+"_"+nameEnds[hNo];

	THStack *hAnaStack = new THStack(sName,sName);
	TH1D* hIntPlot[nType];
	for(unsigned int p=0; p < nType; p++) {
		hIntPlot[p] = (TH1D*) pdfs[p]->projHists[hNo]->Clone();
		if( hNo < 3 ) { hIntPlot[p]->Rebin(5); }
		else if( hNo >= 3 && hNo < 6 ) { hIntPlot[p]->Rebin(2);}
		hIntPlot[p]->GetXaxis()->SetTitle(xTitles[hNo]);
		hIntPlot[p]->GetXaxis()->SetTitleOffset(0.8);
		hIntPlot[p]->Scale(inPars[p]);
		hIntPlot[p]->SetFillColor(pdfs[p]->pdfColor);
		hIntPlot[p]->SetFillStyle(1001);
		hIntPlot[p]->SetLineColor(pdfs[p]->pdfColor);
		hAnaStack->Add(hIntPlot[p]);
	} // end of loop over pdfs to be stacked
	return hAnaStack;
} // end of makeStackedHist

void FitTool::printFitRes(TString outPdfName) {
	
	unsigned int nProjHists = pdfs[0]->projHists.size();
	unsigned int nFit = hFitProj.size();
	if(nFit < 1 || nProjHists < 7) {
		cout << "Not Enough Hists or Fits to Print, nFit: " << nFit << " nProjHist: "<<nProjHists<<endl;
		return;
	}
	TCanvas* canny = new TCanvas("canny","canny",1200,2700);
	unsigned int pType = 3;
	unsigned int gType = 9;
	canny->Divide(pType,gType);
	canny->Print(outPdfName+"[");
	TText* texty = new TText();
	texty->SetTextFont(10);
	float tX = 0.6, tY = 0.6, dY = 0.05;
	for(unsigned int f=0; f < nFit; f++) {
		for(unsigned int h=0; h < nProjHists; h++) {
			hDataProj[f][h]->SetLineColor(6);
			if( h < 3 ) hDataProj[f][h]->Rebin(5);
			else if( h >= 3 && h < 6 ) hDataProj[f][h]->Rebin(2);
			
			//hDataProj[f][h]->SetLineWidth(2);
			float maxY = 1.5*hDataProj[f][h]->GetMaximum();
			hDataProj[f][h]->GetYaxis()->SetRangeUser(0,maxY);
			// Fit Stacks
			canny->cd(1+3*h);
			hDataProj[f][h]->Draw();
			hFitProj[f][h]->Draw("same hist");
			hDataProj[f][h]->Draw("same");
			texty->DrawTextNDC(tX,tY,Form("Fit Min. Lh: %.1f",fitMinLhVec[f]));
			for(unsigned int p=0; p < nType; p++) {
				texty->DrawTextNDC(tX,tY-(p+1)*dY,
						(TString)pdfs[p]->pdfName+Form(" : %.1f",fitParsVec[f][p]));
			}			
			// Sim Stacks
			canny->cd(2+3*h);
			hDataProj[f][h]->Draw();
			hSimProj[f][h]->Draw("same hist");
			hDataProj[f][h]->Draw("same");
			texty->DrawTextNDC(tX,tY,Form("Sim True Lh: %.1f",simTrueLhVec[f]));
			for(unsigned int p=0; p < nType; p++) {
				texty->DrawTextNDC(tX,tY-(p+1)*dY,
						(TString)pdfs[p]->pdfName+Form(" : %.1f",simTrueParsVec[f][p]));
			}			
			// Pois Stacks
			canny->cd(3+3*h);
			hDataProj[f][h]->Draw();
			hPoisProj[f][h]->Draw("same hist");
			hDataProj[f][h]->Draw("same");
			texty->DrawTextNDC(tX,tY,Form("Sim Pois. Lh: %.1f",simPoisLhVec[f]));
			for(unsigned int p=0; p < nType; p++) {
				texty->DrawTextNDC(tX,tY-(p+1)*dY,
						(TString)pdfs[p]->pdfName+Form(" : %.1f",simPoisParsVec[f][p]));
			}			

		} //end of loop over projection hists
		canny->Print(outPdfName);
		cout << " Printing Fit No: " << f << "/" << nFit << endl;
	} //end of loop over # of fit
	
	canny->Print(outPdfName+"]");
	delete canny;
}


// ** Member function variation of likelihood/gradient functions

// calculate likelihood given parameters
double FitTool::getLikelihood(const double *inPars){
	//cout << nType << endl;
	//cout << inPars[0] << endl;
	//cout << "In getLikelihood!" << endl;
	if(dataWasSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	else { // only if data is set
		double totalRate = 0;
		double lh=0, term = 0;
		for(unsigned int d=0; d<dataToFit->nPoints; d++) { // loop over data points
			term = 0; // reset term at each data point
			for(unsigned int p=0; p<nType; p++) { // loop over pdfs
				//cout << "DataPoint of : " << d << " of: " << dataToFit->nPoints << endl;
				if(d==0) totalRate += inPars[p]; // calculate total rate with given pars
				term += ( inPars[p]*pdfs[p]->getProb(dataToFit->dataPoint[d]) ); 
			}
			lh -= log(term); // convert term to NLL contribution from a single data point 
		}
		lh += totalRate; // overall rate term in NLL
		
		// Temporary modification to add restriction, improve later!
		//float p3Exp = 386.3;
		//float p4Exp = 140.7;
		//float p3Err = 0.6 * p3Exp;
		//float p4Err = 0.1 * p4Exp;
		//if( nType > 3 ) lh += 0.5*pow((inPars[3] - p3Exp),2)/pow(p3Err,2);
		//if( nType > 4 ) lh += 0.5*pow((inPars[4] - p4Exp),2)/pow(p4Err,2);

		return lh;
	} //end of dataWasSet else
} // end of getLikelihood

// calculate likelihood gradient given parameters and direction
double FitTool::getLhGradient(const double *inPars, unsigned int coord) {
	//cout << "In getLhGradient" << endl;
	if(dataWasSet == false) {
		cout << "No data was set for the fit, returning 0" << endl;
		return 0;
	}
	else { // only if data is set
		double grad=1.0; // init gradient from 1
		double top=0, bot=0;
		for(unsigned int d=0; d<dataToFit->nPoints; d++) { // loop over data points
			top = pdfs[coord]->getProb(dataToFit->dataPoint[d]); // numerator
			bot = 0; //reset denominator
			for(unsigned int p=0; p<nType; p++) { // loop over pdfs
				bot += ( inPars[p]*pdfs[p]->getProb(dataToFit->dataPoint[d]) ); //up denom
			}
			grad -= top/bot; // subtract contribution of each data point from NLLs derivative
		}
	
		// Temporary modification to add restriction, improve later!
		//float p3Exp = 386.3;
		//float p4Exp = 140.7;
		//float p3Err = 0.6 * p3Exp;
		//float p4Err = 0.1 * p4Exp;
		//if( coord == 3 ) grad += (inPars[3] - p3Exp)/pow(p3Err,2); 
		//if( coord == 4 ) grad += (inPars[4] - p4Exp)/pow(p4Err,2); 

		return grad;	
	} //end of dataWasSet else

} // end of getLhGradient

