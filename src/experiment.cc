/********
*Author: Baran Bodur
*Date: 2022-07-26
*Description: Holds interaction PDFs and mock (or real) data for the given experiment
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
#include "organizer.h" 
#include "experiment.h" 
#include "mockData.h" 
#include "fitKdePdf.h"
#include "vectorUtils.h"

using namespace std;

// Functions go here

// Constructor
Experiment::Experiment(Organizer *orgIn, unsigned int expIndexIn) {
	org = orgIn;
	expIndex = expIndexIn;
	cout << "Creating Experiment Named: " << org->expNames[expIndex] << endl;
	// Create interaction pdfs for this experiment
	pdfs.resize(org->nInt);
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i] = new FitKdePdf(org,expIndex,i);
	}
	// Create mock data for this experiment (only after interactions are initialized)
	mocks.resize(org->nFit);
	for(unsigned int m=0; m<org->nFit; m++) {
		mocks[m] = new MockData(org,expIndex,m);
	}

}

// Reweights 1 systematic with given index and sigma to all systematics at 0 case
void Experiment::genMockData(vector< vector<double> > usedSysVecIn) {
	for(unsigned int m=0; m<org->nFit; m++) {
		mocks[m]->genMockData(usedSysVecIn[m]);
		mocks[m]->saveProjHists();
	}	
}

void Experiment::reweightOneSys(unsigned int sysIndex, double nSigma, unsigned int mockIndexIn) {
	//for(unsigned int m=0; m<org->nFit; m++) {
	mocks[mockIndexIn]->reweightOneSys(sysIndex, nSigma);
	//mocks[mockIndexIn]->saveProjHists();
	//}

}

void Experiment::updateAllPdfs() {
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->updateAllPdfs();
	}
}

void Experiment::saveCurPdfs(TString saveName) {
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->saveCurPdfs(saveName);
	}
}

void Experiment::updateProjHists(double scale) {
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->updateProjHists(scale);
	}
}

void Experiment::saveCurProj(TString saveName) {
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->saveCurProj(saveName);
	}
}

// Returns desired stacked hist
THStack* Experiment::makeStackedHist(vector<double> intFracs, vector<double> sysVals, 
		TString type, unsigned int nX, unsigned int eX, TString nameStart) {

	// Prepare projection histograms as required:
	org->updateCurSys(sysVals);
	updateAllPdfs();
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->updateProjHists( intFracs[i]*pdfs[i]->expEvents );
	}
	vector<TH1D*> hTempProj;
	hTempProj.resize(org->nInt);

	TString neutNames[3] = { "Nall", "N0", "N1" };
	
	if(type == "nTag") {
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_nTagProj";
		TString sAxis = "Tagged Neutrons"; 
		THStack* hAnaStack = new THStack(sName,sName);
		for(unsigned int i=0; i<org->nInt; i++) {
			hTempProj[i] = (TH1D*) pdfs[i]->hNeutProj->Clone();
			hTempProj[i]->SetFillColor(org->intColors[i]);
			hTempProj[i]->SetFillStyle(1001);
			hTempProj[i]->SetLineWidth(0);
			hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hAnaStack->Add( hTempProj[i] );
		}
		return hAnaStack;
	}
	else if(type == "eVis") {
		if(org->neutronInfo[expIndex]==false && nX!=0) {
			cout << "WARNING! In makeStackedHist, asked eVis hist with specific neutron count,"
				   << " but no neutron info available, returning NULL!"<<endl;
			return NULL;
		}
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_eVisProj_" + neutNames[nX];
		TString sAxis = "Visible Energy (MeV) for " + neutNames[nX]; 
		THStack* hAnaStack = new THStack(sName,sName);
		for(unsigned int i=0; i<org->nInt; i++) {
			hTempProj[i] = (TH1D*) pdfs[i]->evisKdes[nX].histProj->Clone();
			hTempProj[i]->SetFillColor(org->intColors[i]);
			hTempProj[i]->SetFillStyle(1001);
			hTempProj[i]->SetLineWidth(0);
			hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hAnaStack->Add( hTempProj[i] );
		}
		return hAnaStack;	
	}
	else if(type == "gTag") {
		if(org->gammaInfo[expIndex]==false) {
			cout << "WARNING! In makeStackedHist, asked gamma hist but no gamma info available," 
					 << " returning NULL!"<<endl;
			return NULL;
		}
		if(org->neutronInfo[expIndex]==false && nX!=0) {
			cout << "WARNING! In makeStackedHist, asked gamma hist with specific neutron count,"
				   << " but no neutron info available, returning NULL!"<<endl;
			return NULL;
		}
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_gTagProj_" + neutNames[nX] + "_" 
				+ pdfs[0]->getGtagSubbinRange(eX);
		TString sAxis = "Gamma Tag for " + neutNames[nX] + " and " + pdfs[0]->getGtagSubbinRange(eX); 
		THStack* hAnaStack = new THStack(sName,sName);
		for(unsigned int i=0; i<org->nInt; i++) {
			hTempProj[i] = (TH1D*) pdfs[i]->gtagKdes[nX][eX].histProj->Clone();
			hTempProj[i]->SetFillColor(org->intColors[i]);
			hTempProj[i]->SetFillStyle(1001);
			hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			hTempProj[i]->SetLineWidth(0);
			hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hAnaStack->Add( hTempProj[i] );
		}
		return hAnaStack;	
	}
	else {
		cout << "WARNING! In makeStackedHist, Type: "<<type<<" not available, returning NULL!"<<endl;
		return NULL;
	}
}

// --- This is under construction
// Returns desired total hist (instead of stacked 1 TH1D that is the sum
TH1D* Experiment::makeTotalHist(vector<double> intFracs, vector<double> sysVals, 
		TString type, unsigned int nX, unsigned int eX, TString nameStart) {

	// Prepare projection histograms as required:
	org->updateCurSys(sysVals);
	updateAllPdfs();
	for(unsigned int i=0; i<org->nInt; i++) {
		pdfs[i]->updateProjHists( intFracs[i]*pdfs[i]->expEvents );
	}
	vector<TH1D*> hTempProj;
	hTempProj.resize(org->nInt);

	TString neutNames[3] = { "Nall", "N0", "N1" };
	
	if(type == "nTag") {
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_nTagProj";
		TString sAxis = "Tagged Neutrons"; 
		TH1D* hTotal;
		for(unsigned int i=0; i<org->nInt; i++) {
			if(i==0) {
				hTotal = (TH1D*) pdfs[i]->hNeutProj->Clone(sName);
				hTotal->Reset();
			}
			hTempProj[i] = (TH1D*) pdfs[i]->hNeutProj->Clone();
			//hTempProj[i]->SetFillColor(org->intColors[i]);
			//hTempProj[i]->SetFillStyle(1001);
			//hTempProj[i]->SetLineWidth(0);
			//hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			//hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hTotal->Add( hTempProj[i] );
		}
		hTotal->GetXaxis()->SetTitle(sAxis);
		hTotal->GetXaxis()->SetTitleOffset(0.8);
		return hTotal;
	}
	else if(type == "eVis") {
		if(org->neutronInfo[expIndex]==false && nX!=0) {
			cout << "WARNING! In makeStackedHist, asked eVis hist with specific neutron count,"
				   << " but no neutron info available, returning NULL!"<<endl;
			return NULL;
		}
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_eVisProj_" + neutNames[nX];
		TString sAxis = "Visible Energy (MeV) for " + neutNames[nX]; 
		TH1D* hTotal;
		for(unsigned int i=0; i<org->nInt; i++) {
			if(i==0) {
				hTotal = (TH1D*) pdfs[i]->evisKdes[nX].histProj->Clone(sName);
				hTotal->Reset();
			}
			hTempProj[i] = (TH1D*) pdfs[i]->evisKdes[nX].histProj->Clone();
			//hTempProj[i]->SetFillColor(org->intColors[i]);
			//hTempProj[i]->SetFillStyle(1001);
			//hTempProj[i]->SetLineWidth(0);
			//hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			//hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hTotal->Add( hTempProj[i] );
		}
		hTotal->GetXaxis()->SetTitle(sAxis);
		hTotal->GetXaxis()->SetTitleOffset(0.8);
		return hTotal;
	}
	else if(type == "gTag") {
		if(org->gammaInfo[expIndex]==false) {
			cout << "WARNING! In makeStackedHist, asked gamma hist but no gamma info available," 
					 << " returning NULL!"<<endl;
			return NULL;
		}
		if(org->neutronInfo[expIndex]==false && nX!=0) {
			cout << "WARNING! In makeStackedHist, asked gamma hist with specific neutron count,"
				   << " but no neutron info available, returning NULL!"<<endl;
			return NULL;
		}
		TString sName = nameStart + "_" + org->expNames[expIndex] + "_gTagProj_" + neutNames[nX] + "_" 
				+ pdfs[0]->getGtagSubbinRange(eX);
		TString sAxis = "Gamma Tag for " + neutNames[nX] + " and " + pdfs[0]->getGtagSubbinRange(eX); 
		TH1D* hTotal;
		for(unsigned int i=0; i<org->nInt; i++) {
			if(i==0) {
				hTotal = (TH1D*) pdfs[i]->gtagKdes[nX][eX].histProj->Clone(sName);
				hTotal->Reset();
			}
			hTempProj[i] = (TH1D*) pdfs[i]->gtagKdes[nX][eX].histProj->Clone();
			//hTempProj[i]->SetFillColor(org->intColors[i]);
			//hTempProj[i]->SetFillStyle(1001);
			//hTempProj[i]->GetXaxis()->SetTitle(sAxis);
			//hTempProj[i]->SetLineWidth(0);
			//hTempProj[i]->GetXaxis()->SetTitleOffset(0.8);
			hTotal->Add( hTempProj[i] );
		}
		hTotal->GetXaxis()->SetTitle(sAxis);
		hTotal->GetXaxis()->SetTitleOffset(0.8);
		return hTotal;
	}
	else {
		cout << "WARNING! In makeStackedHist, Type: "<<type<<" not available, returning NULL!"<<endl;
		return NULL;
	}
}

// Prints 1 sigma (default) systematic histograms and apply KS test between 0 and 1 sigma
// -1 sigma will also be implemented eventually
void Experiment::printSysHists(double sysShift) {
	
	// number of columns and and rows
	unsigned int nRowsSys = 4;
	unsigned int nColsSys = 2;
	
	// neutron
	unsigned int nInd = 0;
	unsigned int eInd = 0;
	TString sAxisEvis = "Visible Energy (MeV)";
	TString sAxisGtag = "Gamma Tag Output";
	TString sAxisNtag = "Tagged Neutrons";

	// init all sys to 0
	unsigned int nSysLoop = org->nSys;
	//unsigned int nSysLoop = 10;
	vector<double> inSys;
	inSys.resize(org->nSys);
	for(unsigned int s=0; s<nSysLoop;s++) inSys[s] = 0;
	org->updateCurSys(inSys);
	updateAllPdfs();
		
	// To print into text file
	TString outSysTextStr = org->expNames[expIndex] + "_systOut.txt"; 
	ofstream outSys(outSysTextStr);
	
	// First loop to fill first row with eVis histograms
	unsigned int cRow;
	unsigned int cCol;
	double lrMargin = 0.05;
	
	// Loop over systematics to update them to 1 in turn
	for(unsigned int s=0; s<nSysLoop;s++) {
		// Update pdfs according to 1 sigma syst
		inSys[s] = sysShift * org->sysSigmas[s];
		if(s>0) inSys[s-1] = 0; // change previous hist to 0
		org->updateCurSys(inSys);
		updateAllPdfs();

		// Init  and resize temporary histogram holders
		vector<TH1D*> hTEvis0s, hTEvis1s, hTEvism1s;
		vector<TH1D*> hREvis0s, hREvis1s, hREvism1s;
		vector<TH1D*> hTGtag0s, hTGtag1s, hTGtagm1s;
		vector<TH1D*> hRGtag0s, hRGtag1s, hRGtagm1s;
		vector<TH1D*> hTNtag0s, hTNtag1s, hTNtagm1s;
		vector<TH1D*> hRNtag0s, hRNtag1s, hRNtagm1s;
		
		hTEvis0s.resize(org->nInt); hTEvis1s.resize(org->nInt); hTEvism1s.resize(org->nInt);
		hREvis0s.resize(org->nInt); hREvis1s.resize(org->nInt); hREvism1s.resize(org->nInt);
		hTGtag0s.resize(org->nInt); hTGtag1s.resize(org->nInt); hTGtagm1s.resize(org->nInt);
		hRGtag0s.resize(org->nInt); hRGtag1s.resize(org->nInt); hRGtagm1s.resize(org->nInt);
		hTNtag0s.resize(org->nInt); hTNtag1s.resize(org->nInt); hTNtagm1s.resize(org->nInt);
		hRNtag0s.resize(org->nInt); hRNtag1s.resize(org->nInt); hRNtagm1s.resize(org->nInt);

		// Init and resize ks test result holders
		vector<double> ksEvis1s, ksEvism1s;
		vector<double> ksGtag1s, ksGtagm1s;
		vector<double> ksNtag1s, ksNtagm1s;

		ksEvis1s.resize(org->nInt); ksEvism1s.resize(org->nInt);
		ksGtag1s.resize(org->nInt); ksGtagm1s.resize(org->nInt);
		ksNtag1s.resize(org->nInt); ksNtagm1s.resize(org->nInt);
			
		for(unsigned int i=0; i<org->nInt; i++) { // loop over interactions to plot
			// Make 0 and 1 sigma PDF copy histograms
			// -Energy
			hTEvis0s[i] = (TH1D*) pdfs[i]->evisKdes[nInd].histPdf0->Clone();
			hTEvis0s[i]->SetLineColor(org->intColors[i]);
			hTEvis0s[i]->SetMarkerColor(org->intColors[i]);
			hTEvis0s[i]->SetLineWidth(1);
			hTEvis0s[i]->SetStats(0);
			hTEvis0s[i]->GetXaxis()->SetTitle(sAxisEvis);
			hTEvis0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hTEvis1s[i] = (TH1D*) pdfs[i]->evisKdes[nInd].histPdf->Clone();
			hTEvis1s[i]->SetLineColor(6);
			hTEvis1s[i]->SetMarkerColor(6);
			hTEvis1s[i]->SetLineStyle(7);
			hTEvis1s[i]->SetLineWidth(1);
			hTEvis1s[i]->SetStats(0);
			hTEvis1s[i]->GetXaxis()->SetTitle(sAxisEvis);
			hTEvis1s[i]->GetXaxis()->SetTitleOffset(0.8);
			// -Gtag
			hTGtag0s[i] = (TH1D*) pdfs[i]->gtagKdes[nInd][eInd].histPdf0->Clone();
			hTGtag0s[i]->SetLineColor(org->intColors[i]);
			hTGtag0s[i]->SetMarkerColor(org->intColors[i]);
			hTGtag0s[i]->SetLineWidth(1);
			hTGtag0s[i]->SetStats(0);
			hTGtag0s[i]->GetXaxis()->SetTitle(sAxisGtag);
			hTGtag0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hTGtag1s[i] = (TH1D*) pdfs[i]->gtagKdes[nInd][eInd].histPdf->Clone();
			hTGtag1s[i]->SetLineColor(6);
			hTGtag1s[i]->SetMarkerColor(6);
			hTGtag1s[i]->SetLineStyle(7);
			hTGtag1s[i]->SetLineWidth(1);
			hTGtag1s[i]->SetStats(0);
			hTGtag1s[i]->GetXaxis()->SetTitle(sAxisGtag);
			hTGtag1s[i]->GetXaxis()->SetTitleOffset(0.8);
			// -Ntag
			hTNtag0s[i] = (TH1D*) pdfs[i]->hNeutProb0->Clone();
			hTNtag0s[i]->SetLineColor(org->intColors[i]);
			hTNtag0s[i]->SetMarkerColor(org->intColors[i]);
			hTNtag0s[i]->SetLineWidth(1);
			hTNtag0s[i]->SetStats(0);
			hTNtag0s[i]->GetXaxis()->SetTitle(sAxisNtag);
			hTNtag0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hTNtag1s[i] = (TH1D*) pdfs[i]->hNeutProb->Clone();
			hTNtag1s[i]->SetLineColor(6);
			hTNtag1s[i]->SetMarkerColor(6);
			hTNtag1s[i]->SetLineStyle(7);
			hTNtag1s[i]->SetLineWidth(1);
			hTNtag1s[i]->SetStats(0);
			hTNtag1s[i]->GetXaxis()->SetTitle(sAxisNtag);
			hTNtag1s[i]->GetXaxis()->SetTitleOffset(0.8);
			
			// Make 0 and 1 sigma sampled histograms
			// -Energy
			hREvis0s[i] = (TH1D*) pdfs[i]->evisKdes[nInd].histPdf0->Clone();
			hREvis0s[i]->Reset();
			hREvis0s[i]->SetName((TString)hTEvis0s[i]->GetName()+"_sampled");
			hREvis0s[i]->FillRandom(hTEvis0s[i],pdfs[i]->expEvents0);
			hREvis0s[i]->SetLineColor(org->intColors[i]);
			hREvis0s[i]->SetMarkerColor(org->intColors[i]);
			hREvis0s[i]->SetLineWidth(1);
			hREvis0s[i]->SetStats(0);
			hREvis0s[i]->GetXaxis()->SetTitle(sAxisEvis);
			hREvis0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hREvis1s[i] = (TH1D*) pdfs[i]->evisKdes[nInd].histPdf->Clone();
			hREvis1s[i]->Reset();
			hREvis1s[i]->SetName((TString)hTEvis1s[i]->GetName()+"_sampled");
			hREvis1s[i]->FillRandom(hTEvis1s[i],pdfs[i]->expEvents0);
			hREvis1s[i]->SetLineColor(6);
			hREvis1s[i]->SetMarkerColor(6);
			hREvis1s[i]->SetLineStyle(1);
			hREvis1s[i]->SetLineWidth(1);
			hREvis1s[i]->SetStats(0);
			hREvis1s[i]->GetXaxis()->SetTitle(sAxisEvis);
			hREvis1s[i]->GetXaxis()->SetTitleOffset(0.8);
			// -Gtag
			hRGtag0s[i] = (TH1D*) pdfs[i]->gtagKdes[nInd][eInd].histPdf0->Clone();
			hRGtag0s[i]->Reset();
			hRGtag0s[i]->SetName((TString)hTGtag0s[i]->GetName()+"_sampled");
			hRGtag0s[i]->FillRandom(hTGtag0s[i],pdfs[i]->expEvents0);
			hRGtag0s[i]->SetLineColor(org->intColors[i]);
			hRGtag0s[i]->SetMarkerColor(org->intColors[i]);
			hRGtag0s[i]->SetLineWidth(1);
			hRGtag0s[i]->SetStats(0);
			hRGtag0s[i]->GetXaxis()->SetTitle(sAxisGtag);
			hRGtag0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hRGtag1s[i] = (TH1D*) pdfs[i]->gtagKdes[nInd][eInd].histPdf->Clone();
			hRGtag1s[i]->Reset();
			hRGtag1s[i]->SetName((TString)hTGtag1s[i]->GetName()+"_sampled");
			hRGtag1s[i]->FillRandom(hTGtag1s[i],pdfs[i]->expEvents0);
			hRGtag1s[i]->SetLineColor(6);
			hRGtag1s[i]->SetMarkerColor(6);
			hRGtag1s[i]->SetLineStyle(1);
			hRGtag1s[i]->SetLineWidth(1);
			hRGtag1s[i]->SetStats(0);
			hRGtag1s[i]->GetXaxis()->SetTitle(sAxisGtag);
			hRGtag1s[i]->GetXaxis()->SetTitleOffset(0.8);
			// -Ntag
			hRNtag0s[i] = (TH1D*) pdfs[i]->hNeutProb0->Clone();
			hRNtag0s[i]->Reset();
			hRNtag0s[i]->SetName((TString)hTNtag0s[i]->GetName()+"_sampled");
			hRNtag0s[i]->FillRandom(hTNtag0s[i],pdfs[i]->expEvents0);
			hRNtag0s[i]->SetLineColor(org->intColors[i]);
			hRNtag0s[i]->SetMarkerColor(org->intColors[i]);
			hRNtag0s[i]->SetLineWidth(1);
			hRNtag0s[i]->SetStats(0);
			hRNtag0s[i]->GetXaxis()->SetTitle(sAxisNtag);
			hRNtag0s[i]->GetXaxis()->SetTitleOffset(0.8);
			hRNtag1s[i] = (TH1D*) pdfs[i]->hNeutProb->Clone();
			hRNtag1s[i]->Reset();
			hRNtag1s[i]->SetName((TString)hTNtag1s[i]->GetName()+"_sampled");
			hRNtag1s[i]->FillRandom(hTNtag1s[i],pdfs[i]->expEvents0);
			hRNtag1s[i]->SetLineColor(6);
			hRNtag1s[i]->SetMarkerColor(6);
			hRNtag1s[i]->SetLineStyle(1);
			hRNtag1s[i]->SetLineWidth(1);
			hRNtag1s[i]->SetStats(0);
			hRNtag1s[i]->GetXaxis()->SetTitle(sAxisNtag);
			hRNtag1s[i]->GetXaxis()->SetTitleOffset(0.8);

			// Arrange maximums
			double max;
			max = hTEvis0s[i]->GetMaximum();
			if( hTEvis1s[i]->GetMaximum() > max ) max = hTEvis1s[i]->GetMaximum();
			max *= 1.1;
			hTEvis0s[i]->GetYaxis()->SetRangeUser(0,max);
			hTEvis1s[i]->GetYaxis()->SetRangeUser(0,max);
			max = hREvis0s[i]->GetMaximum();
			if( hREvis1s[i]->GetMaximum() > max ) max = hREvis1s[i]->GetMaximum();
			max *= 1.5;
			hREvis0s[i]->GetYaxis()->SetRangeUser(0,max);
			hREvis1s[i]->GetYaxis()->SetRangeUser(0,max);
			max = hTGtag0s[i]->GetMaximum();
			if( hTGtag1s[i]->GetMaximum() > max ) max = hTGtag1s[i]->GetMaximum();
			max *= 1.1;
			hTGtag0s[i]->GetYaxis()->SetRangeUser(0,max);
			hTGtag1s[i]->GetYaxis()->SetRangeUser(0,max);
			max = hRGtag0s[i]->GetMaximum();
			if( hRGtag1s[i]->GetMaximum() > max ) max = hRGtag1s[i]->GetMaximum();
			max *= 1.5;
			hRGtag0s[i]->GetYaxis()->SetRangeUser(0,max);
			hRGtag1s[i]->GetYaxis()->SetRangeUser(0,max);
			max = hTNtag0s[i]->GetMaximum();
			if( hTNtag1s[i]->GetMaximum() > max ) max = hTNtag1s[i]->GetMaximum();
			max *= 1.1;
			hTNtag0s[i]->GetYaxis()->SetRangeUser(0,max);
			hTNtag1s[i]->GetYaxis()->SetRangeUser(0,max);
			max = hRNtag0s[i]->GetMaximum();
			if( hRNtag1s[i]->GetMaximum() > max ) max = hRNtag1s[i]->GetMaximum();
			max *= 1.5;
			hRNtag0s[i]->GetYaxis()->SetRangeUser(0,max);
			hRNtag1s[i]->GetYaxis()->SetRangeUser(0,max);

			// Do KS Tests For random samples histograms
			//ksEvis1s[i] = hREvis0s[i]->KolmogorovTest(hREvis1s[i]);
			//ksGtag1s[i] = hREvis0s[i]->KolmogorovTest(hRGtag1s[i]);
			//ksNtag1s[i] = hREvis0s[i]->KolmogorovTest(hRNtag1s[i]);

			unsigned int dataSize = 100000;
			//unsigned int dataSize = 100*pdfs[i]->expEvents0;
			vector<double> dEvis0s, dEvis1s; 
			vector<double> dGtag0s, dGtag1s; 
			vector<double> dNtag0s, dNtag1s; 

			for(unsigned int d=0;d<dataSize;d++) {
				dEvis0s.push_back(hTEvis0s[i]->GetRandom());
				dEvis1s.push_back(hTEvis1s[i]->GetRandom());
				dGtag0s.push_back(hTGtag0s[i]->GetRandom());
				dGtag1s.push_back(hTGtag1s[i]->GetRandom());
				dNtag0s.push_back(hTNtag0s[i]->GetRandom());
				dNtag1s.push_back(hTNtag1s[i]->GetRandom());
			}
			sort(dEvis0s.begin(), dEvis0s.end());
			sort(dEvis1s.begin(), dEvis1s.end());
			sort(dGtag0s.begin(), dGtag0s.end());
			sort(dGtag1s.begin(), dGtag1s.end());
			sort(dNtag0s.begin(), dNtag0s.end());
			sort(dNtag1s.begin(), dNtag1s.end());

			ksEvis1s[i] = TMath::KolmogorovTest( 
					dEvis0s.size(), &dEvis0s[0], dEvis1s.size(), &dEvis1s[0], "" );  
			ksGtag1s[i] = TMath::KolmogorovTest( 
					dGtag0s.size(), &dGtag0s[0], dGtag1s.size(), &dGtag1s[0], "" );  
			ksNtag1s[i] = TMath::KolmogorovTest( 
					dNtag0s.size(), &dNtag0s[0], dNtag1s.size(), &dNtag1s[0], "" );  


			// Prepare legends
			TString legNevStr0s = Form("No sys events: %.1f",pdfs[i]->expEvents0);
			TString legNevStr1s = Form("1 #sigma events: %.1f",pdfs[i]->expEvents);
			double percentDiff = 100*(pdfs[i]->expEvents - pdfs[i]->expEvents0)/pdfs[i]->expEvents0;
			TString legNevPerDiff = Form("Percent Diff.: %.1f",percentDiff);
			double percentStatErr = 100*sqrt(pdfs[i]->expEvents0)/pdfs[i]->expEvents0;
			TString legNevPerStatErr = Form("Stat. Error: %.1f",percentStatErr);
			// -Evis
			TLegend* leggyEvis = new TLegend(0.55,0.35,0.9,0.65,"","nbNDC");
			TString legEvisStr0s = "allSysAt0";
			TString legEvisStr1s = org->sysNames[s] + Form("_at%.1f",sysShift);
			TString legEvisKS = Form("KS_result: %.3f",ksEvis1s[i]); 
			leggyEvis->AddEntry(hTEvis0s[i],legEvisStr0s,"l");
			leggyEvis->AddEntry(hTEvis1s[i],legEvisStr1s,"l");
			leggyEvis->AddEntry(hTEvis1s[i],legEvisKS,"l");
			leggyEvis->AddEntry(hTEvis0s[i],legNevStr0s,"l");
			leggyEvis->AddEntry(hTEvis1s[i],legNevStr1s,"l");
			leggyEvis->AddEntry(hTEvis1s[i],legNevPerDiff,"l");
			leggyEvis->AddEntry(hTEvis0s[i],legNevPerStatErr,"l");
			// -Gtag
			TLegend* leggyGtag = new TLegend(0.55,0.6,0.9,0.9,"","nbNDC");
			TString legGtagStr0s = "allSysAt0";
			TString legGtagStr1s = org->sysNames[s] + Form("_at%.1f",sysShift);
			TString legGtagKS = Form("KS_result: %.3f",ksGtag1s[i]); 
			leggyGtag->AddEntry(hTGtag0s[i],legGtagStr0s,"l");
			leggyGtag->AddEntry(hTGtag1s[i],legGtagStr1s,"l");
			leggyGtag->AddEntry(hTGtag1s[i],legGtagKS,"l");
			// -Ntag
			TLegend* leggyNtag = new TLegend(0.55,0.6,0.9,0.9,"","nbNDC");
			TString legNtagStr0s = "allSysAt0";
			TString legNtagStr1s = org->sysNames[s] + Form("_at%.1f",sysShift);
			TString legNtagKS = Form("KS_result: %.3f",ksNtag1s[i]); 
			leggyNtag->AddEntry(hTNtag0s[i],legNtagStr0s,"l");
			leggyNtag->AddEntry(hTNtag1s[i],legNtagStr1s,"l");
			leggyNtag->AddEntry(hTNtag1s[i],legNtagKS,"l");

			// Draw to canvases
			// -Evis
			org->cannyEvisSys->cd( i+1 );
			hTEvis0s[i]->Draw("HIST");
			hTEvis1s[i]->Draw("HIST SAME");
			leggyEvis->Draw();
			org->cannyEvisSys->cd( i+1+4 );
			hREvis0s[i]->Draw("E");
			hREvis1s[i]->Draw("E SAME");
			leggyEvis->Draw();
			// -Gtag
			org->cannyGtagSys->cd( i+1 );
			hTGtag0s[i]->Draw("HIST");
			hTGtag1s[i]->Draw("HIST SAME");
			leggyGtag->Draw();
			org->cannyGtagSys->cd( i+1+4 );
			hRGtag0s[i]->Draw("E");
			hRGtag1s[i]->Draw("E SAME");
			leggyGtag->Draw();
			// -Ntag
			org->cannyNtagSys->cd( i+1 );
			hTNtag0s[i]->Draw("HIST");
			hTNtag1s[i]->Draw("HIST SAME");
			leggyNtag->Draw();
			org->cannyNtagSys->cd( i+1+4 );
			hRNtag0s[i]->Draw("E");
			hRNtag1s[i]->Draw("E SAME");
			leggyNtag->Draw();
			
		}
		org->cannyEvisSys->Print(org->outPdfFileEvisSys);
		org->cannyGtagSys->Print(org->outPdfFileGtagSys);
		org->cannyNtagSys->Print(org->outPdfFileNtagSys);

		// To print into text file
		outSys << org->sysNames[s];
		for(unsigned int i=0; i<org->nInt; i++) { // loop over interactions to plot
			double expNev = pdfs[i]->expEvents0;
			double percentDiff = 100*(pdfs[i]->expEvents - pdfs[i]->expEvents0)/pdfs[i]->expEvents0;
			outSys << " " << expNev << " " << ksEvis1s[i] << " " << ksGtag1s[i] << " " << ksNtag1s[i]
				     << " " << percentDiff;
		}
		outSys << endl;

	} // end of loop over systematics
	outSys.close();

} // end of printSysHists function

// To study some details of systematics effects in given circumstances
void Experiment::checkSysEffect(unsigned int sysIndex, unsigned int pdfIndex) {

	vector<double> inValues;
	inValues = fillLinear( -3.0*org->sysSigmas[sysIndex],
				3.0*org->sysSigmas[sysIndex], 0.1*org->sysSigmas[sysIndex]);
	
	// init all sys to 0
	unsigned int nSysLoop = org->nSys;
	//unsigned int nSysLoop = 10;
	vector<double> inSys;
	inSys.resize(org->nSys);
	for(unsigned int s=0; s<nSysLoop;s++) inSys[s] = 0;
	
	cout << endl << "***** INFO ABOUT SYSTEMATIC: *****" << endl << endl;
	cout << " Experiment: " << org->expNames[expIndex] << endl;
	cout << " Effect of systematic: " << org->sysNames[sysIndex]
		   << " on the PDF: " << org->intNames[pdfIndex] << endl << endl;

	unsigned int nInd = 0;
	vector<TH1D*> evisPdfsNall; evisPdfsNall.resize(inValues.size()); 
	vector<TH1D*> evisRawsNall; evisRawsNall.resize(inValues.size()); 
	TH1D* evisPdfsNall0;
	evisPdfsNall0 = (TH1D*) pdfs[pdfIndex]->evisKdes[nInd].histPdf0->Clone();
	TH1D* evisRawsNall0;
	evisRawsNall0 = (TH1D*) pdfs[pdfIndex]->evisKdes[nInd].histRaw->Clone();
	evisRawsNall0->Reset();
	unsigned int evCount0 = pdfs[pdfIndex]->evisKdes[nInd].evCount;
	for(unsigned int e=0; e<evCount0;e++) {
		evisRawsNall0->Fill(
				pdfs[pdfIndex]->evisKdes[nInd].evX[e], pdfs[pdfIndex]->evisKdes[nInd].evWeight[e]);
	}
	
	for(unsigned int v=0; v<inValues.size(); v++) {
		inSys[sysIndex] = inValues[v];
		org->updateCurSys(inSys);
		updateAllPdfs();
		TString hName = pdfs[pdfIndex]->evisKdes[nInd].histPdf->GetName();
		hName = hName + Form("_sysVal%d",inValues[v]);
		evisPdfsNall[v] = (TH1D*) pdfs[pdfIndex]->evisKdes[nInd].histPdf->Clone(hName);
		
		hName = pdfs[pdfIndex]->evisKdes[nInd].histPdf->GetName();
		hName = hName + Form("_raw_sysVal%d",inValues[v]);
		evisRawsNall[v] = (TH1D*) pdfs[pdfIndex]->evisKdes[nInd].histPdf->Clone(hName);
		evisRawsNall[v]->Reset();
		unsigned int evCount = pdfs[pdfIndex]->evisKdes[nInd].evCount;
		for(unsigned int e=0; e<evCount;e++) {
			evisRawsNall[v]->Fill(
					pdfs[pdfIndex]->evisKdes[nInd].evX[e], pdfs[pdfIndex]->evisKdes[nInd].evWeight[e]);
		}
		
		unsigned int nBins=0;
		double chi2Val = 0;
		double chi2ValRaw = 0;
		double evisBeginDiff = 0; 
		
		nBins = evisPdfsNall0->GetNbinsX(); 
		for(unsigned int b=0;b<nBins;b++) {
			double theoVal = evisPdfsNall0->GetBinContent(b+1);
			double poisVal = evisPdfsNall[v]->GetBinContent(b+1);
			double poisErrSqr = theoVal;
			if(poisVal>0 && theoVal>0) {
				chi2Val += pow((poisVal-theoVal),2)/poisErrSqr; 
			}	
			if(b<5) {
				evisBeginDiff += (poisVal - theoVal); 
			}
			theoVal = evisRawsNall0->GetBinContent(b+1);
			poisVal = evisRawsNall[v]->GetBinContent(b+1);
			poisErrSqr = theoVal;
			if(poisVal>0 && theoVal>0) {
				chi2ValRaw += pow((poisVal-theoVal),2)/poisErrSqr; 
			}	
		}

		double curIntegral = evisPdfsNall[v]->Integral();
		double curNevents = pdfs[pdfIndex]->expEvents;

		cout << " Value of systematic: " << inValues[v] << endl;
		cout << " Integral: " << curIntegral << " Events: " << curNevents << endl;
		cout << " Chi2 To systematic at 0: " << chi2Val << endl;
		cout << " RAW Chi2 To systematic at 0: " << chi2ValRaw << endl;
		cout << " Diff of first 5 energy bins: " << evisBeginDiff << endl; 
			
	
	}
	
	org->outFile->cd();
	evisPdfsNall0->Write();
	delete evisPdfsNall0;
	for(unsigned int v=0; v<inValues.size(); v++) {
		evisPdfsNall[v]->Write();
		evisRawsNall[v]->Scale(1.0/((double)evisRawsNall[v]->Integral()));
		evisRawsNall[v]->Write();
		delete evisPdfsNall[v];
		delete evisRawsNall[v];
	}
	
}
