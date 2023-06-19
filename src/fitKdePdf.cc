/********
*Author: Baran Bodur
*Date: 2021-07-08
*Description: Prepare and store PDFs from input tree with the use of KDEs, systematics are accounted
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
#include "fitKdePdf.h"
#include "organizer.h"
#include "mockData.h"
#include "experiment.h"
#include "vectorUtils.h"

using namespace std;

//Functions go here

// Constructor
FitKdePdf::FitKdePdf(Organizer *orgIn, unsigned int expIndexIn, unsigned int intIndexIn) {
	org = orgIn;
	expIndex = expIndexIn;
	intIndex = intIndexIn;
	pdfName = org->expNames[expIndex] + "_" + org->intNames[intIndex];
	cout << "Creating Interaction Pdf: " << org->intNames[intIndex] << " for Experiment Named: "
			 << org->expNames[expIndex] << endl;
	fillKdeHolders();
	initKdeHolders();
	fillWhichKdeInd();
	updateAllPdfs("forceUpdate");
	saveCurPdfs("withInitSys");
	cout << "For PDF: " << pdfName << " expected events at 0 systematics is: "
			 << expEvents0 << endl;


	//printVector(evisKdes[0].evSysWeight,"All Neutron Visible Energy Sys Weights");
	//printVector(gtagKdes[0][0].evSysWeight,"All Neutron Gamma Tag Sys Weights");
	//if(expIndex==4 && intIndex==0) {
	//	printVector(gtagKdes[0][12].evX0,"All Neutron Gamma Tag Vals");
	//	printVector(gtagKdes[0][12].evWeight0,"All Neutron Gamma Tag Weights");
	//}
}

// Fill Kde Holders From Tree
void FitKdePdf::fillKdeHolders() {
	
	// Read experiment tree and set branches:
	TFile* inFile = new TFile(org->expFiles[expIndex],"READ");
	TTree* inTree = (TTree*) inFile->Get(org->expTrees[expIndex]);
	
	// Tree reading variables (non systematic) and setting of related branches:
	int intName;
	float eVis;
	float tNuEne;
	float tElEne;
	float gTag;
	float tGammaEne;
	int nGamPmt;
	int nGamUnPmt;
	int nTag;
	float weight;
	int nMueTrue;
	int nTagTrue;
	float dWall;
	float openAng3Hit;
	float openAngRatio;
	float openAngPiRatio;
	float zenith;
	float zenithTrue;
	float zenithLepton;
	float azi;
	float aziTrue;
	float aziLepton;

	inTree->SetBranchAddress("intName",&intName);          
	inTree->SetBranchAddress("eVis",&eVis);          
	inTree->SetBranchAddress("tNuEne",&tNuEne);          
	inTree->SetBranchAddress("tElEne",&tElEne);          
	inTree->SetBranchAddress("gTag",&gTag);          
	inTree->SetBranchAddress("tGammaEne",&tGammaEne);          
	inTree->SetBranchAddress("nGamPmt",&nGamPmt);          
	inTree->SetBranchAddress("nGamUnPmt",&nGamUnPmt);          
	inTree->SetBranchAddress("nTag",&nTag);          
	inTree->SetBranchAddress("weight",&weight);          
	inTree->SetBranchAddress("nMueTrue",&nMueTrue);          
	inTree->SetBranchAddress("nTagTrue",&nTagTrue);          
	inTree->SetBranchAddress("dWall",&dWall);          
	inTree->SetBranchAddress("openAng3Hit",&openAng3Hit);          
	inTree->SetBranchAddress("openAngRatio",&openAngRatio);          
	inTree->SetBranchAddress("openAngPiRatio",&openAngPiRatio);          
	inTree->SetBranchAddress("zenith",&zenith);          
	inTree->SetBranchAddress("zenithTrue",&zenithTrue);          
	inTree->SetBranchAddress("zenithLepton",&zenithLepton);          
	inTree->SetBranchAddress("azi",&azi);          
	inTree->SetBranchAddress("aziTrue",&aziTrue);          
	inTree->SetBranchAddress("aziLepton",&aziLepton);          
	
	// Tree read variables (systematic) and setting of related branches
	relScaleSys = getRelSys(nRelScaleSys,"scale","any");
	relShiftSys = getRelSys(nRelShiftSys,"shift","any");
	vector<float> scaleSysWeights(org->nSys,0.0);
	
	for(unsigned int s = 0; s<org->nSys; s++) {
		if( relScaleSys[s]==true ) {
			cout << "Related Scale Sys. Identified: " << org->sysNames[s] << "!" << endl;
			inTree->SetBranchAddress(org->weightNames[s],&scaleSysWeights[s]);
		}
	}
	
	TString name;
	bool evisOrGtag;
	int neutronType;
	int evisSubbin;
	unsigned int nRelScaleSys;
	unsigned int nRelShiftSys;

	// Resize kde holders
	TString neutNames[3] = { "Nall", "N0", "N1" };
	if(org->neutronInfo[expIndex] == false) {
		evisKdes.resize(1); 
		ntagInds.resize(1); 
		evisInds.resize(1); 
		evisKdes[0].name = pdfName + "_eVisPdf_" + neutNames[0];
		//cout << "Will initialize: " << evisKdes[0].name << endl;
		evisKdes[0].evisOrGtag = false;
		evisKdes[0].neutronType = -1;
		evisKdes[0].evisSubbin = -1;
		evisKdes[0].relScaleSys = getRelSys(evisKdes[0].nRelScaleSys,"scale","evis");
		evisKdes[0].relShiftSys = getRelSys(evisKdes[0].nRelShiftSys,"shift","evis");
		gtagKdes.resize(1);
		if(org->gammaInfo[expIndex] == false) {
			gtagKdes[0].resize(0); 
			evisInds[0].resize(0); 
		}
		if(org->gammaInfo[expIndex] == true) {
			gtagKdes[0].resize(org->evisSubbinsForGtag[intIndex]+1); 
			evisInds[0].resize(org->evisSubbinsForGtag[intIndex]+1); 
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				gtagKdes[0][e].name = pdfName + "_gTagPdf_" + neutNames[0] + "_" + getGtagSubbinRange(e);
				//cout << "Will initialize: " << gtagKdes[0][e].name << endl;
				gtagKdes[0][e].evisOrGtag = true;
				gtagKdes[0][e].neutronType = -1;
				gtagKdes[0][e].evisSubbin = e;
				gtagKdes[0][e].relScaleSys = getRelSys(gtagKdes[0][e].nRelScaleSys,"scale","gtag");
				gtagKdes[0][e].relShiftSys = getRelSys(gtagKdes[0][e].nRelShiftSys,"shift","gtag");
			}
		}
	}
	if(org->neutronInfo[expIndex] == true) {
		evisKdes.resize(3); 
		ntagInds.resize(3); 
		gtagKdes.resize(3);
		evisInds.resize(3); 
		for(unsigned int n=0; n<3; n++) {
			evisKdes[n].name = pdfName + "_eVisPdf_" + neutNames[n];
			//cout << "Will initialize: " << evisKdes[n].name << endl;
			evisKdes[n].evisOrGtag = false;
			evisKdes[n].neutronType = -1 + n; // all -1, 0 neutrons 0 , 1 neutron 1
			evisKdes[n].evisSubbin = -1;
			evisKdes[n].relScaleSys = getRelSys(evisKdes[n].nRelScaleSys,"scale","evis");
			evisKdes[n].relShiftSys = getRelSys(evisKdes[n].nRelShiftSys,"shift","evis");
			if(org->gammaInfo[expIndex] == false) {
				gtagKdes[n].resize(0); 
				evisInds[n].resize(0); 
			}
			if(org->gammaInfo[expIndex] == true) {
				gtagKdes[n].resize(org->evisSubbinsForGtag[intIndex]+1); 
				evisInds[n].resize(org->evisSubbinsForGtag[intIndex]+1); 
				for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
					gtagKdes[n][e].name = pdfName + "_gTagPdf_" + neutNames[n] + "_" + getGtagSubbinRange(e);
					//cout << "Will initialize: " << gtagKdes[n][e].name << endl;
					gtagKdes[n][e].evisOrGtag = true;
					gtagKdes[n][e].neutronType = -1 + n;
					gtagKdes[n][e].evisSubbin = e;
					gtagKdes[n][e].relScaleSys = getRelSys(nRelScaleSys,"scale","gtag");
					gtagKdes[n][e].relShiftSys = getRelSys(nRelShiftSys,"shift","gtag");
				}
			}
		}
	}
	
	// For neutron and event count info
	relScaleSysNeut = getRelSys(nRelScaleSysNeut,"scale","ntag");
	relShiftSysNeut = getRelSys(nRelShiftSysNeut,"shift","ntag");
	curSysNeut.resize(org->nSys);
	relScaleSysNev = getRelSys(nRelScaleSysNev,"scale","nevents");
	relShiftSysNev = getRelSys(nRelShiftSysNev,"shift","nevents");
	curSysNev.resize(org->nSys);
	
	
	// Loop over tree
	unsigned int intEvCtr = 0;
	unsigned int treeEvents = inTree->GetEntries();
	cout << "Tree has: " << treeEvents << " total events!" << endl;
	for(unsigned int eCtr = 0; eCtr<treeEvents; eCtr++) {
		inTree->GetEntry(eCtr);
	
		// Experiment and interaction cuts!
		if( !(org->intNameCutLow[intIndex]<=intName && intName<org->intNameCutUp[intIndex]) ) continue;
		if( !(org->nuEneCutLow[intIndex]<=tNuEne && tNuEne<org->nuEneCutUp[intIndex]) ) continue;
		if( !(org->nmueCutLow[intIndex]<=nMueTrue && nMueTrue<org->nmueCutUp[intIndex]) ) continue;
		if( !(org->wallCutLow[intIndex]<=dWall && dWall<org->wallCutUp[intIndex]) ) continue;
		if( !(org->gtagLowBuf[expIndex]<=gTag && gTag<org->gtagUpBuf[expIndex]) ) continue;
		if( !(org->evisLowBuf[expIndex]<=eVis && eVis<org->evisUpBuf[expIndex]) ) continue;


		bool pushGtag = false;
		bool pushEvis = false;
		if( (org->evisLow[expIndex]<=eVis && eVis<org->evisUp[expIndex]) ) pushGtag = true;
		if( (org->gtagLow[expIndex]<=gTag && gTag<org->gtagUp[expIndex]) ) pushEvis = true;

		//more cuts needed
	
		// Fill KDE/Hist Init Vectors From Tree Events
		unsigned int gTagIndex = getEvisSubbinForGtag((double)eVis);
		//cout << "Evis: " << eVis << " gTagIndex: " << gTagIndex << endl;
		unsigned int nTagIndex = 0;	


		vector<double> tempEvisSysVector;
		vector<double> tempGtagSysVector;
		vector<double> tempNtagSysVector;
		vector<double> tempNevSysVector;
		for(unsigned int s = 0; s<org->nSys; s++) {
			// Simple NaN protection, if systematic weight is NaN assign it 1 (default)
			if( scaleSysWeights[s] != scaleSysWeights[s] ) scaleSysWeights[s] = 1.0;

			if( evisKdes[0].relScaleSys[s]==true ) {
				tempEvisSysVector.push_back((double)scaleSysWeights[s]);
			}
			if( org->gammaInfo[expIndex] == true && gtagKdes[0][0].relScaleSys[s] == true) {
				tempGtagSysVector.push_back((double)scaleSysWeights[s]);
			}
			if( relScaleSysNeut[s] == true ) {
				tempNtagSysVector.push_back((double)scaleSysWeights[s]);	
			}
			if( relScaleSysNev[s] == true ) {
				tempNevSysVector.push_back((double)scaleSysWeights[s]);	
			}	
		}

		if(pushEvis == true) {
			evisKdes[nTagIndex].evX0.push_back((double)eVis);
			evisKdes[nTagIndex].evWeight0.push_back((double)weight);
			evisKdes[nTagIndex].evSysWeight.push_back(tempEvisSysVector);
			// Also push ntag and event count vectors here:
			neutEvX0.push_back((double)eVis);
			neutEvWeight0.push_back( (double)( weight * org->reweightFactors[expIndex] ) );
			ntagVec0.push_back((int)nTag);
			neutSysWeights.push_back(tempNtagSysVector);
			nevSysWeights.push_back(tempNevSysVector);
		}
		if(pushGtag == true && org->gammaInfo[expIndex] == true) {
			gtagKdes[nTagIndex][0].evX0.push_back((double)gTag);
			gtagKdes[nTagIndex][0].evWeight0.push_back((double)weight);
			gtagKdes[nTagIndex][0].evSysWeight.push_back(tempGtagSysVector);
			gtagKdes[nTagIndex][gTagIndex].evX0.push_back((double)gTag);
			gtagKdes[nTagIndex][gTagIndex].evWeight0.push_back((double)weight);
			gtagKdes[nTagIndex][gTagIndex].evSysWeight.push_back(tempGtagSysVector);
		}

		if( org->neutronInfo[expIndex]==true ) { // if neutron info exists use it
			if( nTag==0 ) nTagIndex = 1;
			if( nTag==1 ) nTagIndex = 2;
			if(pushEvis == true ) {
				evisKdes[nTagIndex].evX0.push_back((double)eVis);
				evisKdes[nTagIndex].evWeight0.push_back((double)weight);
				evisKdes[nTagIndex].evSysWeight.push_back(tempEvisSysVector);
			}
			if(pushGtag == true && org->gammaInfo[expIndex] == true) {
				gtagKdes[nTagIndex][0].evX0.push_back((double)gTag);
				gtagKdes[nTagIndex][0].evWeight0.push_back((double)weight);
				gtagKdes[nTagIndex][0].evSysWeight.push_back(tempGtagSysVector);
				gtagKdes[nTagIndex][gTagIndex].evX0.push_back((double)gTag);
				gtagKdes[nTagIndex][gTagIndex].evWeight0.push_back((double)weight);
				gtagKdes[nTagIndex][gTagIndex].evSysWeight.push_back(tempGtagSysVector);
			}
		}

	} // end of loop over tree events 
	inFile->Close();
	
	// Initialize vector sizes of vectors that are not filled yet.
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	unsigned int kdeSize;
	for(unsigned int n=0; n<nMax; n++) {
		kdeSize = evisKdes[n].evX0.size();
		evisKdes[n].evCount = kdeSize;	
		evisKdes[n].evX.resize(kdeSize);	
		evisKdes[n].evWeight.resize(kdeSize);	
		evisKdes[n].curSys.resize(org->nSys);
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				kdeSize = gtagKdes[n][e].evX0.size();
				gtagKdes[n][e].evCount = kdeSize;	
				gtagKdes[n][e].evX.resize(kdeSize);	
				gtagKdes[n][e].evWeight.resize(kdeSize);	
				gtagKdes[n][e].curSys.resize(org->nSys);
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	return;
} // end of fillKdeHolders function

// Initialize histograms and kdes after reading vectors from trees
void FitKdePdf::initKdeHolders() {
	
	org->outFile->cd();

	// Temporary variables
	TString tName;
	TString kdeOpt;
	double kdeRho;
	double hMin;
	double hMax;
	double hMinBuf;
	double hMaxBuf;
	double hBinW;
	unsigned int kdeSize;
	unsigned int hBinsNorm;
	unsigned int hColor;
	TString xName;

	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	
	for(unsigned int n=0; n<nMax; n++) {
		// Here set evis kde and histograms
		kdeOpt = org->evisKdeOptStr[intIndex];
		kdeRho = (double) org->evisKdeRho[intIndex];
		hMin = (double) org->evisLow[expIndex];
		hMax = (double) org->evisUp[expIndex];
		hMinBuf = (double) org->evisLowBuf[expIndex];
		hMaxBuf = (double) org->evisUpBuf[expIndex];
		kdeSize = evisKdes[n].evX0.size();
		hBinsNorm = org->evisBinsNorm[expIndex];
		hBinW = (hMax-hMin)/((double)hBinsNorm);
		hColor = org->intColors[intIndex];
		xName = "Visible Energy (MeV)";


		// These are KDE and histogram PDFs assuming no systematics (initialized now)
		tName = "kdePdf0_" + evisKdes[n].name;
		//evisKdes[n].kdePdf0 = new TKDE(kdeSize,&evisKdes[n].evX0[0],&evisKdes[n].evWeight0[0],
		//		hMin,hMax,kdeOpt,kdeRho);	
		evisKdes[n].kdePdf0 = new TKDE(kdeSize,&evisKdes[n].evX0[0],&evisKdes[n].evWeight0[0],
				hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
		
		tName = "histPdf0_" + evisKdes[n].name;
		evisKdes[n].histPdf0 = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
		evisKdes[n].histPdf0->SetLineColor(hColor);
		evisKdes[n].histPdf0->SetLineWidth(2);
		evisKdes[n].histPdf0->GetXaxis()->SetTitle(xName);

		for(unsigned int b=0; b<hBinsNorm; b++) {
			if(kdeSize>1) {
		 		evisKdes[n].histPdf0->SetBinContent(b+1, 
		 				evisKdes[n].kdePdf0->GetValue(evisKdes[n].histPdf0->GetBinCenter(b+1))); 
			}
			else {
				evisKdes[n].histPdf0->SetBinContent(b+1,0);
			}
		 	evisKdes[n].histPdf0->SetBinError(b+1, 0);	
		}
		if(evisKdes[n].histPdf0->Integral()>0) {
			evisKdes[n].histPdf0->Scale(1.0/(evisKdes[n].histPdf0->Integral()*hBinW));
		}
		evisKdes[n].histPdf0->Write();
		cout << "Histogram: " << evisKdes[n].histPdf0->GetName() << " is created!" << endl;

		tName = "histRaw_" + evisKdes[n].name;
		evisKdes[n].histRaw = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
		evisKdes[n].histRaw->SetLineColor(hColor);
		evisKdes[n].histRaw->SetLineWidth(2);
		evisKdes[n].histRaw->GetXaxis()->SetTitle(xName);
		
		evisKdes[n].evTotalWeight0 = 0;
		evisKdes[n].evTotalWeight = 0;
		for(unsigned int ev=0; ev<kdeSize; ev++) {
			evisKdes[n].histRaw->Fill(evisKdes[n].evX0[ev],evisKdes[n].evWeight0[ev]);
			if( org->evisLow[expIndex]<=evisKdes[n].evX0[ev] && 
					evisKdes[n].evX0[ev]<org->evisUp[expIndex] ) {
				evisKdes[n].evTotalWeight0 += evisKdes[n].evWeight0[ev];
			}
		}
		if(evisKdes[n].histRaw->Integral()>0) {
			evisKdes[n].histRaw->Scale(1.0/(evisKdes[n].histRaw->Integral()*hBinW));
		}
		evisKdes[n].histRaw->Write();
		cout << "Histogram: " << evisKdes[n].histRaw->GetName() << " is created!" << endl;
		cout << "Total weight of events for above histogram: " << evisKdes[n].evTotalWeight0 << endl;
		
		// These need to be updated based on systematics
		tName = "kdePdf_" + evisKdes[n].name;
		//evisKdes[n].kdePdf = new TKDE(0,NULL,NULL,hMin,hMax,kdeOpt,kdeRho);	
		evisKdes[n].kdePdf = new TKDE(0,NULL,NULL,hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
		
		tName = "histPdf_" + evisKdes[n].name;
		evisKdes[n].histPdf = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
		evisKdes[n].histPdf->SetLineColor(hColor);
		evisKdes[n].histPdf->SetLineWidth(2);
		evisKdes[n].histPdf->GetXaxis()->SetTitle(xName);

		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				// Here set gtag kde and histograms
				kdeOpt = org->gtagKdeOptStr[intIndex];
				kdeRho = (double) org->gtagKdeRho[intIndex];
				hMin = (double) org->gtagLow[expIndex];
				hMax = (double) org->gtagUp[expIndex];
				hMinBuf = (double) org->gtagLowBuf[expIndex];
				hMaxBuf = (double) org->gtagUpBuf[expIndex];
				kdeSize = gtagKdes[n][e].evX0.size();
				hBinsNorm = org->gtagBinsNorm[expIndex];
				hBinW = (hMax-hMin)/((double)hBinsNorm);
				hColor = org->intColors[intIndex];
				xName = "Gamma Tagging Output";

				// These are KDE and histogram PDFs assuming no systematics (initialized now)
				tName = "kdePdf0_" + gtagKdes[n][e].name;
				gtagKdes[n][e].kdePdf0 = new TKDE(kdeSize,&gtagKdes[n][e].evX0[0],
						&gtagKdes[n][e].evWeight0[0],hMin,hMax,kdeOpt,kdeRho);	
				//gtagKdes[n][e].kdePdf0 = new TKDE(kdeSize,&gtagKdes[n][e].evX0[0],
				//		&gtagKdes[n][e].evWeight0[0],hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
				
				tName = "histPdf0_" + gtagKdes[n][e].name;
				gtagKdes[n][e].histPdf0 = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
				gtagKdes[n][e].histPdf0->SetLineColor(hColor);
				gtagKdes[n][e].histPdf0->SetLineWidth(2);
				gtagKdes[n][e].histPdf0->GetXaxis()->SetTitle(xName);

				for(unsigned int b=0; b<hBinsNorm; b++) {
					if(kdeSize>1) {
						gtagKdes[n][e].histPdf0->SetBinContent(b+1, 
								gtagKdes[n][e].kdePdf0->GetValue(gtagKdes[n][e].histPdf0->GetBinCenter(b+1))); 
					}
					else{
						gtagKdes[n][e].histPdf0->SetBinContent(b+1,0);
					}
					gtagKdes[n][e].histPdf0->SetBinError(b+1, 0);	
				}
				if(gtagKdes[n][e].histPdf0->Integral()>0) {
					gtagKdes[n][e].histPdf0->Scale(1.0/(gtagKdes[n][e].histPdf0->Integral()*hBinW));
				}
				gtagKdes[n][e].histPdf0->Write();
				cout << "Histogram: " << gtagKdes[n][e].histPdf0->GetName() << " is created!" << endl;

				tName = "histRaw_" + gtagKdes[n][e].name;
				gtagKdes[n][e].histRaw = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
				gtagKdes[n][e].histRaw->SetLineColor(hColor);
				gtagKdes[n][e].histRaw->SetLineWidth(2);
				gtagKdes[n][e].histRaw->GetXaxis()->SetTitle(xName);
				
				gtagKdes[n][e].evTotalWeight0 = 0;
				gtagKdes[n][e].evTotalWeight = 0;
				for(unsigned int ev=0; ev<kdeSize; ev++) {
					gtagKdes[n][e].histRaw->Fill(gtagKdes[n][e].evX0[ev],gtagKdes[n][e].evWeight0[ev]);
					if( org->gtagLow[expIndex]<=gtagKdes[n][e].evX0[ev] && 
							gtagKdes[n][e].evX0[ev]<org->gtagUp[expIndex] ) {
						gtagKdes[n][e].evTotalWeight0 += gtagKdes[n][e].evWeight0[ev];
					}
				}
				if(gtagKdes[n][e].histRaw->Integral()>0) {
					gtagKdes[n][e].histRaw->Scale(1.0/(gtagKdes[n][e].histRaw->Integral()*hBinW));
				}
				gtagKdes[n][e].histRaw->Write();
				cout << "Histogram: " << gtagKdes[n][e].histRaw->GetName() << " is created!" << endl;
				cout << "Total weight of events for above histogram: " 
						 << gtagKdes[n][e].evTotalWeight0 << endl;
				
				// These need to be updated based on systematics
				tName = "kdePdf_" + gtagKdes[n][e].name;
				gtagKdes[n][e].kdePdf = new TKDE(0,NULL,NULL,hMin,hMax,kdeOpt,kdeRho);	
				//gtagKdes[n][e].kdePdf = new TKDE(0,NULL,NULL,hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
				
				tName = "histPdf_" + gtagKdes[n][e].name;
				gtagKdes[n][e].histPdf = new TH1D(tName,tName,hBinsNorm,hMin,hMax);
				gtagKdes[n][e].histPdf->SetLineColor(hColor);
				gtagKdes[n][e].histPdf->SetLineWidth(2);
				gtagKdes[n][e].histPdf->GetXaxis()->SetTitle(xName);
				
			} // evis subbin for gamma loop
		} // gamma info if end
	} // neutron loop
	
	// Handle neutron probability, expected events and save neutron prob hist without systematics
	neutEvCount = neutEvX0.size(); 
	TString hNameNeutron = "hNeutProb0_" + pdfName;
	hNeutProb0 = new TH1D(hNameNeutron, hNameNeutron, 2,0,2);
	hNeutProb0->SetLineColor(hColor);
	hNeutProb0->SetLineWidth(2);
	hNeutProb0->GetXaxis()->SetTitle("Tagged Neutrons");
	
	hNameNeutron = "hNeutProb_" + pdfName;
	hNeutProb = new TH1D(hNameNeutron, hNameNeutron, 2,0,2);
	hNeutProb->SetLineColor(hColor);
	hNeutProb->SetLineWidth(2);
	hNeutProb->GetXaxis()->SetTitle("Tagged Neutrons");

	if(org->neutronInfo[expIndex] == false) {
		expEvents0 = 0;
		for(unsigned int ev=0; ev<neutEvCount; ev++) {
			if( org->evisLow[expIndex]<=neutEvX0[ev] && neutEvX0[ev]<org->evisUp[expIndex] ) {
				expEvents0 += neutEvWeight0[ev];
			}	
		}
		neutProb0[0] = 1.0;
		neutProb0[1] = 0.0;
		cout << "For PDF: " << pdfName 
			<< " there is no neutron info requested, so do not trust these neutron probs!" << endl;
	}
	if(org->neutronInfo[expIndex] == true) {
		neutProb0[0] = 0;
		neutProb0[1] = 0;
		for(unsigned int ev=0; ev<neutEvCount; ev++) {
			if( org->evisLow[expIndex]<=neutEvX0[ev] && neutEvX0[ev]<org->evisUp[expIndex] ) {
				if(ntagVec0[ev] == 0 ) neutProb0[0] += neutEvWeight0[ev];
				else if(ntagVec0[ev] == 1) neutProb0[1] += neutEvWeight0[ev];
			}	
		}
		expEvents0 = neutProb0[0] + neutProb0[1];	
		neutProb0[0] /= expEvents0;
		neutProb0[1] /= expEvents0;
	}
	hNeutProb0->SetBinContent(1,neutProb0[0]);
	hNeutProb0->SetBinContent(2,neutProb0[1]);
	hNeutProb0->Print();
	hNeutProb0->Write();

	return;	
} // initKdeHolders function end

// Update all kde holders according to current systematics
void FitKdePdf::updateAllPdfs(TString upOpt) {
	
	int isUpdated = 0;
	// Initialize KDEs
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	for(unsigned int n=0; n<nMax; n++) {
		isUpdated = updateKdeHolder(evisKdes[n],upOpt);
		if(isUpdated>0) {
			//cout << "KDE: " << evisKdes[n].name << " is updated with current systematics!" << endl;  
		}
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				isUpdated = updateKdeHolder(gtagKdes[n][e],upOpt);
				if(isUpdated>0) {
					//cout<< "KDE: "<<gtagKdes[n][e].name<< " is updated with current systematics!" << endl; 
				}
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	// Update ntag and nev
	isUpdated = updateNtagAndNev(upOpt);
	if(isUpdated>0) {
		//cout << "Ntag prob and exp events of PDF: " << pdfName 
		//		 << " is updated with current systematics!" << endl;  
	}

} // end of update all Pdfs with current systematics function

// Save systematic updated kde pdf histograms
void FitKdePdf::saveCurPdfs(TString saveName) {
	TString tempStore;
	org->outFile->cd();
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	for(unsigned int n=0; n<nMax; n++) {
		tempStore = evisKdes[n].histPdf->GetName();
		evisKdes[n].histPdf->SetName(tempStore+"_"+saveName);
		evisKdes[n].histPdf->Write();
		evisKdes[n].histPdf->SetName(tempStore);
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
			tempStore = gtagKdes[n][e].histPdf->GetName();
			gtagKdes[n][e].histPdf->SetName(tempStore+"_"+saveName);
			gtagKdes[n][e].histPdf->Write();
			gtagKdes[n][e].histPdf->SetName(tempStore);
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	tempStore = hNeutProb->GetName();
	hNeutProb->SetName(tempStore+"_"+saveName);
	hNeutProb->Write();
	hNeutProb->SetName(tempStore);

	cout << "Saved PDF hists generated with current values of the systematics!" << endl;
	return;
}

// Update 1 kde holder according to current systematics
int FitKdePdf::updateKdeHolder(KdeHolder &khIn, TString upOpt) {
	// Decide if update is necessary from the last used one
	bool updateKde = false;
 	for(unsigned int s=0; s<org->nSys; s++) {
 		if( (khIn.relScaleSys[s]==true || khIn.relShiftSys[s]==true) 
 				&& khIn.curSys[s]!=org->sysCurValues[s] ) updateKde = true;	
	}
	if(upOpt=="forceUpdate") {
		updateKde = true;
	}
	
	int updateReturn = 0;
	if( updateKde==true ) { // need to loop over if any is updated
		updateReturn = 1;
		khIn.evTotalWeight = 0;
		for(unsigned int ev=0; ev<khIn.evCount; ev++) {	
			khIn.evX[ev] = 1.0;
			khIn.evWeight[ev] = 1.0;
			unsigned int sInd = 0;
			for(unsigned int s=0; s<org->nSys; s++) {
				if( khIn.relScaleSys[s]==true ) {
					if(org->sysCurValues[s]>0) { 
						khIn.evWeight[ev] += org->sysCurValues[s] * org->sysPosCorr[s] *
							(khIn.evSysWeight[ev][sInd]-1.0);
					}
					else { // else if(org->sysCurValues[s]<0) ... 
						khIn.evWeight[ev] += org->sysCurValues[s] * org->sysNegCorr[s] *
							(khIn.evSysWeight[ev][sInd]-1.0);
					}
					sInd++;
					//cout << pdfName << " Ev: " << ev << " Sys: " << org->sysNames[s] << " sysWeight : "
					//	   << khIn.evSysWeight[ev][sInd-1] << " evWeight: " << khIn.evWeight[ev] << endl;
				}
				if( khIn.relShiftSys[s]==true ) {
					if( khIn.evisOrGtag==false ) khIn.evX[ev] += org->sysCurValues[s]/100.0; 
					if( khIn.evisOrGtag==true ) khIn.evX[ev] += org->sysCurValues[s]; // additive shift
				}
			} // end of loop over systematic weights for current event
			if( khIn.evWeight[ev]<org->sysWeightMinScale ) khIn.evWeight[ev] = org->sysWeightMinScale;
			if( khIn.evWeight[ev]>org->sysWeightMaxScale ) khIn.evWeight[ev] = org->sysWeightMaxScale;
			khIn.evWeight[ev] *= khIn.evWeight0[ev];
			if( khIn.evisOrGtag==false ) { // % shift for energy
				if( khIn.evX[ev]<org->sysWeightMinShift ) khIn.evX[ev] = org->sysWeightMinShift;
				if( khIn.evX[ev]>org->sysWeightMaxShift ) khIn.evX[ev] = org->sysWeightMaxShift;
				khIn.evX[ev] *= khIn.evX0[ev];
			}
			if( khIn.evisOrGtag==true ) { // if gtag make shift additive (% does not well around 0)
				khIn.evX[ev] += (khIn.evX0[ev]-1.0);
			}
			if(khIn.evisOrGtag==false) {
				//cout<< pdfName << " Ev: " << ev << endl;
				//cout<< "Ori Loc: " << khIn.evX0[ev] << " Shifted Loc: " << khIn.evX[ev] << endl;
				//cout<< "Ori Wei: " << khIn.evWeight0[ev] << " Scaled Wei: " << khIn.evWeight[ev] << endl;
				if( org->evisLow[expIndex]<=khIn.evX[ev] && khIn.evX[ev]<org->evisUp[expIndex] ) {
					khIn.evTotalWeight += khIn.evWeight[ev];
				}		
			}
			else if(khIn.evisOrGtag==true) {
				if( org->gtagLow[expIndex]<=khIn.evX[ev] && khIn.evX[ev]<org->gtagUp[expIndex] ) {
					khIn.evTotalWeight += khIn.evWeight[ev];
				}			
			}
	
		} // end of loop over events

		// Below block results in no change as expected...
		// New try to make systematics more symmetric
		//for(unsigned int ev=0; ev<khIn.evCount; ev++) {	
			//cout << "evWeightBef: " << khIn.evWeight[ev];
			//khIn.evWeight[ev] *= (khIn.evTotalWeight0/khIn.evTotalWeight);
			//cout << " evWeightAft: " << khIn.evWeight[ev] << endl;
		//}

		// Now update the kde
		if(khIn.kdePdf!=NULL) delete khIn.kdePdf;
		TString tName = "kdePdf_" + khIn.name;
		TString kdeOpt;
		double kdeRho;
		double hMin;
		double hMax;	
		double hMinBuf;
		double hMaxBuf;
		double hBinsNorm;
		double hBinW;
		if(khIn.evisOrGtag==false) {
			// Here set evis kde and histograms
			kdeOpt = org->evisKdeOptStr[intIndex];
			kdeRho = (double) org->evisKdeRho[intIndex];
			hMin = (double) org->evisLow[expIndex];
			hMax = (double) org->evisUp[expIndex];
			hMinBuf = (double) org->evisLowBuf[expIndex];
			hMaxBuf = (double) org->evisUpBuf[expIndex];
			hBinsNorm = org->evisBinsNorm[expIndex];
			hBinW = (hMax-hMin)/((double)hBinsNorm);
			//khIn.kdePdf = new TKDE(khIn.evCount,&khIn.evX[0],&khIn.evWeight[0],hMin,hMax,kdeOpt,kdeRho);	
			khIn.kdePdf = new TKDE(khIn.evCount,&khIn.evX[0],&khIn.evWeight[0],hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
		}
		else if(khIn.evisOrGtag==true) {
			// Here set gtag kde and histograms
			kdeOpt = org->gtagKdeOptStr[intIndex];
			kdeRho = (double) org->gtagKdeRho[intIndex];
			hMin = (double) org->gtagLow[expIndex];
			hMax = (double) org->gtagUp[expIndex];
			hMinBuf = (double) org->gtagLowBuf[expIndex];
			hMaxBuf = (double) org->gtagUpBuf[expIndex];
			hBinsNorm = org->gtagBinsNorm[expIndex];
			hBinW = (hMax-hMin)/((double)hBinsNorm);
			khIn.kdePdf = new TKDE(khIn.evCount,&khIn.evX[0],&khIn.evWeight[0],hMin,hMax,kdeOpt,kdeRho);	
			//khIn.kdePdf = new TKDE(khIn.evCount,&khIn.evX[0],&khIn.evWeight[0],hMinBuf,hMaxBuf,kdeOpt,kdeRho);	
		}

		// Now update histograms with the KDE
		for(unsigned int b=0; b<hBinsNorm; b++) {
			if(khIn.evCount>1) {
				khIn.histPdf->SetBinContent(b+1, 
						khIn.kdePdf->GetValue(khIn.histPdf->GetBinCenter(b+1))); 
			}
			else{
				khIn.histPdf->SetBinContent(b+1,0);
			}
			khIn.histPdf->SetBinError(b+1, 0);	
		}
		if(khIn.histPdf->Integral()>0) {
			khIn.histPdf->Scale(1.0/(khIn.histPdf->Integral()*hBinW));
		}	
	
	} // end of update
		
	// Finally update currently used systematics for the PDF
	for(unsigned int s=0; s<org->nSys; s++) {
		khIn.curSys[s] = org->sysCurValues[s];
	}

	return updateReturn;
} // end of kde holder updating function

// Update ntag prob and expEvent from systematics
int FitKdePdf::updateNtagAndNev(TString upOpt) {
	// Decide if update is necessary from the last used one
	bool updateNeut = false;
	bool updateNev = false;
 	for(unsigned int s=0; s<org->nSys; s++) {
 		if( (relScaleSysNeut[s]==true || relShiftSysNeut[s]==true) 
 				&& curSysNeut[s]!=org->sysCurValues[s] ) updateNeut = true;	
 		if( (relScaleSysNev[s]==true || relShiftSysNev[s]==true) 
 				&& curSysNev[s]!=org->sysCurValues[s] ) updateNev = true;	
	}
	if( org->neutronInfo[expIndex] == false ) updateNeut = false;
	// Check whether an update is forced!
	if(upOpt=="forceUpdate") {
		updateNeut = true;
		updateNev = true;
	}
	
	// If update required (or forced) do it:
	int updateReturn = 0;
	if( updateNev==true || updateNeut==true ) { // need to loop over if any is updated
		updateReturn = 1;
		expEvents = 0;
		neutProb[0] = 0.0;
		neutProb[1] = 0.0;
		if(org->neutronInfo[expIndex] == false) {
			for(unsigned int ev=0; ev<neutEvCount; ev++) {
				double newWeight = 1.0;
				double newX = 1.0;
				unsigned int sInd = 0;
 				for(unsigned int s=0; s<org->nSys; s++) {
					 if( relScaleSysNev[s]==true ) {
						 newWeight += org->sysCurValues[s] * (nevSysWeights[ev][sInd]-1.0);

						 //cout << "Updating Nev! " << org->sysNames[s] << ": " << org->sysCurValues[s]
						 //    << " Sys Weight: " << nevSysWeights[ev][sInd] << " Final W: " 
						 //	    << newWeight << endl;
						 sInd++;
					 }
					 if( relShiftSysNev[s]==true ) {
						 newX += org->sysCurValues[s]/100.0; 
					 }
				} // loop over sys
				if( newWeight<org->sysWeightMinScale ) newWeight = org->sysWeightMinScale;
				if( newWeight>org->sysWeightMaxScale ) newWeight = org->sysWeightMaxScale;
				newWeight *= neutEvWeight0[ev];	 
				if( newX<org->sysWeightMinShift ) newX = org->sysWeightMinShift;
				if( newX>org->sysWeightMaxShift ) newX = org->sysWeightMaxShift;
				newX *= neutEvX0[ev];
				// Add up event weights
				if( org->evisLow[expIndex]<=newX && newX<org->evisUp[expIndex] ) {
					expEvents += newWeight;
				}		
			} // loop over events in all neutron evis dist
			neutProb[0] = 1.0;
			neutProb[1] = 0.0;
		} // end of neutron info if ( no neutron)

		else if(org->neutronInfo[expIndex] == true) {
			double expEventsNev = 0;
			for(unsigned int ev=0; ev<neutEvCount; ev++) {
				double newWeight = 1.0;
				double newX = 1.0;
				unsigned int sInd = 0;
				double newWeightNev = 1.0;
				double newXNev = 1.0;
				unsigned int sIndNev = 0;
 				for(unsigned int s=0; s<org->nSys; s++) {
					 if( relScaleSysNeut[s]==true ) {
						 newWeight += org->sysCurValues[s] * (neutSysWeights[ev][sInd]-1.0);
						 sInd++;
					 }
					 if( relShiftSysNeut[s]==true ) {
						 newX += org->sysCurValues[s]/100.0; 
					 }
					 if( relScaleSysNev[s]==true ) {
					 	 // 20221010 -Baran fixed sInd to sIndNev
						 newWeightNev += org->sysCurValues[s] * (nevSysWeights[ev][sIndNev]-1.0);
						 sIndNev++;
					 }
					 if( relShiftSysNev[s]==true ) {
						 newXNev += org->sysCurValues[s]/100.0; 
					 }
				} // loop over sys
				if( newWeight<org->sysWeightMinScale ) newWeight = org->sysWeightMinScale;
				if( newWeight>org->sysWeightMaxScale ) newWeight = org->sysWeightMaxScale;
				newWeight *= neutEvWeight0[ev];	 
				if( newX<org->sysWeightMinShift ) newX = org->sysWeightMinShift;
				if( newX>org->sysWeightMaxShift ) newX = org->sysWeightMaxShift;
				newX *= neutEvX0[ev];
				if( newWeightNev<org->sysWeightMinScale ) newWeightNev = org->sysWeightMinScale;
				if( newWeightNev>org->sysWeightMaxScale ) newWeightNev = org->sysWeightMaxScale;
				newWeightNev *= neutEvWeight0[ev];	 
				if( newXNev<org->sysWeightMinShift ) newXNev = org->sysWeightMinShift;
				if( newXNev>org->sysWeightMaxShift ) newXNev = org->sysWeightMaxShift;
				newXNev *= neutEvX0[ev];
				// Add up event weights
				if( org->evisLow[expIndex]<=newX && newX<org->evisUp[expIndex] ) {
					if(ntagVec0[ev] == 0 ) neutProb[0] += newWeight;
					else if(ntagVec0[ev] == 1) neutProb[1] += newWeight;
				}		
				if( org->evisLow[expIndex]<=newXNev && newXNev<org->evisUp[expIndex] ) {
					if(ntagVec0[ev] == 0 || ntagVec0[ev] == 1) expEventsNev += newWeightNev;
				}		
			} // loop over events in neutron dist
			expEvents = neutProb[0] + neutProb[1];	
			neutProb[0] /= expEvents;
			neutProb[1] /= expEvents;
			expEvents = expEventsNev;
		} // end of neutron info if ( yes neutron )
		
		// Simulated # of events:
		simEvents = expEvents * org->intVarFraction[intIndex]; 
	
		hNeutProb->SetBinContent(1,neutProb[0]);
		hNeutProb->SetBinContent(2,neutProb[1]);

		// Finally update currently used systematics for the PDF
		for(unsigned int s=0; s<org->nSys; s++) {
			curSysNeut[s] = org->sysCurValues[s];
			curSysNev[s] = org->sysCurValues[s];
		}
	
	} // end of update

	return updateReturn;
} // end of update neutron and nevents function

// Determine which systematics to this interaction and experiment, 
// Arguments are scaleOrshift (any,scale,shift) and sysType (any, ntag, gtag, evis, nevents)
vector<bool> FitKdePdf::getRelSys(unsigned int &nRelSys, TString scaleOrShift, TString sysType) {
 vector<bool> outVector(org->nSys,false);
 nRelSys = 0;
 for(unsigned int s=0; s<org->nSys; s++) {
 	 bool temp = org->effExp[s][expIndex] && org->effInt[s][intIndex];
	if( scaleOrShift=="any" && sysType=="any")
 		outVector[s] = temp;
	else if( scaleOrShift=="shift" && sysType=="any")
 		outVector[s] = temp && !org->shiftOrScale[s];
	else if( scaleOrShift=="scale" && sysType=="any")
 		outVector[s] = temp && org->shiftOrScale[s];
	else if( scaleOrShift=="any" && sysType=="ntag")
 		outVector[s] = temp && org->effNeutron[s];
	else if( scaleOrShift=="any" && sysType=="gtag")
 		outVector[s] = temp && org->effGamma[s];
	else if( scaleOrShift=="any" && sysType=="evis")
 		outVector[s] = temp && org->effEnergy[s];
	else if( scaleOrShift=="any" && sysType=="nevents")
 		outVector[s] = temp && org->effNevents[s];
	else if( scaleOrShift=="shift" && sysType=="ntag")
 		outVector[s] = temp && !org->shiftOrScale[s] && org->effNeutron[s];
	else if( scaleOrShift=="shift" && sysType=="gtag")
 		outVector[s] = temp && !org->shiftOrScale[s] && org->effGamma[s];
	else if( scaleOrShift=="shift" && sysType=="evis")
 		outVector[s] = temp && !org->shiftOrScale[s] && org->effEnergy[s];
	else if( scaleOrShift=="shift" && sysType=="nevents")
 		outVector[s] = temp && !org->shiftOrScale[s] && org->effNevents[s];
	else if( scaleOrShift=="scale" && sysType=="ntag")
 		outVector[s] = temp && org->shiftOrScale[s] && org->effNeutron[s];
	else if( scaleOrShift=="scale" && sysType=="gtag")
 		outVector[s] = temp && org->shiftOrScale[s] && org->effGamma[s];
	else if( scaleOrShift=="scale" && sysType=="evis")
 		outVector[s] = temp && org->shiftOrScale[s] && org->effEnergy[s];
	else if( scaleOrShift=="scale" && sysType=="nevents")
 		outVector[s] = temp && org->shiftOrScale[s] && org->effNevents[s];
 	else {
 		cout << "Unavailable option passed to FitKdePdf::getRelSys! Either: " << endl;
 		cout << "scaleOrShift: " << scaleOrShift << " is an unavailable option or" << endl;
 		cout << "sysType: " << sysType << " is an unavailable option!" << endl;
 		cout << "Returning all false boolean vector!";
 	}
 	if(outVector[s] == true) nRelSys++;
 }
 return outVector;
}

// Determine which gtag KDE from evis subbins
unsigned int FitKdePdf::getEvisSubbinForGtag(double evis) {
	double eRange = org->evisUp[expIndex] - org->evisLow[expIndex];
	double eLow = org->evisLow[expIndex];
	unsigned int nSubbins = org->evisSubbinsForGtag[intIndex];
	double eBinW = eRange/(double)nSubbins;
	unsigned int curSubbin = 0;
	for(unsigned int ns=0; ns<nSubbins; ns++) {
		if( evis >= eLow + eBinW*((double)ns) ) {
			curSubbin = ns;
		}
	}
	return curSubbin + 1; // 0th element is gtag projection over all evis, so return +1
}

// Return a descriptive text corresponding to evis subbin for gtag (ie: evis40-50 )
TString FitKdePdf::getGtagSubbinRange(unsigned int subbin) {
	double eRange = org->evisUp[expIndex] - org->evisLow[expIndex];
	double eLow = org->evisLow[expIndex];
	unsigned int nSubbins = org->evisSubbinsForGtag[intIndex];
	double eBinW = eRange/(double)nSubbins;
	double ns = (double) subbin;
	TString range;
	if(subbin==0) range = Form("evis%.0f-%.0f",eLow,eRange+eLow);
	else range = Form("evis%.0f-%.0f",eLow+eBinW*(ns-1),eLow+eBinW*ns);
	return range;
}

// For KDEs without enough events determine which more general PDF to use
void FitKdePdf::fillWhichKdeInd() {
	unsigned int nMax = 0;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;

	for(unsigned int n=0; n<nMax; n++) {
		ntagInds[n] = n;
		if(evisKdes[n].evCount < org->evisKdeLow[intIndex]) ntagInds[n] = 0;
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				evisInds[n][e] = e;
				if(gtagKdes[n][e].evCount < org->gtagKdeLow[intIndex]) evisInds[n][e]=0;
			}
		}
	}
}

// For getting projection histograms from the most recent KDE histogram
void FitKdePdf::updateProjHists(double scale) {

	TString tName;
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;	
	if(org->neutronInfo[expIndex] == true) nMax = 3;

	for(unsigned int n=0; n<nMax; n++) {	
		double scaleTrueEne = scale;
		if(evisKdes[n].histProj!=NULL) delete evisKdes[n].histProj;
		tName = "histProj_" + evisKdes[n].name;
		evisKdes[n].histProj = (TH1D*) evisKdes[n].histPdf->Clone(tName);
		evisKdes[n].histProj->SetTitle(tName);
		evisKdes[n].histProj->Rebin(org->evisProjRebinFactor[expIndex]);

		if(n==1) scaleTrueEne *= neutProb[0]; // If projection of only ntag=0
		if(n==2) scaleTrueEne *= neutProb[1];
		double projBinW = ( (double) org->evisUp[expIndex] - (double) org->evisLow[expIndex] ) / 
			( (double) org->evisBinsNorm[expIndex] );
		scaleTrueEne *= projBinW;
		evisKdes[n].histProj->Scale(scaleTrueEne);
		
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
				double scaleTrueGamma = scale;
				if(gtagKdes[n][e].histProj!=NULL) delete gtagKdes[n][e].histProj;
				tName = "histProj_" + gtagKdes[n][e].name;
				gtagKdes[n][e].histProj = (TH1D*) gtagKdes[n][e].histPdf->Clone(tName);
				gtagKdes[n][e].histProj->SetTitle(tName);
				gtagKdes[n][e].histProj->Rebin(org->gtagProjRebinFactor[expIndex]);
				if(n==1) scaleTrueGamma *= neutProb[0]; // If projection of only ntag=0
				if(n==2) scaleTrueGamma *= neutProb[1];
				double projBinW = ( (double) org->gtagUp[expIndex] - (double) org->gtagLow[expIndex] ) / 
					( (double) org->gtagBinsNorm[expIndex] );
				scaleTrueGamma *= projBinW;
				gtagKdes[n][e].histProj->Scale(scaleTrueGamma);
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	double scaleTrueNeut = scale;
	//if(hNeutProj!=NULL) delete hNeutProj;
	tName = "histProj_" + pdfName + "_nTag";
	hNeutProj = (TH1D*) hNeutProb->Clone(tName);
	hNeutProj->SetTitle(tName);
	hNeutProj->Scale(scaleTrueNeut);

	cout << " Updated projection histograms based on the most recent KDE histograms" << endl;
	return;
}


// For saving projection histograms from the most recent KDE histogram
void FitKdePdf::saveCurProj(TString saveName) {
	TString tempStore;
	org->outFile->cd();
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	for(unsigned int n=0; n<nMax; n++) {
		tempStore = evisKdes[n].histProj->GetName();
		evisKdes[n].histProj->SetName(tempStore+"_"+saveName);
		evisKdes[n].histProj->Write();
		evisKdes[n].histProj->SetName(tempStore);
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[intIndex]+1; e++) {
			tempStore = gtagKdes[n][e].histProj->GetName();
			gtagKdes[n][e].histProj->SetName(tempStore+"_"+saveName);
			gtagKdes[n][e].histProj->Write();
			gtagKdes[n][e].histProj->SetName(tempStore);
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	tempStore = hNeutProj->GetName();
	hNeutProj->SetName(tempStore+"_"+saveName);
	hNeutProj->Write();
	hNeutProj->SetName(tempStore);

	cout << "Saved Projection hists generated with current values of the systematics!" << endl;
	return;
}



// !!*** Functions for obtaining probability and random numbers of this PDF

// Most fundamental PDF prob function: P(n)P(E|n)P(g|n,E)
double FitKdePdf::getProb(double eVis, double gTag, int nTag) {
	double prob = 0; // return variable

	// Make sure inputs are within boundaries: 
	if(eVis<org->evisLow[expIndex]) eVis = org->evisLow[expIndex];
	else if(eVis>org->evisUp[expIndex]) eVis = org->evisUp[expIndex];
	if(gTag<org->gtagLow[expIndex]) gTag = org->gtagLow[expIndex];
	else if(gTag>org->gtagUp[expIndex]) gTag = org->gtagUp[expIndex];
	if(nTag<0) nTag = 0;
	else if(nTag>1) nTag = 1;
	
	// Calculate probability (depends on which info is available and which not...)
	if(org->neutronInfo[expIndex]==true) { // if neutron info avaiable
		prob = neutProb[nTag];
		double tInterp = evisKdes[ntagInds[nTag+1]].histPdf->Interpolate(eVis);
		if(tInterp<=0) tInterp = evisKdes[ntagInds[nTag+1]].histPdf->GetMinimum(0);
		prob *= tInterp; 
		if(org->gammaInfo[expIndex]==true) { // if gamma info available
			unsigned int gtagInd = evisInds[ntagInds[nTag+1]][getEvisSubbinForGtag(eVis)];
			tInterp = gtagKdes[ntagInds[nTag+1]][gtagInd].histPdf->Interpolate(gTag);
			if(tInterp<=0) tInterp = gtagKdes[ntagInds[nTag+1]][gtagInd].histPdf->GetMinimum(0);
			prob *= tInterp;	
		}
	} // end of neutron if
	else if(org->neutronInfo[expIndex]==false) { // no neutron info 
		double tInterp = evisKdes[0].histPdf->Interpolate(eVis);
		if(tInterp<=0) tInterp = evisKdes[0].histPdf->GetMinimum(0);
		prob = tInterp; 	
		if(org->gammaInfo[expIndex]==true) { // if gamma info available
			unsigned int gtagInd = evisInds[0][getEvisSubbinForGtag(eVis)];
			tInterp = gtagKdes[0][gtagInd].histPdf->Interpolate(gTag);
			if(tInterp<=0) tInterp = gtagKdes[0][gtagInd].histPdf->GetMinimum(0);
			prob *= tInterp;		
		}
	} // end of neutron if

	if( prob!=prob || prob<=0) {
		cout << "WARNING!! Probability is NaN or negative, prob: " << prob << endl;
	}
	return prob;
}

// Probability with direction info (not yet implemented, though tested before in the old
// fitter and found not helpful... )
double FitKdePdf::getProb(double eVis, double gTag, int nTag, double cosZen, double azi) {
	if(org->dirInfo[expIndex] == false) {
		return getProb(eVis, gTag, nTag);
	}
	else if(org->dirInfo[expIndex] == true) {
		cout << "*** !!! Direction info is not yet implemented with systematics!" << endl;
		cout << "*** !!! Returning probability without direction info:" << endl; 
		return getProb(eVis, gTag, nTag);
	}
}

// Get probability (direction info not yet implemented so just use evis,gtag,ntag 
double FitKdePdf::getProb(DataPoint inPoint) {
	return getProb(inPoint.eVis, inPoint.gTag, inPoint.nTag);
} 

// Get random ntag:
int FitKdePdf::getRandNtag() {
	if(org->neutronInfo[expIndex]==false) return 0; // nothing to do
	else if(org->neutronInfo[expIndex]==true) { // if neutron info is available
		double tRand = gRandom->Rndm();
		if(tRand < neutProb[0]) return 0;
		else return 1;
	}
}

// Get random energy at given ntag:
double FitKdePdf::getRandEneAtNtag(int nTag) {
	double eT = 0;
	if(org->neutronInfo[expIndex]==false) eT = evisKdes[0].histPdf->GetRandom(); 
	else if(org->neutronInfo[expIndex]==true) { // if neutron info is available
		unsigned int ntagInd = ntagInds[nTag+1];
		eT = evisKdes[ntagInd].histPdf->GetRandom();
	}
	return eT;	
}

// Get random gtag at ntag and ene
double FitKdePdf::getRandGtagAtEneAndNtag(double eVis, int nTag) { 
	double gT = 0;
	unsigned int ntagInd = 0;
	if(org->neutronInfo[expIndex]==true) ntagInd = ntagInds[nTag+1];
	if(org->gammaInfo[expIndex]==true) {
		unsigned int gtagInd = evisInds[ntagInd][getEvisSubbinForGtag(eVis)];
		gT = gtagKdes[ntagInd][gtagInd].histPdf->GetRandom();	
	}
	return gT;
}

// get a random data point from the pdf using more specicif functions 
// Direction is always 0 for now
DataPoint FitKdePdf::getRandPoint() {
	DataPoint rPoint;
	rPoint.nTag = getRandNtag();
	rPoint.eVis = getRandEneAtNtag(rPoint.nTag);
	rPoint.gTag = getRandGtagAtEneAndNtag(rPoint.eVis, rPoint.nTag);
	rPoint.weight = 1.0;
	rPoint.cosZen = 0;
	rPoint.azi = 0;
	return rPoint;
}


