/********
*Author: Baran Bodur
*Date: 2022-07-26
*Description: Generate mock data based on exps/PDF to be used in fitting studies
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
#include "fitKdePdf.h"
#include "mockData.h"

using namespace std;

// Simple function to print data points
void printDataPoint(DataPoint dP) {
	cout << "Visible Energy: " << dP.eVis << endl;
	cout << "Gamma Taggin Output: " << dP.gTag << endl;
	cout << "Neutrons Tagged: " << dP.nTag << endl;
	cout << "Cos(Zenith) Angle: " << dP.cosZen << endl;
	cout << "Azimuth Angle: " << dP.azi << endl;
	cout << "Weight: " << dP.weight << endl;
}

// Class functions go here
MockData::MockData(Organizer *orgIn, unsigned int expIndexIn, unsigned int mockIndexIn) {
	org = orgIn;
	expIndex = expIndexIn;
	mockIndex = mockIndexIn;
	mockName = "Mock_" + org->expNames[expIndex] + Form("_%d", mockIndex);	

	//Resize vectors
	simTrueEvents.resize(org->nInt);
	simPoisEvents.resize(org->nInt);
	simExpEvents.resize(org->nInt);
	simSystematics.resize(org->nSys);

	// Initialize mock data histogram vector sizes
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	evisProjs.resize(nMax);
	gtagProjs.resize(nMax);
	for(unsigned int n=0; n<nMax; n++) {
		if(org->gammaInfo[expIndex] == false) gtagProjs[n].resize(0); 
		if(org->gammaInfo[expIndex] == true) gtagProjs[n].resize(org->evisSubbinsForGtag[0]+1); 
	}
	
/*	
	TString h2Name = "hSysVsChiProb" + mockName;
	hSysVsChiProb = new TH2D(h2Name,h2Name,100,-5,5,100,0,10);
	hSysVsChiProb->GetXaxis()->SetTitle(org->sysNames[0]);
	hSysVsChiProb->GetYaxis()->SetTitle("Chi2 Prob");
*/

} // end of constructor

// Generate mock data (fills DataPoint vectors)
unsigned int MockData::genMockData(vector<double> usedSysVecIn) {

	cout << "Generating Mock Data: " << mockName << endl;
	cout << "Gen Rand Data Option: " << org->genRandData << endl;
	if( org->genRandData==0 ) { // use read data:
		nPoints = readData();	
		totalWeight = nPoints;
		// Fill true events as expected events: True events will not be known
		org->updateCurSys(usedSysVecIn);
		simSystematics = org->sysCurValues;
		for(unsigned int p=0; p<org->nInt; p++) {
			org->exps[expIndex]->pdfs[p]->updateAllPdfs();
			org->exps[expIndex]->pdfs[p]->updateProjHists(1.0);
			simExpEvents[p] = round(org->exps[expIndex]->pdfs[p]->expEvents);
			simTrueEvents[p] = simExpEvents[p]; 
			simPoisEvents[p] = simExpEvents[p]; 
		}
	} // end of read data
	else if(org->genRandData>0) { // random data
		// First set systematics and update PDFs accordingly:
		org->updateCurSys(usedSysVecIn);
		simSystematics = org->sysCurValues;
		for(unsigned int p=0; p<org->nInt; p++) {
			org->exps[expIndex]->pdfs[p]->updateAllPdfs();
			org->exps[expIndex]->pdfs[p]->updateProjHists(1.0);
		}

		// Decide how many events from which pdf, based on simEvents of each pdf 
		nPoints = 0;
		totalWeight = 0;
		byTypeEvents.resize(org->nInt);
		for(unsigned int p=0; p<org->nInt; p++) {
			if(org->nEventsPoisson==false) {
				byTypeEvents[p] = round(org->exps[expIndex]->pdfs[p]->simEvents); 
			}
			else if(org->nEventsPoisson==true) {
				byTypeEvents[p] = gRandom->Poisson(org->exps[expIndex]->pdfs[p]->simEvents); 
			}
			// For book keeping later	
			simTrueEvents[p] = round(org->exps[expIndex]->pdfs[p]->simEvents);
			simExpEvents[p] = round(org->exps[expIndex]->pdfs[p]->expEvents);
			simPoisEvents[p] = (double) byTypeEvents[p];
			
			totalWeight += byTypeEvents[p];
			nPoints += byTypeEvents[p];
			cout << org->exps[expIndex]->pdfs[p]->pdfName << " Exp Ctr: " 
					<< org->exps[expIndex]->pdfs[p]->simEvents << " Gen: " << byTypeEvents[p] << endl;
		} // end of loop over interaction pdfs

		// Create data points according to pdfs
		if(org->genRandData==1) { // random data
			for(unsigned int p=0; p<org->nInt; p++) {
				for(unsigned int d=0; d<byTypeEvents[p]; d++) {
					//cout << "Pdf: " << org->exps[expIndex]->pdfs[p]->pdfName 
					//		 << " Event: " << d+1 << "/" << byTypeEvents[p] << endl;
					DataPoint temp = org->exps[expIndex]->pdfs[p]->getRandPoint(); 	
					dataPoint.push_back(temp);			
				} // end of loop over # of events in a pdf
			} // end of loop over pdfs
		}
		else if(org->genRandData==2) { // generate asimovish data instead of random
			for(unsigned int p=0; p<org->nInt; p++) {
				vector<DataPoint> temp = makeAsimovDataset(p, byTypeEvents[p]);
				for(unsigned int d=0; d<byTypeEvents[p]; d++) {
					dataPoint.push_back(temp[d]);
				}
			}
		}
	} // end of if randData option is>0

	// Make Projection Histograms
	makeProjHists();

	// -- This is under construction for now
  /*
	// --- Begin of just checking data generation
	// Get Ideal Evis Histogram
	vector<double> usedIntVecIn;
	for(unsigned int p=0;p<org->nInt;p++) {
		usedIntVecIn.push_back(1.0);
	}
	if(org->neutronInfo[expIndex]==false) {
		hIdealEvis = org->exps[expIndex]->makeTotalHist(usedIntVecIn, usedSysVecIn, 
			"eVis", 0, 0, "hIdealEvis");
	}
	else {
		hIdealEvis = org->exps[expIndex]->makeTotalHist(usedIntVecIn, usedSysVecIn, 
			"eVis", 1, 0, "hIdealEvis");
		TH1D* hIdealEvisN1 = org->exps[expIndex]->makeTotalHist(usedIntVecIn, usedSysVecIn, 
			"eVis", 2, 0, "hIdealEvis");
		hIdealEvis->Add(hIdealEvisN1);
		delete hIdealEvisN1;
	}
	hIdealEvis->Print();
	evisProjs[0]->Print();
	//double ksTestRes = evisProjs[0]->KolmogorovTest(hIdealEvis,"UW");
	double ksTestRes = 0;
	double totalDiff = 0;
	unsigned int maxBin=5;
	for(unsigned int b=0;b<maxBin;b++) {
		double theoVal = hIdealEvis->GetBinContent(b+1);
		double poisVal = evisProjs[0]->GetBinContent(b+1);
		double poisErrSqr = poisVal;
		ksTestRes += pow((poisVal-theoVal),2)/poisErrSqr; 
		totalDiff += (poisVal-theoVal);
	}
	ksTestRes/=(double)maxBin;
	cout << org->sysNames[0] << ": " << usedSysVecIn[0] << endl; 
	cout << "KS Test To Ideal: " << ksTestRes << " Total Diff: " << totalDiff << endl;
	hSysVsChiProb->Fill(usedSysVecIn[0], ksTestRes);
	// end of checking data generation
	*/	
	return nPoints;
} // end of gen mock data

// Reweight an already generated dataset (original dataset has weights 1.0 exact)
// This is done to compare two different systematic configuration with same random sample
// Assume fi = 1.0 when generating data
void MockData::reweightMockData(vector<double> oriSysVec, vector<double> newSysVec) {

	vector<double> top; top.resize(nPoints);
	vector<double> bottom; bottom.resize(nPoints);
	// Calculate semi-PDF at old point (bottom)
	org->updateCurSys(oriSysVec);	
	for(unsigned int p=0; p<org->nInt; p++) {
		org->exps[expIndex]->pdfs[p]->updateAllPdfs();
		org->exps[expIndex]->pdfs[p]->updateProjHists(1.0);
	}
	for(unsigned int d=0; d<nPoints; d++) {
		//cout << "Data Point: " << d << endl; 
		//printDataPoint(dataPoint[d]);
		bottom[d] = 0;
		for(unsigned int p=0; p<org->nInt; p++) {
			bottom[d] += ( org->exps[expIndex]->pdfs[p]->expEvents * 
					           org->exps[expIndex]->pdfs[p]->getProb(dataPoint[d]) );
			//cout << "PDF: " << org->intNames[p] << " Prob: " << 
			//	org->exps[expIndex]->pdfs[p]->getProb(dataPoint[d]) << endl;
		}
	}
	// Calculate semi-PDF at new point (top)
	org->updateCurSys(newSysVec);	
	for(unsigned int p=0; p<org->nInt; p++) {
		org->exps[expIndex]->pdfs[p]->updateAllPdfs();
		org->exps[expIndex]->pdfs[p]->updateProjHists(1.0);
	}
	for(unsigned int d=0; d<nPoints; d++) {
		top[d] = 0;
		for(unsigned int p=0; p<org->nInt; p++) {
			top[d] += ( org->exps[expIndex]->pdfs[p]->expEvents * 
					        org->exps[expIndex]->pdfs[p]->getProb(dataPoint[d]) );
		}
	}
	// Apply reweighting
	totalWeight = 0;
	for(unsigned int d=0; d<nPoints; d++) {
		dataPoint[d].weight = top[d]/bottom[d];
		if( bottom[d]==0 || dataPoint[d].weight != dataPoint[d].weight ) dataPoint[d].weight = 1.0;
		totalWeight += dataPoint[d].weight;
	}

} // end of mock data reweight

// Reweights 1 systematic with given index and sigma to all systematics at 0 case
// Assumes f_i are 1.0
void MockData::reweightOneSys(unsigned int sysIndex, double nSigma) {

	// Make appropriate systematic vectors
	vector<double> zeroSys; zeroSys.resize(org->nSys);
	vector<double> oneSys; oneSys.resize(org->nSys);
	for(unsigned int s=0; s<org->nSys; s++) {
		zeroSys[s] = 0;
		oneSys[s] = 0;
	}
	oneSys[sysIndex] = nSigma * org->sysSigmas[sysIndex];

	// Call generic reweightMockData
	reweightMockData(zeroSys, oneSys);

	// Print what happened with reweighing
	cout << "Reweighting systematic: " << org->sysNames[sysIndex] << " from 0 to 1 sigma!" << endl;
	cout << "This is mock data for experiment: " << org->expNames[expIndex] << endl;
	cout << "Original weights were all 1.0, new weights are: " << endl;
	for(unsigned int d=0; d<nPoints; d++) {
		cout << " DP: " << d << " Weight: " << dataPoint[d].weight << endl;
	}

	//Temporary, comment later
	//makeProjHists();

} // end of mock data reweight of 1 systematic

// Get Asimov Dataset as a vector of DataPoints for a given interaction / # of events
vector<DataPoint> MockData::makeAsimovDataset(unsigned int intIndexIn, unsigned int nEvents) {
	// Determine # of 0 and 1 neutron events:
	unsigned int nNeutron0, nNeutron1;
	vector<double> evisN0, evisN1;
	vector<double> gtagN0, gtagN1;
	vector<double> gtagN0_t, gtagN1_t;
	if( org->neutronInfo[expIndex] == false ) {
		nNeutron0 = nEvents;
		nNeutron1 = 0;
		//get evis and gtag samples
		evisN0 = getAsimovSample(org->exps[expIndex]->pdfs[intIndexIn]->evisKdes[0].histPdf, 
				nNeutron0);
		gtagN0_t = getAsimovSample(org->exps[expIndex]->pdfs[intIndexIn]->gtagKdes[0][0].histPdf, 
				nNeutron0);
		gtagN0 =  shuffleOrder(gtagN0_t);

	
	}
	else {
		nNeutron0 = round( org->exps[expIndex]->pdfs[intIndexIn]->neutProb[0]*(double)nEvents );
		nNeutron1 = nEvents - nNeutron0;
		// get evis and gtag samples
		evisN0 = getAsimovSample( org->exps[expIndex]->pdfs[intIndexIn]->evisKdes[ 
				org->exps[expIndex]->pdfs[intIndexIn]->ntagInds[1] ].histPdf , nNeutron0);  
		evisN1 = getAsimovSample( org->exps[expIndex]->pdfs[intIndexIn]->evisKdes[ 
				org->exps[expIndex]->pdfs[intIndexIn]->ntagInds[2] ].histPdf , nNeutron1);  
		gtagN0_t = getAsimovSample( org->exps[expIndex]->pdfs[intIndexIn]->gtagKdes[ 
				org->exps[expIndex]->pdfs[intIndexIn]->ntagInds[1] ][0].histPdf , nNeutron0);  
		gtagN1_t = getAsimovSample( org->exps[expIndex]->pdfs[intIndexIn]->gtagKdes[ 
				org->exps[expIndex]->pdfs[intIndexIn]->ntagInds[2] ][0].histPdf , nNeutron1);  
		gtagN0 =  shuffleOrder(gtagN0_t);
		gtagN1 =  shuffleOrder(gtagN1_t);
	}

	vector<DataPoint> dataVector; dataVector.resize(nEvents);
	DataPoint rPoint;
	for(unsigned int n0 = 0; n0<nNeutron0; n0++) {
		rPoint.nTag = 0;
		rPoint.eVis = evisN0[n0];
		rPoint.gTag = gtagN0[n0];
		rPoint.weight = 1.0;
		rPoint.cosZen = 0;
		rPoint.azi = 0;
		dataVector[n0] = rPoint;
	}
	for(unsigned int n1 = 0; n1<nNeutron1; n1++) {
		rPoint.nTag = 1;
		rPoint.eVis = evisN1[n1];
		rPoint.gTag = gtagN1[n1];
		rPoint.weight = 1.0;
		rPoint.cosZen = 0;
		rPoint.azi = 0;
		dataVector[nNeutron0+n1] = rPoint;
	}
	return dataVector;
}

// To make projection histograms for data:
void MockData::makeProjHists() {
	
	// Initialize mock data histograms
	TString neutNames[3] = { "Nall", "N0", "N1" };
	TString tempName;
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	

	for(unsigned int n=0; n<nMax; n++) {
		if(evisProjs[n] != NULL) delete evisProjs[n];
		tempName = "hist" + mockName + "_eVisProj_" + neutNames[n];
		//cout << "Initializing mock data hist: " << tempName << endl;
		evisProjs[n] = (TH1D*) org->exps[expIndex]->pdfs[0]->evisKdes[n].histProj->Clone(tempName);
		evisProjs[n]->SetTitle(tempName);
		evisProjs[n]->Reset();
		evisProjs[n]->SetLineColor(6);
		evisProjs[n]->SetMarkerColor(6);
		evisProjs[n]->SetLineWidth(1);
		
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[0]+1; e++) {
				if(gtagProjs[n][e] != NULL) delete gtagProjs[n][e];
				tempName = "hist" + mockName + "_gTagProj_" + neutNames[n] + "_" + 
					org->exps[expIndex]->pdfs[0]->getGtagSubbinRange(e);
				gtagProjs[n][e] = (TH1D*) 
					org->exps[expIndex]->pdfs[0]->gtagKdes[n][e].histProj->Clone(tempName);
				gtagProjs[n][e]->SetTitle(tempName);
				gtagProjs[n][e]->Reset();
				gtagProjs[n][e]->SetLineColor(6);
				gtagProjs[n][e]->SetMarkerColor(6);
				gtagProjs[n][e]->SetLineWidth(1);
			}
		}
	}
	

	tempName = "hist" + mockName + "_nTagProj";
	hNeutProj = (TH1D*) org->exps[expIndex]->pdfs[0]->hNeutProj->Clone(tempName);
	hNeutProj->SetTitle(tempName);
	hNeutProj->Reset();
	hNeutProj->SetLineColor(6);
	hNeutProj->SetMarkerColor(6);
	hNeutProj->SetLineWidth(1);
	hNeutProj->SetDrawOption("E");
	

	// Loop over histograms
	for(unsigned int d=0;d<nPoints;d++) {
		//cout << "Point: " << d+1 << " of " << nPoints << endl;
		evisProjs[0]->Fill(dataPoint[d].eVis, dataPoint[d].weight);
		hNeutProj->Fill(dataPoint[d].nTag, dataPoint[d].weight);
		if(org->neutronInfo[expIndex] == true) {
			if(dataPoint[d].nTag==0) evisProjs[1]->Fill(dataPoint[d].eVis, dataPoint[d].weight);
			if(dataPoint[d].nTag==1) evisProjs[2]->Fill(dataPoint[d].eVis, dataPoint[d].weight);	
		}
		if(org->gammaInfo[expIndex] == true) {
			gtagProjs[0][0]->Fill(dataPoint[d].gTag, dataPoint[d].weight);		
			int subbin = org->exps[expIndex]->pdfs[0]->getEvisSubbinForGtag(dataPoint[d].eVis);
			gtagProjs[0][subbin]->Fill(dataPoint[d].gTag, dataPoint[d].weight);		
			
			if(org->neutronInfo[expIndex] == true) {
				if(dataPoint[d].nTag==0) {
					gtagProjs[1][0]->Fill(dataPoint[d].gTag, dataPoint[d].weight);
					gtagProjs[1][subbin]->Fill(dataPoint[d].gTag, dataPoint[d].weight);		
				}
				if(dataPoint[d].nTag==1) {
					gtagProjs[2][0]->Fill(dataPoint[d].gTag, dataPoint[d].weight);	
					gtagProjs[2][subbin]->Fill(dataPoint[d].gTag, dataPoint[d].weight);		
				}
			}
		}
	}	

	cout << "Prepared projection histograms for mock data!" << endl;
	return;
} // end of amke proj huists for mock data

// Save Projection Histograms
void MockData::saveProjHists() {
	
	org->outFile->cd();
	unsigned int nMax;
	if(org->neutronInfo[expIndex] == false) nMax = 1;
	if(org->neutronInfo[expIndex] == true) nMax = 3;
	
	for(unsigned int n=0; n<nMax; n++) {
		evisProjs[n]->Write();
		if(org->gammaInfo[expIndex] == true) {
			for(unsigned int e=0; e<org->evisSubbinsForGtag[0]+1; e++) {
				gtagProjs[n][e]->Write();
			} // end of evis subbin for gtag loop
		} // end of gamma info if
	} // end of neutron loop

	hNeutProj->Write();
	//hSysVsChiProb->Write();
	//hIdealEvis->Write();

} // end of save projection hists

// Shuffles the order of elements in a vector: Fisher-Yates
vector<double> MockData::shuffleOrder(vector<double> inVec) {

	vector<double> outVec;
	unsigned int inSize = inVec.size();
	// Loop over # of elements
	for(unsigned int e=0;e<inSize;e++) {
		unsigned int diffSize = inVec.size(); // remaining elements
		unsigned int remIndex = gRandom->Integer(diffSize);
		outVec.push_back(inVec[remIndex]);
		inVec.erase(inVec.begin()+remIndex); // remove the element
	}
	if( inSize!=outVec.size() ) cout << "ERROR in shuffleOrder! Input/Output sizes differ" << endl;
	return outVec;
} // end of shuffle

// Read data to fit!
unsigned int MockData::readData() {
	
	// Read experiment data tree and set branches:
	TFile* inDataFile = new TFile(org->dataFiles[expIndex],"READ");
	TTree* inDataTree = (TTree*) inDataFile->Get(org->dataTrees[expIndex]);
	
	float eVis;
	float gTag;
	int nTag;
	float dWall;
	float zenith;
	float azi;
	
	inDataTree->SetBranchAddress("eVis",&eVis);          
	inDataTree->SetBranchAddress("gTag",&gTag);          
	inDataTree->SetBranchAddress("nTag",&nTag);          
	inDataTree->SetBranchAddress("dWall",&dWall);          
	inDataTree->SetBranchAddress("zenith",&zenith);          
	inDataTree->SetBranchAddress("azi",&azi);          
	
	int passCtr = 0;
	unsigned int treeEvents = inDataTree->GetEntries();
	cout << "Data Tree has: " << treeEvents << " total events!" << endl;
	for(unsigned int eCtr = 0; eCtr<treeEvents; eCtr++) {
		
		inDataTree->GetEntry(eCtr);
		// dWall, eVis and gTag limits from card!
		if( !(org->wallCutLow[0]<=dWall && dWall<org->wallCutUp[0]) ) continue;
		if( !(org->evisLow[expIndex]<=eVis && eVis<org->evisUp[expIndex]) ) continue;
		if( !(org->gtagLow[expIndex]<=gTag && gTag<org->gtagUp[expIndex]) ) continue;
		if( !(nTag<2) ) continue;
		
		passCtr++;
		//cout << "Evis: " << eVis << " gTag: " << gTag << " nTag: " << nTag
		//	   << " dWall: " << dWall << " passCtr: " << passCtr << endl;
		// Put the values into a new dataPoint struct
		DataPoint temp;
		temp.eVis = eVis;
		temp.gTag = gTag;
		temp.nTag = nTag;
		temp.weight = 1.0;
		temp.cosZen = zenith;
		temp.azi = azi;
		dataPoint.push_back(temp);			
	
	} // end of loop over data events
	
	//cout << "Final Passed: " << passCtr << " Returned: " << dataPoint.size() << endl;
	return dataPoint.size();
} // end of read data
