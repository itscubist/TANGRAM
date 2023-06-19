/********
*Author: Baran Bodur
*Date: 2021-07-08
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
#include "fitPdf.h"

using namespace std;

// Constructor from file and histogram templates ( based on kde from tree)
FitPdf::FitPdf(TString inName, TFile* inFile, TString intName, int wallOpt, 
		int histOpt, int inNeutType, int inGammaType, int inDirType, int inCol) {

	// Prepare rng
	randAlThor = new TRandom3();
	randAlThor->SetSeed(0);

	// Set name
	pdfName = inName;
	neutInfType = inNeutType;
	gammaInfType = inGammaType;
	dirInfType = inDirType;
	pdfColor = inCol;

	// Get wall and hist strings based on options
	TString wallStr, histStr;
	cout << "PDF: " << inName << " intName: " << intName << " Wall Opt: " << wallOpt 
			 << " histOpt: " << histOpt << " neutronInfo: " << inNeutType 
			 << " gammaInfo: " << inGammaType << " dirInfo: " << inDirType 
			 << "Color: " << inCol << endl; 
	if( wallOpt == 0 ) wallStr = "wall100-200";
	else if( wallOpt == 1 ) wallStr = "Wall200+";
	else if( wallOpt == 2 ) wallStr = "wall100+";
	else cout << " Bad Wall Option, will result in problems later!!!" << endl;
	if( histOpt == 0 ) histStr = "h2dEneGamRaw";
	else if( histOpt == 1 ) histStr = "h2dEneGamKde";
	else if( histOpt == 2 ) histStr = "h2dEneGamSemi";
	else cout << " Bad Hist Option, will result in problems later!!!" << endl;

	// Read basic neutron info
	TString tempName;
	tempName = "hNeut_" + intName + "_" + wallStr;
	cout << tempName << endl;
	TH1F* tempNeutH = (TH1F*) inFile->Get(tempName);  
	cout << "PDF: " << inName << " Prepared Histogram: " << tempNeutH->GetName() << endl;

	// Read direction info if direction info is used in the fit, if not assume isotropic
	// If inDirType == 2 read 4 histograms at 4 energy ranges as below
	tempName = "hDirProb_" + intName + "_" + wallStr;
	TH2F* tempHistDir;
	if( inDirType == 0 ) {
		tempHistDir = new TH2F(tempName, tempName, 2,-1,1,2,0,360);	
		tempHistDir->SetBinContent(1,1,0.25);
		tempHistDir->SetBinContent(1,2,0.25);
		tempHistDir->SetBinContent(2,1,0.25);
		tempHistDir->SetBinContent(2,2,0.25);
	}
	else if( inDirType != 0 ) {
		cout << tempName << endl;
		tempHistDir = (TH2F*) inFile->Get(tempName);	
	}
	tempName = "hDirProb_" + inName;
	hDir2d = (TH2F*) tempHistDir->Clone(tempName);
	delete tempHistDir;
	cout << "PDF: " << inName << " Prepared Histogram: " << hDir2d->GetName() << endl;

	if( inDirType==2 ) { // If inDirType == 2 read 4 histograms at 4 energy ranges as below
		TString dirEvisRanges[4] = {"_eVisR30-70","_eVisR70-100","_eVisR100-150","_eVisR150-200"};
		for(unsigned int e=0;e<4;e++) {	
			tempName = "hDirProb_" + intName + "_" + wallStr + dirEvisRanges[e];
			tempHistDir = (TH2F*) inFile->Get(tempName);
			tempName = "hDirProb_" + inName + dirEvisRanges[e];
			hDir2dEvisR[e] = (TH2F*) tempHistDir->Clone(tempName);
			delete tempHistDir;
			cout << "PDF: " << inName << " Prepared Histogram: " << hDir2dEvisR[e]->GetName() << endl;
		}
	} // if inDirType 2 (direction can depend on eVis)

	if( inNeutType==0 ) {
		
		expEvents = tempNeutH->Integral();
		simEvents = expEvents;
		//cout << "Exp events: " << expEvents << endl;
		neutProb[0] = 1.0;
		neutProb[1] = 0.0;
		
		// Smoothed (or KDEd) hists used in fitting studies
		tempName = histStr + "_" + intName + "_" + wallStr + "_neutAll";
		TH2F* tempHist = (TH2F*) inFile->Get(tempName);
		tempName = "hsEgNall_" + inName;
		hsEgNall = (TH2F*) tempHist->Clone(tempName);
		delete tempHist;
		hsEgNall->Scale(1.0/hsEgNall->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hsEgNall->GetName() << endl;
	
		// Raw hists for reference
		tempName = "h2dEneGamRaw_" + intName + "_" + wallStr + "_neutAll";
		tempHist = (TH2F*) inFile->Get(tempName);
		tempName = "hrEgNall_" + inName;
		hrEgNall = (TH2F*) tempHist->Clone(tempName);
		delete tempHist;
		hrEgNall->Scale(1.0/hrEgNall->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hrEgNall->GetName() << endl;
		
		// Set Raw Ranges
		hrEvisBins = hrEgNall->GetXaxis()->GetNbins();
		hrEvisLow = 1.5*hrEgNall->GetXaxis()->GetBinCenter(1)-
			0.5*hrEgNall->GetXaxis()->GetBinCenter(2);
		hrEvisUp = 1.5*hrEgNall->GetXaxis()->GetBinCenter(hrEvisBins)-
			0.5*hrEgNall->GetXaxis()->GetBinCenter(hrEvisBins-1);
		hrGtagBins = hrEgNall->GetYaxis()->GetNbins();
		hrGtagLow = 1.5*hrEgNall->GetYaxis()->GetBinCenter(1)-
			0.5*hrEgNall->GetYaxis()->GetBinCenter(2);
		hrGtagUp = 1.5*hrEgNall->GetYaxis()->GetBinCenter(hrGtagBins)-
			0.5*hrEgNall->GetYaxis()->GetBinCenter(hrGtagBins-1);
		
		// Set Final (Smoothed) Ranges
		hsEvisBins = hrEvisBins;
		hsEvisLow = hsEvisBins;
		hsEvisUp = hrEvisUp;
		hsGtagBins = hrGtagBins;
		hsGtagLow = hrGtagLow;
		hsGtagUp = hrGtagUp;

		// This is a PDF. no need for any errorbars (which are oversized and wrong anyway)
		for(unsigned int eBin = 0; eBin<hsEvisBins; eBin++) {
			for(unsigned int gBin = 0; gBin<hsGtagBins; gBin++) {
				hsEgNall->SetBinError(eBin+1,gBin+1,0);
				hrEgNall->SetBinError(eBin+1,gBin+1,0);
			}
		}

	}
	if(neutInfType == 1) { // If neutron info - sk4 onwards
		
		expEvents = tempNeutH->GetBinContent(1) + tempNeutH->GetBinContent(2); // expected events
		simEvents = expEvents;
		neutProb[0] = tempNeutH->GetBinContent(1)/expEvents;
		neutProb[1] = tempNeutH->GetBinContent(2)/expEvents;

		// Smoothed (or KDEd) hists used in fitting studies
		tempName = histStr + "_" + intName + "_" + wallStr + "_neut0";
		TH2F* tempHist0 = (TH2F*) inFile->Get(tempName);
		tempName = "hsEgN0_" + inName;
		hsEgN0 = (TH2F*) tempHist0->Clone(tempName);
		delete tempHist0;
		hsEgN0->Scale(1.0/hsEgN0->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hsEgN0->GetName() << endl;
		
		tempName = histStr + "_" + intName + "_" + wallStr + "_neut1";
		TH2F* tempHist1 = (TH2F*) inFile->Get(tempName);
		tempName = "hsEgN1_" + inName;
		hsEgN1 = (TH2F*) tempHist1->Clone(tempName);
		delete tempHist1;
		hsEgN1->Scale(1.0/hsEgN1->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hsEgN1->GetName() << endl;
		
		// Raw histograms for reference
		tempName = "h2dEneGamRaw_" + intName + "_" + wallStr + "_neut0";
		tempHist0 = (TH2F*) inFile->Get(tempName);
		tempName = "hrEgN0_" + inName;
		hrEgN0 = (TH2F*) tempHist0->Clone(tempName);
		delete tempHist0;
		hrEgN0->Scale(1.0/hrEgN0->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hrEgN0->GetName() << endl;
		
		tempName = "h2dEneGamRaw_" + intName + "_" + wallStr + "_neut1";
		tempHist1 = (TH2F*) inFile->Get(tempName);
		tempName = "hrEgN1_" + inName;
		hrEgN1 = (TH2F*) tempHist1->Clone(tempName);
		delete tempHist1;
		hrEgN1->Scale(1.0/hrEgN1->Integral());
		cout << "PDF: " << inName << " Prepared Histogram: " << hrEgN1->GetName() << endl;
		
		// Set Raw Ranges
		hrEvisBins = hrEgN0->GetXaxis()->GetNbins();
		hrEvisLow = 1.5*hrEgN0->GetXaxis()->GetBinCenter(1)-
			0.5*hrEgN0->GetXaxis()->GetBinCenter(2);
		hrEvisUp = 1.5*hrEgN0->GetXaxis()->GetBinCenter(hrEvisBins)-
			0.5*hrEgN0->GetXaxis()->GetBinCenter(hrEvisBins-1);
		hrGtagBins = hrEgN0->GetYaxis()->GetNbins();
		hrGtagLow = 1.5*hrEgN0->GetYaxis()->GetBinCenter(1)-
			0.5*hrEgN0->GetYaxis()->GetBinCenter(2);
		hrGtagUp = 1.5*hrEgN0->GetYaxis()->GetBinCenter(hrGtagBins)-
			0.5*hrEgN0->GetYaxis()->GetBinCenter(hrGtagBins-1);
		
		// Set Final (Smoothed) Ranges
		hsEvisBins = hrEvisBins;
		hsEvisLow = hsEvisBins;
		hsEvisUp = hrEvisUp;
		hsGtagBins = hrGtagBins;
		hsGtagLow = hrGtagLow;
		hsGtagUp = hrGtagUp;
		
		// This is a PDF. no need for any errorbars (which are oversized and wrong anyway)
		for(unsigned int eBin = 0; eBin<hsEvisBins; eBin++) {
			for(unsigned int gBin = 0; gBin<hsGtagBins; gBin++) {
				hsEgN0->SetBinError(eBin+1,gBin+1,0);
				hrEgN0->SetBinError(eBin+1,gBin+1,0);
				hsEgN1->SetBinError(eBin+1,gBin+1,0);
				hrEgN1->SetBinError(eBin+1,gBin+1,0);
			}
		}
		
		cout << "N0 Prob: " << neutProb[0] << " N1 Prob: " << neutProb[1] << endl;
		tempName = "hNeutProb_"+inName;
		hNeutProb = new TH1F(tempName,tempName,2,0,2); 
		hNeutProb->SetBinContent(1,neutProb[0]);
		hNeutProb->SetBinContent(2,neutProb[1]);
		hNeutProb->Print("v");
	
	}
	delete tempNeutH;
	
	// Make projection histograms with no weight
	cout << " Now making projection histograms! " << endl;
	makeProjHists(1.0);

}

// Constructor from tree and cut
FitPdf::FitPdf(TString inName, TTree* inTree, TString inCut, int inNeutType, 
		int inGammaType, int inCol) {

	// Prepare rng
	randAlThor = new TRandom3();
	randAlThor->SetSeed(0);

	// Set name
	pdfName = inName;
	neutInfType = inNeutType;
	gammaInfType = inGammaType;
	pdfColor = inCol;

/*
	// Set Raw Ranges
	hrEvisBins = 50;
	hrEvisLow = 20.0;
	hrEvisUp = 120.0;
	//hrEvisLow = 100.0;
	//hrEvisUp = 200.0;
	hrGtagBins = 22;
	//hrGtagBins = 1;
	hrGtagLow = -1.1;
	hrGtagUp = 1.1; 
	
	// Set Final (Smoothed) Ranges
	hsEvisBins = 35;
	hsEvisLow = 30.0;
	hsEvisUp = 100.0;
	//hsEvisBins = 50;
	//hsEvisLow = 100.0;
	//hsEvisUp = 200.0;
	hsGtagBins = 20;
	//hsGtagBins = 1;
	hsGtagLow = -1.0;
	hsGtagUp = 1.0; 
*/
	
	// Set Raw Ranges
	hrEvisBins = 90;
	hrEvisLow = 20.0;
	hrEvisUp = 200.0;
	hrGtagBins = 82;
	hrGtagLow = -2.2;
	hrGtagUp = 2.2; 
	
	// Set Final (Smoothed) Ranges
	hsEvisBins = 85;
	hsEvisLow = 30.0;
	hsEvisUp = 200.0;
	hsGtagBins = 80;
	hsGtagLow = -2.0;
	hsGtagUp = 2.0;


	// Prepare histograms for projection
	TString tempCut, tempName;
	if(neutInfType == 0) { // If  no neutron info - sk123
		tempCut = "weight*(" + inCut + ")";
		tempName = "hrEgNall_"+inName;
		hrEgNall = new TH2F(tempName,tempName,
				hrEvisBins,hrEvisLow,hrEvisUp,hrGtagBins,hrGtagLow,hrGtagUp);
		inTree->Project(tempName,"gTag:eVis",tempCut);
		hrEgNall->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hrEgNall->GetYaxis()->SetTitle("Gamma Tagging Parameter");
	} // end of if no neutron info 

	if(neutInfType == 1) { // If neutron info - sk4 onwards
		// Below is 0 neutron
		tempCut = "weight*(" + inCut + " && nTag==0)";
		tempName = "hrEgN0_"+inName;
		hrEgN0 = new TH2F(tempName,tempName,
				hrEvisBins,hrEvisLow,hrEvisUp,hrGtagBins,hrGtagLow,hrGtagUp);
		inTree->Project(tempName,"gTag:eVis",tempCut);
		hrEgN0->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hrEgN0->GetYaxis()->SetTitle("Gamma Tagging Parameter");
		// Below is 1 neutron
		tempCut = "weight*(" + inCut + " && nTag==1)";
		tempName = "hrEgN1_"+inName;
		hrEgN1 = new TH2F(tempName,tempName,
				hrEvisBins,hrEvisLow,hrEvisUp,hrGtagBins,hrGtagLow,hrGtagUp);
		inTree->Project(tempName,"gTag:eVis",tempCut);	
		hrEgN1->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hrEgN1->GetYaxis()->SetTitle("Gamma Tagging Parameter");
	} // end of if neutron info 

	makeSmoothDists("k3a",1); // Make smoothed histograms in the actual interest range of the PDFs

	// Get Expected Events and Normalize
	if(neutInfType == 0) { // If no neutron info - sk123
		expEvents = hsEgNall->Integral();
		simEvents = expEvents;
		hsEgNall->Scale(1.0/expEvents);
		neutProb[0] = 1; neutProb[1] = 0; // not used for no neutron info anyway
	}
	if(neutInfType == 1) { // If neutron info - sk4 onwards
		float expEventsN0 = hsEgN0->Integral();
		hsEgN0->Scale(1.0/expEventsN0);
		float expEventsN1 = hsEgN1->Integral();
		hsEgN1->Scale(1.0/expEventsN1);
		cout << "PDF: " << pdfName << " exp N0: " << expEventsN0 << " exp N1: " << expEventsN1 << endl;
		expEvents = expEventsN0 + expEventsN1;
		simEvents = expEvents;
		neutProb[0] = (float) expEventsN0/expEvents;
		neutProb[1] = (float) expEventsN1/expEvents;
		cout << "N0 Prob: " << neutProb[0] << " N1 Prob: " << neutProb[1] << endl;
		tempName = "hNeutProb_"+inName;
		hNeutProb = new TH1F(tempName,tempName,2,0,2); 
		hNeutProb->SetBinContent(1,neutProb[0]);
		hNeutProb->SetBinContent(2,neutProb[1]);
		hNeutProb->Print("v");
	}

	// Make projection histograms with no weight
	makeProjHists(1.0);
		
} // end of constructor!
	
FitPdf::~FitPdf() { // Destructor
/*
	if(neutInfType == 0) {
		delete hrEgNall;
		delete hsEgNall;
	}
	if(neutInfType == 1) {
		delete hNeutProb;
		delete hrEgN1;
		delete hrEgN0;
		delete hsEgN1;
		delete hsEgN0;
	}
	
	// Delete projection hists:
	delete hEvisProjN0;
	delete hEvisProjN1;
	delete hEvisProjNall;
	delete hGtagProjN0;
	delete hGtagProjN1;
	delete hGtagProjNall;
	delete hNtagProj;	

	delete hDir2d;
	delete hCosZenProj;
	delete hAziProj;
*/
} // end of destructor

float FitPdf::getProb(float eVis, float gTag, int nTag) {// get value of pdf at point (interpolates)
	if(gTag<hsGtagLow) gTag = hsGtagLow;
	if(gTag>hsGtagUp) gTag = hsGtagUp;
	float prob;
	//if(gammaInfType == 0) gTag = -10; // no need just to test gamma info is indeed not used

	if(neutInfType == 0) { // no neutron info
		if(gammaInfType != 0 ) {
			float tInterp = hsEgNall->Interpolate(eVis,gTag);
			//float tInterp = hsEgNall->GetBinContent(hsEgNall->FindBin(eVis,gTag));
			if(tInterp<=0) tInterp = (float) hsEgNall->GetMinimum(0); // if pdf is 0, put a very small
			prob = tInterp;
		}
		else if(gammaInfType == 0 ){ // do not use gamma tagging info
			float tInterp = hEvisProjNall->Interpolate(eVis); 
			if(tInterp<=0) tInterp = (float) hEvisProjNall->GetMinimum(0); //if pdf is 0, put a very small
			prob = tInterp;
		}
	}

	if(neutInfType == 1) { // if neutron info, use law of total probability
		if(nTag == 0) {
			if(gammaInfType != 0 ) {
				float tInterp = hsEgN0->Interpolate(eVis,gTag);
				//float tInterp = hsEgN0->GetBinContent(hsEgN0->FindBin(eVis,gTag));
				if(tInterp<=0) tInterp = (float) hsEgN0->GetMinimum(0);//if pdf is 0 , put a very small
				prob = neutProb[0]*tInterp;
			}
			else if(gammaInfType == 0 ) { // do not use gamma tagging info
				float tInterp = hEvisProjN0->Interpolate(eVis); 
				if(tInterp<=0) tInterp = (float) hEvisProjN0->GetMinimum(0); //if pdf is 0...
				prob = neutProb[0]*tInterp;
			}

		}
		if(nTag == 1) {
			if(gammaInfType != 0 ) {
				float tInterp = hsEgN1->Interpolate(eVis,gTag);
				//tInterp = hsEgN1->GetBinContent(hsEgN1->FindBin(eVis,gTag));
				if(tInterp<=0) tInterp = (float) hsEgN1->GetMinimum(0);//if pdf is 0 , put a very small
				prob = neutProb[1]*tInterp;
			}
			else if(gammaInfType == 0 ) { // do not use gamma tagging info
				float tInterp = hEvisProjN1->Interpolate(eVis); 
				if(tInterp<=0) tInterp = (float) hEvisProjN1->GetMinimum(0); //if pdf is 0...
				prob = neutProb[1]*tInterp;
			}
		}
	}
	return prob;
} // end of getProb


float FitPdf::getProb(float eVis, float gTag, int nTag, float cosZen, float azi) { // get value of pdf at point
	if( dirInfType == 0 ) return getProb(eVis,gTag,nTag); // if dirInfo is not used
	else if( dirInfType !=0 && dirInfType !=2 ) {
		float dirProb = hDir2d->GetBinContent( hDir2d->FindBin(cosZen, azi) );
		return dirProb*getProb(eVis,gTag,nTag);
	}
	else if(dirInfType == 2 ) {
		unsigned int eBin = 0;
		if( eVis>=30 && eVis<70 ) eBin=0;
		else if( eVis>=70 && eVis<100 ) eBin=1;
		else if( eVis>=100 && eVis<150 ) eBin=2;
		else if( eVis>=150 && eVis<=200 ) eBin=3;
		float dirProb = hDir2dEvisR[eBin]->GetBinContent( hDir2dEvisR[eBin]->FindBin(cosZen, azi) );
		return dirProb*getProb(eVis,gTag,nTag);
	}

}

float FitPdf::getProb(DataPoint inPoint) { // get value of pdf at point (interpolates)
	float eT = inPoint.eVis;
	float gT = inPoint.gTag;
	float nT = inPoint.nTag;
	float zT = inPoint.cosZen;
	float aT = inPoint.azi;
	return getProb(eT,gT,nT,zT,aT);
} // end of getProb

DataPoint FitPdf::getRandPoint() { // get a random data point from the pdf 
	DataPoint rPoint;
	rPoint.nTag = getRandNtag();	
	getRandEneAtNtag(rPoint.eVis, rPoint.gTag, rPoint.nTag);
	getRandDirAtEne(rPoint.cosZen, rPoint.azi, rPoint.eVis);
	rPoint.weight = 1.0;
	return rPoint; 
} // end of getRandPoint

int FitPdf::getRandNtag() { // get a random ntag (not the whole data point)
	if(neutInfType == 0 ) return -1; // nothing to do if no neutron info
	// If there are neutrons then roll the dice
	//float tRand = randAlThor->Rndm(); 	
	float tRand = gRandom->Rndm(); 	
	if (tRand <= neutProb[0]) return 0;
	else return 1;
} // end of GetRandNtag

void FitPdf::getRandEneAtNtag(float &eVis, float &gTag, int nTag) {//random eVis/gTag at given nTag	
	double eT=0, gT=0;
	if(neutInfType == 0) { // no neutron info
		if(gammaInfType != 0 ) hsEgNall->GetRandom2(eT,gT);	
		if(gammaInfType == 0 ) eT = hEvisProjNall->GetRandom();	
	}
	if(neutInfType == 1) { // no neutron info
		if(nTag == 0) {
			if( gammaInfType != 0 ) hsEgN0->GetRandom2(eT,gT);
			if( gammaInfType == 0 ) eT = hEvisProjN0->GetRandom();
		}
		else if(nTag == 1) {
			if( gammaInfType != 0 ) hsEgN1->GetRandom2(eT,gT);
			if( gammaInfType == 0 ) eT = hEvisProjN1->GetRandom();
		}
	}

	eVis = (float) eT;
	gTag = (float) gT;
} // end of getRandEneAtNtag
	
void FitPdf::getRandDirAtEne(float &cosZen, float &azi, float eVis) { // random zenith and azimuth angle at given eVis
	double cosZenT = 0, aziT = 0;
	//float tRand, binContent;
	//int binCx, binCy;

	// for now not evis dependent at all!!
	if( dirInfType != 0 && dirInfType !=2 ) { // only if direction info is to be used
		hDir2d->GetRandom2(cosZenT, aziT);		
		// For now just allow random values to be on of the 4 bin centers
		if( cosZenT < 0 ) cosZenT = -0.5; 
		else if( cosZenT >= 0 ) cosZenT = 0.5; 
		if( aziT < 180 ) aziT = 90; 
		else if( aziT >= 180 ) aziT = 270; 
	} 
	if( dirInfType == 2 ) { // only if direction info is to be used energy dependently
		unsigned int eBin = 0;
		if( eVis<70 ) eBin=0;
		else if( eVis>=70 && eVis<100 ) eBin=1;
		else if( eVis>=100 && eVis<150 ) eBin=2;
		else if( eVis>=150 ) eBin=3;
		hDir2dEvisR[eBin]->GetRandom2(cosZenT, aziT);		
		if( cosZenT < 0 ) cosZenT = -0.5; 
		else if( cosZenT >= 0 ) cosZenT = 0.5; 
		if( aziT < 180 ) aziT = 90; 
		else if( aziT >= 180 ) aziT = 270; 	
	}

	cosZen = (float) cosZenT;
	azi = (float) aziT;
}
	
void FitPdf::makeSmoothDists(TString smoothOpt, unsigned int nTimes) { // Provide smoothing for pdfs
	if(neutInfType == 0) {
		hTemp = (TH2F*) hrEgNall->Clone("hTemp"); // Clone raw hist  
		for(unsigned int i = 0; i< nTimes; i++) hTemp->Smooth(1,smoothOpt); // smooth clone
		// Manually fill the smoothed hist: 
		// The issue here is that raw range is wider than smoothed range
		TString tempName = "hsEgNall_"+pdfName;
		hsEgNall = new TH2F(tempName,tempName,
				hsEvisBins,hsEvisLow,hsEvisUp,hsGtagBins,hsGtagLow,hsGtagUp);
		hsEgNall->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hsEgNall->GetYaxis()->SetTitle("Gamma Tagging Parameter");
		hsEgNall->Sumw2();
		for(unsigned int e = 1; e < (hsEvisBins+1); e++) {
			float eCenter = hsEgNall->GetXaxis()->GetBinCenter(e);
			for(unsigned int g = 1; g < (hsGtagBins+1); g++) {
				float gCenter = hsEgNall->GetYaxis()->GetBinCenter(g);
				float binVal = hTemp->GetBinContent(hTemp->FindBin(eCenter,gCenter));
				hsEgNall->SetBinContent(e,g,binVal);	
			}
		}
		delete hTemp; 
	} // end of if no neutron info
	
	if(neutInfType == 1) {
		// Do 0 neutron case first:
		hTemp = (TH2F*) hrEgN0->Clone("hTemp"); // Clone raw hist  
		for(unsigned int i = 0; i< nTimes; i++) hTemp->Smooth(1,smoothOpt); // smooth clone
		// Manually fill the smoothed hist: 
		// The issue here is that raw range is wider than smoothed range
		TString tempName = "hsEgN0_"+pdfName;
		hsEgN0 = new TH2F(tempName,tempName,
				hsEvisBins,hsEvisLow,hsEvisUp,hsGtagBins,hsGtagLow,hsGtagUp);
		hsEgN0->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hsEgN0->GetYaxis()->SetTitle("Gamma Tagging Parameter");
		hsEgN0->Sumw2();
		for(unsigned int e = 1; e < (hsEvisBins+1); e++) {
			float eCenter = hsEgN0->GetXaxis()->GetBinCenter(e);
			for(unsigned int g = 1; g < (hsGtagBins+1); g++) {
				float gCenter = hsEgN0->GetYaxis()->GetBinCenter(g);
				float binVal = hTemp->GetBinContent(hTemp->FindBin(eCenter,gCenter));
				hsEgN0->SetBinContent(e,g,binVal);	
			}
		}
		delete hTemp; 
		
		// Do 1 neutron case second:
		hTemp = (TH2F*) hrEgN1->Clone("hTemp"); // Clone raw hist  
		for(unsigned int i = 0; i< nTimes; i++) hTemp->Smooth(1,smoothOpt); // smooth clone
		// Manually fill the smoothed hist: 
		// The issue here is that raw range is wider than smoothed range
		tempName = "hsEgN1_"+pdfName;
		hsEgN1 = new TH2F(tempName,tempName,
				hsEvisBins,hsEvisLow,hsEvisUp,hsGtagBins,hsGtagLow,hsGtagUp);
		hsEgN1->GetXaxis()->SetTitle("Visible Energy (MeV)");
		hsEgN1->GetYaxis()->SetTitle("Gamma Tagging Parameter");
		hsEgN1->Sumw2();
		for(unsigned int e = 1; e < (hsEvisBins+1); e++) {
			float eCenter = hsEgN1->GetXaxis()->GetBinCenter(e);
			for(unsigned int g = 1; g < (hsGtagBins+1); g++) {
				float gCenter = hsEgN1->GetYaxis()->GetBinCenter(g);
				float binVal = hTemp->GetBinContent(hTemp->FindBin(eCenter,gCenter));
				hsEgN1->SetBinContent(e,g,binVal);	
			}
		}
		delete hTemp;

	} // end if there are neutrons

} // end of makeSmoothDists

void FitPdf::savePlots(TFile* outFile) { // save the histograms
	outFile->cd();
	if(neutInfType == 0) {
		// First Main Hists
		hrEgNall->Write();
		hsEgNall->Write();
		
		// Now Projection Hists
		hEvisProjNall->Write();
		hGtagProjNall->Write();
	}
	if(neutInfType == 1) {
		//First Main Hists
		hNeutProb->Write();
		hrEgN0->Write();
		hsEgN0->Write();
		hrEgN1->Write();
		hsEgN1->Write();
	
		// Now Projection Hists
		hEvisProjN0->Write();
		hEvisProjN1->Write();
		hEvisProjNall->Write();
		hGtagProjN0->Write();
		hGtagProjN1->Write();
		hGtagProjNall->Write();
		hNtagProj->Write();
	}
	if(dirInfType !=0 ) {
		hDir2d->Write();
		hCosZenProj->Write();
		hAziProj->Write();
	}

} // end of savePlots

float FitPdf::scaleExpEvents(float scale) { // scale expected events
	expEvents *= scale; 
	return expEvents;
} // end of scaleExpEvents

float FitPdf::scaleSimEvents(float scale) { // scale expected events
	simEvents *= scale; 
	return simEvents;
} // end of scaleExpEvents

TH1D* FitPdf::getEneProj(int nTag,float scale) { // get energy projection of histograms
	TH1D *ht0, *ht1, *hReturn;
	int fBin = 1; // first bin
	int lBin = hsGtagBins; // last bin
	if(neutInfType == 0) { // no neutron info
		int nInf = -1;
		TString name = "hpEvis_" + pdfName + Form("_N%d",nInf);
		hReturn = hsEgNall->ProjectionX(name,fBin,lBin,"e");
	}
	else if(neutInfType == 1) { 
		int nInf = nTag;
		if(nTag==0) { // 0 neutrons
			TString name = "hpEvis_" + pdfName + Form("_N%d",nInf);
			hReturn = hsEgN0->ProjectionX(name,fBin,lBin,"e"); 
			hReturn->Scale(neutProb[0]);
		}	
		else if (nTag==1) { // 1 neutron
			TString name = "hpEvis_" + pdfName + Form("_N%d",nInf);
			hReturn = hsEgN1->ProjectionX(name,fBin,lBin,"e"); 	
			hReturn->Scale(neutProb[1]);
		}
		else { // total neutron info
			TString name = "hpEvis_" + pdfName + Form("_N%d",nInf);
			ht0 = hsEgN0->ProjectionX("ht0",fBin,lBin,"e");
			ht1 = hsEgN1->ProjectionX("ht1",fBin,lBin,"e");
			hReturn = hsEgN1->ProjectionX(name,fBin,lBin,"e");
			hReturn->Reset();
			hReturn->Add(ht0,ht1,neutProb[0],neutProb[1]);
			delete ht0;
			delete ht1;	
		}	// end of nTag input if
	} // end of neutInfo if

	// Scale the histogram if asked
	//hReturn->Sumw2();
	hReturn->Scale(scale);
	hReturn->GetXaxis()->SetTitle("Visible Energy (MeV)");
	return hReturn;
} // end of getEneProj

TH1D* FitPdf::getGtagProj(int nTag, float scale) { // get gTag projection of histograms
	TH1D *ht0, *ht1, *hReturn;
	int fBin = 1; // first bin
	int lBin = hsEvisBins; // last bin
	if(neutInfType == 0) { // no neutron info
		int nInf = -1;
		TString name = "hpGtag_" + pdfName + Form("_N%d",nInf);
		hReturn = hsEgNall->ProjectionY(name,fBin,lBin,"e"); 
	}
	else if(neutInfType == 1) { 
		int nInf = nTag;
		if(nTag==0) { // 0 neutrons
			TString name = "hpGtag_" + pdfName + Form("_N%d",nInf);
			hReturn = hsEgN0->ProjectionY(name,fBin,lBin,"e"); 
			hReturn->Scale(neutProb[0]);
		}	
		else if (nTag==1) { // 1 neutron
			TString name = "hpGtag_" + pdfName + Form("_N%d",nInf);
			hReturn = hsEgN1->ProjectionY(name,fBin,lBin,"e"); 	
			hReturn->Scale(neutProb[1]);
		}
		else { // total neutron info
			TString name = "hpGtag_" + pdfName + Form("_N%d",nInf);
			ht0 = hsEgN0->ProjectionY("ht0",fBin,lBin,"e");
			ht1 = hsEgN1->ProjectionY("ht1",fBin,lBin,"e");
			hReturn = hsEgN1->ProjectionY(name,fBin,lBin,"e");
			hReturn->Reset();
			hReturn->Add(ht0,ht1,neutProb[0],neutProb[1]);
			delete ht0;
			delete ht1;	
		}	// end of nTag input if
	} // end of neutInfo if

	// Scale the histogram if asked
	//hReturn->Sumw2();
	hReturn->Scale(scale);
	hReturn->GetXaxis()->SetTitle("Gamma Tagging Output");
	return hReturn;
} // end of getGtagProj

TH1D* FitPdf::getNtagProj(float scale) { // get gTag projection of histograms
	TString name = "hpNtag_" + pdfName;
	TH1D* hpNeut = new TH1D(name,name,2,0,2); 
	hpNeut->SetBinContent(1,neutProb[0]);
	hpNeut->SetBinContent(2,neutProb[1]);
	if(neutInfType == 0) { // no neutron info
		cout << " Neutron Info is not Available, Do not rely on this hist" << endl;
	}
	hpNeut->Scale(scale);
	hpNeut->GetXaxis()->SetTitle("Neutrons Tagged");
	return hpNeut;
}
	
TH1D* FitPdf::getCosZenProj(float scale) { // get cosZen projection of histograms
	TH1D *hReturn;
	TString name = "hpCosZen_" + pdfName;
	int fBin = 1;
	int lBin = hDir2d->GetYaxis()->GetNbins();
	hReturn = hDir2d->ProjectionX(name,fBin,lBin,"e"); 
	if(dirInfType == 0) { // no neutron info
		cout << " Direction Info is not Available, Do not rely on this hist" << endl;
	}
	hReturn->Scale(scale);
	hReturn->GetXaxis()->SetTitle("Cos(Zenith)");
	return hReturn;
}

TH1D* FitPdf::getAziProj(float scale) { // get azimuth projection of histograms
	TH1D *hReturn;
	TString name = "hpAzi_" + pdfName;
	int fBin = 1;
	int lBin = hDir2d->GetXaxis()->GetNbins();
	hReturn = hDir2d->ProjectionY(name,fBin,lBin,"e"); 
	if(dirInfType == 0) { // no neutron info
		cout << " Direction Info is not Available, Do not rely on this hist" << endl;
	}
	hReturn->Scale(scale);
	hReturn->GetXaxis()->SetTitle("Azimuth Angle (Degrees)");
	return hReturn;
}

void FitPdf::makeProjHists(float scale) { // stores projection histograms of this pdf
	hEvisProjN0 = getEneProj(0,scale);
	hEvisProjN1 = getEneProj(1,scale);
	hEvisProjNall = getEneProj(-1,scale);
	hGtagProjN0 = getGtagProj(0,scale);
	hGtagProjN1 = getGtagProj(1,scale);
	hGtagProjNall = getGtagProj(-1,scale);
	hNtagProj = getNtagProj(scale);
	hCosZenProj = getCosZenProj(scale);
	hAziProj = getAziProj(scale);
	
	// Also push them into vector:
	projHists.push_back(hEvisProjN0);
	projHists.push_back(hEvisProjN1);
	projHists.push_back(hEvisProjNall);
	projHists.push_back(hGtagProjN0);
	projHists.push_back(hGtagProjN1);
	projHists.push_back(hGtagProjNall);
	projHists.push_back(hNtagProj);
	projHists.push_back(hCosZenProj);
	projHists.push_back(hAziProj);
	
} // end of makeProjHists
