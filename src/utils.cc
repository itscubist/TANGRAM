/********
 *
 *Author: Baran Bodur
 *Date: 2022-06-28
 *Description: Some utility functions for unbinned fitter with systematics 
 *
 ********/

// C++ libs
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

// ROOT Libraries
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TRandom3.h"

// Project Libraries
#include "utils.h"

// namespaces and typedefs
using namespace std;
using namespace TMath;

// Set root style
void setRootStyle() {
	gRandom->SetSeed(0);
	gStyle->SetOptFit(1111); 
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatColor(10);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.07);
	gStyle->SetOptStat(11);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFrameFillColor(10);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetHatchesSpacing(2);
}

// Given a ROOT histogram return a vector with most representative sample of requested size
vector<double> getAsimovSample(TH1D* hIn, unsigned int nPoints) {
	
	std::vector<double> cdfEnds; cdfEnds.resize(nPoints);
	std::vector<double> bins; bins.resize(nPoints);
	std::vector<double> samples; samples.resize(nPoints);
	if(nPoints == 0) return samples;
	
	unsigned int nBins = hIn->GetNbinsX();
	double integral = hIn->Integral();
	double *cumIntegral = hIn->GetIntegral();
	double start = 0.5/((double)nPoints);
	double step = 1.0/((double)nPoints);
	
	for(unsigned int i=0; i<nPoints; i++) {
		cdfEnds[i] = start + step*(double)i;
		bins[i] = TMath::BinarySearch(nBins,cumIntegral,cdfEnds[i]);
		samples[i] = hIn->GetBinCenter(bins[i]+1);
		//cout << "Point: " << i+1 << " cdfEnd: " << cdfEnds[i] << " Position: " << samples[i] << endl;
	}

	return samples;
}

