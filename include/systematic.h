/********
*Author: Baran Bodur
*Date: 2022-06-16
*Description: A simple class for storing systematic error related info 
*
********/
#ifndef SYSTEMATIC_INCLUDED
#define SYSTEMATIC_INCLUDED

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
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom3.h"
// Project Headers
#include "mockData.h"

struct Systematic

struct Systematic{
	float value;
	bool evisSys;
	bool gtagSys
	bool ntagSys;
	bool ctrSys;
	bool shiftOrWeight;
	float sigma;
};

class FitKdePdf{
public:

private:

};
#endif
