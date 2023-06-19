/********
*Author: Baran Bodur
*Date: 2021-07-13
*Description: Create and store a set of (mock-)data
*
********/
#ifndef UTILS_INCLUDED
#define UTILS_INCLUDED

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
#include "TH2F.h"
#include "TString.h"
#include "TRandom3.h"
// Project Headers

// Function declarations
void setRootStyle();

std::vector<double> getAsimovSample(TH1D* hIn, unsigned int nPoints);

#endif
