/********
*Author: Baran Bodur
*Date: 	2020-09-08
*Description: Common operations on vectors, such as getting mean, elementwise operations...
*
********/
#ifndef VECTORUTILS_INCLUDED
#define VECTORUTILS_INCLUDED

// c++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
// ROOT headers
#include "TString.h"
#include "TMath.h"


// Function declarations
template <typename T>
void printVector(std::vector<T> v1, std::string m1="", unsigned int nEl=0, unsigned int nCol=1);

template <typename T>
void printVector(std::vector< std::vector<T> > v1, std::string m1="", unsigned int nCol=0, 
		unsigned int nRow=0);

template <typename T>
T getSum(std::vector<T> inVec);

template <typename T>
T getProduct(std::vector<T> inVec);

template <typename T>
float getMean(std::vector<T> inVec);

template <typename T>
float getMomN(std::vector<T> inVec, float center, float momN);

template <typename T>
float getMomN(std::vector<T> inVec, float momN);

template <typename T>
float getStdDev(std::vector<T> inVec);

template <typename T>
float getStdErr(std::vector<T> inVec);

template <typename T>
float getCorrCoeff(std::vector<T> v1, std::vector<T> v2, bool opt=true);

template <typename T>
std::vector<T> checkNans(std::vector<T> inVec, bool& nanFlag); 

template <typename T>
std::vector<T> checkNans(std::vector<T> inVec, bool& nanFlag, std::vector<unsigned int>& index); 

template <typename T>
void equateSize(std::vector<T> v1, std::vector<T>& v2, T fillVal);

template <typename T>
std::vector<T> checkNans(std::vector<T> inVec, bool& nanFlag);

template <typename T>
std::vector<T> operator+(std::vector<T> v1, std::vector<T> v2);

template <typename T>
std::vector<T> operator-(std::vector<T> v1, std::vector<T> v2);

template <typename T>
std::vector<T> operator*(std::vector<T> v1, std::vector<T> v2);

template <typename T>
std::vector<T> operator/(std::vector<T> v1, std::vector<T> v2);

template <typename T>
std::vector<T> fillLinear(T start, T end, T step);

template <typename T>
std::vector<T> fillLinearN(T start, T end, unsigned int entries);

template <typename T>
std::vector < std::vector<T> > fillConst(unsigned int nRow, unsigned int nCol, T fillVal=0);

template <typename T>
unsigned int getMaxIndex(std::vector<T> v1);

template <typename T>
unsigned int getMinIndex(std::vector<T> v1);

template <typename T>
std::vector<T> cutVector(std::vector<T> v1, T cUp, T cLow, std::vector<T> cVec, int opt=0);

template <typename T>
std::vector<unsigned int> getRemIndex(std::vector<T> v1, T cUp, T cLow, int opt=0);

template <typename T>
std::vector<T> cutVector(std::vector<T> v1, T cUp, T cLow, int opt=0);

template <typename T>
std::vector<T> cutVector(std::vector<T> v1, std::vector<unsigned int> ind, bool opt=false);

template <typename T>
std::vector<T> rearrVector(std::vector<T> v1, std::vector<unsigned int> ind);

// Low->High coparator
inline bool pCompLow(const std::pair<unsigned int,float>& l,const std::pair<unsigned int,float>& r);
// High->Low Comparators
inline bool pCompGre(const std::pair<unsigned int,float>& l,const std::pair<unsigned int,float>& r);

template <typename T>
std::vector<unsigned int> getSortedIndex(std::vector<T> v1, bool sortLow=true);

template <typename T> 
std::vector< std::vector<T> > sortVector2D(std::vector< std::vector<T> > v1, 
		unsigned int colInd, bool sortLow=true);

template <typename T>
std::vector<T> getColVector(std::vector< std::vector<T> > v1, unsigned int colInd);

template <typename T>
std::vector< std::vector<T> > cutVecCols(std::vector< std::vector<T> > v1, 
		std::vector<unsigned int> colInd, bool opt=false);

template <typename T> 
std::vector< std::vector<T> > cutVector2D(std::vector< std::vector<T> > v1, unsigned int colInd, 
		T cUp, T cLow, int opt=0); 

template <typename T> 
std::vector< std::vector<T> > transposeVector(std::vector< std::vector<T> > v1);

template <typename T> 
std::vector< std::vector<float> > getCorrCoeffMatrix(std::vector< std::vector<T> > v1,bool opt=true);

#endif
