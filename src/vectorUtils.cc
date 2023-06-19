/********
*Author: Baran Bodur
*Date: 2020-09-09
*Description: Common operations on vectors, such as getting mean, elementwise operations...
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
#include "vectorUtils.h"

using namespace std;

// Print 1D vector, useful for debugging
template <typename T>
void printVector(vector<T> v1, string m1, unsigned int nEl, unsigned int nCol) {
	if( nEl > v1.size() || nEl==0 ) nEl = v1.size();
	unsigned int nRow = nEl/nCol + 1;
	cout << m1 << endl;
	if(nEl == 0 ) { cout << "Empty Vector!" << endl; return; }
	for(unsigned int r=0; r < nRow; r++) {
		for(unsigned int c=0; c < nCol; c++) {
			if( (r*nCol+c) >= nEl ) break;
			if(c==0) cout << "Indices " << r*nCol+c << "-" << (r+1)*nCol-1 << " : ";
			cout << v1[r*nCol+c] << " ";	
		}
		cout << endl;
	}
	cout << "------------ End of 1D vector printing --------" << endl;
	return;
}
template void printVector<float>(vector<float> v1, string m1, unsigned int nEl, unsigned int nCol);
template void printVector<double>(vector<double> v1, string m1, unsigned int nEl, unsigned int nCol);
template void printVector<int>(vector<int> v1, string m1, unsigned int nEl, unsigned int nCol);
template void printVector<bool>(vector<bool> v1, string m1, unsigned int nEl, unsigned int nCol);
template void printVector<unsigned int>(vector<unsigned int> v1, string m1, 
		unsigned int nEl, unsigned int nCol);

// Print 2D vector, usful for debugging
template < typename T>
void printVector(vector< vector<T> > v1, string m1, unsigned int nCol, unsigned int nRow) {
	cout << m1 << endl;
	if(v1.size() == 0) { cout << "Empty Vector" << endl; return; }
	if(v1[0].size() == 0) { cout << "Empty Vector" << endl; return; }
	if( nRow > v1.size() || nRow==0 ) nRow = v1.size();
	if( nCol > v1[0].size() || nCol==0 ) nCol = v1[0].size();
	
	for(unsigned int r=0; r < nRow; r++) {
		for(unsigned int c=0; c < nCol; c++) {
			if(c==0) cout << "Row: " << r << " Columns 0 to " << nCol << " : ";
			cout << v1[r][c] << " ";  
		}
		cout << endl;
	}
	cout << "------------ End of 2D vector printing --------" << endl;
	return;
}
template void printVector<float>(vector< vector<float> > v1, string m1, 
		unsigned int nCol, unsigned int nRow);
template void printVector<double>(vector< vector<double> > v1, string m1, 
		unsigned int nCol, unsigned int nRow);
template void printVector<int>(vector< vector<int> > v1, string m1, 
		unsigned int nCol, unsigned int nRow);
template void printVector<bool>(vector< vector<bool> > v1, string m1, 
		unsigned int nCol, unsigned int nRow);
template void printVector<unsigned int>(vector< vector<unsigned int> > v1, string m1, 
		unsigned int nCol, unsigned int nRow);


// sum the elelements of a vector
template <typename T>
T getSum(vector<T> inVec) {
	return accumulate(inVec.begin(),inVec.end(),0.0);	
}
template float getSum<float> (vector<float> inVec);
template double getSum<double> (vector<double> inVec);
template int getSum<int> (vector<int> inVec);

// Multiply the vector
template <typename T>
T getProduct(vector<T> inVec) {
	return accumulate(inVec.begin(),inVec.end(),1.0, multiplies<T>() );	
}
template float getProduct<float> (vector<float> inVec);
template double getProduct<double> (vector<double> inVec);
template int getProduct<int> (vector<int> inVec);

// Return vector mean
template <typename T>
float getMean(vector<T> inVec) {
	if(inVec.size()==0) {
		cout << "This is getMean function, input vector size is 0" << endl;
		return 0;
	}
	T mean = getSum(inVec);
	mean/=(float)inVec.size();
	return mean;
}
template float getMean<float> (vector<float> inVec);
template float getMean<double> (vector<double> inVec);
template float getMean<int> (vector<int> inVec);

// calc Nth moment around a given point
template <typename T>
float getMomN(vector<T> inVec, float center, float momN) {
	if(inVec.size()==0) return 0;
	float res=0;
	for (unsigned int i=0; i < inVec.size(); i++) res+=pow((inVec[i]-center),momN);
	res/=((float)inVec.size());
	return res;
}
template float getMomN<float> (vector<float> inVec, float center, float momN);
template float getMomN<double> (vector<double> inVec, float center, float momN);
template float getMomN<int> (vector<int> inVec, float center, float momN);

// calc Nth moment around mean
template <typename T>
float getMomN(vector<T> inVec, float momN) {
	if(inVec.size()==0) return 0;
	float mean = getMean(inVec);
	return getMomN(inVec, mean, momN);
}
template float getMomN<float> (vector<float> inVec, float momN);
template float getMomN<double> (vector<double> inVec, float momN);
template float getMomN<int> (vector<int> inVec, float momN);

// return std dev
template <typename T>
float getStdDev(vector<T> inVec) {
	if(inVec.size()==0) {
		cout << "This is getStdDev function, input vector size is 0 "<< endl;
		return 0;
	}
	float var = getMomN(inVec, 2);
	return sqrt(var);	
}
template float getStdDev<float>(vector<float> inVec);
template float getStdDev<double>(vector<double> inVec);
template float getStdDev<int>(vector<int> inVec);

// return std error
template <typename T>
float getStdErr(vector<T> inVec) {
	if(inVec.size()==0) {
		cout << "This is getStdErr function, input vector size is 0 "<< endl;
		return 0;
	}
	return getStdDev(inVec)/( sqrt( (float)inVec.size() ) );
}
template float getStdErr<float>(vector<float> inVec);
template float getStdErr<double>(vector<double> inVec);
template float getStdErr<int>(vector<int> inVec);

// return correlation coeff
template <typename T>
float getCorrCoeff(vector<T> v1, vector<T> v2, bool opt) {
	if(v1.size()==0) {cout << "getCorrCoeff: First vector is empty!" << endl; return 0; }
	if(v2.size()==0) {cout << "getCorrCoeff: Second vector is empty!" << endl; return 0; }
	if(v1.size()!=v2.size()) {cout << "getCorrCoeff: Vectors are different size!" << endl; return 0; }
	float m1 = getMean(v1), m2 = getMean(v2);
	float res = 0;
	for (unsigned int i=0; i < v1.size(); i++) res+=(v1[i]-m1)*(v2[i]-m2);
	res /= (float)v1.size(); // covariance
	if(opt == true) res /= sqrt( getMomN(v1,m1,2)*getMomN(v2,m2,2) ); //corr coeff
	return res;
}
template float getCorrCoeff<float>(vector<float> v1, vector<float> v2, bool opt);
template float getCorrCoeff<double>(vector<double> v1, vector<double> v2, bool opt);
template float getCorrCoeff<int>(vector<int> v1, vector<int> v2, bool opt);

// Detect NaNs in a vector and return the vector without those NaNs
template <typename T>
vector<T> checkNans(vector<T> inVec, bool& nanFlag) { 
	vector<T> outVec;
	nanFlag = false;
	for( unsigned int i = 0; i < inVec.size(); i++) {
		if(isnan(inVec[i])) {
			nanFlag = true;
			cout << "NaN detected at index: " << i << endl;
		}
		else {
			outVec.push_back(inVec[i]);
		} 
	}
	return outVec;
}
template vector<float> checkNans<float>(vector<float> inVec, bool& nanFlag);
template vector<double> checkNans<double>(vector<double> inVec, bool& nanFlag);
template vector<int> checkNans<int>(vector<int> inVec, bool& nanFlag);

// Detect NaNs in a vector and return the vector without those NaNs,
// also push back NaN indexes to a vector
template <typename T>
vector<T> checkNans(vector<T> inVec, bool& nanFlag, vector<unsigned int>& index) { 
	vector<T> outVec;
	nanFlag = false;
	for( unsigned int i = 0; i < inVec.size(); i++) {
		if(isnan(inVec[i])) {
			nanFlag = true;
			cout << "NaN detected at index: " << i << endl;
			index.push_back(i);
		}
		else {
			outVec.push_back(inVec[i]);
		} 
	}
	return outVec;
}
template vector<float> checkNans<float>(vector<float> inVec, bool& nanFlag, 
		vector<unsigned int>& index);
template vector<double> checkNans<double>(vector<double> inVec, bool& nanFlag, 
		vector<unsigned int>& index);
template vector<int> checkNans<int>(vector<int> inVec, bool& nanFlag, 
		vector<unsigned int>& index);

// Modifies size of v2 to the size of v1, if longer remove last elements
// If shorter fill the newly added elements with fillVal
template <typename T>
void equateSize(vector<T> v1, vector<T>& v2, T fillVal) {
	unsigned int outSize = v1.size();
	unsigned int v2Size = v2.size();
	if(outSize == v2Size) return; // if same size do nothing
	else if(v2Size < outSize) { // if v2 is shorter
		cout << "2nd vector is shorter by: " << outSize-v2Size
			<< " elements, filling " << fillVal << " to the 2nd vector!" << endl;
		v2.resize(outSize); // resize v2 to v1
		fill(v2.begin()+v2Size,v2.end(),fillVal); // fill new elements with 0s
	}
	else if (v2Size > outSize) { // if v2 is longer
		cout << "2nd vector is longer by: " << v2Size-outSize
			<< " removing last elements" << endl;
		v2.resize(outSize); // resize v2 to v1	
	} 	
}
template void equateSize(vector<float> v1, vector<float>& v2, float fillVal);
template void equateSize(vector<double> v1, vector<double>& v2, double fillVal);
template void equateSize(vector<int> v1, vector<int>& v2, int fillVal);


// Vector addition operator
template <typename T>
vector<T> operator+(vector<T> v1, vector<T> v2) {
	T fVal = 0;
	equateSize(v1,v2,fVal);
	vector<T> vRes; vRes.reserve(v1.size());
	transform(v1.begin(), v1.end(), v2.begin(), back_inserter(vRes), plus<T>() );
	return vRes;
}
template vector<float> operator+(vector<float> v1, vector<float> v2);
template vector<double> operator+(vector<double> v1, vector<double> v2);
template vector<int> operator+(vector<int> v1, vector<int> v2);

// Vector subtraction operator
template <typename T>
vector<T> operator-(vector<T> v1, vector<T> v2) {
	T fVal = 0;
	equateSize(v1,v2,fVal);
	vector<T> vRes; vRes.reserve(v1.size());
	transform(v1.begin(), v1.end(), v2.begin(), back_inserter(vRes), minus<T>() );
	return vRes;
}
template vector<float> operator-(vector<float> v1, vector<float> v2);
template vector<double> operator-(vector<double> v1, vector<double> v2);
template vector<int> operator-(vector<int> v1, vector<int> v2);

// Vector multiplication operator
template <typename T>
vector<T> operator*(vector<T> v1, vector<T> v2) {
	T fVal = 0;
	equateSize(v1,v2,fVal);
	vector<T> vRes; vRes.reserve(v1.size());
	transform(v1.begin(), v1.end(), v2.begin(), back_inserter(vRes), multiplies<T>() );
	return vRes;
}
template vector<float> operator*(vector<float> v1, vector<float> v2);
template vector<double> operator*(vector<double> v1, vector<double> v2);
template vector<int> operator*(vector<int> v1, vector<int> v2);

// Vector division operator
template <typename T>
vector<T> operator/(vector<T> v1, vector<T> v2) {
	T fVal = 0;
	equateSize(v1,v2,fVal);
	vector<T> vRes; vRes.reserve(v1.size());
	transform(v1.begin(), v1.end(), v2.begin(), back_inserter(vRes), divides<T>() );
	return vRes;
}
template vector<float> operator/(vector<float> v1, vector<float> v2);
template vector<double> operator/(vector<double> v1, vector<double> v2);
template vector<int> operator/(vector<int> v1, vector<int> v2);

// Make a linear vector, start, end, increment
template <typename T>
vector<T> fillLinear(T start, T end, T step) {
	vector<T> vRes;
	vRes.push_back(start);
	unsigned int rCtr = 0;
	while( ( (vRes[rCtr]+step) <= end && start < end ) ||
			( (vRes[rCtr]+step) >= end && start > end ) ) {
		vRes.push_back(vRes[rCtr] + step);
		rCtr++;
	}
	return vRes;
}
template vector<float> fillLinear<float>(float start, float end, float step);
template vector<double> fillLinear<double>(double start, double end, double step);
template vector<int> fillLinear<int>(int start, int end, int step);
template vector<unsigned int> fillLinear<unsigned int>(unsigned int start, unsigned int end, 
		unsigned int step);

// Make a linear vector, start, end, number of entries
template <typename T>
vector<T> fillLinearN(T start, T end, unsigned int entries) {
	T step = (start - end)/((float)entries);
	vector<T> vRes; vRes.resize(entries); 
	vRes[0] = start;
 	for(unsigned int i = 1; i < entries; i++) {
		vRes[i] = vRes[i-1] + step; 
	}
	return vRes;
}
template vector<float> fillLinearN<float>(float start, float end, unsigned int entries);
template vector<double> fillLinearN<double>(double start, double end, unsigned int step);
template vector<int> fillLinearN<int>(int start, int end, unsigned int step);
template vector<unsigned int> fillLinearN<unsigned int>(unsigned int start, unsigned int end, 
		unsigned int step);

template <typename T>
vector < vector<T> > fillConst(unsigned int nRow, unsigned int nCol, T fillVal) {
	vector < vector<T> > vRes;
	vector<T> vTemp; vTemp.resize(nCol); 
	fill(vTemp.begin(),vTemp.end(),fillVal);
 	for(unsigned int i = 0; i < nRow; i++) vRes.push_back(vTemp);
 	return vRes;
}
template vector< vector<float> > fillConst<float>(unsigned int nRow, unsigned int nCol,
		float fillVal);
template vector< vector<double> > fillConst<double>(unsigned int nRow, unsigned int nCol,
		double fillVal);
template vector< vector<int> > fillConst<int>(unsigned int nRow, unsigned int nCol,
		int fillVal);
template vector< vector<bool> > fillConst<bool>(unsigned int nRow, unsigned int nCol,
		bool fillVal);


// Cuts a vector according to another vector's values and returns to resultant vector 
// Need to be same type of vector
// Inputs are vector to be cutten, upper limit, lower limit, reference vector, and option
// opt = 0 means take events between low and up and is the default value
// opt = 1 means do not use upper limit
// opt = 2 means do not use lower limit
// opt = 3 means cut between events low and up
// Convention is lower limit is inclusive and upper limit is exclusive
template <typename T>
vector<T> cutVector(vector<T> v1, T cUp, T cLow, vector<T> cVec, int opt) {
	vector<T> vRes;
	T fVal = 0;
	equateSize(cVec,v1,fVal);
	for(unsigned int i = 0; i < cVec.size(); i++) {
		if(opt==0) if(cVec[i] < cLow || cVec[i] >= cUp ) continue;
		if(opt==1) if(cVec[i] < cLow) continue;
		if(opt==2) if(cVec[i] >= cUp ) continue;
		if(opt==3) if(cVec[i] >= cLow && cVec[i] < cUp ) continue;
		vRes.push_back(v1[i]);
	}
	return vRes;
}
template vector<float> cutVector<float>(vector<float> v1, float cUp, float cLow, 
		vector<float> cVec, int opt);
template vector<double> cutVector<double>(vector<double> v1, double cUp, double cLow, 
		vector<double> cVec, int opt);
template vector<int> cutVector<int>(vector<int> v1, int cUp, int cLow, vector<int> cVec, int opt);


// Get a vector of indexes passing the cuts
// opt = 0 means take events between low and up and is the default value
// opt = 1 means do not use upper limit
// opt = 2 means do not use lower limit
// opt = 3 means cut between events low and up
// Convention is lower limit is inclusive and upper limit is exclusive
template <typename T>
vector<unsigned int> getRemIndex(vector<T> v1, T cUp, T cLow, int opt) {
	vector<unsigned int> vRes;
	for(unsigned int i = 0; i < v1.size(); i++) {
		if(opt==0) if(v1[i] < cLow || v1[i] >= cUp ) continue;
		if(opt==1) if(v1[i] < cLow) continue;
		if(opt==2) if(v1[i] >= cUp ) continue;
		if(opt==3) if(v1[i] >= cLow && v1[i] < cUp ) continue;
		vRes.push_back(i);
	}
	return vRes;
}
template vector<unsigned int> getRemIndex<float>(vector<float> v1, float cUp, float cLow, int opt);
template vector<unsigned int> getRemIndex<double>(vector<double> v1, double cUp,double cLow,int opt);
template vector<unsigned int> getRemIndex<int>(vector<int> v1, int cUp, int cLow, int opt);

// Cuts a vector based on its own values vector's values and returns to resultant vector 
// opt = 0 means take events between low and up and is the default value
// opt = 1 means do not use upper limit
// opt = 2 means do not use lower limit
// opt = 3 means cut between events low and up
// Convention is lower limit is inclusive and upper limit is exclusive
template <typename T>
vector<T> cutVector(vector<T> v1, T cUp, T cLow, int opt) {
	return cutVector(v1,cUp,cLow,v1,opt);
}
template vector<float> cutVector<float>(vector<float> v1, float cUp, float cLow, int opt);
template vector<double> cutVector<double>(vector<double> v1, double cUp, double cLow, int opt);
template vector<int> cutVector<int>(vector<int> v1, int cUp, int cLow, int opt);

// Cut vector based on input index vector
// Based on opt, it will either select or deselect the inputed indices
// Default is throwing out indices
// Preserves original vector's order
template <typename T>
vector<T> cutVector(vector<T> v1, vector<unsigned int> ind, bool opt) {
	if(v1.size() == 0 ) { cout << "Empty input vector!!" << endl; return v1;}
	if(ind.size() == 0 ) { cout << "Empty index vector, returning input vector!!" << endl; return v1;}
	sort(ind.begin(),ind.end(),greater<int>() );  // sort from high to low
	if(opt==false) {
		for(unsigned int i = 0; i < ind.size(); i++) {
 			if( ind[i] >= v1.size() ) continue; // Do not operate for indices larger than v1's size
 			v1.erase(v1.begin()+ind[i]);
		}
	} // end of if removing indices
	else if(opt==true) {
		while(ind[0] >= v1.size()) ind.erase(ind.begin()); // get rid of indices larger than v1 size
		unsigned int sV1 = v1.size(), sInd = ind.size();
		if(ind[0] < sV1-1 ) v1.erase(v1.begin()+ind[0]+1,v1.end()); // cut from last index to end
		for(unsigned int i = 1; i< sInd;i++ ) {
			int rLast = ind[i-1], rFirst = ind[i];
			//if( rLast >= sV1 ) rLast = sV1 - 1;
			if( (rLast - rFirst) <= 1 || rFirst >= sV1 ) continue; // if consecutive no cut needed
			v1.erase(v1.begin() + rFirst + 1, v1.begin() + rLast); // if not cut in between
		}
		if(ind[sInd-1] != 0) v1.erase(v1.begin(),v1.begin()+ind[sInd-1]); // cut beginning to last ind
	} // end of keepn indices
 	return v1;
}
template vector<float> cutVector<float>(vector<float>, vector<unsigned int>, bool opt);
template vector<double> cutVector<double>(vector<double>, vector<unsigned int>, bool opt);
template vector<int> cutVector<int>(vector<int>, vector<unsigned int>, bool opt);
template vector<TString> cutVector<TString>(vector<TString>, vector<unsigned int>, bool opt);

template vector< vector<float> > cutVector< vector<float> >(vector< vector<float> >, 
		vector<unsigned int>, bool opt);
template vector< vector<double> > cutVector< vector<double> >(vector< vector<double> >,
		vector<unsigned int>, bool opt);
template vector< vector<int> > cutVector< vector<int> >(vector< vector<int> >, 
		vector<unsigned int>, bool opt);

// A function that returns index of max value (if multiple occurances then the first occurance)
template <typename T>
unsigned int getMaxIndex(vector<T> v1) {
	if(v1.size()==0) {cout << " Empty Vector!!! Check Input" << endl; return 0;}
	T maxVal= v1[0];
	unsigned int maxInd=0;
	for(unsigned int i = 0; i< v1.size();i++ ) {
		if( maxVal < v1[i] ) {
			maxInd = i;
			maxVal = v1[i];	
		}
	}
	return maxInd;
} 
template unsigned int getMaxIndex<float>(vector<float> v1);
template unsigned int getMaxIndex<double>(vector<double> v1);
template unsigned int getMaxIndex<int>(vector<int> v1);
template unsigned int getMaxIndex<unsigned int>(vector<unsigned int> v1);

// A function that returns index of min value (if multiple occurances then the first occurance)
template <typename T>
unsigned int getMinIndex(vector<T> v1) {
	if(v1.size()==0) {cout << " Empty Vector!!! Check Input" << endl; return 0;}
	T minVal= v1[0];
	unsigned int minInd=0;
	for(unsigned int i = 0; i< v1.size();i++ ) {
		if( minVal > v1[i] ) {
			minInd = i;
			minVal = v1[i];	
		}
	}
	return minInd;
} 
template unsigned int getMinIndex<float>(vector<float> v1);
template unsigned int getMinIndex<double>(vector<double> v1);
template unsigned int getMinIndex<int>(vector<int> v1);
template unsigned int getMinIndex<unsigned int>(vector<unsigned int> v1);

// rearrange vector in the order of given indices
// In principle the rearranged vector can be longer than the original one (duplicate indices) 
template <typename T>
vector<T> rearrVector(vector<T> v1, vector<unsigned int> ind) {
	vector<T> vRes; vRes.reserve(ind.size());
	for(unsigned int i = 0; i < ind.size(); i++) {
		if(ind[i] < 0 || ind[i] >= v1.size() ) continue;
		vRes.push_back(v1[ind[i]]);
	}
	return vRes;
}
template vector<float> rearrVector<float>(vector<float>, vector<unsigned int>);
template vector<double> rearrVector<double>(vector<double>, vector<unsigned int>);
template vector<int> rearrVector<int>(vector<int>, vector<unsigned int>);
template vector< vector<float> > rearrVector< vector<float> >(vector< vector<float> >, 
		vector<unsigned int>);
template vector< vector<double> > rearrVector< vector<double> >(vector< vector<double> >, 
		vector<unsigned int>);
template vector< vector<int> > rearrVector< vector<int> >(vector< vector<int> >, 
		vector<unsigned int>);

// Comparator for pairs
// Low->High coparator
inline bool pCompLow(const std::pair<unsigned int,float>& l,const std::pair<unsigned int,float>& r) {
	return l.second < r.second; 
}
// High->Low Comparators
inline bool pCompGre(const std::pair<unsigned int,float>& l,const std::pair<unsigned int,float>& r) {
	return l.second > r.second; 
}

// function to get sorted index of a vector
template <typename T>
vector<unsigned int> getSortedIndex(vector<T> v1, bool sortLow) {
	vector< pair<unsigned int, T> > vSort; 
	vSort.resize(v1.size());
	for(unsigned int i = 0; i<v1.size(); i++) vSort[i] = make_pair(i,v1[i]);
	if(sortLow == true) sort(vSort.begin(),vSort.end(),pCompLow);
	else if(sortLow == false) sort(vSort.begin(),vSort.end(),pCompGre);
	vector<unsigned int> vInd; vInd.reserve(v1.size());
	for(unsigned int i = 0; i<vSort.size();i++) vInd.push_back(vSort[i].first);	
	return vInd;
}
template vector<unsigned int> getSortedIndex<float>(vector<float> v1, bool sortLow);
template vector<unsigned int> getSortedIndex<double>(vector<double> v1, bool sortLow);
template vector<unsigned int> getSortedIndex<int>(vector<int> v1, bool sortLow);

// Sort a 2d vector based on 1 of its columns
// Input the 2d vector, index of the column used as sort reference, and sorting asc. or desc.
template <typename T> 
vector< vector<T> > sortVector2D(vector< vector<T> > v1, unsigned int colInd, bool sortLow) {
	vector<T> sortCol = getColVector(v1,colInd); // get the vector on its own
	vector<unsigned int> sortInd = getSortedIndex(sortCol,sortLow); // apply cuts and get passed ind
	vector< vector<T> > vRes = rearrVector(v1,sortInd); // select those indices
	return vRes;
}
template vector< vector<float> > sortVector2D<float>(vector< vector<float> > v1, 
		unsigned int colInd, bool sortLow);
template vector< vector<double> > sortVector2D<double>(vector< vector<double> > v1, 
		unsigned int colInd, bool sortLow);
template vector< vector<int> > sortVector2D<int>(vector< vector<int> > v1, 
		unsigned int colInd, bool sortLow);

// Get a column of a 2d vector (a row is easily accessible via v[rowInd], but a column is not)
template <typename T>
vector<T> getColVector(vector< vector<T> > v1, unsigned int colInd) {
	vector<T> vRes;	
	for(unsigned int i=0;i<v1.size();i++) vRes.push_back(v1[i][colInd]);
	return vRes;
}
template vector<float> getColVector<float>(vector< vector<float> > v1, unsigned int colInd);
template vector<double> getColVector<double>(vector< vector<double> > v1, unsigned int colInd);
template vector<int> getColVector<int>(vector< vector<int> > v1, unsigned int colInd);

// Select/deselt columns of a 2d vector
// opt==true keep selected indices
// opt==false remove selected indices
template <typename T>
vector< vector<T> > cutVecCols(vector< vector<T> > v1, vector<unsigned int> colInd, bool opt) {
	if(v1.size() == 0) { cout << "Empty Vector" << endl; return v1; }
	if(v1[0].size() == 0) { cout << "Empty Vector" << endl; return v1; }
	if(colInd.size()==0) { cout << "Mepty Index Vector" << endl; return v1;}
	if(colInd.size()>=v1[0].size()) { cout << "#Indices >= #columns" << endl; return v1;}
	vector<unsigned int> usedInd;
	if(opt==true) usedInd = colInd;
	if(opt==false) { 
		unsigned int start = 0, end = v1[0].size()-1, inc = 1;
		usedInd = fillLinear(start,end,inc);
		usedInd = cutVector(usedInd, colInd, false);			
	}
	vector< vector<T> > vRes;	
	T fillVal = 0;
	vRes = fillConst(v1.size(),usedInd.size(),fillVal);
	for(unsigned int r = 0; r<v1.size(); r++) {
		for(unsigned int c = 0; c<usedInd.size(); c++) { 
			vRes[r][c] = v1[r][ usedInd[c] ];
		}
	}
	return vRes;
}
template vector< vector<float> > cutVecCols<float>(vector< vector<float> > v1, 
		vector<unsigned int> colInd, bool opt);
template vector< vector<double> > cutVecCols<double>(vector< vector<double> > v1, 
		vector<unsigned int> colInd, bool opt);
template vector< vector<int> > cutVecCols<int>(vector< vector<int> > v1, 
		vector<unsigned int> colInd, bool opt);


// Cut a 2d vector based on 1 of its columns
template <typename T> 
vector< vector<T> > cutVector2D(vector< vector<T> > v1, unsigned int colInd, T cUp, T 
		cLow, int opt) {
	vector<T> cutCol = getColVector(v1,colInd); // get the vector on its own
	vector<unsigned int> remInd = getRemIndex(cutCol, cUp, cLow, opt); // apply cuts and get passed ind
	vector< vector<T> > vRes = cutVector(v1,remInd,true); // select those indices
	return vRes;
}
template vector< vector<float> > cutVector2D<float>(vector< vector<float> > v1, 
		unsigned int colInd, float cUp, float cLow, int opt);
template vector< vector<double> > cutVector2D<double>(vector< vector<double> > v1, 
		unsigned int colInd, double cUp, double cLow, int opt);
template vector< vector<int> > cutVector2D<int>(vector< vector<int> > v1, 
		unsigned int colInd, int cUp, int cLow, int opt);



// Transposes a vector (Do not use it for linear algebra, not meant to be fast)
template <typename T> 
vector< vector<T> > transposeVector(vector< vector<T> > v1) {
	if(v1.size() == 0) { cout << "Empty Vector" << endl; return v1; }
	if(v1[0].size() == 0) { cout << "Empty Vector" << endl; return v1; }
	T fillVal = 0;
	vector< vector<T> > vRes = fillConst(v1[0].size(),v1.size(),fillVal);
	for(unsigned int r = 0; r<v1.size(); r++) {
		for(unsigned int c = 0; c<v1[0].size(); c++) { 
			vRes[c][r] = v1[r][c];	
		}
	}
	return vRes;
}
template vector< vector<float> > transposeVector<float>(vector< vector<float> > v1);
template vector< vector<double> > transposeVector<double>(vector< vector<double> > v1);
template vector< vector<int> > transposeVector<int>(vector< vector<int> > v1);
template vector< vector<unsigned int> > transposeVector<unsigned int>(
		vector< vector<unsigned int> > v1);
template vector< vector<bool> > transposeVector<bool>(vector< vector<bool> > v1);

// get cor coeff matrix of a 2d vector (columns are different datasets)
template <typename T>
vector< vector<float> > getCorrCoeffMatrix(vector< vector<T> > v1, bool opt) {
	unsigned int nCol = v1[0].size();
	float fillVal = 0.0;
	vector< vector <float> > vRes = fillConst(nCol,nCol,fillVal);
	for(unsigned int r = 0; r < nCol; r++) {
		vector<T> vTemp1 = getColVector(v1,r);
		if(opt==false) vRes[r][r] = getStdDev(vTemp1);
		else vRes[r][r] = 1.0; // cor coeff with itself is 1
		for(unsigned int c = 0; c < r; c++) {
			vector<T> vTemp2 = getColVector(v1,c);
			vRes[r][c] = getCorrCoeff(vTemp1,vTemp2,opt);
			vRes[c][r] = vRes[r][c]; //symmetric
		}
	}	
	return vRes;
}
template vector< vector<float> > getCorrCoeffMatrix<float>(vector< vector<float> > v1, bool opt);
template vector< vector<float> > getCorrCoeffMatrix<double>(vector< vector<double> > v1, bool opt);
template vector< vector<float> > getCorrCoeffMatrix<int>(vector< vector<int> > v1, bool opt);
