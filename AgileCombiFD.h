/*

Agile Combinatorial Factor Decomposition (AgileFD)
Modified 11/30/2016
Version 1.2

================================================
Copyright 2016 Institute for Computational Sustainability

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
================================================

*/




#include <cstdio>
#include <time.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
#include <armadillo>
#include <fstream>
#include "assert.h"
#include "spline.h"
#include <set>
#include <cctype>

using namespace std;
using namespace arma;

class Coordinate {
public:
    string s;
    vector<double> v;
};

class Instance {
public:
    int N, L;
    mat XrdMat;
    vector<double> qvalues;
    vector<Coordinate> comp_xyz;
    vector<Coordinate> decomp_xy;

    double qmin, qmax;
    int length;
};

class Phase {
public:
    vector<double> xrd;
};

class PhaseMixture {
public:
    bool isInvolved;
    double proportion;
    double shift;
};

class AgileCombiFDParameters {
public:
    int M, K;
    mat initW, initH;

    AgileCombiFDParameters(int _M, int _K) {
        M = _M;
        K = _K;
    }
    AgileCombiFDParameters() {
        M = 10;
        K = 8;
    }
};

class SeedFreeze {
public:
    bool whetherValueInit;
    vector<int> valueFreeze;
    vector<double> Q;
    vector<vector<double> > valueSeeds;
    vector<set<int> > samples;

	SeedFreeze(bool w1, vector<double> q, vector<int> f1, vector<vector<double> > s1, vector<set<int> > s) {
		whetherValueInit = w1;
        Q = q;
		valueFreeze = f1;
		valueSeeds = s1;
		samples = s;
	}
};


class AgileCombiFD {
private:
    const static double eps = 0.0000000001; // threshold
    AgileCombiFDParameters param;

public:
    // constructor
    AgileCombiFD(AgileCombiFDParameters _param);

    // print shifts information into an output file
    void printShiftsInfo(vector<vector<PhaseMixture> > phase_mixtures, vector<Coordinate> xyz, int M, vector<double> Qlogsteps, vector<mat> H, char*);

    // when we shift a matrix, we use zeros to fill the blank rows or column
    mat InsertRows(mat, int, double, bool);

    // reconstruct signals
    mat Rec(mat, vector< mat >, vector<int>);

    // calculate the KL divergence
    double KLdiv(mat, mat, bool);

    // apply updating rules to W and H
    mat KLupdate(mat, mat &, vector<mat > &, vector<int>,mat &, mat, vector<int>, bool);

    // initialize W
    void InitW(mat &);

    // initialize H
    void InitH(vector<mat> &);

    // nomalize W
    mat NormalizeW(mat W);

    // calculate the L1 loss of two respective columns of two matrices
    double diff(mat, int, mat, int);

    // at the end of this algorithm, normalize H and adjust W
    void NormalizeWH(mat &W, vector<mat > &H);
    
    // solver (main function)
    void Solve(Instance, double, int, char*, mat, SeedFreeze, double, mat &, vector<mat> &, mat &, bool, time_t, char*, int, char*);

    // print reconstructed signals (mainly for debugging)
    void printRec(mat D, vector<double> Q, mat Dlog, vector<double> Qlogsteps, mat DRec, char *);

    // print all the necessary information to the output file (phases, reconstructed signals, concentrations, shifts)
    void generateSolTxt(int N, int M, int K, int L, int len, vector<double> Q, vector<double> Qlogsteps, char* sol_file, vector<mat> H, mat, mat, mat, int, char*, int);
};
