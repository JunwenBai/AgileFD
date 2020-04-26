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





#include <ilcplex/ilocplex.h>
#include <iostream>
#include <cstdio>
#include "tclap/CmdLine.h"
#include "AgileCombiFD.h"
ILOSTLBEGIN
#define INF 1

using namespace std;
using namespace TCLAP;

int N, n_threads;

class arg_struct {
public:
    int K;
    int M;
    int sample;
    double sumcost;
    double mipgap;
    vector<colvec> h;
};


class Slice {
public:
    set<int> samples;
    double lowcost;
    double highcost;
    int phaseIndex;
};


arg_struct args[400];


mat D, W;
vector<mat> H;


int readInt(string s) {
    char* str = (char*)s.c_str();
    int v;
    sscanf(str, "%d", &v);
    return v;
}

void readIntArray(vector<int> &vec, string buf) {
    vec.clear();
    size_t found = buf.find(',');
    while (found != string::npos) {
        vec.push_back(readInt(buf.substr(0, found)));
        buf.erase(0, found+1);
		found = buf.find(',');
    }
	vec.push_back(readInt(buf));
}


double readDouble(string s) {
    char* str = (char*)s.c_str();
    double v;
    sscanf(str, "%lf", &v);
    return v;
}

void readDoubleArray(vector<double> &vec, string buf) {
    vec.clear();
    size_t found = buf.find(',');
    while (found != string::npos) {
        vec.push_back(readDouble(buf.substr(0, found)));
        buf.erase(0, found+1);
		found = buf.find(',');
    }
	vec.push_back(readDouble(buf));
}

string readString(string s) {
	while (s.size() && !isalpha(s[0])) s.erase(0, 1);
	while (s.size() && !isalpha(s[s.size()-1])) s.erase(s.size()-1, 1);
	return s;
}

void readStringArray(set<string> &vec, string buf) {
    vec.clear();
    size_t found = buf.find(',');
    while (found != string::npos) {
        vec.insert(readString(buf.substr(0, found)));
        buf.erase(0, found+1);
        found = buf.find(',');
    }
    vec.insert(readString(buf));
}


// read instance files
Instance read_inst(char* filename, string &str) { // this function is mainly a parser

    freopen(filename, "r", stdin);
    Instance inst; // inst comprises all the information that an instance file has
    char c;
    string buf;
    int sampleNo = 0; // the index of a sample
	str = "";
	set<string> setComp;
	set<string> setDecomp;

    while ((c = getchar()) != EOF) {
        if (c == 10) continue; // reach the end of one line
        if (c == 13) { // reach the end of one line
            getchar();continue;
        }

		getline(cin, buf);
		buf = c+buf;
		
		if (buf.find("Composition=") != string::npos) {
			size_t found = buf.find('=');
			if (found == string::npos) continue;
			buf.erase(0, found+1);
			readStringArray(setComp, buf);
			continue;
		}
		if (buf.find("Decomposition=") != string::npos) {
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);
            readStringArray(setDecomp, buf);
			continue;
        }
		size_t found = buf.find('=');
		string subbuf = buf.substr(0, found);
		bool isComp = false;
		for (set<string>::iterator it = setComp.begin(); it != setComp.end(); it++) if (subbuf.find(*it) != string::npos) isComp = true;
		if (isComp) {
			Coordinate v1;
			size_t found = buf.find('=');
			if (found == string::npos) break;
			v1.s = buf.substr(0, found);
			buf.erase(0, found+1);
			readDoubleArray(v1.v, buf);
			inst.comp_xyz.push_back(v1);
			continue;
		}
		bool isDecomp = false;
        for (set<string>::iterator it = setDecomp.begin(); it != setDecomp.end(); it++) if (subbuf.find(*it) != string::npos) isDecomp = true;
        if (isDecomp) {
            Coordinate v1;
            size_t found = buf.find('=');
            if (found == string::npos) break;
            v1.s = buf.substr(0, found);
            buf.erase(0, found+1);
            readDoubleArray(v1.v, buf);
            inst.decomp_xy.push_back(v1);
            continue;
        }
		

        if (c == 'N') { // the number of sample points

			size_t found = buf.find('=');
			if (found == string::npos) continue;
			buf.erase(0, found+1);
			
			vector<int> temp;

			readIntArray(temp, buf);
			inst.N = temp[0];

        } else if (c == 'D') {
			
			if (buf.find("Description") != string::npos) str += buf+"\n";

		} else if (c == 'U') {
			
			if (buf.find("UUID") != string::npos) str += buf+"\n";

		} else if (c == 'Q') {
			
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            readDoubleArray(inst.qvalues, buf);
//			for (int i = 0; i < inst.qvalues.size(); i++) printf("%lf ", inst.qvalues[i]);printf("\n");
			inst.L = inst.qvalues.size();
			inst.XrdMat = zeros<mat>(inst.L, inst.N);

        } else if (c == 'I') { // XRD patterns

            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

			vector<double> temp;
            readDoubleArray(temp, buf);

			for (int i = 0; i < temp.size(); i++) inst.XrdMat(i, sampleNo) = temp[i];
			sampleNo++;

        }

    }

    // print some parameters to make sure the parser processes the data correctly
    cout << "L: " << inst.L << endl; // the length of one XRD pattern
    cout << "N: " << inst.N << endl; // the number of sample points
    cout << "size of qvalues: " << inst.qvalues.size() << endl; // the size of qvalues array (should == L)
    cout << "XrdMat: (" << inst.XrdMat.n_rows << " " << inst.XrdMat.n_cols << ")\n" << endl; // D matrix: the collection of all the XRD patterns at each sample point
    fclose(stdin);
    return inst;
}

mat getHumanIntelligence(char* filename, int N, int K, double beta_pow, vector<int> seeds, Slice &slice, bool whetherSlice) { // read human input(# of phases)

    mat beta = ones<mat>(K, N);
    beta = pow(beta, beta_pow);

    if (strlen(filename) == 0) return beta*beta_pow; // if there is no human input, then beta_pow is just a multiplier

    int specific_k = 0;
    for (int i = 0; i < seeds.size(); i++) {
        if (seeds[i] >= 186) specific_k = i;
    }

    FILE* file = freopen(filename, "r", stdin);
    if (file == NULL) return beta*beta_pow;

    double v;
    for (int i = 0; i < N; i++) {
        scanf("%lf", &v);
        for (int k = 0; k < K; k++) {
            if (v < 1e-6) {
                beta(k, i) = beta_pow*6; // if v < eps, the value of beta(k, i) is set to be the upper bound
                continue;
            }
            if (whetherSlice)
                if ( slice.samples.find(i) != slice.samples.end() ) v = 1.0;
            beta(k, i) = 1.0/v;
        }
    }

    // to enlarge the difference
    beta = pow(beta, beta_pow); // beta_pow now is the index

    fclose(stdin);
    return beta;

}

void initial_from_value(char * filename, vector<int> &freeze, vector<vector<double> > &seeds, bool &whetherInit, vector<set<int> > &samples, int N, Instance &data, vector<double> &Q) { // initialization based on given initial values(vectors)

    if (strlen(filename) == 0) { // users do not want to use this function
        whetherInit = false;
        return;
    }


    FILE* file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherInit = false;
        return;
    }


//	fprintf(stderr, "----------------------------------------------------\n");
    
	samples.clear();
    freeze.clear();
	seeds.clear();
	Q.clear();
    whetherInit = true;
    char c;
    int K;
    string buf;
	double shifting;

    // get values
    while ((c=getchar()) != EOF) {
        if (c == 10) continue;
        if (c == 13) {getchar();continue;}
        if (c == 'K') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);
            K = temp[0];

        } else if (c == 'V') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

			vector<double> temp;
            readDoubleArray(temp, buf);
			shifting = temp[0];

		} else if (c == 'Q') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            readDoubleArray(Q, buf);

		} else if (c == 'B') { // seeds: initial values for seeding

			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<double> phase;
            readDoubleArray(phase, buf);
            seeds.push_back(phase);

        } else if (c == 'F') { // initial values for freezing

			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);
            freeze.push_back(temp[0]);
            
        } else if (c == 'S') { // seeds: initial values for seeding

           	getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);

			set<int> s;
			s.clear();
			for (vector<int>::iterator it = temp.begin(); it != temp.end(); it++) s.insert(*it);

            samples.push_back(s);

        } else {
            getline(cin, buf);
            continue;
        }
    }

	for (vector<double>::iterator it = Q.begin(); it != Q.end(); it++) {
		*it *= shifting;
	}
}

void printVector(vector<int> v) { // for debugging
    int len = v.size();
    for (int i = 0; i < len; i++) printf("%d ", v[i]); printf("\n");
}

void initial_from_sample(char* filename, vector<int> &freeze, vector<vector<double> > &seeds, bool &whetherSampleInit, vector<set<int> > &samples, int N, Instance &data, vector<double> &Q) { // initialize phases based on the XRD patterns at some sample point
    if (strlen(filename) == 0) { // users do not use this functionality
        whetherSampleInit = false;
        return;
    }
    FILE * file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherSampleInit = false;
        return;
    }

	samples.clear();
	freeze.clear();
	seeds.clear();
	Q.clear();
    whetherSampleInit = true;
	char c;
    int K;
    string buf;
	double shifting;

    // get values
	while ((c=getchar()) != EOF) {
        if (c == 10) continue;
        if (c == 13) {getchar();continue;}
        if (c == 'K') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);
            K = temp[0];

        } else if (c == 'V') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

			vector<double> temp;
            readDoubleArray(temp, buf);
			shifting = temp[0];

		} else if (c == 'Q') {
			
			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            readDoubleArray(Q, buf);

		} else if (c == 'B') { // seeds: initial values for seeding

			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);
			int sample = temp[0];
			int L = data.XrdMat.n_rows;

			vector<double> phase;
			for (int l = 0; l < L; l++) phase.push_back(data.XrdMat(l,sample));
            seeds.push_back(phase);

        } else if (c == 'F') { // initial values for freezing

			getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);
            freeze.push_back(temp[0]);
            
        } else if (c == 'S') { // seeds: initial values for seeding

           	getline(cin, buf);
            buf = c+buf;
            size_t found = buf.find('=');
            if (found == string::npos) continue;
            buf.erase(0, found+1);

            vector<int> temp;
            readIntArray(temp, buf);

			set<int> s;
			s.clear();
			for (vector<int>::iterator it = temp.begin(); it != temp.end(); it++) s.insert(*it);

            samples.push_back(s);

        } else {
            getline(cin, buf);
            continue;
        }
    }

	for (vector<double>::iterator it = Q.begin(); it != Q.end(); it++) {
		*it *= shifting;
	}

}

void addGibbs3(mat &W, vector<mat> &H, int N, int K, int M) { // the simplest way to add Gibbs phase rule: choose three phases with highest coefficients in H tensor
    for (int n = 0; n < N; n++) {
        vector<double> vec; // vec[k] represents the sum of H[m](k,n) which means the coefficient of the m-th shifted version of phase k at sample n for all m. Notice that k and n are given
        for (int k = 0; k < K; k++) {
            double sumc = 0.0;
            for (int m = 0; m < M; m++) {
                sumc += H[m](k,n); // calculate the sum
            }
            vec.push_back(sumc);
        }
        int k1, k2, k3; // top 3 phases
        double max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max) {
                max = vec[k];
                k1 = k;
            }
        }
        max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max && k != k1) {
                max = vec[k];
                k2 = k;
            }
        }
        max = 0.0;
        for (int k = 0; k < K; k++) {
            if (vec[k] > max && k != k1 && k != k2) {
                max = vec[k];
                k3 = k;
            }
        }
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                if (k != k1 && k != k2 && k != k3) {
                    H[m](k,n) = 0.0;
                }
            }
        }
    }
}


void addGibbs(int K, int M, int sample0, double mipgap, vector<Coordinate> &xyz, bool onev) { // use mip to enforce Gibbs phase rule

	int sample = sample0;
	bool permitprint = (sample0 % 10 == 0);

    IloEnv env; // set up the environment
    int L = W.n_rows; // the length of one XRD pattern

    try {
        
        IloModel model(env);
        IloIntVarArray I(env);
		IloIntVarArray Im(env);
        for (int i = 0; i < K; i++) {
            I.add(IloIntVar(env,0,1)); // K indicator variables indicating whether to select the k-th phase or not
        }
		for (int i = 0; i < M*K; i++) {
			Im.add(IloIntVar(env,0,1));
		}
        IloNumVarArray Hc(env);
        for (int i = 0; i < M*K; i++) {
            Hc.add(IloNumVar(env, 0.0, IloInfinity)); // K phases, each of which has M shifted versions
        }

        IloNumVarArray t(env); // the L1 loss of the differences between reconstructed signals (of length L) and original signals (one column of D matrix)
        for (int i = 0; i < L; i++) t.add(IloNumVar(env, 0.0, IloInfinity));
        
	    IloExpr obj(env); // objective function      
      
	    for (int l = 0; l < L; l++) { // t[l] = |reconstructed_signals[l]-original_signals[l]|
                IloExpr e(env); // reconstructed signals
                for (int m = 0; m < M; m++) {
                    for (int k = 0; k < K; k++) {
                        if (l-m < 0) continue;
                        e += W(l-m,k)*Hc[m*K+k];
                    }
                }

                // t[l] = |reconstructed_signals[l]-original_signals[l]|

                model.add(t[l]-e>=-D(l,sample));
                model.add(t[l]+e>=D(l,sample));

				obj += t[l];
        }
        
		model.add(IloMinimize(env, obj)); // minimize the objective function


        if (onev) {
		// ---------------new------------------
		for (int k = 0; k < K; k++) {
			IloExpr e(env);
			for (int m = 0; m < M; m++) {
				e += Im[m*K+k];
			}
			model.add(e<=1);
		}
		for (int k = 0; k < K; k++) {
			model.add(Im[0*K+k]+Im[1*K+k]>=Hc[0*K+k]);
			model.add(Im[(M-1)*K+k]+Im[(M-2)*K+k]>=Hc[(M-1)*K+k]);
			for (int m = 1; m < M-1; m++) {
				model.add(Im[(m-1)*K+k]+Im[m*K+k]+Im[(m+1)*K+k]>=Hc[m*K+k]);
			}
		}
		for (int k = 0; k < K; k++) {
			for (int m = 1; m < M; m++)	model.add(Hc[m*K+k]-Hc[(m-1)*K+k]>=Im[m*K+k]-1);
			for (int m = 0; m < M-1; m++) model.add(Hc[m*K+k]-Hc[(m+1)*K+k]>=Im[m*K+k]-1);
		}
		// ---------------end------------------
        }


        for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    model.add(I[k]-Hc[m*K+k]>=0); // if I[k] == 0, the coefficient of every shifted version of the k-th phase should be 0. Otherwise, there is no boundaries for the coefficients
                }
        }
        
        IloExpr e(env);
        for (int k = 0; k < K; k++) {
            e += I[k]; // calculate the number of used phases
        }

        // enforce Gibbs phase rule
		double x = xyz[0].v[sample], y = xyz[1].v[sample], z = xyz[2].v[sample];
		int numlim = 3;
		if (abs(x+y+z-1) < 1e-3) {
			if (x < 1e-6) numlim--;
			if (y < 1e-6) numlim--;
			if (z < 1e-6) numlim--;
		}

        model.add(e<=numlim);
        model.add(e>=1);
        
        // -------------------------------------

        IloCplex cplex(model);

        // -------------------------------------

        cplex.setParam(IloCplex::EpGap, mipgap); // set the mipgap
        cplex.setParam(IloCplex::Threads, 1);
        
        
        if (!cplex.solve()) { // start to solve the MIP
            env.error() << "Failed to optimize LP." << endl;
            throw(-1);
        }


//        sumcost += cplex.getObjValue(); // the minimum value of the objective function
        if (permitprint) fprintf(stderr, "Solution value = %lf\n", cplex.getObjValue());


        // get values of the variables
        IloNumArray Hvalues(env);
        cplex.getValues(Hvalues, Hc);
        
        vector<int> intvalues;
		vector<int> Imvalues;
        for (int i = 0; i < K; i++) {
            intvalues.push_back(cplex.getValue(I[i])); // indicators
        }


        int cnt = 0;
        for (int k = 0; k < K; k++) {
            double tmpsum = 0.0;
            for (int m = 0; m < M; m++) {
//                H[m](k,sample) = Hvalues[m*K+k];
                args[sample].h[m](k) = Hvalues[m*K+k]; // H tensor
                tmpsum += Hvalues[m*K+k];
            }
            if (tmpsum > 1e-6) cnt++;
        }
        if (permitprint) fprintf(stderr, "sample %d: %d\n", sample, cnt);

    }
    catch (IloException &e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();      

}

void readSlice(char* filename, Slice &slice, bool &whetherSlice) {
    if (strlen(filename) == 0) { // users do not use this functionality
        whetherSlice = false;
        return;
    }
    FILE * file = freopen(filename, "r", stdin);
    if (file == NULL) { // no such file
        whetherSlice = false;
        return;
    }

    whetherSlice = true;
    slice.samples.clear();
    char c = getchar();
    while (c != '=') c = getchar();
    while (true) {
        int sample;
        scanf("%d", &sample);
        slice.samples.insert(sample-1);

        char c = getchar();
        if (c == 10 || c == EOF) break;
        if (c == 13) {getchar();break;}
    }

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%lf", &slice.lowcost);

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%lf", &slice.highcost);

    c = getchar();
    while (c != '=') c = getchar();
    scanf("%d", &slice.phaseIndex);
    slice.phaseIndex--;

}


void addSlice(char* filename, char* filename2, int N, int K, double beta_pow, Slice &slice, bool whetherSlice, mat& beta) { // read human input(# of phases)

    if (strlen(filename) == 0) return;
    if (strlen(filename2) == 0) {
        beta.fill(slice.lowcost);
        beta = pow(beta, beta_pow);
    }
   
    int specific_k = slice.phaseIndex; 
    std::set<int>::iterator it;
    if (whetherSlice)
    for (it=slice.samples.begin(); it!=slice.samples.end(); ++it) {
        int i = *it;
        for (int k = 0; k < K; k++) {
            if (k == specific_k) beta(k, i) *= slice.lowcost;
            else beta(k, i) *= slice.highcost;
        }
    }

}


int counter = 0;

int main(int argc, char **argv)
{

    CmdLine cmd("AgileFD", ' ', "1.2");


    int FLAGS_k, FLAGS_m, FLAGS_time, FLAGS_seed, FLAGS_threads, FLAGS_snapshot, FLAGS_length;
    string FLAGS_inst, FLAGS_sol, FLAGS_humanInput, FLAGS_valueInit, FLAGS_sampleInit, FLAGS_slice;
    double FLAGS_c, FLAGS_beta, FLAGS_stepsize, FLAGS_sparsity, FLAGS_mipgap, FLAGS_qmin, FLAGS_qmax;
    bool FLAGS_shiftInfo, FLAGS_rec, FLAGS_gibbs, FLAGS_oneVersion;


    try {
//        instArg, solArg, kArg, mArg, timeArg, cArg, seedArg, betaArg, shiftInfoArg, humanInputArg, recArg, initialArg, sampleInit, stepsizeArg, sparsityArg, mipgapArg, gibbsArg;
        
        ValueArg<string> instArg("", "inst", "Input instance file", true, "", "string", cmd);
        ValueArg<string> solArg("", "sol", "The output file name of the solution", true, "output.txt", "string", cmd);
        ValueArg<int> kArg("", "k", "The number of phases", true, 8, "int", cmd);
        ValueArg<int> mArg("", "m", "The number of possible different shifts", true, 10, "int", cmd);
        ValueArg<int> timeArg("", "time", "The maximum time(seconds) spent to train the model that you could accept", false, 10000, "int", cmd);
        ValueArg<double> cArg("", "c", "Related to termination criterion: In one iteration, if (old_cost-new_cost)<c*old_cost, then the loop terminates", false, 0.00008, "double", cmd);
        ValueArg<int> seedArg("", "seed", "random seed for the random number generator.(The default value -1 means time(0))", false, -1, "int", cmd);
        ValueArg<string> valueInitArg("", "valueInit", "Initialization file containing seeds for phases and phase freezing", false, "", "string", cmd);
        ValueArg<string> sampleInitArg("", "sampleInit", "The filename containing initialzation from single-phase sample points", false, "", "string", cmd);
        ValueArg<double> stepsizeArg("", "stepsize", "Customized stepsize for resampling. default shift value 0 means the user just want to use the std stepsize", false, 0.00000, "dobule", cmd);
        ValueArg<double> sparsityArg("", "sparsity", "The overall sparsity coefficient", false, 1.0, "double", cmd);
        ValueArg<double> mipgapArg("", "mipgap", "mipgap for MIP. dafault: 0.1", false, 0.1, "double", cmd);
        SwitchArg gibbsArg("", "gibbs", "whether to enforce Gibbs phase rule", cmd, false);
		ValueArg<int> snapshotArg("", "snapshot", "Every ? iterations, output one snapshot of the lastest solution", false, 0, "int", cmd);
		ValueArg<double> qminArg("", "qmin", "The customized minimum q value (could be different from the minimum qvalue in the instance file)", false, -1, "double", cmd);
		ValueArg<double> qmaxArg("", "qmax", "The customized maximum q value (could be different from the maximum qvalue in the instance file)", false, -1, "double", cmd);
		ValueArg<int> lengthArg("", "length", "The customized length of the array of qvalues (could be different from the length of the given qvalues in the instance file)", false, -1, "int", cmd);
        
    cmd.parse(argc, argv);
    
    FLAGS_inst = instArg.getValue();
    FLAGS_sol = solArg.getValue();
    FLAGS_k = kArg.getValue();
    FLAGS_m = mArg.getValue();
    FLAGS_time = timeArg.getValue();
    FLAGS_c = cArg.getValue();
    FLAGS_seed = seedArg.getValue();
    FLAGS_valueInit = valueInitArg.getValue();
    FLAGS_sampleInit = sampleInitArg.getValue();
    FLAGS_stepsize = stepsizeArg.getValue();
    FLAGS_sparsity = sparsityArg.getValue();
    FLAGS_mipgap = mipgapArg.getValue();
    FLAGS_gibbs = gibbsArg.getValue();
	FLAGS_snapshot = snapshotArg.getValue();
	FLAGS_length = lengthArg.getValue();
	FLAGS_qmin = qminArg.getValue();
	FLAGS_qmax = qmaxArg.getValue();

    } catch (ArgException &e) {
        cerr << "error: " << e.error() << "for arg" << e.argId() << endl;
    }

    cout << "inst: " << FLAGS_inst << endl;
    cout << "sol: " << FLAGS_sol << endl;
    printf("k: %d\n", FLAGS_k);
    printf("m: %d\n", FLAGS_m);
    printf("time: %d\n", FLAGS_time);
    printf("c: %lf\n", FLAGS_c);
    printf("seed: %d\n", FLAGS_seed);
    cout << "initial: " << FLAGS_valueInit << endl;
    cout << "sampleInit: " << FLAGS_sampleInit << endl;
    printf("stepsize: %lf\n", FLAGS_stepsize);
    printf("sparsity: %lf\n", FLAGS_sparsity);
    printf("mipgap: %lf\n", FLAGS_mipgap);
    cout << "gibbs: " << (FLAGS_gibbs?"true":"false") << endl;
	printf("snapshot: %d\n", FLAGS_snapshot);
	printf("qmin: %lf\n", FLAGS_qmin);
	printf("qmax: %lf\n", FLAGS_qmax);
	printf("length: %d\n", FLAGS_length);

	
    FLAGS_oneVersion = false;
	if (FLAGS_seed == -1) {FLAGS_seed = time(0); srand(FLAGS_seed);}
    else srand(FLAGS_seed);


	char buffer[500];

    AgileCombiFDParameters param(FLAGS_m, FLAGS_k);
	AgileCombiFD conv_nmf(param);

    time_t init_time = time(NULL); // record the starting time

	string desp;
    Instance inst = read_inst((char*)(FLAGS_inst.c_str()), desp); // read the instance file

	
	sprintf(buffer, "%sParamters=solver: \"AgileFD 1.2\", inst: %s, sol: %s, k: %d, m: %d, time: %d, c: %lf, seed: %d, initial: %s, sampleInit: %s, stepsize: %lf, sparsity: %lf, mipgap: %lf, gibbs: %s, snapshot: %d, qmin: %lf, qmax: %lf, length: %d\n", (char*)desp.c_str(), FLAGS_inst.c_str(), FLAGS_sol.c_str(), FLAGS_k, FLAGS_m, FLAGS_time, FLAGS_c, FLAGS_seed, FLAGS_valueInit.c_str(), FLAGS_sampleInit.c_str(), FLAGS_stepsize, FLAGS_sparsity, FLAGS_mipgap, (FLAGS_gibbs?"true":"false"), FLAGS_snapshot, FLAGS_qmin, FLAGS_qmax, FLAGS_length);


	inst.qmin = FLAGS_qmin;
	inst.qmax = FLAGS_qmax;
	inst.length = FLAGS_length;

	// N: # of sample points, K: # of phases, M: # of shifted versions
	N = inst.N;
    int K = FLAGS_k, M = param.M;
	   
   	vector<set<int> > samples;
    
    // value initialization
    vector<int> freeze;
    vector<vector<double> > seeds;
    bool whetherInit;
	vector<double> Q;
    initial_from_sample((char*)(FLAGS_sampleInit.c_str()), freeze, seeds, whetherInit, samples, N, inst, Q);
    initial_from_value((char*)(FLAGS_valueInit.c_str()), freeze, seeds, whetherInit, samples, N, inst, Q);

    mat beta = ones<mat>(K, N);
    beta = beta * FLAGS_sparsity;

    H.clear();
 
	SeedFreeze sf(whetherInit, Q, freeze, seeds, samples);

    // the first-round AgileFD

    conv_nmf.Solve(inst, FLAGS_c, FLAGS_time, (char*)FLAGS_sol.c_str(), beta, sf, FLAGS_stepsize, W, H, D, true, init_time, buffer, FLAGS_snapshot, (char*)"first");


    if (FLAGS_gibbs == false) { // if enforcing Gibbs phase rule is not required, then the program terminates
        return 0;
    }

    fprintf(stderr, "Solution without Gibbs rule added is done.\n\n\n\n");

    fprintf(stderr, "Start to add Gibbs phase rule......\n");

    double sumcost = 0.0;

    fprintf(stderr, "Start to generate args...\n");

    for (int i = 0; i < N; i++) {
        args[i].K = K;
        args[i].M = M;
        args[i].sample = i;
        args[i].mipgap = FLAGS_mipgap;
        args[i].h.clear();
        for (int j = 0; j < M; j++) {
            colvec colh = zeros<colvec>(K);
            args[i].h.push_back(colh);
        }
    }

	// n_threads

    fprintf(stderr, "preparation for MIP is ready!\n");

    time_t time0 = time(NULL);

    // #pragma omp parallel num_threads(n_threads) if (n_threads > 0)
    #pragma omp parallel
    {

        #pragma omp for
        for (int i = 0; i < N; i++) {
            addGibbs(K, M, i, FLAGS_mipgap, inst.comp_xyz, FLAGS_oneVersion);
        }

    }


    H.clear();
    for (int m = 0; m < M; m++) {
        mat tmp;
        for (int i = 0; i < N; i++) {
            tmp.insert_cols(tmp.n_cols, args[i].h[m]);
        }
        H.push_back(tmp);
    }

    time_t time1 = time(NULL);

    fprintf(stderr, "MIP time cost: %d seconds\n", time1-time0);

    fprintf(stderr, "Finish!!!\n\n\n\n");
    AgileCombiFD gibbs(param);


    // the second-round AgileFD (enforcing Gibbs phase rule)
    gibbs.Solve(inst, FLAGS_c, FLAGS_time, (char*)FLAGS_sol.c_str(), beta, sf, FLAGS_stepsize, W, H, D, false, init_time, buffer, FLAGS_snapshot, (char*)"second");

    fprintf(stderr, "Solution with Gibbs rule added is done.\n");

    return 0;
}
