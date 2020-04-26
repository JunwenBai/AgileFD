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



#include "AgileCombiFD.h"
#define Max_N 999

using namespace std;
using namespace arma;

void AgileCombiFD::printShiftsInfo(vector<vector<PhaseMixture> > phase_mixtures, vector<Coordinate> xyz, int M, vector<double> Qlogsteps, vector<mat> H, char* sol_file) { // print shifts information
    char str[1000];
    int length = strlen(sol_file);
    // generate the output file name of "shifts_info" output file
    for (int i = 0; i < length; i++) {
        if (sol_file[i] == '.' && sol_file[i+1] == 't' && sol_file[i+2] == 'x' && sol_file[i+3] == 't') {str[i]='\0';break;}
        str[i] = sol_file[i];
    }
    strcat(str, "_Shifts_Information.txt");
    // open the output file
    freopen(str, "w", stdout);

    int ndim = xyz.size();

    int K = H[0].n_rows; // K: the number of phases
    int N = H[0].n_cols; // N: the number of sample points
    printf("M=%d\n", M); // M: the nubmer of shifted versions
    printf("K=%d\n", K);
    printf("N=%d\n", N);
    if (ndim == 3) {
        printf("\n// (x = %s, y = %s, z = %s)   %s,%s,%s are elements from instance file,", xyz[0].s.c_str(), xyz[1].s.c_str(), xyz[2].s.c_str(), xyz[0].s.c_str(), xyz[1].s.c_str(), xyz[2].s.c_str()); // xyz[i].s is the name of the i-th element
        printf(" or (x, y)   x, y are coordinates of decomposition space\n");
    } else {
        printf("\n// (x, y)   x, y are coordinates of decomposition space\n");
    }
 
    // ternary system coordinates
    for (int i = 0; i < ndim; i++) {
        if (i == 0) printf("x=");
        if (i == 1) printf("y=");
        if (i == 2) printf("z=");
        for (int j = 0; j < N; j++) {
            printf("%lf", xyz[i].v[j]);
            if (j == N-1) printf("\n");
            else printf(",");
        }
    }

    // print Qvalues
    printf("\n// Qvalues\n");
    printf("q=");
    int L = Qlogsteps.size();
    for (int i = 0; i < L; i++) {
        printf("%lf", Qlogsteps[i]);
        if (i == L-1) printf("\n");
        else printf(",");
    }

    // print H tensor
    printf("\n// Activation matrix H :  H[m](i,j) means the activation value for phase i at sample point j with shift m(or say, the m-th layer of the tensor)\n");
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            printf("H[%d](%d,*)=", m, k);
            for (int n = 0; n < N; n++) {
                printf("%lf", H[m](k,n));
                if (n == N-1) printf("\n");
                else printf(",");
            }
        }
    }

    fclose(stdout);
}

// constructor
AgileCombiFD::AgileCombiFD(AgileCombiFDParameters _param) {
    param = _param;
}

mat AgileCombiFD::InsertRows(mat W, int rowcount, double val, bool top = true) { // if top==true, row0 is inserted into the beginning. Otherwise, row0 is inserted into the end.
    rowvec row0(W.n_cols); // one row filled with zeros of length K
    row0.fill(val); // row0=[val,val,...,val]. Usually, val == 0.
    for (int r = 0; r < rowcount; r++) {
        if (top) {
            W.insert_rows(0, row0); // insert row0 to the beginning 
        } else {
            W.insert_rows(W.n_rows, row0); // insert row0 to the end
        }
    }
    return W;
}

// reconstruct signals
mat AgileCombiFD::Rec(mat W, vector<mat> H, vector<int> Phi) { // Phi is an array containing shift values w.r.t # of positions to shift(not the real shift values)
    int L = W.n_rows; // the length of one XRD pattern
    int K = W.n_cols; // # of phases
    int N = H[0].n_cols; // # of sample points
    int M = H.size(); // # of shifted version
    assert((int)Phi.size() == (int)M); // check whether Phi is valid or not

    mat Rec = zeros<mat>(L, N); // reconstructed signals (reconstructed D matrix)

    for (int i = 0; i < M; i++) {
        int p = Phi[i]; // shift down p positions
        mat WW = W.rows(0,L-p-1);
        WW = InsertRows(WW, p, 0.0); // top p lines are filled with zeros

        Rec += WW * H[i]; // one shfited version of W multiplies the corresponding layer of H tensor
    }
    return Rec;
}

// calculate KL-divergence
double AgileCombiFD::KLdiv(mat D, mat Rec, bool show = false) { // "show" indicates whether to show some intermediate results (mainly for debugging)
    int L = D.n_rows;
    int N = D.n_cols;


    double ckl = 0.0; // KL divergence(cost)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < N; j++) {
            double Dij = D(i,j); // the original value at position (i,j)
            double Rij = Rec(i,j); // the reconstructed value at position (i,j)
            if (show) printf("(%.10lf, %.10lf, %.10lf)\n", Dij, Rij, log((Dij+eps)/(Rij+eps)));
            if (show) printf("(%.10lf, %.10lf, %.10lf)\n",  (Dij+eps)/(Rij+eps), log((Dij+eps)/(Rij+eps)), Dij*log((Dij+eps)/(Rij+eps)));

            ckl += Dij*log((Dij+eps)/(Rij+eps)) - Dij + Rij; // KL divergence

            if (show) printf("%d: %.10lf\n", i, ckl);
        }
    }

    return ckl;
}

mat AgileCombiFD::KLupdate(mat D, mat &W, vector<mat> &H, vector<int> Phi, mat &VRec, mat beta, vector<int> freeze, bool whetherInit) { // updating W and H
    // default meanings
    int L = W.n_rows;
    int K = W.n_cols;
    int N = H[0].n_cols;
    int M = H.size();

    double cost_old = KLdiv(D, VRec); // old cost

    vector<mat> H_old; // a copy of H tensor
    for (int i = 0; i < M; i++) {
        mat Hi = H[i];
        H_old.push_back(Hi);
    }
    mat W_old = W; // a copy of W
    W = NormalizeW(W);

    /*
        START to update! Please refer to the technical report for more details.
    */

    mat VR = D / (VRec+eps); // an auxiliary matrix
    mat O = ones<mat>(L, N); // an all-1 matrix

    vector<mat> H_x;
    vector<mat> H_y;
    vector<mat> GradH;


    for (int i = 0; i < M; i++) {
        int p = Phi[i];
        mat Hi_x = W.rows(0, L-p-1).t()*VR.rows(p, L-1);
        H_x.push_back(Hi_x);
        mat Hi_y = W.rows(0, L-p-1).t()*O.rows(p, L-1);
        H_y.push_back(Hi_y);

        mat Gradi = Hi_x / (Hi_y+beta);
//        mat Gradi = Hi_x / (Hi_y+0.5); // do not erase
//        mat Gradi = Hi_x / (Hi_y); // do not erase
        GradH.push_back(Gradi);
    }

    double nuH = 1.0;
    double accel = 1.0;
    int pindex = 0;
    for (int i = 0; i < M; i++) {
        H[pindex] = H_old[pindex] % (pow(GradH[pindex], nuH));
        pindex++;
    }
    VRec = Rec(W, H, Phi);
    double cost = KLdiv(D, VRec);
    if (cost - cost_old > 0.0001 * cost) {
            nuH = std::max(nuH/2.0, 1.0);
    } else {
            nuH *= accel;
            cost_old = cost;
    }

    
    // update W
    VR = D / (VRec + eps); // auxilliary matrix
    mat W_x = zeros<mat>(L, K); // auxiliary matrix (numerator in the updating rule)
    mat W_y = zeros<mat>(L, K); // auxiliary matrix (denominator)

    double* h = new double[K];
    mat H1, W1;
    for (int i = 0; i < M; i++) {
        int p = Phi[i];

        // numerator
        mat Wxp = VR.rows(p, L-1) * H[i].t();
        Wxp = InsertRows(Wxp, p, 0.0, false);
        mat Wyp = O.rows(p, L-1) * H[i].t();
        Wyp = InsertRows(Wyp, p, 0.0, false); 
        
        H1 = zeros<mat>(L, K);
        H1 = Wyp % W;
        for (int k = 0; k < K; k++) {
            h[k] = 0;
            for (int j = 0; j < L; j++) {
                h[k] += H1(j, k);
            }
        }
        W1 = zeros<mat>(L, K);
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K; k++) {
                W1(l, k) = W(l,k) * h[k];
            }
        }
        W_x = W_x+Wxp+W1;
        
        // denominator
        W1 = zeros<mat>(L, K);
        W1 = Wxp % W;
        for (int k = 0; k < K; k++) {
            h[k] = 0;
            for (int l = 0; l < L; l++) {
                h[k] += W1(l, k);
            }
        }
        H1 = zeros<mat>(L,K);
        for (int l = 0; l < L; l++) {
            for (int k = 0; k < K; k++) {
                H1(l, k) = W(l, k) * h[k];
            }
        }

        W_y = W_y+Wyp+H1;
    }
    delete [] h;

    mat GradW = W_x / (W_y+eps);
    for (int i = 0; i < freeze.size(); i++) {
        if (whetherInit && freeze[i] == 1) { // freezing phases
            for (int l = 0; l < L; l++) GradW(l,i) = 1;
        }
    }


    W = W % (pow(GradW, nuH));
    VRec = Rec(W, H, Phi);
    cost = KLdiv(D, VRec);


    if (cost - cost_old > 0.0001 * cost) {
        nuH = max(nuH/2.0, 1.0);
    } else {
        nuH *= accel;
        cost_old = cost;
    }


    return VRec;
}

void AgileCombiFD::InitW(mat &W) { // initialize W
    int L = W.n_rows;
    int K = W.n_cols;

    for (int j = 0; j < K; j++) {
        for (int i = 0; i < L; i++)
            W(i,j) = (rand() % (Max_N+1)+1)/float(Max_N+1); // to avoid initial value 0
    }

}

void AgileCombiFD::InitH(vector<mat> &H) {
    int K = H[0].n_rows;
    int N = H[0].n_cols;
    int M = H.size();

    for (int m = 0; m < M; m++) {

        for (int k = 0; k < K; k++) {
            for (int n = 0; n < N; n++) {
                H[m](k, n) = (rand() % (Max_N+1)+1)/float(Max_N+1); // to avoid initial value 0
            }
        }

    }
}

mat AgileCombiFD::NormalizeW(mat W) { // normalize W phase by phase
    int L = W.n_rows;
    int K = W.n_cols;
    for (int j = 0; j < K; j++) {
        double W2norm = norm(W.cols(j, j), "fro"); // calculate frobenius norm of the j-th column
        for (int i = 0; i < L; i++) {
            W(i,j) /= W2norm;
        }
    }
    return W;
}

void AgileCombiFD::NormalizeWH(mat &W, vector<mat> &H) { // normalize H
    int L = W.n_rows;
    int K = H[0].n_rows;
    int N = H[0].n_cols;
    int M = H.size();

    for (int k = 0; k < K; k++) {
        double maxsumH = 0.0;
        for (int n = 0; n < N; n++) {
            double sumH = 0.0;
            for (int m = 0; m < M; m++) {
                sumH += H[m](k, n);
            }
            maxsumH = max(maxsumH, sumH);
        }

        // normalize maxsumH to be 1
        for (int m = 0; m < M; m++) {
            for (int i = 0; i < N; i++) {
                double valH = H[m](k,i);
                if (maxsumH < 1e-6) H[m](k,i) = 0;
                else H[m](k,i) = valH / maxsumH;
            }
        }

        // to maintain the resulting matrix of W*H, multiply W(i,k) by some maxsumH
        for (int i = 0; i < L; i++) {
            double valW = W(i,k);
            W(i,k) = valW * maxsumH;
        }
    }
}

// print reconstructed signals to a file
void AgileCombiFD::printRec(mat D, vector<double> Q, mat Dlog, vector<double> Qlogsteps, mat DRec, char* sol_file) {
    char str[1000];
    int length = strlen(sol_file);
    for (int i = 0; i < length; i++) {
        if (sol_file[i] == '.' && sol_file[i+1] == 't' && sol_file[i+2] == 'x' && sol_file[i+3] == 't') {str[i] = '\0'; break;} // remove suffix
        str[i] = sol_file[i];
    }
    strcat(str, "_Ori_Rec.txt"); // add new suffix and file type
    freopen(str, "w", stdout);

    int L = Dlog.n_rows; // the length of one XRD pattern
    int N = Dlog.n_cols; // # of sample points
    vector<double> xrd;
    printf("N=%d\n", N);
    printf("L=%d\n", L);
    printf("Q=");
    for (int i = 0; i < L; i++) {
        printf("%lf", Qlogsteps[i]); // qvalues
        if (i == L-1) printf("\n");
        else printf(",");
    }

    for (int n = 0; n < N; n++) {
        printf("O%d=", n);
        for (int j = 0; j < L; j++) { // the original XRD pattern at sample point n
            printf("%lf", Dlog(j, n));
            if (j == L-1) printf("\n");
            else printf(",");
        }
    }
    
    for (int n = 0; n < N; n++) {
        printf("R%d=", n);
        for (int j = 0; j < L; j++) { // the reconstructed XRD pattern at sample point n
            printf("%lf", DRec(j, n));
            if (j == L-1) printf("\n");
            else printf(",");
        }
    }
    
    fclose(stdout);
}

// calculate the L1 loss of two respective columns of two matrices
double AgileCombiFD::diff(mat W, int k, mat Dlog, int n) {
    int L = W.n_rows;
    k = 186, n = 186;
    double cost = 0.0;
    for (int l = 0; l < L; l++) {
        cost += abs(W(l,k)-Dlog(l,n));
    }
    return cost;
}

void AgileCombiFD::Solve(Instance data, double convergenceRate, int time_limit, char* sol_file, mat beta, SeedFreeze sf, double givenstepsize, mat &W, vector<mat> &H, mat &Dlog, bool randini, time_t init_time, char * buffer, int snapshot, char* stageSuffix) {
    time_t now_time; // current time
    time_t local_init_time;

    local_init_time = time(NULL);

    // default meanings
    int N = data.N;
    int M = param.M;
    int K = param.K;
    int L = data.L;
    double qmin = data.qvalues[0];
    double qmax = data.qvalues[L-1];

    fprintf(stderr, "original qmin: %lf, original qmax: %lf\n", qmin, qmax);
    if (data.qmin > 0.0) qmin = data.qmin;
    if (data.qmax > 0.0) qmax = data.qmax;

    // D matrix: the collection of XRD patterns at each sample point
    mat D = data.XrdMat;

    // convert Q into log space
    double Qlogmin = log(qmin);
    double Qlogmax = log(qmax);
    double Qstepsize = (Qlogmax-Qlogmin)/(L-1);
    fprintf(stderr, "std Qstepsize = %.50lf\n", Qstepsize);

    int len;
    if (givenstepsize > 1e-8) { // reset stepsize using a given stepsize
        len = (int)((Qlogmax-Qlogmin)/log(1+givenstepsize)+1.0);
        Qstepsize = (Qlogmax-Qlogmin)/(len-1);
    } else {
        len = L;
    }
    if (data.length > 0) {
        len = data.length;
        Qstepsize = (Qlogmax-Qlogmin)/(len-1);
    } else {
        len = L;
    }

    fprintf(stderr, "Qstepsize = %.50lf\n", Qstepsize);
    fprintf(stderr, "Qstepsize in normal space = %.50lf\n", exp(Qstepsize));
    fprintf(stderr, "Length of qvalues array:\nOld len = %d, New len = %d\n", L, len);

    vector<double> Q;
    vector<double> Qlog;
    vector<double> Qlogsteps;

    for (int i = 0; i < L; i++) {
        Q.push_back(data.qvalues[i]);
    }
    for (int i = 0; i < len; i++) {
        Qlog.push_back(Qlogmin+i*Qstepsize);
        Qlogsteps.push_back(exp(Qlogmin+i*Qstepsize)); // Qlogsteps is a geometric sequence such that shifting Qlogsteps 1 position to the right means multiplying each qvalue by a constant
    }


    Dlog = zeros<mat>(len, N);

    // since we have new qvalues(Qlogsteps), we need to generate a new D matrix(Dlog) according to these new qvalues
    for (int i = 0; i < N; i++) {
        vector<double> xrd;
        xrd.clear();
        for (int j = 0; j < L; j++) {
            xrd.push_back(D(j, i));
        }
        vector<double> X(Q.begin(), Q.end());
        vector<double> Y(xrd.begin(), xrd.end());
//        Spline s(Q, xrd); // 1-d interpolation
        tk::spline s;
        s.set_points(X, Y);

        for (int j = 0; j < len; j++) {
            double pos = Qlogsteps[j];
            if (pos < X[0]) pos = X[0];
            if (pos > X[L-1]) pos = X[L-1];
            Dlog(j, i) = s(pos);
            if (Dlog(j, i) < 0.0) Dlog(j, i) = 0.0; // ensure non-negativity
        }
    }

    // Init matrix W
    W = zeros<mat>(len,K);
    InitW(W); // random initialization

    // convert seeding values into log space
    int rows = sf.valueSeeds.size();
    mat Wseed = zeros<mat>(len,rows);
    for (int r = 0; r < rows; r++) {
        vector<double> xrd;
        xrd.clear();
        for (int l = 0; l < L; l++) {
            xrd.push_back(sf.valueSeeds[r][l]);
        }
        vector<double> X(sf.Q.begin(), sf.Q.end());
        vector<double> Y(xrd.begin(), xrd.end());
//        Spline s(X, Y);
        tk::spline s;
        s.set_points(X,Y);

        for (int l = 0; l < len; l++) {
            double pos = Qlogsteps[l];
            if (pos < X[0]) pos = X[l];
            if (pos > X[L-1]) pos = X[L-1];
            Wseed(l,r) = s(pos);
            if (Wseed(l,r) < 0.0) Wseed(l,r) = 0.0;
        }
    }
    

    vector<int> freeze;
    bool whetherInit = false;

    // if the user specifies to seed with certain values,
    if (sf.whetherValueInit) {
        int rows = Wseed.n_cols;
        for (int r = 0; r < min(rows, K); r++) {
            for (int l = 0; l < len; l++) {
                W(l,r) = Wseed(l,r);
            }
        }
        freeze = sf.valueFreeze;
        whetherInit = sf.whetherValueInit;
    }


    W = NormalizeW(W); // step 1: W becomes W_tilde

    // Init matrices H
    if (randini == true) { // randini indicates whether you want a random initialization or not
        // if yes,
        H.clear();
        for (int m = 0; m < M; m++) {
            mat Hm = zeros<mat>(K,N);
            H.push_back(Hm);
        }
        InitH(H); // initialize H randomly
    } else {
        // if no, we initialize H based on the given H

        for (int n = 0; n < N; n++) {
            int cnt = 0;
            for (int k = 0; k < K; k++) {
                double tmp = 0.0;
                for (int m = 0; m < M; m++) {
                    tmp += H[m](k,n);
                }
                if (tmp > 1e-6) { // this means phase k plays a role at sample point n
                    for (int m = 0; m < M; m++) H[m](k,n) = (rand() % (Max_N+1)+1.0)/(float)(Max_N+1);
                } else {
                    for (int m = 0; m < M; m++) H[m](k,n) = 0.0;
                }
            }
        }

    }

    if (whetherInit) {
            int cols = Wseed.n_cols;

            int kk = sf.samples.size();
            
            for (int k = 0; k < min(kk, K); k++) {
                for (int n = 0; n < N; n++) {
                    if (sf.samples[k].find(n) != sf.samples[k].end()) continue;
                    for (int m = 0; m < M; m++) {
                        H[m](k,n) = 0.0;
                    }
                }
            }
            
    }

    //Init Phi
    vector<int> Phi;
    for (int i = 0; i < M; i++) {
        Phi.push_back(i); // generate shift values with regard to how many positions to shift
    }


    mat DRec = Rec(W, H, Phi); // reconstructed D matrix
    double cost = KLdiv(Dlog, DRec); // calculate cost

    fprintf(stderr, "cost: %lf\n", cost);

    now_time = time(NULL); // current time
    int time_cost = now_time - local_init_time; // duration time
    fprintf(stderr, "Before KLupdate, we have spent %d seconds.\n", time_cost);


    int iter_cnt = 0;
    while (true) {
        DRec = KLupdate(Dlog, W, H, Phi, DRec, beta, freeze, whetherInit); // up-to-date reconstructed D matrix
        double new_cost = KLdiv(Dlog, DRec); // new cost
        
        if (snapshot >= 1 && iter_cnt > 0 && iter_cnt % snapshot == 0) {
            string snapshot_file(sol_file);
            char buf[50];
            sprintf(buf, "_iter_%d_%s.txt", iter_cnt, stageSuffix);
            string suffix(buf);
            string filename = snapshot_file.substr(0, snapshot_file.size()-4)+suffix;
            now_time = time(NULL);
            generateSolTxt(N, M, K, L, len, Q, Qlogsteps, (char*)filename.c_str(), H, Dlog, DRec, W, now_time-init_time, buffer, iter_cnt);
        }

        if (iter_cnt % 25 == 0)
            fprintf(stderr, "%d: old: %lf, new: %lf", iter_cnt, cost, new_cost);
        now_time = time(NULL);
        time_cost = now_time - local_init_time;
        
        // based on cost and new_cost, we determine whether to stop or not
        // termination criterion: 1. reach the convergence point  2. running time expires
        if ((cost-new_cost) < cost*convergenceRate+eps || time_cost > time_limit) {
            fprintf(stderr, "\n");
            fprintf(stderr, "%d: old: %lf, new: %lf\n", iter_cnt, cost, new_cost);
            break;
        }
        else {
            cost = new_cost; // new cost becomes old cost
            if (iter_cnt % 50 == 0) fprintf(stderr, "   Total time cost until now: %d s.", time_cost);
            if (iter_cnt % 25 == 0) fprintf(stderr, "\n");
        }
        iter_cnt++;
    }

    NormalizeWH(W, H); // normalize H
    DRec = Rec(W, H, Phi);
    
//    if (showRec) { // whether to print reconstructed signals
//        printRec(D, Q, Dlog, Qlogsteps, DRec, sol_file);
//    }
    
    
    vector<Phase> phases; // XRD patterns of phases
    vector<vector<PhaseMixture> > phase_mixtures; // shifts, proportion

    for (int k = 0; k < K; k++) {
        vector<double> xrd;
        for (int j = 0; j < len; j++) {
            xrd.push_back(W(j, k));
        }

        // convert qvalues from Qlogsteps array back to original Q array using interpolation
        vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
        vector<double> Y(xrd.begin(), xrd.end());

//        Spline s(X, Y);
        tk::spline s;
        s.set_points(X,Y);

        Phase phase;
        for (int j = 0; j < L; j++) {
            double pos = Q[j];
            if (pos < X[0]) pos = X[0];
            if (pos > X[len-1]) pos = X[len-1];
            phase.xrd.push_back(s(pos));
        }

        phases.push_back(phase);

        vector<PhaseMixture> mixtures;
        for (int i = 0; i < N; i++) {
            double sumH = 0.0;
            double shiftH = 0.0;

            for (int m = 0; m < M; m++) {
                // Notice that Phi[m] == m
                sumH += H[m](k,i); // total proportion
                double shift = Qlogsteps[m]/Qlogsteps[0]; // generate the real shift values
                shiftH += shift * H[m](k,i);
            }
            if (sumH < 1e-6) shiftH = 1.0; // phase k doesn't appear at sample point i
            else shiftH /= sumH; // weighted shift value

            PhaseMixture mixture;
            mixture.isInvolved = (sumH > 0.0001); // whether phase k involves at sample point i
            mixture.proportion = (sumH > 0.0001)?sumH:0.0; // proportion
            mixture.shift = shiftH; // shift
            mixtures.push_back(mixture);
        }
        phase_mixtures.push_back(mixtures);
    }
    
    // generate *_output.txt
    now_time = time(NULL);
    generateSolTxt(N, M, K, L, len, Q, Qlogsteps, sol_file, H, Dlog, DRec, W, now_time-init_time, buffer, iter_cnt);

    // generate *_Ori_Rec.txt
//    if (showShift) {
//        printShiftsInfo(phase_mixtures, data.comp_xyz, M, Qlogsteps, H, sol_file);
//    }


}

void AgileCombiFD::generateSolTxt(int N, int M, int K, int L, int len, vector<double> Q, vector<double> Qlogsteps, char* sol_file, vector<mat> H, mat D, mat DRec, mat W, int cost_time, char * buffer, int iter_cnt) {
    
    vector<Phase> phases; // XRD patterns of phases
    vector<vector<PhaseMixture> > phase_mixtures; // shifts, proportion

    for (int k = 0; k < K; k++) {
        vector<double> xrd;
        for (int j = 0; j < len; j++) {
            xrd.push_back(W(j, k));
        }

        // convert qvalues from Qlogsteps array back to original Q array using interpolation
        vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
        vector<double> Y(xrd.begin(), xrd.end());

//        Spline s(X, Y);
        tk::spline s;
        s.set_points(X,Y);

        Phase phase;
        for (int j = 0; j < L; j++) {
            double pos = Q[j];
            if (pos < X[0]) pos = X[0];
            if (pos > X[len-1]) pos = X[len-1];
            phase.xrd.push_back(s(pos));
        }

        phases.push_back(phase);

        vector<PhaseMixture> mixtures;
        for (int i = 0; i < N; i++) {
            double sumH = 0.0;
            double shiftH = 0.0;

            for (int m = 0; m < M; m++) {
                // Notice that Phi[m] == m
                sumH += H[m](k,i); // total proportion
                double shift = Qlogsteps[m]/Qlogsteps[0]; // generate the real shift values
                shiftH += shift * H[m](k,i);
            }
            if (sumH < 1e-6) shiftH = 1.0; // phase k doesn't appear at sample point i
            else shiftH /= sumH; // weighted shift value

            PhaseMixture mixture;
            mixture.isInvolved = (sumH > 0.0001); // whether phase k involves at sample point i
            mixture.proportion = (sumH > 0.0001)?sumH:0.0; // proportion
            mixture.shift = shiftH; // shift
            mixtures.push_back(mixture);
        }
        phase_mixtures.push_back(mixtures);
    }

    freopen(sol_file, "w", stdout);

    printf("%s", buffer); 
    printf("time cost=%d\n", cost_time); // running time
    printf("n_iters=%d\n", iter_cnt);
    printf("Params=[Q, B, C, R, HS, H, L]\n");

    printf("\nK=%d\n", K);

    double xrd[L+1];

    printf("\n// Phase pattern (basis)\n");
    // print Q values
    printf("Q=");
    for (int i = 0; i < L; i++) {
        double q = Q[i];
        printf("%.3lf", q);
        if (i == L-1) printf("\n");
        else printf(",");
    }
    // print phases patterns
    for (int k = 0; k < K; k++) {
        printf("B%d=", k+1);
        for (int i = 0; i < L; i++) {
            printf("%.6lf", phases[k].xrd[i]);
            if (i == L-1) printf("\n");
            else printf(",");
        }
    }

    // print the proportion each phase accounts for at every sample point
    printf("\n// Phase concentrations at each sample\n");
    for (int  n = 0; n < N; n++) {
        double sum = 0.0;
        for (int k = 0; k < K; k++) {
            sum += phase_mixtures[k][n].proportion;
        }
        sum = 1.0;
        printf("C%d=", n+1);
        for (int k = 0; k < K; k++) {
            double c;
            if (sum > 0.0) c = phase_mixtures[k][n].proportion/sum;
            else c = 0.0;
            printf("%.6lf", c);
            if (k == K-1) printf("\n");
            else printf(",");
        }
    }

    // print the reconstructed signals of each phase at each sample point
    printf("\n// Per-phase model for each sample\n");
    double logxrd[len+1];
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            printf("R%d_%d=", n+1, k+1); // Ri_j denotes the signals of phase j at sample point i
            for (int l = 0; l < len; l++) logxrd[l] = 0.0;

            for (int m = 0; m < M; m++) {
                for (int l = 0; l < len; l++) {
                    int tl = l-m;
                    if (tl < 0) continue;
                    if (tl > len-1) continue;
                    logxrd[l] += W(tl, k) * H[m](k,n);
                }
            }
            vector<double> X(Qlogsteps.begin(), Qlogsteps.end());
            vector<double> Y;
            for (int l = 0; l < len; l++) Y.push_back(logxrd[l]);

//            Spline s(X, Y);
            tk::spline s;
            s.set_points(X,Y);

            vector<double> rlist;
            for (int j = 0; j < L; j++) {
                double pos = Q[j];
                if (pos < X[0]) pos = X[0];
                if (pos > X[len-1]) pos = X[len-1];
                rlist.push_back(s(pos));
            }
            for (int l = 0; l < L; l++) {
                printf("%.6lf", rlist[l]);
                if (l == L-1) printf("\n");
                else printf(",");
            }
 
        }
    }

    // shifts
    printf("\n// Per-phase shift for each sample\n");
    for (int n = 0; n < N; n++) {
        printf("S%d=", n+1);
        for (int k = 0; k < K; k++) {
            printf("%lf", phase_mixtures[k][n].shift);
            if (k == K-1) printf("\n");
            else printf(",");
        }
    }

    // Real shifts
    printf("\n// Real shifts\n");
    printf("HS=");
    for (int m = 0; m < M; m++) {
        printf("%lf", Qlogsteps[m]/Qlogsteps[0]);
        if (m == M-1) printf("\n");
        else printf(",");
    }

    // H tensor
    printf("\n// H tensor\n");
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < K; k++) {
            printf("H[*](%d,%d)=", k+1, n+1);
            for (int m = 0; m < M; m++) {
                printf("%lf", H[m](k,n));
                if (m != M-1) printf(",");
                else printf("\n");
            }
        }
    }

    // sample contribution
    printf("\n// Per-sample contribution (L1 loss)\n");
    printf("L=");
    for (int n = 0; n < N; n++) {
        double loss = 0.0;
        for (int l = 0; l < len; l++) {
            double d = D(l,n), dr = DRec(l,n);
//            loss += d*log((d+eps)/(dr+eps))-d+dr;
            loss += abs(d-dr);
        }
        printf("%lf", loss);
        if (n == N-1) printf("\n");
        else printf(",");
    }

    fclose(stdout);

}
