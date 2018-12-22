#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include "MersenneTwister.h"

using namespace std;

#define MAX_TIME 1000
#define J_MAX 100 //MAX number of jobs
#define O_MAX 10 //MAX number of operations per jobs
#define M_MAX 100 //MAX number of machines
#define L_MAX 10 //MAX number of maintenance tasks per machine
#define NB_ITERS 1000
#define SEED 2018

struct Parameters{
    int n=0; //number of jobs
    int m=0; //number of machines
    int n_o[J_MAX]; //number of operations of jobs i 
    int L[M_MAX]; //number of maintenance tasks on machine j
    bool A[J_MAX][O_MAX][M_MAX]; //availables machines for task k of job i
    double t[J_MAX][O_MAX][M_MAX]; //processing time for task k of job i on machine j
    double p[L_MAX][M_MAX]; //processing time for maintenance task l on machine j
    double t_e[L_MAX][M_MAX]; //maintenance completion time windows begin (earliest)
    double t_l[L_MAX][M_MAX]; //maintenance completion time windows end (lateness)
};

struct Solution{
    //affectation of the operation ik to machine j with position 
    int x[J_MAX][O_MAX][M_MAX];
    double W_tot=0; //Total workload (of all machines)
    double W_max=0; //Maximum Workload (of machines)
    double C_max=0; //Maximum completion time (of machines)
};


//static variables
static MTRand rand_gen(SEED);
static Parameters params;
static Solution solutions[NB_ITERS];

int read_parameters(string filename)
{
    ifstream file(filename);
    if(file.fail()) {
        cout << "Error : file \" "<< filename <<" \" not found" << endl;
        return 1;
    }
    string word;
    file >> word;
    cout << word << endl;
    file >> params.n;
    cout << "Got n : " << params.n << endl;
    file >> params.m;
    cout << "Got m : " << params.m << endl;
    for(int i = 0; i<params.n; i++) {
       file >> params.n_o[i];
       cout << "Got n_o["<<i<<"] : " << params.n_o[i] << endl;
    }
    for(int j = 0; j<params.m; j++) {
        file >> params.L[j];
        cout << "Got L["<<j<<"] : " << params.L[j] << endl;
    }
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            for(int j=0; j<params.m; j++) {
                string temp_word;
                file >> temp_word;
                if(temp_word.find('T') != std::string::npos) {
                    params.A[i][k][j]=true;
                } else {
                    params.A[i][k][j]=false;
                }
                cout << "Got A["<<i<<"]["<<k<<"]["<<j<<"] : " << temp_word << endl;
            }
        }
    }
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            for(int j=0; j<params.m; j++)
            {
                file >> params.t[i][k][j];
                cout << "Got t["<<i<<"]["<<k<<"]["<<j<<"] : " << params.t[i][k][j] << endl;
            }
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.p[j][l];
            cout << "Got p["<<j<<"]["<<l<<"] : " << params.p[j][l] << endl;
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_e[j][l];
            cout << "Got t_e["<<j<<"]["<<l<<"] : " << params.t_e[j][l] << endl;
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_l[j][l];
            cout << "Got t_l["<<j<<"]["<<l<<"] : " << params.t_l[j][l] << endl;
        }
    }
    return 0;
}

void initialize_x(int x[J_MAX][O_MAX][M_MAX])
{ 
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            for(int j=0; j<params.m; j++) {
                    x[i][k][j]=0;
            }
        }
    }
}

void GRASP_routing()
{
    //double W_tot=numeric_limits<double>::max(); // OFV : total workload
    for(int t=0; t<NB_ITERS; t++) {
        initialize_x(solutions[t].x);
        for(int i=0; i<params.n; i++) {
            for(int k=0; k<params.n_o[i]; k++) {
                //Restreined Candidates List of processing times (will be sorted)
                vector<double> RCL;
                vector<int> index_RCL;
                vector<double> proba;
                double totalBias=0;
                double mach_runif=0;
                double minW=numeric_limits<double>::max(); //to change ?
                double maxW=0;
                double range=0;
                double alpha=0;
                double width=0;
                for(int j=0; j<params.m; j++) {
                    if(params.A[i][k][j] == true && params.t[i][k][j] < minW) {
                        minW = params.t[i][k][j];
                    }
                    if(params.A[i][k][j] == true && params.t[i][k][j] > maxW) {
                        maxW = params.t[i][k][j];
                    }
                }
                range = maxW - minW;
                alpha = rand_gen.randDblExc(); //NO ?
                width = range * alpha;
                for(int j=0; j<params.m; j++) {
                    if(params.A[i][k][j] == true && params.t[i][k][j] <= minW + width) {
                        RCL.push_back(params.t[i][k][j]); // maybe costly
                        index_RCL.push_back(j);
                    }
                }
                sort(index_RCL.begin(),index_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
                sort(RCL.begin(), RCL.end()); //to test
                proba.resize(RCL.size(),0);
                for(unsigned int r=0; r<proba.size(); r++) { //r for rank
                    proba[r]=1/(r+1);
                    totalBias+=1/(r+1);
                }
                //~ for(unsigned int r=0; r<proba.size(); r++) {
                    //~ proba[r]/=totalBias;
                //~ }
                mach_runif = rand_gen.randDblExc(totalBias);
                for(unsigned int r=0; r<proba.size(); r++) {
                    if(mach_runif <= proba[r]) {
                        solutions[t].x[i][k][index_RCL[r]] = -1;
                        break;
                    } else {
                        mach_runif -= proba[r];
                    }
                }
            }
        }
        //Calculations of the OFV
        for(int j=0; j<params.m; j++) {
            double W_temp = 0;
            for(int i=0; i<params.n; i++) {
                for(int k=0; k<params.n_o[i]; k++) {
                    if(solutions[t].x[i][k][j] == -1) {
                        W_temp += params.t[i][k][j];
                        solutions[t].W_tot += params.t[i][k][j];
                    }
                }
            }
            for(int l=0; l<params.L[j]; l++) {
                W_temp += params.p[l][j];
                solutions[t].W_tot += params.p[l][j];
            }
            if(W_temp > solutions[t].W_max) {
                solutions[t].W_max = W_temp;
            }
        }
        
        
        cout << "At step " << t << " : " << " ";
        cout << "W total : " << solutions[t].W_tot << " ";
        cout << "W max : " << solutions[t].W_max << endl;
    }
}

int main()
{
    read_parameters("problem8x8.in");
    GRASP_routing();
    //~ for(int i=0; i<params.n; i++) {
        //~ for(int k=0; k<params.n_o[i]; k++) {
            //~ cout << "j" << i << " ";
            //~ cout << "t"<< k << " ";
            //~ for(int j=0; j<params.m; j++) {
                //~ if(vars.x[i][k][j] == true) {
                    //~ cout << "m" << j << " : " << params.t[i][k][j];
                //~ }
            //~ }
            //~ cout << endl;
        //~ }
    //~ }
    // decision varaibles -> c, y, x
    return 0;
}
