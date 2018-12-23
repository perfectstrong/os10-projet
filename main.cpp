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
#define NB_ITERS 10000
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
    //affectation of the operation ik to machine j
    int x[J_MAX][O_MAX];
    double W[M_MAX]; //Workload of each machine
    double W_tot=0; //Total workload (of all machines)
    double W_max=0; //Maximum Workload (of machines)
    double C_max=0; //Maximum completion time (of machines)
};

//index of an operation
struct Indexop{
    int i;
    int j;
    int k;
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

void initialize_x(int x[J_MAX][O_MAX])
{ 
    for(int i=0; i<params.n; i++) {
            for(int k=0; k<params.n_o[i]; k++) {
                    x[i][k]=-1;
            }
    }
}

void GRASP_routing()
{
    double W_tot_best = numeric_limits<double>::max();
    double W_max_best = numeric_limits<double>::max();
    //double W_tot=numeric_limits<double>::max(); // OFV : total workload
    for(int t=0; t<NB_ITERS; t++) {
        int nb_operations_to_assign = 0;
        for(int j=0; j<params.m; j++) {
            for(int l=0; l<params.L[j]; l++) {
                solutions[t].W[j] += params.p[l][j];
                solutions[t].W_tot += params.p[l][j];
            }
        }
        for(int i=0; i<params.n; i++) {
            nb_operations_to_assign += params.n_o[i];
        }
        initialize_x(solutions[t].x);
        while(nb_operations_to_assign > 0) {
            //Restreined Candidates List of processing times (will be sorted)
            vector<double> RCL;
            vector<int> index_i_RCL;
            vector<int> index_j_RCL;
            vector<int> index_k_RCL;
            vector<double> proba;
            double totalBias=0;
            double candidate_runif=0; // uniform varaible to choose the candidate
            double minW=numeric_limits<double>::max();
            double maxW=0;
            double range=0;
            double alpha=0;
            double width=0;
            
            //Determine the range
            for(int i=0; i<params.n; i++) {
                for(int k=0; k<params.n_o[i]; k++) {
                    if(solutions[t].x[i][k] == -1) {
                           for(int j=0; j<params.m; j++) {
                            if(params.A[i][k][j] == true && params.t[i][k][j] < minW) {
                                minW = params.t[i][k][j];
                            }
                            if(params.A[i][k][j] == true && params.t[i][k][j] > maxW) {
                                maxW = params.t[i][k][j];
                            }
                        } 
                    }
                }
            }
            range = maxW - minW;
            alpha = rand_gen.randDblExc(); //NO ?
            width = range * alpha;
            //Determine the candidates
            for(int i=0; i<params.n; i++) {
                for(int k=0; k<params.n_o[i]; k++) {
                    if(solutions[t].x[i][k] == -1) {
                        for(int j=0; j<params.m; j++) {
                            if(params.A[i][k][j] == true && params.t[i][k][j] <= minW + width) {
                                RCL.push_back(params.t[i][k][j]); // maybe costly
                                index_i_RCL.push_back(i);
                                index_j_RCL.push_back(j);
                                index_k_RCL.push_back(k);
                            }
                        }
                    }
                }
            } 
            sort(index_i_RCL.begin(),index_i_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(index_j_RCL.begin(),index_j_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(index_k_RCL.begin(),index_k_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(RCL.begin(), RCL.end()); //to test
            proba.resize(RCL.size(),0);
            for(unsigned int r=0; r<proba.size(); r++) { //r for rank
                    proba[r]=1/(r+1);
                    totalBias+=1/(r+1);
            }
            candidate_runif = rand_gen.randDblExc(totalBias);
            for(unsigned int r=0; r<proba.size(); r++) {
                if(candidate_runif <= proba[r]) {
                    int ind_i = index_i_RCL[r];
                    int ind_j = index_j_RCL[r];
                    int ind_k = index_k_RCL[r];
                    solutions[t].x[ind_i][ind_k] = ind_j;
                    solutions[t].W[ind_j] += params.t[ind_i][ind_k][ind_j];
                    solutions[t].W_tot += params.t[ind_i][ind_k][ind_j];
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            nb_operations_to_assign--;
        }
        //Calculations of the OFV
        for(int j=0; j<params.m; j++) {
            if(solutions[t].W[j] > solutions[t].W_max) {
                solutions[t].W_max = solutions[t].W[j];
            }
        }
        
        if(solutions[t].W_tot < W_tot_best) {
            W_tot_best = solutions[t].W_tot;
        }
        if(solutions[t].W_max < W_max_best) {
            W_max_best = solutions[t].W_max;
        }
        //~ cout << "At step " << t << " : " << " ";
        //~ cout << "W total : " << solutions[t].W_tot << " ";
        //~ cout << "W max : " << solutions[t].W_max << endl;
    }
    cout << "Conclusion :" << endl;
    cout << "W total best : " << W_tot_best << " ";
    cout << "W max best : " << W_max_best << endl;
}

void GRASP_scheduling() {
    for(int t=0; t<NB_ITERS; t++) {
        double machine_completion_time[M_MAX];
        double job_completion_time[J_MAX];
        int last_operation[J_MAX];
        int nb_operations_to_schedule = 0;
        for(int j=0; j<params.m; j++) {
            machine_completion_time[j] = 0;
        }
        for(int i=0; i<params.n; i++) {
            job_completion_time[i] = 0;
            last_operation[i] = 0;
            nb_operations_to_schedule += params.n_o[i];
        }
        while(nb_operations_to_schedule > 0) {
            
            
        }
    }
}

int main()
{
    read_parameters("problem8x8.in");
    GRASP_routing();
    // decision varaibles -> c, y, x
    return 0;
}
