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
    double p[M_MAX][L_MAX]; //processing time for maintenance task l on machine j
    double t_e[M_MAX][L_MAX]; //maintenance completion time windows begin (earliest)
    double t_l[M_MAX][L_MAX]; //maintenance completion time windows end (lateness)
};

struct Solution{
    //affectation of the operation ik to a machine
    int x[J_MAX][O_MAX];
    //completion time of maintenance task l on machine j
    double y[M_MAX][L_MAX];
    //completion time of operation ik on machine x[i][k]
    double c[J_MAX][O_MAX]; 
    double W[M_MAX]; //Workload of each machine
    double W_tot=0; //Total workload (of all machines)
    double W_max=0; //Maximum Workload (of machines)
    double C_j[M_MAX]; //Completion time of each machine
    double C_i[J_MAX]; //Completion time of each job
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
    //~ cout << "Got n : " << params.n << endl;
    file >> params.m;
    //~ cout << "Got m : " << params.m << endl;
    for(int i = 0; i<params.n; i++) {
       file >> params.n_o[i];
       //~ cout << "Got n_o["<<i<<"] : " << params.n_o[i] << endl;
    }
    for(int j = 0; j<params.m; j++) {
        file >> params.L[j];
        //~ cout << "Got L["<<j<<"] : " << params.L[j] << endl;
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
                //~ cout << "Got A["<<i<<"]["<<k<<"]["<<j<<"] : " << temp_word << endl;
            }
        }
    }
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            for(int j=0; j<params.m; j++)
            {
                file >> params.t[i][k][j];
                //~ cout << "Got t["<<i<<"]["<<k<<"]["<<j<<"] : " << params.t[i][k][j] << endl;
            }
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.p[j][l];
            //~ cout << "Got p["<<j<<"]["<<l<<"] : " << params.p[j][l] << endl;
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_e[j][l];
            //~ cout << "Got t_e["<<j<<"]["<<l<<"] : " << params.t_e[j][l] << endl;
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_l[j][l];
            //~ cout << "Got t_l["<<j<<"]["<<l<<"] : " << params.t_l[j][l] << endl;
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
                solutions[t].W[j] += params.p[j][l];
                solutions[t].W_tot += params.p[j][l];
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
            double candidate_runif=0; // uniform variable to choose the candidate
            double W_min=numeric_limits<double>::max();
            double W_max=0;
            double range=0;
            double alpha=0;
            double width=0;
            
            //Determine the range
            for(int i=0; i<params.n; i++) {
                for(int k=0; k<params.n_o[i]; k++) {
                    if(solutions[t].x[i][k] == -1) {
                           for(int j=0; j<params.m; j++) {
                            if(params.A[i][k][j] == true && params.t[i][k][j] < W_min) {
                                W_min = params.t[i][k][j];
                            }
                            if(params.A[i][k][j] == true && params.t[i][k][j] > W_max) {
                                W_max = params.t[i][k][j];
                            }
                        } 
                    }
                }
            }
            range = W_max - W_min;
            alpha = rand_gen.randDblExc(); //NO ?
            width = range * alpha;
            //Determine the candidates
            for(int i=0; i<params.n; i++) {
                for(int k=0; k<params.n_o[i]; k++) {
                    if(solutions[t].x[i][k] == -1) {
                        for(int j=0; j<params.m; j++) {
                            if(params.A[i][k][j] == true && params.t[i][k][j] <= W_min + width) {
                                RCL.push_back(params.t[i][k][j]); // maybe costly
                                index_i_RCL.push_back(i);
                                index_j_RCL.push_back(j);
                                index_k_RCL.push_back(k);
                            }
                        }
                    }
                }
            }
            //Sort the Candidate list 
            sort(index_i_RCL.begin(),index_i_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(index_j_RCL.begin(),index_j_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(index_k_RCL.begin(),index_k_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(RCL.begin(), RCL.end()); //to test
            proba.resize(RCL.size(),0);
            //Determine the probability of selection
            for(unsigned int r=0; r<proba.size(); r++) { //r for rank
                    proba[r]=1/(r+1);
                    totalBias+=1/(r+1);
            }
            //Select a candidate
            candidate_runif = rand_gen.randDblExc(totalBias);
            for(unsigned int r=0; r<proba.size(); r++) {
                if(candidate_runif <= proba[r]) {
                    int i = index_i_RCL[r];
                    int j = index_j_RCL[r];
                    int k = index_k_RCL[r];
                    solutions[t].x[i][k] = j;
                    solutions[t].W[j] += params.t[i][k][j];
                    solutions[t].W_tot += params.t[i][k][j];
                    if(solutions[t].W[j] > solutions[t].W_max) {
                        solutions[t].W_max = solutions[t].W[j];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            nb_operations_to_assign--;
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
    double C_max_best = numeric_limits<double>::max();
    for(int t=0; t<NB_ITERS; t++) {
        int last_operation[J_MAX];
        int nb_operations_to_schedule = 0;
        for(int j=0; j<params.m; j++) {
            solutions[t].C_j[j] = 0;
            for(int l=0; l<params.L[j]; l++) {
                // to indicate that is not scheduled yet
                solutions[t].y[j][l] = -1; 
            }
        }
        for(int i=0; i<params.n; i++) {
            for(int k=0; k<params.n_o[i]; k++) {
                // to indicate that is not scheduled yet
                solutions[t].c[i][k] = -1;
            }
            solutions[t].C_i[i] = 0;
            last_operation[i] = 0;
            nb_operations_to_schedule += params.n_o[i];
        }
        while(nb_operations_to_schedule > 0) {
            //Restreined Candidates List of completion times (will be sorted)
            vector<double> RCL;
            vector<int> index_i_RCL;
            vector<double> proba;
            double totalBias=0;
            double candidate_runif=0; // uniform variable to choose the candidate
            double C_max_temp[J_MAX]; //Temporary Completion times for range and candidates
            double C_max_min = numeric_limits<double>::max();
            double C_max_max = 0;
            double range=0;
            double alpha=0;
            double width=0;
            //Determine the range
            for(int i=0; i<params.n; i++) { 
                int k = last_operation[i]+1;
                if(k < params.n_o[i]) {
                    int j = solutions[t].x[i][k];
                    double begin_time=0;
                    //determine the most recent completion time
                    if(solutions[t].C_i[i] > solutions[t].C_j[j]) {
                        begin_time = solutions[t].C_i[i];
                    } else {
                        begin_time = solutions[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for(int l=0; l<params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if(solutions[t].y[j][l] < 0 && 
                        begin_time+params.t[i][k][j] > params.t_l[j][l]-params.p[j][l]) {
                            if(begin_time < params.t_e[j][l]) {
                                begin_time = params.t_e[j][l];
                            }
                            begin_time = begin_time + params.p[j][l];
                        }
                    }
                    if(begin_time + params.t[i][k][j] > solutions[t].C_max) {
                        C_max_temp[i] = begin_time + params.t[i][k][j];
                    } else {
                        C_max_temp[i] = solutions[t].C_max;
                    }
                    if(C_max_temp[i] < C_max_min) {
                        C_max_min = C_max_temp[i];
                    } else if (C_max_temp[i] > C_max_max) {
                        C_max_max = C_max_temp[i];
                    }
                }
            }
            range = C_max_max - C_max_min;
            alpha = rand_gen.randDblExc(); //NO ?
            width = range * alpha;
            //Determine the candidates
            for(int i=0; i<params.n; i++) {
                int k = last_operation[i]+1;
                int j = solutions[t].x[i][k];
                if(k < params.n_o[i] && C_max_temp[i] <= C_max_min + width) {
                    RCL.push_back(params.t[i][k][j]);
                    index_i_RCL.push_back(i);
                }
            } 
            //Sort the Candidate list 
            sort(index_i_RCL.begin(),index_i_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
            sort(RCL.begin(), RCL.end()); //to test
            proba.resize(RCL.size(),0);
            //Determine the probability of selection
            for(unsigned int r=0; r<proba.size(); r++) { //r for rank
                    proba[r]=1/(r+1);
                    totalBias+=1/(r+1);
            }
            //Select a candidate
            candidate_runif = rand_gen.randDblExc(totalBias);
            for(unsigned int r=0; r<proba.size(); r++) {
                if(candidate_runif <= proba[r]) {
                    int i = index_i_RCL[r];
                    int k = last_operation[i]+1;
                    int j = solutions[t].x[i][k];
                    double begin_time=0;
                    //determine the most recent completion time
                    if(solutions[t].C_i[i] > solutions[t].C_j[j]) {
                        begin_time = solutions[t].C_i[i];
                    } else {
                        begin_time = solutions[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for(int l=0; l<params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if(solutions[t].y[j][l] < 0 && //it can only happen on non-scheduled maintenance task
                        begin_time+params.t[i][k][j] > params.t_l[j][l]-params.p[j][l]) {
                            if(begin_time < params.t_e[j][l]) {
                                begin_time = params.t_e[j][l];
                            }
                            // and schedule it first
                            solutions[t].y[j][l] = begin_time + params.p[j][l];
                            begin_time = solutions[t].y[j][l];
                        }
                    }
                    //schedule the selected operation
                    solutions[t].c[i][k] = begin_time + params.t[i][k][j];
                    solutions[t].C_i[i] = solutions[t].c[i][k];
                    solutions[t].C_j[j] = solutions[t].c[i][k];
                    if(solutions[t].c[i][k] > solutions[t].C_max) {
                        solutions[t].C_max = solutions[t].c[i][k];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            if(solutions[t].C_max < C_max_best) {
                C_max_best = solutions[t].C_max;
            }
            nb_operations_to_schedule--;
        }
    }
    cout << " C max best : " << C_max_best << endl;
}

int main()
{
    read_parameters("problem8x8.in");
    GRASP_routing();
    GRASP_scheduling();
    // decision varaibles -> c, y, x
    return 0;
}
