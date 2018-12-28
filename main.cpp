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
#define NB_ITERS 100
#define SEED 2018

struct Parameters{
    int n=0; //number of jobs
    int m=0; //number of machines
    int n_o[J_MAX]; //number of operations of jobs i 
    int L[M_MAX]; //number of maintenance tasks on machine j
    bool A[J_MAX][O_MAX][M_MAX]; //availables machines for task k of job i
    double t[J_MAX][O_MAX][M_MAX]; //processing time for task k of job i on machine j
    //CAUTION : Maintenance task have tobe given in chronological order
    double p[M_MAX][L_MAX]; //processing time for maintenance task l on machine j
    double t_e[M_MAX][L_MAX]; //maintenance completion time windows begin (earliest)
    double t_l[M_MAX][L_MAX]; //maintenance completion time windows end (lateness)
};

struct Solution{
    // decision varaibles -> x, y, c
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

struct candidate 
{
    int i; //job
    int k; //operation
    int j; //machine
    double obj_value; //workload or completion time
};

//static variables
static MTRand rand_gen(SEED);
static Parameters params;
static Solution solution[NB_ITERS];

bool comparisonCandidates(const candidate& candidate1, const candidate& candidate2)
{
    if (candidate1.obj_value < candidate2.obj_value) {
        return true;
    } else if (candidate1.obj_value == candidate2.obj_value) {
        if(candidate1.i < candidate2.i) {
            return true;
        } else if(candidate1.i == candidate2.i) {
            if(candidate1.k < candidate2.k) {
                return true;
            } else if(candidate1.k == candidate2.k) {
                if(candidate1.j < candidate2.j) {
                    return true;
                }
            }
        }
    }
    return false;
}

int read_parameters(string filename)
{
    ifstream file(filename);
    if(file.fail()) {
        cout << "Error reading : file \" "<< filename <<" \" not found" << endl;
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
        //~ cout << "iteration " << t << endl;
        int nb_operations_to_assign = 0;
        for(int j=0; j<params.m; j++) {
            for(int l=0; l<params.L[j]; l++) {
                solution[t].W[j] += params.p[j][l];
                solution[t].W_tot += params.p[j][l];
            }
        }
        for(int i=0; i<params.n; i++) {
            nb_operations_to_assign += params.n_o[i];
        }
        initialize_x(solution[t].x);
        while(nb_operations_to_assign > 0) {
            //Restreined Candidates List of processing times (will be sorted)
            vector<candidate> RCL;
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
                    if(solution[t].x[i][k] == -1) {
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
                    if(solution[t].x[i][k] == -1) {
                        for(int j=0; j<params.m; j++) {
                            if(params.A[i][k][j] == true && params.t[i][k][j] <= W_min + width) {
                                //~ candidate temp_candidate = 
                                //~ {i, k, j, params.t[i][k][j]};
                                RCL.push_back(candidate());
                                RCL.back().i=i;
                                RCL.back().k=k;
                                RCL.back().j=j;
                                RCL.back().obj_value=params.t[i][k][j];
                            }
                        }
                    }
                }
            }
            //Sort the Candidate list 
            sort(RCL.begin(), RCL.end(), comparisonCandidates); //to test
            //~ if(t == 10) {
                //~ cout << "candidates after :" << endl;
                //~ for(unsigned int r=0; r<RCL.size(); r++) {
                    //~ cout << RCL[r].i << " ";
                    //~ cout << RCL[r].k << " ";
                    //~ cout << RCL[r].j << " ";
                    //~ cout << RCL[r].obj_value << " ";
                    //~ cout << endl;
                //~ }
            //~ }
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
                    int i = RCL[r].i;
                    int k = RCL[r].k;
                    int j = RCL[r].j;
                    //~ if(t==10) {
                        //~ cout << "sélection : " << i << " " << k ;
                        //~ cout << " pour : " << j << endl;
                    //~ }
                    solution[t].x[i][k] = j;
                    solution[t].W[j] += params.t[i][k][j];
                    solution[t].W_tot += params.t[i][k][j];
                    if(solution[t].W[j] > solution[t].W_max) {
                        solution[t].W_max = solution[t].W[j];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            nb_operations_to_assign--;
        }
        
        if(solution[t].W_tot < W_tot_best) {
            W_tot_best = solution[t].W_tot;
        }
        if(solution[t].W_max < W_max_best) {
            W_max_best = solution[t].W_max;
        }
        //~ cout << "At step " << t << " : " << " ";
        //~ cout << "W total : " << solution[t].W_tot << " ";
        //~ cout << "W max : " << solution[t].W_max << endl;
    }
    cout << "Conclusion :" << endl;
    cout << "W total best : " << W_tot_best << " ";
    cout << "W max best : " << W_max_best << endl;
}

void GRASP_scheduling() {
    double C_max_best = numeric_limits<double>::max();
    for(int t=0; t<NB_ITERS; t++) {
        int current_operation[J_MAX];
        int nb_operations_to_schedule = 0;
        for(int j=0; j<params.m; j++) {
            solution[t].C_j[j] = 0;
            for(int l=0; l<params.L[j]; l++) {
                // to indicate that is not scheduled yet
                solution[t].y[j][l] = -1; 
            }
        }
        for(int i=0; i<params.n; i++) {
            for(int k=0; k<params.n_o[i]; k++) {
                // to indicate that is not scheduled yet
                solution[t].c[i][k] = -1;
            }
            solution[t].C_i[i] = 0;
            current_operation[i] = 0;
            nb_operations_to_schedule += params.n_o[i];
        }
        while(nb_operations_to_schedule > 0) {
            //Restreined Candidates List of completion times (will be sorted)
            vector<candidate> RCL;
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
                int k = current_operation[i];
                if(k < params.n_o[i]) {
                    int j = solution[t].x[i][k];
                    double begin_time=0;
                    //determine the most recent completion time
                    if(solution[t].C_i[i] > solution[t].C_j[j]) {
                        begin_time = solution[t].C_i[i];
                    } else {
                        begin_time = solution[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for(int l=0; l<params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if(solution[t].y[j][l] < 0 && 
                        begin_time+params.t[i][k][j] > params.t_l[j][l]-params.p[j][l]) {
                            if(begin_time < params.t_e[j][l]) {
                                begin_time = params.t_e[j][l];
                            }
                            begin_time = begin_time + params.p[j][l];
                        }
                    }
                    if(begin_time + params.t[i][k][j] > solution[t].C_max) {
                        C_max_temp[i] = begin_time + params.t[i][k][j];
                    } else {
                        C_max_temp[i] = solution[t].C_max;
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
                int k = current_operation[i];
                int j = solution[t].x[i][k];
                if(k < params.n_o[i] && C_max_temp[i] <= C_max_min + width) {
                    RCL.push_back(candidate());
                    RCL.back().i=i;
                    RCL.back().k=k;
                    RCL.back().j=j;
                    RCL.back().obj_value=params.t[i][k][j];
                }
            } 
            //Sort the Candidate list 
            sort(RCL.begin(), RCL.end(), comparisonCandidates); 
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
                    int i = RCL[r].i;
                    int k = RCL[r].k;
                    int j = RCL[r].j;
                    double begin_time=0;
                    //determine the most recent completion time
                    if(solution[t].C_i[i] > solution[t].C_j[j]) {
                        begin_time = solution[t].C_i[i];
                    } else {
                        begin_time = solution[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for(int l=0; l<params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if(solution[t].y[j][l] < 0 && //it can only happen on non-scheduled maintenance task
                        begin_time+params.t[i][k][j] > params.t_l[j][l]-params.p[j][l]) {
                            if(begin_time < params.t_e[j][l]) {
                                begin_time = params.t_e[j][l];
                            }
                            // and schedule it first
                            solution[t].y[j][l] = begin_time + params.p[j][l];
                            begin_time = solution[t].y[j][l];
                        }
                    }
                    //schedule the selected operation
                    solution[t].c[i][k] = begin_time + params.t[i][k][j];
                    current_operation[i]++;
                    solution[t].C_i[i] = solution[t].c[i][k];
                    solution[t].C_j[j] = solution[t].c[i][k];
                    if(solution[t].c[i][k] > solution[t].C_max) {
                        solution[t].C_max = solution[t].c[i][k];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            //CAUTION : Maintenance task have tobe given in chronological order
            //Scheduling the last maintenance task
            for(int j=0; j<params.m; j++) {
                double begin_time = solution[t].C_j[j];
                for(int l=0; l<params.L[j]; l++) {
                    if(solution[t].y[j][l] < 0) {
                        if(begin_time < params.t_e[j][l]) {
                            begin_time = params.t_e[j][l];
                        }
                        // and schedule it first
                        solution[t].y[j][l] = begin_time + params.p[j][l];
                        begin_time = solution[t].y[j][l];
                    }
                }
                solution[t].C_j[j] = begin_time;
            }
            
            if(solution[t].C_max < C_max_best) {
                C_max_best = solution[t].C_max;
            }
            nb_operations_to_schedule--;
        }
    }
    cout << " C max best : " << C_max_best << endl;
}

int writeOutSolution(int t, string filename) {
    ofstream file(filename);
    if(file.fail()) {
        cout << "Error writing : file \" "<< filename <<" \" not found" << endl;
        return 1;
    }
    //add more stats ?
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
                int j = solution[t].x[i][k];
                file << "job " << i << " " << k << " " << j << " ";
                file << solution[t].c[i][k]-params.t[i][k][j] << " ";
                file << solution[t].c[i][k] << endl;
        }
    }
    for(int j=0; j<params.m; j++) {
        for(int l=0; l<params.L[j]; l++) {
            file << "maintenance " << j << " " << l << " ";
            file << solution[t].y[j][l]-params.p[j][l] << " ";
            file << solution[t].y[j][l] << endl;
        }
    }
    return 0;
}

int main()
{
    read_parameters("./input/problem8x8.in");
    GRASP_routing();
    GRASP_scheduling();
    writeOutSolution(10, "./output/results8x8.out");
    //~ {
        //~ vector<candidate> toto;
        //~ toto.push_back(candidate());
        //~ toto.back().i=1;
    //~ }
    return 0;
}
