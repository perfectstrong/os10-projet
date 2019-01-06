#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include "MersenneTwister.h"
#include <ctime>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define MAX_TIME 1000
#define J_MAX 100 //MAX number of jobs
#define O_MAX 10 //MAX number of operations per jobs
#define M_MAX 100 //MAX number of machines
#define L_MAX 10 //MAX number of maintenance tasks per machine
#define NB_ITERS 1000
#define SEED 2019

struct Parameters {
    int n = 0; //number of jobs
    int m = 0; //number of machines
    int n_o[J_MAX]; //number of operations of jobs i 
    int L[M_MAX]; //number of maintenance tasks on machine j
    bool A[J_MAX][O_MAX][M_MAX]; //availables machines for task k of job i
    double t[J_MAX][O_MAX][M_MAX]; //processing time for task k of job i on machine j
    //CAUTION : Maintenance task have tobe given in chronological order
    double p[M_MAX][L_MAX]; //processing time for maintenance task l on machine j
    double t_e[M_MAX][L_MAX]; //maintenance completion time windows begin (earliest)
    double t_l[M_MAX][L_MAX]; //maintenance completion time windows end (lateness)
};

struct Solution {
    // decision varaibles -> x, y, c
    //affectation of the operation ik to a machine
    int x[J_MAX][O_MAX];
    //completion time of maintenance task l on machine j
    double y[M_MAX][L_MAX];
    //completion time of operation ik on machine x[i][k]
    double c[J_MAX][O_MAX];
    double W[M_MAX]; //Workload of each machine
    double W_tot = 0; //Total workload (of all machines)
    double W_max = 0; //Maximum Workload (of machines)
    double C_j[M_MAX]; //Completion time of each machine
    double C_i[J_MAX]; //Completion time of each job
    double C_max = 0; //Maximum completion time (of machines)
};

struct candidate {
    int i; //job
    int k; //operation
    int j; //machine
    double obj_value; //workload or completion time
};

//static variables
static MTRand rand_gen(SEED);
static Parameters params;
static Solution solution[NB_ITERS]; //indexed by t
static bool parentoFront[NB_ITERS]; // FALSE = dominated

//Give  the current time of the clock
high_resolution_clock::time_point now() {
    return high_resolution_clock::now();
}

//to compare 2 candidates in the RCL list for sorting
bool comparisonCandidates(const candidate& candidate1, const candidate& candidate2) {
    if (candidate1.obj_value < candidate2.obj_value) {
        return true;
    } else if (candidate1.obj_value == candidate2.obj_value) {
        if (candidate1.i < candidate2.i) {
            return true;
        } else if (candidate1.i == candidate2.i) {
            if (candidate1.k < candidate2.k) {
                return true;
            } else if (candidate1.k == candidate2.k) {
                if (candidate1.j < candidate2.j) {
                    return true;
                }
            }
        }
    }
    return false;
}

//read the FJSSP-nfa parameters from a file
void read_parameters(string filename) {
    ifstream file(filename);
    if (file.fail()) {
        cout << "Error reading : file \" " << filename << " \" not found" << endl;
        getchar();
        exit(1);
    }
    string word;
    file >> word;
    //~ cout << word << endl;
    file >> params.n;
    //~ cout << "Got n : " << params.n << endl;
    file >> params.m;
    //~ cout << "Got m : " << params.m << endl;
    for (int i = 0; i < params.n; i++) {
        file >> params.n_o[i];
        //~ cout << "Got n_o["<<i<<"] : " << params.n_o[i] << endl;
    }
    for (int j = 0; j < params.m; j++) {
        file >> params.L[j];
        //~ cout << "Got L["<<j<<"] : " << params.L[j] << endl;
    }
    for (int i = 0; i < params.n; i++) {
        for (int k = 0; k < params.n_o[i]; k++) {
            for (int j = 0; j < params.m; j++) {
                string temp_word;
                file >> temp_word;
                if (temp_word.find('T') != std::string::npos) {
                    params.A[i][k][j] = true;
                } else {
                    params.A[i][k][j] = false;
                }
                //~ cout << "Got A["<<i<<"]["<<k<<"]["<<j<<"] : " << temp_word << endl;
            }
        }
    }
    for (int i = 0; i < params.n; i++) {
        for (int k = 0; k < params.n_o[i]; k++) {
            for (int j = 0; j < params.m; j++) {
                file >> params.t[i][k][j];
                //~ cout << "Got t["<<i<<"]["<<k<<"]["<<j<<"] : " << params.t[i][k][j] << endl;
            }
        }
    }
    for (int j = 0; j < params.m; j++) {
        for (int l = 0; l < params.L[j]; l++) {
            file >> params.p[j][l];
            //~ cout << "Got p["<<j<<"]["<<l<<"] : " << params.p[j][l] << endl;
        }
    }
    for (int j = 0; j < params.m; j++) {
        for (int l = 0; l < params.L[j]; l++) {
            file >> params.t_e[j][l];
            //~ cout << "Got t_e["<<j<<"]["<<l<<"] : " << params.t_e[j][l] << endl;
        }
    }
    for (int j = 0; j < params.m; j++) {
        for (int l = 0; l < params.L[j]; l++) {
            file >> params.t_l[j][l];
            //~ cout << "Got t_l["<<j<<"]["<<l<<"] : " << params.t_l[j][l] << endl;
        }
    }
}

//initialize the affectation
void initialize_x(int x[J_MAX][O_MAX]) {
    for (int i = 0; i < params.n; i++) {
        for (int k = 0; k < params.n_o[i]; k++) {
            x[i][k] = -1;
        }
    }
}

//Routing heuristic
void GRASP_routing() {
    //construct NB_ITERS solutions indexed by t
    for (int t = 0; t < NB_ITERS; t++) {
        //Initializations
        int nb_operations_to_assign = 0;
        for (int j = 0; j < params.m; j++) {
            for (int l = 0; l < params.L[j]; l++) {
                solution[t].W[j] += params.p[j][l];
                solution[t].W_tot += params.p[j][l];
            }
        }
        for (int i = 0; i < params.n; i++) {
            nb_operations_to_assign += params.n_o[i];
        }
        initialize_x(solution[t].x);
        //Restreined Candidates List of processing times (will be sorted)
        vector<candidate> RCL;
        vector<double> proba;
        double totalBias;
        double candidate_runif; // uniform variable to choose the candidate
        double W_min;
        double W_max;
        double range;
        double alpha = rand_gen.randDblExc(); // Fixed for each solution (iteration)
        double width;
        //Affectations
        while (nb_operations_to_assign > 0) {
            RCL.clear();
            proba.clear();
            totalBias = 0;
            candidate_runif = 0;
            W_min = numeric_limits<double>::max();
            W_max = 0;
            range = 0;
            width = 0;
            //Determine the range
            for (int i = 0; i < params.n; i++) {
                for (int k = 0; k < params.n_o[i]; k++) {
                    if (solution[t].x[i][k] == -1) {
                        for (int j = 0; j < params.m; j++) {
                            if (params.A[i][k][j] && params.t[i][k][j] < W_min) {
                                W_min = params.t[i][k][j];
                            }
                            if (params.A[i][k][j] && params.t[i][k][j] > W_max) {
                                W_max = params.t[i][k][j];
                            }
                        }
                    }
                }
            }
            range = W_max - W_min;
            width = range * alpha;
            //Determine the candidates
            for (int i = 0; i < params.n; i++) {
                for (int k = 0; k < params.n_o[i]; k++) {
                    if (solution[t].x[i][k] == -1) {
                        for (int j = 0; j < params.m; j++) {
                            if (params.A[i][k][j]
                                // && params.t[i][k][j] >= W_min // always true
                                && params.t[i][k][j] <= W_min + width) {
                                RCL.push_back(candidate());
                                RCL.back().i = i;
                                RCL.back().k = k;
                                RCL.back().j = j;
                                RCL.back().obj_value = params.t[i][k][j];
                            }
                        }
                    }
                }
            }
            //Sort the Candidate list 
            sort(RCL.begin(), RCL.end(), comparisonCandidates);
            proba.resize(RCL.size(), 0);
            //Determine the probability of selection
            for (unsigned int r = 0; r < proba.size(); r++) { //r for rank
                proba.at(r) = (double) 1 / (r + 1);
                totalBias += (double) 1 / (r + 1);
            }
            //Select and assign a candidate
            candidate_runif = rand_gen.randDblExc(totalBias);
            for (unsigned int r = 0; r < proba.size(); r++) {
                if (candidate_runif <= proba[r]) {
                    int i = RCL[r].i;
                    int k = RCL[r].k;
                    int j = RCL[r].j;
                    solution[t].x[i][k] = j;
                    solution[t].W[j] += params.t[i][k][j];
                    solution[t].W_tot += params.t[i][k][j];
                    if (solution[t].W[j] > solution[t].W_max) {
                        solution[t].W_max = solution[t].W[j];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            nb_operations_to_assign--;
        }
    }
}

void GRASP_scheduling() {
    //construct NB_ITERS solutions indexed by t
    for (int t = 0; t < NB_ITERS; t++) {
        //Initializations
        int current_operation[J_MAX];
        int nb_operations_to_schedule = 0;
        for (int j = 0; j < params.m; j++) {
            solution[t].C_j[j] = 0;
            for (int l = 0; l < params.L[j]; l++) {
                // to indicate that is not scheduled yet
                solution[t].y[j][l] = -1;
            }
        }
        for (int i = 0; i < params.n; i++) {
            for (int k = 0; k < params.n_o[i]; k++) {
                // to indicate that is not scheduled yet
                solution[t].c[i][k] = -1;
            }
            solution[t].C_i[i] = 0;
            current_operation[i] = 0;
            nb_operations_to_schedule += params.n_o[i];
        }
        solution[t].C_max = 0;
        //Restreined Candidates List of completion times (will be sorted)
        vector<candidate> RCL;
        vector<double> proba;
        double totalBias;
        double candidate_runif; // uniform variable to choose the candidate
        double C_max_temp[J_MAX]; //Temporary Completion times for range and candidates
        double C_max_min;
        double C_max_max;
        double range;
        double alpha = rand_gen.randDblExc(); // Fixed for each solution (iteration)
        double width;
        //Scheduling
        while (nb_operations_to_schedule > 0) {
            RCL.clear();
            proba.clear();
            totalBias = 0;
            candidate_runif = 0;
            C_max_min = numeric_limits<double>::max();
            C_max_max = 0;
            range = 0;
            width = 0;
            //Determine the range
            for (int i = 0; i < params.n; i++) {
                int k = current_operation[i];
                if (k < params.n_o[i]) {
                    int j = solution[t].x[i][k];
                    double begin_time = 0;
                    //determine the most recent completion time
                    if (solution[t].C_i[i] > solution[t].C_j[j]) {
                        begin_time = solution[t].C_i[i];
                    } else {
                        begin_time = solution[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for (int l = 0; l < params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if (solution[t].y[j][l] < 0 &&
                            begin_time + params.t[i][k][j] > params.t_l[j][l] - params.p[j][l]) {
                            double PM_begin_time = solution[t].C_j[j];
                            if (PM_begin_time < params.t_e[j][l] - params.p[j][l]) {
                                PM_begin_time = params.t_e[j][l] - params.p[j][l];
                            }
                            if (begin_time < PM_begin_time + params.p[j][l]) {
                                begin_time = PM_begin_time + params.p[j][l];
                            }
                        }
                    }
                    if (begin_time + params.t[i][k][j] > solution[t].C_max) {
                        C_max_temp[i] = begin_time + params.t[i][k][j];
                    } else {
                        C_max_temp[i] = solution[t].C_max;
                    }
                    if (C_max_temp[i] < C_max_min) {
                        C_max_min = C_max_temp[i];
                    }
                    if (C_max_temp[i] > C_max_max) {
                        C_max_max = C_max_temp[i];
                    }
                }
            }
            range = C_max_max - C_max_min;
            width = 1 + range * alpha;
            //Determine the candidates
            for (int i = 0; i < params.n; i++) {
                int k = current_operation[i];
                int j = solution[t].x[i][k];
                if (k < params.n_o[i] 
                    // && C_max_temp[i] >= C_max_min // always true
                    && C_max_temp[i] <= C_max_min + width) {
                    RCL.push_back(candidate());
                    RCL.back().i = i;
                    RCL.back().k = k;
                    RCL.back().j = j;
                    RCL.back().obj_value = params.t[i][k][j];
                }
            }
            //Sort the Candidate list 
            sort(RCL.begin(), RCL.end(), comparisonCandidates);
            proba.resize(RCL.size(), 0);
            //Determine the probability of selection
            for (unsigned int r = 0; r < proba.size(); r++) { //r for rank
                proba.at(r) = (double) 1 / (r + 1);
                totalBias += (double) 1 / (r + 1);
            }
            //Select a candidate
            candidate_runif = rand_gen.randDblExc(totalBias);
            for (unsigned int r = 0; r < proba.size(); r++) {
                if (candidate_runif <= proba[r]) {
                    int i = RCL[r].i;
                    int k = RCL[r].k;
                    int j = RCL[r].j;
                    double begin_time = 0;
                    //determine the most recent completion time
                    if (solution[t].C_i[i] > solution[t].C_j[j]) {
                        begin_time = solution[t].C_i[i];
                    } else {
                        begin_time = solution[t].C_j[j];
                    }
                    //checking if no conflict with a maintenance task
                    for (int l = 0; l < params.L[j]; l++) {
                        //if conflict, take into account in begin time
                        if (solution[t].y[j][l] < 0 && //it can only happen on non-scheduled maintenance task
                            begin_time + params.t[i][k][j] > params.t_l[j][l] - params.p[j][l]) {
                            double PM_begin_time = solution[t].C_j[j];
                            if (PM_begin_time < params.t_e[j][l] - params.p[j][l]) {
                                PM_begin_time = params.t_e[j][l] - params.p[j][l];
                            }
                            // and schedule it first
                            solution[t].y[j][l] = PM_begin_time + params.p[j][l];
                            if (begin_time < solution[t].y[j][l]) {
                                begin_time = solution[t].y[j][l];
                            }
                        }
                    }
                    //schedule the selected operation
                    solution[t].c[i][k] = begin_time + params.t[i][k][j];
                    current_operation[i]++;
                    solution[t].C_i[i] = solution[t].c[i][k];
                    solution[t].C_j[j] = solution[t].c[i][k];
                    if (solution[t].c[i][k] > solution[t].C_max) {
                        solution[t].C_max = solution[t].c[i][k];
                    }
                    break;
                } else {
                    candidate_runif -= proba[r];
                }
            }
            nb_operations_to_schedule--;
        }
        //CAUTION : Maintenance task have tobe given in chronological order
        //Scheduling the last maintenance task
        for (int j = 0; j < params.m; j++) {
            double begin_time = solution[t].C_j[j];
            for (int l = 0; l < params.L[j]; l++) {
                if (solution[t].y[j][l] < 0) {
                    if (begin_time < params.t_e[j][l] - params.p[j][l]) {
                        begin_time = params.t_e[j][l] - params.p[j][l];
                    }
                    // and schedule it first
                    solution[t].y[j][l] = begin_time + params.p[j][l];
                    begin_time = solution[t].y[j][l];
                }
            }
            solution[t].C_j[j] = begin_time;
            if (solution[t].C_j[j] > solution[t].C_max) {
                solution[t].C_max = solution[t].C_j[j];
            }
        }
    }
}

//Return the best solution according to the minmization of the
//objective function which is a linear combination of W_tot, W_max and C_max
//parameters are the factors of the linear combination and
//return is the index t
int getBestSolution(double F_W_tot, double F_W_max, double F_C_max) {
    int t_best = 0;
    //best OFV : Objective Function Value
    double OFV_best = numeric_limits<double>::max();
    for (int t = 0; t < NB_ITERS; t++) {
        double OFV = F_W_tot * solution[t].W_tot + F_W_max * solution[t].W_max
            + F_C_max * solution[t].C_max;
        if (OFV < OFV_best) {
            t_best = t;
            OFV_best = OFV;
        }
    }
    return t_best;
}

//Write the best solution according to F_W_tot, F_W_max and F_C_max
//on a file
void writeOutSolution(double F_W_tot, double F_W_max, double F_C_max, string filename) {
    //Declaration
    int t = getBestSolution(F_W_tot, F_W_max, F_C_max);
    double OFV = F_W_tot * solution[t].W_tot + F_W_max * solution[t].W_max
        + F_C_max * solution[t].C_max;
    int nb_operations = 0;
    int nb_maintenances = 0;
    ofstream file(filename);
    if (file.fail()) {
        cout << "Error writing : file \" " << filename << " \" not found" << endl;
        getchar();
        exit(1);
    }
    //init
    for (int i = 0; i < params.n; i++) {
        nb_operations += params.n_o[i];
    }
    for (int j = 0; j < params.m; j++) {
        nb_maintenances += params.L[j];
    }
    file << nb_operations << " " << nb_maintenances << endl;
    file << F_W_tot << " " << F_W_max << " "
        << F_C_max << " " << OFV << endl;
    file << " " << solution[t].W_tot << " " << solution[t].W_max << " "
        << solution[t].C_max << endl;
    for (int i = 0; i < params.n; i++) {
        for (int k = 0; k < params.n_o[i]; k++) {
            int j = solution[t].x[i][k];
            file << i << " " << k << " " << j << " ";
            file << solution[t].c[i][k] - params.t[i][k][j] << " ";
            file << solution[t].c[i][k] << endl;
        }
    }
    for (int j = 0; j < params.m; j++) {
        for (int l = 0; l < params.L[j]; l++) {
            file << j << " " << l << " ";
            file << solution[t].y[j][l] - params.p[j][l] << " ";
            file << solution[t].y[j][l] << endl;
        }
    }
}

bool dominate(Solution s1, Solution s2) {
    // Check if s1 dominates s2
    return (s1.C_max <= s2.C_max
        && s1.W_max <= s2.W_max
        && s1.W_tot <= s2.W_tot
        && (s1.C_max < s2.C_max
            || s1.W_max < s2.W_max
            || s1.W_tot < s2.W_tot));
}

bool idem(Solution s1, Solution s2) {
    return s1.C_max == s2.C_max
        && s1.W_max == s2.W_max
        && s1.W_tot == s2.W_tot;
}

void writeOutParetoFront(string output) {
    ofstream f(output);
    if (f.fail()) {
        cout << "Error writing : file \" " << output << " \" not found" << endl;
        getchar();
        exit(1);
    }
    // Brutal check indomination
    for (int i = 0; i < NB_ITERS; i++) {
        parentoFront[i] = true;
        for (int j = 0; j < NB_ITERS; j++) {
            if (j == i) continue;
            if (idem(solution[j], solution[i])
                && parentoFront[j]) {
                parentoFront[i] = false;
                break;
            }
            if (dominate(solution[j], solution[i])) {
                // If solution[j] dominates solution[i]
                parentoFront[i] = false;
                break;
            }
        }
        if (parentoFront[i]) {
            // Write
            f << solution[i].W_tot << " " << solution[i].W_max << " "
                << solution[i].C_max << "\n";
        }
    }
}

void doJob(string input) {
    read_parameters(input);
    high_resolution_clock::time_point start, finish;
    duration<double, std::milli> timespan;
    start = now();
    GRASP_routing();
    GRASP_scheduling();
    finish = now();
    string prefix = string("./output/results") + to_string(params.n) + "x" + to_string(params.m);
    writeOutSolution(1, 0, 0, prefix + "_Wtot.out");
    writeOutSolution(0, 1, 0, prefix + "_Wmax.out");
    writeOutSolution(0, 0, 1, prefix + "_Cmax.out");
    writeOutSolution(0.5, 0.3, 0.2, prefix + "_F050302.out");
    writeOutSolution(0.5, 0.2, 0.3, prefix + "_F050203.out");
    writeOutParetoFront(prefix + "_ParetoFront.out");
    timespan = finish - start;
    cout << "Execution time = " << timespan.count() << " ms\n";
    //~ getchar();
}

int main() {
    cout << "Input=";
    string input;
    getline(cin, input);
    doJob(input);
    return 0;
}
