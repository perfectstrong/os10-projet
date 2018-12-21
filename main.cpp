#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include "MersenneTwister.h"

using namespace std;

#define MAX_TIME 1000

struct parameters{
    int n=0; //number of jobs
    int m=0; //number of machines
    int* n_o; //number of operations of jobs i 
    int* L; //number of maintenance tasks on machine j
    bool*** A; //availables machines for task k of job i
    double*** t; //processing time for task k of job i on machine j
    double** p; //processing time for maintenance task l on machine j
    double** t_e; //maintenance completion time windows begin (earliest)
    double** t_l; //maintenance completion time windows end (lateness)
};

//~ struct variables{
    //~ double** y; //completion time of of maintenance task jl
    //~ bool*** x; //affectation of the operation ik to machine j
//~ };

int read_parameters(string filename, parameters& params)
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
    params.n_o = new int[params.n];
    for(int i = 0; i<params.n; i++) {
       file >> params.n_o[i];
       cout << "Got n_o["<<i<<"] : " << params.n_o[i] << endl;
    }
    params.L = new int[params.m];
    for(int j = 0; j<params.m; j++) {
        file >> params.L[j];
        cout << "Got L["<<j<<"] : " << params.L[j] << endl;
    }
    params.A = new bool**[params.n];
    for(int i=0; i<params.n; i++) {
        params.A[i] = new bool*[params.n_o[i]];
        for(int k=0; k<params.n_o[i]; k++) {
            params.A[i][k] = new bool[params.m];
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
    params.t = new double**[params.n];
    for(int i=0; i<params.n; i++) {
        params.t[i] = new double*[params.n_o[i]];
        for(int k=0; k<params.n_o[i]; k++) {
            params.t[i][k] = new double[params.m];
            for(int j=0; j<params.m; j++)
            {
                file >> params.t[i][k][j];
                cout << "Got t["<<i<<"]["<<k<<"]["<<j<<"] : " << params.t[i][k][j] << endl;
            }
        }
    }
    params.p = new double*[params.m];
    for(int j=0; j<params.m; j++) {
        params.p[j] = new double[params.L[j]];
        for(int l=0; l<params.L[j]; l++) {
            file >> params.p[j][l];
            cout << "Got p["<<j<<"]["<<l<<"] : " << params.p[j][l] << endl;
        }
    }
    params.t_e = new double*[params.m];
    for(int j=0; j<params.m; j++) {
        params.t_e[j] = new double[params.L[j]];
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_e[j][l];
            cout << "Got t_e["<<j<<"]["<<l<<"] : " << params.t_e[j][l] << endl;
        }
    }
    params.t_l = new double*[params.m];
    for(int j=0; j<params.m; j++) {
        params.t_l[j] = new double[params.L[j]];
        for(int l=0; l<params.L[j]; l++) {
            file >> params.t_l[j][l];
            cout << "Got t_l["<<j<<"]["<<l<<"] : " << params.t_l[j][l] << endl;
        }
    }
    return 0;
}

void delete_parameters(parameters& params)
{
    free(params.n_o);
    free(params.L);
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            free(params.A[i][k]);
        }
        free(params.A[i]);
    }
    free(params.A);
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            free(params.t[i][k]);
        }
        free(params.t[i]);
    }
    free(params.t);
    for(int j=0; j<params.m; j++) {
        free(params.p[j]);
    }
    free(params.p);
    for(int j=0; j<params.m; j++) {
        free(params.t_e[j]);
    }
    free(params.t_e);
    for(int j=0; j<params.m; j++) {
        free(params.t_l[j]);
    }
    free(params.t_l);
}

    //~ vars.y = new double*[params.m];
    //~ for(int j=0; j<params.n; j++) {
        //~ vars.y[j] = new double[params.L[j]];
        //~ for(int l=0; l<params.L[j]; l++) {
            //~ vars.y[j][l] = 0;
        //~ }
    //~ }

void initialize_x(bool*** x, parameters& params)
{
    x = new bool**[params.n]; 
    for(int i=0; i<params.n; i++) {
        x[i] = new bool*[params.n_o[i]];
        for(int k=0; k<params.n_o[i]; k++) {
            x[i][k] = new bool[params.m];
            for(int j=0; j<params.m; j++) {
                    x[i][k][j]=false;
            }
        }
    }
}

//~ for(int j=0; j<params.n; j++) {
        //~ free(vars.y[j]);
    //~ }
    //~ free(vars.y);

void delete_x(bool*** x, parameters& params)
{
    for(int i=0; i<params.n; i++) {
        for(int k=0; k<params.n_o[i]; k++) {
            free(x[i][k]);
        }
        free(x[i]);
    }
    free(x);
}

bool*** GRASP_routing(MTRand rand_gen, int* alpha, int nb_iters, parameters& params)
{
    bool*** x_sol=NULL;
    double W_tot=numeric_limits<double>::max(); // OFV : total workload
    //bool
    for(int t=0; t<nb_iters; t++) {
        bool*** x_temp;
        double W_tot_temp=0;
        initialize_x(x_temp,params);
        for(int i=0; i<params.n; i++) {
            for(int k=0; k<params.n_o[i]; k++) {
                //Restreined Candidates List (will be sorted)
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
                    if(params.t[i][k][j] < minW) {
                        minW = params.t[i][k][j];
                    }
                    if(params.t[i][k][j] < maxW) {
                        maxW = params.t[i][k][j];
                    }
                }
                range = maxW - minW;
                alpha = rand_gen.randDblExc(); //NO
                width = range * alpha;
                for(int j=0; j<params.m; j++) {
                    if(params.t[i][k][j] <= minW + width) {
                        RCL.push_back(params.t[i][k][j]); // maybe costly
                        index_RCL.push_back(j);
                    }
                }
                sort(index_RCL.begin(),index_RCL.end(), [&](int i,int j){return RCL[i]<RCL[j];} );
                sort(RCL.begin(), RCL.end());
                proba.resize(RCL.size(),0);
                for(unsigned int r=0; r<proba.size(); r++) { //r for rank
                    proba[r]=1/r;
                    totalBias+=1/r;
                }
                for(unsigned int r=0; r<proba.size(); r++) {
                    proba[r]/=totalBias;
                }
                mach_runif = rand_gen.randDblExc();
                for(unsigned int r=0; r<proba.size(); r++) {
                    if(mach_runif <= proba[r]) {
                        x_temp[i][k][index_RCL[r]] = true;
                        W_tot_temp+=params.t[i][k][index_RCL[r]];
                        break;
                    } else {
                        mach_runif -= proba[r];
                    }
                }
            }
        }
        if(t==0) {
            x_sol = x_temp;
        } else if (W_tot_temp < W_tot) {
            delete_x(x_sol,params);
            x_sol = x_temp;
        } else {
            delete_x(x_temp,params);
        }
    }
    return x_sol;
}

int main()
{
    parameters params;
    read_parameters("problem8x8.in",params);
    bool*** x;
    delete_parameters(params);
    // decision varaibles -> c, y, x
    return 0;
}
