#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

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
    for(int i = 0; i<params.m; i++) {
        file >> params.L[i];
        cout << "Got L["<<i<<"] : " << params.L[i] << endl;
    }
    params.A = new bool**[params.n];
    for(int i=0; i<params.n; i++) {
        params.A[i] = new bool*[params.n_o[i]];
        for(int k=0; k<params.n_o[i]; k++) {
            params.A[i][k] = new bool[params.m];
            for(int j=0; j<params.m; j++)
            {
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

int main()
{
    //~ int n=0; //number of jobs
    //~ int m=0; //number of machines
    //~ int* n_o = new int[n]; //number of operations of jobs i 
    //~ int* L = new int[m]; //number of maintenance tasks on machine j
    //~ bool*** A = new bool** [n]; //availables machines for task k of job i
    //~ double*** t = new double**[n]; //processing time for task k of job i on machine j
    //~ for(int i=0; i<n; i++) {
        //~ A[i] = new bool*[n_o[i]];
        //~ t[i] = new double*[n_o[i]];
        //~ for(int k=0; k<n_o[i]; i++) {
            //~ A[i][k] = new bool[m];
            //~ t[i][k] = new double[m];
        //~ }
    //~ }
    //~ double** p = new double*[m]; //processing time for maintenance task l on machine j
    //~ double** t_e = new double*[m]; //time windows begin (earliest)
    //~ double** t_l = new double*[m]; //time windows end (lateness)
    //~ for(int j=0; j<m; j++) {
        //~ p[j] = new double[L[m]];
        //~ t_e[j] = new double[L[m]];
        //~ t_l[j] = new double[L[m]];
    //~ }
    
    parameters params;
    read_parameters("data_example.in",params);
    
    // decision varaibles -> c, y, x
    return 0;
}
