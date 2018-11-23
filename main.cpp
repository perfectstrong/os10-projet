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
    while(!file.eof()) {
        string id_word;
        file >> id_word;
        if(id_word.find('n') != std::string::npos) {
            file >> params.n;
            cout << "Got n : " << params.n << endl;
        } else if(id_word == "m") {
            file >> params.m;
            cout << "Got m : " << params.n << endl;
        } else if(id_word == "n_o") {
            int* n_o = new int[params.n];
            for(int i = 0; i<params.n; i++) {
                file >> n_o[i];
            }
        }
        else
        {
            cout << id_word << endl;
            if(id_word != "n") {
                cout << "NO" << endl;
            }
        }
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
