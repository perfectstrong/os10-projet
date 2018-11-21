#include <vector>
#include <fstream>

using namespace std;

int read_problem()
{
    return 0;
}

int main()
{
    int n=0; //number of jobs
    int m=0; //number of machines
    int* n_o = new int[n]; //number of operations of jobs i 
    int* L = new int[m]; //number of maintenance tasks on machine j
    bool*** A = new bool** [n]; //availables machines for task k of job i
    double*** t = new double**[n]; //processing time for task k of job i on machine j
    for(int i=0; i<n; i++) {
        A[i] = new bool*[n_o[i]];
        t[i] = new double*[n_o[i]];
        for(int k=0; k<n_o[i]; i++) {
            A[i][k] = new bool[m];
            t[i][k] = new double[m];
        }
    }
    double** p = new double*[m]; //processing time for maintenance task l on machine j
    double** t_e = new double*[m]; //time windows begin (earliest)
    double** t_l = new double*[m]; //time windows end (lateness)
    for(int j=0; j<m; j++) {
        p[j] = new double[L[m]];
        t_e[j] = new double[L[m]];
        t_l[j] = new double[L[m]];
    }
    
    // decision varaibles -> c, y, x
    return 0;
}
