#include<iostream>
using std::cout; using std::endl;
#include<math.h>
#include<omp.h>
#include <fstream>
using namespace std;

#define PI 3.14159265
#define N 100000
#define left_b 0
#define right_b 0

int main(int argc, char **argv){
	int n_tct = atoi(argv[1]);
	int m;
	int n;
	ifstream f1("A.txt");
	ifstream f2("b.txt");
	f1 >> m >> n;
	float mass_matrix[n][m];
	float left_vector[n];
    float u[m];
    float u_new[m];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++){
	  f1 >> mass_matrix[i][j];
	}
	for (int i = 0; i < m; i++){
	  f2 >> left_vector[i] ;
	}
	omp_set_nested(1);
	omp_set_num_threads(n_tct);
	double t1 = omp_get_wtime();
	for (int t=0; t<10000; t++){
	    #pragma omp parallel for shared(u_new,u)
		for (int i=0; i<m;i++){
		    u_new[i] = left_vector[i];
		    for (int j=0; j <m;j++){
		        if(i!=j)
			    u_new[i] = u_new[i] - mass_matrix[i][j]*u[j];
		    }
		    u_new[i] = u_new[i]/mass_matrix[i][i];
		}
	    #pragma omp parallel for shared(u,u_new)
		for (int i=0; i<m;i++){
		    u[i] =u_new[i];
		}
	}
	double t2 = omp_get_wtime();
    ofstream myfile;
    myfile.open ("result.txt");
	for (int i = 0; i<m;i++){
		myfile << u[i] << " ";
	}
  myfile.close();
    cout << (t2 - t1)/100000 << endl;
}
