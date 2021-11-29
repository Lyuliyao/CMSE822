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

int main()
{
	int m;
	int n;
	ifstream f1("A.txt");
	ifstream f2("b.txt");
	f1 >> m >> n;
	float mass_matrix[n][m];
	float left_vector[n];
    float u[m];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++){
	  f1 >> mass_matrix[i][j];
	}
	for (int i = 0; i < m; i++){
	  f2 >> left_vector[i] ;
	}
	for (int t=0; t<100; t++){
		for (int i=0; i<m;i++){
		    u[i] = left_vector[i];
		    for (int j=0; j <m;j++){
		        if(i!=j)
			    u[i] = u[i] - mass_matrix[i][j]*u[j];
		    }
		    u[i] = u[i]/mass_matrix[i][i];
		    cout << u[i] << endl;
		}
	}
    ofstream myfile;
    myfile.open ("result.txt");
	for (int i = 0; i<=m;i++){
		myfile << u[i] << " ";
	}
	myfile << u[m+1];
  myfile.close();
}
