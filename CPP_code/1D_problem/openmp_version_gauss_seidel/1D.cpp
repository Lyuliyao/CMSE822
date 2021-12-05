#include<iostream>
using std::cout; using std::endl;
#include<math.h>
#include<omp.h>
#include <fstream>
using namespace std;

#define PI 3.14159265
#define left_b 0
#define right_b 0

float fcn_source(float x){
	return pow(PI,2)*sin(PI*x);
}
//f_l(x,n) = (x+1)/N/4*fcn_left((x+1)/2/N+(n-1)/N)
float f_l(float x, int n, int N){
	return (x+1)/float(N)/4*fcn_source((x+1)/2/float(N)+(float(n)-1.0)/float(N));
}
float f_r(float x, int n, int N){
	return (1-x)/float(N)/4*fcn_source((x+1)/2/float(N)+float(n)/float(N));
}
int main(int argc, char **argv){
    int N = atoi(argv[2]);
	float mass_matrix[N+2][N+2];
	float left_vector[N+2];
	float u[N+2];
    int i;
	int n_tct = atoi(argv[1]);
	omp_set_nested(1);
	omp_set_num_threads(n_tct);
    #pragma omp parallel for shared(mass_matrix,left_vector)
	for(int i = 1; i <= N; i++){
		mass_matrix[i][i] = 2*float(N);
		mass_matrix[i][i-1] = -float(N);
		mass_matrix[i][i+1] = -float(N);
		left_vector[i]= f_l(-1/sqrt(3),i-1,N) + f_l(1/sqrt(3),i-1,N);
		left_vector[i]= left_vector[i] + f_r(-1/sqrt(3),i-1,N) + f_r(1/sqrt(3),i-1,N);
	}
	u[0]= left_b;
	u[N] = right_b;


	double t1 = omp_get_wtime();
	
	for (int t=0; t<100000; t++){
	    #pragma omp parallel for shared(u)
		for (int i=1; i<=N;i++){
			u[i] = (left_vector[i] - mass_matrix[i][i-1]*u[i-1] - mass_matrix[i][i+1] *u[i+1])/mass_matrix[i][i];
		}
	}
	double t2 = omp_get_wtime();
	ofstream myfile;
    myfile.open ("result.txt");
	for (int i = 0; i<N+1;i++){
		myfile << u[i] << " ";
	}
  myfile << u[N+1];
  myfile.close();
  cout << (t2 - t1)/100000 << endl;
}
