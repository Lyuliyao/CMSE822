#include<iostream>
using std::cout; using std::endl;
#include<math.h>
#include<omp.h>
#include <fstream>
using namespace std;

#define PI 3.14159265
#define N 10
#define left_b 0
#define right_b 0

float fcn_source(float x){
	return pow(PI,2)*sin(PI*x);
}
//f_l(x,n) = (x+1)/N/4*fcn_left((x+1)/2/N+(n-1)/N)
float f_l(float x, int n){
	return (x+1)/float(N)/4*fcn_source((x+1)/2/float(N)+(float(n)-1.0)/float(N));
}
float f_r(float x, int n){
	return (1-x)/float(N)/4*fcn_source((x+1)/2/float(N)+float(n)/float(N));
}
int main(){
	float mass_matrix[N+2][N+2];
	float left_vector[N+2];
	float u[N+2];
	for(int i = 1; i <= N; i++){
		mass_matrix[i][i] = 2*float(N);
		mass_matrix[i][i-1] = -float(N);
		mass_matrix[i][i+1] = -float(N);
		left_vector[i]= f_l(-1/sqrt(3),i-1) + f_l(1/sqrt(3),i-1);
		left_vector[i]= left_vector[i] + f_r(-1/sqrt(3),i-1) + f_r(1/sqrt(3),i-1);
	}

	u[0]= left_b;
	u[N] = right_b;
	for (int t=0; t<100; t++){
		for (int i=1; i<=N;i++){
			u[i] = (left_vector[i] - mass_matrix[i][i-1]*u[i-1] - mass_matrix[i][i+1] *u[i+1])/mass_matrix[i][i];
		}
	}
	ofstream myfile;
  myfile.open ("result.txt");
	for (int i = 0; i<N+1;i++){
		myfile << u[i] << " ";
	}
	myfile << u[N+1];
  myfile.close();
}
