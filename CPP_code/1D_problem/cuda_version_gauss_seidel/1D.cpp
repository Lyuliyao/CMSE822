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

__global__ void gauss_seidel(float* dev_mass_matrix, float* dev_left_vector, float* dev_u, float* dev_u1, int width ,left_b,right_b) {
  //calculate the row and column for this element of the matrix
  int row = threadIdx.x + (blockDim.x * blockIdx.x);
  if ((row == 1) ){
      dev_u1[row] = left_b;
  }
  else if ((row == width +1) ){
      dev_u1[row] = right_b;
  }
  else if ((row <= width)) {
    dev_u1[row] = (dev_left_vector[i] - dev_mass_matrix[i][i-1]*u[i-1] - dev_mass_matrix[i][i+1] *dev_u[i+1])/dev_mass_matrix[i][i];
  }
}


int main(int argc, char **argv){
    int N = atoi(argv[1]);
    int thread_cnt = atoi(argv[2]);
	float mass_matrix[N+2][N+2];
	float left_vector[N+2];
	float dev_mass_matrix[N+2][N+2];
	float dev_left_vecotr[N+2][N+2];
	float dev_u[N+2];
	float dev_u1[N+2];
	float u[N+2];
    int block_cnt;
	for(int i = 1; i <= N; i++){
		mass_matrix[i][i] = 2*float(N);
		mass_matrix[i][i-1] = -float(N);
		mass_matrix[i][i+1] = -float(N);
		left_vector[i]= f_l(-1/sqrt(3),i-1) + f_l(1/sqrt(3),i-1);
		left_vector[i]= left_vector[i] + f_r(-1/sqrt(3),i-1) + f_r(1/sqrt(3),i-1);
	}

	u[0]= left_b;
	u[N] = right_b;
	int bytes1 = (N+2)*(N+2) * sizeof(float);
    int bytes2 = (N+2) * sizeof(float);
    cudaMalloc((void **) &dev_mass_matrix, bytes1);
      
    cudaMalloc((void **) &dev_left_vecotr, bytes2);
    cudaMalloc((void **) &dev_u, bytes2);
    cudaMalloc((void **) &dev_u1, bytes2);
    cudaMemcpy(dev_mass_matrix, mass_matrix, bytes1, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_left_vecotr, left_vector, bytes2, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_u, u, bytes2, cudaMemcpyHostToDevice);
    block_cnt = n/thread_cnt + (n % thread_cnt > 0);
	for (int t=0; t<10000; t++){
	    gauss_seidel<<<block_cnt, thread_cnt>>>(dev_mass_matrix, dev_left_vecotr, dev_u, dev_u1, N,left_b,right_b);
	    dev_u = dev_u1;
	}
	ofstream myfile;
  myfile.open ("result.txt");
	for (int i = 0; i<N+1;i++){
		myfile << u[i] << " ";
	}
	myfile << u[N+1];
  myfile.close();
}
