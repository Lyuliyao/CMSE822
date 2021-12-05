#include<iostream>
using std::cout; using std::endl;
#include<math.h>
#include <mpi.h>
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
    int N = atoi(argv[1]);
	float mass_matrix[N+2][N+2];
	float left_vector[N+2];
	float u[N+2];
    int myrank,from,to,i,P;
    int tag = 666;		/* any value will do */
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	/* who am i */
    MPI_Comm_size(MPI_COMM_WORLD, &P); /* number of processors */
    double starttime = MPI_Wtime();
    float work = floor((float)N/(float)P);
    from = myrank * work +1;
    to = (myrank+1) * work ;
    if (myrank == P-1)
    {
    to = N;
    }
	for(int i = from; i <= to; i++){
		mass_matrix[i][i] = 2*float(N);
		mass_matrix[i][i-1] = -float(N);
		mass_matrix[i][i+1] = -float(N);
		left_vector[i]= f_l(-1/sqrt(3),i-1,N) + f_l(1/sqrt(3),i-1,N);
		left_vector[i]= left_vector[i] + f_r(-1/sqrt(3),i-1,N) + f_r(1/sqrt(3),i-1,N);
	}

    for (int i=1; i<=N;i++){
			u[i] = 0;
			}
	u[0]= left_b;
	u[N+1] = right_b;
	float global_u[N+2];
	MPI_Allreduce(&u,&global_u,(N+2),MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);


	for (int t=0; t<100000; t++){
		for (int i=from; i<=to;i++){
			u[i] = (left_vector[i] - mass_matrix[i][i-1]*global_u[i-1] - mass_matrix[i][i+1] *global_u[i+1])/mass_matrix[i][i];
		}
		MPI_Allreduce(&u,&global_u,(N+2),MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	}
    double endtime = MPI_Wtime();
    double full_time  = endtime-starttime;
	if (myrank==0)
	{
	ofstream myfile;
    myfile.open ("result.txt");
	for (int i = 0; i<N+1;i++){
		myfile << u[i] << " ";
	}
	myfile << u[N+1];
    myfile.close();
	cout  << full_time/100000 <<endl ;
	}

    // Output time taken by this rank
    //std::

  //cout << (t2 - t1)/100000 << endl;
}
