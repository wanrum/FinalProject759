#include <iostream>
#include "omp.h"
#include "Jacobi.h"
using namespace std;
void create_matrix(float *A, std::size_t n){
  fill(A,A+n*n,0.0);
  for(size_t i=0; i<n; i++){
    A[i*n+i] = 4;
    if(i>0){
    A[i*n+i-1] = -1;
    }
    A[i*n+i+1] = -1;
    if(i>2){
      A[i*n+i-3] = -1;
    }
    if(i+3<n){
      A[i*n+i+3] = -1;
    }
  }
}
int main(int argc, char* argv[]){
  size_t n = atoi(argv[1]);
  // float arr[] = {4,-1,0,-2,0,0,-1,4,-1,0,-2,0,0,-1,4,0,0,-2,-1,0,-1,4,-1,0,0,-1,0,-1,4,-1,0,0,-1,0,-1,4};
  // float array[] = {-1,0,1,-2,1,2};
  //float arr[] = {4,-1,0,-2,0,0,0,0,-1,4,-1,0,-2,0,0,0,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,0,0,0,-1,0,-1,4,-1,0,0,0,0,-1,0,-1,4};
  float array[] = {-1,0,1,-2,1,2,-1,0};
  float *A = new float[n*n];
  create_matrix(A, n);
  float *L_U = new float[n*n];
  float *x0 = new float[n];
  //float *T= new float[n*n];
  float *b = new float[n];
  float *x = new float[n];
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      //A[i*n+j] = arr[i*n+j];
      L_U[i*n+j] = -A[i*n+j];
      //[i*n+j] =0.0;
    }
    b[i] =array[i];
    ////x[i] = 0.0; ///start iteration from x[i] = 1;
  }
  double start = omp_get_wtime();
  jacobi(A, b, L_U, x0, x,n);
  double end = omp_get_wtime();
  for(size_t i=0; i<n; i++){
    cout<<"x:"<<x[i]<<endl;
  }
  cout << " time : " << end - start << endl;
  // for(size_t i=0; i<n; i++){
  //   for(size_t j=0; j<n; j++){
  //     cout<<T[i*n+j] <<", ";
  //   }
  // }
  cout<<endl;
  delete[] A;
  delete[] L_U;
  delete[] x0;
  ///delete[] T;
  delete[] b;
  delete[] x;
  return 0;
}
