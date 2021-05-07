#include <iostream>
#include "omp.h"
#include <chrono>
#include <random>
#include "LU_decom.h"
#include "solve_LU.h"
using namespace std;
// Solving linear equations by using LU_decomposition
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
int main(int argc, char *argv[]){
    size_t n = atoi(argv[1]);
    // float arr[] = {4,-1,0,-2,0,0,-1,4,-1,0,-2,0,0,-1,4,0,0,-2,-1,0,-1,4,-1,0,0,-1,0,-1,4,-1,0,0,-1,0,-1,4};
    // float array[] = {-1,0,1,-2,1,2};
    //float arr[] = {4,-1,0,-2,0,0,0,0,-1,4,-1,0,-2,0,0,0,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,-2,0,0,-1,0,-1,4,-1,0,0,0,0,-1,0,-1,4,-1,0,0,0,0,-1,0,-1,4};
    float array[] = {-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2,1,-1,-2,-1,0,1,2,4,3,2,1,-1,1,-1,-2,-1,0,1,2,4,3,2,1,-1,1,-1,-2,-1,0,1,2,4,3,2,1,-1,1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,1,2,4,3,2,1,-1,-2,-1,0,1,-2,1,2,-1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1,-2-1,0,1,-2,1,2,-1,0,1,2,4,3,2,1,-1,-2-1,0,1};
    //const int n_threads = 4;
    // chrono::high_resolution_clock::time_point start;
    // chrono::high_resolution_clock::time_point end;
    // chrono::duration<double,std::micro> duration_sec;
    float *A = new float[n*n];
    create_matrix(A, n);
    float *L = new float[n*n];
    float *U = new float[n*n];
    float *vector1 = new float[n];
    float *vector2 = new float[n];
    float *b = new float[n];
    for(size_t i=0; i<n; i++){
      for(size_t j=0; j<n; j++){
        //A[i*n+j] = arr[i*n+j];
        U[i*n+j] = A[i*n+j];
        L[i*n+j] = 0.0;
        if(i==j){
          L[i*n+j] = 1;
        }
      }
      vector1[i] =0.0;
      vector2[i] =0.0;
      b[i] =array[i];
    }
    int num = (int)n;
    double start = omp_get_wtime();
    LU_decom(A,L,U,n);
    // for(size_t i=0; i<n; i++){
    //   for(size_t j=0; j<n; j++){
    //     cout<<"L"<<L[i*n+j]<<endl;
    //     cout<<"U"<<U[i*n+j]<<endl;
    //   }
    // }
    // start = chrono::high_resolution_clock::now();
    solve_LU(L,U,b,vector1, vector2,num);
    double end = omp_get_wtime();
    // end = chrono::high_resolution_clock::now();
    //duration_sec = chrono::duration_cast<std::chrono::duration<double, std::micro>>(end - start);
    for(size_t i=0; i<n; i++){
      //cout<<"vec1:"<<vector1[i]<<endl;
      cout<<"vec2:"<<vector2[i]<<endl;
    }
    //cout <<"Time:   "<< duration_sec.count()<< endl;
    cout << " LU decomposition without openmp time : " << end - start << endl;
    delete[] A;
    delete[] L;
    delete[] U;
    delete[] vector1;
    delete[] vector2;
    delete[] b;
    return 0;
}
