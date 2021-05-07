#include <iostream>
#include "omp.h"
#include <cmath>
#include "Jacobi.h"

using namespace std;

void jacobi(const float *A, const float *b, float *L_U, float *x0, float *x,std::size_t n){
  float sum,diff,dummy = 1.0,temp;
  while (dummy>0.000001){
    #pragma omp parallel num_threads(4)
    {
    for(size_t i=0; i<n; i++){
      temp = 1/A[i*n+i]*b[i];
      L_U[i*n+i] = 0.0;
      sum=0.0;
      #pragma omp simd reduction(+:sum)
      for(size_t j=0; j<n; j++){
         sum += 1/A[i*n+i]*L_U[i*n+j]*x[j];
      }
      x0[i] = sum+temp;
    }

      diff = 0.0;
      #pragma omp simd reduction(+:diff)
      for(size_t i=0; i<n; i++){
        diff += (x0[i]-x[i])*(x0[i]-x[i]);
        x[i] = x0[i];
      }
     dummy = sqrt(diff);
  // cout<<"difference:"<<dummy<<endl;
}
}
}
