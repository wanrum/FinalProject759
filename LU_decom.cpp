#include <iostream>
#include "omp.h"
#include "LU_decom.h"
using namespace std;

void LU_decom(const float *A, float *L, float *U, std::size_t n){
  // A,L and U are stored in row directions.
  for(size_t i=0; i<n-1; i++){
    #pragma omp parallel for schedule (dynamic)
    for(size_t j=i+1; j<n; j++){
    L[j*n+i]= U[j*n+i]/U[i*n+i];
    for(size_t k=0; k<n; k++){
    U[j*n+k] = U[j*n+k]-L[j*n+i]*U[i*n+k];
  }
  }
}
}
