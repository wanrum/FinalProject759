#include <iostream>
#include "omp.h"
#include "solve_LU.h"

using namespace std;

void solve_LU(const float *L, const float *U, const float *b,
              float *vec_1, float *vec_2, int n){
                float dummy,temp;
                //L*vec_1 =b
                #pragma omp parallel num_threads(4) shared(vec_1,vec_2)
                  {
                for(int i=0; i<n; i++){//row
                  //vec_1[0] = b[0] because L[0,0] =  1
                  temp =0.0;
                   #pragma omp simd reduction(+:temp)
                   for(int j=0; j<i; j++){//columns
                     temp += L[i*n+j]*vec_1[j];
                   }
                   vec_1[i] = b[i]-temp;
                }
                //U*vec_2 = vec_1
                for(int i=n-1; i>=0; i--){
                  dummy =0.0;
                  #pragma omp simd reduction(+:dummy)
                  for(int j=i+1;j<n; j++){
                    dummy += U[i*n+j]*vec_2[j];
                  }
                  vec_2[i] = (vec_1[i]-dummy)/U[i*n+i];
                }
                }
              }
