#ifndef JACOBI_H
#define JACOBI_H

//iteration method for solving Ax=b
void jacobi(const float *A, const float *b, float *L_U, float *x0, float *x,std::size_t n); 

#endif
