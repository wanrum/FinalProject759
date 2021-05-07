#ifndef LU_DECOM_H
#define LU_DECOM_H
#include <cstddef>
#include <omp.h>
///Thomas method
///matrix A
///matrix L
///matrix U


void LU_decom(const float *A, float *L, float *U, std::size_t n);

#endif
