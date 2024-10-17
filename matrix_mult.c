#include <stdio.h>
#include <immintrin.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

#define N 2048 

void initialize_matrices(double *A, double *B, double *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = 0.001;
            B[i * N + j] = 0.002;
            C[i * N + j] = 0.0;
        }
    }
}

void matrix_multiply_serial(double *A, double *B, double *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}


void matrix_multiply_vectorized(double *A, double *B, double *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256d sum_vec = _mm256_setzero_pd(); 
            int k;
            for (k = 0; k + 3 < N; k += 4) {
                // Загрузка 4 значений из A
                __m256d a_vec = _mm256_loadu_pd(&A[i * N + k]);
                __m256d b_vec = _mm256_broadcast_sd(&B[k * N + j]);
                sum_vec = _mm256_add_pd(sum_vec, _mm256_mul_pd(a_vec, b_vec));
            }
            
            // double temp[4];
            __attribute__ ((aligned (32))) double temp[4];

            _mm256_store_pd(temp, sum_vec);
            C[i * N + j] += temp[0] + temp[1] + temp[2] + temp[3];

            // Обычная обработка оставшихся элементов (если есть)
            for (; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

// Пока не нужно 
// void matrix_multiply_vectorized(double *A, double *B, double *C) {
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             __m256d sum_vec = _mm256_setzero_pd();
//             for (int k = 0; k < N; k += 4) {
               
//                 __m256d a_vec = _mm256_load_pd(&A[i * N + k]);
               
//                 __m256d b_vec = _mm256_broadcast_sd(&B[k * N + j]);
              
//                 sum_vec = _mm256_add_pd(sum_vec, _mm256_mul_pd(a_vec, b_vec));
//             }
         
//             double temp[4];
//             _mm256_store_pd(temp, sum_vec);
//             C[i * N + j] += temp[0] + temp[1] + temp[2] + temp[3];
//         }
//     }
// }

int main() {
    // double *A = (double *)_mm_malloc(N * N * sizeof(double), 32);
    // double *B = (double *)_mm_malloc(N * N * sizeof(double), 32);
    // double *C_serial = (double *)_mm_malloc(N * N * sizeof(double), 32);
    // double *C_vectorized = (double *)_mm_malloc(N * N * sizeof(double), 32);
    double *A = aligned_alloc(32, N * N * sizeof(double));
    double *B = aligned_alloc(32, N * N * sizeof(double));
    double *C_serial = aligned_alloc(32, N * N * sizeof(double));
    double *C_vectorized = aligned_alloc(32, N * N * sizeof(double));
   

    initialize_matrices(A, B, C_serial);
    
    double start_time, end_time;
   
    start_time = omp_get_wtime();
    matrix_multiply_serial(A, B, C_serial);
    end_time = omp_get_wtime();
    printf("Serial Time: %lf secondsn \n ", end_time - start_time);


    initialize_matrices(A, B, C_vectorized); 
    start_time = omp_get_wtime();
    matrix_multiply_vectorized(A, B, C_vectorized);
    end_time = omp_get_wtime();
    printf("Vectorized Time: %lf secondsn \n", end_time - start_time);

   
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         if (fabs(C_serial[i * N + j] - C_vectorized[i * N + j]) > 1e-6) {
    //             printf("Results do not match at (%d, %d)n", i, j);
    //             return -1;
    //         }
    //     }
    // }

    // printf("Results match!n");

    // _mm_free(A);
    // _mm_free(B);
    // _mm_free(C_serial);
    // _mm_free(C_vectorized);
    free(A);
    free(B);
    free(C_serial);
    free(C_vectorized);

    return 0;
}
