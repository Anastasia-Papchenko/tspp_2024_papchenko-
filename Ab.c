#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define N 1000

void multiply(float *A, float *b, float *local_c, int local_rows) {
    for (int i = 0; i < local_rows; i++) {
        local_c[i] = 0.0f;
        for (int j = 0; j < N; j++) {
            local_c[i] += A[i * N + j] * b[j];
        }
    }
}

void print_matrix(float *A, int rows, int cols) {
    printf("Matrix A:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%5.1f ", A[i * cols + j]);
        }
        printf("\n");
    }
}

void print_vector(float *v, int size) {
    for (int i = 0; i < size; i++) {
        printf("%5.1f ", v[i]);
    }
    printf("\n");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);  

    int periods[2] = {0, 0};  
    MPI_Comm comm_2d;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);  // reorder=1
    int coords[2];
    MPI_Cart_coords(comm_2d, rank, 2, coords);

    int rows_size = N / size;
    int remains = N % size;
    int local_rows = rows_size + (rank < remains ? 1 : 0);

    int row_start = rank * rows_size + (rank < remains ? rank : remains);

    float *local_A = (float*)malloc(local_rows * N * sizeof(float));
    float *local_c = (float*)malloc(local_rows * sizeof(float));
    float *b = (float*)malloc(N * sizeof(float));

    for(int i = 0; i < local_rows; i++) {
        int global_row = row_start + i;
        for(int j = 0; j < N; j++) {
            local_A[i * N + j] = (float)(global_row * N + j); 
        }
    }
 
    float *A_serial = NULL; 
    float *c_final = NULL;
    if(rank == 0) {
        for(int i = 0; i < N; i++) {
            b[i] = (float)(i + 1);  
        }

        A_serial = (float*)malloc(N * N * sizeof(float));
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                A_serial[i * N + j] = (float)(i * N + j);
            }
        }
        c_final = (float*)calloc(N, sizeof(float));

        // print_matrix(A_serial, N, N);
        // printf("Vector b:\n");
        // print_vector(b, N);
        // printf("\n");
    }

    MPI_Win win_b, win_c;
    if(rank == 0) {
        MPI_Win_create(b, N * sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_b);
        MPI_Win_create(c_final, N * sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c);
    }
    else {
        MPI_Win_create(NULL, 0, sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_b);
        MPI_Win_create(NULL, 0, sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c);
    }

    MPI_Win_fence(0, win_b);
    MPI_Win_fence(0, win_c);

    if(rank != 0) {
        MPI_Get(b, N, MPI_FLOAT, 0, 0, N, MPI_FLOAT, win_b);
    }

    MPI_Win_fence(0, win_b);
    MPI_Win_fence(0, win_c);

    double start_time = MPI_Wtime();

    multiply(local_A, b, local_c, local_rows);

    MPI_Win_fence(0, win_c);  

    MPI_Put(local_c, local_rows, MPI_FLOAT, 0, row_start, local_rows, MPI_FLOAT, win_c);

    MPI_Win_fence(0, win_c);  

    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed;

    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        // printf("Result Vector c:\n");
        // print_vector(c_final, N);
        // printf("\n");
        printf("Execution Time: %f seconds\n", max_elapsed);
        printf("\n");
    }

    MPI_Win_free(&win_b);
    MPI_Win_free(&win_c);
    free(local_A);
    free(local_c);
    free(b);
    if(rank == 0) {
        free(c_final);
        free(A_serial);
    }

    MPI_Comm_free(&comm_2d);
    MPI_Finalize();
    return 0;
}
