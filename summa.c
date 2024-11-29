#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


void initialize_matrix(float *matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = rand() % 10; 
    }
}


void sequential_matrix_multiplication(float *A, float *B, float *C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i * N + j] = 0;
            for (int k = 0; k < N; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}


void print_matrix(float *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.1f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int rank, size;
    int N, b; 
    MPI_Comm comm_2d;

    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) {
            printf("Usage: %s <matrix_size> <block_size>\n", argv[0]);
        }
        MPI_Finalize();
        return -1;
    }

    N = atoi(argv[1]); 
    b = atoi(argv[2]); 

    int sqrt_p = (int)sqrt(size); 
    if (sqrt_p * sqrt_p != size || N % sqrt_p != 0 || (N / sqrt_p) % b != 0) {
        if (rank == 0) {
            printf("Error: incorrect parameters. Make sure that P is a complete square, N is divisible by sqrt(P), and b is a divisor of N/sqrt(P).\n");
        }
        MPI_Finalize();
        return -1;
    }


    int dims[2] = {sqrt_p, sqrt_p};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);

    int coords[2];
    MPI_Cart_coords(comm_2d, rank, 2, coords);

    int local_size = N / sqrt_p; 
    float *A_local = malloc(local_size * local_size * sizeof(float));
    float *B_local = malloc(local_size * local_size * sizeof(float));
    float *C_local = calloc(local_size * local_size, sizeof(float));

    float *A = NULL, *B = NULL, *C_seq = NULL;
    if (rank == 0) {
        A = malloc(N * N * sizeof(float));
        B = malloc(N * N * sizeof(float));
        C_seq = malloc(N * N * sizeof(float));
        initialize_matrix(A, N, N);
        initialize_matrix(B, N, N);
    }

    
    int *recvcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));

        for (int i = 0; i < sqrt_p; i++) {
            for (int j = 0; j < sqrt_p; j++) {
                recvcounts[i * sqrt_p + j] = local_size * local_size;
                displs[i * sqrt_p + j] = (i * N * local_size) + (j * local_size);
            }
        }
    }

    MPI_Scatterv(A, recvcounts, displs, MPI_FLOAT, A_local, local_size * local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, recvcounts, displs, MPI_FLOAT, B_local, local_size * local_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    
    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int k = 0; k < sqrt_p; k++) {
        float *A_k = malloc(local_size * local_size * sizeof(float));
        float *B_k = malloc(local_size * local_size * sizeof(float));

        
        MPI_Bcast(A_k, local_size * local_size, MPI_FLOAT, k, comm_2d);
        MPI_Bcast(B_k, local_size * local_size, MPI_FLOAT, k, comm_2d);

        
        for (int i = 0; i < local_size; i++) {
            for (int j = 0; j < local_size; j++) {
                for (int l = 0; l < local_size; l++) {
                    C_local[i * local_size + j] += A_k[i * local_size + l] * B_k[l * local_size + j];
                }
            }
        }

        free(A_k);
        free(B_k);
    }

    gettimeofday(&end, NULL);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;

    if (rank == 0) {
        printf("Execution time: %f seconds\n", elapsed);
    }

  
    MPI_Gatherv(C_local, local_size * local_size, MPI_FLOAT, C_seq, recvcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Матрица A:\n");
        print_matrix(A, N, N);

        printf("Матрица B:\n");
        print_matrix(B, N, N);


        printf("Матрица C (параллельная):\n");
        print_matrix(C_seq, N, N);

        float *C_check = malloc(N * N * sizeof(float));
        sequential_matrix_multiplication(A, B, C_check, N);
        printf("Матрица C (последовательная):\n");
        print_matrix(C_check, N, N);
        free(C_check);

        free(recvcounts);
        free(displs);
    }

    free(A_local);
    free(B_local);
    free(C_local);
    if (rank == 0) {
        free(A);
        free(B);
        free(C_seq);
    }

    MPI_Finalize();
    return 0;
}