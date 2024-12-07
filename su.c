#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <sys/time.h> 

void sequential_multiply(float *A, float *B, float *C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < N; k++) {
                sum += A[i*N + k] * B[k*N + j];
            }
            C[i*N + j] = sum;
        }
    }
}

void print_matrix(const char *name, float *M, int N) {
    printf("%s:\n", name);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%5.1f ", M[i*N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s N b\n", argv[0]);
            fprintf(stderr, "N - size of the NxN matrix\n");
            fprintf(stderr, "b - block size for SUMMA steps\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int N = atoi(argv[1]);
    int b = atoi(argv[2]);

    int q = (int) sqrt(size);
    if (q * q != size) {
        if (rank == 0)
            fprintf(stderr, "Number of processes must be a perfect square.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (N % q != 0) {
        if (rank == 0)
            fprintf(stderr, "N must be divisible by sqrt(P).\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int local_n = N / q;
    if (local_n % b != 0) {
        if (rank == 0)
            fprintf(stderr, "local_n must be divisible by b.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int dims[2] = {q, q};
    int periods[2] = {0, 0};
    MPI_Comm grid_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

    int coords[2];
    MPI_Cart_coords(grid_comm, rank, 2, coords);
    int row_coord = coords[0];
    int col_coord = coords[1];

    MPI_Comm row_comm, col_comm;
    MPI_Comm_split(grid_comm, row_coord, rank, &row_comm);
    MPI_Comm_split(grid_comm, col_coord, rank, &col_comm);

    float *A_local = (float*)malloc(local_n * local_n * sizeof(float));
    float *B_local = (float*)malloc(local_n * local_n * sizeof(float));
    float *C_local = (float*)calloc(local_n * local_n, sizeof(float));

    float *A_global = NULL;
    float *B_global = NULL;

    srand(time(NULL) + rank);

    if (rank == 0) {
        A_global = (float*)malloc(N * N * sizeof(float));
        B_global = (float*)malloc(N * N * sizeof(float));

        for (int i = 0; i < N*N; i++) {
            A_global[i] = (float)(rand() % 10);
        }

        srand(time(NULL) + size); 
        for (int i = 0; i < N*N; i++) {
            B_global[i] = (float)(rand() % 10);
        }

        // printf("Matrices on root (rank=0):\n");
        // print_matrix("A", A_global, N);
        // print_matrix("B", B_global, N);
    }

    MPI_Datatype block_type, block_resized, global_block;
    MPI_Type_contiguous(local_n, MPI_FLOAT, &block_type);
    MPI_Type_commit(&block_type);

    MPI_Type_vector(local_n, local_n, N, MPI_FLOAT, &global_block);
    MPI_Type_create_resized(global_block, 0, sizeof(float), &block_resized);
    MPI_Type_commit(&block_resized);

    int *sendcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        sendcounts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                sendcounts[i*q + j] = 1;
                displs[i*q + j] = i*(N*local_n) + j*local_n;
            }
        }
    }

    MPI_Scatterv(A_global, sendcounts, displs, block_resized,
                 A_local, local_n*local_n, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(B_global, sendcounts, displs, block_resized,
                 B_local, local_n*local_n, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

    int nb = local_n / b;
    float *A_panel = (float*)malloc(local_n * b * sizeof(float));
    float *B_panel = (float*)malloc(b * local_n * sizeof(float));


    MPI_Barrier(MPI_COMM_WORLD);
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // SUMMA
    for (int k = 0; k < N; k += b) {
        int owner_col = (k / b) / nb; 
        int owner_row = (k / b) / nb;

        int local_kA = (k % N) % local_n;
        int local_kB = (k % N) % local_n;


        if (col_coord == owner_col) {
            for (int ii = 0; ii < local_n; ii++) {
                for (int jj = 0; jj < b; jj++) {
                    A_panel[ii*b + jj] = A_local[ii*local_n + (local_kA + jj)];
                }
            }
        }
        MPI_Bcast(A_panel, local_n*b, MPI_FLOAT, owner_col, row_comm);

    
        if (row_coord == owner_row) {
            for (int ii = 0; ii < b; ii++) {
                for (int jj = 0; jj < local_n; jj++) {
                    B_panel[ii*local_n + jj] = B_local[(local_kB + ii)*local_n + jj];
                }
            }
        }
        MPI_Bcast(B_panel, b*local_n, MPI_FLOAT, owner_row, col_comm);


        for (int ii = 0; ii < local_n; ii++) {
            for (int jj = 0; jj < local_n; jj++) {
                float sum = 0.0f;
                for (int kk = 0; kk < b; kk++) {
                    sum += A_panel[ii*b + kk] * B_panel[kk*local_n + jj];
                }
                C_local[ii*local_n + jj] += sum;
            }
        }
    }

  
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&end, NULL);

    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1.0e-6;

    free(A_panel);
    free(B_panel);

    float *C_global = NULL;
    int *recvcounts = NULL;
    int *displs_c = NULL;
    if (rank == 0) {
        C_global = (float*)malloc(N*N*sizeof(float));
        recvcounts = (int*)malloc(size * sizeof(int));
        displs_c = (int*)malloc(size * sizeof(int));
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                recvcounts[i*q + j] = 1;
                displs_c[i*q + j] = i*(N*local_n) + j*local_n;
            }
        }
    }

    MPI_Gatherv(C_local, local_n*local_n, MPI_FLOAT,
                C_global, recvcounts, displs_c, block_resized,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        float *C_seq = (float*)malloc(N*N*sizeof(float));
        sequential_multiply(A_global, B_global, C_seq, N);

        // print_matrix("C (sequential)", C_seq, N);
        // print_matrix("C (parallel)", C_global, N);

        int correct = 1;
        for (int i = 0; i < N*N; i++) {
            if (fabs(C_global[i] - C_seq[i]) > 1e-5) {
                correct = 0;
                break;
            }
        }

        if (correct) {
            printf("Result is correct.\n");
        } else {
            printf("Result is incorrect.\n");
        }

        printf("Time for SUMMA: %f seconds\n", elapsed);

        free(C_seq);
        free(C_global);
        free(A_global);
        free(B_global);
        free(recvcounts);
        free(displs_c);
    }

    free(A_local);
    free(B_local);
    free(C_local);

    MPI_Type_free(&block_type);
    MPI_Type_free(&block_resized);
    MPI_Type_free(&global_block);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&grid_comm);

    MPI_Finalize();
    return 0;
}
