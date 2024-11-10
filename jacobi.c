#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define N_ITER 100  

void initialize(int *subgrid, int rows, int N) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < N; j++) {
            subgrid[i * N + j] = rand() % 100; 
        }
    }
}

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

double calculate_norm(int *subgrid, int *new_subgrid, int rows, int N) {
    double norm = 0.0;
    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            double diff = new_subgrid[i * N + j] - subgrid[i * N + j];
            norm += diff * diff;
        }
    }
    return sqrt(norm);
}

int is_power_of_two(int n) {
    return (n && !(n & (n - 1)));
}

int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    
    return MPI_Send(buf, count, datatype, dest, tag, comm);  
}

int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
   
    return MPI_Recv(buf, count, datatype, source, tag, comm, status);  
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N;
    
    if (argc < 2) {
        if (rank == 0) {
            printf("Error: Please provide the grid size N as a command line argument.\n");
            printf("N must be a power of two (e.g., 2, 4, 8, 16, ..., up to the limits of your system resources).\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    N = atoi(argv[1]); 

    
    if (!is_power_of_two(N)) {
        if (rank == 0) {
            printf("Error: The grid size N (%d) must be a power of two.\n", N);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    
    if (N % size != 0) {
        if (rank == 0) {
            printf("Error: The grid size N must be a multiple of the number of processes.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    int rows_per_process = N / size;
    int *subgrid = (int *)malloc(rows_per_process * N * sizeof(int));
    int *new_subgrid = (int *)malloc(rows_per_process * N * sizeof(int));

  
    initialize(subgrid, rows_per_process, N);

    double start_time = get_time();

    for (int iter = 0; iter < N_ITER; iter++) {
        if (rank > 0) {
            PMPI_Send(subgrid, N, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            PMPI_Recv(&subgrid[0], N, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            PMPI_Send(&subgrid[(rows_per_process - 1) * N], N, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            PMPI_Recv(&subgrid[rows_per_process * N], N, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

       
        for (int i = 1; i < rows_per_process - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                new_subgrid[i * N + j] = 0.25 * (subgrid[(i - 1) * N + j] + 
                                                   subgrid[(i + 1) * N + j] + 
                                                   subgrid[i * N + (j - 1)] + 
                                                   subgrid[i * N + (j + 1)]);
            }
        }
        
        for (int i = 1; i < rows_per_process - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                subgrid[i * N + j] = new_subgrid[i * N + j];
            }
        }

        if (iter == N_ITER - 1) {
            double norm = calculate_norm(subgrid, new_subgrid, rows_per_process, N);
            printf("Process %d, Norm of difference: %f\n", rank, norm);
        }
    }

    double end_time = get_time();

    if (rank == 0) {
        printf("Total execution time: %f seconds\n", end_time - start_time);
    }

    free(subgrid);
    free(new_subgrid);
    MPI_Finalize();
    return 0;
}
