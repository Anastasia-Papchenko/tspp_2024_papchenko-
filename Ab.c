#include <stdio.h> 
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

void fill_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = rand() % 10;
    }
}

void fill_vector(double *vector, int size) {
    for (int i = 0; i < size; i++) {
        vector[i] = rand() % 10;
    }
}

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void print_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
}

void print_vector(double *vector, int size) {
    for (int i = 0; i < size; i++) {
        printf("%.2f ", vector[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[2];
    int periods[2] = {0, 0};
    MPI_Comm grid_comm;

    MPI_Dims_create(size, 2, dims);

    int create_result = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &grid_comm);
    if (create_result != MPI_SUCCESS) {
        if (rank == 0) {
            printf("Error: Failed to create Cartesian topology.\n");
        }
        MPI_Finalize();
        return -1;
    }

    int coords[2];
    MPI_Cart_coords(grid_comm, rank, 2, coords);

    int global_rows = 8; 
    int global_cols = 8;

    if (global_rows % dims[0] != 0 || global_cols % dims[1] != 0) {
        if (rank == 0) {
            printf("Error: Matrix dimensions must be divisible by process grid dimensions.\n");
        }
        MPI_Finalize();
        return -1;
    }

    int local_rows = global_rows / dims[0];
    int local_cols = global_cols / dims[1];

    double *local_matrix = (double *)malloc(local_rows * local_cols * sizeof(double));
    if (local_matrix == NULL) {
        printf("Error: Unable to allocate memory for local_matrix.\n");
        MPI_Finalize();
        return -1;
    }

    double *vector = NULL;
    double *local_result = (double *)calloc(local_rows, sizeof(double));
    if (local_result == NULL) {
        printf("Error: Unable to allocate memory for local_result.\n");
        free(local_matrix);
        MPI_Finalize();
        return -1;
    }
    
    double *result = NULL;

    if (rank == 0) {
        vector = (double *)malloc(global_cols * sizeof(double));
        if (vector == NULL) {
            printf("Error: Unable to allocate memory for vector.\n");
            free(local_matrix);
            free(local_result);
            MPI_Finalize();
            return -1;
        }
        fill_vector(vector, global_cols);

        printf("Vector b:\n");
        print_vector(vector, global_cols);

        result = (double *)malloc(global_rows * sizeof(double));
        if (result == NULL) {
            printf("Error: Unable to allocate memory for result.\n");
            free(local_matrix);
            free(local_result);
            free(vector);
            MPI_Finalize();
            return -1;
        }
    }

    fill_matrix(local_matrix, local_rows, local_cols);

    double *full_matrix = NULL;
    if (rank == 0) {
        full_matrix = (double *)malloc(global_rows * global_cols * sizeof(double));
    }

    MPI_Gather(local_matrix, local_rows * local_cols, MPI_DOUBLE, full_matrix, local_rows * local_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Matrix A:\n");
        print_matrix(full_matrix, global_rows, global_cols);
        free(full_matrix);
    }

    MPI_Win win;
    double *shared_vector;
    MPI_Win_allocate(global_cols * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &shared_vector, &win);

    if (rank == 0) {
        for (int i = 0; i < global_cols; i++) {
            shared_vector[i] = vector[i];
        }
    }

    MPI_Win_fence(0, win);

    double *local_vector = (double *)malloc(local_cols * sizeof(double));
    if (local_vector == NULL) {
        printf("Error: Unable to allocate memory for local_vector.\n");
        MPI_Win_free(&win);
        free(local_matrix);
        free(local_result);
        MPI_Finalize();
        return -1;
    }

    MPI_Get(local_vector, local_cols, MPI_DOUBLE, 0, 0, local_cols, MPI_DOUBLE, win);

    MPI_Win_fence(0, win);

    double start_time = get_time();

    for (int i = 0; i < local_rows; i++) {
        for (int j = 0; j < local_cols; j++) {
            local_result[i] += local_matrix[i * local_cols + j] * local_vector[j];
        }
    }

    double end_time = get_time();

    MPI_Gather(local_result, local_rows, MPI_DOUBLE, result, local_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Result vector:\n");
        for (int i = 0; i < global_rows; i++) {
            printf("%f\n", result[i]);
        }

        printf("Execution time: %.6f seconds\n", end_time - start_time);

        free(result);
        free(vector);
    }

    MPI_Win_free(&win);
    free(local_matrix);
    free(local_vector);
    free(local_result);

    MPI_Finalize();
    return 0;
}
