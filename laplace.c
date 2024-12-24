#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define N_ITER 100

void initialize(double *subgrid, int nx, int ny, int nz) {
    for (int i = 0; i < nx * ny * nz; i++) {
        subgrid[i] = (double)(rand() % 100) / 100.0;
    }
}

double calculate_norm(double *grid, double *new_grid, int nx, int ny, int nz) {
    double norm = 0.0;
    for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++)
            for (int k = 1; k < nz - 1; k++) {
                int idx = i * ny * nz + j * nz + k;
                double diff = new_grid[idx] - grid[idx];
                norm += diff * diff;
            }
    return sqrt(norm);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct timeval start, end;

    int N;
    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: %s <N>\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    N = atoi(argv[1]);
    if ((N & (N - 1)) != 0) {
        if (rank == 0)
            printf("Error: N must be a power of two.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);

    int periods[3] = {0, 0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);

    int coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    int nx = N / dims[0];
    int ny = N / dims[1];
    int nz = N / dims[2];

    double *subgrid = (double *)malloc((nx + 2) * (ny + 2) * (nz + 2) * sizeof(double));
    double *new_subgrid = (double *)malloc((nx + 2) * (ny + 2) * (nz + 2) * sizeof(double));
    initialize(subgrid, nx + 2, ny + 2, nz + 2);

    gettimeofday(&start, NULL);

    MPI_Datatype yz_face_type, xz_face_type, xy_face_type;
    

    MPI_Type_vector(ny, nz, nz + 2, MPI_DOUBLE, &yz_face_type);
    MPI_Type_commit(&yz_face_type);

    MPI_Type_vector(nx, nz, (ny + 2) * (nz + 2), MPI_DOUBLE, &xz_face_type);
    MPI_Type_commit(&xz_face_type);

    MPI_Type_vector(nx * ny, 1, nz + 2, MPI_DOUBLE, &xy_face_type);
    MPI_Type_commit(&xy_face_type);

    for (int iter = 0; iter < N_ITER; iter++) {
        int north, south, east, west, up, down;
        MPI_Cart_shift(cart_comm, 0, 1, &west, &east);
        MPI_Cart_shift(cart_comm, 1, 1, &north, &south);
        MPI_Cart_shift(cart_comm, 2, 1, &down, &up);

        if (west != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[1 * (ny + 2) * (nz + 2)], 1, yz_face_type, west, 0,
                         &subgrid[0], 1, yz_face_type, west, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        if (east != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[nx * (ny + 2) * (nz + 2)], 1, yz_face_type, east, 0,
                         &subgrid[(nx + 1) * (ny + 2) * (nz + 2)], 1, yz_face_type, east, 0, 
                         cart_comm, MPI_STATUS_IGNORE);
        }

        if (north != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[1 * (nz + 2)], 1, xz_face_type, north, 0,
                         &subgrid[0], 1, xz_face_type, north, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        if (south != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[ny * (nz + 2)], 1, xz_face_type, south, 0,
                         &subgrid[(ny + 1) * (nz + 2)], 1, xz_face_type, south, 0,
                         cart_comm, MPI_STATUS_IGNORE);
        }

        if (down != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[1], 1, xy_face_type, down, 0,
                         &subgrid[0], 1, xy_face_type, down, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        if (up != MPI_PROC_NULL) {
            MPI_Sendrecv(&subgrid[nz], 1, xy_face_type, up, 0,
                         &subgrid[nz + 1], 1, xy_face_type, up, 0, cart_comm, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i <= nx; i++) {
            for (int j = 1; j <= ny; j++) {
                for (int k = 1; k <= nz; k++) {
                    int idx = i * (ny + 2) * (nz + 2) + j * (nz + 2) + k;
                    new_subgrid[idx] = (subgrid[idx - (ny + 2) * (nz + 2)] +
                                        subgrid[idx + (ny + 2) * (nz + 2)] +
                                        subgrid[idx - (nz + 2)] +
                                        subgrid[idx + (nz + 2)] +
                                        subgrid[idx - 1] +
                                        subgrid[idx + 1]) / 6.0;
                }
            }
        }

        for (int i = 1; i <= nx; i++)
            for (int j = 1; j <= ny; j++)
                for (int k = 1; k <= nz; k++) {
                    int idx = i * (ny + 2) * (nz + 2) + j * (nz + 2) + k;
                    subgrid[idx] = new_subgrid[idx];
                }
    }

    gettimeofday(&end, NULL);

    double local_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    double global_time;

    MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Execution time: %f seconds\n", global_time);
    }

    free(subgrid);
    free(new_subgrid);

    MPI_Type_free(&yz_face_type);
    MPI_Type_free(&xz_face_type);
    MPI_Type_free(&xy_face_type);

    MPI_Finalize();
    return 0;
}
