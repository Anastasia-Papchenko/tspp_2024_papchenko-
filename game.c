#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdbool.h>

#define ALIVE 1
#define DEAD 0

void initialize_grid(int *grid, int rows, int cols, int rank, int size) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            grid[i * cols + j] = (i == rank && j == 2) ? ALIVE : DEAD; 
        }
    }
}

void exchange_borders(int *local_grid, int local_rows, int cols, int rank, int size, MPI_Request *requests) {
    int top_rank = (rank == 0) ? size - 1 : rank - 1;
    int bottom_rank = (rank == size - 1) ? 0 : rank + 1;

    MPI_Isend(&local_grid[0], cols, MPI_INT, top_rank, 0, MPI_COMM_WORLD, &requests[0]);
    MPI_Isend(&local_grid[(local_rows - 1) * cols], cols, MPI_INT, bottom_rank, 1, MPI_COMM_WORLD, &requests[1]);
    MPI_Irecv(&local_grid[-cols], cols, MPI_INT, top_rank, 1, MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(&local_grid[local_rows * cols], cols, MPI_INT, bottom_rank, 0, MPI_COMM_WORLD, &requests[3]);
}

int count_neighbors(int *grid, int rows, int cols, int x, int y) {
    int count = 0;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            if (dx == 0 && dy == 0) continue;
            int nx = (x + dx + rows) % rows;
            int ny = (y + dy + cols) % cols;
            count += grid[nx * cols + ny];
        }
    }
    return count;
}

void update_grid(int *grid, int *new_grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int neighbors = count_neighbors(grid, rows, cols, i, j);
            if (grid[i * cols + j] == ALIVE) {
                new_grid[i * cols + j] = (neighbors == 2 || neighbors == 3) ? ALIVE : DEAD;
            } else {
                new_grid[i * cols + j] = (neighbors == 3) ? ALIVE : DEAD;
            }
        }
    }
}

bool check_stability(int current_alive, int prev_alive) {
    return current_alive == prev_alive;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <rows> <cols> <iterations>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    int max_iterations = atoi(argv[3]);

    int local_rows = rows / size + (rank < rows % size ? 1 : 0);
    int *local_grid = (int *)malloc((local_rows + 2) * cols * sizeof(int)); 
    int *new_grid = (int *)malloc(local_rows * cols * sizeof(int));

    initialize_grid(&local_grid[cols], local_rows, cols, rank, size);

    MPI_Barrier(MPI_COMM_WORLD);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    MPI_Request requests[4];
    int global_alive = 0, local_alive, prev_local_alive = -1;

    for (int iter = 0; iter < max_iterations; ++iter) {
      
        exchange_borders(&local_grid[cols], local_rows, cols, rank, size, requests);
        MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

  
        update_grid(&local_grid[cols], new_grid, local_rows, cols);

      
        local_alive = 0;
        for (int i = 0; i < local_rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (new_grid[i * cols + j] == ALIVE) {
                    local_alive++;
                }
            }
        }

       
        if (check_stability(local_alive, prev_local_alive)) {
            break;
        }
        prev_local_alive = local_alive;

        int *temp = local_grid;
        local_grid = new_grid;
        new_grid = temp;

        MPI_Allreduce(&local_alive, &global_alive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (global_alive == 0) break;
    }

    gettimeofday(&end, NULL);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    if (rank == 0) {
        printf("Total number of living cells: %d\n", global_alive);
        printf("Execution time: %f seconds\n", elapsed_time);
    }

    free(local_grid);
    free(new_grid);

    MPI_Finalize();
    return 0;
}

