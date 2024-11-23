#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#define ALIVE 1
#define DEAD 0


void initialize_grid(int *grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            grid[i * cols + j] = rand() % 2; 
        }
    }
}


int count_neighbors(int *grid, int rows, int cols, int x, int y) {
    int count = 0;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            if (!(i == 0 && j == 0)) {
                int ni = (x + i + rows) % rows; 
                int nj = (y + j + cols) % cols;
                count += grid[ni * cols + nj];
            }
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
    int iterations = atoi(argv[3]); 
    int local_rows = rows / size; 

    int *grid = NULL;
    int *new_grid = (int *)malloc(local_rows * cols * sizeof(int));
    int *local_grid = (int *)malloc(local_rows * cols * sizeof(int));

   
    if (rank == 0) {
        grid = (int *)malloc(rows * cols * sizeof(int));
        initialize_grid(grid, rows, cols);
    }

    
    MPI_Scatter(grid, local_rows * cols, MPI_INT, local_grid, local_rows * cols, MPI_INT, 0, MPI_COMM_WORLD);

    int game_active = 1;
    int global_alive_count = 0;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    for (int iteration = 0; iteration < iterations && game_active; ++iteration) {
        
        update_grid(local_grid, new_grid, local_rows, cols);

        
        int local_alive_count = 0;
        for (int i = 0; i < local_rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (local_grid[i * cols + j] == ALIVE) {
                    local_alive_count++;
                }
            }
        }

        
        MPI_Allreduce(&local_alive_count, &global_alive_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        
        static int previous_alive_count = -1;
        if (previous_alive_count == global_alive_count) {
            game_active = 0;
        }
        previous_alive_count = global_alive_count;

        
        int *temp = local_grid;
        local_grid = new_grid;
        new_grid = temp;
    }

    gettimeofday(&end, NULL);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    if (rank == 0) {
        printf("Total number of living cells: %d\n", global_alive_count);
        printf("Execution time: %f sec\n", elapsed_time);
    }

    free(new_grid);
    free(local_grid);
    if (rank == 0) {
        free(grid);
    }

    MPI_Finalize();
    return 0;
}
