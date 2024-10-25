#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void random_walk(int a, int b, double p, int x, int* count_a, int* count_b) {
    int position = x;
 

    while (position > a && position < b) {
        double rand_val = (double) rand() / RAND_MAX; 
        if (rand_val < p) { 
            position++;
        } else {
            position--;
        }
    }

    if (position == a) {
        (*count_a)++;
    } else if (position == b) {
        (*count_b)++;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 7) {
        printf("Usage: %s a b p x N P\n", argv[0]);
        return 1;
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    double p = atof(argv[3]);
    int x = atoi(argv[4]);
    int N = atoi(argv[5]);
    int P = atoi(argv[6]);

    if (a >= x || x >= b) {
        printf("Убедитесь, что a < x < b.\n");
        return 1;
    }

    omp_set_num_threads(P);
    
    int count_a = 0;
    int count_b = 0;

    double start_time = omp_get_wtime();

   
    #pragma omp parallel
    {
       
        int private_count_a = 0;
        int private_count_b = 0;


        #pragma omp for
        for (int i = 0; i < N; i++) {
            random_walk(a, b, p, x, &private_count_a, &private_count_b);
        }


        #pragma omp atomic
        count_a += private_count_a;

        #pragma omp atomic
        count_b += private_count_b;
    }

    double end_time = omp_get_wtime();

    double prob_a = (double) count_a / N; 
    double prob_b = (double) count_b / N;
    double avg_time = (end_time - start_time) / N; 


    printf("number of particles in A: %d\n", count_a);
    printf("number of particles in B: %d\n", count_b);
    printf("probability of hitting A: %f\n", prob_a);
    printf("probability of hitting B: %f\n", prob_b);
    printf("average lifetime of one particle: %f sec\n", avg_time);
    printf("total working time: %f sec\n", end_time - start_time);

    return 0;
}