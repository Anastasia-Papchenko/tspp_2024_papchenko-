#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <vector>

void merge(int* arr, int* temp, int left, int mid, int right) {
    int i = left; 
    int j = mid + 1; 
    int k = left;   

    while (i <= mid && j <= right) {
        if (arr[i] <= arr[j]) {
            temp[k++] = arr[i++];
        } else {
            temp[k++] = arr[j++];
        }
    }

    while (i <= mid) {
        temp[k++] = arr[i++];
    }

    while (j <= right) {
        temp[k++] = arr[j++];
    }

    for (i = left; i <= right; i++) {
        arr[i] = temp[i];
    }
}

void merge_sort(int* arr, int* temp, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        #pragma omp task shared(arr, temp) if (right - left > 1000)
        merge_sort(arr, temp, left, mid);

        #pragma omp task shared(arr, temp) if (right - left > 1000)
        merge_sort(arr, temp, mid + 1, right);

        #pragma omp taskwait
        merge(arr, temp, left, mid, right);
    }
}

int cmpfunc(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}


bool check(int* arr, int N) {
    for (int i = 0; i < N - 1; ++i) {
        if (arr[i] > arr[i + 1]) {
            return false; 
        }
    }
    return true; 
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <N> <p>n", argv[0]);
        return -1;
    }

    int N = atoi(argv[1]);
    int p = atoi(argv[2]);
    int* arr = (int*)malloc(N * sizeof(int));
    int* arr_copy = (int*)malloc(N * sizeof(int));
    int* temp = (int*)malloc(N * sizeof(int)); 

    srand(omp_get_wtime());
    for (int i = 0; i < N; i++) {
        arr[i] = rand() % 1000000; 
    }


    for (int i = 0; i < N; i++) {
        arr_copy[i] = arr[i];
    }

    
    struct timeval start_qsort, end_qsort;
    gettimeofday(&start_qsort, NULL);
    qsort(arr_copy, N, sizeof(int), cmpfunc);
    gettimeofday(&end_qsort, NULL);
    
    double time_qsort = (end_qsort.tv_sec - start_qsort.tv_sec) +
                        (end_qsort.tv_usec - start_qsort.tv_usec) / 1000000.0;

    
    struct timeval start_merge, end_merge;
    gettimeofday(&start_merge, NULL);
    
    
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            merge_sort(arr, temp, 0, N - 1);
        }
    }

    gettimeofday(&end_merge, NULL);
    double time_merge = (end_merge.tv_sec - start_merge.tv_sec) +
                        (end_merge.tv_usec - start_merge.tv_usec) / 1000000.0;


    // for (int i = 0; i < N; i++) {
    //     printf("%d ", arr[i]);
    // }
    // printf("\n");

    if (check(arr, N)) {
        std::cout << "Массив не убывает (arr)." << std::endl;
    } else {
        std::cout << "Массив убывает (arr)." << std::endl;
    }

    if (check(arr_copy, N)) {
        std::cout << "Массив не убывает (arr_copy)." << std::endl;
    } else {
        std::cout << "Массив убывает (arr_copy)." << std::endl;
    }


    printf("Time qsort: %lf sec\n", time_qsort);
    printf("Multithreaded sorting time: %lf sec\n", time_merge);

    free(arr);
    free(arr_copy);
    free(temp); 
    return 0;
}
