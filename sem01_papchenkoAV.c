#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
// #include <time.h>
#include <sys/time.h>

double total_area = 0.0; // Общая площадь
pthread_mutex_t mutex;   // Мьютекс для защиты общей переменной

typedef struct {
    int thread_id;
    int num_threads;
    int num_rectangles;
    double width;
} thread_data_t;

void* calculate_area(void* threadarg) {
    thread_data_t* my_data = (thread_data_t*) threadarg;
    int thread_id = my_data->thread_id;
    int num_threads = my_data->num_threads;
    int num_rectangles = my_data->num_rectangles;
    double width = my_data->width;

    double local_area = 0.0;
    for (int i = thread_id; i < num_rectangles; i += num_threads) {
        double x = (i + 0.5) * width; // Центр прямоугольника
        local_area += 4.0 / (1.0 + x * x); // Высота прямоугольника
    }
    local_area *= width; // Умножаем на ширину прямоугольника

    pthread_mutex_lock(&mutex);
    total_area += local_area; 
    pthread_mutex_unlock(&mutex);

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Использование: %s <число отрезков> <количество нитей>\n", argv[0]);
        return 1;
    }

    int num_rectangles = atoi(argv[1]);
    int num_threads = atoi(argv[2]);

    if (num_rectangles <= 0 || num_threads <= 0) {
        fprintf(stderr, "Число отрезков и количество нитей должны быть положительными.\n");
        return 1;
    }

    pthread_t threads[num_threads];
    thread_data_t thread_data[num_threads];
    double width = 1.0 / num_rectangles;


    pthread_mutex_init(&mutex, NULL);


    // clock_t start_time = clock();
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Создание потоков
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].num_threads = num_threads;
        thread_data[i].num_rectangles = num_rectangles;
        thread_data[i].width = width;
        pthread_create(&threads[i], NULL, calculate_area, (void*)&thread_data[i]);
    }

    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

  
    pthread_mutex_destroy(&mutex);


    // clock_t end_time = clock();
    // double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC; // измеряем время и прииводим к секундам 
    gettimeofday(&end, NULL);

    
    double time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;


    // Вывод результата
    printf("Приблизительное значение π : %.10f\n", total_area);
    printf("Время работы: %.6f секунд\n", time_spent);

    return 0;
}