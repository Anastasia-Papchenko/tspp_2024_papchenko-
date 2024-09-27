#include "MyConcurrentQueue.h"
#include "pthread.h"
#include <queue>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <sys/time.h>
#include <semaphore.h>

MyConcurrentQueue my_queue;

int numberOfProducers;
int numberOfConsumers;

void* producer_func(void* params) {
    int thread_id = *((int*)params);
    int numberOfElements = *(int*)params + 100000 * numberOfConsumers;
    for (int i = 0; i < numberOfElements; i++) {
        my_queue.put(thread_id * 100000 + i); 
    }
    return NULL;
}

void *consumer_func(void *params) {
    int numberOfElements = *(int*)params + 100000 * numberOfProducers;
    for (int i = 0; i < numberOfElements; i++) {
        int a = my_queue.get();
        assert(a >= 0);
    }
    return NULL;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <number_of_producers> <number_of_consumers>" << std::endl;
        return 1;
    }

    numberOfProducers = std::atoi(argv[1]);
    numberOfConsumers = std::atoi(argv[2]);

    my_queue.reserve(1000);

    pthread_t* producers = new pthread_t[numberOfProducers];
    pthread_t* consumers = new pthread_t[numberOfConsumers];

    int* producer_ids = new int[numberOfProducers];
    int* consumer_ids = new int[numberOfConsumers];

    //фиксируем время начала
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);


    for (int i = 0; i < numberOfProducers; ++i) {
        producer_ids[i] = i; 
        pthread_create(&producers[i], NULL, producer_func, &producer_ids[i]);
    }


    for (int i = 0; i < numberOfConsumers; ++i) {
        pthread_create(&consumers[i], NULL, consumer_func, &consumer_ids[i]);
    }

    
    for (int i = 0; i < numberOfProducers; ++i) {
        pthread_join(producers[i], NULL);
    }

  
    for (int i = 0; i < numberOfConsumers; ++i) {
        pthread_join(consumers[i], NULL);
    }

    // фиксируем время конца
    gettimeofday(&end_time, NULL);

    
    long seconds = end_time.tv_sec - start_time.tv_sec;
    long micros = end_time.tv_usec - start_time.tv_usec;
    long elapsed_time = seconds * 1000000 + micros; 

    std::cout << "Elapsed time: " << elapsed_time << " microseconds" << std::endl;

   
   
    delete[] producer_ids;
    delete[] consumer_ids;

    return 0;
}



