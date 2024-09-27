#include "pthread.h"
#include <queue>
#include <cassert>
#include <iostream>
#include <semaphore.h>

class MyConcurrentQueue {
public:
    MyConcurrentQueue() : queue_limit(0) {}

    void reserve(int n) {
        queue_limit = n;
        sem_init(&empty_slots, 0, n); // семафор для свободных мест
        sem_init(&filled_slots, 0, 0); // семафор для занятых мест
        pthread_mutex_init(&mutex, NULL); 
    }

    void put(int value) {
        sem_wait(&empty_slots); // ждем когда появятся свободные места

        pthread_mutex_lock(&mutex); 
        queue.push(value); 
        pthread_mutex_unlock(&mutex);

        sem_post(&filled_slots); // Увеличиваем счетчик занятых мест
    }

    int get() {
        sem_wait(&filled_slots); // ждем когда появятся занятые места 

        pthread_mutex_lock(&mutex); 
        int value = queue.front(); 
        queue.pop(); 
        pthread_mutex_unlock(&mutex);

        sem_post(&empty_slots); // Увеличиваем счетчик свободных мест
        return value; 
    }

private:
    int queue_limit;
    std::queue<int> queue;
    
    sem_t empty_slots; // Семафор для свободных мест в очереди
    sem_t filled_slots; // Семафор для занятых мест в очереди
    pthread_mutex_t mutex; // Мьютекс для защиты доступа к очереди
};
