#include <stdio.h>
#include <stdlib.h>
#include <papi.h>

typedef struct {
    int num_vertices;// вершины
    int *row_ptr;// Указатели на начало ряда для каждой вершины
    int *col_ind;// Индексы инцидентных вершин
    float *edge_weights;//массив весов ребер  
} Graph;

//рассчер ранга вершин
float calculate_rank(int vertex, Graph *g) {
    float rank = 0.0;
    
    for (int i = g->row_ptr[vertex]; i < g->row_ptr[vertex + 1]; ++i) {
        int neighbor = g->col_ind[i];//получаем индекс соседней вершины
        float W_neigh = 0.0;
        
        // Считаем W для соседней вершины
        for (int j = g->row_ptr[neighbor]; j < g->row_ptr[neighbor + 1]; ++j) {
            int second_neighbor = g->col_ind[j];
            // Обновляем значение W для соседней вершины, умножая вес ребра на количество инцидентных рёбер второго соседа.
            W_neigh += g->edge_weights[j] * (g->row_ptr[second_neighbor + 1] - g->row_ptr[second_neighbor]);
        }
        
        // Расчет ранга
        rank += g->edge_weights[i] * W_neigh;
    }
    
    return rank;
}

// Функция для нахождения вершины с наибольшим весом инцидентных рёбер
// делаем указатель на граф 
int find_max_weight_even(Graph *g) { 
    float max_weight = 0.0;
    int max_vertex = -1;

    for (int v = 0; v < g->num_vertices; v++) {
        if (v % 2 == 0) {  // Если вершина четная
            float weight = 0.0;

            for (int i = g->row_ptr[v]; i < g->row_ptr[v + 1]; ++i) {
                weight += g->edge_weights[i]; // Суммируем веса инцидентных рёбер
            }

            if (weight > max_weight) {
                max_weight = weight;
                max_vertex = v;
            }
        }
    }
    
    return max_vertex;
}

// Основная функция
int main() {

    int PAPI_library_init(int 6.0.0.0);
    // Загрузка графа из CSR
    Graph g;
    // здесь нужно заполнить g.row_ptr, g.col_ind, g.edge_weights... :(

    // Начало измерений PAPI
    int retval;
    long long values[2];
    
    if (PAPI_start((int[]){PAPI_L1_TCM, PAPI_L2_TCM}, 2) != PAPI_OK) {
        fprintf(stderr, "Ошибка инициализации PAPI\n");
        return 1;
    }

    // Находим вершину с наибольшим весом и вершину с наибольшим рангом
    int max_weight_vertex = find_max_weight_even(&g);
    float max_rank = 0.0;
    int max_rank_vertex = -1;

    for (int v = 0; v < g.num_vertices; v++) {
        float rank = calculate_rank(v, &g);
        if (rank > max_rank) {
            max_rank = rank;
            max_rank_vertex = v;
        }
    }

    // Остановка счётчиков PAPI и вывод результатов
    PAPI_stop(values, 2);
    
    printf("Вершина с наибольшим весом (чётная вершина): %d\n", max_weight_vertex);
    printf("Вершина с наибольшим рангом: %d\n", max_rank_vertex);
    printf("PAPI_L1_TCM: %lld\n", values[0]);
    printf("PAPI_L2_TCM: %lld\n", values[1]);
    

    return 0;
}