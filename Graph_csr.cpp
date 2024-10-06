#include <iostream>
#include <vector>
#include "papi.h"

class CSR_graph {
    int row_count; //количество вершин
    unsigned int col_count; //количество ребер

    std::vector<unsigned int> row_ptr;
    std::vector<int> col_ids;
    std::vector<double> vals;

public:
    void read(const char* filename) {
        FILE *graph_file = fopen(filename, "rb");
        fread(reinterpret_cast<char*>(&row_count), sizeof(int), 1, graph_file);
        fread(reinterpret_cast<char*>(&col_count), sizeof(unsigned int), 1, graph_file);

        std::cout << "Количество вершин = " << row_count << ", Количество ребер = " << col_count << std::endl;

        row_ptr.resize(row_count + 1);
        col_ids.resize(col_count);
        vals.resize(col_count);

        fread(reinterpret_cast<char*>(row_ptr.data()), sizeof(unsigned int), row_count + 1, graph_file);
        fread(reinterpret_cast<char*>(col_ids.data()), sizeof(int), col_count, graph_file);
        fread(reinterpret_cast<char*>(vals.data()), sizeof(double), col_count, graph_file);
        fclose(graph_file);
    }

    //Метод, который принимает индекс вершины и возвращает ее вес
    double calculate_weight(int vertex) {
        double weight = 0.0;
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) {
            weight += vals[i];
        }
        return weight;
    }

    //Метод, который принимает индекс вершины и возвращает ее ранг
    double calculate_rank(int vertex) {
        double rank = 0.0;
        double W_vertex = calculate_weight(vertex); 

        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) {
            int neighbor = col_ids[i];
            double w_edge = vals[i];
            double W_neighbor = calculate_weight(neighbor); 
            rank += w_edge * W_neighbor;
        }
        return rank;
    }

    //Метод, который находит и возвращает индекс вершины с наибольшим рангом
    int find_highest_rank_vertex() {
        int max_vertex = -1;
        double max_rank = -1.0;

        for (int i = 0; i < row_count; i++) {
            double rank = calculate_rank(i);
            if (rank > max_rank) {
                max_rank = rank;
                max_vertex = i;
            }
        }
        return max_vertex;
    }

    //Сбрасываем состояние графа
    void reset() {
        row_count = 0;
        col_count = 0;
        row_ptr.clear();
        col_ids.clear();
        vals.clear();
    }
};

#define N_TESTS 5

int main() {
    //объеденила в один запуск
    const char* filenames[N_TESTS] = {
        "synt",
        "road_graph",
        "stanford",
        "youtube",
        "syn_rmat"
    };

    for (int n_test = 0; n_test < N_TESTS; n_test++) {
        CSR_graph a;
        a.read(filenames[n_test]);

        // Инициализация PAPI
        int retval;
        if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
            std::cerr << "PAPI library init error!" << std::endl;
            return -1;
        }

        long long values[3];
        int event_set = PAPI_NULL;
        
        PAPI_create_eventset(&event_set);
        
        // Добавление событий PAPI_L1_TCM и PAPI_L2_TCM
        PAPI_add_event(event_set, PAPI_L1_TCM);
       
        
        // Если поддерживается, добавляем PAPI_L2_TCM
        if (PAPI_add_event(event_set, PAPI_L2_TCM) != PAPI_OK) {
            std::cerr << "L2 Cache event not available." << std::endl;
        }
        
        // Добавление события PAPI_TOT_CYC
        PAPI_add_event(event_set, PAPI_TOT_CYC);


        // Начало отслеживания событий
        PAPI_start(event_set);

        // Алгоритм для нахождения вершины с наибольшим рангом
        int highest_rank_vertex = a.find_highest_rank_vertex();

        // Остановка отслеживания событий
        PAPI_stop(event_set, values);

        std::cout << "Название теста: " << filenames[n_test] << ", Индекс вершины с наибольшим рангом: " << highest_rank_vertex 
                  << ", Кэш-промах L1: " << values[0] << ", Кэш-промах L2: " << (PAPI_add_event(event_set, PAPI_L2_TCM) == PAPI_OK ? values[1] : 0)  
                  << ", Общее количество циклов: " << values[2] << std::endl;

        // Сброс графа для следующего теста
        a.reset();
    }

    return 0;
}
