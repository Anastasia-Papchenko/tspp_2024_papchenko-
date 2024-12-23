#include <iostream>
#include <vector>
#include <cstdio> // Для работы с fopen, fread, fclose
#include "papi.h"
#include <algorithm>

class CSR_graph {
    int row_count; // количество вершин
    unsigned int col_count; // количество ребер

    std::vector<unsigned int> row_ptr;
    std::vector<int> col_ids;
    std::vector<double> vals;

public:
    void read(const char* filename) {
        FILE *graph_file = fopen(filename, "rb");
        if (!graph_file) {
            std::cerr << "Ошибка открытия файла: " << filename << std::endl;
            return;
        }
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

    double calculate_weight(int vertex) {
        double weight = 0.0;
        for (int i = row_ptr[vertex]; i < row_ptr[vertex + 1]; i++) {
            weight += vals[i];
        }
        return weight;
    }

    void reset() {
        row_count = 0;
        col_count = 0;
        row_ptr.clear();
        col_ids.clear();
        vals.clear();
    }

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

    int find_highest_rank_vertex() {
        int max_vertex = -1;
        float max_weight = -1, sum_weight;

        for (int i = 0; i < row_count; i++) {
            sum_weight = 0;
            for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
                if (col_ids[j] % 2 == 0) sum_weight += vals[j];
            }
            if (sum_weight > max_weight) {
                max_weight = sum_weight;
                max_vertex = i;
            }
        }
        return max_vertex;
    }

    int find_highest_weight_vertex() {
        int max_rank_vertex = -1;
        float max_rank = -1, rank, W_vert;

        for (int i = 0; i < row_count; i++) {
            rank = 0;
            for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
                W_vert = 0;
                for (int k = row_ptr[col_ids[j]]; k < row_ptr[col_ids[j] + 1]; k++) {
                    W_vert += vals[k] * (row_ptr[col_ids[j] + 1] - row_ptr[col_ids[j]]);
                }
                rank += vals[j] * W_vert;
            }
            if (rank > max_rank) {
                max_rank = rank;
                max_rank_vertex = i;
            }
        }
        return max_rank_vertex;
    }
};

#define N_TESTS 5

int main() {
    const char* filenames[N_TESTS] = {
        "synt",
        "road_graph",
        "stanford",
        "youtube",
        "syn_rmat"
    };

    int Eventset = PAPI_NULL, code, retval;
    long long values[3];

    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        std::cerr << "Ошибка инициализации PAPI!" << std::endl;
        return -1;
    }
    PAPI_create_eventset(&Eventset);

    PAPI_add_event(Eventset, PAPI_L1_TCM);
    PAPI_add_event(Eventset, PAPI_L2_TCM);
    char event_name[] = "perf::PERF_COUNT_HW_CACHE_REFERENCES";
    PAPI_event_name_to_code(event_name, &code);
    PAPI_add_event(Eventset, code);

    for (int n_test = 0; n_test < N_TESTS; n_test++) {
        CSR_graph a;
        a.read(filenames[n_test]);

        PAPI_start(Eventset);
        std::cout << "Название теста: " << filenames[n_test] << std::endl;
        std::cout << "Алгоритм 1 - Индекс вершины с наибольшим рангом: " 
                  << a.find_highest_rank_vertex() + 1 << std::endl;

        PAPI_stop(Eventset, values);
        std::cout << "Кэш-промах L1: " << values[0] 
                  << ", Кэш-промах L2: " << values[1] 
                  << ", PERF_COUNT_HW_CACHE_REFERENCES: " << values[2] 
                  << std::endl;

        PAPI_reset(Eventset);
        PAPI_start(Eventset);
        std::cout << "Алгоритм 2 - Индекс вершины с наибольшим рангом: " 
                  << a.find_highest_weight_vertex() + 1 << std::endl;

        PAPI_stop(Eventset, values);
        std::cout << "Кэш-промах L1: " << values[0] 
                  << ", Кэш-промах L2: " << values[1] 
                  << ", PERF_COUNT_HW_CACHE_REFERENCES: " << values[2] 
                  << std::endl;
        std::cout << std::endl;
        a.reset();
    }

    PAPI_destroy_eventset(&Eventset);
    PAPI_shutdown();
    return 0;
}
