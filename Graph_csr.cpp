#include <iostream>
#include <vector>
#include "papi.h"
#include <algorithm>

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

//Поиск вершины с наибольшим рангом с помощью сортировка
    int find_highest_weight_vertex() {
        std::vector<std::pair<int, double>> weights(row_count);

        for (int i = 0; i < row_count; i++) {
            weights[i] = {i, calculate_weight(i)};
        }

        // Сортировка по весу в порядке убывания
        std::sort(weights.begin(), weights.end(), [](const auto& a, const auto& b) {
            return a.second > b.second; 
        });

        return weights.front().first; // Возвращаем индекс вершины с наибольшим рангом
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

    int retval;
    int EventCode;
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
        
        if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
            std::cerr << "PAPI library init error!" << std::endl;
            return -1;
        }




        long long values[2];
        int event_set = PAPI_NULL;

        retval = PAPI_event_name_to_code("PAPI_TOT_INS", &event_set);
        if ( retval != PAPI_OK) {
            printf("PAPI_event_name_to_code error %d n", retval);
        }

        if (PAPI_query_event(EventCode) != PAPI_OK) {
            printf("Can not measure %d event ", EventCode);
        } else {
            printf("Will be succeeded in measuring %d event ", EventCode);
        }
        
        if (PAPI_create_eventset(&event_set) != PAPI_OK) {
            printf("PAPI_create_eventset error!n");
            // exit(1);
        }

        // values[0]
        
        PAPI_add_event(event_set, PAPI_L1_TCM);
       
        // values[1]

        if (PAPI_add_event(event_set, PAPI_L2_DCM) != PAPI_OK) {
            std::cerr << "L2 Cache event not available." << std::endl;
        }

        // values[3]
    
        if (PAPI_add_event(event_set, EventCode) != PAPI_OK) {
            printf("PAPI_add_event error !n");
            // exit(1);
        }

        unsigned int native = 0x0;
        PAPI_event_info_t info;

        if (PAPI_get_event_info(EventCode, &info) != PAPI_OK) {
            printf("Error in get event_infon");
            // exit(1);
        }

        printf("n%d, %s, count: %sn     ", info.event_code, info.symbol, info.short_descr);


        // Начало отслеживания событий
        PAPI_start(event_set);

        // Алгоритм для нахождения вершины с наибольшим рангом
        int highest_rank_vertex = a.find_highest_rank_vertex();
        // int highest_weight_vertex = a.find_highest_weight_vertex();

        // Остановка отслеживания событий
        PAPI_stop(event_set, values);

        std::cout << "Название теста: " << filenames[n_test] << ", Алгоритм 1 - Индекс вершины с наибольшим рангом: " << highest_rank_vertex 
                  << ", Кэш-промах L1: " << values[0] << ", Кэш-промах L2: " << (PAPI_add_event(event_set, PAPI_L2_DCM) == PAPI_OK ? values[1] : 0)  
                  << std::endl;
                //    << ", Общее количество циклов: " << values[2]
        printf("values[2] = %lldn", values[2]);
        // std::cout << "Тест: " << filenames[n_test] 
        //            << ", Алгоритм 2 - Индекс вершины с наибольшим рангом: " << highest_weight_vertex 
        //            << ", Кэш-промах L1: " << values[0] 
        //            << ", Кэш-промах L2: " << (PAPI_add_event(event_set, PAPI_L2_ICA) == PAPI_OK ? values[1] : 0)
        //            << ", Общее количество циклов: " << values[2] << std::endl;

        // Сброс графа для следующего теста
        // a.reset();

         PAPI_start(event_set);

        // Алгоритм для нахождения вершины с наибольшим рангом
        // int highest_rank_vertex = a.find_highest_rank_vertex();
        int highest_weight_vertex = a.find_highest_weight_vertex();

        // Остановка отслеживания событий
        PAPI_stop(event_set, values);

        // std::cout << "Название теста: " << filenames[n_test] << ", Алгоритм 1 - Индекс вершины с наибольшим рангом: " << highest_rank_vertex 
        //           << ", Кэш-промах L1: " << values[0] << ", Кэш-промах L2: " << (PAPI_add_event(event_set, PAPI_L2_ICA) == PAPI_OK ? values[1] : 0)  
        //           << ", Общее количество циклов: " << values[2] << std::endl;

        std::cout << "Тест: " << filenames[n_test] 
                   << ", Алгоритм 2 - Индекс вершины с наибольшим рангом: " << highest_weight_vertex 
                   << ", Кэш-промах L1: " << values[0] 
                   << ", Кэш-промах L2: " << (PAPI_add_event(event_set, PAPI_L2_DCM) == PAPI_OK ? values[1] : 0)
                    << std::endl;
                    // << ", Общее количество циклов: " << values[2]
        printf("values[2] = %lldn", values[2]);
        // Сброс графа для следующего теста
        a.reset();
    }

    return 0;
}