#include <cstdlib>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <chrono>
#include "fitting_tree.hpp"

#define Filename "../data/SOSD_Dataset/books_800M_uint64/books_800M_uint64"

typedef uint32_t K;

void load_data(std::vector<K> &data) {
    std::ifstream in(Filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "Error opening file: " << Filename << std::endl;
        return;
    }

    K n_keys;
    in.read(reinterpret_cast<char*>(&n_keys), sizeof(K));  
    data.resize(n_keys);
    in.read(reinterpret_cast<char*>(data.data()), n_keys * sizeof(K));  
    in.close();
}

template<typename T>
void testForModel(std::vector<K> data, T &fitting_tree) {
    double time_sum = 0;
    double avg_time_sum = 0;
    
    for(int k = 0; k < 3; k++) {
        time_sum = 0;
        for(int i = 0; i < 1000; i++) {
            int index = rand() % data.size();
            auto start_search = std::chrono::high_resolution_clock::now();
            fitting_tree.search(data[index]);
            auto end_search = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> search_duration = end_search - start_search;
            double time = search_duration.count();
            
            if(time > 1e-5) {
                printf("搜索 %u 用时 %f 秒\n", data[index], time);
            }
            time_sum += search_duration.count();
        }
        
        auto avg_time = time_sum / 1000;
        avg_time_sum += avg_time;
        std::cout << "平均搜索时间: " << avg_time << " 秒" << std::endl;
        std::cout << std::endl;
    }
    
    auto avg_time = avg_time_sum / 3;
    std::cout << "最终平均搜索时间: " << avg_time << " 秒" << std::endl;
    std::cout << "索引大小: " << fitting_tree.size_in_bytes() << " 字节" << std::endl;
    std::cout << "树高度: " << fitting_tree.height() << std::endl;
    
    // 获取线段信息
    std::cout << "线段数量: " << fitting_tree.get_segment_count() << std::endl;
    std::cout << "线段占用字节数: " << fitting_tree.get_segment_size_bytes() << " 字节" << std::endl;
}

template<size_t epsilon = 2>
void experiment_fitting(std::vector<K> data, int threads_num = 16) {
    using FITTING_Index = fitting::FITTING_Tree<K, epsilon>;
    
    auto start = std::chrono::high_resolution_clock::now();
    FITTING_Index fitting_tree;
    fitting_tree.build(data);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    
    std::cout << "构建索引用时: " << duration.count() << " 秒" << std::endl;
    testForModel<FITTING_Index>(data, fitting_tree);
}

int main() {
    std::vector<K> data;
    load_data(data);
    printf("数据大小: %ld\n", data.size());
    
    constexpr size_t epi = 64;
    experiment_fitting<epi>(data);
    
    return 0;
} 