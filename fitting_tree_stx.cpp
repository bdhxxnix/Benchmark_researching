#include "fitting_Tree_stx.hpp"
#include "greedy.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

typedef double K; // Define the type of the data, can be changed to any other type
#define FILE_NAME "../data/fbdata.txt" // Define the file name to load the data from

template<typename K>
void load_data(std::vector<K> &data) {

	// Load data from a file or generate it
	std::ifstream file(FILE_NAME);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << FILE_NAME << std::endl;
        
		return;
	}
	K value;
	while(file >> value) { // Limit to 1 million entries for testing
		data.push_back(value);
	}
	file.close();
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
                printf("搜索 %f 用时 %f 秒\n", data[index], time);
            }
            time_sum += search_duration.count();
            if(i<10){
                printf("The predicted position is %d\n", fitting_tree.search(data[index]));
                printf("The real position is %d\n", index);
            }
        }
        
        auto avg_time = time_sum / 1000;
        avg_time_sum += avg_time;
        std::cout << "平均搜索时间: " << avg_time << " 秒" << std::endl;
        std::cout << std::endl;
    }
    
    auto avg_time = avg_time_sum / 3;
    std::cout << "最终平均搜索时间: " << avg_time << " 秒" << std::endl;
    std::cout << "索引大小: " << fitting_tree.size_in_bytes() << " 字节" << std::endl;
    
    // 获取线段信息
    std::cout << "线段数量: " << fitting_tree.get_segment_count() << std::endl;
    std::cout << "线段占用字节数: " << fitting_tree.get_segment_size_bytes() << " 字节" << std::endl;
}

void experiment_fitting(std::vector<K> data, int threads_num = 16) {
    using FITTING_Index = FittingTree<K>;
    
    auto start = std::chrono::high_resolution_clock::now();
    FITTING_Index fitting_tree;
    fitting_tree.build(data);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    
    std::cout << "构建索引用时: " << duration.count() << " 秒" << std::endl;
    testForModel<FITTING_Index>(data, fitting_tree);
}

int main() {

    // Test the FittingTree class
    std::vector<K> data; 
    load_data(data); 
    
    experiment_fitting(data); 
    
    return 0;
}