
#include <cstdlib>
#include <fstream>
#include "pgm_index.hpp"
#include "greedyplrBasedPGM.hpp"
#include "optimalPLR_C_BasedPGM_V2.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

#define Filename "../data/osmdata.txt"
template<typename K>
void load_data(std::vector<K> &data) {

	// Load data from a file or generate it
	std::ifstream file(Filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << Filename << std::endl;
		return;
	}
	K value;
	while(file >> value) {
		data.push_back(value);
	}
	file.close();
}

template<typename T>
void testForModel(std::vector<double> data, T &pgm_index) {
	// Randomly pick 1000 elements from the data and search for them
	double time_sum = 0;
	double avg_time_sum = 0;
	for(int k = 0; k< 3 ; k++){
		for(int i = 0;i<1000;i++){
			int index = rand() % data.size();
			auto start_search = std::chrono::high_resolution_clock::now();
			auto result = pgm_index.search(data[index]);
			auto end_search = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> search_duration = end_search - start_search;
			double time = search_duration.count();	
			if(time > 1e-5){
				printf("The time taken to search for the %f is %f seconds\n",data[index],time);
			}
			time_sum += search_duration.count();
		}
		auto avg_time  = time_sum / 1000;
		avg_time_sum += avg_time;
		std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
        std::cout << std::endl;
	}
	auto avg_time = avg_time_sum / 3;
	std::cout << "Final Average time taken to search: " << avg_time << " seconds" << std::endl;

	// Print the size of the index
	std::cout << "Size of the index: " << pgm_index.size_in_bytes() << " bytes" << std::endl;
	// Print the number of segments of each level
	auto levels_offsets = pgm_index.get_levels_offsets();
	std::cout << "Number of segments: ";
	for (size_t i = 0; i < levels_offsets.size() - 1; ++i) {
		std::cout << levels_offsets[i + 1] - levels_offsets[i] << " ";
	}
	std::cout << std::endl;
	// Print the height of the index
	std::cout << "Height of the index: " << pgm_index.height() << std::endl;
}

template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_linear(std::vector<double> data, int threads_num = 16) {
	// Build the whole PGM-Index
	using PGM_Index = pgm::PGMIndex<double, epsilon, epsilon_recursive>;
	// Get the time of constructing the whole PGM-Index
	auto start = std::chrono::high_resolution_clock::now();
	PGM_Index pgm_index(data.begin(), data.end());
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "Time taken to build the index(linear): " << duration.count() << " seconds" << std::endl;

	// Test the model
	testForModel<PGM_Index>(data, pgm_index);
}

template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_greedy(std::vector<double> data, int threads_num = 16) {
	// Build the whole PGM-Index
	using PGM_Index = Greedy::PGMIndex<double, epsilon, epsilon_recursive>;
	// Get the time of constructing the whole PGM-Index
	auto start = std::chrono::high_resolution_clock::now();
	PGM_Index greedy_model(data.begin(), data.end());
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "Time taken to build the index(greedy): " << duration.count() << " seconds" << std::endl;


	// Test the model
	testForModel<PGM_Index>(data, greedy_model);
}


template<size_t epsilon = 2, size_t epsilon_recursive = 4>
void experiment_optimal_v2(std::vector<double> data, int threads_num = 16) {
	// Build the whole PGM-Index
	using PGM_Index = PGM_C2::PGMIndex<double, epsilon, epsilon_recursive>;
	// Get the time of constructing the whole PGM-Index
	auto start = std::chrono::high_resolution_clock::now();
	PGM_Index opt(data);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "Time taken to build the index(OptimalPLR_C_BasedPGM_V2): " << duration.count() << " seconds" << std::endl;

	// Test the model
	testForModel<PGM_Index>(data, opt);
}

int main() {
	// Load the data
	std::vector<double> data;
	load_data(data);
	printf("The size of the data is %ld\n",data.size());
	constexpr size_t epi = 64;
	constexpr size_t epi_recursive = 8;
	experiment_greedy<epi, epi_recursive>(data);
	return 0;
}
