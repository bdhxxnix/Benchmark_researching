
#include <cstdlib>
#include <fstream>
#include "pgm_index.hpp"
#include "greedyplrBasedPGM.hpp"
#include "swingBasedPGM.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

#define Filename "../data/osmdata.txt"
typedef uint32_t K;

struct Result {
	double build_time;
	double search_time;
	size_t index_size;
};

struct Results{
	std::vector<double> build_time;
	std::vector<double> search_time;
	std::vector<size_t> index_size;
	
};

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

template <typename T>
double searchForKey(std::vector<K> data, T model, int parallelism = 1) {

	// Randomly pick 1000 elements from the data and search for them
	
	double time_sum = 0;
	double avg_time_sum = 0;


	for(int k = 0; k< 3 ; k++){
		time_sum = 0;
		# pragma omp parallel for reduction(+:time_sum) num_threads(parallelism)
		for(int i = 0;i<1000;i++){
			int index = rand() % data.size();
			auto start_search = std::chrono::high_resolution_clock::now();
			model.search(data[index]);
			auto end_search = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> search_duration = end_search - start_search;
			double time = search_duration.count();	
			
			time_sum += time;
		}
		auto avg_time  = time_sum / 1000;
		avg_time_sum += avg_time;
		// std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
	}
	auto avg_time = avg_time_sum / 3;

    return avg_time;
}


template<typename T>
Result testForModel(std::vector<K> data) {
    Result result;	
	// Get the time of constructing the whole PGM-Index
	auto start = std::chrono::high_resolution_clock::now();
	T model(data.begin(), data.end());
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	double time = duration.count();
	result.build_time = time;

	result.search_time = searchForKey(data, model);
	result.index_size = model.size_in_bytes();
	
	return result;
}

void printout(Results results){
	auto build_time_vector = results.build_time;
	auto search_time_vector = results.search_time;
	auto index_size_vector = results.index_size;
	printf("The build time is :");
	for(size_t i = 0 ;i<build_time_vector.size();i++){
		std::cout<< build_time_vector[i] <<",";
	}
	printf("\nThe search time is: ");
	for(size_t i = 0 ;i < search_time_vector.size();i++){
		std::cout << search_time_vector[i] << ",";
	}
	printf("\nThe index size is :");
	for(size_t i = 0 ; i <index_size_vector.size();i++){
		std::cout << index_size_vector[i]<< ",";
	}
	printf("\n");
}

void experiment_linear(std::vector<K> data, int threads_num = 16) {
	
	
	Results results;
	Result result;

	using PGM_Index1 = pgm::PGMIndex<K, 4, 4>;
	result = testForModel<PGM_Index1>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index2 = pgm::PGMIndex<K, 8, 8>;
	result = testForModel<PGM_Index2>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index3 = pgm::PGMIndex<K, 16, 16>;
	result = testForModel<PGM_Index3>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index4 = pgm::PGMIndex<K, 32, 32>;
	result = testForModel<PGM_Index4>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index5 = pgm::PGMIndex<K, 64, 64>;
	result = testForModel<PGM_Index5>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index6 = pgm::PGMIndex<K, 128, 128>;
	result = testForModel<PGM_Index6>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index7 = pgm::PGMIndex<K, 256, 256>;
	result = testForModel<PGM_Index7>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index8 = pgm::PGMIndex<K, 512, 512>;
	result = testForModel<PGM_Index8>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index9 = pgm::PGMIndex<K, 1024, 1024>;
	result = testForModel<PGM_Index9>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index10 = pgm::PGMIndex<K, 2048, 2048>;
	result = testForModel<PGM_Index10>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index11 = pgm::PGMIndex<K, 4096, 4096>;
	result = testForModel<PGM_Index11>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index12 = pgm::PGMIndex<K, 8192, 8192>;
	result = testForModel<PGM_Index12>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);



	printf("The results for OptimalPLR:\n");
	printout(results);
}

void experiment_greedy(std::vector<K> data, int threads_num = 16) {

	
	Results results;
	Result result;

	using PGM_Index1 = Greedy::PGMIndex<K, 4, 4>;
	result = testForModel<PGM_Index1>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index2 = Greedy::PGMIndex<K, 8, 8>;
	result = testForModel<PGM_Index2>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index3 = Greedy::PGMIndex<K, 16, 16>;
	result = testForModel<PGM_Index3>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index4 = Greedy::PGMIndex<K, 32, 32>;
	result = testForModel<PGM_Index4>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index5 = Greedy::PGMIndex<K, 64, 64>;
	result = testForModel<PGM_Index5>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index6 = Greedy::PGMIndex<K, 128, 128>;
	result = testForModel<PGM_Index6>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index7 = Greedy::PGMIndex<K, 256, 256>;
	result = testForModel<PGM_Index7>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index8 = Greedy::PGMIndex<K, 512, 512>;
	result = testForModel<PGM_Index8>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index10 = Greedy::PGMIndex<K, 1024, 1024>;
	result = testForModel<PGM_Index10>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index11 = Greedy::PGMIndex<K, 2048, 2048>;
	result = testForModel<PGM_Index11>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index12 = Greedy::PGMIndex<K, 4096, 4096>;
	result = testForModel<PGM_Index12>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index13 = Greedy::PGMIndex<K, 8192, 8192>;
	result = testForModel<PGM_Index13>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);


	printf("The results for greedyPLR:\n");
	printout(results);
}

void experiment_Swing(std::vector<K> data, int threads_num = 16) {
	
		
	Results results;
	Result result;

	using PGM_Index1 = PGM_S::PGMIndex<K, 4, 4>;
	result = testForModel<PGM_Index1>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index2 = PGM_S::PGMIndex<K, 8, 8>;
	result = testForModel<PGM_Index2>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index3 = PGM_S::PGMIndex<K, 16, 16>;
	result = testForModel<PGM_Index3>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index4 = PGM_S::PGMIndex<K, 32, 32>;
	result = testForModel<PGM_Index4>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index5 = PGM_S::PGMIndex<K, 64, 64>;
	result = testForModel<PGM_Index5>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index6 = PGM_S::PGMIndex<K, 128, 128>;
	result = testForModel<PGM_Index6>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index7 = PGM_S::PGMIndex<K, 256, 256>;
	result = testForModel<PGM_Index7>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index8 = PGM_S::PGMIndex<K, 512, 512>;
	result = testForModel<PGM_Index8>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index9 = PGM_S::PGMIndex<K, 1024, 1024>;
	result = testForModel<PGM_Index9>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index10 =PGM_S::PGMIndex<K, 2048, 2048>;
	result = testForModel<PGM_Index10>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index11 =PGM_S::PGMIndex<K, 4096, 4096>;
	result = testForModel<PGM_Index11>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);

	using PGM_Index12 =PGM_S::PGMIndex<K, 8192, 8192>;
	result = testForModel<PGM_Index12>(data);
	results.build_time.push_back(result.build_time);
	results.search_time.push_back(result.search_time);
	results.index_size.push_back(result.index_size);



	printf("The results for SwingFilter:\n");
	printout(results);
}



int main() {
	// Load the data
	std::vector<K> data;
	load_data(data);
	printf("The size of the data is %ld\n",data.size());
	
	experiment_Swing(data);
	// experiment_greedy(data);
	// experiment_linear(data);
	return 0;
}
