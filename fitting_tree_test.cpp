#include "swingFTree.hpp"
#include "greedyFTree.hpp"
#include "optimalFTree.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

typedef uint64_t K; // Define the type of the data, can be changed to any other type

struct Result {
	double build_time;
	double search_time;
	size_t index_size;
};

/**
 * @brief Translate the data from the file to the vector
 *
 * @param filename
 * @return std::vector<K>
 */
template<typename K>
std::vector<K> load_data(std::string FILE_NAME)
{

    /* Open file. */
    std::ifstream in(FILE_NAME, std::ios::binary);
    if (!in.is_open())
        exit(EXIT_FAILURE);

    /* Read number of keys. */
    K n_keys;
    in.read(reinterpret_cast<char *>(&n_keys), sizeof(K));

    /* Initialize vector. */
    std::vector<K> data;
    data.resize(n_keys);

    /* Read keys. */
    in.read(reinterpret_cast<char *>(data.data()), n_keys * sizeof(K));
    in.close();

    return data;
}

/**
 *  Print out the results of a certain algorithm based on epi variation
 */
void printout(std::vector<Result> results){
	
	std::vector<double> build_times;
	std::vector<double> search_times;
	std::vector<size_t> index_size;

	for(auto result : results){
		build_times.push_back(result.build_time);
		search_times.push_back(result.search_time);
		index_size.push_back(result.index_size);
		
	}

	printf("\nThe build time is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		printf("%f,",build_times[i]);
	}
	printf("\nThe search time(last) is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		std::cout<<search_times[i] << ",";
	}

	printf("\nThe index size is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		printf("%ld,",index_size[i]);
	}

	printf("\n");
}


template<typename T>
Result testForModel(std::vector<K> data, int epi) {
	
    Result result;	
	// Get the time of constructing the whole PGM-Index
	auto start = std::chrono::high_resolution_clock::now();
	T model(epi);
    model.build(data);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	double time = duration.count();
	result.build_time = time;
		
	
	// Randomly pick 1000 elements from the data and search for them
	std::vector<double> times;
	std::random_device rd;  // Seed source
    std::mt19937 gen(rd()); // Mersenne Twister RNG
    std::uniform_int_distribution<> distrib(0, data.size()-1);
	for(int k = 0; k< 10 ; k++){
		auto start = std::chrono::high_resolution_clock::now();
		for(int i = 0;i<1000;i++){
			srand(42);
			size_t index = distrib(gen);
			
			size_t pos = model.search(data[index]);
        
            size_t lo = pos - epi - 1 > 0 ? (pos - epi - 1) : 0;
            size_t hi = pos + epi + 1 < data.size() ? (pos + epi + 1) : data.size();
			
			auto it = std::prev(std::upper_bound(data.begin()+lo,data.begin()+hi,data[index]));
			if(it == data.end() || it == data.begin()) {
				std::cout << "Error: No valid segment found for key " << data[index] << std::endl;
				continue;
			}
			
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = end -start;
		time =duration.count();
		time = time/1000.0;
		times.push_back(time);
		// std::cout << "Average time taken to search: " << avg_time << " seconds" << std::endl;
	}
	
	std::sort(times.begin(),times.end());

	// Remove the extremes
	double time_sum  = 0;
	for(int i  = 1; i < 9 ; i++){
		time_sum += times[i];
	}

	// Calculate the average time
	result.search_time = time_sum / 8.0;
	result.index_size = model.size_in_bytes();
	
	return result;
}



void experiment_swing(std::vector<K> data, int threads_num = 16) {
    using FITTING_Index = Swing::FittingTree<K>;

    std::vector<Result> results;
	printf("The result for Swing Filter: \n");
    for(int i = 2 ;i<14;i++){
        int epi = 1<<i;
        
        auto result = testForModel<FITTING_Index>(data, epi);
        results.push_back(result);
    }
    printout(results);

}

void experiment_Greedy(std::vector<K> data, int threads_num = 16) {
    using FITTING_Index = Greedy::FittingTree<K>;

    std::vector<Result> results;
	printf("The result for GreedyPLA: \n");
    for(int i = 2 ;i<14;i++){
        int epi = 1<<i;
        
        auto result = testForModel<FITTING_Index>(data, epi);
        results.push_back(result);
    }
    printout(results);

}

void experiment_Optimal(std::vector<K> data, int threads_num = 16) {
    using FITTING_Index = Optimal::FittingTree<K>;

	std::vector<Result> results;	
	printf("The result for OptimalPLR: \n");
    for(int i = 2 ;i<14;i++){
        int epi = 1<<i;
        
        auto result = testForModel<FITTING_Index>(data, epi);
        results.push_back(result);
    }
    printout(results);

}

int main(int argc, char* argv[]) {
	// Check if the dataset path is provided
	if (argc < 2) {
        std::cerr << "Usage: ./fitting_tree_test <dataset_path>\n";
        return 1;
    }

    std::string file_path = argv[1];

    // Test the FittingTree class
    std::vector<K> data; 
    data = load_data<K>(file_path); 

	printf("We are running on FIT test\n");
	printf("The size of the data, %ld\n",data.size());
    
	experiment_swing(data); 
	experiment_Greedy(data);
	experiment_Optimal(data);

	return 0;
}