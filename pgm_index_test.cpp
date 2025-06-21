
#include <cstdlib>
#include <fstream>
#include "optimalPGM.hpp"
#include "greedyPGM.hpp"
#include "swingPGM.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <omp.h>
#include <random>
#include <chrono>

#define FILE_NAME "../data/SOSD_Dataset/books"

struct Result {
	double build_time;
	double search_time;
	size_t index_size;
	size_t height;
	std::vector<int> segOfLayer;
};


/**
 * @brief Translate the data from the file to the vector
 *
 * @param filename
 * @return std::vector<K>
 */
template<typename K>
std::vector<K> load_data()
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

    /* Sort the data in increasing order. */
    return data;
}



/**
 *  Print out the results of a certain algorithm based on epi variation
 */
void printout(std::vector<Result> results){
	
	std::vector<double> build_times;
	std::vector<double> search_times;
	std::vector<size_t> index_size;
	std::vector<size_t> heights;
	std::vector<std::vector<int>> segOfLayers;

	for(auto result : results){
		build_times.push_back(result.build_time);
		search_times.push_back(result.search_time);
		index_size.push_back(result.index_size);
		heights.push_back(result.height);
		segOfLayers.push_back(result.segOfLayer);
	}

	printf("\nThe build time is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		printf("%f,",build_times[i]);
		if((i+1)%12 == 0) printf("\n");
	}
	printf("\nThe search time is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		if(i%12 == 0) printf("[");
		std::cout<<search_times[i];
		if((i+1)%12 == 0) printf("],\n");
		else printf(",");
	}

	printf("\nThe index size is :");
	for(size_t i = 0 ; i< build_times.size();i++){
		
		printf("%ld,",index_size[i]);
		if((i+1)%12 == 0) printf("\n");
	}

	printf("\nThe height is :");
	for(size_t i = 0 ; i < build_times.size();i++){
		printf("%ld,",heights[i]);
		if((i+1)%12 == 0) printf("\n");
	}
	printf("\nThe #seg of each layers is :");
	for(size_t i = 0 ; i<  build_times.size();i++){
		printf("{");
		for(size_t j = 0 ; j< segOfLayers[i].size();j++){
			printf("%d,",segOfLayers[i][j]);
		}
		printf("} ");
		if((i+1)%12 == 0) printf("\n");
	}

	printf("\n");
}


template<typename K,typename T>
void testForModel(std::vector<K> data) {
	
	
    Result result;	
	std::vector<Result> results;

	int begin = 2;
	int ed = 13;
	for(int i = begin; i< ed;i++){
		for(int j = begin; j< ed;j++){
			// Begin an evaluation with different eps and eps_rec
			size_t epsilon = 1<<i;
			size_t epsilon_rec = 1<<j;

			printf("The epsilon is:%ld epsilon_i is %ld\n",epsilon,epsilon_rec);
			std::vector<std::vector<long>> residuals;
			// Get the time of constructing the whole PGM-Index
			auto start = std::chrono::high_resolution_clock::now();
			T model(data.begin(),data.end(),epsilon , epsilon_rec);
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> duration = end - start;
			double time = duration.count();
			result.build_time = time;
			result.height = model.height();
			result.segOfLayer = model.get_levels_offsets();
	
			// Randomly pick 1000 elements from the data and search for them
			std::vector<double> times;
			for(int k = 0; k< 10 ; k++){
				
				auto start = std::chrono::high_resolution_clock::now(); 	
				for(int i = 0;i<1000;i++){
					int index = rand() % data.size();
					
					auto app = model.search(data[index]);
					residuals.push_back(app.residuals);
					auto it = std::prev(std::upper_bound(data.begin()+app.lo,data.begin()+app.hi,data[index]));
				
					if(i>1000){
						printf("%ld",*it);
					}
				}
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> duration = end - start;
				double time = duration.count();
				time = time/1000.0;
				times.push_back(time);
			}
			
			std::sort(times.begin(),times.end());
			double avg_time_sum = 0;
			for(int i  = 1; i < 9 ; i++){
				avg_time_sum += times[i];
			}
			
			result.search_time = avg_time_sum / 8.0;
			result.index_size = model.size_in_bytes();
			
			results.push_back(result);
			for(size_t i = 0;i<residuals[0].size();i++){
				long res = 0;
				for(size_t j = 0;j<residuals.size();j++){
					res +=residuals[j][i];
				}
				double res_d = res/10000.0;
				printf("The average residual for layer %ld is %f\n",i,res_d);
			}	
			printf("\n");
		}
	}
	printout(results);

}


template<typename K>
void experiment_Optimal(std::vector<K> data) {
	
	printf("The results for OptimalPLA:\n");
	using PGM_Index1 = Optimal::PGMIndex<K>;
	testForModel<K,PGM_Index1>(data);

}

template<typename K>
void experiment_Greedy(std::vector<K> data) {

	
	using PGM_Index = Greedy::PGMIndex<K>;
	
	printf("The results for GreedyPLA:\n");
	testForModel<K,PGM_Index>(data);
}

template<typename K>
void experiment_Swing(std::vector<K> data) {
	
		

	using PGM_Index = Swing::PGMIndex<K>;
	printf("The results for SwingFilter:\n");
	testForModel<K,PGM_Index>(data);
}



int main() {
	// Load the data
	typedef uint64_t K;
	std::vector<K> data;
	data = load_data<K>();

	data.pop_back();
	std::cout<<data[data.size()-1]<<std::endl;
	printf("The size of the data is %ld\n",data.size());
	
	experiment_Swing(data);
	experiment_Greedy(data);
	experiment_Optimal(data);

	return 0;
}
