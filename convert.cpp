#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>


/**
 * @brief Translate the data from the file to the vector
 *
 * @param filename
 * @return std::vector<K>
 */
typedef int K;
std::vector<K> load_data(const std::string &filename)
{
    std::cout << "dataset:" << filename.substr(8) << std::endl;

    /* Open file. */
    std::ifstream in(filename, std::ios::binary);
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
    std::sort(data.begin(), data.end());

    return data;
}