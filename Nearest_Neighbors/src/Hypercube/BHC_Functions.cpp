#include "BHC_Functions.h"

using namespace std;

void show_bhc_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k <int>         Dimension of Hypercube\n"
              << "\t-M <int>         Max allowed points to check\n"
              << "\t-probes <int>    Max allowed vertices to check\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

void fill_dictionary (map<int,int>* dictionary,vector<vector<int>> amplified_g) {
    /* Map every amplified function g to 0 or 1
     * Keep them in a map as a dictionary so as not to calculate them every time */

    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,1);
    map<int,int>::iterator map_it;

    for (int i = 0; i < amplified_g.size(); i++) {
        for (int j=0; j < amplified_g[0].size(); j++) {
            map_it = dictionary->find(amplified_g[i][j]);
            if (map_it == dictionary->end()) {
                dictionary->insert(pair<int, int>(amplified_g[i][j], distribution(generator)));
            }
        }
    }
}

int calculate_vertex (vector<vector<int>> amplified_g, map<int,int> dictionary, int index) {
    /* Concatenate amplified g mappings to find the corresponding vertex of hupercube */

    int vertex = 0;
    map<int,int>::iterator map_it;

    for (int i=0; i < amplified_g.size(); i++){
        map_it = dictionary.find(amplified_g[i][index]);
        if (map_it != dictionary.end()) {
            vertex = vertex << 1 | map_it->second;
        }
    }

    return vertex;
}

int hammingDistance(int n1, int n2){
    /* Calculate hamming distance of two vertices of hyoercube */

    int x = n1 ^ n2;
    int setBits = 0;

    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }

    return setBits;
}