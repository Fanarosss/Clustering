#include "LSH_Functions.h"
#include "Helper_Functions.h"

using namespace std;

template double compute_window<int>(vector<vector<int>>*);
template double compute_window<double>(vector<vector<double>>*);
template void projections<int>(vector<vector<int>>*, vector<vector<int>>*, vector<double>*, double, int);
template void projections<double>(vector<vector<int>>*, vector<vector<double>>*, vector<double>*, double, int);

void show_lsh_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-d <input file>  Path to data file\n"
              << "\t-q <query file>  Path to search file\n"
              << "\t-k <int>         Number of hi function for construction of g function\n"
              << "\t-L <int>         Number of Hash tables\n"
              << "\t-o <output file> Path to file of results\n"
              << endl;
}

template <typename Point>
double compute_window(vector<vector<Point>>* dataset) {
    /* 1. Take all points in dataset
     * 2. Find their nearest neighbor using L1 metric
     * 3. Average all distances between points
     * L1 = sum(|P1i - P2i|) for i in 0,d-1
     * Note: dist is a function scalable for all metrics(L1,L2,L3 etc.) */

    vector<Point> distances;
    Point L1, min_distance;
    int size = dataset->size();
    int d = dataset->at(0).size();
    vector<Point> P1;
    vector<Point> P2;

    for (int i = 0; i < size; i++) {
        P1 = dataset->at(i);
        L1 = 0;
        min_distance = -1;
        for (int j = 0; j < size; j++) {
            if (i != j) {
                P2 = dataset->at(j);
                /* Default is L1 metric, for Lk metric, add a 4th argument*/
                L1 = dist(&P1, &P2, d);
                if (min_distance == -1)
                    min_distance = L1;
                if (L1 < min_distance)
                    min_distance = L1;
            }
            L1 = 0;
        }
        distances.push_back(min_distance);
    }
    double w = accumulate(distances.begin(), distances.end(), 0) / size;

    return w;
}

void generate_shifts(vector<vector<vector<double>>>* s, double w, int d, int k, int L){
    /* Generate K * Si for every dimension
     * At the end, s will be a vector of size (k,d) */

    unsigned seed;
    uniform_real_distribution<double> distribution (0, w);
    vector<double> Sj;
    vector<vector<double>> Sl;

    for (int l = 0; l < L; l++) {
        for (int i = 0; i < k; i++) {
            seed = chrono::system_clock::now().time_since_epoch().count();
            default_random_engine generator(seed);
            for (int j = 1; j < d; j++) {
                Sj.push_back(distribution(generator));
            }
            Sl.push_back(Sj);
            vector<double>().swap(Sj);
        }
        s->push_back(Sl);
        vector<vector<double>>().swap(Sl);
    }
}

template <typename Point>
void projections(vector<vector<int>>* a_projects, vector<vector<Point>>* x, vector<double>* s, double w, int d) {
    /* Ai = (Xi - Si) / W
     * Project every X to A in d-dimensional grid shifted by S, where every cell size = W */

    int ai;
    vector<int> a;

    for (int i = 0; i < x->size(); i++) {
        for (int dim = 1; dim < d; dim++) {
            ai = floor((double)((*x)[i][dim] - (*s)[dim - 1]) / w); // used to be + w
            a.push_back(ai);
        }
        a_projects->push_back(a);
        a.clear();
        a.shrink_to_fit();
    }
}

void compute_hash(vector<int>* H, vector<vector<int>> *a, int** power, int d, int k, double w){
    /* Compute K of hash functions hi for every point
     * Vector H at the end will have a size of (dataset.size(), k)
     * Usage of moduloMultiplication and modulo functions was necessary to avoid overflow */

    int M = 0, h = 0, term = 0, dim = d - 1;
    M = pow(2, 32/k);

    for (int i = 0; i < a->size(); i++){
        h=0;
        term = 0;
        for (int j = dim - 1; j >= 0; j--) {
            term += moduloMultiplication((*a)[i][j], (*power)[(dim - 1) - j], M);
        }
        h = modulo(term, M);
        H->push_back(h);
    }
}

void amplify_hash(vector<int>* amplified_g, vector<vector<int>>* hash_functions, int k){
    /* For every item it amplifies the hash from K dimensions to 1
     * g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)] */

    int g;
    int concat_dist = 31/k;

    /* for all points in dataset */
    for (int i = 0; i < (*hash_functions)[0].size(); i++) {
        g=0;
        for (int j = 0; j < k; j++) {
            g += (g << concat_dist) | ((*hash_functions)[j][i]);
        }
        amplified_g->push_back(g);
    }
}
