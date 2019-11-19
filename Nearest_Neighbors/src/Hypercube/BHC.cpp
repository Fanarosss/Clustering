#include "LSH.h"
#include "LSH_Functions.h"
#include "BHC.h"
#include "BHC_Functions.h"
#include "Helper_Functions.h"
#include "HashTable.h"

using namespace std;

template void BHC <int>(vector<vector<int>>*, vector<vector<int>>*, int, int, int, int, int, double, vector<vector<int>>*, int**, double**, int**);
template void BHC <double>(vector<vector<double>>*, vector<vector<double>>*, int, int, int, int, double, double, vector<vector<int>>*, double**, double**, int**);

template <typename Point>
void BHC (vector<vector<Point>>* dataset, vector<vector<Point>>* searchset, int k, int dim, int M, int probes, Point w, double R, vector<vector<int>>* R_Neighbors, Point** min_distance, double** time, int** nearest_neighbor){

    /* Size of dataset and searchset */
    int d_size = dataset->size();
    int s_size = searchset->size();
    /* Size of Vectors (d-dimensional vectors) */
    int d = (*dataset)[0].size();
    /* Vector containing shifts of size(k,d) */
    vector<vector<double>> s;
    /* Vector containing H of size (k, d_size) */
    vector<vector<int>> hash_functions;
    /* Vector containing projections of data */
    vector<vector<int>> a_projects;
    /* internal vector H for pushing */
    vector<int> H;
    /* Amplified hash for dataset*/
    vector<vector<int>> data_amplified_g;
    /* Amplified hash for searchset */
    vector<vector<int>> query_amplified_g;
    /* Temporary g vector for push back */
    vector<int> temp_g;
    /* ANN results */
    vector<vector<vector<vector<Point>>>> ANN;
    /* Neighbors from probes */
    vector<vector<vector<Point>>> Neighbors;
    /* Map to assign every g function to 0 or 1 */
    map <int, int> dictionary;
    /* Vector for R-Neighbors (BONUS) */
    vector<int> Curr_R_Neighbors;

    /* Computing big numbers */
    int capital_M = pow(2, 32/k);
    int m = 3;                              //chosen after testing because of the results that it gave us
    int * power = new int [d-1];
    for (int j = 0; j < d-1; j++)
        power[j] = moduloPower(m, j, capital_M);

    /* Loop for creation of dim(d') amplified functions g */
    for (int l = 0; l < dim; l++) {

        /* Shifts Generation */
        generate_shifts(&s, w, d, k);

        /* ----------------------- DATA SET -------------------------------*/

        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, dataset, &(s[i]), w, d);
            /* Hash Computation */
            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            /* clear only deletes data, but does not free the underlying storage
            so we use swap */
            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        }

        /* Amplified hash computation */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* ------------------------ SEARCH SET ----------------------------*/

        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, searchset, &(s[i]), w, d);
            /* Hash Computation */
            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);

            vector<int>().swap(H);
            vector<vector<int>>().swap(a_projects);
        }

        /* Amplified hash computation */
        amplify_hash(&temp_g, &hash_functions, k);
        query_amplified_g.push_back(temp_g);

        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        vector<vector<double>>().swap(s);
    }

    /* clean memory */
    delete[] power;

    /* Mapping every amplified g function to 0 or 1 */
    fill_dictionary(&dictionary, data_amplified_g);

    int hypercube_size = pow(2,dim);
    int vertex = 0;

    /* Mapping every dataset point to a hypercube's vertex */
    vector<vector<Point>> MyVerticesTable[hypercube_size];
    for (int i = 0; i < dataset->size(); i++) {
        vertex = calculate_vertex(data_amplified_g, dictionary, i);
        MyVerticesTable[vertex].push_back((*dataset)[i]);
    }

    /* Finding for every searchset point its neighbors in their vertex and vertices of hamming distance = 1 */
    int ham_dist = 0;
    for (int i = 0; i < searchset->size(); i++) {
        vertex = calculate_vertex(query_amplified_g, dictionary, i);
        Neighbors.push_back(MyVerticesTable[vertex]);
        for (int j = 0; j < hypercube_size && Neighbors.size() < probes; j++) {
            ham_dist = hammingDistance(vertex, j);
            if(ham_dist == 1){
                Neighbors.push_back(MyVerticesTable[j]);
            }
        }
        ANN.push_back(Neighbors);
        vector<vector<vector<Point>>>().swap(Neighbors);
    }

    double distance;
    /* default metric L1 Manhattan */
    int Metric = 1;
    int computations = 0;
    /* For every query in searchset */
    for (int q = 0; q < searchset->size(); q++) {
        auto start = chrono::high_resolution_clock::now();
        /* For every probe */
        for (int n = 0; n < ANN[q].size(); n++) {
            computations = 0;
            /* For every vector in the same probe (max M calculations) */
            for (int j = 0; j < ANN[q][n].size(); j++) {
                if (computations == M && R == 0) break;
                distance = dist(&ANN[q][n][j], &searchset->at(q), dataset->at(0).size(), Metric);
                /* Find R-Neighbors */
                if(distance <= R){
                    Curr_R_Neighbors.push_back(ANN[q][n][j][0]);
                }
                /* Find Nearest Neighbor */
                if (((distance < (*min_distance)[q]) || (*min_distance)[q] == -1) && computations < M) {
                    (*min_distance)[q] = distance;
                    (*nearest_neighbor)[q] = ANN[q][n][j][0];
                }
                computations++;
            }
        }
        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();
        (*time)[q] = time_elapsed;

        sort( Curr_R_Neighbors.begin(), Curr_R_Neighbors.end() );
        Curr_R_Neighbors.erase( unique( Curr_R_Neighbors.begin(), Curr_R_Neighbors.end() ), Curr_R_Neighbors.end() );
        vector<int>::iterator position = find(Curr_R_Neighbors.begin(), Curr_R_Neighbors.end(), (*nearest_neighbor)[q]);
        if (position != Curr_R_Neighbors.end())
            Curr_R_Neighbors.erase(position);
        R_Neighbors->push_back(Curr_R_Neighbors);
        vector<int>().swap(Curr_R_Neighbors);
    }
}
