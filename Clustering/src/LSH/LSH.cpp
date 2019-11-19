#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
LSH<Point>::LSH(int k, int L, Point w) {
    this->k = k;
    this->L = L;
    this->w = w;
    this->power = NULL;
    MyHashTable = new HashTable <Point>* [L];
    for (int i = 0; i < L; i++)
        MyHashTable[i] = NULL;
}

template <class Point>
void LSH<Point>::fit(vector<vector<Point>>* dataset) {
    /* ----------------------- DATA SET FIT -------------------------------*/

    this->dataset = dataset;
    /* Size of dataset */
    int d_size = dataset->size();
    /* Size of Vectors (d-dimensional vectors) */
    this->d = (*dataset)[0].size();
    /* Vector containing H of size (k, d_size) */
    vector<vector<int>> hash_functions;
    /* Vector containing projections of data */
    vector<vector<int>> a_projects;
    /* Internal vector H for pushing */
    vector<int> H;
    /* Temporary g vector for push back */
    vector<int> temp_g;
    /* hash table size */
    this->TableSize = d_size/8;

    /* Computing big numbers */
    int M = pow(2, 32/k);
    int m = 3;                              //chosen after testing because of the results that it gave us
    power = new int [d-1];
    for (int j = 0; j < d-1; j++)
        power[j] = moduloPower(m, j, M);

    /* Shifts Generation */
    generate_shifts(&s, w, d, k, L);

    /* Loop for creation of L hash tables */
    for (int l = 0; l < L; l++) {
        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, dataset, &(s[l][i]), w, d);
            /* Hash Computation */
            compute_hash(&H, &a_projects, &power, d, k, w);
            hash_functions.push_back(H);
            /* clear only deletes data, but does not free the underlying storage
            so we use swap */
            vector<int>().swap(H);
            vector < vector < int >> ().swap(a_projects);
        }

        /* Amplified hash computation */
        amplify_hash(&temp_g, &hash_functions, k);
        data_amplified_g.push_back(temp_g);

        /* clean vectors */
        vector<vector<int>>().swap(hash_functions);
        vector<int>().swap(temp_g);

        /* Insert all items inside the Hash Table */
        MyHashTable[l] = new HashTable<Point>(TableSize);
        for (int i = 0; i < dataset->size(); i++) {
            MyHashTable[l]->Insert(data_amplified_g[l][i], (*dataset)[i]);
        }
    }
}

template <class Point>
void LSH<Point>::evaluate(vector<vector<Point>>* searchset, double R, vector<vector<int>>* R_Neighbors, Point** min_distance, double** time, int** nearest_neighbor) {
    /* ------------------------ QUERY SEARCH ----------------------------*/

    /* Vector containing H of size (k, d_size) */
    vector <vector<int>> hash_functions;
    /* Vector containing projections of data */
    vector <vector<int>> a_projects;
    /* Internal vector H for pushing */
    vector<int> H;
    /* Amplified hash for searchset */
    vector <vector<int>> query_amplified_g;
    /* Temporary g vector for push back */
    vector<int> temp_g;
    /* ANN results */
    vector<vector<vector<vector<Point>>>> ANN;
    /* ANN for each L iteration */
    vector<vector<vector<Point>>> ANNi;
    /* Vector for R-Neighbors (BONUS) */
    vector<int> Curr_R_Neighbors;

    /* Loop for creation of L hash tables */
    for (int l = 0; l < L; l++) {
        /* Loop for creation of K hi to create amplified g*/
        for (int i = 0; i < k; i++) {
            /* Projections Computation */
            projections(&a_projects, searchset, &(s[l][i]), w, d);
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

        /* ANN Calculation */
        for (int i = 0; i < searchset->size(); i++) {
            ANNi.push_back(*MyHashTable[l]->Search_Neighbors(query_amplified_g[l][i]));
        }
        ANN.push_back(ANNi);

        /* clear hash functions and s for next iteration */
        vector<vector<vector<Point>>>().swap(ANNi);
    }

    double distance;
    /* default metric L1 Manhattan */
    int Metric = 1;
    int computations = 0;
    /* For every query in searchset */
    for (int q = 0; q < searchset->size(); q++) {
        auto start = chrono::high_resolution_clock::now();
        /* For every hash table L */
        for (int i = 0; i < ANN.size(); i++) {
            computations = 0;
            /* For every vector in the same bucket (max 4*L calculations - we use 25 as limit) */
            for (int j = 0; j < ANN[i][q].size(); j++) {
                if (computations == 25 && R == 0) break;
                /* Check points that have the same amplified g */
                if (query_amplified_g[i][q] == data_amplified_g[i][(int)ANN[i][q][j][0]]) {
                    distance = dist(&ANN[i][q][j], &searchset->at(q), dataset->at(0).size(), Metric);
                    /* Find R-Neighbors */
                    if(distance <= R){
                        Curr_R_Neighbors.push_back(ANN[i][q][j][0]);
                    }
                    /* Find Nearest Neighbor */
                    if (((distance < (*min_distance)[q]) || (*min_distance)[q] == -1) && computations < 25) {
                        (*min_distance)[q] = distance;
                        (*nearest_neighbor)[q] = ANN[i][q][j][0];
                    }
                    computations++;
                }
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

        /* clean underlying memory of the vector */
        vector<int>().swap(Curr_R_Neighbors);
    }
    /* clean memory */
    vector<vector<vector<vector<Point>>>>().swap(ANN);
}

template <class Point>
LSH<Point>::~LSH() {
    /* clean memory */
    if (power != NULL)  delete[] power;
    for (int l = 0; l < L; l++) {
        if (MyHashTable[l] != NULL) delete MyHashTable[l];
    }
    delete[] MyHashTable;
}

template class LSH<int>;
template class LSH<double>;