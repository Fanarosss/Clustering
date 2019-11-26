#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "LSH_DataTypes.h"


using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<pair<vector<Point>*, int>> centroids) {
    cout << '\t' << "Assigning with Lloyd's Assignment" << endl;
    int num_of_centroids = centroids.size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    /* vector containing only the centroid id */
    vector<int> centroid_ids;
    for (auto it : centroids) {
        centroid_ids.push_back(it.second);
    }
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
        clusters[i] = new vector<int>;
    int centroid;
    double min_dist, curr_dist;
    for (int i = 0; i < data_size; i++) {
        if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
        centroid = -1;
        min_dist = DBL_MAX;
        for (int j = 0; j < num_of_centroids; j++) {
            curr_dist = dist(centroids[j].first, &(*dataset)[i], dimension);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                centroid = j;
            }
        }
        clusters[centroid]->push_back(i);
    }
    return clusters;
}

template <class Point>
string Lloyd_assignment<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<int>** Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<pair<vector<Point>*, int>> centroids) {
    cout << '\t' << "Assigning with Inverse Assignment" << endl;
    int num_of_centroids = centroids.size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    /* vector containing only the centroid id */
    vector<int> centroid_ids;
    for (auto it : centroids) {
        centroid_ids.push_back(it.second);
    }

   vector<vector<Point>> lsh_searchset;                 //centroids
//    vector<vector<Point>> lsh_searchset;               //queries
    for (int i = 0; i < num_of_centroids; i++){
        lsh_searchset.push_back((*centroids[i].first));
    }
//    for (int j = 0; j < data_size; j++){
//        if (find(centroid_ids.begin(), centroid_ids.end(), j) != centroid_ids.end()) continue;
//        lsh_searchset.push_back((*dataset)[j]);
//    }
//    cout << "num of centroids " << num_of_centroids << endl;
//    cout << "data size " << data_size <<  endl;
//    cout <<  "dataset size " << lsh_dataset.size() <<endl;
//    cout <<  "searchset size " << lsh_searchset.size() << endl;

    /* LSH Call for Vectors or Curves */
    lsh_datatype(dataset, &lsh_searchset, this->Grids, this->k, this->L, &centroid_ids, clusters);

    return clusters;
}

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<double>;
template class Lloyd_assignment<double*>;
template class Inverse_assignment<double>;
template class Inverse_assignment<double*>;