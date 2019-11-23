#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "LSH_DataTypes.h"


using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<vector<Point>*> centroids) {
    cout << '\t' << "Assigning with Lloyd's Assignment" << endl;
    int num_of_centroids = centroids.size();
    cout << "NUMBER OF CENTROIDS : " << num_of_centroids << endl;
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
        clusters[i] = new vector<int>;
    int centroid;
    double min_dist, curr_dist;
    for (int i = 0; i < data_size; i++) {
        if (find(centroids.begin(), centroids.end(), &(*dataset)[i]) != centroids.end()) continue;
        centroid = -1;
        min_dist = DBL_MAX;
        for (int j = 0; j < num_of_centroids; j++) {
            curr_dist = dist(centroids[j], &(*dataset)[i], dimension);
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
vector<int>** Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<vector<Point>*> centroids) {
    cout << '\t' << "Assigning with Inverse Assignment" << endl;
    int num_of_centroids = centroids.size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];

    vector<vector<Point>> lsh_dataset;                 //centroids
    vector<vector<Point>> lsh_searchset;               //queries
    for (int i = 0; i < num_of_centroids; i++){
        lsh_dataset.push_back((*centroids[i]));
    }
    for (int j = 0; j < data_size; j++){
        if (find(centroids.begin(), centroids.end(), &(*dataset)[j]) != centroids.end()) continue;
        lsh_searchset.push_back((*dataset)[j]);
    }

    /* LSH Call for Vectors or Curves */
    lsh_datatype(&lsh_dataset, &lsh_searchset, this->Grids, this->k, this->L, centroids, clusters);

    return clusters;
}

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<int>;
template class Lloyd_assignment<double*>;
template class Inverse_assignment<int>;
template class Inverse_assignment<double*>;