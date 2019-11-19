#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {
    cout << '\t' << "Assigning with Lloyd's Assignment" << endl;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
        clusters[i] = new vector<int>;
    int centroid;
    double min_dist, curr_dist;
    cout << num_of_centroids << endl << data_size << endl << dimension <<endl;
    for (int i = 0; i < data_size; i++) {
        if (find(centroids->begin(), centroids->end(), i) != centroids->end()) continue;
        centroid = -1;
        min_dist = DBL_MAX;
        for (int j = 0; j < num_of_centroids; j++) {
            curr_dist = dist(&(*dataset)[(*centroids)[j]], &(*dataset)[i], dimension);
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
vector<int>** Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {

}

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<int>;
template class Lloyd_assignment<double*>;
template class Inverse_assignment<int>;
template class Inverse_assignment<double*>;