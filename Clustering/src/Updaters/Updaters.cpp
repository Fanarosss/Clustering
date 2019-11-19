#include "Updaters.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
void PAM<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<int>* centroids) {
    /* minimize Sum(dist(i,t)) over all objects t in cluster C */
    /* OPTIMIZATIONS: 1) compute cluster size only once
     *                2) Keep distances between points in a triangular array
     *                because of the symmetricity dist(i,j) == dist(j,i) */
    int t, row, col, cluster_size;
    double sum, min;
    double **distances;
    for (int i = 0; i < this->get_K(); i++) {
        cout << "PAM for cluster <" << i << "> " << endl;
        min = DBL_MAX;
        /* size */
        cluster_size = clusters[i]->size();
        /* initialize array for distances */
        distances = new double* [cluster_size];
        for (int j = 0; j < cluster_size; j++) {
            /* triangular storage to avoid duplicate distances ->
             * distance from 1 to 10 is the same as from 10 to 1 */
            distances[j] = new double[cluster_size - j];
            for (int l = 0; l < (cluster_size - j); l++) {
                distances[j][l] = -1;
            }
        }
        /* iterate over cluster data */
        for (int j = 0; j < cluster_size; j++) {
            /* find medoid t to minimize the distances in this cluster */
            sum = 0.0;
            for (int l = 0; l < cluster_size; l++) {
                /* find indexes */
                row = (j > l) ? j : l;
                col = (j < l) ? j : l;
                /* break point */
                if (j == l) continue;
                /* sum */
                if (distances[row][col] == -1)
                    distances[row][col] = dist(&(*dataset)[(*clusters[i])[j]], &(*dataset)[(*clusters[i])[l]]);
                sum += distances[row][col];
            }
            /* find min and the id of min, make it centroid for this cluster */
            if (sum < min) {
                min = sum;
                t = (*clusters[i])[j];
            }
        }
        /* new centroid for this cluster */
        centroids->at(i) = t;
        /* clear memory */
        for (int j = 0; j < cluster_size; j++) {
            delete[] distances[j];
        }
        delete[] distances;
    }
    return;
}

template <class Point>
string PAM<Point>::get_name() {
    return this->name;
}

template <class Point>
void MV_DTW<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<int>* centroids) {
    return;
}

template <class Point>
string MV_DTW<Point>::get_name() {
    return this->name;
}

template class PAM<int>;
template class PAM<double*>;
template class MV_DTW<int>;
template class MV_DTW<double*>;