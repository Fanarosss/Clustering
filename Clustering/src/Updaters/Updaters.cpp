#include "Updaters.h"
#include "Library.h"
#include "Helper_Functions.h"

using namespace std;

template <class Point>
int PAM<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<pair<vector<Point>*, int>>* centroids, DistanceDatabase<Point>* db) {
    /* minimize Sum(dist(i,t)) over all objects t in cluster C */
    /* OPTIMIZATIONS: 1) compute cluster size only once
     *                2) Keep distances between points in a triangular array
     *                because of the symmetry dist(i,j) == dist(j,i) */
    int convergence = 0;
    int t, row, col, cluster_size;
    double sum, min;
    double **distances;
//    cout << '\t' << "Updating with PAM: " << endl;
    for (int i = 0; i < this->get_K(); i++) {
//        cout << '\t' << "Cluster<" << i << "> " << endl;
        min = DBL_MAX;
        /* size */
        cluster_size = clusters[i]->size();
        /* initialize array for distances */
        distances = new double* [cluster_size];
        for (int j = 0; j < cluster_size; j++) {
            /* triangular storage to avoid duplicate distances ->
             * distance from 1 to 10 is the same as from 10 to 1 */
            distances[j] = new double[cluster_size];
            for (int l = 0; l < cluster_size; l++) {
                distances[j][l] = -1;
            }
        }
        /* iterate over cluster data */
        for (int j = 0; j < cluster_size; j++) {
            /* find medoid t to minimize the distances in this cluster */
            sum = 0.0;
            for (int l = 0; l < cluster_size; l++) {
                /* find indexes to minimize computations
                 * dist from 985 to 4751 is the same as 4751 to 985 */
                row = (j > l) ? j : l;
                col = (j < l) ? j : l;
                /* break point */
                if (j == l) continue;
                /* sum */
                if (distances[row][col] == -1)
                    distances[row][col] = db->get_distance((*clusters[i])[j], (*clusters[i])[l]);
                sum += distances[row][col];
            }
            /* find min and the id of min, make it centroid for this cluster */
            if (sum < min) {
                min = sum;
                t = (*clusters[i])[j];
            }
        }
        /* new centroid for this cluster */
        if (db->get_distance((*centroids)[i].second, t) > CONVERGENCE_DISTANCE) convergence++;
        (*centroids)[i].first = &(*dataset)[t];
        (*centroids)[i].second = t;
        /* clear memory */
        for (int j = 0; j < cluster_size; j++) {
            delete[] distances[j];
        }
        delete[] distances;
    }
//    cout << endl;
    return (convergence != 0) ? 0 : 1;
}

template <class Point>
string PAM<Point>::get_name() {
    return this->name;
}

template <class Point>
int MV_DTW<Point>::update(vector<vector<Point>>* dataset, vector<int>** clusters, vector<pair<vector<Point>*, int>>* centroids, DistanceDatabase<Point>* db) {

//    cout << '\t' << "Updating with MV_DTW: " << endl;

    int convergence = mv_dtw_datatype(dataset, clusters, centroids);
//    cout << "Convergence : " << convergence << endl;
    return convergence;
}

template <class Point>
string MV_DTW<Point>::get_name() {
    return this->name;
}


template class PAM<double>;
template class PAM<double*>;
template class MV_DTW<double>;
template class MV_DTW<double*>;