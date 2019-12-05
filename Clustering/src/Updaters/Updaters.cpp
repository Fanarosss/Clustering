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
    for (int i = 0; i < this->get_K(); i++) {
        min = DBL_MAX;
        /* size */
        cluster_size = clusters[i]->size();
        /* distance for the current centroid */
        sum = 0.0;
        for (int j = 0; j < cluster_size; j++) {
            sum += db->get_distance((*clusters[i])[j], (*centroids)[i].second);
        }
        min = sum;
        t = (*centroids)[i].second;
        /* iterate over cluster data */
        for (int j = 0; j < cluster_size; j++) {
            /* find medoid t to minimize the distances in this cluster */
            sum = 0.0;
            for (int l = 0; l < cluster_size; l++) {
                /* break point */
                if (j == l) continue;
                /* sum */
                sum += db->get_distance((*clusters[i])[j], (*clusters[i])[l]);
            }
            /* increment sum because centroid belongs to cluster, even though its not in the vector */
            sum += db->get_distance((*clusters[i])[j], (*centroids)[i].second);
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
    }
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