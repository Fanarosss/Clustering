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
    double sum, min, objective_function;
    for (int i = 0; i < this->get_K(); i++) {
        /* size */
        cluster_size = clusters[i]->size();
        /* distance for the current centroid */
        sum = 0.0;
        for (int j = 0; j < cluster_size; j++) {
            sum += db->get_distance((*clusters[i])[j], (*centroids)[i].second);
        }
        objective_function = sum;
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
        if ( objective_function == 0 || (fabs(objective_function - min) / objective_function) > RATE_OF_CHANGE_LOWERBOUND) convergence++;
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

    /* Choose Between Mean Vector for Vectors and DTW for Curves and return convergence value */
    int convergence = mv_dtw_datatype(dataset, clusters, centroids);

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

/* MV_DTW Update Operations */

/* Mean Vector */
int mv_dtw_datatype(vector<vector<double>>* dataset, vector<int>** clusters, vector<pair<vector<double>*, int>>* centroids){

    int convergence = 0;
    int cluster_size;
    int num_of_centroids = centroids->size();
    int dimension = (*dataset)[0].size() - 1;
    double sum, mean, global_mean;

    /* Loop for every centroid */
    for (int i = 0; i < num_of_centroids; i++) {
        global_mean = 0.0;

        /* Cluster size */
        cluster_size = clusters[i]->size();

        /* Mean Centroid */
        vector<double>* new_centroid = new vector<double>;

        /* Push back -1 as ID */
        new_centroid->push_back(-1);

        /* Calculate mean value for every dimension of vector */
        for (int j = 1; j <= dimension; j++) {
            sum = 0;
            mean = 0;
            for (int k = 0; k < cluster_size; k++) {
                sum += (*dataset)[(*clusters[i])[k]][j];
            }
            mean = sum / cluster_size;
            global_mean += mean;
            new_centroid->push_back(mean);
        }

        /* global mean for rate of centroid change */
        global_mean /= dimension;

        /* Calculate distance of previous and current centroid */
        if (dist((*centroids)[i].first, new_centroid, dimension) / global_mean > MEAN_PERCENTAGE_RATE_LOWERBOUND) convergence++;

        /* Update centroid pointer and ID */
        (*centroids)[i].first = new_centroid;
        (*centroids)[i].second = -1;
    }

    /* If all new centroids and previous centroids had distance < CONVERGENCE_DISTANCE, then convergence reached */
    return (convergence != 0) ? 0 : 1;
}

/* DTW centroid Curve */
int mv_dtw_datatype(vector<vector<double*>>* dataset, vector<int>** clusters, vector<pair<vector<double*>*, int>>* centroids){

    int convergence = 0;
    int cluster_size;
    int num_of_centroids = centroids->size();
    int lamda[num_of_centroids];
    double sum, mean, global_mean;
    double* point;
    vector<vector<double*>*> Initial_C;

    global_mean = 0;
    /* DBA Initialization - Loop to find lamda and Initial_C for every cluster  */
    for (int i = 0; i < num_of_centroids; i++) {

        /* Cluster size */
        cluster_size = clusters[i]->size();

        /* Lamda Calculation */
        sum = 0;
        mean = 0;
        for (int k = 0; k < cluster_size; k++) {
            sum += (*dataset)[(*clusters[i])[k]][0][1];
        }
        mean = sum / cluster_size;
        lamda[i] = mean;

        /* Initial_C Calculation */
        for (int k = 0; k < cluster_size; k++) {
            vector<double*>* curr_c = new vector<double*>;
            if((*dataset)[(*clusters[i])[k]][0][1] >= lamda[i]){
                point = new double [2];
                point[0] = -1;
                point[1] = lamda[i]+1;
                curr_c->push_back(point);
                for ( int l = 1; l <= lamda[i]; l++){
                    curr_c->push_back((*dataset)[(*clusters[i])[k]][l]);
                }
                Initial_C.push_back(curr_c);
                break;
            }
        }
    }

    /* DBA Algorithm - Loop to find Index-Pairs of Best-Traversal and calculate DTW centroid curve */
    for (int i = 0; i < num_of_centroids; i++) {
        global_mean = 0;
        /* DTW_pairs Calculation */
        vector<double*>* C = new vector<double*>;
        cluster_size = clusters[i]->size();
        vector<double*> A[lamda[i]];
        for (int k = 0; k < cluster_size; k++) {
            vector<pair<int,int>> pairs;
            DTW_pairs(Initial_C[i], &(*dataset)[(*clusters[i])[k]], &pairs);
            for (auto iter : pairs) {
                A[iter.first].push_back((*dataset)[(*clusters[i])[k]][iter.second + 1]);
            }
        }

        /* DTW centroid curve calculation */
        double sumx = 0, meanx = 0;
        double sumy = 0, meany = 0;
        point = new double [2];
        point[0] = -1;
        point[1] = lamda[i]+1;
        C->push_back(point);
        for (int l = 0; l < lamda[i]; l++){
            sumx = 0;
            sumy = 0;
            for(unsigned int m = 0; m < A[l].size(); m++){
                sumx += A[l][m][0];
                sumy += A[l][m][1];
            }
            meanx = sumx / A[l].size();
            meany = sumy / A[l].size();
            point = new double [2];
            point[0] = meanx;
            point[1] = meany;
            C->push_back(point);
            global_mean += (meanx + meany) / 2.0;
        }

        global_mean /= lamda[i];

        /* Calculate distance of previous and current centroid */
        if ((dist((*centroids)[i].first, C, C->size()) / global_mean) > MEAN_PERCENTAGE_RATE_LOWERBOUND) {
            convergence++;
        }

        /* Update centroid pointer and ID */
        (*centroids)[i].first = C;
        (*centroids)[i].second = -1;
    }

    /* If all new centroids and previous centroids had distance < CONVERGENCE_DISTANCE, then convergence reached */
    return (convergence != 0) ? 0 : 1;
} //todo : check convergence