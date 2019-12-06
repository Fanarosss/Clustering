#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "LSH_DataTypes.h"


using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<pair<vector<Point>*, int>>* centroids, DistanceDatabase<Point>* db) {

    int all_clusters_non_empty = 1;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size() - 1;

    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ ) clusters[i] = new vector<int>;

    /* vector with ids of all centroids */
    vector<int> centroid_ids;
    /* vector with all ids a certain centroid has changed in this current iteration to avoid being reassigned to the same point over and over again */
    vector<int> centroid_pool[num_of_centroids];

    /* Loop to secure that no cluster is left with size 0 */
    while(all_clusters_non_empty == 1){
        all_clusters_non_empty = 0;

        /* Clear vectors of previous loops */
        for(int i = 0 ; i < num_of_centroids ; i++ ) {
            clusters[i]->clear();
        }
        centroid_ids.clear();
        for (auto it : (*centroids)) {
            centroid_ids.push_back(it.second);
        }

        int centroid;
        double min_dist, max_dist,  curr_dist;
        /* Loop to assign all queries to closest centroid */
        for (int i = 0; i < data_size; i++) {
            if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
            centroid = -1;
            min_dist = DBL_MAX;
            for (int j = 0; j < num_of_centroids; j++) {
                if((*centroids)[j].second == -1){
                    curr_dist = curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                }else{
                    curr_dist = db->get_distance((*centroids)[j].second, i);
                }
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    centroid = j;
                }
            }
            clusters[centroid]->push_back(i);
        }

        /* Loop to find if there is any empty cluster and change its centroid to the most distant point */
        for (int j = 0; j < num_of_centroids; j++){
            if(clusters[j]->size() == 0){
                centroid_pool[j].push_back((*centroids)[j].second);
                all_clusters_non_empty = 1;
                centroid = -1;
                max_dist = 0;
                for (int i = 0; i < data_size; i++) {
                    /* If i is already a centroid of different clusters, then continue */
                    if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
                    /* If i is already calculated as centroid for this cluster in previous iteration, then continue */
                    if (find(centroid_pool[j].begin(), centroid_pool[j].end(), i) != centroid_pool[j].end()) continue;
                    if((*centroids)[j].second == -1){
                        /* Case of Medoid   -> Calculate Distance */
                        curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                    }else{
                        /* Case of Centroid -> Find Already Calculated Distance in Database */
                        curr_dist = db->get_distance((*centroids)[j].second, i);
                    }
                    if (curr_dist > max_dist) {
                        max_dist = curr_dist;
                        centroid = i;
                    }
                }
                centroid_ids.erase(remove(centroid_ids.begin(), centroid_ids.end(), (*centroids)[j].second ), centroid_ids.end());
                (*centroids)[j].first = &(*dataset)[centroid];
                (*centroids)[j].second = centroid;
                centroid_ids.push_back(centroid);
            }
        }
    }
    return clusters;
}

template <class Point>
string Lloyd_assignment<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<int>** Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<pair<vector<Point>*, int>>* centroids, DistanceDatabase<Point>* db) {

    int all_clusters_non_empty = 1;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size() - 1;

    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ ) clusters[i] = new vector<int>;

    /* vector with ids of all centroids */
    vector<int> centroid_ids;
    /* vector with all ids a certain centroid has changed in this current iteration to avoid being reassigned to the same point over and over again */
    vector<int> centroid_pool[num_of_centroids];
    /* vector with centroids to be sent as searchset to LSH */
    vector<vector<Point>> lsh_searchset;

    /* Loop to secure that no cluster is left with size 0 */
    while(all_clusters_non_empty == 1){
        all_clusters_non_empty = 0;

        /* Clear vectors of previous loops */
        for(int i = 0 ; i < num_of_centroids ; i++ ){
            clusters[i]->clear();
        }
        centroid_ids.clear();
        for (auto it : (*centroids)) {
            centroid_ids.push_back(it.second);
        }
        lsh_searchset.clear();
        for (int i = 0; i < num_of_centroids; i++){
            lsh_searchset.push_back((*(*centroids)[i].first));
        }


        /* LSH Call for Vectors or Curves */
        lsh_datatype(dataset, &lsh_searchset, &centroid_ids, this->Grids, this->k, this->L, this->w, clusters, db);

        int centroid;
        double max_dist, curr_dist;
        /* Loop to find if there is any empty cluster and change its centroid to the most distant point */
        for (int j = 0; j < num_of_centroids; j++){
            if(clusters[j]->size() == 0){
                centroid_pool[j].push_back((*centroids)[j].second);
                all_clusters_non_empty = 1;
                centroid = -1;
                max_dist = 0;
                for (int i = 0; i < data_size; i++) {
                    /* If i is already a centroid of different clusters, then continue */
                    if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
                    /* If i is already calculated as centroid for this cluster in previous iteration, then continue */
                    if (find(centroid_pool[j].begin(), centroid_pool[j].end(), i) != centroid_pool[j].end()) continue;
                    if((*centroids)[j].second == -1){
                        /* Case of Medoid   -> Calculate Distance */
                        curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                    }else{
                        /* Case of Centroid -> Find Already Calculated Distance in Database */
                        curr_dist = db->get_distance((*centroids)[j].second, i);
                    }
                    if (curr_dist > max_dist) {
                        max_dist = curr_dist;
                        centroid = i;
                    }
                }
                centroid_ids.erase(remove(centroid_ids.begin(), centroid_ids.end(), (*centroids)[j].second ), centroid_ids.end());
                (*centroids)[j].first = &(*dataset)[centroid];
                (*centroids)[j].second = centroid;
                centroid_ids.push_back(centroid);
            }
        }
    }
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