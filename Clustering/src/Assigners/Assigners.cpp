#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"
#include "LSH.h"
#include "LSH_Functions.h"
#include "LSH_DataTypes.h"


using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<pair<vector<Point>*, int>>* centroids, DistanceDatabase<Point>* db) {
    cout << '\t' << "Assigning with Lloyd's Assignment" << endl;
    int all_clusters_non_empty = 1;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size() - 1;
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
    clusters[i] = new vector<int>;
    /* vector containing only the centroid id */
    vector<int> centroid_ids;
    vector<int> centroid_pool[num_of_centroids];

    while(all_clusters_non_empty == 1){
        for(int i = 0 ; i < num_of_centroids ; i++ )
            clusters[i]->clear();
        all_clusters_non_empty = 0;
        /* vector containing only the centroid id */
        centroid_ids.clear();
        for (auto it : (*centroids)) {
            centroid_ids.push_back(it.second);
        //            cout << "CENTROID  : " << it.second << endl;
        }
        int centroid;
        double min_dist, max_dist,  curr_dist;
        for (int i = 0; i < data_size; i++) {
            if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
            centroid = -1;
            min_dist = DBL_MAX;
            //            cout << "DATA " << i << " : " << endl;
            for (int j = 0; j < num_of_centroids; j++) {
//                curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);                //todo: already calculated
                if((*centroids)[j].second == -1){
                    curr_dist = curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                }else{
                    curr_dist = db->get_distance((*centroids)[j].second, i);
                }
                //                cout << curr_dist << endl;
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    centroid = j;
                }
            }
            //            cout << "Centroid : " << centroid << endl;
            clusters[centroid]->push_back(i);
        }
        for (int j = 0; j < num_of_centroids; j++){
            cout << clusters[j]->size() << endl;
        }
        for (int j = 0; j < num_of_centroids; j++){
            if(clusters[j]->size() == 0){
                //                cout << "HERE : " << j << endl;
                centroid_pool[j].push_back((*centroids)[j].second);
                all_clusters_non_empty = 1;
                centroid = -1;
                max_dist = 0;
                for (int i = 0; i < data_size; i++) {
                    if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
                    if (find(centroid_pool[j].begin(), centroid_pool[j].end(), i) != centroid_pool[j].end()) continue;
//                    curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);           //todo: already calculated
                    if((*centroids)[j].second == -1){
                        curr_dist = curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                    }else{
                        curr_dist = db->get_distance((*centroids)[j].second, i);
                    }
                    if (curr_dist > max_dist) {
                        max_dist = curr_dist;
                        centroid = i;
                    }
                }
                cout << "OLD : " << (*centroids)[j].second << endl;
                cout << "NEW : " << centroid << endl;
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
    cout << '\t' << "Assigning with Inverse Assignment" << endl;
    int all_clusters_non_empty = 1;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size() - 1;
//
    for( int i = 0; i < num_of_centroids; i++){
        cout << "Centroid : " << (*centroids)[i].first << " , ID: " << (*centroids)[i].second << endl;
    }

    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
        clusters[i] = new vector<int>;
    vector<int> centroid_ids;
    vector<int> centroid_pool[num_of_centroids];
    vector<vector<Point>> lsh_searchset;                 //centroids

    while(all_clusters_non_empty == 1){
        for(int i = 0 ; i < num_of_centroids ; i++ ){
            clusters[i]->clear();
        }
        all_clusters_non_empty = 0;
        /* vector containing only the centroid id */
        centroid_ids.clear();
        for (auto it : (*centroids)) {
            centroid_ids.push_back(it.second);
//            cout << "CENTROID  : " << it.second << endl;
        }
        lsh_searchset.clear();
        for (int i = 0; i < num_of_centroids; i++){
            lsh_searchset.push_back((*(*centroids)[i].first));
        }


        /* LSH Call for Vectors or Curves */
        lsh_datatype(dataset, &lsh_searchset, &centroid_ids, this->Grids, this->k, this->L, clusters);

        for (int j = 0; j < num_of_centroids; j++){
            cout << clusters[j]->size() << endl;
        }
        int centroid;
        double max_dist, curr_dist;
        for (int j = 0; j < num_of_centroids; j++){
            if(clusters[j]->size() == 0){
//                cout << "HERE : " << j << endl;
                centroid_pool[j].push_back((*centroids)[j].second);
                all_clusters_non_empty = 1;
                centroid = -1;
                max_dist = 0;
                for (int i = 0; i < data_size; i++) {
                    if (find(centroid_ids.begin(), centroid_ids.end(), i) != centroid_ids.end()) continue;
                    if (find(centroid_pool[j].begin(), centroid_pool[j].end(), i) != centroid_pool[j].end()) continue;
//                    curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                    if((*centroids)[j].second == -1){
                        curr_dist = curr_dist = dist((*centroids)[j].first, &(*dataset)[i], dimension);
                    }else{
                        curr_dist = db->get_distance((*centroids)[j].second, i);
                    }
                    if (curr_dist > max_dist) {
                        max_dist = curr_dist;
                        centroid = i;
                    }
                }
                cout << "OLD : " << (*centroids)[j].second << endl;
                cout << "NEW : " << centroid << endl;
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