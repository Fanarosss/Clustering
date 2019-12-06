#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "LSH_DataTypes.h"

using namespace std;

template void compute_unassigned<double>(vector<vector<double>>* ,vector<vector<double>>* , vector<int>*, int, double**, int**, DistanceDatabase<double>* db);
template void compute_unassigned<double*>(vector<vector<double*>>* ,vector<vector<double*>>* , vector<int>*, int, double**, int**, DistanceDatabase<double*>* db);

void lsh_datatype(vector<vector<double>>* lsh_dataset, vector<vector<double>>* lsh_searchset, vector<int>* centroid_ids, int Grids, int k, int L, double w, vector<int>** clusters, DistanceDatabase<double>* db){

    int data_size = lsh_dataset->size();

    /* Arrays for results */
    double *min_distance = new double[data_size];
    int *nearest_centroid = new int[data_size];
    int unassigned_vectors = data_size;

    /* Initialize arrays */
    for (int i = 0; i < data_size; i++) {
        min_distance[i] = DBL_MAX;
        nearest_centroid[i] = -1;
    }

    /* ---- LSH model ---- */
    LSH <double>* model = new LSH <double> (k, L, w);
    model->fit(lsh_dataset);

    int iterations = 0;
    /* LSH computations until all vectors are assigned to centroid or LSH_ITERATIONS reached */
    while((unassigned_vectors > k) && (iterations < LSH_ITERATIONS)){
        model->evaluate_clusters(lsh_searchset, centroid_ids, &min_distance, &nearest_centroid, &unassigned_vectors);
        iterations++;
    }

    /* If there are any unassigned vectors, assign them using brute force */
    if(unassigned_vectors > 0){
        compute_unassigned(lsh_dataset, lsh_searchset, centroid_ids, data_size, &min_distance, &nearest_centroid, db);
    }

    delete (model);

    /* Fill every cluster with its assigned vectors */
    for (int i = 0; i < data_size; i++){
        if (min_distance[i] == 0) continue;
        clusters[nearest_centroid[i]]->push_back(i);
    }

    /* Clean Pointers */
    delete[] min_distance;
    delete[] nearest_centroid;

}

void lsh_datatype(vector<vector<double*>>* lsh_dataset, vector<vector<double*>>* lsh_searchset, vector<int>* centroid_ids, int Grids, int k, int L, double w, vector<int>** clusters, DistanceDatabase<double*>* db){

    /* default 2D curves */
    int d = 2;
    double delta = 0.00006;

    /* Vectors to be used in Grid Vectorization */
    vector<vector<double>> data_vectored_curves;
    vector<vector<double>> search_vectored_curves;

    int data_size = lsh_dataset->size();
    int centroid_num = lsh_searchset->size();

    /* Arrays for results */
    double *min_distance = new double[centroid_num];
    double *min_distance_data = new double[data_size];
    int *nearest_neighbor = new int[centroid_num];
    int *nearest_centroid = new int[data_size];
    double *time = new double[centroid_num];
    int unassigned_curves = data_size;

    /* Initialize arrays */
    for (int i = 0; i < centroid_num; i++) {
        min_distance[i] = DBL_MAX;
        nearest_neighbor[i] = -1;
        time[i] = 0.0;
    }
    for (int i = 0; i < data_size; i++) {
        min_distance_data[i] = DBL_MAX;
        nearest_centroid[i] = -1;
    }

    /* Loop for all Grids */
    for (int i = 0; i < Grids; i++) {
        /* Vectorization */
        Grid_Vectorization(delta, d, lsh_dataset, lsh_searchset, &data_vectored_curves, &search_vectored_curves);

        /* ---- LSH model ---- */
        LSH<double>* model = new LSH<double>(k, L, w);
        model->fit(&data_vectored_curves);

        int iterations = 0;
        vector<vector<int>> R_Neighbors;
        int R = (int)(data_size / (centroid_num)) - 1;
        /* LSH computations until all vectors are assigned to centroid or LSH_ITERATIONS reached */
        model->evaluate(&search_vectored_curves, R, &R_Neighbors, &min_distance, &time, &nearest_neighbor);

        for (int j = 0; j < centroid_num; j++) {
            for (auto centroid_neighbors : R_Neighbors) {
                for (auto neighbor : centroid_neighbors) {
                    if (nearest_centroid[neighbor] == -1) {
                        nearest_centroid[neighbor] = j;
                        if((*centroid_ids)[j] == -1){
                            /* Case of Centroid   -> Calculate Distance */
                            min_distance_data[neighbor] = dist(&lsh_dataset->at(neighbor), &lsh_searchset->at(j), lsh_dataset->at(neighbor).size());
                        }else{
                            /* Case of Medoid -> Find Already Calculated Distance in Database */
                            min_distance_data[neighbor] = db->get_distance((*centroid_ids)[j], neighbor);
                        }

                    } else {
                        /* measure distances and take the smaller one */
                        double distance;
                        if((*centroid_ids)[j] == -1){
                            /* Case of Centroid   -> Calculate Distance */
                            distance = dist(&lsh_dataset->at(neighbor), &lsh_searchset->at(j), lsh_dataset->at(neighbor).size());
                        }else{
                            /* Case of Medoid -> Find Already Calculated Distance in Database */
                            distance = db->get_distance((*centroid_ids)[j], neighbor);
                        }
                        if (distance < min_distance_data[neighbor]) {
                            nearest_centroid[neighbor] = j;
                            min_distance_data[neighbor] = distance;
                        }
                    }
                }
            }
        }
        delete (model);

        /* Clean Vectors */
        vector<vector<double>>().swap(data_vectored_curves);
        vector<vector<double>>().swap(search_vectored_curves);
    }

    /* If there are any unassigned curves, assign them using brute force */
    compute_unassigned(lsh_dataset, lsh_searchset, centroid_ids, data_size, &min_distance_data, &nearest_centroid, db);

    /* Fill every cluster with its assigned curves */
    for (int i = 0; i < data_size; i++){
        if (find(centroid_ids->begin(), centroid_ids->end(), nearest_centroid[i]) != centroid_ids->end()) continue;
        clusters[nearest_centroid[i]]->push_back(i);
    }

    /* Clean Pointers */
    delete[] min_distance;
    delete[] nearest_neighbor;
    delete[] nearest_centroid;
    delete[] time;

}


/* Function to assign unassigned elements using Brute Force */

template <typename Point>
void compute_unassigned(vector<vector<Point>>* lsh_dataset,vector<vector<Point>>* lsh_searchset, vector<int>* centroid_ids, int data_size, double** min_distance, int** nearest_centroid, DistanceDatabase<Point>* db){

    /* default metric L1 Manhattan */
    int Metric = 1;
    double curr_dist;
    /* For every element in dataset */
    for(int i = 0; i < data_size; i++){
        /* If element is not assigned to centroid */
        if((*nearest_centroid)[i] == -1){
            /* For every centroid calculate distance to element and assign to closest */
            for(int j = 0; j < lsh_searchset->size(); j++){
                if((*centroid_ids)[j] == -1){
                    /* Case of Centroid   -> Calculate Distance */
                    curr_dist = dist(&lsh_dataset->at(i) ,&lsh_searchset->at(j), lsh_dataset->at(0).size(), Metric);
                }else{
                    /* Case of Medoid -> Find Already Calculated Distance in Database */
                    curr_dist = db->get_distance((*centroid_ids)[j], i);
                }
                if (curr_dist < (*min_distance)[i]) {
                    (*min_distance)[i] = curr_dist;
                    (*nearest_centroid)[i] = j;
                }
            }
        }
    }
}