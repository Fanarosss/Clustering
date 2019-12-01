#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "LSH_DataTypes.h"

#define MAX_ITERATIONS 1000

using namespace std;

template void compute_unassigned<double>(vector<vector<double>>* ,vector<vector<double>>* , int, double**, int**);
template void compute_unassigned<double*>(vector<vector<double*>>* ,vector<vector<double*>>* , int, double**, int**);

void lsh_datatype(vector<vector<double>>* lsh_dataset, vector<vector<double>>* lsh_searchset, int Grids, int k, int L, vector<int>* centroid_ids, vector<int>** clusters){

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
    double w = 4*compute_window(lsh_dataset);
    LSH <double>* model = new LSH <double> (k, L, w);
    model->fit(lsh_dataset);

    int iterations = 0;
    while((unassigned_vectors > k) && (iterations < MAX_ITERATIONS)){
        model->evaluate_clusters(lsh_searchset, &min_distance, &nearest_centroid, &unassigned_vectors);
        iterations++;
    }
    if(unassigned_vectors > 0){
        compute_unassigned(lsh_dataset, lsh_searchset, data_size, &min_distance, &nearest_centroid);
    }

    delete (model);

    for (int i = 0; i < data_size; i++){
        if (min_distance[i] == 0) continue;
        clusters[nearest_centroid[i]]->push_back(i);
    }

    /* Clean Pointers */
    delete[] min_distance;
    delete[] nearest_centroid;

}

void lsh_datatype(vector<vector<double*>>* lsh_dataset, vector<vector<double*>>* lsh_searchset, int Grids, int k, int L, vector<int>* centroid_ids, vector<int>** clusters){

    double delta = 0.00006;
    int d = 2;                                                      /* default 2D curves */
    vector<vector<double>> data_vectored_curves;
    vector<vector<double>> search_vectored_curves;

    int data_size = lsh_dataset->size();
    /* Arrays for results */
    double *min_distance = new double[data_size];
    int *nearest_centroid = new int[data_size];
    int unassigned_curves = data_size;

    /* Initialize arrays */
    for (int i = 0; i < data_size; i++) {
        min_distance[i] = DBL_MAX;
        nearest_centroid[i] = -1;
    }

    for (int i = 0; i < Grids; i++) {

        /* Vectorization */
        Grid_Vectorization(delta, d, lsh_dataset, lsh_searchset, &data_vectored_curves, &search_vectored_curves);

        /* ---- LSH model ---- */
        double w = 4 * compute_window(&data_vectored_curves);
        LSH<double> *model = new LSH<double>(k, L, w);
        model->fit(&data_vectored_curves);

        int iterations = 0;
        while((unassigned_curves > k) && (iterations < MAX_ITERATIONS)){
            model->evaluate_clusters(&search_vectored_curves, &min_distance, &nearest_centroid, &unassigned_curves);
            iterations++;
        }

        delete (model);

        /* Clean Vectors */
        vector<vector<double>>().swap(data_vectored_curves);
        vector<vector<double>>().swap(search_vectored_curves);
    }
    if(unassigned_curves > 0){
        compute_unassigned(lsh_dataset, lsh_searchset, data_size, &min_distance, &nearest_centroid);
    }

    for (int i = 0; i < data_size; i++){
        if (min_distance[i] == 0) continue;
        clusters[nearest_centroid[i]]->push_back(i);
    }

    /* Clean Pointers */
    delete[] min_distance;
    delete[] nearest_centroid;

}

template <typename Point>
void compute_unassigned(vector<vector<Point>>* lsh_dataset,vector<vector<Point>>* lsh_searchset, int data_size, double** min_distance, int** nearest_centroid){
    double curr_dist;
    /* default metric L1 Manhattan */
    int Metric = 1;
    for(int i = 0; i < data_size; i++){
        if((*nearest_centroid)[i] == -1){
            for(int j = 0; j < lsh_searchset->size(); j++){
                curr_dist = dist(&lsh_dataset->at(i) ,&lsh_searchset->at(j), lsh_dataset->at(0).size(), Metric);
                //todo: use already computed distances
                if (curr_dist < (*min_distance)[i]) {
                    (*min_distance)[i] = curr_dist;
                    (*nearest_centroid)[i] = j;
                }
            }
        }
    }
}