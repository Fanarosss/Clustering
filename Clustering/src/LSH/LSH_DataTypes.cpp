#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "LSH_DataTypes.h"

using namespace std;

void lsh_datatype(vector<vector<int>>* lsh_dataset, vector<vector<int>>* lsh_searchset, int Grids, int k, int L, vector<vector<int>*> centroids, vector<int>** clusters){

    double R;
    vector<vector<int>> R_neighbors;

    /* Arrays for results */
    int *min_distance = new int[lsh_searchset->size()];
    int *nearest_centroid = new int[lsh_searchset->size()];
    double *time = new double[lsh_searchset->size()];

    /* Initialize arrays */
    for (int i = 0; i < lsh_searchset->size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_centroid[i] = -1;
        time[i] = 0;
    }

    /* ---- LSH model ---- */
    int w = 4*compute_window(lsh_dataset);
    LSH <int>* model = new LSH <int> (k, L, w);
    model->fit(lsh_dataset);
    model->evaluate(lsh_searchset, R, &R_neighbors, &min_distance, &time, &nearest_centroid);
    delete (model);

    int counter_for_queries = 0;                            // skip ids of centroids
    for (int i = 0; i < lsh_searchset->size(); i++){
        if (find(centroids.begin(), centroids.end(), &(*lsh_dataset)[i]) != centroids.end()) {
            counter_for_queries++;
            continue;
        }
        clusters[nearest_centroid[i]]->push_back(counter_for_queries);
        counter_for_queries++;
    }

}

void lsh_datatype(vector<vector<double*>>* lsh_dataset, vector<vector<double*>>* lsh_searchset, int Grids, int k, int L, vector<vector<double*>*> centroids, vector<int>** clusters){

    double delta = 0.00006;
    int d = 2;                                                      /* default 2D curves */
    vector<vector<double>> data_vectored_curves;
    vector<vector<double>> search_vectored_curves;
    vector<int*> hashed_neighbors;

    for (int i = 0; i < Grids; i++) {

        /* Vectorization */
        Grid_Vectorization(delta, d, lsh_dataset, lsh_searchset, &data_vectored_curves, &search_vectored_curves);

        double R;
        vector<vector<int>> R_neighbors;

        /* Arrays for results */
        double *min_distance = new double[lsh_searchset->size()];
        int *nearest_centroid = new int[lsh_searchset->size()];
        double *time = new double[lsh_searchset->size()];

        /* Initialize arrays */
        for (int i = 0; i < lsh_searchset->size(); i++) {
            min_distance[i] = -1;
            nearest_centroid[i] = -1;
            time[i] = 0;
        }

        /* ---- LSH model ---- */
        double w = 4 * compute_window(&data_vectored_curves);
        LSH<double> *model = new LSH<double>(k, L, w);
        model->fit(&data_vectored_curves);
        model->evaluate(&search_vectored_curves, R, &R_neighbors, &min_distance, &time, &nearest_centroid);
        delete (model);

        /* Store Results for all Grids */
        hashed_neighbors.push_back(nearest_centroid);

        /* Clean Vectors */
        vector<vector<double>>().swap(data_vectored_curves);
        vector<vector<double>>().swap(search_vectored_curves);
        vector<vector<int>>().swap(R_neighbors);
        /* Clean Pointers */
        delete[] min_distance;
        delete[] time;

    }

    /* Find Nearest Centroid from results of LSH */
    int counter_for_queries = 0;                            // skip ids of centroids
    double distance = 0.0;
    double min_distance = -1;
    int nearest_centroid = -1;
    for (int i = 0; i < lsh_searchset->size(); i++) {
        if (find(centroids.begin(), centroids.end(), &(*lsh_dataset)[i]) != centroids.end()) {
            counter_for_queries++;
            continue;
        }
        for (int j = 0; j < Grids; j++) {
            if (hashed_neighbors[j][i] == -1) continue;
            distance = DTW(&(*lsh_searchset)[i], &(*lsh_dataset)[hashed_neighbors[j][i]]);
            if (j == 0) {
                min_distance = distance;
                nearest_centroid = hashed_neighbors[j][i];
            } else if (distance < min_distance) {
                min_distance = distance;
                nearest_centroid = hashed_neighbors[j][i];
            }
        }
        clusters[i]->push_back(counter_for_queries);
        counter_for_queries++;
    }
}