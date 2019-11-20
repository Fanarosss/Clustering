#include "LSH.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"
#include "LSH_DataTypes.h"

using namespace std;

void lsh_datatype(vector<vector<int>>* lsh_dataset, vector<vector<int>>* lsh_searchset, int k, int L, vector<int>* centroids, vector<int>** clusters){

    int data_size = lsh_dataset->size() + lsh_searchset->size();

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

    int centroids_found = 0;
    for (int i = 0; i < data_size; i++){
        if (find(centroids->begin(), centroids->end(), i) != centroids->end()) {
            centroids_found++;
            continue;
        }
        clusters[nearest_centroid[i-centroids_found]]->push_back(i);
    }

}

void lsh_datatype(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, int k, int L, vector<int>* centroids, vector<int>** clusters){

}