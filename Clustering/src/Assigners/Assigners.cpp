#include "Assigners.h"
#include "Library.h"
#include "Helper_Functions.h"
#include "LSH.h"
#include "LSH_Functions.h"


using namespace std;

template <class Point>
vector<int>** Lloyd_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {
    cout << '\t' << "Assigning with Lloyd's Assignment" << endl;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];
    for(int i = 0 ; i < num_of_centroids ; i++ )
        clusters[i] = new vector<int>;
    int centroid;
    double min_dist, curr_dist;
    for (int i = 0; i < data_size; i++) {
        if (find(centroids->begin(), centroids->end(), i) != centroids->end()) continue;
        centroid = -1;
        min_dist = DBL_MAX;
        for (int j = 0; j < num_of_centroids; j++) {
            curr_dist = dist(&(*dataset)[(*centroids)[j]], &(*dataset)[i], dimension);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                centroid = j;
            }
        }
        clusters[centroid]->push_back(i);
    }
    return clusters;
}

template <class Point>
string Lloyd_assignment<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<int>** Inverse_assignment<Point>::assign(vector<vector<Point>>* dataset, vector<int>* centroids) {
    cout << '\t' << "Assigning with Inverse Assignment" << endl;
    int num_of_centroids = centroids->size();
    int data_size = dataset->size();
    int dimension = (*dataset)[0].size();
    vector<int>** clusters;
    clusters = new vector<int>*[num_of_centroids];

    vector<vector<Point>> lsh_dataset;                 //centroids
    vector<vector<Point>> lsh_searchset;               //queries
    for (int i = 0; i < num_of_centroids; i++){
        lsh_dataset.push_back((*dataset)[(*centroids)[i]]);
    }
    for (int j = 0; j < data_size; j++){
        if (find(centroids->begin(), centroids->end(), j) != centroids->end()) continue;
        lsh_searchset.push_back((*dataset)[(*centroids)[j]]);
    }

    double R;
    vector<vector<int>> R_neighbors;

    /* Arrays for results */
    int *min_distance = new int[lsh_searchset.size()];
    int *nearest_centroid = new int[lsh_searchset.size()];
    double *time = new double[lsh_searchset.size()];

    /* Initialize arrays */
    for (int i = 0; i < lsh_searchset.size(); i++) {
        min_distance[i] = INT_MAX;
        nearest_centroid[i] = -1;
        time[i] = 0;
    }

    /* ---- LSH model ---- */
    Point w = 4*compute_window(&lsh_dataset);
    LSH <Point>* model = new LSH <Point> (this->k, this->L, w);
    model->fit(&lsh_dataset);
    model->evaluate(&lsh_searchset, R, &R_neighbors, &min_distance, &time, &nearest_centroid);
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

template <class Point>
string Inverse_assignment<Point>::get_name() {
    return this->name;
}

template class Lloyd_assignment<int>;
template class Lloyd_assignment<double*>;
template class Inverse_assignment<int>;
template class Inverse_assignment<double*>;