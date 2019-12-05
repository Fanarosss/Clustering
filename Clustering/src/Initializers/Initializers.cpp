#include "Initializers.h"
#include "Library.h"
#include "Helper_Functions.h"
#include <algorithm>
#include <bits/stdc++.h>

using namespace std;

template <class Point>
vector<pair<vector<Point>*, int>> Random_Selection<Point>::init(vector<vector<Point>>* dataset) {


    int id;
    unsigned seed;
    vector<int> ids;
    vector<pair<vector<Point>*, int>> centroids;

    uniform_int_distribution<int> distribution (0, dataset->size()-1);
    seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    for (int i = 0; i < this->get_K(); i++) {
        id = distribution(generator);
        while (find(ids.begin(), ids.end(), id) != ids.end()) {
            id = distribution(generator);
        }
        centroids.push_back(make_pair(&(*dataset)[id], id));
        ids.push_back(id);
    }
    return centroids;
}

template <class Point>
string Random_Selection<Point>::get_name() {
    return this->name;
}

template <class Point>
vector<pair<vector<Point>*, int>> KMeans_plusplus<Point>::init(vector<vector<Point>>* dataset) {

    int r, id;
    /* # of centroids */
    int t = 1;
    /* # of data */
    int n = dataset->size();
    /* point, ids of points in dataset */
    vector<pair<vector<Point>*, int>> centroids;
    /* ids of points*/
    vector<int> centroids_id;
    /* array with min distances of points from centroids */
    vector<double> D;
    /* vector for partial sums */
    vector< pair <int,double> > P;
    /* ids map */
    vector<int> id_map;
    /* uniform distribution */
    unsigned seed;
    uniform_int_distribution<int> distribution (0, dataset->size()-1);
    seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    /* Step 1: Choose a centroid uniformly at random; t←1 */
    id = distribution(generator);
    centroids.push_back(make_pair(&(*dataset)[id], id));
    centroids_id.push_back(id);
    while (t < this->get_K()) {
        /* Step 2: for all non-centroid point i=1,...,n−t, letD(i)←min distance to some centroid,
         * among t chosen centroids. */
        for (int i = 0; i < n; i++) {
            // for all non centroids
            if (find(centroids_id.begin(), centroids_id.end(), i) != centroids_id.end()) continue;
            D.push_back(min_distance(i, &centroids_id, dataset));
            id_map.push_back(i);
        }

        /* Normalization */
        normalize(&D);

        /* Step 3: Choose new centroid: r chosen with probability proportional to D(r)^2 */
        double value, temp;
        r = 0;
        for (int i = 0; i < n; i++) {
            if (find(centroids_id.begin(), centroids_id.end(), i) != centroids_id.end()) continue;
            value = Sum(0, r, &D, 2);
            /* keep pair of value and id -> we want to know after the sorting
             * the id that corresponds to each value */
            P.push_back(make_pair(i, value));
            r++;
        }

        /* sort P */
        sort(P.begin(), P.end(), sortbysec);

        /* uniform distribution */
        unsigned seed;
        uniform_real_distribution<double> distribution(0, P[n - t - 1].second);
        seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        double x = distribution(generator);

        r = 0;
        for (r = 0; r < P.size(); r++) {
            if (x <= P[r].second) break;
        }
        t++;

        /* careful, r is the index of the sorted data, i need to keep index of the initial ones */
        centroids.push_back(make_pair(&(*dataset)[P[r].first], P[r].first));
        centroids_id.push_back(P[r].first);
        /* Step 4: Go to (2) until t = k = given #centroids. */
        vector<int>().swap(id_map);
        vector<double>().swap(D);
        vector<pair<int,double>>().swap(P);
    }
    return centroids;
}

template <class Point>
string KMeans_plusplus<Point>::get_name() {
    return this->name;
}

template class KMeans_plusplus<double>;
template class KMeans_plusplus<double*>;
template class Random_Selection<double>;
template class Random_Selection<double*>;

// Driver function to sort the vector elements
// by second element of pairs
bool sortbysec(const pair<int, double> &a, const pair<int, double> &b) {
    return (a.second < b.second);
}