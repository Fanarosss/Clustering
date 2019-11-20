#include "Library.h"
#include "LSH_Functions.h"
#include "Helper_Functions.h"

using namespace std;

template double compute_window<int>(vector<vector<int>>*);
template double compute_window<double>(vector<vector<double>>*);
template void projections<int>(vector<vector<int>>*, vector<vector<int>>*, vector<double>*, double, int);
template void projections<double>(vector<vector<int>>*, vector<vector<double>>*, vector<double>*, double, int);

template <typename Point>
double compute_window(vector<vector<Point>>* dataset) {
    /* 1. Take all points in dataset
     * 2. Find their nearest neighbor using L1 metric
     * 3. Average all distances between points
     * L1 = sum(|P1i - P2i|) for i in 0,d-1
     * Note: dist is a function scalable for all metrics(L1,L2,L3 etc.) */

    vector<Point> distances;
    Point L1, min_distance;
    int size = dataset->size();
    int d = dataset->at(0).size();
    vector<Point> P1;
    vector<Point> P2;

    for (int i = 0; i < size; i++) {
        P1 = dataset->at(i);
        L1 = 0;
        min_distance = -1;
        for (int j = 0; j < size; j++) {
            if (i != j) {
                P2 = dataset->at(j);
                /* Default is L1 metric, for Lk metric, add a 4th argument*/
                L1 = dist(&P1, &P2, d);
                if (min_distance == -1)
                    min_distance = L1;
                if (L1 < min_distance)
                    min_distance = L1;
            }
            L1 = 0;
        }
        distances.push_back(min_distance);
    }
    double w = accumulate(distances.begin(), distances.end(), 0) / size;

    return w;
}

void generate_shifts(vector<vector<vector<double>>>* s, double w, int d, int k, int L){
    /* Generate K * Si for every dimension
     * At the end, s will be a vector of size (k,d) */

    unsigned seed;
    uniform_real_distribution<double> distribution (0, w);
    vector<double> Sj;
    vector<vector<double>> Sl;

    for (int l = 0; l < L; l++) {
        for (int i = 0; i < k; i++) {
            seed = chrono::system_clock::now().time_since_epoch().count();
            default_random_engine generator(seed);
            for (int j = 1; j < d; j++) {
                Sj.push_back(distribution(generator));
            }
            Sl.push_back(Sj);
            vector<double>().swap(Sj);
        }
        s->push_back(Sl);
        vector<vector<double>>().swap(Sl);
    }
}

template <typename Point>
void projections(vector<vector<int>>* a_projects, vector<vector<Point>>* x, vector<double>* s, double w, int d) {
    /* Ai = (Xi - Si) / W
     * Project every X to A in d-dimensional grid shifted by S, where every cell size = W */

    int ai;
    vector<int> a;

    for (int i = 0; i < x->size(); i++) {
        for (int dim = 1; dim < d; dim++) {
            ai = floor((double)((*x)[i][dim] - (*s)[dim - 1]) / w); // used to be + w
            a.push_back(ai);
        }
        a_projects->push_back(a);
        a.clear();
        a.shrink_to_fit();
    }
}

void compute_hash(vector<int>* H, vector<vector<int>> *a, int** power, int d, int k, double w){
    /* Compute K of hash functions hi for every point
     * Vector H at the end will have a size of (dataset.size(), k)
     * Usage of moduloMultiplication and modulo functions was necessary to avoid overflow */

    int M = 0, h = 0, term = 0, dim = d - 1;
    M = pow(2, 32/k);

    for (int i = 0; i < a->size(); i++){
        h=0;
        term = 0;
        for (int j = dim - 1; j >= 0; j--) {
            term += moduloMultiplication((*a)[i][j], (*power)[(dim - 1) - j], M);
        }
        h = modulo(term, M);
        H->push_back(h);
    }
}

void amplify_hash(vector<int>* amplified_g, vector<vector<int>>* hash_functions, int k){
    /* For every item it amplifies the hash from K dimensions to 1
     * g(x) = [h1(x)|h2(x)|h3(x)....|hk(x)] */

    int g;
    int concat_dist = 31/k;

    /* for all points in dataset */
    for (int i = 0; i < (*hash_functions)[0].size(); i++) {
        g=0;
        for (int j = 0; j < k; j++) {
            g += (g << concat_dist) | ((*hash_functions)[j][i]);
        }
        amplified_g->push_back(g);
    }
}


/* Grid Vectorization Functions*/

void Grid_Vectorization(double delta, int d, vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, vector<vector<double>>* data_vectored_curves, vector<vector<double>>* search_vectored_curves) {

    /* variable decl */
    double max_element = 0.0;
    int max_points = 0, elements = 0;
    /* orthogonal grid of size d */
    vector<double> orthogonal_grid;
    /* vector for hashed curves */
    vector<vector<vector<double>>> hashed_curves;
    vector<vector<double>> temp_hash;
    /* temp */
    vector<double> curve;

    /* ----------------------- HASHING with ORTHOGONAL GRID ---------------------- */
    shift_grid(&orthogonal_grid, delta, d);

    /* ------------------ DATA SET hashing ----------------- */
    /* hash all dataset curves */
    for (int i = 0; i < dataset->size(); i++) {
        hash_curve(&temp_hash, &(*dataset)[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        /* clean temp hash */
        vector<vector<double>>().swap(temp_hash);
    }

    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */
    elements = 0;
    max_points = 0;
    max_element = 0.0;
    /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
    for (int i = 0; i < hashed_curves.size(); i++) {
        for (int j = 0; j < hashed_curves[i].size(); j++) {
            curve.push_back(hashed_curves[i][j][0]);
            /* find max element */
            /* j != 0 because at index 0 there is the id and length of curve */
            if (j != 0) {
                curve.push_back(hashed_curves[i][j][1]);
                if (hashed_curves[i][j][0] > max_element) {
                    max_element = hashed_curves[i][j][0];
                }
                if (hashed_curves[i][j][1] > max_element) {
                    max_element = hashed_curves[i][j][1];
                }
            }
            elements++;
        }

        if (elements > max_points) {
            max_points = elements;
        }
        elements = 0;
        data_vectored_curves->push_back(curve);
        vector<double>().swap(curve);
    }

    /* ------------------ SEARCH SET hashing ----------------- */
    /* clean hashed curves */
    vector<vector<vector<double>>>().swap(hashed_curves);
    /* end of cleaning */

    /* hash all curves */
    for (int i = 0; i < searchset->size(); i++) {
        hash_curve(&temp_hash, &(*searchset)[i], &orthogonal_grid, delta, d);
        hashed_curves.push_back(temp_hash);
        /* clean temp hash */
        vector<vector<double>>().swap(temp_hash);
    }

    /* now that we have each hash, we can find by adding the orthogonal grid to the hash points
     * the equivalent points in our new grid, that the polygonal curve projects */

    /* concat the 2d points in every h to make it from (x1,y1)(x2,y2) to x1,y1,x2,y2 */
    elements = 0;
    for (int i = 0; i < hashed_curves.size(); i++) {
        for (int j = 0; j < hashed_curves[i].size(); j++) {
            /* we always push back the first element - even the id on (0,0)
             * we won't push back the dimensions of a curve in (0,1) */
            curve.push_back(hashed_curves[i][j][0]);
            /* j != 0 because at index 0 there is the id and length of curve */
            if (j != 0) {
                /* we dont push back the length size for the lsh */
                curve.push_back(hashed_curves[i][j][1]);
                if (hashed_curves[i][j][0] > max_element) {
                    max_element = hashed_curves[i][j][0];
                }
                if (hashed_curves[i][j][1] > max_element) {
                    max_element = hashed_curves[i][j][1];
                }
            }
            elements++;
        }
        if (elements > max_points) {
            max_points = elements;
        }
        elements = 0;
        search_vectored_curves->push_back(curve);
        vector<double>().swap(curve);
    }

    /* clear no longer used vectors for memory optimizations */
    vector<double>().swap(orthogonal_grid);
    /* clean hashed curves */
    vector<vector<vector<double>>>().swap(hashed_curves);
    /* end of cleaning */

    /* ----------------- PADDING for both sets ------------ */
    /* pad special number > max coord */
    for (int i = 0; i < data_vectored_curves->size(); i++) {
        while ((*data_vectored_curves)[i].size() < max_points) {
            (*data_vectored_curves)[i].push_back(2 * max_element);
        }
    }
    for (int i = 0; i < search_vectored_curves->size(); i++) {
        while ((*search_vectored_curves)[i].size() < max_points) {
            (*search_vectored_curves)[i].push_back(2 * max_element);
        }
    }
    /* ----------- end of padding ------------ */
    return;
}


void shift_grid(vector<double>* orthogonal_grid, int delta, int d) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);

    /* t uniformly in [0,d) */
    uniform_real_distribution<double> distribution (0, d);
    for (int i = 0; i < d; i++) {
        orthogonal_grid->push_back(distribution(generator));
    }
}

void hash_curve(vector<vector<double>>* hashed_curve, vector<double*>* curve, vector<double>* orthogonal_grid, double delta, int d) {
    double* pi;
    vector<double> pi_new;
    vector<double> pi_old;
    /* first index has the id and the dimensions of the curve */
    pi = (*curve)[0];
    pi_new.push_back(pi[0]);
    pi_new.push_back(pi[1]);
    hashed_curve->push_back(pi_new);
    pi_old = pi_new;
    for (int i = 1; i < (*curve)[0][1]; i++) {
        pi = (*curve)[i];
        pi_new = arg_min(&pi, orthogonal_grid, delta, d);
    /* remove consecutive duplicates pi' from the hashed_curve */
    if ((pi_new[0] != pi_old[0]) || (pi_new[1] != pi_old[1])) {
        hashed_curve->push_back(pi_new);
    }
    pi_old = pi_new;
    }
}

/* function for snapping a point onto the grid */
vector<double> arg_min(double** pi, vector<double>* orthogonal_grid, double delta, int d) {
    double min, num, shift;
    int q;
    vector<double> argmin;
    /* Point is to minimize the ||pi-q|| for all q */
    for (int i = 0; i < d; i++) {
        num = (*pi)[i];
        shift = (*orthogonal_grid)[i];
        /* snapping */
        argmin.push_back(num + (delta + shift)/2);
        argmin[i] -= fmod(argmin[i], delta);
    }
    return argmin;
}
