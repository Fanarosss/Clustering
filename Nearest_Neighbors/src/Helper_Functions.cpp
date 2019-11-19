#include "Helper_Functions.h"
#include "HashTable.h"
#include "Traversals.h"

using namespace std;

template double dist<int>(vector<int>*, vector<int>*, int, int=1);
template double dist<double>(vector<double>*, vector<double>*, int, int=1);

int Read_point_files(vector<vector<int>>* dataset, vector<vector<int>>* searchset, double* R, string data_filename, string query_filename) {
    /* Create dataset and searchset vectors from files */

    string line;
    int id;
    int number;
    vector<int> v;

    ifstream input_file(data_filename);
    ifstream query_file(query_filename);

    id = 0;
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* discard id */
        ss >> number;
        v.push_back(id);
        while (ss >> number) {
            v.push_back(number);
        }
        dataset->push_back(v);
        v.clear();
        id++;
    }

    id = -1;
    while (getline(query_file, line)) {
        stringstream ss(line);
        if(id == -1){
            /* Radius: */
            ss >> line;
            ss >> (*R);
            id++;
        }else {
            /* discard id */
            ss >> number;
            v.push_back(id);
            while (ss >> number) {
                v.push_back(number);
            }
            searchset->push_back(v);
            v.clear();
            id++;
        }
    }

    if (dataset->size() == 0) {
        cerr << "Error: dataset file is empty!" << endl;
        return -1;
    }
    /* <= 1 because first line is radius */
    if (searchset->size() <= 1) {
        cerr << "Error: searchset file is empty!" << endl;
        return -1;
    }

    return 1;
}

int Read_curve_files(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, string data_filename, string query_filename) {
    string line;
    int id, trash;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

    ifstream input_file(data_filename);
    ifstream query_file(query_filename);

    id = 0;
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
        /* length */
        ss >> point[1];
        v.push_back(point);
        /* for all data points */
        while (ss >> bracket) {
            point = new double [2];
            ss >> point[0];
            ss >> comma;
            ss >> point[1];
            ss >> bracket;
            v.push_back(point);
        }
        dataset->push_back(v);
        v.clear();
        id++;
    }

    id = 0;
    while (getline(query_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
        /* length */
        ss >> point[1];
        v.push_back(point);
        /* for all data points */
        while (ss >> bracket) {
            point = new double [2];
            ss >> point[0];
            ss >> comma;
            ss >> point[1];
            ss >> bracket;
            v.push_back(point);
        }
        searchset->push_back(v);
        v.clear();
        id++;
    }

    if (dataset->size() == 0) {
        cerr << "Error: dataset file is empty!" << endl;
        return -1;
    }
    if (searchset->size() == 0) {
        cerr << "Error: searchset file is empty!" << endl;
        return -1;
    }
    return 1;
}

int Read_curve_files_max_dim(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, string data_filename, string query_filename, double max_dim) {
    string line;
    int id, trash;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

    ifstream input_file(data_filename);
    ifstream query_file(query_filename);

    id = 0;
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
        /* length */
        ss >> point[1];
        if (point[1] <= max_dim) {
            v.push_back(point);
            /* for all data points */
            while (ss >> bracket) {
                point = new double[2];
                ss >> point[0];
                ss >> comma;
                ss >> point[1];
                ss >> bracket;
                v.push_back(point);
            }
            dataset->push_back(v);
            v.clear();
            id++;
        }
    }

    id = 0;
    while (getline(query_file, line)) {
        stringstream ss(line);
        /* id */
        ss >> trash;
        point = new double [2];
        point[0] = id;
        /* length */
        ss >> point[1];
        if (point[1] <= max_dim) {
            v.push_back(point);
            /* for all data points */
            while (ss >> bracket) {
                point = new double[2];
                ss >> point[0];
                ss >> comma;
                ss >> point[1];
                ss >> bracket;
                v.push_back(point);
            }
            searchset->push_back(v);
            v.clear();
            id++;
        }
    }

    if (dataset->size() == 0) {
        cerr << "Error: dataset file is empty!" << endl;
        return -1;
    }
    if (searchset->size() == 0) {
        cerr << "Error: searchset file is empty!" << endl;
        return -1;
    }
    return 1;
}

template <typename Point>
double dist(vector<Point>* P1, vector<Point>* P2, int d, int Metric) {
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    double dist = 0;
    for (int dim = 1; dim < d; dim++)
        dist += pow(fabs((*P1)[dim] - (*P2)[dim]),Metric);
    return pow(dist,1/(double)Metric);
}

double point_dist(double* p, double* q, int Metric){
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    double d1, d2;
    d1 = pow(fabs(p[0] - q[0]), Metric);
    d2 = pow(fabs(p[1] - q[1]), Metric);
    return pow(d1+d2,1/(double)Metric);
}

double DTW(vector<double*>* P, vector<double*>* Q) {
    /* Initialize c(1,1) = ||p1-q1||
    * * if j > 1, then c(1,j) = c(1,j-1) + ||pi - qj||
   * if i > 1, then c(i,1) = c(i-1,1) + ||pi - qj||
   * if i > 1, j > 1, then c(i,j) = min{c(i-1,j), c(i-1,j-1), c(i,j-1)} + ||pi - qj|| */
    int m1 = P->size() - 1;
    int m2 = Q->size() - 1;

    /* allocate space */
    double ** c = new double* [m1];
    for (int i = 0; i < m1; i++) {
        c[i] = new double [m2];
    }

    c[0][0] = point_dist((*P)[1], (*Q)[1], 2);
    for (int i = 0; i < m1; i++) {
        for (int j = 0; j < m2; j++) {
            if (i == 0 && j == 0) continue;
            if (j > 0 && i == 0) {
                c[i][j] = c[i][j - 1] + point_dist((*P)[i+1], (*Q)[j+1], 2);
            } else if (i > 0 && j == 0) {
                c[i][j] = c[i - 1][j] + point_dist((*P)[i+1], (*Q)[j+1], 2);
            } else {
                c[i][j] = min(c[i - 1][j], c[i - 1][j - 1], c[i][j - 1]) + point_dist((*P)[i+1], (*Q)[j+1], 2);
            }
        }
    }
    double res = c[m1-1][m2-1];

    /* Free allocated space */
    for (int i = 0; i < m1; i++) {
        delete[] c[i];
    }
    delete[] c;

    return res;
}

int modulo (int a, int b){
    int m = a % b;
    if (m < 0){
        m = (b < 0) ? m - b : m + b;
    }
    return m;
}

// Returns (a * b) % mod
int moduloMultiplication(int a, int b, int mod) {
    int res = 0;
    a %= mod;
    while (b)
    {
        // If b is odd, add a with result
        if (b & 1)
            res = (res + a) % mod;
        // Here we assume that doing 2*a
        // doesn't cause overflow
        a = (2 * a) % mod;
        b >>= 1; // b = b / 2
    }
    return res;
}

long moduloPower(long base,long exp,long div) {
    if (exp == 0)   return 1;
    if (exp == 1)   return base % div;
    if (exp % 2 == 0)   return (moduloPower(base, exp / 2, div) * moduloPower(base, exp / 2, div)) % div;
    return (moduloPower(base, exp - 1, div) * moduloPower(base, 1, div)) % div;
}

void brute_force(vector<vector<int>>* dataset, vector<vector<int>>* searchset, vector<int>* TrueDistances, vector<double>* TrueTimes) {
    /* vectors init */
    vector<int> n_neighbors;
    vector<int> distances;
    vector<int> P1;
    vector<int> P2;
    /* variable declaration */
    int L1, min_distance, n_neighbor;
    int d_size = dataset->size();
    int s_size = searchset->size();
    int d = (*dataset)[0].size();
    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/nneighbors_brute_force.txt");
    for (int i = 0; i < s_size; i++) {
        P1 = (*searchset)[i];
        L1 = 0;
        min_distance = -1;
        auto start = chrono::high_resolution_clock::now();
        for (int j = 1; j < d_size; j++) {
            P2 = (*dataset)[j];
            /* default is L1 metric, for Lk metric, add a 4th argument, k */
            L1 = dist(&P1, &P2, d);
            if (min_distance == -1) {
                min_distance = L1;
                n_neighbor = j;
            }
            if (L1 < min_distance) {
                min_distance = L1;
                n_neighbor = j;
            }
        }
        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();
        TrueDistances->push_back(min_distance);
        TrueTimes->push_back(time_elapsed);
        distances.push_back(min_distance);
        n_neighbors.push_back(n_neighbor);
        neighbors_file << "Item:" << setw(floor(log10(s_size) + 1)) << setfill('0') << i + 1 << ", Neighbor: " << setw(floor(log10(d_size) + 1)) << setfill('0') << n_neighbor + 1 << " | Distance: " << setw(4) << setfill('0') << min_distance <<  " | Duration: " << time_elapsed << endl;
    }
    neighbors_file.close();
}

void read_vectors_brute_force_file(string bffile, vector<int>* TrueDistances, vector<double>* TrueTimes){
    ifstream brute_force_file(bffile);
    string line, trash;
    double tdistance, ttime;
    while (getline(brute_force_file, line)) {
        stringstream ss(line);
        /* Item */
        ss >> trash;
        /* Neighbor: */
        ss >> trash;
        /* True Neighbor */
        ss >> trash;
        /* |*/
        ss >> trash;
        /* Distance: */
        ss >> trash;
        /* True Distance */
        ss >> tdistance;
        TrueDistances->push_back(tdistance);
        /* |*/
        ss >> trash;
        /* Duration: */
        ss >> trash;
        /* True Time */
        ss >> ttime;
        TrueTimes->push_back(ttime);
    }
}

void curves_brute_force(vector<vector<double*>>* dataset, vector<vector<double*>>* searchset, vector<double>* TrueDistances, vector<double>* TrueTimes, vector<int>* TrueNeighbors) {
    /* vectors init */
    vector<double*> P1;
    vector<double*> P2;
    /* variable declaration */
    double L1 = 0.0, min_distance = -1.0;
    int n_neighbor = -1;
    int d_size = dataset->size();
    int s_size = searchset->size();
    int d = (*dataset)[0].size();

    /* open file to write results */
    ofstream neighbors_file;
    neighbors_file.open ("./output/curves_brute_force.txt");
    for (int i = 0; i < s_size; i++) {
        P1 = (*searchset)[i];
        min_distance = -1.0;
        auto start = chrono::high_resolution_clock::now();
        for (int j = 0; j < d_size; j++) {
            P2 = (*dataset)[j];
            /* DTW metric for curves */
            L1 = DTW(&P1, &P2);
            if (min_distance == -1) {
                min_distance = L1;
                n_neighbor = j;
            }
            if (L1 < min_distance) {
                min_distance = L1;
                n_neighbor = j;
            }
        }

        auto finish = chrono::high_resolution_clock::now();
        auto elapsed = finish - start;
        double time_elapsed = chrono::duration<double>(elapsed).count();

        TrueDistances->push_back(min_distance);
        TrueTimes->push_back(time_elapsed);
        TrueNeighbors->push_back(n_neighbor);

        neighbors_file << "Item:" << setw(floor(log10(s_size) + 1)) << setfill('0') << i + 1 << ", Neighbor: " << setw(floor(log10(d_size) + 1)) << setfill('0') << n_neighbor + 1 << " | Distance: " << setw(10) << setprecision(9) << setfill('0') << min_distance <<  " | Duration: " << time_elapsed << endl;
    }
    neighbors_file.close();
}

void read_curves_brute_force_file(string bffile, vector<double>* TrueDistances, vector<double>* TrueTimes, vector<int>* TrueNeighbors){
    ifstream brute_force_file(bffile);
    string line, trash;
    double tdistance, ttime;
    int tneighbor;
    while (getline(brute_force_file, line)) {
        stringstream ss(line);
        /* Item */
        ss >> trash;
        /* Neighbor: */
        ss >> trash;
        /* True Neighbor */
        ss >> tneighbor;
        TrueNeighbors->push_back(tneighbor);
        /* |*/
        ss >> trash;
        /* Distance: */
        ss >> trash;
        /* True Distance */
        ss >> tdistance;
        TrueDistances->push_back(tdistance);
        /* |*/
        ss >> trash;
        /* Duration: */
        ss >> trash;
        /* True Time */
        ss >> ttime;
        TrueTimes->push_back(ttime);
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

double min(double x, double y, double z) {
    /* get the min out of 3 real numbers */
    double temp = (x < y) ? x : y;
    return (z < temp) ? z : temp;
}
