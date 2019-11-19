#include "Helper_Functions.h"

using namespace std;

template double min_distance<int>(int, vector<int>*, vector<vector<int>>*);
template double min_distance<double*>(int, vector<int>*, vector<vector<double*>>*);

void show_cluster_usage(string name)
{
    cerr      << "Usage:   " << name << " -letter(s) <option(s)>\n"
              << "Options:\n"
              << "\t-i <input file>             Path to input.dat\n"
              << "\t-c <configuration file>     Path to cluster.conf\n"
              << "\t-o <output file>            Path to file of results\n"
              << "\t-complete <optional>        Thorough Results\n"
              << endl;
}

int Read_input_file(string input){
    string line;
    ifstream input_file(input);
    getline(input_file, line);
    if (!line.empty() && line[line.size() - 1] == '\r')
        line.erase(line.size() - 1);
    if (!line.empty() && line[line.size() - 1] == '\n')
        line.erase(line.size() - 1);
    if (line == "vectors"){
        return 1;
    }else if (line == "curves"){
        return 2;
    }else{
        return 0;
    }
}

int Read_files(vector<vector<int>>* cluster_data, int* cluster_config, string input_file_name, string config_file_name) {

    string line;
    int id;
    int number;
    vector<int> v;

    ifstream input_file(input_file_name);
    ifstream config_file(config_file_name);

    id = 0;
    /* discard "vectors" */
    getline(input_file, line);
    while (getline(input_file, line)) {
        stringstream ss(line);
        /* discard id */
        ss >> number;
        v.push_back(id);
        while (ss >> number) {
            v.push_back(number);
        }
        cluster_data->push_back(v);
        v.clear();
        id++;
    }

    id = 0;
    while (getline(config_file, line)) {
        stringstream ss(line);
        while (line != ":"){
            ss >> line;
        }
        ss >> number;
        cluster_config[id] = number;
        id++;
    }

    if (cluster_data->size() == 0) {
        cerr << "Error: input.dat file is empty!" << endl;
        return -1;
    }

    return 1;
}

int Read_files(vector<vector<double*>>* cluster_data, int* cluster_config, string input_file_name, string config_file_name) {

    string line;
    int id, trash;
    char bracket, comma;
    double number;
    double* point;
    vector<double*> v;

    ifstream input_file(input_file_name);
    ifstream config_file(config_file_name);

    id = 0;
    /* discard "curves" */
    getline(input_file, line);
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
        cluster_data->push_back(v);
        v.clear();
        id++;
    }

    id = 0;
    while (getline(config_file, line)) {
        stringstream ss(line);
        while (line != ":"){
            ss >> line;
        }
        ss >> number;
        cluster_config[id] = number;
        id++;
    }

    if (cluster_data->size() == 0) {
        cerr << "Error: input.dat file is empty!" << endl;
        return -1;
    }
    return 1;
}

/* distance of vectors-curves */
double dist(vector<int>* P1, vector<int>* P2, int Metric) {
    /* Lk metric
     * for metric = 1 we have L1 metric
     * for metric = 2 we have L2 metric etc.
     * (default value = L1 Metric) -> Manhattan distance */
    double dist = 0;
    for (int dim = 1; dim < P1->size(); dim++)
        dist += pow(fabs((*P1)[dim] - (*P2)[dim]),Metric);
    return pow(dist,1/(double)Metric);
}

/* distance of vectors-curves */
double dist(vector<double*>* P1, vector<double*>* P2, int Metric) {
    return DTW(P1, P2);
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

void normalize(vector<double>* D) {
    auto it = max_element(D->begin(), D->end());
    double max_D = D->at(distance(D->begin(), it));
    for (int i = 0; i < D->size(); i++) {
        (*D)[i] = (*D)[i] / max_D;
    }
    return;
}

double Sum(int start, int end, vector<double>* D, int power) {
    double sum = 0.0;
    for (int i = start; i <= end; i++) {
        sum += pow((*D)[i], power);
    }
    return sum;
}

template <typename Point>
double min_distance(int index, vector<int>* centroids, vector<vector<Point>>* dataset) {
    double distance = 0.0;
    double min_distance = DBL_MAX;
    for (int i = 0; i < centroids->size(); i++) {
        distance = dist(&dataset->at(centroids->at(i)), &dataset->at(index));
        if (distance < min_distance)
            min_distance = distance;
    }
    return min_distance;
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

double min(double x, double y, double z) {
    /* get the min out of 3 real numbers */
    double temp = (x < y) ? x : y;
    return (z < temp) ? z : temp;
}