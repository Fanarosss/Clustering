#include "Library.h"

using namespace std;

template <typename Point> double min_distance(int, vector<int>*, vector<vector<Point>>*);

void show_cluster_usage(string name);
int Read_input_file(string input);
int Read_files(vector<vector<int>>*, int*, string, string );
int Read_files(vector<vector<double*>>*, int*, string, string);
double point_dist(double*, double*, int=2);
double DTW(vector<double*>*, vector<double*>*);
double dist(vector<int>* P1, vector<int>* P2, int=1);
double dist(vector<double*>* P1, vector<double*>* P2, int=1);
double min(double, double, double);
void normalize(vector<double>*);
double Sum(int, int, vector<double>*, int);
